//! Scoring primitives: affine gap costs, substitution models and precomputed
//! per-allele substitution profiles.
//!
//! **Trick — reference profiles.** The DP inner loop needs `s(q_i, r_j)` for every
//! cell. Instead of branching on match / mismatch / ambiguous per cell, each
//! allele gets a `5 × |r|` table `profile[c][j] = s(c, r_j)` built once at database
//! load. With the query as the outer DP loop, the inner loop streams one
//! contiguous row: branch-free, auto-vectorisable, and the layout an inter-sequence
//! or striped SIMD kernel wants anyway.
//!
//! Because profiles are built through [`SubstitutionModel`], position-specific
//! scoring (SHM hotspot-aware mismatch costs, softer penalties at known
//! polymorphic positions, codon-aware scores) plugs in without touching any
//! aligner, clipper or adjudicator — they all read the same profile.

use crate::alphabet::{A, ALPHABET_SIZE, C, G, N, T};

pub type Score = i32;

#[inline]
fn is_w(c: u8) -> bool {
    c == A || c == T
}
#[inline]
fn is_r(c: u8) -> bool {
    c == A || c == G
}
#[inline]
fn is_y(c: u8) -> bool {
    c == C || c == T
}
#[inline]
fn is_s(c: u8) -> bool {
    c == C || c == G
}

/// Classic AID hotspot target at `position` of `reference` (WRCY / RGYW).
#[inline]
pub fn is_hotspot_target(reference: &[u8], position: usize) -> bool {
    let n = reference.len();
    let b = reference[position];
    // WRCY — mutable C
    if b == C
        && position >= 2
        && position + 1 < n
        && is_w(reference[position - 2])
        && is_r(reference[position - 1])
        && is_y(reference[position + 1])
    {
        return true;
    }
    // RGYW — mutable G
    if b == G
        && position >= 1
        && position + 2 < n
        && is_r(reference[position - 1])
        && is_y(reference[position + 1])
        && is_w(reference[position + 2])
    {
        return true;
    }
    false
}

/// Coldspot target at `position` (SYC / GRS).
#[inline]
pub fn is_coldspot_target(reference: &[u8], position: usize) -> bool {
    let n = reference.len();
    let b = reference[position];
    // SYC — mutable C
    if b == C && position >= 2 && is_s(reference[position - 2]) && is_y(reference[position - 1]) {
        return true;
    }
    // GRS — mutable G
    if b == G && position + 2 < n && is_r(reference[position + 1]) && is_s(reference[position + 2])
    {
        return true;
    }
    false
}

/// "Minus infinity" with headroom so that adding a few penalties never overflows.
pub const NEG_INF: Score = i32::MIN / 4;

/// Affine gap costs (both negative). A gap of length `L` costs
/// `open + (L - 1) × extend`, matching the recurrence in the method note.
///
/// Requires `open ≤ extend ≤ 0`. If opening were cheaper than extending, the
/// Gotoh recurrence would price a gap run as repeated opens and "consecutive
/// gap columns extend" — the convention used by clipping and reporting —
/// would no longer describe the optimum. The DP score itself is exact for any
/// costs (see the aligner tests); the constructor enforces the precondition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GapCosts {
    pub open: Score,
    pub extend: Score,
}

impl GapCosts {
    /// Panics unless `open <= extend <= 0` (a configuration programming error).
    pub const fn new(open: Score, extend: Score) -> Self {
        assert!(
            open <= extend && extend <= 0,
            "affine gaps need open <= extend <= 0"
        );
        Self { open, extend }
    }

    pub fn cost(&self, length: u32) -> Score {
        match length {
            0 => 0,
            l => self.open + (l as Score - 1) * self.extend,
        }
    }
}

/// Score of aligning a query code against one reference position.
pub trait SubstitutionModel {
    fn score(&self, query_code: u8, reference: &[u8], position: usize) -> Score;
}

/// The `(m, μ, g_o, g_e)` tuple of the method note.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AffineScoring {
    pub matched: Score,
    pub mismatch: Score,
    pub gaps: GapCosts,
}

impl AffineScoring {
    pub const fn new(matched: Score, mismatch: Score, gap_open: Score, gap_extend: Score) -> Self {
        Self {
            matched,
            mismatch,
            gaps: GapCosts::new(gap_open, gap_extend),
        }
    }

    /// Uniform substitution score; ambiguous bases score 0.
    #[inline]
    pub fn substitution(&self, a: u8, b: u8) -> Score {
        if a >= N || b >= N {
            0
        } else if a == b {
            self.matched
        } else {
            self.mismatch
        }
    }
}

impl SubstitutionModel for AffineScoring {
    #[inline]
    fn score(&self, query_code: u8, reference: &[u8], position: usize) -> Score {
        self.substitution(query_code, reference[position])
    }
}

/// SHM hotspot-aware substitution model.
///
/// Matches and gaps follow [`AffineScoring`]. Mismatches at germline AID
/// hotspot positions (WRCY / RGYW) cost `hotspot_mismatch`; coldspot positions
/// (SYC / GRS) cost `coldspot_mismatch`; everything else uses `base.mismatch`.
/// Ambiguous bases still score 0.
///
/// Built into J profiles by default: under heavy SHM the true allele's hotspot
/// mismatches are discounted relative to decoys whose disagreements sit on
/// cold or neutral positions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HotspotScoring {
    pub base: AffineScoring,
    pub hotspot_mismatch: Score,
    pub coldspot_mismatch: Score,
}

impl HotspotScoring {
    /// Soft hotspot (−2); coldspot stays at `base.mismatch` (no extra harden).
    /// A −1 hotspot / −4 coldspot pair over-shrunk high-SHM J scores and let
    /// decoys win; this milder discount is enough to tip close gene ties.
    pub const fn for_shm(base: AffineScoring) -> Self {
        Self {
            base,
            hotspot_mismatch: -2,
            coldspot_mismatch: base.mismatch,
        }
    }

    #[inline]
    pub fn substitution(&self, query_code: u8, reference: &[u8], position: usize) -> Score {
        let r = reference[position];
        if query_code >= N || r >= N {
            0
        } else if query_code == r {
            self.base.matched
        } else if is_hotspot_target(reference, position) {
            self.hotspot_mismatch
        } else if is_coldspot_target(reference, position) {
            self.coldspot_mismatch
        } else {
            self.base.mismatch
        }
    }
}

impl SubstitutionModel for HotspotScoring {
    #[inline]
    fn score(&self, query_code: u8, reference: &[u8], position: usize) -> Score {
        self.substitution(query_code, reference, position)
    }
}

/// Match and mismatch scores of a position-independent profile.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UniformScores {
    pub matched: Score,
    pub mismatch: Score,
}

/// `ALPHABET_SIZE` rows × `len` columns of substitution scores for one allele.
#[derive(Clone, Debug)]
pub struct SubstitutionProfile {
    len: usize,
    rows: Vec<Score>,
    max: Score,
    min: Score,
    /// Reference codes (ambiguity collapsed to `N`).
    codes: Vec<u8>,
    /// `Some` when every entry is `matched` / `mismatch` / 0 (for `N`) by
    /// identity alone, as for [`AffineScoring`]; kernels may then compare
    /// codes instead of reading the profile.
    uniform: Option<UniformScores>,
}

impl SubstitutionProfile {
    pub fn build<M: SubstitutionModel + ?Sized>(model: &M, reference: &[u8]) -> Self {
        let len = reference.len();
        let mut rows = vec![0; ALPHABET_SIZE * len];
        for (code, row) in rows
            .chunks_exact_mut(len.max(1))
            .enumerate()
            .take(ALPHABET_SIZE)
        {
            let code = u8::try_from(code).expect("alphabet is tiny");
            for (j, slot) in row.iter_mut().enumerate() {
                *slot = model.score(code, reference, j);
            }
        }
        let max = rows.iter().copied().max().unwrap_or(0);
        let min = rows.iter().copied().min().unwrap_or(0);
        let codes: Vec<u8> = reference.iter().map(|&c| c.min(N)).collect();
        let uniform = uniform_scores(&rows, &codes);
        Self {
            len,
            rows,
            max,
            min,
            codes,
            uniform,
        }
    }

    /// Match/mismatch scores when the profile depends on identity alone.
    pub fn uniform(&self) -> Option<UniformScores> {
        self.uniform
    }

    /// Reference codes, ambiguity collapsed to `N`.
    pub fn codes(&self) -> &[u8] {
        &self.codes
    }

    /// Largest single-cell score; bounds any alignment score by `max × min(n, m)`.
    pub fn max_score(&self) -> Score {
        self.max
    }

    /// Most negative single-cell score.
    pub fn min_score(&self) -> Score {
        self.min
    }

    /// Reference length.
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Scores of `query_code` against every reference position.
    #[inline]
    pub fn row(&self, query_code: u8) -> &[Score] {
        let c = usize::from(query_code.min(N));
        &self.rows[c * self.len..(c + 1) * self.len]
    }

    #[inline]
    pub fn at(&self, query_code: u8, position: usize) -> Score {
        self.row(query_code)[position]
    }
}

/// The uniform scores reproducing `rows` exactly, if any: `matched` where
/// the query code equals a canonical reference code, 0 where either is `N`,
/// `mismatch` elsewhere.
fn uniform_scores(rows: &[Score], codes: &[u8]) -> Option<UniformScores> {
    let len = codes.len();
    let at = |q: usize, j: usize| rows[q * len + j];
    let canonical = codes.iter().position(|&c| c < N);
    let scores = match canonical {
        Some(j) => UniformScores {
            matched: at(usize::from(codes[j]), j),
            mismatch: at(usize::from((codes[j] + 1) % N), j),
        },
        None => UniformScores {
            matched: 0,
            mismatch: 0,
        },
    };
    let expected = |q: u8, r: u8| {
        if q >= N || r >= N {
            0
        } else if q == r {
            scores.matched
        } else {
            scores.mismatch
        }
    };
    let exact = (0..ALPHABET_SIZE as u8).all(|q| {
        codes
            .iter()
            .enumerate()
            .all(|(j, &r)| at(usize::from(q), j) == expected(q, r))
    });
    exact.then_some(scores)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;

    #[test]
    fn uniform_profiles_are_recognised_and_hotspot_profiles_are_not() {
        let scoring = AffineScoring::new(2, -3, -5, -1);
        let reference = encode(b"AGCTGTCNNA");
        let uniform = SubstitutionProfile::build(&scoring, &reference);
        assert_eq!(
            uniform.uniform(),
            Some(UniformScores {
                matched: 2,
                mismatch: -3
            })
        );
        assert_eq!(uniform.codes(), &reference[..]);
        let hotspot = HotspotScoring::for_shm(AffineScoring::new(2, -3, -17, -2));
        assert_eq!(
            SubstitutionProfile::build(&hotspot, &reference).uniform(),
            None
        );
    }

    #[test]
    fn profile_matches_uniform_scoring() {
        let scoring = AffineScoring::new(2, -3, -5, -1);
        let reference = encode(b"ACGN");
        let profile = SubstitutionProfile::build(&scoring, &reference);
        for q in 0..5u8 {
            for (j, &r) in reference.iter().enumerate() {
                assert_eq!(profile.at(q, j), scoring.substitution(q, r));
            }
        }
        assert_eq!(scoring.gaps.cost(3), -7);
    }

    #[test]
    fn hotspot_softens_wrcy_and_hardens_syc() {
        // AGCT = WRCY (mutable C at 2); GTC = SYC (mutable C at 6).
        let reference = encode(b"AGCTGTC");
        let model = HotspotScoring::for_shm(AffineScoring::new(2, -3, -17, -2));
        assert!(is_hotspot_target(&reference, 2));
        assert_eq!(model.substitution(T, &reference, 2), -2);
        assert_eq!(model.substitution(A, &reference, 2), -2);
        assert_eq!(model.substitution(C, &reference, 0), -3);
        assert!(is_coldspot_target(&reference, 6));
        // coldspot penalty equals base mismatch in for_shm
        assert_eq!(model.substitution(T, &reference, 6), -3);
        assert_eq!(model.substitution(C, &reference, 6), 2);
    }
}
