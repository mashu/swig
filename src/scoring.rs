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
/// Not used by the shipped profiles: as a J profile it gained nothing on
/// simulated reads, whose SHM has no hotspot targeting (docs/ASSESSMENT.md).
/// It also exercises the aligners' position-specific profile path.
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

/// An anchored alignment as a local one: reference position 0 is an *anchor
/// column* that scores `bonus` against every query code, `N` included; every
/// other position scores `inner`. The reference carries the anchor as its
/// first position (`C′ = [N] + C`, see [`crate::germline::GermlineSet::with_anchor_column`]).
///
/// With `bonus = matched · ℓ + 1`, `ℓ` the longest `C`, any alignment that
/// avoids column 0 scores at most `matched · ℓ < bonus`, and a local
/// alignment entering column 0 can only start there (it is the first column;
/// a leading gap only lowers a score). So every optimal local alignment
/// starts with the anchor column, and
///
/// ```text
/// max local score = bonus + max over starts of the anchored score
/// ```
///
/// where a start is a band diagonal through column 0 and the anchored score
/// is the best alignment of `C` from its first base onward with a free end.
/// Any local aligner — every SIMD backend included — therefore computes the
/// exact anchored optimum, start set and drift corridor being the band.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AnchorColumn {
    pub inner: AffineScoring,
    pub bonus: Score,
}

impl AnchorColumn {
    /// The anchor model for references up to `longest` bases (anchor excluded).
    pub fn new(inner: AffineScoring, longest: usize) -> Self {
        let longest = Score::try_from(longest).expect("reference length fits a score");
        Self {
            inner,
            bonus: inner.matched.max(0) * longest + 1,
        }
    }
}

impl SubstitutionModel for AnchorColumn {
    #[inline]
    fn score(&self, query_code: u8, reference: &[u8], position: usize) -> Score {
        match position {
            0 => self.bonus,
            _ => self.inner.substitution(query_code, reference[position]),
        }
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
    /// Largest column maximum outside the first column attaining `max`.
    rest: Score,
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
        let mut columns = rows[..len].to_vec();
        for row in rows.chunks_exact(len.max(1)).skip(1) {
            for (column, &score) in columns.iter_mut().zip(row) {
                *column = (*column).max(score);
            }
        }
        let top = columns.iter().position(|&c| c == max);
        let rest = columns
            .iter()
            .enumerate()
            .filter(|&(j, _)| Some(j) != top)
            .map(|(_, &c)| c)
            .max()
            .unwrap_or(max);
        let codes: Vec<u8> = reference.iter().map(|&c| c.min(N)).collect();
        let uniform = uniform_scores(&rows, &codes);
        Self {
            len,
            rows,
            max,
            rest,
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

    /// Upper bound on any local alignment score with at most `k` substitution
    /// columns (`k = min(n, m)` for an `n`-base query):
    ///
    /// ```text
    /// S ≤ max⁺ + rest⁺ · (k − 1)
    /// ```
    ///
    /// with `max` the largest entry and `rest` the largest column maximum
    /// outside the first column attaining it, `x⁺ = max(x, 0)`. Each reference
    /// column is aligned at most once and gaps cost ≤ 0, so the alignment
    /// collects at most one entry from the top column and at most `rest` from
    /// each other. For [`AffineScoring`] and [`HotspotScoring`] every canonical
    /// column's maximum is `matched`, so this is `matched · k`; a profile with
    /// one outsized column (an anchor bonus) is bounded by that column once.
    pub fn local_bound(&self, k: usize) -> i64 {
        match k {
            0 => 0,
            k => i64::from(self.max.max(0)) + i64::from(self.rest.max(0)) * (k as i64 - 1),
        }
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

    /// `local_bound` is `matched · k` for the shipped models on references
    /// with two or more canonical bases, so lane admission is unchanged for
    /// every V/D/J job, and it bounds every local alignment score.
    #[test]
    fn local_bound_is_matched_per_column_and_bounds_local_scores() {
        use crate::align::{Aligner, BandedGotoh, GotohWorkspace};
        use crate::band::Band;
        use crate::testing::{random_codes, rng};
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let hotspot = HotspotScoring::for_shm(scoring);
        let mut next = rng(11);
        let mut ws = GotohWorkspace::default();
        for _ in 0..200 {
            let (m, n) = (3 + (next() % 60) as usize, 1 + (next() % 80) as usize);
            let mut reference = random_codes(&mut next, m);
            reference[(next() as usize) % m] = N;
            let query = random_codes(&mut next, n);
            let models: [&dyn SubstitutionModel; 2] = [&scoring, &hotspot];
            for model in models {
                let profile = SubstitutionProfile::build(model, &reference);
                for k in 0..=reference.len() {
                    assert_eq!(profile.local_bound(k), 2 * k as i64);
                }
                let k = query.len().min(reference.len());
                let hit = BandedGotoh.score(&query, &profile, scoring.gaps, Band::Full, &mut ws);
                let best = hit.map_or(0, |hit| i64::from(hit.score));
                assert!(best <= profile.local_bound(k));
            }
        }
    }

    /// An anchor profile is bounded by its bonus once plus `matched` per
    /// further column, and that bounds its local scores.
    #[test]
    fn anchor_profiles_count_the_bonus_once() {
        use crate::align::{Aligner, BandedGotoh, GotohWorkspace};
        use crate::band::Band;
        use crate::testing::{random_codes, rng};
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let model = AnchorColumn::new(scoring, 50);
        assert_eq!(model.bonus, 101);
        let mut next = rng(12);
        let mut ws = GotohWorkspace::default();
        for _ in 0..200 {
            let (m, n) = (2 + (next() % 50) as usize, 1 + (next() % 80) as usize);
            let reference = random_codes(&mut next, m);
            let query = random_codes(&mut next, n);
            let profile = SubstitutionProfile::build(&model, &reference);
            let k = n.min(m);
            assert_eq!(profile.local_bound(k), 101 + 2 * (k as i64 - 1));
            let hit = BandedGotoh.score(&query, &profile, scoring.gaps, Band::Full, &mut ws);
            assert!(i64::from(hit.expect("the anchor scores").score) <= profile.local_bound(k));
        }
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
