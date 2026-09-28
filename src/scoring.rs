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

use crate::alphabet::{ALPHABET_SIZE, N};

pub type Score = i32;

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

/// `ALPHABET_SIZE` rows × `len` columns of substitution scores for one allele.
#[derive(Clone, Debug)]
pub struct SubstitutionProfile {
    len: usize,
    rows: Vec<Score>,
    max: Score,
    min: Score,
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
        Self {
            len,
            rows,
            max,
            min,
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;

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
}
