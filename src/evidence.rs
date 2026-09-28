//! Column-level evidence summaries of an alignment (or a clipped span of one).
//!
//! Admission gates (exact-run floor, distributed-D), D-presence priors and
//! reports all read the same [`Evidence`], computed in one pass over columns.

use crate::align::Op;
use crate::alphabet::N;
use crate::scoring::{AffineScoring, Score};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Evidence {
    /// All alignment columns.
    pub columns: u32,
    /// Columns consuming both query and reference (matches, mismatches, ambiguous).
    pub match_columns: u32,
    /// Identical canonical bases.
    pub matches: u32,
    /// Differing canonical bases.
    pub mismatches: u32,
    /// Number of gap runs.
    pub gap_opens: u32,
    /// Gap columns (insertions + deletions).
    pub gap_columns: u32,
    /// Longest run of consecutive identical canonical bases (`N` and gaps break runs).
    pub longest_exact_run: u32,
}

impl Evidence {
    pub fn from_columns(
        columns: impl IntoIterator<Item = (Op, u32, u32)>,
        query: &[u8],
        reference: &[u8],
    ) -> Self {
        let mut evidence = Evidence::default();
        let mut run = 0u32;
        let mut previous: Option<Op> = None;
        for (op, q, r) in columns {
            evidence.columns += 1;
            match op {
                Op::Match => {
                    evidence.match_columns += 1;
                    let (a, b) = (query[q as usize], reference[r as usize]);
                    if a < N && a == b {
                        evidence.matches += 1;
                        run += 1;
                        evidence.longest_exact_run = evidence.longest_exact_run.max(run);
                    } else {
                        if a < N && b < N {
                            evidence.mismatches += 1;
                        }
                        run = 0;
                    }
                }
                gap => {
                    evidence.gap_columns += 1;
                    if previous != Some(gap) {
                        evidence.gap_opens += 1;
                    }
                    run = 0;
                }
            }
            previous = Some(op);
        }
        evidence
    }

    /// Identical bases per alignment column (the AIRR `*_identity` convention).
    pub fn identity(&self) -> f64 {
        match self.columns {
            0 => 0.0,
            c => f64::from(self.matches) / f64::from(c),
        }
    }
}

/// Affine score of `columns` under a (possibly different) scoring tuple.
/// Columns must begin with a `Match` for gap accounting to be exact.
pub fn rescore(
    columns: impl IntoIterator<Item = (Op, u32, u32)>,
    query: &[u8],
    reference: &[u8],
    scoring: &AffineScoring,
) -> Score {
    let mut score = 0;
    let mut previous: Option<Op> = None;
    for (op, q, r) in columns {
        score += match op {
            Op::Match => scoring.substitution(query[q as usize], reference[r as usize]),
            gap if previous == Some(gap) => scoring.gaps.extend,
            _ => scoring.gaps.open,
        };
        previous = Some(op);
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;

    #[test]
    fn counts_runs_mismatches_and_gaps() {
        let query = encode(b"ACGTAACGG");
        let reference = encode(b"ACGTTACG");
        // 4 matches, 1 mismatch, 3 matches, 1 insertion
        let mut cols = Vec::new();
        for k in 0..8u32 {
            cols.push((Op::Match, k, k));
        }
        cols.push((Op::Ins, 8, 8));
        let e = Evidence::from_columns(cols.iter().copied(), &query, &reference);
        assert_eq!(
            (e.matches, e.mismatches, e.gap_columns, e.gap_opens),
            (7, 1, 1, 1)
        );
        assert_eq!(e.longest_exact_run, 4);
        let s = rescore(
            cols.iter().copied(),
            &query,
            &reference,
            &AffineScoring::new(2, -3, -13, -1),
        );
        assert_eq!(s, 7 * 2 - 3 - 13);
    }
}
