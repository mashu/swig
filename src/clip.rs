//! O(1) clip-and-rescore of an alignment restricted to a query interval.
//!
//! Adjudication constantly asks "what is V worth if it must stop before query
//! position `c`?" — for every candidate cut `c`, every V×D×J combination. The
//! method note answers by clipping the alignment and rescoring it. Done naively
//! that is `O(alignment length)` per question.
//!
//! **Trick — prefix sums over alignment columns.** At placement time each
//! alignment is flattened into columns with a prefix-score array, plus two
//! lookup tables over query offsets: the first `Match` column at or after a
//! position, and the last `Match` column before a position. A clipped span is
//! then two table lookups and one subtraction.
//!
//! Exactness: a clipped span always starts and ends on a `Match` column, so no
//! gap is cut in half and the affine open/extend accounting of every retained
//! column is unchanged. The span score is therefore *exactly* the affine score
//! of the retained sub-alignment — identical to re-running the rescoring loop.
//!
//! Separability: `score(lo, hi) = exit(hi).prefix − entry(lo).prefix`. The
//! entry term depends only on the left cut and the exit term only on the right
//! cut. [`crate::adjudicate::PlacementDag`] exploits this to evaluate all cuts
//! of a transition with prefix maxima instead of a double loop.

use crate::align::{Alignment, Op};
use crate::scoring::{GapCosts, Score, SubstitutionProfile};

const NONE: u32 = u32::MAX;

/// Column-indexed view of one alignment supporting O(1) clipping.
#[derive(Clone, Debug)]
pub struct ClipIndex {
    q_start: u32,
    q_end: u32,
    col_op: Vec<Op>,
    col_q: Vec<u32>,
    col_r: Vec<u32>,
    /// `prefix[k]` = score of columns `0..k`.
    prefix: Vec<Score>,
    /// For query offset `t` (relative to `q_start`), first `Match` column whose
    /// query position is `≥ q_start + t`.
    first_m_from: Vec<u32>,
    /// For query offset `t`, last `Match` column whose query position is `< q_start + t`.
    last_m_before: Vec<u32>,
}

/// One side of a clipped span.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Boundary {
    /// Alignment column of the boundary `Match`.
    pub col: u32,
    /// Query coordinate of the boundary: first retained base for an entry,
    /// one past the last retained base for an exit.
    pub q: u32,
    /// Prefix score at the boundary (`prefix[col]` for an entry, `prefix[col + 1]` for an exit).
    pub prefix: Score,
}

/// A retained sub-alignment. Query and reference coordinates are 0-based half-open.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ClipSpan {
    pub first_col: u32,
    /// Inclusive.
    pub last_col: u32,
    pub q_lo: u32,
    pub q_hi: u32,
    pub r_lo: u32,
    pub r_hi: u32,
    pub score: Score,
}

impl ClipSpan {
    pub fn q_len(&self) -> u32 {
        self.q_hi - self.q_lo
    }

    pub fn r_len(&self) -> u32 {
        self.r_hi - self.r_lo
    }
}

impl ClipIndex {
    /// Flatten `alignment` (query coordinates absolute in `query`) against the
    /// profile and gap costs it was computed with.
    pub fn build(
        alignment: &Alignment,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
    ) -> Self {
        let columns = alignment.ops.len();
        let q_len = (alignment.q_end - alignment.q_start) as usize;
        let mut col_op = Vec::with_capacity(columns);
        let mut col_q = Vec::with_capacity(columns);
        let mut col_r = Vec::with_capacity(columns);
        let mut prefix = Vec::with_capacity(columns + 1);
        let mut match_col_of_q = vec![NONE; q_len];
        prefix.push(0);
        let mut previous: Option<Op> = None;
        let mut total: Score = 0;
        for (k, (op, q, r)) in alignment.columns().enumerate() {
            total += match op {
                Op::Match => {
                    match_col_of_q[(q - alignment.q_start) as usize] = k as u32;
                    profile.at(query[q as usize], r as usize)
                }
                gap if previous == Some(gap) => gaps.extend,
                _ => gaps.open,
            };
            col_op.push(op);
            col_q.push(q);
            col_r.push(r);
            prefix.push(total);
            previous = Some(op);
        }
        debug_assert_eq!(
            total, alignment.score,
            "column rescoring must reproduce the DP score"
        );

        let mut first_m_from = vec![NONE; q_len + 1];
        for t in (0..q_len).rev() {
            first_m_from[t] = match match_col_of_q[t] {
                NONE => first_m_from[t + 1],
                col => col,
            };
        }
        let mut last_m_before = vec![NONE; q_len + 1];
        for t in 0..q_len {
            last_m_before[t + 1] = match match_col_of_q[t] {
                NONE => last_m_before[t],
                col => col,
            };
        }
        Self {
            q_start: alignment.q_start,
            q_end: alignment.q_end,
            col_op,
            col_q,
            col_r,
            prefix,
            first_m_from,
            last_m_before,
        }
    }

    pub fn q_start(&self) -> u32 {
        self.q_start
    }

    pub fn q_end(&self) -> u32 {
        self.q_end
    }

    pub fn q_len(&self) -> u32 {
        self.q_end - self.q_start
    }

    /// Left boundary when the retained part must start at or after query position `lo`.
    #[inline]
    pub fn entry(&self, lo: u32) -> Option<Boundary> {
        let lo = lo.max(self.q_start);
        if lo >= self.q_end {
            return None;
        }
        match self.first_m_from[(lo - self.q_start) as usize] {
            NONE => None,
            col => Some(Boundary {
                col,
                q: self.col_q[col as usize],
                prefix: self.prefix[col as usize],
            }),
        }
    }

    /// Right boundary when the retained part must end before query position `hi`.
    #[inline]
    pub fn exit(&self, hi: u32) -> Option<Boundary> {
        let hi = hi.min(self.q_end);
        if hi <= self.q_start {
            return None;
        }
        match self.last_m_before[(hi - self.q_start) as usize] {
            NONE => None,
            col => Some(Boundary {
                col,
                q: self.col_q[col as usize] + 1,
                prefix: self.prefix[col as usize + 1],
            }),
        }
    }

    /// The sub-alignment restricted to query interval `[lo, hi)`, or `None` if
    /// no `Match` column survives.
    #[inline]
    pub fn span(&self, lo: u32, hi: u32) -> Option<ClipSpan> {
        let entry = self.entry(lo)?;
        let exit = self.exit(hi)?;
        (entry.col <= exit.col).then(|| ClipSpan {
            first_col: entry.col,
            last_col: exit.col,
            q_lo: entry.q,
            q_hi: exit.q,
            r_lo: self.col_r[entry.col as usize],
            r_hi: self.col_r[exit.col as usize] + 1,
            score: exit.prefix - entry.prefix,
        })
    }

    /// The whole alignment as a span.
    pub fn full(&self) -> ClipSpan {
        self.span(self.q_start, self.q_end)
            .expect("aligner contract: alignments start and end with Match columns")
    }

    /// `(op, query_pos, ref_pos)` of the columns of `span`.
    pub fn columns(&self, span: &ClipSpan) -> impl Iterator<Item = (Op, u32, u32)> + '_ {
        (span.first_col as usize..=span.last_col as usize)
            .map(move |k| (self.col_op[k], self.col_q[k], self.col_r[k]))
    }

    /// Query position aligned to reference position `r`, extrapolated along the
    /// alignment's end diagonals when `r` lies outside the aligned reference range.
    /// Used to project conserved anchors (V Cys, J Trp/Phe) onto the query.
    pub fn project_reference(&self, r: u32) -> i64 {
        let r_start = self.col_r.first().copied().unwrap_or(0);
        let r_end = self.col_r.last().map_or(0, |&last| last + 1);
        if r < r_start {
            return i64::from(self.q_start) - i64::from(r_start - r);
        }
        if r >= r_end {
            return i64::from(self.q_end) + i64::from(r - r_end);
        }
        let mut k = self.col_r.partition_point(|&x| x < r);
        while k < self.col_op.len() && self.col_op[k] == Op::Ins {
            k += 1;
        }
        i64::from(self.col_q[k.min(self.col_q.len() - 1)])
    }

    /// `true` when both indices describe the same query columns with the same
    /// per-column scores, so every clip of one scores exactly like the other.
    /// Reference coordinates may differ (e.g. alleles with different 5' lengths).
    pub fn equivalent(&self, other: &Self) -> bool {
        self.q_start == other.q_start && self.col_op == other.col_op && self.prefix == other.prefix
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Aligner, BandedGotoh, GotohWorkspace};
    use crate::alphabet::encode;
    use crate::band::Band;
    use crate::scoring::AffineScoring;

    const SCORING: AffineScoring = AffineScoring::new(2, -3, -5, -1);

    fn naive_clip_score(
        aln: &Alignment,
        query: &[u8],
        reference: &[u8],
        lo: u32,
        hi: u32,
    ) -> Option<Score> {
        let cols: Vec<_> = aln.columns().collect();
        let first = cols
            .iter()
            .position(|&(op, q, _)| op == Op::Match && q >= lo)?;
        let last = cols
            .iter()
            .rposition(|&(op, q, _)| op == Op::Match && q < hi)?;
        if first > last {
            return None;
        }
        let mut score = 0;
        let mut previous = None;
        for &(op, q, r) in &cols[first..=last] {
            score += match op {
                Op::Match => SCORING.substitution(query[q as usize], reference[r as usize]),
                gap if previous == Some(gap) => SCORING.gaps.extend,
                _ => SCORING.gaps.open,
            };
            previous = Some(op);
        }
        Some(score)
    }

    #[test]
    fn clip_scores_equal_naive_rescoring() {
        let reference = encode(b"ACGTTGCAAGCTTGACCGTAGGCTAACGTTAGCCAT");
        let query = encode(b"TTACGTTGCAAGCGGGTTGACCGTAGCTAACGTTAGCCATGG");
        let profile = SubstitutionProfile::build(&SCORING, &reference);
        let mut ws = GotohWorkspace::default();
        let aln = BandedGotoh
            .align(&query, &profile, SCORING.gaps, Band::Full, &mut ws)
            .unwrap();
        assert!(
            aln.ops.iter().any(|&op| op != Op::Match),
            "test wants a gapped alignment"
        );
        let clip = ClipIndex::build(&aln, &query, &profile, SCORING.gaps);
        assert_eq!(clip.full().score, aln.score);
        for lo in 0..query.len() as u32 {
            for hi in lo + 1..=query.len() as u32 {
                let fast = clip.span(lo, hi).map(|s| s.score);
                assert_eq!(
                    fast,
                    naive_clip_score(&aln, &query, &reference, lo, hi),
                    "[{lo},{hi})"
                );
            }
        }
    }

    #[test]
    fn projects_reference_positions() {
        let reference = encode(b"ACGTTGCAAGCTTGACCGTA");
        let query = encode(b"GGGACGTTGCAAGCTTGACC");
        let profile = SubstitutionProfile::build(&SCORING, &reference);
        let mut ws = GotohWorkspace::default();
        let aln = BandedGotoh
            .align(&query, &profile, SCORING.gaps, Band::Full, &mut ws)
            .unwrap();
        let clip = ClipIndex::build(&aln, &query, &profile, SCORING.gaps);
        assert_eq!(clip.project_reference(0), 3);
        assert_eq!(clip.project_reference(10), 13);
        assert_eq!(clip.project_reference(19), 22);
    }
}
