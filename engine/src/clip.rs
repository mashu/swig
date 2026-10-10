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

/// One alignment column: its operation and the query and reference
/// positions before it consumes its bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Column {
    op: Op,
    q: u32,
    r: u32,
}

/// Column-indexed view of one alignment supporting O(1) clipping.
///
/// Three allocations per alignment: the columns, their prefix scores, and
/// both query-offset lookup tables in one vector.
#[derive(Clone, Debug)]
pub struct ClipIndex {
    q_start: u32,
    q_end: u32,
    columns: Vec<Column>,
    /// `prefix[k]` = score of columns `0..k`.
    prefix: Vec<Score>,
    /// `lookup[t]` for `t ≤ q_len`: first `Match` column whose query position
    /// is `≥ q_start + t`. `lookup[q_len + 1 + t]`: last `Match` column whose
    /// query position is `< q_start + t`.
    lookup: Vec<u32>,
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
        let q_len = (alignment.q_end - alignment.q_start) as usize;
        let mut columns = Vec::with_capacity(alignment.ops.len());
        let mut prefix = Vec::with_capacity(alignment.ops.len() + 1);
        // First the Match column of each query offset, then the two scans.
        let mut lookup = vec![NONE; 2 * (q_len + 1)];
        prefix.push(0);
        let mut previous: Option<Op> = None;
        let mut total: Score = 0;
        for (k, (op, q, r)) in alignment.columns().enumerate() {
            total += match op {
                Op::Match => {
                    lookup[(q - alignment.q_start) as usize] = k as u32;
                    profile.at(query[q as usize], r as usize)
                }
                gap if previous == Some(gap) => gaps.extend,
                _ => gaps.open,
            };
            columns.push(Column { op, q, r });
            prefix.push(total);
            previous = Some(op);
        }
        debug_assert_eq!(
            total, alignment.score,
            "column rescoring must reproduce the DP score"
        );
        let (first_from, last_before) = lookup.split_at_mut(q_len + 1);
        for t in 0..q_len {
            last_before[t + 1] = match first_from[t] {
                NONE => last_before[t],
                col => col,
            };
        }
        for t in (0..q_len).rev() {
            if first_from[t] == NONE {
                first_from[t] = first_from[t + 1];
            }
        }
        Self {
            q_start: alignment.q_start,
            q_end: alignment.q_end,
            columns,
            prefix,
            lookup,
        }
    }

    /// Canonical cuts, one per `Match` column in query order: its query
    /// position `q`, the prefix score of an entry cut at `q` (the column is
    /// kept first) and of an exit cut at `q + 1` (kept last). Every distinct
    /// clipped span has exactly one such pair of cuts.
    pub fn matches(&self) -> impl Iterator<Item = (u32, Score, Score)> + '_ {
        self.columns
            .iter()
            .enumerate()
            .filter(|(_, c)| c.op == Op::Match)
            .map(|(k, c)| (c.q, self.prefix[k], self.prefix[k + 1]))
    }

    #[inline]
    fn first_match_from(&self, offset: u32) -> u32 {
        self.lookup[offset as usize]
    }

    #[inline]
    fn last_match_before(&self, offset: u32) -> u32 {
        self.lookup[(self.q_len() + 1 + offset) as usize]
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
        match self.first_match_from(lo - self.q_start) {
            NONE => None,
            col => Some(Boundary {
                col,
                q: self.columns[col as usize].q,
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
        match self.last_match_before(hi - self.q_start) {
            NONE => None,
            col => Some(Boundary {
                col,
                q: self.columns[col as usize].q + 1,
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
            r_lo: self.columns[entry.col as usize].r,
            r_hi: self.columns[exit.col as usize].r + 1,
            score: exit.prefix - entry.prefix,
        })
    }

    /// The whole alignment as a span.
    pub fn full(&self) -> ClipSpan {
        self.span(self.q_start, self.q_end)
            .expect("aligner contract: alignments start and end with Match columns")
    }

    /// `(op, query_pos, ref_pos)` of the columns of `span`.
    pub fn columns(&self, span: &ClipSpan) -> impl Iterator<Item = (Op, u32, u32)> + use<'_> {
        (span.first_col as usize..=span.last_col as usize).map(move |k| {
            let Column { op, q, r } = self.columns[k];
            (op, q, r)
        })
    }

    /// The sub-alignment of `span`: its columns, bounds and score. It starts
    /// and ends with a `Match`, as a span does.
    pub fn alignment(&self, span: &ClipSpan) -> Alignment {
        Alignment {
            q_start: span.q_lo,
            q_end: span.q_hi,
            r_start: span.r_lo,
            r_end: span.r_hi,
            score: span.score,
            ops: self.columns(span).map(|(op, ..)| op).collect(),
        }
    }

    /// Query position aligned to reference position `r`, extrapolated along the
    /// alignment's end diagonals when `r` lies outside the aligned reference range.
    /// Used to project conserved anchors (V Cys, J Trp/Phe) onto the query.
    pub fn project_reference(&self, r: u32) -> i64 {
        let r_start = self.columns.first().map_or(0, |c| c.r);
        let r_end = self.columns.last().map_or(0, |c| c.r + 1);
        if r < r_start {
            return i64::from(self.q_start) - i64::from(r_start - r);
        }
        if r >= r_end {
            return i64::from(self.q_end) + i64::from(r - r_end);
        }
        let mut k = self.columns.partition_point(|c| c.r < r);
        while k < self.columns.len() && self.columns[k].op == Op::Ins {
            k += 1;
        }
        i64::from(self.columns[k.min(self.columns.len() - 1)].q)
    }

    /// `true` when both indices describe the same query columns with the same
    /// per-column scores, so every clip of one scores exactly like the other.
    /// Reference coordinates may differ (e.g. alleles with different 5' lengths).
    pub fn equivalent(&self, other: &Self) -> bool {
        self.q_start == other.q_start
            && self.prefix == other.prefix
            && self.columns.len() == other.columns.len()
            && self
                .columns
                .iter()
                .zip(&other.columns)
                .all(|(a, b)| a.op == b.op)
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
                // The span's sub-alignment is a whole alignment worth the span.
                if let Some(span) = clip.span(lo, hi) {
                    let sub = clip.alignment(&span);
                    let whole = naive_clip_score(&sub, &query, &reference, 0, u32::MAX);
                    assert_eq!(whole, Some(span.score));
                    assert_eq!((sub.q_start, sub.q_end), (span.q_lo, span.q_hi));
                }
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
