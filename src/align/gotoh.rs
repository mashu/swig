//! Banded affine local alignment (Gotoh recurrence, Smith–Waterman boundary).
//!
//! Exactly the recurrence of the method note:
//!
//! ```text
//! I(i,j) = max(S(i-1,j) + g_o, I(i-1,j) + g_e)   gap in reference (query insertion)
//! D(i,j) = max(S(i,j-1) + g_o, D(i,j-1) + g_e)   gap in query (reference deletion)
//! S(i,j) = max(0, S(i-1,j-1) + s(q_i, r_j), I(i,j), D(i,j))
//! ```
//!
//! Implementation tricks:
//!
//! * **Rolling rows.** One `H` and one `I` array of length `m + 1` are updated
//!   in place, with the diagonal carried in a register. Out-of-band cells read as
//!   `H = 0`, `I = -∞`; for a *local* recurrence `H = 0` is equivalent to a fresh
//!   start, so the result is exact inside the band.
//! * **Band-only traceback.** Traceback is stored for in-band cells only, with
//!   per-row offsets: `O(n·w)` bytes rather than `O(n·m)`, allocated once per
//!   call from the band's cell count.
//! * **Compact traceback codes.** Each cell stores its `S` source in 2 bits plus
//!   one "extends" bit per gap state — one byte per cell, no `I`/`D` matrices.
//! * **Const-generic pass.** Score-only and traceback passes are one function
//!   specialised by `const TRACE: bool`; after monomorphisation the score-only
//!   loop carries no traceback code at all.
//! * **Branch-free, bounds-check-free inner loop.** LLVM turns `i32::max` into
//!   data-dependent branches here, which mispredict on real alignments (3×
//!   slower than on the synthetic, perfectly predictable inputs typical of
//!   micro-benchmarks). Maxima are therefore computed with bit arithmetic,
//!   traceback codes with branch-free selects, and rows are pre-sliced so the
//!   zipped iteration carries no bounds checks. The best cell is tracked as a
//!   row maximum and located once per improving row.
//! * **Short loop-carried chain.** The horizontal gap state obeys
//!   `D(i,j+1) = max(S(i,j)+g_o, D(i,j)+g_e)` with `S = max(t, D)`, where
//!   `t = max(0, diag, I)` does not depend on `D`. Hence
//!   `D(i,j+1) = max(t + g_o, D(i,j) + max(g_o, g_e))`: the serial dependency
//!   per cell shrinks from three maxima to one add and one max, and the
//!   result is identical for any gap costs.

use super::{Aligner, Alignment, Op, ScoreHit};
use crate::band::Band;
use crate::scoring::{GapCosts, NEG_INF, Score, SubstitutionProfile};

// Traceback codes, shared with the SIMD kernel (`super::inter_seq`).
pub(super) const SRC_ZERO: u8 = 0;
pub(super) const SRC_DIAG: u8 = 1;
pub(super) const SRC_INS: u8 = 2;
pub(super) const SRC_DEL: u8 = 3;
const SRC_MASK: u8 = 0b0011;
pub(super) const INS_EXTENDS: u8 = 0b0100; // == 1 << 2, set arithmetically in `fill`
pub(super) const DEL_EXTENDS: u8 = 0b1000; // == 1 << 3

/// Scalar banded Gotoh aligner; the reference implementation of [`Aligner`].
#[derive(Clone, Copy, Debug, Default)]
pub struct BandedGotoh;

/// Reusable buffers for [`BandedGotoh`].
#[derive(Debug, Default)]
pub struct GotohWorkspace {
    h: Vec<Score>,
    ins: Vec<Score>,
    trace: Vec<u8>,
    row_first_col: Vec<u32>,
    row_offset: Vec<u32>,
}

impl GotohWorkspace {
    #[inline]
    fn trace_cell(&self, row: usize, col: usize) -> Option<u8> {
        if row == 0 || col == 0 || row + 1 >= self.row_offset.len() {
            return None;
        }
        let start = self.row_offset[row] as usize;
        let width = (self.row_offset[row + 1] as usize).checked_sub(start)?;
        let first = self.row_first_col[row] as usize;
        (col >= first && col - first < width).then(|| self.trace[start + col - first])
    }
}

#[derive(Clone, Copy, Debug)]
struct Best {
    score: Score,
    row: usize,
    col: usize,
}

/// `max` that always compiles to straight-line code. LLVM lowers `i32::max`
/// to a branch in this loop, which mispredicts on real alignments.
#[inline(always)]
fn max_branchless(a: Score, b: Score) -> Score {
    a ^ ((a ^ b) & -Score::from(a < b))
}

fn fill<const TRACE: bool>(
    query: &[u8],
    profile: &SubstitutionProfile,
    gaps: GapCosts,
    band: Band,
    ws: &mut GotohWorkspace,
) -> Option<Best> {
    let (n, m) = (query.len(), profile.len());
    if n == 0 || m == 0 {
        return None;
    }
    let GotohWorkspace {
        h,
        ins,
        trace,
        row_first_col,
        row_offset,
    } = ws;
    h.clear();
    h.resize(m + 1, 0);
    ins.clear();
    ins.resize(m + 1, NEG_INF);
    if TRACE {
        // Every in-band cell is written before traceback reads it: grow only.
        let cells = band.cells(n, m);
        if trace.len() < cells {
            trace.resize(cells, 0);
        }
        row_first_col.clear();
        row_first_col.resize(n + 1, 0);
        row_offset.clear();
        row_offset.resize(n + 2, 0);
    }
    let (open, extend) = (gaps.open, gaps.extend);
    // D(i,j+1) = max(S(i,j) + go, D(i,j) + ge) with S = max(t, D), t = max(0, diag, I):
    //          = max(t + go, D + max(go, ge)).  `t` is off the loop-carried chain.
    let del_step = open.max(extend);
    let mut best = Best {
        score: 0,
        row: 0,
        col: 0,
    };

    for i in 1..=n {
        let Some((lo, hi)) = band.columns(i, m) else {
            if TRACE {
                row_offset[i + 1] = row_offset[i];
            }
            if band.exhausted(i, m) {
                break;
            }
            continue;
        };
        let width = hi - lo + 1;
        let trace_row: &mut [u8] = if TRACE {
            let offset = row_offset[i] as usize;
            row_first_col[i] = lo as u32;
            row_offset[i + 1] = (offset + width) as u32;
            &mut trace[offset..offset + width]
        } else {
            &mut []
        };

        // Pre-sliced rows: the zip below compiles without bounds checks.
        let substitution = &profile.row(query[i - 1])[lo - 1..hi];
        let mut h_diag = h[lo - 1];
        let h_row = &mut h[lo..=hi];
        let ins_row = &mut ins[lo..=hi];
        let mut del: Score = NEG_INF;
        let mut del_extends = false;
        let mut row_max: Score = 0;

        for (k, ((h_cell, ins_cell), &s)) in h_row
            .iter_mut()
            .zip(ins_row.iter_mut())
            .zip(substitution)
            .enumerate()
        {
            let h_up = *h_cell;
            let ins_open = h_up + open;
            let ins_extend = *ins_cell + extend;
            let ins_score = max_branchless(ins_open, ins_extend);
            let diag = h_diag + s;
            let t0 = max_branchless(diag, 0);
            let t = max_branchless(t0, ins_score);
            let score = max_branchless(t, del);

            if TRACE {
                // Codes are ordered like the preference (0 < diag < ins < del on strict
                // improvement), so "the last condition that holds" is plain arithmetic.
                let source = u8::from(diag > 0) * SRC_DIAG;
                let source = source + u8::from(ins_score > t0) * (SRC_INS - source);
                let source = source + u8::from(del > t) * (SRC_DEL - source);
                trace_row[k] =
                    source | (u8::from(ins_extend > ins_open) << 2) | (u8::from(del_extends) << 3);
            }

            let del_open = t + open;
            let del_extend = del + del_step;
            del_extends = del + extend > del_open;
            del = max_branchless(del_open, del_extend);
            h_diag = h_up;
            *h_cell = score;
            *ins_cell = ins_score;
            row_max = max_branchless(row_max, score);
        }
        if row_max > best.score {
            let k = h_row
                .iter()
                .position(|&x| x == row_max)
                .expect("row maximum is attained in the row");
            best = Best {
                score: row_max,
                row: i,
                col: lo + k,
            };
        }
    }
    (best.score > 0).then_some(best)
}

/// Walk traceback codes from the best cell `(row, col)`, 1-based.
///
/// `cell(i, j)` returns the code of an in-band cell and `None` for row 0,
/// column 0 and out-of-band cells. Every kernel that produces the same codes
/// therefore produces the same alignment.
pub(super) fn traceback_with(
    cell: impl Fn(usize, usize) -> Option<u8>,
    row: usize,
    col: usize,
) -> (usize, usize, Vec<Op>) {
    enum State {
        Diagonal,
        Insertion,
        Deletion,
    }
    let (mut i, mut j) = (row, col);
    let mut ops = Vec::with_capacity(row.min(col) + 8);
    let mut state = State::Diagonal;
    loop {
        match state {
            State::Diagonal => {
                let Some(code) = cell(i, j) else { break };
                match code & SRC_MASK {
                    SRC_DIAG => {
                        ops.push(Op::Match);
                        i -= 1;
                        j -= 1;
                    }
                    SRC_INS => state = State::Insertion,
                    SRC_DEL => state = State::Deletion,
                    SRC_ZERO => break,
                    _ => unreachable!("source codes are two bits"),
                }
            }
            State::Insertion => {
                let code = cell(i, j).expect("an insertion run stays inside the band");
                ops.push(Op::Ins);
                i -= 1;
                if code & INS_EXTENDS == 0 {
                    state = State::Diagonal;
                }
            }
            State::Deletion => {
                let code = cell(i, j).expect("a deletion run stays inside the band");
                ops.push(Op::Del);
                j -= 1;
                if code & DEL_EXTENDS == 0 {
                    state = State::Diagonal;
                }
            }
        }
    }
    ops.reverse();
    (i, j, ops)
}

impl Aligner for BandedGotoh {
    type Workspace = GotohWorkspace;

    fn score(
        &self,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
        band: Band,
        workspace: &mut GotohWorkspace,
    ) -> Option<ScoreHit> {
        let best = fill::<false>(query, profile, gaps, band, workspace)?;
        Some(ScoreHit {
            score: best.score,
            q_end: best.row as u32,
            r_end: best.col as u32,
        })
    }

    fn align(
        &self,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
        band: Band,
        workspace: &mut GotohWorkspace,
    ) -> Option<Alignment> {
        let best = fill::<true>(query, profile, gaps, band, workspace)?;
        let (q_start, r_start, ops) =
            traceback_with(|i, j| workspace.trace_cell(i, j), best.row, best.col);
        Some(Alignment {
            q_start: q_start as u32,
            q_end: best.row as u32,
            r_start: r_start as u32,
            r_end: best.col as u32,
            score: best.score,
            ops,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::scoring::AffineScoring;

    const SCORING: AffineScoring = AffineScoring::new(2, -3, -5, -1);

    use crate::testing::{random_codes, rng};

    fn rescore(aln: &Alignment, query: &[u8], reference: &[u8]) -> Score {
        let mut score = 0;
        let mut previous = None;
        for (op, q, r) in aln.columns() {
            score += match op {
                Op::Match => SCORING.substitution(query[q as usize], reference[r as usize]),
                gap if previous == Some(gap) => SCORING.gaps.extend,
                _ => SCORING.gaps.open,
            };
            previous = Some(op);
        }
        score
    }

    #[test]
    fn identical_sequences_align_fully() {
        let reference = encode(b"ACGTACGGTACCA");
        let profile = SubstitutionProfile::build(&SCORING, &reference);
        let mut ws = GotohWorkspace::default();
        let aln = BandedGotoh
            .align(&reference, &profile, SCORING.gaps, Band::Full, &mut ws)
            .unwrap();
        assert_eq!(aln.score, 26);
        assert_eq!(
            (aln.q_start, aln.q_end, aln.r_start, aln.r_end),
            (0, 13, 0, 13)
        );
        assert!(aln.ops.iter().all(|&op| op == Op::Match));
    }

    #[test]
    fn finds_affine_insertion() {
        let reference = encode(b"ACGTACGGTACCATTGACCA");
        let query = encode(b"ACGTACGGTAGGGCCATTGACCA");
        let profile = SubstitutionProfile::build(&SCORING, &reference);
        let mut ws = GotohWorkspace::default();
        let aln = BandedGotoh
            .align(&query, &profile, SCORING.gaps, Band::Full, &mut ws)
            .unwrap();
        assert_eq!(aln.score, 20 * 2 + SCORING.gaps.cost(3));
        assert_eq!(aln.ops.iter().filter(|&&op| op == Op::Ins).count(), 3);
        assert_eq!(rescore(&aln, &query, &reference), aln.score);
    }

    #[test]
    fn wide_band_equals_full_matrix_and_tracebacks_rescore() {
        let mut next = rng(7);
        let mut ws = GotohWorkspace::default();
        for _ in 0..200 {
            let len = 40 + (next() % 40) as usize;
            let reference = random_codes(&mut next, len);
            let mut query = reference.clone();
            for _ in 0..(next() % 6) {
                let p = (next() as usize) % query.len();
                query[p] = (next() % 4) as u8;
            }
            if next().is_multiple_of(2) {
                let p = (next() as usize) % query.len();
                query.insert(p, (next() % 4) as u8);
            }
            let profile = SubstitutionProfile::build(&SCORING, &reference);
            let full = BandedGotoh.score(&query, &profile, SCORING.gaps, Band::Full, &mut ws);
            let banded = BandedGotoh.score(
                &query,
                &profile,
                SCORING.gaps,
                Band::around(0, 100),
                &mut ws,
            );
            assert_eq!(full, banded);
            let aln =
                BandedGotoh.align(&query, &profile, SCORING.gaps, Band::around(0, 8), &mut ws);
            if let Some(aln) = aln {
                assert_eq!(rescore(&aln, &query, &reference), aln.score);
                assert_eq!(aln.ops.first(), Some(&Op::Match));
                assert_eq!(aln.ops.last(), Some(&Op::Match));
            }
        }
    }

    /// Textbook three-matrix Gotoh, written for clarity, as an oracle.
    fn naive_local_score(query: &[u8], reference: &[u8], scoring: &AffineScoring) -> Score {
        let (n, m) = (query.len(), reference.len());
        let (go, ge) = (scoring.gaps.open, scoring.gaps.extend);
        let mut s = vec![vec![0; m + 1]; n + 1];
        let mut ins = vec![vec![NEG_INF; m + 1]; n + 1];
        let mut del = vec![vec![NEG_INF; m + 1]; n + 1];
        let mut best = 0;
        for i in 1..=n {
            for j in 1..=m {
                ins[i][j] = (s[i - 1][j] + go).max(ins[i - 1][j] + ge);
                del[i][j] = (s[i][j - 1] + go).max(del[i][j - 1] + ge);
                let diag = s[i - 1][j - 1] + scoring.substitution(query[i - 1], reference[j - 1]);
                s[i][j] = 0.max(diag).max(ins[i][j]).max(del[i][j]);
                best = best.max(s[i][j]);
            }
        }
        best
    }

    #[test]
    fn matches_textbook_gotoh_for_any_gap_costs() {
        let mut next = rng(19);
        let mut ws = GotohWorkspace::default();
        for case in 0..300 {
            let (open, extend) = match case % 3 {
                0 => (-5, -1),
                1 => (-13, -2),
                _ => (-1, -4), // opening cheaper than extending: score still exact
            };
            let affine = open <= extend;
            let scoring = AffineScoring {
                matched: 2,
                mismatch: -3,
                gaps: GapCosts { open, extend },
            };
            let len = 20 + (next() % 30) as usize;
            let reference = random_codes(&mut next, len);
            let mut query = reference.clone();
            for _ in 0..(next() % 5) {
                let p = (next() as usize) % query.len();
                query.remove(p);
            }
            for _ in 0..(next() % 5) {
                let p = (next() as usize) % query.len();
                query.insert(p, (next() % 4) as u8);
            }
            let profile = SubstitutionProfile::build(&scoring, &reference);
            let fast = BandedGotoh
                .score(&query, &profile, scoring.gaps, Band::Full, &mut ws)
                .map_or(0, |h| h.score);
            assert_eq!(
                fast,
                naive_local_score(&query, &reference, &scoring),
                "case {case}"
            );
            // Column accounting ("consecutive gaps extend") is only defined for open <= extend.
            if let Some(aln) = BandedGotoh
                .align(&query, &profile, scoring.gaps, Band::Full, &mut ws)
                .filter(|_| affine)
            {
                assert_eq!(aln.score, fast);
                let mut rescored = 0;
                let mut previous = None;
                for (op, q, r) in aln.columns() {
                    rescored += match op {
                        Op::Match => scoring.substitution(query[q as usize], reference[r as usize]),
                        gap if previous == Some(gap) => scoring.gaps.extend,
                        _ => scoring.gaps.open,
                    };
                    previous = Some(op);
                }
                assert_eq!(
                    rescored, fast,
                    "traceback must realise the optimum (case {case})"
                );
            }
        }
    }

    #[test]
    fn narrow_band_can_only_lose_score() {
        let mut next = rng(11);
        let mut ws = GotohWorkspace::default();
        for _ in 0..100 {
            let reference = random_codes(&mut next, 60);
            let query = random_codes(&mut next, 60);
            let profile = SubstitutionProfile::build(&SCORING, &reference);
            let full = BandedGotoh.score(&query, &profile, SCORING.gaps, Band::Full, &mut ws);
            let narrow =
                BandedGotoh.score(&query, &profile, SCORING.gaps, Band::around(5, 3), &mut ws);
            assert!(narrow.map_or(0, |h| h.score) <= full.map_or(0, |h| h.score));
        }
    }
}
