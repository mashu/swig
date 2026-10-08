//! Constant-region (C) assignment: an exact splice-anchored alignment right
//! of the called J.
//!
//! **Model.** J and C are joined by an RNA splice at a fixed phase: no N/P
//! bases, no hypermutation, and C's germline base 0 lies next to the J's
//! germline 3′ end. With `x` the called J's exit in the read, hypothesis
//! `H(a, δ)` puts base 0 of C allele `a` at read position `x + δ`, `|δ| ≤ s`
//! (reference conventions differ by a base — IMGT CH1 records start after
//! the J's last base, codon-aligned ones repeat it, 3′-partial J records end
//! a base early — and SHM indels near the J's 3′ end move `x`).
//!
//! **Trick — the anchor column.** An anchored alignment (start pinned to
//! `C[0]` on one of `2s + 1` diagonals, free end) is a local alignment of
//! `C′ = [N] + C` whose column 0 scores a bonus larger than any alignment
//! without it ([`AnchorColumn`]). Every optimal local alignment then starts
//! with column 0, so the existing banded kernels — scalar and every SIMD
//! backend — return the exact anchored optimum plus the bonus, and
//! [`strip_anchor`] turns the traceback into the alignment on `C`.
//! [`SpliceWindow`] sets the query window and the band.

use std::ops::Range;

use crate::align::{Alignment, Op};
use crate::band::Band;
use crate::scoring::{AnchorColumn, GapCosts};

/// Where one read's anchored hypotheses are searched.
///
/// The anchor column takes the read base before `C[0]`, at `x + δ − 1`, so
/// the window starts at the earliest one, `max(0, x − s − 1)`, and runs to
/// the read's end. In window coordinates (`x_rel = x − start`) the anchor of
/// `H(a, δ)` lies on diagonal `x_rel − 1 + δ` of `C′`, so the band is
///
/// ```text
/// lo = x_rel − 1 − s,   hi = x_rel − 1 + s
/// ```
///
/// exactly the `2s + 1` start diagonals (fewer when `x ≤ s`: an anchor
/// cannot precede the read). The band is also the corridor indels drift in:
/// a C starting at `x + δ` admits `s + δ` net deletions and `s − δ` net
/// insertions; beyond that the free end truncates the alignment. Widening
/// the band would add start hypotheses, which the significance counts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpliceWindow {
    /// Read position of the splice: the called J's exit.
    pub x: u32,
    /// Query window of the anchored alignment.
    pub window: Range<usize>,
    /// Band of the anchored alignment on `C′`, relative to the window.
    pub band: Band,
}

impl SpliceWindow {
    /// The window for a splice at `x` in a read of `read_len` bases with
    /// spread `s`; `None` when nothing follows `x`.
    pub fn new(x: u32, read_len: usize, spread: u32) -> Option<Self> {
        if x as usize >= read_len {
            return None;
        }
        let start = x.saturating_sub(spread + 1);
        let centre = (x - start) as i32 - 1;
        let spread = spread as i32;
        Some(Self {
            x,
            window: start as usize..read_len,
            band: Band::Diagonals {
                lo: centre - spread,
                hi: centre + spread,
            },
        })
    }
}

/// The alignment on `C` of a local alignment on `C′ = [N] + C` that starts
/// with the anchor column, as every optimal one does ([`AnchorColumn`]).
///
/// Column 0 is dropped (query start `+1`, reference coordinates `−1`, score
/// `− bonus`), and so are gap columns right after it — an alignment starts
/// with a `Match` ([`crate::align::Aligner`]'s contract, which clipping
/// relies on) — with their cost added back to the score. `None` when no
/// `Match` column remains.
pub fn strip_anchor(alignment: &Alignment, model: &AnchorColumn) -> Option<Alignment> {
    debug_assert!(
        alignment.r_start == 0 && alignment.ops.first() == Some(&Op::Match),
        "an optimal alignment starts with the anchor column"
    );
    let GapCosts { open, extend } = model.inner.gaps;
    let lead = alignment.ops[1..].iter().position(|&op| op == Op::Match)? + 1;
    let (mut q_start, mut r_start) = (alignment.q_start + 1, 0);
    let mut score = alignment.score - model.bonus;
    let mut previous = None;
    for &op in &alignment.ops[1..lead] {
        match op {
            Op::Ins => q_start += 1,
            _ => r_start += 1,
        }
        score -= if previous == Some(op) { extend } else { open };
        previous = Some(op);
    }
    Some(Alignment {
        q_start,
        q_end: alignment.q_end,
        r_start,
        r_end: alignment.r_end - 1,
        score,
        ops: alignment.ops[lead..].to_vec(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{
        AlignJob, Aligner, Backend, BandedGotoh, GotohWorkspace, InterSequenceGotoh,
    };
    use crate::alphabet::N;
    use crate::evidence::rescore;
    use crate::scoring::{AffineScoring, NEG_INF, Score, SubstitutionProfile};
    use crate::testing::{random_codes, rng};

    const SCORING: AffineScoring = AffineScoring::new(2, -3, -9, -1);

    /// Textbook anchored Gotoh, the oracle: `C[0]` enters at window position
    /// `a` for every start `a = x_rel + δ`, `|δ| ≤ s`, `1 ≤ a ≤ |query|`
    /// (`H(a, 0) = 0`, every other cell of column 0 −∞); cells keep
    /// `|(i − j) − x_rel| ≤ s`; the end is free. Returns the best score.
    fn anchored_oracle(query: &[u8], c: &[u8], x_rel: i64, s: i64) -> Score {
        let (n, m) = (query.len(), c.len());
        let in_band = |i: usize, j: usize| ((i as i64 - j as i64) - x_rel).abs() <= s;
        let GapCosts { open, extend } = SCORING.gaps;
        let mut h = vec![vec![NEG_INF; m + 1]; n + 1];
        let mut ins = h.clone();
        let mut del = h.clone();
        let mut best = NEG_INF;
        for i in 0..=n {
            for j in 0..=m {
                if !in_band(i, j) {
                    continue;
                }
                if j == 0 && i >= 1 {
                    // A start: the anchor took query base i − 1.
                    h[i][0] = 0;
                }
                if i > 0 && in_band(i - 1, j) {
                    ins[i][j] = (h[i - 1][j] + open).max(ins[i - 1][j] + extend);
                }
                if j > 0 && in_band(i, j - 1) {
                    del[i][j] = (h[i][j - 1] + open).max(del[i][j - 1] + extend);
                }
                let diagonal = if i > 0 && j > 0 {
                    h[i - 1][j - 1] + SCORING.substitution(query[i - 1], c[j - 1])
                } else {
                    NEG_INF
                };
                h[i][j] = h[i][j].max(diagonal).max(ins[i][j]).max(del[i][j]);
                best = best.max(h[i][j]);
            }
        }
        best
    }

    /// A read: `x` random J-side bases, then `C[..len]` from `δ` bases after
    /// the splice, mutated with substitutions, `N`s and indels, then a random
    /// tail; returns the read and `x`.
    fn read(next: &mut impl FnMut() -> u64, c: &[u8], s: i64) -> (Vec<u8>, u32) {
        let x = 1 + (next() % 12) as usize;
        let delta = (next() % (2 * s as u64 + 5)) as i64 - s - 2;
        let mut read = random_codes(next, x);
        if delta > 0 {
            read.extend(random_codes(next, delta as usize));
        } else {
            read.truncate((x as i64 + delta).max(0) as usize);
        }
        let len = (next() % (c.len() as u64 + 1)) as usize;
        let mut body = c[..len].to_vec();
        for _ in 0..(next() % 4) {
            if let Some(p) = (!body.is_empty()).then(|| (next() as usize) % body.len()) {
                body[p] = (next() % 5) as u8;
            }
        }
        for _ in 0..(next() % 3) {
            let p = (next() as usize) % (body.len() + 1);
            match next() % 2 {
                0 if p < body.len() => {
                    body.remove(p);
                }
                _ => body.insert(p, (next() % 4) as u8),
            }
        }
        read.extend(body);
        let tail = (next() % 20) as usize;
        read.extend(random_codes(next, tail));
        read.resize(read.len().max(x + 1), N);
        (read, x as u32)
    }

    /// The anchor column computes the textbook anchored optimum exactly on
    /// the scalar kernel and every SIMD backend, and every traced optimum
    /// starts with the anchor column on a start diagonal.
    #[test]
    fn the_anchor_column_is_the_exact_anchored_optimum_on_every_backend() {
        let mut next = rng(6);
        let mut cases = Vec::new();
        for case in 0..500 {
            let s = [1, 3, 6][case % 3];
            let len = 5 + (next() % 200) as usize;
            let c = random_codes(&mut next, len);
            let (query, x) = read(&mut next, &c, s);
            cases.push((s, c, query, x));
        }
        let longest = cases.iter().map(|(_, c, ..)| c.len()).max().unwrap();
        let model = AnchorColumn::new(SCORING, longest);
        let profiles: Vec<_> = cases
            .iter()
            .map(|(_, c, ..)| {
                let anchored: Vec<u8> = std::iter::once(N).chain(c.iter().copied()).collect();
                SubstitutionProfile::build(&model, &anchored)
            })
            .collect();
        let windows: Vec<_> = cases
            .iter()
            .map(|(s, _, query, x)| SpliceWindow::new(*x, query.len(), *s as u32).unwrap())
            .collect();
        let jobs: Vec<AlignJob> = cases
            .iter()
            .zip(&profiles)
            .zip(&windows)
            .map(|(((_, _, query, _), profile), w)| AlignJob {
                query: &query[w.window.clone()],
                profile,
                gaps: SCORING.gaps,
                band: w.band,
            })
            .collect();
        let mut expected = Vec::new();
        for ((s, c, query, x), w) in cases.iter().zip(&windows) {
            let x_rel = i64::from(*x) - w.window.start as i64;
            expected.push(anchored_oracle(&query[w.window.clone()], c, x_rel, *s));
        }
        let mut scalar = GotohWorkspace::default();
        let backends = Backend::ALL
            .into_iter()
            .filter_map(InterSequenceGotoh::with_backend);
        let (mut hits, mut alignments) = (Vec::new(), Vec::new());
        let mut ws = Default::default();
        for simd in backends {
            simd.score_batch(&jobs, &mut ws, &mut hits);
            simd.align_batch(&jobs, &mut ws, &mut alignments);
            for (k, job) in jobs.iter().enumerate() {
                let backend = simd.backend();
                let score = hits[k].map(|hit| hit.score - model.bonus);
                assert_eq!(score, Some(expected[k]), "{backend:?} case {k}");
                let local =
                    BandedGotoh.align(job.query, job.profile, job.gaps, job.band, &mut scalar);
                assert_eq!(alignments[k], local, "{backend:?} case {k}");
                let aligned = alignments[k].as_ref().unwrap();
                // Starts with the anchor column, on a start diagonal.
                assert_eq!((aligned.r_start, aligned.ops[0]), (0, Op::Match));
                let Band::Diagonals { lo, hi } = job.band else {
                    unreachable!()
                };
                assert!((lo..=hi).contains(&(aligned.q_start as i32)));
                // The traceback after the anchor rescores to the anchored
                // optimum (a gap right after the anchor opens).
                let c = &cases[k].1;
                let after = aligned.columns().skip(1).map(|(op, q, r)| (op, q, r - 1));
                let rescored = rescore(after, job.query, c, &SCORING);
                assert_eq!(rescored, expected[k], "{backend:?} rescoring, case {k}");
                assert_eq!(aligned.score, model.bonus + expected[k]);
                // Stripping leaves an alignment on C worth its own columns.
                if let Some(stripped) = strip_anchor(aligned, &model) {
                    assert_eq!(stripped.ops[0], Op::Match);
                    let own = rescore(stripped.columns(), job.query, c, &SCORING);
                    assert_eq!(stripped.score, own, "case {k}");
                    assert!(stripped.score >= expected[k]);
                }
            }
        }
    }

    #[test]
    fn the_window_holds_every_anchor_and_clamps_at_the_read_start() {
        let w = SpliceWindow::new(40, 100, 3).unwrap();
        assert_eq!(w.window, 36..100);
        // x_rel = 4: anchors at window positions 0..=6, diagonals 0..=6.
        assert_eq!(w.band, Band::Diagonals { lo: 0, hi: 6 });
        let w = SpliceWindow::new(2, 100, 3).unwrap();
        assert_eq!(w.window, 0..100);
        assert_eq!(w.band, Band::Diagonals { lo: -2, hi: 4 });
        assert_eq!(SpliceWindow::new(100, 100, 3), None);
    }

    #[test]
    fn stripping_drops_the_anchor_and_leading_gaps() {
        let model = AnchorColumn::new(SCORING, 10);
        let ops = [Op::Match, Op::Del, Op::Del, Op::Ins, Op::Match, Op::Match];
        let alignment = Alignment {
            q_start: 5,
            q_end: 9,
            r_start: 0,
            r_end: 5,
            score: model.bonus - 9 - 1 - 9 + 4,
            ops: ops.to_vec(),
        };
        let stripped = strip_anchor(&alignment, &model).unwrap();
        assert_eq!(
            (
                stripped.q_start,
                stripped.q_end,
                stripped.r_start,
                stripped.r_end
            ),
            (7, 9, 2, 4)
        );
        assert_eq!((stripped.score, stripped.ops.len()), (4, 2));
        let only_anchor = Alignment {
            q_start: 5,
            q_end: 6,
            r_start: 0,
            r_end: 1,
            score: model.bonus,
            ops: vec![Op::Match],
        };
        assert_eq!(strip_anchor(&only_anchor, &model), None);
    }
}
