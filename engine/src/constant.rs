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
//! [`SpliceWindow`] sets the query window and the band, and
//! [`ReachClasses`] aligns one allele per class of alleles the window cannot
//! tell apart.
//!
//! **A stage after adjudication.** [`ConstantStage`] reads finished calls
//! only — a path whose last segment is J, the J's exit `x` — so V, D and J
//! calls are unchanged by construction, and a read without the stage costs
//! nothing. C never enters the placement lattice: a C kind there would let
//! C evidence move the J exit or the J allele, and every read would pay for
//! it. Its placements live in [`Call::constant`], clipped to `[x, n)` so C
//! explains only bases the J does not.
//!
//! **Significance.** Under random bases after `x`, the walk of one start
//! `(gene, δ)` reaches score `S` with probability at most `e^(−λS)`
//! (Lundberg's bound for an ungapped walk with `E[e^(λX)] = 1`; gaps are
//! mildly anti-conservative, as for every E-value of the crate). The union
//! over the `N_C` genes of the database and the `2s + 1` start diagonals
//! gives the anchored E-value
//!
//! ```text
//! E_C(S) = N_C · (2s + 1) · e^(−λS)
//! ```
//!
//! of the clipped score `S`, and C is called iff `E_C` is at most the gate,
//! by default `10⁻³`, the gate of a hit that must stand on its own
//! ([`crate::segment::EvalueGates::anchor`]): C is not implied by V(D)J, it
//! depends on the library. For `(2, −3)` scoring
//! (`λ = 0.634`) that is `S ≥ 18` for the 10 human IGHC genes. `N_C`
//! counts every gene of the C database, so one database per locus keeps the
//! gate tight. This is not the K·m·n E-value of V, D and J.
//!
//! **Posterior.** Among the read's C placements (the best `top` classes and
//! ties), a placement's weight is the log-odds of its alignment at C's own
//! divergence — C is not hypermutated, so the read's SHM rate `p̂` does not
//! apply — estimated from the called C's columns,
//! `p_C = (mismatches + ½)/(informative columns + 1)` clamped to
//! `[10⁻⁴, ½]`:
//!
//! ```text
//! w_p = exp(Σ_columns (λ_C·s + μ_C) + λ·S_gap)      (λ_C, μ_C) at p_C
//! c_posterior = Σ_{p tied or carried} w_p / Σ_p w_p
//! ```
//!
//! with ties and carried placements as for V, D and J
//! ([`crate::report`]): the probability that the reported names include
//! the read's C allele, given that the read reaches C. No C-presence
//! probability is reported.

use std::ops::Range;
use std::sync::Arc;

use crate::adjudicate::{CalledSegment, EvidenceScale, co_optimal};
use crate::align::{Aligner, Alignment, Op};
use crate::band::Band;
use crate::clip::ClipIndex;
use crate::config::CSegmentConfig;
use crate::germline::{GermlineSet, SegmentKind};
use crate::locus::{Loci, Locus};
use crate::pipeline::Call;
use crate::placement::{Placement, merge_identical};
use crate::retrieve::{Candidate, ReachClasses, RetrievalScratch};
use crate::scoring::{AffineScoring, AnchorColumn, GapCosts, Score};
use crate::segment::{BatchScratch, KeepPolicy, PreparedSegment, SegmentRequest, SegmentStage};
use crate::stats::KarlinAltschul;

/// The anchored E-value of a clipped C score (see the module documentation):
/// `E_C(S) = N_C · (2s + 1) · e^(−λS)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnchoredSignificance {
    /// Karlin–Altschul `λ` of the C scoring.
    pub lambda: f64,
    /// `ln(N_C · (2s + 1))`: genes × start diagonals.
    pub ln_hypotheses: f64,
    /// Largest E-value called.
    pub gate: f64,
}

impl AnchoredSignificance {
    /// `None` when the scoring has no `λ` (a non-negative expected score).
    pub fn new(set: &GermlineSet, scoring: &AffineScoring, spread: u32, gate: f64) -> Option<Self> {
        let hypotheses = set.gene_count() as f64 * f64::from(2 * spread + 1);
        Some(Self {
            lambda: KarlinAltschul::ungapped(scoring)?.lambda,
            ln_hypotheses: hypotheses.ln(),
            gate,
        })
    }

    pub fn evalue(&self, score: Score) -> f64 {
        (self.ln_hypotheses - self.lambda * f64::from(score)).exp()
    }

    /// Smallest score called: `⌈(ln(N_C·(2s + 1)) − ln gate) / λ⌉`.
    pub fn min_score(&self) -> Score {
        ((self.ln_hypotheses - self.gate.ln()) / self.lambda)
            .ceil()
            .max(1.0) as Score
    }
}

/// The C call of a read.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConstantCall {
    /// The called C: `placement` indexes [`Call::constant`], cut to `[x, n)`.
    pub segment: CalledSegment,
    /// Anchored E-value of its score, `E_C`.
    pub evalue: f64,
    /// Probability that the reported names include the read's C allele,
    /// given that the read reaches C.
    pub posterior: f64,
}

/// Reusable buffers of [`ConstantStage::assign`].
#[derive(Debug, Default)]
pub struct CScratch {
    windows: Vec<Option<SpliceWindow>>,
    candidates: Vec<Vec<Candidate>>,
    placed: Vec<Vec<Placement>>,
}

/// Constant-region assignment after adjudication (see the module
/// documentation).
pub struct ConstantStage {
    /// The anchored search over `C′ = [N] + C`: reach classes, banded DP
    /// with the anchor column, the best classes traced.
    pub anchored: SegmentStage<ReachClasses>,
    /// The C database, in whose coordinates C placements are reported.
    pub plain: Arc<PreparedSegment>,
    /// The anchor model of `anchored`.
    pub model: AnchorColumn,
    /// `s`: C's base 0 lies within `s` bases of the J exit.
    pub spread: u32,
    pub significance: AnchoredSignificance,
    /// For a C database of several loci, the gate of each: `N_C` counts the
    /// C genes of the read's locus only ([`crate::locus`]).
    pub locus_significance: Vec<(Locus, AnchoredSignificance)>,
}

impl ConstantStage {
    /// The stage for C database `set`; `None` when `config.scoring` has no
    /// Karlin–Altschul `λ`.
    pub fn new(set: GermlineSet, config: &CSegmentConfig) -> Option<Self> {
        let scoring = config.scoring;
        let significance = AnchoredSignificance::new(&set, &scoring, config.spread, config.evalue)?;
        let longest = set
            .alleles()
            .iter()
            .map(|a| a.codes.len())
            .max()
            .unwrap_or(0);
        let model = AnchorColumn::new(scoring, longest);
        let anchored = PreparedSegment::new(set.with_anchor_column(), &model, scoring.gaps);
        let keep = KeepPolicy {
            min_score: model.bonus + config.keep.min_score,
            ..config.keep
        };
        let retriever = ReachClasses::new(&set, config.spread);
        let loci: Vec<Locus> = Locus::ALL
            .into_iter()
            .filter(|&l| set.alleles().iter().any(|a| a.loci == Loci::of(l)))
            .collect();
        let locus_significance = if loci.len() > 1 {
            loci.into_iter()
                .filter_map(|locus| {
                    let (part, _) = set.subset(|a| a.loci.contains(locus));
                    let gate =
                        AnchoredSignificance::new(&part, &scoring, config.spread, config.evalue)?;
                    Some((locus, gate))
                })
                .collect()
        } else {
            Vec::new()
        };
        Some(Self {
            anchored: SegmentStage::new(Arc::new(anchored), retriever, keep),
            plain: Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps)),
            model,
            spread: config.spread,
            significance,
            locus_significance,
        })
    }

    /// The gate of a read whose J belongs to `loci`.
    fn significance(&self, loci: Loci) -> &AnchoredSignificance {
        loci.single()
            .and_then(|l| self.locus_significance.iter().find(|(at, _)| *at == l))
            .map_or(&self.significance, |(_, gate)| gate)
    }

    /// Whether some allele of the reach class of `representative` belongs
    /// to `loci`.
    fn class_pairs(&self, representative: u32, reach: usize, loci: Loci) -> bool {
        let set = self.plain.set();
        self.anchored
            .retriever
            .members(representative, reach)
            .iter()
            .any(|&a| set.allele(a).loci.pairs_with(loci))
    }

    /// The splice window of a call whose path ends with a J, if C can follow.
    pub fn window(&self, call: &Call) -> Option<SpliceWindow> {
        let last = call.path.as_ref()?.segments.last()?;
        if call.placements[last.placement].kind != SegmentKind::J {
            return None;
        }
        SpliceWindow::new(last.span.q_hi, call.query.len(), self.spread)
    }

    /// Place and call C for a chunk of finished calls, the J of call `i`
    /// belonging to `loci[i]`: one score batch and traceback waves for the
    /// whole chunk, then per read the C placements ([`Call::constant`]) of
    /// the J's locus and the call ([`Call::c`]).
    #[allow(clippy::too_many_arguments)]
    pub fn assign<A: Aligner>(
        &self,
        aligner: &A,
        calls: &mut [Call],
        loci: &[Loci],
        scratch: &mut CScratch,
        retrieval: &mut RetrievalScratch,
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        let CScratch {
            windows,
            candidates,
            placed,
        } = scratch;
        windows.clear();
        windows.extend(calls.iter().map(|call| self.window(call)));
        if placed.len() < calls.len() {
            candidates.resize_with(calls.len(), Vec::new);
            placed.resize_with(calls.len(), Vec::new);
        }
        let placed = &mut placed[..calls.len()];
        placed.iter_mut().for_each(Vec::clear);
        for (((window, found), call), &loci) in windows
            .iter()
            .zip(candidates.iter_mut())
            .zip(&*calls)
            .zip(loci)
        {
            if let Some(w) = window {
                self.anchored
                    .retrieve(&call.query, w.window.clone(), retrieval, found);
                let reach = self.anchored.retriever.reach(w.window.len());
                found.retain(|c| self.class_pairs(c.allele, reach, loci));
            }
        }
        let requests: Vec<SegmentRequest<'_>> = windows
            .iter()
            .zip(calls.iter())
            .enumerate()
            .filter_map(|(read, (window, call))| {
                window.as_ref().map(|w| SegmentRequest {
                    read,
                    query: &call.query,
                    window: w.window.clone(),
                    candidates: &candidates[read],
                })
            })
            .collect();
        self.anchored
            .place_many(aligner, &requests, dp, batch, placed);
        drop(requests);
        for (((call, window), anchored), &loci) in calls
            .iter_mut()
            .zip(windows.iter())
            .zip(placed.iter())
            .zip(loci)
        {
            if let Some(w) = window {
                self.settle(call, w, anchored, loci);
            }
        }
    }

    /// The read's C placements from its anchored ones: anchor stripped,
    /// clipped to `[x, n)`, class members joined, ranked by (score, longest
    /// exact run, matches, allele id) and merged; then the call.
    fn settle(&self, call: &mut Call, w: &SpliceWindow, anchored: &[Placement], loci: Loci) {
        let n = call.query.len() as u32;
        let reach = self.anchored.retriever.reach(w.window.len());
        for p in anchored {
            let Some(stripped) = strip_anchor(&p.alignment, &self.model) else {
                continue;
            };
            let representative = p.primary_allele();
            let profile = self.plain.profile(representative);
            let clip = ClipIndex::build(&stripped, &call.query, profile, self.plain.gaps());
            let Some(span) = clip.span(w.x, n) else {
                continue;
            };
            let mut placement = Placement::new(
                SegmentKind::C,
                representative,
                clip.alignment(&span),
                &call.query,
                &self.plain,
            );
            let set = self.plain.set();
            placement.alleles = p
                .alleles
                .iter()
                .flat_map(|&a| self.anchored.retriever.members(a, reach))
                .copied()
                .filter(|&a| set.allele(a).loci.pairs_with(loci))
                .collect();
            placement.alleles.sort_unstable();
            placement.alleles.dedup();
            // Members of a reach class align identically, so the first of
            // the locus reports for the class.
            if placement.alleles.is_empty() {
                continue;
            }
            call.constant.push(placement);
        }
        call.constant.sort_by_key(|p| {
            (
                std::cmp::Reverse(p.score()),
                std::cmp::Reverse(p.evidence.longest_exact_run),
                std::cmp::Reverse(p.evidence.matches),
                p.primary_allele(),
            )
        });
        merge_identical(&mut call.constant, 0);
        call.c = self.decide(call, w.x, self.significance(loci));
    }

    /// The call from the ranked C placements: the best, if significant
    /// under `gate`.
    fn decide(&self, call: &Call, x: u32, gate: &AnchoredSignificance) -> Option<ConstantCall> {
        let best = call.constant.first()?;
        if best.score() < gate.min_score() {
            return None;
        }
        let n = call.query.len() as u32;
        let segment = CalledSegment {
            placement: 0,
            entry: x,
            exit: n,
            span: best.clip.full(),
            allele: None,
        };
        let co = co_optimal(
            &call.constant,
            &call.query,
            self.plain.set(),
            &segment,
            false,
            x..n,
        );
        // Log-odds at C's own divergence, from the called columns.
        let AffineScoring {
            matched, mismatch, ..
        } = self.model.inner;
        let e = &best.evidence;
        let p = ((f64::from(e.mismatches) + 0.5) / f64::from(e.matches + e.mismatches + 1))
            .clamp(1e-4, 0.5);
        let scale = EvidenceScale::at_divergence(p, matched, mismatch);
        let nats = |q: &Placement| {
            let e = &q.evidence;
            let substitution = matched * e.matches as Score + mismatch * e.mismatches as Score;
            scale.lambda * f64::from(substitution)
                + scale.mu * f64::from(e.matches + e.mismatches)
                + gate.lambda * f64::from(q.score() - substitution)
        };
        let top = call
            .constant
            .iter()
            .map(nats)
            .fold(f64::NEG_INFINITY, f64::max);
        let weight = |q: &Placement| (nats(q) - top).exp();
        let total: f64 = call.constant.iter().map(weight).sum();
        let counted: f64 = co.counted.iter().map(|&i| weight(&call.constant[i])).sum();
        Some(ConstantCall {
            segment,
            evalue: gate.evalue(best.score()),
            posterior: (counted / total).min(1.0),
        })
    }
}

/// Where one read's anchored hypotheses are searched.
///
/// The anchor column takes the read base before `C[0]`, at `x + δ − 1`, so
/// the window starts at the earliest one, `x − s − 1`, and runs to the
/// read's end. The splice lies `s + 1` bases into every window, so the
/// anchor of `H(a, δ)` lies on diagonal `s + δ` of `C′` and the band is
///
/// ```text
/// lo = 0,   hi = 2s
/// ```
///
/// exactly the `2s + 1` start diagonals, and an alignment reaches at most
/// `|window| − 1 = (n − x) + s` bases of `C` ([`crate::retrieve::ReachClasses`]).
/// The band is also the corridor indels drift in: a C starting at `x + δ`
/// admits `s + δ` net deletions and `s − δ` net insertions; beyond that the
/// free end truncates the alignment. Widening the band would add start
/// hypotheses, which the significance counts.
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
    /// spread `s`; `None` when nothing follows `x`, or when `x ≤ s` leaves no
    /// room for every anchor (a called J is far longer than `s`).
    pub fn new(x: u32, read_len: usize, spread: u32) -> Option<Self> {
        let start = x.checked_sub(spread + 1)?;
        (read_len > x as usize).then(|| Self {
            x,
            window: start as usize..read_len,
            band: Band::Diagonals {
                lo: 0,
                hi: 2 * spread as i32,
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
        let x = s as usize + 1 + (next() % 8) as usize;
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

    /// The members of a reach class align exactly like their representative
    /// in the window the class was formed for.
    #[test]
    fn reach_class_members_align_identically() {
        use crate::alphabet::decode;
        use crate::germline::{GermlineSet, SegmentKind};
        use crate::retrieve::{ReachClasses, RetrievalScratch, Retriever};
        let mut next = rng(21);
        let root = random_codes(&mut next, 300);
        let mut set = GermlineSet::new(SegmentKind::C);
        for k in 0..24 {
            let mut allele = root[..150 + (next() % 150) as usize].to_vec();
            let p = (next() as usize) % allele.len();
            allele[p] = (allele[p] + 1) % 4;
            let name = format!("IGHC{}*{k:02}", k % 4);
            set.push(&name, decode(&allele).as_bytes(), None).unwrap();
        }
        let model = AnchorColumn::new(SCORING, 300);
        let anchored = set.with_anchor_column();
        let profiles: Vec<_> = anchored
            .alleles()
            .iter()
            .map(|a| SubstitutionProfile::build(&model, &a.codes))
            .collect();
        let classes = ReachClasses::new(&set, 3);
        let (mut scratch, mut candidates) = (RetrievalScratch::default(), Vec::new());
        let mut ws = GotohWorkspace::default();
        for _ in 0..200 {
            let (query, x) = read(&mut next, &set.allele(0).codes, 3);
            let w = SpliceWindow::new(x, query.len(), 3).unwrap();
            let window = &query[w.window.clone()];
            classes.retrieve(window, &mut scratch, &mut candidates);
            for candidate in &candidates {
                let align = |allele: u32, ws: &mut GotohWorkspace| {
                    let profile = &profiles[allele as usize];
                    BandedGotoh.align(window, profile, SCORING.gaps, w.band, ws)
                };
                let expected = align(candidate.allele, &mut ws);
                assert_eq!(candidate.band, w.band);
                for &member in classes.members(candidate.allele, classes.reach(window.len())) {
                    assert_eq!(align(member, &mut ws), expected);
                }
            }
        }
    }

    /// `S ≥ ⌈(ln(N_C·(2s + 1)) + ln 10³)/λ⌉`: 18 for the 10 human IGHC genes,
    /// 19 with IGKC and IGLC (16 genes), and the E-value crosses the gate there.
    #[test]
    fn the_anchored_gate_counts_genes_and_start_diagonals() {
        use crate::germline::{GermlineSet, SegmentKind};
        let genes = |count: usize| {
            let mut set = GermlineSet::new(SegmentKind::C);
            for g in 0..count {
                set.push(&format!("IGHC{g}*01"), b"ACGT", None).unwrap();
                set.push(&format!("IGHC{g}*02"), b"ACGA", None).unwrap();
            }
            set
        };
        let gate = |count: usize, spread: u32| {
            AnchoredSignificance::new(&genes(count), &SCORING, spread, 1e-3).unwrap()
        };
        assert_eq!(gate(10, 3).min_score(), 18);
        assert_eq!(gate(16, 3).min_score(), 19);
        assert!(gate(10, 6).min_score() > gate(10, 1).min_score());
        let human = gate(10, 3);
        assert!(human.evalue(18) <= 1e-3 && human.evalue(17) > 1e-3);
        assert!((human.lambda - 0.634).abs() < 1e-3);
    }

    #[test]
    fn the_window_holds_every_anchor() {
        let w = SpliceWindow::new(40, 100, 3).unwrap();
        assert_eq!(w.window, 36..100);
        // The splice 4 bases in: anchors at window positions 0..=6.
        assert_eq!(w.band, Band::Diagonals { lo: 0, hi: 6 });
        assert_eq!(SpliceWindow::new(4, 5, 3).unwrap().window, 0..5);
        assert_eq!(SpliceWindow::new(3, 100, 3), None);
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
