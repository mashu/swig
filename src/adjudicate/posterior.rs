//! Posterior decoding: forward–backward (sum-product) over the placement DAG.
//!
//! The path objective is a log-likelihood ratio in units of `1/λ` nats (see
//! [`super::EngineeringObjective`]), so
//!
//! ```text
//! P(path, cuts | read) = exp(λ · score(path, cuts)) / Z
//! ```
//!
//! is a distribution over every grammatical path through the placements and
//! every pair of cuts between neighbours. The Viterbi path of
//! [`PlacementDag`] is its mode. Summing instead of maximising gives the
//! marginal probability of each placement (an allele call) and of each D
//! count (D presence).
//!
//! **Domain.** As in the Viterbi DAG, the first segment keeps its start and
//! the last segment its end: V and J are extended to the read's ends on
//! purpose ([`crate::refine::extend_outer`]). Between neighbours, every pair
//! of cuts is summed, including entries beyond the overlap that the Viterbi
//! search skips because they cannot win: they lose nothing for the maximum
//! but carry mass. Each segmentation is counted once by taking *canonical*
//! cuts: an entry cut on the first retained base, an exit cut right after
//! the last.
//!
//! **Separability.** The DAG's two shortcuts carry over with sums for maxima:
//! prefix sums over a node's entries replace running maxima over its states,
//! and running sums over its exits replace the running best exit. The
//! backward pass mirrors both with suffix sums. A read costs
//! `O(Σ_node len + Σ_edges entries)` multiply–adds.
//!
//! **Arithmetic.** Values are linear, not logarithmic. Each `(node, group)`
//! slot carries one log-scale (as in scaled HMM recursions) and is
//! renormalised before use, and `exp(λk)` comes from a table, so no
//! transcendental function is evaluated per cut. Within a placement, weights
//! are centred on the midpoint of its prefix-score range; a placement whose
//! range exceeds [`MAX_SPREAD`] nats would leave the `f64` range, and the
//! read then gets no marginals.

use super::model::{Grammar, PathModel};
use super::PlacementDag;
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// Largest prefix-score range of one placement, in nats (`λ · score`).
pub const MAX_SPREAD: f64 = 600.0;

/// Scores whose weight `exp(λk)` is tabulated (`|k| ≤ TABLE_RANGE`).
const TABLE_RANGE: Score = 1024;

/// Posterior marginals of one read's placements.
#[derive(Clone, Debug, PartialEq)]
pub struct Marginals {
    /// `P(placement i is on the path)`, indexed like the placement slice.
    pub placement: Vec<f64>,
    /// `P(the path has exactly k D segments)`, for `k = 0, 1, …`.
    pub d_count: Vec<f64>,
    /// `ln Z`, with `Z = Σ_{path, cuts} exp(λ · score)`.
    pub log_partition: f64,
}

impl Marginals {
    /// Probability that the path contains at least one D.
    pub fn d_present(&self) -> f64 {
        self.d_count.iter().skip(1).sum()
    }
}

/// `exp(λk)` for integer scores `k`.
#[derive(Clone, Debug)]
struct Weights {
    lambda: f64,
    table: Vec<f64>,
}

impl Weights {
    fn new(lambda: f64) -> Self {
        Self {
            lambda,
            table: (-TABLE_RANGE..=TABLE_RANGE)
                .map(|k| (lambda * f64::from(k)).exp())
                .collect(),
        }
    }

    #[inline]
    fn of(&self, score: Score) -> f64 {
        match usize::try_from(score + TABLE_RANGE) {
            Ok(i) if i < self.table.len() => self.table[i],
            _ => (self.lambda * f64::from(score)).exp(),
        }
    }
}

/// One placement in query order, with its canonical cut weights.
#[derive(Clone, Copy, Debug)]
struct Node {
    placement: usize,
    kind: SegmentKind,
    q_start: u32,
    /// `q_end − q_start`; cut offsets run over `0..=len`.
    len: u32,
    /// Start of this node's dense weight arrays in `entry_w` / `exit_w`.
    cuts: usize,
    /// This node's canonical entry offsets in `entry_cuts`.
    entries: (usize, usize),
}

impl Node {
    fn q_end(&self) -> u32 {
        self.q_start + self.len
    }

    fn width(&self) -> usize {
        self.len as usize + 1
    }
}

/// Forward and backward values of one `(node, D count, first kind)` group,
/// each `scale + ln(values)`.
#[derive(Clone, Copy, Debug)]
struct Slot {
    d_used: u8,
    first: SegmentKind,
    /// Start of this slot's values in `forward` / `backward`.
    values: usize,
    forward_scale: f64,
    backward_scale: f64,
}

/// Per-thread buffers for [`PlacementDag::posterior`].
#[derive(Debug, Default)]
pub struct PosteriorScratch {
    weights: Option<Weights>,
    order: Vec<usize>,
    nodes: Vec<Node>,
    /// `exp(λ(centre − entry prefix))` at canonical entries, else 0.
    entry_w: Vec<f64>,
    /// `exp(λ(exit prefix − centre))` at canonical exits, else 0.
    exit_w: Vec<f64>,
    entry_cuts: Vec<u32>,
    node_slots: Vec<Vec<usize>>,
    slots: Vec<Slot>,
    forward: Vec<f64>,
    backward: Vec<f64>,
    sums: Vec<f64>,
    exits: Vec<f64>,
    running: Vec<f64>,
    inner: Vec<f64>,
}

/// `ln(e^a + e^b)`.
fn log_add(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    if lo == f64::NEG_INFINITY {
        hi
    } else {
        hi + (lo - hi).exp().ln_1p()
    }
}

/// Divide `values` by their maximum; returns its log (`−∞` when all are 0).
fn normalise(values: &mut [f64]) -> f64 {
    let max = values.iter().copied().fold(0.0, f64::max);
    if max <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let inverse = 1.0 / max;
    values.iter_mut().for_each(|v| *v *= inverse);
    max.ln()
}

/// Bring a slot to a common scale with an incoming batch at `incoming`;
/// returns the factor the batch must be multiplied by.
fn align_scale(scale: &mut f64, values: &mut [f64], incoming: f64) -> f64 {
    if *scale == f64::NEG_INFINITY {
        *scale = incoming;
        1.0
    } else if incoming > *scale {
        let shrink = (*scale - incoming).exp();
        values.iter_mut().for_each(|v| *v *= shrink);
        *scale = incoming;
        1.0
    } else {
        (incoming - *scale).exp()
    }
}

fn find_slot(
    node_slots: &[usize],
    slots: &[Slot],
    d_used: u8,
    first: SegmentKind,
) -> Option<usize> {
    node_slots
        .iter()
        .copied()
        .find(|&s| slots[s].d_used == d_used && slots[s].first == first)
}

impl<G: Grammar, M: PathModel> PlacementDag<G, M> {
    /// Posterior marginals of `placements` under `P(path, cuts) ∝ exp(λ·score)`.
    ///
    /// `None` when no grammatical path exists, or a placement's score range
    /// is too wide for linear arithmetic ([`MAX_SPREAD`]).
    pub fn posterior(
        &self,
        placements: &[Placement],
        lambda: f64,
        scratch: &mut PosteriorScratch,
    ) -> Option<Marginals> {
        let PosteriorScratch {
            weights,
            order,
            nodes,
            entry_w,
            exit_w,
            entry_cuts,
            node_slots,
            slots,
            forward,
            backward,
            sums,
            exits,
            running,
            inner,
        } = scratch;
        if weights.as_ref().map_or(true, |w| w.lambda != lambda) {
            *weights = Some(Weights::new(lambda));
        }
        let w = weights.as_ref().expect("set above");
        let n = placements.len();
        if n == 0 {
            return None;
        }

        // Geometry: nodes in DAG order with dense canonical-cut weights.
        order.clear();
        order.extend(0..n);
        order.sort_by_key(|&i| (placements[i].q_start(), placements[i].q_end(), i));
        nodes.clear();
        entry_w.clear();
        exit_w.clear();
        entry_cuts.clear();
        for &i in order.iter() {
            let p = &placements[i];
            let (q_start, len) = (p.q_start(), p.q_end() - p.q_start());
            let entry = |o: u32| {
                let y = q_start + o;
                p.clip.entry(y).filter(|b| o < len && b.q == y)
            };
            let exit = |o: u32| {
                let x = q_start + o;
                p.clip.exit(x).filter(|b| o > 0 && b.q == x)
            };
            let (mut lo, mut hi) = (Score::MAX, Score::MIN);
            for o in 0..=len {
                for b in [entry(o), exit(o)].into_iter().flatten() {
                    lo = lo.min(b.prefix);
                    hi = hi.max(b.prefix);
                }
            }
            if lo <= hi && lambda * f64::from(hi - lo) > MAX_SPREAD {
                return None;
            }
            let centre = if lo <= hi { lo + (hi - lo) / 2 } else { 0 };
            let (cuts, first_entry) = (entry_w.len(), entry_cuts.len());
            for o in 0..=len {
                let ew = entry(o).map_or(0.0, |b| w.of(centre - b.prefix));
                if ew > 0.0 {
                    entry_cuts.push(o);
                }
                entry_w.push(ew);
                exit_w.push(exit(o).map_or(0.0, |b| w.of(b.prefix - centre)));
            }
            nodes.push(Node {
                placement: i,
                kind: p.kind,
                q_start,
                len,
                cuts,
                entries: (first_entry, entry_cuts.len()),
            });
        }

        // Slots are created on first use; starts get weight 1 per entry.
        if node_slots.len() < n {
            node_slots.resize_with(n, Vec::new);
        }
        node_slots.iter_mut().for_each(Vec::clear);
        slots.clear();
        forward.clear();
        backward.clear();
        let slot_of = |node: usize,
                       d_used: u8,
                       first: SegmentKind,
                       node_slots: &mut Vec<Vec<usize>>,
                       slots: &mut Vec<Slot>,
                       forward: &mut Vec<f64>,
                       backward: &mut Vec<f64>| {
            if let Some(s) = find_slot(&node_slots[node], slots, d_used, first) {
                return s;
            }
            let values = forward.len();
            forward.resize(values + nodes[node].width(), 0.0);
            backward.resize(values + nodes[node].width(), 0.0);
            slots.push(Slot {
                d_used,
                first,
                values,
                forward_scale: f64::NEG_INFINITY,
                backward_scale: f64::NEG_INFINITY,
            });
            node_slots[node].push(slots.len() - 1);
            slots.len() - 1
        };
        for (j, node) in nodes.iter().enumerate() {
            if !self.grammar.can_start(node.kind) {
                continue;
            }
            let d0 = u8::from(node.kind == SegmentKind::D);
            let s = slot_of(j, d0, node.kind, node_slots, slots, forward, backward);
            let range = slots[s].values..slots[s].values + node.width();
            if entry_w[node.cuts] > 0.0 {
                let factor = align_scale(
                    &mut slots[s].forward_scale,
                    &mut forward[range.clone()],
                    0.0,
                );
                forward[range.start] += factor;
            }
        }

        // Forward: F(e, y) = Σ_b Σ_{x ≤ y} X_b(x) · exp(λι(y − x)), with
        // X_b(x) = exp(λ(exit prefix + bonus)) · Σ_{y' ≤ x − min} F(b, y') · exp(−λ entry prefix).
        let mut log_z: Vec<f64> = Vec::new();
        for i in 0..n {
            let b = nodes[i];
            let min_b = self.min_bases(b.kind) as usize;
            let width = b.width();
            let successors = nodes.partition_point(|e| e.q_start <= b.q_start);
            for k in 0..node_slots[i].len() {
                let s = node_slots[i][k];
                let Slot {
                    d_used,
                    first,
                    values,
                    ..
                } = slots[s];
                let f = &mut forward[values..values + width];
                slots[s].forward_scale += normalise(f);
                let forward_scale = slots[s].forward_scale;
                if forward_scale == f64::NEG_INFINITY {
                    continue;
                }
                sums.clear();
                let mut acc = 0.0;
                for (o, &value) in f.iter().enumerate() {
                    acc += value * entry_w[b.cuts + o];
                    sums.push(acc);
                }
                let bonus = w.of(self.model.segment_bonus(&placements[b.placement], d_used));
                exits.clear();
                exits.extend((0..width).map(|o| {
                    let xw = exit_w[b.cuts + o];
                    if o >= min_b && xw > 0.0 {
                        sums[o - min_b] * xw * bonus
                    } else {
                        0.0
                    }
                }));
                let exit_scale = forward_scale + normalise(exits);
                if exit_scale == f64::NEG_INFINITY {
                    continue;
                }
                if self.grammar.can_end(first, b.kind) && exits[width - 1] > 0.0 {
                    let end = w.of(self.model.end(first, b.kind, d_used));
                    let d = usize::from(d_used);
                    if log_z.len() <= d {
                        log_z.resize(d + 1, f64::NEG_INFINITY);
                    }
                    log_z[d] = log_add(log_z[d], exit_scale + (exits[width - 1] * end).ln());
                }
                for kind in SegmentKind::ALL {
                    if !self.grammar.can_follow(b.kind, kind, d_used) {
                        continue;
                    }
                    let iota = self.model.junction_per_base(b.kind, kind);
                    running.clear();
                    let mut acc = 0.0;
                    for (o, &x) in exits.iter().enumerate() {
                        acc += x * w.of(-iota * o as Score);
                        running.push(acc);
                    }
                    let d_next = d_used + u8::from(kind == SegmentKind::D);
                    for j in successors..n {
                        let e = nodes[j];
                        if e.kind != kind {
                            continue;
                        }
                        let t = slot_of(j, d_next, first, node_slots, slots, forward, backward);
                        let range = slots[t].values..slots[t].values + e.width();
                        let factor = align_scale(
                            &mut slots[t].forward_scale,
                            &mut forward[range.clone()],
                            exit_scale,
                        );
                        for &o in &entry_cuts[e.entries.0..e.entries.1] {
                            let y = e.q_start + o;
                            let z = (y.min(b.q_end()) - b.q_start) as usize;
                            let gap = w.of(iota * (y - b.q_start) as Score);
                            forward[range.start + o as usize] += running[z] * gap * factor;
                        }
                    }
                }
            }
        }
        let log_partition = log_z.iter().copied().fold(f64::NEG_INFINITY, log_add);
        if !log_partition.is_finite() {
            return None;
        }

        // Backward: B(b, y) = exp(−λ entry prefix) · Σ_{x ≥ y + min} Out_b(x), with
        // Out_b(x) = exp(λ(exit prefix + bonus)) · ([x = end of b]·end + Σ_e Σ_{y' ≥ x} exp(λι(y' − x)) · B(e, y')).
        for i in (0..n).rev() {
            let b = nodes[i];
            let min_b = self.min_bases(b.kind) as usize;
            let width = b.width();
            let successors = nodes.partition_point(|e| e.q_start <= b.q_start);
            for k in 0..node_slots[i].len() {
                let s = node_slots[i][k];
                let Slot {
                    d_used,
                    first,
                    values,
                    forward_scale,
                    ..
                } = slots[s];
                if forward_scale == f64::NEG_INFINITY {
                    continue;
                }
                let can_end = self.grammar.can_end(first, b.kind);
                let target = |kind: SegmentKind, j: usize| {
                    let d_next = d_used + u8::from(kind == SegmentKind::D);
                    find_slot(&node_slots[j], slots, d_next, first)
                        .filter(|&t| slots[t].backward_scale > f64::NEG_INFINITY)
                };
                let mut reference = if can_end { 0.0 } else { f64::NEG_INFINITY };
                for kind in SegmentKind::ALL {
                    if !self.grammar.can_follow(b.kind, kind, d_used) {
                        continue;
                    }
                    for j in (successors..n).filter(|&j| nodes[j].kind == kind) {
                        if let Some(t) = target(kind, j) {
                            reference = reference.max(slots[t].backward_scale);
                        }
                    }
                }
                if reference == f64::NEG_INFINITY {
                    continue;
                }
                let end = if can_end {
                    w.of(self.model.end(first, b.kind, d_used)) * (-reference).exp()
                } else {
                    0.0
                };
                inner.clear();
                inner.resize(width, 0.0);
                inner[width - 1] = end;
                for kind in SegmentKind::ALL {
                    if !self.grammar.can_follow(b.kind, kind, d_used) {
                        continue;
                    }
                    let iota = self.model.junction_per_base(b.kind, kind);
                    running.clear();
                    running.resize(width, 0.0);
                    for j in (successors..n).filter(|&j| nodes[j].kind == kind) {
                        let Some(t) = target(kind, j) else {
                            continue;
                        };
                        let e = nodes[j];
                        let factor = (slots[t].backward_scale - reference).exp();
                        for &o in &entry_cuts[e.entries.0..e.entries.1] {
                            let y = e.q_start + o;
                            let z = (y.min(b.q_end()) - b.q_start) as usize;
                            let gap = w.of(iota * (y - b.q_start) as Score);
                            running[z] += backward[slots[t].values + o as usize] * gap * factor;
                        }
                    }
                    let mut acc = 0.0;
                    for o in (0..width).rev() {
                        acc += running[o];
                        inner[o] += w.of(-iota * o as Score) * acc;
                    }
                }
                let bonus = w.of(self.model.segment_bonus(&placements[b.placement], d_used));
                sums.clear();
                sums.resize(width + 1, 0.0);
                for o in (0..width).rev() {
                    sums[o] = sums[o + 1] + exit_w[b.cuts + o] * bonus * inner[o];
                }
                let out = &mut backward[values..values + width];
                for (o, value) in out.iter_mut().enumerate() {
                    let ew = entry_w[b.cuts + o];
                    *value = if ew > 0.0 && o + min_b < width {
                        ew * sums[o + min_b]
                    } else {
                        0.0
                    };
                }
                slots[s].backward_scale = reference + normalise(out);
            }
        }

        let mut placement = vec![0.0; n];
        for (i, node) in nodes.iter().enumerate() {
            for &s in &node_slots[i] {
                let slot = slots[s];
                let scale = slot.forward_scale + slot.backward_scale;
                if scale == f64::NEG_INFINITY {
                    continue;
                }
                let range = slot.values..slot.values + node.width();
                let mass: f64 = forward[range.clone()]
                    .iter()
                    .zip(&backward[range])
                    .map(|(f, b)| f * b)
                    .sum();
                if mass > 0.0 {
                    placement[node.placement] += (scale + mass.ln() - log_partition).exp();
                }
            }
        }
        Some(Marginals {
            placement,
            d_count: log_z.iter().map(|z| (z - log_partition).exp()).collect(),
            log_partition,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adjudicate::{Adjudicator, DagScratch, EngineeringObjective, VdjGrammar};
    use crate::align::{Alignment, Op};
    use crate::alphabet::encode;
    use crate::clip::ClipIndex;
    use crate::evidence::Evidence;
    use crate::scoring::{AffineScoring, SubstitutionProfile};

    const LAMBDA: f64 = 0.6337;

    fn placement(kind: SegmentKind, query: &[u8], reference: &[u8], q_start: u32) -> Placement {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let profile = SubstitutionProfile::build(&scoring, reference);
        let len = reference.len() as u32;
        let score = (0..reference.len())
            .map(|k| profile.at(query[q_start as usize + k], k))
            .sum();
        let alignment = Alignment {
            q_start,
            q_end: q_start + len,
            r_start: 0,
            r_end: len,
            score,
            ops: vec![Op::Match; reference.len()],
        };
        Placement {
            kind,
            alleles: vec![0],
            clip: ClipIndex::build(&alignment, query, &profile, scoring.gaps),
            evidence: Evidence::default(),
            anchor: true,
            alignment,
        }
    }

    /// Exhaustive sum over grammatical paths and canonical cut pairs.
    struct Brute<'a> {
        dag: &'a PlacementDag,
        placements: &'a [Placement],
        z: f64,
        best: Score,
        placement: Vec<f64>,
        d_count: Vec<f64>,
    }

    impl<'a> Brute<'a> {
        fn entries(p: &Placement) -> Vec<u32> {
            (p.q_start()..p.q_end())
                .filter(|&y| p.clip.entry(y).is_some_and(|b| b.q == y))
                .collect()
        }

        fn exits(p: &Placement) -> Vec<u32> {
            (p.q_start() + 1..=p.q_end())
                .filter(|&x| p.clip.exit(x).is_some_and(|b| b.q == x))
                .collect()
        }

        fn visit(
            &mut self,
            path: &mut Vec<usize>,
            entry: u32,
            d_used: u8,
            first: SegmentKind,
            acc: Score,
        ) {
            let (dag, placements) = (self.dag, self.placements);
            let b = &placements[*path.last().expect("non-empty")];
            let min_b = dag.min_retained[b.kind.index()].max(1);
            for x in Self::exits(b).into_iter().filter(|&x| x >= entry + min_b) {
                let span = b.clip.span(entry, x).expect("canonical cuts");
                let segment = span.score + dag.model.segment_bonus(b, d_used);
                if dag.grammar.can_end(first, b.kind) && x == b.q_end() {
                    let total = acc + segment + dag.model.end(first, b.kind, d_used);
                    let weight = (LAMBDA * f64::from(total)).exp();
                    self.z += weight;
                    self.best = self.best.max(total);
                    for &p in path.iter() {
                        self.placement[p] += weight;
                    }
                    self.d_count[usize::from(d_used)] += weight;
                }
                for (j, e) in placements.iter().enumerate() {
                    if e.q_start() <= b.q_start() || !dag.grammar.can_follow(b.kind, e.kind, d_used)
                    {
                        continue;
                    }
                    let iota = dag.model.junction_per_base(b.kind, e.kind);
                    let d_next = d_used + u8::from(e.kind == SegmentKind::D);
                    for y in Self::entries(e).into_iter().filter(|&y| y >= x) {
                        path.push(j);
                        self.visit(
                            path,
                            y,
                            d_next,
                            first,
                            acc + segment + iota * (y - x) as Score,
                        );
                        path.pop();
                    }
                }
            }
        }

        fn run(dag: &'a PlacementDag, placements: &'a [Placement]) -> Self {
            let mut brute = Brute {
                dag,
                placements,
                z: 0.0,
                best: Score::MIN,
                placement: vec![0.0; placements.len()],
                d_count: vec![0.0; 3],
            };
            for (s, p) in placements.iter().enumerate() {
                if dag.grammar.can_start(p.kind) {
                    let d0 = u8::from(p.kind == SegmentKind::D);
                    brute.visit(&mut vec![s], p.q_start(), d0, p.kind, 0);
                }
            }
            brute
        }
    }

    fn assert_close(a: f64, b: f64, what: &str) {
        assert!(
            (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs())),
            "{what}: {a} vs {b}"
        );
    }

    fn check(dag: &PlacementDag, placements: &[Placement]) {
        let brute = Brute::run(dag, placements);
        let m = dag
            .posterior(placements, LAMBDA, &mut PosteriorScratch::default())
            .expect("a path exists");
        assert_close(m.log_partition, brute.z.ln(), "ln Z");
        for (i, (&got, &want)) in m.placement.iter().zip(&brute.placement).enumerate() {
            assert_close(got, want / brute.z, &format!("placement {i}"));
        }
        for (k, &want) in brute.d_count.iter().enumerate() {
            let got = m.d_count.get(k).copied().unwrap_or(0.0);
            assert_close(got, want / brute.z, &format!("D count {k}"));
        }
        // The Viterbi path is the mode of the same distribution.
        let path = dag
            .adjudicate(0, placements, &mut DagScratch::default())
            .expect("a path exists");
        assert_eq!(path.score.score, brute.best);
    }

    /// Competing V alleles, two D candidates (one spurious, overlapping J)
    /// and two J alleles, with overlaps on every side.
    fn junction() -> (Vec<u8>, Vec<Placement>) {
        let query = encode(b"ACGTTGCAAGTCGGATTACAGGCTTAACGGTACCATGGA");
        let p = |kind, reference: &[u8], at: u32| placement(kind, &query, &encode(reference), at);
        let placements = vec![
            p(SegmentKind::V, b"ACGTTGCAAGTCGG", 0),
            p(SegmentKind::V, b"ACGTAGCAAGTCGA", 0),
            p(SegmentKind::D, b"GATTACAGG", 13),
            p(SegmentKind::D, b"TACAGCCTT", 16),
            p(SegmentKind::J, b"GCTTAACGGTAC", 21),
            p(SegmentKind::J, b"CTTTACGGTACC", 22),
        ];
        (query, placements)
    }

    #[test]
    fn marginals_equal_exhaustive_summation() {
        let (_, placements) = junction();
        let model = EngineeringObjective {
            d_presence: 4,
            ..EngineeringObjective::default()
        };
        for grammar in [VdjGrammar::VDJ, VdjGrammar::TANDEM_D, VdjGrammar::VJ] {
            check(&PlacementDag::new(grammar, model, [4, 3, 4]), &placements);
        }
        let per_base = EngineeringObjective {
            insertion_per_base: -1,
            ..model
        };
        check(
            &PlacementDag::new(VdjGrammar::TANDEM_D, per_base, [4, 3, 4]),
            &placements,
        );
    }

    #[test]
    fn marginals_are_probabilities() {
        let (_, placements) = junction();
        let dag = PlacementDag::new(VdjGrammar::VDJ, EngineeringObjective::default(), [4, 3, 4]);
        let m = dag
            .posterior(&placements, LAMBDA, &mut PosteriorScratch::default())
            .expect("a path exists");
        assert_close(m.d_count.iter().sum(), 1.0, "D counts sum to 1");
        for kind in [SegmentKind::V, SegmentKind::J] {
            let total: f64 = placements
                .iter()
                .zip(&m.placement)
                .filter(|(p, _)| p.kind == kind)
                .map(|(_, &q)| q)
                .sum();
            assert!(total <= 1.0 + 1e-12, "{kind:?} marginals sum to {total}");
        }
        // The better V allele dominates.
        assert!(m.placement[0] > 0.9, "V marginal {}", m.placement[0]);
    }
}
