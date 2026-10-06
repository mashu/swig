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
//! **Evidence at the read's own divergence.** `(2, −3)` scores are log-odds
//! only at 11.2 % divergence, so one λ would make unmutated reads
//! underconfident and heavily mutated ones overconfident. The prior terms
//! (`P_D`, pair bonus, N/P) stay at the λ they were derived with, but the
//! alignment evidence is rescored as log-odds at the read's divergence `p̂`
//! ([`divergence`]): a match is worth `ln 4(1 − p̂)` nats and a mismatch
//! `ln(4p̂/3)` ([`EvidenceScale::at_divergence`]). For a match/mismatch
//! scoring that is an affine map of each `Match` column's score,
//! `λ_r·s + μ_r`, which keeps every clip score separable. At `p̂ = 11.2 %` it
//! is the fixed-λ model. The posterior's mode can then differ from the
//! Viterbi call, which stays on the fixed scores; the reported probability
//! is that of the reported call under the per-read model.
//!
//! **Domain.** As in the Viterbi DAG, the first segment keeps its start and
//! the last segment its end: V and J are extended to the read's ends on
//! purpose ([`crate::refine::extend_outer`]). Between neighbours, every pair
//! of cuts `x ≤ y` is summed, including entries beyond the overlap that the
//! Viterbi search skips because they cannot win: they lose nothing for the
//! maximum but carry mass. Each segmentation is counted once by taking
//! *canonical* cuts ([`crate::clip::ClipIndex::matches`]): an entry cut on
//! the first retained base, an exit cut right after the last.
//!
//! **Recursions.** For a node `b` in group `g = (D count, first kind)`, with
//! `F(b, y)` the mass of path prefixes entering `b` at `y` and `ι` the
//! per-base N/P score,
//!
//! ```text
//! X(b, x)  = e^{λ(exit(x) + bonus)} · Σ_{y ≤ x − min} F(b, y) · e^{−λ entry(y)}
//! F(e, y') = Σ_b Σ_{x ≤ y'} X(b, x) · e^{λι(y' − x)}            (+ 1 at a path start)
//! Z        = Σ_b X(b, end of b) · e^{λ end(g, b)}
//! ```
//!
//! and the backward pass mirrors them. Every inner sum is a prefix or suffix
//! sum, as the DAG's maxima are running maxima. On an edge `b → e`, every
//! entry of `e` beyond `b`'s end sees the same exit mass `Σ_x X(b, x)`, so
//! forward adds it once to a difference array and backward reads it from a
//! suffix sum: an edge costs O(overlap), a read
//! `O(Σ_nodes len + Σ_edges (overlap + 1))` multiply–adds.
//!
//! **Arithmetic.** Values are linear, not logarithmic. Each `(node, group)`
//! slot carries one log-scale (as in scaled HMM recursions) and is
//! renormalised once complete. Prior weights `exp(λk)` come from a table and
//! cut weights from running products along the alignment, so no
//! transcendental function is evaluated per cut. Within a placement, weights
//! are centred on the midpoint of its evidence range; a placement whose
//! range exceeds [`MAX_SPREAD`] nats would leave the `f64` range, and the
//! read then gets no marginals.

use super::model::{Grammar, PathModel};
use super::{CallPath, PlacementDag};
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// Largest prefix-score range of one placement, in nats (`λ · score`).
pub const MAX_SPREAD: f64 = 600.0;

/// Scores whose weight `exp(λk)` is tabulated (`|k| ≤ TABLE_RANGE`).
const TABLE_RANGE: Score = 1024;

/// Segment kinds, one array per kind in [`PosteriorScratch::running`] and
/// [`PosteriorScratch::suffix`].
const KINDS: usize = SegmentKind::ALL.len();

/// How [`PlacementDag`] turns scores into probabilities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PosteriorModel {
    /// Nats per unit of the path model's prior terms (`P_D`, pair bonus, N/P).
    pub lambda: f64,
    /// Match and mismatch score of each segment kind's substitution scoring
    /// (indexed by [`SegmentKind::index`]), to rescore its alignments as
    /// log-odds at the read's divergence.
    pub substitution: [(Score, Score); 3],
}

impl PosteriorModel {
    /// Evidence scales of the three segment kinds at divergence `p`.
    pub fn evidence(&self, p: f64) -> [EvidenceScale; 3] {
        self.substitution
            .map(|(matched, mismatch)| EvidenceScale::at_divergence(p, matched, mismatch))
    }
}

/// Nats per unit of alignment evidence: a `Match` column scoring `s` is
/// worth `lambda · s + mu`, a gap scoring `g` is worth `lambda · g`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvidenceScale {
    pub lambda: f64,
    pub mu: f64,
}

impl EvidenceScale {
    /// Scores read as `lambda` nats per unit, with no per-column term.
    pub fn fixed(lambda: f64) -> Self {
        Self { lambda, mu: 0.0 }
    }

    /// Log-odds of "templated at divergence `p`" against "random" (uniform
    /// background) for a match/mismatch scoring: a match is worth
    /// `ln 4(1 − p)` nats and a mismatch `ln(4p/3)`.
    pub fn at_divergence(p: f64, matched: Score, mismatch: Score) -> Self {
        let (a, b) = ((4.0 * (1.0 - p)).ln(), (4.0 * p / 3.0).ln());
        let lambda = (a - b) / f64::from(matched - mismatch);
        Self {
            lambda,
            mu: a - lambda * f64::from(matched),
        }
    }
}

/// Divergence of a read from its called germlines: the fraction of
/// mismatched `Match` columns in the called V and J segments (a D is too
/// short to tell), with half a pseudo-mismatch so that an unmutated read gets
/// a small positive rate. Clamped to `[10⁻⁴, 0.5]`.
pub fn divergence(placements: &[Placement], path: &CallPath) -> f64 {
    let (mut columns, mut mismatches) = (0usize, 0usize);
    for segment in &path.segments {
        let p = &placements[segment.placement];
        if p.kind == SegmentKind::D {
            continue;
        }
        let span = segment.span.q_lo..segment.span.q_hi;
        for (_, entry, exit) in p.clip.matches().filter(|(q, ..)| span.contains(q)) {
            columns += 1;
            mismatches += usize::from(exit < entry);
        }
    }
    ((mismatches as f64 + 0.5) / (columns as f64 + 1.0)).clamp(1e-4, 0.5)
}

/// `exp(±(λd + μ))` for small integer score steps `d`, to walk cut weights
/// along an alignment by multiplication.
#[derive(Clone, Debug)]
struct Steps {
    scale: EvidenceScale,
    /// `exp(λd + μ)` and `exp(−(λd + μ))` for `|d| ≤ STEP_RANGE`.
    up: Vec<f64>,
    down: Vec<f64>,
}

/// Steps between consecutive `Match` columns that are tabulated.
const STEP_RANGE: Score = 24;

impl Steps {
    fn new(scale: EvidenceScale) -> Self {
        let nats = |d: Score| scale.lambda * f64::from(d) + scale.mu;
        Self {
            scale,
            up: (-STEP_RANGE..=STEP_RANGE).map(|d| nats(d).exp()).collect(),
            down: (-STEP_RANGE..=STEP_RANGE)
                .map(|d| (-nats(d)).exp())
                .collect(),
        }
    }

    fn nats(&self, d: Score) -> f64 {
        self.scale.lambda * f64::from(d) + self.scale.mu
    }

    #[inline]
    fn up(&self, d: Score) -> f64 {
        match usize::try_from(d + STEP_RANGE) {
            Ok(i) if i < self.up.len() => self.up[i],
            _ => self.nats(d).exp(),
        }
    }

    #[inline]
    fn down(&self, d: Score) -> f64 {
        match usize::try_from(d + STEP_RANGE) {
            Ok(i) if i < self.down.len() => self.down[i],
            _ => (-self.nats(d)).exp(),
        }
    }
}

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

/// One placement in DAG order (query start, then end).
#[derive(Clone, Copy, Debug)]
struct Node {
    placement: usize,
    kind: SegmentKind,
    q_start: u32,
    /// `q_end − q_start`; cut offsets run over `0..=len`.
    len: u32,
    /// Start of the node's dense cut weights in `entry_w` / `exit_w`.
    cuts: usize,
    /// The node's canonical entry offsets (ascending) in `entries`.
    entries: (usize, usize),
    /// The node's slots in `slots`.
    slots: (usize, usize),
}

impl Node {
    fn q_end(&self) -> u32 {
        self.q_start + self.len
    }

    fn width(&self) -> usize {
        self.len as usize + 1
    }
}

/// The paths through one node with one `(D count, first kind)`: forward
/// mass per entry, backward mass per entry, each `e^{scale} · values`.
#[derive(Clone, Copy, Debug)]
struct Slot {
    d_used: u8,
    first: SegmentKind,
    /// Start of the slot's `width` values in `forward` / `backward`;
    /// `KINDS · values` starts its per-kind arrays in `running` / `suffix`.
    values: usize,
    forward_scale: f64,
    /// Scale of the slot's exit mass; `−∞` when no path leaves it.
    exit_scale: f64,
    backward_scale: f64,
}

impl Slot {
    /// The slot's array for `kind` in `running` / `suffix`.
    fn per_kind(&self, kind: SegmentKind, width: usize) -> std::ops::Range<usize> {
        let start = KINDS * self.values + kind.index() * width;
        start..start + width
    }
}

/// Per-thread buffers for [`PlacementDag::posterior`].
#[derive(Debug, Default)]
pub struct PosteriorScratch {
    weights: Option<Weights>,
    order: Vec<usize>,
    nodes: Vec<Node>,
    /// `e^{λ(centre − entry prefix)}` at canonical entries, else 0.
    entry_w: Vec<f64>,
    /// `e^{λ(exit prefix − centre)}` at canonical exits, else 0.
    exit_w: Vec<f64>,
    entries: Vec<u32>,
    slots: Vec<Slot>,
    forward: Vec<f64>,
    backward: Vec<f64>,
    /// Per slot and successor kind: `R(z) = Σ_{x ≤ z} X(x) · e^{−λιx}`.
    running: Vec<f64>,
    /// Per slot and predecessor kind: `S(o) = Σ_{o' ≥ o} B(o') · e^{λιo'}`.
    suffix: Vec<f64>,
    groups: Vec<(u8, SegmentKind)>,
    work: Vec<f64>,
    tails: Vec<f64>,
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

/// Slots of `predecessors` whose paths continue into the group
/// `(d_used, first)` of a node of kind `kind`.
fn incoming<'s, G: Grammar>(
    s: &'s PosteriorScratch,
    grammar: &'s G,
    predecessors: &'s [Node],
    kind: SegmentKind,
    d_used: u8,
    first: SegmentKind,
) -> impl Iterator<Item = (Node, Slot)> + 's {
    let d_before = d_used.checked_sub(u8::from(kind == SegmentKind::D));
    predecessors.iter().flat_map(move |b| {
        s.slots[b.slots.0..b.slots.1]
            .iter()
            .filter(move |u| {
                u.first == first
                    && Some(u.d_used) == d_before
                    && u.exit_scale > f64::NEG_INFINITY
                    && grammar.can_follow(b.kind, kind, u.d_used)
            })
            .map(move |u| (*b, *u))
    })
}

/// The slot of each node of kind `kind` in `successors` that continues the
/// paths of `slot`, when it has completions.
fn outgoing<'s>(
    s: &'s PosteriorScratch,
    successors: &'s [Node],
    kind: SegmentKind,
    slot: &'s Slot,
) -> impl Iterator<Item = (Node, Slot)> + 's {
    let d_next = slot.d_used + u8::from(kind == SegmentKind::D);
    successors
        .iter()
        .filter(move |e| e.kind == kind)
        .filter_map(move |e| {
            s.slots[e.slots.0..e.slots.1]
                .iter()
                .find(|u| u.d_used == d_next && u.first == slot.first)
                .filter(|u| u.backward_scale > f64::NEG_INFINITY)
                .map(|u| (*e, *u))
        })
}

impl<G: Grammar, M: PathModel> PlacementDag<G, M> {
    /// Posterior marginals of `placements` under `P(path, cuts) ∝ exp(nats)`,
    /// with each segment's alignment evidence on `evidence[kind]` and the
    /// model's prior terms at `lambda` nats per unit.
    ///
    /// `None` when no grammatical path exists, or a placement's evidence
    /// range is too wide for linear arithmetic ([`MAX_SPREAD`]).
    pub fn posterior(
        &self,
        placements: &[Placement],
        lambda: f64,
        evidence: &[EvidenceScale; 3],
        scratch: &mut PosteriorScratch,
    ) -> Option<Marginals> {
        let weights = match scratch.weights.take() {
            Some(w) if w.lambda == lambda => w,
            _ => Weights::new(lambda),
        };
        let marginals = Pass {
            dag: self,
            placements,
            w: &weights,
            steps: evidence.map(Steps::new),
        }
        .run(scratch);
        scratch.weights = Some(weights);
        marginals
    }
}

/// One read's forward–backward.
struct Pass<'a, G, M> {
    dag: &'a PlacementDag<G, M>,
    placements: &'a [Placement],
    /// Prior weights.
    w: &'a Weights,
    /// Evidence steps per segment kind.
    steps: [Steps; 3],
}

impl<G: Grammar, M: PathModel> Pass<'_, G, M> {
    fn run(&self, s: &mut PosteriorScratch) -> Option<Marginals> {
        if self.placements.is_empty() || !self.layout(s) {
            return None;
        }
        s.slots.clear();
        s.forward.clear();
        s.backward.clear();
        s.running.clear();
        s.suffix.clear();
        let mut log_z: Vec<f64> = Vec::new();
        for e in 0..s.nodes.len() {
            self.gather(s, e);
            let (first, end) = s.nodes[e].slots;
            for t in first..end {
                if let Some((d, mass)) = self.emit(s, e, t) {
                    if log_z.len() <= d {
                        log_z.resize(d + 1, f64::NEG_INFINITY);
                    }
                    log_z[d] = log_add(log_z[d], mass);
                }
            }
        }
        let log_partition = log_z.iter().copied().fold(f64::NEG_INFINITY, log_add);
        if !log_partition.is_finite() {
            return None;
        }
        for b in (0..s.nodes.len()).rev() {
            let (first, end) = s.nodes[b].slots;
            for t in first..end {
                self.pull(s, b, t);
            }
        }
        Some(self.collect(s, &log_z, log_partition))
    }

    /// Nodes in DAG order with dense canonical-cut weights. `false` when a
    /// placement's prefix-score range exceeds [`MAX_SPREAD`].
    fn layout(&self, s: &mut PosteriorScratch) -> bool {
        let placements = self.placements;
        s.order.clear();
        s.order.extend(0..placements.len());
        s.order
            .sort_by_key(|&i| (placements[i].q_start(), placements[i].q_end(), i));
        s.nodes.clear();
        s.entry_w.clear();
        s.exit_w.clear();
        s.entries.clear();
        for &i in &s.order {
            let p = &placements[i];
            let steps = &self.steps[p.kind.index()];
            let EvidenceScale { lambda, mu } = steps.scale;
            // Evidence before Match column j (entry) is λ·prefix + μ·j nats,
            // through it (exit) λ·prefix + μ·(j + 1).
            let (lo, hi) = p.clip.matches().enumerate().fold(
                (f64::INFINITY, f64::NEG_INFINITY),
                |(lo, hi), (j, (_, entry, exit))| {
                    let a = lambda * f64::from(entry) + mu * j as f64;
                    let b = lambda * f64::from(exit) + mu * (j + 1) as f64;
                    (lo.min(a.min(b)), hi.max(a.max(b)))
                },
            );
            if lo <= hi && hi - lo > MAX_SPREAD {
                return false;
            }
            let centre = 0.5 * (lo + hi);
            let (q_start, len) = (p.q_start(), p.q_end() - p.q_start());
            let cuts = s.entry_w.len();
            s.entry_w.resize(cuts + len as usize + 1, 0.0);
            s.exit_w.resize(cuts + len as usize + 1, 0.0);
            let first_entry = s.entries.len();
            // Running products: `before` = e^{centre − entry nats},
            // `after` = its inverse, stepped column by column.
            let (mut before, mut after, mut previous) = (0.0, 0.0, None);
            for (q, entry, exit) in p.clip.matches() {
                match previous {
                    None => {
                        let a = lambda * f64::from(entry);
                        before = (centre - a).exp();
                        after = (a - centre).exp();
                    }
                    Some(previous) => {
                        before *= steps.down(entry - previous);
                        after *= steps.up(entry - previous);
                    }
                }
                previous = Some(entry);
                let o = q - q_start;
                s.entry_w[cuts + o as usize] = before;
                s.exit_w[cuts + o as usize + 1] = after * steps.up(exit - entry);
                s.entries.push(o);
            }
            s.nodes.push(Node {
                placement: i,
                kind: p.kind,
                q_start,
                len,
                cuts,
                entries: (first_entry, s.entries.len()),
                slots: (0, 0),
            });
        }
        true
    }

    /// Forward: create node `e`'s slots and their entry masses `F(e, ·)`
    /// from path starts and from every predecessor's exits.
    fn gather(&self, s: &mut PosteriorScratch, e: usize) {
        let grammar = &self.dag.grammar;
        let node = s.nodes[e];
        let width = node.width();
        let predecessors = s.nodes.partition_point(|b| b.q_start < node.q_start);
        let start = grammar
            .can_start(node.kind)
            .then_some((u8::from(node.kind == SegmentKind::D), node.kind))
            .filter(|_| s.entry_w[node.cuts] > 0.0);

        s.groups.clear();
        s.groups.extend(start);
        for b in &s.nodes[..predecessors] {
            for u in &s.slots[b.slots.0..b.slots.1] {
                let group = (u.d_used + u8::from(node.kind == SegmentKind::D), u.first);
                if u.exit_scale > f64::NEG_INFINITY
                    && grammar.can_follow(b.kind, node.kind, u.d_used)
                    && !s.groups.contains(&group)
                {
                    s.groups.push(group);
                }
            }
        }
        let first_slot = s.slots.len();
        for &(d_used, first) in &s.groups {
            let values = s.forward.len();
            s.forward.resize(values + width, 0.0);
            s.backward.resize(values + width, 0.0);
            s.running.resize(KINDS * (values + width), 0.0);
            s.suffix.resize(KINDS * (values + width), 0.0);
            s.slots.push(Slot {
                d_used,
                first,
                values,
                forward_scale: f64::NEG_INFINITY,
                exit_scale: f64::NEG_INFINITY,
                backward_scale: f64::NEG_INFINITY,
            });
        }
        s.nodes[e].slots = (first_slot, s.slots.len());

        for t in first_slot..s.slots.len() {
            let Slot {
                d_used,
                first,
                values,
                ..
            } = s.slots[t];
            let starts_here = start == Some((d_used, first));
            let reference = incoming(
                s,
                grammar,
                &s.nodes[..predecessors],
                node.kind,
                d_used,
                first,
            )
            .map(|(_, u)| u.exit_scale)
            .fold(if starts_here { 0.0 } else { f64::NEG_INFINITY }, f64::max);

            let mut tails = std::mem::take(&mut s.tails);
            tails.clear();
            tails.resize(KINDS * width, 0.0);
            let mut f = std::mem::take(&mut s.work);
            f.clear();
            f.resize(width, 0.0);
            if starts_here {
                f[0] = (-reference).exp();
            }
            let entries = &s.entries[node.entries.0..node.entries.1];
            for (b, u) in incoming(
                s,
                grammar,
                &s.nodes[..predecessors],
                node.kind,
                d_used,
                first,
            ) {
                let factor = (u.exit_scale - reference).exp();
                let iota = self.dag.model.junction_per_base(b.kind, node.kind);
                let r = &s.running[u.per_kind(node.kind, b.width())];
                // Entries inside b: exits up to the entry.
                for &o in entries {
                    let y = node.q_start + o;
                    if y >= b.q_end() {
                        break;
                    }
                    let z = y - b.q_start;
                    f[o as usize] += factor * r[z as usize] * self.w.of(iota * z as Score);
                }
                // Entries at or beyond b's end: all exits, added once.
                let beyond = b.q_end().saturating_sub(node.q_start) as usize;
                if beyond < width {
                    let gap = self.w.of(iota * (node.q_start - b.q_start) as Score);
                    tails[b.kind.index() * width + beyond] += factor * r[b.len as usize] * gap;
                }
            }
            for kind in SegmentKind::ALL {
                let tail = &tails[kind.index() * width..(kind.index() + 1) * width];
                if tail.iter().all(|&v| v == 0.0) {
                    continue;
                }
                let iota = self.dag.model.junction_per_base(kind, node.kind);
                let mut acc = 0.0;
                for ((o, &v), gap) in tail.iter().enumerate().zip(self.geometric(iota)) {
                    acc += v;
                    if s.entry_w[node.cuts + o] > 0.0 {
                        f[o] += acc * gap;
                    }
                }
            }
            let scale = normalise(&mut f);
            s.forward[values..values + width].copy_from_slice(&f);
            s.slots[t].forward_scale = reference + scale;
            s.work = f;
            s.tails = tails;
        }
    }

    /// Forward, continued: the exit masses `X(e, ·)` of slot `t`, their
    /// running sums per successor kind, and the slot's path-end mass
    /// `(D count, ln mass)` if paths may end here.
    fn emit(&self, s: &mut PosteriorScratch, e: usize, t: usize) -> Option<(usize, f64)> {
        let node = s.nodes[e];
        let slot = s.slots[t];
        if slot.forward_scale == f64::NEG_INFINITY {
            return None;
        }
        let width = node.width();
        let min = self.dag.min_bases(node.kind) as usize;
        let bonus = self.w.of(self
            .dag
            .model
            .segment_bonus(&self.placements[node.placement], slot.d_used));
        let mut x = std::mem::take(&mut s.work);
        x.clear();
        let mut acc = 0.0;
        // x[o] holds Σ_{o' ≤ o} F(o') e^{−λ entry} until it is overwritten
        // from the right with the exit masses.
        for o in 0..width {
            acc += s.forward[slot.values + o] * s.entry_w[node.cuts + o];
            x.push(acc);
        }
        for o in (0..width).rev() {
            let exit = s.exit_w[node.cuts + o];
            x[o] = if o >= min && exit > 0.0 {
                x[o - min] * exit * bonus
            } else {
                0.0
            };
        }
        let exit_scale = slot.forward_scale + normalise(&mut x);
        s.slots[t].exit_scale = exit_scale;
        let mut end = None;
        if exit_scale > f64::NEG_INFINITY {
            let mut done: Option<(Score, usize)> = None;
            for kind in SegmentKind::ALL {
                if !self.dag.grammar.can_follow(node.kind, kind, slot.d_used) {
                    continue;
                }
                let iota = self.dag.model.junction_per_base(node.kind, kind);
                let range = slot.per_kind(kind, width);
                match done {
                    Some((same, from)) if same == iota => {
                        s.running.copy_within(from..from + width, range.start)
                    }
                    _ => {
                        let mut acc = 0.0;
                        for ((r, &v), gap) in s.running[range.clone()]
                            .iter_mut()
                            .zip(&x)
                            .zip(self.geometric(-iota))
                        {
                            acc += v * gap;
                            *r = acc;
                        }
                        done = Some((iota, range.start));
                    }
                }
            }
            if self.dag.grammar.can_end(slot.first, node.kind) && x[width - 1] > 0.0 {
                let bonus = self
                    .w
                    .of(self.dag.model.end(slot.first, node.kind, slot.d_used));
                end = Some((
                    usize::from(slot.d_used),
                    exit_scale + (x[width - 1] * bonus).ln(),
                ));
            }
        }
        s.work = x;
        end
    }

    /// Backward: the completion masses `B(b, ·)` of slot `t` from path ends
    /// and every successor's entries, then its suffix sums per predecessor kind.
    fn pull(&self, s: &mut PosteriorScratch, b: usize, t: usize) {
        let (grammar, model) = (&self.dag.grammar, &self.dag.model);
        let node = s.nodes[b];
        let slot = s.slots[t];
        if slot.exit_scale == f64::NEG_INFINITY {
            return;
        }
        let width = node.width();
        let successors = s.nodes.partition_point(|e| e.q_start <= node.q_start);
        let can_end = grammar.can_end(slot.first, node.kind);
        let follows = |kind: SegmentKind| grammar.can_follow(node.kind, kind, slot.d_used);
        let mut reference = if can_end { 0.0 } else { f64::NEG_INFINITY };
        for kind in SegmentKind::ALL.into_iter().filter(|&k| follows(k)) {
            for (_, u) in outgoing(s, &s.nodes[successors..], kind, &slot) {
                reference = reference.max(u.backward_scale);
            }
        }
        if reference == f64::NEG_INFINITY {
            return;
        }

        // inner[o]: completions after an exit at offset o.
        let mut inner = std::mem::take(&mut s.work);
        inner.clear();
        inner.resize(width, 0.0);
        if can_end {
            inner[width - 1] =
                self.w.of(model.end(slot.first, node.kind, slot.d_used)) * (-reference).exp();
        }
        let mut tails = std::mem::take(&mut s.tails);
        for kind in SegmentKind::ALL.into_iter().filter(|&k| follows(k)) {
            let iota = model.junction_per_base(node.kind, kind);
            tails.clear();
            tails.resize(width, 0.0);
            for (e, u) in outgoing(s, &s.nodes[successors..], kind, &slot) {
                let factor = (u.backward_scale - reference).exp();
                let backward = &s.backward[u.values..u.values + e.width()];
                for &o in &s.entries[e.entries.0..e.entries.1] {
                    let y = e.q_start + o;
                    if y >= node.q_end() {
                        break;
                    }
                    let z = y - node.q_start;
                    tails[z as usize] +=
                        factor * backward[o as usize] * self.w.of(iota * z as Score);
                }
                let beyond = node.q_end().saturating_sub(e.q_start) as usize;
                if beyond < e.width() {
                    let gap = self.w.of(iota * (e.q_start - node.q_start) as Score);
                    let suffix = &s.suffix[u.per_kind(node.kind, e.width())];
                    tails[width - 1] += factor * gap * suffix[beyond];
                }
            }
            let mut acc = 0.0;
            for v in tails.iter_mut().rev() {
                acc += *v;
                *v = acc;
            }
            for ((i, &t), gap) in inner.iter_mut().zip(&tails).zip(self.geometric(-iota)) {
                *i += gap * t;
            }
        }
        s.tails = tails;

        // B(y) = e^{−λ entry(y)} · Σ_{x ≥ y + min} X-weight(x) · inner(x).
        let min = self.dag.min_bases(node.kind) as usize;
        let bonus = self
            .w
            .of(model.segment_bonus(&self.placements[node.placement], slot.d_used));
        let mut acc = 0.0;
        for o in (0..width).rev() {
            acc += s.exit_w[node.cuts + o] * bonus * inner[o];
            inner[o] = acc;
        }
        let out = &mut s.backward[slot.values..slot.values + width];
        for (o, value) in out.iter_mut().enumerate() {
            let entry = s.entry_w[node.cuts + o];
            *value = if entry > 0.0 && o + min < width {
                entry * inner[o + min]
            } else {
                0.0
            };
        }
        let scale = normalise(out);
        s.slots[t].backward_scale = reference + scale;
        s.work = inner;

        // Suffix sums for predecessors of each kind that can precede this slot.
        let mut done: Option<(Score, usize)> = None;
        for kind in SegmentKind::ALL {
            if !(0..=slot.d_used).any(|d| grammar.can_follow(kind, node.kind, d)) {
                continue;
            }
            let iota = model.junction_per_base(kind, node.kind);
            let range = slot.per_kind(kind, width);
            match done {
                Some((same, from)) if same == iota => {
                    s.suffix.copy_within(from..from + width, range.start)
                }
                _ => {
                    let backward = &s.backward[slot.values..slot.values + width];
                    let suffix = &mut s.suffix[range.clone()];
                    for ((v, &b), gap) in suffix.iter_mut().zip(backward).zip(self.geometric(iota))
                    {
                        *v = b * gap;
                    }
                    let mut acc = 0.0;
                    for v in suffix.iter_mut().rev() {
                        acc += *v;
                        *v = acc;
                    }
                    done = Some((iota, range.start));
                }
            }
        }
    }

    /// `e^{λιo}` for `o = 0, 1, …` as a running product: exactly 1 for `ι = 0`.
    fn geometric(&self, iota: Score) -> impl Iterator<Item = f64> {
        let step = self.w.of(iota);
        std::iter::successors(Some(1.0), move |g| Some(g * step))
    }

    fn collect(&self, s: &PosteriorScratch, log_z: &[f64], log_partition: f64) -> Marginals {
        let mut placement = vec![0.0; self.placements.len()];
        for node in &s.nodes {
            for slot in &s.slots[node.slots.0..node.slots.1] {
                let scale = slot.forward_scale + slot.backward_scale;
                if scale == f64::NEG_INFINITY {
                    continue;
                }
                let range = slot.values..slot.values + node.width();
                let mass: f64 = s.forward[range.clone()]
                    .iter()
                    .zip(&s.backward[range])
                    .map(|(f, b)| f * b)
                    .sum();
                if mass > 0.0 {
                    placement[node.placement] += (scale + mass.ln() - log_partition).exp();
                }
            }
        }
        Marginals {
            placement,
            d_count: log_z.iter().map(|z| (z - log_partition).exp()).collect(),
            log_partition,
        }
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
        evidence: EvidenceScale,
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

        /// Nats of `b`'s evidence on `[entry, exit)`.
        fn evidence(&self, b: &Placement, entry: u32, exit: u32) -> f64 {
            let span = b.clip.span(entry, exit).expect("canonical cuts");
            let columns = b
                .clip
                .matches()
                .filter(|(q, ..)| (entry..exit).contains(q))
                .count();
            self.evidence.lambda * f64::from(span.score) + self.evidence.mu * columns as f64
        }

        /// `acc`: path score so far; `nats`: its log weight.
        #[allow(clippy::too_many_arguments)]
        fn visit(
            &mut self,
            path: &mut Vec<usize>,
            entry: u32,
            d_used: u8,
            first: SegmentKind,
            acc: Score,
            nats: f64,
        ) {
            let (dag, placements) = (self.dag, self.placements);
            let b = &placements[*path.last().expect("non-empty")];
            let min_b = dag.min_retained[b.kind.index()].max(1);
            for x in Self::exits(b).into_iter().filter(|&x| x >= entry + min_b) {
                let span = b.clip.span(entry, x).expect("canonical cuts");
                let bonus = dag.model.segment_bonus(b, d_used);
                let segment = span.score + bonus;
                let segment_nats = self.evidence(b, entry, x) + LAMBDA * f64::from(bonus);
                if dag.grammar.can_end(first, b.kind) && x == b.q_end() {
                    let end = dag.model.end(first, b.kind, d_used);
                    let total = acc + segment + end;
                    let weight = (nats + segment_nats + LAMBDA * f64::from(end)).exp();
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
                        let junction = iota * (y - x) as Score;
                        path.push(j);
                        self.visit(
                            path,
                            y,
                            d_next,
                            first,
                            acc + segment + junction,
                            nats + segment_nats + LAMBDA * f64::from(junction),
                        );
                        path.pop();
                    }
                }
            }
        }

        fn run(
            dag: &'a PlacementDag,
            placements: &'a [Placement],
            evidence: EvidenceScale,
        ) -> Self {
            let mut brute = Brute {
                dag,
                placements,
                evidence,
                z: 0.0,
                best: Score::MIN,
                placement: vec![0.0; placements.len()],
                d_count: vec![0.0; 3],
            };
            for (s, p) in placements.iter().enumerate() {
                if dag.grammar.can_start(p.kind) {
                    let d0 = u8::from(p.kind == SegmentKind::D);
                    brute.visit(&mut vec![s], p.q_start(), d0, p.kind, 0, 0.0);
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

    fn check(dag: &PlacementDag, placements: &[Placement], evidence: EvidenceScale) {
        let brute = Brute::run(dag, placements, evidence);
        let m = dag
            .posterior(
                placements,
                LAMBDA,
                &[evidence; 3],
                &mut PosteriorScratch::default(),
            )
            .expect("a path exists");
        assert_close(m.log_partition, brute.z.ln(), "ln Z");
        for (i, (&got, &want)) in m.placement.iter().zip(&brute.placement).enumerate() {
            assert_close(got, want / brute.z, &format!("placement {i}"));
        }
        for (k, &want) in brute.d_count.iter().enumerate() {
            let got = m.d_count.get(k).copied().unwrap_or(0.0);
            assert_close(got, want / brute.z, &format!("D count {k}"));
        }
        // The Viterbi path is the best path under the same scores.
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
        let fixed = EvidenceScale::fixed(LAMBDA);
        for grammar in [VdjGrammar::VDJ, VdjGrammar::TANDEM_D, VdjGrammar::VJ] {
            check(
                &PlacementDag::new(grammar, model, [4, 3, 4]),
                &placements,
                fixed,
            );
        }
        let per_base = EngineeringObjective {
            insertion_per_base: -1,
            ..model
        };
        let dag = PlacementDag::new(VdjGrammar::TANDEM_D, per_base, [4, 3, 4]);
        check(&dag, &placements, fixed);
        // Evidence at a read's own divergence: a per-column term μ ≠ 0.
        check(&dag, &placements, EvidenceScale::at_divergence(0.02, 2, -3));
    }

    #[test]
    fn evidence_scale_is_log_odds_at_the_divergence() {
        let at = |p: f64| EvidenceScale::at_divergence(p, 2, -3);
        let nats = |e: EvidenceScale, s: Score| e.lambda * f64::from(s) + e.mu;
        for p in [0.001, 0.05, 0.2] {
            assert_close(nats(at(p), 2), (4.0 * (1.0 - p)).ln(), "match");
            assert_close(nats(at(p), -3), (4.0 * p / 3.0).ln(), "mismatch");
        }
        // (2, −3) is log-odds at 11.2 % divergence: the fixed-λ model.
        let base = at(0.112);
        assert!(
            (base.lambda - 0.634).abs() < 1e-3 && base.mu.abs() < 1e-3,
            "{base:?}"
        );
    }

    #[test]
    fn marginals_are_probabilities() {
        let (_, placements) = junction();
        let dag = PlacementDag::new(VdjGrammar::VDJ, EngineeringObjective::default(), [4, 3, 4]);
        let evidence = [EvidenceScale::fixed(LAMBDA); 3];
        let m = dag
            .posterior(
                &placements,
                LAMBDA,
                &evidence,
                &mut PosteriorScratch::default(),
            )
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
