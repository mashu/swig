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
//! every pair of cuts between neighbours. With the evidence read at the fixed
//! λ, the Viterbi path of [`PlacementDag`] is its mode (both run over the same
//! domain). Summing instead of maximising gives the marginal probability of
//! each placement (an allele call) and of each D count (D presence).
//!
//! **Evidence at the read's own divergence.** `(2, −3)` scores are log-odds
//! only at 11.2 % divergence, so one λ would make unmutated reads
//! underconfident and heavily mutated ones overconfident. The alignment
//! evidence is therefore read as log-odds at the read's divergence `p̂`, and
//! the prior terms at the λ they were derived with ([`super::units`]). The
//! weights come from a [`PathNats`]: the call's own [`ReadModel`] in nats
//! under [`CallModel::PerRead`], whose Viterbi path is then a mode of this
//! distribution; under [`CallModel::FixedScale`] the evidence unquantised and
//! the presence costs `λ·P*_k` ([`super::EvidenceModel::nats`]), while the call stays
//! on the fixed scores, so the reported probability is that of a call made
//! under another model.
//!
//! **Domain.** As in the Viterbi DAG, the first segment keeps its start and
//! the last segment its end: V and J are extended to the read's ends on
//! purpose ([`crate::flank::extend_outer`]). Between neighbours, every pair
//! of cuts `x ≤ y` is summed, the domain the Viterbi DAG maximises over.
//! Each segmentation is counted once by taking *canonical* cuts
//! ([`crate::clip::ClipIndex::matches`]): an entry cut on the first retained
//! base, an exit cut right after the last.
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
//! transcendental function is evaluated per cut. Within a placement, entry
//! and exit weights are centred separately, each on the midpoint of the
//! evidence range over the cuts a path can use on that side, and the
//! difference is carried as a per-node `shift`. Junction recursions keep only
//! relative factors (`≤ 1`) and renormalise before multiplying weights; the
//! N/P cost of a gap between non-overlapping neighbours joins the slot's
//! log-scale. A placement whose range on one side exceeds [`MAX_SPREAD`]
//! nats, or an overlap across which both neighbours gain more than
//! [`OVERLAP_SPREAD`], would leave the `f64` range, and the read then gets no
//! marginals.

use std::cell::Cell;

use super::PlacementDag;
use super::lattice::{Group, KINDS, Lattice, Node, column_steps};
use super::model::{Grammar, PathModel};
#[cfg(doc)]
use super::units::{CallModel, ReadModel};
use super::units::{EvidenceScale, PathNats};
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// Largest range of evidence, in nats, over a placement's usable entry cuts
/// or over its usable exit cuts. Each weight then stays within `e^{±600}`.
pub const MAX_SPREAD: f64 = 1200.0;

/// Largest evidence a predecessor's exits and a successor's entries may
/// both gain across their overlap. A junction product `F(y) · B(y)` there
/// is about `e^{−min of the two}` relative to its factors' maxima, so this
/// keeps it inside the `f64` range. Real overlaps (a junction's length) gain
/// a few dozen nats.
pub const OVERLAP_SPREAD: f64 = 600.0;

/// D counts the grammar is asked about when deciding which cuts are usable.
const MAX_D: u8 = 3;

/// Scores whose weight `exp(λk)` is tabulated (`|k| ≤ TABLE_RANGE`).
const TABLE_RANGE: Score = 1024;

/// `exp(±(λd + μ))` for small integer score steps `d`, to walk cut weights
/// along an alignment by multiplication.
#[derive(Clone, Debug)]
struct Steps {
    scale: EvidenceScale,
    /// `exp(λd + μ)` and `exp(−(λd + μ))` for `|d| ≤ STEP_RANGE`, each
    /// computed when first needed (NaN until then): a read's columns score
    /// a match or a mismatch, so a few of them are.
    up: [Cell<f64>; STEPS],
    down: [Cell<f64>; STEPS],
}

/// Steps between consecutive `Match` columns that are tabulated.
const STEP_RANGE: Score = 24;

/// Tabulated steps, `−STEP_RANGE..=STEP_RANGE`.
const STEPS: usize = 2 * STEP_RANGE as usize + 1;

/// The value of `cell`, computing it with `value` the first time.
#[inline]
fn tabulated(cell: &Cell<f64>, value: impl FnOnce() -> f64) -> f64 {
    let known = cell.get();
    if known.is_nan() {
        let value = value();
        cell.set(value);
        value
    } else {
        known
    }
}

impl Steps {
    fn new(scale: EvidenceScale) -> Self {
        Self {
            scale,
            up: std::array::from_fn(|_| Cell::new(f64::NAN)),
            down: std::array::from_fn(|_| Cell::new(f64::NAN)),
        }
    }

    fn nats(&self, d: Score) -> f64 {
        self.scale.lambda * f64::from(d) + self.scale.mu
    }

    #[inline]
    fn up(&self, d: Score) -> f64 {
        match usize::try_from(d + STEP_RANGE) {
            Ok(i) if i < STEPS => tabulated(&self.up[i], || self.nats(d).exp()),
            _ => self.nats(d).exp(),
        }
    }

    #[inline]
    fn down(&self, d: Score) -> f64 {
        match usize::try_from(d + STEP_RANGE) {
            Ok(i) if i < STEPS => tabulated(&self.down[i], || (-self.nats(d)).exp()),
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
    /// Probability that the path contains at least one D: `+0.0` when no
    /// path can hold one (an empty `f64` sum is `−0.0`). Summed over the D
    /// counts, not `1 − P(no D)`, which loses a small probability.
    pub fn d_present(&self) -> f64 {
        self.d_count.iter().skip(1).fold(0.0, |sum, p| sum + p)
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

/// Log-scales of one slot: its values are `e^{scale} · values`.
#[derive(Clone, Copy, Debug)]
struct Scales {
    forward: f64,
    /// Scale of the slot's exit mass; `−∞` when no path leaves it.
    exit: f64,
    backward: f64,
}

impl Scales {
    const EMPTY: Self = Self {
        forward: f64::NEG_INFINITY,
        exit: f64::NEG_INFINITY,
        backward: f64::NEG_INFINITY,
    };
}

/// Per-thread buffers for [`PlacementDag::posterior`].
#[derive(Debug, Default)]
pub struct PosteriorScratch {
    weights: Option<Weights>,
    lattice: Lattice,
    /// Per node: `c_out − c_in`. Entry and exit weights are centred
    /// separately, so a span's weight is `entry_w · exit_w · e^{shift}`.
    shift: Vec<f64>,
    /// Per slot.
    scales: Vec<Scales>,
    /// `e^{λ(centre − entry prefix)}` at canonical entries, else 0.
    entry_w: Vec<f64>,
    /// `e^{λ(exit prefix − centre)}` at canonical exits, else 0.
    exit_w: Vec<f64>,
    /// Per slot value: forward mass per entry, backward mass per entry.
    forward: Vec<f64>,
    backward: Vec<f64>,
    /// Per slot and successor kind: `R(z) = Σ_{x ≤ z} X(x) · e^{−λιx}`.
    running: Vec<f64>,
    /// Per slot and predecessor kind: `S(o) = Σ_{o' ≥ o} B(o') · e^{λιo'}`.
    suffix: Vec<f64>,
    work: Vec<f64>,
    tails: Vec<f64>,
    /// Per Match column: query position, the gap run before it, its
    /// substitution score ([`column_steps`]), and the evidence (nats) before and
    /// through it.
    columns: Vec<(u32, Score, Score, f64, f64)>,
    /// Uncentred evidence at each usable entry / exit offset (NaN elsewhere).
    entry_nats: Vec<f64>,
    exit_nats: Vec<f64>,
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

/// Slots of the `successors` of kind `kind` that continue `group` and have
/// completions.
fn outgoing<'a>(
    lattice: &'a Lattice,
    scales: &'a [Scales],
    successors: &'a [Node],
    kind: SegmentKind,
    group: Group,
) -> impl Iterator<Item = (Node, usize)> + 'a {
    let next = group.then(kind);
    successors
        .iter()
        .filter(move |e| e.kind == kind)
        .filter_map(move |e| {
            lattice
                .slot_of(e, next)
                .filter(|&u| scales[u].backward > f64::NEG_INFINITY)
                .map(|u| (*e, u))
        })
}

impl<G: Grammar, M: PathModel> PlacementDag<G, M> {
    /// Posterior marginals of `placements` under `P(path, cuts) ∝ exp(nats)`,
    /// with each path term read by `nats`.
    ///
    /// `None` when no grammatical path exists, or when linear arithmetic
    /// cannot hold the read: a placement's evidence range ([`MAX_SPREAD`]) or
    /// an overlap's ([`OVERLAP_SPREAD`]) is too wide.
    pub fn posterior(
        &self,
        placements: &[Placement],
        nats: &PathNats,
        scratch: &mut PosteriorScratch,
    ) -> Option<Marginals> {
        let weights = match scratch.weights.take() {
            Some(w) if w.lambda == nats.lambda => w,
            _ => Weights::new(nats.lambda),
        };
        let marginals = Pass {
            dag: self,
            placements,
            w: &weights,
            presence: nats.presence.map(|cost| (-cost).exp()),
            steps: nats.evidence.map(Steps::new),
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
    /// `e^{−c_k}`, the weight of a segment's presence cost, per kind.
    presence: [f64; 3],
    /// Evidence steps per rearranged segment kind.
    steps: [Steps; 3],
}

impl<G: Grammar, M: PathModel> Pass<'_, G, M> {
    /// Weight of `node`'s cut-independent term, its kind's presence cost.
    fn presence(&self, node: Node) -> f64 {
        self.presence[node.kind.index()]
    }

    fn run(&self, s: &mut PosteriorScratch) -> Option<Marginals> {
        if self.placements.is_empty() || !self.layout(s) {
            return None;
        }
        s.scales.clear();
        s.forward.clear();
        s.backward.clear();
        // `running` and `suffix` are not cleared: a slot's sums are written
        // before any other slot reads them, `running` in `emit` for every
        // kind that can follow it and `suffix` in `pull` for every kind that
        // can precede it, by the predicates that record the edges reading
        // them.
        let mut log_z: Vec<f64> = Vec::new();
        for e in 0..s.lattice.nodes.len() {
            self.gather(s, e);
            for t in s.lattice.nodes[e].slots() {
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
        for b in (0..s.lattice.nodes.len()).rev() {
            for t in s.lattice.nodes[b].slots() {
                self.pull(s, b, t);
            }
        }
        Some(self.collect(s, &log_z, log_partition))
    }

    /// The lattice with dense canonical-cut weights. `false` when a
    /// placement's usable entries or exits span more than [`MAX_SPREAD`] nats,
    /// or an overlap more than [`OVERLAP_SPREAD`].
    fn layout(&self, s: &mut PosteriorScratch) -> bool {
        let placements = self.placements;
        let grammar = &self.dag.grammar;
        let lambda0 = self.w.lambda;
        s.lattice.layout(placements, |_| true);
        let cells = s.lattice.cells();
        for (weights, fill) in [(&mut s.entry_w, 0.0), (&mut s.exit_w, 0.0)]
            .into_iter()
            .chain([(&mut s.entry_nats, f64::NAN), (&mut s.exit_nats, f64::NAN)])
        {
            weights.clear();
            weights.resize(cells, fill);
        }
        s.shift.clear();
        let any_d = |f: &dyn Fn(u8) -> bool| (0..=MAX_D).any(f);
        for node in &s.lattice.nodes {
            let p = &placements[node.placement];
            let steps = &self.steps[p.kind.index()];
            let EvidenceScale { lambda, mu } = steps.scale;
            let (q_start, len, cuts) = (node.q_start, node.len, node.cuts);
            // Cuts a path can use: entries past the start only when a kind
            // can precede this one, exits before the end only when one can
            // follow it. Only those set the centres and the spread.
            let starts = grammar.can_start(p.kind);
            let ends = SegmentKind::REARRANGED
                .iter()
                .any(|&f| any_d(&|d| grammar.can_end(f, p.kind, d)));
            let preceded = SegmentKind::REARRANGED
                .iter()
                .any(|&k| any_d(&|d| grammar.can_follow(k, p.kind, d)));
            let followed = SegmentKind::REARRANGED
                .iter()
                .any(|&k| any_d(&|d| grammar.can_follow(p.kind, k, d)));
            let entry_used = |o: u32| preceded || (starts && o == 0);
            let exit_used = |o: u32| followed || (ends && o == len);

            // Evidence nats: a substitution column scoring s adds λ_r·s + μ_r,
            // an uninformative one (an N, s = 0) nothing; gaps keep λ.
            s.columns.clear();
            let mut previous: Option<f64> = None;
            let (mut in_lo, mut in_hi) = (f64::INFINITY, f64::NEG_INFINITY);
            let (mut out_lo, mut out_hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for (q, gap, column) in column_steps(p) {
                let gap_nats = lambda0 * f64::from(gap);
                let before = previous.map_or(gap_nats, |nats| nats + gap_nats);
                let through = before
                    + if column == 0 {
                        0.0
                    } else {
                        lambda * f64::from(column) + mu
                    };
                if entry_used(q - q_start) {
                    (in_lo, in_hi) = (in_lo.min(before), in_hi.max(before));
                }
                if exit_used(q + 1 - q_start) {
                    (out_lo, out_hi) = (out_lo.min(through), out_hi.max(through));
                }
                s.columns.push((q, gap, column, before, through));
                previous = Some(through);
            }
            if in_hi - in_lo > MAX_SPREAD || out_hi - out_lo > MAX_SPREAD {
                return false;
            }
            let centre = |lo: f64, hi: f64| if lo <= hi { 0.5 * (lo + hi) } else { 0.0 };
            let (c_in, c_out) = (centre(in_lo, in_hi), centre(out_lo, out_hi));

            // Running products along the alignment, `before` = e^{c_in − entry
            // nats} and `inv` = e^{entry nats − c_out}, reseeded exactly
            // wherever a used cut finds them out of range.
            let substitution = |column: Score, up: bool| match (column, up) {
                (0, _) => 1.0,
                (c, true) => steps.up(c),
                (c, false) => steps.down(c),
            };
            let (mut before, mut inv) = (f64::NAN, f64::NAN);
            let mut last: Option<Score> = None;
            for &(q, gap, column, nats, through) in &s.columns {
                if let Some(previous) = last {
                    before *= substitution(previous, false) * self.w.of(-gap);
                    inv *= substitution(previous, true) * self.w.of(gap);
                }
                last = Some(column);
                let o = q - q_start;
                let usable = |v: f64| v.is_finite() && v > 0.0;
                if entry_used(o) {
                    if !usable(before) {
                        before = (c_in - nats).exp();
                    }
                    s.entry_w[cuts + o as usize] = before;
                    s.entry_nats[cuts + o as usize] = nats;
                }
                if exit_used(o + 1) {
                    if !usable(inv) {
                        inv = (nats - c_out).exp();
                    }
                    s.exit_w[cuts + o as usize + 1] = inv * substitution(column, true);
                    s.exit_nats[cuts + o as usize + 1] = through;
                }
            }
            s.shift.push(c_out - c_in);
        }
        self.overlaps_in_range(s)
    }

    /// Whether every overlap between a node and a possible successor keeps
    /// its junction products in range ([`OVERLAP_SPREAD`]).
    fn overlaps_in_range(&self, s: &PosteriorScratch) -> bool {
        let grammar = &self.dag.grammar;
        // Range of the finite values of `nats` at query positions `lo..=hi`.
        let range = |nats: &[f64], node: &Node, lo: u32, hi: u32| {
            let at = |q: u32| nats[node.cuts + (q - node.q_start) as usize];
            let (min, max) = (lo..=hi)
                .map(at)
                .filter(|v| v.is_finite())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
                    (a.min(v), b.max(v))
                });
            (max - min).max(0.0)
        };
        let nodes = &s.lattice.nodes;
        for (k, e) in nodes.iter().enumerate() {
            for b in &nodes[..k] {
                let follows = (0..=MAX_D).any(|d| grammar.can_follow(b.kind, e.kind, d));
                if !(b.precedes(e) && follows && b.q_end() > e.q_start && e.q_end() > b.q_start) {
                    continue;
                }
                let (lo, hi) = (b.q_start.max(e.q_start), b.q_end().min(e.q_end()));
                let gained = range(&s.exit_nats, b, lo, hi).min(range(&s.entry_nats, e, lo, hi));
                if gained > OVERLAP_SPREAD {
                    return false;
                }
            }
        }
        true
    }

    /// Forward: open node `e`'s slots and compute their entry masses
    /// `F(e, ·)` from path starts and from every predecessor's exits.
    fn gather(&self, s: &mut PosteriorScratch, e: usize) {
        let grammar = &self.dag.grammar;
        let node = s.lattice.nodes[e];
        let width = node.width();
        let start = grammar
            .can_start(node.kind)
            .then(|| Group::start(node.kind))
            .filter(|_| s.entry_w[node.cuts] > 0.0);
        let scales = &s.scales;
        let slots = s
            .lattice
            .open(e, grammar, start, |u| scales[u].exit > f64::NEG_INFINITY);
        s.scales.resize(s.lattice.slots.len(), Scales::EMPTY);
        let values = s.lattice.values();
        s.forward.resize(values, 0.0);
        s.backward.resize(values, 0.0);
        // Grown, never cleared (see `run`).
        if s.running.len() < KINDS * values {
            s.running.resize(KINDS * values, 0.0);
            s.suffix.resize(KINDS * values, 0.0);
        }

        // A predecessor ending before this node reaches it across an N/P
        // gap; its weight joins the log scale, so a long, costly gap cannot
        // underflow when every path crosses it. 0 for an overlapping one.
        let gap = |b: &Node| {
            let iota = self.dag.model.junction_per_base(b.kind, node.kind);
            self.w.lambda * f64::from(iota * node.q_start.saturating_sub(b.q_end()) as Score)
        };
        if s.tails.len() < KINDS * width {
            s.tails.resize(KINDS * width, 0.0);
        }
        for t in slots {
            let slot = s.lattice.slots[t];
            let starts_here = start == Some(slot.group);
            let scales = &s.scales;
            let reference = s
                .lattice
                .incoming(t)
                .map(|(b, u)| scales[u].exit + gap(&b))
                .fold(if starts_here { 0.0 } else { f64::NEG_INFINITY }, f64::max);

            let mut tails = std::mem::take(&mut s.tails);
            // Per predecessor kind, cleared when an edge first reaches it.
            let mut touched = [false; KINDS];
            let mut f = std::mem::take(&mut s.work);
            f.clear();
            f.resize(width, 0.0);
            if starts_here {
                f[0] = (-reference).exp();
            }
            for (b, u) in s.lattice.incoming(t) {
                let factor = (scales[u].exit + gap(&b) - reference).exp();
                let r = &s.running[s.lattice.slots[u].per_kind(node.kind, b.width())];
                // Entries inside b: exits up to the entry. r[z] already
                // carries e^{λι(z − x)} for each exit x.
                for (o, z) in s.lattice.inside(&b, &node) {
                    f[o] += factor * r[z];
                }
                // Entries at or beyond b's end: all exits, added once.
                if let Some(beyond) = b.beyond(&node) {
                    let k = b.kind.index();
                    let tail = &mut tails[k * width..(k + 1) * width];
                    if !std::mem::replace(&mut touched[k], true) {
                        tail.fill(0.0);
                    }
                    tail[beyond] += factor * r[b.len as usize];
                }
            }
            for kind in SegmentKind::REARRANGED {
                if !touched[kind.index()] {
                    continue;
                }
                let tail = &tails[kind.index() * width..(kind.index() + 1) * width];
                // acc(o) = Σ_{o' ≤ o} tail(o') · e^{λι(o − o')}.
                let step = self.w.of(self.dag.model.junction_per_base(kind, node.kind));
                let mut acc = 0.0;
                for (o, &v) in tail.iter().enumerate() {
                    acc = acc * step + v;
                    if s.entry_w[node.cuts + o] > 0.0 {
                        f[o] += acc;
                    }
                }
            }
            let scale = normalise(&mut f);
            s.forward[slot.values..slot.values + width].copy_from_slice(&f);
            s.scales[t].forward = reference + scale;
            s.work = f;
            s.tails = tails;
        }
    }

    /// Forward, continued: the exit masses `X(e, ·)` of slot `t`, their
    /// running sums per successor kind, and the slot's path-end mass
    /// `(D count, ln mass)` if paths may end here.
    fn emit(&self, s: &mut PosteriorScratch, e: usize, t: usize) -> Option<(usize, f64)> {
        let node = s.lattice.nodes[e];
        let slot = s.lattice.slots[t];
        let Group { d_used, first } = slot.group;
        let forward_scale = s.scales[t].forward;
        if forward_scale == f64::NEG_INFINITY {
            return None;
        }
        let width = node.width();
        let min = self.dag.min_bases(node.kind) as usize;
        let bonus = self.presence(node);
        let mut x = std::mem::take(&mut s.work);
        x.clear();
        let mut acc = 0.0;
        // x[o] holds Σ_{o' ≤ o} F(o') e^{−λ entry} (normalised) until it is
        // overwritten from the right with the exit masses.
        for o in 0..width {
            acc += s.forward[slot.values + o] * s.entry_w[node.cuts + o];
            x.push(acc);
        }
        let entered = normalise(&mut x);
        for o in (0..width).rev() {
            let exit = s.exit_w[node.cuts + o];
            x[o] = if o >= min && exit > 0.0 {
                x[o - min] * exit * bonus
            } else {
                0.0
            };
        }
        let exit_scale = forward_scale + entered + normalise(&mut x) + s.shift[e];
        s.scales[t].exit = exit_scale;
        let mut end = None;
        if exit_scale > f64::NEG_INFINITY {
            let mut done: Option<(Score, usize)> = None;
            for kind in SegmentKind::REARRANGED {
                if !self.dag.grammar.can_follow(node.kind, kind, d_used) {
                    continue;
                }
                let iota = self.dag.model.junction_per_base(node.kind, kind);
                let range = slot.per_kind(kind, width);
                match done {
                    Some((same, from)) if same == iota => {
                        s.running.copy_within(from..from + width, range.start)
                    }
                    _ => {
                        // R(z) = Σ_{x ≤ z} X(x) · e^{λι(z − x)}: the N/P weight
                        // is carried relative to z, so it never exceeds 1.
                        let step = self.w.of(iota);
                        let mut acc = 0.0;
                        for (r, &v) in s.running[range.clone()].iter_mut().zip(&x) {
                            acc = acc * step + v;
                            *r = acc;
                        }
                        done = Some((iota, range.start));
                    }
                }
            }
            if self.dag.grammar.can_end(first, node.kind, d_used) && x[width - 1] > 0.0 {
                let bonus = self.w.of(self.dag.model.end(first, node.kind, d_used));
                end = Some((
                    usize::from(d_used),
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
        let node = s.lattice.nodes[b];
        let slot = s.lattice.slots[t];
        let Group { d_used, first } = slot.group;
        if s.scales[t].exit == f64::NEG_INFINITY {
            return;
        }
        let width = node.width();
        let successors = &s.lattice.nodes[s.lattice.successors(b)..];
        let can_end = grammar.can_end(first, node.kind, d_used);
        let follows = |kind: SegmentKind| grammar.can_follow(node.kind, kind, d_used);
        // As in the forward pass, an N/P gap to a successor that starts after
        // this node's end joins the log scale.
        let gap = |e: &Node| {
            let iota = model.junction_per_base(node.kind, e.kind);
            self.w.lambda * f64::from(iota * e.q_start.saturating_sub(node.q_end()) as Score)
        };
        let mut reference = if can_end { 0.0 } else { f64::NEG_INFINITY };
        for kind in SegmentKind::REARRANGED.into_iter().filter(|&k| follows(k)) {
            for (e, u) in outgoing(&s.lattice, &s.scales, successors, kind, slot.group) {
                reference = reference.max(s.scales[u].backward + gap(&e));
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
            inner[width - 1] = self.w.of(model.end(first, node.kind, d_used)) * (-reference).exp();
        }
        let mut tails = std::mem::take(&mut s.tails);
        for kind in SegmentKind::REARRANGED.into_iter().filter(|&k| follows(k)) {
            let iota = model.junction_per_base(node.kind, kind);
            tails.clear();
            tails.resize(width, 0.0);
            for (e, u) in outgoing(&s.lattice, &s.scales, successors, kind, slot.group) {
                let successor = s.lattice.slots[u];
                let factor = (s.scales[u].backward + gap(&e) - reference).exp();
                let backward = &s.backward[successor.values..successor.values + e.width()];
                for (o, z) in s.lattice.inside(&node, &e) {
                    tails[z] += factor * backward[o];
                }
                if let Some(beyond) = node.beyond(&e) {
                    let suffix = &s.suffix[successor.per_kind(node.kind, e.width())];
                    tails[width - 1] += factor * suffix[beyond];
                }
            }
            // inner(o) += Σ_{z ≥ o} tails(z) · e^{λι(z − o)}.
            let step = self.w.of(iota);
            let mut acc = 0.0;
            for o in (0..width).rev() {
                acc = acc * step + tails[o];
                inner[o] += acc;
            }
        }
        s.tails = tails;

        // B(y) = e^{−λ entry(y)} · Σ_{x ≥ y + min} X-weight(x) · inner(x).
        let min = self.dag.min_bases(node.kind) as usize;
        let bonus = self.presence(node);
        let mut acc = 0.0;
        for o in (0..width).rev() {
            acc += s.exit_w[node.cuts + o] * bonus * inner[o];
            inner[o] = acc;
        }
        let left = normalise(&mut inner);
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
        s.scales[t].backward = reference + left + scale + s.shift[b];
        s.work = inner;

        // Suffix sums for predecessors of each kind that can precede this slot.
        let mut done: Option<(Score, usize)> = None;
        for kind in SegmentKind::REARRANGED {
            if !(0..=d_used).any(|d| grammar.can_follow(kind, node.kind, d)) {
                continue;
            }
            let iota = model.junction_per_base(kind, node.kind);
            let range = slot.per_kind(kind, width);
            match done {
                Some((same, from)) if same == iota => {
                    s.suffix.copy_within(from..from + width, range.start)
                }
                _ => {
                    // S(o) = Σ_{o' ≥ o} B(o') · e^{λι(o' − o)}.
                    let step = self.w.of(iota);
                    let mut acc = 0.0;
                    for o in (0..width).rev() {
                        acc = acc * step + s.backward[slot.values + o];
                        s.suffix[range.start + o] = acc;
                    }
                    done = Some((iota, range.start));
                }
            }
        }
    }

    fn collect(&self, s: &PosteriorScratch, log_z: &[f64], log_partition: f64) -> Marginals {
        let mut placement = vec![0.0; self.placements.len()];
        for node in &s.lattice.nodes {
            for t in node.slots() {
                let scale = s.scales[t].forward + s.scales[t].backward;
                if scale == f64::NEG_INFINITY {
                    continue;
                }
                let values = s.lattice.slots[t].values;
                let range = values..values + node.width();
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
    use crate::adjudicate::units::{EvidenceModel, ReadModel, Units};
    use crate::adjudicate::{
        CalledSegment, DagScratch, EngineeringObjective, Mode, PathScore, VdjGrammar,
    };
    use crate::align::{Alignment, Op};
    use crate::alphabet::encode;
    use crate::clip::ClipIndex;
    use crate::evidence::Evidence;
    use crate::scoring::{AffineScoring, SubstitutionProfile};
    use crate::testing::placement;

    const LAMBDA: f64 = 0.6337;

    /// VJ, VDJ and VDDJ paths: the engine builds at most one D, the DAG
    /// handles any number.
    const TANDEM_D: VdjGrammar = VdjGrammar {
        max_d: 2,
        ..VdjGrammar::VDJ
    };

    /// The fixtures' scorings: `(2, −3)` for every kind, priors at `LAMBDA`.
    const EVIDENCE: EvidenceModel = EvidenceModel {
        lambda: LAMBDA,
        substitution: [(2, -3); 3],
        presence_lambda: [LAMBDA; 3],
    };

    /// Exhaustive search over grammatical paths and canonical cut pairs: the
    /// sum of `e^{nats}` (marginals) and the best key under `read`.
    struct Brute<'a> {
        dag: &'a PlacementDag,
        placements: &'a [Placement],
        nats: PathNats,
        read: ReadModel,
        z: f64,
        best: Option<PathScore<Units>>,
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

        /// `b`'s alignment evidence on `[entry, exit)`, in nats and in units:
        /// substitution columns on the evidence scale (N columns, scoring 0,
        /// are worth nothing), gaps at the prior λ.
        fn evidence(&self, b: &Placement, entry: u32, exit: u32) -> (f64, Units) {
            let scale = self.nats.evidence[b.kind.index()];
            let (mut nats, mut units) = (0.0, 0);
            let mut last: Option<Score> = None;
            for (_, before, after) in b.clip.matches().filter(|(q, ..)| (entry..exit).contains(q)) {
                if let Some(previous) = last {
                    nats += self.nats.lambda * f64::from(before - previous);
                    units += self.read.prior(before - previous);
                }
                let column = after - before;
                if column != 0 {
                    nats += scale.lambda * f64::from(column) + scale.mu;
                }
                units += self.read.column(b.kind, column);
                last = Some(after);
            }
            (nats, units)
        }

        /// `b`'s cut-independent term, its kind's presence cost.
        fn bonus(&self, b: &Placement) -> (f64, Units) {
            (
                -self.nats.presence[b.kind.index()],
                -self.read.presence(b.kind),
            )
        }

        /// A prior term scoring `k`.
        fn prior(&self, k: Score) -> (f64, Units) {
            (self.nats.lambda * f64::from(k), self.read.prior(k))
        }

        /// `path` so far is worth `nats` and `units`, its last segment entered at `entry`.
        #[allow(clippy::too_many_arguments)]
        fn visit(
            &mut self,
            path: &mut Vec<usize>,
            entry: u32,
            d_used: u8,
            first: SegmentKind,
            nats: f64,
            units: Units,
        ) {
            let (dag, placements) = (self.dag, self.placements);
            let b = &placements[*path.last().expect("non-empty")];
            let min_b = dag.min_retained[b.kind.index()].max(1);
            let bonus = self.bonus(b);
            for x in Self::exits(b).into_iter().filter(|&x| x >= entry + min_b) {
                let evidence = self.evidence(b, entry, x);
                let (nats, units) = (nats + evidence.0 + bonus.0, units + evidence.1 + bonus.1);
                if dag.grammar.can_end(first, b.kind, d_used) && x == b.q_end() {
                    let end = self.prior(dag.model.end(first, b.kind, d_used));
                    let weight = (nats + end.0).exp();
                    self.z += weight;
                    let key = PathScore {
                        score: units + end.1,
                        d_count: d_used,
                        support: path.iter().map(|&p| placements[p].score()).sum(),
                    };
                    self.best = self.best.max(Some(key));
                    for &p in path.iter() {
                        self.placement[p] += weight;
                    }
                    self.d_count[usize::from(d_used)] += weight;
                }
                for (j, e) in placements.iter().enumerate() {
                    // Kinds in grammar order wherever they start; one kind
                    // (tandem D) in order of start.
                    let same_kind_before = e.kind == b.kind && e.q_start() <= b.q_start();
                    if same_kind_before || !dag.grammar.can_follow(b.kind, e.kind, d_used) {
                        continue;
                    }
                    let iota = dag.model.junction_per_base(b.kind, e.kind);
                    let d_next = d_used + u8::from(e.kind == SegmentKind::D);
                    for y in Self::entries(e).into_iter().filter(|&y| y >= x) {
                        let junction = self.prior(iota * (y - x) as Score);
                        path.push(j);
                        self.visit(
                            path,
                            y,
                            d_next,
                            first,
                            nats + junction.0,
                            units + junction.1,
                        );
                        path.pop();
                    }
                }
            }
        }

        fn run(
            dag: &'a PlacementDag,
            placements: &'a [Placement],
            nats: PathNats,
            read: ReadModel,
        ) -> Self {
            let mut brute = Brute {
                dag,
                placements,
                nats,
                read,
                z: 0.0,
                best: None,
                placement: vec![0.0; placements.len()],
                d_count: vec![0.0; 3],
            };
            for (s, p) in placements.iter().enumerate() {
                if dag.grammar.can_start(p.kind) {
                    let d0 = u8::from(p.kind == SegmentKind::D);
                    brute.visit(&mut vec![s], p.q_start(), d0, p.kind, 0.0, 0);
                }
            }
            brute
        }

        /// The value of a called path under `read`, term by term.
        fn rescore(&self, segments: &[CalledSegment]) -> Units {
            let (dag, placements) = (self.dag, self.placements);
            let mut d_used = 0;
            let mut units = 0;
            for (k, segment) in segments.iter().enumerate() {
                let b = &placements[segment.placement];
                d_used += u8::from(b.kind == SegmentKind::D);
                units += self.evidence(b, segment.entry, segment.exit).1 + self.bonus(b).1;
                if let Some(before) = k.checked_sub(1).map(|k| &segments[k]) {
                    let kind = placements[before.placement].kind;
                    let inserted = segment.entry - before.exit;
                    units += self.read.prior(dag.model.junction(kind, b.kind, inserted));
                }
            }
            let kind = |s: &CalledSegment| placements[s.placement].kind;
            let (first, last) = (kind(&segments[0]), kind(&segments[segments.len() - 1]));
            units + self.read.prior(dag.model.end(first, last, d_used))
        }
    }

    fn assert_close(a: f64, b: f64, what: &str) {
        assert!(
            (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs())),
            "{what}: {a} vs {b}"
        );
    }

    /// Marginals under `nats` and the Viterbi key under `read` against
    /// exhaustive search; the traced path re-scores to its key.
    fn check(dag: &PlacementDag, placements: &[Placement], nats: PathNats, read: ReadModel) {
        let brute = Brute::run(dag, placements, nats, read);
        let m = dag
            .posterior(placements, &nats, &mut PosteriorScratch::default())
            .expect("a path exists");
        assert_close(m.log_partition, brute.z.ln(), "ln Z");
        for (i, (&got, &want)) in m.placement.iter().zip(&brute.placement).enumerate() {
            assert_close(got, want / brute.z, &format!("placement {i}"));
        }
        for (k, &want) in brute.d_count.iter().enumerate() {
            let got = m.d_count.get(k).copied().unwrap_or(0.0);
            assert_close(got, want / brute.z, &format!("D count {k}"));
        }
        let (key, segments) = dag
            .viterbi(placements, &read, &mut DagScratch::default())
            .expect("a path exists");
        assert_eq!(Some(key), brute.best);
        assert_eq!(brute.rescore(&segments), key.score);
    }

    /// The fixed-scale call with the posterior on `evidence` (the
    /// [`CallModel::FixedScale`] pair): a segment costs `⌊P*_k⌋` and
    /// `λ·P*_k`.
    fn check_fixed(dag: &PlacementDag, placements: &[Placement], evidence: EvidenceScale) {
        let presence = dag.model.presence;
        let nats = PathNats {
            lambda: LAMBDA,
            evidence: [evidence; 3],
            presence: presence.map(|cost| LAMBDA * cost),
        };
        check(dag, placements, nats, ReadModel::fixed(LAMBDA, presence));
    }

    /// One read model at divergence `p` for both passes: the Viterbi path is
    /// a mode of the posterior (the brute-force maximum of `e^{units/F}`).
    fn check_per_read(dag: &PlacementDag, placements: &[Placement], p: f64) {
        let read = ReadModel::at(&EVIDENCE, p, dag.model.presence);
        check(dag, placements, read.nats(), read);
    }

    /// Competing V alleles, two D candidates (one spurious, overlapping J)
    /// and two J alleles, with overlaps on every side.
    fn junction() -> (Vec<u8>, Vec<Placement>) {
        let query = encode(b"ACGTTGCAAGTCGGATTACAGGCTTAACGGTACCATGGA");
        let placements = junction_placements(&query);
        (query, placements)
    }

    fn junction_placements(query: &[u8]) -> Vec<Placement> {
        let p = |kind, reference: &[u8], at: u32| placement(kind, query, &encode(reference), at);
        vec![
            p(SegmentKind::V, b"ACGTTGCAAGTCGG", 0),
            p(SegmentKind::V, b"ACGTAGCAAGTCGA", 0),
            p(SegmentKind::D, b"GATTACAGG", 13),
            p(SegmentKind::D, b"TACAGCCTT", 16),
            p(SegmentKind::J, b"GCTTAACGGTAC", 21),
            p(SegmentKind::J, b"CTTTACGGTACC", 22),
        ]
    }

    /// The junction read with an N inside J1 and a J whose alignment has a
    /// one-base insertion: uninformative columns and gaps.
    fn gapped_junction() -> (Vec<u8>, Vec<Placement>) {
        let (mut query, _) = junction();
        query[27] = crate::alphabet::N;
        let mut placements = junction_placements(&query);
        let scoring = AffineScoring::new(2, -3, -9, -1);
        // J aligned to query 21..34 with query base 26 inserted.
        let reference = encode(b"GCTTACGGTACC");
        let profile = SubstitutionProfile::build(&scoring, &reference);
        let mut ops = vec![Op::Match; 5];
        ops.push(Op::Ins);
        ops.extend(std::iter::repeat_n(Op::Match, 7));
        let (mut q, mut r, mut score) = (21usize, 0usize, 0);
        for &op in &ops {
            match op {
                Op::Match => {
                    score += profile.at(query[q], r);
                    q += 1;
                    r += 1;
                }
                _ => {
                    score += scoring.gaps.open;
                    q += 1;
                }
            }
        }
        let alignment = Alignment {
            q_start: 21,
            q_end: q as u32,
            r_start: 0,
            r_end: r as u32,
            score,
            ops,
        };
        let clip = ClipIndex::build(&alignment, &query, &profile, scoring.gaps);
        placements.push(Placement {
            kind: SegmentKind::J,
            alleles: vec![1],
            clip,
            evidence: Evidence::default(),
            anchor: Some(Score::MIN),
            alignment,
        });
        (query, placements)
    }

    /// A J alignment that starts before the D it follows on the best path
    /// (its first bases mismatch the D region).
    fn j_reaching_over_d() -> (Vec<u8>, Vec<Placement>) {
        let ascii = b"TTTTTTTTTTTTTTTTTTTTACGACTTGCAGTCAAAAAAAAAAAAAAAAAAAAAAAA";
        let query = encode(ascii);
        let other = |b: &u8| match b {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            _ => b'G',
        };
        let mut j: Vec<u8> = ascii[21..34].iter().map(other).collect();
        j.extend_from_slice(&ascii[34..]);
        let placements = vec![
            placement(SegmentKind::V, &query, &encode(&ascii[..20]), 0),
            placement(SegmentKind::D, &query, &encode(&ascii[22..34]), 22),
            placement(SegmentKind::J, &query, &encode(&j), 21),
        ];
        (query, placements)
    }

    #[test]
    fn marginals_equal_exhaustive_summation() {
        let (_, placements) = junction();
        // Unrounded presence priors: the posterior uses 4.5 for a D and 3.5
        // for a J, the Viterbi call 4 and 3.
        let model = EngineeringObjective::priced([0.0, 4.5, 3.5]);
        let fixed = EvidenceScale::fixed(LAMBDA);
        for grammar in [VdjGrammar::VDJ, TANDEM_D, VdjGrammar::VJ] {
            check_fixed(
                &PlacementDag::new(grammar, model, [4, 3, 4], Mode::Fixed),
                &placements,
                fixed,
            );
        }
        let per_base = EngineeringObjective {
            insertion_per_base: -1,
            ..model
        };
        let dag = PlacementDag::new(TANDEM_D, per_base, [4, 3, 4], Mode::Fixed);
        check_fixed(&dag, &placements, fixed);
        // Evidence at a read's own divergence: a per-column term μ ≠ 0.
        check_fixed(&dag, &placements, EvidenceScale::at_divergence(0.02, 2, -3));
        // A D before a J whose alignment starts earlier.
        let (_, over) = j_reaching_over_d();
        for grammar in [VdjGrammar::VDJ, TANDEM_D] {
            check_fixed(
                &PlacementDag::new(grammar, model, [4, 3, 4], Mode::Fixed),
                &over,
                fixed,
            );
        }
        // An N column and an insertion, at the read's divergence and fixed.
        let (_, gapped) = gapped_junction();
        for evidence in [fixed, EvidenceScale::at_divergence(0.005, 2, -3)] {
            check_fixed(&dag, &gapped, evidence);
            check_fixed(
                &PlacementDag::new(VdjGrammar::VDJ, model, [4, 3, 4], Mode::Fixed),
                &gapped,
                evidence,
            );
        }
    }

    /// One integer model for both passes, from nearly unmutated to heavily
    /// mutated reads: the marginals are those of `e^{units/F}` and the
    /// Viterbi key is its exhaustive maximum, so the call is a mode.
    #[test]
    fn per_read_call_is_a_mode_of_the_posterior() {
        let model = EngineeringObjective::priced([0.0, 4.5, 3.5]);
        let per_base = EngineeringObjective {
            insertion_per_base: -1,
            ..model
        };
        let delta0 = EVIDENCE.log_odds_divergence(SegmentKind::V);
        let fixtures = [junction().1, gapped_junction().1, j_reaching_over_d().1];
        for p in [0.005, delta0, 0.3] {
            for placements in &fixtures {
                for grammar in [VdjGrammar::VDJ, TANDEM_D, VdjGrammar::VJ] {
                    for objective in [model, per_base] {
                        let dag = PlacementDag::new(grammar, objective, [4, 3, 4], Mode::Fixed);
                        check_per_read(&dag, placements, p);
                    }
                }
            }
        }
    }

    /// A strong per-base N/P cost across a 300-base V: the junction weights
    /// are carried relative to each cut, so nothing overflows.
    #[test]
    fn strong_junction_costs_stay_in_range() {
        let query = crate::sim::Rng::new(7).codes(340);
        let mutate = |range: std::ops::Range<usize>, every: usize| -> Vec<u8> {
            range
                .enumerate()
                .map(|(k, i)| {
                    if k % every == every - 1 {
                        (query[i] + 1) & 3
                    } else {
                        query[i]
                    }
                })
                .collect()
        };
        let placements = vec![
            placement(SegmentKind::V, &query, &mutate(0..300, 7), 0),
            placement(SegmentKind::J, &query, &mutate(290..330, 6), 290),
        ];
        for iota in [-4, -8] {
            let model = EngineeringObjective {
                insertion_per_base: iota,
                ..EngineeringObjective::priced([0.0, 11.0, 0.0])
            };
            let dag = PlacementDag::new(VdjGrammar::VDJ, model, [24, 4, 10], Mode::Fixed);
            check_fixed(&dag, &placements, EvidenceScale::fixed(LAMBDA));
            check_fixed(
                &dag,
                &placements,
                EvidenceScale::at_divergence(0.001, 2, -3),
            );
            check_per_read(&dag, &placements, 0.001);
        }
    }

    /// Every path crosses a 300-base N/P gap at ι = −4 (−760 nats), which
    /// is below the `f64` range as a factor; it must go into the log scale.
    #[test]
    fn long_costly_gaps_keep_their_mass() {
        let query = crate::sim::Rng::new(11).codes(360);
        let placements = vec![
            placement(SegmentKind::V, &query, &query[..30], 0),
            placement(SegmentKind::J, &query, &query[330..], 330),
        ];
        let model = EngineeringObjective {
            insertion_per_base: -4,
            ..EngineeringObjective::priced([0.0, 11.0, 0.0])
        };
        let grammar = VdjGrammar {
            require_v: true,
            require_j: true,
            ..VdjGrammar::VDJ
        };
        let dag = PlacementDag::new(grammar, model, [4, 3, 4], Mode::Fixed);
        check_fixed(&dag, &placements, EvidenceScale::fixed(LAMBDA));
        check_per_read(&dag, &placements, 0.01);
    }

    /// V and J overlap on 540 matched bases: a junction product would fall
    /// out of the `f64` range, so the read gets no marginals rather than
    /// wrong ones.
    #[test]
    fn overlaps_beyond_the_f64_range_are_refused() {
        let query = crate::sim::Rng::new(5).codes(700);
        let placements = vec![
            placement(SegmentKind::V, &query, &query[..600], 0),
            placement(SegmentKind::J, &query, &query[60..], 60),
        ];
        let dag = PlacementDag::new(
            VdjGrammar::VDJ,
            EngineeringObjective::priced([0.0, 11.0, 0.0]),
            [4, 3, 4],
            Mode::Fixed,
        );
        let nats = EVIDENCE.nats(0.01, dag.model.presence);
        let m = dag.posterior(&placements, &nats, &mut PosteriorScratch::default());
        assert!(m.is_none(), "{m:?}");
        // A junction-sized overlap is fine.
        let short = vec![
            placement(SegmentKind::V, &query, &query[..300], 0),
            placement(SegmentKind::J, &query, &query[260..360], 260),
        ];
        check_fixed(&dag, &short, EvidenceScale::at_divergence(0.01, 2, -3));
    }

    /// The passes' buffers carry over between reads uncleared: a slot's
    /// running maxima and sums are written before any slot reads them. A
    /// large lattice, then a small one, then large ones again give, through
    /// one scratch each, exactly what fresh buffers give.
    #[test]
    fn reused_scratch_gives_what_fresh_scratch_gives() {
        let model = EngineeringObjective::priced([0.0, 4.5, 3.5]);
        let big = PlacementDag::new(TANDEM_D, model, [4, 3, 4], Mode::Fixed);
        let small = PlacementDag::new(VdjGrammar::VJ, model, [1, 1, 1], Mode::Fixed);
        let (_, junction) = junction();
        let (_, gapped) = gapped_junction();
        let (_, over) = j_reaching_over_d();
        let vj = [junction[0].clone(), junction[5].clone()];
        let read = ReadModel::at(&EVIDENCE, 0.02, model.presence);
        let nats = read.nats();
        let (mut dag_scratch, mut scratch) = (DagScratch::default(), PosteriorScratch::default());
        for (dag, placements) in [
            (&big, &junction[..]),
            (&small, &vj[..]),
            (&big, &gapped[..]),
            (&big, &over[..]),
            (&big, &junction[..]),
        ] {
            let fresh = (
                dag.viterbi(placements, &read, &mut DagScratch::default()),
                dag.posterior(placements, &nats, &mut PosteriorScratch::default()),
            );
            let reused = (
                dag.viterbi(placements, &read, &mut dag_scratch),
                dag.posterior(placements, &nats, &mut scratch),
            );
            assert!(fresh.0.is_some() && fresh.1.is_some());
            assert_eq!(fresh, reused);
        }
    }

    /// A lattice without a D (a light chain, or a D-free database) has one
    /// D count: D presence is `+0.0`, which prints as `0.000000`.
    #[test]
    fn d_presence_without_a_d_class_is_positive_zero() {
        let (_, placements) = junction();
        let dag = PlacementDag::new(
            VdjGrammar::VJ,
            EngineeringObjective::priced([0.0, 11.0, 0.0]),
            [4, 3, 4],
            Mode::Fixed,
        );
        let nats = PathNats {
            lambda: LAMBDA,
            evidence: [EvidenceScale::fixed(LAMBDA); 3],
            presence: dag.model.presence.map(|cost| LAMBDA * cost),
        };
        let m = dag
            .posterior(&placements, &nats, &mut PosteriorScratch::default())
            .expect("a path exists");
        assert_eq!(m.d_count.len(), 1);
        assert!(m.d_present().is_sign_positive(), "{}", m.d_present());
        assert_eq!(crate::report::probability(m.d_present()), "0.000000");
    }

    #[test]
    fn marginals_are_probabilities() {
        let (_, placements) = junction();
        let dag = PlacementDag::new(
            VdjGrammar::VDJ,
            EngineeringObjective::priced([0.0, 11.0, 0.0]),
            [4, 3, 4],
            Mode::Fixed,
        );
        let nats = PathNats {
            lambda: LAMBDA,
            evidence: [EvidenceScale::fixed(LAMBDA); 3],
            presence: dag.model.presence.map(|cost| LAMBDA * cost),
        };
        let m = dag
            .posterior(&placements, &nats, &mut PosteriorScratch::default())
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
