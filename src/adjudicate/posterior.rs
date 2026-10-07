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
//! underconfident and heavily mutated ones overconfident. The prior terms
//! (`P_D`, pair bonus, N/P) stay at the λ they were derived with, but the
//! alignment evidence is rescored as log-odds at the read's divergence `p̂`
//! ([`divergence`]): a match is worth `ln 4(1 − p̂)` nats and a mismatch
//! `ln(4p̂/3)` ([`EvidenceScale::at_divergence`]). For a match/mismatch
//! scoring that is an affine map of each substitution column's score,
//! `λ_r·s + μ_r`, which keeps every clip score separable. A column scoring 0
//! (an `N`) is worth 0 either way, and gaps keep the prior λ: `p̂` says
//! nothing about indel rates. At `p̂ = 11.2 %` it is the fixed-λ model. The posterior's mode can then differ from the
//! Viterbi call, which stays on the fixed scores; the reported probability
//! is that of the reported call under the per-read model.
//!
//! **Domain.** As in the Viterbi DAG, the first segment keeps its start and
//! the last segment its end: V and J are extended to the read's ends on
//! purpose ([`crate::refine::extend_outer`]). Between neighbours, every pair
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

use super::lattice::{Group, KINDS, Lattice, Node};
use super::model::{Grammar, PathModel};
use super::{CallPath, PlacementDag};
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

/// Nats per unit of alignment evidence: a substitution column scoring
/// `s ≠ 0` is worth `lambda · s + mu`, an `N` column (`s = 0`) nothing. Gaps
/// are read at the prior λ ([`PosteriorModel::lambda`]).
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
            // N columns (score 0) carry no evidence either way.
            columns += usize::from(exit != entry);
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
    /// Per Match column: query position, entry and exit prefix scores, and
    /// the evidence (nats) before and through it.
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
    /// with each segment's alignment evidence on `evidence[kind]` and the
    /// model's prior terms at `lambda` nats per unit.
    ///
    /// `None` when no grammatical path exists, or when linear arithmetic
    /// cannot hold the read: a placement's evidence range ([`MAX_SPREAD`]) or
    /// an overlap's ([`OVERLAP_SPREAD`]) is too wide.
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
        s.scales.clear();
        s.forward.clear();
        s.backward.clear();
        s.running.clear();
        s.suffix.clear();
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
        s.lattice.layout(placements);
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
            let ends = SegmentKind::ALL.iter().any(|&f| grammar.can_end(f, p.kind));
            let preceded = SegmentKind::ALL
                .iter()
                .any(|&k| any_d(&|d| grammar.can_follow(k, p.kind, d)));
            let followed = SegmentKind::ALL
                .iter()
                .any(|&k| any_d(&|d| grammar.can_follow(p.kind, k, d)));
            let entry_used = |o: u32| preceded || (starts && o == 0);
            let exit_used = |o: u32| followed || (ends && o == len);

            // Evidence nats: a substitution column scoring s adds λ_r·s + μ_r,
            // an uninformative one (an N, s = 0) nothing; gaps keep λ.
            s.columns.clear();
            let mut previous: Option<(Score, f64)> = None;
            let (mut in_lo, mut in_hi) = (f64::INFINITY, f64::NEG_INFINITY);
            let (mut out_lo, mut out_hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for (q, entry, exit) in p.clip.matches() {
                let before = match previous {
                    None => lambda0 * f64::from(entry),
                    Some((exit_prefix, nats)) => nats + lambda0 * f64::from(entry - exit_prefix),
                };
                let column = exit - entry;
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
                s.columns.push((q, entry, exit, before, through));
                previous = Some((exit, through));
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
            let mut last: Option<(Score, Score)> = None;
            for &(q, entry, exit, nats, through) in &s.columns {
                if let Some((last_entry, last_exit)) = last {
                    let (column, gap) = (last_exit - last_entry, entry - last_exit);
                    before *= substitution(column, false) * self.w.of(-gap);
                    inv *= substitution(column, true) * self.w.of(gap);
                }
                last = Some((entry, exit));
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
                    s.exit_w[cuts + o as usize + 1] = inv * substitution(exit - entry, true);
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
        s.running.resize(KINDS * values, 0.0);
        s.suffix.resize(KINDS * values, 0.0);

        // A predecessor ending before this node reaches it across an N/P
        // gap; its weight joins the log scale, so a long, costly gap cannot
        // underflow when every path crosses it. 0 for an overlapping one.
        let gap = |b: &Node| {
            let iota = self.dag.model.junction_per_base(b.kind, node.kind);
            self.w.lambda * f64::from(iota * node.q_start.saturating_sub(b.q_end()) as Score)
        };
        for t in slots {
            let slot = s.lattice.slots[t];
            let starts_here = start == Some(slot.group);
            let scales = &s.scales;
            let live = |u: usize| scales[u].exit > f64::NEG_INFINITY;
            let reference = s
                .lattice
                .incoming(e, grammar, slot.group, &live)
                .map(|(b, u)| scales[u].exit + gap(&b))
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
            let entries = s.lattice.entries(&node);
            for (b, u) in s.lattice.incoming(e, grammar, slot.group, &live) {
                let factor = (scales[u].exit + gap(&b) - reference).exp();
                let r = &s.running[s.lattice.slots[u].per_kind(node.kind, b.width())];
                // Entries inside b: exits up to the entry (none before b's
                // start, which a predecessor of another kind may have).
                for &o in entries {
                    let y = node.q_start + o;
                    if y >= b.q_end() {
                        break;
                    }
                    let Some(z) = y.checked_sub(b.q_start) else {
                        continue;
                    };
                    // r[z] already carries e^{λι(z − x)} for each exit x.
                    f[o as usize] += factor * r[z as usize];
                }
                // Entries at or beyond b's end: all exits, added once.
                let beyond = b.q_end().saturating_sub(node.q_start) as usize;
                if beyond < width {
                    tails[b.kind.index() * width + beyond] += factor * r[b.len as usize];
                }
            }
            for kind in SegmentKind::ALL {
                let tail = &tails[kind.index() * width..(kind.index() + 1) * width];
                if tail.iter().all(|&v| v == 0.0) {
                    continue;
                }
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
        let bonus = self.w.of(self
            .dag
            .model
            .segment_bonus(&self.placements[node.placement], d_used));
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
            for kind in SegmentKind::ALL {
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
            if self.dag.grammar.can_end(first, node.kind) && x[width - 1] > 0.0 {
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
        let can_end = grammar.can_end(first, node.kind);
        let follows = |kind: SegmentKind| grammar.can_follow(node.kind, kind, d_used);
        // As in the forward pass, an N/P gap to a successor that starts after
        // this node's end joins the log scale.
        let gap = |e: &Node| {
            let iota = model.junction_per_base(node.kind, e.kind);
            self.w.lambda * f64::from(iota * e.q_start.saturating_sub(node.q_end()) as Score)
        };
        let mut reference = if can_end { 0.0 } else { f64::NEG_INFINITY };
        for kind in SegmentKind::ALL.into_iter().filter(|&k| follows(k)) {
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
        for kind in SegmentKind::ALL.into_iter().filter(|&k| follows(k)) {
            let iota = model.junction_per_base(node.kind, kind);
            tails.clear();
            tails.resize(width, 0.0);
            for (e, u) in outgoing(&s.lattice, &s.scales, successors, kind, slot.group) {
                let successor = s.lattice.slots[u];
                let factor = (s.scales[u].backward + gap(&e) - reference).exp();
                let backward = &s.backward[successor.values..successor.values + e.width()];
                for &o in s.lattice.entries(&e) {
                    let y = e.q_start + o;
                    if y >= node.q_end() {
                        break;
                    }
                    let Some(z) = y.checked_sub(node.q_start) else {
                        continue;
                    };
                    tails[z as usize] += factor * backward[o as usize];
                }
                let beyond = node.q_end().saturating_sub(e.q_start) as usize;
                if beyond < e.width() {
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
        let bonus = self
            .w
            .of(model.segment_bonus(&self.placements[node.placement], d_used));
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
        for kind in SegmentKind::ALL {
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
        /// Substitution columns at the evidence scale (N columns, scoring
        /// 0, are worth nothing), gaps at the prior λ.
        fn evidence(&self, b: &Placement, entry: u32, exit: u32) -> f64 {
            let span = b.clip.span(entry, exit).expect("canonical cuts");
            let (mut nats, mut substitutions) = (0.0, 0);
            for (_, before, after) in b.clip.matches().filter(|(q, ..)| (entry..exit).contains(q)) {
                let column = after - before;
                substitutions += column;
                if column != 0 {
                    nats += self.evidence.lambda * f64::from(column) + self.evidence.mu;
                }
            }
            nats + LAMBDA * f64::from(span.score - substitutions)
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
                    // Kinds in grammar order wherever they start; one kind
                    // (tandem D) in order of start.
                    let same_kind_before = e.kind == b.kind && e.q_start() <= b.q_start();
                    if same_kind_before || !dag.grammar.can_follow(b.kind, e.kind, d_used) {
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
            anchor: true,
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
        // A D before a J whose alignment starts earlier.
        let (_, over) = j_reaching_over_d();
        for grammar in [VdjGrammar::VDJ, VdjGrammar::TANDEM_D] {
            check(&PlacementDag::new(grammar, model, [4, 3, 4]), &over, fixed);
        }
        // An N column and an insertion, at the read's divergence and fixed.
        let (_, gapped) = gapped_junction();
        for evidence in [fixed, EvidenceScale::at_divergence(0.005, 2, -3)] {
            check(&dag, &gapped, evidence);
            check(
                &PlacementDag::new(VdjGrammar::VDJ, model, [4, 3, 4]),
                &gapped,
                evidence,
            );
        }
    }

    /// A strong per-base N/P cost across a 300-base V: the junction weights
    /// are carried relative to each cut, so nothing overflows.
    #[test]
    fn strong_junction_costs_stay_in_range() {
        let mut x = 7u64;
        let query: Vec<u8> = (0..340)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((x >> 33) & 3) as u8
            })
            .collect();
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
                ..EngineeringObjective::default()
            };
            let dag = PlacementDag::new(VdjGrammar::VDJ, model, [24, 4, 10]);
            check(&dag, &placements, EvidenceScale::fixed(LAMBDA));
            check(
                &dag,
                &placements,
                EvidenceScale::at_divergence(0.001, 2, -3),
            );
        }
    }

    /// Every path crosses a 300-base N/P gap at ι = −4 (−760 nats), which
    /// is below the `f64` range as a factor; it must go into the log scale.
    #[test]
    fn long_costly_gaps_keep_their_mass() {
        let mut x = 11u64;
        let query: Vec<u8> = (0..360)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((x >> 33) & 3) as u8
            })
            .collect();
        let placements = vec![
            placement(SegmentKind::V, &query, &query[..30], 0),
            placement(SegmentKind::J, &query, &query[330..], 330),
        ];
        let model = EngineeringObjective {
            insertion_per_base: -4,
            ..EngineeringObjective::default()
        };
        let grammar = VdjGrammar {
            require_v: true,
            require_j: true,
            ..VdjGrammar::VDJ
        };
        check(
            &PlacementDag::new(grammar, model, [4, 3, 4]),
            &placements,
            EvidenceScale::fixed(LAMBDA),
        );
    }

    /// V and J overlap on 540 matched bases: a junction product would fall
    /// out of the `f64` range, so the read gets no marginals rather than
    /// wrong ones.
    #[test]
    fn overlaps_beyond_the_f64_range_are_refused() {
        let mut x = 5u64;
        let query: Vec<u8> = (0..700)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((x >> 33) & 3) as u8
            })
            .collect();
        let placements = vec![
            placement(SegmentKind::V, &query, &query[..600], 0),
            placement(SegmentKind::J, &query, &query[60..], 60),
        ];
        let dag = PlacementDag::new(VdjGrammar::VDJ, EngineeringObjective::default(), [4, 3, 4]);
        let evidence = EvidenceScale::at_divergence(0.01, 2, -3);
        let m = dag.posterior(
            &placements,
            LAMBDA,
            &[evidence; 3],
            &mut PosteriorScratch::default(),
        );
        assert!(m.is_none(), "{m:?}");
        // A junction-sized overlap is fine.
        let short = vec![
            placement(SegmentKind::V, &query, &query[..300], 0),
            placement(SegmentKind::J, &query, &query[260..360], 260),
        ];
        check(&dag, &short, evidence);
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
