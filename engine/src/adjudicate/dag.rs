//! Placement DAG: bounded-topology semi-Markov Viterbi over clipped placements.
//!
//! The lattice ([`super::lattice`]) orders placements by segment kind, then
//! query start. A path visits placements in read order of their retained
//! bases: a transition between kinds needs only an exit cut `x` and an entry
//! cut `y ≥ x`, wherever the two alignments start, so a D can precede a J
//! whose alignment reaches back over it. Each placement on the path is
//! clipped independently: the earlier of two consecutive placements keeps
//! query bases before its **exit cut** `x`, the later one keeps bases from
//! its **entry cut** `y ≥ x` on, and the bases in between are N/P. This is
//! the full segmentation model: when two placements overlap, the cut where
//! D wins the overlap, the cut where V/J win it, every allocation in
//! between, and leaving any part of the overlap or of either
//! alignment to N/P are all evaluated exactly. The first segment keeps its
//! start and the last its end (V and J are extended to the read's ends on
//! purpose, see [`crate::flank::extend_outer`]).
//!
//! Cuts are *canonical* ([`crate::clip::ClipIndex::matches`]): an entry on
//! the first retained base, an exit right after the last, so every distinct
//! clipped span is one `(entry, exit)` pair.
//!
//! **Values.** Every term of a path is read in one [`ReadModel`]'s integer
//! units: the fixed scores themselves ([`ReadModel::fixed`]), or the read's
//! own log-odds at its divergence ([`ReadModel::at`], the units the posterior
//! reads too). A clip's value is then `exit(x) − entry(y)` for prefix values
//! accumulated once per node along its columns: a substitution column adds
//! [`ReadModel::column`], the gap run before it [`ReadModel::prior`] of its
//! score. Canonical cuts never split a gap run, so this is exact, as for
//! scores ([`crate::clip`]).
//!
//! **Slots.** Each slot of the lattice holds, per entry cut, the best key
//! (value, then germline support) of a path prefix entering there. The
//! recursion is the max-plus twin of the posterior's sum-product
//! ([`super::posterior`]):
//!
//! ```text
//! X(b, x)  = exit(x) + bonus + max_{y ≤ x − min} [F(b, y) − entry(y)]
//! F(e, y') = max_b max_{x ≤ y'} [X(b, x) + ι·(y' − x)]          (0 at a path start)
//! best     = max_b X(b, end of b) + end(g, b)
//! ```
//!
//! **Trick — separable values → running maxima.** A clipped span is worth
//! `exit(x) − entry(y)`, the model adds a cut-independent bonus, and the N/P
//! value is linear in the inserted length.
//! So the inner maximum over entries is a running maximum along the node, and
//! the maximum over exits `x ≤ y'` another, read in O(1) per entry. On an edge
//! `b → e`, every entry of `e` at or beyond `b`'s end sees the same best exit,
//! so it is added once to a running maximum over `e`'s entries: an edge costs
//! `O(overlap + 1)`, a read `O(Σ_nodes len + Σ_edges (overlap + 1))`, plus
//! one scan of its predecessors' slots per node, which opens the node's
//! slots and records their edges ([`super::lattice`]).
//!
//! Ties: final paths compare by [`PathScore`] (value, then fewer D, then
//! germline support, the summed unclipped fixed scores), so a with-D path must
//! strictly beat the best no-D path. Integer values make every tie exact.
//! Among equal keys: a node keeps its earliest entry, so a junction base that
//! scores the same on either side goes to the later segment; for that entry,
//! the predecessor's latest exit; and across predecessors, a path start, then
//! the first predecessor in DAG order.

use super::decode::{DAlleleDecoding, Topology, Within, max_marginal_d};
use super::lattice::{Group, KINDS, Lattice, column_steps};
use super::model::{Grammar, PathModel, VdjGrammar};
use super::posterior::{Marginals, PosteriorScratch};
use super::units::{CallModel, EvidenceModel, PathNats, ReadModel, Units, divergence};
use super::{Adjudicator, CallPath, CalledSegment, EngineeringObjective, PathScore};
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// No back pointer: a path start.
const START: u32 = u32::MAX;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlacementDag<G = VdjGrammar, M = EngineeringObjective> {
    pub grammar: G,
    pub model: M,
    /// Minimum query bases a clipped segment must keep, per rearranged kind
    /// ([`SegmentKind::REARRANGED`]), indexed by [`SegmentKind::index`].
    pub min_retained: [u32; 3],
    mode: Mode,
}

/// What a [`PlacementDag`] computes for a read: the model its call
/// maximises and whether it reports posteriors. Each variant holds what it
/// needs, so no mode asks for a per-read call or posteriors without the
/// evidence model they read. The engine's DAG comes from
/// [`crate::config::EngineConfig`] ([`crate::config::DefaultAssigner::from_config`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mode {
    /// Calls on the fixed scores, without posteriors.
    Fixed,
    /// Alignment evidence read as log-odds at each read's divergence.
    Read {
        evidence: EvidenceModel,
        /// The model the call maximises.
        calls: CallModel,
        /// Posterior marginals ([`Adjudicator::marginals`]), with the rule
        /// that chooses the reported D allele from them
        /// ([`Adjudicator::decode`]); `None`: no posteriors.
        posterior: Option<DAlleleDecoding>,
    },
}

/// Value, then germline support (the summed unclipped scores of the path's
/// placements), compared lexicographically.
type Key = (Units, Score);

/// A best value with where it came from.
#[derive(Clone, Copy, Debug)]
struct Cell {
    key: Key,
    /// Entry cells: the predecessor slot ([`START`] at a path start).
    /// Exit cells: unused. Running cells: unused.
    slot: u32,
    /// Entry cells: the predecessor's exit offset. Exit cells: this slot's
    /// entry offset. Running cells: this slot's exit offset.
    at: u32,
}

impl Cell {
    #[inline]
    fn shift(self, by: Units) -> Self {
        Self {
            key: (self.key.0 + by, self.key.1),
            ..self
        }
    }
}

/// Keep `candidate` if it is better (strictly, or also when equal with `ties`).
#[inline]
fn improve(best: &mut Option<Cell>, candidate: Cell, ties: bool) {
    let better = match best {
        None => true,
        Some(current) => candidate.key > current.key || (ties && candidate.key == current.key),
    };
    if better {
        *best = Some(candidate);
    }
}

/// Keep `candidate` if its key is better, or equal and its predecessor comes
/// first (a path start, then predecessors in node order). Candidates for one
/// entry may be folded in any order: the result does not depend on it.
#[inline]
fn prefer(best: &mut Option<Cell>, candidate: Cell) {
    // START (u32::MAX) wraps to 0, slot u to u + 1; slots follow node order.
    let rank = |cell: &Cell| cell.slot.wrapping_add(1);
    let better = match best {
        None => true,
        Some(current) => {
            candidate.key > current.key
                || (candidate.key == current.key && rank(&candidate) < rank(current))
        }
    };
    if better {
        *best = Some(candidate);
    }
}

/// Per-thread buffers for [`PlacementDag`].
#[derive(Debug, Default)]
pub struct DagScratch {
    lattice: Lattice,
    /// Prefix value at each canonical entry / exit cut (`None` elsewhere).
    entry: Vec<Option<Units>>,
    exit: Vec<Option<Units>>,
    /// Per slot: whether any path leaves it.
    open: Vec<bool>,
    /// Per slot value: best path prefix entering at each offset.
    forward: Vec<Option<Cell>>,
    /// Per slot value: best path prefix through each exit offset (segment included).
    exits: Vec<Option<Cell>>,
    /// Per slot and successor kind: best `X(x) − ι·x` over exits `≤ z`.
    running: Vec<Option<Cell>>,
    tails: Vec<Option<Cell>>,
    posterior: PosteriorScratch,
}

impl<G: Grammar, M: PathModel> PlacementDag<G, M> {
    pub fn new(grammar: G, model: M, min_retained: [u32; 3], mode: Mode) -> Self {
        Self {
            grammar,
            model,
            min_retained,
            mode,
        }
    }

    pub fn mode(&self) -> Mode {
        self.mode
    }

    #[inline]
    pub(super) fn min_bases(&self, kind: SegmentKind) -> u32 {
        self.min_retained[kind.index()].max(1)
    }

    /// The lattice of the kept placements with the prefix value of every
    /// canonical cut under `read`.
    fn layout(
        &self,
        placements: &[Placement],
        keep: impl Fn(usize) -> bool,
        read: &ReadModel,
        s: &mut DagScratch,
    ) {
        s.lattice.layout(placements, keep);
        s.entry.clear();
        s.entry.resize(s.lattice.cells(), None);
        s.exit.clear();
        s.exit.resize(s.lattice.cells(), None);
        for node in &s.lattice.nodes {
            let mut value: Units = 0;
            for (q, gap, column) in column_steps(&placements[node.placement]) {
                let o = (q - node.q_start) as usize;
                value += read.prior(gap);
                s.entry[node.cuts + o] = Some(value);
                value += read.column(node.kind, column);
                s.exit[node.cuts + o + 1] = Some(value);
            }
        }
        s.open.clear();
        s.forward.clear();
        s.exits.clear();
        // `running` is not cleared: `emit` writes a slot's running maxima
        // for every kind that can follow it before `gather` reads them, by
        // the same `can_follow`.
    }

    /// Open node `e`'s slots and fill their entry values from path starts
    /// and from every predecessor's running exit maxima.
    fn gather(&self, read: &ReadModel, s: &mut DagScratch, e: usize) {
        let grammar = &self.grammar;
        let node = s.lattice.nodes[e];
        let width = node.width();
        let start = grammar
            .can_start(node.kind)
            .then(|| Group::start(node.kind))
            .filter(|_| s.entry[node.cuts].is_some());
        let open = &s.open;
        let slots = s.lattice.open(e, grammar, start, |t| open[t]);
        s.open.resize(s.lattice.slots.len(), false);
        s.forward.resize(s.lattice.values(), None);
        s.exits.resize(s.lattice.values(), None);
        // Grown, never cleared (see `layout`).
        if s.running.len() < KINDS * s.lattice.values() {
            s.running.resize(KINDS * s.lattice.values(), None);
        }
        if s.tails.len() < KINDS * width {
            s.tails.resize(KINDS * width, None);
        }

        for t in slots {
            let slot = s.lattice.slots[t];
            let forward = &mut s.forward[slot.values..slot.values + width];
            if start == Some(slot.group) {
                prefer(
                    &mut forward[0],
                    Cell {
                        key: (0, 0),
                        slot: START,
                        at: 0,
                    },
                );
            }
            // Per predecessor kind, the entries at or beyond its nodes' ends:
            // cleared when an edge of that kind first reaches them.
            let mut tails = [false; KINDS];
            for (b, u) in s.lattice.incoming(t) {
                let iota = read.prior(self.model.junction_per_base(b.kind, node.kind));
                let running = &s.running[s.lattice.slots[u].per_kind(node.kind, b.width())];
                // Entries inside b: the best exit at or before the entry.
                for (o, z) in s.lattice.inside(&b, &node) {
                    if let Some(cell) = running[z] {
                        let candidate = Cell {
                            key: (cell.key.0 + iota * z as Units, cell.key.1),
                            slot: u as u32,
                            at: cell.at,
                        };
                        prefer(&mut forward[o], candidate);
                    }
                }
                // Entries at or beyond b's end: b's best exit overall, once.
                if let (Some(beyond), Some(cell)) = (b.beyond(&node), running[b.len as usize]) {
                    let candidate = Cell {
                        key: (
                            cell.key.0
                                + iota * (Units::from(node.q_start) - Units::from(b.q_start)),
                            cell.key.1,
                        ),
                        slot: u as u32,
                        at: cell.at,
                    };
                    let k = b.kind.index();
                    let tail = &mut s.tails[k * width..(k + 1) * width];
                    if !std::mem::replace(&mut tails[k], true) {
                        tail.fill(None);
                    }
                    prefer(&mut tail[beyond], candidate);
                }
            }
            for kind in SegmentKind::REARRANGED {
                if !tails[kind.index()] {
                    continue;
                }
                let tail = &s.tails[kind.index() * width..(kind.index() + 1) * width];
                let iota = read.prior(self.model.junction_per_base(kind, node.kind));
                let mut best: Option<Cell> = None;
                for (o, cell) in tail.iter().enumerate() {
                    if let Some(cell) = *cell {
                        prefer(&mut best, cell);
                    }
                    if let (Some(cell), true) = (best, s.entry[node.cuts + o].is_some()) {
                        prefer(&mut forward[o], cell.shift(iota * o as Units));
                    }
                }
            }
        }
    }

    /// Exit values `X(e, ·)` of slot `t`, their running maxima per successor
    /// kind, and the path ending here if the grammar allows one.
    fn emit(
        &self,
        placements: &[Placement],
        read: &ReadModel,
        s: &mut DagScratch,
        t: usize,
    ) -> Option<(PathScore<Units>, Cell)> {
        let slot = s.lattice.slots[t];
        let Group { d_used, first } = slot.group;
        let node = s.lattice.nodes[slot.node as usize];
        let width = node.width();
        let placement = &placements[node.placement];
        let min = self.min_bases(node.kind) as usize;
        let bonus = -read.presence(node.kind);
        let support = placement.score();

        // Running best of F(y) − entry(y) over entries, read min bases later.
        let mut inner: Option<Cell> = None;
        let mut open = false;
        for o in 0..width {
            if o >= min {
                let entered = o - min;
                if let (Some(cell), Some(prefix)) = (
                    s.forward[slot.values + entered],
                    s.entry[node.cuts + entered],
                ) {
                    improve(
                        &mut inner,
                        Cell {
                            key: (cell.key.0 - prefix, cell.key.1),
                            slot: 0,
                            at: entered as u32,
                        },
                        false,
                    );
                }
            }
            let exit = match (inner, s.exit[node.cuts + o]) {
                (Some(cell), Some(prefix)) => {
                    open = true;
                    Some(Cell {
                        key: (cell.key.0 + prefix + bonus, cell.key.1 + support),
                        slot: 0,
                        at: cell.at,
                    })
                }
                _ => None,
            };
            s.exits[slot.values + o] = exit;
        }
        s.open[t] = open;
        if !open {
            return None;
        }
        for kind in SegmentKind::REARRANGED {
            if !self.grammar.can_follow(node.kind, kind, d_used) {
                continue;
            }
            let iota = read.prior(self.model.junction_per_base(node.kind, kind));
            let range = slot.per_kind(kind, width);
            let mut best: Option<Cell> = None;
            for o in 0..width {
                if let Some(cell) = s.exits[slot.values + o] {
                    // `ties`: the later exit retains more of the placement.
                    improve(
                        &mut best,
                        Cell {
                            key: (cell.key.0 - iota * o as Units, cell.key.1),
                            slot: 0,
                            at: o as u32,
                        },
                        true,
                    );
                }
                s.running[range.start + o] = best;
            }
        }
        let end = s.exits[slot.values + width - 1]?;
        self.grammar.can_end(first, node.kind, d_used).then(|| {
            let score = PathScore {
                score: end.key.0 + read.prior(self.model.end(first, node.kind, d_used)),
                d_count: d_used,
                support: end.key.1,
            };
            (
                score,
                Cell {
                    key: end.key,
                    slot: t as u32,
                    at: (width - 1) as u32,
                },
            )
        })
    }

    /// Walk back from the best path end.
    fn trace(&self, placements: &[Placement], s: &DagScratch, end: Cell) -> Vec<CalledSegment> {
        let mut segments = Vec::with_capacity(4);
        let (mut t, mut exit) = (end.slot as usize, end.at);
        loop {
            let slot = s.lattice.slots[t];
            let node = s.lattice.nodes[slot.node as usize];
            let entry = s.exits[slot.values + exit as usize]
                .expect("traced exits are reachable")
                .at;
            let placement = node.placement;
            let (y, x) = (node.q_start + entry, node.q_start + exit);
            segments.push(CalledSegment {
                placement,
                entry: y,
                exit: x,
                span: placements[placement]
                    .clip
                    .span(y, x)
                    .expect("canonical cuts keep at least min_retained bases"),
                allele: None,
            });
            let back =
                s.forward[slot.values + entry as usize].expect("traced entries are reachable");
            if back.slot == START {
                break;
            }
            (t, exit) = (back.slot as usize, back.at);
        }
        segments.reverse();
        segments
    }
}

impl<G: Grammar, M: PathModel> PlacementDag<G, M> {
    /// The best path under `read`, with its key.
    pub(crate) fn viterbi(
        &self,
        placements: &[Placement],
        read: &ReadModel,
        s: &mut DagScratch,
    ) -> Option<(PathScore<Units>, Vec<CalledSegment>)> {
        self.best_path(placements, |_| true, read, s)
    }

    /// The best path under `read` with `topology` through the placements
    /// that `keep`.
    pub(crate) fn conditioned(
        &self,
        placements: &[Placement],
        keep: impl Fn(usize) -> bool,
        topology: Topology,
        read: &ReadModel,
        s: &mut DagScratch,
    ) -> Option<(PathScore<Units>, Vec<CalledSegment>)> {
        let within = PlacementDag {
            grammar: Within {
                grammar: &self.grammar,
                topology,
            },
            model: &self.model,
            min_retained: self.min_retained,
            mode: self.mode,
        };
        within.best_path(placements, keep, read, s)
    }

    fn best_path(
        &self,
        placements: &[Placement],
        keep: impl Fn(usize) -> bool,
        read: &ReadModel,
        s: &mut DagScratch,
    ) -> Option<(PathScore<Units>, Vec<CalledSegment>)> {
        self.layout(placements, keep, read, s);
        let mut best: Option<(PathScore<Units>, Cell)> = None;
        for e in 0..s.lattice.nodes.len() {
            self.gather(read, s, e);
            for t in s.lattice.nodes[e].slots() {
                if let Some((score, cell)) = self.emit(placements, read, s, t)
                    && best.is_none_or(|(current, _)| score > current)
                {
                    best = Some((score, cell));
                }
            }
        }
        let (score, end) = best?;
        Some((score, self.trace(placements, s, end)))
    }

    /// The fixed-scale objective of a path's `segments` (non-empty, presence
    /// costs `⌊P*_k⌋`), with its D count and germline support.
    pub fn path_score(&self, placements: &[Placement], segments: &[CalledSegment]) -> PathScore {
        let kind = |segment: &CalledSegment| placements[segment.placement].kind;
        let (first, last) = (kind(&segments[0]), kind(&segments[segments.len() - 1]));
        let (mut score, mut d_count, mut support) = (0, 0u8, 0);
        for (k, segment) in segments.iter().enumerate() {
            let placement = &placements[segment.placement];
            d_count += u8::from(placement.kind == SegmentKind::D);
            score += segment.span.score - self.model.presence_score(placement.kind);
            support += placement.score();
            if let Some(before) = k.checked_sub(1).map(|b| &segments[b]) {
                score +=
                    self.model
                        .junction(kind(before), placement.kind, segment.entry - before.exit);
            }
        }
        PathScore {
            score: score + self.model.end(first, last, d_count),
            d_count,
            support,
        }
    }

    /// The model the call maximises.
    pub fn calls(&self) -> CallModel {
        match self.mode {
            Mode::Read { calls, .. } => calls,
            Mode::Fixed => CallModel::FixedScale,
        }
    }

    /// `P*_k` of each rearranged kind.
    fn presence(&self) -> [f64; 3] {
        SegmentKind::REARRANGED.map(|kind| self.model.presence(kind))
    }

    /// The model a call at divergence `p` maximises.
    fn read_model(&self, p: f64) -> ReadModel {
        match self.mode {
            Mode::Read {
                evidence,
                calls: CallModel::PerRead,
                ..
            } => ReadModel::at(&evidence, p, self.presence()),
            // One unit per score unit, λ nats each (1 without an evidence
            // model; nothing then reads them as nats).
            Mode::Read { evidence, .. } => ReadModel::fixed(evidence.lambda, self.presence()),
            Mode::Fixed => ReadModel::fixed(1.0, self.presence()),
        }
    }

    /// The posterior's model at divergence `p`, when the DAG reports
    /// posteriors: the call's own under [`CallModel::PerRead`].
    fn nats(&self, p: f64) -> Option<PathNats> {
        let Mode::Read {
            evidence,
            calls,
            posterior: Some(_),
        } = self.mode
        else {
            return None;
        };
        Some(match calls {
            CallModel::PerRead => self.read_model(p).nats(),
            CallModel::FixedScale => evidence.nats(p, self.presence()),
        })
    }
}

impl<G: Grammar, M: PathModel> Adjudicator for PlacementDag<G, M> {
    type Scratch = DagScratch;

    fn adjudicate(&self, placements: &[Placement], scratch: &mut DagScratch) -> Option<CallPath> {
        let divergence = divergence(placements);
        let (_, segments) = self.viterbi(placements, &self.read_model(divergence), scratch)?;
        Some(CallPath {
            score: self.path_score(placements, &segments),
            segments,
            divergence,
            model: self.calls(),
        })
    }

    fn marginals(
        &self,
        placements: &[Placement],
        path: &CallPath,
        scratch: &mut DagScratch,
    ) -> Option<Marginals> {
        let nats = self.nats(path.divergence)?;
        self.posterior(placements, &nats, &mut scratch.posterior)
    }

    /// A fixed-scale call without posteriors never takes a D scoring at
    /// most `⌊P*_D⌋` ([`PathModel::d_floor`]); a per-read call can, and a
    /// posterior gives every D its probability.
    fn d_score_floor(&self) -> Score {
        match self.mode {
            Mode::Read {
                calls: CallModel::PerRead,
                ..
            }
            | Mode::Read {
                posterior: Some(_), ..
            } => Score::MIN,
            _ => self.model.d_floor(),
        }
    }

    /// Under [`DAlleleDecoding::MaxMarginal`], a call with one D whose
    /// allele is not the max-marginal one is re-traced with the call's
    /// topology through the D placements of that allele.
    fn decode(
        &self,
        placements: &[Placement],
        path: &CallPath,
        marginals: &Marginals,
        scratch: &mut DagScratch,
    ) -> Option<CallPath> {
        let Mode::Read {
            posterior: Some(DAlleleDecoding::MaxMarginal),
            ..
        } = self.mode
        else {
            return None;
        };
        let kind = |segment: &CalledSegment| placements[segment.placement].kind;
        let mut ds = path.segments.iter().filter(|s| kind(s) == SegmentKind::D);
        let (Some(called), None) = (ds.next(), ds.next()) else {
            return None;
        };
        let allele = max_marginal_d(placements, &placements[called.placement], marginals)?;
        let topology = Topology {
            first: kind(&path.segments[0]),
            last: kind(&path.segments[path.segments.len() - 1]),
            d_count: 1,
        };
        let keep = |i: usize| {
            let p = &placements[i];
            p.kind != SegmentKind::D || p.alleles.contains(&allele)
        };
        let read = self.read_model(path.divergence);
        let (_, mut segments) = self.conditioned(placements, keep, topology, &read, scratch)?;
        for segment in &mut segments {
            segment.allele = (kind(segment) == SegmentKind::D).then_some(allele);
        }
        Some(CallPath {
            score: self.path_score(placements, &segments),
            segments,
            ..*path
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::germline::AlleleId;
    use crate::testing::placement;

    /// V's tail and J's head overlap on six bases, V mismatching at the 2nd and
    /// 3rd and J at the 4th and 5th. One shared cut must give one side the
    /// mismatches; two cuts leave bases 2–5 to N, worth 2 points more.
    #[test]
    fn two_cuts_leave_a_shared_mismatch_stretch_to_n() {
        let query = encode(b"TTTTTTTTTTACGTACCCCCCCCCCC");
        let v = encode(b"TTTTTTTTTTAGTTAC"); //      overlap ACGTAC: M X X M M M
        let j = encode(b"ACGAGCCCCCCCCCCC"); // overlap ACGTAC: M M M X X M
        let placements = [
            placement(SegmentKind::V, &query, &v, 0),
            placement(SegmentKind::J, &query, &j, 10),
        ];
        assert_eq!((placements[0].score(), placements[1].score()), (22, 22));
        let dag = PlacementDag::new(
            VdjGrammar::VJ,
            EngineeringObjective::priced([0.0, 11.0, 0.0]),
            [1, 1, 1],
            Mode::Fixed,
        );
        let path = dag
            .adjudicate(&placements, &mut DagScratch::default())
            .expect("a VJ path");
        // Single cut: 40 + 2 = 42. Two cuts (V to 11, J from 15): 44.
        assert_eq!(path.score.score, 44);
        let (v_seg, j_seg) = (path.segments[0], path.segments[1]);
        assert_eq!((v_seg.span.q_hi, j_seg.span.q_lo), (11, 15));
    }

    /// C is assigned after adjudication and never enters a lattice, whose
    /// per-kind arrays cover the rearranged kinds only.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "a C placement reached the lattice")]
    fn a_c_placement_never_enters_the_lattice() {
        let query = encode(b"TTTTTTTTTTACGTACCCCCCCCCCCGGGGGGGGGG");
        let placements = [
            placement(SegmentKind::V, &query, &query[..16], 0),
            placement(SegmentKind::J, &query, &query[10..26], 10),
            placement(SegmentKind::C, &query, &query[26..], 26),
        ];
        let dag = PlacementDag::new(
            VdjGrammar::VJ,
            EngineeringObjective::priced([0.0, 11.0, 0.0]),
            [1, 1, 1],
            Mode::Fixed,
        );
        dag.adjudicate(&placements, &mut DagScratch::default());
    }

    /// The J alignment starts one base before the D, its first 13 bases
    /// mismatching the D region. V → D → J, with J entered after the D,
    /// beats V → J: the D follows a J that starts earlier.
    #[test]
    fn a_d_can_precede_a_j_that_starts_earlier() {
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
        let placements = [
            placement(SegmentKind::V, &query, &encode(&ascii[..20]), 0),
            placement(SegmentKind::D, &query, &encode(&ascii[22..34]), 22),
            placement(SegmentKind::J, &query, &encode(&j), 21),
        ];
        assert_eq!(
            placements.iter().map(Placement::score).collect::<Vec<_>>(),
            [40, 24, 7]
        );
        let dag = PlacementDag::new(
            VdjGrammar::VDJ,
            EngineeringObjective::priced([0.0, 11.0, 0.0]),
            [1, 1, 1],
            Mode::Fixed,
        );
        let path = dag
            .adjudicate(&placements, &mut DagScratch::default())
            .expect("a path");
        // V 40 + D (24 − 11) + J from base 34 (46) = 99; V → J alone is 86.
        assert_eq!(path.score.score, 99);
        let spans: Vec<(u32, u32)> = path
            .segments
            .iter()
            .map(|s| (s.span.q_lo, s.span.q_hi))
            .collect();
        assert_eq!(spans, vec![(0, 20), (22, 34), (34, 57)]);
    }

    /// J starts inside V, and its alignment has two mismatches right after
    /// V's end. Keeping V whole and entering J after the mismatches (N/P
    /// between) beats every cut inside the overlap: the entry lies beyond
    /// the overlap, in the middle of J's alignment.
    #[test]
    fn entries_beyond_the_overlap_are_tried() {
        let query = encode(b"TTTTTTTTTTTTTTTTACGTCCAAAAAAAAAAAAAAAAAA");
        let v = encode(b"TTTTTTTTTTTTTTTTACGT");
        let j = encode(b"ACGTGGAAAAAAAAAAAAAAAAAA");
        let placements = [
            placement(SegmentKind::V, &query, &v, 0),
            placement(SegmentKind::J, &query, &j, 16),
        ];
        assert_eq!((placements[0].score(), placements[1].score()), (40, 38));
        let dag = PlacementDag::new(
            VdjGrammar::VJ,
            EngineeringObjective::priced([0.0, 11.0, 0.0]),
            [1, 1, 1],
            Mode::Fixed,
        );
        let path = dag
            .adjudicate(&placements, &mut DagScratch::default())
            .expect("a VJ path");
        // Inside the overlap the best is V[0,16) + J whole = 32 + 38 = 70.
        assert_eq!(path.score.score, 76);
        let spans: Vec<(u32, u32)> = path
            .segments
            .iter()
            .map(|s| (s.span.q_lo, s.span.q_hi))
            .collect();
        assert_eq!(spans, vec![(0, 20), (22, 40)]);
    }

    /// A V of 20 exact bases, then a J whose alignment matches the read only
    /// in its last 3 bases. The J ends the path, so its clip keeps its end:
    /// those 3 bases, worth 6 points or `3·a(p̂)`. Free, they complete the
    /// path; at `P*_J = 7.91` (exact runs of 4 win) they do not pay for a J,
    /// on the fixed scale or at any divergence.
    #[test]
    fn a_j_pays_its_presence_cost() {
        let ascii = b"TTTTTTTTTTTTTTTTTTTTACGTACGTACGTCCA";
        let query = encode(ascii);
        let mut j = encode(&ascii[20..]);
        for base in &mut j[..12] {
            *base = (*base + 1) & 3;
        }
        let placements = [
            placement(SegmentKind::V, &query, &query[..20], 0),
            placement(SegmentKind::J, &query, &j, 20),
        ];
        let kinds = |presence: [f64; 3], read: &ReadModel| {
            let model = EngineeringObjective::priced(presence);
            let dag = PlacementDag::new(VdjGrammar::VJ, model, [1, 1, 1], Mode::Fixed);
            let (_, segments) = dag
                .viterbi(&placements, read, &mut DagScratch::default())
                .expect("a path");
            segments.len()
        };
        let (free, priced) = ([0.0; 3], [0.0, 0.0, 7.91]);
        for p in [None, Some(0.001), Some(0.1), Some(0.4)] {
            let read = |presence| match p {
                None => ReadModel::fixed(1.0, presence),
                Some(p) => ReadModel::at(&evidence(), p, presence),
            };
            assert_eq!(kinds(free, &read(free)), 2, "{p:?}");
            assert_eq!(kinds(priced, &read(priced)), 1, "{p:?}");
        }
    }

    /// `(2, −3)` scored at the KA λ0 for every kind.
    fn evidence() -> EvidenceModel {
        let scoring = crate::scoring::AffineScoring::new(2, -3, -9, -1);
        EvidenceModel::new([&scoring; 3]).expect("λ")
    }

    /// A read with two D candidates in its N region: an exact 7-mer (score
    /// 14, placement 1) and a 10-mer with a mismatch in its middle (score 15,
    /// placement 2), whose best clip is 5 bases.
    fn two_ds() -> (Vec<u8>, [Placement; 4]) {
        let v = b"TTTTTTTTTTTTTTTTTTTT";
        let d1 = b"GACCTAG";
        let d2 = b"CAGTCTTGCA";
        let j = b"AAAAAAAAAAAAAAAAAAAA";
        let mut ascii = v.to_vec();
        ascii.extend_from_slice(b"CC");
        ascii.extend_from_slice(d1);
        ascii.extend_from_slice(b"GG");
        ascii.extend_from_slice(d2);
        ascii.extend_from_slice(b"CC");
        ascii.extend_from_slice(j);
        let query = encode(&ascii);
        let mut d2_reference = encode(d2);
        d2_reference[5] = (d2_reference[5] + 1) & 3;
        let mut placements = [
            placement(SegmentKind::V, &query, &encode(v), 0),
            placement(SegmentKind::D, &query, &encode(d1), 22),
            placement(SegmentKind::D, &query, &d2_reference, 31),
            placement(SegmentKind::J, &query, &encode(j), 43),
        ];
        placements[2].alleles = vec![1];
        assert_eq!((placements[1].score(), placements[2].score()), (14, 15));
        (query, placements)
    }

    /// The fixed scale takes the 10-mer of [`two_ds`]. At low divergence a
    /// mismatch costs `ln(4p̂/3)` (−4.1 nats at `p̂ = 1.2 %`, not −1.9), so the
    /// read's own model takes the exact 7-mer: 9.6 nats against 8.3.
    #[test]
    fn a_per_read_call_takes_the_exact_short_d_at_low_divergence() {
        let (_, placements) = two_ds();
        let model = EngineeringObjective::priced([0.0, 10.0, 0.0]);
        let dag = |calls| {
            let mode = Mode::Read {
                evidence: evidence(),
                calls,
                posterior: None,
            };
            PlacementDag::new(VdjGrammar::VDJ, model, [1, 4, 1], mode)
        };
        let (fixed, per_read) = (dag(CallModel::FixedScale), dag(CallModel::PerRead));
        let mut scratch = DagScratch::default();
        let d_of = |path: &CallPath| path.segments[1].placement;
        let call = fixed.adjudicate(&placements, &mut scratch).unwrap();
        assert_eq!((d_of(&call), call.model), (2, CallModel::FixedScale));
        let call = per_read.adjudicate(&placements, &mut scratch).unwrap();
        assert_eq!(call.divergence, 0.5 / 41.0);
        assert_eq!((d_of(&call), call.model), (1, CallModel::PerRead));
        // The reported path score stays the fixed scale's: 40 + 14 − 10 + 40.
        assert_eq!(call.score.score, 84);
        // At δ0 the per-read model is the fixed one.
        let delta0 = evidence().log_odds_divergence(SegmentKind::D);
        let at = ReadModel::at(&evidence(), delta0, model.presence);
        let (_, segments) = per_read.viterbi(&placements, &at, &mut scratch).unwrap();
        assert_eq!(segments[1].placement, 2);
    }

    /// On the fixed scale the call takes the 10-mer of [`two_ds`], but the
    /// posterior at the read's divergence puts more mass on the exact 7-mer:
    /// decoding its max-marginal allele re-traces the path through it, with
    /// the call's topology, V and J.
    #[test]
    fn the_max_marginal_d_allele_is_retraced() {
        let (_, placements) = two_ds();
        let dag = |d_allele| {
            let mode = Mode::Read {
                evidence: evidence(),
                calls: CallModel::FixedScale,
                posterior: Some(d_allele),
            };
            let model = EngineeringObjective::priced([0.0, 10.0, 0.0]);
            PlacementDag::new(VdjGrammar::VDJ, model, [1, 4, 1], mode)
        };
        let (map, dag) = (dag(DAlleleDecoding::Map), dag(DAlleleDecoding::MaxMarginal));
        let mut scratch = DagScratch::default();
        let call = dag.adjudicate(&placements, &mut scratch).unwrap();
        assert_eq!(call.segments[1].placement, 2);
        let m = dag.marginals(&placements, &call, &mut scratch).unwrap();
        assert!(m.placement[1] > m.placement[2], "{:?}", m.placement);
        let decoded = dag.decode(&placements, &call, &m, &mut scratch).unwrap();
        let picks = |path: &CallPath| -> Vec<(usize, Option<AlleleId>)> {
            path.segments
                .iter()
                .map(|s| (s.placement, s.allele))
                .collect()
        };
        assert_eq!(picks(&decoded), [(0, None), (1, Some(0)), (3, None)]);
        assert_eq!(decoded.score.score, 84);
        // The call's own allele is kept when it is the max-marginal one.
        assert!(map.decode(&placements, &call, &m, &mut scratch).is_none());
        assert!(
            dag.decode(&placements, &decoded, &m, &mut scratch)
                .is_none()
        );
    }

    /// At high divergence a D scoring `P_D` can win a per-read call: 8
    /// matches and 2 mismatches score 10, no better than `P_D` = 10, but are
    /// worth 6.40 nats at `p̂ = 0.3`, above `c_D` = 5.15. So a per-read call
    /// keeps every admissible D, below the fixed-scale floor too.
    #[test]
    fn a_d_at_the_fixed_floor_can_win_a_per_read_call() {
        let v = b"TTTTTTTTTTTTTTTTTTTT";
        let d = b"GACCTAGCAG";
        let j = b"AAAAAAAAAAAAAAAAAAAA";
        let mut ascii = v.to_vec();
        ascii.extend_from_slice(b"CC");
        ascii.extend_from_slice(d);
        ascii.extend_from_slice(b"CC");
        ascii.extend_from_slice(j);
        let query = encode(&ascii);
        let mut reference = encode(d);
        for k in [3, 6] {
            reference[k] = (reference[k] + 1) & 3;
        }
        let placements = [
            placement(SegmentKind::V, &query, &encode(v), 0),
            placement(SegmentKind::D, &query, &reference, 22),
            placement(SegmentKind::J, &query, &encode(j), 34),
        ];
        let model = EngineeringObjective::priced([0.0, 10.0, 0.0]);
        assert_eq!(placements[1].score(), model.d_floor() - 1);
        let dag = PlacementDag::new(VdjGrammar::VDJ, model, [1, 4, 1], Mode::Fixed);
        let mut scratch = DagScratch::default();
        let d_count = |read: &ReadModel, scratch: &mut DagScratch| {
            dag.viterbi(&placements, read, scratch).unwrap().0.d_count
        };
        assert_eq!(
            d_count(&ReadModel::fixed(1.0, model.presence), &mut scratch),
            0
        );
        let high = ReadModel::at(&evidence(), 0.3, model.presence);
        assert_eq!(d_count(&high, &mut scratch), 1);
    }
}
