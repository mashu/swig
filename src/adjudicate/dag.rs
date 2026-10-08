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
//! the full segmentation model: when two placements overlap, AER-R's
//! candidate A (D wins the overlap), candidate B (V/J win it), every
//! allocation in between, and leaving any part of the overlap or of either
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
//! `O(overlap + 1)`, a read `O(Σ_nodes len + Σ_edges (overlap + 1))`.
//!
//! Ties: final paths compare by [`PathScore`] (value, then fewer D, then
//! germline support, the summed unclipped fixed scores), so a with-D path must
//! strictly beat the best no-D path. Integer values make every tie exact.
//! Among equal keys: a node keeps its earliest entry, so a junction base that
//! scores the same on either side goes to the later segment; for that entry,
//! the predecessor's latest exit; and across predecessors, a path start, then
//! the first predecessor in DAG order.

use super::lattice::{Group, KINDS, Lattice};
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
    /// Minimum query bases a clipped segment must keep, indexed by [`SegmentKind::index`].
    pub min_retained: [u32; 3],
    /// How alignments read as log-odds at a read's divergence. Without it
    /// everything stays on the fixed scale and there are no posteriors.
    pub evidence: Option<EvidenceModel>,
    /// The model the call maximises.
    pub calls: CallModel,
    /// Whether [`Adjudicator::marginals`] returns posterior marginals (see
    /// [`PlacementDag::posterior`]); needs `evidence`.
    pub posterior: bool,
}

impl Default for PlacementDag {
    fn default() -> Self {
        Self::new(VdjGrammar::VDJ, EngineeringObjective::default(), [10, 3, 6])
    }
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
    /// A DAG on the fixed scale, without posteriors.
    pub fn new(grammar: G, model: M, min_retained: [u32; 3]) -> Self {
        Self {
            grammar,
            model,
            min_retained,
            evidence: None,
            calls: CallModel::FixedScale,
            posterior: false,
        }
    }

    #[inline]
    pub(super) fn min_bases(&self, kind: SegmentKind) -> u32 {
        self.min_retained[kind.index()].max(1)
    }

    /// The lattice with the prefix value of every canonical cut under `read`.
    fn layout(&self, placements: &[Placement], read: &ReadModel, s: &mut DagScratch) {
        s.lattice.layout(placements);
        s.entry.clear();
        s.entry.resize(s.lattice.cells(), None);
        s.exit.clear();
        s.exit.resize(s.lattice.cells(), None);
        for node in &s.lattice.nodes {
            // (exit prefix score, exit prefix value) of the previous column.
            let mut previous: Option<(Score, Units)> = None;
            for (q, entry, exit) in placements[node.placement].clip.matches() {
                let before = match previous {
                    None => read.prior(entry),
                    Some((score, value)) => value + read.prior(entry - score),
                };
                let through = before + read.column(node.kind, exit - entry);
                let o = (q - node.q_start) as usize;
                s.entry[node.cuts + o] = Some(before);
                s.exit[node.cuts + o + 1] = Some(through);
                previous = Some((exit, through));
            }
        }
        s.open.clear();
        s.forward.clear();
        s.exits.clear();
        s.running.clear();
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
        s.running.resize(KINDS * s.lattice.values(), None);

        let entries = s.lattice.entries(&node);
        let open = &s.open;
        let live = |u: usize| open[u];
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
            s.tails.clear();
            s.tails.resize(KINDS * width, None);
            for (b, u) in s.lattice.incoming(e, grammar, slot.group, &live) {
                let iota = read.prior(self.model.junction_per_base(b.kind, node.kind));
                let running = &s.running[s.lattice.slots[u].per_kind(node.kind, b.width())];
                // Entries inside b: the best exit at or before the entry. A
                // predecessor of another kind may start after this node
                // (its entries before b's start have no exit to follow).
                for &o in entries {
                    let y = node.q_start + o;
                    if y >= b.q_end() {
                        break;
                    }
                    let Some(z) = y.checked_sub(b.q_start) else {
                        continue;
                    };
                    if let Some(cell) = running[z as usize] {
                        let candidate = Cell {
                            key: (cell.key.0 + iota * Units::from(z), cell.key.1),
                            slot: u as u32,
                            at: cell.at,
                        };
                        prefer(&mut forward[o as usize], candidate);
                    }
                }
                // Entries at or beyond b's end: b's best exit overall, once.
                let beyond = b.q_end().saturating_sub(node.q_start) as usize;
                if let (true, Some(cell)) = (beyond < width, running[b.len as usize]) {
                    let candidate = Cell {
                        key: (
                            cell.key.0
                                + iota * (Units::from(node.q_start) - Units::from(b.q_start)),
                            cell.key.1,
                        ),
                        slot: u as u32,
                        at: cell.at,
                    };
                    prefer(&mut s.tails[b.kind.index() * width + beyond], candidate);
                }
            }
            for kind in SegmentKind::ALL {
                let tail = &s.tails[kind.index() * width..(kind.index() + 1) * width];
                if tail.iter().all(Option::is_none) {
                    continue;
                }
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
        let presence = if node.kind == SegmentKind::D {
            read.d_cost()
        } else {
            0
        };
        let bonus = read.prior(self.model.segment_bonus(placement, d_used)) - presence;
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
        for kind in SegmentKind::ALL {
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
        self.grammar.can_end(first, node.kind).then(|| {
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
    pub fn viterbi(
        &self,
        placements: &[Placement],
        read: &ReadModel,
        s: &mut DagScratch,
    ) -> Option<(PathScore<Units>, Vec<CalledSegment>)> {
        self.layout(placements, read, s);
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

    /// The fixed-scale objective of a path's `segments` (non-empty, D cost
    /// `⌊P*_D⌋`), with its D count and germline support.
    pub fn path_score(&self, placements: &[Placement], segments: &[CalledSegment]) -> PathScore {
        let kind = |segment: &CalledSegment| placements[segment.placement].kind;
        let (first, last) = (kind(&segments[0]), kind(&segments[segments.len() - 1]));
        let (mut score, mut d_count, mut support) = (0, 0u8, 0);
        for (k, segment) in segments.iter().enumerate() {
            let placement = &placements[segment.placement];
            if placement.kind == SegmentKind::D {
                d_count += 1;
                score -= self.model.d_presence_score();
            }
            score += segment.span.score + self.model.segment_bonus(placement, d_count);
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

    /// The model a call at divergence `p` maximises.
    fn read_model(&self, p: f64) -> ReadModel {
        let d_presence = self.model.d_presence();
        match (self.calls, &self.evidence) {
            (CallModel::PerRead, Some(evidence)) => ReadModel::at(evidence, p, d_presence),
            // Without an evidence model the fixed scores are not read as nats.
            (_, evidence) => ReadModel::fixed(evidence.map_or(1.0, |e| e.lambda), d_presence),
        }
    }

    /// The posterior's model at divergence `p`: the call's own under
    /// [`CallModel::PerRead`].
    fn nats(&self, p: f64) -> Option<PathNats> {
        let evidence = self.evidence.as_ref()?;
        Some(match self.calls {
            CallModel::PerRead => self.read_model(p).nats(),
            CallModel::FixedScale => evidence.nats(p, self.model.d_presence()),
        })
    }
}

impl<G: Grammar, M: PathModel> Adjudicator for PlacementDag<G, M> {
    type Scratch = DagScratch;

    fn adjudicate(
        &self,
        _query_len: usize,
        placements: &[Placement],
        scratch: &mut DagScratch,
    ) -> Option<CallPath> {
        let divergence = divergence(placements);
        let (_, segments) = self.viterbi(placements, &self.read_model(divergence), scratch)?;
        Some(CallPath {
            score: self.path_score(placements, &segments),
            segments,
            divergence,
        })
    }

    fn marginals(
        &self,
        placements: &[Placement],
        path: &CallPath,
        scratch: &mut DagScratch,
    ) -> Option<Marginals> {
        let nats = self.nats(path.divergence).filter(|_| self.posterior)?;
        self.posterior(placements, &nats, &mut scratch.posterior)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
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
        let dag = PlacementDag::new(VdjGrammar::VJ, EngineeringObjective::default(), [1, 1, 1]);
        let path = dag
            .adjudicate(query.len(), &placements, &mut DagScratch::default())
            .expect("a VJ path");
        // Single cut: 40 + 2 = 42. Two cuts (V to 11, J from 15): 44.
        assert_eq!(path.score.score, 44);
        let (v_seg, j_seg) = (path.segments[0], path.segments[1]);
        assert_eq!((v_seg.span.q_hi, j_seg.span.q_lo), (11, 15));
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
        let dag = PlacementDag::new(VdjGrammar::VDJ, EngineeringObjective::default(), [1, 1, 1]);
        let path = dag
            .adjudicate(query.len(), &placements, &mut DagScratch::default())
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
        let dag = PlacementDag::new(VdjGrammar::VJ, EngineeringObjective::default(), [1, 1, 1]);
        let path = dag
            .adjudicate(query.len(), &placements, &mut DagScratch::default())
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
}
