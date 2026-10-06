//! Placement DAG: bounded-topology semi-Markov Viterbi over clipped placements.
//!
//! Nodes are placements sorted by segment kind (V, D, J), then query start.
//! A path visits placements in read order of their retained bases: a
//! transition between kinds needs only an exit cut `x` and an entry cut
//! `y ≥ x`, wherever the two alignments start, so a D can precede a J whose
//! alignment reaches back over it. Two segments of one kind (tandem D) also
//! follow the order of their alignment starts. Each placement on the path is
//! clipped independently: the
//! earlier of two consecutive placements keeps query bases before its **exit
//! cut** `x`, the later one keeps bases from its **entry cut** `y ≥ x` on,
//! and the bases in between are N/P. This is the full segmentation model:
//! when two placements overlap, AER-R's candidate A (D wins the overlap),
//! candidate B (V/J win it), every allocation in between, and leaving any part
//! of the overlap or of either alignment to N/P are all evaluated exactly.
//! The first segment keeps its start and the last its end (V and J are
//! extended to the read's ends on purpose, see [`crate::refine::extend_outer`]).
//!
//! Cuts are *canonical* ([`crate::clip::ClipIndex::matches`]): an entry on
//! the first retained base, an exit right after the last, so every distinct
//! clipped span is one `(entry, exit)` pair.
//!
//! **Slots.** The paths through a node with a given D count and first kind
//! form a slot, with one value per entry cut: the best key (score, then
//! germline support) of a path prefix entering there. The recursion is the
//! max-plus twin of the posterior's sum-product ([`super::posterior`]):
//!
//! ```text
//! X(b, x)  = exit(x) + bonus + max_{y ≤ x − min} [F(b, y) − entry(y)]
//! F(e, y') = max_b max_{x ≤ y'} [X(b, x) + ι·(y' − x)]          (0 at a path start)
//! best     = max_b X(b, end of b) + end(g, b)
//! ```
//!
//! **Trick — separable scores → running maxima.** A clipped span scores
//! `exit(x).prefix − entry(y).prefix` (see [`crate::clip`]), the model adds a
//! cut-independent bonus, and the N/P score is linear in the inserted length.
//! So the inner maximum over entries is a running maximum along the node, and
//! the maximum over exits `x ≤ y'` another, read in O(1) per entry. On an edge
//! `b → e`, every entry of `e` at or beyond `b`'s end sees the same best exit,
//! so it is added once to a running maximum over `e`'s entries: an edge costs
//! `O(overlap + 1)`, a read `O(Σ_nodes len + Σ_edges (overlap + 1))`.
//!
//! Ties: final paths compare by [`PathScore`] (score, then fewer D, then
//! germline support), so a with-D path must strictly beat the best no-D path.
//! Among equal keys: a node keeps its earliest entry, so a junction base that
//! scores the same on either side goes to the later segment; for that entry,
//! the predecessor's latest exit; and across predecessors, a path start, then
//! the first predecessor in DAG order.

use super::model::{Grammar, PathModel, VdjGrammar};
use super::posterior::{divergence, Marginals, PosteriorModel, PosteriorScratch};
use super::{Adjudicator, CallPath, CalledSegment, EngineeringObjective, PathScore};
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// Segment kinds: one running-maximum array per kind and slot.
const KINDS: usize = SegmentKind::ALL.len();

/// No back pointer: a path start.
const START: u32 = u32::MAX;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlacementDag<G = VdjGrammar, M = EngineeringObjective> {
    pub grammar: G,
    pub model: M,
    /// Minimum query bases a clipped segment must keep, indexed by [`SegmentKind::index`].
    pub min_retained: [u32; 3],
    /// When set, [`Adjudicator::marginals`] returns posterior marginals
    /// (see [`PlacementDag::posterior`]).
    pub posterior: Option<PosteriorModel>,
}

impl Default for PlacementDag {
    fn default() -> Self {
        Self {
            grammar: VdjGrammar::VDJ,
            model: EngineeringObjective::default(),
            min_retained: [10, 3, 6],
            posterior: None,
        }
    }
}

/// Score, then germline support (the summed unclipped scores of the path's
/// placements), compared lexicographically.
type Key = (Score, Score);

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
    fn shift(self, by: Score) -> Self {
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

/// DAG order: segment kind (V, D, J), then query start and end. Every
/// transition between kinds goes forward in this order whatever the
/// alignments' starts, so a D can precede a J whose alignment starts earlier
/// (the cuts decide whether the path is feasible); two segments of one kind
/// (tandem D) follow the order of their starts.
pub(super) fn order_key(p: &Placement, i: usize) -> (usize, u32, u32, usize) {
    (p.kind.index(), p.q_start(), p.q_end(), i)
}

/// One placement in DAG order (kind, then query start and end).
#[derive(Clone, Copy, Debug)]
struct Node {
    placement: usize,
    kind: SegmentKind,
    q_start: u32,
    /// `q_end − q_start`; cut offsets run over `0..=len`.
    len: u32,
    /// Start of the node's dense cut prefixes in `entry` / `exit`.
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

    /// Whether `self` comes before `other` in DAG order and may precede it
    /// on a path: an earlier kind, or the same kind starting earlier.
    fn precedes(&self, other: &Node) -> bool {
        (self.kind.index(), self.q_start) < (other.kind.index(), other.q_start)
    }

    fn width(&self) -> usize {
        self.len as usize + 1
    }
}

/// The paths through one node with one `(D count, first kind)`.
#[derive(Clone, Copy, Debug)]
struct Slot {
    node: u32,
    d_used: u8,
    first: SegmentKind,
    /// Start of the slot's `width` cells in `forward` / `exits`;
    /// `KINDS · values` starts its per-successor-kind arrays in `running`.
    values: usize,
    /// Whether any path leaves the slot.
    open: bool,
}

impl Slot {
    fn per_kind(&self, kind: SegmentKind, width: usize) -> std::ops::Range<usize> {
        let start = KINDS * self.values + kind.index() * width;
        start..start + width
    }
}

/// Per-thread buffers for [`PlacementDag`].
#[derive(Debug, Default)]
pub struct DagScratch {
    order: Vec<usize>,
    nodes: Vec<Node>,
    /// Prefix score at each canonical entry / exit offset (`None` elsewhere).
    entry: Vec<Option<Score>>,
    exit: Vec<Option<Score>>,
    entries: Vec<u32>,
    slots: Vec<Slot>,
    /// Best path prefix entering at each offset.
    forward: Vec<Option<Cell>>,
    /// Best path prefix through each exit offset (segment included).
    exits: Vec<Option<Cell>>,
    /// Per slot and successor kind: best `X(x) − ι·x` over exits `≤ z`.
    running: Vec<Option<Cell>>,
    groups: Vec<(u8, SegmentKind)>,
    tails: Vec<Option<Cell>>,
    posterior: PosteriorScratch,
}

/// Slots of `predecessors` whose paths continue into the group
/// `(d_used, first)` of a node of kind `kind`.
fn incoming<'s, G: Grammar>(
    slots: &'s [Slot],
    grammar: &'s G,
    predecessors: &'s [Node],
    kind: SegmentKind,
    d_used: u8,
    first: SegmentKind,
) -> impl Iterator<Item = (Node, usize)> + 's {
    let d_before = d_used.checked_sub(u8::from(kind == SegmentKind::D));
    predecessors.iter().flat_map(move |b| {
        (b.slots.0..b.slots.1)
            .filter(move |&u| {
                let slot = &slots[u];
                slot.open
                    && slot.first == first
                    && Some(slot.d_used) == d_before
                    && grammar.can_follow(b.kind, kind, slot.d_used)
            })
            .map(move |u| (*b, u))
    })
}

impl<G: Grammar, M: PathModel> PlacementDag<G, M> {
    pub fn new(grammar: G, model: M, min_retained: [u32; 3]) -> Self {
        Self {
            grammar,
            model,
            min_retained,
            posterior: None,
        }
    }

    /// Also compute posterior marginals under `model`.
    pub fn with_posterior(mut self, model: PosteriorModel) -> Self {
        self.posterior = Some(model);
        self
    }

    #[inline]
    pub(super) fn min_bases(&self, kind: SegmentKind) -> u32 {
        self.min_retained[kind.index()].max(1)
    }

    /// Nodes in DAG order with dense canonical-cut prefixes.
    fn layout(&self, placements: &[Placement], s: &mut DagScratch) {
        s.order.clear();
        s.order.extend(0..placements.len());
        s.order.sort_by_key(|&i| order_key(&placements[i], i));
        s.nodes.clear();
        s.entry.clear();
        s.exit.clear();
        s.entries.clear();
        for &i in &s.order {
            let p = &placements[i];
            let (q_start, len) = (p.q_start(), p.q_end() - p.q_start());
            let cuts = s.entry.len();
            s.entry.resize(cuts + len as usize + 1, None);
            s.exit.resize(cuts + len as usize + 1, None);
            let first_entry = s.entries.len();
            for (q, entry, exit) in p.clip.matches() {
                let o = q - q_start;
                s.entry[cuts + o as usize] = Some(entry);
                s.exit[cuts + o as usize + 1] = Some(exit);
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
    }

    /// Create node `e`'s slots and their entry values from path starts and
    /// from every predecessor's running exit maxima.
    fn gather(&self, s: &mut DagScratch, e: usize) {
        let grammar = &self.grammar;
        let node = s.nodes[e];
        let width = node.width();
        let predecessors = s.nodes.partition_point(|b| b.precedes(&node));
        let start = grammar
            .can_start(node.kind)
            .then_some((u8::from(node.kind == SegmentKind::D), node.kind))
            .filter(|_| s.entry[node.cuts].is_some());

        s.groups.clear();
        s.groups.extend(start);
        for b in &s.nodes[..predecessors] {
            for slot in &s.slots[b.slots.0..b.slots.1] {
                let group = (
                    slot.d_used + u8::from(node.kind == SegmentKind::D),
                    slot.first,
                );
                if slot.open
                    && grammar.can_follow(b.kind, node.kind, slot.d_used)
                    && !s.groups.contains(&group)
                {
                    s.groups.push(group);
                }
            }
        }
        let first_slot = s.slots.len();
        for &(d_used, first) in &s.groups {
            let values = s.forward.len();
            s.forward.resize(values + width, None);
            s.exits.resize(values + width, None);
            s.running.resize(KINDS * (values + width), None);
            s.slots.push(Slot {
                node: e as u32,
                d_used,
                first,
                values,
                open: false,
            });
        }
        s.nodes[e].slots = (first_slot, s.slots.len());

        let entries = &s.entries[node.entries.0..node.entries.1];
        for t in first_slot..s.slots.len() {
            let Slot {
                d_used,
                first,
                values,
                ..
            } = s.slots[t];
            let forward = &mut s.forward[values..values + width];
            if start == Some((d_used, first)) {
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
            for (b, u) in incoming(
                &s.slots,
                grammar,
                &s.nodes[..predecessors],
                node.kind,
                d_used,
                first,
            ) {
                let iota = self.model.junction_per_base(b.kind, node.kind);
                let running = &s.running[s.slots[u].per_kind(node.kind, b.width())];
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
                            key: (cell.key.0 + iota * z as Score, cell.key.1),
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
                            cell.key.0 + iota * (node.q_start as Score - b.q_start as Score),
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
                let iota = self.model.junction_per_base(kind, node.kind);
                let mut best: Option<Cell> = None;
                for (o, cell) in tail.iter().enumerate() {
                    if let Some(cell) = *cell {
                        prefer(&mut best, cell);
                    }
                    if let (Some(cell), true) = (best, s.entry[node.cuts + o].is_some()) {
                        prefer(&mut forward[o], cell.shift(iota * o as Score));
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
        s: &mut DagScratch,
        t: usize,
    ) -> Option<(PathScore, Cell)> {
        let slot = s.slots[t];
        let node = s.nodes[slot.node as usize];
        let width = node.width();
        let placement = &placements[node.placement];
        let min = self.min_bases(node.kind) as usize;
        let bonus = self.model.segment_bonus(placement, slot.d_used);
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
        s.slots[t].open = open;
        if !open {
            return None;
        }
        for kind in SegmentKind::ALL {
            if !self.grammar.can_follow(node.kind, kind, slot.d_used) {
                continue;
            }
            let iota = self.model.junction_per_base(node.kind, kind);
            let range = slot.per_kind(kind, width);
            let mut best: Option<Cell> = None;
            for o in 0..width {
                if let Some(cell) = s.exits[slot.values + o] {
                    // `ties`: the later exit retains more of the placement.
                    improve(
                        &mut best,
                        Cell {
                            key: (cell.key.0 - iota * o as Score, cell.key.1),
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
        self.grammar.can_end(slot.first, node.kind).then(|| {
            let score = PathScore {
                score: end.key.0 + self.model.end(slot.first, node.kind, slot.d_used),
                d_count: slot.d_used,
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
            let slot = s.slots[t];
            let node = s.nodes[slot.node as usize];
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

impl<G: Grammar, M: PathModel> Adjudicator for PlacementDag<G, M> {
    type Scratch = DagScratch;

    fn adjudicate(
        &self,
        _query_len: usize,
        placements: &[Placement],
        scratch: &mut DagScratch,
    ) -> Option<CallPath> {
        let s = scratch;
        self.layout(placements, s);
        s.slots.clear();
        s.forward.clear();
        s.exits.clear();
        s.running.clear();
        let mut best: Option<(PathScore, Cell)> = None;
        for e in 0..s.nodes.len() {
            self.gather(s, e);
            let (first, end) = s.nodes[e].slots;
            for t in first..end {
                if let Some((score, cell)) = self.emit(placements, s, t) {
                    if best.map_or(true, |(current, _)| score > current) {
                        best = Some((score, cell));
                    }
                }
            }
        }
        let (score, end) = best?;
        Some(CallPath {
            segments: self.trace(placements, s, end),
            score,
        })
    }

    fn marginals(
        &self,
        placements: &[Placement],
        path: &CallPath,
        scratch: &mut DagScratch,
    ) -> Option<Marginals> {
        let model = self.posterior?;
        let evidence = model.evidence(divergence(placements, path));
        self.posterior(placements, model.lambda, &evidence, &mut scratch.posterior)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Alignment, Op};
    use crate::alphabet::encode;
    use crate::clip::ClipIndex;
    use crate::evidence::Evidence;
    use crate::scoring::{AffineScoring, SubstitutionProfile};

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
