//! The placement lattice shared by the Viterbi pass ([`super::PlacementDag`])
//! and the forward–backward pass ([`super::posterior`]).
//!
//! **Nodes** are placements in DAG order: segment kind (V, D, J), then query
//! start and end. Every transition between kinds goes forward in this order
//! whatever the alignments' starts, so a D can precede a J whose alignment
//! starts earlier (the cuts decide whether a path is feasible); two segments
//! of one kind (tandem D) follow the order of their starts. A node has one
//! cell per cut offset `0..=len` and a list of canonical entry offsets
//! ([`crate::clip::ClipIndex::matches`]).
//!
//! **Slots.** The paths through a node with one [`Group`] — D count and
//! first segment kind — form a slot, which holds one value per cut offset.
//! A pass opens a node's slots when it reaches the node: the group of paths
//! starting there, then every group that a *live* slot of a predecessor
//! continues into. Liveness is the pass's own (some exit reachable), so
//! both passes build their slots by the same rule over the same domain.
//! Opening a slot records its incoming edges, the live predecessor slots
//! whose paths continue into it ([`Lattice::incoming`]).
//!
//! **Shared index arithmetic.** Every pass values a placement's cuts along
//! its columns ([`column_steps`]), and joins an edge `b → e` the same way: each
//! entry of `e` inside `b` meets `b`'s exits up to the same read position
//! ([`Lattice::inside`]), and every entry at or beyond `b`'s end meets all
//! of them, so the edge adds them once at its first such entry
//! ([`Node::beyond`]).

use std::ops::Range;

use super::model::Grammar;
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// Rearranged segment kinds: passes keep one array per successor or
/// predecessor kind. C never enters a lattice ([`SegmentKind::REARRANGED`]).
pub(super) const KINDS: usize = SegmentKind::REARRANGED.len();

/// One placement in DAG order.
#[derive(Clone, Copy, Debug)]
pub(super) struct Node {
    /// Index into the placement slice.
    pub placement: usize,
    pub kind: SegmentKind,
    pub q_start: u32,
    /// `q_end − q_start`; cut offsets run over `0..=len`.
    pub len: u32,
    /// Start of the node's [`Node::width`] cells in per-cut arrays.
    pub cuts: usize,
    /// The node's canonical entry offsets in the lattice's entry list.
    entries: (usize, usize),
    /// The node's slots in [`Lattice::slots`], once opened.
    pub slots: (usize, usize),
}

impl Node {
    pub fn q_end(&self) -> u32 {
        self.q_start + self.len
    }

    /// Whether `self` comes before `other` in DAG order and may precede it
    /// on a path: an earlier kind, or the same kind starting earlier.
    pub fn precedes(&self, other: &Node) -> bool {
        (self.kind.index(), self.q_start) < (other.kind.index(), other.q_start)
    }

    pub fn width(&self) -> usize {
        self.len as usize + 1
    }

    pub fn slots(&self) -> Range<usize> {
        self.slots.0..self.slots.1
    }

    /// On an edge from this node to `e`: the first offset of `e` at or
    /// beyond this node's end, from which every entry of `e` follows all of
    /// this node's exits; `None` when `e` ends before it.
    pub fn beyond(&self, e: &Node) -> Option<usize> {
        let beyond = self.q_end().saturating_sub(e.q_start) as usize;
        (beyond < e.width()).then_some(beyond)
    }
}

/// A placement's canonical cuts along its alignment: per `Match` column,
/// its read position `q`, its reference position `r`, the score of the gap
/// run since the previous `Match` column (or the alignment's start) and its
/// substitution score. A pass accumulates their values once per node: the
/// entry cut at `q` is worth the prefix through the gap run, the exit cut at
/// `q + 1` the prefix through the column. Canonical cuts never split a gap
/// run.
pub(super) fn column_steps(
    placement: &Placement,
) -> impl Iterator<Item = (u32, u32, Score, Score)> + '_ {
    let mut previous = 0;
    placement
        .clip
        .match_columns()
        .map(move |(q, r, entry, exit)| {
            let gap = entry - previous;
            previous = exit;
            (q, r, gap, exit - entry)
        })
}

/// Whether a path may step from `b` to `e`: a V that cannot open a D's
/// flank does not precede a D, and a J that cannot close one does not
/// follow a D (`closes`, per placement; empty: every one can).
#[inline]
pub(super) fn joins(closes: &[bool], b: &Node, e: &Node) -> bool {
    match (b.kind, e.kind) {
        _ if closes.is_empty() => true,
        (SegmentKind::V, SegmentKind::D) => closes[b.placement],
        (SegmentKind::D, SegmentKind::J) => closes[e.placement],
        _ => true,
    }
}

/// D segments up to and including a node, and the path's first kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Group {
    pub d_used: u8,
    pub first: SegmentKind,
}

impl Group {
    /// The paths that start at a node of `kind`.
    pub fn start(kind: SegmentKind) -> Self {
        Self {
            d_used: u8::from(kind == SegmentKind::D),
            first: kind,
        }
    }

    /// The group these paths join at a successor of `kind`.
    pub fn then(self, kind: SegmentKind) -> Self {
        Self {
            d_used: self.d_used + u8::from(kind == SegmentKind::D),
            ..self
        }
    }
}

/// The paths through one node with one [`Group`].
#[derive(Clone, Copy, Debug)]
pub(super) struct Slot {
    pub node: u32,
    pub group: Group,
    /// Start of the slot's `width` values in per-slot arrays; per-kind
    /// arrays start at `KINDS · values` ([`Slot::per_kind`]).
    pub values: usize,
    /// The slot's incoming edges in [`Lattice::edges`].
    edges: (u32, u32),
}

impl Slot {
    /// The slot's array for `kind` in a per-slot, per-kind buffer.
    pub fn per_kind(&self, kind: SegmentKind, width: usize) -> Range<usize> {
        let start = KINDS * self.values + kind.index() * width;
        start..start + width
    }
}

#[derive(Debug, Default)]
pub(super) struct Lattice {
    order: Vec<usize>,
    pub nodes: Vec<Node>,
    entries: Vec<u32>,
    pub slots: Vec<Slot>,
    /// Cells of all nodes: the length of per-cut arrays.
    cells: usize,
    /// Values of all slots opened so far: the length of per-slot arrays.
    values: usize,
    groups: Vec<Group>,
    /// Incoming edges of every slot opened so far, `(predecessor node,
    /// predecessor slot)`, each slot's in DAG order.
    edges: Vec<(u32, u32)>,
    /// `(group, predecessor node, predecessor slot)` while a node opens.
    pending: Vec<(u32, u32, u32)>,
    /// Per placement: whether a V or J can open or close a D's flank
    /// ([`super::units::FlankTable::closes`]); empty: every one can.
    closes: Vec<bool>,
    /// Whether some placement cannot.
    blocks: bool,
}

impl Lattice {
    /// Lay out the `placements` that `keep` (by index) as nodes in DAG order,
    /// with no slots yet. Nodes index the whole slice.
    pub fn layout(&mut self, placements: &[Placement], keep: impl Fn(usize) -> bool) {
        self.order.clear();
        self.order
            .extend((0..placements.len()).filter(|&i| keep(i)));
        self.order.sort_by_key(|&i| {
            let p = &placements[i];
            (p.kind.index(), p.q_start(), p.q_end(), i)
        });
        self.nodes.clear();
        self.entries.clear();
        self.cells = 0;
        for &i in &self.order {
            let p = &placements[i];
            debug_assert!(p.kind.is_rearranged(), "a C placement reached the lattice");
            let first = self.entries.len();
            self.entries
                .extend(p.clip.matches().map(|(q, ..)| q - p.q_start()));
            let node = Node {
                placement: i,
                kind: p.kind,
                q_start: p.q_start(),
                len: p.q_end() - p.q_start(),
                cuts: self.cells,
                entries: (first, self.entries.len()),
                slots: (0, 0),
            };
            self.cells += node.width();
            self.nodes.push(node);
        }
        self.slots.clear();
        self.edges.clear();
        self.values = 0;
    }

    pub fn cells(&self) -> usize {
        self.cells
    }

    pub fn values(&self) -> usize {
        self.values
    }

    /// Canonical entry offsets of `node`, ascending.
    pub fn entries(&self, node: &Node) -> &[u32] {
        &self.entries[node.entries.0..node.entries.1]
    }

    /// The nodes that may precede node `e` are `nodes[..predecessors(e)]`.
    pub fn predecessors(&self, e: usize) -> usize {
        let node = self.nodes[e];
        self.nodes.partition_point(|b| b.precedes(&node))
    }

    /// The nodes that node `b` may precede are `nodes[successors(b)..]`.
    pub fn successors(&self, b: usize) -> usize {
        let node = self.nodes[b];
        self.nodes.partition_point(|e| !node.precedes(e))
    }

    /// Open node `e`'s slots: `start` (paths starting here), then every group
    /// a live predecessor slot continues into, in DAG order, each with its
    /// incoming edges. Per-slot arrays must then hold [`Lattice::values`]
    /// values.
    pub fn open<G: Grammar>(
        &mut self,
        e: usize,
        grammar: &G,
        start: Option<Group>,
        live: impl Fn(usize) -> bool,
    ) -> Range<usize> {
        let node = self.nodes[e];
        self.groups.clear();
        self.groups.extend(start);
        self.pending.clear();
        for (b, before) in self.nodes[..self.predecessors(e)].iter().enumerate() {
            for t in before.slots() {
                let group = self.slots[t].group;
                if !(live(t)
                    && grammar.can_follow(before.kind, node.kind, group.d_used)
                    && (!self.blocks || joins(&self.closes, before, &node)))
                {
                    continue;
                }
                let next = group.then(node.kind);
                let g = match self.groups.iter().position(|&known| known == next) {
                    Some(g) => g,
                    None => {
                        self.groups.push(next);
                        self.groups.len() - 1
                    }
                };
                self.pending.push((g as u32, b as u32, t as u32));
            }
        }
        let first = self.slots.len();
        for (g, &group) in self.groups.iter().enumerate() {
            let from = self.edges.len() as u32;
            self.edges.extend(
                self.pending
                    .iter()
                    .filter(|&&(of, ..)| of == g as u32)
                    .map(|&(_, b, t)| (b, t)),
            );
            self.slots.push(Slot {
                node: e as u32,
                group,
                values: self.values,
                edges: (from, self.edges.len() as u32),
            });
            self.values += node.width();
        }
        self.nodes[e].slots = (first, self.slots.len());
        first..self.slots.len()
    }

    /// The live predecessor slots whose paths continue into slot `t`, with
    /// their node, in DAG order (recorded when `t` was opened).
    pub fn incoming(&self, t: usize) -> impl Iterator<Item = (Node, usize)> + '_ {
        let (from, to) = self.slots[t].edges;
        self.edges[from as usize..to as usize]
            .iter()
            .map(|&(b, u)| (self.nodes[b as usize], u as usize))
    }

    /// On an edge `b → e`: every canonical entry offset `o` of `e` inside
    /// `b`, with `z`, `b`'s exit offset at the same read position (the exits
    /// `b`'s paths can leave by for that entry are its offsets `≤ z`). A
    /// predecessor of another kind may start after `e`: entries before its
    /// start have no exit to follow.
    pub fn inside(&self, b: &Node, e: &Node) -> impl Iterator<Item = (usize, usize)> + '_ {
        let (start, end, at) = (b.q_start, b.q_end(), e.q_start);
        self.entries(e)
            .iter()
            .map(move |&o| (o, at + o))
            .take_while(move |&(_, y)| y < end)
            .filter_map(move |(o, y)| y.checked_sub(start).map(|z| (o as usize, z as usize)))
    }

    /// Which placements can open or close a D's flank (empty: all).
    pub fn set_closes(&mut self, closes: &[bool]) {
        self.closes.clear();
        self.closes.extend_from_slice(closes);
        self.blocks = closes.contains(&false);
    }

    /// Whether a D may follow `node` on a path: false only for a V that
    /// cannot open a D's flank ([`joins`]).
    #[inline]
    pub fn opens(&self, node: &Node) -> bool {
        !self.blocks || node.kind != SegmentKind::V || self.closes[node.placement]
    }

    /// The slot of `node` holding `group`, if opened.
    pub fn slot_of(&self, node: &Node, group: Group) -> Option<usize> {
        node.slots().find(|&u| self.slots[u].group == group)
    }
}
