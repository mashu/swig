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

use std::ops::Range;

use super::model::Grammar;
use crate::germline::SegmentKind;
use crate::placement::Placement;

/// Segment kinds: passes keep one array per successor or predecessor kind.
pub(super) const KINDS: usize = SegmentKind::ALL.len();

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
    /// a live predecessor slot continues into, in DAG order. Per-slot arrays
    /// must then hold [`Lattice::values`] values.
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
        for b in &self.nodes[..self.predecessors(e)] {
            for t in b.slots() {
                let group = self.slots[t].group;
                let next = group.then(node.kind);
                if live(t)
                    && grammar.can_follow(b.kind, node.kind, group.d_used)
                    && !self.groups.contains(&next)
                {
                    self.groups.push(next);
                }
            }
        }
        let first = self.slots.len();
        for &group in &self.groups {
            self.slots.push(Slot {
                node: e as u32,
                group,
                values: self.values,
            });
            self.values += node.width();
        }
        self.nodes[e].slots = (first, self.slots.len());
        first..self.slots.len()
    }

    /// Live slots of node `e`'s predecessors whose paths continue into
    /// `group` at `e`, with their node.
    pub fn incoming<'a, G: Grammar, L: Fn(usize) -> bool>(
        &'a self,
        e: usize,
        grammar: &'a G,
        group: Group,
        live: &'a L,
    ) -> impl Iterator<Item = (Node, usize)> + 'a {
        let kind = self.nodes[e].kind;
        self.nodes[..self.predecessors(e)]
            .iter()
            .flat_map(move |b| {
                b.slots()
                    .filter(move |&u| {
                        let before = self.slots[u].group;
                        live(u)
                            && before.then(kind) == group
                            && grammar.can_follow(b.kind, kind, before.d_used)
                    })
                    .map(move |u| (*b, u))
            })
    }

    /// The slot of `node` holding `group`, if opened.
    pub fn slot_of(&self, node: &Node, group: Group) -> Option<usize> {
        node.slots().find(|&u| self.slots[u].group == group)
    }
}
