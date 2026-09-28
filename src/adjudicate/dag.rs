//! Placement DAG: bounded-topology semi-Markov Viterbi over clipped placements.
//!
//! Nodes are placements sorted by query start. A path visits placements in
//! query order; consecutive placements either do not overlap (the gap between
//! them is an N/P region) or overlap, in which case a **cut** `c` splits the
//! overlap: the earlier placement keeps query bases `< c`, the later one keeps
//! bases `≥ c`. Every cut in the overlap is considered, so AER-R's candidate A
//! (D wins the overlap), candidate B (V/J win it) and every allocation in
//! between are all evaluated, exactly, and the "one-base V restore" patch is
//! unnecessary.
//!
//! State: `(node, entry cut, D count, first kind)`. The score of a node's own
//! segment depends on its entry *and* exit cut, which is what makes this
//! semi-Markov. Naively each transition costs `#states × #cuts` clip
//! evaluations.
//!
//! **Trick — separable segment scores → prefix maxima.** A clipped span scores
//! `exit(c).prefix − entry(e).prefix` (see [`crate::clip`]), and the model adds
//! a cut-independent bonus. So for a fixed exit cut the best predecessor state
//! is the one maximising `score_before − entry.prefix` among states whose
//! entry leaves at least `min_retained` bases — a *threshold* on the entry
//! position. Sorting a node's states by entry position once and keeping a
//! running maximum turns every (transition, cut) into one binary search:
//! `O(S log S + T·C·log S)` instead of `O(S·T·C)`.
//!
//! Ties: final paths compare by [`PathScore`] (score, then fewer D), so a
//! with-D path must strictly beat the best no-D path.

use super::model::{Grammar, PathModel, VdjGrammar};
use super::{Adjudicator, CallPath, CalledSegment, EngineeringObjective, PathScore};
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlacementDag<G = VdjGrammar, M = EngineeringObjective> {
    pub grammar: G,
    pub model: M,
    /// Minimum query bases a clipped segment must keep, indexed by [`SegmentKind::index`].
    pub min_retained: [u32; 3],
}

impl Default for PlacementDag {
    fn default() -> Self {
        Self {
            grammar: VdjGrammar::VDJ,
            model: EngineeringObjective::default(),
            min_retained: [10, 3, 6],
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Back {
    node: u32,
    state: u32,
    exit: u32,
}

#[derive(Clone, Copy, Debug)]
struct State {
    /// Query cut this node's placement is clipped to on the left.
    entry: u32,
    /// First retained query base for that cut, and the prefix score there.
    q_lo: u32,
    entry_prefix: Score,
    d_used: u8,
    first: SegmentKind,
    /// Path score accumulated *before* this node's own segment.
    score: Score,
    back: Option<Back>,
}

impl State {
    #[inline]
    fn key(&self) -> Score {
        self.score - self.entry_prefix
    }
}

/// Running maximum of `State::key` over states sorted by `q_lo`.
#[derive(Clone, Copy, Debug)]
struct Ranked {
    q_lo: u32,
    best_key: Score,
    state: u32,
}

#[derive(Clone, Copy, Debug)]
struct Group {
    d_used: u8,
    first: SegmentKind,
    start: usize,
    end: usize,
}

/// Per-thread buffers for [`PlacementDag`].
#[derive(Debug, Default)]
pub struct DagScratch {
    order: Vec<u32>,
    states: Vec<Vec<State>>,
    perm: Vec<u32>,
    ranked: Vec<Ranked>,
    groups: Vec<Group>,
}

fn upsert(states: &mut Vec<State>, candidate: State) {
    match states.iter_mut().find(|s| {
        s.entry == candidate.entry && s.d_used == candidate.d_used && s.first == candidate.first
    }) {
        Some(existing) if candidate.score > existing.score => *existing = candidate,
        Some(_) => {}
        None => states.push(candidate),
    }
}

/// Group a node's states by `(d_used, first)` and build prefix maxima by entry position.
fn rank_states(
    states: &[State],
    perm: &mut Vec<u32>,
    ranked: &mut Vec<Ranked>,
    groups: &mut Vec<Group>,
) {
    perm.clear();
    perm.extend(0..states.len() as u32);
    perm.sort_by_key(|&i| {
        let s = &states[i as usize];
        (s.d_used, s.first, s.q_lo)
    });
    ranked.clear();
    groups.clear();
    for &i in perm.iter() {
        let s = &states[i as usize];
        let new_group = groups
            .last()
            .map_or(true, |g| g.d_used != s.d_used || g.first != s.first);
        if new_group {
            groups.push(Group {
                d_used: s.d_used,
                first: s.first,
                start: ranked.len(),
                end: ranked.len(),
            });
        }
        let running = match ranked.last() {
            Some(prev) if !new_group && prev.best_key >= s.key() => *prev,
            _ => Ranked {
                q_lo: s.q_lo,
                best_key: s.key(),
                state: i,
            },
        };
        ranked.push(Ranked {
            q_lo: s.q_lo,
            ..running
        });
        groups.last_mut().expect("pushed above").end = ranked.len();
    }
}

/// Best ranked state whose first retained base is `≤ threshold`.
#[inline]
fn best_up_to(ranked: &[Ranked], threshold: i64) -> Option<Ranked> {
    let k = ranked.partition_point(|r| i64::from(r.q_lo) <= threshold);
    k.checked_sub(1).map(|i| ranked[i])
}

impl<G: Grammar, M: PathModel> PlacementDag<G, M> {
    pub fn new(grammar: G, model: M, min_retained: [u32; 3]) -> Self {
        Self {
            grammar,
            model,
            min_retained,
        }
    }

    #[inline]
    fn min_bases(&self, kind: SegmentKind) -> u32 {
        self.min_retained[kind.index()].max(1)
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
        let n = placements.len();
        let DagScratch {
            order,
            states,
            perm,
            ranked,
            groups,
        } = scratch;
        order.clear();
        order.extend(0..n as u32);
        order.sort_by_key(|&i| {
            let p = &placements[i as usize];
            (p.q_start(), p.q_end(), i)
        });
        if states.len() < n {
            states.resize_with(n, Vec::new);
        }
        for node in states.iter_mut().take(n) {
            node.clear();
        }

        for (pos, &idx) in order.iter().enumerate() {
            let p = &placements[idx as usize];
            if !self.grammar.can_start(p.kind) {
                continue;
            }
            if let Some(entry) = p.clip.entry(p.q_start()) {
                states[pos].push(State {
                    entry: p.q_start(),
                    q_lo: entry.q,
                    entry_prefix: entry.prefix,
                    d_used: u8::from(p.kind == SegmentKind::D),
                    first: p.kind,
                    score: 0,
                    back: None,
                });
            }
        }

        let mut best: Option<(PathScore, Back)> = None;
        for pos in 0..n {
            if states[pos].is_empty() {
                continue;
            }
            let b = &placements[order[pos] as usize];
            let min_b = i64::from(self.min_bases(b.kind));
            rank_states(&states[pos], perm, ranked, groups);

            // Paths ending at b (b unclipped on the right).
            if let Some(exit) = b.clip.exit(b.q_end()) {
                for g in groups
                    .iter()
                    .filter(|g| self.grammar.can_end(g.first, b.kind))
                {
                    let Some(r) = best_up_to(&ranked[g.start..g.end], i64::from(exit.q) - min_b)
                    else {
                        continue;
                    };
                    let total = r.best_key
                        + exit.prefix
                        + self.model.segment_bonus(b, g.d_used)
                        + self.model.end(g.first, b.kind, g.d_used);
                    let candidate = PathScore {
                        score: total,
                        d_count: g.d_used,
                    };
                    if best.map_or(true, |(current, _)| candidate > current) {
                        best = Some((
                            candidate,
                            Back {
                                node: pos as u32,
                                state: r.state,
                                exit: b.q_end(),
                            },
                        ));
                    }
                }
            }

            // Transitions b → e.
            for epos in pos + 1..n {
                let e = &placements[order[epos] as usize];
                if e.q_start() <= b.q_start() {
                    continue;
                }
                let min_e = self.min_bases(e.kind);
                let overlapping = e.q_start() < b.q_end();
                let (first_cut, last_cut) = if overlapping {
                    (e.q_start(), b.q_end().min(e.q_end().saturating_sub(1)))
                } else {
                    (0, 0)
                };
                for g in groups.iter() {
                    if !self.grammar.can_follow(b.kind, e.kind, g.d_used) {
                        continue;
                    }
                    let d_next = g.d_used + u8::from(e.kind == SegmentKind::D);
                    let bonus = self.model.segment_bonus(b, g.d_used);
                    let group = &ranked[g.start..g.end];
                    let cuts = if overlapping {
                        first_cut..=last_cut
                    } else {
                        0..=0
                    };
                    for c in cuts {
                        let (exit_at, entry_at) = if overlapping {
                            (c, c)
                        } else {
                            (b.q_end(), e.q_start())
                        };
                        let (Some(xb), Some(en)) = (b.clip.exit(exit_at), e.clip.entry(entry_at))
                        else {
                            continue;
                        };
                        if e.q_end() < en.q + min_e {
                            continue;
                        }
                        let Some(r) = best_up_to(group, i64::from(xb.q) - min_b) else {
                            continue;
                        };
                        let score = r.best_key
                            + xb.prefix
                            + bonus
                            + self
                                .model
                                .junction(b.kind, e.kind, en.q.saturating_sub(xb.q));
                        upsert(
                            &mut states[epos],
                            State {
                                entry: entry_at,
                                q_lo: en.q,
                                entry_prefix: en.prefix,
                                d_used: d_next,
                                first: g.first,
                                score,
                                back: Some(Back {
                                    node: pos as u32,
                                    state: r.state,
                                    exit: exit_at,
                                }),
                            },
                        );
                    }
                }
            }
        }

        let (score, mut at) = best?;
        let mut segments = Vec::with_capacity(4);
        loop {
            let state = states[at.node as usize][at.state as usize];
            let placement = order[at.node as usize] as usize;
            let span = placements[placement]
                .clip
                .span(state.entry, at.exit)
                .expect("DAG only records spans that keep min_retained bases");
            segments.push(CalledSegment {
                placement,
                entry: state.entry,
                exit: at.exit,
                span,
            });
            match state.back {
                Some(back) => at = back,
                None => break,
            }
        }
        segments.reverse();
        Some(CallPath { segments, score })
    }
}
