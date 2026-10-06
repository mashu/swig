//! Placement DAG: bounded-topology semi-Markov Viterbi over clipped placements.
//!
//! Nodes are placements sorted by query start. A path visits placements in
//! query order. Each placement on the path is clipped independently: the
//! earlier of two consecutive placements keeps query bases before its **exit
//! cut** `x`, the later one keeps bases from its **entry cut** `y ≥ x` on,
//! and the bases in between are N/P. This is the full segmentation model:
//! when two placements overlap, AER-R's candidate A (D wins the overlap),
//! candidate B (V/J win it), every allocation in between, and leaving part of
//! the overlap to neither are all evaluated exactly.
//!
//! State: `(node, entry cut, D count, first kind)`. The score of a node's own
//! segment depends on its entry *and* exit cut, which is what makes this
//! semi-Markov. Naively each transition costs `#states × #exits × #entries`
//! clip evaluations.
//!
//! **Trick — separable scores → running maxima.** A clipped span scores
//! `exit(x).prefix − entry(e).prefix` (see [`crate::clip`]), the model adds a
//! cut-independent bonus, and the N/P score is linear in the inserted length.
//! So:
//!
//! 1. For a fixed exit `x`, the best predecessor state maximises
//!    `score_before − entry.prefix` among states whose entry leaves at least
//!    `min_retained` bases — a threshold on the entry position, answered by a
//!    binary search in the node's states sorted by entry with running maxima.
//! 2. The best exit for an entry `y` of the next placement is the maximum of
//!    those values over `x ≤ y` — a running maximum over the node's exits,
//!    built once per node and state group, then read in O(1) per entry.
//!
//! Entries are only tried inside the overlap (or at the next placement's start
//! when there is none): moving an entry further right clips a prefix of an
//! optimal local alignment, which never scores below zero, so it never helps.
//! Cost: `O(S log S + L log S + T·C)` per node, for `S` states, `L` exit cuts,
//! `T` successors and `C` entries — instead of `O(S·T·L·C)`.
//!
//! Ties: final paths compare by [`PathScore`] (score, then fewer D, then
//! germline support), so a with-D path must strictly beat the best no-D path.

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
    /// Unclipped placement scores accumulated before this node (tie-break).
    support: Score,
    back: Option<Back>,
}

impl State {
    /// `(score − entry prefix, support)`: lexicographic, and both parts are
    /// separable from the exit cut, so prefix maxima stay exact.
    #[inline]
    fn key(&self) -> (Score, Score) {
        (self.score - self.entry_prefix, self.support)
    }
}

/// Running maximum of `State::key` over states sorted by `q_lo`.
#[derive(Clone, Copy, Debug)]
struct Ranked {
    q_lo: u32,
    best_key: (Score, Score),
    state: u32,
}

#[derive(Clone, Copy, Debug)]
struct Group {
    d_used: u8,
    first: SegmentKind,
    start: usize,
    end: usize,
}

/// Best exit cut `≤ x` of a node for one state group: the predecessor state
/// and exit maximising `(score before − entry prefix + exit prefix − per_base ·
/// exit position, support)`.
#[derive(Clone, Copy, Debug)]
struct ExitBest {
    key: (Score, Score),
    exit: u32,
    state: u32,
}

/// Per-thread buffers for [`PlacementDag`].
#[derive(Debug, Default)]
pub struct DagScratch {
    order: Vec<u32>,
    states: Vec<Vec<State>>,
    perm: Vec<u32>,
    ranked: Vec<Ranked>,
    groups: Vec<Group>,
    /// Running best exits of the current node and group, `exits[kind][x − q_start]`
    /// per successor kind (the junction score may depend on it).
    exits: [Vec<Option<ExitBest>>; 3],
}

/// Fill `out[x − b.q_start]` with the best exit cut `≤ x` of `b` for states
/// `group`, folding in the linear junction score `per_base`.
fn exit_table(
    b: &Placement,
    group: &[Ranked],
    min_b: i64,
    per_base: Score,
    out: &mut Vec<Option<ExitBest>>,
) {
    out.clear();
    let mut running: Option<ExitBest> = None;
    for x in b.q_start()..=b.q_end() {
        if let Some(exit) = b.clip.exit(x) {
            if let Some(r) = best_up_to(group, i64::from(exit.q) - min_b) {
                let key = (
                    r.best_key.0 + exit.prefix - per_base * exit.q as Score,
                    r.best_key.1,
                );
                // `>=`: on ties keep the later exit, which retains more of `b`.
                if running.map_or(true, |best| key >= best.key) {
                    running = Some(ExitBest {
                        key,
                        exit: x,
                        state: r.state,
                    });
                }
            }
        }
        out.push(running);
    }
}

fn upsert(states: &mut Vec<State>, candidate: State) {
    match states.iter_mut().find(|s| {
        s.entry == candidate.entry && s.d_used == candidate.d_used && s.first == candidate.first
    }) {
        Some(existing)
            if (candidate.score, candidate.support) > (existing.score, existing.support) =>
        {
            *existing = candidate
        }
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
            exits,
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
                    support: 0,
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
                    let total = r.best_key.0
                        + exit.prefix
                        + self.model.segment_bonus(b, g.d_used)
                        + self.model.end(g.first, b.kind, g.d_used);
                    let candidate = PathScore {
                        score: total,
                        d_count: g.d_used,
                        support: r.best_key.1 + b.score(),
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

            // Transitions b → e, one state group and successor kind at a time.
            for g in groups.iter() {
                let group = &ranked[g.start..g.end];
                let bonus = self.model.segment_bonus(b, g.d_used);
                let mut built = [false; 3];
                for epos in pos + 1..n {
                    let e = &placements[order[epos] as usize];
                    if e.q_start() <= b.q_start()
                        || !self.grammar.can_follow(b.kind, e.kind, g.d_used)
                    {
                        continue;
                    }
                    let per_base = self.model.junction_per_base(b.kind, e.kind);
                    let table = &mut exits[e.kind.index()];
                    if !built[e.kind.index()] {
                        exit_table(b, group, min_b, per_base, table);
                        built[e.kind.index()] = true;
                    }
                    let min_e = self.min_bases(e.kind);
                    let d_next = g.d_used + u8::from(e.kind == SegmentKind::D);
                    let last_entry = if e.q_start() < b.q_end() {
                        b.q_end().min(e.q_end().saturating_sub(1))
                    } else {
                        e.q_start()
                    };
                    for y in e.q_start()..=last_entry {
                        let Some(en) = e.clip.entry(y) else {
                            continue;
                        };
                        if e.q_end() < en.q + min_e {
                            continue;
                        }
                        let Some(best_exit) = table[(y.min(b.q_end()) - b.q_start()) as usize]
                        else {
                            continue;
                        };
                        upsert(
                            &mut states[epos],
                            State {
                                entry: y,
                                q_lo: en.q,
                                entry_prefix: en.prefix,
                                d_used: d_next,
                                first: g.first,
                                score: best_exit.key.0 + bonus + per_base * en.q as Score,
                                support: best_exit.key.1 + b.score(),
                                back: Some(Back {
                                    node: pos as u32,
                                    state: best_exit.state,
                                    exit: best_exit.exit,
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
        // Single cut: 40 + 2 + pair bonus 24 = 66. Two cuts (V to 11, J from 15): 68.
        assert_eq!(path.score.score, 68);
        let (v_seg, j_seg) = (path.segments[0], path.segments[1]);
        assert_eq!((v_seg.span.q_hi, j_seg.span.q_lo), (11, 15));
    }
}
