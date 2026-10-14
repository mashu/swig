//! Path grammar (which topologies are legal) and path model (how they score).
//!
//! Separating the two keeps "which topologies?" (a grammar question) apart
//! from "what does a segment cost?" (a scoring question). The DAG needs one
//! structural property from a [`PathModel`]: a segment contributes
//! *clip score − a presence cost that does not depend on the cut
//! positions*, and a junction contributes a score linear in its inserted
//! length. That
//! separability is what lets running maxima (or sums) evaluate every pair of
//! cuts in O(1) amortised per cut (see [`super::PlacementDag`]).

use super::units::ResolvedOpenFlank;
use crate::germline::SegmentKind;
use crate::scoring::Score;

/// Which segment sequences form a path. Transitions between kinds must go
/// forward in [`SegmentKind::REARRANGED`] order (V, D, J): the DAG visits kinds in
/// that order ([`super::PlacementDag`]).
pub trait Grammar: Send + Sync {
    fn can_start(&self, kind: SegmentKind) -> bool;
    /// `d_used` counts D segments up to and including `prev`.
    fn can_follow(&self, prev: SegmentKind, next: SegmentKind, d_used: u8) -> bool;
    /// Whether a path from `first` to `last` with `d_count` D segments is complete.
    fn can_end(&self, first: SegmentKind, last: SegmentKind, d_count: u8) -> bool;
}

/// V → D{0..max_d} → J with optional truncation at either end. The engine
/// builds [`VdjGrammar::VDJ`], or [`VdjGrammar::VJ`] without a D database:
/// an AIRR row holds one D, and Swig's double-D screen reports tandem D
/// beside the call ([`crate::double_d`]). The DAG itself handles any
/// `max_d`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VdjGrammar {
    pub max_d: u8,
    /// Paths must start with V.
    pub require_v: bool,
    /// Paths must end with J.
    pub require_j: bool,
}

impl VdjGrammar {
    /// Heavy-chain default: VJ or VDJ, partial reads allowed.
    pub const VDJ: Self = Self {
        max_d: 1,
        require_v: false,
        require_j: false,
    };
    /// Light chains: no D.
    pub const VJ: Self = Self {
        max_d: 0,
        require_v: false,
        require_j: false,
    };
}

impl Default for VdjGrammar {
    fn default() -> Self {
        Self::VDJ
    }
}

impl Grammar for VdjGrammar {
    fn can_start(&self, kind: SegmentKind) -> bool {
        match kind {
            SegmentKind::V => true,
            SegmentKind::D => !self.require_v && self.max_d > 0,
            SegmentKind::J => !self.require_v,
            SegmentKind::C => false,
        }
    }

    fn can_follow(&self, prev: SegmentKind, next: SegmentKind, d_used: u8) -> bool {
        use SegmentKind::{D, J, V};
        match (prev, next) {
            (V, D) | (D, D) => d_used < self.max_d,
            (V, J) | (D, J) => true,
            _ => false,
        }
    }

    fn can_end(&self, first: SegmentKind, last: SegmentKind, _d_count: u8) -> bool {
        let has_flank = first != SegmentKind::D || last != SegmentKind::D;
        has_flank && (!self.require_j || last == SegmentKind::J)
    }
}

pub trait PathModel: Send + Sync {
    /// Prior cost of including a segment of `kind` (`P*_k`), in score units
    /// of that kind's scoring and unrounded. Every segment of the kind on a
    /// path pays it. On the fixed scale an integer path score beats it
    /// exactly when it beats `⌊P*_k⌋` ([`PathModel::presence_score`]).
    fn presence(&self, kind: SegmentKind) -> f64;

    /// `⌊P*_k⌋`: for integer scores `S`, `S > P*_k ⟺ S > ⌊P*_k⌋`, so this is
    /// the fixed-scale cost of a segment of `kind`, not an approximation of
    /// it.
    fn presence_score(&self, kind: SegmentKind) -> Score {
        self.presence(kind).floor() as Score
    }

    /// Smallest D local score that can improve a fixed-scale path; lower
    /// ones can be pruned before adjudication. `Score::MIN` (no floor)
    /// unless the model proves one.
    fn d_floor(&self) -> Score {
        Score::MIN
    }

    /// Score per unassigned query base (N/P region) between two segments.
    /// Junction scores are linear in the inserted length, which keeps both
    /// cuts of a transition separable (see [`super::PlacementDag`]).
    fn junction_per_base(&self, _prev: SegmentKind, _next: SegmentKind) -> Score {
        0
    }

    /// Edge score for `inserted` unassigned query bases between two segments.
    fn junction(&self, prev: SegmentKind, next: SegmentKind, inserted: u32) -> Score {
        self.junction_per_base(prev, next) * inserted as Score
    }

    /// Whole-path term (e.g. the V+J pair bonus).
    fn end(&self, _first: SegmentKind, _last: SegmentKind, _d_count: u8) -> Score {
        0
    }

    /// The cost at the open end of a path whose first or last segment is a
    /// D ([`super::units::FlankTable`]), read per read at the path's two
    /// ends beside [`PathModel::end`]; `None`: no such cost.
    fn open_flank(&self) -> Option<&ResolvedOpenFlank> {
        None
    }
}

impl<M: PathModel> PathModel for &M {
    fn presence(&self, kind: SegmentKind) -> f64 {
        (**self).presence(kind)
    }

    fn d_floor(&self) -> Score {
        (**self).d_floor()
    }

    fn junction_per_base(&self, prev: SegmentKind, next: SegmentKind) -> Score {
        (**self).junction_per_base(prev, next)
    }

    fn end(&self, first: SegmentKind, last: SegmentKind, d_count: u8) -> Score {
        (**self).end(first, last, d_count)
    }

    fn open_flank(&self) -> Option<&ResolvedOpenFlank> {
        (**self).open_flank()
    }
}

/// The path objective, `Σ S_clip − P_V·#V − P_D·#D − P_J·#J +
/// pair_bonus·[V…J] + per-base N/P cost`, read as a log-likelihood ratio.
///
/// Substitution scores are log-odds (in units of `1/λ` nats) of "templated by
/// this germline" against "random sequence"; N/P bases score 0 under both.
/// The remaining terms are log prior odds of the path's hypotheses:
///
/// * `presence` (`P*_k` per kind) — a D is one of `N` genes at one of `L`
///   positions of the junction, so its prior is spread over `N·L`
///   hypotheses: including the best one is worth it when
///   `λ·S_D > ln(N·L) − ln(prior odds of D)`. A J is one of the J genes at
///   one of `L` positions after the V, and a V one of the V genes at one of
///   `L` positions before the J, so the same rule prices them: without that
///   cost any J (V) clip of positive value would complete (start) a path
///   for free, and a read that stops before its J (starts after its V)
///   would get one from a few chance bases of an outer flank or of a
///   random end. [`crate::config::DPrior`] derives all three from the
///   databases at build time. They are kept unrounded: the posterior uses
///   them as they are, and the fixed-scale call as `⌊P*_k⌋`
///   ([`PathModel::presence_score`]).
/// * `pair_bonus` — prior odds of a complete V…J rearrangement over a partial
///   read. 0 by default: on full-length reads the J's own clip score already
///   completes the path, and 0–48 changed no call.
/// * `open_flank` — a path that ends on a D asserts a junction whose J was
///   missed in the read bases after it (and one that starts on a D, a V
///   missed before it). It pays the log-likelihood ratio of that miss
///   against "no junction there" ([`super::units::FlankTable`]), which
///   removes chance Ds far from the anchor of reads without a D.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EngineeringObjective {
    pub pair_bonus: Score,
    /// `P*_k` per rearranged kind ([`SegmentKind::REARRANGED`]), indexed by
    /// [`SegmentKind::index`].
    pub presence: [f64; 3],
    /// Score per unassigned junction base (0 = no N/P length prior).
    pub insertion_per_base: Score,
    /// The cost at the open end of a D-first or D-last path (`None` = off).
    pub open_flank: Option<ResolvedOpenFlank>,
}

#[cfg(test)]
impl EngineeringObjective {
    /// The objective of the unit tests: presence costs `presence`, no pair
    /// bonus, no N/P cost.
    pub(crate) fn priced(presence: [f64; 3]) -> Self {
        Self {
            pair_bonus: 0,
            presence,
            insertion_per_base: 0,
            open_flank: None,
        }
    }
}

impl PathModel for EngineeringObjective {
    fn presence(&self, kind: SegmentKind) -> f64 {
        self.presence[kind.index()]
    }

    /// Smallest D local score that can improve a fixed-scale path: `P_D + 1`,
    /// with `P_D = ⌊P*_D⌋`. It bounds nothing for a per-read call: there a D
    /// scoring 10 ≤ `P_D` (8 matches, 2 mismatches) is worth 6.40 nats at
    /// `p̂ = 0.3`, above its cost `c_D` of 5.38 for `P*_D = 10.45`
    /// (`ReadModel`).
    ///
    /// A clipped span of an optimal local alignment never outscores the
    /// whole alignment (its prefixes and suffixes score ≥ 0), so a D scoring
    /// at most `P_D` adds at most 0. Dropping it leaves its neighbours' cuts
    /// valid (the exit before it is ≤ the entry after it), so the no-D path
    /// scores at least as much and "with-D must strictly beat no-D" rejects
    /// it. With a negative per-base junction score a D could pay for itself
    /// by replacing N bases, so no floor is derived then.
    fn d_floor(&self) -> Score {
        if self.insertion_per_base < 0 {
            Score::MIN
        } else {
            self.presence_score(SegmentKind::D) + 1
        }
    }

    fn junction_per_base(&self, _prev: SegmentKind, _next: SegmentKind) -> Score {
        self.insertion_per_base
    }

    fn end(&self, first: SegmentKind, last: SegmentKind, _d_count: u8) -> Score {
        if first == SegmentKind::V && last == SegmentKind::J {
            self.pair_bonus
        } else {
            0
        }
    }

    fn open_flank(&self) -> Option<&ResolvedOpenFlank> {
        self.open_flank.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adjudicate::units::{CallModel, FlankTable, ReadContext, ReadModel};
    use crate::adjudicate::{CalledSegment, DagScratch, Mode, PlacementDag};
    use crate::alphabet::encode;
    use crate::config::{DefaultAssigner, EngineConfig};
    use crate::pipeline::VdjAssigner;
    use crate::placement::Placement;
    use crate::sim::{Rng, SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

    /// The open-flank cost is non-negative and charges D only, so the D
    /// floor `P_D + 1` still bounds a fixed-scale call: on mid and
    /// hypermutated reads cut 0–80 bases around a D, with the cost on,
    /// pruning every D below the floor leaves the call (path and key)
    /// unchanged.
    #[test]
    fn d_floor_bounds_pruning_with_the_open_flank_cost() {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
        let config = EngineConfig::default();
        let assigner =
            DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
                .unwrap();
        let objective =
            config
                .objective
                .resolve(&db.v, &config.v, Some(&db.d), &config.d, &db.j, &config.j);
        assert!(objective.open_flank.is_some());
        let dag = PlacementDag::new(VdjGrammar::VDJ, objective, config.min_retained, Mode::Fixed);
        let floor = dag.model.d_floor();
        let read_model = ReadModel::fixed(1.0, objective.presence);
        let d_presence = objective.presence[SegmentKind::D.index()];
        let (mut ws, mut scratch, mut table) = (
            assigner.workspace(),
            DagScratch::default(),
            FlankTable::default(),
        );
        let mut rng = Rng::new(31);
        let (mut checked, mut pruned, mut charged) = (0, 0, 0);
        for preset in [SimParams::mid(), SimParams::hard()] {
            let sim = Simulator::new(&db.v, Some(&db.d), &db.j, preset).unwrap();
            for read in sim.reads(40, 17).unwrap() {
                let call = assigner.assign(&encode(&read.sequence), &mut ws);
                let Some(d) = call.placements.iter().find(|p| p.kind == SegmentKind::D) else {
                    continue;
                };
                let lo = d.q_start().saturating_sub(rng.below(81) as u32) as usize;
                let hi = (d.q_end() + rng.below(81) as u32).min(call.query.len() as u32) as usize;
                let cut = assigner.assign(&call.query[lo..hi], &mut ws);
                for dropped in [None, Some(SegmentKind::J), Some(SegmentKind::V)] {
                    let all: Vec<Placement> = cut
                        .placements
                        .iter()
                        .filter(|p| Some(p.kind) != dropped)
                        .cloned()
                        .collect();
                    let kept: Vec<Placement> = all
                        .iter()
                        .filter(|p| p.kind != SegmentKind::D || p.score() >= floor)
                        .cloned()
                        .collect();
                    pruned += all.len() - kept.len();
                    let context = ReadContext::of(&all, cut.query.len() as u32);
                    let flank = objective.open_flank.as_ref();
                    table.fill(flank, d_presence, CallModel::FixedScale, &context, &all);
                    let full = dag.viterbi(&all, &read_model, &table, &mut scratch);
                    table.fill(flank, d_presence, CallModel::FixedScale, &context, &kept);
                    let fewer = dag.viterbi(&kept, &read_model, &table, &mut scratch);
                    // The calls by what they hold, as the placement lists differ.
                    let held = |placements: &[Placement], path: &Option<(_, Vec<_>)>| {
                        path.as_ref()
                            .map(|(key, segments): &(_, Vec<CalledSegment>)| {
                                let spans: Vec<_> = segments
                                    .iter()
                                    .map(|s| {
                                        let p = &placements[s.placement];
                                        (p.kind, p.alleles.clone(), p.q_start(), s.entry, s.exit)
                                    })
                                    .collect();
                                (*key, spans)
                            })
                    };
                    if let Some((_, path)) = &full {
                        charged += usize::from(path.iter().any(|s| {
                            let p = &all[s.placement];
                            table.units(super::super::Flank::End, p.kind, p.q_start(), p.q_end())
                                > 0
                        }));
                    }
                    assert_eq!(held(&all, &full), held(&kept, &fewer));
                    checked += 1;
                }
            }
        }
        assert!(
            checked >= 150 && pruned >= 100 && charged >= 10,
            "{checked} reads, {pruned} D pruned, {charged}"
        );
    }
}
