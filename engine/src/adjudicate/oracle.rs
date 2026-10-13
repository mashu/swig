//! Exhaustive oracles for the placement DAG on the pipeline's own
//! placements: the Viterbi key under the fixed and the per-read model, and
//! the conditioned pass behind max-marginal D decoding, against a memoised
//! search over every grammatical path and every pair of canonical cuts.

use std::collections::{HashMap, HashSet};

use super::decode::{Topology, Within};
use super::units::ReadModel;
use super::{
    Adjudicator, CalledSegment, DagScratch, EngineeringObjective, EvidenceModel, Grammar, Mode,
    PathModel, PathScore, PlacementDag, Units, VdjGrammar, divergence,
};
use crate::alphabet::encode;
use crate::config::{DefaultAssigner, EngineConfig};
use crate::germline::SegmentKind;
use crate::pipeline::VdjAssigner;
use crate::placement::Placement;
use crate::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

/// VJ, VDJ and VDDJ paths: the engine builds at most one D, the DAG
/// handles any number.
const TANDEM_D: VdjGrammar = VdjGrammar {
    max_d: 2,
    ..VdjGrammar::VDJ
};

/// Exhaustive search of the segmentation model: every path of `grammar`, and
/// for every transition every canonical exit cut `x` and every canonical
/// entry cut `y ≥ x` of the next placement, with no shortcut, valued term by
/// term under `read`. Memoised on `(placement, entry, D count, first)`.
struct Brute<'a, G = VdjGrammar> {
    dag: &'a PlacementDag,
    grammar: &'a G,
    placements: &'a [Placement],
    read: ReadModel,
    memo: HashMap<(usize, u32, u8, SegmentKind), Option<PathScore<Units>>>,
}

/// Canonical entry and exit cuts of a placement (one per Match column).
fn cuts(p: &Placement) -> (Vec<u32>, Vec<u32>) {
    p.clip.matches().map(|(q, ..)| (q, q + 1)).unzip()
}

impl<'a> Brute<'a> {
    fn new(dag: &'a PlacementDag, placements: &'a [Placement], read: ReadModel) -> Self {
        Self::within(dag, &dag.grammar, placements, read)
    }
}

impl<'a, G: Grammar> Brute<'a, G> {
    /// The search under another grammar than `dag`'s.
    fn within(
        dag: &'a PlacementDag,
        grammar: &'a G,
        placements: &'a [Placement],
        read: ReadModel,
    ) -> Self {
        Self {
            dag,
            grammar,
            placements,
            read,
            memo: HashMap::new(),
        }
    }

    fn min(&self, kind: SegmentKind) -> u32 {
        self.dag.min_retained[kind.index()].max(1)
    }

    /// `p`'s alignment on `[entry, exit)`: its substitution columns and the
    /// gap runs between them.
    fn clip(&self, p: &Placement, entry: u32, exit: u32) -> Units {
        let mut value = 0;
        let mut last: Option<i32> = None;
        for (_, before, after) in p.clip.matches().filter(|(q, ..)| (entry..exit).contains(q)) {
            if let Some(previous) = last {
                value += self.read.prior(before - previous);
            }
            value += self.read.column(p.kind, after - before);
            last = Some(after);
        }
        value
    }

    /// The cut-independent term of `p`, its kind's presence cost.
    fn bonus(&self, p: &Placement) -> Units {
        -self.read.presence(p.kind)
    }

    /// Best rest of a path whose current segment is `b`, entered at `entry`.
    fn best_from(
        &mut self,
        b: usize,
        entry: u32,
        d_used: u8,
        first: SegmentKind,
    ) -> Option<PathScore<Units>> {
        let key = (b, entry, d_used, first);
        if let Some(known) = self.memo.get(&key) {
            return *known;
        }
        let (grammar, model) = (self.grammar, &self.dag.model);
        let pb = &self.placements[b];
        let bonus = self.bonus(pb);
        let mut best: Option<PathScore<Units>> = None;
        let mut consider = |candidate: PathScore<Units>| {
            if best.is_none_or(|current| candidate > current) {
                best = Some(candidate);
            }
        };
        if grammar.can_end(first, pb.kind, d_used)
            && pb
                .clip
                .span(entry, pb.q_end())
                .is_some_and(|s| s.q_len() >= self.min(pb.kind))
        {
            consider(PathScore {
                score: self.clip(pb, entry, pb.q_end())
                    + bonus
                    + self.read.prior(model.end(first, pb.kind, d_used)),
                d_count: d_used,
                support: pb.score(),
            });
        }
        for e in 0..self.placements.len() {
            let pe = &self.placements[e];
            // Kinds in grammar order wherever they start; one kind (tandem D)
            // in order of start.
            let same_kind_before = pe.kind == pb.kind && pe.q_start() <= pb.q_start();
            if same_kind_before || !grammar.can_follow(pb.kind, pe.kind, d_used) {
                continue;
            }
            let d_next = d_used + u8::from(pe.kind == SegmentKind::D);
            let (_, exits) = cuts(pb);
            let (entries, _) = cuts(pe);
            for x in exits.into_iter().filter(|&x| x > entry) {
                let Some(span) = pb
                    .clip
                    .span(entry, x)
                    .filter(|s| s.q_len() >= self.min(pb.kind))
                else {
                    continue;
                };
                let here = self.clip(pb, entry, x) + bonus;
                for &y in entries.iter().filter(|&&y| y >= x) {
                    let Some(next) = pe.clip.entry(y) else {
                        continue;
                    };
                    let inserted = next.q.saturating_sub(span.q_hi);
                    let junction = self.read.prior(model.junction(pb.kind, pe.kind, inserted));
                    let Some(rest) = self.best_from(e, y, d_next, first) else {
                        continue;
                    };
                    consider(PathScore {
                        score: here + junction + rest.score,
                        d_count: rest.d_count,
                        support: pb.score() + rest.support,
                    });
                }
            }
        }
        self.memo.insert(key, best);
        best
    }

    /// The best path over every start.
    fn best(&mut self) -> Option<PathScore<Units>> {
        (0..self.placements.len())
            .filter(|&s| self.grammar.can_start(self.placements[s].kind))
            .filter_map(|s| {
                let p = &self.placements[s];
                self.best_from(s, p.q_start(), u8::from(p.kind == SegmentKind::D), p.kind)
            })
            .max()
    }

    /// A called path's value, term by term.
    fn rescore(&self, segments: &[CalledSegment]) -> Units {
        let kind = |s: &CalledSegment| self.placements[s.placement].kind;
        let mut d_used = 0;
        let mut value = 0;
        for (k, segment) in segments.iter().enumerate() {
            let p = &self.placements[segment.placement];
            d_used += u8::from(p.kind == SegmentKind::D);
            value += self.clip(p, segment.entry, segment.exit) + self.bonus(p);
            if let Some(before) = k.checked_sub(1).map(|k| &segments[k]) {
                let inserted = segment.entry - before.exit;
                value += self
                    .read
                    .prior(self.dag.model.junction(kind(before), p.kind, inserted));
            }
        }
        let (first, last) = (kind(&segments[0]), kind(&segments[segments.len() - 1]));
        value + self.read.prior(self.dag.model.end(first, last, d_used))
    }
}

/// The fixed scale of `dag`'s model: scores are units.
fn fixed(dag: &PlacementDag) -> ReadModel {
    ReadModel::fixed(1.0, dag.model.presence)
}

/// The DAG's best key and the exhaustive one on `placements` under `read`;
/// the DAG's traced path re-scores to its key.
fn dag_and_brute(
    dag: &PlacementDag,
    placements: &[Placement],
    read: ReadModel,
    scratch: &mut DagScratch,
) -> (Option<PathScore<Units>>, Option<PathScore<Units>>) {
    let fast = dag.viterbi(placements, &read, scratch);
    let mut brute = Brute::new(dag, placements, read);
    if let Some((key, segments)) = &fast {
        assert_eq!(brute.rescore(segments), key.score, "traced path re-scored");
    }
    (fast.map(|(key, _)| key), brute.best())
}

/// IgSim reads mutate V, D and J, so an alignment can carry mismatches right
/// after its neighbour's end, where entering later and leaving them to N/P
/// wins (an entry beyond the overlap). The pipeline's own placements, with and
/// without a per-base N/P cost.
#[test]
fn dag_equals_brute_force_on_hypermutated_igsim_reads() {
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let base = config
        .objective
        .resolve(Some(&db.d), &config.d, &db.j, &config.j.scoring);
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::hard()).unwrap();
    let mut ws = assigner.workspace();
    let mut scratch = DagScratch::default();
    let mut beyond = 0;
    for read in sim.reads(40, 5).unwrap() {
        let call = assigner.assign(&encode(&read.sequence), &mut ws);
        let mut placements = top_by_kind(&call.placements, SegmentKind::V, 3);
        placements.extend(top_by_kind(&call.placements, SegmentKind::D, 6));
        placements.extend(top_by_kind(&call.placements, SegmentKind::J, 3));
        for iota in [0, -1] {
            let model = EngineeringObjective {
                insertion_per_base: iota,
                ..base
            };
            let dag = PlacementDag::new(VdjGrammar::VDJ, model, config.min_retained, Mode::Fixed);
            let (fast, exhaustive) = dag_and_brute(&dag, &placements, fixed(&dag), &mut scratch);
            assert_eq!(
                fast,
                exhaustive,
                "read {} at iota {iota}",
                read.sequence.len()
            );
            // Count reads whose optimum enters some placement beyond the overlap.
            if let Some(path) = dag.adjudicate(&placements, &mut scratch) {
                beyond += usize::from(path.segments.windows(2).any(|w| {
                    let (b, e) = (&placements[w[0].placement], &placements[w[1].placement]);
                    e.q_start() < b.q_end() && w[1].entry > b.q_end()
                }));
            }
        }
    }
    assert!(
        beyond > 0,
        "the reads should exercise entries beyond the overlap"
    );
}

/// The per-read model in integer units is an exact max-plus problem: on the
/// pipeline's placements for mid and hypermutated IgSim reads, under a tandem
/// grammar with and without a per-base N/P cost, the DAG's key equals the
/// exhaustive maximum at the read's own divergence and across the range, and
/// its traced path re-scores to it. At the log-odds divergence δ0 (with an
/// integer `P_D`) every value is the fixed score times `u0`, so the path is
/// the fixed-scale call's.
#[test]
fn per_read_dag_equals_exhaustive_enumeration() {
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let base = config
        .objective
        .resolve(Some(&db.d), &config.d, &db.j, &config.j.scoring);
    let evidence =
        EvidenceModel::new([&config.v.scoring, &config.d.scoring, &config.j.scoring]).unwrap();
    let delta0 = evidence.log_odds_divergence(SegmentKind::D);
    let mut ws = assigner.workspace();
    let mut scratch = DagScratch::default();
    let (mut checked, mut differs) = (0, 0);
    for preset in [SimParams::hard(), SimParams::mid()] {
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, preset).unwrap();
        for read in sim.reads(40, 5).unwrap() {
            let call = assigner.assign(&encode(&read.sequence), &mut ws);
            let mut placements = top_by_kind(&call.placements, SegmentKind::V, 3);
            placements.extend(top_by_kind(&call.placements, SegmentKind::D, 6));
            placements.extend(top_by_kind(&call.placements, SegmentKind::J, 3));
            for iota in [0, -1] {
                let objective = EngineeringObjective {
                    insertion_per_base: iota,
                    ..base
                };
                let dag = PlacementDag::new(TANDEM_D, objective, config.min_retained, Mode::Fixed);
                let Some((fixed_key, fixed_path)) =
                    dag.viterbi(&placements, &fixed(&dag), &mut scratch)
                else {
                    continue;
                };
                let own = divergence(&placements);
                for p in [own, 1e-3, 0.02, 0.3] {
                    let model = ReadModel::at(&evidence, p, objective.presence);
                    let (fast, exhaustive) = dag_and_brute(&dag, &placements, model, &mut scratch);
                    assert_eq!(fast, exhaustive, "p {p}, iota {iota}");
                    let (_, path) = dag.viterbi(&placements, &model, &mut scratch).unwrap();
                    differs += usize::from(path != fixed_path);
                    checked += 1;
                }
                // At δ0 the per-read model is the fixed one times u0.
                let integer = PlacementDag::new(
                    TANDEM_D,
                    EngineeringObjective {
                        presence: objective.presence.map(f64::floor),
                        ..objective
                    },
                    config.min_retained,
                    Mode::Fixed,
                );
                let model = ReadModel::at(&evidence, delta0, integer.model.presence);
                let (fast, exhaustive) = dag_and_brute(&integer, &placements, model, &mut scratch);
                assert_eq!(fast, exhaustive, "δ0, iota {iota}");
                let (key, path) = integer.viterbi(&placements, &model, &mut scratch).unwrap();
                assert_eq!(path, fixed_path, "the fixed-scale path at δ0");
                assert_eq!(
                    (key.score, key.d_count, key.support),
                    (
                        Units::from(fixed_key.score) * model.prior(1),
                        fixed_key.d_count,
                        fixed_key.support
                    )
                );
                checked += 1;
            }
        }
    }
    // The per-read optimum is often not the fixed-scale one, so the test
    // discriminates.
    assert!(
        checked >= 600 && differs >= 50,
        "{checked} optima, {differs} differ"
    );
}

/// The conditioned pass behind max-marginal D decoding is the same DAG on a
/// masked lattice under a restricted grammar: for every D placement of
/// hypermutated IgSim reads, the best path with the call's topology through
/// V, J and that D equals the exhaustive search over those placements under
/// [`Within`]; with every placement kept it is the unconditioned call.
#[test]
fn conditioned_dag_equals_exhaustive_enumeration() {
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let objective = config
        .objective
        .resolve(Some(&db.d), &config.d, &db.j, &config.j.scoring);
    let evidence =
        EvidenceModel::new([&config.v.scoring, &config.d.scoring, &config.j.scoring]).unwrap();
    let dag = PlacementDag::new(VdjGrammar::VDJ, objective, config.min_retained, Mode::Fixed);
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::hard()).unwrap();
    let mut ws = assigner.workspace();
    let mut scratch = DagScratch::default();
    let mut checked = 0;
    for read in sim.reads(40, 7).unwrap() {
        let call = assigner.assign(&encode(&read.sequence), &mut ws);
        let mut placements = top_by_kind(&call.placements, SegmentKind::V, 3);
        placements.extend(top_by_kind(&call.placements, SegmentKind::D, 6));
        placements.extend(top_by_kind(&call.placements, SegmentKind::J, 3));
        let model = ReadModel::at(&evidence, divergence(&placements), objective.presence);
        let Some((key, path)) = dag.viterbi(&placements, &model, &mut scratch) else {
            continue;
        };
        let kind = |s: &CalledSegment| placements[s.placement].kind;
        let topology = Topology {
            first: kind(&path[0]),
            last: kind(&path[path.len() - 1]),
            d_count: key.d_count,
        };
        let all = dag.conditioned(&placements, |_| true, topology, &model, &mut scratch);
        assert_eq!(all, Some((key, path.clone())));
        let within = Within {
            grammar: &dag.grammar,
            topology,
        };
        for d in (0..placements.len()).filter(|&i| placements[i].kind == SegmentKind::D) {
            let keep = |i: usize| placements[i].kind != SegmentKind::D || i == d;
            let fast = dag.conditioned(&placements, keep, topology, &model, &mut scratch);
            let kept: Vec<Placement> = (0..placements.len())
                .filter(|&i| keep(i))
                .map(|i| placements[i].clone())
                .collect();
            let exhaustive = Brute::within(&dag, &within, &kept, model).best();
            assert_eq!(fast.as_ref().map(|(key, _)| *key), exhaustive);
            if let Some((key, segments)) = fast {
                assert!(segments.iter().all(|s| keep(s.placement)));
                assert_eq!(key.d_count, topology.d_count);
                checked += 1;
            }
        }
    }
    assert!(checked >= 100, "{checked} conditioned optima");
}

fn top_by_kind(placements: &[Placement], kind: SegmentKind, n: usize) -> Vec<Placement> {
    let mut chosen: Vec<Placement> = placements
        .iter()
        .filter(|p| p.kind == kind)
        .cloned()
        .collect();
    chosen.sort_by_key(|p| std::cmp::Reverse(p.score()));
    chosen.truncate(n);
    chosen
}

/// The DAG under a two-D grammar, with a per-base N/P cost, against the
/// exhaustive search on the pipeline's placements of reads with no D, one D
/// and two (VJ, VDJ and VDDJ optima all occur).
#[test]
fn dag_equals_brute_force_over_all_paths_and_cuts() {
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let model = EngineeringObjective {
        insertion_per_base: -1,
        ..config
            .objective
            .resolve(Some(&db.d), &config.d, &db.j, &config.j.scoring)
    };
    let dag = PlacementDag::new(TANDEM_D, model, [10, 3, 6], Mode::Fixed);
    let params = SimParams {
        include_d: 0.7,
        include_dd: 0.5,
        ..SimParams::easy()
    }
    .with_shm(0.04);
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params).unwrap();
    let mut ws = assigner.workspace();
    let mut scratch = DagScratch::default();
    let mut distinct_topologies = HashSet::new();
    for read in sim.reads(40, 23).unwrap() {
        let call = assigner.assign(&encode(&read.sequence), &mut ws);
        let mut placements = top_by_kind(&call.placements, SegmentKind::V, 4);
        placements.extend(top_by_kind(&call.placements, SegmentKind::D, 6));
        placements.extend(top_by_kind(&call.placements, SegmentKind::J, 3));

        let (fast, exhaustive) = dag_and_brute(&dag, &placements, fixed(&dag), &mut scratch);
        assert_eq!(fast, exhaustive);
        if let Some(path) = dag.adjudicate(&placements, &mut scratch) {
            // On the fixed scale the reported path score is the DAG's key.
            let key = fast.expect("a path");
            assert_eq!(
                (
                    Units::from(path.score.score),
                    path.score.d_count,
                    path.score.support
                ),
                (key.score, key.d_count, key.support)
            );
            let kinds: Vec<SegmentKind> = path
                .segments
                .iter()
                .map(|s| placements[s.placement].kind)
                .collect();
            distinct_topologies.insert(kinds);
            let clip_sum: i32 = path.segments.iter().map(|s| s.span.score).sum();
            assert!(
                clip_sum >= path.score.score - 24,
                "path score is built from clip scores"
            );
        }
    }
    assert!(
        distinct_topologies.len() >= 3,
        "test should exercise VJ, VDJ and VDDJ paths"
    );
}
