//! Exhaustive oracles for the placement DAG on the pipeline's own
//! placements: the Viterbi key under the fixed and the per-read model, and
//! the conditioned pass behind max-marginal D decoding, against a memoised
//! search over every grammatical path and every pair of canonical cuts.

use std::collections::{HashMap, HashSet};

use super::decode::{Topology, Within};
use super::units::{Columns, Flank, FlankTable, ReadContext, ReadModel};
use super::{
    Adjudicator, CallModel, CalledSegment, DagScratch, EngineeringObjective, EvidenceModel,
    Grammar, Mode, PathModel, PathScore, PlacementDag, Units, VdjGrammar, divergence,
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
/// term under `read` and the open-flank costs `flank`. Memoised on
/// `(placement, entry, D count, first)`.
struct Brute<'a, G = VdjGrammar> {
    dag: &'a PlacementDag,
    grammar: &'a G,
    placements: &'a [Placement],
    read: ReadModel,
    flank: &'a FlankTable,
    memo: HashMap<(usize, u32, u8, SegmentKind), Option<PathScore<Units>>>,
}

/// Canonical entry and exit cuts of a placement (one per Match column).
fn cuts(p: &Placement) -> (Vec<u32>, Vec<u32>) {
    p.clip.matches().map(|(q, ..)| (q, q + 1)).unzip()
}

impl<'a> Brute<'a> {
    fn new(
        dag: &'a PlacementDag,
        placements: &'a [Placement],
        read: ReadModel,
        flank: &'a FlankTable,
    ) -> Self {
        Self::within(dag, &dag.grammar, placements, read, flank)
    }
}

impl<'a, G: Grammar> Brute<'a, G> {
    /// The search under another grammar than `dag`'s.
    fn within(
        dag: &'a PlacementDag,
        grammar: &'a G,
        placements: &'a [Placement],
        read: ReadModel,
        flank: &'a FlankTable,
    ) -> Self {
        Self {
            dag,
            grammar,
            placements,
            read,
            flank,
            memo: HashMap::new(),
        }
    }

    /// The open-flank cost of a path that starts or ends at `p`.
    fn open(&self, side: Flank, p: &Placement) -> Units {
        self.flank.units(side, p.kind, p.q_start(), p.q_end())
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
                    + self.read.prior(model.end(first, pb.kind, d_used))
                    - self.open(Flank::End, pb),
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
            // A V that cannot open a D's flank, a J that cannot close one.
            let closes = self.flank.closes();
            let blocked = !closes.is_empty()
                && match (pb.kind, pe.kind) {
                    (SegmentKind::V, SegmentKind::D) => !closes[b],
                    (SegmentKind::D, SegmentKind::J) => !closes[e],
                    _ => false,
                };
            if blocked {
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
                let start = self.open(Flank::Start, p);
                let rest =
                    self.best_from(s, p.q_start(), u8::from(p.kind == SegmentKind::D), p.kind)?;
                Some(PathScore {
                    score: rest.score - start,
                    ..rest
                })
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
        let start = self.open(Flank::Start, &self.placements[segments[0].placement]);
        let end = self.open(
            Flank::End,
            &self.placements[segments[segments.len() - 1].placement],
        );
        value + self.read.prior(self.dag.model.end(first, last, d_used)) - start - end
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
    dag_and_brute_with(dag, placements, read, &FlankTable::default(), scratch)
}

/// [`dag_and_brute`] with the open-flank costs of `flank`.
fn dag_and_brute_with(
    dag: &PlacementDag,
    placements: &[Placement],
    read: ReadModel,
    flank: &FlankTable,
    scratch: &mut DagScratch,
) -> (Option<PathScore<Units>>, Option<PathScore<Units>>) {
    let fast = dag.viterbi(placements, &read, flank, scratch);
    let mut brute = Brute::new(dag, placements, read, flank);
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
        .resolve(&db.v, &config.v, Some(&db.d), &config.d, &db.j, &config.j);
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
        .resolve(&db.v, &config.v, Some(&db.d), &config.d, &db.j, &config.j);
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
                let Some((fixed_key, fixed_path)) = dag.viterbi(
                    &placements,
                    &fixed(&dag),
                    &FlankTable::default(),
                    &mut scratch,
                ) else {
                    continue;
                };
                let own = divergence(&placements);
                for p in [own, 1e-3, 0.02, 0.3] {
                    let model = ReadModel::at(&evidence, p, objective.presence);
                    let (fast, exhaustive) = dag_and_brute(&dag, &placements, model, &mut scratch);
                    assert_eq!(fast, exhaustive, "p {p}, iota {iota}");
                    let (_, path) = dag
                        .viterbi(&placements, &model, &FlankTable::default(), &mut scratch)
                        .unwrap();
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
                let (key, path) = integer
                    .viterbi(&placements, &model, &FlankTable::default(), &mut scratch)
                    .unwrap();
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
    let objective =
        config
            .objective
            .resolve(&db.v, &config.v, Some(&db.d), &config.d, &db.j, &config.j);
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
        let none = FlankTable::default();
        let Some((key, path)) = dag.viterbi(&placements, &model, &none, &mut scratch) else {
            continue;
        };
        let kind = |s: &CalledSegment| placements[s.placement].kind;
        let topology = Topology {
            first: kind(&path[0]),
            last: kind(&path[path.len() - 1]),
            d_count: key.d_count,
        };
        let all = dag.conditioned(&placements, |_| true, topology, &model, &none, &mut scratch);
        assert_eq!(all, Some((key, path.clone())));
        let within = Within {
            grammar: &dag.grammar,
            topology,
        };
        for d in (0..placements.len()).filter(|&i| placements[i].kind == SegmentKind::D) {
            let keep = |i: usize| placements[i].kind != SegmentKind::D || i == d;
            let fast = dag.conditioned(&placements, keep, topology, &model, &none, &mut scratch);
            let kept: Vec<Placement> = (0..placements.len())
                .filter(|&i| keep(i))
                .map(|i| placements[i].clone())
                .collect();
            let exhaustive = Brute::within(&dag, &within, &kept, model, &none).best();
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
            .resolve(&db.v, &config.v, Some(&db.d), &config.d, &db.j, &config.j)
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

/// The open-flank cost at a path's ends ([`FlankTable`]) on the pipeline's
/// placements of mid and hypermutated reads cut 0–80 bases before and after
/// a D, so that paths start or end on a D with or without a V or J, at
/// random anchored columns, at two allowances, with the cost on and
/// off: the DAG's key equals the exhaustive search under the per-read and
/// the fixed model, its traced path re-scores to it, and `path_score`
/// reports the oracle's fixed-scale value of that path. Some reads hold a
/// V or J that cannot pay for itself, which then neither opens nor closes
/// a D's flank ([`FlankTable::closes`]).
#[test]
fn open_flank_dag_equals_exhaustive_enumeration() {
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let base = config
        .objective
        .resolve(&db.v, &config.v, Some(&db.d), &config.d, &db.j, &config.j);
    let resolved = base
        .open_flank
        .expect("the open-flank cost is on by default");
    let evidence =
        EvidenceModel::new([&config.v.scoring, &config.d.scoring, &config.j.scoring]).unwrap();
    let mut rng = crate::sim::Rng::new(29);
    let mut ws = assigner.workspace();
    let mut scratch = DagScratch::default();
    let (mut table, mut fixed_table) = (FlankTable::default(), FlankTable::default());
    let (mut checked, mut charged, mut differs, mut blocked) = (0, 0, 0, 0);
    for preset in [SimParams::mid(), SimParams::hard()] {
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, preset).unwrap();
        for read in sim.reads(30, 13).unwrap() {
            let call = assigner.assign(&encode(&read.sequence), &mut ws);
            let Some(d) = top_by_kind(&call.placements, SegmentKind::D, 1).pop() else {
                continue;
            };
            let (lead, trail) = (rng.below(81) as u32, rng.below(81) as u32);
            let lo = d.q_start().saturating_sub(lead) as usize;
            let hi = (d.q_end() + trail).min(call.query.len() as u32) as usize;
            let cut = assigner.assign(&call.query[lo..hi], &mut ws);
            let mut all = top_by_kind(&cut.placements, SegmentKind::V, 3);
            all.extend(top_by_kind(&cut.placements, SegmentKind::D, 6));
            all.extend(top_by_kind(&cut.placements, SegmentKind::J, 3));
            // With and without the V or the J: paths that start or end on a D.
            for dropped in [None, Some(SegmentKind::V), Some(SegmentKind::J)] {
                let placements: Vec<Placement> = all
                    .iter()
                    .filter(|p| Some(p.kind) != dropped)
                    .cloned()
                    .collect();
                if !placements.iter().any(|p| p.kind == SegmentKind::D) {
                    continue;
                }
                let mut context = ReadContext::of(&placements, cut.query.len() as u32);
                let columns = 30 + rng.below(120);
                context.anchored = Columns {
                    columns,
                    mismatches: rng.below(columns * 2 / 5),
                };
                for n_open in [24, 12] {
                    let flank = super::ResolvedOpenFlank { n_open, ..resolved };
                    let objective = EngineeringObjective {
                        open_flank: Some(flank),
                        ..base
                    };
                    let dag =
                        PlacementDag::new(TANDEM_D, objective, config.min_retained, Mode::Fixed);
                    for calls in [CallModel::PerRead, CallModel::FixedScale] {
                        let read = match calls {
                            CallModel::PerRead => {
                                ReadModel::at(&evidence, context.divergence, objective.presence)
                            }
                            CallModel::FixedScale => fixed(&dag),
                        };
                        let mut keys = [None, None];
                        for (k, on) in [true, false].into_iter().enumerate() {
                            let model = on.then_some(&flank);
                            table.fill(
                                model,
                                objective.presence,
                                &read,
                                calls,
                                &context,
                                &placements,
                            );
                            blocked += usize::from(table.closes().iter().any(|&c| !c));
                            let (fast, exhaustive) =
                                dag_and_brute_with(&dag, &placements, read, &table, &mut scratch);
                            assert_eq!(fast, exhaustive, "{n_open} {calls:?} on {on}");
                            keys[k] = fast;
                            let Some((_, segments)) =
                                dag.viterbi(&placements, &read, &table, &mut scratch)
                            else {
                                continue;
                            };
                            checked += 1;
                            // The reported fixed score is the oracle's fixed value.
                            fixed_table.fill(
                                model,
                                objective.presence,
                                &fixed(&dag),
                                CallModel::FixedScale,
                                &context,
                                &placements,
                            );
                            let oracle = Brute::new(&dag, &placements, fixed(&dag), &fixed_table)
                                .rescore(&segments);
                            let reported = dag.path_score(&placements, &segments, &table);
                            assert_eq!(Units::from(reported.score), oracle);
                            let ends = [
                                (Flank::Start, &segments[0]),
                                (Flank::End, &segments[segments.len() - 1]),
                            ];
                            charged += usize::from(ends.iter().any(|(side, s)| {
                                let p = &placements[s.placement];
                                table.units(*side, p.kind, p.q_start(), p.q_end()) > 0
                            }));
                        }
                        differs += usize::from(keys[0] != keys[1]);
                    }
                }
            }
        }
    }
    assert!(
        checked >= 600 && charged >= 50 && differs >= 20 && blocked >= 20,
        "{checked} optima, {charged} charged, {differs} changed by the cost, {blocked} blocking"
    );
}
