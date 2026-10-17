//! Exhaustive oracles for the placement DAG on the pipeline's own
//! placements: the Viterbi key under the fixed and the per-read model, and
//! the conditioned pass behind max-marginal D decoding, against a memoised
//! search over every grammatical path and every pair of canonical cuts.

use std::collections::{HashMap, HashSet};

use super::decode::{Topology, Within};
use super::units::{Columns, DTerms, Flank, FlankTable, ReadContext, ReadModel};
use super::{
    Adjudicator, CallModel, CallPath, CalledSegment, DagScratch, EngineeringObjective,
    EvidenceModel, Grammar, Mode, PathModel, PathScore, PlacementDag, Units, VdjGrammar,
    divergence,
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
        let span = p
            .clip
            .span(entry, exit)
            .expect("canonical cuts keep a column");
        value + self.read.d_span(p, &span)
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
            let on = DTerms::geometric(Some((5.0, 7.0)), Some(8.0));
            for (iota, d_terms) in [
                (0, DTerms::default()),
                (-1, DTerms::default()),
                (0, on),
                (-1, on),
            ] {
                let objective = EngineeringObjective {
                    insertion_per_base: iota,
                    d_terms,
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
                    let model = ReadModel::at_with(&evidence, p, objective.presence, d_terms);
                    let (fast, exhaustive) = dag_and_brute(&dag, &placements, model, &mut scratch);
                    assert_eq!(fast, exhaustive, "p {p}, iota {iota}, {d_terms:?}");
                    let (_, path) = dag
                        .viterbi(&placements, &model, &FlankTable::default(), &mut scratch)
                        .unwrap();
                    differs += usize::from(path != fixed_path);
                    checked += 1;
                }
                // At δ0 the per-read model is the fixed one times u0
                // (without the D terms, which are fractions of a unit).
                if !d_terms.is_off() {
                    continue;
                }
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
    // The placements of the unguarded engine: the pair guards' hit floors
    // would drop the short V and J that cannot pay for themselves, which
    // this test feeds the DAG ([`FlankTable::closes`]).
    let mut config = EngineConfig::default();
    config.options.guards = None;
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

/// The Combine rule (Swig 0.39.5), stated here on its own rather than
/// taken from [`crate::pipeline`]: a path keeps the guards when
/// * each of its V and J keeps SwiftIG's score and identity floor where the
///   path clips it (`floor`), and
/// * it holds a V and a J that pair (`pairs`, SwiftIG's test) with at most
///   `span` read bases between the V's exit and the J's entry; or
/// * it holds no D; or
/// * it holds a V and no J, and every D placement ends at most `reach`
///   bases after the V placement's end (the D window of a read without a
///   J); or a J and no V, and every D placement starts at most `reach`
///   bases before the J placement's start.
///
/// A D with neither a V nor a J breaks it.
struct CombineRule<'a, F, P> {
    floor: F,
    pairs: P,
    span: u32,
    reach: u32,
    query: &'a [u8],
}

impl<F, P> CombineRule<'_, F, P>
where
    F: Fn(&[u8], &[Placement], &CalledSegment) -> bool,
    P: Fn(&Placement, &Placement) -> bool,
{
    fn keeps(&self, placed: &[Placement], segments: &[CalledSegment]) -> bool {
        let of = |kind| segments.iter().find(|s| placed[s.placement].kind == kind);
        let flanks = segments
            .iter()
            .filter(|s| matches!(placed[s.placement].kind, SegmentKind::V | SegmentKind::J));
        if !flanks.clone().all(|s| (self.floor)(self.query, placed, s)) {
            return false;
        }
        let mut ds = segments
            .iter()
            .map(|s| &placed[s.placement])
            .filter(|p| p.kind == SegmentKind::D);
        match (of(SegmentKind::V), of(SegmentKind::J)) {
            (Some(v), Some(j)) => {
                j.span.q_lo.saturating_sub(v.span.q_hi) <= self.span
                    && (self.pairs)(&placed[v.placement], &placed[j.placement])
            }
            (Some(v), None) => {
                let end = placed[v.placement].q_end().saturating_add(self.reach);
                ds.all(|d| d.q_end() <= end)
            }
            (None, Some(j)) => {
                let start = placed[j.placement].q_start().saturating_sub(self.reach);
                ds.all(|d| d.q_start() >= start)
            }
            (None, None) => ds.next().is_none(),
        }
    }
}

/// Reads of the shapes the pair guards and the open-flank cost act on:
/// whole, cut two bases after their D (V–D), from ten bases before their D
/// (D–J), a V followed by 100 random bases, 100 random bases followed by a
/// J, and a V and D with their J beyond the junction span.
fn guarded_fixtures(db: &SyntheticDatabase) -> Vec<Vec<u8>> {
    use crate::sim::Region;
    let mut rng = crate::sim::Rng::new(41);
    let mut reads = Vec::new();
    for (preset, seed) in [(SimParams::mid(), 43), (SimParams::hard(), 44)] {
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, preset).unwrap();
        for read in sim.reads(40, seed).unwrap() {
            let f = &read.forward;
            let (Some(v), Some(d), Some(j)) = (
                read.span(Region::V),
                read.span(Region::D),
                read.span(Region::J),
            ) else {
                continue;
            };
            reads.push(f.clone());
            reads.push(f[..(d.end + 2).min(f.len())].to_vec());
            reads.push(f[d.start.saturating_sub(10)..].to_vec());
            reads.push([&f[..v.end], &rng.bases(100)[..]].concat());
            reads.push([&rng.bases(100)[..], &f[j.start..]].concat());
            let gap = 230usize.saturating_sub(d.end - v.end);
            reads.push([&f[..d.end], &rng.bases(gap)[..], &f[j.start..]].concat());
        }
    }
    reads
}

/// The pair guards with the open-flank cost, under both call models with
/// posteriors: on every read of [`guarded_fixtures`],
/// * the call keeps the guards (a D beside a V or a J at the read's open
///   end is allowed and charged);
/// * the call's key re-scores, on all the read's placements and in the
///   read's context, to the exhaustive oracle's value of its path, and its
///   fixed score is the oracle's fixed-scale value: the DAG, the oracle and
///   `path_score` agree on the guarded call;
/// * it never beats the unguarded optimum, and equals it when that keeps
///   the guards; otherwise it is the best of the guard subsets' optima
///   (each the DAG's and the exhaustive search's alike) that hold an
///   anchor and keep the guards, the earlier subset on an exact tie;
/// * its posterior is the forward–backward pass over the subset that gave
///   the call, in the same context, so that under the per-read model the
///   call is a mode of the reported distribution;
/// * the engine's guard verdict on the call, the unguarded optimum and
///   every subset optimum is the verdict of [`CombineRule`], the rule
///   written out independently;
/// * a path the rule allows rarely beats the call: for every choice of at
///   most one V, one D and one J among the read's best placements (4 V,
///   8 D, 4 J), the optimum over that choice, when it holds an anchor and
///   keeps [`CombineRule`], is compared with the call. The subsets are a
///   search heuristic (a subset whose optimum breaks the guards is passed
///   over whole), so this measures what they miss: at most 2 calls of the
///   fixtures' per model (measured: 0 per-read, 1 fixed-scale).
#[test]
fn guarded_calls_with_open_flanks_equal_the_exhaustive_subset_optima() {
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let reads = guarded_fixtures(&db);
    let mut posterior = super::PosteriorScratch::default();
    let mut scratch = DagScratch::default();
    let (mut table, mut fixed_table) = (FlankTable::default(), FlankTable::default());
    for calls in [CallModel::PerRead, CallModel::FixedScale] {
        let config = EngineConfig {
            calls,
            ..EngineConfig::default()
        };
        let assigner =
            DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
                .unwrap();
        let dag = &assigner.adjudicator;
        let flank = dag.model.open_flank;
        assert!(flank.is_some() && dag.charges_open_flank());
        let mut ws = assigner.workspace();
        let (mut fallbacks, mut open_ends, mut open_subsets) = (0, 0, 0);
        let (mut skeletons, mut allowed, mut missed, mut calls_checked) = (0, 0, 0, 0);
        let floor =
            |q: &[u8], p: &[Placement], s: &CalledSegment| assigner.segment_keeps_floor(q, p, s);
        let pairs = |v: &Placement, j: &Placement| assigner.pair_check(v, j);
        for read in &reads {
            let call = assigner.assign(&encode(read), &mut ws);
            let placed = &call.placements;
            let (query, path) = (&call.query, call.path.as_ref());
            let rule = CombineRule {
                floor: &floor,
                pairs: &pairs,
                span: config.options.guards.unwrap().max_junction_span,
                reach: config.d_window.reach,
                query,
            };
            let (subsets, keeps) = assigner.guard_check(query, placed, path);
            assert!(keeps, "the call breaks the guards");
            let Some(path) = path else {
                continue;
            };
            assert!(
                rule.keeps(placed, &path.segments),
                "the call breaks the rule"
            );
            let context = ReadContext::of(placed, query.len() as u32);
            assert_eq!(path.context, context);
            let read_model = dag.read_model(context.divergence);
            let fill = |table: &mut FlankTable,
                        model: &ReadModel,
                        calls: CallModel,
                        placements: &[Placement]| {
                table.fill(
                    flank.as_ref(),
                    dag.model.presence,
                    model,
                    calls,
                    &context,
                    placements,
                );
            };
            // The call's value and fixed score, re-scored on every placement.
            fill(&mut table, &read_model, calls, placed);
            let mut all = Brute::new(dag, placed, read_model, &table);
            assert_eq!(all.rescore(&path.segments), path.value.score);
            let unguarded = all.best().expect("a call has a path");
            assert!(path.value <= unguarded);
            fill(&mut fixed_table, &fixed(dag), CallModel::FixedScale, placed);
            let oracle = Brute::new(dag, placed, fixed(dag), &fixed_table).rescore(&path.segments);
            assert_eq!(Units::from(path.score.score), oracle);
            let kind = |s: &CalledSegment| placed[s.placement].kind;
            let (first, last) = (kind(&path.segments[0]), kind(path.segments.last().unwrap()));
            let open = (first == SegmentKind::D) || (last == SegmentKind::D);
            open_ends += usize::from(open);
            // Which subset the call came from: all placements when their
            // optimum keeps the guards, else the best of the guard subsets.
            let traced =
                dag.viterbi(placed, &read_model, &table, &mut scratch)
                    .map(|(value, segments)| CallPath {
                        score: dag.path_score(placed, &segments, &table),
                        segments,
                        value,
                        context,
                        model: calls,
                    });
            let (_, all_keep) = assigner.guard_check(query, placed, traced.as_ref());
            if let Some(t) = &traced {
                assert_eq!(
                    all_keep,
                    rule.keeps(placed, &t.segments),
                    "unguarded optimum"
                );
            }
            let anchored = |placements: &[Placement], segments: &[CalledSegment]| {
                segments
                    .iter()
                    .any(|s| placements[s.placement].anchors(&s.span))
            };
            let chosen: Vec<usize> = if all_keep
                && traced
                    .as_ref()
                    .is_some_and(|t| anchored(placed, &t.segments))
            {
                assert_eq!(traced.as_ref().unwrap().value, path.value);
                (0..placed.len()).collect()
            } else {
                fallbacks += 1;
                let mut best: Option<(usize, PathScore<Units>)> = None;
                for (k, ids) in subsets.iter().enumerate() {
                    let subset: Vec<Placement> = ids.iter().map(|&i| placed[i].clone()).collect();
                    fill(&mut table, &read_model, calls, &subset);
                    let (fast, exhaustive) =
                        dag_and_brute_with(dag, &subset, read_model, &table, &mut scratch);
                    assert_eq!(fast, exhaustive, "subset {ids:?}");
                    let Some((value, segments)) =
                        dag.viterbi(&subset, &read_model, &table, &mut scratch)
                    else {
                        continue;
                    };
                    let candidate = CallPath {
                        score: dag.path_score(&subset, &segments, &table),
                        segments,
                        value,
                        context,
                        model: calls,
                    };
                    let verdict = assigner.guard_check(query, &subset, Some(&candidate)).1;
                    assert_eq!(
                        verdict,
                        rule.keeps(&subset, &candidate.segments),
                        "subset {ids:?}"
                    );
                    if !anchored(&subset, &candidate.segments) || !verdict {
                        continue;
                    }
                    if best.is_none_or(|(_, b)| value > b) {
                        best = Some((k, value));
                    }
                }
                let (k, value) = best.expect("the call is a subset optimum");
                assert_eq!(value, path.value);
                // Only the open-end subsets hold a D without both a V and a J.
                let has = |kind| subsets[k].iter().any(|&i| placed[i].kind == kind);
                open_subsets += usize::from(
                    has(SegmentKind::D) && !(has(SegmentKind::V) && has(SegmentKind::J)),
                );
                subsets[k].clone()
            };
            // No path the rule allows beats the call, over every choice of
            // at most one V, D and J among the best placements.
            let best_ids = |kind: SegmentKind, n: usize| {
                let mut ids: Vec<usize> = (0..placed.len())
                    .filter(|&i| placed[i].kind == kind)
                    .collect();
                ids.sort_by_key(|&i| std::cmp::Reverse(placed[i].score()));
                ids.truncate(n);
                let mut options: Vec<Option<usize>> = ids.into_iter().map(Some).collect();
                options.push(None);
                options
            };
            let (vs, ds, js) = (
                best_ids(SegmentKind::V, 4),
                best_ids(SegmentKind::D, 8),
                best_ids(SegmentKind::J, 4),
            );
            let mut beaten = false;
            for &v in &vs {
                for &d in &ds {
                    for &j in &js {
                        if v.is_none() && j.is_none() {
                            continue;
                        }
                        let ids: Vec<usize> = [v, d, j].into_iter().flatten().collect();
                        let choice: Vec<Placement> =
                            ids.iter().map(|&i| placed[i].clone()).collect();
                        fill(&mut table, &read_model, calls, &choice);
                        let Some((value, segments)) =
                            dag.viterbi(&choice, &read_model, &table, &mut scratch)
                        else {
                            continue;
                        };
                        skeletons += 1;
                        if !anchored(&choice, &segments) || !rule.keeps(&choice, &segments) {
                            continue;
                        }
                        allowed += 1;
                        beaten |= value > path.value;
                    }
                }
            }
            missed += usize::from(beaten);
            calls_checked += 1;
            // The reported posterior is the chosen subset's.
            let subset: Vec<Placement> = chosen.iter().map(|&i| placed[i].clone()).collect();
            fill(&mut table, &read_model, calls, &subset);
            let nats = dag.nats(context.divergence).expect("posteriors are on");
            let expected = dag
                .posterior(&subset, &nats, &table, &mut posterior)
                .expect("the subset has a path");
            let reported = call.marginals.as_ref().expect("posteriors are on");
            let mut mapped = vec![0.0; placed.len()];
            for (&id, &p) in chosen.iter().zip(&expected.placement) {
                mapped[id] = p;
            }
            for (a, b) in reported.placement.iter().zip(&mapped) {
                assert!((a - b).abs() <= 1e-12, "{a} != {b}");
            }
            assert_eq!(reported.d_count.len(), expected.d_count.len());
            for (a, b) in reported.d_count.iter().zip(&expected.d_count) {
                assert!((a - b).abs() <= 1e-12, "{a} != {b}");
            }
        }
        assert!(
            fallbacks >= 40 && open_ends >= 100 && open_subsets >= 5,
            "{calls:?}: {fallbacks} fallbacks, {open_ends} open-ended calls, \
             {open_subsets} from an open-end subset"
        );
        eprintln!(
            "{calls:?}: {skeletons} choices, {allowed} allowed by the rule; \
             on {missed} of {calls_checked} calls one beats the call"
        );
        // Measured: none under the per-read model, one call under the fixed
        // scale, where the optimum of the V–J pair with every D clips the V
        // so that the pair's span exceeds 225 bases, and the subset is
        // passed over although a V–D–J path with another D and the V's
        // full exit keeps it (Swig 0.39.4's subsets, the D-free pair wins).
        assert!(
            allowed >= 1000 && missed <= 2,
            "{calls:?}: on {missed} of {calls_checked} calls an allowed path beats the call"
        );
    }
}

/// The two switches of [`guarded_calls_with_open_flanks_equal_the_exhaustive_subset_optima`]:
/// without the open-flank cost a D is called only between a V and a J
/// (Swig 0.39.4's rule), and without the pair guards the call is the
/// exhaustive optimum over every placement, the cost included.
#[test]
fn without_the_cost_or_the_guards() {
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let reads = guarded_fixtures(&db);
    let mut table = FlankTable::default();
    for calls in [CallModel::PerRead, CallModel::FixedScale] {
        let mut closed = EngineConfig {
            calls,
            ..EngineConfig::default()
        };
        closed.objective.open_flank = None;
        let mut free = EngineConfig {
            calls,
            ..EngineConfig::default()
        };
        free.options.guards = None;
        let build = |config: &EngineConfig| {
            DefaultAssigner::from_config(config, db.v.clone(), Some(db.d.clone()), db.j.clone())
                .unwrap()
        };
        let (closed, free) = (build(&closed), build(&free));
        assert!(!closed.adjudicator.charges_open_flank());
        let (mut ws, mut free_ws) = (closed.workspace(), free.workspace());
        let (mut with_d, mut open_ends) = (0, 0);
        for read in &reads {
            let call = closed.assign(&encode(read), &mut ws);
            if let Some(path) = &call.path {
                let kinds: Vec<SegmentKind> = path
                    .segments
                    .iter()
                    .map(|s| call.placements[s.placement].kind)
                    .collect();
                if kinds.contains(&SegmentKind::D) {
                    with_d += 1;
                    assert_eq!(kinds.first(), Some(&SegmentKind::V));
                    assert_eq!(kinds.last(), Some(&SegmentKind::J));
                }
            }
            let call = free.assign(&encode(read), &mut free_ws);
            let Some(path) = &call.path else {
                continue;
            };
            let dag = &free.adjudicator;
            let context = ReadContext::of(&call.placements, call.query.len() as u32);
            let model = dag.read_model(context.divergence);
            table.fill(
                dag.model.open_flank.as_ref(),
                dag.model.presence,
                &model,
                calls,
                &context,
                &call.placements,
            );
            let best = Brute::new(dag, &call.placements, model, &table).best();
            assert_eq!(best, Some(path.value));
            let kind = |s: &CalledSegment| call.placements[s.placement].kind;
            open_ends += usize::from(
                kind(&path.segments[0]) == SegmentKind::D
                    || kind(path.segments.last().unwrap()) == SegmentKind::D,
            );
        }
        assert!(
            with_d >= 50 && open_ends >= 100,
            "{with_d} D calls, {open_ends} open ends"
        );
    }
}
