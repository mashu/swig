//! End-to-end tests on a synthetic germline database.
//!
//! The database is random but deterministic (SplitMix64): 6 V genes × 2
//! alleles (SNPs at 50/150/250), 8 D genes of 14–21 nt, 4 J genes of 48 nt.
//! Queries are simulated rearrangements with trimming and random N regions.

use std::collections::HashSet;

use vdj_assign::adjudicate::{Adjudicator, DagScratch, PathScore};
use vdj_assign::placement::Placement;
use vdj_assign::prelude::*;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn range(&mut self, lo: usize, hi: usize) -> usize {
        lo + self.below(hi - lo + 1)
    }

    fn bases(&mut self, len: usize) -> String {
        (0..len).map(|_| char::from(b"ACGT"[self.below(4)])).collect()
    }
}

fn substitute(base: char, rng: &mut Rng) -> char {
    loop {
        let b = char::from(b"ACGT"[rng.below(4)]);
        if b != base {
            return b;
        }
    }
}

struct Db {
    v: Vec<(String, String)>,
    d: Vec<(String, String)>,
    j: Vec<(String, String)>,
}

const V_LEN: usize = 290;
const V_ANCHOR: usize = 281;
const J_ANCHOR: usize = 15;

impl Db {
    fn new() -> Self {
        let mut rng = Rng(2024);
        let mut v = Vec::new();
        for g in 1..=6 {
            let first = rng.bases(V_LEN);
            let mut second: Vec<char> = first.chars().collect();
            for p in [50, 150, 250] {
                second[p] = substitute(second[p], &mut rng);
            }
            v.push((format!("IGHV{g}-1*01"), first));
            v.push((format!("IGHV{g}-1*02"), second.into_iter().collect()));
        }
        let d = (1..=8).map(|g| (format!("IGHD{g}-1*01"), rng.bases(13 + g))).collect();
        let j = (1..=4).map(|g| (format!("IGHJ{g}*01"), rng.bases(48))).collect();
        Self { v, d, j }
    }

    fn set(kind: SegmentKind, entries: &[(String, String)], anchor: Option<u32>) -> GermlineSet {
        let mut set = GermlineSet::new(kind);
        for (name, seq) in entries {
            set.push(name, seq.as_bytes(), anchor).unwrap();
        }
        set
    }

    fn assigner(&self, config: &EngineConfig) -> DefaultAssigner {
        DefaultAssigner::from_config(
            config,
            Self::set(SegmentKind::V, &self.v, Some(V_ANCHOR as u32)),
            Some(Self::set(SegmentKind::D, &self.d, None)),
            Self::set(SegmentKind::J, &self.j, Some(J_ANCHOR as u32)),
        )
        .unwrap()
    }
}

struct Sim {
    query: String,
    v: String,
    d: Vec<String>,
    j: String,
}

/// V[10..] − trim, then for each chosen D: N + trimmed D, then N + trimmed J.
fn simulate(db: &Db, rng: &mut Rng, d_choice: &[usize], shm_rate: f64) -> Sim {
    let (v_name, v_seq) = &db.v[rng.below(db.v.len())];
    let (j_name, j_seq) = &db.j[rng.below(db.j.len())];
    let v_part: String = v_seq[10..V_LEN - rng.range(0, 3)]
        .chars()
        .map(|b| if (rng.below(10_000) as f64) < shm_rate * 10_000.0 { substitute(b, rng) } else { b })
        .collect();
    let mut query = v_part;
    let mut d_names = Vec::new();
    for &d in d_choice {
        let (d_name, d_seq) = &db.d[d];
        let (left, right) = (rng.range(0, 2), rng.range(0, 2));
        let n_len = rng.range(3, 6);
        query += &rng.bases(n_len);
        query += &d_seq[left..d_seq.len() - right];
        d_names.push(d_name.clone());
    }
    let n_len = rng.range(3, 6);
    query += &rng.bases(n_len);
    query += &j_seq[rng.range(0, 3)..];
    Sim { query, v: v_name.clone(), d: d_names, j: j_name.clone() }
}

fn first_label(report: &Option<vdj_assign::report::SegmentReport>) -> Option<String> {
    report.as_ref().map(|r| r.call.split(',').next().unwrap().to_string())
}

fn gene(label: &str) -> &str {
    label.split('*').next().unwrap()
}

fn record<A: VdjAssigner + AlleleCatalog>(assigner: &A, query: &str) -> AirrRecord {
    let mut ws = assigner.workspace();
    let call = assigner.assign(&encode(query.as_bytes()), &mut ws);
    AirrRecord::from_call("q", query, &call, assigner, &[None, None, None])
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

#[test]
fn calls_vdj_alleles_on_clean_rearrangements() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let mut rng = Rng(1);
    let (mut d_correct, trials) = (0, 60);
    for _ in 0..trials {
        let d = rng.below(db.d.len());
        let sim = simulate(&db, &mut rng, &[d], 0.0);
        let rec = record(&assigner, &sim.query);
        assert!(!rec.rev_comp);
        assert_eq!(first_label(&rec.v).as_deref(), Some(sim.v.as_str()), "V allele for {}", sim.query);
        assert_eq!(first_label(&rec.j).as_deref(), Some(sim.j.as_str()), "J allele");
        if first_label(&rec.d).as_deref() == Some(sim.d[0].as_str()) {
            d_correct += 1;
        }
        assert!(rec.d2.is_none(), "default grammar allows one D");
    }
    assert!(d_correct * 100 >= trials * 95, "D accuracy {d_correct}/{trials}");
}

#[test]
fn reverse_complement_reads_are_reoriented() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let mut rng = Rng(7);
    for _ in 0..20 {
        let d = rng.below(8);
        let sim = simulate(&db, &mut rng, &[d], 0.0);
        let forward = record(&assigner, &sim.query);
        let reverse = record(&assigner, &reverse_complement(&sim.query));
        assert!(reverse.rev_comp);
        assert_eq!(forward.v.as_ref().map(|s| &s.call), reverse.v.as_ref().map(|s| &s.call));
        assert_eq!(forward.j.as_ref().map(|s| &s.call), reverse.j.as_ref().map(|s| &s.call));
        assert_eq!(forward.v.as_ref().map(|s| s.cigar.clone()), reverse.v.as_ref().map(|s| s.cigar.clone()));
    }
}

#[test]
fn hypermutated_reads_keep_their_v_gene() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let mut rng = Rng(99);
    for _ in 0..40 {
        let d = rng.below(8);
        let sim = simulate(&db, &mut rng, &[d], 0.08);
        let rec = record(&assigner, &sim.query);
        let called = first_label(&rec.v).expect("V must be called at 8% SHM");
        assert_eq!(gene(&called), gene(&sim.v));
        assert_eq!(first_label(&rec.j).as_deref(), Some(sim.j.as_str()));
    }
}

#[test]
fn junctions_without_d_do_not_invent_one() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let mut rng = Rng(3);
    let trials = 40;
    let with_d = (0..trials)
        .filter(|_| {
            let sim = simulate(&db, &mut rng, &[], 0.0);
            let rec = record(&assigner, &sim.query);
            assert_eq!(first_label(&rec.v).as_deref(), Some(sim.v.as_str()));
            rec.d.is_some()
        })
        .count();
    assert!(with_d * 10 <= trials, "{with_d}/{trials} D-free junctions received a D");
}

#[test]
fn tandem_d_is_a_first_class_path() {
    let db = Db::new();
    let tandem = db.assigner(&EngineConfig { grammar: VdjGrammar::TANDEM_D, ..EngineConfig::default() });
    let single = db.assigner(&EngineConfig::default());
    let mut rng = Rng(5);
    let mut both_correct = 0;
    let trials = 30;
    for _ in 0..trials {
        let first = 5 + rng.below(3);
        let second = (first - 5 + 1 + rng.below(2)) % 3 + 5;
        let sim = simulate(&db, &mut rng, &[first, second], 0.0);
        let rec = record(&tandem, &sim.query);
        let calls = (first_label(&rec.d), first_label(&rec.d2));
        if calls == (Some(sim.d[0].clone()), Some(sim.d[1].clone())) {
            both_correct += 1;
        }
        assert!(record(&single, &sim.query).d2.is_none());
        assert!(rec.np3.is_some(), "VDDJ reports three N/P regions");
    }
    assert!(both_correct * 100 >= trials * 90, "tandem D {both_correct}/{trials}");
}

#[test]
fn residual_adjudicator_is_a_drop_in_replacement() {
    let db = Db::new();
    let residual = db.assigner(&EngineConfig::default()).with_adjudicator(ResidualAdjudicator::default());
    let mut rng = Rng(11);
    for _ in 0..20 {
        let d = rng.below(8);
        let sim = simulate(&db, &mut rng, &[d], 0.0);
        let rec = record(&residual, &sim.query);
        assert_eq!(first_label(&rec.v).as_deref(), Some(sim.v.as_str()));
        assert_eq!(first_label(&rec.j).as_deref(), Some(sim.j.as_str()));
    }
}

#[test]
fn called_segments_never_share_query_bases() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig { grammar: VdjGrammar::TANDEM_D, ..EngineConfig::default() });
    let mut rng = Rng(13);
    let mut ws = assigner.workspace();
    for _ in 0..50 {
        let n_d = rng.below(3);
        let ds: Vec<usize> = (0..n_d).map(|_| rng.below(8)).collect();
        let sim = simulate(&db, &mut rng, &ds, 0.03);
        let call = assigner.assign(&encode(sim.query.as_bytes()), &mut ws);
        let path = call.path.expect("a path exists");
        for pair in path.segments.windows(2) {
            assert!(pair[0].span.q_hi <= pair[1].span.q_lo);
            assert!(pair[0].exit <= pair[1].entry);
        }
    }
}

#[test]
fn batch_matches_sequential() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let mut rng = Rng(17);
    let queries: Vec<Vec<u8>> = (0..300)
        .map(|_| {
            let d = rng.below(8);
            let sim = simulate(&db, &mut rng, &[d], 0.02);
            let q = if rng.below(2) == 0 { sim.query } else { reverse_complement(&sim.query) };
            encode(q.as_bytes())
        })
        .collect();
    let parallel = assign_all(&assigner, &queries, 4);
    let mut ws = assigner.workspace();
    for (query, call) in queries.iter().zip(&parallel) {
        let sequential = assigner.assign(query, &mut ws);
        assert_eq!(sequential.path, call.path);
        assert_eq!(sequential.orientation, call.orientation);
    }
}

// ---------------------------------------------------------------------------
// Brute-force check of the DAG's running-maximum shortcuts.
// ---------------------------------------------------------------------------

/// Exhaustive search of the segmentation model: every grammatical path, and
/// for every transition every exit cut `x` and entry cut `y ≥ x` within
/// `MARGIN` bases of the placements' meeting point — wider than the DAG's own
/// search (entries only inside the overlap), so it also checks that the DAG's
/// restriction loses nothing. Memoised on `(placement, entry, D count, first)`.
struct Brute<'a> {
    dag: &'a PlacementDag,
    placements: &'a [Placement],
    memo: std::collections::HashMap<(usize, u32, u8, SegmentKind), Option<PathScore>>,
}

const MARGIN: u32 = 12;

impl Brute<'_> {
    fn min(&self, kind: SegmentKind) -> u32 {
        self.dag.min_retained[kind.index()].max(1)
    }

    /// Best rest of a path whose current segment is `b`, entered at `entry`.
    fn best_from(&mut self, b: usize, entry: u32, d_used: u8, first: SegmentKind) -> Option<PathScore> {
        let key = (b, entry, d_used, first);
        if let Some(known) = self.memo.get(&key) {
            return *known;
        }
        let (grammar, model) = (&self.dag.grammar, &self.dag.model);
        let pb = &self.placements[b];
        let bonus = model.segment_bonus(pb, d_used);
        let mut best: Option<PathScore> = None;
        let mut consider = |candidate: PathScore| {
            if best.map_or(true, |current| candidate > current) {
                best = Some(candidate);
            }
        };
        if grammar.can_end(first, pb.kind) {
            if let Some(span) = pb.clip.span(entry, pb.q_end()).filter(|s| s.q_len() >= self.min(pb.kind)) {
                let score = span.score + bonus + model.end(first, pb.kind, d_used);
                consider(PathScore { score, d_count: d_used, support: pb.score() });
            }
        }
        for e in 0..self.placements.len() {
            let pe = &self.placements[e];
            if pe.q_start() <= pb.q_start() || !grammar.can_follow(pb.kind, pe.kind, d_used) {
                continue;
            }
            let meet = pb.q_end().min(pe.q_start());
            let exits = meet.saturating_sub(MARGIN).max(entry + 1)..=pb.q_end();
            let last_entry = (pb.q_end().max(pe.q_start()) + MARGIN).min(pe.q_end() - 1);
            let d_next = d_used + u8::from(pe.kind == SegmentKind::D);
            for x in exits {
                let Some(span) = pb.clip.span(entry, x).filter(|s| s.q_len() >= self.min(pb.kind)) else {
                    continue;
                };
                for y in pe.q_start().max(x)..=last_entry {
                    let Some(next) = pe.clip.entry(y) else {
                        continue;
                    };
                    let junction = model.junction(pb.kind, pe.kind, next.q.saturating_sub(span.q_hi));
                    let Some(rest) = self.best_from(e, y, d_next, first) else {
                        continue;
                    };
                    consider(PathScore {
                        score: span.score + bonus + junction + rest.score,
                        d_count: rest.d_count,
                        support: pb.score() + rest.support,
                    });
                }
            }
        }
        self.memo.insert(key, best);
        best
    }
}

fn top_by_kind(placements: &[Placement], kind: SegmentKind, n: usize) -> Vec<Placement> {
    let mut chosen: Vec<Placement> = placements.iter().filter(|p| p.kind == kind).cloned().collect();
    chosen.sort_by_key(|p| std::cmp::Reverse(p.score()));
    chosen.truncate(n);
    chosen
}

#[test]
fn dag_equals_brute_force_over_all_paths_and_cuts() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig { grammar: VdjGrammar::TANDEM_D, ..EngineConfig::default() });
    let config = EngineConfig::default();
    let model = EngineeringObjective {
        insertion_per_base: -1,
        ..config.objective.resolve(Some(&Db::set(SegmentKind::D, &db.d, None)), &config.d.scoring)
    };
    let dag = PlacementDag::new(VdjGrammar::TANDEM_D, model, [10, 3, 6]);
    let mut rng = Rng(23);
    let mut ws = assigner.workspace();
    let mut scratch = DagScratch::default();
    let mut distinct_topologies = HashSet::new();
    for _ in 0..40 {
        let ds: Vec<usize> = (0..rng.below(3)).map(|_| rng.below(8)).collect();
        let sim = simulate(&db, &mut rng, &ds, 0.04);
        let call = assigner.assign(&encode(sim.query.as_bytes()), &mut ws);
        let mut placements = top_by_kind(&call.placements, SegmentKind::V, 4);
        placements.extend(top_by_kind(&call.placements, SegmentKind::D, 6));
        placements.extend(top_by_kind(&call.placements, SegmentKind::J, 3));

        let fast = dag.adjudicate(call.query.len(), &placements, &mut scratch);
        let mut brute = Brute { dag: &dag, placements: &placements, memo: Default::default() };
        let exhaustive = (0..placements.len())
            .filter(|&s| dag.grammar.can_start(placements[s].kind))
            .filter_map(|s| {
                let p = &placements[s];
                brute.best_from(s, p.q_start(), u8::from(p.kind == SegmentKind::D), p.kind)
            })
            .max();
        assert_eq!(fast.as_ref().map(|path| path.score), exhaustive);
        if let Some(path) = fast {
            let kinds: Vec<SegmentKind> = path.segments.iter().map(|s| placements[s.placement].kind).collect();
            distinct_topologies.insert(kinds);
            let clip_sum: i32 = path.segments.iter().map(|s| s.span.score).sum();
            assert!(clip_sum >= path.score.score - 24, "path score is built from clip scores");
        }
    }
    assert!(distinct_topologies.len() >= 3, "test should exercise VJ, VDJ and VDDJ paths");
}

#[test]
fn simd_batched_pipeline_equals_scalar_per_read() {
    let db = Db::new();
    let config = EngineConfig { grammar: VdjGrammar::TANDEM_D, ..EngineConfig::default() };
    let simd = db.assigner(&config);
    let scalar = db.assigner(&config).with_aligner(BandedGotoh);
    let mut rng = Rng(29);
    let queries: Vec<Vec<u8>> = (0..400)
        .map(|_| {
            let ds: Vec<usize> = (0..rng.below(3)).map(|_| rng.below(8)).collect();
            let sim = simulate(&db, &mut rng, &ds, 0.05);
            let q = if rng.below(2) == 0 { sim.query } else { reverse_complement(&sim.query) };
            encode(q.as_bytes())
        })
        .collect();
    let batched = assign_all(&simd, &queries, 3);
    let mut ws = scalar.workspace();
    for (query, call) in queries.iter().zip(&batched) {
        let reference = scalar.assign(query, &mut ws);
        assert_eq!(call.orientation, reference.orientation);
        assert_eq!(call.path, reference.path);
        assert_eq!(call.placements.len(), reference.placements.len());
        for (a, b) in call.placements.iter().zip(&reference.placements) {
            assert_eq!(
                (a.kind, &a.alleles, a.score(), a.q_start(), a.q_end()),
                (b.kind, &b.alleles, b.score(), b.q_start(), b.q_end())
            );
        }
    }
}
