//! End-to-end tests on a synthetic germline database.
//!
//! The database is random but deterministic (SplitMix64): 6 V genes × 2
//! alleles (SNPs at 50/150/250), 8 D genes of 14–21 nt, 4 J genes of 48 nt,
//! and 2 C isotypes × 2 subclasses of 330 nt sharing their first 24 nt.
//! Queries are simulated rearrangements with trimming and random N regions.

use std::collections::{HashMap, HashSet};

use swig::adjudicate::{
    Adjudicator, DagScratch, EvidenceModel, PathScore, ReadModel, Topology, Units, Within,
    divergence,
};
use swig::placement::Placement;
use swig::prelude::*;
use swig::sim::Rng;

/// Random ASCII bases as a string.
fn bases(rng: &mut Rng, len: usize) -> String {
    String::from_utf8(rng.bases(len)).expect("ASCII bases")
}

/// A different base than `base`.
fn substitute(base: char, rng: &mut Rng) -> char {
    char::from(rng.substitute(base as u8))
}

struct Db {
    v: Vec<(String, String)>,
    d: Vec<(String, String)>,
    j: Vec<(String, String)>,
    c: Vec<(String, String)>,
}

const V_LEN: usize = 290;
const V_ANCHOR: usize = 281;
const J_ANCHOR: usize = 15;

impl Db {
    fn new() -> Self {
        let mut rng = Rng::new(2024);
        let mut v = Vec::new();
        for g in 1..=6 {
            let first = bases(&mut rng, V_LEN);
            let mut second: Vec<char> = first.chars().collect();
            for p in [50, 150, 250] {
                second[p] = substitute(second[p], &mut rng);
            }
            v.push((format!("IGHV{g}-1*01"), first));
            v.push((format!("IGHV{g}-1*02"), second.into_iter().collect()));
        }
        let d = (1..=8)
            .map(|g| (format!("IGHD{g}-1*01"), bases(&mut rng, 13 + g)))
            .collect();
        let j = (1..=4)
            .map(|g| (format!("IGHJ{g}*01"), bases(&mut rng, 48)))
            .collect();
        let mut c = Vec::new();
        for isotype in 1..=2 {
            let root = bases(&mut rng, 330);
            for subclass in 1..=2 {
                let sequence: String = root
                    .chars()
                    .enumerate()
                    .map(|(k, b)| {
                        if k >= 24 && rng.below(10) == 0 {
                            substitute(b, &mut rng)
                        } else {
                            b
                        }
                    })
                    .collect();
                c.push((format!("IGHG{isotype}{subclass}*01"), sequence));
            }
        }
        Self { v, d, j, c }
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
    let v_part: String = v_seq[10..V_LEN - rng.range(0, 3) as usize]
        .chars()
        .map(|b| {
            if (rng.below(10_000) as f64) < shm_rate * 10_000.0 {
                substitute(b, rng)
            } else {
                b
            }
        })
        .collect();
    let mut query = v_part;
    let mut d_names = Vec::new();
    for &d in d_choice {
        let (d_name, d_seq) = &db.d[d];
        let (left, right) = (rng.range(0, 2) as usize, rng.range(0, 2) as usize);
        let n_len = rng.range(3, 6) as usize;
        query += &bases(rng, n_len);
        query += &d_seq[left..d_seq.len() - right];
        d_names.push(d_name.clone());
    }
    let n_len = rng.range(3, 6) as usize;
    query += &bases(rng, n_len);
    query += &j_seq[rng.range(0, 3) as usize..];
    Sim {
        query,
        v: v_name.clone(),
        d: d_names,
        j: j_name.clone(),
    }
}

/// `query` with `c[..len]` spliced after its J and `tail` random bases.
fn with_constant(query: &str, c: &str, len: usize, tail: usize, rng: &mut Rng) -> String {
    format!("{query}{}{}", &c[..len], bases(rng, tail))
}

/// The assigner with the database's C stage.
fn with_c(db: &Db, config: &EngineConfig) -> DefaultAssigner {
    db.assigner(config)
        .with_constant(Db::set(SegmentKind::C, &db.c, None), &config.c)
        .unwrap()
}

fn first_label(report: &Option<swig::report::SegmentReport>) -> Option<String> {
    report
        .as_ref()
        .map(|r| r.call.split(',').next().unwrap().to_string())
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
    let mut rng = Rng::new(1);
    let (mut d_correct, trials) = (0, 60);
    for _ in 0..trials {
        let d = rng.below(db.d.len());
        let sim = simulate(&db, &mut rng, &[d], 0.0);
        let rec = record(&assigner, &sim.query);
        assert!(!rec.rev_comp);
        assert_eq!(
            first_label(&rec.v).as_deref(),
            Some(sim.v.as_str()),
            "V allele for {}",
            sim.query
        );
        assert_eq!(
            first_label(&rec.j).as_deref(),
            Some(sim.j.as_str()),
            "J allele"
        );
        if first_label(&rec.d).as_deref() == Some(sim.d[0].as_str()) {
            d_correct += 1;
        }
    }
    assert!(
        d_correct * 100 >= trials * 95,
        "D accuracy {d_correct}/{trials}"
    );
}

#[test]
fn reverse_complement_reads_are_reoriented() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let mut rng = Rng::new(7);
    for _ in 0..20 {
        let d = rng.below(8);
        let sim = simulate(&db, &mut rng, &[d], 0.0);
        let forward = record(&assigner, &sim.query);
        let reverse = record(&assigner, &reverse_complement(&sim.query));
        assert!(reverse.rev_comp);
        assert_eq!(
            forward.v.as_ref().map(|s| &s.call),
            reverse.v.as_ref().map(|s| &s.call)
        );
        assert_eq!(
            forward.j.as_ref().map(|s| &s.call),
            reverse.j.as_ref().map(|s| &s.call)
        );
        assert_eq!(
            forward.v.as_ref().map(|s| s.cigar.clone()),
            reverse.v.as_ref().map(|s| s.cigar.clone())
        );
    }
}

#[test]
fn hypermutated_reads_keep_their_v_gene() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let mut rng = Rng::new(99);
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
    let mut rng = Rng::new(3);
    let trials = 40;
    let with_d = (0..trials)
        .filter(|_| {
            let sim = simulate(&db, &mut rng, &[], 0.0);
            let rec = record(&assigner, &sim.query);
            assert_eq!(first_label(&rec.v).as_deref(), Some(sim.v.as_str()));
            rec.d.is_some()
        })
        .count();
    assert!(
        with_d * 10 <= trials,
        "{with_d}/{trials} D-free junctions received a D"
    );
}

/// AIRR names N/P regions by their neighbours: `np1` follows V, `np2` the
/// first D and `np3` the second, so a read without V reports its D–J
/// bases as `np2`.
#[test]
fn np_regions_are_named_by_their_neighbours() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
    let (d_name, d) = &db.d[7];
    let (j_name, j) = &db.j[0];
    let rec = record(&assigner, &format!("{d}GTA{j}"));
    assert!(rec.v.is_none());
    assert_eq!(first_label(&rec.d).as_deref(), Some(d_name.as_str()));
    assert_eq!(first_label(&rec.j).as_deref(), Some(j_name.as_str()));
    assert_eq!(
        (rec.np1.as_deref(), rec.np2.as_deref()),
        (None, Some("GTA"))
    );
}

#[test]
fn tandem_d_is_a_first_class_path() {
    let db = Db::new();
    let tandem = db.assigner(&EngineConfig {
        grammar: VdjGrammar::TANDEM_D,
        ..EngineConfig::default()
    });
    let single = db.assigner(&EngineConfig::default());
    let mut rng = Rng::new(5);
    let mut both_correct = 0;
    let trials = 30;
    for _ in 0..trials {
        let first = 5 + rng.below(3);
        let second = (first - 5 + 1 + rng.below(2)) % 3 + 5;
        let sim = simulate(&db, &mut rng, &[first, second], 0.0);
        let d_labels = |assigner: &DefaultAssigner| {
            let call = assigner.assign(&encode(sim.query.as_bytes()), &mut assigner.workspace());
            let set = assigner.segment(SegmentKind::D).unwrap().set();
            call.path
                .iter()
                .flat_map(|p| &p.segments)
                .map(|s| &call.placements[s.placement])
                .filter(|p| p.kind == SegmentKind::D)
                .map(|p| set.allele(p.primary_allele()).name.clone())
                .collect::<Vec<_>>()
        };
        if d_labels(&tandem) == sim.d[..2] {
            both_correct += 1;
        }
        assert!(
            d_labels(&single).len() <= 1,
            "the default grammar allows one D"
        );
    }
    assert!(
        both_correct * 100 >= trials * 90,
        "tandem D {both_correct}/{trials}"
    );
}

#[test]
fn residual_adjudicator_is_a_drop_in_replacement() {
    let db = Db::new();
    let config = EngineConfig {
        posterior: false,
        ..EngineConfig::default()
    };
    let residual = db
        .assigner(&config)
        .with_adjudicator(ResidualAdjudicator::default());
    let mut rng = Rng::new(11);
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
    let assigner = db.assigner(&EngineConfig {
        grammar: VdjGrammar::TANDEM_D,
        ..EngineConfig::default()
    });
    let mut rng = Rng::new(13);
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
    let mut rng = Rng::new(17);
    let queries: Vec<Vec<u8>> = (0..300)
        .map(|_| {
            let d = rng.below(8);
            let sim = simulate(&db, &mut rng, &[d], 0.02);
            let q = if rng.below(2) == 0 {
                sim.query
            } else {
                reverse_complement(&sim.query)
            };
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

/// The C stage reads finished calls only: placements, path and marginals are
/// those of the assigner without it, and a read without C gets none.
#[test]
fn the_constant_stage_leaves_vdj_calls_unchanged() {
    let db = Db::new();
    let config = EngineConfig::default();
    let (plain, constant) = (db.assigner(&config), with_c(&db, &config));
    let mut rng = Rng::new(31);
    let queries: Vec<Vec<u8>> = (0..300)
        .map(|_| {
            let d = rng.below(8);
            let sim = simulate(&db, &mut rng, &[d], 0.05);
            let (c, len) = (&db.c[rng.below(4)].1, rng.below(80));
            let q = with_constant(&sim.query, c, len, rng.below(20), &mut rng);
            let q = if rng.below(2) == 0 {
                q
            } else {
                reverse_complement(&q)
            };
            encode(q.as_bytes())
        })
        .collect();
    let without = assign_all(&plain, &queries, 2);
    let with = assign_all(&constant, &queries, 3);
    let mut called = 0;
    for (a, b) in without.iter().zip(&with) {
        assert_eq!(
            (a.orientation, &a.path, &a.marginals),
            (b.orientation, &b.path, &b.marginals)
        );
        assert_eq!(a.placements.len(), b.placements.len());
        for (p, q) in a.placements.iter().zip(&b.placements) {
            assert_eq!(
                (p.kind, &p.alleles, &p.alignment),
                (q.kind, &q.alleles, &q.alignment)
            );
        }
        assert!(a.constant.is_empty() && a.c.is_none());
        assert!(b.placements.iter().all(|p| p.kind != SegmentKind::C));
        assert!(b.constant.iter().all(|p| p.kind == SegmentKind::C));
        called += usize::from(b.c.is_some());
    }
    assert!(called > 100, "{called} C calls");
}

/// A spliced C of 30 or more bases is called, as an allele equal to the
/// read's on the spliced bases; random tails are not called; one read alone,
/// a batch and the scalar kernel agree.
#[test]
fn the_constant_region_is_called_from_the_j_exit() {
    let db = Db::new();
    let config = EngineConfig::default();
    let simd = with_c(&db, &config);
    let scalar = with_c(&db, &config).with_aligner(BandedGotoh);
    let mut rng = Rng::new(37);
    let mut reads = Vec::new();
    for k in 0..400 {
        let d = rng.below(8);
        let sim = simulate(&db, &mut rng, &[d], 0.03);
        let c = &db.c[rng.below(4)].1;
        let len = if k % 4 == 0 { 0 } else { 30 + rng.below(150) };
        let tail = if len == 0 { 40 } else { rng.below(10) };
        let query = with_constant(&sim.query, c, len, tail, &mut rng);
        let equivalent: Vec<&str> =
            db.c.iter()
                .filter(|(_, other)| other[..len] == c[..len])
                .map(|(name, _)| name.as_str())
                .collect();
        reads.push((encode(query.as_bytes()), (len > 0).then_some(equivalent)));
    }
    let queries: Vec<&[u8]> = reads.iter().map(|(q, _)| q.as_slice()).collect();
    let batched = assign_all(&simd, &queries, 3);
    let mut ws = scalar.workspace();
    let mut false_calls = 0;
    for ((query, truth), call) in reads.iter().zip(&batched) {
        let alone = scalar.assign(query, &mut ws);
        assert_eq!(alone.c, call.c);
        assert_eq!(alone.constant.len(), call.constant.len());
        for (a, b) in alone.constant.iter().zip(&call.constant) {
            assert_eq!((&a.alleles, &a.alignment), (&b.alleles, &b.alignment));
        }
        let record = AirrRecord::from_call("q", "", call, &simd, &[None, None, None]);
        let primary = call.c.map(|c| {
            let allele = call.constant[c.segment.placement].primary_allele();
            simd.segment(SegmentKind::C)
                .unwrap()
                .set()
                .allele(allele)
                .name
                .clone()
        });
        match truth {
            Some(equivalent) => {
                assert!(equivalent.contains(&primary.unwrap().as_str()));
                let c = call.c.unwrap();
                assert!(c.evalue <= 1e-3 && (0.0..=1.0).contains(&c.posterior));
                assert!(c.segment.span.q_lo >= record.j.as_ref().unwrap().sequence_end);
            }
            None => false_calls += usize::from(primary.is_some()),
        }
    }
    assert!(false_calls <= 1, "{false_calls} false C calls");
}

/// The `c_` columns are empty without a C stage; with one, every other
/// column is unchanged.
#[test]
fn c_columns_are_filled_only_with_a_c_stage() {
    let db = Db::new();
    let config = EngineConfig::default();
    let (plain, constant) = (db.assigner(&config), with_c(&db, &config));
    let header = AirrRecord::tsv_header();
    let constant_column: Vec<bool> = header.split('\t').map(|c| c.starts_with("c_")).collect();
    assert_eq!(constant_column.iter().filter(|&&c| c).count(), 14);
    let mut rng = Rng::new(41);
    let (mut ws_plain, mut ws_c) = (plain.workspace(), constant.workspace());
    for k in 0..60 {
        let d = rng.below(8);
        let sim = simulate(&db, &mut rng, &[d], 0.02);
        let (name, c) = &db.c[rng.below(4)];
        let len = if k % 3 == 0 { 0 } else { 60 + rng.below(100) };
        let query = with_constant(&sim.query, c, len, 5, &mut rng);
        let encoded = encode(query.as_bytes());
        let row = |assigner: &DefaultAssigner, ws| {
            let call = assigner.assign(&encoded, ws);
            AirrRecord::from_call("q", &query, &call, assigner, &[None, None, None])
        };
        let (without, with) = (row(&plain, &mut ws_plain), row(&constant, &mut ws_c));
        assert!(without.c.is_none());
        let (a, b) = (without.to_tsv(), with.to_tsv());
        for ((x, y), &is_c) in a.split('\t').zip(b.split('\t')).zip(&constant_column) {
            if is_c {
                assert!(x.is_empty());
            } else {
                assert_eq!(x, y);
            }
        }
        assert_eq!(b.split('\t').count(), constant_column.len());
        match &with.c {
            Some(report) if len > 0 => {
                assert!(report.call.split(',').any(|n| n == name));
                let j_end = with.j.as_ref().unwrap().sequence_end;
                assert_eq!(report.sequence_start, j_end + 1);
                assert_eq!(report.germline_start, 1);
                assert!(report.sequence_end as usize >= j_end as usize + len);
                assert!(report.support.unwrap() <= 1e-3);
                assert!(report.prefix_identity.is_some_and(|p| p > 0.9));
            }
            None => assert_eq!(len, 0),
            Some(report) => panic!("a C call without C: {report:?}"),
        }
    }
}

/// The V(D)J report refuses a C placement on the path.
#[test]
#[should_panic(expected = "V(D)J")]
fn a_c_placement_never_reaches_the_vdj_report() {
    let db = Db::new();
    let config = EngineConfig::default();
    let assigner = with_c(&db, &config);
    let mut rng = Rng::new(43);
    let sim = simulate(&db, &mut rng, &[0], 0.0);
    let query = with_constant(&sim.query, &db.c[0].1, 100, 0, &mut rng);
    let mut call = assigner.assign(&encode(query.as_bytes()), &mut assigner.workspace());
    let c = call.c.expect("a C call");
    call.placements
        .push(call.constant[c.segment.placement].clone());
    let path = call.path.as_mut().unwrap();
    path.segments.last_mut().unwrap().placement = call.placements.len() - 1;
    AirrRecord::from_call("q", &query, &call, &assigner, &[None, None, None]);
}

// ---------------------------------------------------------------------------
// Brute-force check of the DAG's running-maximum shortcuts.
// ---------------------------------------------------------------------------

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

    /// The cut-independent terms of `p` as a segment with `d_used` D so far.
    fn bonus(&self, p: &Placement, d_used: u8) -> Units {
        let presence = if p.kind == SegmentKind::D {
            self.read.d_cost()
        } else {
            0
        };
        self.read.prior(self.dag.model.segment_bonus(p, d_used)) - presence
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
        let bonus = self.bonus(pb, d_used);
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
            value += self.clip(p, segment.entry, segment.exit) + self.bonus(p, d_used);
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
    ReadModel::fixed(1.0, dag.model.d_presence)
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
    use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let base = config.objective.resolve(Some(&db.d), &config.d);
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
            let dag = PlacementDag::new(VdjGrammar::VDJ, model, config.min_retained);
            let (fast, exhaustive) = dag_and_brute(&dag, &placements, fixed(&dag), &mut scratch);
            assert_eq!(
                fast,
                exhaustive,
                "read {} at iota {iota}",
                read.sequence.len()
            );
            // Count reads whose optimum enters some placement beyond the overlap.
            if let Some(path) = dag.adjudicate(0, &placements, &mut scratch) {
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
    use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let base = config.objective.resolve(Some(&db.d), &config.d);
    let lambda = KarlinAltschul::ungapped(&config.v.scoring).unwrap().lambda;
    let substitution = |s: &AffineScoring| (s.matched, s.mismatch);
    let evidence = EvidenceModel {
        lambda,
        substitution: [&config.v.scoring, &config.d.scoring, &config.j.scoring].map(substitution),
    };
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
                let dag = PlacementDag::new(VdjGrammar::TANDEM_D, objective, config.min_retained);
                let Some((fixed_key, fixed_path)) =
                    dag.viterbi(&placements, &fixed(&dag), &mut scratch)
                else {
                    continue;
                };
                let own = divergence(&placements);
                for p in [own, 1e-3, 0.02, 0.3] {
                    let model = ReadModel::at(&evidence, p, objective.d_presence);
                    let (fast, exhaustive) = dag_and_brute(&dag, &placements, model, &mut scratch);
                    assert_eq!(fast, exhaustive, "p {p}, iota {iota}");
                    let (_, path) = dag.viterbi(&placements, &model, &mut scratch).unwrap();
                    differs += usize::from(path != fixed_path);
                    checked += 1;
                }
                // At δ0 the per-read model is the fixed one times u0.
                let integer = PlacementDag::new(
                    VdjGrammar::TANDEM_D,
                    EngineeringObjective {
                        d_presence: objective.d_presence.floor(),
                        ..objective
                    },
                    config.min_retained,
                );
                let model = ReadModel::at(&evidence, delta0, integer.model.d_presence);
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
    use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let config = EngineConfig::default();
    let assigner =
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap();
    let objective = config.objective.resolve(Some(&db.d), &config.d);
    let lambda = KarlinAltschul::ungapped(&config.v.scoring).unwrap().lambda;
    let evidence = EvidenceModel {
        lambda,
        substitution: [(2, -3); 3],
    };
    let dag = PlacementDag::new(VdjGrammar::VDJ, objective, config.min_retained);
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::hard()).unwrap();
    let mut ws = assigner.workspace();
    let mut scratch = DagScratch::default();
    let mut checked = 0;
    for read in sim.reads(40, 7).unwrap() {
        let call = assigner.assign(&encode(&read.sequence), &mut ws);
        let mut placements = top_by_kind(&call.placements, SegmentKind::V, 3);
        placements.extend(top_by_kind(&call.placements, SegmentKind::D, 6));
        placements.extend(top_by_kind(&call.placements, SegmentKind::J, 3));
        let model = ReadModel::at(&evidence, divergence(&placements), objective.d_presence);
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

/// Max-marginal D decoding keeps the call's topology and D presence on every
/// read, and changes a D only to an allele of larger posterior mass.
#[test]
fn max_marginal_d_keeps_the_topology_and_raises_the_d_mass() {
    use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 11);
    let assigner = |d_allele| {
        let config = EngineConfig {
            d_allele,
            ..EngineConfig::default()
        };
        DefaultAssigner::from_config(&config, db.v.clone(), Some(db.d.clone()), db.j.clone())
            .unwrap()
    };
    let (map, decoded) = (
        assigner(DAlleleDecoding::Map),
        assigner(DAlleleDecoding::MaxMarginal),
    );
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::hard()).unwrap();
    let (mut ws, mut ws2) = (map.workspace(), decoded.workspace());
    let mut changed = 0;
    for read in sim.reads(600, 9).unwrap() {
        let query = encode(&read.sequence);
        let (a, b) = (
            map.assign(&query, &mut ws),
            decoded.assign(&query, &mut ws2),
        );
        assert_eq!(a.marginals, b.marginals, "the posterior does not move");
        let (Some(pa), Some(pb)) = (&a.path, &b.path) else {
            assert_eq!(a.path.is_some(), b.path.is_some());
            continue;
        };
        let kinds = |c: &Call, p: &CallPath| -> Vec<SegmentKind> {
            p.segments
                .iter()
                .map(|s| c.placements[s.placement].kind)
                .collect()
        };
        assert_eq!(kinds(&a, pa), kinds(&b, pb));
        let d = |c: &Call, p: &CallPath| {
            p.segments
                .iter()
                .find(|s| c.placements[s.placement].kind == SegmentKind::D)
                .map(|s| c.placements[s.placement].primary_allele())
        };
        let (da, db_) = (d(&a, pa), d(&b, pb));
        if da != db_ {
            let m = b.marginals.as_ref().unwrap();
            let mass = |allele| -> f64 {
                b.placements
                    .iter()
                    .zip(&m.placement)
                    .filter(|(p, _)| p.kind == SegmentKind::D && p.alleles.contains(&allele))
                    .map(|(_, &q)| q)
                    .sum()
            };
            assert!(mass(db_.unwrap()) > mass(da.unwrap()));
            changed += 1;
        }
    }
    assert!(changed > 0, "some D calls should change ({changed})");
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

#[test]
fn dag_equals_brute_force_over_all_paths_and_cuts() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig {
        grammar: VdjGrammar::TANDEM_D,
        ..EngineConfig::default()
    });
    let config = EngineConfig::default();
    let model = EngineeringObjective {
        insertion_per_base: -1,
        ..config
            .objective
            .resolve(Some(&Db::set(SegmentKind::D, &db.d, None)), &config.d)
    };
    let dag = PlacementDag::new(VdjGrammar::TANDEM_D, model, [10, 3, 6]);
    let mut rng = Rng::new(23);
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

        let (fast, exhaustive) = dag_and_brute(&dag, &placements, fixed(&dag), &mut scratch);
        assert_eq!(fast, exhaustive);
        if let Some(path) = dag.adjudicate(call.query.len(), &placements, &mut scratch) {
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

#[test]
fn simd_batched_pipeline_equals_scalar_per_read() {
    let db = Db::new();
    let config = EngineConfig {
        grammar: VdjGrammar::TANDEM_D,
        ..EngineConfig::default()
    };
    let simd = db.assigner(&config);
    let scalar = db.assigner(&config).with_aligner(BandedGotoh);
    let mut rng = Rng::new(29);
    let queries: Vec<Vec<u8>> = (0..400)
        .map(|_| {
            let ds: Vec<usize> = (0..rng.below(3)).map(|_| rng.below(8)).collect();
            let sim = simulate(&db, &mut rng, &ds, 0.05);
            let q = if rng.below(2) == 0 {
                sim.query
            } else {
                reverse_complement(&sim.query)
            };
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
