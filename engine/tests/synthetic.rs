//! End-to-end tests on a synthetic germline database.
//!
//! The database is random but deterministic (SplitMix64): 6 V genes × 2
//! alleles (SNPs at 50/150/250), 8 D genes of 14–21 nt, 4 J genes of 48 nt,
//! and 2 C isotypes × 2 subclasses of 330 nt sharing their first 24 nt.
//! Queries are simulated rearrangements with trimming and random N regions.

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

/// AIRR names N/P regions by their neighbours: `np1` follows V and `np2`
/// the D, so a read without V reports its D–J bases as `np2`.
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
fn called_segments_never_share_query_bases() {
    let db = Db::new();
    let assigner = db.assigner(&EngineConfig::default());
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
        let record = AirrRecord::from_call("q", &decode(query), call, &simd, &[None, None, None]);
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

/// Max-marginal D decoding keeps the call's topology and D presence on every
/// read, and changes a D only to an allele of larger posterior mass that
/// the report rates at least as probable as the call's. The row names the
/// decoded allele first, also when it shares its placement's alignment with
/// another allele, and reports that allele's own germline: its bases
/// between the reported germline coordinates.
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
    let set = decoded.segment(SegmentKind::D).unwrap().set();
    let (mut changed, mut merged) = (0, 0);
    for read in sim.reads(1200, 9).unwrap() {
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
        let d = |c: &'_ Call, p: &CallPath| {
            p.segments
                .iter()
                .find(|s| c.placements[s.placement].kind == SegmentKind::D)
                .map(|s| (s.allele, c.placements[s.placement].alleles.clone()))
        };
        let (Some((None, called)), Some((chosen, alleles))) = (d(&a, pa), d(&b, pb)) else {
            continue;
        };
        let sequence = String::from_utf8_lossy(&read.sequence);
        let record = AirrRecord::from_call("q", &sequence, &b, &decoded, &[None, None, None]);
        let report = record.d.expect("a D call");
        // The reported first allele's own bases lie between its coordinates.
        let first = report.call.split(',').next().unwrap();
        let allele = set.alleles().iter().find(|a| a.names().any(|n| n == first));
        let codes = &allele.expect("a known name").codes;
        let germline: String = report
            .germline_alignment
            .chars()
            .filter(|&c| c != '-')
            .collect();
        let range = report.germline_start as usize - 1..report.germline_end as usize;
        assert_eq!(germline, decode(&codes[range]), "{first}");
        let Some(chosen) = chosen else {
            continue;
        };
        // The decoded allele is named first, outweighs the call's, and is
        // reported at least as probable.
        assert!(set.allele(chosen).names().any(|n| n == first), "{first}");
        let map_record = AirrRecord::from_call("q", &sequence, &a, &map, &[None, None, None]);
        let posterior = |r: &swig::report::SegmentReport| r.posterior.expect("posteriors");
        assert!(posterior(&report) >= posterior(map_record.d.as_ref().unwrap()));
        let m = b.marginals.as_ref().unwrap();
        let mass = |allele| -> f64 {
            b.placements
                .iter()
                .zip(&m.placement)
                .filter(|(p, _)| p.kind == SegmentKind::D && p.alleles.contains(&allele))
                .map(|(_, &q)| q)
                .sum()
        };
        assert!(called.iter().all(|&c| mass(chosen) > mass(c)));
        changed += 1;
        merged += usize::from(alleles[0] != chosen);
    }
    assert!(changed > 0, "some D calls should change ({changed})");
    assert!(
        merged > 0,
        "some decoded alleles share a placement ({merged})"
    );
}

#[test]
fn simd_batched_pipeline_equals_scalar_per_read() {
    let db = Db::new();
    let config = EngineConfig::default();
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
