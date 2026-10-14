//! Chance D on reads without one: a D that ends its path (no J after it)
//! or starts it (no V before it) pays for the J (V) its read would show
//! beyond it ([`swig::config::OpenFlank`]). Reads of a full V followed by
//! random sequence, or random sequence followed by a J, get no D; reads
//! whose real D sits near their end keep it.

use swig::germline::GermlineSet;
use swig::prelude::*;
use swig::reads::ReadRecord;
use swig::run::{Run, RunOptions};
use swig::sim::SyntheticSpec;
use swig::sim::{Region, Rng, SimParams, SimulatedRead, Simulator, SyntheticDatabase};

fn fasta(set: &GermlineSet) -> String {
    let mut text = Vec::new();
    set.write_fasta(&mut text, |_| None).unwrap();
    String::from_utf8(text).unwrap()
}

/// The heavy-chain database of the tests.
fn heavy() -> SyntheticDatabase {
    SyntheticDatabase::generate(&SyntheticSpec::default(), 7)
}

/// A heavy-chain database and, built from another seed, a light-chain one
/// renamed to IGK without D: a BCR database of two loci.
fn bcr() -> (SyntheticDatabase, String, String, String) {
    let heavy = heavy();
    let light = SyntheticDatabase::generate(&SyntheticSpec::default(), 8);
    let rename = |set: &GermlineSet, from: &str, to: &str| {
        fasta(set).replace(from, to).replace(" locus=IGH", "")
    };
    let v = fasta(&heavy.v) + &rename(&light.v, "IGHV", "IGKV");
    let j = fasta(&heavy.j) + &rename(&light.j, "IGHJ", "IGKJ");
    let d = fasta(&heavy.d);
    (heavy, v, d, j)
}

/// Germline-sense simulated reads with a D, from `seed`.
fn simulated(
    db: &SyntheticDatabase,
    params: SimParams,
    count: usize,
    seed: u64,
) -> Vec<SimulatedRead> {
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params).unwrap();
    sim.reads(count * 4, seed)
        .unwrap()
        .into_iter()
        .filter(|r| !r.reverse_complemented && r.span(Region::D).is_some())
        .take(count)
        .collect()
}

fn record(id: String, sequence: Vec<u8>) -> ReadRecord {
    ReadRecord {
        id,
        sequence,
        quality: String::new(),
    }
}

/// Each read's V (its 5′ flank and V) followed by 100 random bases.
fn v_then_random(reads: &[SimulatedRead], seed: u64) -> Vec<ReadRecord> {
    let mut rng = Rng::new(seed);
    reads
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let v_end = r.span(Region::V).unwrap().end;
            let mut sequence = r.forward[..v_end].to_vec();
            sequence.extend(rng.bases(100));
            record(format!("vcut{i}"), sequence)
        })
        .collect()
}

/// 100 random bases followed by each read's J and what follows it.
fn random_then_j(reads: &[SimulatedRead], seed: u64) -> Vec<ReadRecord> {
    let mut rng = Rng::new(seed);
    reads
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let j_start = r.span(Region::J).unwrap().start;
            let mut sequence = rng.bases(100);
            sequence.extend_from_slice(&r.forward[j_start..]);
            record(format!("jcut{i}"), sequence)
        })
        .collect()
}

/// `d_call` and `d_presence_posterior` of every row.
fn d_calls(
    v: &str,
    d: &str,
    j: &str,
    options: &RunOptions,
    reads: &[ReadRecord],
) -> Vec<(String, f64)> {
    let texts = GermlineTexts {
        v,
        d: Some(d),
        j,
        c: None,
    };
    let database = PreparedDatabase::from_fasta(texts, options.annotate()).unwrap();
    let mut run = Run::new(options, database).unwrap();
    let header = run.airr_header();
    let columns: Vec<&str> = header.split('\t').collect();
    let column = |name: &str| columns.iter().position(|&c| c == name).unwrap();
    let (call, presence) = (column("d_call"), column("d_presence_posterior"));
    let rows = run.assign(reads);
    rows.airr
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let p = fields[presence].parse().unwrap_or(0.0);
            (fields[call].to_string(), p)
        })
        .collect()
}

/// The D calls and the mean D presence of `rows`.
fn summary(rows: &[(String, f64)]) -> (usize, f64) {
    let called = rows.iter().filter(|(d, _)| !d.is_empty()).count();
    let mean = rows.iter().map(|(_, p)| p).sum::<f64>() / rows.len() as f64;
    (called, mean)
}

fn heavy_texts(db: &SyntheticDatabase) -> (String, String, String) {
    (fasta(&db.v), fasta(&db.d), fasta(&db.j))
}

/// (a) A full V followed by 100 random bases: a chance D in the random
/// bases would assert a junction whose J the read does not show. Swig
/// 0.39.3 calls a D on 56 of these 100 reads at mean presence 0.60 (this
/// database has 35 D alleles; 93 % on human IGH); the cost on 5, at 0.14.
/// The bound is the acceptance gate's 10 %.
#[test]
fn a_full_v_followed_by_random_bases_gets_no_d() {
    let db = heavy();
    let (v, d, j) = heavy_texts(&db);
    let reads = v_then_random(&simulated(&db, SimParams::mid(), 100, 11), 101);
    let (called, presence) = summary(&d_calls(&v, &d, &j, &RunOptions::default(), &reads));
    assert!(
        called <= 10 && presence < 0.2,
        "{called} D calls, mean presence {presence:.3}"
    );
}

/// (b) 100 random bases followed by a J (and its C): a chance D first on
/// its path would assert a junction whose V the read does not show. Swig
/// 0.39.3 calls a D on all 100 reads at mean presence 0.99; the cost on 14,
/// at 0.41. The V side prices less than the J side: the read's only anchor
/// is its J, whose few columns leave the V's divergence uncertain, and a V
/// stays hidden in more bases than a J (on human IGH 12 against 8: accept
/// scores 24 and 16, as measured on reads cut inside them). This bound is
/// what the cost reaches, so that it catches any weakening of the V side.
#[test]
fn random_bases_followed_by_a_j_get_few_ds() {
    let db = heavy();
    let (v, d, j) = heavy_texts(&db);
    let reads = random_then_j(&simulated(&db, SimParams::mid(), 100, 12), 102);
    let (called, presence) = summary(&d_calls(&v, &d, &j, &RunOptions::default(), &reads));
    assert!(
        called <= 15 && presence < 0.45,
        "{called} D calls, mean presence {presence:.3}"
    );
}

/// (c) Fixture (a) against a BCR database of two loci: the heavy locus is
/// adjudicated alone and pays the same cost. Swig 0.39.3: 86 D calls at
/// mean presence 0.86; the cost: 4, at 0.18.
#[test]
fn a_full_v_followed_by_random_bases_gets_no_d_in_a_bcr_database() {
    let (db, v, d, j) = bcr();
    let reads = v_then_random(&simulated(&db, SimParams::mid(), 100, 11), 101);
    let (called, presence) = summary(&d_calls(&v, &d, &j, &RunOptions::default(), &reads));
    assert!(
        called <= 10 && presence < 0.2,
        "{called} D calls, mean presence {presence:.3}"
    );
}

/// (d) Fixture (a) under `--calls fixed`: the fixed scale charges the cost
/// as `⌊P*_D + c/λ_D⌋ − ⌊P*_D⌋` score units, its posterior `c` nats. Swig
/// 0.39.3: 55 D calls at mean presence 0.62; the cost: 8, at 0.17.
#[test]
fn a_full_v_followed_by_random_bases_gets_no_d_on_the_fixed_scale() {
    let db = heavy();
    let (v, d, j) = heavy_texts(&db);
    let reads = v_then_random(&simulated(&db, SimParams::mid(), 100, 11), 101);
    let mut options = RunOptions::default();
    options.set("calls", Some("fixed")).unwrap();
    let (called, presence) = summary(&d_calls(&v, &d, &j, &options, &reads));
    assert!(
        called <= 10 && presence < 0.2,
        "{called} D calls, mean presence {presence:.3}"
    );
}

/// Guard: a read that ends 2 bases after its real D (V, N1, D) keeps it;
/// the cost is 0 within the N/P allowance.
#[test]
fn guard_a_read_ending_just_after_its_d_keeps_it() {
    let db = heavy();
    let (v, d, j) = heavy_texts(&db);
    let reads: Vec<ReadRecord> = simulated(&db, SimParams::mid(), 60, 13)
        .iter()
        .filter(|r| r.span(Region::D).unwrap().len() >= 12)
        .enumerate()
        .map(|(i, r)| {
            let end = r.span(Region::D).unwrap().end + 2;
            record(
                format!("dend{i}"),
                r.forward[..end.min(r.forward.len())].to_vec(),
            )
        })
        .collect();
    assert!(reads.len() >= 40, "{}", reads.len());
    // Swig 0.39.3 calls 29 of these 48 Ds: a D at a read's end must beat
    // the V-only path on its own, which the cost leaves as it was.
    let (called, _) = summary(&d_calls(&v, &d, &j, &RunOptions::default(), &reads));
    assert!(called >= 29, "{called} D calls of {}", reads.len());
}

/// Guard: a read that starts 10 bases before its real D keeps it.
#[test]
fn guard_a_read_starting_just_before_its_d_keeps_it() {
    let db = heavy();
    let (v, d, j) = heavy_texts(&db);
    let reads: Vec<ReadRecord> = simulated(&db, SimParams::mid(), 60, 14)
        .iter()
        .filter(|r| r.span(Region::D).unwrap().len() >= 12)
        .enumerate()
        .map(|(i, r)| {
            let start = r.span(Region::D).unwrap().start.saturating_sub(10);
            record(format!("dstart{i}"), r.forward[start..].to_vec())
        })
        .collect();
    assert!(reads.len() >= 40, "{}", reads.len());
    let (called, _) = summary(&d_calls(&v, &d, &j, &RunOptions::default(), &reads));
    assert_eq!(called, reads.len());
}

/// Guard: at 30 % SHM, a read that ends 20–40 bases into its J keeps its
/// D. The J stage misses such a short, mutated J on about a fifth of these
/// reads, whose D then ends its path some 20–60 bases before the read's
/// end and pays for the J it implies there. The miss model reads the
/// divergence of the read's V, at which those bases could well hide a J,
/// so the D pays little. A model that took the read for unmutated drops 5
/// of these D calls. The reference is the same run with
/// `--d-open-flank off`, which reproduces Swig 0.39.3.
#[test]
fn guard_a_mutated_read_ending_in_a_short_j_keeps_its_d() {
    let db = heavy();
    let (v, d, j) = heavy_texts(&db);
    let simulated = simulated(&db, SimParams::mid().with_shm(0.3), 120, 15);
    let reads: Vec<ReadRecord> = simulated
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let end = r.span(Region::J).unwrap().start + 20 + i % 21;
            record(
                format!("shortj{i}"),
                r.forward[..end.min(r.forward.len())].to_vec(),
            )
        })
        .collect();
    let truth: Vec<&str> = simulated
        .iter()
        .map(|r| db.d.alleles()[r.d.unwrap()].name.as_str())
        .collect();
    let calls = |options: &RunOptions| {
        let rows = d_calls(&v, &d, &j, options, &reads);
        let correct = rows
            .iter()
            .zip(&truth)
            .filter(|((call, _), name)| call.split(',').any(|c| c == **name))
            .count();
        (summary(&rows).0, correct)
    };
    let mut off = RunOptions::default();
    off.set("d-open-flank", Some("off")).unwrap();
    // 114 D calls, 51 of them right, as with the cost (52).
    let (called_off, correct_off) = calls(&off);
    let (called, correct) = calls(&RunOptions::default());
    assert!(
        called + 1 >= called_off && correct + 1 >= correct_off,
        "{called} D calls ({correct} right), off: {called_off} ({correct_off})"
    );
}
