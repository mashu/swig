//! Accuracy regression tests on simulated repertoires.
//!
//! Reads come from the IgSim model ([`swig::sim`]) on a synthetic
//! IGH-like germline database with a family → gene → allele hierarchy, half
//! of them reverse-complemented. Seeds are fixed, so every run sees the same
//! reads and reproduces the measured accuracy exactly. Each floor sits three
//! standard errors below the measured rate ([`at_least`]): a change that
//! moves a few reads passes, one that loses accuracy fails. Accuracy is
//! germline-equivalence aware (see `swig::sim::evaluate`).

use swig::prelude::*;
use swig::sim::{
    Assignment, Evaluation, Rng, SimParams, Simulator, SyntheticDatabase, SyntheticSpec, Truth,
};

struct Setup {
    db: SyntheticDatabase,
}

impl Setup {
    fn new() -> Self {
        Self {
            db: SyntheticDatabase::generate(&SyntheticSpec::default(), 11),
        }
    }

    fn assigner(&self, config: &EngineConfig) -> DefaultAssigner {
        self.db.assigner(config).expect("synthetic database builds")
    }

    fn evaluate(
        &self,
        params: SimParams,
        reads: usize,
        seed: u64,
        config: &EngineConfig,
    ) -> Evaluation {
        self.evaluate_with(params, reads, seed, config, false)
    }

    /// With `constant`, reads reach into the synthetic C regions and the
    /// assigner calls them.
    fn evaluate_with(
        &self,
        params: SimParams,
        reads: usize,
        seed: u64,
        config: &EngineConfig,
        constant: bool,
    ) -> Evaluation {
        let mut sim = Simulator::new(
            &self.db.v,
            Some(&self.db.d),
            &self.db.j,
            params.with_reverse_complement(0.5),
        )
        .expect("simulator");
        let mut assigner = self.assigner(config);
        if constant {
            sim = sim.with_constant(&self.db.c).expect("C simulator");
            assigner = assigner
                .with_constant(self.db.c.clone(), &config.c)
                .expect("C stage");
        }
        let simulated = sim.reads(reads, seed).expect("reads");
        let queries: Vec<Vec<u8>> = simulated.iter().map(|r| encode(&r.sequence)).collect();
        let calls = assign_all(&assigner, &queries, 0);
        let evaluation = Evaluation::of_calls(&sim, &simulated, &calls, &assigner);
        println!("{evaluation}");
        evaluation
    }
}

fn percent(part: usize, whole: usize) -> f64 {
    100.0 * part as f64 / whole as f64
}

/// Assert that `part` of `whole` reads is no more than three standard
/// errors below the `measured` percentage: the standard error of a rate
/// near it on `whole` reads (Agresti–Coull, so that a measured 100 % keeps
/// a margin too).
fn at_least(label: &str, part: usize, whole: usize, measured: f64) {
    let n = whole as f64;
    let p = (measured / 100.0 * n + 2.0) / (n + 4.0);
    let floor = measured - 300.0 * (p * (1.0 - p) / (n + 4.0)).sqrt();
    let rate = percent(part, whole);
    println!("{label}: {rate:.2} % of {whole} (floor {floor:.2} %)");
    assert!(
        rate >= floor,
        "{label}: {rate:.2} % ({part} of {whole}) is below {floor:.2} % (measured {measured} %)"
    );
}

/// [`at_least`] for the V, V gene, D, J and junction rates of `e`, measured
/// at `measured` (in that order).
fn segments_at_least(e: &Evaluation, measured: [f64; 5]) {
    let [v, v_gene, d, j, junction] = measured;
    at_least("V", e.v.equivalent, e.v.reads, v);
    at_least("V gene", e.v.gene, e.v.reads, v_gene);
    at_least("D", e.d.equivalent, e.d.reads, d);
    at_least("J", e.j.equivalent, e.j.reads, j);
    at_least("junction", e.junction_exact, e.junction_reads, junction);
}

#[test]
fn unmutated_reads_are_assigned_exactly() {
    let e = Setup::new().evaluate(SimParams::easy(), 300, 1, &EngineConfig::default());
    assert_eq!(e.oriented, e.reads, "orientation");
    assert_eq!(e.unassigned, 0);
    segments_at_least(&e, [100.00, 100.00, 99.66, 100.00, 100.00]);
}

#[test]
fn igm_like_reads_keep_their_accuracy() {
    let e = Setup::new().evaluate(SimParams::mid(), 500, 2, &EngineConfig::default());
    at_least("orientation", e.oriented, e.reads, 100.00);
    segments_at_least(&e, [99.20, 100.00, 83.18, 99.60, 100.00]);
}

#[test]
fn igg_like_hypermutated_reads_keep_their_accuracy() {
    let e = Setup::new().evaluate(SimParams::hard(), 500, 3, &EngineConfig::default());
    at_least("orientation", e.oriented, e.reads, 100.00);
    segments_at_least(&e, [97.60, 100.00, 76.04, 97.20, 100.00]);
}

#[test]
fn long_d_remnants_are_found() {
    let e = Setup::new().evaluate(SimParams::mid(), 500, 4, &EngineConfig::default());
    // Bins 12–15 and 16+ (see `D_LENGTH_BINS`): a D this long is identifiable.
    let long: Vec<_> = e.d_presence[4..].iter().collect();
    let reads: usize = long.iter().map(|p| p.reads).sum();
    let correct: usize = long.iter().map(|p| p.correct).sum();
    at_least("long D", correct, reads, 98.53);
}

/// Reported V and J posteriors match observed accuracy (expected calibration
/// error and Brier score), on IgM-like and on IgG-like hypermutated reads.
#[test]
fn posteriors_are_calibrated() {
    let setup = Setup::new();
    for (params, seed) in [(SimParams::mid(), 11), (SimParams::hard(), 12)] {
        let e = setup.evaluate(params, 500, seed, &EngineConfig::default());
        for (label, c) in [("V", e.calibration.v), ("J", e.calibration.j)] {
            eprintln!(
                "{label}: n {} ECE {:.3} % Brier {:.4}",
                c.count(),
                100.0 * c.ece(),
                c.brier()
            );
            assert!(c.count() >= 495, "{label}: {} posteriors", c.count());
            assert!(c.ece() <= 0.03, "{label} ECE {:.3} %", 100.0 * c.ece());
            assert!(c.brier() <= 0.04, "{label} Brier {:.4}", c.brier());
        }
        let d = e.calibration.d_presence;
        assert!(d.ece() <= 0.08, "D presence ECE {:.3} %", 100.0 * d.ece());
    }
}

/// The constant region: no call on reads without C, the primary name
/// anchored-equivalent once the C part passes the gate's reach, few names
/// outside the equivalent set, and a calibrated `c_posterior`.
#[test]
fn constant_regions_are_called_from_the_splice() {
    let setup = Setup::new();
    let config = EngineConfig::default();
    // (preset, seed, measured 8-15 nt C, 16+ nt C, 16+ nt first gene)
    let presets = [
        (SimParams::mid(), 21, [78.26, 99.86, 95.88]),
        (SimParams::hard(), 22, [91.67, 100.00, 96.03]),
    ];
    for (params, seed, [medium_c, longer_c, longer_gene]) in presets {
        let e = setup.evaluate_with(params, 1000, seed, &config, true);
        let [none, short, medium, longer @ ..] = &e.c;
        assert!(none.reads >= 150 && none.called <= 1, "{none:?}");
        assert!(short.called * 10 <= short.reads, "1-7 nt {short:?}");
        let label = |what: &str| format!("seed {seed} {what}");
        at_least(
            &label("C 8-15 nt"),
            medium.equivalent,
            medium.reads,
            medium_c,
        );
        let reads: usize = longer.iter().map(|c| c.reads).sum();
        let equivalent: usize = longer.iter().map(|c| c.equivalent).sum();
        let first_gene: usize = longer.iter().map(|c| c.first_gene).sum();
        let calls: usize = longer.iter().map(|c| c.called).sum();
        let over: usize = longer.iter().map(|c| c.over_call).sum();
        at_least(&label("C 16+ nt"), equivalent, reads, longer_c);
        at_least(
            &label("C 16+ nt first gene"),
            first_gene,
            reads,
            longer_gene,
        );
        assert!(percent(over, calls) <= 1.0, "{over}/{calls}");
        let c = e.calibration.c;
        assert!(c.ece() <= 0.02, "C ECE {:.3} %", 100.0 * c.ece());
    }
}

/// Sibling bounds (`family_delta`, V and J) are lossless under every stage
/// configuration they interact with: banding, the ungapped prefilter, `top`,
/// `cap`, the band width and V acceptance. Calls with and without them are
/// identical, down to every placement and posterior.
#[test]
fn family_bounds_agree_with_full_scoring_under_config_variants() {
    use swig::retrieve::UngappedParams;
    use swig::segment::HitAcceptance;

    let setup = Setup::new();
    let sim = Simulator::new(
        &setup.db.v,
        Some(&setup.db.d),
        &setup.db.j,
        SimParams::hard().with_reverse_complement(0.5),
    )
    .expect("simulator");
    let queries: Vec<Vec<u8>> = sim
        .reads(300, 21)
        .expect("reads")
        .iter()
        .map(|read| encode(&read.sequence))
        .collect();
    type Variant = (&'static str, fn(&mut EngineConfig));
    let variants: [Variant; 12] = [
        ("default", |_| {}),
        ("unbanded", |c| c.v.seeds.unbanded = true),
        ("ungapped", |c| {
            c.v.ungapped = Some(UngappedParams::default())
        }),
        ("V top 1", |c| c.v.keep.top = 1),
        ("V top 2", |c| c.v.keep.top = 2),
        ("V top 5", |c| c.v.keep.top = 5),
        ("V cap 3", |c| c.v.keep.cap = 3),
        ("band_max 3", |c| {
            c.v.seeds.band_max = 3;
            c.v.seeds.band_pad = 0;
        }),
        ("J top 1", |c| c.j.keep.top = 1),
        ("J top 5", |c| c.j.keep.top = 5),
        ("J cap 3", |c| c.j.keep.cap = 3),
        ("V acceptance", |c| {
            c.v.accept = HitAcceptance {
                min_subst: 60,
                min_identity: 0.93,
            }
        }),
    ];
    for (name, apply) in variants {
        let calls = |family_delta: bool| -> Vec<String> {
            let mut config = EngineConfig::default();
            apply(&mut config);
            config.v.family_delta = family_delta;
            config.j.family_delta = family_delta;
            assign_all(&setup.assigner(&config), &queries, 0)
                .iter()
                .map(|call| format!("{call:?}"))
                .collect()
        };
        let (bounded, full) = (calls(true), calls(false));
        let differ = bounded.iter().zip(&full).filter(|(a, b)| a != b).count();
        assert_eq!(differ, 0, "{name}: {differ} of {} calls differ", full.len());
    }
}

#[test]
fn random_reads_receive_no_call() {
    let setup = Setup::new();
    let assigner = setup.assigner(&EngineConfig::default());
    let mut rng = Rng::new(5);
    let queries: Vec<Vec<u8>> = (0..300)
        .map(|_| {
            let len = 300 + rng.below(300);
            (0..len).map(|_| encode(&[rng.base()])[0]).collect()
        })
        .collect();
    let with_path = assign_all(&assigner, &queries, 0)
        .iter()
        .filter(|call| call.path.is_some())
        .count();
    assert!(with_path <= 1, "{with_path}/300 random reads were assigned");
}

#[test]
fn airr_tsv_rows_score_like_the_records_they_were_written_from() {
    let setup = Setup::new();
    let sim = Simulator::new(
        &setup.db.v,
        Some(&setup.db.d),
        &setup.db.j,
        SimParams::mid().with_reverse_complement(0.5),
    )
    .expect("simulator");
    let assigner = setup.assigner(&EngineConfig::default());
    let header = AirrRecord::tsv_header();
    let columns: Vec<&str> = header.split('\t').collect();
    let mut ws = assigner.workspace();
    for (i, read) in sim.reads(100, 6).expect("reads").iter().enumerate() {
        let call = assigner.assign(&encode(&read.sequence), &mut ws);
        let sequence = String::from_utf8_lossy(&read.sequence);
        let id = format!("r{i}");
        let record = AirrRecord::from_call(&id, &sequence, &call, &assigner, &[None, None, None]);
        let (parsed_id, mut parsed) =
            Assignment::from_airr(&columns, &record.to_tsv()).expect("AIRR row");
        assert_eq!(parsed_id, id);
        let mut direct = Assignment::from_record(&record);
        assert_eq!(
            record.d_presence_posterior.is_some(),
            call.path.is_some(),
            "posteriors by default"
        );
        // Posteriors are written with six decimals (more near 0 and 1).
        let probabilities = |a: &mut Assignment| {
            let mut out = vec![a.d_presence.take()];
            out.extend(a.segments.iter_mut().flatten().map(|s| s.posterior.take()));
            out
        };
        for (p, q) in probabilities(&mut parsed)
            .into_iter()
            .zip(probabilities(&mut direct))
        {
            assert!(
                p.zip(q).map_or(p == q, |(p, q)| (p - q).abs() <= 5e-7),
                "{p:?} vs {q:?}"
            );
        }
        assert_eq!(parsed, direct);
        let truth = Truth::from_read(&sim, read, &id);
        let header = Truth::tsv_header();
        let truth_columns: Vec<&str> = header.split('\t').collect();
        assert_eq!(
            Truth::from_tsv(&truth_columns, &truth.to_tsv()).expect("truth row"),
            truth
        );
    }
}

/// Accuracy at scale: 20,000 reads per IgSim preset. A rate measured on 500
/// reads has a 95 % interval of ±1–4 pp, which hides most regressions; on
/// 20,000 reads three standard errors are 0.1–0.9 pp. Run with
/// `cargo test --release -- --ignored`.
#[test]
#[ignore = "20,000 reads per preset; run with `cargo test --release -- --ignored`"]
fn accuracy_at_scale() {
    let setup = Setup::new();
    // (preset, seed, measured V, V gene, D, J, junction; V and J ECE bounds)
    let presets = [
        (
            SimParams::easy(),
            21,
            [100.00, 100.00, 99.50, 99.95, 100.00],
            0.2,
            0.5,
        ),
        (
            SimParams::mid(),
            22,
            [99.53, 100.00, 84.08, 99.41, 99.99],
            0.3,
            0.5,
        ),
        (
            SimParams::hard(),
            23,
            [97.62, 99.92, 75.50, 95.95, 99.73],
            0.5,
            0.8,
        ),
    ];
    for (params, seed, measured, v_ece, j_ece) in presets {
        let e = setup.evaluate(params, 20_000, seed, &EngineConfig::default());
        segments_at_least(&e, measured);
        assert!(100.0 * e.calibration.v.ece() <= v_ece, "seed {seed} V ECE");
        assert!(100.0 * e.calibration.j.ece() <= j_ece, "seed {seed} J ECE");
        let random_calls = e.unassigned;
        assert!(
            random_calls <= e.reads / 1000,
            "seed {seed}: {random_calls} unassigned"
        );
    }
}
