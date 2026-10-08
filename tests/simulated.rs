//! Accuracy regression tests on simulated repertoires.
//!
//! Reads come from the IgSim model ([`vdj_assign::sim`]) on a synthetic
//! IGH-like germline database with a family → gene → allele hierarchy, half
//! of them reverse-complemented. Seeds are fixed, so every run sees the same
//! reads; the floors sit a little below the measured accuracy, so a change
//! that loses accuracy fails here. Accuracy is germline-equivalence aware
//! (see `vdj_assign::sim::evaluate`).

use vdj_assign::prelude::*;
use vdj_assign::sim::{
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
        DefaultAssigner::from_config(
            config,
            self.db.v.clone(),
            Some(self.db.d.clone()),
            self.db.j.clone(),
        )
        .expect("synthetic database builds")
    }

    fn evaluate(
        &self,
        params: SimParams,
        reads: usize,
        seed: u64,
        config: &EngineConfig,
    ) -> Evaluation {
        let sim = Simulator::new(
            &self.db.v,
            Some(&self.db.d),
            &self.db.j,
            params.with_reverse_complement(0.5),
        )
        .expect("simulator");
        let simulated = sim.reads(reads, seed).expect("reads");
        let assigner = self.assigner(config);
        let queries: Vec<Vec<u8>> = simulated.iter().map(|r| encode(&r.sequence)).collect();
        let calls = assign_all(&assigner, &queries, 0);
        let mut evaluation = Evaluation::default();
        for (read, call) in simulated.iter().zip(&calls) {
            let record = AirrRecord::from_call("r", "", call, &assigner, &[None, None, None]);
            let truth = Truth::from_read(&sim, read, "r");
            evaluation.add(&truth, &Assignment::from_record(&record));
        }
        println!("{evaluation}");
        evaluation
    }
}

fn percent(part: usize, whole: usize) -> f64 {
    100.0 * part as f64 / whole as f64
}

#[test]
fn unmutated_reads_are_assigned_exactly() {
    let e = Setup::new().evaluate(SimParams::easy(), 300, 1, &EngineConfig::default());
    assert_eq!(e.oriented, e.reads, "orientation");
    assert_eq!(e.unassigned, 0);
    assert!(
        e.v.equivalent_rate() >= 99.5,
        "V {:.2}",
        e.v.equivalent_rate()
    );
    assert!(
        e.j.equivalent_rate() >= 99.5,
        "J {:.2}",
        e.j.equivalent_rate()
    );
    assert!(
        e.d.equivalent_rate() >= 97.0,
        "D {:.2}",
        e.d.equivalent_rate()
    );
    assert!(
        e.junction_rate() >= 99.5,
        "junction {:.2}",
        e.junction_rate()
    );
}

#[test]
fn igm_like_reads_keep_their_accuracy() {
    let e = Setup::new().evaluate(SimParams::mid(), 500, 2, &EngineConfig::default());
    assert!(percent(e.oriented, e.reads) >= 99.8);
    assert!(
        e.v.equivalent_rate() >= 98.5,
        "V {:.2}",
        e.v.equivalent_rate()
    );
    assert!(e.v.gene_rate() >= 99.5, "V gene {:.2}", e.v.gene_rate());
    assert!(
        e.j.equivalent_rate() >= 98.0,
        "J {:.2}",
        e.j.equivalent_rate()
    );
    assert!(
        e.d.equivalent_rate() >= 82.5,
        "D {:.2}",
        e.d.equivalent_rate()
    );
    assert!(
        e.junction_rate() >= 99.0,
        "junction {:.2}",
        e.junction_rate()
    );
}

#[test]
fn igg_like_hypermutated_reads_keep_their_accuracy() {
    let e = Setup::new().evaluate(SimParams::hard(), 500, 3, &EngineConfig::default());
    assert!(percent(e.oriented, e.reads) >= 99.6);
    assert!(
        e.v.equivalent_rate() >= 96.0,
        "V {:.2}",
        e.v.equivalent_rate()
    );
    assert!(e.v.gene_rate() >= 99.0, "V gene {:.2}", e.v.gene_rate());
    assert!(
        e.j.equivalent_rate() >= 95.5,
        "J {:.2}",
        e.j.equivalent_rate()
    );
    assert!(
        e.d.equivalent_rate() >= 74.0,
        "D {:.2}",
        e.d.equivalent_rate()
    );
    assert!(
        e.junction_rate() >= 98.5,
        "junction {:.2}",
        e.junction_rate()
    );
}

#[test]
fn long_d_remnants_are_found() {
    let e = Setup::new().evaluate(SimParams::mid(), 500, 4, &EngineConfig::default());
    // Bins 12–15 and 16+ (see `D_LENGTH_BINS`): a D this long is identifiable.
    let long: Vec<_> = e.d_presence[4..].iter().collect();
    let reads: usize = long.iter().map(|p| p.reads).sum();
    let correct: usize = long.iter().map(|p| p.correct).sum();
    assert!(percent(correct, reads) >= 98.0, "long D {correct}/{reads}");
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

/// Sibling bounds (`family_delta`, V and J) are lossless under every stage
/// configuration they interact with: banding, the ungapped prefilter, `top`,
/// `cap`, the band width and V acceptance. Calls with and without them are
/// identical, down to every placement and posterior.
#[test]
fn family_bounds_agree_with_full_scoring_under_config_variants() {
    use vdj_assign::retrieve::UngappedParams;
    use vdj_assign::segment::HitAcceptance;

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
            record.path_score.is_some(),
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
/// 20,000 reads it is ±0.1–0.6 pp, so these floors sit within ~1 pp of the
/// measured values. Run with `cargo test --release -- --ignored`.
#[test]
#[ignore = "20,000 reads per preset; run with `cargo test --release -- --ignored`"]
fn accuracy_at_scale() {
    let setup = Setup::new();
    // (preset, seed, V, V gene, D, J, junction, V ECE, J ECE)
    let presets = [
        // Measured: 99.99 100 99.53 99.95 100, ECE 0.02 0.25.
        (
            SimParams::easy(),
            21,
            99.8,
            99.9,
            99.2,
            99.8,
            99.9,
            0.2,
            0.5,
        ),
        // Measured: 99.53 100 84.50 99.40 99.99, ECE 0.04 0.20.
        (SimParams::mid(), 22, 99.3, 99.9, 83.5, 99.1, 99.9, 0.3, 0.5),
        // Measured: 97.57 99.85 75.86 95.93 99.71, ECE 0.13 0.40.
        (
            SimParams::hard(),
            23,
            97.0,
            99.6,
            74.8,
            95.3,
            99.4,
            0.5,
            0.8,
        ),
    ];
    for (params, seed, v, v_gene, d, j, junction, v_ece, j_ece) in presets {
        let e = setup.evaluate(params, 20_000, seed, &EngineConfig::default());
        let checks = [
            ("V", e.v.equivalent_rate(), v),
            ("V gene", e.v.gene_rate(), v_gene),
            ("D", e.d.equivalent_rate(), d),
            ("J", e.j.equivalent_rate(), j),
            ("junction", e.junction_rate(), junction),
        ];
        for (label, measured, floor) in checks {
            assert!(
                measured >= floor,
                "seed {seed} {label} {measured:.2} < {floor}"
            );
        }
        assert!(100.0 * e.calibration.v.ece() <= v_ece, "seed {seed} V ECE");
        assert!(100.0 * e.calibration.j.ece() <= j_ece, "seed {seed} J ECE");
        let random_calls = e.unassigned;
        assert!(
            random_calls <= e.reads / 1000,
            "seed {seed}: {random_calls} unassigned"
        );
    }
}
