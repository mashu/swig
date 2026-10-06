//! Accuracy regression tests on simulated repertoires.
//!
//! Reads come from the IgSim model ([`vdj_assign::sim`]) on a synthetic
//! IGH-like germline database with a family → gene → allele hierarchy, half
//! of them reverse-complemented. Seeds are fixed, so every run sees the same
//! reads; the floors sit a little below the measured accuracy, so a change
//! that loses accuracy fails here. Accuracy is germline-equivalence aware
//! (see `vdj_assign::sim::evaluate`).

use vdj_assign::prelude::*;
use vdj_assign::sim::{Evaluation, Rng, SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

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
        DefaultAssigner::from_config(config, self.db.v.clone(), Some(self.db.d.clone()), self.db.j.clone())
            .expect("synthetic database builds")
    }

    fn evaluate(&self, params: SimParams, reads: usize, seed: u64, config: &EngineConfig) -> Evaluation {
        let sim = Simulator::new(&self.db.v, Some(&self.db.d), &self.db.j, params.with_reverse_complement(0.5))
            .expect("simulator");
        let simulated = sim.reads(reads, seed).expect("reads");
        let assigner = self.assigner(config);
        let queries: Vec<Vec<u8>> = simulated.iter().map(|r| encode(&r.sequence)).collect();
        let calls = assign_all(&assigner, &queries, 0);
        let mut evaluation = Evaluation::default();
        for (read, call) in simulated.iter().zip(&calls) {
            let record = AirrRecord::from_call("r", "", call, &assigner, &[None, None, None]);
            evaluation.add(&sim, read, &record);
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
    assert!(e.v.equivalent_rate() >= 99.5, "V {:.2}", e.v.equivalent_rate());
    assert!(e.j.equivalent_rate() >= 99.5, "J {:.2}", e.j.equivalent_rate());
    assert!(e.d.equivalent_rate() >= 97.0, "D {:.2}", e.d.equivalent_rate());
    assert!(e.junction_rate() >= 99.5, "junction {:.2}", e.junction_rate());
}

#[test]
fn igm_like_reads_keep_their_accuracy() {
    let e = Setup::new().evaluate(SimParams::mid(), 500, 2, &EngineConfig::default());
    assert!(percent(e.oriented, e.reads) >= 99.8);
    assert!(e.v.equivalent_rate() >= 98.5, "V {:.2}", e.v.equivalent_rate());
    assert!(e.v.gene_rate() >= 99.5, "V gene {:.2}", e.v.gene_rate());
    assert!(e.j.equivalent_rate() >= 98.0, "J {:.2}", e.j.equivalent_rate());
    assert!(e.d.equivalent_rate() >= 82.5, "D {:.2}", e.d.equivalent_rate());
    assert!(e.junction_rate() >= 99.0, "junction {:.2}", e.junction_rate());
}

#[test]
fn igg_like_hypermutated_reads_keep_their_accuracy() {
    let e = Setup::new().evaluate(SimParams::hard(), 500, 3, &EngineConfig::default());
    assert!(percent(e.oriented, e.reads) >= 99.6);
    assert!(e.v.equivalent_rate() >= 96.0, "V {:.2}", e.v.equivalent_rate());
    assert!(e.v.gene_rate() >= 99.0, "V gene {:.2}", e.v.gene_rate());
    assert!(e.j.equivalent_rate() >= 95.5, "J {:.2}", e.j.equivalent_rate());
    assert!(e.d.equivalent_rate() >= 74.0, "D {:.2}", e.d.equivalent_rate());
    assert!(e.junction_rate() >= 98.5, "junction {:.2}", e.junction_rate());
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
