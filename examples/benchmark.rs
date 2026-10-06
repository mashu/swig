//! Accuracy and throughput on simulated reads (IgSim model).
//!
//! ```text
//! cargo run --release --example benchmark -- [--v V.fasta --d D.fasta --j J.fasta]
//!     [--preset easy|mid|hard] [--shm RATE] [--reads N] [--seed S] [--rc P]
//!     [--profile r|truth|igblast|sensitive-d] [--tandem-d] [--threads N]
//! ```
//!
//! Without germline files a synthetic IGH-like database is generated
//! ([`vdj_assign::sim::SyntheticDatabase`]). Reads are reverse-complemented
//! with probability `--rc` (default 0.5). Throughput is measured on the
//! given number of threads (default 1).

use std::error::Error;
use std::time::Instant;

use vdj_assign::annotate::AnnotateOptions;
use vdj_assign::prelude::*;
use vdj_assign::prepare::prepare;
use vdj_assign::sim::{Evaluation, SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

fn load(kind: SegmentKind, path: &str) -> Result<GermlineSet, Box<dyn Error>> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    Ok(GermlineSet::from_fasta(kind, &text)?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let value = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let flag = |name: &str| args.iter().any(|a| a == name);
    let reads: usize = value("--reads").map_or(Ok(5000), |v| v.parse())?;
    let seed: u64 = value("--seed").map_or(Ok(1), |v| v.parse())?;
    let threads: usize = value("--threads").map_or(Ok(1), |v| v.parse())?;
    let rc: f64 = value("--rc").map_or(Ok(0.5), |v| v.parse())?;
    let preset = value("--preset").unwrap_or_else(|| "mid".to_string());
    let mut params = SimParams::preset(&preset)
        .ok_or("unknown preset")?
        .with_reverse_complement(rc);
    if let Some(rate) = value("--shm") {
        params = params.with_shm(rate.parse()?);
    }
    let profile = match value("--profile") {
        Some(name) => CallingProfile::from_name(&name).ok_or("unknown profile")?,
        None => CallingProfile::ROptimized,
    };
    let mut config = profile.config();
    if flag("--tandem-d") {
        config.grammar = VdjGrammar::TANDEM_D;
    }
    if let Some(penalty) = value("--p-d") {
        config.objective.d_prior = vdj_assign::config::DPrior::Fixed(penalty.parse()?);
    }
    if let Some(e) = value("--anchor-evalue") {
        let e: f64 = e.parse()?;
        let gates = (e > 0.0).then_some(vdj_assign::segment::EvalueGates {
            anchor: e,
            ..Default::default()
        });
        config.v.evalue = gates;
        config.j.evalue = gates;
    }
    if let Some(limit) = value("--max-alleles") {
        config.v.seeds.max_alleles = limit.parse()?;
    }
    if let Some(limit) = value("--max-occ") {
        config.v.seeds.max_occurrences = limit.parse()?;
    }
    if let Some(reclaim) = value("--reclaim") {
        let reclaim: u32 = reclaim.parse()?;
        config.d_window.reclaim = reclaim;
        config.j.window.reclaim = reclaim;
    }
    if let Some(pad) = value("--band-pad") {
        config.v.seeds.band_pad = pad.parse()?;
    }
    for ablation in value("--ablate")
        .unwrap_or_default()
        .split(',')
        .filter(|a| !a.is_empty())
    {
        match ablation {
            "family-delta" => config.v.family_delta = false,
            "anchor-window" => config.anchor_window = false,
            "distributed-d" => config.d.distributed = None,
            other => return Err(format!("unknown ablation '{other}'").into()),
        }
    }

    let (v, d, j) = match (value("--v"), value("--j")) {
        (Some(v), Some(j)) => (
            load(SegmentKind::V, &v)?,
            value("--d").map(|d| load(SegmentKind::D, &d)).transpose()?,
            load(SegmentKind::J, &j)?,
        ),
        _ => {
            let db = SyntheticDatabase::generate(&SyntheticSpec::default(), seed);
            (db.v, Some(db.d), db.j)
        }
    };
    let db = prepare(v, d, j, AnnotateOptions::default())?;
    eprintln!("database: {}", db.summary());
    let sim = Simulator::new(&db.v, db.d.as_ref(), &db.j, params)?;
    let simulated = sim.reads(reads, seed)?;
    let assigner = DefaultAssigner::from_config(&config, db.v.clone(), db.d.clone(), db.j.clone())?;
    let queries: Vec<Vec<u8>> = simulated.iter().map(|r| encode(&r.sequence)).collect();

    let started = Instant::now();
    let calls = assign_all(&assigner, &queries, threads);
    let seconds = started.elapsed().as_secs_f64();

    let mut evaluation = Evaluation::default();
    for (read, call) in simulated.iter().zip(&calls) {
        let sequence = String::from_utf8_lossy(&read.sequence);
        let record = AirrRecord::from_call("r", &sequence, call, &assigner, &[None, None, None]);
        evaluation.add(&sim, read, &record);
    }
    println!(
        "preset {preset}  profile {}  {} reads  {:.0} reads/s on {threads} thread(s)",
        profile.name(),
        reads,
        reads as f64 / seconds
    );
    println!("{evaluation}");

    // Negative control: random reads with the simulated length distribution.
    if let Some(n) = value("--random") {
        let n: usize = n.parse()?;
        let mut rng = vdj_assign::sim::Rng::new(seed ^ 0xA5A5);
        let random: Vec<Vec<u8>> = simulated
            .iter()
            .cycle()
            .take(n)
            .map(|r| (0..r.sequence.len()).map(|_| rng.base()).collect())
            .collect();
        let queries: Vec<Vec<u8>> = random.iter().map(|r| encode(r)).collect();
        let calls = assign_all(&assigner, &queries, threads);
        let (mut v, mut j, mut path) = (0, 0, 0);
        for call in &calls {
            let record = AirrRecord::from_call("r", "", call, &assigner, &[None, None, None]);
            v += usize::from(record.v.is_some());
            j += usize::from(record.j.is_some());
            path += usize::from(record.path_score.is_some());
        }
        println!(
            "random reads {n}: with a path {:.2}%  V called {:.2}%  J called {:.2}%",
            100.0 * path as f64 / n as f64,
            100.0 * v as f64 / n as f64,
            100.0 * j as f64 / n as f64
        );
    }
    Ok(())
}
