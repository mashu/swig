//! Accuracy and throughput on simulated reads (IgSim model).
//!
//! ```text
//! cargo run --release --example benchmark -- [--v V.fasta --d D.fasta --j J.fasta]
//!     [--preset easy|mid|hard] [--shm RATE] [--reads N] [--seed S] [--rc P]
//!     [--profile r|truth|igblast|sensitive-d] [--threads N]
//!     [--simd avx512|avx2|neon|portable|scalar]
//! ```
//!
//! Without germline files a synthetic IGH-like database is generated
//! ([`swig::sim::SyntheticDatabase`]). Reads are reverse-complemented
//! with probability `--rc` (default 0.5). Throughput is measured on the
//! given number of threads (default 1).

mod common;

use std::error::Error;
use std::time::Instant;

use common::Args;
use swig::prelude::*;
use swig::sim::{Assignment, Evaluation, SimParams, Simulator, Truth};

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::from_env();
    let value = |name: &str| args.value(name);
    let reads: usize = args.number("--reads", 20000)?;
    let seed: u64 = args.number("--seed", 1)?;
    let threads: usize = args.number("--threads", 1)?;
    let rc: f64 = args.number("--rc", 0.5)?;
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
    if let Some(name) = value("--simd") {
        config.simd = Some(swig::align::Backend::from_name(&name).ok_or("unknown --simd backend")?);
    }
    if let Some(penalty) = value("--p-d") {
        config.objective.d_prior = swig::config::DPrior::Fixed(penalty.parse()?);
    }
    if let Some(e) = value("--anchor-evalue") {
        let e: f64 = e.parse()?;
        let gates = (e > 0.0).then_some(swig::segment::EvalueGates {
            anchor: e,
            ..Default::default()
        });
        config.v.evalue = gates;
        config.j.evalue = gates;
    }
    if let Some(base) = value("--budget-base") {
        config.v.budget.base = base.parse()?;
    }
    if let Some(max) = value("--budget-max") {
        config.v.budget.max = max.parse()?;
    }
    // V gene margin: a gene is kept while its votes are within z standard
    // deviations of the leader's (AmbiguityBudget).
    if let Some(z) = value("--margin-z") {
        config.v.budget.z = z.parse()?;
    }
    if let Some(limit) = value("--max-occ") {
        config.v.seeds.max_occurrences = limit.parse()?;
    }
    if let Some(reclaim) = value("--reclaim") {
        let reclaim: u32 = reclaim.parse()?;
        config.d_window.reclaim = reclaim;
        config.j.window.reclaim = reclaim;
    }
    if let Some(sampling) = value("--sampling") {
        // all | closed:W | open:W | min:W
        use swig::seed::SeedSampling;
        let (kind, w) = match sampling.split_once(':') {
            Some((kind, w)) => (kind, w.parse::<u8>()?),
            None if sampling == "all" => ("all", 0),
            None => {
                return Err(
                    format!("--sampling {sampling} needs a window, e.g. {sampling}:4").into(),
                );
            }
        };
        config.v.seeds.sampling = match kind {
            "all" => SeedSampling::All,
            "closed" => SeedSampling::ClosedSyncmer { w },
            "open" => SeedSampling::OpenSyncmer { w },
            "min" => SeedSampling::Minimizer { w },
            other => return Err(format!("unknown sampling '{other}'").into()),
        };
    }
    if let Some(seeds) = value("--weak-seeds") {
        config.v.seeds.weak_seeds = seeds.parse()?;
    }
    if let Some(reclaim) = value("--j-reclaim") {
        // How far before the earliest V end J is searched (huge = whole read).
        config.j.window.reclaim = reclaim.parse()?;
    }
    if let Some(pad) = value("--band-pad") {
        config.v.seeds.band_pad = pad.parse()?;
    }
    // Sensitivity analysis of the remaining constants (comma-separated values).
    let list = |name: &str| -> Option<Vec<String>> {
        value(name).map(|v| v.split(',').map(str::to_string).collect())
    };
    if let Some(bonus) = value("--pair-bonus") {
        config.objective.pair_bonus = bonus.parse()?;
    }
    for (name, accept) in [
        ("--v-accept", &mut config.v.accept),
        ("--j-accept", &mut config.j.accept),
    ] {
        if let Some(v) = list(name) {
            accept.min_subst = v[0].parse()?;
            accept.min_identity = v[1].parse()?;
        }
    }
    if let Some(v) = list("--seeds") {
        config.v.primary = swig::seed::SeedShape::contiguous(v[0].parse()?);
        let fallback: usize = v[1].parse()?;
        config.v.fallback = (fallback > 0).then(|| swig::seed::SeedShape::contiguous(fallback));
    }
    // Primary V seed as a 0/1 pattern, e.g. a spaced seed `11101101111`.
    if let Some(pattern) = value("--v-shape") {
        config.v.primary = swig::seed::SeedShape::from_pattern(&pattern)?;
    }
    if let Some(v) = list("--budget") {
        if v.len() != 3 {
            return Err("--budget expects base,max,z".into());
        }
        config.v.budget.base = v[0].parse()?;
        config.v.budget.max = v[1].parse()?;
        config.v.budget.z = v[2].parse()?;
    }
    if let Some(run) = value("--min-d-run") {
        config.d.min_exact_run = run.parse()?;
    }
    if let Some(v) = list("--dist-d")
        && let Some(gate) = config.d.distributed.as_mut()
    {
        gate.min_columns = v[0].parse()?;
        gate.min_matches = v[1].parse()?;
        gate.min_identity = v[2].parse()?;
        gate.min_rescore = v[3].parse()?;
    }
    if let Some(v) = list("--min-retained") {
        config.min_retained = [v[0].parse()?, v[1].parse()?, v[2].parse()?];
    }
    for ablation in value("--ablate")
        .unwrap_or_default()
        .split(',')
        .filter(|a| !a.is_empty())
    {
        match ablation {
            "family-delta" => {
                config.v.family_delta = false;
                config.j.family_delta = false;
            }
            "anchor-window" => config.anchor_window = false,
            "distributed-d" => config.d.distributed = None,
            "posterior" => config.posterior = false,
            other => return Err(format!("unknown ablation '{other}'").into()),
        }
    }

    let db = args.database(seed)?;
    let sim = Simulator::new(&db.v, db.d.as_ref(), &db.j, params)?;
    let simulated = sim.reads(reads, seed)?;
    let assigner = DefaultAssigner::from_config(&config, db.v.clone(), db.d.clone(), db.j.clone())?;
    // `--tail N`: keep the last N bases of each read (germline sense only),
    // as a short read covering part of V.
    let tail: Option<usize> = value("--tail").map(|v| v.parse()).transpose()?;
    if tail.is_some() && rc != 0.0 {
        return Err("--tail needs --rc 0".into());
    }
    let queries: Vec<Vec<u8>> = simulated
        .iter()
        .map(|r| {
            let keep = tail.map_or(0, |n| r.sequence.len().saturating_sub(n));
            encode(&r.sequence[keep..])
        })
        .collect();

    let started = Instant::now();
    let calls = assign_all(&assigner, &queries, threads);
    let seconds = started.elapsed().as_secs_f64();

    let mut evaluation = Evaluation::default();
    for (read, call) in simulated.iter().zip(&calls) {
        let sequence = String::from_utf8_lossy(&read.sequence);
        let record = AirrRecord::from_call("r", &sequence, call, &assigner, &[None, None, None]);
        let truth = Truth::from_read(&sim, read, "r");
        evaluation.add(&truth, &Assignment::from_record(&record));
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
        let mut rng = swig::sim::Rng::new(seed ^ 0xA5A5);
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
            path += usize::from(call.path.is_some());
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
