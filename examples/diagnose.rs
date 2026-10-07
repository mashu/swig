//! Split V/J allele errors into search errors and model errors.
//!
//! For every read the pipeline's call is compared with an exhaustive search
//! under the same objective: every allele aligned unbanded (local Gotoh plus
//! the co-linear outer flank, as the V/J stages rank them). If the exhaustive
//! optimum is germline-equivalent to the truth but the pipeline's call is
//! not, the pipeline lost it (search error). If the exhaustive optimum is
//! itself wrong, the objective cannot separate the alleles (model error).
//!
//! ```text
//! cargo run --release --example diagnose -- --v V.fasta --d D.fasta --j J.fasta --preset hard [--reads N] [--show K]
//! ```

use std::error::Error;

use vdj_assign::align::{Aligner, InterSequenceWorkspace};
use vdj_assign::annotate::AnnotateOptions;
use vdj_assign::prelude::*;
use vdj_assign::prepare::prepare;
use vdj_assign::refine::{OuterExtend, extend_outer};
use vdj_assign::scoring::GapCosts;
use vdj_assign::scoring::SubstitutionProfile;
use vdj_assign::sim::{Region, SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

fn load(kind: SegmentKind, path: &str) -> Result<GermlineSet, Box<dyn Error>> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    Ok(GermlineSet::from_fasta(kind, &text)?)
}

/// Exhaustive ranking of every allele: (score, allele), best first.
fn exhaustive(
    query: &[u8],
    set: &GermlineSet,
    profiles: &[SubstitutionProfile],
    gaps: GapCosts,
    outer: OuterExtend,
    ws: &mut InterSequenceWorkspace,
) -> Vec<(i32, usize)> {
    let jobs: Vec<AlignJob<'_>> = profiles
        .iter()
        .map(|profile| AlignJob {
            query,
            profile,
            gaps,
            band: Band::Full,
        })
        .collect();
    let mut out = Vec::new();
    InterSequenceGotoh::new().align_batch(&jobs, ws, &mut out);
    let mut ranked: Vec<(i32, usize)> = out
        .into_iter()
        .enumerate()
        .filter_map(|(i, alignment)| {
            let mut alignment = alignment?;
            extend_outer(
                &mut alignment,
                query,
                &set.alleles()[i].codes,
                &profiles[i],
                outer,
            );
            Some((alignment.score, i))
        })
        .collect();
    ranked.sort_by_key(|&(s, i)| (std::cmp::Reverse(s), i));
    ranked
}

fn equivalent(set: &GermlineSet, body: &[u8]) -> Vec<usize> {
    let codes = encode(body);
    set.alleles()
        .iter()
        .enumerate()
        .filter(|(_, a)| a.codes.windows(codes.len()).any(|w| w == codes.as_slice()))
        .map(|(i, _)| i)
        .collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let value = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let reads: usize = value("--reads").map_or(Ok(2000), |v| v.parse())?;
    let show: usize = value("--show").map_or(Ok(5), |v| v.parse())?;
    let seed: u64 = value("--seed").map_or(Ok(1), |v| v.parse())?;
    let preset = value("--preset").unwrap_or_else(|| "hard".to_string());
    let params = SimParams::preset(&preset)
        .ok_or("unknown preset")?
        .with_reverse_complement(0.5);
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
    let config = CallingProfile::ROptimized.config();
    let sim = Simulator::new(&db.v, db.d.as_ref(), &db.j, params)?;
    let simulated = sim.reads(reads, seed)?;
    let assigner = DefaultAssigner::from_config(&config, db.v.clone(), db.d.clone(), db.j.clone())?;
    let mut ws = assigner.workspace();
    let mut simd = InterSequenceWorkspace::default();

    for (kind, region, set, scoring, outer) in [
        (
            SegmentKind::V,
            Region::V,
            &db.v,
            config.v.scoring,
            OuterExtend::Prefix,
        ),
        (
            SegmentKind::J,
            Region::J,
            &db.j,
            config.j.scoring,
            OuterExtend::Suffix,
        ),
    ] {
        let profiles: Vec<SubstitutionProfile> = set
            .alleles()
            .iter()
            .map(|a| SubstitutionProfile::build(&scoring, &a.codes))
            .collect();
        let (mut search, mut model, mut both_ok, mut orientation, mut ties) = (0, 0, 0, 0, 0);
        let mut shown = 0;
        for read in &simulated {
            let call = assigner.assign(&encode(&read.sequence), &mut ws);
            if call.orientation.is_reverse() != read.reverse_complemented {
                orientation += 1;
                continue;
            }
            let truth_allele = if kind == SegmentKind::V {
                read.v
            } else {
                read.j
            };
            let body = sim
                .germline_body(read, region)
                .expect("V and J always present");
            let ok = equivalent(set, body);
            let forward = encode(&read.forward);
            let ranked = exhaustive(&forward, set, &profiles, scoring.gaps, outer, &mut simd);
            let best = ranked[0];
            let called = call.path.as_ref().and_then(|p| {
                p.segments
                    .iter()
                    .map(|s| &call.placements[s.placement])
                    .find(|pl| pl.kind == kind)
                    .map(|pl| pl.primary_allele() as usize)
            });
            let exhaustive_ok = ranked
                .iter()
                .take_while(|r| r.0 == best.0)
                .any(|r| ok.contains(&r.1));
            let pipeline_ok = called.is_some_and(|c| ok.contains(&c));
            let called_score = called.and_then(|c| ranked.iter().find(|r| r.1 == c).map(|r| r.0));
            // The called allele ties the exhaustive optimum: no objective-based
            // caller can do better than report the tie.
            let tie = called_score == Some(best.0) && exhaustive_ok;
            match (pipeline_ok, exhaustive_ok) {
                (true, _) => both_ok += 1,
                (false, true) if tie => ties += 1,
                (false, true) => {
                    search += 1;
                    if shown < show {
                        shown += 1;
                        let score_of = |a: usize| ranked.iter().find(|r| r.1 == a).map(|r| r.0);
                        let chosen = call.path.as_ref().and_then(|p| {
                            p.segments
                                .iter()
                                .find(|s| call.placements[s.placement].kind == kind)
                                .map(|s| {
                                    (
                                        set.allele(call.placements[s.placement].primary_allele())
                                            .name
                                            .clone(),
                                        s.span.score,
                                        s.span.q_lo,
                                        s.span.q_hi,
                                    )
                                })
                        });
                        let truth_span = read.span(region);
                        eprintln!("    DAG chose {chosen:?}; truth span {truth_span:?}");
                        let placed: Vec<String> = call
                            .placements
                            .iter()
                            .filter(|p| p.kind == kind)
                            .map(|p| {
                                format!("{}:{}", set.allele(p.primary_allele()).name, p.score())
                            })
                            .collect();
                        eprintln!(
                            "{} search error: truth {} (exhaustive {:?}), called {:?} (exhaustive {:?}); shm {:.3}, placements {:?}",
                            kind.label(),
                            set.allele(truth_allele as u32).name,
                            score_of(truth_allele),
                            called.map(|c| &set.allele(c as u32).name),
                            called.and_then(score_of),
                            read.shm_rate,
                            placed
                        );
                    }
                }
                (false, false) => model += 1,
            }
        }
        println!(
            "{}: correct {both_ok}, ties {ties}, search errors {search}, model errors {model}, mis-oriented {orientation} (of {})",
            kind.label(),
            simulated.len()
        );
    }
    Ok(())
}
