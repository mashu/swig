//! Attribute every D error on simulated reads to its cause.
//!
//! For each read with a simulated D, the call is classified as correct
//! (primary allele germline-equivalent) or wrong / missing, and each error by
//! where it arose:
//!
//! * **outside windows** — the true D is not inside the D search windows
//!   (the windows are recomputed from the call's own V and J placements);
//! * **not found** — inside the windows, but no D placement of an equivalent
//!   allele overlaps the true D (prefilter, admission, score floor, keep);
//! * **adjudication** — such a placement exists, but the path took another D
//!   or none.
//!
//! Missing D calls are further split by which called segment covers the true
//! D's bases: the J (it "ate" the D), the V, or neither (left to N by the D
//! prior). A J-covered miss whose window was cut short by the J placement is
//! the failure mode of a J alignment eating into the D.
//!
//! ```text
//! cargo run --release --example d_errors -- --v V.fasta --d D.fasta --j J.fasta [--preset hard] [--reads N]
//!     [--no-anchor-window] [--j-reclaim N]
//! ```

use std::error::Error;
use std::ops::Range;

use vdj_assign::annotate::AnnotateOptions;
use vdj_assign::placement::Placement;
use vdj_assign::prelude::*;
use vdj_assign::prepare::prepare;
use vdj_assign::sim::{SimParams, Simulator, Slot, SyntheticDatabase, SyntheticSpec, Truth};
use vdj_assign::window::DWindowPolicy;

fn load(kind: SegmentKind, path: &str) -> Result<GermlineSet, Box<dyn Error>> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    Ok(GermlineSet::from_fasta(kind, &text)?)
}

fn overlaps(a: &Range<usize>, b: &Range<usize>) -> bool {
    a.start < b.end && b.start < a.end
}

#[derive(Default)]
struct Tally {
    reads: usize,
    correct: usize,
    wrong_adjudication: usize,
    wrong_not_found: usize,
    wrong_outside: usize,
    missing_adjudication: [usize; 3],
    missing_not_found: [usize; 3],
    missing_outside: [usize; 3],
    /// Correct or not: the called J segment covers ≥ 3 bases of the true D.
    j_covers_d: usize,
    v_covers_d: usize,
}

const COVER: [&str; 3] = ["J covers", "V covers", "neither"];

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let value = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let reads: usize = value("--reads").map_or(Ok(5000), |v| v.parse())?;
    let seed: u64 = value("--seed").map_or(Ok(1), |v| v.parse())?;
    let preset = value("--preset").unwrap_or_else(|| "hard".to_string());
    let params = SimParams::preset(&preset)
        .ok_or("unknown preset")?
        .with_reverse_complement(0.0);
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
    let mut config = CallingProfile::ROptimized.config();
    // Counterfactuals: D windows from the V/J placements alone, and J
    // searched over the whole read.
    config.anchor_window = !args.iter().any(|a| a == "--no-anchor-window");
    if let Some(reclaim) = value("--j-reclaim") {
        config.j.window.reclaim = reclaim.parse()?;
    }
    let sim = Simulator::new(&db.v, db.d.as_ref(), &db.j, params)?;
    let assigner = DefaultAssigner::from_config(&config, db.v.clone(), db.d.clone(), db.j.clone())?;
    let d_set = db.d.as_ref().ok_or("D database required")?;
    let mut ws = assigner.workspace();

    // Detectable (≥ 8 nt) and short remnants separately.
    let mut tallies = [Tally::default(), Tally::default()];
    let mut windows = Vec::new();
    for (i, read) in sim.reads(reads, seed)?.iter().enumerate() {
        let truth = Truth::from_read(&sim, read, &i.to_string());
        let Some(true_d) = truth.segment(Slot::D) else {
            continue;
        };
        let call = assigner.assign(&encode(&read.sequence), &mut ws);
        if call.orientation.is_reverse() {
            continue;
        }
        let tally = &mut tallies[usize::from(true_d.span.len() < 8)];
        tally.reads += 1;
        let record = AirrRecord::from_call("x", "", &call, &assigner, &[None, None, None]);
        let span = |r: &Option<vdj_assign::report::SegmentReport>| {
            r.as_ref()
                .map(|s| s.sequence_start as usize - 1..s.sequence_end as usize)
        };
        let covers = |segment: &Option<Range<usize>>| {
            segment
                .as_ref()
                .is_some_and(|s| s.start.max(true_d.span.start) + 3 <= s.end.min(true_d.span.end))
        };
        let (j_covers, v_covers) = (covers(&span(&record.j)), covers(&span(&record.v)));
        tally.j_covers_d += usize::from(j_covers);
        tally.v_covers_d += usize::from(v_covers);
        let primary = record
            .d
            .as_ref()
            .and_then(|d| d.call.split(',').next().map(str::to_string));
        if primary
            .as_ref()
            .is_some_and(|p| true_d.equivalent.contains(p))
        {
            tally.correct += 1;
            continue;
        }

        let found = call.placements.iter().any(|p| {
            p.kind == SegmentKind::D
                && overlaps(&(p.q_start() as usize..p.q_end() as usize), &true_d.span)
                && p.alleles
                    .iter()
                    .any(|&a| p_names(d_set, a).any(|n| true_d.equivalent.iter().any(|e| e == n)))
        });
        let kind_of = |k: SegmentKind| -> Vec<Placement> {
            call.placements
                .iter()
                .filter(|p| p.kind == k)
                .cloned()
                .collect()
        };
        windows.clear();
        assigner.d_window.windows(
            call.query.len(),
            &kind_of(SegmentKind::V),
            &kind_of(SegmentKind::J),
            &mut windows,
        );
        let inside = windows
            .iter()
            .any(|w| w.start <= true_d.span.start && true_d.span.end <= w.end);
        let cover = if j_covers {
            0
        } else if v_covers {
            1
        } else {
            2
        };
        match (primary.is_some(), found, inside) {
            (true, true, _) => tally.wrong_adjudication += 1,
            (true, false, true) => tally.wrong_not_found += 1,
            (true, false, false) => tally.wrong_outside += 1,
            (false, true, _) => tally.missing_adjudication[cover] += 1,
            (false, false, true) => tally.missing_not_found[cover] += 1,
            (false, false, false) => tally.missing_outside[cover] += 1,
        }
    }

    for (label, t) in ["D remnant ≥ 8 nt", "D remnant 1–7 nt"]
        .iter()
        .zip(&tallies)
    {
        let pct = |x: usize| 100.0 * x as f64 / t.reads.max(1) as f64;
        println!(
            "{label}: {} reads, correct {:.2} %",
            t.reads,
            pct(t.correct)
        );
        println!(
            "  wrong D:   adjudication {:.2} %  not found {:.2} %  outside windows {:.2} %",
            pct(t.wrong_adjudication),
            pct(t.wrong_not_found),
            pct(t.wrong_outside)
        );
        for (name, counts) in [
            ("adjudication", &t.missing_adjudication),
            ("not found", &t.missing_not_found),
            ("outside windows", &t.missing_outside),
        ] {
            println!(
                "  missing D, {name:<16} {}",
                COVER
                    .iter()
                    .zip(counts.iter())
                    .map(|(c, &n)| format!("{c} {:.2} %", pct(n)))
                    .collect::<Vec<_>>()
                    .join("  ")
            );
        }
        println!(
            "  called J covers ≥ 3 bases of the true D: {:.2} %; called V: {:.2} %",
            pct(t.j_covers_d),
            pct(t.v_covers_d)
        );
    }
    Ok(())
}

fn p_names(set: &GermlineSet, allele: u32) -> impl Iterator<Item = &str> {
    set.allele(allele).names()
}
