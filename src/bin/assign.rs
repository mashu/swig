//! `vdj-assign` — assign V(D)J germline alleles to FASTA or FASTQ reads,
//! AIRR TSV on stdout; `vdj-assign prepare` writes an annotated database.
//!
//! ```text
//! vdj-assign --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta [options] reads.fastq > calls.tsv
//! vdj-assign prepare --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --out prepared/
//! ```

use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use vdj_assign::annotate::{AnnotateOptions, ChainType};
use vdj_assign::prelude::*;
use vdj_assign::prepare::{prepare, PreparedDatabase};
use vdj_assign::reads::ReadReader;

const USAGE: &str = "\
usage: vdj-assign --v V.fasta --j J.fasta [--d D.fasta] [options] <reads.fasta|reads.fastq|->
       vdj-assign prepare --v V.fasta --j J.fasta [--d D.fasta] [--chain heavy|light] --out DIR

Reads may be FASTA or FASTQ (detected from the first character); '-' reads
stdin, so compressed input can be streamed: zcat reads.fastq.gz | vdj-assign ... -

Germline databases are FASTA. Headers may carry 0-based positions:
  anchor=POS            first base of the V Cys104 / J Trp-Phe118 codon (fixes the frame)
  fr1= cdr1= fr2= cdr2= fr3=   V region starts (IMGT definitions)
Missing values are detected on every start (motifs, database consensus,
transfer from the closest sibling allele; milliseconds). Identical sequences
are merged and reported under all their names. `prepare` writes the result as
annotated FASTA plus annotation.tsv for review; header values always win.

options:
  --profile NAME       r (default) | truth | igblast | sensitive-d
                       r is Swig's r_optimized; also accepts r_optimized, truth_optimized, sensitive_d
  --threads N          worker threads (default: all cores)
  --chain TYPE         heavy | light (default: voted from the database)
  --tandem-d           allow VDDJ paths (second D as a first-class node)
  --forward-only       do not try the reverse complement
  --no-anchor-window   do not add the conserved-anchor D window
  --residual           ordinary commit-V/J-then-D adjudication instead of the DAG
  -h, --help           this message";

/// Reads assigned per streaming batch (bounds memory on large files).
const BATCH: usize = 16_384;

#[derive(Default)]
struct Args {
    prepare: bool,
    v: Option<String>,
    d: Option<String>,
    j: Option<String>,
    out: Option<String>,
    reads: Option<String>,
    profile: Option<String>,
    threads: Option<String>,
    chain: Option<String>,
    tandem_d: bool,
    forward_only: bool,
    no_anchor_window: bool,
    residual: bool,
}

/// Store `value` in `slot`, rejecting a second occurrence of the option.
fn set_once(slot: &mut Option<String>, option: &str, value: String) -> Result<(), String> {
    match slot.replace(value) {
        Some(previous) => Err(format!(
            "{option} given more than once ('{previous}' and '{}')",
            slot.as_deref().unwrap_or_default()
        )),
        None => Ok(()),
    }
}

fn parse_args() -> Result<Option<Args>, String> {
    let mut args = Args::default();
    let mut input = std::env::args().skip(1).peekable();
    if input.peek().map(String::as_str) == Some("prepare") {
        input.next();
        args.prepare = true;
    }
    while let Some(arg) = input.next() {
        let mut value = |name: &str| input.next().ok_or_else(|| format!("{name} needs a value"));
        match arg.as_str() {
            "-h" | "--help" => return Ok(None),
            "--v" => set_once(&mut args.v, "--v", value("--v")?)?,
            "--d" => set_once(&mut args.d, "--d", value("--d")?)?,
            "--j" => set_once(&mut args.j, "--j", value("--j")?)?,
            "--out" if args.prepare => set_once(&mut args.out, "--out", value("--out")?)?,
            "--chain" => set_once(&mut args.chain, "--chain", value("--chain")?)?,
            "--profile" if !args.prepare => {
                set_once(&mut args.profile, "--profile", value("--profile")?)?
            }
            "--threads" if !args.prepare => {
                set_once(&mut args.threads, "--threads", value("--threads")?)?
            }
            "--tandem-d" if !args.prepare => args.tandem_d = true,
            "--forward-only" if !args.prepare => args.forward_only = true,
            "--no-anchor-window" if !args.prepare => args.no_anchor_window = true,
            "--residual" if !args.prepare => args.residual = true,
            flag if flag.starts_with("--") => {
                let mode = if args.prepare {
                    "prepare"
                } else {
                    "assignment"
                };
                return Err(format!("unknown option '{flag}' for {mode}"));
            }
            path if !args.prepare => set_once(&mut args.reads, "the reads file", path.to_string())?,
            extra => return Err(format!("prepare takes no positional argument ('{extra}')")),
        }
    }
    Ok(Some(args))
}

fn load_set(kind: SegmentKind, option: &str, path: &str) -> Result<GermlineSet, Box<dyn Error>> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    if text.trim_start().starts_with('@') {
        return Err(format!(
            "{path}: looks like FASTQ, but {option} expects a germline FASTA database (reads go last, as a positional argument)"
        )
        .into());
    }
    Ok(GermlineSet::from_fasta(kind, &text).map_err(|e| format!("{path}: {e}"))?)
}

fn load_database(args: &Args) -> Result<PreparedDatabase, Box<dyn Error>> {
    let chain = match args.chain.as_deref() {
        Some(name) => Some(
            ChainType::from_name(name)
                .ok_or_else(|| format!("--chain expects heavy or light, got '{name}'"))?,
        ),
        None => None,
    };
    let v = load_set(
        SegmentKind::V,
        "--v",
        args.v.as_deref().ok_or("--v is required")?,
    )?;
    let d = args
        .d
        .as_deref()
        .map(|path| load_set(SegmentKind::D, "--d", path))
        .transpose()?;
    let j = load_set(
        SegmentKind::J,
        "--j",
        args.j.as_deref().ok_or("--j is required")?,
    )?;
    let started = Instant::now();
    let database = prepare(v, d, j, AnnotateOptions { chain })?;
    eprintln!(
        "vdj-assign: database prepared in {:.1} ms: {}",
        started.elapsed().as_secs_f64() * 1e3,
        database.summary()
    );
    if database.annotation.uncertain() > 0 {
        eprintln!(
            "vdj-assign: {} alleles have transferred, low-confidence or missing anchors; see `vdj-assign prepare` annotation.tsv",
            database.annotation.uncertain()
        );
    }
    Ok(database)
}

/// Stream reads in batches through any assigner (static dispatch per variant).
fn stream<R: BufRead, T: VdjAssigner + AlleleCatalog + Sync>(
    input: R,
    source: &str,
    assigner: &T,
    threads: usize,
    support: &SupportModels,
) -> Result<(usize, usize), Box<dyn Error>> {
    let mut reads = ReadReader::new(input);
    let mut out = BufWriter::new(io::stdout().lock());
    writeln!(out, "{}", AirrRecord::tsv_header())?;
    let (mut total, mut with_path) = (0, 0);
    let mut batch = Vec::with_capacity(BATCH);
    loop {
        batch.clear();
        if reads
            .read_batch(BATCH, &mut batch)
            .map_err(|e| format!("{source}: {e}"))?
            == 0
        {
            break;
        }
        let encoded: Vec<Vec<u8>> = batch.iter().map(|r| encode(&r.sequence)).collect();
        let rows = assign_map(assigner, &encoded, threads, |i, call| {
            let sequence = String::from_utf8_lossy(&batch[i].sequence);
            let has_path = call.path.is_some();
            (
                AirrRecord::from_call(batch[i].id(), &sequence, &call, assigner, support).to_tsv(),
                has_path,
            )
        });
        for (row, has_path) in rows {
            writeln!(out, "{row}")?;
            with_path += usize::from(has_path);
        }
        total += batch.len();
    }
    out.flush()?;
    Ok((total, with_path))
}

fn stream_input<T: VdjAssigner + AlleleCatalog + Sync>(
    path: &str,
    assigner: &T,
    threads: usize,
    support: &SupportModels,
) -> Result<(usize, usize), Box<dyn Error>> {
    match path {
        "-" => stream(io::stdin().lock(), "stdin", assigner, threads, support),
        _ => {
            let file = File::open(path).map_err(|e| format!("{path}: {e}"))?;
            stream(
                BufReader::with_capacity(1 << 20, file),
                path,
                assigner,
                threads,
                support,
            )
        }
    }
}

fn run_prepare(args: &Args) -> Result<(), Box<dyn Error>> {
    let out = args.out.as_deref().ok_or("prepare needs --out DIR")?;
    let database = load_database(args)?;
    database.write(Path::new(out))?;
    eprintln!(
        "vdj-assign: wrote {out}/V.fasta, J.fasta{}, annotation.tsv",
        if database.d.is_some() {
            ", D.fasta"
        } else {
            ""
        }
    );
    Ok(())
}

fn run_assign(args: &Args) -> Result<(), Box<dyn Error>> {
    if cfg!(debug_assertions) {
        eprintln!("vdj-assign: warning: debug build, typically 10-50x slower; build with `cargo build --release`");
    }
    let reads = args
        .reads
        .as_deref()
        .ok_or("a reads file (FASTA, FASTQ or '-') is required")?;
    let profile = match args.profile.as_deref() {
        Some(name) => {
            CallingProfile::from_name(name).ok_or_else(|| format!("unknown profile '{name}'"))?
        }
        None => CallingProfile::ROptimized,
    };
    let threads = match args.threads.as_deref() {
        Some(n) => n
            .parse()
            .map_err(|_| format!("--threads expects a number, got '{n}'"))?,
        None => 0,
    };
    let mut config = profile.config();
    if args.tandem_d {
        config.grammar = VdjGrammar::TANDEM_D;
    }
    if args.forward_only {
        config.options.orientation = OrientationPolicy::FORWARD_ONLY;
    }
    config.anchor_window = !args.no_anchor_window;
    let support: SupportModels = [
        KarlinAltschul::ungapped(&config.v.scoring),
        KarlinAltschul::ungapped(&config.d.scoring),
        KarlinAltschul::ungapped(&config.j.scoring),
    ];
    let database = load_database(args)?;
    let assigner = DefaultAssigner::from_config(&config, database.v, database.d, database.j)?;

    let started = Instant::now();
    let (total, with_path) = if args.residual {
        let residual = ResidualAdjudicator {
            min_retained: config.min_retained,
            ..ResidualAdjudicator::default()
        };
        stream_input(
            reads,
            &assigner.with_adjudicator(residual),
            threads,
            &support,
        )?
    } else {
        stream_input(reads, &assigner, threads, &support)?
    };
    let seconds = started.elapsed().as_secs_f64();
    eprintln!(
        "vdj-assign: {total} reads ({with_path} with a V(D)J path) in {seconds:.1} s ({:.0} reads/s)",
        total as f64 / seconds.max(1e-9)
    );
    Ok(())
}

fn main() -> ExitCode {
    match parse_args() {
        Ok(None) => {
            println!("{USAGE}");
            ExitCode::SUCCESS
        }
        Ok(Some(args)) => {
            let result = if args.prepare {
                run_prepare(&args)
            } else {
                run_assign(&args)
            };
            match result {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("vdj-assign: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Err(e) => {
            eprintln!("vdj-assign: {e}\n\n{USAGE}");
            ExitCode::from(2)
        }
    }
}
