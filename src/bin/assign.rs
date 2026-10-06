//! `vdj-assign` — assign V(D)J germline alleles to FASTA or FASTQ reads,
//! AIRR TSV on stdout; `vdj-assign prepare` writes an annotated database;
//! `simulate` and `evaluate` score any tool on identical simulated reads.
//!
//! ```text
//! vdj-assign --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta [options] reads.fastq > calls.tsv
//! vdj-assign prepare --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --out prepared/
//! vdj-assign simulate --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --preset hard --out sim
//! vdj-assign evaluate --truth sim.truth.tsv calls.tsv
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
use vdj_assign::sim::{Assignment, Evaluation, Rng, SimParams, Simulator, Truth};

const USAGE: &str = "\
usage: vdj-assign --v V.fasta --j J.fasta [--d D.fasta] [options] <reads.fasta|reads.fastq|->
       vdj-assign prepare --v V.fasta --j J.fasta [--d D.fasta] [--chain heavy|light] --out DIR
       vdj-assign simulate --v V.fasta --j J.fasta [--d D.fasta] --out PREFIX
                  [--preset easy|mid|hard] [--reads N] [--seed S] [--rc P] [--shm RATE]
       vdj-assign evaluate --truth PREFIX.truth.tsv <calls.tsv>

simulate writes PREFIX.fasta and PREFIX.truth.tsv: reads from the IgSim model
(presets as in IgSim.jl; --rc reverse-complements that fraction, default 0.5;
--shm fixes the SHM rate of every read) with their true alleles, the
germline-equivalent alleles, segment intervals and junctions. evaluate scores
an AIRR Rearrangement TSV from any tool against that truth (germline-
equivalence-aware allele, gene, boundary, D-presence and junction accuracy).

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
  --no-posterior       skip posterior probabilities (*_posterior columns; ~30 % faster)
  --residual           ordinary commit-V/J-then-D adjudication instead of the DAG
  -h, --help           this message";

/// Reads assigned per streaming batch (bounds memory on large files).
const BATCH: usize = 16_384;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Mode {
    #[default]
    Assign,
    Prepare,
    Simulate,
    Evaluate,
}

impl Mode {
    fn name(self) -> &'static str {
        match self {
            Mode::Assign => "assignment",
            Mode::Prepare => "prepare",
            Mode::Simulate => "simulate",
            Mode::Evaluate => "evaluate",
        }
    }
}

#[derive(Default)]
struct Args {
    mode: Mode,
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
    no_posterior: bool,
    residual: bool,
    preset: Option<String>,
    count: Option<String>,
    seed: Option<String>,
    rc: Option<String>,
    shm: Option<String>,
    truth: Option<String>,
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
    args.mode = match input.peek().map(String::as_str) {
        Some("prepare") => Mode::Prepare,
        Some("simulate") => Mode::Simulate,
        Some("evaluate") => Mode::Evaluate,
        _ => Mode::Assign,
    };
    if args.mode != Mode::Assign {
        input.next();
    }
    let mode = args.mode;
    let (assign, prepare, simulate, evaluate) = (
        mode == Mode::Assign,
        mode == Mode::Prepare,
        mode == Mode::Simulate,
        mode == Mode::Evaluate,
    );
    while let Some(arg) = input.next() {
        let mut value = |name: &str| input.next().ok_or_else(|| format!("{name} needs a value"));
        match arg.as_str() {
            "-h" | "--help" => return Ok(None),
            "--v" if !evaluate => set_once(&mut args.v, "--v", value("--v")?)?,
            "--d" if !evaluate => set_once(&mut args.d, "--d", value("--d")?)?,
            "--j" if !evaluate => set_once(&mut args.j, "--j", value("--j")?)?,
            "--out" if prepare || simulate => set_once(&mut args.out, "--out", value("--out")?)?,
            "--chain" if !evaluate => set_once(&mut args.chain, "--chain", value("--chain")?)?,
            "--profile" if assign => set_once(&mut args.profile, "--profile", value("--profile")?)?,
            "--threads" if assign => set_once(&mut args.threads, "--threads", value("--threads")?)?,
            "--tandem-d" if assign => args.tandem_d = true,
            "--forward-only" if assign => args.forward_only = true,
            "--no-anchor-window" if assign => args.no_anchor_window = true,
            "--no-posterior" if assign => args.no_posterior = true,
            "--residual" if assign => args.residual = true,
            "--preset" if simulate => set_once(&mut args.preset, "--preset", value("--preset")?)?,
            "--reads" if simulate => set_once(&mut args.count, "--reads", value("--reads")?)?,
            "--seed" if simulate => set_once(&mut args.seed, "--seed", value("--seed")?)?,
            "--rc" if simulate => set_once(&mut args.rc, "--rc", value("--rc")?)?,
            "--shm" if simulate => set_once(&mut args.shm, "--shm", value("--shm")?)?,
            "--truth" if evaluate => set_once(&mut args.truth, "--truth", value("--truth")?)?,
            flag if flag.starts_with("--") => {
                return Err(format!("unknown option '{flag}' for {}", mode.name()));
            }
            path if assign || evaluate => {
                set_once(&mut args.reads, "the reads or calls file", path.to_string())?
            }
            extra => {
                return Err(format!(
                    "{} takes no positional argument ('{extra}')",
                    mode.name()
                ))
            }
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
    // The residual adjudicator computes no posteriors.
    config.posterior = !args.no_posterior && !args.residual;
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

/// Parse an optional numeric option.
fn number<T: std::str::FromStr>(
    value: &Option<String>,
    option: &str,
    default: T,
) -> Result<T, String> {
    match value {
        Some(text) => text
            .parse()
            .map_err(|_| format!("{option} expects a number, got '{text}'")),
        None => Ok(default),
    }
}

fn run_simulate(args: &Args) -> Result<(), Box<dyn Error>> {
    let out = args.out.as_deref().ok_or("simulate needs --out PREFIX")?;
    let preset = args.preset.as_deref().unwrap_or("mid");
    let mut params = SimParams::preset(preset)
        .ok_or_else(|| format!("unknown preset '{preset}' (easy, mid, train, hard)"))?
        .with_reverse_complement(number(&args.rc, "--rc", 0.5)?);
    if args.shm.is_some() {
        params = params.with_shm(number(&args.shm, "--shm", 0.0)?);
    }
    let count: usize = number(&args.count, "--reads", 10_000)?;
    let seed: u64 = number(&args.seed, "--seed", 1)?;
    let database = load_database(args)?;
    let sim = Simulator::new(&database.v, database.d.as_ref(), &database.j, params)?;
    let mut reads = BufWriter::new(File::create(format!("{out}.fasta"))?);
    let mut truth = BufWriter::new(File::create(format!("{out}.truth.tsv"))?);
    writeln!(truth, "{}", Truth::tsv_header())?;
    let mut rng = Rng::new(seed);
    for index in 0..count {
        let read = sim.read(&mut rng)?;
        let id = format!("sim{index}");
        writeln!(reads, ">{id}\n{}", String::from_utf8_lossy(&read.sequence))?;
        writeln!(truth, "{}", Truth::from_read(&sim, &read, &id).to_tsv())?;
    }
    reads.flush()?;
    truth.flush()?;
    eprintln!("vdj-assign: wrote {count} reads to {out}.fasta and their truth to {out}.truth.tsv (preset {preset}, seed {seed})");
    Ok(())
}

fn run_evaluate(args: &Args) -> Result<(), Box<dyn Error>> {
    let truth_path = args
        .truth
        .as_deref()
        .ok_or("evaluate needs --truth PREFIX.truth.tsv")?;
    let calls_path = args
        .reads
        .as_deref()
        .ok_or("evaluate needs an AIRR TSV of calls")?;
    let truth_text = fs::read_to_string(truth_path).map_err(|e| format!("{truth_path}: {e}"))?;
    let calls_text = fs::read_to_string(calls_path).map_err(|e| format!("{calls_path}: {e}"))?;
    let mut lines = calls_text.lines();
    let header: Vec<&str> = lines
        .next()
        .ok_or("empty calls file")?
        .split('\t')
        .collect();
    let mut calls = std::collections::HashMap::new();
    for line in lines.filter(|l| !l.is_empty()) {
        let (id, assignment) = Assignment::from_airr(&header, line)?;
        calls.insert(id, assignment);
    }
    let mut lines = truth_text.lines();
    let header: Vec<&str> = lines
        .next()
        .ok_or("empty truth file")?
        .split('\t')
        .collect();
    let mut evaluation = Evaluation::default();
    let mut missing = 0;
    for line in lines.filter(|l| !l.is_empty()) {
        let truth = Truth::from_tsv(&header, line).map_err(|e| format!("{truth_path}: {e}"))?;
        let call = calls.get(&truth.id).cloned().unwrap_or_else(|| {
            missing += 1;
            Assignment::default()
        });
        evaluation.add(&truth, &call);
    }
    if missing > 0 {
        eprintln!("vdj-assign: {missing} simulated reads have no row in {calls_path}; counted as unassigned");
    }
    println!("{evaluation}");
    Ok(())
}

fn main() -> ExitCode {
    match parse_args() {
        Ok(None) => {
            println!("{USAGE}");
            ExitCode::SUCCESS
        }
        Ok(Some(args)) => {
            let result = match args.mode {
                Mode::Assign => run_assign(&args),
                Mode::Prepare => run_prepare(&args),
                Mode::Simulate => run_simulate(&args),
                Mode::Evaluate => run_evaluate(&args),
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
