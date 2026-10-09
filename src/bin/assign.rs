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

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use vdj_assign::prelude::*;
use vdj_assign::prepare::DatabaseError;
use vdj_assign::sim::{Assignment, C_STREAM, Evaluation, Paired, Rng, SimParams, Simulator, Truth};

const USAGE: &str = "\
usage: vdj-assign --v V.fasta --j J.fasta [--d D.fasta] [--c C.fasta] [options] <reads.fasta|reads.fastq|->
       vdj-assign prepare --v V.fasta --j J.fasta [--d D.fasta] [--c C.fasta] [--chain heavy|light] --out DIR
       vdj-assign simulate --v V.fasta --j J.fasta [--d D.fasta] [--c C.fasta] --out PREFIX
                  [--preset easy|mid|hard] [--reads N] [--seed S] [--rc P] [--shm RATE]
       vdj-assign evaluate --truth PREFIX.truth.tsv <calls.tsv> [other.tsv]

simulate writes PREFIX.fasta and PREFIX.truth.tsv: reads from the IgSim model
(presets as in IgSim.jl; --rc reverse-complements that fraction, default 0.5;
--shm fixes the SHM rate of every read) with their true alleles, the
germline-equivalent alleles, segment intervals and junctions. evaluate scores
an AIRR Rearrangement TSV from any tool against that truth (germline-
equivalence-aware allele, gene, boundary, D-presence and junction accuracy).
Given a second calls file, evaluate also compares the two read by read:
accuracy difference with its 95 % confidence interval and an exact McNemar
test on the reads only one of them gets right (D also by remnant length),
and the D allele and D-presence Brier scores with a paired t.

Reads may be FASTA or FASTQ (detected from the first character); '-' reads
stdin, so compressed input can be streamed: zcat reads.fastq.gz | vdj-assign ... -

Germline databases are FASTA. Headers may carry 0-based positions:
  anchor=POS            first base of the V Cys104 / J Trp-Phe118 codon (fixes the frame)
  fr1= cdr1= fr2= cdr2= fr3=   V region starts (IMGT definitions)
Missing values are detected on every start (motifs, database consensus,
transfer from the closest sibling allele; milliseconds). Identical sequences
are merged and reported under all their names. `prepare` writes the result as
annotated FASTA plus annotation.tsv for review; header values always win.

--c adds constant-region (isotype) calls: C base 0 is anchored within 3 bases
of the called J's end and aligned to the read's end, and the c_ columns follow
productive. Give one locus' C regions (coding CH1 onward, e.g. human IGHC):
the gate counts every C gene in the file. prepare warns about C records that
cannot be anchored (5'-truncated, genomic) and J records whose 3' end strays.
simulate --c reads 80 % of reads through the J-C splice into a C allele and
adds c_call, c_equivalent, c_start and c_end to the truth; evaluate then
scores C by its length against anchored equivalence.

options:
  --profile NAME       r (default) | truth | igblast | sensitive-d
                       r is Swig's r_optimized; also accepts r_optimized, truth_optimized, sensitive_d
  --threads N          worker threads (default: all cores)
  --chain TYPE         heavy | light (default: voted from the database)
  --tandem-d           allow VDDJ paths (second D as a first-class node)
  --forward-only       do not try the reverse complement
  --no-anchor-window   do not add the conserved-anchor D window
  --no-posterior       skip posterior probabilities (*_posterior columns; ~30 % faster)
  --calls MODEL        the model the call maximises: per-read (each read's own, the
                       model its posteriors come from) | fixed (the (2, -3) scores)
  --d-allele RULE      the reported D allele: map (the call's) | max-marginal (largest
                       posterior marginal given the call's D presence and topology;
                       needs posteriors, so --no-posterior and --residual report map)
  --d-presence-rate R  D prior from the D search's own chance-hit count (alleles x
                       offsets x positions) at D-presence rate R, 0 < R < 1, instead
                       of genes x positions at even odds (docs/ALGORITHM.md 2.3)
  --c-evalue E         call C when its anchored E-value, C genes x 7 start
                       diagonals x e^(-lambda S), is at most E (default 1e-3)
  --simd NAME          alignment kernel: avx512 | avx2 | neon | portable | scalar
                       (default: the fastest this CPU supports, avx512 excepted)
  --residual           ordinary commit-V/J-then-D adjudication instead of the DAG
  -h, --help           this message";

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

    /// Whether option `name` of this mode takes a value; `None` when the
    /// mode has no such option. Assignment takes every [`RunOptions`] name.
    fn takes_value(self, name: &str) -> Option<bool> {
        let germline = matches!(name, "v" | "d" | "j" | "c");
        let own = match self {
            Mode::Assign => germline,
            Mode::Prepare => germline || matches!(name, "out" | "chain"),
            Mode::Simulate => {
                germline
                    || matches!(
                        name,
                        "out" | "chain" | "preset" | "reads" | "seed" | "rc" | "shm"
                    )
            }
            Mode::Evaluate => name == "truth",
        };
        match self {
            _ if own => Some(true),
            Mode::Assign => RunOptions::takes_value(name),
            _ => None,
        }
    }
}

/// Germline FASTA paths.
#[derive(Default)]
struct Files {
    v: Option<String>,
    d: Option<String>,
    j: Option<String>,
    c: Option<String>,
}

impl Files {
    /// The option and path of the `kind` database.
    fn of(&self, kind: SegmentKind) -> (&'static str, &str) {
        let (option, path) = match kind {
            SegmentKind::V => ("--v", &self.v),
            SegmentKind::D => ("--d", &self.d),
            SegmentKind::J => ("--j", &self.j),
            SegmentKind::C => ("--c", &self.c),
            kind => unreachable!("no {} database is loaded", kind.label()),
        };
        (option, path.as_deref().unwrap_or_default())
    }
}

#[derive(Default)]
struct Args {
    mode: Mode,
    files: Files,
    out: Option<String>,
    /// Assignment: the reads; evaluate: the calls.
    input: Option<String>,
    /// Evaluate: a second calls file to compare against the first.
    other: Option<String>,
    truth: Option<String>,
    run: RunOptions,
    preset: Option<String>,
    count: Option<String>,
    seed: Option<String>,
    rc: Option<String>,
    shm: Option<String>,
}

fn parse_args() -> Result<Option<Args>, String> {
    let mut input = std::env::args().skip(1).peekable();
    let mode = match input.peek().map(String::as_str) {
        Some("prepare") => Mode::Prepare,
        Some("simulate") => Mode::Simulate,
        Some("evaluate") => Mode::Evaluate,
        _ => Mode::Assign,
    };
    if mode != Mode::Assign {
        input.next();
    }
    let mut args = Args {
        mode,
        ..Args::default()
    };
    let mut seen = HashSet::new();
    while let Some(arg) = input.next() {
        if arg == "-h" || arg == "--help" {
            return Ok(None);
        }
        let Some(name) = arg.strip_prefix("--") else {
            let slot = match mode {
                Mode::Assign | Mode::Evaluate if args.input.is_none() => &mut args.input,
                Mode::Evaluate if args.other.is_none() => &mut args.other,
                Mode::Assign | Mode::Evaluate => {
                    return Err(format!("unexpected argument '{arg}' for {}", mode.name()));
                }
                Mode::Prepare | Mode::Simulate => {
                    return Err(format!(
                        "{} takes no positional argument ('{arg}')",
                        mode.name()
                    ));
                }
            };
            *slot = Some(arg);
            continue;
        };
        let takes_value = mode
            .takes_value(name)
            .ok_or_else(|| format!("unknown option '{arg}' for {}", mode.name()))?;
        if !seen.insert(name.to_string()) {
            return Err(format!("{arg} given more than once"));
        }
        let value = match takes_value {
            true => Some(input.next().ok_or_else(|| format!("{arg} needs a value"))?),
            false => None,
        };
        let slot = match name {
            "v" => &mut args.files.v,
            "d" => &mut args.files.d,
            "j" => &mut args.files.j,
            "c" => &mut args.files.c,
            "out" => &mut args.out,
            "truth" => &mut args.truth,
            "preset" => &mut args.preset,
            "reads" => &mut args.count,
            "seed" => &mut args.seed,
            "rc" => &mut args.rc,
            "shm" => &mut args.shm,
            _ => {
                args.run
                    .set(name, value.as_deref())
                    .map_err(|e| e.to_string())?;
                continue;
            }
        };
        *slot = value;
    }
    Ok(Some(args))
}

/// Read the germline FASTA files and prepare the database, logging its
/// summary and warnings.
fn load_database(files: &Files, options: &RunOptions) -> Result<PreparedDatabase, Box<dyn Error>> {
    let read = |path: &str| fs::read_to_string(path).map_err(|e| format!("{path}: {e}"));
    let v = read(files.v.as_deref().ok_or("--v is required")?)?;
    let j = read(files.j.as_deref().ok_or("--j is required")?)?;
    let d = files.d.as_deref().map(read).transpose()?;
    let c = files.c.as_deref().map(read).transpose()?;
    let texts = GermlineTexts {
        v: &v,
        d: d.as_deref(),
        j: &j,
        c: c.as_deref(),
    };
    let started = Instant::now();
    let database = PreparedDatabase::from_fasta(texts, options.annotate()).map_err(|e| {
        let (option, path) = files.of(e.kind());
        match e {
            DatabaseError::Fastq(_) => format!(
                "{path}: looks like FASTQ, but {option} expects a germline FASTA database (reads go last, as a positional argument)"
            ),
            e => format!("{path}: {e}"),
        }
    })?;
    eprintln!(
        "vdj-assign: database prepared in {:.1} ms: {}",
        started.elapsed().as_secs_f64() * 1e3,
        database.summary()
    );
    for warning in database.warnings(CSegmentConfig::default().spread) {
        eprintln!("vdj-assign: warning: {warning}");
    }
    Ok(database)
}

fn run_prepare(args: &Args) -> Result<(), Box<dyn Error>> {
    let out = args.out.as_deref().ok_or("prepare needs --out DIR")?;
    let database = load_database(&args.files, &args.run)?;
    database.write(Path::new(out))?;
    let names: Vec<&str> = database.files().iter().map(|&(name, _)| name).collect();
    eprintln!("vdj-assign: wrote {} into {out}", names.join(", "));
    Ok(())
}

fn run_assign(args: &Args) -> Result<(), Box<dyn Error>> {
    if cfg!(debug_assertions) {
        eprintln!(
            "vdj-assign: warning: debug build, typically 10-50x slower; build with `cargo build --release`"
        );
    }
    let reads = args
        .input
        .as_deref()
        .ok_or("a reads file (FASTA, FASTQ or '-') is required")?;
    let database = load_database(&args.files, &args.run)?;
    let run = Run::new(&args.run, database)?;
    let started = Instant::now();
    let mut out = BufWriter::new(io::stdout().lock());
    let (source, written) = match reads {
        "-" => ("stdin", run.write_airr(io::stdin().lock(), &mut out)),
        path => {
            let file = File::open(path).map_err(|e| format!("{path}: {e}"))?;
            let input = BufReader::with_capacity(1 << 20, file);
            (path, run.write_airr(input, &mut out))
        }
    };
    let tally = written.map_err(|e| match e {
        RunError::Reads(e) => format!("{source}: {e}"),
        e => e.to_string(),
    })?;
    out.flush()?;
    let seconds = started.elapsed().as_secs_f64();
    eprintln!(
        "vdj-assign: {} reads ({} with a V(D)J path) in {seconds:.1} s ({:.0} reads/s)",
        tally.reads,
        tally.called,
        tally.reads as f64 / seconds.max(1e-9)
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
    let database = load_database(&args.files, &args.run)?;
    let mut sim = Simulator::new(&database.v, database.d.as_ref(), &database.j, params)?;
    if let Some(c) = &database.c {
        sim = sim.with_constant(c)?;
    }
    let mut reads = BufWriter::new(File::create(format!("{out}.fasta"))?);
    let mut truth = BufWriter::new(File::create(format!("{out}.truth.tsv"))?);
    writeln!(truth, "{}", Truth::tsv_header_with(database.c.is_some()))?;
    let (mut rng, mut c_rng) = (Rng::new(seed), Rng::new(seed ^ C_STREAM));
    for index in 0..count {
        let mut read = sim.read(&mut rng)?;
        sim.append_constant(&mut read, &mut c_rng);
        let id = format!("sim{index}");
        writeln!(reads, ">{id}\n{}", String::from_utf8_lossy(&read.sequence))?;
        writeln!(truth, "{}", Truth::from_read(&sim, &read, &id).to_tsv())?;
    }
    reads.flush()?;
    truth.flush()?;
    eprintln!(
        "vdj-assign: wrote {count} reads to {out}.fasta and their truth to {out}.truth.tsv (preset {preset}, seed {seed})"
    );
    Ok(())
}

/// The rows of an AIRR TSV of calls, by `sequence_id`.
fn read_calls(path: &str) -> Result<HashMap<String, Assignment>, Box<dyn Error>> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines
        .next()
        .ok_or_else(|| format!("{path}: empty calls file"))?
        .split('\t')
        .collect();
    lines
        .filter(|l| !l.is_empty())
        .map(|line| Ok(Assignment::from_airr(&header, line).map_err(|e| format!("{path}: {e}"))?))
        .collect()
}

fn run_evaluate(args: &Args) -> Result<(), Box<dyn Error>> {
    let truth_path = args
        .truth
        .as_deref()
        .ok_or("evaluate needs --truth PREFIX.truth.tsv")?;
    let calls_path = args
        .input
        .as_deref()
        .ok_or("evaluate needs an AIRR TSV of calls")?;
    let calls = read_calls(calls_path)?;
    let other = args.other.as_deref().map(read_calls).transpose()?;
    let truth_text = fs::read_to_string(truth_path).map_err(|e| format!("{truth_path}: {e}"))?;
    let mut lines = truth_text.lines();
    let header: Vec<&str> = lines
        .next()
        .ok_or("empty truth file")?
        .split('\t')
        .collect();
    let mut evaluation = Evaluation::default();
    let mut paired = Paired::default();
    let mut missing = [0usize; 2];
    let unassigned = Assignment::default();
    for line in lines.filter(|l| !l.is_empty()) {
        let truth = Truth::from_tsv(&header, line).map_err(|e| format!("{truth_path}: {e}"))?;
        let lookup = |calls: &HashMap<String, Assignment>, missing: &mut usize| {
            calls.get(&truth.id).cloned().unwrap_or_else(|| {
                *missing += 1;
                unassigned.clone()
            })
        };
        let call = lookup(&calls, &mut missing[0]);
        evaluation.add(&truth, &call);
        if let Some(other) = &other {
            paired.add(&truth, &call, &lookup(other, &mut missing[1]));
        }
    }
    for (path, missing) in [Some(calls_path), args.other.as_deref()]
        .into_iter()
        .zip(missing)
    {
        if let (Some(path), 1..) = (path, missing) {
            eprintln!(
                "vdj-assign: {missing} simulated reads have no row in {path}; counted as unassigned"
            );
        }
    }
    println!("{evaluation}");
    if let Some(other_path) = &args.other {
        println!("\nfirst: {calls_path}\nsecond: {other_path}\n{paired}");
    }
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
