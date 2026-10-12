//! `swig` — assign V(D)J germline alleles to FASTA or FASTQ reads,
//! AIRR TSV on stdout; `swig prepare` writes an annotated database;
//! `simulate` and `evaluate` score any tool on identical simulated reads.
//!
//! ```text
//! swig --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta [options] reads.fastq > calls.tsv
//! swig prepare --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --out prepared/
//! swig simulate --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --preset hard --out sim
//! swig evaluate --truth sim.truth.tsv calls.tsv
//! ```

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use swig::double_d::DoubleDMode;
use swig::prelude::*;
use swig::prepare::DatabaseError;
use swig::sim::{Assignment, C_STREAM, Evaluation, Paired, Rng, SimParams, Simulator, Truth};

const USAGE: &str = "\
usage: swig --v V.fasta --j J.fasta [--d D.fasta] [--c C.fasta] [options] <reads|->
       swig prepare --v V.fasta --j J.fasta [--d D.fasta] [--c C.fasta] [--chain heavy|light] --out DIR
       swig simulate --v V.fasta --j J.fasta [--d D.fasta] [--c C.fasta] [--chain heavy|light]
                  --out PREFIX [--preset easy|mid|train|hard] [--reads N] [--seed S]
                  [--rc P] [--shm RATE] [--dd P]
       swig evaluate --truth PREFIX.truth.tsv <calls.tsv> [other.tsv]

Assigns V(D)J (and C) germline alleles to reads and writes Swig's AIRR
Rearrangement TSV on stdout. Reads may be FASTA, FASTQ or an AIRR table
with a 'sequence' column (detected from the first character; TSV, or CSV
without quoting); '-' reads stdin, so compressed input can be streamed:
zcat reads.fastq.gz | swig ... -

simulate writes PREFIX.fasta and PREFIX.truth.tsv: reads from the IgSim model
(presets as in IgSim.jl; --rc reverse-complements that fraction, default 0.5;
--shm fixes the SHM rate of every read; --dd sets the tandem-DD fraction) with
their true alleles, the germline-equivalent alleles, segment intervals and
junctions; each read's V, D and C come from its J's locus, and V and J records
without a trusted anchor and the C records prepare warns about are not drawn.
evaluate scores an AIRR Rearrangement TSV from any tool against that truth
(germline-equivalence-aware allele, gene, boundary, D-presence and junction
accuracy). Given a second calls file, evaluate also compares the two read by
read: accuracy difference with its 95 % confidence interval and an exact
McNemar test on the reads only one of them gets right (D also by remnant
length), and the D allele and D-presence Brier scores with a paired t.

Germline databases are FASTA and may hold several loci (IGH, IGK and IGL;
or the TR loci): each allele's locus is read from its name, and every
rearrangement, D and C call stays within one locus. Headers may carry
0-based positions:
  anchor=POS            first base of the V Cys104 / J Trp-Phe118 codon (fixes the frame)
  fr1= cdr1= fr2= cdr2= fr3=   V region starts (IMGT definitions)
  frame=F               offset of the first complete codon (D: d_frame)
  SWIGMETA=...          Swig's prepared coordinates (filling what the keys leave open)
Missing values are detected on every start, per locus (motifs, database
consensus, transfer from the closest sibling allele; milliseconds).
Identical sequences are merged and reported under all their names.
`prepare` writes the result as annotated FASTA plus annotation.tsv for
review; header values always win.

--c adds constant-region (isotype) calls: C base 0 is anchored within 3 bases
of the called J's end, in the J's locus, and aligned to the read's end. The
gate counts the C genes of that locus. prepare warns about C records that
cannot be anchored (5'-truncated, genomic) and J records whose 3' end strays.
simulate --c reads 80 % of reads through the J-C splice into a C allele and
adds c_call, c_equivalent, c_start and c_end to the truth; evaluate then
scores C by its length against anchored equivalence.

options:
  --profile NAME       r (default) | truth | igblast | sensitive-d
                       (also r_optimized, truth_optimized, igblast_compatible, sensitive_d)
  --threads N          worker threads (default: all cores)
  --chain TYPE         heavy | light (default: voted per locus from the database)
  --strand S           both (default) | forward | reverse
  --min-identity F     drop V, D and J hits below identity F (0-1; default none)
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
  --min-d-match N      shortest exact run of an admissible D (default: the profile's)
  --d-mismatch S       D mismatch score, negative (default: the profile's)
  --j-mismatch S       J mismatch score, negative (default: the profile's)
  --min-j-length N     fewest aligned J columns of a J hit (default none)
  --c-evalue E         call C when its anchored E-value, C genes x 7 start
                       diagonals x e^(-lambda S), is at most E (default 1e-3)
  --c-prefix-identity F  report C only where the germline from its first base
                       reaches identity F (c_prefix_identity; default 0, off)
  --double-d MODE      off (default) | all | long-span: Swig's double-D screen, a
                       seeded search for two D between every call's V and J,
                       written to --double-d-out
  --double-d-out FILE  the screen's sidecar TSV
  --double-d-min-span N, --double-d-seed-length N, --double-d-pseudo-trim N,
  --double-d-pseudo-mismatches N, --double-d-min-gain S
                       the screen's evidence rules (defaults 40, 11, 5, 3, 8)
  --simd NAME          alignment kernel: avx512 | avx2 | neon | portable | scalar
                       (default: the fastest this CPU supports, avx512 excepted)
  --residual           ordinary commit-V/J-then-D adjudication instead of the DAG
  -V, --version        the version
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
            Mode::Assign => germline || name == "double-d-out",
            Mode::Prepare => germline || matches!(name, "out" | "chain"),
            Mode::Simulate => {
                germline
                    || matches!(
                        name,
                        "out" | "chain" | "preset" | "reads" | "seed" | "rc" | "shm" | "dd"
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
    dd: Option<String>,
    /// Assignment: the double-D sidecar.
    double_d_out: Option<String>,
}

/// What the command line asks for.
enum Command {
    Help,
    Version,
    Run(Box<Args>),
}

fn parse_args() -> Result<Command, String> {
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
            return Ok(Command::Help);
        }
        if arg == "-V" || arg == "--version" {
            return Ok(Command::Version);
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
        // A repeated value is ambiguous; a repeated flag changes nothing.
        if takes_value && !seen.insert(name.to_string()) {
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
            "dd" => &mut args.dd,
            "double-d-out" => &mut args.double_d_out,
            _ => {
                args.run
                    .set(name, value.as_deref())
                    .map_err(|e| e.to_string())?;
                continue;
            }
        };
        *slot = value;
    }
    args.check()?;
    Ok(Command::Run(Box::new(args)))
}

impl Args {
    /// Refuse a command line that lacks an input its mode needs, or gives
    /// the double-D screen without its sidecar (or the reverse): a usage
    /// error like any other bad command line.
    fn check(&self) -> Result<(), String> {
        let mode = self.mode;
        let germline = matches!(mode, Mode::Assign | Mode::Prepare | Mode::Simulate);
        let missing = [
            (germline && self.files.v.is_none(), "--v V.fasta"),
            (germline && self.files.j.is_none(), "--j J.fasta"),
            (mode == Mode::Prepare && self.out.is_none(), "--out DIR"),
            (mode == Mode::Simulate && self.out.is_none(), "--out PREFIX"),
            (
                mode == Mode::Evaluate && self.truth.is_none(),
                "--truth PREFIX.truth.tsv",
            ),
            (
                mode == Mode::Evaluate && self.input.is_none(),
                "an AIRR TSV of calls",
            ),
            (
                mode == Mode::Assign && self.input.is_none(),
                "a reads file (FASTA, FASTQ, AIRR TSV or '-')",
            ),
        ];
        if let Some((_, what)) = missing.iter().find(|(missing, _)| *missing) {
            return Err(format!("{} needs {what}", mode.name()));
        }
        let screening = self.run.double_d.mode != DoubleDMode::Off;
        if screening != self.double_d_out.is_some() {
            return Err("--double-d and --double-d-out go together".into());
        }
        Ok(())
    }
}

/// Read the germline FASTA files and prepare the database, logging its
/// summary and warnings.
fn load_database(files: &Files, options: &RunOptions) -> Result<PreparedDatabase, Box<dyn Error>> {
    let read = |path: &str| fs::read_to_string(path).map_err(|e| format!("{path}: {e}"));
    let v = read(files.v.as_deref().expect("checked by Args::check"))?;
    let j = read(files.j.as_deref().expect("checked by Args::check"))?;
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
        "swig: database prepared in {:.1} ms: {}",
        started.elapsed().as_secs_f64() * 1e3,
        database.summary()
    );
    for warning in database.warnings(CSegmentConfig::default().spread) {
        eprintln!("swig: warning: {warning}");
    }
    Ok(database)
}

fn run_prepare(args: &Args) -> Result<(), Box<dyn Error>> {
    let out = args.out.as_deref().expect("checked by Args::check");
    let database = load_database(&args.files, &args.run)?;
    database.write(Path::new(out))?;
    let names: Vec<&str> = database.files().iter().map(|&(name, _)| name).collect();
    eprintln!("swig: wrote {} into {out}", names.join(", "));
    Ok(())
}

fn run_assign(args: &Args) -> Result<(), Box<dyn Error>> {
    if cfg!(debug_assertions) {
        eprintln!(
            "swig: warning: debug build, typically 10-50x slower; build with `cargo build --release`"
        );
    }
    let reads = args.input.as_deref().expect("checked by Args::check");
    let screening = args.run.double_d.mode != DoubleDMode::Off;
    let database = load_database(&args.files, &args.run)?;
    let mut run = Run::new(&args.run, database)?;
    let mut sidecar = args
        .double_d_out
        .as_deref()
        .map(|path| {
            File::create(path)
                .map(BufWriter::new)
                .map_err(|e| format!("{path}: {e}"))
        })
        .transpose()?;
    let started = Instant::now();
    let mut out = BufWriter::new(io::stdout().lock());
    let double_d = sidecar.as_mut().map(|s| s as &mut dyn Write);
    let (source, written) = match reads {
        "-" => (
            "stdin",
            run.write_airr(io::stdin().lock(), &mut out, double_d),
        ),
        path => {
            let file = File::open(path).map_err(|e| format!("{path}: {e}"))?;
            let input = BufReader::with_capacity(1 << 20, file);
            (path, run.write_airr(input, &mut out, double_d))
        }
    };
    let tally = written.map_err(|e| match e {
        RunError::Reads(e) => format!("{source}: {e}"),
        e => e.to_string(),
    })?;
    out.flush()?;
    if let Some(sidecar) = sidecar.as_mut() {
        sidecar.flush()?;
    }
    let seconds = started.elapsed().as_secs_f64();
    let screened = match screening {
        true => format!(", {} double-D", tally.double_d),
        false => String::new(),
    };
    eprintln!(
        "swig: {} reads ({} with a V(D)J path{screened}) in {seconds:.1} s ({:.0} reads/s)",
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

/// Parse an optional probability or per-base rate, which must lie in
/// `[0, 1]` (not NaN).
fn probability(value: &Option<String>, option: &str, default: f64) -> Result<f64, String> {
    let p = number(value, option, default)?;
    match (0.0..=1.0).contains(&p) {
        true => Ok(p),
        false => Err(format!(
            "{option} expects a value in [0, 1], got '{}'",
            value.as_deref().unwrap_or_default()
        )),
    }
}

fn run_simulate(args: &Args) -> Result<(), Box<dyn Error>> {
    let out = args.out.as_deref().expect("checked by Args::check");
    let preset = args.preset.as_deref().unwrap_or("mid");
    let mut params = SimParams::preset(preset)
        .ok_or_else(|| format!("unknown preset '{preset}' (easy, mid, train, hard)"))?
        .with_reverse_complement(probability(&args.rc, "--rc", 0.5)?);
    if args.shm.is_some() {
        params = params.with_shm(probability(&args.shm, "--shm", 0.0)?);
    }
    if args.dd.is_some() {
        params.include_dd = probability(&args.dd, "--dd", 0.0)?;
    }
    let count: usize = number(&args.count, "--reads", 10_000)?;
    if count == 0 {
        return Err("--reads expects a positive count, got '0'".into());
    }
    let seed: u64 = number(&args.seed, "--seed", 1)?;
    let database = load_database(&args.files, &args.run)?;
    if params.include_dd > 0.0 && database.d.as_ref().map_or(0, GermlineSet::len) < 2 {
        eprintln!("swig: warning: --dd has no effect without at least two D alleles");
    }
    let sim = Simulator::from_database(&database, params)?;
    let sets = [
        (SegmentKind::V, Some(&database.v)),
        (SegmentKind::J, Some(&database.j)),
        (SegmentKind::C, database.c.as_ref()),
    ];
    for (kind, set) in sets {
        let (all, drawn) = (set.map_or(0, GermlineSet::len), sim.drawn(kind));
        if drawn < all {
            eprintln!(
                "swig: drawing {} alleles from {drawn} of {all}: the others lack a trusted anchor or a complete coding region",
                kind.label()
            );
        }
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
        "swig: wrote {count} reads to {out}.fasta and their truth to {out}.truth.tsv (preset {preset}, seed {seed})"
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
    let truth_path = args.truth.as_deref().expect("checked by Args::check");
    let calls_path = args.input.as_deref().expect("checked by Args::check");
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
                "swig: {missing} simulated reads have no row in {path}; counted as unassigned"
            );
        }
    }
    // Written, not printed: a closed or full stdout is an error, not a panic.
    let mut out = io::stdout().lock();
    writeln!(out, "{evaluation}")?;
    if let Some(other_path) = &args.other {
        writeln!(out, "\nfirst: {calls_path}\nsecond: {other_path}\n{paired}")?;
    }
    out.flush()?;
    Ok(())
}

fn main() -> ExitCode {
    match parse_args() {
        // A closed stdout leaves no one to tell: neither is an error.
        Ok(Command::Help) => {
            let _ = writeln!(io::stdout(), "{USAGE}");
            ExitCode::SUCCESS
        }
        Ok(Command::Version) => {
            let _ = writeln!(io::stdout(), "swig {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        Ok(Command::Run(args)) => {
            let result = match args.mode {
                Mode::Assign => run_assign(&args),
                Mode::Prepare => run_prepare(&args),
                Mode::Simulate => run_simulate(&args),
                Mode::Evaluate => run_evaluate(&args),
            };
            match result {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("swig: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        Err(e) => {
            eprintln!("swig: {e}\n\n{USAGE}");
            ExitCode::from(2)
        }
    }
}
