//! An assignment run as the command line and the WebAssembly module perform
//! it: options by name, the assigner they build, and reads to AIRR rows.
//!
//! [`RunOptions`] holds the user-facing switches. [`RunOptions::set`] reads
//! them by the names the command line spells (`profile`, `strand`, …), so
//! every front end shares one vocabulary, one validation and one mapping to
//! [`EngineConfig`]. [`Run`] builds the assigner from the options and a
//! [`PreparedDatabase`]; [`Run::assign`] turns a batch of reads into AIRR
//! rows (and double-D sidecar rows, [`crate::double_d`]), and
//! [`Run::write_airr`] streams any [`BufRead`] to any [`Write`].
//!
//! Identical sequences are assigned once: within a batch, and across
//! batches through a bounded cache of rows ([`CACHE`]). A read's row
//! depends on its sequence alone, besides its identifier and qualities, so
//! this changes no output.
//!
//! ```ignore
//! use swig::prelude::*;
//! let texts = GermlineTexts { v: &v_fasta, d: Some(&d_fasta), j: &j_fasta, c: None };
//! let mut options = RunOptions::default();
//! options.set("profile", Some("r"))?;
//! let database = PreparedDatabase::from_fasta(texts, options.annotate())?;
//! let mut run = Run::new(&options, database)?;
//! let tally = run.write_airr(reads_fastq.as_bytes(), &mut std::io::stdout().lock(), None)?;
//! ```

use std::collections::{HashMap, VecDeque};
use std::fmt;
use std::io::{self, BufRead, Write};

use crate::adjudicate::{CallModel, DAlleleDecoding, ResidualAdjudicator};
use crate::admission::DefaultAdmission;
use crate::align::{Backend, InterSequenceGotoh};
use crate::alphabet::encode;
use crate::annotate::{AnnotateOptions, ChainType};
use crate::batch::assign_map;
use crate::config::{BuildError, DPrior, DefaultAssigner, DefaultVRetriever, EngineConfig};
use crate::double_d::{self, DoubleDMode, DoubleDOptions, DoubleDScreen};
use crate::germline::SegmentKind;
use crate::orient::OrientationPolicy;
use crate::pipeline::{AlleleCatalog, Assigner, VdjAssigner};
use crate::prepare::PreparedDatabase;
use crate::profiles::CallingProfile;
use crate::reads::{ReadError, ReadReader, ReadRecord};
use crate::report::{AirrRecord, SupportModels};
use crate::retrieve::{ExhaustiveRetriever, LmerMaskRetriever};
use crate::scoring::{AffineScoring, Score};
use crate::segment::HitAcceptance;
use crate::stats::KarlinAltschul;
use crate::window::{ExpandedResidual, WithAnchorWindow};

/// Reads assigned per batch by [`Run::write_airr`] (bounds memory on large
/// inputs). swig-cli keeps more records than this in flight when it streams
/// through the binary (`ENGINE_BATCH` in cli-src/native-aligner.mjs).
pub const BATCH: usize = 16_384;

/// Rows [`Run`] keeps for sequences it may see again.
pub const CACHE: usize = 4096;

/// An option that is unknown, lacks its value, or cannot take the one given.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptionError(pub String);

impl fmt::Display for OptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for OptionError {}

/// The strands reads are assigned on.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Strand {
    /// The read or its reverse complement, whichever germline k-mers support.
    #[default]
    Both,
    /// The read as given.
    Forward,
    /// The reverse complement only.
    Reverse,
}

impl Strand {
    fn policy(self) -> OrientationPolicy {
        match self {
            Self::Both => OrientationPolicy::default(),
            Self::Forward => OrientationPolicy::FORWARD_ONLY,
            Self::Reverse => OrientationPolicy::REVERSE_ONLY,
        }
    }
}

/// The switches of an assignment run.
#[derive(Clone, Debug, PartialEq)]
pub struct RunOptions {
    pub profile: CallingProfile,
    /// Chain type of the V database; `None` votes it per locus.
    pub chain: Option<ChainType>,
    pub strand: Strand,
    /// Add the conserved-anchor D window.
    pub anchor_window: bool,
    /// Report posterior probabilities (the residual adjudicator has none).
    pub posterior: bool,
    /// Commit V and J, then search D in the residual, instead of the DAG.
    pub residual: bool,
    /// The model the call maximises; `None` keeps the profile's.
    pub calls: Option<CallModel>,
    /// How the reported D allele is chosen; `None` keeps the profile's.
    pub d_allele: Option<DAlleleDecoding>,
    /// D prior from the D search at this D-presence rate ([`DPrior::Search`]).
    pub d_presence_rate: Option<f64>,
    /// Call C at an anchored E-value at most this (needs a C database).
    pub c_evalue: Option<f64>,
    /// Report C only where the germline from its first base reaches this
    /// identity (`c_prefix_identity`).
    pub c_prefix_identity: Option<f64>,
    /// Drop V, D and J hits below this identity.
    pub min_identity: Option<f64>,
    /// Shortest exact run of an admissible D; `None` keeps the profile's.
    pub min_d_match: Option<u32>,
    /// D and J mismatch scores; `None` keeps the profile's.
    pub d_mismatch: Option<Score>,
    pub j_mismatch: Option<Score>,
    /// Fewest aligned J columns (matches and mismatches) of a J hit.
    pub min_j_length: Option<u32>,
    /// Swig's double-D screen ([`crate::double_d`]).
    pub double_d: DoubleDOptions,
    /// Alignment kernel; `None` picks the fastest this CPU runs.
    pub simd: Option<Backend>,
    /// Worker threads; 0 uses every core.
    pub threads: usize,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self {
            profile: CallingProfile::ROptimized,
            chain: None,
            strand: Strand::Both,
            anchor_window: true,
            posterior: true,
            residual: false,
            calls: None,
            d_allele: None,
            d_presence_rate: None,
            c_evalue: None,
            c_prefix_identity: None,
            min_identity: None,
            min_d_match: None,
            d_mismatch: None,
            j_mismatch: None,
            min_j_length: None,
            double_d: DoubleDOptions::default(),
            simd: None,
            threads: 0,
        }
    }
}

/// Parse `value` of option `name` with `parse`, or say what it expects.
fn parsed<T>(
    name: &str,
    value: &str,
    expects: &str,
    parse: impl Fn(&str) -> Option<T>,
) -> Result<T, OptionError> {
    parse(value).ok_or_else(|| OptionError(format!("--{name} expects {expects}, got '{value}'")))
}

/// A fraction in [0, 1].
fn fraction(value: &str) -> Option<f64> {
    value.parse().ok().filter(|f: &f64| (0.0..=1.0).contains(f))
}

/// `scoring` with mismatch score `mismatch`.
fn with_mismatch(scoring: AffineScoring, mismatch: Score) -> AffineScoring {
    AffineScoring::new(
        scoring.matched,
        mismatch,
        scoring.gaps.open,
        scoring.gaps.extend,
    )
}

impl RunOptions {
    /// Every option [`Self::set`] reads, and whether it takes a value.
    pub const NAMES: [(&'static str, bool); 24] = [
        ("profile", true),
        ("chain", true),
        ("strand", true),
        ("threads", true),
        ("simd", true),
        ("calls", true),
        ("d-allele", true),
        ("d-presence-rate", true),
        ("c-evalue", true),
        ("c-prefix-identity", true),
        ("min-identity", true),
        ("min-d-match", true),
        ("d-mismatch", true),
        ("j-mismatch", true),
        ("min-j-length", true),
        ("double-d", true),
        ("double-d-min-span", true),
        ("double-d-seed-length", true),
        ("double-d-pseudo-trim", true),
        ("double-d-pseudo-mismatches", true),
        ("double-d-min-gain", true),
        ("no-anchor-window", false),
        ("no-posterior", false),
        ("residual", false),
    ];

    /// Whether option `name` takes a value; `None` for an unknown option.
    pub fn takes_value(name: &str) -> Option<bool> {
        Self::NAMES
            .iter()
            .find(|(known, _)| *known == name)
            .map(|&(_, value)| value)
    }

    /// Set option `name` (spelled as on the command line, without `--`)
    /// from `value`, which a flag must not have and any other option needs.
    pub fn set(&mut self, name: &str, value: Option<&str>) -> Result<(), OptionError> {
        let Some(takes_value) = Self::takes_value(name) else {
            return Err(OptionError(format!("unknown option '--{name}'")));
        };
        let value = match (takes_value, value) {
            (true, Some(value)) => value,
            (true, None) => return Err(OptionError(format!("--{name} needs a value"))),
            (false, Some(value)) => {
                return Err(OptionError(format!(
                    "--{name} takes no value, got '{value}'"
                )));
            }
            (false, None) => "",
        };
        let count = |v: &str| v.parse::<u32>().ok();
        let positive = |v: &str| v.parse::<u32>().ok().filter(|&n| n > 0);
        let penalty = |v: &str| v.parse::<Score>().ok().filter(|s| *s < 0);
        match name {
            "profile" => {
                self.profile = parsed(
                    name,
                    value,
                    "r, truth, igblast or sensitive-d",
                    CallingProfile::from_name,
                )?;
            }
            "chain" => {
                self.chain = Some(parsed(name, value, "heavy or light", ChainType::from_name)?)
            }
            "strand" => {
                self.strand = parsed(name, value, "both, forward or reverse", |v| match v {
                    "both" => Some(Strand::Both),
                    "forward" => Some(Strand::Forward),
                    "reverse" => Some(Strand::Reverse),
                    _ => None,
                })?;
            }
            "threads" => self.threads = parsed(name, value, "a number", count)? as usize,
            "simd" => {
                let names: Vec<&str> = Backend::ALL.iter().map(|b| b.name()).collect();
                self.simd = Some(parsed(name, value, &names.join(", "), Backend::from_name)?);
            }
            "calls" => {
                self.calls = Some(parsed(name, value, "per-read or fixed", |v| match v {
                    "per-read" => Some(CallModel::PerRead),
                    "fixed" => Some(CallModel::FixedScale),
                    _ => None,
                })?);
            }
            "d-allele" => {
                self.d_allele = Some(parsed(name, value, "map or max-marginal", |v| match v {
                    "map" => Some(DAlleleDecoding::Map),
                    "max-marginal" => Some(DAlleleDecoding::MaxMarginal),
                    _ => None,
                })?);
            }
            "d-presence-rate" => {
                let rate = |v: &str| v.parse().ok().filter(|r: &f64| *r > 0.0 && *r < 1.0);
                self.d_presence_rate = Some(parsed(name, value, "a rate in (0, 1)", rate)?);
            }
            "c-evalue" => {
                let evalue = |v: &str| v.parse().ok().filter(|e: &f64| *e > 0.0);
                self.c_evalue = Some(parsed(name, value, "a positive number", evalue)?);
            }
            "c-prefix-identity" => {
                // 0 turns the gate off.
                let gate = parsed(name, value, "a fraction in [0, 1]", fraction)?;
                self.c_prefix_identity = (gate > 0.0).then_some(gate);
            }
            "min-identity" => {
                self.min_identity = Some(parsed(name, value, "a fraction in [0, 1]", fraction)?);
            }
            "min-d-match" => {
                self.min_d_match = Some(parsed(name, value, "a positive number", positive)?);
            }
            "d-mismatch" => {
                self.d_mismatch = Some(parsed(name, value, "a negative score", penalty)?);
            }
            "j-mismatch" => {
                self.j_mismatch = Some(parsed(name, value, "a negative score", penalty)?);
            }
            "min-j-length" => self.min_j_length = Some(parsed(name, value, "a number", count)?),
            "double-d" => {
                self.double_d.mode =
                    parsed(name, value, "off, all or long-span", DoubleDMode::from_name)?;
            }
            "double-d-min-span" => {
                self.double_d.min_vj_span = parsed(name, value, "a number", count)?;
            }
            "double-d-seed-length" => {
                let seed = |v: &str| {
                    v.parse()
                        .ok()
                        .filter(|l| double_d::SEED_LENGTHS.contains(l))
                };
                self.double_d.seed_length = parsed(name, value, "a number from 6 to 24", seed)?;
            }
            "double-d-pseudo-trim" => {
                self.double_d.pseudo_trim = parsed(name, value, "a number", count)?;
            }
            "double-d-pseudo-mismatches" => {
                self.double_d.pseudo_mismatches = parsed(name, value, "a number", count)?;
            }
            "double-d-min-gain" => {
                let gain = |v: &str| v.parse::<Score>().ok();
                self.double_d.min_score_gain = parsed(name, value, "a score", gain)?;
            }
            "no-anchor-window" => self.anchor_window = false,
            "no-posterior" => self.posterior = false,
            "residual" => self.residual = true,
            _ => unreachable!("every name in NAMES is handled"),
        }
        Ok(())
    }

    /// How germline databases are annotated.
    pub fn annotate(&self) -> AnnotateOptions {
        AnnotateOptions { chain: self.chain }
    }

    /// The engine configuration: the profile's, with these options applied.
    pub fn config(&self) -> EngineConfig {
        let mut config = self.profile.config();
        config.options.orientation = self.strand.policy();
        config.options.min_identity = self.min_identity.unwrap_or(0.0);
        config.anchor_window = self.anchor_window;
        config.simd = self.simd;
        // The residual adjudicator computes no posteriors.
        config.posterior = self.posterior && !self.residual;
        if let Some(calls) = self.calls {
            config.calls = calls;
        }
        if let Some(d_allele) = self.d_allele {
            config.d_allele = d_allele;
        }
        if let Some(presence) = self.d_presence_rate {
            config.objective.d_prior = DPrior::search(presence);
        }
        if let Some(evalue) = self.c_evalue {
            config.c.evalue = evalue;
        }
        if let Some(run) = self.min_d_match {
            config.d.min_exact_run = run;
            config.min_retained[SegmentKind::D.index()] = run;
        }
        if let Some(mismatch) = self.d_mismatch {
            config.d.scoring = with_mismatch(config.d.scoring, mismatch);
        }
        if let Some(mismatch) = self.j_mismatch {
            config.j.scoring = with_mismatch(config.j.scoring, mismatch);
        }
        if let Some(columns) = self.min_j_length {
            config.j.accept = HitAcceptance {
                min_subst: columns,
                ..config.j.accept
            };
        }
        config
    }
}

/// Why a run could not be built or finished.
#[derive(Debug)]
pub enum RunError {
    Option(OptionError),
    Build(BuildError),
    /// The reads could not be read.
    Reads(ReadError),
    /// The rows could not be written.
    Output(io::Error),
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Option(e) => write!(f, "{e}"),
            Self::Build(e) => write!(f, "{e}"),
            Self::Reads(e) => write!(f, "{e}"),
            Self::Output(e) => write!(f, "writing output: {e}"),
        }
    }
}

impl std::error::Error for RunError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Option(e) => Some(e),
            Self::Build(e) => Some(e),
            Self::Reads(e) => Some(e),
            Self::Output(e) => Some(e),
        }
    }
}

impl From<OptionError> for RunError {
    fn from(e: OptionError) -> Self {
        Self::Option(e)
    }
}

impl From<BuildError> for RunError {
    fn from(e: BuildError) -> Self {
        Self::Build(e)
    }
}

/// Reads assigned, how many of them have a V(D)J path, and how many
/// double-D rows the screen wrote.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Tally {
    pub reads: usize,
    pub called: usize,
    pub double_d: usize,
}

impl std::ops::AddAssign for Tally {
    fn add_assign(&mut self, other: Self) {
        self.reads += other.reads;
        self.called += other.called;
        self.double_d += other.double_d;
    }
}

/// The rows of one batch: AIRR and double-D sidecar lines, each ending in
/// a line break, without headers.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Rows {
    pub airr: String,
    pub double_d: String,
    pub tally: Tally,
}

/// The default assigner with the residual adjudicator swapped in.
type ResidualAssigner = Assigner<
    InterSequenceGotoh,
    DefaultVRetriever,
    ExhaustiveRetriever,
    LmerMaskRetriever,
    WithAnchorWindow<ExpandedResidual>,
    DefaultAdmission,
    ResidualAdjudicator,
>;

/// The assigner of a run; each variant monomorphises the pipeline once.
enum Engine {
    Dag(DefaultAssigner),
    Residual(ResidualAssigner),
}

/// What a sequence's rows are, besides the identifier and qualities.
#[derive(Clone, Debug)]
struct Assigned {
    /// The AIRR row from the tab before `sequence_aa`.
    airr: String,
    /// The sidecar row from the tab before `locus`.
    double_d: Option<String>,
    called: bool,
}

/// Rows of recently assigned sequences, oldest dropped first.
#[derive(Debug, Default)]
struct Cache {
    rows: HashMap<Vec<u8>, Assigned>,
    order: VecDeque<Vec<u8>>,
}

impl Cache {
    fn insert(&mut self, sequence: Vec<u8>, assigned: Assigned) {
        if self.order.len() == CACHE
            && let Some(oldest) = self.order.pop_front()
        {
            self.rows.remove(&oldest);
        }
        self.order.push_back(sequence.clone());
        self.rows.insert(sequence, assigned);
    }
}

/// Free text as one TSV field: tabs and line breaks become spaces.
fn field(text: &str) -> String {
    text.replace(['\t', '\r', '\n'], " ")
}

/// `row` from the tab after its first `fields` fields.
fn tail(row: &str, fields: usize) -> &str {
    row.match_indices('\t')
        .nth(fields - 1)
        .map_or("", |(at, _)| &row[at..])
}

/// The [`Assigned`] rows of each of `sequences`, formatted on the worker
/// threads.
fn assigned<T: VdjAssigner + AlleleCatalog + Sync>(
    assigner: &T,
    support: &SupportModels,
    sequences: &[&[u8]],
    threads: usize,
    c_prefix_identity: Option<f64>,
    screen: Option<&DoubleDScreen>,
) -> Vec<Assigned> {
    let encoded: Vec<Vec<u8>> = sequences.iter().map(|s| encode(s)).collect();
    assign_map(assigner, &encoded, threads, |i, call| {
        let sequence = String::from_utf8_lossy(sequences[i]);
        let mut record = AirrRecord::from_call("", &sequence, &call, assigner, support);
        if let Some(gate) = c_prefix_identity
            && record
                .c
                .as_ref()
                .is_some_and(|c| c.prefix_identity.is_none_or(|p| p + 1e-9 < gate))
        {
            record.c = None;
        }
        let double_d = screen
            .zip(assigner.segment(SegmentKind::D))
            .and_then(|(screen, d)| {
                screen.screen(
                    d.set(),
                    &call.query,
                    record.locus,
                    [record.v.as_ref(), record.d.as_ref(), record.j.as_ref()],
                )
            });
        Assigned {
            // `sequence_id`, `sequence` and `quality` are the read's own.
            airr: tail(&record.to_tsv(), 3).to_string(),
            double_d,
            called: call.path.is_some(),
        }
    })
}

/// An assigner built from [`RunOptions`] and a database, writing AIRR rows.
pub struct Run {
    engine: Engine,
    support: SupportModels,
    threads: usize,
    c_prefix_identity: Option<f64>,
    /// The double-D screen, when it is on and there is a D database.
    screen: Option<DoubleDScreen>,
    cache: Cache,
    /// Reads seen so far: a read without an identifier is `sequence_N`.
    ordinal: usize,
}

impl Run {
    /// Build the assigner `options` describe over `database`.
    pub fn new(options: &RunOptions, database: PreparedDatabase) -> Result<Self, RunError> {
        if options.c_evalue.is_some() && database.c.is_none() {
            return Err(OptionError("--c-evalue needs a C database (--c)".into()).into());
        }
        let config = options.config();
        let support =
            [&config.v.scoring, &config.d.scoring, &config.j.scoring].map(KarlinAltschul::ungapped);
        let mut assigner =
            DefaultAssigner::from_config(&config, database.v, database.d, database.j)?;
        if let Some(c) = database.c {
            assigner = assigner.with_constant(c, &config.c)?;
        }
        let screen = assigner
            .segment(SegmentKind::D)
            .filter(|_| options.double_d.mode != DoubleDMode::Off)
            .map(|d| DoubleDScreen::new(d.set(), options.double_d));
        let engine = if options.residual {
            Engine::Residual(assigner.with_adjudicator(ResidualAdjudicator {
                min_retained: config.min_retained,
                ..ResidualAdjudicator::default()
            }))
        } else {
            Engine::Dag(assigner)
        };
        Ok(Self {
            engine,
            support,
            threads: options.threads,
            c_prefix_identity: options.c_prefix_identity,
            screen,
            cache: Cache::default(),
            ordinal: 0,
        })
    }

    /// The alignment kernel the run uses.
    pub fn backend(&self) -> Backend {
        match &self.engine {
            Engine::Dag(assigner) => assigner.aligner.backend(),
            Engine::Residual(assigner) => assigner.aligner.backend(),
        }
    }

    /// The AIRR TSV header, without a line break.
    pub fn airr_header(&self) -> String {
        AirrRecord::tsv_header()
    }

    /// The double-D sidecar's header, without a line break.
    pub fn double_d_header(&self) -> String {
        double_d::header()
    }

    /// The AIRR and double-D rows of `records`, in order: a sidecar row
    /// starts with its read's index in `records`.
    pub fn assign(&mut self, records: &[ReadRecord]) -> Rows {
        // The sequences to assign: each once, and none the cache holds.
        let mut fresh: Vec<&[u8]> = Vec::new();
        let mut queued: HashMap<&[u8], ()> = HashMap::new();
        for record in records {
            let sequence = record.sequence.as_slice();
            if !self.cache.rows.contains_key(sequence) && queued.insert(sequence, ()).is_none() {
                fresh.push(sequence);
            }
        }
        let found = match &self.engine {
            Engine::Dag(a) => assigned(
                a,
                &self.support,
                &fresh,
                self.threads,
                self.c_prefix_identity,
                self.screen.as_ref(),
            ),
            Engine::Residual(a) => assigned(
                a,
                &self.support,
                &fresh,
                self.threads,
                self.c_prefix_identity,
                self.screen.as_ref(),
            ),
        };
        let batch: HashMap<&[u8], Assigned> = fresh.into_iter().zip(found).collect();
        let mut rows = Rows::default();
        for (index, record) in records.iter().enumerate() {
            self.ordinal += 1;
            let sequence = record.sequence.as_slice();
            let assigned = batch
                .get(sequence)
                .or_else(|| self.cache.rows.get(sequence))
                .expect("every sequence is assigned or cached");
            let id = if record.id.is_empty() {
                format!("sequence_{}", self.ordinal)
            } else {
                field(&record.id)
            };
            rows.airr += &id;
            rows.airr.push('\t');
            rows.airr += &AirrRecord::normalized(&String::from_utf8_lossy(sequence));
            rows.airr.push('\t');
            rows.airr += &field(&record.quality);
            rows.airr += &assigned.airr;
            rows.airr.push('\n');
            if let Some(rest) = &assigned.double_d {
                rows.double_d += &format!("{index}\t{id}\t{rest}\n");
                rows.tally.double_d += 1;
            }
            rows.tally.reads += 1;
            rows.tally.called += usize::from(assigned.called);
        }
        for (sequence, assigned) in batch {
            self.cache.insert(sequence.to_vec(), assigned);
        }
        rows
    }

    /// Write the header line, then a line per read of the FASTA, FASTQ or
    /// AIRR `reads`, [`BATCH`] reads at a time; with `double_d`, the
    /// screen's sidecar there (header first, `swig_batch_record_index` the
    /// read's 0-based index in the whole input).
    pub fn write_airr<R: BufRead, W: Write>(
        &mut self,
        reads: R,
        out: &mut W,
        mut double_d: Option<&mut dyn Write>,
    ) -> Result<Tally, RunError> {
        writeln!(out, "{}", self.airr_header()).map_err(RunError::Output)?;
        if let Some(sidecar) = double_d.as_mut() {
            writeln!(sidecar, "{}", self.double_d_header()).map_err(RunError::Output)?;
        }
        let mut reader = ReadReader::new(reads);
        let mut batch = Vec::with_capacity(BATCH);
        let mut tally = Tally::default();
        loop {
            batch.clear();
            if reader
                .read_batch(BATCH, &mut batch)
                .map_err(RunError::Reads)?
                == 0
            {
                return Ok(tally);
            }
            let rows = self.assign(&batch);
            out.write_all(rows.airr.as_bytes())
                .map_err(RunError::Output)?;
            if let Some(sidecar) = double_d.as_mut() {
                // Batch indices become indices in the whole input.
                for line in rows.double_d.lines() {
                    let (index, rest) = line.split_once('\t').expect("an indexed row");
                    let index: usize = index.parse().expect("a batch index");
                    writeln!(sidecar, "{}\t{rest}", tally.reads + index)
                        .map_err(RunError::Output)?;
                }
            }
            tally += rows.tally;
            // A batch's rows reach a reader as soon as they are assigned:
            // swig-cli streams records through a pipe and waits for them.
            out.flush().map_err(RunError::Output)?;
            if let Some(sidecar) = double_d.as_mut() {
                sidecar.flush().map_err(RunError::Output)?;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prepare::GermlineTexts;
    use crate::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

    #[test]
    fn options_are_read_by_their_command_line_names() {
        let mut options = RunOptions::default();
        for (name, value) in [
            ("profile", Some("igblast")),
            ("chain", Some("light")),
            ("strand", Some("reverse")),
            ("threads", Some("3")),
            ("simd", Some("scalar")),
            ("calls", Some("fixed")),
            ("d-allele", Some("max-marginal")),
            ("d-presence-rate", Some("0.9")),
            ("c-evalue", Some("0.01")),
            ("c-prefix-identity", Some("0.97")),
            ("min-identity", Some("0.7")),
            ("min-d-match", Some("5")),
            ("d-mismatch", Some("-4")),
            ("j-mismatch", Some("-4")),
            ("min-j-length", Some("12")),
            ("double-d", Some("long-span")),
            ("double-d-min-span", Some("30")),
            ("double-d-seed-length", Some("9")),
            ("double-d-pseudo-trim", Some("4")),
            ("double-d-pseudo-mismatches", Some("2")),
            ("double-d-min-gain", Some("6")),
            ("no-anchor-window", None),
            ("no-posterior", None),
            ("residual", None),
        ] {
            options.set(name, value).unwrap();
        }
        let expected = RunOptions {
            profile: CallingProfile::IgBlastAgreement,
            chain: Some(ChainType::Light),
            strand: Strand::Reverse,
            anchor_window: false,
            posterior: false,
            residual: true,
            calls: Some(CallModel::FixedScale),
            d_allele: Some(DAlleleDecoding::MaxMarginal),
            d_presence_rate: Some(0.9),
            c_evalue: Some(0.01),
            c_prefix_identity: Some(0.97),
            min_identity: Some(0.7),
            min_d_match: Some(5),
            d_mismatch: Some(-4),
            j_mismatch: Some(-4),
            min_j_length: Some(12),
            double_d: DoubleDOptions {
                mode: DoubleDMode::LongSpan,
                min_vj_span: 30,
                seed_length: 9,
                pseudo_trim: 4,
                pseudo_mismatches: 2,
                min_score_gain: 6,
            },
            simd: Some(Backend::Scalar),
            threads: 3,
        };
        assert_eq!(options, expected);
        let config = options.config();
        assert_eq!(config.options.orientation, OrientationPolicy::REVERSE_ONLY);
        assert_eq!(config.options.min_identity, 0.7);
        assert_eq!(config.objective.d_prior, DPrior::search(0.9));
        assert!(!config.posterior && !config.anchor_window);
        assert_eq!(config.c.evalue, 0.01);
        assert_eq!(config.d.min_exact_run, 5);
        assert_eq!(config.d.scoring, AffineScoring::new(2, -4, -11, -1));
        assert_eq!(config.j.accept.min_subst, 12);
        assert_eq!(
            RunOptions::default().config(),
            CallingProfile::ROptimized.config()
        );
        for (name, value) in [
            ("profile", Some("best")),
            ("strand", Some("plus")),
            ("threads", Some("many")),
            ("d-presence-rate", Some("1")),
            ("c-evalue", Some("0")),
            ("c-prefix-identity", Some("1.5")),
            ("d-mismatch", Some("3")),
            ("double-d", Some("some")),
            ("simd", Some("mmx")),
            ("profile", None),
            ("residual", Some("yes")),
            ("unknown", None),
        ] {
            assert!(options.set(name, value).is_err(), "--{name} {value:?}");
        }
    }

    fn texts_and_reads(count: usize) -> ([String; 4], String) {
        texts_and_simulated(count, SimParams::mid())
    }

    fn texts_and_simulated(count: usize, params: SimParams) -> ([String; 4], String) {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 5);
        let fasta = |set: &crate::germline::GermlineSet| {
            let mut text = Vec::new();
            set.write_fasta(&mut text, |_| None).unwrap();
            String::from_utf8(text).unwrap()
        };
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params).unwrap();
        let reads: String = sim
            .reads(count, 3)
            .unwrap()
            .iter()
            .enumerate()
            .map(|(i, read)| format!(">r{i}\n{}\n", String::from_utf8_lossy(&read.sequence)))
            .collect();
        (
            [fasta(&db.v), fasta(&db.d), fasta(&db.j), fasta(&db.c)],
            reads,
        )
    }

    /// The streamed TSV is the header and one row per read, the rows equal
    /// to [`AirrRecord::from_call`] on the assigner the options build.
    #[test]
    fn a_run_writes_one_airr_row_per_read() {
        let ([v, d, j, c], reads) = texts_and_reads(40);
        let texts = GermlineTexts {
            v: &v,
            d: Some(&d),
            j: &j,
            c: Some(&c),
        };
        let mut options = RunOptions::default();
        options.set("threads", Some("2")).unwrap();
        let database = PreparedDatabase::from_fasta(texts, options.annotate()).unwrap();
        let mut run = Run::new(&options, database).unwrap();
        let mut out = Vec::new();
        let tally = run.write_airr(reads.as_bytes(), &mut out, None).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(tally.reads, 40);
        assert!(tally.called >= 38, "{tally:?}");
        assert_eq!(lines.len(), 41);
        assert_eq!(lines[0], AirrRecord::tsv_header());
        assert_eq!(lines[0].split('\t').nth(12), Some("j_call"));
        let columns = lines[0].split('\t').count();
        assert!(lines.iter().all(|l| l.split('\t').count() == columns));
        assert!(lines[1].starts_with("r0\t"));

        let database = PreparedDatabase::from_fasta(texts, options.annotate()).unwrap();
        let config = options.config();
        let assigner = DefaultAssigner::from_config(&config, database.v, database.d, database.j)
            .unwrap()
            .with_constant(database.c.unwrap(), &config.c)
            .unwrap();
        let mut ws = assigner.workspace();
        let first = ReadReader::new(reads.as_bytes()).next().unwrap().unwrap();
        let call = assigner.assign(&encode(&first.sequence), &mut ws);
        let sequence = String::from_utf8_lossy(&first.sequence);
        let support =
            [&config.v.scoring, &config.d.scoring, &config.j.scoring].map(KarlinAltschul::ungapped);
        let record = AirrRecord::from_call("r0", &sequence, &call, &assigner, &support);
        assert_eq!(lines[1], record.to_tsv());
    }

    /// Repeated sequences, within a batch and across batches, get the rows
    /// a fresh run gives them, under their own identifiers and qualities.
    #[test]
    fn repeated_sequences_reuse_their_rows() {
        let ([v, d, j, _], reads) = texts_and_reads(12);
        let texts = GermlineTexts {
            v: &v,
            d: Some(&d),
            j: &j,
            c: None,
        };
        let options = RunOptions::default();
        let records: Vec<ReadRecord> = ReadReader::new(reads.as_bytes())
            .map(Result::unwrap)
            .collect();
        let mut repeated = records.clone();
        for (k, record) in records.iter().enumerate().take(6) {
            repeated.push(ReadRecord {
                id: format!("copy{k}"),
                quality: "I".repeat(record.sequence.len()),
                ..record.clone()
            });
        }
        let database = || PreparedDatabase::from_fasta(texts, options.annotate()).unwrap();
        let mut run = Run::new(&options, database()).unwrap();
        let first = run.assign(&repeated[..12]);
        let second = run.assign(&repeated[12..]);
        let mut fresh = Run::new(&options, database()).unwrap();
        let alone = fresh.assign(&repeated);
        assert_eq!(first.airr.clone() + &second.airr, alone.airr);
        let lines: Vec<&str> = alone.airr.lines().collect();
        let fields: Vec<&str> = lines[12].split('\t').collect();
        assert_eq!(fields[0], "copy0");
        assert_eq!(fields[2].len(), records[0].sequence.len());
        let original: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields[3..], original[3..]);
    }

    /// The double-D screen leaves every AIRR row as it is and writes
    /// sidecar rows for tandem-D reads.
    #[test]
    fn the_double_d_screen_adds_a_sidecar_only() {
        let params = SimParams {
            include_dd: 0.5,
            ..SimParams::preset("easy").unwrap()
        };
        let ([v, d, j, _], reads) = texts_and_simulated(120, params);
        let texts = GermlineTexts {
            v: &v,
            d: Some(&d),
            j: &j,
            c: None,
        };
        let records: Vec<ReadRecord> = ReadReader::new(reads.as_bytes())
            .map(Result::unwrap)
            .collect();
        let database = || PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
        let plain = Run::new(&RunOptions::default(), database())
            .unwrap()
            .assign(&records);
        let mut options = RunOptions::default();
        options.set("double-d", Some("all")).unwrap();
        options.set("double-d-seed-length", Some("6")).unwrap();
        let screened = Run::new(&options, database()).unwrap().assign(&records);
        assert_eq!(plain.airr, screened.airr);
        assert!(plain.double_d.is_empty());
        let rows: Vec<&str> = screened.double_d.lines().collect();
        assert_eq!(rows.len(), screened.tally.double_d);
        assert!(rows.len() >= 10, "{} double-D rows", rows.len());
        let columns = double_d::COLUMNS;
        for row in rows {
            let fields: Vec<&str> = row.split('\t').collect();
            assert_eq!(fields.len(), columns.len());
            let index: usize = fields[0].parse().unwrap();
            assert_eq!(fields[1], records[index].id);
            let at = |name: &str| fields[columns.iter().position(|c| *c == name).unwrap()];
            assert!(!at("d_call").is_empty() && !at("d2_call").is_empty());
            let d_end: u32 = at("d_sequence_end").parse().unwrap();
            let d2_start: u32 = at("d2_sequence_start").parse().unwrap();
            assert!(d_end < d2_start);
            assert!(at("swig_double_d_score_gain").parse::<i32>().unwrap() >= 8);
        }
    }

    #[test]
    fn a_c_gate_needs_a_c_database() {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 5);
        let database =
            crate::prepare::prepare(db.v, Some(db.d), db.j, AnnotateOptions::default()).unwrap();
        let mut options = RunOptions::default();
        options.set("c-evalue", Some("0.01")).unwrap();
        assert!(matches!(
            Run::new(&options, database),
            Err(RunError::Option(_))
        ));
    }
}
