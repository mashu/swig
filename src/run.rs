//! An assignment run as the command line and the WebAssembly module perform
//! it: options by name, the assigner they build, and reads streamed to AIRR
//! TSV.
//!
//! [`RunOptions`] holds the user-facing switches. [`RunOptions::set`] reads
//! them by the names the command line spells (`profile`, `tandem-d`, …), so
//! every front end shares one vocabulary, one validation and one mapping to
//! [`EngineConfig`]. [`Run`] builds the assigner from the options and a
//! [`PreparedDatabase`] and writes AIRR rows from any [`BufRead`] to any
//! [`Write`]; [`Run::write_rows`] takes parsed records, for front ends that
//! schedule batches themselves.
//!
//! ```ignore
//! use vdj_assign::prelude::*;
//! let texts = GermlineTexts { v: &v_fasta, d: Some(&d_fasta), j: &j_fasta, c: None };
//! let mut options = RunOptions::default();
//! options.set("profile", Some("r"))?;
//! let database = PreparedDatabase::from_fasta(texts, options.annotate())?;
//! let run = Run::new(&options, database)?;
//! let tally = run.write_airr(reads_fastq.as_bytes(), &mut std::io::stdout().lock())?;
//! ```

use std::fmt;
use std::io::{self, BufRead, Write};

use crate::adjudicate::{CallModel, DAlleleDecoding, ResidualAdjudicator, VdjGrammar};
use crate::admission::DefaultAdmission;
use crate::align::{Backend, InterSequenceGotoh};
use crate::alphabet::encode;
use crate::annotate::{AnnotateOptions, ChainType};
use crate::batch::assign_map;
use crate::config::{BuildError, DPrior, DefaultAssigner, DefaultVRetriever, EngineConfig};
use crate::germline::SegmentKind;
use crate::orient::OrientationPolicy;
use crate::pipeline::{AlleleCatalog, Assigner, VdjAssigner};
use crate::prepare::PreparedDatabase;
use crate::profiles::CallingProfile;
use crate::reads::{ReadError, ReadReader, ReadRecord};
use crate::report::{AirrRecord, SupportModels};
use crate::retrieve::{ExhaustiveRetriever, LmerMaskRetriever};
use crate::stats::KarlinAltschul;
use crate::window::{ExpandedResidual, WithAnchorWindow};

/// Reads assigned per batch by [`Run::write_airr`] (bounds memory on large
/// inputs).
pub const BATCH: usize = 16_384;

/// An option that is unknown, lacks its value, or cannot take the one given.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptionError(pub String);

impl fmt::Display for OptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for OptionError {}

/// The switches of an assignment run.
#[derive(Clone, Debug, PartialEq)]
pub struct RunOptions {
    pub profile: CallingProfile,
    /// Chain type of the V database; `None` votes it from the database.
    pub chain: Option<ChainType>,
    /// Allow VDDJ paths.
    pub tandem_d: bool,
    /// Do not try the reverse complement.
    pub forward_only: bool,
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
            tandem_d: false,
            forward_only: false,
            anchor_window: true,
            posterior: true,
            residual: false,
            calls: None,
            d_allele: None,
            d_presence_rate: None,
            c_evalue: None,
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

impl RunOptions {
    /// Every option [`Self::set`] reads, and whether it takes a value.
    pub const NAMES: [(&'static str, bool); 13] = [
        ("profile", true),
        ("chain", true),
        ("threads", true),
        ("simd", true),
        ("calls", true),
        ("d-allele", true),
        ("d-presence-rate", true),
        ("c-evalue", true),
        ("tandem-d", false),
        ("forward-only", false),
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
            "threads" => self.threads = parsed(name, value, "a number", |v| v.parse().ok())?,
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
            "tandem-d" => self.tandem_d = true,
            "forward-only" => self.forward_only = true,
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
        if self.tandem_d {
            config.grammar = VdjGrammar::TANDEM_D;
        }
        if self.forward_only {
            config.options.orientation = OrientationPolicy::FORWARD_ONLY;
        }
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

/// Reads assigned, and how many of them have a V(D)J path.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Tally {
    pub reads: usize,
    pub called: usize,
}

impl std::ops::AddAssign for Tally {
    fn add_assign(&mut self, other: Self) {
        self.reads += other.reads;
        self.called += other.called;
    }
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

/// An assigner built from [`RunOptions`] and a database, writing AIRR rows.
pub struct Run {
    engine: Engine,
    support: SupportModels,
    threads: usize,
}

/// AIRR rows of `records`, each with whether its read has a path, formatted
/// on the worker threads.
fn rows<T: VdjAssigner + AlleleCatalog + Sync>(
    assigner: &T,
    support: &SupportModels,
    records: &[ReadRecord],
    threads: usize,
) -> Vec<(String, bool)> {
    let encoded: Vec<Vec<u8>> = records.iter().map(|r| encode(&r.sequence)).collect();
    assign_map(assigner, &encoded, threads, |i, call| {
        let record = &records[i];
        let sequence = String::from_utf8_lossy(&record.sequence);
        let row = AirrRecord::from_call(record.id(), &sequence, &call, assigner, support);
        (row.to_tsv(), call.path.is_some())
    })
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
        })
    }

    /// Whether rows carry the `c_` columns (a C database was given).
    pub fn has_constant(&self) -> bool {
        match &self.engine {
            Engine::Dag(assigner) => assigner.segment(SegmentKind::C).is_some(),
            Engine::Residual(assigner) => assigner.segment(SegmentKind::C).is_some(),
        }
    }

    /// The AIRR TSV header of this run's rows, without a line break.
    pub fn airr_header(&self) -> String {
        AirrRecord::tsv_header_with(self.has_constant())
    }

    /// Assign `records` and write one AIRR TSV line per record, in order.
    pub fn write_rows<W: Write>(&self, records: &[ReadRecord], out: &mut W) -> io::Result<Tally> {
        let rows = match &self.engine {
            Engine::Dag(assigner) => rows(assigner, &self.support, records, self.threads),
            Engine::Residual(assigner) => rows(assigner, &self.support, records, self.threads),
        };
        let mut tally = Tally::default();
        for (row, called) in rows {
            writeln!(out, "{row}")?;
            tally.reads += 1;
            tally.called += usize::from(called);
        }
        Ok(tally)
    }

    /// Write the header line, then a line per read of the FASTA or FASTQ
    /// `reads`, [`BATCH`] reads at a time.
    pub fn write_airr<R: BufRead, W: Write>(
        &self,
        reads: R,
        out: &mut W,
    ) -> Result<Tally, RunError> {
        writeln!(out, "{}", self.airr_header()).map_err(RunError::Output)?;
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
            tally += self.write_rows(&batch, out).map_err(RunError::Output)?;
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
            ("threads", Some("3")),
            ("simd", Some("scalar")),
            ("calls", Some("fixed")),
            ("d-allele", Some("max-marginal")),
            ("d-presence-rate", Some("0.9")),
            ("c-evalue", Some("0.01")),
            ("tandem-d", None),
            ("forward-only", None),
            ("no-anchor-window", None),
            ("no-posterior", None),
            ("residual", None),
        ] {
            options.set(name, value).unwrap();
        }
        let expected = RunOptions {
            profile: CallingProfile::IgBlastAgreement,
            chain: Some(ChainType::Light),
            tandem_d: true,
            forward_only: true,
            anchor_window: false,
            posterior: false,
            residual: true,
            calls: Some(CallModel::FixedScale),
            d_allele: Some(DAlleleDecoding::MaxMarginal),
            d_presence_rate: Some(0.9),
            c_evalue: Some(0.01),
            simd: Some(Backend::Scalar),
            threads: 3,
        };
        assert_eq!(options, expected);
        let config = options.config();
        assert_eq!(config.grammar, VdjGrammar::TANDEM_D);
        assert_eq!(config.objective.d_prior, DPrior::search(0.9));
        assert!(!config.posterior && !config.anchor_window);
        assert_eq!(config.c.evalue, 0.01);
        assert_eq!(
            RunOptions::default().config(),
            CallingProfile::ROptimized.config()
        );
        for (name, value) in [
            ("profile", Some("best")),
            ("threads", Some("many")),
            ("d-presence-rate", Some("1")),
            ("c-evalue", Some("0")),
            ("simd", Some("mmx")),
            ("profile", None),
            ("residual", Some("yes")),
            ("unknown", None),
        ] {
            assert!(options.set(name, value).is_err(), "--{name} {value:?}");
        }
    }

    /// The streamed TSV is the header and one row per read, the rows equal
    /// to [`AirrRecord::from_call`] on the assigner the options build.
    #[test]
    fn a_run_writes_one_airr_row_per_read() {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 5);
        let fasta = |set: &crate::germline::GermlineSet| {
            let mut text = Vec::new();
            set.write_fasta(&mut text, |_| None).unwrap();
            String::from_utf8(text).unwrap()
        };
        let (v, d, j, c) = (fasta(&db.v), fasta(&db.d), fasta(&db.j), fasta(&db.c));
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::mid()).unwrap();
        let reads: String = sim
            .reads(40, 3)
            .unwrap()
            .iter()
            .enumerate()
            .map(|(i, read)| format!(">r{i}\n{}\n", String::from_utf8_lossy(&read.sequence)))
            .collect();
        let texts = GermlineTexts {
            v: &v,
            d: Some(&d),
            j: &j,
            c: Some(&c),
        };
        let mut options = RunOptions::default();
        options.set("threads", Some("2")).unwrap();
        let database = PreparedDatabase::from_fasta(texts, options.annotate()).unwrap();
        let run = Run::new(&options, database).unwrap();
        let mut out = Vec::new();
        let tally = run.write_airr(reads.as_bytes(), &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(tally.reads, 40);
        assert!(tally.called >= 38, "{tally:?}");
        assert_eq!(lines.len(), 41);
        assert_eq!(lines[0], AirrRecord::tsv_header_with(true));
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
