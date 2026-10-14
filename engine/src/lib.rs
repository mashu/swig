//! # swig
//!
//! Seeded, banded, affine V(D)J germline assignment with placement-DAG
//! adjudication: Swig's aligner. It began as a Rust implementation of the
//! SwiftIG method note, organised so that every heuristic is a replaceable
//! component.
//!
//! ```text
//! alphabet / fasta / germline      germline data
//! locus / metadata                 loci of a database; SWIGMETA header values
//! translate / annotate / prepare   codon table, anchors and region starts,
//!                                  deduplicated and annotated databases
//! reads                            streaming FASTA / FASTQ / AIRR-table reads
//! orient                           strand from germline k-mer support
//! scoring                          substitution profiles and affine gaps
//! seed → index → retrieve          candidate generation (k-mers as integers)
//! band → align                     banded affine Smith–Waterman (Gotoh):
//!                                  scalar kernel + SIMD inter-sequence batches
//!                                  (AVX-512, AVX2, NEON, WebAssembly simd128)
//! segment / family_bound / flank   retrieve → score-only DP → keep → traceback
//! clip / evidence / placement      O(1) clip-and-rescore, admission evidence
//! window / admission               where D is searched, which D hits count
//! adjudicate                       placement DAG (VJ / VDJ)
//! constant                         C anchored at the J exit, after adjudication
//! pipeline / config / profiles     stage-wise, chunk-batched composition and the
//!                                  note's calling profiles
//! report / read_regions / stats    AIRR output, FR/CDR projection, E-values
//! double_d                         Swig's double-D screen and its sidecar
//! batch / run                      parallel batches; a whole run from named
//!                                  options to AIRR TSV (CLI and WebAssembly)
//! sim                              IgSim simulator and accuracy evaluation
//! ```
//!
//! Quick start: [`prepare`](prepare::prepare) deduplicates the germline sets
//! and annotates their anchors and regions, which the assigner reads.
//!
//! ```ignore
//! use swig::prelude::*;
//! let database = prepare(v, Some(d), j, AnnotateOptions::default())?;
//! let config = CallingProfile::ROptimized.config();
//! let assigner = DefaultAssigner::from_config(&config, database.v, database.d, database.j)?;
//! let mut ws = assigner.workspace();
//! let call = assigner.assign(&encode(b"CAGGTGCAG..."), &mut ws);
//! ```
//!
//! A whole run from FASTA text to AIRR TSV, as the command line performs
//! it, is [`run::Run`].

pub mod adjudicate;
pub mod admission;
pub mod align;
pub mod alphabet;
pub mod annotate;
pub mod band;
pub mod batch;
pub mod clip;
pub mod config;
pub mod constant;
pub mod double_d;
pub mod evidence;
pub mod family_bound;
pub mod fasta;
pub mod flank;
pub mod germline;
pub mod index;
pub mod locus;
pub mod locus_infer;
pub mod metadata;
pub mod orient;
pub mod pipeline;
pub mod placement;
pub mod prepare;
pub mod profiles;
pub mod read_regions;
pub mod reads;
pub mod report;
pub mod retrieve;
pub mod run;
pub mod scoring;
pub mod seed;
pub mod segment;
pub mod sim;
pub mod stats;
#[cfg(test)]
mod testing;
pub mod translate;
pub mod window;

/// The names most programs need.
pub mod prelude {
    pub use crate::adjudicate::{
        Adjudicator, CallModel, CallPath, CalledSegment, DAlleleDecoding, EngineeringObjective,
        Grammar, PathModel, PathScore, PlacementDag, VdjGrammar,
    };
    pub use crate::admission::{
        Admission, AdmitAll, AllOf, AnyOf, DistributedEvidence, ExactRunFloor,
    };
    pub use crate::align::{AlignJob, Aligner, BandedGotoh, InterSequenceGotoh};
    pub use crate::alphabet::{decode, encode};
    pub use crate::annotate::{AnnotateOptions, ChainType};
    pub use crate::band::Band;
    pub use crate::batch::{assign_all, assign_map};
    pub use crate::config::{
        BuildError, CSegmentConfig, DPrior, DSearchSpace, DSegmentConfig, DefaultAssigner,
        EngineConfig, JSegmentConfig, ObjectiveConfig, VSegmentConfig,
    };
    pub use crate::germline::{GermlineSet, SegmentKind};
    pub use crate::orient::{Orientation, OrientationPolicy, StrandIndex};
    pub use crate::pipeline::{AlleleCatalog, Assigner, AssignerOptions, Call, VdjAssigner};
    pub use crate::prepare::{GermlineTexts, PreparedDatabase, prepare};
    pub use crate::profiles::CallingProfile;
    pub use crate::reads::{ReadFormat, ReadReader, ReadRecord};
    pub use crate::report::{AirrRecord, SupportModels};
    pub use crate::retrieve::{Candidate, ExhaustiveRetriever, Retriever};
    pub use crate::run::{Rows, Run, RunError, RunOptions, Strand, Tally};
    pub use crate::scoring::{AffineScoring, HotspotScoring};
    pub use crate::stats::KarlinAltschul;
    pub use crate::window::{
        DWindowPolicy, DownstreamWindow, ExpandedResidual, ResidualWindow, WithAnchorWindow,
    };
}
