//! # vdj_assign
//!
//! Seeded, banded, affine V(D)J germline assignment with placement-DAG
//! adjudication — a Rust implementation of the SwiftIG method note, organised
//! so that every heuristic is a replaceable component.
//!
//! ```text
//! alphabet / fasta / germline      germline data
//! reads                            streaming FASTA / FASTQ reads
//! seed → index → retrieve          candidate generation (k-mers as integers)
//! band → align                     banded affine Smith–Waterman (Gotoh):
//!                                  scalar kernel + AVX2 inter-sequence batches
//! clip / evidence / placement      O(1) clip-and-rescore, admission evidence
//! segment                          retrieve → score-only DP → keep → traceback
//! window / admission               where D is searched, which D hits count
//! adjudicate                       placement DAG (VJ / VDJ / VDDJ) or residual
//! refine                           winning-path polish
//! pipeline / config / profiles     stage-wise, chunk-batched composition and the
//!                                  note's calling profiles
//! report / stats / batch           AIRR output, E-values, parallel batches
//! ```
//!
//! Quick start:
//!
//! ```ignore
//! use vdj_assign::prelude::*;
//! let assigner = DefaultAssigner::from_config(&CallingProfile::ROptimized.config(), v, Some(d), j)?;
//! let mut ws = assigner.workspace();
//! let call = assigner.assign(&encode(b"CAGGTGCAG..."), &mut ws);
//! ```

pub mod adjudicate;
pub mod admission;
pub mod align;
pub mod alphabet;
pub mod annotate;
pub mod band;
pub mod batch;
pub mod clip;
pub mod config;
pub mod evidence;
pub mod family_bound;
pub mod fasta;
pub mod flank;
pub mod germline;
pub mod index;
pub mod orient;
pub mod pipeline;
pub mod placement;
pub mod prepare;
pub mod profiles;
pub mod read_regions;
pub mod reads;
pub mod report;
pub mod retrieve;
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
        Grammar, PathModel, PathScore, PlacementDag, ResidualAdjudicator, VdjGrammar,
    };
    pub use crate::admission::{
        Admission, AdmitAll, AllOf, AnyOf, DistributedEvidence, ExactRunFloor,
    };
    pub use crate::align::{AlignJob, Aligner, BandedGotoh, InterSequenceGotoh};
    pub use crate::alphabet::{decode, encode};
    pub use crate::band::Band;
    pub use crate::batch::{assign_all, assign_map};
    pub use crate::config::{
        BuildError, DPrior, DSegmentConfig, DefaultAssigner, EngineConfig, JSegmentConfig,
        ObjectiveConfig, VSegmentConfig,
    };
    pub use crate::germline::{GermlineSet, SegmentKind};
    pub use crate::orient::{Orientation, OrientationPolicy, StrandIndex};
    pub use crate::pipeline::{AlleleCatalog, Assigner, AssignerOptions, Call, VdjAssigner};
    pub use crate::profiles::CallingProfile;
    pub use crate::reads::{ReadFormat, ReadReader, ReadRecord};
    pub use crate::report::{AirrRecord, SupportModels};
    pub use crate::retrieve::{Candidate, ExhaustiveRetriever, Retriever};
    pub use crate::scoring::{AffineScoring, HotspotScoring};
    pub use crate::stats::KarlinAltschul;
    pub use crate::window::{
        DWindowPolicy, DownstreamWindow, ExpandedResidual, ResidualWindow, WithAnchorWindow,
    };
}
