//! Pairwise local alignment: the `Aligner` abstraction and its data types.
//!
//! The pipeline talks to aligners only through [`Aligner`], split in two calls:
//!
//! * [`Aligner::score`] — score-only pass, `O(m)` memory, no traceback.
//!   Used to rank every retrieved candidate.
//! * [`Aligner::align`] — full pass with traceback, used only for survivors.
//!
//! **Trick — score first, trace back later.** Traceback storage and the walk
//! back dominate memory traffic; most candidates are discarded after scoring,
//! so they never pay for it.
//!
//! Both calls also exist in batch form ([`Aligner::score_batch`],
//! [`Aligner::align_batch`]) over independent [`AlignJob`]s. The default
//! implementations loop over the single-job calls; [`InterSequenceGotoh`]
//! overrides them with an inter-sequence SIMD kernel that aligns 16 or 32
//! jobs at once on AVX-512, AVX2, NEON, WebAssembly `simd128` or portable
//! lanes ([`Backend`]).

mod gotoh;
mod inter_seq;
mod lanes;

pub use gotoh::{BandedGotoh, GotohWorkspace};
pub use inter_seq::{Backend, InterSequenceGotoh, InterSequenceWorkspace};

use crate::band::Band;
use crate::scoring::{GapCosts, Score, SubstitutionProfile};

/// One alignment column.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Op {
    /// Query and reference both consumed (match or substitution).
    Match,
    /// Query base with no reference base (insertion in the query).
    Ins,
    /// Reference base with no query base (deletion from the query).
    Del,
}

/// A local alignment island. Coordinates are 0-based, half-open.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub q_start: u32,
    pub q_end: u32,
    pub r_start: u32,
    pub r_end: u32,
    pub score: Score,
    pub ops: Vec<Op>,
}

impl Alignment {
    /// Move query coordinates by `offset` (alignments computed on a window slice).
    pub fn shift_query(&mut self, offset: u32) {
        self.q_start += offset;
        self.q_end += offset;
    }

    /// `(op, query_pos, ref_pos)` per column; positions are those *before*
    /// the column consumes its bases.
    pub fn columns(&self) -> impl Iterator<Item = (Op, u32, u32)> + '_ {
        let (mut q, mut r) = (self.q_start, self.r_start);
        self.ops.iter().map(move |&op| {
            let column = (op, q, r);
            match op {
                Op::Match => {
                    q += 1;
                    r += 1;
                }
                Op::Ins => q += 1,
                Op::Del => r += 1,
            }
            column
        })
    }
}

/// Result of a score-only pass. `q_end`/`r_end` are exclusive end coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ScoreHit {
    pub score: Score,
    pub q_end: u32,
    pub r_end: u32,
}

/// One independent alignment problem of a batch.
#[derive(Clone, Copy)]
pub struct AlignJob<'a> {
    pub query: &'a [u8],
    pub profile: &'a SubstitutionProfile,
    pub gaps: GapCosts,
    pub band: Band,
}

/// Local affine alignment of a query slice against one allele profile.
///
/// Implementations must return `None` when the best local score is not positive,
/// and alignments that begin and end with [`Op::Match`] columns (true for any
/// optimal local alignment; [`crate::clip`] relies on it). An alignment's
/// score is the score-only optimum of the same problem,
/// `align(job).score == score(job).score` (both `None` together):
/// [`crate::segment`] ranks a candidate it aligned instead of scoring by that
/// score.
pub trait Aligner {
    /// Per-thread scratch memory, reused across calls.
    type Workspace: Default + Send;

    fn score(
        &self,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
        band: Band,
        workspace: &mut Self::Workspace,
    ) -> Option<ScoreHit>;

    fn align(
        &self,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
        band: Band,
        workspace: &mut Self::Workspace,
    ) -> Option<Alignment>;

    /// Score every job; `out[k]` answers `jobs[k]`. Results must equal
    /// [`Aligner::score`] job by job.
    fn score_batch(
        &self,
        jobs: &[AlignJob<'_>],
        workspace: &mut Self::Workspace,
        out: &mut Vec<Option<ScoreHit>>,
    ) {
        out.clear();
        out.extend(
            jobs.iter()
                .map(|job| self.score(job.query, job.profile, job.gaps, job.band, workspace)),
        );
    }

    /// Align every job; `out[k]` answers `jobs[k]`. Results must equal
    /// [`Aligner::align`] job by job.
    fn align_batch(
        &self,
        jobs: &[AlignJob<'_>],
        workspace: &mut Self::Workspace,
        out: &mut Vec<Option<Alignment>>,
    ) {
        out.clear();
        out.extend(
            jobs.iter()
                .map(|job| self.align(job.query, job.profile, job.gaps, job.band, workspace)),
        );
    }
}
