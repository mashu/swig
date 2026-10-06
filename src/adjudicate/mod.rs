//! Adjudication: choose a topology, alleles and cut points from placements.
//!
//! Input: all admitted V, D and J placements of one oriented query.
//! Output: a [`CallPath`] — an ordered list of clipped placements.
//!
//! Two implementations:
//! * [`PlacementDag`] — the note's target architecture: semi-Markov Viterbi
//!   over placements with a pluggable [`Grammar`] (VJ / VDJ / VDDJ …) and
//!   [`PathModel`] (bonuses, D-presence prior, N/P edge scores). It subsumes
//!   AER-R: "candidate A" and "candidate B" are two particular cuts of an
//!   overlap, and the DAG evaluates *every* cut exactly.
//! * [`ResidualAdjudicator`] — the ordinary commit-V/J-then-D path, kept as a
//!   baseline and for regression comparisons.
//!
//! [`PlacementDag`] also reads the objective as a distribution over paths and
//! cuts, `P ∝ exp(λ·score)`, and returns posterior [`Marginals`] of the
//! placements and of the D count by forward–backward over the same DAG.

mod dag;
mod model;
mod posterior;
mod residual;

pub use dag::{DagScratch, PlacementDag};
pub use model::{EngineeringObjective, Grammar, PathModel, VdjGrammar};
pub use posterior::{Marginals, PosteriorScratch, MAX_SPREAD};
pub use residual::{ResidualAdjudicator, ResidualScratch};

use std::cmp::Ordering;

use crate::clip::ClipSpan;
use crate::placement::Placement;
use crate::scoring::Score;

/// Path objective with its tie rules built in, compared lexicographically:
///
/// 1. `score` — the objective;
/// 2. fewer D segments ("with-D must strictly beat no-D");
/// 3. `support` — the summed *unclipped* scores of the path's placements.
///
/// Rule 3 separates paths that explain the read equally well but differ in
/// how much germline sequence supports them, typically alleles of different
/// 3′ length at a junction where a base can belong to either side. The
/// allele with more support admits more cuts that reach the same score, so
/// it has the larger likelihood summed over cuts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PathScore {
    pub score: Score,
    pub d_count: u8,
    pub support: Score,
}

impl Ord for PathScore {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .cmp(&other.score)
            .then(other.d_count.cmp(&self.d_count))
            .then(self.support.cmp(&other.support))
    }
}

impl PartialOrd for PathScore {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// One placement on the called path, clipped to `[entry, exit)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CalledSegment {
    /// Index into the placement slice given to the adjudicator.
    pub placement: usize,
    /// Query cut bounds the placement was clipped to.
    pub entry: u32,
    pub exit: u32,
    /// The retained sub-alignment (its bounds may be tighter than the cuts).
    pub span: ClipSpan,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CallPath {
    /// Segments in query order.
    pub segments: Vec<CalledSegment>,
    pub score: PathScore,
}

pub trait Adjudicator: Send + Sync {
    type Scratch: Default + Send;

    fn adjudicate(
        &self,
        query_len: usize,
        placements: &[Placement],
        scratch: &mut Self::Scratch,
    ) -> Option<CallPath>;

    /// Posterior marginals of `placements`, for adjudicators with a
    /// probabilistic reading that are configured to compute them.
    fn marginals(
        &self,
        _placements: &[Placement],
        _scratch: &mut Self::Scratch,
    ) -> Option<Marginals> {
        None
    }
}
