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

mod dag;
mod model;
mod residual;

pub use dag::{DagScratch, PlacementDag};
pub use model::{DPresencePenalty, EngineeringObjective, Grammar, PathModel, VdjGrammar};
pub use residual::{ResidualAdjudicator, ResidualScratch};

use std::cmp::Ordering;

use crate::clip::ClipSpan;
use crate::placement::Placement;
use crate::scoring::Score;

/// Path objective with the note's tie rule built in: at equal score the path
/// with fewer D segments wins ("with-D must strictly beat no-D").
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PathScore {
    pub score: Score,
    pub d_count: u8,
}

impl Ord for PathScore {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .cmp(&other.score)
            .then(other.d_count.cmp(&self.d_count))
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
}
