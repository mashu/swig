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
//! cuts, `P ∝ exp(nats)`, and returns posterior [`Marginals`] of the
//! placements and of the D count by forward–backward over the same DAG, with
//! alignment evidence rescored at the read's own divergence. Both passes read
//! one [`ReadModel`]: under [`CallModel::PerRead`] the call is a mode of the
//! reported posterior.

mod dag;
mod decode;
mod lattice;
mod model;
mod posterior;
mod residual;
mod units;

pub use dag::{DagScratch, PlacementDag};
pub use decode::{DAlleleDecoding, Topology, Within, max_marginal_d};
pub use model::{EngineeringObjective, Grammar, PathModel, VdjGrammar};
pub use posterior::{MAX_SPREAD, Marginals, OVERLAP_SPREAD, PosteriorScratch};
pub use residual::{ResidualAdjudicator, ResidualScratch};
pub use units::{
    CallModel, EvidenceModel, EvidenceScale, PathNats, ReadModel, UNITS_PER_NAT, Units, divergence,
};

use std::cmp::Ordering;

use crate::clip::ClipSpan;
use crate::placement::Placement;
use crate::scoring::Score;

/// Path objective with its tie rules built in, compared lexicographically:
///
/// 1. `score` — the objective, in fixed scores or in a [`ReadModel`]'s units;
/// 2. fewer D segments ("with-D must strictly beat no-D");
/// 3. `support` — the summed *unclipped* scores of the path's placements.
///
/// Rule 3 separates paths that explain the read equally well but differ in
/// how much germline sequence supports them, typically alleles of different
/// 3′ length at a junction where a base can belong to either side. The
/// allele with more support admits more cuts that reach the same score, so
/// it has the larger likelihood summed over cuts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PathScore<S = Score> {
    pub score: S,
    pub d_count: u8,
    pub support: Score,
}

impl<S: Ord> Ord for PathScore<S> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .cmp(&other.score)
            .then(other.d_count.cmp(&self.d_count))
            .then(self.support.cmp(&other.support))
    }
}

impl<S: Ord> PartialOrd for PathScore<S> {
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
    /// Chosen by its allele's posterior marginal
    /// ([`DAlleleDecoding::MaxMarginal`]), so that allele is reported first;
    /// otherwise the call's tie rules order co-optimal names.
    pub by_marginal: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CallPath {
    /// Segments in query order.
    pub segments: Vec<CalledSegment>,
    /// The path's fixed-scale objective, whichever model chose it.
    pub score: PathScore,
    /// The read's divergence `p̂` ([`divergence`]): the posterior reads
    /// alignment evidence at it, and so does a [`CallModel::PerRead`] call.
    pub divergence: f64,
    /// The model the call maximised.
    pub model: CallModel,
}

pub trait Adjudicator: Send + Sync {
    type Scratch: Default + Send;

    fn adjudicate(
        &self,
        query_len: usize,
        placements: &[Placement],
        scratch: &mut Self::Scratch,
    ) -> Option<CallPath>;

    /// Posterior marginals of `placements` given the adjudicated `path`, for
    /// adjudicators with a probabilistic reading that are configured to
    /// compute them.
    fn marginals(
        &self,
        _placements: &[Placement],
        _path: &CallPath,
        _scratch: &mut Self::Scratch,
    ) -> Option<Marginals> {
        None
    }

    /// Another path to report, when the adjudicator decodes a segment from
    /// the `marginals` instead of keeping the call's (see
    /// [`DAlleleDecoding`]).
    fn decode(
        &self,
        _placements: &[Placement],
        _path: &CallPath,
        _marginals: &Marginals,
        _scratch: &mut Self::Scratch,
    ) -> Option<CallPath> {
        None
    }
}
