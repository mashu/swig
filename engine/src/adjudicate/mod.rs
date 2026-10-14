//! Adjudication: choose a topology, alleles and cut points from placements.
//!
//! Input: all admitted V, D and J placements of one oriented query.
//! Output: a [`CallPath`] — an ordered list of clipped placements.
//!
//! [`PlacementDag`] is the note's target architecture: semi-Markov Viterbi
//! over placements with a pluggable [`Grammar`] and [`PathModel`]
//! (presence priors, N/P edge scores). It evaluates every cut of an
//! overlap exactly, including the two extreme cuts (the D takes the
//! overlap, or the V or J does). The note's ordinary path, which commits V
//! and J first and searches a D in the residual between them, was a
//! baseline exhibiting "V/J eat D"; the DAG replaced it.
//!
//! [`PlacementDag`] also reads the objective as a distribution over paths and
//! cuts, `P ∝ exp(nats)`, and returns posterior [`Marginals`] of the
//! placements and of the D count by forward–backward over the same DAG, with
//! alignment evidence rescored at the read's own divergence. Both passes read
//! one per-read model, `ReadModel`: under [`CallModel::PerRead`] the call is a mode of the
//! reported posterior.

mod dag;
mod decode;
mod lattice;
mod model;
#[cfg(test)]
mod oracle;
mod posterior;
mod ties;
mod units;

pub use dag::{DagScratch, Mode, PlacementDag};
pub use decode::DAlleleDecoding;
pub use model::{EngineeringObjective, Grammar, PathModel, VdjGrammar};
pub use posterior::{MAX_SPREAD, Marginals, OVERLAP_SPREAD, PosteriorScratch};
pub(crate) use ties::{CoOptimal, co_optimal};
pub use units::{
    CallModel, Columns, DivergenceCounts, EvidenceModel, EvidenceScale, Flank, PathNats,
    ReadContext, ResolvedOpenFlank, UNITS_PER_NAT, Units, beta_binomial_tail, divergence,
    divergence_counts,
};

use std::cmp::Ordering;

use crate::clip::ClipSpan;
use crate::germline::AlleleId;
use crate::placement::Placement;
use crate::scoring::Score;

/// Path objective with its tie rules built in, compared lexicographically:
///
/// 1. `score` — the objective, in fixed scores or in a read model's units;
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
    /// The allele chosen by its posterior marginal
    /// ([`DAlleleDecoding::MaxMarginal`]), one of the placement's, which is
    /// reported first; `None`: the call's tie rules order co-optimal names.
    pub allele: Option<AlleleId>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CallPath {
    /// Segments in query order.
    pub segments: Vec<CalledSegment>,
    /// The path's fixed-scale objective, whichever model chose it.
    pub score: PathScore,
    /// The key the call maximised: its value in `model`'s units for
    /// `context` (the fixed scores under [`CallModel::FixedScale`]), then
    /// fewer D, then support. A path decoded from the posterior keeps its
    /// call's.
    pub value: PathScore<Units>,
    /// The read the call was made for: its divergence `p̂` ([`divergence`]),
    /// at which the posterior reads alignment evidence, and so does a
    /// [`CallModel::PerRead`] call, and what the open-flank cost reads.
    pub context: ReadContext,
    /// The model the call maximised.
    pub model: CallModel,
}

pub trait Adjudicator: Send + Sync {
    type Scratch: Default + Send;

    /// The call of `placements` in the context of their own
    /// ([`ReadContext::of_placements`]: the read ends where they end). A D
    /// that ends its path is then charged for the bases after it up to the
    /// last placement's end only, not to the read's end: a caller that
    /// knows the read's length passes it with
    /// [`Adjudicator::adjudicate_at`].
    fn adjudicate(&self, placements: &[Placement], scratch: &mut Self::Scratch)
    -> Option<CallPath>;

    /// [`Adjudicator::adjudicate`] for a read with context `context` (of
    /// all its placements, of which `placements` may be a part, and its
    /// length): paths of one read's placement subsets adjudicated in one
    /// context are hypotheses of one model, and compare by
    /// [`CallPath::value`]. Adjudicators without a per-read model ignore it.
    fn adjudicate_at(
        &self,
        placements: &[Placement],
        context: &ReadContext,
        scratch: &mut Self::Scratch,
    ) -> Option<CallPath> {
        let _ = context;
        self.adjudicate(placements, scratch)
    }

    /// The smallest local score a D placement needs to matter: one that
    /// scores less can be neither called nor given probability, so the D
    /// stage need not keep it. `Score::MIN` keeps every D.
    fn d_score_floor(&self) -> Score {
        Score::MIN
    }

    /// Whether a path that ends on a D (no J after it) or starts on one (no
    /// V before it) pays the open-flank cost for the J or V it lacks
    /// ([`ResolvedOpenFlank`]). The pair guards then let such a path be
    /// called ([`crate::pipeline::PairGuards`]); without the cost a D is
    /// called only between a V and a J.
    fn charges_open_flank(&self) -> bool {
        false
    }

    /// Posterior marginals of `placements` given the adjudicated `path`, for
    /// adjudicators with a probabilistic reading that are configured to
    /// compute them. With the scratch of the adjudication that returned
    /// `path` (and of no other read's since), what that adjudication
    /// computed for the read is reused.
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
