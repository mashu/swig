//! The assignment pipeline.
//!
//! ```text
//!  chunk of reads ─► per read: orient (germline k-mer support), V candidates
//!                 ─► V stage over the whole chunk     score batch + traceback waves
//!                 ─► per read: J window right of V; J stage over the whole chunk
//!                 ─► per read: D windows from V/J geometry, D candidates
//!                 ─► D stage over the whole chunk ─► per read: admission gates
//!                 ─► per read: adjudicator (placement DAG) ─► refiner ─► Call
//! ```
//!
//! **Trick — stage-wise batching across reads.** After prefiltering a read has
//! only a few V, J and D candidates, far too few to fill a 16-lane SIMD
//! register. Running each stage over a chunk of reads hands the aligner one
//! batch per stage for the whole chunk. Every per-read decision (keep policy,
//! D windows, admission, adjudication) is still taken per read, and batch
//! aligners must answer each job exactly as a single call would, so a read's
//! result does not depend on its chunk: `assign(q)` is literally a chunk of one.
//!
//! [`Assigner`] is generic over every replaceable component; each one is a
//! trait with a single responsibility, and the whole pipeline monomorphises
//! into straight-line code (no trait objects, no runtime type checks).
//! Swapping a component changes the type:
//!
//! ```ignore
//! let experimental = assigner
//!     .with_adjudicator(ResidualAdjudicator::default())
//!     .with_d_window(ResidualWindow { overlap_margin: 2 });
//! ```

use std::ops::Range;

use crate::adjudicate::{Adjudicator, CallPath, DagScratch, Marginals, PlacementDag};
use crate::admission::{Admission, DefaultAdmission};
use crate::align::{Aligner, InterSequenceGotoh, InterSequenceWorkspace};
use crate::alphabet::reverse_complement_into;
use crate::config::DefaultVRetriever;
use crate::germline::SegmentKind;
use crate::orient::{Orientation, OrientationPolicy, StrandIndex};
use crate::placement::Placement;
use crate::refine::{NoRefine, Refiner};
use crate::retrieve::{
    Candidate, ExhaustiveRetriever, LmerMaskRetriever, RetrievalScratch, Retriever,
};
use crate::segment::{BatchScratch, PreparedSegment, SegmentRequest, SegmentStage};
use crate::window::{
    merge_windows, DWindowPolicy, DownstreamWindow, ExpandedResidual, WithAnchorWindow,
};

/// Result of assigning one query.
#[derive(Clone, Debug)]
pub struct Call {
    pub orientation: Orientation,
    /// Oriented query codes; all placement coordinates refer to these.
    pub query: Vec<u8>,
    /// Every retained placement (V, J, admitted D), including non-winners.
    pub placements: Vec<Placement>,
    /// The adjudicated path, if any segment could be placed.
    pub path: Option<CallPath>,
    /// Posterior marginals of `placements` (indexed alike), when the
    /// adjudicator computes them and the read has a call.
    pub marginals: Option<Marginals>,
}

/// Anything that assigns V(D)J to encoded queries with a reusable workspace.
pub trait VdjAssigner: Send + Sync {
    type Workspace: Send;

    fn workspace(&self) -> Self::Workspace;

    /// `query` holds codes from [`crate::alphabet::encode`].
    fn assign(&self, query: &[u8], workspace: &mut Self::Workspace) -> Call;

    /// Assign a chunk of queries, appending one call per query to `out` in
    /// order. Results must equal [`VdjAssigner::assign`] query by query;
    /// implementations may batch work across the chunk.
    fn assign_many<Q: AsRef<[u8]>>(
        &self,
        queries: &[Q],
        workspace: &mut Self::Workspace,
        out: &mut Vec<Call>,
    ) {
        out.extend(
            queries
                .iter()
                .map(|query| self.assign(query.as_ref(), workspace)),
        );
    }
}

/// Access to the prepared germline databases (reporting needs names and sequences).
pub trait AlleleCatalog {
    fn segment(&self, kind: SegmentKind) -> Option<&PreparedSegment>;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AssignerOptions {
    pub orientation: OrientationPolicy,
}

pub struct Assigner<
    A = InterSequenceGotoh,
    RV = DefaultVRetriever,
    RJ = ExhaustiveRetriever,
    RD = LmerMaskRetriever,
    W = WithAnchorWindow<ExpandedResidual>,
    Adm = DefaultAdmission,
    Adj = PlacementDag,
    Rf = NoRefine,
> {
    pub v: SegmentStage<RV>,
    pub j: SegmentStage<RJ>,
    /// `None` for loci without D (light chains).
    pub d: Option<SegmentStage<RD>>,
    /// Germline k-mer presence for orientation.
    pub strand: StrandIndex,
    /// Where J is searched relative to V.
    pub j_window: DownstreamWindow,
    pub aligner: A,
    pub d_window: W,
    pub d_admission: Adm,
    pub adjudicator: Adj,
    pub refiner: Rf,
    pub options: AssignerOptions,
}

/// Per-read state while a chunk moves through the stages.
#[derive(Debug, Default)]
struct ReadState {
    orientation: Option<Orientation>,
    /// Oriented query codes.
    query: Vec<u8>,
    v: Vec<Candidate>,
    j: Vec<Candidate>,
    j_window: Range<usize>,
}

/// Buffers of the D stage.
#[derive(Debug, Default)]
struct DScratch {
    /// `(read, window)` of every D search of the chunk.
    windows: Vec<(usize, Range<usize>)>,
    read_windows: Vec<Range<usize>>,
    candidates: Vec<Vec<Candidate>>,
    /// Per read, the index of its first D placement.
    start: Vec<usize>,
}

/// Per-thread buffers of an [`Assigner`]; every buffer is reused across chunks.
#[derive(Debug, Default)]
pub struct Workspace<W = InterSequenceWorkspace, S = DagScratch> {
    retrieval: RetrievalScratch,
    dp: W,
    batch: BatchScratch,
    reads: Vec<ReadState>,
    placements: Vec<Vec<Placement>>,
    d: DScratch,
    adjudication: S,
    calls: Vec<Call>,
}

/// Keep only elements of `items[from..]` satisfying `keep`, preserving order.
fn retain_from<T>(items: &mut Vec<T>, from: usize, mut keep: impl FnMut(&T) -> bool) {
    let mut kept = from;
    for i in from..items.len() {
        if keep(&items[i]) {
            items.swap(kept, i);
            kept += 1;
        }
    }
    items.truncate(kept);
}

impl<A, RV, RJ, RD, W, Adm, Adj, Rf> Assigner<A, RV, RJ, RD, W, Adm, Adj, Rf> {
    pub fn with_aligner<A2>(self, aligner: A2) -> Assigner<A2, RV, RJ, RD, W, Adm, Adj, Rf> {
        Assigner {
            v: self.v,
            j: self.j,
            d: self.d,
            strand: self.strand,
            j_window: self.j_window,
            aligner,
            d_window: self.d_window,
            d_admission: self.d_admission,
            adjudicator: self.adjudicator,
            refiner: self.refiner,
            options: self.options,
        }
    }

    pub fn with_d_window<W2>(self, d_window: W2) -> Assigner<A, RV, RJ, RD, W2, Adm, Adj, Rf> {
        Assigner {
            v: self.v,
            j: self.j,
            d: self.d,
            strand: self.strand,
            j_window: self.j_window,
            aligner: self.aligner,
            d_window,
            d_admission: self.d_admission,
            adjudicator: self.adjudicator,
            refiner: self.refiner,
            options: self.options,
        }
    }

    pub fn with_d_admission<Adm2>(
        self,
        d_admission: Adm2,
    ) -> Assigner<A, RV, RJ, RD, W, Adm2, Adj, Rf> {
        Assigner {
            v: self.v,
            j: self.j,
            d: self.d,
            strand: self.strand,
            j_window: self.j_window,
            aligner: self.aligner,
            d_window: self.d_window,
            d_admission,
            adjudicator: self.adjudicator,
            refiner: self.refiner,
            options: self.options,
        }
    }

    pub fn with_adjudicator<Adj2>(
        self,
        adjudicator: Adj2,
    ) -> Assigner<A, RV, RJ, RD, W, Adm, Adj2, Rf> {
        Assigner {
            v: self.v,
            j: self.j,
            d: self.d,
            strand: self.strand,
            j_window: self.j_window,
            aligner: self.aligner,
            d_window: self.d_window,
            d_admission: self.d_admission,
            adjudicator,
            refiner: self.refiner,
            options: self.options,
        }
    }

    pub fn with_refiner<Rf2>(self, refiner: Rf2) -> Assigner<A, RV, RJ, RD, W, Adm, Adj, Rf2> {
        Assigner {
            v: self.v,
            j: self.j,
            d: self.d,
            strand: self.strand,
            j_window: self.j_window,
            aligner: self.aligner,
            d_window: self.d_window,
            d_admission: self.d_admission,
            adjudicator: self.adjudicator,
            refiner,
            options: self.options,
        }
    }

    pub fn with_options(mut self, options: AssignerOptions) -> Self {
        self.options = options;
        self
    }
}

/// Requests of one stage: one per read, over `window` with `candidates`.
fn read_requests<'a>(
    reads: &'a [ReadState],
    window: impl Fn(&ReadState) -> Range<usize>,
    candidates: impl Fn(&'a ReadState) -> &'a [Candidate],
) -> Vec<SegmentRequest<'a>> {
    reads
        .iter()
        .enumerate()
        .map(|(read, state)| SegmentRequest {
            read,
            query: &state.query,
            window: window(state),
            candidates: candidates(state),
        })
        .collect()
}

/// The stages of [`VdjAssigner::assign_many`], in order. Per-read stages
/// loop over the chunk; chunk stages hand all reads' jobs to the aligner at once.
impl<A, RV, RJ, RD, W, Adm, Adj, Rf> Assigner<A, RV, RJ, RD, W, Adm, Adj, Rf>
where
    A: Aligner + Send + Sync,
    RV: Retriever,
    RJ: Retriever,
    RD: Retriever,
    W: DWindowPolicy,
    Adm: Admission,
    Adj: Adjudicator,
    Rf: Refiner,
{
    /// Strand from germline k-mer support; when that is not significant,
    /// from the strongest V seed cluster on each strand (`state.query` and
    /// `state.v` serve as scratch).
    fn orient(
        &self,
        query: &[u8],
        state: &mut ReadState,
        retrieval: &mut RetrievalScratch,
    ) -> Orientation {
        let policy = self.options.orientation;
        if !policy.try_reverse {
            return Orientation::Forward;
        }
        let support = self.strand.support(query);
        if support.decisive() {
            return policy.choose(support.forward, support.reverse);
        }
        let strongest =
            |candidates: &[Candidate]| candidates.iter().map(|c| c.votes).max().unwrap_or(0);
        self.v
            .retrieve(query, 0..query.len(), retrieval, &mut state.v);
        let forward = strongest(&state.v);
        reverse_complement_into(query, &mut state.query);
        self.v
            .retrieve(&state.query, 0..query.len(), retrieval, &mut state.v);
        policy.choose(forward, strongest(&state.v))
    }

    /// Stage 1, per read: orientation, oriented query, V candidates.
    fn stage_orient<Q: AsRef<[u8]>>(
        &self,
        queries: &[Q],
        reads: &mut [ReadState],
        retrieval: &mut RetrievalScratch,
    ) {
        for (state, query) in reads.iter_mut().zip(queries) {
            let query = query.as_ref();
            let orientation = self.orient(query, state, retrieval);
            state.orientation = Some(orientation);
            if orientation.is_reverse() {
                reverse_complement_into(query, &mut state.query);
            } else {
                state.query.clear();
                state.query.extend_from_slice(query);
            }
            self.v
                .retrieve(&state.query, 0..state.query.len(), retrieval, &mut state.v);
        }
    }

    /// Stage 2, whole chunk: V placements.
    fn stage_v(
        &self,
        reads: &[ReadState],
        placements: &mut [Vec<Placement>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        let requests = read_requests(reads, |r| 0..r.query.len(), |r| &r.v);
        self.v
            .place_many(&self.aligner, &requests, dp, batch, placements);
    }

    /// Stage 3: per read, the J window right of V and its candidates; then J
    /// placements over the whole chunk.
    fn stage_j(
        &self,
        reads: &mut [ReadState],
        placements: &mut [Vec<Placement>],
        retrieval: &mut RetrievalScratch,
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        for (state, placed) in reads.iter_mut().zip(placements.iter()) {
            state.j_window = self.j_window.window(state.query.len(), placed);
            self.j.retrieve(
                &state.query,
                state.j_window.clone(),
                retrieval,
                &mut state.j,
            );
        }
        let requests = read_requests(reads, |r| r.j_window.clone(), |r| &r.j);
        self.j
            .place_many(&self.aligner, &requests, dp, batch, placements);
    }

    /// Stage 4: D windows from V/J geometry and their candidates per read, D
    /// placements over the whole chunk, then the admission gates per read.
    #[allow(clippy::too_many_arguments)]
    fn stage_d(
        &self,
        d: &SegmentStage<RD>,
        reads: &[ReadState],
        placements: &mut [Vec<Placement>],
        scratch: &mut DScratch,
        retrieval: &mut RetrievalScratch,
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        let DScratch {
            windows,
            read_windows,
            candidates,
            start,
        } = scratch;
        windows.clear();
        for (k, (state, placed)) in reads.iter().zip(placements.iter()).enumerate() {
            read_windows.clear();
            let v_count = placed
                .iter()
                .take_while(|p| p.kind == SegmentKind::V)
                .count();
            let (v, j) = placed.split_at(v_count);
            self.d_window.windows(state.query.len(), v, j, read_windows);
            merge_windows(read_windows);
            windows.extend(read_windows.drain(..).map(|window| (k, window)));
        }
        if candidates.len() < windows.len() {
            candidates.resize_with(windows.len(), Vec::new);
        }
        for ((k, window), found) in windows.iter().zip(candidates.iter_mut()) {
            d.retrieve(&reads[*k].query, window.clone(), retrieval, found);
        }
        let requests: Vec<SegmentRequest<'_>> = windows
            .iter()
            .zip(candidates.iter())
            .map(|((k, window), candidates)| SegmentRequest {
                read: *k,
                query: &reads[*k].query,
                window: window.clone(),
                candidates,
            })
            .collect();
        start.clear();
        start.extend(placements.iter().map(Vec::len));
        d.place_many(&self.aligner, &requests, dp, batch, placements);
        for ((state, placed), &first) in reads.iter().zip(placements.iter_mut()).zip(start.iter()) {
            retain_from(placed, first, |p| {
                self.d_admission
                    .admit(p, &state.query, d.segment.reference(p.primary_allele()))
            });
        }
    }

    /// Whether calls need an anchor: some V or J hit significant on its own.
    fn requires_anchor(&self) -> bool {
        self.v.significance.is_some() || self.j.significance.is_some()
    }

    /// Stage 5, per read: adjudication and refinement into calls. A path
    /// without an anchor ([`Placement::anchor`]) is not a call.
    fn stage_adjudicate(
        &self,
        reads: &[ReadState],
        placements: &mut [Vec<Placement>],
        adjudication: &mut Adj::Scratch,
        out: &mut Vec<Call>,
    ) {
        out.reserve(reads.len());
        for (state, placed) in reads.iter().zip(placements.iter_mut()) {
            let mut path = self
                .adjudicator
                .adjudicate(state.query.len(), placed, adjudication)
                .filter(|path| {
                    !self.requires_anchor()
                        || path.segments.iter().any(|s| placed[s.placement].anchor)
                });
            let marginals = path
                .as_ref()
                .and_then(|path| self.adjudicator.marginals(placed, path, adjudication));
            if let Some(path) = path.as_mut() {
                self.refiner.refine(&state.query, placed, path);
            }
            out.push(Call {
                orientation: state.orientation.expect("stage 1 orients every read"),
                query: state.query.clone(),
                placements: std::mem::take(placed),
                path,
                marginals,
            });
        }
    }
}

impl<A, RV, RJ, RD, W, Adm, Adj, Rf> VdjAssigner for Assigner<A, RV, RJ, RD, W, Adm, Adj, Rf>
where
    A: Aligner + Send + Sync,
    RV: Retriever,
    RJ: Retriever,
    RD: Retriever,
    W: DWindowPolicy,
    Adm: Admission,
    Adj: Adjudicator,
    Rf: Refiner,
{
    type Workspace = Workspace<A::Workspace, Adj::Scratch>;

    fn workspace(&self) -> Self::Workspace {
        Workspace::default()
    }

    fn assign(&self, query: &[u8], ws: &mut Self::Workspace) -> Call {
        let mut calls = std::mem::take(&mut ws.calls);
        calls.clear();
        self.assign_many(&[query], ws, &mut calls);
        let call = calls.pop().expect("one call per query");
        ws.calls = calls;
        call
    }

    fn assign_many<Q: AsRef<[u8]>>(
        &self,
        queries: &[Q],
        ws: &mut Self::Workspace,
        out: &mut Vec<Call>,
    ) {
        let count = queries.len();
        if ws.reads.len() < count {
            ws.reads.resize_with(count, ReadState::default);
            ws.placements.resize_with(count, Vec::new);
        }
        let Workspace {
            retrieval,
            dp,
            batch,
            reads,
            placements,
            d,
            adjudication,
            ..
        } = ws;
        let (reads, placements) = (&mut reads[..count], &mut placements[..count]);
        placements.iter_mut().for_each(Vec::clear);

        self.stage_orient(queries, reads, retrieval);
        self.stage_v(reads, placements, dp, batch);
        self.stage_j(reads, placements, retrieval, dp, batch);
        if let Some(stage) = &self.d {
            self.stage_d(stage, reads, placements, d, retrieval, dp, batch);
        }
        self.stage_adjudicate(reads, placements, adjudication, out);
    }
}

impl<A, RV, RJ, RD, W, Adm, Adj, Rf> AlleleCatalog for Assigner<A, RV, RJ, RD, W, Adm, Adj, Rf> {
    fn segment(&self, kind: SegmentKind) -> Option<&PreparedSegment> {
        match kind {
            SegmentKind::V => Some(&self.v.segment),
            SegmentKind::D => self.d.as_ref().map(|stage| &*stage.segment),
            SegmentKind::J => Some(&self.j.segment),
        }
    }
}
