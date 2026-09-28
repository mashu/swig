//! The assignment pipeline.
//!
//! ```text
//!  chunk of reads ─► per read: orient (V votes, fwd vs rc), retrieve V and J candidates
//!                 ─► V stage over the whole chunk ─┐  one score batch + one traceback batch
//!                 ─► J stage over the whole chunk ─┤
//!                 ─► per read: D windows from V/J geometry, D candidates
//!                 ─► D stage over the whole chunk ─► per read: admission gates
//!                 ─► per read: adjudicator (placement DAG) ─► refiner ─► Call
//! ```
//!
//! **Trick — stage-wise batching across reads.** After prefiltering a read has
//! only ~3 V, ~3 J and ~9 D candidates, far too few to fill a 16-lane SIMD
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

use crate::adjudicate::{Adjudicator, CallPath, DagScratch, PlacementDag};
use crate::admission::{Admission, DefaultAdmission};
use crate::align::{Aligner, InterSequenceGotoh, InterSequenceWorkspace};
use crate::alphabet::reverse_complement_into;
use crate::config::DefaultVjRetriever;
use crate::evidence::Evidence;
use crate::germline::SegmentKind;
use crate::orient::{Orientation, OrientationPolicy};
use crate::placement::{mark_repeated_d_tracts, Placement};
use crate::refine::{extend_junction_tie, NoRefine, Refiner, JUNCTION_TIE_IDENTITY};
use crate::retrieve::{strongest, Candidate, LmerMaskRetriever, RetrievalScratch, Retriever};
use crate::segment::{BatchScratch, PreparedSegment, SegmentRequest, SegmentStage};
use crate::window::{merge_windows, DWindowPolicy, ExpandedResidual, WithAnchorWindow};

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
    /// R-optimized / Sensitive-D: write a zero-score five-base block onto a
    /// high-SHM V 3′ end or J 5′ start after the joint partition is chosen.
    pub junction_tie: bool,
}

pub struct Assigner<
    A = InterSequenceGotoh,
    RV = DefaultVjRetriever,
    RJ = RV,
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
}

/// Per-thread buffers of an [`Assigner`]; every buffer is reused across chunks.
#[derive(Debug, Default)]
pub struct Workspace<W = InterSequenceWorkspace, S = DagScratch> {
    retrieval: RetrievalScratch,
    dp: W,
    batch: BatchScratch,
    reverse: Vec<u8>,
    chosen: Vec<Candidate>,
    reverse_candidates: Vec<Candidate>,
    spare_candidates: Vec<Candidate>,
    reads: Vec<ReadState>,
    placements: Vec<Vec<Placement>>,
    v_count: Vec<usize>,
    d_start: Vec<usize>,
    read_windows: Vec<Range<usize>>,
    d_windows: Vec<(usize, Range<usize>)>,
    d_candidates: Vec<Vec<Candidate>>,
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

    /// High-SHM junction repair. V identity is measured on the called span,
    /// after terminal extension and clipping. A zero-score five-base block is
    /// written only into the gap before the selected D, or before the opposite
    /// segment when D is absent.
    fn apply_junction_tie(
        &self,
        query: &[u8],
        placements: &mut [Placement],
        path: &mut crate::adjudicate::CallPath,
    ) {
        use crate::germline::SegmentKind;
        let Some(v_at) = path
            .segments
            .iter()
            .position(|s| placements[s.placement].kind == SegmentKind::V)
        else {
            return;
        };
        let Some(j_at) = path
            .segments
            .iter()
            .position(|s| placements[s.placement].kind == SegmentKind::J)
        else {
            return;
        };
        let v_index = path.segments[v_at].placement;
        let v_reference = self
            .v
            .segment
            .set()
            .allele(placements[v_index].primary_allele())
            .codes
            .clone();
        let v_profile = self
            .v
            .segment
            .profile(placements[v_index].primary_allele())
            .clone();
        let identity = Evidence::from_columns(
            placements[v_index].clip.columns(&path.segments[v_at].span),
            query,
            &v_reference,
        )
        .identity();
        if identity > JUNCTION_TIE_IDENTITY {
            return;
        }
        let d_span = path
            .segments
            .iter()
            .find(|s| placements[s.placement].kind == SegmentKind::D)
            .map(|s| s.span);
        let v_boundary = d_span
            .map(|s| s.q_lo)
            .unwrap_or(path.segments[j_at].span.q_lo);
        extend_junction_tie(
            &mut placements[v_index],
            &mut path.segments[v_at],
            query,
            &v_reference,
            &v_profile,
            self.v.segment.gaps(),
            true,
            v_boundary,
        );
        let j_index = path.segments[j_at].placement;
        let j_reference = self
            .j
            .segment
            .set()
            .allele(placements[j_index].primary_allele())
            .codes
            .clone();
        let j_profile = self
            .j
            .segment
            .profile(placements[j_index].primary_allele())
            .clone();
        let j_boundary = d_span
            .map(|s| s.q_hi)
            .unwrap_or(path.segments[v_at].span.q_hi);
        extend_junction_tie(
            &mut placements[j_index],
            &mut path.segments[j_at],
            query,
            &j_reference,
            &j_profile,
            self.j.segment.gaps(),
            false,
            j_boundary,
        );
    }
}

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
    /// Decide the strand from V retrieval strength (J breaks exact ties).
    /// On return `chosen` holds the chosen strand's V candidates and, for
    /// `Reverse`, `reverse` holds the reverse complement.
    fn orient(
        &self,
        query: &[u8],
        retrieval: &mut RetrievalScratch,
        chosen: &mut Vec<Candidate>,
        reverse: &mut Vec<u8>,
        reverse_candidates: &mut Vec<Candidate>,
        spare: &mut Vec<Candidate>,
    ) -> Orientation {
        let n = query.len();
        self.v.retrieve(query, 0..n, retrieval, chosen);
        let forward = strongest(chosen);
        if !self.options.orientation.wants_reverse(forward) {
            return Orientation::Forward;
        }
        reverse_complement_into(query, reverse);
        self.v
            .retrieve(reverse, 0..n, retrieval, reverse_candidates);
        let backward = strongest(reverse_candidates);
        let orientation = if forward == backward {
            self.j.retrieve(query, 0..n, retrieval, spare);
            let j_forward = strongest(spare);
            self.j.retrieve(reverse, 0..n, retrieval, spare);
            OrientationPolicy::choose(j_forward, strongest(spare))
        } else {
            OrientationPolicy::choose(forward, backward)
        };
        if orientation.is_reverse() {
            std::mem::swap(chosen, reverse_candidates);
        }
        orientation
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
        let Workspace {
            retrieval,
            dp,
            batch,
            reverse,
            chosen,
            reverse_candidates,
            spare_candidates,
            reads,
            placements,
            v_count,
            d_start,
            read_windows,
            d_windows,
            d_candidates,
            adjudication,
            ..
        } = ws;
        if reads.len() < count {
            reads.resize_with(count, ReadState::default);
            placements.resize_with(count, Vec::new);
        }
        v_count.resize(count, 0);
        d_start.resize(count, 0);
        let (reads, placements) = (&mut reads[..count], &mut placements[..count]);

        // Stage 1, per read: orientation, V and J candidates.
        for (state, query) in reads.iter_mut().zip(queries) {
            let query = query.as_ref();
            let orientation = self.orient(
                query,
                retrieval,
                chosen,
                reverse,
                reverse_candidates,
                spare_candidates,
            );
            state.orientation = Some(orientation);
            state.query.clear();
            state.query.extend_from_slice(if orientation.is_reverse() {
                reverse
            } else {
                query
            });
            std::mem::swap(&mut state.v, chosen);
            self.j
                .retrieve(&state.query, 0..state.query.len(), retrieval, &mut state.j);
        }
        placements.iter_mut().for_each(Vec::clear);

        // Stages 2 and 3, whole chunk: V, then J.
        let requests: Vec<SegmentRequest<'_>> = reads
            .iter()
            .enumerate()
            .map(|(k, r)| SegmentRequest {
                read: k,
                query: &r.query,
                window: 0..r.query.len(),
                candidates: &r.v,
            })
            .collect();
        self.v
            .place_many(&self.aligner, &requests, dp, batch, placements);
        for (count, placed) in v_count.iter_mut().zip(placements.iter()) {
            *count = placed.len();
        }
        let requests: Vec<SegmentRequest<'_>> = reads
            .iter()
            .enumerate()
            .map(|(k, r)| SegmentRequest {
                read: k,
                query: &r.query,
                window: 0..r.query.len(),
                candidates: &r.j,
            })
            .collect();
        self.j
            .place_many(&self.aligner, &requests, dp, batch, placements);

        // Stage 4: D windows per read, D over the whole chunk, admission per read.
        if let Some(d) = &self.d {
            d_windows.clear();
            for (k, (state, placed)) in reads.iter().zip(placements.iter()).enumerate() {
                read_windows.clear();
                let (v, j) = placed.split_at(v_count[k]);
                self.d_window.windows(state.query.len(), v, j, read_windows);
                merge_windows(read_windows);
                d_windows.extend(read_windows.drain(..).map(|window| (k, window)));
            }
            if d_candidates.len() < d_windows.len() {
                d_candidates.resize_with(d_windows.len(), Vec::new);
            }
            for ((k, window), candidates) in d_windows.iter().zip(d_candidates.iter_mut()) {
                d.retrieve(&reads[*k].query, window.clone(), retrieval, candidates);
            }
            let requests: Vec<SegmentRequest<'_>> = d_windows
                .iter()
                .zip(d_candidates.iter())
                .map(|((k, window), candidates)| SegmentRequest {
                    read: *k,
                    query: &reads[*k].query,
                    window: window.clone(),
                    candidates,
                })
                .collect();
            for (start, placed) in d_start.iter_mut().zip(placements.iter()) {
                *start = placed.len();
            }
            d.place_many(&self.aligner, &requests, dp, batch, placements);
            let d_references: Vec<&[u8]> = d
                .segment
                .set()
                .alleles()
                .iter()
                .map(|allele| allele.codes.as_slice())
                .collect();
            for (k, placed) in placements.iter_mut().enumerate() {
                let query = &reads[k].query;
                retain_from(placed, d_start[k], |p| {
                    let reference = &d.segment.set().allele(p.primary_allele()).codes;
                    self.d_admission.admit(p, query, reference)
                });
                mark_repeated_d_tracts(placed, query, &d_references);
            }
        }

        // Stage 5, per read: adjudication and refinement.
        out.reserve(count);
        for (state, placed) in reads.iter().zip(placements.iter_mut()) {
            let mut path = self
                .adjudicator
                .adjudicate(state.query.len(), placed, adjudication);
            if let Some(path) = path.as_mut() {
                if self.options.junction_tie {
                    self.apply_junction_tie(&state.query, placed, path);
                }
                self.refiner.refine(&state.query, placed, path);
            }
            out.push(Call {
                orientation: state.orientation.expect("stage 1 orients every read"),
                query: state.query.clone(),
                placements: std::mem::take(placed),
                path,
            });
        }
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
