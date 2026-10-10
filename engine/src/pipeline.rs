//! The assignment pipeline.
//!
//! ```text
//!  chunk of reads ─► per read: orient (germline k-mer support), V candidates
//!                 ─► V stage over the whole chunk     score batch + traceback waves
//!                 ─► per read: J window right of V; J stage over the whole chunk
//!                 ─► per read: D windows from V/J geometry, D candidates
//!                 ─► D stage over the whole chunk ─► per read: admission gates
//!                 ─► per read: adjudicator (placement DAG) ─► Call
//!                 ─► C stage over the whole chunk (optional; reads finished
//!                    calls only, see [`crate::constant`])
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
use crate::constant::{CScratch, ConstantCall, ConstantStage};
use crate::germline::SegmentKind;
use crate::locus::{Loci, Locus};
use crate::orient::{Orientation, OrientationPolicy, StrandIndex};
use crate::placement::Placement;
use crate::retrieve::{
    Candidate, ExhaustiveRetriever, LmerMaskRetriever, RetrievalScratch, Retriever,
};
use crate::segment::{BatchScratch, PreparedSegment, SegmentRequest, SegmentStage};
use crate::window::{
    DWindowPolicy, DownstreamWindow, ExpandedResidual, WithAnchorWindow, merge_windows,
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
    /// C placements right of a J-final path, clipped to the J exit; empty
    /// without a C stage. Never part of `placements`, so C cannot change
    /// the V(D)J call.
    pub constant: Vec<Placement>,
    /// The C call, when a C stage finds a significant one.
    pub c: Option<ConstantCall>,
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
    /// Adjudicators for reads of one locus, where it differs from
    /// `adjudicator`: a D database of several loci prices a D by its own
    /// locus' D genes ([`crate::config::DPrior`]).
    pub locus_adjudicators: Vec<(Locus, Adj)>,
    /// Constant-region assignment after adjudication; `None` without a C
    /// database ([`Assigner::with_constant`]).
    pub c: Option<ConstantStage>,
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
    c: CScratch,
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

/// The assigner with one component replaced: the other fields move across.
macro_rules! replace {
    ($assigner:ident { $($kept:ident),* } $field:ident: $value:expr) => {
        Assigner {
            $field: $value,
            $($kept: $assigner.$kept),*
        }
    };
}

/// Swapping a component changes the assigner's type, so the whole pipeline
/// monomorphises again around it.
impl<A, RV, RJ, RD, W, Adm, Adj> Assigner<A, RV, RJ, RD, W, Adm, Adj> {
    pub fn with_aligner<A2>(self, aligner: A2) -> Assigner<A2, RV, RJ, RD, W, Adm, Adj> {
        replace!(self { v, j, d, strand, j_window, d_window, d_admission, adjudicator, locus_adjudicators, c, options } aligner: aligner)
    }

    pub fn with_d_window<W2>(self, d_window: W2) -> Assigner<A, RV, RJ, RD, W2, Adm, Adj> {
        replace!(self { v, j, d, strand, j_window, aligner, d_admission, adjudicator, locus_adjudicators, c, options } d_window: d_window)
    }

    pub fn with_d_admission<Adm2>(
        self,
        d_admission: Adm2,
    ) -> Assigner<A, RV, RJ, RD, W, Adm2, Adj> {
        replace!(self { v, j, d, strand, j_window, aligner, d_window, adjudicator, locus_adjudicators, c, options } d_admission: d_admission)
    }

    /// Swap the adjudicator, for reads of every locus. The D stage keeps
    /// the D placements the original one needed: built with posteriors, it
    /// keeps those below the DAG's D floor, so pair an adjudicator without
    /// its own D significance floor (such as
    /// [`crate::adjudicate::ResidualAdjudicator`]) with `posterior: false`.
    pub fn with_adjudicator<Adj2>(
        self,
        adjudicator: Adj2,
    ) -> Assigner<A, RV, RJ, RD, W, Adm, Adj2> {
        let assigner = self;
        Assigner {
            adjudicator,
            locus_adjudicators: Vec::new(),
            v: assigner.v,
            j: assigner.j,
            d: assigner.d,
            strand: assigner.strand,
            j_window: assigner.j_window,
            aligner: assigner.aligner,
            d_window: assigner.d_window,
            d_admission: assigner.d_admission,
            c: assigner.c,
            options: assigner.options,
        }
    }

    pub fn with_options(mut self, options: AssignerOptions) -> Self {
        self.options = options;
        self
    }
}

impl<A, RV, RJ, RD, W, Adm, Adj> Assigner<A, RV, RJ, RD, W, Adm, Adj> {
    /// The loci of a placement: those of its alleles together.
    fn placement_loci(&self, placement: &Placement) -> Loci {
        let Some(segment) = self.segment(placement.kind) else {
            return Loci::ANY;
        };
        let set = segment.set();
        placement
            .alleles
            .iter()
            .fold(Loci::NONE, |loci, &a| loci.union(set.allele(a).loci))
    }

    /// The loci of `placements` together; [`Loci::ANY`] when there are none.
    fn placed_loci(&self, placements: &[Placement]) -> Loci {
        if placements.is_empty() {
            return Loci::ANY;
        }
        placements
            .iter()
            .fold(Loci::NONE, |loci, p| loci.union(self.placement_loci(p)))
    }

    /// The adjudicator for reads of `locus` (`None`: unknown or several).
    fn adjudicator_for(&self, locus: Option<Locus>) -> &Adj {
        locus
            .and_then(|l| self.locus_adjudicators.iter().find(|(at, _)| *at == l))
            .map_or(&self.adjudicator, |(_, adjudicator)| adjudicator)
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
impl<A, RV, RJ, RD, W, Adm, Adj> Assigner<A, RV, RJ, RD, W, Adm, Adj>
where
    A: Aligner + Send + Sync,
    RV: Retriever,
    RJ: Retriever,
    RD: Retriever,
    W: DWindowPolicy,
    Adm: Admission,
    Adj: Adjudicator,
{
    /// Stage 1, per read: orientation, oriented query, V candidates.
    ///
    /// The strand comes from germline k-mer support. When that is not
    /// significant, V is retrieved on both strands and the strongest seed
    /// cluster decides; the chosen strand's candidates are kept.
    fn stage_orient<Q: AsRef<[u8]>>(
        &self,
        queries: &[Q],
        reads: &mut [ReadState],
        retrieval: &mut RetrievalScratch,
    ) {
        let policy = self.options.orientation;
        let strongest =
            |candidates: &[Candidate]| candidates.iter().map(|c| c.votes).max().unwrap_or(0);
        for (state, query) in reads.iter_mut().zip(queries) {
            let query = query.as_ref();
            let whole = 0..query.len();
            let decided = if policy.try_reverse {
                let support = self.strand.support(query);
                support
                    .decisive()
                    .then(|| policy.choose(support.forward, support.reverse))
            } else {
                Some(Orientation::Forward)
            };
            let orientation = match decided {
                Some(orientation) => {
                    if orientation.is_reverse() {
                        reverse_complement_into(query, &mut state.query);
                    } else {
                        state.query.clear();
                        state.query.extend_from_slice(query);
                    }
                    self.v
                        .retrieve(&state.query, whole, retrieval, &mut state.v);
                    orientation
                }
                None => {
                    // `state.j` is free until stage 3: it holds the reverse strand's candidates.
                    self.v
                        .retrieve(query, whole.clone(), retrieval, &mut state.v);
                    reverse_complement_into(query, &mut state.query);
                    self.v
                        .retrieve(&state.query, whole, retrieval, &mut state.j);
                    let orientation = policy.choose(strongest(&state.v), strongest(&state.j));
                    if orientation.is_reverse() {
                        std::mem::swap(&mut state.v, &mut state.j);
                    } else {
                        state.query.clear();
                        state.query.extend_from_slice(query);
                    }
                    orientation
                }
            };
            state.orientation = Some(orientation);
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
        let set = self.j.segment.set();
        for (state, placed) in reads.iter_mut().zip(placements.iter()) {
            state.j_window = self.j_window.window(state.query.len(), placed);
            self.j.retrieve(
                &state.query,
                state.j_window.clone(),
                retrieval,
                &mut state.j,
            );
            // A J joins a V of its own locus.
            let loci = self.placed_loci(placed);
            state
                .j
                .retain(|c| set.allele(c.allele).loci.pairs_with(loci));
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
        let set = d.segment.set();
        for ((k, window), found) in windows.iter().zip(candidates.iter_mut()) {
            d.retrieve(&reads[*k].query, window.clone(), retrieval, found);
            // A D joins the V and J of its own locus.
            let loci = self.placed_loci(&placements[*k]);
            found.retain(|c| set.allele(c.allele).loci.pairs_with(loci));
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

    /// Stage 5, per read: adjudication into calls (see
    /// [`Self::adjudicate_loci`] for reads whose placements span loci).
    fn stage_adjudicate(
        &self,
        reads: &[ReadState],
        placements: &mut [Vec<Placement>],
        adjudication: &mut Adj::Scratch,
        out: &mut Vec<Call>,
    ) {
        out.reserve(reads.len());
        for (state, placed) in reads.iter().zip(placements.iter_mut()) {
            let (path, marginals) = self.adjudicate_loci(state.query.len(), placed, adjudication);
            out.push(Call {
                orientation: state.orientation.expect("stage 1 orients every read"),
                query: state.query.clone(),
                placements: std::mem::take(placed),
                path,
                marginals,
                constant: Vec::new(),
                c: None,
            });
        }
    }

    /// One read's path and marginals, each path within one locus. When
    /// every placement shares a locus, the read is adjudicated once (by
    /// that locus' adjudicator). Otherwise each locus is adjudicated over
    /// the placements that belong to it, and the call is the path of the
    /// locus with the largest likelihood `Z` (or, without posteriors, the
    /// best path): the posterior is the mixture over loci, weighted by `Z`,
    /// as each locus' paths are hypotheses of one model.
    fn adjudicate_loci(
        &self,
        query_len: usize,
        placed: &[Placement],
        scratch: &mut Adj::Scratch,
    ) -> (Option<CallPath>, Option<Marginals>) {
        let loci: Vec<Loci> = placed.iter().map(|p| self.placement_loci(p)).collect();
        let shared = loci.iter().fold(Loci::ANY, |all, &l| all.intersection(l));
        if !shared.is_empty() {
            let adjudicator = self.adjudicator_for(shared.single());
            return self.adjudicate_read(adjudicator, query_len, placed, scratch);
        }
        let all = loci.iter().fold(Loci::NONE, |all, &l| all.union(l));
        // Per locus: the placement ids, the path and the marginals.
        let mut found: Vec<(Vec<usize>, CallPath, Option<Marginals>)> = Vec::new();
        for locus in all.iter() {
            let ids: Vec<usize> = (0..placed.len())
                .filter(|&i| loci[i].contains(locus))
                .collect();
            if ids.iter().all(|&i| placed[i].kind == SegmentKind::D) {
                continue;
            }
            let subset: Vec<Placement> = ids.iter().map(|&i| placed[i].clone()).collect();
            let adjudicator = self.adjudicator_for(Some(locus));
            if let (Some(path), marginals) =
                self.adjudicate_read(adjudicator, query_len, &subset, scratch)
            {
                found.push((ids, path, marginals));
            }
        }
        let posterior = found.iter().all(|(_, _, m)| m.is_some());
        let best = if posterior {
            found.iter().enumerate().max_by(|(_, a), (_, b)| {
                let z = |m: &Option<Marginals>| m.as_ref().map_or(0.0, |m| m.log_partition);
                z(&a.2).total_cmp(&z(&b.2))
            })
        } else {
            found
                .iter()
                .enumerate()
                .max_by_key(|(_, (_, path, _))| path.score)
        };
        let Some((best, _)) = best else {
            return (None, None);
        };
        let marginals = posterior.then(|| {
            let log_z: Vec<f64> = found
                .iter()
                .map(|(_, _, m)| m.as_ref().expect("posterior").log_partition)
                .collect();
            let top = log_z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let total = log_z.iter().map(|z| (z - top).exp()).sum::<f64>();
            let mut mixed = Marginals {
                placement: vec![0.0; placed.len()],
                d_count: Vec::new(),
                log_partition: top + total.ln(),
            };
            for ((ids, _, m), z) in found.iter().zip(&log_z) {
                let m = m.as_ref().expect("posterior");
                let weight = (z - top).exp() / total;
                for (&id, p) in ids.iter().zip(&m.placement) {
                    mixed.placement[id] += weight * p;
                }
                if mixed.d_count.len() < m.d_count.len() {
                    mixed.d_count.resize(m.d_count.len(), 0.0);
                }
                for (sum, p) in mixed.d_count.iter_mut().zip(&m.d_count) {
                    *sum += weight * p;
                }
            }
            mixed
        });
        let (ids, mut path, _) = found.swap_remove(best);
        for segment in &mut path.segments {
            segment.placement = ids[segment.placement];
        }
        (Some(path), marginals)
    }

    /// A path without an anchor ([`Placement::anchor`]) is not a call. A
    /// path the adjudicator decodes from the posterior replaces the call
    /// when it keeps an anchor.
    fn adjudicate_read(
        &self,
        adjudicator: &Adj,
        query_len: usize,
        placed: &[Placement],
        scratch: &mut Adj::Scratch,
    ) -> (Option<CallPath>, Option<Marginals>) {
        let anchored = |path: &CallPath| {
            !self.requires_anchor() || path.segments.iter().any(|s| placed[s.placement].anchor)
        };
        let mut path = adjudicator
            .adjudicate(query_len, placed, scratch)
            .filter(anchored);
        let marginals = path
            .as_ref()
            .and_then(|path| adjudicator.marginals(placed, path, scratch));
        if let (Some(called), Some(m)) = (&path, &marginals)
            && let Some(decoded) = adjudicator.decode(placed, called, m, scratch)
            && anchored(&decoded)
        {
            path = Some(decoded);
        }
        (path, marginals)
    }
}

impl<A, RV, RJ, RD, W, Adm, Adj> VdjAssigner for Assigner<A, RV, RJ, RD, W, Adm, Adj>
where
    A: Aligner + Send + Sync,
    RV: Retriever,
    RJ: Retriever,
    RD: Retriever,
    W: DWindowPolicy,
    Adm: Admission,
    Adj: Adjudicator,
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
            c,
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
        let first = out.len();
        self.stage_adjudicate(reads, placements, adjudication, out);
        if let Some(stage) = &self.c {
            // C follows the J of its own locus.
            let loci: Vec<Loci> = out[first..]
                .iter()
                .map(|call| {
                    let last = call.path.as_ref().and_then(|p| p.segments.last());
                    last.map_or(Loci::ANY, |s| {
                        self.placement_loci(&call.placements[s.placement])
                    })
                })
                .collect();
            stage.assign(
                &self.aligner,
                &mut out[first..],
                &loci,
                c,
                retrieval,
                dp,
                batch,
            );
        }
    }
}

impl<A, RV, RJ, RD, W, Adm, Adj> AlleleCatalog for Assigner<A, RV, RJ, RD, W, Adm, Adj> {
    fn segment(&self, kind: SegmentKind) -> Option<&PreparedSegment> {
        match kind {
            SegmentKind::V => Some(&self.v.segment),
            SegmentKind::D => self.d.as_ref().map(|stage| &*stage.segment),
            SegmentKind::J => Some(&self.j.segment),
            SegmentKind::C => self.c.as_ref().map(|stage| &*stage.plain),
        }
    }
}
