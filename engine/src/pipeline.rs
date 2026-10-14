//! The assignment pipeline.
//!
//! ```text
//!  chunk of reads ─► per read: orient (germline k-mer support), V candidates
//!                 ─► V stage over the whole chunk     score batch + traceback waves
//!                 ─► per read: J window right of V; J stage over the whole chunk
//!                 ─► per read: D windows from V/J geometry, D candidates
//!                 ─► D stage over the whole chunk (admission gates as acceptance)
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
//!     .with_d_window(ExpandedResidual::default())
//!     .with_d_admission(ExactRunFloor { min_run: 5 });
//! ```

use std::ops::Range;

use crate::adjudicate::{Adjudicator, CallPath, DagScratch, Marginals, PlacementDag, ReadContext};
use crate::admission::{Admission, DefaultAdmission};
use crate::align::{Aligner, InterSequenceGotoh, InterSequenceWorkspace};
use crate::alphabet::reverse_complement_into;
use crate::config::DefaultVRetriever;
use crate::constant::{CScratch, ConstantCall, ConstantStage};
use crate::germline::SegmentKind;
use crate::locus::{Loci, Locus};
use crate::orient::{Orientation, OrientationPolicy, StrandIndex};
use crate::placement::Placement;
use crate::report::reported_posterior;
use crate::retrieve::{
    Candidate, ExhaustiveRetriever, LmerMaskRetriever, RetrievalScratch, Retriever,
};
use crate::segment::{self, BatchScratch, PreparedSegment, SegmentRequest, SegmentStage};
use crate::window::{
    DWindowPolicy, DownstreamWindow, ExpandedResidual, WithAnchorWindow, bounding_v, merge_windows,
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

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AssignerOptions {
    pub orientation: OrientationPolicy,
    /// V, D and J hits whose local alignment has an identity (matches over
    /// alignment columns) below this are dropped before they shape windows
    /// or enter adjudication; 0 keeps every hit. The outer flank a stage
    /// stitches on afterwards ([`crate::flank`]) does not count: it scores
    /// at most 0, so at `(2, −3)` its identity is at most 3/5, and it would
    /// drag a hypermutated V or J below a floor near 0.6 that its own
    /// alignment clears.
    pub min_identity: f64,
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
    /// `None` for loci without D (light chains). Its acceptance predicate
    /// is the D admission gate, so an inadmissible hit takes no slot.
    pub d: Option<SegmentStage<RD, Adm>>,
    /// Germline k-mer presence for orientation.
    pub strand: StrandIndex,
    /// Where J is searched relative to V.
    pub j_window: DownstreamWindow,
    pub aligner: A,
    pub d_window: W,
    pub adjudicator: Adj,
    /// Adjudicators for reads of one locus, where it differs from
    /// `adjudicator`: a D or J database of several loci prices a D and a J
    /// by its own locus' genes ([`crate::config::DPrior`]).
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
    /// The candidates of each J request, as ranges of `j`.
    j_requests: Vec<Range<usize>>,
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
        replace!(self { v, j, d, strand, j_window, d_window, adjudicator, locus_adjudicators, c, options } aligner: aligner)
    }

    pub fn with_d_window<W2>(self, d_window: W2) -> Assigner<A, RV, RJ, RD, W2, Adm, Adj> {
        replace!(self { v, j, d, strand, j_window, aligner, adjudicator, locus_adjudicators, c, options } d_window: d_window)
    }

    /// Swap the D stage's admission gate (its acceptance predicate).
    pub fn with_d_admission<Adm2>(
        self,
        d_admission: Adm2,
    ) -> Assigner<A, RV, RJ, RD, W, Adm2, Adj> {
        let d = self.d.map(|stage| stage.with_accept(d_admission));
        replace!(self { v, j, strand, j_window, aligner, d_window, adjudicator, locus_adjudicators, c, options } d: d)
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

    /// Drop the hits of `placed[from..]` whose local alignment is below the
    /// identity floor.
    fn drop_dissimilar(&self, placed: &mut Vec<Placement>, from: usize) {
        let floor = self.options.min_identity;
        if floor > 0.0 {
            retain_from(placed, from, |p| {
                p.local_evidence.identity() + 1e-9 >= floor
            });
        }
    }

    /// The loci of the anchors among `placements` together; [`Loci::ANY`]
    /// when none is an anchor. A chance hit does not choose the locus of
    /// the next segment, as it does not bound the J window
    /// ([`DownstreamWindow`]).
    fn anchored_loci(&self, placements: &[Placement]) -> Loci {
        placements
            .iter()
            .filter(|p| p.anchor.is_some())
            .map(|p| self.placement_loci(p))
            .reduce(Loci::union)
            .unwrap_or(Loci::ANY)
    }

    /// The adjudicator for reads of `locus` (`None`: unknown or several).
    fn adjudicator_for(&self, locus: Option<Locus>) -> &Adj {
        locus
            .and_then(|l| self.locus_adjudicators.iter().find(|(at, _)| *at == l))
            .map_or(&self.adjudicator, |(_, adjudicator)| adjudicator)
    }
}

impl<A, RV, RJ, RD, W: DWindowPolicy, Adm, Adj> Assigner<A, RV, RJ, RD, W, Adm, Adj> {
    /// Set `out` to the merged D windows of a read of `query_len` bases
    /// with `placements`, its V placements followed by its J placements
    /// (each kind best first; any D placements after them are ignored).
    /// The windows lie between the J and the V that bound the D
    /// ([`bounding_v`]).
    ///
    /// On a read holding two chains (an scFv's heavy and light chain), the
    /// best V and the first J may belong to different chains, and the
    /// window between them holds neither chain's D. Each locus of the
    /// read's anchors that has D genes then bounds a window of its own, by
    /// its own V and J.
    pub fn d_windows(
        &self,
        query_len: usize,
        placements: &[Placement],
        out: &mut Vec<Range<usize>>,
    ) {
        let v_count = placements
            .iter()
            .take_while(|p| p.kind == SegmentKind::V)
            .count();
        let (v, rest) = placements.split_at(v_count);
        let j_count = rest.iter().take_while(|p| p.kind == SegmentKind::J).count();
        let j = &rest[..j_count];
        let anchoring = self.v.significance.is_some();
        out.clear();
        match self.chains(v.iter().chain(j)) {
            None => self
                .d_window
                .windows(query_len, bounding_v(v, anchoring), j.first(), out),
            Some(loci) => {
                for locus in loci.iter().filter(|l| l.has_d()) {
                    let of = |p: &&Placement| self.placement_loci(p).contains(locus);
                    let v = bounding_v(v.iter().filter(of), anchoring);
                    self.d_window.windows(query_len, v, j.iter().find(of), out);
                }
            }
        }
        merge_windows(out);
    }

    /// The loci of a read whose anchors belong to more than one chain: no
    /// locus is shared by all of them, as an scFv's heavy and light V share
    /// none. `None` for a read of one chain (TRAV/DV alleles and a TRAJ
    /// share TRA) or with no anchor.
    fn chains<'a>(&self, placements: impl Iterator<Item = &'a Placement>) -> Option<Loci> {
        let mut anchors = placements
            .filter(|p| p.anchor.is_some())
            .map(|p| self.placement_loci(p));
        let first = anchors.next()?;
        let (shared, all) = anchors.fold((first, first), |(shared, all), loci| {
            (shared.intersection(loci), all.union(loci))
        });
        shared.is_empty().then_some(all)
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
            let decided = match policy.only() {
                Some(only) => Some(only),
                None => {
                    let support = self.strand.support(query);
                    support
                        .decisive()
                        .then(|| policy.choose(support.forward, support.reverse))
                }
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
        for placed in placements.iter_mut() {
            self.drop_dissimilar(placed, 0);
        }
    }

    /// Stage 3: per read, the J window right of V and its candidates; then J
    /// placements over the whole chunk.
    ///
    /// A J joins an anchoring V of its own locus. When the anchors choose
    /// several loci (a read holding two chains, as an scFv does), each locus
    /// is a request of its own: in one request, the J alleles of one chain
    /// could fill the stage's keep and crowd out the other chain's J. An
    /// allele of several of these loci is searched with the first.
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
            let loci = self.anchored_loci(placed);
            state
                .j
                .retain(|c| set.allele(c.allele).loci.pairs_with(loci));
            state.j_requests.clear();
            if loci == Loci::ANY || loci.single().is_some() {
                state.j_requests.push(0..state.j.len());
                continue;
            }
            let locus = |c: &Candidate| set.allele(c.allele).loci.intersection(loci).first();
            state.j.sort_by_key(locus);
            let mut at = 0;
            for group in state.j.chunk_by(|a, b| locus(a) == locus(b)) {
                state.j_requests.push(at..at + group.len());
                at += group.len();
            }
        }
        let first: Vec<usize> = placements.iter().map(Vec::len).collect();
        let requests: Vec<SegmentRequest<'_>> = reads
            .iter()
            .enumerate()
            .flat_map(|(read, state)| {
                state.j_requests.iter().map(move |range| SegmentRequest {
                    read,
                    query: &state.query,
                    window: state.j_window.clone(),
                    candidates: &state.j[range.clone()],
                })
            })
            .collect();
        self.j
            .place_many(&self.aligner, &requests, dp, batch, placements);
        for ((placed, from), state) in placements.iter_mut().zip(first).zip(reads.iter()) {
            self.drop_dissimilar(placed, from);
            // Each request's placements are ranked best first, as the stage
            // ranks them; so are the read's.
            if state.j_requests.len() > 1 {
                placed[from..].sort_by_key(segment::rank);
            }
        }
    }

    /// Stage 4: D windows from V/J geometry and their candidates per read,
    /// then D placements over the whole chunk (the stage applies the
    /// admission gates).
    #[allow(clippy::too_many_arguments)]
    fn stage_d(
        &self,
        d: &SegmentStage<RD, Adm>,
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
            self.d_windows(state.query.len(), placed, read_windows);
            windows.extend(read_windows.drain(..).map(|window| (k, window)));
        }
        if candidates.len() < windows.len() {
            candidates.resize_with(windows.len(), Vec::new);
        }
        let set = d.segment.set();
        for ((k, window), found) in windows.iter().zip(candidates.iter_mut()) {
            d.retrieve(&reads[*k].query, window.clone(), retrieval, found);
            // A D joins the anchoring V and J of its own locus.
            let loci = self.anchored_loci(&placements[*k]);
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
        for (placed, &first) in placements.iter_mut().zip(start.iter()) {
            self.drop_dissimilar(placed, first);
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
            let (path, marginals) = self.adjudicate_loci(&state.query, placed, adjudication);
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
    /// the placements that belong to it, in the context of all the read's
    /// placements (their divergence and anchored columns, and the read's
    /// length: [`ReadContext`]), so that every locus' paths are hypotheses
    /// of one model. The call is the best of their paths in the model the call
    /// maximises ([`CallPath::value`]), with or without posteriors; an exact
    /// tie goes to the path whose last segment comes first in the DAG's
    /// order (by kind, then query start and end), as it would within one
    /// locus. The posterior is the mixture over loci, weighted by their
    /// likelihoods `Z`.
    fn adjudicate_loci(
        &self,
        query: &[u8],
        placed: &[Placement],
        scratch: &mut Adj::Scratch,
    ) -> (Option<CallPath>, Option<Marginals>) {
        let loci: Vec<Loci> = placed.iter().map(|p| self.placement_loci(p)).collect();
        let shared = loci.iter().fold(Loci::ANY, |all, &l| all.intersection(l));
        let context = ReadContext::of(placed, query.len() as u32);
        if !shared.is_empty() {
            let adjudicator = self.adjudicator_for(shared.single());
            return self.adjudicate_read(adjudicator, query, placed, &context, scratch);
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
                self.adjudicate_read(adjudicator, query, &subset, &context, scratch)
            {
                found.push((ids, path, marginals));
            }
        }
        let last = |(ids, path, _): &(Vec<usize>, CallPath, _)| {
            let segment = path.segments.last().expect("a path has a segment");
            let placement = &placed[ids[segment.placement]];
            let order = (
                placement.kind.index(),
                placement.q_start(),
                placement.q_end(),
            );
            std::cmp::Reverse(order)
        };
        let best = (0..found.len()).reduce(|best, k| {
            let key = |k: usize| (found[k].1.value, last(&found[k]));
            if key(k) > key(best) { k } else { best }
        });
        let Some(best) = best else {
            return (None, None);
        };
        let posterior = found.iter().all(|(_, _, m)| m.is_some());
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

    /// The call of `placed` in the read's context `context`
    /// ([`Adjudicator::adjudicate_at`]), and its marginals. A path without
    /// an anchor, a segment whose retained clip is
    /// significant on its own ([`Placement::anchors`]), is not a call. A
    /// path the adjudicator decodes from the posterior replaces the call
    /// when it keeps an anchor and its D is reported at least as probable
    /// ([`Self::keeps_the_d_posterior`]).
    fn adjudicate_read(
        &self,
        adjudicator: &Adj,
        query: &[u8],
        placed: &[Placement],
        context: &ReadContext,
        scratch: &mut Adj::Scratch,
    ) -> (Option<CallPath>, Option<Marginals>) {
        let anchored = |path: &CallPath| {
            !self.requires_anchor()
                || path
                    .segments
                    .iter()
                    .any(|s| placed[s.placement].anchors(&s.span))
        };
        let mut path = adjudicator
            .adjudicate_at(placed, context, scratch)
            .filter(anchored);
        let marginals = path
            .as_ref()
            .and_then(|path| adjudicator.marginals(placed, path, scratch));
        if let (Some(called), Some(m)) = (&path, &marginals)
            && let Some(decoded) = adjudicator.decode(placed, called, m, scratch)
            && anchored(&decoded)
            && self.keeps_the_d_posterior(query, placed, called, &decoded, m)
        {
            path = Some(decoded);
        }
        (path, marginals)
    }

    /// Whether the report rates `decoded`'s D at least as probable as
    /// `called`'s. The decoder maximises one allele's posterior marginal;
    /// the report's probability is that its names include the read's allele
    /// ([`crate::report::reported_posterior`]), which also counts shorter
    /// placements whose germline the first allele carries. A switch the
    /// report rates less probable would contradict the reason for it.
    fn keeps_the_d_posterior(
        &self,
        query: &[u8],
        placed: &[Placement],
        called: &CallPath,
        decoded: &CallPath,
        marginals: &Marginals,
    ) -> bool {
        let Some(d) = self.segment(SegmentKind::D) else {
            return false;
        };
        let posterior = |path: &CallPath| {
            path.segments
                .iter()
                .find(|s| placed[s.placement].kind == SegmentKind::D)
                .map_or(0.0, |s| {
                    reported_posterior(placed, query, d.set(), path, s, marginals)
                })
        };
        posterior(decoded) >= posterior(called)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::admission::AdmitAll;
    use crate::alphabet::encode;
    use crate::annotate::AnnotateOptions;
    use crate::batch::assign_all;
    use crate::config::DefaultAssigner;
    use crate::prepare::prepare;
    use crate::profiles::CallingProfile;
    use crate::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

    /// Swapping a component rebuilds the assigner around it and moves the
    /// others across: the D window and admission of an identical assigner,
    /// swapped in, give the same calls; admitting every D hit keeps at
    /// least as many D placements per read.
    #[test]
    fn swapped_components_keep_the_others() {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 5);
        let assigner = || {
            let database = prepare(
                db.v.clone(),
                Some(db.d.clone()),
                db.j.clone(),
                AnnotateOptions::default(),
            )
            .unwrap();
            let config = CallingProfile::ROptimized.config();
            DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap()
        };
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::mid()).unwrap();
        let queries: Vec<Vec<u8>> = sim
            .reads(60, 4)
            .unwrap()
            .iter()
            .map(|read| encode(&read.sequence))
            .collect();
        let shown = |calls: &[Call]| -> Vec<String> {
            calls
                .iter()
                .map(|c| format!("{:?} {:?}", c.path, c.placements.len()))
                .collect()
        };
        let base = assign_all(&assigner(), &queries, 1);
        let donor = assigner();
        let swapped = assigner()
            .with_d_window(donor.d_window)
            .with_d_admission(donor.d.expect("the test database has D").accept);
        assert_eq!(shown(&assign_all(&swapped, &queries, 1)), shown(&base));
        let d_count = |call: &Call| {
            call.placements
                .iter()
                .filter(|p| p.kind == SegmentKind::D)
                .count()
        };
        let admitted = assign_all(&assigner().with_d_admission(AdmitAll), &queries, 1);
        assert!(base.iter().any(|call| d_count(call) > 0));
        for (all, gated) in admitted.iter().zip(&base) {
            assert!(d_count(all) >= d_count(gated));
        }
    }
}
