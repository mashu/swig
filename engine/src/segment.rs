//! One segment kind's search stage: prepared database + retriever + keep policy.
//!
//! `SegmentStage` is the unit reused for V, D, J and C (C as
//! `SegmentStage<ReachClasses>`, see [`crate::constant`]). For
//! every request (a query window and its retrieved candidates) it runs:
//!
//! ```text
//! candidates ─► score-only DP ─► (V, J) sibling bounds ─► rank by local score
//!            ─► traceback in waves ─► acceptance ─► outer flank extension
//!               until `top` extended scores are kept (+ ties)
//!            ─► keep the top by extended and by local score
//!            ─► rank by extended score ─► merge identical alignments
//! ```
//!
//! The stage's objective is the *extended* score: a V (J) is valued with
//! its co-linear 5′ (3′) flank stitched on, as adjudication values it. The
//! extension never raises a score, so the local score bounds the extended
//! one, and the waves walk the local ranking until it falls below the
//! `top`-th extended score kept. The `top` best local alignments are kept
//! as well: a hit with no room for a flank pays none, so ranked by extended
//! score alone it could push out a stronger hit with a mutated flank.
//!
//! Every request of a batch is processed in the same waves, so a SIMD
//! aligner sees one score batch and one traceback batch per wave for the
//! whole chunk of reads, while every decision stays per request.

use std::ops::Range;
use std::sync::Arc;

use crate::admission::Admission;
use crate::align::{AlignJob, Aligner, Alignment, ScoreHit};
use crate::clip::ClipIndex;
use crate::evidence::Evidence;
use crate::family_bound::{FamilyBoundTable, FramePath};
use crate::flank::{OuterExtend, extend_outer};
use crate::germline::{AlleleId, GermlineSet, SegmentKind};
use crate::placement::{Placement, merge_identical};
use crate::retrieve::{Candidate, RetrievalScratch, Retriever};
use crate::scoring::{AffineScoring, GapCosts, Score, SubstitutionModel, SubstitutionProfile};
use crate::stats::KarlinAltschul;

/// Length and identity floor applied to a local alignment before outer
/// extension: the acceptance predicate ([`Admission`]) of the V, J and C
/// stages.
///
/// A hit needs at least `min_subst` substitution columns (matches plus
/// mismatches) and an identity of at least `min_identity`. `min_identity` of
/// 0 admits every hit that meets the length floor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HitAcceptance {
    pub min_subst: u32,
    pub min_identity: f64,
}

impl HitAcceptance {
    pub const PERMISSIVE: Self = Self {
        min_subst: 0,
        min_identity: 0.0,
    };

    /// Whether [`Self::allows`] admits every hit: no length and no identity
    /// floor. Only then is the `top`-th best local score also the `top`-th
    /// best *accepted* one, which is what cutoff pruning relies on.
    pub fn is_permissive(&self) -> bool {
        self.min_subst == 0 && self.min_identity <= 0.0
    }

    pub fn allows(&self, evidence: &Evidence) -> bool {
        let subst = evidence.matches + evidence.mismatches;
        let columns = evidence.columns;
        subst >= self.min_subst
            && (columns == 0
                || f64::from(evidence.matches) + 1e-9 >= self.min_identity * f64::from(columns))
    }
}

impl Admission for HitAcceptance {
    fn admit(&self, placement: &Placement, _: &[u8], _: &[u8]) -> bool {
        self.allows(&placement.evidence)
    }

    fn admits_all(&self) -> bool {
        self.is_permissive()
    }
}

/// A germline set with per-allele substitution profiles and gap costs.
#[derive(Clone, Debug)]
pub struct PreparedSegment {
    set: GermlineSet,
    profiles: Vec<SubstitutionProfile>,
    gaps: GapCosts,
}

impl PreparedSegment {
    pub fn new<M: SubstitutionModel + ?Sized>(set: GermlineSet, model: &M, gaps: GapCosts) -> Self {
        let profiles = set
            .alleles()
            .iter()
            .map(|a| SubstitutionProfile::build(model, &a.codes))
            .collect();
        Self {
            set,
            profiles,
            gaps,
        }
    }

    pub fn set(&self) -> &GermlineSet {
        &self.set
    }

    pub fn kind(&self) -> SegmentKind {
        self.set.kind()
    }

    #[inline]
    pub fn profile(&self, allele: AlleleId) -> &SubstitutionProfile {
        &self.profiles[allele as usize]
    }

    #[inline]
    pub fn reference(&self, allele: AlleleId) -> &[u8] {
        &self.set.allele(allele).codes
    }

    pub fn gaps(&self) -> GapCosts {
        self.gaps
    }
}

/// Which alignments a stage keeps per request.
///
/// The stage keeps the `top` best *accepted* alignments by extended score
/// (the local score when the stage extends no flank) and every alignment
/// tied with the last of them, and likewise the `top` best by local score.
/// A high score that the stage's acceptance predicate rejects does not take
/// a slot, so traceback proceeds in waves down the local-score ranking
/// until `top` alignments are accepted and the next local score is below
/// the `top`-th extended score (so below the `top`-th local one too), the
/// list is exhausted, or `cap` candidates were traced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KeepPolicy {
    pub top: usize,
    /// Local scores below this are discarded.
    pub min_score: Score,
    /// Upper bound on tracebacks per request.
    pub cap: usize,
}

impl KeepPolicy {
    /// Exclusive end of the next traceback wave in the score-ranked `scores[from..]`.
    ///
    /// While fewer than `top` hits are accepted, the wave takes as many more
    /// candidates as slots are open, plus the tie group of the last one. Once
    /// `top` are accepted, `cutoff` is the `top`-th best kept score and the
    /// wave takes the remaining candidates whose score reaches it: a score
    /// bounds the kept (extended) score from above.
    pub fn next_wave(
        &self,
        scores: &[Score],
        from: usize,
        accepted: usize,
        cutoff: Option<Score>,
    ) -> usize {
        let limit = self.cap.min(scores.len());
        if from >= limit || self.top == 0 {
            return from;
        }
        let mut end = from;
        if let Some(cutoff) = cutoff {
            while end < limit && scores[end] >= cutoff {
                end += 1;
            }
            return end;
        }
        let need = self.top.saturating_sub(accepted);
        while end < limit && (end - from < need || scores[end] == scores[end - 1]) {
            end += 1;
        }
        end
    }
}

/// Order of equal local scores in a request's ranking.
///
/// Ties decide which alleles take the last traceback slots when `cap` binds,
/// and the order of equally scored placements, whose first allele is the
/// primary allele of a merged placement. The order is therefore part of the
/// stage's definition: a key over the candidates, independent of which
/// scores were computed, bounded or traced first.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TieOrder {
    /// Retrieval order.
    #[default]
    Retrieval,
    /// Every gene's *leader* first, then the other candidates, each group in
    /// retrieval order. A gene's leader is its candidate with most votes (the
    /// lowest allele id among equals).
    GeneLeadersFirst,
}

impl TieOrder {
    /// Sort key of candidate `index` of a request whose candidates have the
    /// gene leaders `leaders` (candidate indices).
    #[inline]
    fn key(self, leaders: &[u32], index: u32) -> (bool, u32) {
        match self {
            Self::Retrieval => (false, index),
            Self::GeneLeadersFirst => (leaders[index as usize] != index, index),
        }
    }
}

/// One search of a stage: candidates of one query window, whose placements
/// are appended to `outs[read]` by [`SegmentStage::place_many`].
#[derive(Clone, Debug)]
pub struct SegmentRequest<'a> {
    pub read: usize,
    pub query: &'a [u8],
    pub window: Range<usize>,
    pub candidates: &'a [Candidate],
}

/// Reusable buffers of [`SegmentStage::place_many`].
#[derive(Debug, Default)]
pub struct BatchScratch {
    hits: Vec<Option<ScoreHit>>,
    alignments: Vec<Option<Alignment>>,
    /// Per request: `(candidate index, local score)` ranked by score, then
    /// by the stage's [`TieOrder`]; `Score::MAX` = not scored.
    ranked: Vec<Vec<(u32, Score)>>,
    /// Per request: extended scores of accepted alignments, descending.
    accepted: Vec<Vec<Score>>,
    /// Per request: siblings whose bound fell below the bounded ranking's
    /// cutoff, `(candidate index, upper bound)`, and that cutoff, the
    /// *barrier* the traceback waves stop at until they are settled (see
    /// [`SegmentStage::rank_with_bounds`]).
    deferred: Vec<Vec<(u32, Score)>>,
    barriers: Vec<Option<Score>>,
    /// Per request: accepted placements before ranking.
    pending: Vec<Vec<Placement>>,
    /// `(request, candidate index)` of each job of the current batch.
    owners: Vec<(usize, u32)>,
    scores: Vec<Score>,
    /// Request `r`'s candidates occupy `starts[r]..starts[r + 1]` of the
    /// per-candidate buffers below.
    starts: Vec<usize>,
    /// Per candidate: the candidate index of its gene's leader (see
    /// [`TieOrder::GeneLeadersFirst`]); filled for requests that are scored.
    leaders: Vec<u32>,
    /// Per candidate: exact local score once known (0 = no positive cell).
    exact: Vec<Score>,
    /// Per sibling: the lower bound `L` from the traced paths (0 = none).
    lower: Vec<Score>,
    /// `(request, candidate index)` of the siblings aligned while ranking.
    early: Vec<(usize, u32)>,
    /// Per request: where its next traceback wave starts, or [`DONE`].
    cursors: Vec<usize>,
    /// Per candidate: its frame's representative for the sibling bounds
    /// (see [`SegmentStage::rank_with_bounds`]).
    representatives: Vec<u32>,
    /// Per gene or frame: the strongest candidate found so far.
    group_best: Vec<u32>,
    /// Alignments computed while ranking, consumed by traceback.
    traced: Vec<Option<Alignment>>,
    /// Per candidate: its entry in `traced`, or [`NONE`].
    cached: Vec<u32>,
    /// Frame paths of the traced representatives; `paths[..path_count]`
    /// are live, the rest keep their buffers.
    paths: Vec<FramePath>,
    path_count: usize,
    /// Per live path: its request and representative candidate index, in
    /// request order.
    path_owners: Vec<(usize, u32)>,
}

impl BatchScratch {
    fn reset(&mut self, requests: &[SegmentRequest<'_>]) {
        fn clear<T>(buffers: &mut Vec<Vec<T>>, len: usize) {
            if buffers.len() < len {
                buffers.resize_with(len, Vec::new);
            }
            buffers.iter_mut().for_each(Vec::clear);
        }
        clear(&mut self.ranked, requests.len());
        clear(&mut self.accepted, requests.len());
        clear(&mut self.deferred, requests.len());
        clear(&mut self.pending, requests.len());
        self.barriers.clear();
        self.barriers.resize(requests.len(), None);
        self.starts.clear();
        self.starts.push(0);
        let mut total = 0;
        for request in requests {
            total += request.candidates.len();
            self.starts.push(total);
        }
        self.leaders.clear();
        self.leaders.resize(total, 0);
        self.representatives.clear();
        self.representatives.resize(total, 0);
        self.exact.clear();
        self.exact.resize(total, 0);
        self.lower.clear();
        self.lower.resize(total, 0);
        self.early.clear();
        self.traced.clear();
        self.cached.clear();
        self.cached.resize(total, NONE);
        self.path_count = 0;
        self.path_owners.clear();
    }

    /// Per-candidate buffer slots of request `index`.
    #[inline]
    fn slots(&self, index: usize) -> Range<usize> {
        self.starts[index]..self.starts[index + 1]
    }
}

/// E-value thresholds of a stage's hits (Karlin–Altschul, against the whole
/// database for the query window searched).
///
/// * `accept`: the largest E-value of a hit kept at all — at most one chance
///   hit of that score is expected in the window.
/// * `anchor`: the largest E-value of a hit that can *anchor* a call. A read
///   is called only when a segment of its path keeps a clip that is this
///   significant ([`Placement::anchors`]), so a random read — whose best V
///   or J is never that significant — gets no call, while a J downstream of
///   a significant V needs only `accept`: given the V, a J is expected
///   there.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvalueGates {
    pub accept: f64,
    pub anchor: f64,
}

impl Default for EvalueGates {
    fn default() -> Self {
        Self {
            accept: 1.0,
            anchor: 1e-3,
        }
    }
}

/// Karlin–Altschul score thresholds for [`EvalueGates`]: a hit has E-value
/// at most `E` when its local score is at least `⌈(ln(K·m·n) − ln E) / λ⌉`
/// (`m` window length, `n` database length).
///
/// Without such a floor every read, random sequence included, has some
/// 10-column J alignment at 60 % identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Significance {
    pub statistics: KarlinAltschul,
    /// Total length of the database (all alleles; conservative for
    /// near-identical alleles).
    pub database_len: usize,
    pub gates: EvalueGates,
}

impl Significance {
    pub fn new(set: &GermlineSet, scoring: &AffineScoring, gates: EvalueGates) -> Option<Self> {
        Some(Self {
            statistics: KarlinAltschul::ungapped(scoring)?,
            database_len: set.total_len(),
            gates,
        })
    }

    /// Smallest local score with E-value ≤ `evalue` for a window of
    /// `query_len` ([`KarlinAltschul::min_score`]).
    pub fn min_score(&self, query_len: usize, evalue: f64) -> Score {
        self.statistics
            .min_score(query_len, self.database_len, evalue)
    }
}

pub struct SegmentStage<R, Acc = HitAcceptance> {
    pub segment: Arc<PreparedSegment>,
    pub retriever: R,
    pub keep: KeepPolicy,
    /// Which traced alignments take a slot, judged before the outer flank
    /// is stitched on: [`HitAcceptance`] for V, J and C, the admission
    /// gates for D. A rejected alignment frees its slot for the next one.
    pub accept: Acc,
    /// `None` accepts hits of any significance, and none anchors a call.
    pub significance: Option<Significance>,
    /// Stitch the co-linear V prefix or J suffix onto every accepted
    /// alignment; the kept alignments are ranked with it. The flank must
    /// lie in the searched window (V windows start at the read's start, J
    /// windows end at its end): a local optimum's flank then scores at most
    /// 0, so the local score bounds the extended one.
    pub outer: OuterExtend,
    /// Order of equal local scores in the ranking.
    pub ties: TieOrder,
    /// Score one allele per gene frame first, then skip the DP of
    /// same-frame siblings whose bounds settle them (see
    /// [`crate::family_bound`]); lossless.
    pub family_delta: Option<FamilyBoundTable>,
}

/// No candidate (a per-gene slot of [`BatchScratch`]'s leader search), or no
/// cached alignment.
const NONE: u32 = u32::MAX;

/// A request whose traceback waves have ended.
const DONE: usize = usize::MAX;

/// A stage's ranking key, best first: the score with the outer flank, then
/// the longest exact run.
pub(crate) fn rank(p: &Placement) -> std::cmp::Reverse<(Score, u32)> {
    std::cmp::Reverse((p.score(), p.evidence.longest_exact_run))
}

/// Whether a local score enters the ranking: positive (an alignment exists)
/// and at least the request's floor.
#[inline]
fn ranks(score: Score, floor: Score) -> bool {
    score > 0 && score >= floor
}

impl<R: Retriever> SegmentStage<R> {
    pub fn new(segment: Arc<PreparedSegment>, retriever: R, keep: KeepPolicy) -> Self {
        Self {
            segment,
            retriever,
            keep,
            accept: HitAcceptance::PERMISSIVE,
            significance: None,
            outer: OuterExtend::None,
            ties: TieOrder::Retrieval,
            family_delta: None,
        }
    }
}

impl<R, Acc> SegmentStage<R, Acc> {
    /// The stage with acceptance predicate `accept`.
    pub fn with_accept<Acc2>(self, accept: Acc2) -> SegmentStage<R, Acc2> {
        SegmentStage {
            segment: self.segment,
            retriever: self.retriever,
            keep: self.keep,
            accept,
            significance: self.significance,
            outer: self.outer,
            ties: self.ties,
            family_delta: self.family_delta,
        }
    }
}

impl<R: Retriever, Acc: Admission> SegmentStage<R, Acc> {
    pub fn kind(&self) -> SegmentKind {
        self.segment.kind()
    }

    /// Retrieve candidates for `query[window]` (bands relative to the window).
    pub fn retrieve(
        &self,
        query: &[u8],
        window: Range<usize>,
        scratch: &mut RetrievalScratch,
        out: &mut Vec<Candidate>,
    ) {
        self.retriever.retrieve(&query[window], scratch, out);
    }

    /// Lowest local score worth keeping for `request`: the keep floor, raised
    /// to the significance floor of its window.
    fn floor(&self, request: &SegmentRequest<'_>) -> Score {
        let significant = self.significance.map_or(Score::MIN, |s| {
            s.min_score(request.window.len(), s.gates.accept)
        });
        self.keep.min_score.max(significant)
    }

    /// The alignment problem of candidate `index` of `request`.
    fn job<'a>(&'a self, request: &SegmentRequest<'a>, index: u32) -> AlignJob<'a> {
        let candidate = &request.candidates[index as usize];
        AlignJob {
            query: &request.query[request.window.clone()],
            profile: self.segment.profile(candidate.allele),
            gaps: self.segment.gaps(),
            band: candidate.band,
        }
    }

    /// Align `candidates` in `query[window]` and append the kept placements to `out`.
    #[allow(clippy::too_many_arguments)]
    pub fn place_candidates<A: Aligner>(
        &self,
        aligner: &A,
        candidates: &[Candidate],
        query: &[u8],
        window: Range<usize>,
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
        out: &mut Vec<Placement>,
    ) {
        let request = SegmentRequest {
            read: 0,
            query,
            window,
            candidates,
        };
        self.place_many(
            aligner,
            std::slice::from_ref(&request),
            dp,
            batch,
            std::slice::from_mut(out),
        );
    }

    /// Place the candidates of many requests (typically many reads) with one
    /// score batch and one traceback batch per wave, so SIMD aligners fill
    /// their lanes. Per request the result is exactly that of placing it alone.
    pub fn place_many<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
        outs: &mut [Vec<Placement>],
    ) {
        batch.reset(requests);
        self.rank(aligner, requests, dp, batch);
        self.trace(aligner, requests, dp, batch);
        for (index, request) in requests.iter().enumerate() {
            let out = &mut outs[request.read];
            let first = out.len();
            out.append(&mut batch.pending[index]);
            self.finish(out, first, &mut batch.scores);
        }
    }

    /// With at most `top + 1` candidates, tracing all of them is cheaper than
    /// scoring all of them and then tracing `top`, unless `cap` would stop
    /// the unranked traceback before the last of them.
    fn needs_scores(&self, request: &SegmentRequest<'_>) -> bool {
        request.candidates.len() > self.keep.top.saturating_add(1).min(self.keep.cap)
    }

    /// Whether a request ranked without sibling bounds aligns every
    /// candidate instead of scoring it. A traceback costs about 1.5
    /// score-only DPs (1.65 on the scalar kernel), so once at least half the
    /// candidates reach traceback anyway, aligning all of them is cheaper
    /// than scoring all and aligning the top again; an alignment's score is
    /// the exact score ([`Aligner`]'s contract), so the ranking is the same.
    /// The D stage (`top = cap = 32`, about 40 candidates per window) is the
    /// case.
    fn aligns_all(&self, request: &SegmentRequest<'_>) -> bool {
        let traced = self.keep.top.min(self.keep.cap);
        request.candidates.len() <= traced.saturating_mul(2)
    }

    /// Each candidate's gene leader (see [`TieOrder::GeneLeadersFirst`]) as
    /// a candidate index.
    fn find_leaders(&self, candidates: &[Candidate], leaders: &mut [u32], best: &mut Vec<u32>) {
        debug_assert!(
            (1..candidates.len()).all(|i| candidates[..i]
                .iter()
                .all(|c| c.allele != candidates[i].allele)),
            "a request lists each allele once"
        );
        let set = self.segment.set();
        let gene = |allele| set.allele(allele).gene as usize;
        let strength = |c: &Candidate| (c.votes, std::cmp::Reverse(c.allele));
        strongest_per_group(candidates, set.gene_count(), gene, strength, leaders, best);
    }

    /// Fill `batch.ranked[r]` with request `r`'s candidates ranked by local
    /// score, ties in the stage's [`TieOrder`]. Requests that skip scoring
    /// keep retrieval order.
    fn rank<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        for (index, request) in requests.iter().enumerate() {
            if self.needs_scores(request) {
                let slots = batch.slots(index);
                self.find_leaders(
                    request.candidates,
                    &mut batch.leaders[slots],
                    &mut batch.group_best,
                );
            } else {
                let count = request.candidates.len() as u32;
                batch.ranked[index].extend((0..count).map(|c| (c, Score::MAX)));
            }
        }
        match &self.family_delta {
            None => {
                batch.owners.clear();
                for (index, request) in requests.iter().enumerate() {
                    if self.needs_scores(request) {
                        let count = request.candidates.len() as u32;
                        let jobs = if self.aligns_all(request) {
                            &mut batch.early
                        } else {
                            &mut batch.owners
                        };
                        jobs.extend((0..count).map(|c| (index, c)));
                    }
                }
                self.trace_early(aligner, requests, dp, batch);
                self.score_owners(aligner, requests, dp, batch);
            }
            Some(table) => self.rank_with_bounds(aligner, requests, table, dp, batch),
        }
        for (index, request) in requests.iter().enumerate() {
            if self.needs_scores(request) {
                self.sort_ranked(batch, index, 0);
            }
        }
    }

    /// Sort `batch.ranked[index][from..]` by local score, ties in the
    /// stage's [`TieOrder`].
    fn sort_ranked(&self, batch: &mut BatchScratch, index: usize, from: usize) {
        let leaders = &batch.leaders[batch.starts[index]..batch.starts[index + 1]];
        batch.ranked[index][from..].sort_unstable_by_key(|&(candidate, score)| {
            (std::cmp::Reverse(score), self.ties.key(leaders, candidate))
        });
    }

    /// Score-only DP for every job in `batch.owners`: record each exact
    /// score and append those that rank to their request's ranking.
    fn score_owners<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        if batch.owners.is_empty() {
            return;
        }
        let jobs: Vec<AlignJob<'a>> = batch
            .owners
            .iter()
            .map(|&(index, candidate)| self.job(&requests[index], candidate))
            .collect();
        aligner.score_batch(&jobs, dp, &mut batch.hits);
        for (&(index, candidate), hit) in batch.owners.iter().zip(&batch.hits) {
            let score = hit.map_or(0, |h| h.score);
            batch.exact[batch.starts[index] + candidate as usize] = score;
            if ranks(score, self.floor(&requests[index])) {
                batch.ranked[index].push((candidate, score));
            }
        }
    }

    /// Family-delta ranking (see [`crate::family_bound`]). Each frame
    /// present in a request has one representative `R` for its siblings:
    /// its strongest candidate by votes, then by the frame's medoid order
    /// ([`FamilyBoundTable::medoid_rank`]), then by lowest id. They serve
    /// the bounds only; ties keep the stage's [`TieOrder`], whatever was
    /// scored first.
    ///
    /// 1. Score the representatives. A score below the floor still bounds.
    ///    The strongest one is aligned instead: it is most often the best
    ///    and traced anyway.
    /// 2. Align each request's best-scoring representatives (every tie) too.
    ///    Their paths and the strongest's bound siblings from below, and
    ///    traceback reuses the alignments.
    /// 3. Bound each sibling `A`: `U = S*(R) + Δ(R, A)` through its
    ///    representative, `L` the best lower bound over the request's traced
    ///    paths. The cutoff `c` is the `top`-th largest of one lower bound
    ///    per candidate that ranks — `S*(R)` of a representative, `L` of a
    ///    sibling — so `c ≤ T`, the `top`-th largest exact score.
    /// 4. Decide each sibling, in this order:
    ///    * *pruned* when `U < floor`: it does not rank;
    ///    * *deferred* when `U < c`: it scores below `T`, so the ranking
    ///      from `c` up is that of full scoring;
    ///    * *known* when `U = L`: that is its exact score, without DP;
    ///    * *traced now* when `L ≥ c`: it likely reaches the traceback, so
    ///      it is aligned instead of scored, its alignment's score is its
    ///      exact score ([`Aligner`]'s contract) and traceback reuses it;
    ///    * *scored* otherwise.
    /// 5. The traceback waves stop at `c`, the request's *barrier*. Waves
    ///    that end above it end as under full scoring: a deferred sibling
    ///    scores below `c`, below every wave's cutoff. When they reach it,
    ///    every deferred sibling that may still be traced is scored and
    ///    ranked — all of them while fewer than `top` alignments are kept,
    ///    else those with `U` at least the `top`-th kept score (the outer
    ///    extension can put it below `c`) — and the waves continue past
    ///    the barrier over a ranking that is, from the cutoff up, that of
    ///    full scoring, positions included (so `cap` binds alike).
    ///
    /// Aligning representatives, the cutoff and early tracing apply only
    /// under permissive acceptance ([`Admission::admits_all`]).
    /// Otherwise a rejected alignment makes the traceback waves refill from
    /// below `T`, where a deferred sibling may lie, so siblings are pruned
    /// only below the floor.
    /// A sibling without an upper bound (another frame, or a band that is
    /// not its representative's frame band) is never pruned or known.
    fn rank_with_bounds<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        table: &FamilyBoundTable,
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        let cut = self.accept.admits_all() && self.keep.top > 0;
        batch.owners.clear();
        for (index, request) in requests.iter().enumerate() {
            if self.needs_scores(request) {
                let slots = batch.slots(index);
                let frame = |allele| table.frame(allele) as usize;
                let strength = |c: &Candidate| {
                    let medoid = table.medoid_rank(c.allele);
                    (c.votes, std::cmp::Reverse((medoid, c.allele)))
                };
                strongest_per_group(
                    request.candidates,
                    table.frame_count(),
                    frame,
                    strength,
                    &mut batch.representatives[slots.clone()],
                    &mut batch.group_best,
                );
                let representatives = (0..)
                    .zip(&batch.representatives[slots])
                    .filter(|&(c, &representative)| c == representative)
                    .map(|(c, _)| c);
                let strongest = cut
                    .then(|| {
                        representatives
                            .clone()
                            .max_by_key(|&c| strength(&request.candidates[c as usize]))
                    })
                    .flatten();
                for c in representatives {
                    if Some(c) == strongest {
                        batch.early.push((index, c));
                    } else {
                        batch.owners.push((index, c));
                    }
                }
            }
        }
        self.score_owners(aligner, requests, dp, batch);
        if cut {
            self.trace_early(aligner, requests, dp, batch);
            self.trace_best_representatives(aligner, requests, table, dp, batch);
            batch.early.clear();
        }
        batch.owners.clear();
        // `batch.path_owners[paths]`: the current request's traced paths.
        let mut paths = 0..0;
        for (index, request) in requests.iter().enumerate() {
            paths.start = paths.end;
            while paths.end < batch.path_count && batch.path_owners[paths.end].0 == index {
                paths.end += 1;
            }
            if !self.needs_scores(request) {
                continue;
            }
            let floor = self.floor(request);
            let start = batch.starts[index];
            let cutoff = if cut {
                // Representatives that rank are in `ranked`; add the siblings' `L`.
                batch.scores.clear();
                batch
                    .scores
                    .extend(batch.ranked[index].iter().map(|&(_, s)| s));
                for (c, sibling) in (0..).zip(request.candidates) {
                    let slot = start + c as usize;
                    if batch.representatives[slot] == c {
                        continue;
                    }
                    let lower = paths
                        .clone()
                        .filter_map(|p| {
                            let representative = batch.path_owners[p].1;
                            table.path_lower_bound(
                                &batch.paths[p],
                                &request.candidates[representative as usize],
                                sibling,
                                &self.segment,
                            )
                        })
                        .max()
                        .unwrap_or(0);
                    batch.lower[slot] = lower;
                    if ranks(lower, floor) {
                        batch.scores.push(lower);
                    }
                }
                nth_largest(&mut batch.scores, self.keep.top)
            } else {
                None
            };
            for (c, sibling) in (0..).zip(request.candidates) {
                let slot = start + c as usize;
                let representative = batch.representatives[slot];
                if representative == c {
                    continue;
                }
                // No positive cell means the representative's optimum is 0.
                let score = batch.exact[start + representative as usize];
                let representative = &request.candidates[representative as usize];
                let upper = table.upper_bound(representative, sibling, score);
                let lower = batch.lower[slot];
                let reaches = |bound: Score| cutoff.is_none_or(|cutoff| bound >= cutoff);
                match upper {
                    Some(upper) if upper < floor => {}
                    Some(upper) if !reaches(upper) => {
                        batch.deferred[index].push((c, upper));
                        batch.barriers[index] = cutoff;
                    }
                    Some(upper) if upper == lower => {
                        batch.exact[slot] = upper;
                        if ranks(upper, floor) {
                            batch.ranked[index].push((c, upper));
                        }
                    }
                    _ if cutoff.is_some() && reaches(lower) => batch.early.push((index, c)),
                    _ => batch.owners.push((index, c)),
                }
            }
        }
        self.trace_early(aligner, requests, dp, batch);
        self.score_owners(aligner, requests, dp, batch);
    }

    /// Align the candidates in `batch.early`: record each exact score, rank
    /// those that rank, and cache the alignments for traceback.
    fn trace_early<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        if batch.early.is_empty() {
            return;
        }
        let jobs: Vec<AlignJob<'a>> = batch
            .early
            .iter()
            .map(|&(index, candidate)| self.job(&requests[index], candidate))
            .collect();
        aligner.align_batch(&jobs, dp, &mut batch.alignments);
        for (&(index, candidate), alignment) in batch.early.iter().zip(batch.alignments.drain(..)) {
            let slot = batch.starts[index] + candidate as usize;
            let score = alignment.as_ref().map_or(0, |a| a.score);
            batch.exact[slot] = score;
            if ranks(score, self.floor(&requests[index])) {
                batch.ranked[index].push((candidate, score));
            }
            batch.cached[slot] = entry(batch.traced.len());
            batch.traced.push(alignment);
        }
    }

    /// Align each scored request's best-scoring representatives (every tie;
    /// `batch.ranked` holds the representatives that rank) not aligned yet,
    /// cache the alignments for traceback, and index their paths and the
    /// strongest representative's (`batch.early`) by frame position.
    fn trace_best_representatives<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        table: &FamilyBoundTable,
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        batch.owners.clear();
        for (index, request) in requests.iter().enumerate() {
            if !self.needs_scores(request) {
                continue;
            }
            let ranked = &batch.ranked[index];
            let start = batch.starts[index];
            if let Some(best) = ranked.iter().map(|&(_, s)| s).max() {
                batch.owners.extend(
                    ranked
                        .iter()
                        .filter(|&&(c, s)| s == best && batch.cached[start + c as usize] == NONE)
                        .map(|&(c, _)| (index, c)),
                );
            }
        }
        if !batch.owners.is_empty() {
            let jobs: Vec<AlignJob<'a>> = batch
                .owners
                .iter()
                .map(|&(index, candidate)| self.job(&requests[index], candidate))
                .collect();
            aligner.align_batch(&jobs, dp, &mut batch.alignments);
            for (&(index, candidate), alignment) in
                batch.owners.iter().zip(batch.alignments.drain(..))
            {
                batch.cached[batch.starts[index] + candidate as usize] = entry(batch.traced.len());
                batch.traced.push(alignment);
            }
        }
        let BatchScratch {
            ranked,
            early,
            starts,
            exact,
            traced,
            cached,
            paths,
            path_count,
            path_owners,
            ..
        } = batch;
        // Paths of the best representatives and of the strongest, if another.
        let mut strongest = early.iter().peekable();
        for (index, request) in requests.iter().enumerate() {
            if !self.needs_scores(request) {
                continue;
            }
            let start = starts[index];
            let best = ranked[index].iter().map(|&(_, s)| s).max();
            let other = strongest
                .next_if(|&&(i, _)| i == index)
                .map(|&(_, c)| c)
                .filter(|&c| best != Some(exact[start + c as usize]));
            let chosen = ranked[index]
                .iter()
                .filter(|&&(_, s)| Some(s) == best)
                .map(|&(c, _)| c)
                .chain(other);
            for candidate in chosen {
                let Some(alignment) = &traced[cached[start + candidate as usize] as usize] else {
                    continue;
                };
                let allele = request.candidates[candidate as usize].allele;
                if paths.len() == *path_count {
                    paths.push(FramePath::default());
                }
                table.frame_path(
                    allele,
                    alignment,
                    &request.query[request.window.clone()],
                    self.segment.profile(allele),
                    self.segment.gaps(),
                    &mut paths[*path_count],
                );
                *path_count += 1;
                path_owners.push((index, candidate));
            }
        }
    }

    /// The next traceback wave of request `index` from `cursor`, stopping
    /// at its barrier: `(end, barrier reached)`.
    fn wave(&self, batch: &mut BatchScratch, index: usize, cursor: usize) -> (usize, bool) {
        let ranked = &batch.ranked[index];
        batch.scores.clear();
        batch.scores.extend(ranked.iter().map(|&(_, s)| s));
        let shared = batch.barriers[index].map_or(ranked.len(), |barrier| {
            batch.scores.partition_point(|&s| s >= barrier)
        });
        // `accepted` is kept in descending order.
        let accepted = &batch.accepted[index];
        let end = self.keep.next_wave(
            &batch.scores[..shared],
            cursor,
            accepted.len(),
            self.cutoff(accepted),
        );
        (end, end == cursor && batch.barriers[index].is_some())
    }

    /// The `top`-th best kept score of a request, once `top` are kept.
    fn cutoff(&self, accepted: &[Score]) -> Option<Score> {
        self.keep
            .top
            .checked_sub(1)
            .and_then(|k| accepted.get(k))
            .copied()
    }

    /// Lift the barriers the traceback waves have reached: score and rank
    /// the deferred siblings that may still be traced (see
    /// [`Self::rank_with_bounds`]).
    fn lift_barriers<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        batch.owners.clear();
        for index in 0..requests.len() {
            let cursor = batch.cursors[index];
            if cursor == DONE
                || batch.barriers[index].is_none()
                || !self.wave(batch, index, cursor).1
            {
                continue;
            }
            batch.barriers[index] = None;
            // At or past `cap` no further candidate is traced.
            if cursor < self.keep.cap {
                let reach = self.cutoff(&batch.accepted[index]).unwrap_or(Score::MIN);
                batch.owners.extend(
                    batch.deferred[index]
                        .iter()
                        .filter(|&&(_, upper)| upper >= reach)
                        .map(|&(c, _)| (index, c)),
                );
            }
        }
        if batch.owners.is_empty() {
            return;
        }
        self.score_owners(aligner, requests, dp, batch);
        // They score below the barrier, below every candidate traced so far.
        let mut owners = std::mem::take(&mut batch.owners);
        owners.dedup_by_key(|&mut (index, _)| index);
        for &(index, _) in &owners {
            self.sort_ranked(batch, index, batch.cursors[index]);
        }
        batch.owners = owners;
    }

    /// Traceback waves until every request has `top` accepted alignments
    /// (plus ties) and the next local score is below the `top`-th extended
    /// one, runs out of candidates, or reaches `cap`. Alignments cached
    /// while ranking are taken, not recomputed.
    fn trace<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        batch.cursors.clear();
        batch.cursors.resize(requests.len(), 0);
        loop {
            self.lift_barriers(aligner, requests, dp, batch);
            batch.owners.clear();
            for index in 0..requests.len() {
                let cursor = batch.cursors[index];
                if cursor == DONE {
                    continue;
                }
                let (end, _) = self.wave(batch, index, cursor);
                batch.owners.extend(
                    batch.ranked[index][cursor..end]
                        .iter()
                        .map(|&(c, _)| (index, c)),
                );
                // An empty wave changes nothing the next one depends on.
                batch.cursors[index] = if end == cursor { DONE } else { end };
            }
            if batch.owners.is_empty() {
                break;
            }
            let jobs: Vec<AlignJob<'a>> = batch
                .owners
                .iter()
                .filter(|&&(index, c)| batch.cached[batch.starts[index] + c as usize] == NONE)
                .map(|&(index, candidate)| self.job(&requests[index], candidate))
                .collect();
            aligner.align_batch(&jobs, dp, &mut batch.alignments);
            let BatchScratch {
                owners,
                alignments,
                starts,
                traced,
                cached,
                accepted,
                pending,
                ..
            } = &mut *batch;
            let mut fresh = alignments.drain(..);
            for &(index, candidate) in owners.iter() {
                let alignment = match cached[starts[index] + candidate as usize] {
                    NONE => fresh.next().expect("one alignment per job"),
                    cached => traced[cached as usize].take(),
                };
                let Some(alignment) = alignment else {
                    continue;
                };
                self.keep_alignment(
                    &requests[index],
                    candidate,
                    alignment,
                    &mut accepted[index],
                    &mut pending[index],
                );
            }
        }
    }

    /// Place a traced alignment of candidate `candidate` of `request` and,
    /// when it clears the floor and acceptance, extend its outer flank and
    /// keep it under its extended score.
    fn keep_alignment(
        &self,
        request: &SegmentRequest<'_>,
        candidate: u32,
        mut alignment: Alignment,
        accepted: &mut Vec<Score>,
        pending: &mut Vec<Placement>,
    ) {
        if alignment.score < self.floor(request) {
            return;
        }
        alignment.shift_query(request.window.start as u32);
        let local = alignment.score;
        let mut placement = Placement::new(
            self.kind(),
            request.candidates[candidate as usize].allele,
            alignment,
            request.query,
            &self.segment,
        );
        placement.anchor = self
            .significance
            .map(|s| s.min_score(request.window.len(), s.gates.anchor))
            .filter(|&floor| local >= floor);
        let reference = self.segment.reference(placement.primary_allele());
        if self.accept.admit(&placement, request.query, reference) {
            // The flank stays in the searched window: a J window ends at
            // the read's end unless the pair guards cap it.
            self.extend_outer(&mut placement, &request.query[..request.window.end]);
            let key = placement.score();
            debug_assert!(key <= local, "the outer flank lies in the searched window");
            let at = accepted.partition_point(|&s| s >= key);
            accepted.insert(at, key);
            pending.push(placement);
        }
    }

    /// Keep the top accepted alignments and their ties by extended score
    /// and those by local score, rank them by extended score, and merge
    /// identical alignments. Equal keys keep the traceback order: local
    /// score, then the stage's [`TieOrder`]. `scores` is scratch.
    ///
    /// The local top guards the extended ranking's blind spot: a hit with no
    /// room for a flank (at the read's or the germline's end) pays none,
    /// and a few alleles sharing one such hit could take every slot from a
    /// hit that is stronger locally but pays for its mutated flank. The
    /// waves have traced the local top: the `top`-th local score bounds the
    /// `top`-th extended one from above, and they stop below the latter.
    fn finish(&self, out: &mut Vec<Placement>, first: usize, scores: &mut Vec<Score>) {
        out[first..].sort_by_key(rank);
        if self.keep.top > 0 && out.len() - first > self.keep.top {
            let extended = out[first + self.keep.top - 1].score();
            scores.clear();
            scores.extend(out[first..].iter().map(|p| p.local_score));
            let local = nth_largest(scores, self.keep.top).expect("more than `top` placements");
            let mut kept = first;
            for i in first..out.len() {
                if out[i].score() >= extended || out[i].local_score >= local {
                    out.swap(kept, i);
                    kept += 1;
                }
            }
            out.truncate(kept);
        }
        merge_identical(out, first);
    }

    fn extend_outer(&self, placement: &mut Placement, query: &[u8]) {
        if self.outer == OuterExtend::None {
            return;
        }
        let allele = placement.primary_allele();
        let (reference, profile) = (self.segment.reference(allele), self.segment.profile(allele));
        if extend_outer(
            &mut placement.alignment,
            query,
            reference,
            profile,
            self.outer,
        ) {
            placement.clip =
                ClipIndex::build(&placement.alignment, query, profile, self.segment.gaps());
            placement.evidence =
                Evidence::from_columns(placement.alignment.columns(), query, reference);
        }
    }
}

/// The `n`-th largest of `values` (reordered), if there are that many.
fn nth_largest(values: &mut [Score], n: usize) -> Option<Score> {
    let k = n.checked_sub(1).filter(|&k| k < values.len())?;
    let (_, &mut nth, _) = values.select_nth_unstable_by_key(k, |&s| std::cmp::Reverse(s));
    Some(nth)
}

/// For each candidate, the strongest candidate of its group (the largest
/// `strength`; unique per candidate) as a candidate index. `best` is
/// per-group scratch, all [`NONE`] on entry and on return.
fn strongest_per_group<K: Ord>(
    candidates: &[Candidate],
    groups: usize,
    group: impl Fn(AlleleId) -> usize,
    strength: impl Fn(&Candidate) -> K,
    out: &mut [u32],
    best: &mut Vec<u32>,
) {
    if best.len() < groups {
        best.resize(groups, NONE);
    }
    for (index, candidate) in (0..).zip(candidates) {
        let slot = &mut best[group(candidate.allele)];
        if *slot == NONE || strength(candidate) > strength(&candidates[*slot as usize]) {
            *slot = index;
        }
    }
    for (strongest, candidate) in out.iter_mut().zip(candidates) {
        *strongest = best[group(candidate.allele)];
    }
    for candidate in candidates {
        best[group(candidate.allele)] = NONE;
    }
}

/// An index into [`BatchScratch::traced`].
fn entry(index: usize) -> u32 {
    u32::try_from(index).expect("traced alignments fit u32")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::BandedGotoh;
    use crate::alphabet::{decode, encode};
    use crate::band::Band;
    use crate::family_bound::FamilyBoundTable;
    use crate::retrieve::{GeneFrames, Retriever};
    use crate::scoring::AffineScoring;

    struct EmptyRetriever;
    impl Retriever for EmptyRetriever {
        fn retrieve(&self, _: &[u8], _: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
            out.clear();
        }
    }

    fn v_stage(family_delta: bool) -> (SegmentStage<EmptyRetriever>, Vec<u8>, Vec<Candidate>) {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let mut set = GermlineSet::new(SegmentKind::V);
        // Gene 1: two alleles, one SNP. Gene 2: distant decoy.
        let a1 = b"CAGGTCACCTTGAAGGAGTCTGGTCCTGCCTGGTGAAACCCACACAGACC";
        let mut a2 = a1.to_vec();
        a2[10] = b'G'; // SNP
        let decoy = b"GACATGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCC";
        set.push("IGHV1-1*01", a1, None).unwrap();
        set.push("IGHV1-1*02", &a2, None).unwrap();
        set.push("IGHV2-1*01", decoy, None).unwrap();
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        let keep = KeepPolicy {
            top: 2,
            min_score: 1,
            cap: 8,
        };
        let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
        stage.accept = HitAcceptance {
            min_subst: 10,
            min_identity: 0.50,
        };
        if family_delta {
            stage.family_delta = Some(FamilyBoundTable::build(&segment));
        }
        let query = encode(a1);
        let candidates = [(0, 40), (1, 30), (2, 5)]
            .map(|(allele, votes)| Candidate {
                allele,
                band: Band::Full,
                votes,
            })
            .to_vec();
        (stage, query, candidates)
    }

    #[test]
    fn family_delta_keeps_the_same_alleles_as_full_scoring() {
        let aligner = BandedGotoh;
        let mut placed = Vec::new();
        for family_delta in [false, true] {
            let (stage, query, candidates) = v_stage(family_delta);
            let (mut dp, mut batch, mut out) =
                (Default::default(), BatchScratch::default(), Vec::new());
            stage.place_candidates(
                &aligner,
                &candidates,
                &query,
                0..query.len(),
                &mut dp,
                &mut batch,
                &mut out,
            );
            let mut ids: Vec<_> = out.iter().flat_map(|p| p.alleles.clone()).collect();
            ids.sort_unstable();
            placed.push(ids);
        }
        assert_eq!(placed[0], placed[1]);
        assert!(placed[0].contains(&0));
    }

    /// `(alleles, score)` of the placements `stage` keeps for `candidates`.
    fn placed<R: Retriever, Acc: Admission>(
        stage: &SegmentStage<R, Acc>,
        candidates: &[Candidate],
        query: &[u8],
    ) -> Vec<(Vec<AlleleId>, Score)> {
        let (mut dp, mut batch, mut out) =
            (Default::default(), BatchScratch::default(), Vec::new());
        stage.place_candidates(
            &BandedGotoh,
            candidates,
            query,
            0..query.len(),
            &mut dp,
            &mut batch,
            &mut out,
        );
        out.iter().map(|p| (p.alleles.clone(), p.score())).collect()
    }

    /// The kept placements are the top ones by the *extended* score, the
    /// stage's objective: a 5′-partial allele whose full-length sibling
    /// scores higher locally but loses more to its stitched-on flank is the
    /// one kept, ranked or not, bounded or not.
    #[test]
    fn top_cut_uses_the_extended_score() {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let mut next = crate::testing::rng(11);
        let full = crate::testing::random_codes(&mut next, 60);
        let decoy = crate::testing::random_codes(&mut next, 60);
        let partial = full[12..].to_vec();
        let mut set = GermlineSet::new(SegmentKind::V);
        for (name, codes) in [
            ("IGHV1-1*01", &full),
            ("IGHV1-1*02", &partial),
            ("IGHV2-1*01", &decoy),
        ] {
            set.push(name, decode(codes).as_bytes(), None).unwrap();
        }
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        // Twelve flank bases that mismatch the full allele, but for two
        // chance matches just before the partial allele's start.
        let mut query: Vec<u8> = full[..12].iter().map(|&b| (b + 1) % 4).collect();
        query[10] = full[10];
        query[11] = full[11];
        query.extend_from_slice(&partial);
        let keep = KeepPolicy {
            top: 1,
            min_score: 1,
            cap: 8,
        };
        for family_delta in [false, true] {
            let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
            stage.outer = OuterExtend::Prefix;
            stage.family_delta = family_delta.then(|| FamilyBoundTable::build(&segment));
            let candidates = [0, 1, 2].map(|allele| Candidate {
                allele,
                band: Band::Full,
                votes: 0,
            });
            // Local: full 100 > partial 96. Extended: full 70 < partial 96.
            assert_eq!(placed(&stage, &candidates, &query)[0], (vec![1], 96));
            // Two candidates are traced unranked.
            assert_eq!(placed(&stage, &candidates[..2], &query)[0], (vec![1], 96));
        }
    }

    /// The top local alignments are kept too. Three J alleles share a hit at
    /// the read's end, which has no room for a flank, so its extended score
    /// is its local one; the real J scores higher locally but pays its
    /// mutated 3′ flank. By extended score alone the copies would fill every
    /// slot and drop it; the local ranking keeps it.
    #[test]
    fn top_cut_keeps_the_best_local_alignments() {
        let scoring = AffineScoring::new(2, -3, -17, -2);
        let mut next = crate::testing::rng(13);
        let real = crate::testing::random_codes(&mut next, 40);
        let shared = crate::testing::random_codes(&mut next, 12);
        let mut set = GermlineSet::new(SegmentKind::J);
        set.push("IGHJ1*01", decode(&real).as_bytes(), None)
            .unwrap();
        for name in ["IGHJ2*01", "IGHJ2*02", "IGHJ2*03"] {
            // The copies differ only past the read's end.
            let tail = crate::testing::random_codes(&mut next, 20);
            let copy = [shared.clone(), tail].concat();
            set.push(name, decode(&copy).as_bytes(), None).unwrap();
        }
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        // The real J's first 25 bases exact (local 50), then 15 mutated
        // flank bases but at offsets 5 and 10: stitched over 13 bases at
        // −29, extended 21. Then the shared hit, 12 bases at 24.
        let mut query = real[..25].to_vec();
        query.extend((25..40).map(|k| {
            let base = real[k];
            if [30, 35].contains(&k) {
                base
            } else {
                (base + 1) % 4
            }
        }));
        query.extend(&shared);
        let keep = KeepPolicy {
            top: 3,
            min_score: 1,
            cap: 64,
        };
        for family_delta in [false, true] {
            let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
            stage.outer = OuterExtend::Suffix;
            stage.family_delta = family_delta.then(|| FamilyBoundTable::build(&segment));
            let candidates = [0, 1, 2, 3].map(|allele| Candidate {
                allele,
                band: Band::Full,
                votes: 0,
            });
            let kept = placed(&stage, &candidates, &query);
            assert_eq!(kept, [(vec![1, 2, 3], 24), (vec![0], 21)]);
        }
    }

    /// The D stage's acceptance predicate is its admission gate: a hit
    /// without an exact run of `min_d_match` bases takes no slot, so the
    /// admissible hit ranked behind it is kept. The default D policy traces
    /// beyond its `top` for that.
    #[test]
    fn inadmissible_hits_free_their_slots() {
        use crate::admission::{AnyOf, DistributedEvidence, ExactRunFloor};
        let scoring = AffineScoring::new(2, -3, -13, -1);
        let mut next = crate::testing::rng(41);
        let (long, short) = (
            crate::testing::random_codes(&mut next, 32),
            crate::testing::random_codes(&mut next, 10),
        );
        let mut set = GermlineSet::new(SegmentKind::D);
        set.push("IGHD1-1*01", decode(&long).as_bytes(), None)
            .unwrap();
        set.push("IGHD2-2*01", decode(&short).as_bytes(), None)
            .unwrap();
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        // The long allele with every 4th base mutated (runs of 3, score
        // about 3 per 4 columns), then the short allele exact (run 10,
        // score 20).
        let mut query = crate::testing::random_codes(&mut next, 6);
        query.extend(
            (0..)
                .zip(&long)
                .map(|(i, &c)| if i % 4 == 3 { (c + 1) % 4 } else { c }),
        );
        query.extend(crate::testing::random_codes(&mut next, 6));
        query.extend(&short);
        query.extend(crate::testing::random_codes(&mut next, 6));
        let keep = KeepPolicy {
            top: 1,
            min_score: 1,
            cap: 2,
        };
        let gate = AnyOf(ExactRunFloor { min_run: 4 }, None::<DistributedEvidence>);
        let stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep).with_accept(gate);
        let candidates = [0, 1].map(|allele| Candidate {
            allele,
            band: Band::Full,
            votes: 0,
        });
        let traced = |stage: &SegmentStage<EmptyRetriever, _>| {
            let score = |allele| {
                BandedGotoh
                    .score(
                        &query,
                        segment.profile(allele),
                        scoring.gaps,
                        Band::Full,
                        &mut Default::default(),
                    )
                    .map_or(0, |hit| hit.score)
            };
            assert!(score(0) > score(1), "the inadmissible hit ranks first");
            placed(stage, &candidates, &query)
        };
        assert_eq!(traced(&stage), vec![(vec![1], 20)]);
        let d = crate::config::DSegmentConfig::default();
        assert!(d.keep.cap > d.keep.top);
    }

    /// A request of `top + 1` candidates is ranked by score when `cap` binds
    /// below it: traced unranked, in retrieval order, the last candidate
    /// would never be aligned (the D stage keeps `top = cap = 32`).
    #[test]
    fn a_binding_cap_still_ranks_top_plus_one_candidates() {
        let scoring = AffineScoring::new(2, -3, -13, -1);
        let mut set = GermlineSet::new(SegmentKind::D);
        set.push("IGHD1-1*01", b"GGTACAACTGGAACGAC", None).unwrap();
        set.push("IGHD2-2*01", b"AGGATATTGTAGTAGTACCAGCTGCTATGCC", None)
            .unwrap();
        set.push("IGHD3-3*01", b"GTATTACGATTTTTGGAGTGGTTATTATACC", None)
            .unwrap();
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        let keep = KeepPolicy {
            top: 2,
            min_score: 1,
            cap: 2,
        };
        let stage = SegmentStage::new(segment, EmptyRetriever, keep);
        let query = encode(b"CCCCGTATTACGATTTTTGGAGTGGTTATTATACCCCCC");
        let candidates: Vec<Candidate> = (0..3)
            .map(|allele| Candidate {
                allele,
                band: Band::Full,
                votes: 0,
            })
            .collect();
        let kept = placed(&stage, &candidates, &query);
        assert_eq!(
            kept.first().map(|(alleles, _)| alleles.as_slice()),
            Some(&[2][..]),
            "{kept:?}"
        );
    }

    /// Three alleles in two genes align identically to the read; the tie
    /// order alone decides the merged placement's allele order. Bounded and
    /// full scoring rank in a different sequence (gene leaders first), so
    /// they agree only because the order is an explicit key.
    #[test]
    fn tie_order_is_independent_of_bounding() {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let mut next = crate::testing::rng(8);
        let shared = crate::testing::random_codes(&mut next, 60);
        let other = crate::testing::random_codes(&mut next, 60);
        // SNPs past the read's end leave the three alignments identical.
        let snp = |at: usize| {
            let mut codes = shared.clone();
            codes[at] = (codes[at] + 1) % 4;
            codes
        };
        let mut set = GermlineSet::new(SegmentKind::V);
        for (name, codes) in [
            ("IGHV1-1*01", shared.clone()),
            ("IGHV1-1*02", snp(55)),
            ("IGHV2-1*01", other),
            ("IGHV2-1*02", snp(58)),
        ] {
            set.push(name, decode(&codes).as_bytes(), None).unwrap();
        }
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        let query = shared[..50].to_vec();
        // Leaders: IGHV1-1*01 (index 2) and IGHV2-1*01 (index 1).
        let candidates = [(1, 10), (2, 30), (0, 20), (3, 5)].map(|(allele, votes)| Candidate {
            allele,
            band: Band::Full,
            votes,
        });
        let keep = KeepPolicy {
            top: 1,
            min_score: 1,
            cap: 8,
        };
        for (ties, primary) in [(TieOrder::Retrieval, 1), (TieOrder::GeneLeadersFirst, 0)] {
            let placements: Vec<_> = [false, true]
                .map(|family_delta| {
                    let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
                    stage.ties = ties;
                    stage.family_delta = family_delta.then(|| FamilyBoundTable::build(&segment));
                    placed(&stage, &candidates, &query)
                })
                .into();
            assert_eq!(placements[0], placements[1], "{ties:?}");
            assert_eq!(placements[0][0].0[0], primary, "{ties:?}");
            assert_eq!(placements[0][0].0.len(), 3, "{ties:?}");
        }
    }

    /// The best representative outscores everything but fails acceptance,
    /// and the allele that replaces it is a sibling whose bound lies below
    /// that representative's score. Cutoff pruning would drop it, and the
    /// traceback waves would settle on its weaker representative instead.
    #[test]
    fn non_permissive_acceptance_disables_cutoff() {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let mut next = crate::testing::rng(5);
        let short = crate::testing::random_codes(&mut next, 20);
        let sibling = crate::testing::random_codes(&mut next, 48);
        let mut representative = sibling.clone();
        representative[25] = (representative[25] + 1) % 4;
        let mut set = GermlineSet::new(SegmentKind::V);
        for (name, codes) in [
            ("IGHV1-1*01", &short),
            ("IGHV2-1*01", &representative),
            ("IGHV2-1*02", &sibling),
        ] {
            set.push(name, decode(codes).as_bytes(), None).unwrap();
        }
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        // An exact copy of the short allele (20 columns), then the sibling
        // with every 4th base mutated (> 40 columns, identity ~ 3/4).
        let mut query = short.clone();
        query.extend(crate::testing::random_codes(&mut next, 10));
        query.extend(
            sibling
                .iter()
                .enumerate()
                .map(|(i, &c)| if i % 4 == 2 { (c + 2) % 4 } else { c }),
        );
        let candidates = [(0, 50), (1, 40), (2, 30)].map(|(allele, votes)| Candidate {
            allele,
            band: Band::Full,
            votes,
        });
        let keep = KeepPolicy {
            top: 1,
            min_score: 1,
            cap: 8,
        };
        let stage = |family_delta: bool| {
            let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
            stage.accept = HitAcceptance {
                min_subst: 30,
                min_identity: 0.0,
            };
            stage.family_delta = family_delta.then(|| FamilyBoundTable::build(&segment));
            stage
        };
        // Preconditions: the short allele is the best representative; the
        // sibling's bound is below its score; it alone fails acceptance.
        let score = |allele: AlleleId| {
            BandedGotoh
                .score(
                    &query,
                    segment.profile(allele),
                    scoring.gaps,
                    Band::Full,
                    &mut Default::default(),
                )
                .map_or(0, |hit| hit.score)
        };
        let table = FamilyBoundTable::build(&segment);
        let bound = table
            .upper_bound(&candidates[1], &candidates[2], score(1))
            .unwrap();
        assert!(score(0) > score(2) && score(2) > score(1) && bound < score(0));
        let full = placed(&stage(false), &candidates, &query);
        assert_eq!(full, vec![(vec![2], score(2))]);
        assert_eq!(placed(&stage(true), &candidates, &query), full);
    }

    /// The `(query, profile)` of a DP and whether it traced back.
    type Run = (*const u8, *const SubstitutionProfile, bool);

    /// [`BandedGotoh`] that records every DP it runs.
    #[derive(Default)]
    struct Recording(std::cell::RefCell<Vec<Run>>);

    impl Aligner for Recording {
        type Workspace = <BandedGotoh as Aligner>::Workspace;

        fn score(
            &self,
            query: &[u8],
            profile: &SubstitutionProfile,
            gaps: GapCosts,
            band: Band,
            workspace: &mut Self::Workspace,
        ) -> Option<ScoreHit> {
            self.0.borrow_mut().push((query.as_ptr(), profile, false));
            BandedGotoh.score(query, profile, gaps, band, workspace)
        }

        fn align(
            &self,
            query: &[u8],
            profile: &SubstitutionProfile,
            gaps: GapCosts,
            band: Band,
            workspace: &mut Self::Workspace,
        ) -> Option<Alignment> {
            self.0.borrow_mut().push((query.as_ptr(), profile, true));
            BandedGotoh.align(query, profile, gaps, band, workspace)
        }
    }

    /// A request ranked without sibling bounds whose candidates mostly reach
    /// traceback anyway (`n ≤ 2·min(top, cap)`, the D stage's case) aligns
    /// each candidate once and scores none, and keeps what ranking by score
    /// keeps: the first `cap` candidates by score, then retrieval order.
    #[test]
    fn few_candidates_are_aligned_not_scored() {
        let scoring = AffineScoring::new(2, -3, -13, -1);
        let mut next = crate::testing::rng(31);
        let mut set = GermlineSet::new(SegmentKind::D);
        let alleles: Vec<Vec<u8>> = (0..6)
            .map(|_| crate::testing::random_codes(&mut next, 20))
            .collect();
        for (k, codes) in alleles.iter().enumerate() {
            let name = format!("IGHD{}-1*01", k + 1);
            set.push(&name, decode(codes).as_bytes(), None).unwrap();
        }
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        // Allele k contributes its last 20 − 2k bases: six different scores.
        let mut query = Vec::new();
        for (k, codes) in alleles.iter().enumerate() {
            query.extend(crate::testing::random_codes(&mut next, 4));
            query.extend(&codes[2 * k..]);
        }
        let candidates: Vec<Candidate> = [3, 0, 5, 1, 4, 2]
            .map(|allele| Candidate {
                allele,
                band: Band::Full,
                votes: 0,
            })
            .into();
        for (bound, scored) in [(4, 0), (2, 6)] {
            let keep = KeepPolicy {
                top: bound,
                min_score: 1,
                cap: bound,
            };
            let stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
            let aligner = Recording::default();
            let (mut dp, mut batch, mut out) =
                (Default::default(), BatchScratch::default(), Vec::new());
            let window = 0..query.len();
            stage.place_candidates(
                &aligner,
                &candidates,
                &query,
                window,
                &mut dp,
                &mut batch,
                &mut out,
            );
            let runs = aligner.0.into_inner();
            assert_eq!(
                runs.iter().filter(|r| !r.2).count(),
                scored,
                "top = cap = {bound}"
            );
            if scored == 0 {
                assert_eq!(runs.len(), candidates.len());
            }
            // Ranking by exact score, ties in retrieval order; `cap` traced.
            let score = |c: &Candidate| {
                let profile = segment.profile(c.allele);
                let hit = BandedGotoh.score(
                    &query,
                    profile,
                    scoring.gaps,
                    c.band,
                    &mut Default::default(),
                );
                hit.map_or(0, |hit| hit.score)
            };
            let mut ranked: Vec<(Score, usize)> = (0..candidates.len())
                .map(|i| (score(&candidates[i]), i))
                .filter(|&(s, _)| s >= keep.min_score)
                .collect();
            ranked.sort_by_key(|&(s, i)| (std::cmp::Reverse(s), i));
            let mut expected: Vec<(AlleleId, Score)> = ranked[..keep.cap.min(ranked.len())]
                .iter()
                .map(|&(s, i)| (candidates[i].allele, s))
                .collect();
            expected.sort_unstable();
            let mut kept: Vec<(AlleleId, Score)> = out
                .iter()
                .flat_map(|p| p.alleles.iter().map(|&a| (a, p.score())))
                .collect();
            kept.sort_unstable();
            assert_eq!(kept, expected, "top = cap = {bound}");
        }
    }

    /// Two V genes in gene frames: SNP alleles, 5′ and 3′ truncations, a 5′
    /// extension, an exact duplicate, which ties with its twin, and an
    /// internal deletion with a SNP variant (a second frame of the gene).
    fn frame_family() -> (Arc<PreparedSegment>, GeneFrames) {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let mut next = crate::testing::rng(17);
        let b = crate::testing::random_codes(&mut next, 90);
        let c = crate::testing::random_codes(&mut next, 90);
        let snp = |codes: &[u8], at: &[usize]| {
            let mut codes = codes.to_vec();
            for &at in at {
                codes[at] = (codes[at] + 1) % 4;
            }
            codes
        };
        let mut deleted = b.clone();
        deleted.drain(45..48);
        let extended = [crate::testing::random_codes(&mut next, 4), b.clone()].concat();
        let mut set = GermlineSet::new(SegmentKind::V);
        for (name, codes) in [
            ("IGHV1-1*01", b.clone()),
            ("IGHV1-1*02", snp(&b, &[20])),
            ("IGHV1-1*03", snp(&b, &[20, 61])),
            ("IGHV1-1*04", b[5..].to_vec()),
            ("IGHV1-1*05", snp(&b[..82], &[40])),
            ("IGHV1-1*06", extended),
            ("IGHV1-1*07", deleted.clone()),
            ("IGHV1-1*08", snp(&b, &[20])),
            ("IGHV2-1*01", c.clone()),
            ("IGHV2-1*02", snp(&c, &[30])),
            ("IGHV2-1*03", c[3..].to_vec()),
            ("IGHV1-1*09", snp(&deleted, &[70])),
        ] {
            set.push(name, decode(&codes).as_bytes(), None).unwrap();
        }
        let frames = GeneFrames::build(&set, &crate::seed::SeedShape::contiguous(7));
        let frame = |allele| frames.frame(allele);
        assert!((1..6).chain([7]).all(|allele| frame(allele) == frame(0)));
        assert!(frame(6) != frame(0) && frame(11) == frame(6));
        let segment = Arc::new(PreparedSegment::new(set, &scoring, scoring.gaps));
        (segment, frames)
    }

    /// Bounded ranking places exactly what scoring every candidate does —
    /// alleles, their order and alignments — for batches of requests with
    /// more candidates than `top + 1` (and some with fewer), shifted frame
    /// bands with either edge moved, unbanded requests, siblings with an
    /// inconsistent band, every `top` and `cap`, both acceptance kinds and
    /// every outer extension over windows that hold its flank; and it skips
    /// the DP of some siblings entirely and aligns others without scoring
    /// them (known exactly, or traced early). Without a binding `cap`, both
    /// keep the `top` accepted alignments and their ties by extended score
    /// and those by local score, as aligning every candidate does.
    #[test]
    fn bounded_ranking_equals_brute_force() {
        let (segment, frames) = frame_family();
        let table = FamilyBoundTable::with_frames(&segment, &frames);
        let alleles = segment.set().len() as AlleleId;
        let mut next = crate::testing::rng(23);
        let (mut pruned, mut unscored) = (0, 0);
        for case in 0..720 {
            let keep = KeepPolicy {
                top: [1, 2, 3, 5][case % 4],
                min_score: 1,
                cap: [2, 3, 64][case / 4 % 3],
            };
            let outer =
                [OuterExtend::None, OuterExtend::Prefix, OuterExtend::Suffix][case / 12 % 3];
            let accept = if case % 7 == 6 {
                HitAcceptance {
                    min_subst: 40,
                    min_identity: 0.9,
                }
            } else {
                HitAcceptance::PERMISSIVE
            };
            let mut queries = Vec::new();
            let mut candidates = Vec::new();
            for _ in 0..3 {
                // A mutated stretch of one allele between random flanks.
                let source = (next() % u64::from(alleles)) as AlleleId;
                let codes = &segment.set().allele(source).codes;
                let from = (next() % 20) as usize;
                let to = codes.len() - (next() % 20) as usize;
                let lead = (next() % 12) as usize;
                let mut query = crate::testing::random_codes(&mut next, lead);
                for &code in &codes[from..to] {
                    match next() % 40 {
                        0..=2 => query.push((code + 1 + (next() % 3) as u8) % 4),
                        3 => query.push(4),
                        4 => {}
                        5 => query.extend([code, (next() % 4) as u8]),
                        _ => query.push(code),
                    }
                }
                let trail = (next() % 12) as usize;
                query.extend(crate::testing::random_codes(&mut next, trail));
                // Frame diagonal of the source's placement, widened unevenly.
                let diagonal = lead as i32 + 1 - from as i32 - frames.offset(source) as i32;
                let (below, above) = ((next() % 6) as i32, (next() % 6) as i32);
                let unbanded = next().is_multiple_of(4);
                let count = if next().is_multiple_of(5) {
                    1 + (next() % (keep.top as u64 + 1)) as usize
                } else {
                    keep.top + 2 + (next() % 8) as usize
                };
                let mut chosen: Vec<AlleleId> = (0..alleles).collect();
                for k in 0..chosen.len() {
                    let other = k + (next() as usize) % (chosen.len() - k);
                    chosen.swap(k, other);
                }
                chosen.truncate(count.min(chosen.len()));
                let request: Vec<Candidate> = chosen
                    .iter()
                    .map(|&allele| {
                        let offset = frames.offset(allele) as i32;
                        let skew = i32::from(next().is_multiple_of(10));
                        let band = if unbanded {
                            Band::Full
                        } else {
                            Band::Diagonals {
                                lo: diagonal - below + offset + skew,
                                hi: diagonal + above + offset,
                            }
                        };
                        Candidate {
                            allele,
                            band,
                            votes: (next() % 4) as u32,
                        }
                    })
                    .collect();
                queries.push(query);
                candidates.push(request);
            }
            // The stage's precondition: the outer flank lies in the window,
            // so a V (prefix) window starts at the read's start.
            let prefix = outer == OuterExtend::Prefix;
            let starts: Vec<usize> = (0..3)
                .map(|_| (next() % 4) as usize)
                .map(|start| if prefix { 0 } else { start })
                .collect();
            let run = |bounded: bool, aligner: &Recording| {
                let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
                stage.accept = accept;
                stage.outer = outer;
                stage.ties = TieOrder::GeneLeadersFirst;
                stage.family_delta = bounded.then(|| table.clone());
                let requests: Vec<SegmentRequest<'_>> = (0..3)
                    .map(|read| SegmentRequest {
                        read,
                        query: &queries[read],
                        window: starts[read]..queries[read].len(),
                        candidates: &candidates[read],
                    })
                    .collect();
                let mut outs = vec![Vec::new(); 3];
                let (mut dp, mut batch) = (Default::default(), BatchScratch::default());
                stage.place_many(aligner, &requests, &mut dp, &mut batch, &mut outs);
                let placed: Vec<Vec<_>> = outs
                    .iter()
                    .map(|out| {
                        out.iter()
                            .map(|p| (p.alleles.clone(), p.alignment.clone()))
                            .collect()
                    })
                    .collect();
                (placed, batch.ranked[..3].to_vec())
            };
            let (full, bounded) = (Recording::default(), Recording::default());
            let (expected, all) = run(false, &full);
            let (placed, ranked) = run(true, &bounded);
            assert_eq!(placed, expected, "case {case}");
            // Brute force: align, accept and extend every candidate.
            for read in (0..3).filter(|_| keep.cap >= 64) {
                let window = starts[read]..queries[read].len();
                let mut kept = Vec::new();
                for candidate in &candidates[read] {
                    let profile = segment.profile(candidate.allele);
                    let query = &queries[read][window.clone()];
                    let mut ws = Default::default();
                    let alignment =
                        BandedGotoh.align(query, profile, segment.gaps(), candidate.band, &mut ws);
                    let Some(mut alignment) = alignment.filter(|a| a.score >= keep.min_score)
                    else {
                        continue;
                    };
                    alignment.shift_query(window.start as u32);
                    let allele = candidate.allele;
                    let placement =
                        Placement::new(SegmentKind::V, allele, alignment, &queries[read], &segment);
                    if accept.allows(&placement.evidence) {
                        let mut alignment = placement.alignment;
                        let (local, reference) = (alignment.score, segment.reference(allele));
                        extend_outer(&mut alignment, &queries[read], reference, profile, outer);
                        kept.push((allele, alignment.score, local));
                    }
                }
                // The top by extended score and the top by local score.
                let cutoff = |score: fn(&(AlleleId, Score, Score)) -> Score| {
                    let mut scores: Vec<Score> = kept.iter().map(score).collect();
                    nth_largest(&mut scores, keep.top).unwrap_or(Score::MIN)
                };
                let (extended, local) = (cutoff(|k| k.1), cutoff(|k| k.2));
                let mut kept: Vec<(AlleleId, Score)> = kept
                    .iter()
                    .filter(|k| k.1 >= extended || k.2 >= local)
                    .map(|k| (k.0, k.1))
                    .collect();
                kept.sort_unstable();
                let mut stage: Vec<(AlleleId, Score)> = expected[read]
                    .iter()
                    .flat_map(|(alleles, alignment)| alleles.iter().map(|&a| (a, alignment.score)))
                    .collect();
                stage.sort_unstable();
                assert_eq!(stage, kept, "case {case} read {read}");
            }
            let ran = bounded.0.into_inner();
            for (read, request) in candidates.iter().enumerate() {
                if request.len() <= keep.top + 1 {
                    assert_eq!(ranked[read], all[read], "case {case}");
                    continue;
                }
                // The bounded ranking is the full one, each entry with its
                // exact score, less some entries below the `top`-th score
                // (none under non-permissive acceptance).
                let mut scores: Vec<Score> = all[read].iter().map(|&(_, s)| s).collect();
                let cutoff = nth_largest(&mut scores, keep.top).unwrap_or(Score::MIN);
                let mut rest = ranked[read].iter().peekable();
                for entry in &all[read] {
                    if rest.next_if_eq(&entry).is_none() {
                        let permissive = accept.is_permissive();
                        assert!(permissive && entry.1 < cutoff, "case {case} read {read}");
                    }
                }
                assert_eq!(rest.next(), None, "case {case} read {read}");
                let query = queries[read][starts[read]..].as_ptr();
                for candidate in request {
                    let profile: *const _ = segment.profile(candidate.allele);
                    let runs = |traced| ran.contains(&(query, profile, traced));
                    pruned += usize::from(!runs(false) && !runs(true));
                    unscored += usize::from(!runs(false) && runs(true));
                }
            }
        }
        assert!(pruned > 100, "{pruned} siblings pruned");
        assert!(unscored > 10, "{unscored} siblings aligned without scoring");
    }

    /// A request's placements do not depend on the schedule: on the other
    /// requests of its batch and their order, or on bounding. Ties keep the
    /// key: gene leaders first in a scored request, retrieval order in one
    /// with at most `top + 1` candidates, which is traced without scoring.
    #[test]
    fn schedule_independence() {
        let (segment, frames) = frame_family();
        let table = FamilyBoundTable::with_frames(&segment, &frames);
        let keep = KeepPolicy {
            top: 2,
            min_score: 1,
            cap: 64,
        };
        // IGHV1-1*02 (allele 1) and its duplicate *08 (allele 7) tie.
        let twin = segment.set().allele(1).codes.clone();
        let other = segment.set().allele(9).codes.clone();
        let queries = [
            twin[10..80].to_vec(),
            twin[..70].to_vec(),
            other[5..85].to_vec(),
        ];
        let full = |alleles: &[(AlleleId, u32)]| -> Vec<Candidate> {
            alleles
                .iter()
                .map(|&(allele, votes)| Candidate {
                    allele,
                    band: Band::Full,
                    votes,
                })
                .collect()
        };
        let candidates = [
            full(&[(7, 1), (0, 3), (1, 5), (2, 2), (8, 4)]),
            full(&[(7, 1), (1, 5)]),
            full(&[(9, 2), (8, 6), (10, 1), (3, 3), (5, 2)]),
        ];
        let place = |reads: &[usize], bounded: bool| {
            let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
            stage.ties = TieOrder::GeneLeadersFirst;
            stage.family_delta = bounded.then(|| table.clone());
            let requests: Vec<SegmentRequest<'_>> = reads
                .iter()
                .map(|&read| SegmentRequest {
                    read,
                    query: &queries[read],
                    window: 0..queries[read].len(),
                    candidates: &candidates[read],
                })
                .collect();
            let mut outs = vec![Vec::new(); 3];
            let (mut dp, mut batch) = (Default::default(), BatchScratch::default());
            stage.place_many(&BandedGotoh, &requests, &mut dp, &mut batch, &mut outs);
            outs.iter()
                .map(|out| {
                    out.iter()
                        .map(|p| (p.alleles.clone(), p.alignment.clone()))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        let expected = place(&[0, 1, 2], false);
        for bounded in [false, true] {
            assert_eq!(place(&[2, 1, 0], bounded), expected, "{bounded}");
            for (read, expected) in expected.iter().enumerate() {
                assert_eq!(&place(&[read], bounded)[read], expected, "{bounded}");
            }
        }
        assert_eq!(expected[0][0].0, vec![1, 7]);
        assert_eq!(expected[1][0].0, vec![7, 1]);
    }

    #[test]
    fn significance_floor_follows_the_search_space() {
        let (stage, _, _) = v_stage(false);
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let significance =
            Significance::new(stage.segment.set(), &scoring, EvalueGates::default()).unwrap();
        // E = K·m·n·e^{−λS}: S grows with ln(m) and with −ln(E) / λ.
        assert!(significance.min_score(600, 1e-3) > significance.min_score(60, 1e-3));
        let gap = significance.min_score(300, 1e-3) - significance.min_score(300, 1.0);
        assert!((gap - (1e3f64.ln() / significance.statistics.lambda) as Score).abs() <= 1);
        // The gate reads the formula of the reported support.
        let strict = significance;
        let floor = strict.min_score(300, 1e-3);
        let evalue = |s: Score| strict.statistics.evalue(s, 300, strict.database_len);
        assert!(evalue(floor) <= 1e-3 && evalue(floor - 1) > 1e-3);
    }

    #[test]
    fn traceback_waves_stop_at_the_top_group_and_refill_after_a_rejection() {
        let policy = KeepPolicy {
            top: 2,
            min_score: 1,
            cap: 8,
        };
        let scores = [40, 30, 30, 12, 4];
        assert_eq!(policy.next_wave(&scores, 0, 0, None), 3);
        // Two accepted with cutoff 30: nothing below the tie group is traced.
        assert_eq!(policy.next_wave(&scores, 3, 2, Some(30)), 3);
        // One of the top scores was rejected: the next wave fills the open slot.
        assert_eq!(policy.next_wave(&scores, 3, 1, None), 4);
        assert_eq!(policy.next_wave(&scores, 4, 2, Some(12)), 4);
        // `cap` bounds the traced prefix.
        let capped = KeepPolicy { cap: 2, ..policy };
        assert_eq!(capped.next_wave(&scores, 0, 0, None), 2);
    }
}
