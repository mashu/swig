//! One segment kind's search stage: prepared database + retriever + keep policy.
//!
//! `SegmentStage` is the unit reused for V, D and J (and would be for C). For
//! every request (a query window and its retrieved candidates) it runs:
//!
//! ```text
//! candidates ─► score-only DP ─► (V) sibling bound pruning ─► rank by score
//!            ─► traceback in waves until `top` alignments pass acceptance (+ ties)
//!            ─► outer flank extension ─► rank ─► merge identical alignments
//! ```
//!
//! Every request of a batch is processed in the same waves, so a SIMD
//! aligner sees one score batch and one traceback batch per wave for the
//! whole chunk of reads, while every decision stays per request.

use std::ops::Range;
use std::sync::Arc;

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

/// Length and identity floor applied to a local alignment before outer extension.
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
/// The stage keeps the `top` best *accepted* local alignments and every
/// alignment tied with the last of them. A high score that fails
/// [`HitAcceptance`] does not take a slot, so traceback proceeds in waves
/// down the score-ranked list until `top` alignments are accepted, the list
/// is exhausted, or `cap` candidates were traced.
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
    /// `top` are accepted, `cutoff` is the worst accepted score and the wave
    /// takes only the remaining candidates tied with it.
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
    /// Per request: local scores of accepted alignments.
    accepted: Vec<Vec<Score>>,
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
        clear(&mut self.pending, requests.len());
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
///   is called only when its path contains an anchor ([`Placement::anchor`]),
///   so a random read — whose best V or J is never that significant — gets
///   no call, while a J downstream of a significant V needs only `accept`:
///   given the V, a J is expected there.
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

    /// Smallest local score with E-value ≤ `evalue` for a window of `query_len`.
    pub fn min_score(&self, query_len: usize, evalue: f64) -> Score {
        let KarlinAltschul { lambda, k, .. } = self.statistics;
        let space = k * query_len.max(1) as f64 * self.database_len.max(1) as f64;
        ((space.ln() - evalue.ln()) / lambda).ceil().max(0.0) as Score
    }
}

pub struct SegmentStage<R> {
    pub segment: Arc<PreparedSegment>,
    pub retriever: R,
    pub keep: KeepPolicy,
    pub accept: HitAcceptance,
    /// `None` accepts hits of any significance, and none anchors a call.
    pub significance: Option<Significance>,
    /// Stitch the co-linear V prefix or J suffix back on after traceback, then re-rank.
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
            self.finish(out, first, request.query);
        }
    }

    /// With at most `top + 1` candidates, tracing all of them is cheaper than
    /// scoring all of them and then tracing `top`.
    fn needs_scores(&self, request: &SegmentRequest<'_>) -> bool {
        request.candidates.len() > self.keep.top.saturating_add(1)
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
                        batch.owners.extend((0..count).map(|c| (index, c)));
                    }
                }
                self.score_owners(aligner, requests, dp, batch);
            }
            Some(table) => self.rank_with_bounds(aligner, requests, table, dp, batch),
        }
        for (index, request) in requests.iter().enumerate() {
            if self.needs_scores(request) {
                let leaders = &batch.leaders[batch.starts[index]..batch.starts[index + 1]];
                batch.ranked[index].sort_unstable_by_key(|&(candidate, score)| {
                    (std::cmp::Reverse(score), self.ties.key(leaders, candidate))
                });
            }
        }
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
    ///    * *pruned* when `U < floor` or `U < c`: it scores below `T`, so
    ///      it can enter neither the top `top` scores nor their ties, and
    ///      the ranking from `T` up is that of full scoring;
    ///    * *known* when `U = L`: that is its exact score, without DP;
    ///    * *traced now* when `L ≥ c`: it likely reaches the traceback, so
    ///      it is aligned instead of scored, its alignment's score is its
    ///      exact score ([`Aligner`]'s contract) and traceback reuses it;
    ///    * *scored* otherwise.
    ///
    /// Aligning representatives, the cutoff and early tracing apply only
    /// under permissive acceptance ([`HitAcceptance::is_permissive`]).
    /// Otherwise a rejected alignment makes the traceback waves refill from
    /// below `T`, where a pruned sibling may lie, so siblings are pruned only
    /// below the floor.
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
        let cut = self.accept.is_permissive() && self.keep.top > 0;
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
                    Some(upper) if upper < floor || !reaches(upper) => {}
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

    /// Align the siblings in `batch.early`: record each exact score, rank
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

    /// Traceback waves until every request has `top` accepted alignments
    /// (plus ties), runs out of candidates, or reaches `cap`. Alignments
    /// cached while ranking are taken, not recomputed.
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
            batch.owners.clear();
            for (index, cursor) in batch.cursors.iter_mut().enumerate() {
                if *cursor == DONE {
                    continue;
                }
                let ranked = &batch.ranked[index];
                batch.scores.clear();
                batch.scores.extend(ranked.iter().map(|&(_, s)| s));
                // `accepted` is kept in descending order.
                let accepted = &batch.accepted[index];
                let cutoff = self
                    .keep
                    .top
                    .checked_sub(1)
                    .and_then(|k| accepted.get(k))
                    .copied();
                let end = self
                    .keep
                    .next_wave(&batch.scores, *cursor, accepted.len(), cutoff);
                batch
                    .owners
                    .extend(ranked[*cursor..end].iter().map(|&(c, _)| (index, c)));
                // An empty wave changes nothing the next one depends on.
                *cursor = if end == *cursor { DONE } else { end };
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

    /// Place a traced alignment of candidate `candidate` of `request` and
    /// keep it when it clears the floor and acceptance.
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
            .is_some_and(|s| local >= s.min_score(request.window.len(), s.gates.anchor));
        if self.accept.allows(&placement.evidence) {
            let at = accepted.partition_point(|&s| s >= local);
            accepted.insert(at, local);
            pending.push(placement);
        }
    }

    /// Keep the top accepted alignments and their ties (by local score),
    /// extend their outer flanks, rank, and merge identical alignments.
    fn finish(&self, out: &mut Vec<Placement>, first: usize, query: &[u8]) {
        let rank = |p: &Placement| std::cmp::Reverse((p.score(), p.evidence.longest_exact_run));
        out[first..].sort_by_key(rank);
        if self.keep.top > 0 && out.len() - first > self.keep.top {
            let cutoff = out[first + self.keep.top - 1].score();
            let kept = out[first..]
                .iter()
                .take_while(|p| p.score() >= cutoff)
                .count();
            out.truncate(first + kept);
        }
        for placement in &mut out[first..] {
            self.extend_outer(placement, query);
        }
        out[first..].sort_by_key(rank);
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
    fn placed<R: Retriever>(
        stage: &SegmentStage<R>,
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
    /// inconsistent band, every `top` and `cap`, and both acceptance kinds;
    /// and it skips the DP of some siblings entirely and aligns others
    /// without scoring them (known exactly, or traced early).
    #[test]
    fn bounded_ranking_equals_brute_force() {
        let (segment, frames) = frame_family();
        let table = FamilyBoundTable::with_frames(&segment, &frames);
        let alleles = segment.set().len() as AlleleId;
        let mut next = crate::testing::rng(23);
        let (mut pruned, mut unscored) = (0, 0);
        for case in 0..240 {
            let keep = KeepPolicy {
                top: [1, 2, 3, 5][case % 4],
                min_score: 1,
                cap: [2, 3, 64][case / 4 % 3],
            };
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
            let starts: Vec<usize> = (0..3).map(|_| (next() % 4) as usize).collect();
            let run = |bounded: bool, aligner: &Recording| {
                let mut stage = SegmentStage::new(Arc::clone(&segment), EmptyRetriever, keep);
                stage.accept = accept;
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
        let strict = significance;
        let floor = strict.min_score(300, 1e-3);
        let evalue = |s: Score| {
            strict.statistics.k
                * 300.0
                * strict.database_len as f64
                * (-strict.statistics.lambda * f64::from(s)).exp()
        };
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
