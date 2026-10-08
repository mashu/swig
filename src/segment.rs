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
use crate::family_bound::FamilyBoundTable;
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
    /// Per gene id: the leader found so far in the current request.
    gene_leader: Vec<u32>,
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
        self.exact.clear();
        self.exact.resize(total, 0);
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
    /// Score one allele per gene first and skip same-length siblings whose
    /// upper bound cannot enter the top `keep.top` scores.
    pub family_delta: Option<FamilyBoundTable>,
}

/// No candidate (a per-gene slot of [`BatchScratch`]'s leader search).
const NONE: u32 = u32::MAX;

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
    /// a candidate index. `best` is per-gene scratch, all [`NONE`] on entry
    /// and on return.
    fn find_leaders(&self, candidates: &[Candidate], leaders: &mut [u32], best: &mut Vec<u32>) {
        debug_assert!(
            (1..candidates.len()).all(|i| candidates[..i]
                .iter()
                .all(|c| c.allele != candidates[i].allele)),
            "a request lists each allele once"
        );
        let set = self.segment.set();
        if best.len() < set.gene_count() {
            best.resize(set.gene_count(), NONE);
        }
        let gene = |c: &Candidate| set.allele(c.allele).gene as usize;
        let strength = |c: &Candidate| (c.votes, std::cmp::Reverse(c.allele));
        for (index, candidate) in (0..).zip(candidates) {
            let slot = &mut best[gene(candidate)];
            if *slot == NONE || strength(candidate) > strength(&candidates[*slot as usize]) {
                *slot = index;
            }
        }
        for (leader, candidate) in leaders.iter_mut().zip(candidates) {
            *leader = best[gene(candidate)];
        }
        for candidate in candidates {
            best[gene(candidate)] = NONE;
        }
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
                    &mut batch.gene_leader,
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

    /// Family-delta ranking: score each gene's leader (its representative),
    /// then only the siblings whose sound upper bound `S(rep) + Δ(rep, sibling)`
    /// reaches the `top`-th best representative score (see [`crate::family_bound`]).
    /// The cutoff can only rise as more alleles are scored, so a sibling below
    /// it can never enter the top scores. A sibling in another band than its
    /// representative has no bound and is always scored.
    ///
    /// The cutoff applies only under permissive acceptance
    /// ([`HitAcceptance::is_permissive`]). Otherwise a rejected alignment makes
    /// the traceback waves refill from below the `top`-th score, where a pruned
    /// sibling may lie, so siblings are pruned only below the floor.
    fn rank_with_bounds<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        table: &FamilyBoundTable,
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        // Representative scores are recorded even below the floor: a weak
        // representative can still bound a sibling that clears it.
        batch.owners.clear();
        for (index, request) in requests.iter().enumerate() {
            if self.needs_scores(request) {
                let leaders = &batch.leaders[batch.starts[index]..batch.starts[index + 1]];
                batch.owners.extend(
                    (0..)
                        .zip(leaders)
                        .filter(|&(c, &leader)| c == leader)
                        .map(|(c, _)| (index, c)),
                );
            }
        }
        self.score_owners(aligner, requests, dp, batch);
        let permissive = self.accept.is_permissive();
        batch.owners.clear();
        for (index, request) in requests.iter().enumerate() {
            if !self.needs_scores(request) {
                continue;
            }
            let cutoff = permissive
                .then(|| top_score(&batch.ranked[index], self.keep.top, &mut batch.scores))
                .flatten();
            let floor = self.floor(request);
            let start = batch.starts[index];
            for (c, candidate) in (0..).zip(request.candidates) {
                let leader = batch.leaders[start + c as usize];
                if leader == c {
                    continue;
                }
                let rep = &request.candidates[leader as usize];
                // No positive cell means the representative's optimum is 0.
                let rep_score = batch.exact[start + leader as usize];
                let bound = (rep.band == candidate.band)
                    .then(|| table.upper_bound(rep.allele, candidate.allele, rep_score))
                    .flatten();
                let prunable =
                    bound.is_some_and(|b| b < floor || cutoff.is_some_and(|cutoff| b < cutoff));
                if !prunable {
                    batch.owners.push((index, c));
                }
            }
        }
        self.score_owners(aligner, requests, dp, batch);
    }

    /// Traceback waves until every request has `top` accepted alignments
    /// (plus ties), runs out of candidates, or reaches `cap`.
    fn trace<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
    ) {
        let mut cursors = vec![0usize; requests.len()];
        loop {
            batch.owners.clear();
            for (index, cursor) in cursors.iter_mut().enumerate() {
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
                *cursor = end;
            }
            if batch.owners.is_empty() {
                break;
            }
            let jobs: Vec<AlignJob<'a>> = batch
                .owners
                .iter()
                .map(|&(index, candidate)| self.job(&requests[index], candidate))
                .collect();
            aligner.align_batch(&jobs, dp, &mut batch.alignments);
            for (&(index, candidate), alignment) in
                batch.owners.iter().zip(batch.alignments.drain(..))
            {
                let Some(mut alignment) = alignment else {
                    continue;
                };
                let request = &requests[index];
                if alignment.score < self.floor(request) {
                    continue;
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
                    let accepted = &mut batch.accepted[index];
                    let at = accepted.partition_point(|&s| s >= local);
                    accepted.insert(at, local);
                    batch.pending[index].push(placement);
                }
            }
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

/// The `top`-th best score of `ranked`, if it has that many entries.
fn top_score(ranked: &[(u32, Score)], top: usize, buffer: &mut Vec<Score>) -> Option<Score> {
    if top == 0 || ranked.len() < top {
        return None;
    }
    buffer.clear();
    buffer.extend(ranked.iter().map(|&(_, s)| s));
    let (_, &mut nth, _) = buffer.select_nth_unstable_by_key(top - 1, |&s| std::cmp::Reverse(s));
    Some(nth)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::BandedGotoh;
    use crate::alphabet::{decode, encode};
    use crate::band::Band;
    use crate::family_bound::FamilyBoundTable;
    use crate::retrieve::Retriever;
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
        let bound = table.upper_bound(1, 2, score(1)).unwrap();
        assert!(score(0) > score(2) && score(2) > score(1) && bound < score(0));
        let full = placed(&stage(false), &candidates, &query);
        assert_eq!(full, vec![(vec![2], score(2))]);
        assert_eq!(placed(&stage(true), &candidates, &query), full);
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
