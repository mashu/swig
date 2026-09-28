//! One segment kind's search stage: prepared database + retriever + keep policy.
//!
//! `SegmentStage` is the unit reused for V, D and J (and would be for C): it
//! retrieves candidates, scores them score-only, keeps the best few, and runs
//! traceback only for those.

use std::ops::Range;
use std::sync::Arc;

use crate::align::{AlignJob, Aligner, Alignment, ScoreHit};
use crate::evidence::Evidence;
use crate::germline::{AlleleId, GermlineSet, SegmentKind};
use crate::placement::{merge_identical, Placement};
use crate::refine::{extend_outer, OuterExtend};
use crate::retrieve::{Candidate, RetrievalScratch, Retriever};
use crate::scoring::{GapCosts, Score, SubstitutionModel, SubstitutionProfile};

/// Length and identity floor applied to a local alignment before outer extension.
///
/// Swig drops a hit below `min_subst` substitution columns (matches plus
/// mismatches) or below `min_identity`. `min_identity` of 0 admits every hit
/// that meets the length floor.
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

    pub fn gaps(&self) -> GapCosts {
        self.gaps
    }
}

/// Which score-only hits proceed to traceback.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KeepPolicy {
    /// Always keep this many best hits (if they reach `min_score`)...
    pub top: usize,
    /// ...plus any hit within `slack` of the best...
    pub slack: Score,
    /// ...never below `min_score`...
    pub min_score: Score,
    /// ...and never more than `cap` in total.
    pub cap: usize,
    /// Run the score-only pre-pass only when more candidates than this were
    /// retrieved. With `top` or fewer candidates the pre-pass can prune only
    /// through `min_score`, so for V/J (where that floor rarely binds) aligning
    /// directly saves a full DP per candidate; D keeps `0` because its
    /// lossless score floor prunes most candidates cheaply.
    pub prescore_above: usize,
    /// Swig's score screen. Score every candidate, then traceback only until
    /// `top` alignments pass acceptance, plus the exact-score tie group at that
    /// cutoff. A high score that fails acceptance does not fill a slot.
    pub screen: bool,
}

impl KeepPolicy {
    /// Whether `candidates` retrieved candidates warrant the score-only pre-pass.
    pub fn prescore(&self, candidates: usize) -> bool {
        candidates > self.prescore_above
    }

    /// Sort `items` by score (descending, stable) and truncate to the kept prefix.
    pub fn select<T>(&self, items: &mut Vec<T>, score: impl Fn(&T) -> Score) {
        items.sort_by_key(|item| std::cmp::Reverse(score(item)));
        let Some(best) = items.first().map(&score) else {
            return;
        };
        let kept = items
            .iter()
            .enumerate()
            .take_while(|&(rank, item)| {
                let s = score(item);
                rank < self.cap
                    && s >= self.min_score
                    && (rank < self.top || s >= best - self.slack)
            })
            .count();
        items.truncate(kept);
    }

    /// Exclusive end index of the next traceback wave in a score-sorted list.
    ///
    /// While fewer than `top` hits have been accepted, take that many more
    /// scores and the tie group of the last one. Once `top` are accepted,
    /// `cutoff` is the worst of those scores and only the remaining tie group
    /// at or above it is traced.
    pub fn screen_end(
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
        if let Some(cutoff) = cutoff {
            let mut index = from;
            while index < limit && scores[index] >= cutoff && scores[index] >= self.min_score {
                index += 1;
            }
            return index;
        }
        let need = self.top.saturating_sub(accepted);
        let mut index = from;
        let mut taken = 0usize;
        let mut last = None;
        while index < limit && scores[index] >= self.min_score {
            if taken >= need && Some(scores[index]) != last {
                break;
            }
            last = Some(scores[index]);
            taken += 1;
            index += 1;
        }
        index
    }

    /// Worst score still inside the top `top` accepted hits.
    pub fn screen_cutoff(&self, accepted: &[Score]) -> Option<Score> {
        if self.top == 0 || accepted.len() < self.top {
            return None;
        }
        let mut ranked = accepted.to_vec();
        ranked.sort_unstable_by_key(|score| std::cmp::Reverse(*score));
        Some(ranked[self.top - 1])
    }
}

/// Per-thread buffers of one stage.
#[derive(Debug, Default)]
pub struct StageScratch<W> {
    pub retrieval: RetrievalScratch,
    pub dp: W,
    pub candidates: Vec<Candidate>,
    pub batch: BatchScratch,
}

/// Result buffers of batched alignment calls.
#[derive(Debug, Default)]
pub struct BatchScratch {
    pub scored: Vec<(Candidate, ScoreHit)>,
    pub hits: Vec<Option<ScoreHit>>,
    pub alignments: Vec<Option<Alignment>>,
    kept: Vec<Candidate>,
    bounds: Vec<usize>,
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

pub struct SegmentStage<R> {
    pub segment: Arc<PreparedSegment>,
    pub retriever: R,
    pub keep: KeepPolicy,
    pub accept: HitAcceptance,
    /// Stitch the co-linear V prefix or J suffix back on after traceback, then re-rank.
    pub outer: OuterExtend,
}

impl<R: Retriever> SegmentStage<R> {
    pub fn new(segment: Arc<PreparedSegment>, retriever: R, keep: KeepPolicy) -> Self {
        Self {
            segment,
            retriever,
            keep,
            accept: HitAcceptance::PERMISSIVE,
            outer: OuterExtend::None,
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

    /// Retrieve, score, keep, trace back; append placements (absolute coordinates) to `out`.
    pub fn place<A: Aligner>(
        &self,
        aligner: &A,
        query: &[u8],
        window: Range<usize>,
        scratch: &mut StageScratch<A::Workspace>,
        out: &mut Vec<Placement>,
    ) {
        let StageScratch {
            retrieval,
            dp,
            candidates,
            batch,
        } = scratch;
        self.retrieve(query, window.clone(), retrieval, candidates);
        self.place_candidates(aligner, candidates, query, window, dp, batch, out);
    }

    /// Score-only pass over `candidates`, keep policy, traceback for survivors.
    #[allow(clippy::too_many_arguments)]
    fn job<'a>(&'a self, request: &SegmentRequest<'a>, candidate: &Candidate) -> AlignJob<'a> {
        AlignJob {
            query: &request.query[request.window.clone()],
            profile: self.segment.profile(candidate.allele),
            gaps: self.segment.gaps(),
            band: candidate.band,
        }
    }

    /// Align `candidates` in `query[window]` and append the kept placements.
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
    /// score batch and one traceback batch, so SIMD aligners fill their lanes.
    ///
    /// Per request the result is exactly that of aligning it alone: the keep
    /// policy is applied per request, and batch aligners must answer each job
    /// as their single-job calls would.
    pub fn place_many<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
        outs: &mut [Vec<Placement>],
    ) {
        if self.keep.screen {
            self.place_screened(aligner, requests, dp, batch, outs);
            return;
        }
        let job =
            |request: &SegmentRequest<'a>, candidate: &Candidate| self.job(request, candidate);
        let BatchScratch {
            scored,
            hits,
            alignments,
            kept,
            bounds,
        } = batch;

        // 1. Score-only pre-pass for requests with more candidates than they keep.
        let mut jobs: Vec<AlignJob<'a>> = Vec::new();
        for request in requests
            .iter()
            .filter(|r| self.keep.prescore(r.candidates.len()))
        {
            jobs.extend(request.candidates.iter().map(|c| job(request, c)));
        }
        hits.clear();
        if !jobs.is_empty() {
            aligner.score_batch(&jobs, dp, hits);
        }
        kept.clear();
        bounds.clear();
        let mut hit_results = hits.iter();
        for request in requests {
            bounds.push(kept.len());
            if self.keep.prescore(request.candidates.len()) {
                scored.clear();
                scored.extend(
                    request
                        .candidates
                        .iter()
                        .zip(hit_results.by_ref())
                        .filter_map(|(c, hit)| hit.map(|hit| (*c, hit))),
                );
                self.keep.select(scored, |(_, hit)| hit.score);
                kept.extend(scored.iter().map(|(c, _)| *c));
            } else {
                kept.extend_from_slice(request.candidates);
            }
        }
        bounds.push(kept.len());

        // 2. One traceback batch over every kept candidate.
        jobs.clear();
        for (k, request) in requests.iter().enumerate() {
            jobs.extend(
                kept[bounds[k]..bounds[k + 1]]
                    .iter()
                    .map(|c| job(request, c)),
            );
        }
        alignments.clear();
        if !jobs.is_empty() {
            aligner.align_batch(&jobs, dp, alignments);
        }

        // 3. Placements per request, in request order.
        let mut results = alignments.drain(..);
        for (k, request) in requests.iter().enumerate() {
            let out = &mut outs[request.read];
            let first_new = out.len();
            for candidate in &kept[bounds[k]..bounds[k + 1]] {
                if let Some(mut alignment) = results.next().expect("one alignment result per job") {
                    alignment.shift_query(request.window.start as u32);
                    if alignment.score < self.keep.min_score {
                        continue;
                    }
                    let mut placement = Placement::new(
                        self.kind(),
                        candidate.allele,
                        alignment,
                        request.query,
                        &self.segment,
                    );
                    if !self.accept.allows(&placement.evidence) {
                        continue;
                    }
                    // Rank on the extended score. Swig traces by the local score,
                    // then stitches the outer flank and re-sorts.
                    self.extend_outer(&mut placement, request.query);
                    out.push(placement);
                }
            }
            out[first_new..].sort_by_key(|p| std::cmp::Reverse(p.score()));
            merge_identical(out, first_new);
        }
    }

    /// Score every candidate, then traceback only the prefix that can enter the
    /// top `keep.top` accepted alignments. One wave is enough when those
    /// alignments pass the length and identity floors; a rejected high score
    /// opens another wave.
    fn place_screened<'a, A: Aligner>(
        &'a self,
        aligner: &A,
        requests: &[SegmentRequest<'a>],
        dp: &mut A::Workspace,
        batch: &mut BatchScratch,
        outs: &mut [Vec<Placement>],
    ) {
        let BatchScratch {
            hits, alignments, ..
        } = batch;
        let mut score_jobs: Vec<AlignJob<'a>> = Vec::new();
        let mut score_owner: Vec<(usize, usize)> = Vec::new();
        for (request_index, request) in requests.iter().enumerate() {
            if request.candidates.len() > self.keep.top.saturating_add(1) {
                for (candidate_index, candidate) in request.candidates.iter().enumerate() {
                    score_jobs.push(self.job(request, candidate));
                    score_owner.push((request_index, candidate_index));
                }
            }
        }
        hits.clear();
        if !score_jobs.is_empty() {
            aligner.score_batch(&score_jobs, dp, hits);
        }

        let mut grouped: Vec<Vec<(Candidate, Score)>> = vec![Vec::new(); requests.len()];
        for ((request_index, candidate_index), hit) in score_owner.iter().zip(hits.iter()) {
            if let Some(hit) = hit {
                if hit.score >= self.keep.min_score {
                    let candidate = requests[*request_index].candidates[*candidate_index];
                    grouped[*request_index].push((candidate, hit.score));
                }
            }
        }
        for (request_index, request) in requests.iter().enumerate() {
            if request.candidates.len() > self.keep.top.saturating_add(1) {
                grouped[request_index].sort_by_key(|(_, score)| std::cmp::Reverse(*score));
            } else {
                grouped[request_index] = request
                    .candidates
                    .iter()
                    .copied()
                    .map(|candidate| (candidate, Score::MAX))
                    .collect();
            }
        }

        let mut scored: Vec<(Candidate, Score)> = Vec::new();
        let mut bounds = Vec::with_capacity(requests.len() + 1);
        bounds.push(0);
        for group in &grouped {
            scored.extend(group.iter().copied());
            bounds.push(scored.len());
        }
        let mut cursors: Vec<usize> = bounds[..requests.len()].to_vec();
        let mut accepted: Vec<Vec<Score>> = vec![Vec::new(); requests.len()];
        let mut pending: Vec<Vec<Placement>> = vec![Vec::new(); requests.len()];
        let mut wave: Vec<(usize, Candidate)> = Vec::new();
        let mut jobs: Vec<AlignJob<'a>> = Vec::new();

        loop {
            wave.clear();
            jobs.clear();
            let mut progress = false;
            for request_index in 0..requests.len() {
                let start = bounds[request_index];
                let end = bounds[request_index + 1];
                let scores: Vec<Score> =
                    scored[start..end].iter().map(|(_, score)| *score).collect();
                let cursor = cursors[request_index] - start;
                let next = start
                    + self.keep.screen_end(
                        &scores,
                        cursor,
                        accepted[request_index].len(),
                        self.keep.screen_cutoff(&accepted[request_index]),
                    );
                if next <= cursors[request_index] {
                    continue;
                }
                progress = true;
                for (candidate, _) in &scored[cursors[request_index]..next] {
                    wave.push((request_index, *candidate));
                    jobs.push(self.job(&requests[request_index], candidate));
                }
                cursors[request_index] = next;
            }
            if !progress {
                break;
            }
            alignments.clear();
            aligner.align_batch(&jobs, dp, alignments);
            for ((request_index, candidate), alignment) in wave.iter().zip(alignments.drain(..)) {
                let Some(mut alignment) = alignment else {
                    continue;
                };
                let request = &requests[*request_index];
                alignment.shift_query(request.window.start as u32);
                if alignment.score < self.keep.min_score {
                    continue;
                }
                let local = alignment.score;
                let placement = Placement::new(
                    self.kind(),
                    candidate.allele,
                    alignment,
                    request.query,
                    &self.segment,
                );
                if !self.accept.allows(&placement.evidence) {
                    continue;
                }
                accepted[*request_index].push(local);
                pending[*request_index].push(placement);
            }
        }

        for (request_index, request) in requests.iter().enumerate() {
            let out = &mut outs[request.read];
            let first_new = out.len();
            out.append(&mut pending[request_index]);
            retain_top_ties(out, first_new, self.keep.top);
            for placement in &mut out[first_new..] {
                self.extend_outer(placement, request.query);
            }
            out[first_new..].sort_by_key(|placement| std::cmp::Reverse(placement.score()));
            merge_identical(out, first_new);
        }
    }

    fn extend_outer(&self, placement: &mut Placement, query: &[u8]) {
        if self.outer == OuterExtend::None {
            return;
        }
        let allele = placement.primary_allele();
        let reference = self.segment.set().allele(allele).codes.clone();
        let profile = self.segment.profile(allele).clone();
        if extend_outer(
            &mut placement.alignment,
            query,
            &reference,
            &profile,
            self.outer,
        ) {
            placement.clip = crate::clip::ClipIndex::build(
                &placement.alignment,
                query,
                &profile,
                self.segment.gaps(),
            );
            placement.evidence =
                Evidence::from_columns(placement.alignment.columns(), query, &reference);
        }
    }
}

/// Keep the first `top` placements and every later one tied with the last of those.
fn retain_top_ties(placements: &mut Vec<Placement>, from: usize, top: usize) {
    let slice = &mut placements[from..];
    slice.sort_by_key(|placement| std::cmp::Reverse(placement.score()));
    if top == 0 || slice.len() <= top {
        return;
    }
    let cutoff = slice[top - 1].score();
    let mut kept = top;
    while kept < slice.len() && slice[kept].score() == cutoff {
        kept += 1;
    }
    placements.truncate(from + kept);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keep_policy_combines_top_slack_floor_and_cap() {
        let policy = KeepPolicy {
            top: 2,
            slack: 10,
            min_score: 20,
            cap: 4,
            prescore_above: 0,
            screen: false,
        };
        let mut scores = vec![50, 15, 45, 30, 42, 41];
        policy.select(&mut scores, |&s| s);
        assert_eq!(scores, vec![50, 45, 42, 41]);
        let mut low = vec![25, 12, 21];
        policy.select(&mut low, |&s| s);
        assert_eq!(low, vec![25, 21]);
    }

    #[test]
    fn score_screen_stops_at_the_top_group_and_continues_after_a_gap() {
        let policy = KeepPolicy {
            top: 2,
            slack: 0,
            min_score: 1,
            cap: 8,
            prescore_above: 3,
            screen: true,
        };
        let scores = [40, 30, 30, 12, 4];
        assert_eq!(policy.screen_end(&scores, 0, 0, None), 3);
        assert_eq!(policy.screen_cutoff(&[40, 30]), Some(30));
        assert_eq!(policy.screen_end(&scores, 3, 2, Some(30)), 3);
        // One of the top scores was rejected, so the next wave fills the open slot
        // and keeps the tie at that new score.
        assert_eq!(policy.screen_end(&scores, 3, 1, None), 4);
        assert_eq!(policy.screen_cutoff(&[40, 12]), Some(12));
        assert_eq!(policy.screen_end(&scores, 4, 2, Some(12)), 4);
    }
}
