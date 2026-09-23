//! One segment kind's search stage: prepared database + retriever + keep policy.
//!
//! `SegmentStage` is the unit reused for V, D and J (and would be for C): it
//! retrieves candidates, scores them score-only, keeps the best few, and runs
//! traceback only for those.

use std::ops::Range;
use std::sync::Arc;

use crate::align::{AlignJob, Aligner, Alignment, ScoreHit};
use crate::germline::{AlleleId, GermlineSet, SegmentKind};
use crate::placement::{merge_identical, Placement};
use crate::retrieve::{Candidate, RetrievalScratch, Retriever};
use crate::scoring::{GapCosts, Score, SubstitutionModel, SubstitutionProfile};

/// A germline set with per-allele substitution profiles and gap costs.
#[derive(Clone, Debug)]
pub struct PreparedSegment {
    set: GermlineSet,
    profiles: Vec<SubstitutionProfile>,
    gaps: GapCosts,
}

impl PreparedSegment {
    pub fn new<M: SubstitutionModel + ?Sized>(set: GermlineSet, model: &M, gaps: GapCosts) -> Self {
        let profiles = set.alleles().iter().map(|a| SubstitutionProfile::build(model, &a.codes)).collect();
        Self { set, profiles, gaps }
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
}

impl KeepPolicy {
    /// Whether `candidates` retrieved candidates warrant the score-only pre-pass.
    pub fn prescore(&self, candidates: usize) -> bool {
        candidates > self.prescore_above
    }

    /// Sort `items` by score (descending, stable) and truncate to the kept prefix.
    pub fn select<T>(&self, items: &mut Vec<T>, score: impl Fn(&T) -> Score) {
        items.sort_by_key(|item| std::cmp::Reverse(score(item)));
        let Some(best) = items.first().map(&score) else { return };
        let kept = items
            .iter()
            .enumerate()
            .take_while(|&(rank, item)| {
                let s = score(item);
                rank < self.cap && s >= self.min_score && (rank < self.top || s >= best - self.slack)
            })
            .count();
        items.truncate(kept);
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
}

impl<R: Retriever> SegmentStage<R> {
    pub fn new(segment: Arc<PreparedSegment>, retriever: R, keep: KeepPolicy) -> Self {
        Self { segment, retriever, keep }
    }

    pub fn kind(&self) -> SegmentKind {
        self.segment.kind()
    }

    /// Retrieve candidates for `query[window]` (bands relative to the window).
    pub fn retrieve(&self, query: &[u8], window: Range<usize>, scratch: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
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
        let StageScratch { retrieval, dp, candidates, batch } = scratch;
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
        let request = SegmentRequest { read: 0, query, window, candidates };
        self.place_many(aligner, std::slice::from_ref(&request), dp, batch, std::slice::from_mut(out));
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
        let job = |request: &SegmentRequest<'a>, candidate: &Candidate| self.job(request, candidate);
        let BatchScratch { scored, hits, alignments, kept, bounds } = batch;

        // 1. Score-only pre-pass for requests with more candidates than they keep.
        let mut jobs: Vec<AlignJob<'a>> = Vec::new();
        for request in requests.iter().filter(|r| self.keep.prescore(r.candidates.len())) {
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
                    request.candidates.iter().zip(hit_results.by_ref()).filter_map(|(c, hit)| hit.map(|hit| (*c, hit))),
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
            jobs.extend(kept[bounds[k]..bounds[k + 1]].iter().map(|c| job(request, c)));
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
                    out.push(Placement::new(self.kind(), candidate.allele, alignment, request.query, &self.segment));
                }
            }
            if !self.keep.prescore(request.candidates.len()) {
                // Every candidate was kept: apply the same ordering and floor
                // that the pre-pass would have applied.
                let fresh = &mut out[first_new..];
                fresh.sort_by_key(|p| std::cmp::Reverse(p.score()));
                let floor = fresh.iter().take_while(|p| p.score() >= self.keep.min_score).count();
                out.truncate(first_new + floor);
            }
            merge_identical(out, first_new);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keep_policy_combines_top_slack_floor_and_cap() {
        let policy = KeepPolicy { top: 2, slack: 10, min_score: 20, cap: 4, prescore_above: 0 };
        let mut scores = vec![50, 15, 45, 30, 42, 41];
        policy.select(&mut scores, |&s| s);
        assert_eq!(scores, vec![50, 45, 42, 41]);
        let mut low = vec![25, 12, 21];
        policy.select(&mut low, |&s| s);
        assert_eq!(low, vec![25, 21]);
    }
}
