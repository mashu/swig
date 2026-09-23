//! Ungapped diagonal pre-ranking: a retriever decorator.
//!
//! **Trick — ungapped before gapped (BLAST staging).** Seed votes say *where*
//! an allele might align; they are a poor proxy for *how well*. Before paying
//! for a banded affine DP (`O(n·w)` cells), each candidate gets an ungapped
//! score: the maximum-sum segment (Kadane) along the few diagonals around its
//! band centre, read straight from the substitution profile — `O(n)` adds and
//! maxes, roughly 1/50 of the banded DP. Only candidates within `slack` of the
//! best ungapped score (and always the `keep_top` best) proceed to DP.
//!
//! Family siblings at ~15 % divergence differ from the true gene by ~200 score
//! units ungapped and are dropped; alleles of the same gene (a few SNPs apart)
//! stay. The gapped optimum is never below the ungapped score on a diagonal
//! inside the band, so a lower bound ranks candidates; the approximation is
//! that an indel present in only *some* candidates could shift their relative
//! ranking — a generous `slack` covers that, and `None` disables the stage.

use std::sync::Arc;

use super::{Candidate, RetrievalScratch, Retriever};
use crate::band::Band;
use crate::scoring::{Score, SubstitutionProfile};
use crate::segment::PreparedSegment;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UngappedParams {
    /// Always keep this many best candidates.
    pub keep_top: usize,
    /// Also keep any candidate within this many score units of the best.
    pub slack: Score,
    /// Diagonals scanned on each side of the band centre.
    pub radius: i32,
}

impl Default for UngappedParams {
    fn default() -> Self {
        Self { keep_top: 3, slack: 30, radius: 2 }
    }
}

pub struct UngappedPrefilter<R> {
    pub inner: R,
    segment: Arc<PreparedSegment>,
    /// `None` passes the inner retriever's candidates through unchanged.
    pub params: Option<UngappedParams>,
}

impl<R> UngappedPrefilter<R> {
    pub fn new(inner: R, segment: Arc<PreparedSegment>, params: Option<UngappedParams>) -> Self {
        Self { inner, segment, params }
    }
}

/// Best ungapped local segment score on diagonal `d = query_pos − ref_pos`.
#[inline]
fn diagonal_score(query: &[u8], profile: &SubstitutionProfile, d: i32) -> Score {
    let r_start = (-d).max(0) as usize;
    let q_start = d.max(0) as usize;
    let len = (query.len().saturating_sub(q_start)).min(profile.len().saturating_sub(r_start));
    let (mut current, mut best) = (0, 0);
    for k in 0..len {
        current = (current + profile.at(query[q_start + k], r_start + k)).max(0);
        best = best.max(current);
    }
    best
}

impl<R: Retriever> Retriever for UngappedPrefilter<R> {
    fn retrieve(&self, query: &[u8], scratch: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
        self.inner.retrieve(query, scratch, out);
        let Some(params) = self.params else { return };
        let ranked = &mut scratch.ranked_candidates;
        ranked.clear();
        ranked.extend(out.iter().map(|&candidate| {
            let profile = self.segment.profile(candidate.allele);
            let score = match candidate.band {
                Band::Diagonals { lo, hi } => {
                    let center = (lo + hi).div_euclid(2);
                    (center - params.radius..=center + params.radius)
                        .map(|d| diagonal_score(query, profile, d))
                        .max()
                        .unwrap_or(0)
                }
                // Unbanded candidates carry no geometry: never filter them.
                Band::Full => Score::MAX,
            };
            (score, candidate)
        }));
        ranked.sort_by(|a, b| b.0.cmp(&a.0).then(b.1.votes.cmp(&a.1.votes)));
        let best = ranked.first().map_or(0, |&(s, _)| s);
        out.clear();
        out.extend(
            ranked
                .iter()
                .enumerate()
                .take_while(|&(rank, &(score, _))| rank < params.keep_top || score >= best.saturating_sub(params.slack))
                .map(|(_, &(_, candidate))| candidate),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::scoring::AffineScoring;

    #[test]
    fn kadane_finds_the_best_ungapped_segment() {
        let scoring = AffineScoring::new(2, -3, -5, -1);
        let reference = encode(b"ACGTACGTAC");
        let profile = SubstitutionProfile::build(&scoring, &reference);
        // query = "GG" + reference with one mismatch in the middle: diagonal +2
        let query = encode(b"GGACGTTCGTAC");
        assert_eq!(diagonal_score(&query, &profile, 2), 9 * 2 - 3);
        assert!(diagonal_score(&query, &profile, 0) < 9 * 2 - 3);
    }
}
