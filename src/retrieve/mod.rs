//! Candidate retrieval: which alleles are worth a DP, and where (band).
//!
//! A [`Retriever`] turns a query window into [`Candidate`]s. Retrieval is the
//! part most likely to be swapped (spaced seeds, minimizers, allele-family
//! trees, learned retrieval), so it knows nothing about alignment and the
//! aligner knows nothing about seeds — the [`Band`] is the whole contract.

mod budget;
mod lmer_mask;
mod seed_vote;
mod ungapped;

pub use budget::{AmbiguityBudget, CandidateBudget, FixedBudget};
pub use lmer_mask::{LmerMaskRetriever, MAX_LMER};
pub use seed_vote::{SeedVoteParams, SeedVoteRetriever};
pub use ungapped::{end_anchored_ungapped, UngappedParams, UngappedPrefilter};

use crate::band::Band;
use crate::germline::{AlleleId, GeneId};
use crate::scoring::Score;

/// An allele worth aligning, with the band to align it in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Candidate {
    pub allele: AlleleId,
    pub band: Band,
    /// Retrieval strength (weighted seed votes or shared L-mers).
    pub votes: u32,
}

/// Best diagonal cluster of one allele.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlleleVote {
    pub allele: AlleleId,
    pub votes: u32,
    /// Smallest and largest seed diagonal (`query_pos − ref_pos`) in the cluster.
    pub d_lo: i32,
    pub d_hi: i32,
}

/// Reusable per-thread retrieval buffers. Retrievers may use any subset.
#[derive(Debug, Default)]
pub struct RetrievalScratch {
    pub keys: Vec<u64>,
    pub votes: Vec<AlleleVote>,
    pub gene_votes: Vec<u32>,
    pub gene_order: Vec<GeneId>,
    pub ranked: Vec<u32>,
    pub counts: Vec<u32>,
    pub ranked_candidates: Vec<(Score, Candidate)>,
}

pub trait Retriever: Send + Sync {
    /// Fill `out` (cleared first) with candidates for `query`, strongest first.
    /// Band diagonals are relative to `query` (the window slice).
    fn retrieve(&self, query: &[u8], scratch: &mut RetrievalScratch, out: &mut Vec<Candidate>);
}

/// Strongest candidate's votes (orientation and weakness decisions).
pub fn strongest(candidates: &[Candidate]) -> u32 {
    candidates.iter().map(|c| c.votes).max().unwrap_or(0)
}
