//! Candidate retrieval: which alleles are worth a DP, and where (band).
//!
//! A [`Retriever`] turns a query window into [`Candidate`]s. Retrieval is the
//! part most likely to be swapped (spaced seeds, minimizers, allele-family
//! trees, learned retrieval), so it knows nothing about alignment and the
//! aligner knows nothing about seeds — the [`Band`] is the whole contract.

mod budget;
mod exhaustive;
mod frames;
mod lmer_mask;
mod reach;
mod seed_vote;
mod ungapped;

pub use budget::{AmbiguityBudget, CandidateBudget, FixedBudget};
pub use exhaustive::ExhaustiveRetriever;
pub use frames::GeneFrames;
pub use lmer_mask::{LmerMaskRetriever, MAX_LMER};
pub use reach::ReachClasses;
pub use seed_vote::{SeedVoteParams, SeedVoteRetriever};
pub use ungapped::{UngappedParams, UngappedPrefilter};

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

/// Best diagonal cluster of one indexed sequence (a gene frame), spanning
/// the far side of an indel when the frame's runner-up cluster is one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ClusterVote {
    pub target: u32,
    pub votes: u32,
    /// Smallest and largest seed diagonal (`query_pos − ref_pos`) in the cluster.
    pub d_lo: i32,
    pub d_hi: i32,
}

/// Reusable per-thread retrieval buffers. Retrievers may use any subset.
#[derive(Debug, Default)]
pub struct RetrievalScratch {
    pub keys: Vec<u64>,
    pub buffer: Vec<u64>,
    pub votes: Vec<ClusterVote>,
    pub gene_votes: Vec<u32>,
    pub gene_order: Vec<GeneId>,
    pub ranked: Vec<u32>,
    pub counts: Vec<u32>,
    pub ranked_candidates: Vec<(Score, Candidate)>,
    /// Sampled query seeds and their s-mer ranks ([`crate::seed::SeedSampling`]).
    pub seeds: Vec<(usize, u64)>,
    pub ranks: Vec<u64>,
}

pub trait Retriever: Send + Sync {
    /// Fill `out` (cleared first) with candidates for `query`, strongest first.
    /// Band diagonals are relative to `query` (the window slice).
    fn retrieve(&self, query: &[u8], scratch: &mut RetrievalScratch, out: &mut Vec<Candidate>);
}
