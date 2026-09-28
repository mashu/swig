//! Seed voting with diagonal clustering (V and J retrieval).
//!
//! Every query seed that hits allele `a` at germline position `p` from query
//! position `i` votes for the diagonal `d = i − p` of `a`. Correct alleles
//! collect many votes on (nearly) one diagonal; the cluster's diagonal range
//! becomes the DP band.
//!
//! **Trick — sort-and-scan voting with packed keys.** Instead of a hash map
//! `allele → diagonal → votes`, each vote is one `u64`
//! `allele << 40 | (d + 2³¹) << 8 | weight`; the batch is sorted once
//! (`sort_unstable` on integers, cache-friendly and branch-light) and a single
//! linear scan groups votes by allele and by adjacent diagonal bins. No
//! allocation per query once the scratch vector has grown.
//!
//! **Trick — frequency cutoff.** Seeds with more than `max_occurrences`
//! postings (low-complexity or family-wide k-mers) carry little information
//! and dominate work; they are skipped, as in the method note.
//!
//! **Two-tier seeds.** Primary seeds vote with weight 4. Only when the best
//! primary cluster is weak are fallback (shorter) seeds added with weight 1.

use super::budget::AmbiguityBudget;
use super::{AlleleVote, Candidate, CandidateBudget, RetrievalScratch, Retriever};
use crate::band::Band;
use crate::germline::{GeneId, GermlineSet};
use crate::index::{IndexError, KmerIndex};
use crate::seed::SeedShape;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SeedVoteParams {
    pub primary_weight: u8,
    pub fallback_weight: u8,
    /// Add fallback seeds when the best primary cluster has fewer weighted votes.
    pub weak_votes: u32,
    /// Skip seeds with more postings than this.
    pub max_occurrences: usize,
    /// Diagonal bin width for clustering.
    pub bin_width: i32,
    /// Empty bins tolerated inside one cluster.
    pub cluster_gap_bins: i32,
    /// Band half-width = clamp(cluster span + pad, min, max).
    pub band_pad: i32,
    pub band_min: i32,
    pub band_max: i32,
    /// Ignore seed geometry and align unbanded (`w < 0` in the note).
    pub unbanded: bool,
    /// Minimum weighted votes for an allele to be a candidate.
    pub min_votes: u32,
    /// Safety pool: never return more alleles than this.
    pub max_alleles: usize,
}

impl Default for SeedVoteParams {
    fn default() -> Self {
        Self {
            primary_weight: 4,
            fallback_weight: 1,
            weak_votes: 16,
            max_occurrences: 256,
            bin_width: 4,
            cluster_gap_bins: 1,
            band_pad: 8,
            band_min: 28,
            band_max: 256,
            unbanded: false,
            min_votes: 1,
            max_alleles: 64,
        }
    }
}

const DIAGONAL_BIAS: i64 = 1 << 31;

#[inline]
fn pack(allele: u32, diagonal: i64, weight: u8) -> u64 {
    (u64::from(allele) << 40)
        | ((((diagonal + DIAGONAL_BIAS) as u64) & 0xFFFF_FFFF) << 8)
        | u64::from(weight)
}

#[inline]
fn unpack(key: u64) -> (u32, i32, u32) {
    let allele = (key >> 40) as u32;
    let diagonal = (((key >> 8) & 0xFFFF_FFFF) as i64 - DIAGONAL_BIAS) as i32;
    (allele, diagonal, (key & 0xFF) as u32)
}

pub struct SeedVoteRetriever<B = AmbiguityBudget> {
    primary: KmerIndex,
    fallback: Option<KmerIndex>,
    allele_gene: Vec<GeneId>,
    gene_count: usize,
    params: SeedVoteParams,
    budget: B,
}

impl<B: CandidateBudget> SeedVoteRetriever<B> {
    pub fn build(
        set: &GermlineSet,
        primary: SeedShape,
        fallback: Option<SeedShape>,
        params: SeedVoteParams,
        budget: B,
    ) -> Result<Self, IndexError> {
        Ok(Self {
            primary: KmerIndex::build(set, primary)?,
            fallback: fallback
                .map(|shape| KmerIndex::build(set, shape))
                .transpose()?,
            allele_gene: set.alleles().iter().map(|a| a.gene).collect(),
            gene_count: set.gene_count(),
            params,
            budget,
        })
    }

    pub fn params(&self) -> &SeedVoteParams {
        &self.params
    }

    fn cast_votes(&self, index: &KmerIndex, weight: u8, query: &[u8], keys: &mut Vec<u64>) {
        for (qpos, key) in index.shape().kmers(query) {
            let postings = index.postings(key);
            if postings.len() > self.params.max_occurrences {
                continue;
            }
            keys.extend(
                postings
                    .iter()
                    .map(|p| pack(p.allele(), qpos as i64 - i64::from(p.position()), weight)),
            );
        }
    }

    /// Sort votes and keep each allele's best diagonal cluster. Returns the best cluster's votes.
    fn cluster(&self, keys: &mut [u64], out: &mut Vec<AlleleVote>) -> u32 {
        keys.sort_unstable();
        out.clear();
        let bin = |d: i32| d.div_euclid(self.params.bin_width);
        let max_step = 1 + self.params.cluster_gap_bins;
        let mut strongest = 0;
        let mut i = 0;
        while i < keys.len() {
            let (allele, d0, w0) = unpack(keys[i]);
            let mut best = AlleleVote {
                allele,
                votes: 0,
                d_lo: d0,
                d_hi: d0,
            };
            let mut current = AlleleVote {
                allele,
                votes: w0,
                d_lo: d0,
                d_hi: d0,
            };
            i += 1;
            while i < keys.len() {
                let (a, d, w) = unpack(keys[i]);
                if a != allele {
                    break;
                }
                if bin(d) - bin(current.d_hi) <= max_step {
                    current.votes += w;
                    current.d_hi = d;
                } else {
                    if current.votes > best.votes {
                        best = current;
                    }
                    current = AlleleVote {
                        allele,
                        votes: w,
                        d_lo: d,
                        d_hi: d,
                    };
                }
                i += 1;
            }
            if current.votes > best.votes {
                best = current;
            }
            strongest = strongest.max(best.votes);
            out.push(best);
        }
        strongest
    }

    fn band_for(&self, vote: &AlleleVote) -> Band {
        if self.params.unbanded {
            return Band::Full;
        }
        let center = (vote.d_lo + vote.d_hi).div_euclid(2);
        let half = (vote.d_hi - vote.d_lo + self.params.band_pad)
            .clamp(self.params.band_min, self.params.band_max);
        Band::around(center, half)
    }
}

impl<B: CandidateBudget> Retriever for SeedVoteRetriever<B> {
    fn retrieve(&self, query: &[u8], scratch: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
        out.clear();
        let RetrievalScratch {
            keys,
            votes,
            gene_votes,
            gene_order,
            ranked,
            ..
        } = scratch;
        keys.clear();
        self.cast_votes(&self.primary, self.params.primary_weight, query, keys);
        let mut strongest = self.cluster(keys, votes);
        if let (true, Some(fallback)) = (strongest < self.params.weak_votes, &self.fallback) {
            self.cast_votes(fallback, self.params.fallback_weight, query, keys);
            strongest = self.cluster(keys, votes);
        }
        if strongest < self.params.min_votes {
            return;
        }

        // Collapse per gene (best allele's votes), rank genes, apply the budget.
        gene_votes.clear();
        gene_votes.resize(self.gene_count, 0);
        for v in votes.iter() {
            let g = self.allele_gene[v.allele as usize] as usize;
            gene_votes[g] = gene_votes[g].max(v.votes);
        }
        gene_order.clear();
        gene_order.extend(
            (0..self.gene_count as GeneId)
                .filter(|&g| gene_votes[g as usize] >= self.params.min_votes),
        );
        gene_order.sort_by(|&a, &b| {
            gene_votes[b as usize]
                .cmp(&gene_votes[a as usize])
                .then(a.cmp(&b))
        });
        ranked.clear();
        ranked.extend(gene_order.iter().map(|&g| gene_votes[g as usize]));
        let keep = self.budget.genes_to_keep(ranked);
        const KEPT: u32 = u32::MAX;
        for &g in &gene_order[..keep] {
            gene_votes[g as usize] = KEPT;
        }

        votes.retain(|v| {
            v.votes >= self.params.min_votes
                && gene_votes[self.allele_gene[v.allele as usize] as usize] == KEPT
        });
        votes.sort_by(|a, b| b.votes.cmp(&a.votes).then(a.allele.cmp(&b.allele)));
        votes.truncate(self.params.max_alleles);
        out.extend(votes.iter().map(|v| Candidate {
            allele: v.allele,
            band: self.band_for(v),
            votes: v.votes,
        }));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::germline::SegmentKind;

    #[test]
    fn packing_round_trips_negative_diagonals() {
        for &(a, d, w) in &[(0u32, -300i64, 4u8), (12345, 0, 1), (1 << 19, 4000, 255)] {
            let (a2, d2, w2) = unpack(pack(a, d, w));
            assert_eq!((a2, i64::from(d2), w2), (a, d, u32::from(w)));
        }
        assert!(pack(1, -5, 1) < pack(1, 3, 1) && pack(1, 3, 1) < pack(2, -100, 1));
    }

    #[test]
    fn finds_the_right_allele_and_diagonal() {
        let mut set = GermlineSet::new(SegmentKind::V);
        set.push("X1*01", b"ACGTTGCAAGCTTGACCGTAGGCTAACGTTAGCCATTGCA", None)
            .unwrap();
        set.push("X2*01", b"TTGACCATGCAGTCAGGATCCAGTACGATCAGGTCAACCT", None)
            .unwrap();
        let retriever = SeedVoteRetriever::build(
            &set,
            SeedShape::contiguous(9),
            Some(SeedShape::contiguous(7)),
            SeedVoteParams::default(),
            AmbiguityBudget::default(),
        )
        .unwrap();
        let query = encode(b"GGGGGTTGACCATGCAGTCAGGATCCAGTACGATCAGG");
        let mut scratch = RetrievalScratch::default();
        let mut out = Vec::new();
        retriever.retrieve(&query, &mut scratch, &mut out);
        assert_eq!(out[0].allele, 1);
        assert!(
            out[0].band.contains(10, 5),
            "diagonal +5 must be inside the band"
        );
    }
}
