//! Seed voting with diagonal clustering (V retrieval).
//!
//! Every query seed that hits allele `a` at germline position `p` from query
//! position `i` votes for the diagonal `d = i − p` of `a`. Correct alleles
//! collect many votes on (nearly) one diagonal; the cluster's diagonal range,
//! widened by `band_pad`, becomes the DP band.
//!
//! **Gene frames.** Seeds are indexed per gene frame ([`GeneFrames`]), not
//! per allele, so a seed shared by a gene's alleles is one posting and one
//! vote. Retrieval ranks genes; every allele of a kept gene becomes a
//! candidate, with the frame's band shifted by the allele's offset.
//!
//! **Trick — radix-sorted packed votes.** Instead of a hash map
//! `allele → diagonal → votes`, each vote is one `u64`
//! `(frame · 2^b + d + max_len) << 8 | weight`, with `b` just wide enough
//! for the diagonals of this query. The sort key is therefore ~21 bits for an
//! IGH-sized database, and two LSD radix passes of 11 bits order the votes
//! by frame and diagonal in linear time; a single scan then groups them into
//! diagonal clusters. No allocation per query once the scratch has grown.
//!
//! **One vote per seed and frame.** A seed occurring twice in one gene frame
//! cannot tell which diagonal it supports, so it does not vote for that
//! gene. Without this rule a period-3 repeat (IGHV3-48's `AGTAGTAGT`)
//! collects votes quadratic in its length.
//!
//! **Trick — frequency cutoff.** Seeds with more than `max_occurrences`
//! postings (low-complexity or family-wide k-mers) carry little information
//! and dominate work; they are skipped, as in the method note.
//!
//! **Two-tier seeds.** Primary seeds vote with weight 4. Only when the best
//! primary cluster is weak are fallback (shorter) seeds added with weight 1.

use super::budget::AmbiguityBudget;
use super::frames::GeneFrames;
use super::{Candidate, CandidateBudget, ClusterVote, RetrievalScratch, Retriever};
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
    /// The band is the cluster's diagonal range widened by `band_pad` on each
    /// side, at most `2 · band_max + 1` diagonals wide.
    pub band_pad: i32,
    pub band_max: i32,
    /// Ignore seed geometry and align unbanded (`w < 0` in the note).
    pub unbanded: bool,
    /// Minimum weighted votes for an allele to be a candidate.
    pub min_votes: u32,
    /// Safety pool: never return more alleles than this (all alleles of the
    /// best genes first).
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
            band_max: 256,
            unbanded: false,
            min_votes: 1,
            max_alleles: 64,
        }
    }
}

/// Bit layout of a packed vote for one query (see the module documentation).
#[derive(Clone, Copy, Debug)]
struct VoteLayout {
    /// Added to a diagonal so it is non-negative (the longest allele).
    bias: i64,
    /// Bits of the biased diagonal.
    diagonal_bits: u32,
}

impl VoteLayout {
    const WEIGHT_BITS: u32 = 8;

    fn new(query_len: usize, max_allele_len: usize) -> Self {
        let span = (query_len + max_allele_len + 1) as u64;
        Self {
            bias: max_allele_len as i64,
            diagonal_bits: 64 - span.leading_zeros(),
        }
    }

    #[inline]
    fn pack(&self, allele: u32, diagonal: i64, weight: u8) -> u64 {
        let field = (u64::from(allele) << self.diagonal_bits) | (diagonal + self.bias) as u64;
        (field << Self::WEIGHT_BITS) | u64::from(weight)
    }

    #[inline]
    fn unpack(&self, key: u64) -> (u32, i32, u32) {
        let field = key >> Self::WEIGHT_BITS;
        let allele = (field >> self.diagonal_bits) as u32;
        let diagonal = (field & ((1 << self.diagonal_bits) - 1)) as i64 - self.bias;
        (allele, diagonal as i32, (key & 0xFF) as u32)
    }

    /// Sort-key bits above the weight.
    fn key_bits(&self, alleles: usize) -> u32 {
        self.diagonal_bits + (64 - (alleles as u64).leading_zeros())
    }
}

/// Stable LSD radix sort of `keys` on bits `shift .. shift + bits`.
fn radix_sort(keys: &mut Vec<u64>, buffer: &mut Vec<u64>, shift: u32, bits: u32) {
    const DIGIT: u32 = 11;
    if keys.len() < 512 {
        keys.sort_unstable_by_key(|&k| (k >> shift) & ((1 << bits) - 1));
        return;
    }
    let mut counts = [0u32; 1 << DIGIT];
    buffer.clear();
    buffer.resize(keys.len(), 0);
    let mut low = shift;
    while low < shift + bits {
        let width = DIGIT.min(shift + bits - low);
        let mask = (1u64 << width) - 1;
        let counts = &mut counts[..1 << width];
        counts.fill(0);
        for &key in keys.iter() {
            counts[((key >> low) & mask) as usize] += 1;
        }
        let mut sum = 0;
        for count in counts.iter_mut() {
            (*count, sum) = (sum, sum + *count);
        }
        for &key in keys.iter() {
            let digit = ((key >> low) & mask) as usize;
            buffer[counts[digit] as usize] = key;
            counts[digit] += 1;
        }
        std::mem::swap(keys, buffer);
        low += width;
    }
}

pub struct SeedVoteRetriever<B = AmbiguityBudget> {
    primary: KmerIndex,
    fallback: Option<KmerIndex>,
    frames: GeneFrames,
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
        let frames = GeneFrames::build(set, &primary);
        Ok(Self {
            primary: KmerIndex::build_from(frames.sequences(set), primary)?,
            fallback: fallback
                .map(|shape| KmerIndex::build_from(frames.sequences(set), shape))
                .transpose()?,
            frames,
            gene_count: set.gene_count(),
            params,
            budget,
        })
    }

    pub fn params(&self) -> &SeedVoteParams {
        &self.params
    }

    fn cast_votes(
        &self,
        index: &KmerIndex,
        weight: u8,
        query: &[u8],
        layout: &VoteLayout,
        keys: &mut Vec<u64>,
    ) {
        for (qpos, key) in index.shape().kmers(query) {
            let postings = index.postings(key);
            if postings.len() > self.params.max_occurrences {
                continue;
            }
            // Postings are sorted by frame: a repeat within one frame is a run.
            for (i, p) in postings.iter().enumerate() {
                let repeated = (i > 0 && postings[i - 1].id() == p.id())
                    || postings.get(i + 1).is_some_and(|n| n.id() == p.id());
                if !repeated {
                    keys.push(layout.pack(p.id(), qpos as i64 - i64::from(p.position()), weight));
                }
            }
        }
    }

    /// Sort votes and keep each frame's best diagonal cluster. Returns the best cluster's votes.
    fn cluster(
        &self,
        layout: &VoteLayout,
        keys: &mut Vec<u64>,
        buffer: &mut Vec<u64>,
        out: &mut Vec<ClusterVote>,
    ) -> u32 {
        let bits = layout.key_bits(self.frames.len());
        radix_sort(keys, buffer, VoteLayout::WEIGHT_BITS, bits);
        out.clear();
        let bin = |d: i32| d.div_euclid(self.params.bin_width);
        let max_step = 1 + self.params.cluster_gap_bins;
        let mut strongest = 0;
        let mut i = 0;
        while i < keys.len() {
            let (target, d0, w0) = layout.unpack(keys[i]);
            let mut best = ClusterVote {
                target,
                votes: 0,
                d_lo: d0,
                d_hi: d0,
            };
            let mut current = ClusterVote {
                target,
                votes: w0,
                d_lo: d0,
                d_hi: d0,
            };
            i += 1;
            while i < keys.len() {
                let (a, d, w) = layout.unpack(keys[i]);
                if a != target {
                    break;
                }
                if bin(d) - bin(current.d_hi) <= max_step {
                    current.votes += w;
                    current.d_hi = d;
                } else {
                    if current.votes > best.votes {
                        best = current;
                    }
                    current = ClusterVote {
                        target,
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

    /// Band of `allele` from its frame's cluster: allele diagonals are frame
    /// diagonals plus the allele's offset in the frame.
    fn band_for(&self, vote: &ClusterVote, allele: u32) -> Band {
        if self.params.unbanded {
            return Band::Full;
        }
        let shift = self.frames.offset(allele) as i32;
        let (lo, hi) = (
            vote.d_lo + shift - self.params.band_pad,
            vote.d_hi + shift + self.params.band_pad,
        );
        if hi - lo > 2 * self.params.band_max {
            Band::around((lo + hi).div_euclid(2), self.params.band_max)
        } else {
            Band::Diagonals { lo, hi }
        }
    }
}

impl<B: CandidateBudget> Retriever for SeedVoteRetriever<B> {
    fn retrieve(&self, query: &[u8], scratch: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
        out.clear();
        let RetrievalScratch {
            keys,
            buffer,
            votes,
            gene_votes,
            gene_order,
            ranked,
            ..
        } = scratch;
        keys.clear();
        let layout = VoteLayout::new(query.len(), self.frames.max_len());
        self.cast_votes(
            &self.primary,
            self.params.primary_weight,
            query,
            &layout,
            keys,
        );
        let mut strongest = self.cluster(&layout, keys, buffer, votes);
        if let (true, Some(fallback)) = (strongest < self.params.weak_votes, &self.fallback) {
            self.cast_votes(fallback, self.params.fallback_weight, query, &layout, keys);
            strongest = self.cluster(&layout, keys, buffer, votes);
        }
        if strongest < self.params.min_votes {
            return;
        }

        // Collapse per gene (best frame's votes), rank genes, apply the budget.
        gene_votes.clear();
        gene_votes.resize(self.gene_count, 0);
        for v in votes.iter() {
            let g = self.frames.gene(v.target) as usize;
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

        // Every allele of a kept frame, strongest frames first.
        votes.retain(|v| {
            v.votes >= self.params.min_votes
                && gene_votes[self.frames.gene(v.target) as usize] == KEPT
        });
        votes.sort_by(|a, b| b.votes.cmp(&a.votes).then(a.target.cmp(&b.target)));
        for vote in votes.iter() {
            for &allele in self.frames.members(vote.target) {
                if out.len() == self.params.max_alleles {
                    return;
                }
                out.push(Candidate {
                    allele,
                    band: self.band_for(vote, allele),
                    votes: vote.votes,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::germline::SegmentKind;

    #[test]
    fn packing_round_trips_negative_diagonals() {
        let layout = VoteLayout::new(600, 400);
        for &(a, d, w) in &[(0u32, -400i64, 4u8), (12345, 0, 1), (1 << 19, 600, 255)] {
            let (a2, d2, w2) = layout.unpack(layout.pack(a, d, w));
            assert_eq!((a2, i64::from(d2), w2), (a, d, u32::from(w)));
        }
        let key = |a, d| layout.pack(a, d, 1);
        assert!(key(1, -5) < key(1, 3) && key(1, 3) < key(2, -100));
    }

    #[test]
    fn radix_sort_orders_by_the_key_bits_only() {
        let mut state = 7u64;
        let mut keys: Vec<u64> = (0..5000)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                state >> 20
            })
            .collect();
        let mut expected = keys.clone();
        expected.sort_by_key(|&k| (k >> 8) & ((1 << 30) - 1));
        radix_sort(&mut keys, &mut Vec::new(), 8, 30);
        let field = |k: &u64| (k >> 8) & ((1 << 30) - 1);
        assert_eq!(
            keys.iter().map(field).collect::<Vec<_>>(),
            expected.iter().map(field).collect::<Vec<_>>()
        );
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
