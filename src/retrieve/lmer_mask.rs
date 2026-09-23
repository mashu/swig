//! Bit-parallel L-mer prefilter for short segments (D).
//!
//! D genes are 10–40 nt and the D window is a few dozen bases, so seed
//! *diagonals* carry little information; the question is only "which alleles
//! can possibly clear the admission gate?".
//!
//! **Trick — allele bitsets per L-mer.** For every L-mer (an integer in
//! `0..4^L`) the index stores a bitset over alleles containing it
//! (`⌈alleles/64⌉` machine words). Scanning the window ORs nothing and
//! branches little: for each window L-mer, iterate the set bits of its row
//! (`trailing_zeros` + `b &= b − 1`) and bump a per-allele counter.
//!
//! **Lossless by pigeonhole.** With `min_shared = 1` the filter never drops an
//! allele that could pass the admission gates of the method note:
//! * exact-run floor `k ≥ L` ⇒ the run itself contains an L-mer;
//! * distributed-D (≥ 16 aligned columns, ≥ 80 % identity, so ≤ 3 mismatches
//!   in an ungapped 16-column stretch) ⇒ 13 matches split into ≤ 4 runs, one of
//!   length ≥ ⌈13/4⌉ = 4 — so `L = 4` is lossless for that gate too.
//!
//! Survivors are aligned unbanded (the window is short, a full matrix is a few
//! hundred cells).

use super::{Candidate, RetrievalScratch, Retriever};
use crate::band::Band;
use crate::germline::{AlleleId, GermlineSet};
use crate::seed::SeedShape;

/// Longest supported L (table has `4^L` rows).
pub const MAX_LMER: usize = 10;

pub struct LmerMaskRetriever {
    shape: SeedShape,
    words: usize,
    masks: Vec<u64>,
    allele_count: usize,
    min_shared: u32,
}

impl LmerMaskRetriever {
    /// Panics unless `1 <= l <= MAX_LMER` (a configuration programming error).
    pub fn build(set: &GermlineSet, l: usize, min_shared: u32) -> Self {
        assert!((1..=MAX_LMER).contains(&l), "L-mer length must be in 1..={MAX_LMER}, got {l}");
        let shape = SeedShape::contiguous(l);
        let allele_count = set.len();
        let words = allele_count.div_ceil(64).max(1);
        let mut masks = vec![0u64; (1usize << (2 * l)) * words];
        for (id, allele) in set.alleles().iter().enumerate() {
            for (_, key) in shape.kmers(&allele.codes) {
                masks[key as usize * words + id / 64] |= 1u64 << (id % 64);
            }
        }
        Self { shape, words, masks, allele_count, min_shared }
    }
}

impl Retriever for LmerMaskRetriever {
    fn retrieve(&self, query: &[u8], scratch: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
        out.clear();
        let counts = &mut scratch.counts;
        counts.clear();
        counts.resize(self.allele_count, 0);
        for (_, key) in self.shape.kmers(query) {
            let row = &self.masks[key as usize * self.words..(key as usize + 1) * self.words];
            for (w, &bits) in row.iter().enumerate() {
                let mut b = bits;
                while b != 0 {
                    counts[w * 64 + b.trailing_zeros() as usize] += 1;
                    b &= b - 1;
                }
            }
        }
        out.extend(
            counts
                .iter()
                .enumerate()
                .filter(|&(_, &c)| c >= self.min_shared.max(1))
                .map(|(id, &c)| Candidate { allele: id as AlleleId, band: Band::Full, votes: c }),
        );
        out.sort_by(|a, b| b.votes.cmp(&a.votes).then(a.allele.cmp(&b.allele)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::germline::SegmentKind;

    #[test]
    fn counts_shared_lmers() {
        let mut set = GermlineSet::new(SegmentKind::D);
        set.push("D1*01", b"GTATTACGATATTTTGACTGG", None).unwrap();
        set.push("D2*01", b"AGGATATTGTAGTAGTACCAGCTGCTATGCC", None).unwrap();
        set.push("D3*01", b"CCCCCCCCCC", None).unwrap();
        let retriever = LmerMaskRetriever::build(&set, 4, 1);
        let mut scratch = RetrievalScratch::default();
        let mut out = Vec::new();
        retriever.retrieve(&encode(b"AAAACGATATTTTGAAAA"), &mut scratch, &mut out);
        assert_eq!(out[0].allele, 0);
        assert!(out.iter().all(|c| c.allele != 2));
    }
}
