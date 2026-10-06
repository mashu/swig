//! Gene frames: alleles of one gene in one shared coordinate system.
//!
//! Alleles of a gene differ by a few substitutions and by 5′/3′ length, so
//! almost every seed of one allele occurs in its siblings at a fixed offset.
//! Indexing every allele separately repeats each shared seed once per
//! allele — for IGHV3-30 and its ~20 alleles, twenty postings and twenty
//! votes per read seed — although retrieval only decides *which genes* (and
//! diagonals) deserve alignment; alleles are separated by the DP.
//!
//! A frame maps each member allele into the coordinates of the gene's
//! longest allele: `frame position = allele position + offset`. An allele
//! joins a frame only when its shared seeds are *position-consistent* with
//! the representative: at least [`MIN_SHARED`] of its seeds occur in the
//! representative, and at least [`MIN_CONSISTENT`] of those at one offset.
//! Substitutions only remove seeds, so SNP alleles and 5′/3′-truncated
//! alleles join; an allele with an internal indel would place half its seeds
//! at a second offset, so it gets a frame of its own, as does any divergent
//! entry. No seed is ever forced into a wrong coordinate system. Indexing frames instead of alleles
//! (with identical postings merged) keeps every allele-specific seed but
//! stores and votes each shared seed once per frame.

use std::collections::HashMap;

use crate::germline::{AlleleId, GeneId, GermlineSet};
use crate::index::IndexedSequence;
use crate::seed::SeedShape;

/// Fraction of an allele's seeds that must occur in the representative.
const MIN_SHARED: f64 = 0.3;
/// Fraction of the shared seeds that must agree on one offset.
const MIN_CONSISTENT: f64 = 0.9;

/// Alleles grouped into gene frames, with each allele's offset.
#[derive(Clone, Debug)]
pub struct GeneFrames {
    /// `(frame, offset)` of each allele.
    allele_frame: Vec<(u32, u32)>,
    /// Members of each frame, by allele id.
    members: Vec<Vec<AlleleId>>,
    frame_gene: Vec<GeneId>,
    /// Frame coordinates span `0..frame_len`.
    frame_len: Vec<usize>,
}

impl GeneFrames {
    pub fn build(set: &GermlineSet, shape: &SeedShape) -> Self {
        let mut by_gene: Vec<Vec<AlleleId>> = vec![Vec::new(); set.gene_count()];
        for (id, allele) in set.alleles().iter().enumerate() {
            by_gene[allele.gene as usize].push(id as AlleleId);
        }
        let mut frames: Vec<(GeneId, Vec<(AlleleId, i64)>)> = Vec::new();
        for (gene, mut alleles) in by_gene.into_iter().enumerate() {
            // Longest allele first: it covers most of the gene's sequence.
            alleles.sort_by_key(|&a| (std::cmp::Reverse(set.allele(a).codes.len()), a));
            let first = frames.len();
            for allele in alleles {
                let codes = &set.allele(allele).codes;
                let joined = (first..frames.len()).find_map(|f| {
                    let representative = frames[f].1[0].0;
                    let offset = shared_offset(shape, codes, &set.allele(representative).codes)?;
                    Some((f, offset + frames[f].1[0].1))
                });
                match joined {
                    Some((f, offset)) => frames[f].1.push((allele, offset)),
                    None => frames.push((gene as GeneId, vec![(allele, 0)])),
                }
            }
        }

        let mut allele_frame = vec![(0, 0); set.len()];
        let mut members = Vec::with_capacity(frames.len());
        let mut frame_gene = Vec::with_capacity(frames.len());
        let mut frame_len = Vec::with_capacity(frames.len());
        for (index, (gene, entries)) in frames.into_iter().enumerate() {
            let base = entries.iter().map(|&(_, o)| o).min().unwrap_or(0);
            let mut len = 0;
            for &(allele, offset) in &entries {
                let offset = (offset - base) as u32;
                allele_frame[allele as usize] = (index as u32, offset);
                len = len.max(offset as usize + set.allele(allele).codes.len());
            }
            let mut ids: Vec<AlleleId> = entries.iter().map(|&(a, _)| a).collect();
            ids.sort_unstable();
            members.push(ids);
            frame_gene.push(gene);
            frame_len.push(len);
        }
        Self {
            allele_frame,
            members,
            frame_gene,
            frame_len,
        }
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    pub fn gene(&self, frame: u32) -> GeneId {
        self.frame_gene[frame as usize]
    }

    /// Member alleles of `frame`, by id.
    pub fn members(&self, frame: u32) -> &[AlleleId] {
        &self.members[frame as usize]
    }

    /// `frame position = allele position + offset`.
    pub fn offset(&self, allele: AlleleId) -> u32 {
        self.allele_frame[allele as usize].1
    }

    pub fn max_len(&self) -> usize {
        self.frame_len.iter().copied().max().unwrap_or(0)
    }

    /// Every allele positioned in its frame, for [`crate::index::KmerIndex::build_from`].
    pub fn sequences<'a>(
        &'a self,
        set: &'a GermlineSet,
    ) -> impl Iterator<Item = IndexedSequence<'a>> {
        set.alleles()
            .iter()
            .zip(&self.allele_frame)
            .map(|(allele, &(frame, offset))| IndexedSequence {
                id: frame,
                name: &allele.name,
                codes: &allele.codes,
                offset,
            })
    }
}

/// Offset `o` with `representative position = allele position + o`, when
/// the seeds `allele` shares with `representative` are position-consistent.
fn shared_offset(shape: &SeedShape, allele: &[u8], representative: &[u8]) -> Option<i64> {
    let mut positions: HashMap<u64, Vec<usize>> = HashMap::new();
    for (pos, key) in shape.kmers(representative) {
        positions.entry(key).or_default().push(pos);
    }
    let mut votes: HashMap<i64, usize> = HashMap::new();
    let (mut seeds, mut shared) = (0usize, 0usize);
    for (pos, key) in shape.kmers(allele) {
        seeds += 1;
        let hits = positions.get(&key).map_or(&[][..], Vec::as_slice);
        shared += usize::from(!hits.is_empty());
        for &rep_pos in hits {
            *votes.entry(rep_pos as i64 - pos as i64).or_default() += 1;
        }
    }
    let (&offset, &count) = votes
        .iter()
        .max_by_key(|&(&offset, &count)| (count, std::cmp::Reverse(offset.abs())))?;
    let consistent = shared as f64 >= MIN_SHARED * seeds as f64
        && count as f64 >= MIN_CONSISTENT * shared as f64;
    consistent.then_some(offset)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::germline::SegmentKind;

    #[test]
    fn partial_and_snp_alleles_share_their_gene_frame() {
        let full = b"CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTC";
        let mut snp = full.to_vec();
        snp[30] = b'T';
        let partial = &full[7..];
        let mut set = GermlineSet::new(SegmentKind::V);
        set.push("IGHV1-2*01", full, None).unwrap();
        set.push("IGHV1-2*02", &snp, None).unwrap();
        set.push("IGHV1-2*03", partial, None).unwrap();
        // An internal 3-nt deletion: seeds after it sit at a second offset.
        let mut deletion = full.to_vec();
        deletion.drain(30..33);
        set.push("IGHV1-2*04", &deletion, None).unwrap();
        set.push(
            "IGHV9-9*01",
            b"TTTTTTTTTTGGGGGGGGGGAAAAAAAAAACCCCCCCCCC",
            None,
        )
        .unwrap();
        let frames = GeneFrames::build(&set, &SeedShape::contiguous(9));
        assert_eq!(frames.len(), 3);
        assert_eq!(frames.members(0), &[0, 1, 2]);
        assert_eq!(frames.members(1), &[3]);
        assert_eq!(
            (frames.offset(0), frames.offset(1), frames.offset(2)),
            (0, 0, 7)
        );
        assert_eq!(frames.max_len(), full.len());
    }
}
