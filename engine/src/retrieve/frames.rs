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
//! representative, and every one of those at one offset. Substitutions only
//! remove seeds, so SNP alleles and 5′/3′-truncated alleles join; an allele
//! with an indel, however near its end, places the seeds past it at a second
//! offset, so it gets a frame of its own, as does any divergent entry. No
//! seed is ever forced into a wrong coordinate system: a seed posted at two
//! positions of one frame could not tell which diagonal it supports, so it
//! would vote for none of the gene's alleles. Indexing frames instead of
//! alleles (with identical postings merged) keeps every allele-specific seed
//! but stores and votes each shared seed once per frame.

use std::collections::HashMap;

use crate::germline::{AlleleId, GeneId, GermlineSet};
use crate::index::IndexedSequence;
use crate::seed::SeedShape;

/// Fraction of an allele's seeds that must occur in the representative.
const MIN_SHARED: f64 = 0.3;

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

    /// Frame of `allele`.
    pub fn frame(&self, allele: AlleleId) -> u32 {
        self.allele_frame[allele as usize].0
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
///
/// A seed votes for each offset at which the representative holds its key,
/// at most once per offset, so an offset with as many votes as there are
/// shared seeds is one that every shared seed supports. Of several such
/// offsets (a repeat in the representative) the one nearest 0 wins, then
/// the smaller: hash order must not decide.
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
    let (&offset, &count) = votes.iter().max_by_key(|&(&offset, &count)| {
        (
            count,
            std::cmp::Reverse(offset.abs()),
            std::cmp::Reverse(offset),
        )
    })?;
    let consistent = shared as f64 >= MIN_SHARED * seeds as f64 && count == shared;
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

    /// Every seed a member shares with its frame's representative is posted
    /// where the representative holds it: no seed is forced into a second
    /// coordinate system of its frame.
    fn assert_shared_seeds_are_consistent(set: &GermlineSet, frames: &GeneFrames, k: usize) {
        let shape = SeedShape::contiguous(k);
        for frame in 0..frames.len() as u32 {
            let members = frames.members(frame);
            let mut posted: HashMap<u64, Vec<usize>> = HashMap::new();
            for &a in members {
                for (pos, key) in shape.kmers(&set.allele(a).codes) {
                    posted
                        .entry(key)
                        .or_default()
                        .push(pos + frames.offset(a) as usize);
                }
            }
            let representative = members
                .iter()
                .copied()
                .max_by_key(|&a| (set.allele(a).codes.len(), std::cmp::Reverse(a)))
                .expect("a frame has members");
            let mut held: HashMap<u64, Vec<usize>> = HashMap::new();
            for (pos, key) in shape.kmers(&set.allele(representative).codes) {
                held.entry(key)
                    .or_default()
                    .push(pos + frames.offset(representative) as usize);
            }
            for (key, at) in &posted {
                if let Some(rep) = held.get(key) {
                    assert!(
                        at.iter().all(|p| rep.contains(p)),
                        "frame {frame}: a shared seed posted at {at:?}, the representative holds it at {rep:?}"
                    );
                }
            }
        }
    }

    /// An indel near an allele's end leaves > 90 % of its shared seeds at
    /// one offset, but the seeds past it would be posted a second time at
    /// shifted frame positions, where the one-vote-per-seed-and-frame rule
    /// silences them for every allele of the gene.
    #[test]
    fn an_allele_with_an_indel_near_its_end_gets_its_own_frame() {
        let mut next = crate::testing::rng(5);
        let full = crate::testing::random_codes(&mut next, 300);
        let mut snp = full.clone();
        snp[150] = (snp[150] + 1) % 4;
        let mut deletion = full.clone();
        deletion.drain(270..273);
        let mut set = GermlineSet::new(SegmentKind::V);
        for (name, codes) in [
            ("IGHV1-2*01", &full),
            ("IGHV1-2*02", &snp),
            ("IGHV1-2*03", &deletion),
        ] {
            set.push(name, crate::alphabet::decode(codes).as_bytes(), None)
                .unwrap();
        }
        let frames = GeneFrames::build(&set, &SeedShape::contiguous(9));
        assert_eq!(frames.members(0), &[0, 1]);
        assert_eq!(frames.members(1), &[2]);
        assert_shared_seeds_are_consistent(&set, &frames, 9);
    }

    #[test]
    fn frames_of_the_synthetic_database_post_shared_seeds_consistently() {
        let db = crate::sim::SyntheticDatabase::generate(&Default::default(), 11);
        let frames = GeneFrames::build(&db.v, &SeedShape::contiguous(9));
        assert_shared_seeds_are_consistent(&db.v, &frames, 9);
    }

    /// Two offsets that every shared seed supports, equally far from 0:
    /// the choice must not depend on hash order.
    #[test]
    fn offset_ties_are_broken_by_the_offset() {
        let mut next = crate::testing::rng(8);
        let p = crate::testing::random_codes(&mut next, 30);
        let q = crate::testing::random_codes(&mut next, 30);
        let representative = [&p[..], &q[..], &p[..]].concat();
        let n = [crate::alphabet::N; 30];
        let allele = [&n[..], &p[..], &n[..]].concat();
        let shape = SeedShape::contiguous(9);
        let offsets: Vec<Option<i64>> = (0..64)
            .map(|_| shared_offset(&shape, &allele, &representative))
            .collect();
        assert!(offsets.iter().all(|&o| o == Some(-30)), "{offsets:?}");
    }
}
