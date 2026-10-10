//! Reach classes: the alleles an anchored alignment cannot tell apart.
//!
//! An anchored alignment of a splice window ([`crate::constant::SpliceWindow`])
//! starts at an allele's base 0 within `s` bases of the splice and reaches at
//! most `R = (n − x) + s` of its bases (`n` the read length, `x` the
//! splice): every column past `R` lies outside the band. Two alleles whose
//! *reached sequences* `a[..min(R, |a|)]` are equal therefore fill identical
//! in-band DP cells — the same score, the same traceback — so one
//! representative per class is aligned and the others join its placement
//! afterwards ([`ReachClasses::members`]). This is exact for the anchored
//! hypothesis set, and cheap where it matters: human IGHC alleles differ
//! mostly in CH2/CH3, far beyond a short read's reach (81 distinct
//! sequences, 7 distinct 32-nt prefixes).
//!
//! **Trick — sorted alleles and an LCP array.** Sorted by sequence, the
//! alleles sharing a reached sequence are one run, and consecutive alleles
//! `a`, `b` share it iff `lcp(a, b) ≥ min(R, max(|a|, |b|))` (equal prefixes
//! of length `R`, or equal sequences shorter than `R`). Classes are found per
//! read in one pass over the alleles, with no hashing.

use super::{Candidate, RetrievalScratch, Retriever};
use crate::band::Band;
use crate::germline::{AlleleId, GermlineSet};

/// Retriever of anchored splice windows: one candidate per reach class.
#[derive(Clone, Debug)]
pub struct ReachClasses {
    /// Allele ids sorted by sequence.
    sorted: Vec<AlleleId>,
    /// `lcp[k]`: common prefix length of `sorted[k − 1]` and `sorted[k]`.
    lcp: Vec<u32>,
    /// Sequence length of `sorted[k]`.
    lens: Vec<u32>,
    /// Position of each allele in `sorted`.
    rank: Vec<u32>,
    /// Splice spread `s`: the splice lies `s + 1` bases into a window.
    spread: u32,
}

impl ReachClasses {
    pub fn new(set: &GermlineSet, spread: u32) -> Self {
        let codes = |id: AlleleId| &set.allele(id).codes;
        let mut sorted: Vec<AlleleId> = (0..set.len() as AlleleId).collect();
        sorted.sort_by(|&a, &b| codes(a).cmp(codes(b)).then(a.cmp(&b)));
        let lcp = (0..sorted.len())
            .map(|k| match k {
                0 => 0,
                k => codes(sorted[k - 1])
                    .iter()
                    .zip(codes(sorted[k]))
                    .take_while(|(a, b)| a == b)
                    .count() as u32,
            })
            .collect();
        let lens = sorted.iter().map(|&a| codes(a).len() as u32).collect();
        let mut rank = vec![0; sorted.len()];
        for (k, &allele) in sorted.iter().enumerate() {
            rank[allele as usize] = k as u32;
        }
        Self {
            sorted,
            lcp,
            lens,
            rank,
            spread,
        }
    }

    /// Allele bases an anchored alignment of a `window_len`-base splice
    /// window can reach: `R = window_len − 1 = (n − x) + s`.
    pub fn reach(&self, window_len: usize) -> usize {
        window_len.saturating_sub(1)
    }

    /// Whether `sorted[k − 1]` and `sorted[k]` differ within `reach`.
    #[inline]
    fn splits(&self, k: usize, reach: usize) -> bool {
        let longer = self.lens[k - 1].max(self.lens[k]) as usize;
        (self.lcp[k] as usize) < reach.min(longer)
    }

    /// Every allele of `representative`'s class at `reach` (itself included),
    /// in sequence order.
    pub fn members(&self, representative: AlleleId, reach: usize) -> &[AlleleId] {
        let k = self.rank[representative as usize] as usize;
        let mut lo = k;
        while lo > 0 && !self.splits(lo, reach) {
            lo -= 1;
        }
        let mut hi = k + 1;
        while hi < self.sorted.len() && !self.splits(hi, reach) {
            hi += 1;
        }
        &self.sorted[lo..hi]
    }
}

impl Retriever for ReachClasses {
    /// One candidate per reach class of the splice window `query`: its
    /// lowest allele id, the class size as votes, and the band of the `2s + 1`
    /// start diagonals (anchors at window positions `0..=2s`). Candidates
    /// are in allele-id order.
    fn retrieve(&self, query: &[u8], _: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
        out.clear();
        let reach = self.reach(query.len());
        let band = Band::Diagonals {
            lo: 0,
            hi: 2 * self.spread as i32,
        };
        let mut start = 0;
        for k in 1..=self.sorted.len() {
            if k == self.sorted.len() || self.splits(k, reach) {
                let class = &self.sorted[start..k];
                out.push(Candidate {
                    allele: *class.iter().min().expect("a class is not empty"),
                    band,
                    votes: class.len() as u32,
                });
                start = k;
            }
        }
        out.sort_unstable_by_key(|candidate| candidate.allele);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::decode;
    use crate::germline::SegmentKind;
    use crate::testing::{random_codes, rng};

    /// Random databases with duplicates, alleles that are prefixes of one
    /// another, shared prefixes and mixed lengths.
    fn database(next: &mut impl FnMut() -> u64) -> GermlineSet {
        let mut set = GermlineSet::new(SegmentKind::C);
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        for k in 0..(1 + next() % 30) {
            let sequence = match (next() % 4, sequences.len()) {
                (_, 0) | (0, _) => {
                    let len = 1 + (next() % 24) as usize;
                    random_codes(next, len)
                }
                (kind, count) => {
                    let mut base = sequences[(next() as usize) % count].clone();
                    match kind {
                        1 => base.truncate(1 + (next() as usize) % base.len()),
                        2 => {
                            let p = (next() as usize) % base.len();
                            base[p] = (base[p] + 1) % 4;
                        }
                        _ => {
                            let tail = (next() % 6) as usize;
                            base.extend(random_codes(next, tail));
                        }
                    }
                    base
                }
            };
            set.push(
                &format!("IGHC{}*01", k % 5),
                decode(&sequence).as_bytes(),
                None,
            )
            .unwrap();
            sequences.push(sequence);
        }
        set
    }

    /// Classes, representatives and members equal brute-force grouping by
    /// the reached sequence `codes[..min(R, len)]`, for every reach.
    #[test]
    fn classes_group_alleles_by_their_reached_sequence() {
        let mut next = rng(17);
        let mut scratch = RetrievalScratch::default();
        let mut out = Vec::new();
        for _ in 0..300 {
            let set = database(&mut next);
            let longest = set.alleles().iter().map(|a| a.codes.len()).max().unwrap();
            let classes = ReachClasses::new(&set, 3);
            for reach in 0..=longest + 2 {
                let reached = |a: AlleleId| {
                    let codes = &set.allele(a).codes;
                    &codes[..reach.min(codes.len())]
                };
                classes.retrieve(&vec![0; reach + 1], &mut scratch, &mut out);
                let ids: Vec<AlleleId> = (0..set.len() as AlleleId).collect();
                let representatives: Vec<AlleleId> = ids
                    .iter()
                    .copied()
                    .filter(|&a| ids.iter().all(|&b| b >= a || reached(b) != reached(a)))
                    .collect();
                let found: Vec<AlleleId> = out.iter().map(|c| c.allele).collect();
                assert_eq!(found, representatives, "reach {reach}");
                for candidate in &out {
                    let mut members = classes.members(candidate.allele, reach).to_vec();
                    members.sort_unstable();
                    let expected: Vec<AlleleId> = ids
                        .iter()
                        .copied()
                        .filter(|&b| reached(b) == reached(candidate.allele))
                        .collect();
                    assert_eq!(members, expected, "reach {reach}");
                    assert_eq!(candidate.votes as usize, expected.len());
                    assert_eq!(candidate.band, Band::Diagonals { lo: 0, hi: 6 });
                }
            }
        }
    }
}
