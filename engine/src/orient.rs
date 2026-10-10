//! Query orientation from germline k-mer support on both strands.
//!
//! A read is in germline sense when many of its k-mers occur in the V and J
//! germline. The statistic is the number of read positions whose k-mer
//! occurs anywhere in the germline, counted for the read and for its
//! reverse complement; the larger count wins and ties stay forward.
//!
//! Each read position counts at most once per strand, so repeats cannot
//! inflate the evidence. (Seed *votes* count postings, and a period-3 repeat
//! such as IGHJ6's `TACTACTAC` matches the antisense of IGHV3-48's
//! `AGTAGTAGT` CDR2 many times over.) With `k = 11` a random k-mer occurs in
//! an IGH-sized germline with probability ~1 %, while a germline-sense read
//! at 15 % SHM still keeps ~17 % of its V/J k-mers, so the two strands
//! separate by an order of magnitude.
//!
//! The count is only trusted when it is significant: under the null (a
//! strand unrelated to the germline) it is Poisson with mean
//! `λ₀ = k-mers × density`, the density being the fraction of all `4^k`
//! k-mers present in the germline (known exactly from the bitset). When
//! neither strand exceeds `λ₀ + 3√λ₀` — a read at 25–35 % SHM keeps almost no
//! exact 11-mer — [`StrandSupport::decisive`] is false and the caller falls
//! back to a more sensitive, more expensive test (seed votes with the short
//! fallback seed on both strands).
//!
//! **Trick — a presence bitset over all `4^k` k-mers** (512 KiB for `k = 11`):
//! one pass over the read computes the forward and reverse-complement keys
//! together, two bit tests per base, no posting lists. Only the chosen
//! strand is then searched.

use crate::alphabet::N;
use crate::germline::GermlineSet;

/// Strand of the query that the call refers to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    Reverse,
}

impl Orientation {
    pub fn is_reverse(self) -> bool {
        self == Orientation::Reverse
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrientationPolicy {
    /// Consider the reverse complement at all.
    pub try_reverse: bool,
}

impl Default for OrientationPolicy {
    fn default() -> Self {
        Self { try_reverse: true }
    }
}

impl OrientationPolicy {
    /// Forward only (pre-oriented input).
    pub const FORWARD_ONLY: Self = Self { try_reverse: false };

    /// Strictly stronger reverse evidence flips the strand; ties stay forward.
    pub fn choose(&self, forward: u32, reverse: u32) -> Orientation {
        if self.try_reverse && reverse > forward {
            Orientation::Reverse
        } else {
            Orientation::Forward
        }
    }
}

/// Read positions whose k-mer occurs in the germline, per strand.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StrandSupport {
    pub forward: u32,
    pub reverse: u32,
    /// Expected count on a strand unrelated to the germline.
    pub null: f64,
}

impl StrandSupport {
    /// The stronger strand is significantly above chance (`> λ₀ + 3√λ₀`).
    pub fn decisive(&self) -> bool {
        f64::from(self.forward.max(self.reverse)) > self.null + 3.0 * self.null.sqrt()
    }
}

/// Presence bitset of every germline k-mer.
#[derive(Clone, Debug)]
pub struct StrandIndex {
    k: usize,
    bits: Vec<u64>,
    /// Fraction of all `4^k` k-mers present.
    density: f64,
}

impl StrandIndex {
    /// Default k-mer length (see the module documentation).
    pub const K: usize = 11;

    /// Panics unless `1 <= k <= 14` (the bitset has `4^k` bits).
    pub fn build<'a>(sets: impl IntoIterator<Item = &'a GermlineSet>, k: usize) -> Self {
        assert!(
            (1..=14).contains(&k),
            "strand k-mer length must be in 1..=14"
        );
        let mut bits = vec![0u64; (1usize << (2 * k)).div_ceil(64)];
        for set in sets {
            for allele in set.alleles() {
                for_each_kmer(&allele.codes, k, |key, _| {
                    bits[(key >> 6) as usize] |= 1 << (key & 63);
                });
            }
        }
        let present: u64 = bits.iter().map(|w| u64::from(w.count_ones())).sum();
        let density = present as f64 / (1u64 << (2 * k)) as f64;
        Self { k, bits, density }
    }

    #[inline]
    fn contains(&self, key: u64) -> bool {
        self.bits[(key >> 6) as usize] >> (key & 63) & 1 == 1
    }

    pub fn support(&self, query: &[u8]) -> StrandSupport {
        let mut support = StrandSupport::default();
        let mut kmers = 0u32;
        for_each_kmer(query, self.k, |forward, reverse| {
            kmers += 1;
            support.forward += u32::from(self.contains(forward));
            support.reverse += u32::from(self.contains(reverse));
        });
        support.null = f64::from(kmers) * self.density;
        support
    }
}

/// Call `f(key, reverse_complement_key)` for every k-mer without `N`.
#[inline]
fn for_each_kmer(codes: &[u8], k: usize, mut f: impl FnMut(u64, u64)) {
    let mask = (1u64 << (2 * k)) - 1;
    let top = 2 * (k as u32 - 1);
    let (mut forward, mut reverse, mut valid) = (0u64, 0u64, 0usize);
    for &code in codes {
        if code >= N {
            valid = 0;
            continue;
        }
        forward = ((forward << 2) | u64::from(code)) & mask;
        reverse = (reverse >> 2) | (u64::from(3 - code) << top);
        valid += 1;
        if valid >= k {
            f(forward, reverse);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{encode, reverse_complement_into};
    use crate::germline::SegmentKind;

    #[test]
    fn reverse_keys_are_the_keys_of_the_reverse_complement() {
        let codes = encode(b"ACGTTGCAAGCTNGGACCA");
        let mut rc = Vec::new();
        reverse_complement_into(&codes, &mut rc);
        let (mut forward, mut reverse) = (Vec::new(), Vec::new());
        for_each_kmer(&codes, 4, |f, r| {
            forward.push(f);
            reverse.push(r);
        });
        let mut rc_forward = Vec::new();
        for_each_kmer(&rc, 4, |f, _| rc_forward.push(f));
        rc_forward.reverse();
        assert_eq!(reverse, rc_forward);
    }

    #[test]
    fn germline_sense_wins_and_repeats_count_once() {
        let mut v = GermlineSet::new(SegmentKind::V);
        v.push(
            "V1*01",
            b"CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCC",
            None,
        )
        .unwrap();
        v.push("V2*01", b"AGTAGTAGTAGTAGTAGTAGTAGTAGT", None)
            .unwrap();
        let index = StrandIndex::build([&v], 11);
        let read = encode(b"TTTTCAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGTTTT");
        let support = index.support(&read);
        assert!(support.forward > 20 && support.reverse == 0);
        assert!(support.decisive());
        let choose = |s: StrandSupport| OrientationPolicy::default().choose(s.forward, s.reverse);
        assert_eq!(choose(support), Orientation::Forward);
        let mut rc = Vec::new();
        reverse_complement_into(&read, &mut rc);
        assert_eq!(choose(index.support(&rc)), Orientation::Reverse);
        // ACTACT… is the antisense of the AGTAGT… repeat: each position counts once.
        let repeat = encode(b"ACTACTACTACTACTACT");
        assert_eq!(index.support(&repeat).reverse, 8);
    }
}
