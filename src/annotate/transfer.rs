//! Sibling transfer: carry positions from the most similar confidently
//! annotated allele through a local alignment.
//!
//! **Trick — inverted 8-mer index.** Distinct 8-mers are 16-bit keys; one CSR
//! table (65,536 offsets) lists the alleles containing each key. Finding the
//! nearest relative of an allele touches only the postings of its own
//! k-mers, i.e. only alleles that share something with it, instead of
//! comparing against every allele. The index is built lazily, only when some
//! allele needs a relative, and a relative must share at least half of the
//! target's k-mers (alleles of one gene share ~90 %, unrelated genes a few
//! per cent), so a chance match can never donate an annotation.

use crate::align::{Aligner, Alignment, BandedGotoh, GotohWorkspace, Op};
use crate::alphabet::N;
use crate::band::Band;
use crate::scoring::{AffineScoring, SubstitutionProfile};

const K: usize = 8;
/// Fewer shared 8-mers than this: no meaningful relative.
pub const MIN_SHARED_KMERS: usize = 12;
/// A relative must share at least this fraction of the target's distinct 8-mers.
pub const MIN_SHARED_FRACTION: f64 = 0.5;
const SCORING: AffineScoring = AffineScoring::new(2, -3, -5, -1);

/// Sorted distinct 8-mers of a sequence (an N resets the window).
pub fn kmers(codes: &[u8]) -> Vec<u16> {
    let mut out = Vec::with_capacity(codes.len());
    let (mut key, mut valid) = (0u16, 0usize);
    for &code in codes {
        if code >= N {
            valid = 0;
            continue;
        }
        key = (key << 2) | u16::from(code);
        valid += 1;
        if valid >= K {
            out.push(key);
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// Inverted index: alleles containing each distinct 8-mer.
#[derive(Debug, Default)]
pub struct RelativeIndex {
    offsets: Vec<u32>,
    postings: Vec<u32>,
    kmers: Vec<Vec<u16>>,
    counts: Vec<u32>,
    touched: Vec<u32>,
}

impl RelativeIndex {
    pub fn build<'a>(sequences: impl Iterator<Item = &'a [u8]>) -> Self {
        let kmers: Vec<Vec<u16>> = sequences.map(kmers).collect();
        let mut offsets = vec![0u32; (1 << 16) + 1];
        for key in kmers.iter().flatten() {
            offsets[usize::from(*key) + 1] += 1;
        }
        for i in 1..offsets.len() {
            offsets[i] += offsets[i - 1];
        }
        let mut fill = offsets.clone();
        let mut postings = vec![0u32; kmers.iter().map(Vec::len).sum()];
        for (allele, set) in kmers.iter().enumerate() {
            for &key in set {
                let slot = &mut fill[usize::from(key)];
                postings[*slot as usize] = allele as u32;
                *slot += 1;
            }
        }
        let counts = vec![0; kmers.len()];
        Self {
            offsets,
            postings,
            kmers,
            counts,
            touched: Vec::new(),
        }
    }

    /// Most similar `eligible` allele to `target` (most shared 8-mers, lowest
    /// index on ties), if it shares enough to count as a relative.
    pub fn nearest(&mut self, target: usize, eligible: &[bool]) -> Option<usize> {
        let Self {
            offsets,
            postings,
            kmers,
            counts,
            touched,
        } = self;
        for &key in &kmers[target] {
            let key = usize::from(key);
            for &allele in &postings[offsets[key] as usize..offsets[key + 1] as usize] {
                if counts[allele as usize] == 0 {
                    touched.push(allele);
                }
                counts[allele as usize] += 1;
            }
        }
        let needed = MIN_SHARED_KMERS
            .max((MIN_SHARED_FRACTION * kmers[target].len() as f64).ceil() as usize);
        let best = touched
            .iter()
            .map(|&a| (a as usize, counts[a as usize] as usize))
            .filter(|&(a, shared)| a != target && eligible[a] && shared >= needed)
            .max_by_key(|&(a, shared)| (shared, std::cmp::Reverse(a)))
            .map(|(a, _)| a);
        for &allele in touched.iter() {
            counts[allele as usize] = 0;
        }
        touched.clear();
        best
    }
}

/// Size of the intersection of two sorted sets.
pub fn shared(a: &[u16], b: &[u16]) -> usize {
    let (mut i, mut j, mut count) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
        }
    }
    count
}

/// Query position aligned to reference position `r`: the aligned base, the
/// next query base across a deletion, or a diagonal extrapolation beyond the
/// aligned part; `None` if that falls outside the query.
fn map_position(alignment: &Alignment, r: u32, query_len: usize) -> Option<u32> {
    let mapped = if r < alignment.r_start {
        i64::from(alignment.q_start) - i64::from(alignment.r_start - r)
    } else if r >= alignment.r_end {
        i64::from(alignment.q_end) + i64::from(r - alignment.r_end)
    } else {
        let (_, q, _) = alignment
            .columns()
            .find(|&(op, _, rc)| op != Op::Ins && rc == r)?;
        i64::from(q)
    };
    u32::try_from(mapped)
        .ok()
        .filter(|&q| (q as usize) < query_len)
}

/// Map annotated `positions` of `reference` onto `query`.
pub fn transfer(
    query: &[u8],
    reference: &[u8],
    positions: &[Option<u32>],
    workspace: &mut GotohWorkspace,
) -> Option<Vec<Option<u32>>> {
    let profile = SubstitutionProfile::build(&SCORING, reference);
    let alignment = BandedGotoh.align(query, &profile, SCORING.gaps, Band::Full, workspace)?;
    Some(
        positions
            .iter()
            .map(|p| p.and_then(|r| map_position(&alignment, r, query.len())))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;

    #[test]
    fn transfers_positions_across_substitutions_and_indels() {
        let reference = encode(b"GATTACACCGGTTAACCGGTACGTAGCTAGCTAGGATCCAGTC");
        // Same sequence with one substitution and a 2-nt deletion before position 30.
        let query = encode(
            b"GATTACACCGGTTTACCGGTACGTAGCTAGATAGGATCCAGTC"
                .as_slice()
                .iter()
                .enumerate()
                .filter(|&(i, _)| i != 28 && i != 29)
                .map(|(_, &b)| b)
                .collect::<Vec<_>>()
                .as_slice(),
        );
        let mapped = transfer(
            &query,
            &reference,
            &[Some(5), Some(35), None],
            &mut GotohWorkspace::default(),
        )
        .unwrap();
        assert_eq!(mapped, vec![Some(5), Some(33), None]);
        assert!(shared(&kmers(&reference), &kmers(&query)) >= MIN_SHARED_KMERS);
    }

    #[test]
    fn index_finds_the_true_sibling_and_rejects_chance_similarity() {
        let a = encode(b"GATTACACCGGTTAACCGGTACGTAGCTAGCTAGGATCCAGTCAAGGCTTACCGATGCA");
        let mut b = a.clone();
        b[20] = (b[20] + 1) % 4;
        let unrelated = encode(b"TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAAT");
        let mut index =
            RelativeIndex::build([a.as_slice(), unrelated.as_slice(), b.as_slice()].into_iter());
        assert_eq!(index.nearest(2, &[true, true, true]), Some(0));
        assert_eq!(index.nearest(1, &[true, true, true]), None);
        assert_eq!(index.nearest(2, &[false, true, true]), None);
    }
}
