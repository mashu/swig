//! Seed index over a germline set: CSR layout, direct addressing, packed postings.
//!
//! **Trick — direct-address CSR ("k-mer array").** Because a seed key is an
//! integer (see [`crate::seed`]), the index for weights up to 12 is simply an
//! offset array of length `4^w + 1` into one flat posting array — a lookup is
//! two loads, no hashing, no pointer chasing, and all postings of a key are
//! contiguous. For larger weights the table would be mostly empty, so keys are
//! stored sorted and found by binary search (same CSR postings). This is the
//! k-mer analogue of a suffix-array interval: every occurrence of a key is one
//! contiguous slice.
//!
//! **Trick — packed postings.** A posting `(allele, position)` is one `u32`
//! (`allele << 12 | position`), halving memory traffic versus a pair of `u32`s.
//! Germline segments are < 4096 nt, so 12 position bits suffice.
//!
//! Why not an FM-index / suffix array here: germline databases are a few hundred
//! kilobases and seeds are fixed-length, so a flat k-mer array is both smaller
//! in constant factors and faster than backward search. Variable-length exact
//! matching (maximal exact matches, matching statistics) is where a suffix
//! array or suffix automaton would pay off — see the README.

use std::fmt;

use crate::germline::GermlineSet;
use crate::seed::{SeedSampling, SeedShape};

const POSITION_BITS: u32 = 12;
/// Largest germline position representable in a posting.
pub const MAX_POSITION: u32 = (1 << POSITION_BITS) - 1;
/// Largest allele id representable in a posting.
pub const MAX_ALLELE: u32 = (1 << (32 - POSITION_BITS)) - 1;
/// Weights up to this many bits use a direct-address offset table.
const DIRECT_KEY_BITS: u32 = 24;

/// One seed occurrence: allele and 0-based start position, packed in a `u32`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Posting(u32);

impl Posting {
    #[inline]
    fn new(id: u32, position: u32) -> Self {
        Posting((id << POSITION_BITS) | position)
    }

    /// Id of the indexed sequence (an allele, or a gene frame).
    #[inline]
    pub fn id(self) -> u32 {
        self.0 >> POSITION_BITS
    }

    #[inline]
    pub fn position(self) -> u32 {
        self.0 & MAX_POSITION
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IndexError {
    SequenceTooLong { allele: String, len: usize },
    TooManyAlleles(usize),
}

impl fmt::Display for IndexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SequenceTooLong { allele, len } => write!(
                f,
                "allele '{allele}' is {len} nt; seed postings support at most {} nt",
                MAX_POSITION + 1
            ),
            Self::TooManyAlleles(n) => write!(
                f,
                "{n} alleles exceed the posting limit of {}",
                MAX_ALLELE + 1
            ),
        }
    }
}

impl std::error::Error for IndexError {}

#[derive(Clone, Debug)]
enum Lookup {
    /// `offsets[key] .. offsets[key + 1]`, length `4^weight + 1`.
    Direct { offsets: Vec<u32> },
    /// Sorted distinct keys with parallel offsets (one extra sentinel offset).
    Sorted { keys: Vec<u64>, offsets: Vec<u32> },
}

#[derive(Clone, Debug)]
pub struct KmerIndex {
    shape: SeedShape,
    lookup: Lookup,
    postings: Vec<Posting>,
}

/// One sequence to index: its k-mers are posted under `id` at `offset + position`.
#[derive(Clone, Copy, Debug)]
pub struct IndexedSequence<'a> {
    pub id: u32,
    pub name: &'a str,
    pub codes: &'a [u8],
    pub offset: u32,
}

impl KmerIndex {
    /// Index every allele of `set` under its own id.
    pub fn build(set: &GermlineSet, shape: SeedShape) -> Result<Self, IndexError> {
        let sequences = set
            .alleles()
            .iter()
            .enumerate()
            .map(|(id, allele)| IndexedSequence {
                id: id as u32,
                name: &allele.name,
                codes: &allele.codes,
                offset: 0,
            });
        Self::build_from(sequences, shape)
    }

    /// Index `sequences`; identical postings (same key, id and position) are
    /// stored once, so sequences sharing a frame share their common seeds.
    pub fn build_from<'a>(
        sequences: impl IntoIterator<Item = IndexedSequence<'a>>,
        shape: SeedShape,
    ) -> Result<Self, IndexError> {
        Self::build_sampled(sequences, shape, SeedSampling::All)
    }

    /// As [`KmerIndex::build_from`], indexing only the seeds `sampling` keeps.
    pub fn build_sampled<'a>(
        sequences: impl IntoIterator<Item = IndexedSequence<'a>>,
        shape: SeedShape,
        sampling: SeedSampling,
    ) -> Result<Self, IndexError> {
        let mut pairs: Vec<(u64, Posting)> = Vec::new();
        let (mut seeds, mut ranks) = (Vec::new(), Vec::new());
        for sequence in sequences {
            if sequence.id > MAX_ALLELE {
                return Err(IndexError::TooManyAlleles(sequence.id as usize + 1));
            }
            let end = sequence.offset as usize + sequence.codes.len();
            if end > MAX_POSITION as usize + 1 {
                return Err(IndexError::SequenceTooLong {
                    allele: sequence.name.to_string(),
                    len: end,
                });
            }
            sampling.seeds(&shape, sequence.codes, &mut seeds, &mut ranks);
            pairs.extend(
                seeds.iter().map(|&(pos, key)| {
                    (key, Posting::new(sequence.id, sequence.offset + pos as u32))
                }),
            );
        }
        pairs.sort_unstable();
        pairs.dedup();
        let postings: Vec<Posting> = pairs.iter().map(|&(_, p)| p).collect();

        let lookup = if shape.key_bits() <= DIRECT_KEY_BITS {
            let table_len = 1usize << shape.key_bits();
            let mut offsets = vec![0u32; table_len + 1];
            for &(key, _) in &pairs {
                offsets[key as usize + 1] += 1;
            }
            for k in 0..table_len {
                offsets[k + 1] += offsets[k];
            }
            Lookup::Direct { offsets }
        } else {
            let mut keys = Vec::new();
            let mut offsets = Vec::new();
            for (i, &(key, _)) in pairs.iter().enumerate() {
                if keys.last() != Some(&key) {
                    keys.push(key);
                    offsets.push(i as u32);
                }
            }
            offsets.push(pairs.len() as u32);
            Lookup::Sorted { keys, offsets }
        };
        Ok(Self {
            shape,
            lookup,
            postings,
        })
    }

    pub fn shape(&self) -> &SeedShape {
        &self.shape
    }

    /// All occurrences of `key`, sorted by allele then position.
    #[inline]
    pub fn postings(&self, key: u64) -> &[Posting] {
        let (start, end) = match &self.lookup {
            Lookup::Direct { offsets } => (offsets[key as usize], offsets[key as usize + 1]),
            Lookup::Sorted { keys, offsets } => match keys.binary_search(&key) {
                Ok(i) => (offsets[i], offsets[i + 1]),
                Err(_) => return &[],
            },
        };
        &self.postings[start as usize..end as usize]
    }

    /// Total number of postings (index size).
    pub fn len(&self) -> usize {
        self.postings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.postings.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::germline::SegmentKind;

    fn set() -> GermlineSet {
        let mut set = GermlineSet::new(SegmentKind::V);
        set.push("A*01", b"ACGTACGTTT", None).unwrap();
        set.push("B*01", b"TTACGTACAA", None).unwrap();
        set
    }

    #[test]
    fn direct_and_sorted_lookups_agree() {
        let mut set = set();
        set.push("C*01", b"GGATCCATTGACCAGTAAGC", None).unwrap();
        let small = KmerIndex::build(&set, SeedShape::contiguous(4)).unwrap();
        let key = SeedShape::contiguous(4)
            .kmers(&encode(b"ACGT"))
            .next()
            .unwrap()
            .1;
        let hits: Vec<_> = small
            .postings(key)
            .iter()
            .map(|p| (p.id(), p.position()))
            .collect();
        assert_eq!(hits, vec![(0, 0), (0, 4), (1, 2)]);

        let shape13 = SeedShape::contiguous(13);
        let large = KmerIndex::build(&set, shape13.clone()).unwrap();
        let key13 = shape13.kmers(&encode(b"ATCCATTGACCAG")).next().unwrap().1;
        let hits: Vec<_> = large
            .postings(key13)
            .iter()
            .map(|p| (p.id(), p.position()))
            .collect();
        assert_eq!(hits, vec![(2, 2)]);
        assert!(large.postings(key13 ^ 1).is_empty());
    }
}
