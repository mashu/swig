//! Seed shapes and rolling 2-bit k-mer extraction.
//!
//! **Trick — k-mers as integers.** With 2 bits per base a k-mer of span ≤ 31 is
//! a `u64`. Sliding the window costs one shift, one OR and one AND per base, and
//! the resulting integer is directly usable as an array index (see
//! [`crate::index::KmerIndex`]) — no hashing, no string comparisons.
//!
//! **Trick — spaced seeds.** Somatic hypermutation scatters substitutions, and a
//! contiguous seed is destroyed by any mismatch inside it. A spaced seed such
//! as `1101101111` ignores the "0" positions, so a substitution there leaves the
//! seed intact (PatternHunter; also used by modern long-read mappers). The
//! shape is data, not a type: switching from contiguous to spaced seeds changes
//! nothing downstream. Contiguous shapes (the default) skip the gather in
//! [`SeedShape::key`] entirely.

use std::fmt;

use crate::alphabet::N;

/// Longest supported seed span (2 bits per base must fit a `u64` with room for
/// the mask arithmetic).
pub const MAX_SPAN: usize = 31;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedShape {
    span: usize,
    /// Bit shifts (inside the rolling window) of the "care" positions, oldest first.
    shifts: Vec<u8>,
    contiguous: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SeedShapeError {
    Empty,
    TooLong(usize),
    InvalidSymbol(char),
    /// Patterns must start and end with a care position (`1`).
    NotAnchored,
}

impl fmt::Display for SeedShapeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "seed pattern is empty"),
            Self::TooLong(len) => write!(f, "seed span {len} exceeds {MAX_SPAN}"),
            Self::InvalidSymbol(c) => write!(f, "seed pattern symbol '{c}' is not 0 or 1"),
            Self::NotAnchored => write!(f, "seed pattern must start and end with '1'"),
        }
    }
}

impl std::error::Error for SeedShapeError {}

impl SeedShape {
    /// Contiguous k-mer. Panics unless `1 <= k <= MAX_SPAN` (a programming error,
    /// like an out-of-range slice index).
    pub fn contiguous(k: usize) -> Self {
        assert!(
            (1..=MAX_SPAN).contains(&k),
            "k must be in 1..={MAX_SPAN}, got {k}"
        );
        Self::from_pattern(&"1".repeat(k)).expect("an all-ones pattern of valid length is valid")
    }

    /// Spaced seed from a `0`/`1` pattern, e.g. `"1101101111"`.
    pub fn from_pattern(pattern: &str) -> Result<Self, SeedShapeError> {
        let bytes = pattern.as_bytes();
        let span = bytes.len();
        if span == 0 {
            return Err(SeedShapeError::Empty);
        }
        if span > MAX_SPAN {
            return Err(SeedShapeError::TooLong(span));
        }
        if let Some(&bad) = bytes.iter().find(|&&b| b != b'0' && b != b'1') {
            return Err(SeedShapeError::InvalidSymbol(char::from(bad)));
        }
        if bytes[0] != b'1' || bytes[span - 1] != b'1' {
            return Err(SeedShapeError::NotAnchored);
        }
        let shifts = bytes
            .iter()
            .enumerate()
            .filter(|&(_, &b)| b == b'1')
            .map(|(p, _)| u8::try_from(2 * (span - 1 - p)).expect("span <= 31 keeps shifts < 64"))
            .collect::<Vec<_>>();
        let contiguous = shifts.len() == span;
        Ok(Self {
            span,
            shifts,
            contiguous,
        })
    }

    /// Number of query bases a seed covers.
    pub fn span(&self) -> usize {
        self.span
    }

    /// Number of care positions.
    pub fn weight(&self) -> usize {
        self.shifts.len()
    }

    /// Bits in a seed key (`2 × weight`).
    pub fn key_bits(&self) -> u32 {
        2 * u32::try_from(self.weight()).expect("weight <= 31")
    }

    /// Key of a rolling window holding the last `span` bases (newest in the low bits).
    #[inline]
    pub fn key(&self, window: u64) -> u64 {
        if self.contiguous {
            return window;
        }
        self.shifts
            .iter()
            .fold(0u64, |key, &shift| (key << 2) | ((window >> shift) & 3))
    }

    /// Iterate `(start, key)` over all seeds of `codes` that contain no `N`.
    pub fn kmers<'a>(&'a self, codes: &'a [u8]) -> Kmers<'a> {
        Kmers {
            shape: self,
            codes,
            next: 0,
            window: 0,
            valid: 0,
            mask: (1u64 << (2 * self.span)) - 1,
        }
    }
}

/// Which k-mers of a sequence are seeds. Applied alike to the germline index
/// and to the query, so a seed both share is kept in both or in neither
/// (exactly so for the context-free syncmers; minimisers depend on
/// neighbouring k-mers). Densities are for random sequence.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SeedSampling {
    /// Every k-mer.
    #[default]
    All,
    /// Closed syncmers: k-mers whose smallest s-mer (`s = k − w`, ordered by
    /// a hash) is first or last. Density `2/(w + 1)`. Any `w + 1`
    /// consecutive k-mers contain one (the smallest s-mer they hold is first
    /// in, or last in, one of them), so every unmutated stretch of `k + w`
    /// bases yields a seed that the germline index holds too.
    ClosedSyncmer { w: u8 },
    /// Open syncmers: the smallest s-mer sits in the middle. Density `1/(w + 1)`.
    OpenSyncmer { w: u8 },
    /// Minimisers: the smallest k-mer (by hash) of every `w` consecutive
    /// k-mers. Density about `2/(w + 1)`.
    Minimizer { w: u8 },
}

/// SplitMix64 finaliser: a fixed random order on keys.
#[inline]
fn order(x: u64) -> u64 {
    let x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    let x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

impl SeedSampling {
    /// Seeds of `codes` under `shape` (contiguous shapes for syncmers);
    /// `ranks` is scratch.
    pub fn seeds(
        &self,
        shape: &SeedShape,
        codes: &[u8],
        out: &mut Vec<(usize, u64)>,
        ranks: &mut Vec<u64>,
    ) {
        out.clear();
        let k = shape.span();
        match *self {
            Self::All => out.extend(shape.kmers(codes)),
            Self::ClosedSyncmer { w } | Self::OpenSyncmer { w } => {
                let w = usize::from(w).min(k - 1);
                let s = k - w;
                // Rank of the s-mer starting at each position, hashed once.
                ranks.clear();
                ranks.resize(codes.len(), u64::MAX);
                let mask = (1u64 << (2 * s)) - 1;
                let (mut window, mut valid) = (0u64, 0usize);
                for (i, &code) in codes.iter().enumerate() {
                    if code >= N {
                        (window, valid) = (0, 0);
                        continue;
                    }
                    window = ((window << 2) | u64::from(code)) & mask;
                    valid += 1;
                    if valid >= s {
                        ranks[i + 1 - s] = order(window);
                    }
                }
                // The k-mer at `pos` holds the s-mers at `pos ..= pos + w`;
                // keep it when the first smallest one sits where required.
                let wanted = |at: usize| match self {
                    Self::ClosedSyncmer { .. } => at == 0 || at == w,
                    _ => at == w / 2,
                };
                out.extend(shape.kmers(codes).filter(|&(pos, _)| {
                    let inner = &ranks[pos..=pos + w];
                    let mut at = 0;
                    for (i, &rank) in inner.iter().enumerate().skip(1) {
                        if rank < inner[at] {
                            at = i;
                        }
                    }
                    wanted(at)
                }));
            }
            Self::Minimizer { w } => {
                let w = usize::from(w).max(1);
                let kmers: Vec<(usize, u64)> = shape.kmers(codes).collect();
                let mut start = 0;
                while start < kmers.len() {
                    // A run of consecutive k-mers (ambiguity codes break runs).
                    let mut end = start + 1;
                    while end < kmers.len() && kmers[end].0 == kmers[end - 1].0 + 1 {
                        end += 1;
                    }
                    let run = &kmers[start..end];
                    for window in run.windows(w.min(run.len())) {
                        let best = *window
                            .iter()
                            .min_by_key(|&&(pos, key)| (order(key), pos))
                            .expect("non-empty window");
                        if out.last() != Some(&best) {
                            out.push(best);
                        }
                    }
                    start = end;
                }
            }
        }
    }
}

/// Rolling seed iterator: O(1) work per base, ambiguity resets the window.
pub struct Kmers<'a> {
    shape: &'a SeedShape,
    codes: &'a [u8],
    next: usize,
    window: u64,
    valid: usize,
    mask: u64,
}

impl Iterator for Kmers<'_> {
    type Item = (usize, u64);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.next < self.codes.len() {
            let i = self.next;
            self.next += 1;
            let code = self.codes[i];
            if code >= N {
                self.valid = 0;
                self.window = 0;
                continue;
            }
            self.window = ((self.window << 2) | u64::from(code)) & self.mask;
            self.valid += 1;
            if self.valid >= self.shape.span {
                return Some((i + 1 - self.shape.span, self.shape.key(self.window)));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;

    #[test]
    fn contiguous_keys_are_two_bit_packed() {
        let shape = SeedShape::contiguous(3);
        let codes = encode(b"ACGTA");
        let kmers: Vec<_> = shape.kmers(&codes).collect();
        // ACG = 00 01 10, CGT = 01 10 11, GTA = 10 11 00
        assert_eq!(kmers, vec![(0, 0b000110), (1, 0b011011), (2, 0b101100)]);
    }

    #[test]
    fn ambiguity_breaks_seeds() {
        let shape = SeedShape::contiguous(3);
        let codes = encode(b"ACNGTAC");
        let starts: Vec<_> = shape.kmers(&codes).map(|(s, _)| s).collect();
        assert_eq!(starts, vec![3, 4]);
    }

    #[test]
    fn spaced_seed_ignores_dont_care_positions() {
        let shape = SeedShape::from_pattern("101").unwrap();
        let a = encode(b"AGC");
        let b = encode(b"ATC");
        let ka = shape.kmers(&a).next().unwrap().1;
        let kb = shape.kmers(&b).next().unwrap().1;
        assert_eq!(ka, kb);
        assert_eq!(shape.key_bits(), 4);
    }

    #[test]
    fn rejects_bad_patterns() {
        assert_eq!(
            SeedShape::from_pattern("0110"),
            Err(SeedShapeError::NotAnchored)
        );
        assert_eq!(
            SeedShape::from_pattern("1x1"),
            Err(SeedShapeError::InvalidSymbol('x'))
        );
    }

    fn random_codes(len: usize, seed: u64) -> Vec<u8> {
        crate::sim::Rng::new(seed).codes(len)
    }

    #[test]
    fn closed_syncmers_are_context_free_dense_and_windowed() {
        let shape = SeedShape::contiguous(9);
        let (mut seeds, mut ranks) = (Vec::new(), Vec::new());
        for w in [2u8, 4, 6] {
            let sampling = SeedSampling::ClosedSyncmer { w };
            let codes = random_codes(20_000, u64::from(w));
            sampling.seeds(&shape, &codes, &mut seeds, &mut ranks);
            let kmers = codes.len() - 8;
            let density = seeds.len() as f64 / kmers as f64;
            let expected = 2.0 / (f64::from(w) + 1.0);
            assert!(
                (density / expected - 1.0).abs() < 0.1,
                "w {w}: density {density}"
            );
            // Every w + 1 consecutive k-mers contain a syncmer.
            let gap = seeds.windows(2).map(|p| p[1].0 - p[0].0).max().unwrap();
            assert!(gap <= usize::from(w) + 1, "w {w}: gap {gap}");
            // Context-free: each selected k-mer is selected on its own.
            let (mut alone, mut alone_ranks) = (Vec::new(), Vec::new());
            for &(pos, key) in seeds.iter().take(200) {
                sampling.seeds(&shape, &codes[pos..pos + 9], &mut alone, &mut alone_ranks);
                assert_eq!(alone, vec![(0, key)]);
            }
        }
    }

    #[test]
    fn sampling_keeps_a_subset_of_the_kmers() {
        let shape = SeedShape::contiguous(9);
        let mut codes = random_codes(2_000, 7);
        codes[500] = N; // ambiguity breaks k-mers and minimiser runs
        let all: Vec<(usize, u64)> = shape.kmers(&codes).collect();
        let (mut seeds, mut ranks) = (Vec::new(), Vec::new());
        for sampling in [
            SeedSampling::All,
            SeedSampling::ClosedSyncmer { w: 4 },
            SeedSampling::OpenSyncmer { w: 4 },
            SeedSampling::Minimizer { w: 4 },
        ] {
            sampling.seeds(&shape, &codes, &mut seeds, &mut ranks);
            assert!(!seeds.is_empty());
            assert!(seeds.iter().all(|s| all.contains(s)), "{sampling:?}");
            assert!(
                seeds.windows(2).all(|p| p[0].0 < p[1].0),
                "{sampling:?}: positions ascend"
            );
        }
    }
}
