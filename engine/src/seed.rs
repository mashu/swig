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

    /// Whether every position of the span is a care position.
    pub fn is_contiguous(&self) -> bool {
        self.contiguous
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
///
/// Syncmers rank the `w + 1` s-mers (`s = k − w`, `k` the shape's weight) of
/// a seed's key by a fixed random order and keep the seed when the first
/// smallest one sits where required (Edgar 2021). A *window guarantee* holds
/// when any `w` consecutive k-mers contain a seed: then every unmutated
/// stretch of `k + w − 1` bases shared by read and germline yields a seed
/// both sides hold. It needs neighbouring seeds to share s-mers, i.e. a
/// contiguous shape; on a spaced shape selection is still per key (so read
/// and germline agree), but seeds can be far apart.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SeedSampling {
    /// Every k-mer.
    #[default]
    All,
    /// Closed syncmers: the smallest s-mer is first or last. Density about
    /// `2/(w + 1)`. Context-free, with the window guarantee on contiguous
    /// shapes: of the `2w` s-mers that `w` consecutive k-mers hold, the first
    /// smallest is first in, or last in, one of them.
    ClosedSyncmer { w: u8 },
    /// Open syncmers: the smallest s-mer is at the middle position `w / 2`,
    /// the offset that maximises conservation (Shaw & Yu 2022). Density
    /// `1/(w + 1)`. Context-free, but without a window guarantee: gaps
    /// between seeds can be long, most of all on low-complexity sequence.
    OpenSyncmer { w: u8 },
    /// Minimisers: the smallest k-mer (by the same order) of every `w`
    /// consecutive k-mers. Density about `2/(w + 1)`, with the window
    /// guarantee, but not context-free: a mutation near a shared k-mer can
    /// unselect it on one side only.
    Minimizer { w: u8 },
}

/// Why a sampling cannot be used with a seed shape.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SamplingError {
    /// Syncmers need `1 ≤ w < k`: an s-mer of at least one base.
    SyncmerWindow { w: u8, k: usize },
    /// Minimisers need a window of at least one k-mer.
    EmptyWindow,
}

impl fmt::Display for SamplingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SyncmerWindow { w, k } => {
                write!(f, "syncmer w = {w} needs 1 <= w < k = {k}")
            }
            Self::EmptyWindow => write!(f, "minimiser w must be at least 1"),
        }
    }
}

impl std::error::Error for SamplingError {}

/// SplitMix64: a fixed random order on keys. The golden-ratio increment
/// keeps the all-zero key (poly-A) from being its own image, which the bare
/// finaliser would make the smallest of all.
#[inline]
fn order(x: u64) -> u64 {
    let x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    let x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// Position (0-based, the first if tied) of the smallest s-mer
/// (`s = k − w`) in a `k`-base key, oldest base in the high bits.
#[inline]
fn smallest_smer(key: u64, k: usize, w: usize) -> usize {
    let mask = (1u64 << (2 * (k - w))) - 1;
    let mut at = 0;
    let mut best = u64::MAX;
    for j in 0..=w {
        let rank = order((key >> (2 * (w - j))) & mask);
        if rank < best {
            (at, best) = (j, rank);
        }
    }
    at
}

impl SeedSampling {
    /// Whether this sampling is defined for seeds of `shape`.
    pub fn check(&self, shape: &SeedShape) -> Result<(), SamplingError> {
        let k = shape.weight();
        match *self {
            Self::All => Ok(()),
            Self::ClosedSyncmer { w } | Self::OpenSyncmer { w } => {
                if w == 0 || usize::from(w) >= k {
                    Err(SamplingError::SyncmerWindow { w, k })
                } else {
                    Ok(())
                }
            }
            Self::Minimizer { w } => {
                if w == 0 {
                    Err(SamplingError::EmptyWindow)
                } else {
                    Ok(())
                }
            }
        }
    }

    /// Seeds of `codes` under `shape`; `ranks` is scratch. Panics if
    /// [`SeedSampling::check`] fails (index construction checks first).
    pub fn seeds(
        &self,
        shape: &SeedShape,
        codes: &[u8],
        out: &mut Vec<(usize, u64)>,
        ranks: &mut Vec<u64>,
    ) {
        if let Err(error) = self.check(shape) {
            panic!("{error}");
        }
        out.clear();
        let k = shape.weight();
        match *self {
            Self::All => out.extend(shape.kmers(codes)),
            Self::ClosedSyncmer { w } | Self::OpenSyncmer { w } => {
                let w = usize::from(w);
                let wanted = |at: usize| match self {
                    Self::ClosedSyncmer { .. } => at == 0 || at == w,
                    _ => at == w / 2,
                };
                if !shape.is_contiguous() {
                    // The s-mers of the key itself, so that two seeds with
                    // one key are kept alike whatever their don't-care bases.
                    out.extend(
                        shape
                            .kmers(codes)
                            .filter(|&(_, key)| wanted(smallest_smer(key, k, w))),
                    );
                    return;
                }
                // Contiguous: the same definition with each s-mer hashed once.
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
                let w = usize::from(w);
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
    fn syncmers_are_context_free_dense_and_windowed() {
        let (mut seeds, mut ranks) = (Vec::new(), Vec::new());
        let (mut alone, mut alone_ranks) = (Vec::new(), Vec::new());
        for k in [7usize, 9] {
            let shape = SeedShape::contiguous(k);
            for w in [2u8, 4, 6] {
                let codes = random_codes(20_000, 31 * k as u64 + u64::from(w));
                let kmers = codes.len() + 1 - k;
                for (sampling, expected) in [
                    (
                        SeedSampling::ClosedSyncmer { w },
                        2.0 / (f64::from(w) + 1.0),
                    ),
                    (SeedSampling::OpenSyncmer { w }, 1.0 / (f64::from(w) + 1.0)),
                ] {
                    sampling.seeds(&shape, &codes, &mut seeds, &mut ranks);
                    let density = seeds.len() as f64 / kmers as f64;
                    // 4^s ≫ k is needed for the nominal density (k 7, w 6: s 1).
                    if 4f64.powi(k as i32 - i32::from(w)) > 4.0 * k as f64 {
                        assert!(
                            (density / expected - 1.0).abs() < 0.1,
                            "{sampling:?} k {k}: density {density}"
                        );
                    }
                    if let SeedSampling::ClosedSyncmer { .. } = sampling {
                        // Any w consecutive k-mers contain a closed syncmer.
                        let gap = seeds.windows(2).map(|p| p[1].0 - p[0].0).max().unwrap();
                        assert!(gap <= usize::from(w), "k {k} w {w}: gap {gap}");
                    }
                    // Context-free both ways: a k-mer is a seed in context
                    // exactly when it is a seed on its own.
                    let mut selected = seeds.iter().map(|&(pos, _)| pos).peekable();
                    for pos in 0..kmers.min(3_000) {
                        sampling.seeds(&shape, &codes[pos..pos + k], &mut alone, &mut alone_ranks);
                        let in_context = selected.next_if_eq(&pos).is_some();
                        assert_eq!(!alone.is_empty(), in_context, "{sampling:?} k {k} at {pos}");
                    }
                }
            }
        }
    }

    /// The window bound is tight: every 12-base string holds a closed
    /// (9, 5)-syncmer, and some 11-base strings hold none.
    #[test]
    fn closed_syncmer_window_bound_is_tight() {
        let (k, w) = (9usize, 4u8);
        let shape = SeedShape::contiguous(k);
        let sampling = SeedSampling::ClosedSyncmer { w };
        let (mut seeds, mut ranks) = (Vec::new(), Vec::new());
        let mut x = 99u64;
        let (mut empty_short, mut empty_full) = (0, 0);
        for _ in 0..200_000 {
            x = order(x);
            let codes: Vec<u8> = (0..k + usize::from(w) - 1)
                .map(|i| ((x >> (2 * i)) & 3) as u8)
                .collect();
            sampling.seeds(&shape, &codes, &mut seeds, &mut ranks);
            empty_full += usize::from(seeds.is_empty());
            sampling.seeds(&shape, &codes[..codes.len() - 1], &mut seeds, &mut ranks);
            empty_short += usize::from(seeds.is_empty());
        }
        assert_eq!(empty_full, 0);
        assert!(empty_short > 0);
    }

    /// Spaced shapes rank the s-mers of the key, so two seeds with one key
    /// are kept alike; on contiguous shapes the rolling ranks agree with it.
    #[test]
    fn syncmers_of_spaced_shapes_follow_the_key() {
        let (mut seeds, mut ranks) = (Vec::new(), Vec::new());
        let spaced = SeedShape::from_pattern("11011011101").unwrap();
        let codes = random_codes(20_000, 5);
        for sampling in [
            SeedSampling::ClosedSyncmer { w: 3 },
            SeedSampling::OpenSyncmer { w: 2 },
        ] {
            sampling.seeds(&spaced, &codes, &mut seeds, &mut ranks);
            let kept: std::collections::HashSet<u64> = seeds.iter().map(|&(_, key)| key).collect();
            for (_, key) in spaced.kmers(&codes) {
                let at = smallest_smer(
                    key,
                    spaced.weight(),
                    match sampling {
                        SeedSampling::ClosedSyncmer { w } | SeedSampling::OpenSyncmer { w } => {
                            usize::from(w)
                        }
                        _ => unreachable!(),
                    },
                );
                let wanted = match sampling {
                    SeedSampling::ClosedSyncmer { w } => at == 0 || at == usize::from(w),
                    SeedSampling::OpenSyncmer { w } => at == usize::from(w) / 2,
                    _ => unreachable!(),
                };
                assert_eq!(kept.contains(&key), wanted, "{sampling:?}");
            }
        }
        let shape = SeedShape::contiguous(9);
        let sampling = SeedSampling::ClosedSyncmer { w: 4 };
        sampling.seeds(&shape, &codes, &mut seeds, &mut ranks);
        let by_key: Vec<(usize, u64)> = shape
            .kmers(&codes)
            .filter(|&(_, key)| matches!(smallest_smer(key, 9, 4), 0 | 4))
            .collect();
        assert_eq!(seeds, by_key);
    }

    /// The all-zero key (poly-A) is not its own image under the order, so it
    /// does not rank first among the s-mers.
    #[test]
    fn poly_a_is_not_the_smallest_smer() {
        assert_ne!(order(0), 0);
        let below = (0u64..1024).filter(|&v| order(v) < order(0)).count();
        assert!(below > 0, "poly-A 5-mer must not rank first");
    }

    #[test]
    fn invalid_windows_are_rejected() {
        let shape = SeedShape::contiguous(9);
        assert_eq!(
            SeedSampling::ClosedSyncmer { w: 9 }.check(&shape),
            Err(SamplingError::SyncmerWindow { w: 9, k: 9 })
        );
        assert!(SeedSampling::OpenSyncmer { w: 0 }.check(&shape).is_err());
        assert_eq!(
            SeedSampling::Minimizer { w: 0 }.check(&shape),
            Err(SamplingError::EmptyWindow)
        );
        assert!(SeedSampling::ClosedSyncmer { w: 8 }.check(&shape).is_ok());
        // A spaced shape's k is its weight.
        let spaced = SeedShape::from_pattern("1101").unwrap();
        assert!(SeedSampling::ClosedSyncmer { w: 3 }.check(&spaced).is_err());
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
