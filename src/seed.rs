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
//! nothing downstream. On x86-64 with BMI2 the gather in [`SeedShape::key`] is a
//! single `pext` instruction; the portable loop here is the reference version.

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
        assert!((1..=MAX_SPAN).contains(&k), "k must be in 1..={MAX_SPAN}, got {k}");
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
        Ok(Self { span, shifts, contiguous })
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
        assert_eq!(SeedShape::from_pattern("0110"), Err(SeedShapeError::NotAnchored));
        assert_eq!(SeedShape::from_pattern("1x1"), Err(SeedShapeError::InvalidSymbol('x')));
    }
}
