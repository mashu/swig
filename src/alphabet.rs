//! Nucleotide alphabet.
//!
//! Two representations are used on purpose:
//!
//! * **byte codes** (`0..=4`, one base per byte) for dynamic programming, where
//!   table lookups dominate and bit-packing would only add shifts and masks;
//! * **2-bit codes packed into a `u64`** for k-mers (see [`crate::seed`]), where a
//!   k-mer becomes an integer that directly indexes an array.
//!
//! Every non-`ACGT` symbol collapses to [`N`]. `N` scores 0 against anything and
//! breaks k-mers, which is exactly the treatment of ambiguity in the method note.

pub const A: u8 = 0;
pub const C: u8 = 1;
pub const G: u8 = 2;
pub const T: u8 = 3;
pub const N: u8 = 4;

/// Number of distinct codes (`A C G T N`).
pub const ALPHABET_SIZE: usize = 5;

const ENCODE: [u8; 256] = build_encode_table();

const fn build_encode_table() -> [u8; 256] {
    let mut table = [N; 256];
    table[b'A' as usize] = A;
    table[b'a' as usize] = A;
    table[b'C' as usize] = C;
    table[b'c' as usize] = C;
    table[b'G' as usize] = G;
    table[b'g' as usize] = G;
    table[b'T' as usize] = T;
    table[b't' as usize] = T;
    table[b'U' as usize] = T;
    table[b'u' as usize] = T;
    table
}

/// Encode one ASCII symbol.
#[inline]
pub fn encode_base(symbol: u8) -> u8 {
    ENCODE[usize::from(symbol)]
}

/// `true` for `A`, `C`, `G`, `T`.
#[inline]
pub fn is_canonical(code: u8) -> bool {
    code < N
}

/// Watson–Crick complement of a code; `N` stays `N`.
#[inline]
pub fn complement(code: u8) -> u8 {
    if code < N {
        3 - code
    } else {
        N
    }
}

/// Decode one code back to an upper-case ASCII symbol.
#[inline]
pub fn decode_base(code: u8) -> u8 {
    b"ACGTN"[usize::from(code.min(N))]
}

/// Encode an ASCII sequence into `out` (cleared first, so buffers can be reused).
pub fn encode_into(sequence: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(sequence.iter().map(|&b| encode_base(b)));
}

/// Encode an ASCII sequence into a fresh vector.
pub fn encode(sequence: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(sequence.len());
    encode_into(sequence, &mut out);
    out
}

/// Reverse complement of encoded `codes` into `out` (cleared first).
pub fn reverse_complement_into(codes: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(codes.iter().rev().map(|&c| complement(c)));
}

/// Decode codes into an upper-case string.
pub fn decode(codes: &[u8]) -> String {
    codes.iter().map(|&c| char::from(decode_base(c))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_and_ambiguity() {
        let codes = encode(b"ACGTuRn");
        assert_eq!(codes, vec![A, C, G, T, T, N, N]);
        assert_eq!(decode(&codes), "ACGTTNN");
    }

    #[test]
    fn reverse_complement_is_involution() {
        let codes = encode(b"AACGTNT");
        let mut rc = Vec::new();
        let mut back = Vec::new();
        reverse_complement_into(&codes, &mut rc);
        reverse_complement_into(&rc, &mut back);
        assert_eq!(decode(&rc), "ANACGTT");
        assert_eq!(back, codes);
    }
}
