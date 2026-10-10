//! Codon translation.
//!
//! **Trick — 2-bit codon index.** Bases are already 2-bit codes, so a codon
//! is the 6-bit number `16a + 4b + c` and translation is one lookup in a
//! 64-entry table built at compile time. Any N (code 4 has bit 2 set) makes
//! the codon unknown, tested with a single OR. All three frames are produced
//! in one pass over the sequence.

use crate::alphabet::N;

pub const STOP: u8 = b'*';
pub const UNKNOWN: u8 = b'X';

/// Standard code in the textbook T, C, A, G order.
const TCAG_CODE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// The table re-indexed by our base codes (A = 0, C = 1, G = 2, T = 3).
const TABLE: [u8; 64] = {
    // Our code -> position in T, C, A, G order.
    const TO_TCAG: [usize; 4] = [2, 1, 3, 0];
    let mut table = [0u8; 64];
    let mut i = 0;
    while i < 64 {
        let (a, b, c) = (i / 16, (i / 4) % 4, i % 4);
        table[i] = TCAG_CODE[16 * TO_TCAG[a] + 4 * TO_TCAG[b] + TO_TCAG[c]];
        i += 1;
    }
    table
};

/// Amino acid of one codon of base codes.
#[inline]
pub fn codon(a: u8, b: u8, c: u8) -> u8 {
    if (a | b | c) & N != 0 {
        UNKNOWN
    } else {
        TABLE[usize::from(a) * 16 + usize::from(b) * 4 + usize::from(c)]
    }
}

/// Translate complete codons starting at `frame` (0, 1 or 2), appending to `out`.
pub fn translate_into(codes: &[u8], frame: usize, out: &mut Vec<u8>) {
    out.extend(
        codes
            .get(frame..)
            .unwrap_or_default()
            .as_chunks::<3>()
            .0
            .iter()
            .map(|&[a, b, c]| codon(a, b, c)),
    );
}

pub fn translate(codes: &[u8], frame: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(codes.len() / 3);
    translate_into(codes, frame, &mut out);
    out
}

/// All three frames in one pass: `frames[f][k]` is the codon starting at `3k + f`.
pub fn translate_frames(codes: &[u8], frames: &mut [Vec<u8>; 3]) {
    frames.iter_mut().for_each(Vec::clear);
    for start in 0..codes.len().saturating_sub(2) {
        frames[start % 3].push(codon(codes[start], codes[start + 1], codes[start + 2]));
    }
}

/// The frames (0, 1, 2) whose translation has no stop codon except as its
/// last codon: the possible reading frames of a coding sequence such as a
/// constant region (`C` ends with its stop codon).
pub fn open_frames(codes: &[u8]) -> impl Iterator<Item = usize> + '_ {
    (0..3).filter(move |&frame| {
        let protein = translate(codes, frame);
        let body = protein.split_last().map_or(&[][..], |(_, body)| body);
        !body.contains(&STOP)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;

    #[test]
    fn translates_the_standard_code() {
        let codes = encode(b"ATGTGGTGTTAATAGTGACAGGGNTTC");
        assert_eq!(translate(&codes, 0), b"MWC***QXF");
        let mut frames = [Vec::new(), Vec::new(), Vec::new()];
        translate_frames(&codes, &mut frames);
        for (frame, translated) in frames.iter().enumerate() {
            assert_eq!(*translated, translate(&codes, frame));
        }
    }

    #[test]
    fn open_frames_allow_a_terminal_stop_only() {
        // Frame 0: M K * (terminal stop); frame 1: stop inside; frame 2: open.
        let codes = encode(b"ATGAAATAA");
        assert_eq!(open_frames(&codes).collect::<Vec<_>>(), vec![0, 2]);
        // Frame 0: * * * (stops inside); frames 1 and 2: N N N and I I.
        let frames: Vec<_> = open_frames(&encode(b"TAATAATAAT")).collect();
        assert_eq!(frames, vec![1, 2]);
    }
}
