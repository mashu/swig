//! J anchors by motif: Trp118 (heavy) / Phe118 (light) in `[WF]G.G`.
//!
//! All three frames are translated in one pass; every `W`/`F` codon is scored
//! by the two conserved glycines, the anchor's distance from the 3′ end (J
//! genes end ~10 codons after it) and an open reading frame through FR4.
//! The work per allele is one translation plus a scan of ~15 codons per frame.

use crate::translate::{translate_frames, STOP};

const STOP_PENALTY: i32 = 15;
const PRIOR_WEIGHT: i32 = 6;
/// Typical distance (nt) from the anchor codon start to the J 3′ end.
const TYPICAL_TAIL_NT: std::ops::RangeInclusive<usize> = 24..=42;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JHit {
    pub frame: usize,
    /// Codon index of Trp/Phe118 within `frame`.
    pub codon: usize,
    pub residue: u8,
    pub score: i32,
    /// Both glycines present and FR4 open.
    pub confident: bool,
}

impl JHit {
    pub fn anchor(&self) -> u32 {
        (3 * self.codon + self.frame) as u32
    }

    /// `+1` heavy (`W`), `−1` light (`F`).
    pub fn chain_vote(&self) -> i8 {
        match self.residue {
            b'W' => 1,
            b'F' => -1,
            _ => 0,
        }
    }
}

/// Best J anchor hit. `tail_prior` is the database-consensus distance from the
/// anchor to the 3′ end; `frames` is scratch space.
pub fn find_j(codes: &[u8], tail_prior: Option<u32>, frames: &mut [Vec<u8>; 3]) -> Option<JHit> {
    translate_frames(codes, frames);
    let n = codes.len();
    let mut best: Option<JHit> = None;
    for (frame, aa) in frames.iter().enumerate() {
        for codon in 0..aa.len() {
            let residue = aa[codon];
            if !matches!(residue, b'W' | b'F') {
                continue;
            }
            let glycine1 = aa.get(codon + 1) == Some(&b'G');
            let glycine3 = aa.get(codon + 3) == Some(&b'G');
            let middle = aa.get(codon + 2).is_some_and(|x| b"QKPAGSTREN".contains(x));
            let tail = n - (3 * codon + frame);
            let prior = match tail_prior {
                Some(mode) => PRIOR_WEIGHT - (tail as i32 - mode as i32).abs().min(PRIOR_WEIGHT),
                None => 3 * i32::from(TYPICAL_TAIL_NT.contains(&tail)),
            };
            let stops = aa[codon..].iter().filter(|&&a| a == STOP).count() as i32;
            let score = 20 + 5 * i32::from(glycine1) + 5 * i32::from(glycine3) + i32::from(middle) + prior
                - STOP_PENALTY * stops;
            let hit = JHit { frame, codon, residue, score, confident: glycine1 && glycine3 && stops == 0 };
            if best.map_or(true, |b| (hit.score, hit.confident) > (b.score, b.confident)) {
                best = Some(hit);
            }
        }
    }
    best
}
