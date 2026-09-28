//! V anchors by motif: Cys23, Trp41 and Cys104 placed jointly.
//!
//! One pass translates all three frames. Cys104 candidates are the few `C`
//! codons near the 3′ end; for each, the best Trp41 is searched in a fixed
//! window upstream and the best Cys23 upstream of that. Every candidate
//! triple is scored from conserved context and spacing, and stop codons
//! upstream of Cys104 are heavily penalised, which also selects the frame.
//! The work per allele is one translation plus a few dozen score lookups.
//!
//! Conserved context used (IMGT numbering):
//!
//! | residue | context | heavy example | light example |
//! |---|---|---|---|
//! | Cys23 | `[ST]` before it | `VSC` | `ITC` |
//! | Trp41 | `P` at +5, `G` at +6, `Q` at +3 | `WVRQAPG` | `WYQQKPG` |
//! | Cys104 | `[YF]` at −2, `[YFH]` at −1 | `AVYYC` | `ATYYC` |
//!
//! Spacing (codons): Trp41 → Cys104 44–70 (CDR2 of 0–12 aa plus FR3),
//! Cys23 → Trp41 10–20 (CDR1 of 4–14 aa), Cys104 → 3′ end 1–15.

use crate::translate::{translate_frames, STOP};

/// Nucleotides from the Cys104 codon start to the V 3′ end.
const TAIL_NT: std::ops::RangeInclusive<usize> = 3..=45;
const TRP41_TO_CYS104: std::ops::RangeInclusive<usize> = 44..=70;
const CYS23_TO_TRP41: std::ops::RangeInclusive<usize> = 10..=20;
const STOP_PENALTY: i32 = 15;
/// Largest bonus from the database-consensus tail prior.
const PRIOR_WEIGHT: i32 = 6;
/// Bonus for a Cys23 at the database-consensus codon index; exceeds the
/// largest context difference between two Cys23 candidates (3).
const CYS23_PRIOR_WEIGHT: i32 = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VHit {
    pub frame: usize,
    /// Codon indices within `frame`.
    pub cys104: usize,
    pub trp41: Option<usize>,
    pub cys23: Option<usize>,
    pub score: i32,
    /// Conserved context around Cys104 and Trp41 and no upstream stop codon.
    pub confident: bool,
    /// `+1` for a heavy-like Trp41 context (`W.R`), `−1` for light-like (`W.Q`).
    pub chain_vote: i8,
    /// Other `C` codons in the Cys23 window (ambiguity worth reporting).
    pub cys23_alternatives: usize,
}

impl VHit {
    fn nt(&self, codon: usize) -> u32 {
        (3 * codon + self.frame) as u32
    }

    pub fn anchor(&self) -> u32 {
        self.nt(self.cys104)
    }

    pub fn trp41_nt(&self) -> Option<u32> {
        self.trp41.map(|k| self.nt(k))
    }

    pub fn cys23_nt(&self) -> Option<u32> {
        self.cys23.map(|k| self.nt(k))
    }
}

/// Search constraints.
#[derive(Clone, Copy, Debug, Default)]
pub struct VSearch {
    /// Database-consensus distance from Cys104 to the 3′ end.
    pub tail_prior: Option<u32>,
    /// Database-consensus `(frame, codon index)` of Cys23. Full-length alleles
    /// share it (IMGT 23 minus the locus' FR1 gaps); it is applied only to
    /// alleles in the same frame, since a 5′-truncated allele's codon indices
    /// are shifted.
    pub cys23_prior: Option<(usize, usize)>,
    /// Trusted anchor (header or earlier pass): fixes frame and Cys104.
    pub fixed_anchor: Option<u32>,
    /// FR3 length in codons once the chain type is known: CDR2 of 0–12 codons
    /// then confines Trp41 → Cys104 to `[14 + L, 26 + L]`.
    pub fr3_codons: Option<u32>,
}

#[inline]
fn at(aa: &[u8], k: isize) -> u8 {
    usize::try_from(k)
        .ok()
        .and_then(|k| aa.get(k).copied())
        .unwrap_or(b'-')
}

/// `YYC` dominates germline Cys104 contexts; `F`/`H` variants score less.
fn cys104_context(aa: &[u8], k: usize) -> (i32, bool) {
    let k = k as isize;
    let flank1 = match at(aa, k - 1) {
        b'Y' => 5,
        b'F' | b'H' => 3,
        _ => 0,
    };
    let flank2 = match at(aa, k - 2) {
        b'Y' => 5,
        b'F' => 3,
        _ => 0,
    };
    let flank3 = matches!(
        at(aa, k - 3),
        b'V' | b'A' | b'T' | b'M' | b'I' | b'L' | b'D' | b'G' | b'E' | b'S'
    );
    (flank1 + flank2 + i32::from(flank3), flank1 + flank2 > 0)
}

fn trp41_context(aa: &[u8], k: usize) -> (i32, bool, i8) {
    let k = k as isize;
    let proline = at(aa, k + 5) == b'P';
    let glycine = at(aa, k + 6) == b'G';
    let glutamine = at(aa, k + 3) == b'Q';
    let second = at(aa, k + 2);
    let vote = match second {
        b'R' => 1,
        b'Q' => -1,
        _ => 0,
    };
    let score = 3 * i32::from(proline)
        + 3 * i32::from(glycine)
        + 2 * i32::from(glutamine)
        + i32::from(matches!(second, b'R' | b'Q'));
    (score, proline || glycine, vote)
}

fn cys23_context(aa: &[u8], k: usize) -> i32 {
    let k = k as isize;
    2 * i32::from(matches!(at(aa, k - 1), b'S' | b'T'))
        + i32::from(matches!(at(aa, k - 2), b'V' | b'I' | b'L' | b'M'))
}

/// Best `(position, score)` of `residue` over codons `range`, scored by `context`.
fn best_in(
    aa: &[u8],
    range: std::ops::RangeInclusive<isize>,
    residue: u8,
    context: impl Fn(usize) -> i32,
) -> Option<(usize, i32)> {
    let lo = (*range.start()).max(0) as usize;
    let hi = (*range.end()).min(aa.len() as isize - 1);
    if hi < 0 {
        return None;
    }
    (lo..=hi as usize)
        .filter(|&k| aa[k] == residue)
        .map(|k| (k, context(k)))
        .max_by_key(|&(k, score)| (score, std::cmp::Reverse(k)))
}

/// Best V anchor hit over the three frames. `frames` is scratch space.
pub fn find_v(codes: &[u8], search: VSearch, frames: &mut [Vec<u8>; 3]) -> Option<VHit> {
    translate_frames(codes, frames);
    let n = codes.len();
    let mut best: Option<VHit> = None;
    for (frame, aa) in frames.iter().enumerate() {
        let fixed = search
            .fixed_anchor
            .filter(|&a| a as usize % 3 == frame)
            .map(|a| a as usize / 3);
        if search.fixed_anchor.is_some() && fixed.is_none() {
            continue;
        }
        for cys104 in 0..aa.len() {
            let tail = n - (3 * cys104 + frame);
            let candidate = match fixed {
                Some(k) => k == cys104,
                None => aa[cys104] == b'C' && TAIL_NT.contains(&tail),
            };
            if !candidate {
                continue;
            }
            let (c104_bonus, c104_ok) = cys104_context(aa, cys104);
            let mut score = 20 + c104_bonus;
            if let Some(mode) = search.tail_prior {
                score += PRIOR_WEIGHT - (tail as i32 - mode as i32).abs().min(PRIOR_WEIGHT);
            }
            let k = cys104 as isize;
            let (near, far) = match search.fr3_codons {
                Some(l) => (14 + l as isize, 26 + l as isize),
                None => (
                    *TRP41_TO_CYS104.start() as isize,
                    *TRP41_TO_CYS104.end() as isize,
                ),
            };
            let trp = best_in(aa, k - far..=k - near, b'W', |w| trp41_context(aa, w).0);
            let (mut trp41_ok, mut chain_vote) = (false, 0);
            if let Some((w, bonus)) = trp {
                let (_, ok, vote) = trp41_context(aa, w);
                score += 20 + bonus;
                trp41_ok = ok;
                chain_vote = vote;
            }
            let prior = search
                .cys23_prior
                .filter(|&(f, _)| f == frame)
                .map(|(_, c)| c);
            let window = trp.map(|(w, _)| {
                let w = w as isize;
                w - *CYS23_TO_TRP41.end() as isize..=w - *CYS23_TO_TRP41.start() as isize
            });
            let cys = window.clone().and_then(|window| {
                best_in(aa, window, b'C', |c| {
                    cys23_context(aa, c) + CYS23_PRIOR_WEIGHT * i32::from(prior == Some(c))
                })
            });
            let cys23_alternatives = window.map_or(0, |window| {
                window
                    .filter(|&c| c >= 0 && at(aa, c) == b'C')
                    .count()
                    .saturating_sub(1)
            });
            if let Some((_, bonus)) = cys {
                score += 12 + bonus;
            }
            let stops = aa[..cys104].iter().filter(|&&a| a == STOP).count() as i32;
            score -= STOP_PENALTY * stops;
            let hit = VHit {
                frame,
                cys104,
                trp41: trp.map(|(w, _)| w),
                cys23: cys.map(|(c, _)| c),
                score,
                confident: c104_ok && trp41_ok && stops == 0,
                chain_vote,
                cys23_alternatives,
            };
            if best.map_or(true, |b| {
                (hit.score, hit.confident) > (b.score, b.confident)
            }) {
                best = Some(hit);
            }
        }
    }
    best
}
