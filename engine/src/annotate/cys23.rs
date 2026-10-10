//! Cys23 learned from the database's unambiguous alleles.
//!
//! When the Cys23 window of an allele holds several `C` codons, or none (a
//! mutated Cys23), and no annotated relative settles it, choosing Cys23 is a
//! small classification problem. The database's alleles with a given CDR1
//! start or a single `C` in the window are labelled examples, and three
//! independent pieces of evidence are learned from them:
//!
//! - `x`: Cys23's nucleotide offset from the sequence start. Full-length
//!   alleles of a locus share a few (IMGT 23 minus the FR1 gaps); counting
//!   nucleotides, not codons of one frame, keeps it meaningful for an allele
//!   whose 5′ end lost one or two bases.
//! - `ℓ`: the CDR1 length the candidate implies, `Trp41 − Cys23 − 6` codons,
//!   which no 5′ truncation changes.
//! - the residues at −1 and −2, as a likelihood ratio against the residue
//!   frequencies upstream of Trp41 (`[ST]C` is an IG feature; TRB Cys23 is
//!   rarely preceded by S or T).
//!
//! ```text
//! P(k is Cys23) ∝ P_x(3k + frame) · P_ℓ(Trp41 − k − 6) · Π_{j=1,2} P_j(aa[k − j]) / P_bg(aa[k − j])
//! ```
//!
//! normalised over the candidates. Every table is a Jeffreys (½) estimate,
//! so an unseen value is unlikely, never impossible: `P_x` over the offsets
//! up to one codon past the largest seen (anything further gets the same
//! floor), `P_ℓ` over 4–14 codons, the residue tables over 20 amino acids,
//! stop and unknown.

use std::ops::RangeInclusive;

use crate::translate::{STOP, UNKNOWN};

/// CDR1 lengths (codons) the motif window admits: Trp41 − Cys23 = 10..=20.
const CDR1: RangeInclusive<usize> = 4..=14;
/// Residues of the context tables: the 20 amino acids, stop and unknown.
const ALPHABET: &[u8; 22] = b"ACDEFGHIKLMNPQRSTVWY*X";

fn residue_index(aa: u8) -> usize {
    let aa = if aa == STOP { b'*' } else { aa };
    ALPHABET
        .iter()
        .position(|&a| a == aa)
        .unwrap_or_else(|| residue_index(UNKNOWN))
}

/// Jeffreys estimate of a count among `total` over `bins` values.
fn jeffreys(count: u32, total: u32, bins: usize) -> f64 {
    (f64::from(count) + 0.5) / (f64::from(total) + 0.5 * bins as f64)
}

/// One labelled Cys23: the translation in the allele's frame and the codon
/// indices of Cys23 and Trp41 in it.
pub struct Label<'a> {
    pub aa: &'a [u8],
    pub frame: usize,
    pub cys23: usize,
    pub trp41: usize,
}

/// Cys23 offset, CDR1 length and context of the labelled alleles.
#[derive(Clone, Debug)]
pub struct Cys23Model {
    /// Labels per nucleotide offset of Cys23.
    offsets: Vec<u32>,
    /// Labels per CDR1 length (codons).
    cdr1: [u32; 15],
    /// Residue counts at −1 and −2 of the labelled Cys23.
    context: [[u32; 22]; 2],
    /// Residue counts upstream of Trp41 in the labelled alleles.
    background: [u32; 22],
    labels: u32,
}

impl Cys23Model {
    /// `None` without a label.
    pub fn learn<'a>(labels: impl IntoIterator<Item = Label<'a>>) -> Option<Self> {
        let mut model = Self {
            offsets: Vec::new(),
            cdr1: [0; 15],
            context: [[0; 22]; 2],
            background: [0; 22],
            labels: 0,
        };
        for label in labels {
            let x = 3 * label.cys23 + label.frame;
            if model.offsets.len() <= x {
                model.offsets.resize(x + 1, 0);
            }
            model.offsets[x] += 1;
            if let Some(l) = label.trp41.checked_sub(label.cys23 + 6)
                && CDR1.contains(&l)
            {
                model.cdr1[l] += 1;
            }
            for (j, table) in model.context.iter_mut().enumerate() {
                if let Some(&aa) = label.cys23.checked_sub(j + 1).and_then(|k| label.aa.get(k)) {
                    table[residue_index(aa)] += 1;
                }
            }
            for &aa in &label.aa[..label.trp41.min(label.aa.len())] {
                model.background[residue_index(aa)] += 1;
            }
            model.labels += 1;
        }
        (model.labels > 0).then_some(model)
    }

    /// Number of labelled alleles.
    pub fn labels(&self) -> u32 {
        self.labels
    }

    fn ln_p_offset(&self, x: usize) -> f64 {
        let bins = self.offsets.len() + 3;
        let count = self.offsets.get(x).copied().unwrap_or(0);
        jeffreys(count, self.labels, bins).ln()
    }

    fn ln_p_cdr1(&self, l: usize) -> f64 {
        if !CDR1.contains(&l) {
            return f64::NEG_INFINITY;
        }
        let total = self.cdr1.iter().sum();
        jeffreys(self.cdr1[l], total, CDR1.count()).ln()
    }

    fn ln_context(&self, aa: &[u8], k: usize) -> f64 {
        let background_total = self.background.iter().sum();
        let mut sum = 0.0;
        for (j, table) in self.context.iter().enumerate() {
            let Some(&residue) = k.checked_sub(j + 1).and_then(|i| aa.get(i)) else {
                continue;
            };
            let r = residue_index(residue);
            let p = jeffreys(table[r], table.iter().sum(), ALPHABET.len());
            let q = jeffreys(self.background[r], background_total, ALPHABET.len());
            sum += (p / q).ln();
        }
        sum
    }

    /// Whether a labelled allele has its Cys23 at nucleotide offset `x`:
    /// evidence that a candidate there sits in a full-length 5′ end.
    pub fn seen_offset(&self, x: usize) -> bool {
        self.offsets.get(x).is_some_and(|&n| n > 0)
    }

    /// Log score of codon `k` of `aa` (frame `frame`) as Cys23 given Trp41
    /// at codon `trp41`.
    pub fn ln_score(&self, aa: &[u8], frame: usize, k: usize, trp41: usize) -> f64 {
        let Some(l) = trp41.checked_sub(k + 6) else {
            return f64::NEG_INFINITY;
        };
        self.ln_p_offset(3 * k + frame) + self.ln_p_cdr1(l) + self.ln_context(aa, k)
    }

    /// The most probable Cys23 among `candidates` and its posterior among
    /// them; ties go to the later codon (nearer Trp41).
    pub fn choose(
        &self,
        aa: &[u8],
        frame: usize,
        trp41: usize,
        candidates: impl Iterator<Item = usize>,
    ) -> Option<(usize, f64)> {
        let scored: Vec<(usize, f64)> = candidates
            .map(|k| (k, self.ln_score(aa, frame, k, trp41)))
            .filter(|&(_, s)| s.is_finite())
            .collect();
        let best = scored
            .iter()
            .copied()
            .max_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)))?;
        let total: f64 = scored.iter().map(|&(_, s)| (s - best.1).exp()).sum();
        Some((best.0, 1.0 / total))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The documented posterior written out again from raw label lists.
    fn oracle(labels: &[Label<'_>], aa: &[u8], frame: usize, trp41: usize, k: usize) -> f64 {
        let n = labels.len() as f64;
        let offsets: Vec<usize> = labels.iter().map(|l| 3 * l.cys23 + l.frame).collect();
        let bins = (*offsets.iter().max().unwrap() + 4) as f64;
        let x = 3 * k + frame;
        let p_x = (offsets.iter().filter(|&&o| o == x).count() as f64 + 0.5) / (n + 0.5 * bins);
        let lengths: Vec<usize> = labels
            .iter()
            .filter_map(|l| l.trp41.checked_sub(l.cys23 + 6))
            .filter(|l| (4..=14).contains(l))
            .collect();
        let l = trp41 - k - 6;
        if !(4..=14).contains(&l) {
            return 0.0;
        }
        let p_l = (lengths.iter().filter(|&&m| m == l).count() as f64 + 0.5)
            / (lengths.len() as f64 + 5.5);
        let upstream: Vec<u8> = labels
            .iter()
            .flat_map(|l| l.aa[..l.trp41].iter().copied())
            .collect();
        let mut ratio = 1.0;
        for j in 1..=2 {
            let Some(r) = k.checked_sub(j).map(|i| aa[i]) else {
                continue;
            };
            let at: Vec<u8> = labels
                .iter()
                .filter_map(|l| l.cys23.checked_sub(j).map(|i| l.aa[i]))
                .collect();
            let p =
                (at.iter().filter(|&&a| a == r).count() as f64 + 0.5) / (at.len() as f64 + 11.0);
            let q = (upstream.iter().filter(|&&a| a == r).count() as f64 + 0.5)
                / (upstream.len() as f64 + 11.0);
            ratio *= p / q;
        }
        p_x * p_l * ratio
    }

    /// Over labelled alleles of two FR1 lengths (Cys23 at codons 21 and
    /// 22) and two contexts, every candidate set of a probe with decoys
    /// gets the argmax and posterior of the brute-force product.
    #[test]
    fn choice_maximises_the_documented_posterior() {
        let heavy = b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPG";
        let light = b"DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPG";
        let tr = b"KAGVTQTPRYLIKTRGQQVTLRCSPISGHRSVSWYQQTPG";
        let mut labels = Vec::new();
        for _ in 0..7 {
            labels.push(Label {
                aa: heavy,
                frame: 0,
                cys23: 21,
                trp41: 35,
            });
        }
        for _ in 0..4 {
            labels.push(Label {
                aa: light,
                frame: 1,
                cys23: 22,
                trp41: 34,
            });
        }
        for _ in 0..2 {
            labels.push(Label {
                aa: tr,
                frame: 2,
                cys23: 22,
                trp41: 33,
            });
        }
        let model = Cys23Model::learn(labels.iter().map(|l| Label { ..*l })).unwrap();
        assert_eq!(model.labels(), 13);
        // Probe: the TR gene 5′-truncated by one codon, with decoy Cys
        // codons at every position of the window in turn.
        let mut probe = tr[3..].to_vec();
        let trp41 = 30;
        for decoy in trp41 - 20..=trp41 - 10 {
            let saved = probe[decoy];
            probe[decoy] = b'C';
            for frame in 0..3 {
                let candidates: Vec<usize> = (trp41 - 20..=trp41 - 10)
                    .filter(|&k| probe[k] == b'C' || k == 10)
                    .collect();
                let (k, p) = model
                    .choose(&probe, frame, trp41, candidates.iter().copied())
                    .unwrap();
                let weights: Vec<(usize, f64)> = candidates
                    .iter()
                    .map(|&c| (c, oracle(&labels, &probe, frame, trp41, c)))
                    .collect();
                let total: f64 = weights.iter().map(|w| w.1).sum();
                let best = weights
                    .iter()
                    .copied()
                    .max_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)))
                    .unwrap();
                assert_eq!(k, best.0, "decoy {decoy} frame {frame}");
                assert!((p - best.1 / total).abs() < 1e-9, "{p} {}", best.1 / total);
            }
            probe[decoy] = saved;
        }
        assert!(model.seen_offset(63) && model.seen_offset(67) && model.seen_offset(68));
        assert!(!model.seen_offset(66));
    }
}
