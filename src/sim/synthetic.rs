//! Synthetic heavy-chain germline databases with realistic similarity structure.
//!
//! Regression tests cannot ship IMGT-derived databases, and random,
//! unrelated sequences make allele calling trivially easy. This generator
//! reproduces what makes real IGH databases hard — a hierarchy of
//! similarity — from a phylogeny-like process:
//!
//! ```text
//! V root (IGHV protein template, back-translated)
//!  └─ family ancestors     root + `family_divergence` substitutions
//!      └─ genes            family + `gene_divergence` substitutions
//!          └─ alleles      gene + 1..=`max_allele_snps` SNPs,
//!                          a few 5′-truncated (partial) alleles
//! ```
//!
//! The conserved codons (V Cys23, Trp41, Cys104 and the J Trp118 of
//! `WGXG`) are never mutated, so every allele has a known anchor. D genes
//! follow the same family → gene → allele scheme at D lengths (11–37 nt), and
//! J genes share one FR4 (`WGQGTLVTVSS`) after a gene-specific CDR3 part,
//! like IGHJ1–6.

use super::random::Rng;
use crate::germline::{GermlineSet, SegmentKind};
use crate::translate::codon;

/// Sizes and divergences of a synthetic database.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SyntheticSpec {
    pub v_families: usize,
    pub v_genes_per_family: usize,
    /// Alleles per gene are uniform in `1..=max_alleles`.
    pub max_alleles: usize,
    pub max_allele_snps: usize,
    /// Fraction of non-first alleles that are 5′-truncated partial sequences.
    pub partial_alleles: f64,
    pub family_divergence: f64,
    pub gene_divergence: f64,
    pub d_families: usize,
    pub d_genes_per_family: usize,
    pub j_genes: usize,
}

impl Default for SyntheticSpec {
    /// About 170 V, 35 D and 12 J alleles: IGH-like structure at a third of
    /// the size, so tests stay fast.
    fn default() -> Self {
        Self {
            v_families: 7,
            v_genes_per_family: 12,
            max_alleles: 3,
            max_allele_snps: 3,
            partial_alleles: 0.1,
            family_divergence: 0.14,
            gene_divergence: 0.05,
            d_families: 7,
            d_genes_per_family: 4,
            j_genes: 6,
        }
    }
}

/// A generated database. Every V and J allele carries its anchor.
pub struct SyntheticDatabase {
    pub v: GermlineSet,
    pub d: GermlineSet,
    pub j: GermlineSet,
}

/// IGHV protein template: FR1 (Cys23 at 21), CDR1, FR2 (Trp41 at 2), CDR2,
/// FR3 ending in Cys104, then the `AR` tail.
const V_FR1: &[u8] = b"QVQLVQSGAEVKKPGASVKVSCKAS";
const V_CDR1: &[u8] = b"GYTFTSYA";
const V_FR2: &[u8] = b"MHWVRQAPGQRLEWMGW";
const V_CDR2: &[u8] = b"INAGNGNT";
const V_FR3: &[u8] = b"KYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYC";
const V_TAIL: &[u8] = b"AR";
const J_FR4: &[u8] = b"WGQGTLVTVSS";

fn back_translate(protein: &[u8], rng: &mut Rng) -> Vec<u8> {
    let mut out = Vec::with_capacity(3 * protein.len());
    for &residue in protein {
        let options: Vec<[u8; 3]> = (0..64u8)
            .filter(|&i| codon(i / 16, (i / 4) % 4, i % 4) == residue)
            .map(|i| {
                [
                    b"ACGT"[usize::from(i / 16)],
                    b"ACGT"[usize::from((i / 4) % 4)],
                    b"ACGT"[usize::from(i % 4)],
                ]
            })
            .collect();
        out.extend(options[rng.below(options.len())]);
    }
    out
}

/// Substitute each unprotected base with probability `rate`.
fn mutate(
    sequence: &[u8],
    rate: f64,
    protected: &[std::ops::Range<usize>],
    rng: &mut Rng,
) -> Vec<u8> {
    sequence
        .iter()
        .enumerate()
        .map(|(i, &b)| {
            if !protected.iter().any(|r| r.contains(&i)) && rng.bernoulli(rate) {
                rng.substitute(b)
            } else {
                b
            }
        })
        .collect()
}

/// Exactly `count` substitutions at distinct unprotected positions.
fn snps(
    sequence: &[u8],
    count: usize,
    protected: &[std::ops::Range<usize>],
    rng: &mut Rng,
) -> Vec<u8> {
    let mut out = sequence.to_vec();
    let mut done = Vec::new();
    while done.len() < count {
        let at = rng.below(out.len());
        if done.contains(&at) || protected.iter().any(|r| r.contains(&at)) {
            continue;
        }
        out[at] = rng.substitute(out[at]);
        done.push(at);
    }
    out
}

impl SyntheticDatabase {
    pub fn generate(spec: &SyntheticSpec, seed: u64) -> Self {
        let mut rng = Rng::new(seed);
        Self {
            v: Self::v_set(spec, &mut rng),
            d: Self::d_set(spec, &mut rng),
            j: Self::j_set(spec, &mut rng),
        }
    }

    fn v_set(spec: &SyntheticSpec, rng: &mut Rng) -> GermlineSet {
        let protein: Vec<u8> = [V_FR1, V_CDR1, V_FR2, V_CDR2, V_FR3, V_TAIL].concat();
        let mut root = back_translate(&protein, rng);
        root.extend(b"GA");
        let codon_at = |index: usize| 3 * index..3 * index + 3;
        let cys23 = codon_at(21);
        let trp41 = codon_at(V_FR1.len() + V_CDR1.len() + 2);
        let cys104 = codon_at(protein.len() - V_TAIL.len() - 1);
        let anchor = cys104.start;
        // Keep Cys23, Trp41 and the YYC motif ending in Cys104.
        let protected = [cys23, trp41, cys104.start - 6..cys104.end];

        let mut set = GermlineSet::new(SegmentKind::V);
        for family in 1..=spec.v_families {
            let ancestor = mutate(&root, spec.family_divergence, &protected, rng);
            for gene in 1..=spec.v_genes_per_family {
                let first = mutate(&ancestor, spec.gene_divergence, &protected, rng);
                let alleles = 1 + rng.below(spec.max_alleles.max(1));
                for allele in 1..=alleles {
                    let name = format!("IGHV{family}-{gene}*{allele:02}");
                    if allele == 1 {
                        set.push(&name, &first, Some(anchor as u32))
                            .expect("non-empty");
                        continue;
                    }
                    let count = 1 + rng.below(spec.max_allele_snps.max(1));
                    let variant = snps(&first, count, &protected, rng);
                    // Partial alleles lack a few 5′ bases, as many IMGT entries do.
                    let cut = if rng.bernoulli(spec.partial_alleles) {
                        1 + rng.below(12)
                    } else {
                        0
                    };
                    set.push(&name, &variant[cut..], Some((anchor - cut) as u32))
                        .expect("non-empty");
                }
            }
        }
        set
    }

    fn d_set(spec: &SyntheticSpec, rng: &mut Rng) -> GermlineSet {
        let mut set = GermlineSet::new(SegmentKind::D);
        for family in 1..=spec.d_families {
            let length = 11 + rng.below(27);
            let ancestor: Vec<u8> = (0..length).map(|_| rng.base()).collect();
            for gene in 1..=spec.d_genes_per_family {
                let first = mutate(&ancestor, 0.15, &[], rng);
                set.push(&format!("IGHD{family}-{gene}*01"), &first, None)
                    .expect("non-empty");
                if rng.bernoulli(0.3) {
                    let second = snps(&first, 1, &[], rng);
                    set.push(&format!("IGHD{family}-{gene}*02"), &second, None)
                        .expect("non-empty");
                }
            }
        }
        set
    }

    fn j_set(spec: &SyntheticSpec, rng: &mut Rng) -> GermlineSet {
        let mut set = GermlineSet::new(SegmentKind::J);
        let fr4 = back_translate(J_FR4, rng);
        // The Trp118 codon opens FR4.
        let protected = std::slice::from_ref(&(0..3usize));
        for gene in 1..=spec.j_genes {
            let lead_length = 6 + rng.below(16);
            let mut first: Vec<u8> = (0..lead_length).map(|_| rng.base()).collect();
            let shift = first.len();
            first.extend(mutate(&fr4, 0.08, protected, rng).iter());
            first.push(b'G');
            let shifted = shift..shift + 3;
            set.push(&format!("IGHJ{gene}*01"), &first, Some(shift as u32))
                .expect("non-empty");
            if rng.bernoulli(0.6) {
                let second = snps(&first, 1, std::slice::from_ref(&shifted), rng);
                set.push(&format!("IGHJ{gene}*02"), &second, Some(shift as u32))
                    .expect("non-empty");
            }
        }
        set
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::translate::translate;

    #[test]
    fn anchors_are_conserved_codons_and_structure_is_hierarchical() {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 5);
        for allele in db.v.alleles() {
            let at = allele.anchor.expect("V anchors are known") as usize;
            let aa = translate(&allele.codes[at..at + 3], 0);
            assert_eq!(aa, b"C", "{}", allele.name);
        }
        for allele in db.j.alleles() {
            let at = allele.anchor.expect("J anchors are known") as usize;
            assert_eq!(
                translate(&allele.codes[at..at + 3], 0),
                b"W",
                "{}",
                allele.name
            );
        }
        let identity = |a: &[u8], b: &[u8]| {
            let n = a.len().min(b.len());
            a[..n].iter().zip(&b[..n]).filter(|(x, y)| x == y).count() as f64 / n as f64
        };
        let named = |prefix: &str| {
            db.v.alleles()
                .iter()
                .find(|a| a.name.starts_with(prefix))
                .expect("allele exists")
                .codes
                .clone()
        };
        let first = named("IGHV1-1*01");
        let same_family = identity(&first, &named("IGHV1-2*01"));
        let other_family = identity(&first, &named("IGHV7-1*01"));
        assert!(
            same_family > other_family,
            "{same_family} vs {other_family}"
        );
        assert!(db.d.len() >= 28 && db.j.len() >= 6);
    }
}
