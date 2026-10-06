//! Sibling upper bounds: skip DP for alleles of the same gene that cannot
//! enter the top-N local scores.
//!
//! After the representative allele \(R\) of a gene is scored, every same-length
//! sibling \(A\) satisfies
//! \[
//! S^*(A) \le S^*(R) + \Delta(R,A)
//! \]
//! where \(\Delta\) is the sum over reference positions of
//! \(\max_q\bigl(s(A[j],q) - s(R[j],q)\bigr)\), taken from the substitution
//! profiles. The proof maps the optimal path of \(A\) onto \(R\) at the same
//! coordinates (valid when \(|A| = |R|\)); gap columns are unchanged and each
//! substitution column gains at most that per-position maximum.
//!
//! Alleles of different length are not bounded and are always scored.

use crate::alphabet::ALPHABET_SIZE;
use crate::germline::{AlleleId, GeneId, GermlineSet};
use crate::scoring::{Score, SubstitutionProfile};
use crate::segment::PreparedSegment;

/// Precomputed per-gene allele lists for family-delta pruning.
#[derive(Clone, Debug, Default)]
pub struct FamilyBoundTable {
    /// Gene of each allele id.
    allele_gene: Vec<GeneId>,
    /// Alleles of each gene, in allele-id order.
    gene_alleles: Vec<Vec<AlleleId>>,
}

impl FamilyBoundTable {
    pub fn build(set: &GermlineSet) -> Self {
        let allele_gene: Vec<GeneId> = set.alleles().iter().map(|a| a.gene).collect();
        let mut gene_alleles = vec![Vec::new(); set.gene_count()];
        for (id, gene) in allele_gene.iter().copied().enumerate() {
            gene_alleles[gene as usize].push(id as AlleleId);
        }
        Self {
            allele_gene,
            gene_alleles,
        }
    }

    pub fn gene(&self, allele: AlleleId) -> GeneId {
        self.allele_gene[allele as usize]
    }

    pub fn alleles_of(&self, gene: GeneId) -> &[AlleleId] {
        &self.gene_alleles[gene as usize]
    }

    /// Upper bound on the sibling's local score given the representative's.
    /// `None` when the alleles differ in length (bound does not apply).
    pub fn upper_bound(
        &self,
        segment: &PreparedSegment,
        representative: AlleleId,
        sibling: AlleleId,
        representative_score: Score,
    ) -> Option<Score> {
        let delta = sibling_delta(
            segment.profile(representative),
            segment.profile(sibling),
        )?;
        Some(representative_score.saturating_add(delta))
    }
}

/// \(\Delta(R,A) = \sum_j \max\bigl(0,\ \max_q (s_A(q,j) - s_R(q,j))\bigr)\).
///
/// `None` when the profiles have different lengths.
pub fn sibling_delta(
    representative: &SubstitutionProfile,
    sibling: &SubstitutionProfile,
) -> Option<Score> {
    if representative.len() != sibling.len() {
        return None;
    }
    let mut delta = 0;
    for j in 0..representative.len() {
        let mut best = 0;
        for code in 0..ALPHABET_SIZE {
            let gain =
                sibling.at(code as u8, j) - representative.at(code as u8, j);
            if gain > best {
                best = gain;
            }
        }
        delta += best;
    }
    Some(delta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Aligner, BandedGotoh};
    use crate::alphabet::encode;
    use crate::band::Band;
    use crate::germline::SegmentKind;
    use crate::scoring::{AffineScoring, GapCosts};

    fn profiles(a: &[u8], b: &[u8]) -> (SubstitutionProfile, SubstitutionProfile) {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        (
            SubstitutionProfile::build(&scoring, &encode(a)),
            SubstitutionProfile::build(&scoring, &encode(b)),
        )
    }

    fn score_against(query: &[u8], reference: &[u8]) -> Score {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let profile = SubstitutionProfile::build(&scoring, reference);
        let mut ws = Default::default();
        BandedGotoh
            .score(
                query,
                &profile,
                GapCosts::new(-9, -1),
                Band::Full,
                &mut ws,
            )
            .map(|hit| hit.score)
            .unwrap_or(0)
    }

    #[test]
    fn identical_alleles_have_zero_delta() {
        let (r, a) = profiles(b"ACGTACGT", b"ACGTACGT");
        assert_eq!(sibling_delta(&r, &a), Some(0));
    }

    #[test]
    fn one_snp_costs_match_minus_mismatch() {
        let (r, a) = profiles(b"ACGTACGT", b"AGGTACGT");
        assert_eq!(sibling_delta(&r, &a), Some(5)); // 2 - (-3)
    }

    #[test]
    fn different_lengths_are_unbounded() {
        let (r, a) = profiles(b"ACGTACGT", b"ACGTACG");
        assert_eq!(sibling_delta(&r, &a), None);
    }

    #[test]
    fn table_groups_alleles_by_gene() {
        let mut set = GermlineSet::new(SegmentKind::V);
        set.push("IGHV1-1*01", b"ACGTACGTACGTACGTACGTACGT", None)
            .unwrap();
        set.push("IGHV1-1*02", b"AGGTACGTACGTACGTACGTACGT", None)
            .unwrap();
        set.push("IGHV1-2*01", b"TTTTACGTACGTACGTACGTACGT", None)
            .unwrap();
        let table = FamilyBoundTable::build(&set);
        assert_eq!(table.alleles_of(0), &[0, 1]);
        assert_eq!(table.alleles_of(1), &[2]);
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let segment = PreparedSegment::new(set, &scoring, scoring.gaps);
        assert_eq!(table.upper_bound(&segment, 0, 1, 100), Some(105));
        // Same length, three SNPs in the prefix (position 3 is T in both):
        // bound is available; the stage only applies it inside one gene.
        assert_eq!(table.upper_bound(&segment, 0, 2, 100), Some(115));
    }

    #[test]
    fn bound_never_underestimates_local_score() {
        // Same length; two substitutions (positions 20 and 23).
        let rep = encode(b"CAGGTCACCTTGAAGGAGTCTGGTCC");
        let sib = encode(b"CAGGTCACCTTGAAGGAGTCCGGCCC");
        let (rep_profile, sib_profile) = {
            let scoring = AffineScoring::new(2, -3, -9, -1);
            (
                SubstitutionProfile::build(&scoring, &rep),
                SubstitutionProfile::build(&scoring, &sib),
            )
        };
        let delta = sibling_delta(&rep_profile, &sib_profile).unwrap();
        assert_eq!(delta, 10);

        // Queries close to the sibling, the representative, and a decoy.
        let queries = [
            sib.clone(),
            rep.clone(),
            encode(b"CAGGTCACCTTGAAGGAGTCTGGTCA"),
            encode(b"NNNNNCACCTTGAAGGAGTCTGGTCC"),
        ];
        for query in &queries {
            let s_rep = score_against(query, &rep);
            let s_sib = score_against(query, &sib);
            assert!(
                s_sib <= s_rep + delta,
                "query score sibling {s_sib} > rep {s_rep} + delta {delta}"
            );
        }
    }
}
