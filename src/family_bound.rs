//! Sibling upper bounds: skip DP for alleles of the same gene that cannot
//! enter the top-N local scores.
//!
//! After the representative allele \(R\) of a gene is scored in band \(B\),
//! every same-length sibling \(A\) aligned in the *same* band satisfies
//! \[
//! S^*_B(A) \le S^*_B(R) + \Delta(R,A)
//! \]
//! where \(\Delta\) is the sum over reference positions of
//! \(\max\bigl(0, \max_q (s(A[j],q) - s(R[j],q))\bigr)\), taken from the
//! substitution profiles. The proof maps the optimal path of \(A\) onto \(R\)
//! cell for cell: the path stays inside \(B\) and inside \(R\)'s matrix
//! (\(|A| = |R|\)), gap columns score alike, and each substitution column
//! gains at most that per-position maximum. A sibling in another band (another
//! gene frame, or another offset in the frame) or of another length is not
//! bounded and is always scored.
//!
//! \(\Delta\) depends only on the two profiles, so it is computed once per
//! same-length allele pair of each gene when the stage is built.

use crate::alphabet::ALPHABET_SIZE;
use crate::germline::{AlleleId, GeneId};
use crate::scoring::{Score, SubstitutionProfile};
use crate::segment::PreparedSegment;

/// \(\Delta(R, A)\) of every same-length allele pair within each gene.
#[derive(Clone, Debug, Default)]
pub struct FamilyBoundTable {
    /// Gene of each allele id.
    allele_gene: Vec<GeneId>,
    /// Position of each allele in its gene's allele list.
    rank: Vec<u32>,
    /// Per gene: `delta[r · n + a]` for representative rank `r` and sibling
    /// rank `a` (`n` alleles in the gene), `None` across lengths.
    delta: Vec<Vec<Option<Score>>>,
}

impl FamilyBoundTable {
    pub fn build(segment: &PreparedSegment) -> Self {
        let set = segment.set();
        let allele_gene: Vec<GeneId> = set.alleles().iter().map(|a| a.gene).collect();
        let mut members: Vec<Vec<AlleleId>> = vec![Vec::new(); set.gene_count()];
        let mut rank = vec![0; allele_gene.len()];
        for (id, &gene) in allele_gene.iter().enumerate() {
            let list = &mut members[gene as usize];
            rank[id] = u32::try_from(list.len()).expect("allele count fits u32");
            list.push(AlleleId::try_from(id).expect("allele count fits u32"));
        }
        let delta = members
            .iter()
            .map(|alleles| {
                alleles
                    .iter()
                    .flat_map(|&r| {
                        alleles
                            .iter()
                            .map(move |&a| sibling_delta(segment.profile(r), segment.profile(a)))
                    })
                    .collect()
            })
            .collect();
        Self {
            allele_gene,
            rank,
            delta,
        }
    }

    #[inline]
    pub fn gene(&self, allele: AlleleId) -> GeneId {
        self.allele_gene[allele as usize]
    }

    /// Upper bound on the sibling's local score given the representative's,
    /// both aligned in the same band. `None` when the alleles belong to
    /// different genes or differ in length (the bound does not apply).
    pub fn upper_bound(
        &self,
        representative: AlleleId,
        sibling: AlleleId,
        representative_score: Score,
    ) -> Option<Score> {
        let gene = self.gene(representative);
        if gene != self.gene(sibling) {
            return None;
        }
        let deltas = &self.delta[gene as usize];
        let n = deltas.len().isqrt();
        let (r, a) = (
            self.rank[representative as usize] as usize,
            self.rank[sibling as usize] as usize,
        );
        deltas[r * n + a].map(|delta| representative_score.saturating_add(delta))
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
    let gain = |code: u8| {
        sibling
            .row(code)
            .iter()
            .zip(representative.row(code))
            .map(|(a, r)| a - r)
    };
    let mut best = vec![0; representative.len()];
    for code in 0..ALPHABET_SIZE as u8 {
        for (slot, g) in best.iter_mut().zip(gain(code)) {
            *slot = (*slot).max(g);
        }
    }
    Some(best.iter().sum())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Aligner, BandedGotoh};
    use crate::alphabet::encode;
    use crate::band::Band;
    use crate::germline::{GermlineSet, SegmentKind};
    use crate::scoring::AffineScoring;

    const SCORING: AffineScoring = AffineScoring::new(2, -3, -9, -1);

    fn profile(sequence: &[u8]) -> SubstitutionProfile {
        SubstitutionProfile::build(&SCORING, &encode(sequence))
    }

    fn score_against(query: &[u8], reference: &[u8], band: Band) -> Score {
        let profile = SubstitutionProfile::build(&SCORING, reference);
        BandedGotoh
            .score(query, &profile, SCORING.gaps, band, &mut Default::default())
            .map_or(0, |hit| hit.score)
    }

    #[test]
    fn deltas_count_profile_gains() {
        let delta = |a: &[u8], b: &[u8]| sibling_delta(&profile(a), &profile(b));
        assert_eq!(delta(b"ACGTACGT", b"ACGTACGT"), Some(0));
        assert_eq!(delta(b"ACGTACGT", b"AGGTACGT"), Some(5)); // 2 − (−3)
        assert_eq!(delta(b"ACGTACGT", b"ACGTACG"), None);
    }

    #[test]
    fn table_bounds_same_gene_same_length_pairs_only() {
        let mut set = GermlineSet::new(SegmentKind::V);
        set.push("IGHV1-1*01", b"ACGTACGTACGTACGTACGTACGT", None)
            .unwrap();
        set.push("IGHV1-2*01", b"TTTTACGTACGTACGTACGTACGT", None)
            .unwrap();
        set.push("IGHV1-1*02", b"AGGTACGTACGTACGTACGTACGT", None)
            .unwrap();
        set.push("IGHV1-1*03", b"ACGTACGTACGTACGTACGTACG", None)
            .unwrap();
        let segment = PreparedSegment::new(set, &SCORING, SCORING.gaps);
        let table = FamilyBoundTable::build(&segment);
        assert_eq!(table.upper_bound(0, 2, 100), Some(105));
        assert_eq!(table.upper_bound(2, 0, 100), Some(105));
        assert_eq!(table.upper_bound(0, 0, 100), Some(100));
        // Another gene, another length.
        assert_eq!(table.upper_bound(0, 1, 100), None);
        assert_eq!(table.upper_bound(0, 3, 100), None);
    }

    /// The bound holds for every query and every band both alleles share.
    #[test]
    fn bound_never_underestimates_local_score() {
        let rep = encode(b"CAGGTCACCTTGAAGGAGTCTGGTCCTGCCTGGTGAAACCCACACAGACC");
        let mut sib = rep.clone();
        sib[20] = 1;
        sib[23] = 1;
        sib[41] = 3;
        let delta = sibling_delta(
            &SubstitutionProfile::build(&SCORING, &rep),
            &SubstitutionProfile::build(&SCORING, &sib),
        )
        .unwrap();
        let mut words = crate::testing::rng(3);
        let mut next = || words() as usize;
        for _ in 0..300 {
            // A mutated, shifted, indel-carrying copy of either allele.
            let source = if next() % 2 == 0 { &rep } else { &sib };
            let mut query: Vec<u8> = (0..next() % 12).map(|_| (next() % 4) as u8).collect();
            let lead = query.len() as i32;
            query.extend(source.iter().map(|&c| {
                if next() % 8 == 0 {
                    (next() % 5) as u8
                } else {
                    c
                }
            }));
            if next() % 2 == 0 {
                let at = next() % query.len();
                query.insert(at, (next() % 4) as u8);
            }
            let band = match next() % 3 {
                0 => Band::Full,
                1 => Band::around(lead, 2),
                _ => Band::around(lead + (next() % 7) as i32 - 3, 8),
            };
            let (s_rep, s_sib) = (
                score_against(&query, &rep, band),
                score_against(&query, &sib, band),
            );
            assert!(
                s_sib <= s_rep + delta,
                "{band:?}: sibling {s_sib} > {s_rep} + {delta}"
            );
        }
    }
}
