//! Decoding the D allele from the posterior instead of the Viterbi path.
//!
//! The call is a mode of `P(path, cuts | read)`: the best single explanation
//! of the read. Its D allele need not be the most probable one. A short exact
//! D remnant and a longer one with a mismatch, or one allele placed at
//! several offsets, can each lose the maximum and still hold most of the
//! summed mass. Under 0/1 loss on the D allele the Bayes act is the allele of
//! largest posterior marginal,
//!
//! ```text
//! a* = argmax_a P(a),   P(a) = Σ_{D placements p ∋ a} P(p on the path)
//! ```
//!
//! ([`super::Marginals::placement`]; an allele merged into several placements
//! sums them). D presence stays the Viterbi decision: the presence marginal
//! sums the look-elsewhere mass of every chance hit, and taken as a decision
//! it gives a D to a third of D-free reads (docs/ASSESSMENT.md §5.3). So only
//! a call with exactly one D is decoded, and a tie within a relative `10⁻⁹`
//! keeps the call's allele.
//!
//! The reported alignment is then the best path, under the call's own model,
//! with the call's [`Topology`] whose D contains `a*`: the same DAG over the
//! V and J placements and the D placements of `a*`, with the grammar
//! restricted by [`Within`].

use super::Marginals;
use super::model::Grammar;
use crate::germline::{AlleleId, SegmentKind};
use crate::placement::Placement;

/// How the reported D allele is chosen.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DAlleleDecoding {
    /// The call's own: the allele of the Viterbi path.
    #[default]
    Map,
    /// The allele of largest posterior marginal, given the call's D presence
    /// and topology. Needs posteriors: without them (`--no-posterior`) the
    /// D allele is the call's.
    MaxMarginal,
}

/// A path's shape: its first and last segment kinds and its D count.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Topology {
    pub first: SegmentKind,
    pub last: SegmentKind,
    pub d_count: u8,
}

/// The paths of `grammar` with one [`Topology`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct Within<'a, G> {
    pub grammar: &'a G,
    pub topology: Topology,
}

impl<G: Grammar> Grammar for Within<'_, G> {
    fn can_start(&self, kind: SegmentKind) -> bool {
        self.grammar.can_start(kind) && kind == self.topology.first
    }

    fn can_follow(&self, prev: SegmentKind, next: SegmentKind, d_used: u8) -> bool {
        self.grammar.can_follow(prev, next, d_used)
            && (next != SegmentKind::D || d_used < self.topology.d_count)
    }

    fn can_end(&self, first: SegmentKind, last: SegmentKind, d_count: u8) -> bool {
        self.grammar.can_end(first, last, d_count)
            && last == self.topology.last
            && d_count == self.topology.d_count
    }
}

/// The allele of largest posterior marginal among the D placements, when
/// it beats every allele of `called`, the call's only D, by more than a
/// relative `10⁻⁹` (ties go to the call, then to the smaller id).
pub(crate) fn max_marginal_d(
    placements: &[Placement],
    called: &Placement,
    marginals: &Marginals,
) -> Option<AlleleId> {
    let d_placements = || {
        placements
            .iter()
            .zip(&marginals.placement)
            .filter(|(p, _)| p.kind == SegmentKind::D)
    };
    let mass = |allele: AlleleId| -> f64 {
        d_placements()
            .filter(|(p, _)| p.alleles.contains(&allele))
            .map(|(_, &m)| m)
            .sum()
    };
    let called_mass = called.alleles.iter().map(|&a| mass(a)).fold(0.0, f64::max);
    let mut alleles: Vec<AlleleId> = d_placements()
        .filter(|&(_, &m)| m > 0.0)
        .flat_map(|(p, _)| p.alleles.iter().copied())
        .collect();
    alleles.sort_unstable();
    alleles.dedup();
    let mut best: Option<(AlleleId, f64)> = None;
    for allele in alleles {
        let m = mass(allele);
        if m > called_mass * (1.0 + 1e-9) && best.is_none_or(|(_, b)| m > b) {
            best = Some((allele, m));
        }
    }
    best.map(|(allele, _)| allele)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adjudicate::VdjGrammar;
    use crate::testing::placement;

    fn marginals(placement: Vec<f64>) -> Marginals {
        Marginals {
            placement,
            d_count: vec![0.0, 1.0],
            log_partition: 0.0,
        }
    }

    fn d(alleles: &[AlleleId]) -> Placement {
        let mut p = placement(SegmentKind::D, &[0, 1, 2, 3], &[0, 1, 2, 3], 0);
        p.alleles = alleles.to_vec();
        p
    }

    #[test]
    fn the_largest_summed_allele_mass_wins() {
        let placements = [d(&[0]), d(&[1]), d(&[1, 2]), d(&[3])];
        let called = &placements[0];
        // Allele 1 sums two placements (0.25 + 0.2) and beats 0 (0.4).
        let m = marginals(vec![0.4, 0.25, 0.2, 0.1]);
        assert_eq!(max_marginal_d(&placements, called, &m), Some(1));
        // The called allele keeps a tie within 10⁻⁹.
        let m = marginals(vec![0.45, 0.25, 0.2, 0.1]);
        assert_eq!(max_marginal_d(&placements, called, &m), None);
        let m = marginals(vec![0.45 * (1.0 - 1e-12), 0.25, 0.2, 0.1]);
        assert_eq!(max_marginal_d(&placements, called, &m), None);
        // Among others, the smaller id wins a tie.
        let m = marginals(vec![0.1, 0.0, 0.3, 0.3]);
        assert_eq!(max_marginal_d(&placements, called, &m), Some(1));
        // A merged placement holding the winner is the call's.
        let m = marginals(vec![0.1, 0.0, 0.5, 0.0]);
        assert_eq!(max_marginal_d(&placements, &placements[2], &m), None);
    }

    #[test]
    fn within_keeps_one_topology() {
        use SegmentKind::{D, J, V};
        let vdj = Within {
            grammar: &VdjGrammar {
                max_d: 2,
                ..VdjGrammar::VDJ
            },
            topology: Topology {
                first: V,
                last: J,
                d_count: 1,
            },
        };
        assert!(vdj.can_start(V) && !vdj.can_start(D) && !vdj.can_start(J));
        assert!(vdj.can_follow(V, D, 0) && !vdj.can_follow(D, D, 1));
        assert!(vdj.can_end(V, J, 1) && !vdj.can_end(V, J, 0) && !vdj.can_end(V, D, 1));
    }
}
