//! How many genes survive retrieval.

/// Decides how many top-ranked genes to keep, given gene votes sorted descending.
pub trait CandidateBudget: Send + Sync {
    fn genes_to_keep(&self, ranked_gene_votes: &[u32]) -> usize;
}

/// Always keep the same number of genes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FixedBudget(pub usize);

impl CandidateBudget for FixedBudget {
    fn genes_to_keep(&self, ranked_gene_votes: &[u32]) -> usize {
        self.0.min(ranked_gene_votes.len())
    }
}

/// The AER rule: refine deeper when the vote ranking is ambiguous.
///
/// If the runner-up gene is within `relative` (fraction of the leader) or
/// `absolute` weighted votes of the leader, keep up to `max` genes, otherwise
/// keep `base`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AmbiguityBudget {
    pub base: usize,
    pub max: usize,
    pub relative: f64,
    pub absolute: u32,
}

impl Default for AmbiguityBudget {
    fn default() -> Self {
        Self { base: 4, max: 16, relative: 0.05, absolute: 8 }
    }
}

impl AmbiguityBudget {
    fn ambiguous(&self, leader: u32, other: u32) -> bool {
        let gap = leader.saturating_sub(other);
        gap <= self.absolute || f64::from(gap) <= self.relative * f64::from(leader)
    }
}

impl CandidateBudget for AmbiguityBudget {
    fn genes_to_keep(&self, ranked_gene_votes: &[u32]) -> usize {
        let keep = match ranked_gene_votes {
            [leader, second, ..] if self.ambiguous(*leader, *second) => self.max,
            _ => self.base,
        };
        keep.min(ranked_gene_votes.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ambiguity_widens_the_pool() {
        let budget = AmbiguityBudget::default();
        let clear: Vec<u32> = [400, 200].into_iter().chain(std::iter::repeat(10).take(20)).collect();
        let tied: Vec<u32> = [400, 395].into_iter().chain(std::iter::repeat(10).take(20)).collect();
        assert_eq!(budget.genes_to_keep(&clear), 4);
        assert_eq!(budget.genes_to_keep(&tied), 16);
        assert_eq!(FixedBudget(3).genes_to_keep(&[5]), 1);
    }
}
