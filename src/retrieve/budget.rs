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

/// Swig's AER V budget.
///
/// Always keep `base` genes (3, which is `max(top_v, 3)`). Then keep each
/// further gene whose votes sit within `max(absolute, relative · leader)` of
/// the leader, up to `max` (16). J uses `base == max == 3`, so it does not expand.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AmbiguityBudget {
    pub base: usize,
    pub max: usize,
    pub relative: f64,
    pub absolute: u32,
}

impl Default for AmbiguityBudget {
    fn default() -> Self {
        Self {
            base: 3,
            max: 16,
            relative: 0.05,
            absolute: 8,
        }
    }
}

impl AmbiguityBudget {
    fn margin(&self, leader: u32) -> u32 {
        let relative = (f64::from(leader) * self.relative) as u32;
        self.absolute.max(relative)
    }
}

impl CandidateBudget for AmbiguityBudget {
    fn genes_to_keep(&self, ranked_gene_votes: &[u32]) -> usize {
        let Some(&leader) = ranked_gene_votes.first() else {
            return 0;
        };
        let mut keep = self.base.min(ranked_gene_votes.len());
        let margin = self.margin(leader);
        while keep < self.max.min(ranked_gene_votes.len())
            && leader.saturating_sub(ranked_gene_votes[keep]) <= margin
        {
            keep += 1;
        }
        keep
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ambiguity_widens_the_pool() {
        let budget = AmbiguityBudget::default();
        let clear: Vec<u32> = [400, 200]
            .into_iter()
            .chain(std::iter::repeat(10).take(20))
            .collect();
        let near: Vec<u32> = [400, 395, 390]
            .into_iter()
            .chain(std::iter::repeat(10).take(20))
            .collect();
        let deep = vec![100u32; 20];
        assert_eq!(budget.genes_to_keep(&clear), 3);
        assert_eq!(budget.genes_to_keep(&near), 3);
        assert_eq!(budget.genes_to_keep(&deep), 16);
        assert_eq!(
            AmbiguityBudget {
                base: 3,
                max: 3,
                ..budget
            }
            .genes_to_keep(&deep),
            3
        );
        assert_eq!(FixedBudget(3).genes_to_keep(&[5]), 1);
    }
}
