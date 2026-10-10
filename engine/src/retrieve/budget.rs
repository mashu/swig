//! How many genes survive retrieval.

/// Decides how many top-ranked genes to keep, given gene votes sorted
/// descending and the votes one seed casts (`vote_unit`).
pub trait CandidateBudget: Send + Sync {
    fn genes_to_keep(&self, ranked_gene_votes: &[u32], vote_unit: u32) -> usize;
}

/// Always keep the same number of genes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FixedBudget(pub usize);

impl CandidateBudget for FixedBudget {
    fn genes_to_keep(&self, ranked_gene_votes: &[u32], _vote_unit: u32) -> usize {
        self.0.min(ranked_gene_votes.len())
    }
}

/// Keep the `base` best genes, then every further gene whose votes are
/// within the leader's sampling noise, up to `max`.
///
/// A gene's votes count the seeds it shares with the read, so a count of `n`
/// seeds fluctuates by about `√n` (Poisson). Two genes near the leader's
/// count `L` (in votes, `unit` votes per seed) differ by chance by about
/// `unit · √(2L / unit) = √(2 · unit · L)` votes. The margin is `z` such
/// standard deviations: relative to the leader, but wider in proportion where
/// counts are small. A proportional margin is too narrow for weak (heavily
/// mutated) reads (docs/ASSESSMENT.md, item 18). When 7-mers vote too
/// (1 vote each), the counts mix weights and the margin, computed at
/// `unit = 4`, is between one and two times `z` actual standard deviations.
/// Used for V; J is searched exhaustively.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AmbiguityBudget {
    pub base: usize,
    pub max: usize,
    /// Standard deviations of the vote difference within which a gene is kept.
    pub z: f64,
}

impl Default for AmbiguityBudget {
    fn default() -> Self {
        Self {
            base: 3,
            max: 16,
            z: 2.0,
        }
    }
}

impl AmbiguityBudget {
    /// Votes a gene may trail the leader by and still be kept.
    pub fn margin(&self, leader: u32, vote_unit: u32) -> u32 {
        (self.z * (2.0 * f64::from(vote_unit.max(1)) * f64::from(leader)).sqrt()) as u32
    }
}

impl CandidateBudget for AmbiguityBudget {
    fn genes_to_keep(&self, ranked_gene_votes: &[u32], vote_unit: u32) -> usize {
        let Some(&leader) = ranked_gene_votes.first() else {
            return 0;
        };
        let mut keep = self.base.min(ranked_gene_votes.len());
        let margin = self.margin(leader, vote_unit);
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
        // Leader 400 votes at 4 votes per seed: margin 2·√3200 = 113.
        assert_eq!(budget.margin(400, 4), 113);
        let clear: Vec<u32> = [400, 200]
            .into_iter()
            .chain(std::iter::repeat_n(10, 20))
            .collect();
        let near: Vec<u32> = [400, 395, 390, 300, 280]
            .into_iter()
            .chain(std::iter::repeat_n(10, 20))
            .collect();
        let deep = vec![100u32; 20];
        assert_eq!(budget.genes_to_keep(&clear, 4), 3);
        assert_eq!(budget.genes_to_keep(&near, 4), 4);
        assert_eq!(budget.genes_to_keep(&deep, 4), 16);
        // The margin grows like the square root of the leader: relatively
        // wider for a weak (mutated) read.
        assert!(budget.margin(25, 4) as f64 / 25.0 > budget.margin(400, 4) as f64 / 400.0);
        assert_eq!(
            AmbiguityBudget {
                base: 3,
                max: 3,
                ..budget
            }
            .genes_to_keep(&deep, 4),
            3
        );
        assert_eq!(FixedBudget(3).genes_to_keep(&[5], 4), 1);
    }
}
