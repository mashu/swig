//! Karlin–Altschul statistics: the E-value gates of V and J hits
//! ([`crate::segment::Significance`]), the λ that turns scores into nats for
//! the D prior and the posteriors, and the reported `*_support` E-values.
//!
//! `E = K m′ n′ e^{−λS}` with the BLAST length adjustment
//! `ℓ = ln(K m n) / H`, `m′ = m − ℓ`, `n′ = n − N ℓ`.
//!
//! [`KarlinAltschul::ungapped`] computes λ, H and K exactly for a match /
//! mismatch tuple under uniform base composition: λ by bisection on
//! `Σ p_s e^{λs} = 1`, H as the tilted mean, and K from the Karlin–Altschul
//! series `σ = Σ_k k⁻¹ E[min(e^{λ S_k}, 1)]`. (It reproduces NCBI's ungapped
//! blastn table, e.g. reward 1 / penalty −3 → λ 1.374, K 0.711, H 1.31.)
//! Gapped alignments score slightly higher than ungapped ones, so these
//! E-values are mildly anti-conservative; for calibrated values fit λ and K
//! by simulation and construct the struct directly.
//!
//! **The σ series.** `k` columns with `a` matches score
//! `S_k = a·r − (k − a)·q` (reward `r`, penalty `−q`), `a ~ Binomial(k, ¼)`,
//! so `S_k ≥ 0` exactly from `a* = ⌈k·q / (r + q)⌉` matches on, and
//!
//! ```text
//! E[min(e^{λS_k}, 1)] = P(a ≥ a*) + E[e^{λS_k}; a < a*] = P(a ≥ a*) + P′(a < a*)
//! ```
//!
//! where `P′` is the binomial with the tilted weights `¼·e^{λr}` and
//! `¾·e^{−λq}` (they sum to 1 at λ). Each term is two binomial tails,
//! advanced one column at a time over the probabilities that are not
//! negligible: no exponential per atom, and a cost independent of the
//! scores' size. The series stops at the first term below 10⁻¹³, or after
//! [`SIGMA_TERMS`] terms: a tuple whose expected score per column is small
//! next to its scores converges slowly (5/−2 takes about 5,800 terms) and
//! one slower still keeps its truncated sum.

use crate::scoring::{AffineScoring, Score};

/// Most terms of the σ series (see the module documentation).
pub const SIGMA_TERMS: usize = 20_000;

/// Binomial probabilities of the match count after some columns, over the
/// counts whose probability is not negligible; advanced one column at a
/// time.
struct Binomial {
    /// `pmf[i]` is the probability of `first + i` matches.
    pmf: Vec<f64>,
    first: usize,
    matched: f64,
    mismatched: f64,
}

impl Binomial {
    /// Probabilities this small are dropped from either end.
    const TINY: f64 = 1e-280;

    /// Zero columns, with per-column weights `matched` and `mismatched`.
    fn new(matched: f64, mismatched: f64) -> Self {
        Self {
            pmf: vec![1.0],
            first: 0,
            matched,
            mismatched,
        }
    }

    /// Add one column.
    fn step(&mut self) {
        let Self {
            pmf,
            first,
            matched,
            mismatched,
        } = self;
        pmf.push(0.0);
        for i in (1..pmf.len()).rev() {
            pmf[i] = pmf[i] * *mismatched + pmf[i - 1] * *matched;
        }
        pmf[0] *= *mismatched;
        let low = pmf
            .iter()
            .position(|&p| p >= Self::TINY)
            .unwrap_or(pmf.len());
        pmf.drain(..low);
        *first += low;
        while pmf.last().is_some_and(|&p| p < Self::TINY) {
            pmf.pop();
        }
    }

    /// Index into `pmf` of `count` matches, clamped to its range.
    fn at(&self, count: usize) -> usize {
        count.saturating_sub(self.first).min(self.pmf.len())
    }

    /// Probability of at least `count` matches.
    fn at_least(&self, count: usize) -> f64 {
        self.pmf[self.at(count)..].iter().sum()
    }

    /// Probability of fewer than `count` matches.
    fn below(&self, count: usize) -> f64 {
        self.pmf[..self.at(count)].iter().sum()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KarlinAltschul {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
}

fn gcd(a: u32, b: u32) -> u32 {
    if b == 0 { a } else { gcd(b, a % b) }
}

impl KarlinAltschul {
    /// Exact ungapped parameters for uniform `A C G T`; `None` when the tuple
    /// has no local-alignment statistics (non-positive match or non-negative
    /// expected score).
    pub fn ungapped(scoring: &AffineScoring) -> Option<Self> {
        let (reward, penalty) = (scoring.matched, scoring.mismatch);
        let (p_match, p_mismatch) = (0.25, 0.75);
        let expected = p_match * f64::from(reward) + p_mismatch * f64::from(penalty);
        if reward <= 0 || expected >= 0.0 {
            return None;
        }
        let mgf = |lambda: f64| {
            p_match * (lambda * f64::from(reward)).exp()
                + p_mismatch * (lambda * f64::from(penalty)).exp()
        };
        let (mut lo, mut hi) = (1e-9, 1.0);
        while mgf(hi) < 1.0 {
            hi *= 2.0;
        }
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if mgf(mid) > 1.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let lambda = 0.5 * (lo + hi);
        let tilted_mean = p_match * f64::from(reward) * (lambda * f64::from(reward)).exp()
            + p_mismatch * f64::from(penalty) * (lambda * f64::from(penalty)).exp();
        let h = lambda * tilted_mean;

        // σ series: E[min(e^{λS_k}, 1)] as two binomial tails (module docs).
        let (r, q) = (
            u64::from(reward.unsigned_abs()),
            u64::from(penalty.unsigned_abs()),
        );
        let mut plain = Binomial::new(p_match, p_mismatch);
        let mut tilted = Binomial::new(
            p_match * (lambda * f64::from(reward)).exp(),
            p_mismatch * (lambda * f64::from(penalty)).exp(),
        );
        let mut sigma = 0.0;
        for k in 1..=SIGMA_TERMS {
            plain.step();
            tilted.step();
            let first = usize::try_from((k as u64 * q).div_ceil(r + q)).unwrap_or(usize::MAX);
            let term = (tilted.below(first) + plain.at_least(first)) / k as f64;
            sigma += term;
            if term < 1e-13 {
                break;
            }
        }
        let delta = f64::from(gcd(reward as u32, penalty.unsigned_abs()));
        let k = lambda * delta * (-2.0 * sigma).exp() / (h * (1.0 - (-lambda * delta).exp()));
        Some(Self { lambda, k, h })
    }

    /// Expectation of chance local alignments scoring at least `score` between a
    /// query of `query_len` and a database of `db_len` residues in `db_seqs` sequences.
    pub fn evalue(&self, score: Score, query_len: usize, db_len: usize, db_seqs: usize) -> f64 {
        let (m, n) = (query_len as f64, db_len as f64);
        let ell = ((self.k * m * n).ln() / self.h).max(0.0);
        let m_eff = (m - ell).max(1.0);
        let n_eff = (n - db_seqs as f64 * ell).max(1.0);
        self.k * m_eff * n_eff * (-self.lambda * f64::from(score)).exp()
    }

    /// Normalised bit score.
    pub fn bit_score(&self, score: Score) -> f64 {
        (self.lambda * f64::from(score) - self.k.ln()) / std::f64::consts::LN_2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn reproduces_known_ungapped_parameters() {
        let pm1 = KarlinAltschul::ungapped(&AffineScoring::new(1, -1, -5, -1)).unwrap();
        assert!(close(pm1.lambda, 3f64.ln(), 1e-9));
        assert!(
            close(pm1.k, 1.0 / 3.0, 1e-6),
            "±1 walk has K = (q−p)²/q = 1/3, got {}",
            pm1.k
        );
        let ncbi = KarlinAltschul::ungapped(&AffineScoring::new(1, -3, -5, -2)).unwrap();
        assert!(
            close(ncbi.lambda, 1.374, 1e-3)
                && close(ncbi.k, 0.711, 1e-3)
                && close(ncbi.h, 1.31, 1e-2)
        );
        assert!(KarlinAltschul::ungapped(&AffineScoring::new(2, 0, -5, -1)).is_none());
    }

    /// The σ series of the dense score distribution, for reference, if it
    /// converges within `terms` terms.
    fn dense_sigma(reward: Score, penalty: Score, lambda: f64, terms: usize) -> Option<f64> {
        let (r, q) = (reward as usize, penalty.unsigned_abs() as usize);
        let mut dist = vec![1.0f64];
        let mut sigma = 0.0;
        for k in 1..=terms {
            let mut next = vec![0.0f64; dist.len() + r + q];
            for (i, &pr) in dist.iter().enumerate() {
                next[i + r + q] += pr * 0.25;
                next[i] += pr * 0.75;
            }
            dist = next;
            let offset = (k * q) as f64;
            let term: f64 = (0..)
                .zip(&dist)
                .map(|(i, &pr): (usize, _)| {
                    let s = i as f64 - offset;
                    pr * if s >= 0.0 { 1.0 } else { (lambda * s).exp() }
                })
                .sum::<f64>()
                / k as f64;
            sigma += term;
            if term < 1e-13 {
                return Some(sigma);
            }
        }
        None
    }

    /// The binomial-tail series sums what the dense score distribution
    /// does, term for term, whatever the tuple.
    #[test]
    fn sigma_series_matches_the_dense_score_distribution() {
        let mut compared = 0;
        for reward in 1..=6 {
            for penalty in -12..=-1 {
                let scoring = AffineScoring::new(reward, penalty, -5, -1);
                let Some(ka) = KarlinAltschul::ungapped(&scoring) else {
                    continue;
                };
                let Some(sigma) = dense_sigma(reward, penalty, ka.lambda, 1500) else {
                    continue;
                };
                compared += 1;
                let delta = f64::from(gcd(reward as u32, penalty.unsigned_abs()));
                let (lambda, h) = (ka.lambda, ka.h);
                let k =
                    lambda * delta * (-2.0 * sigma).exp() / (h * (1.0 - (-lambda * delta).exp()));
                assert!(
                    close(ka.k, k, 1e-12 * k),
                    "{reward}/{penalty}: {} vs {k}",
                    ka.k
                );
            }
        }
        assert!(compared >= 50, "{compared} tuples compared");
    }

    /// A slowly converging tuple runs past 500 terms (K was 5 % high at
    /// 500), and an enormous penalty costs no more than a small one.
    #[test]
    fn sigma_series_converges_for_any_penalty() {
        let slow = KarlinAltschul::ungapped(&AffineScoring::new(5, -2, -5, -1)).unwrap();
        assert!(close(slow.k, 0.011_370_919_8, 1e-10), "{}", slow.k);
        let started = std::time::Instant::now();
        let steep = KarlinAltschul::ungapped(&AffineScoring::new(2, Score::MIN, -5, -1)).unwrap();
        assert!(started.elapsed().as_secs_f64() < 0.5);
        assert!(close(steep.lambda, 4f64.ln() / 2.0, 1e-9));
        assert!(steep.k.is_finite() && steep.k > 0.0);
    }

    #[test]
    fn evalue_decreases_with_score() {
        let ka = KarlinAltschul::ungapped(&AffineScoring::new(2, -3, -5, -1)).unwrap();
        assert!(ka.evalue(60, 400, 100_000, 300) < ka.evalue(40, 400, 100_000, 300));
    }
}
