//! Karlin–Altschul support values for reporting (never used for calling).
//!
//! `E = K m′ n′ e^{−λS}` with the BLAST length adjustment
//! `ℓ = ln(K m n) / H`, `m′ = m − ℓ`, `n′ = n − N ℓ`.
//!
//! [`KarlinAltschul::ungapped`] computes λ, H and K exactly for a match /
//! mismatch tuple under uniform base composition: λ by bisection on
//! `Σ p_s e^{λs} = 1`, H as the tilted mean, and K from the Karlin–Altschul
//! series `σ = Σ_k k⁻¹ E[min(e^{λ S_k}, 1)]` over convolutions of the score
//! distribution. (It reproduces NCBI's ungapped blastn table, e.g. reward 1 /
//! penalty −3 → λ 1.374, K 0.711, H 1.31.) Gapped alignments score slightly
//! higher than ungapped ones, so these E-values are mildly anti-conservative;
//! for calibrated values fit λ and K by simulation and construct the struct
//! directly.

use crate::scoring::{AffineScoring, Score};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KarlinAltschul {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
}

fn gcd(a: u32, b: u32) -> u32 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
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

        // σ series over the distribution of S_k (support k·penalty ..= k·reward).
        let (r, q) = (reward as usize, penalty.unsigned_abs() as usize);
        let mut dist = vec![1.0f64];
        let mut sigma = 0.0;
        for k in 1..=500usize {
            let mut next = vec![0.0f64; dist.len() + r + q];
            for (i, &pr) in dist.iter().enumerate().filter(|&(_, &pr)| pr > 0.0) {
                next[i + r + q] += pr * p_match;
                next[i] += pr * p_mismatch;
            }
            dist = next;
            let offset = (k * q) as f64; // index i holds score i − k·q
            let term: f64 = dist
                .iter()
                .enumerate()
                .map(|(i, &pr)| {
                    let s = i as f64 - offset;
                    pr * if s >= 0.0 { 1.0 } else { (lambda * s).exp() }
                })
                .sum::<f64>()
                / k as f64;
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

    #[test]
    fn evalue_decreases_with_score() {
        let ka = KarlinAltschul::ungapped(&AffineScoring::new(2, -3, -5, -1)).unwrap();
        assert!(ka.evalue(60, 400, 100_000, 300) < ka.evalue(40, 400, 100_000, 300));
    }
}
