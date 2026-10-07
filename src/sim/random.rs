//! Deterministic random numbers and the handful of distributions IgSim uses.
//!
//! The generator is xoshiro256** seeded through SplitMix64, so a seed fixes
//! every simulated read on every platform. The distributions mirror the
//! `Distributions.jl` types IgSim.jl is written against, with the same
//! parameterisations:
//!
//! | Here | Distributions.jl | Support, mean |
//! |---|---|---|
//! | [`Count::Uniform`] | `DiscreteUniform(a, b)` | `a..=b` |
//! | [`Count::Geometric`] | `Geometric(p)` | failures before the first success, `(1−p)/p` |
//! | [`Count::NegativeBinomial`] | `NegativeBinomial(r, p)` | failures before the `r`-th success, `r(1−p)/p` |
//! | [`Count::Mixture`] | IgSim `DomainFlankMix` | two uniforms |
//! | [`Rate::Uniform`] | `Uniform(a, b)` | |
//! | [`Rate::LogNormal`] | `LogNormal(μ, σ)` | |

/// xoshiro256** (Blackman & Vigna), seeded with SplitMix64.
#[derive(Clone, Debug)]
pub struct Rng {
    state: [u64; 4],
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Self {
            state: [next(), next(), next(), next()],
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.state;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Uniform in `[0, 1)` with 53 random bits.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform integer in `0..n` (Lemire's multiply-shift; `n > 0`).
    pub fn below(&mut self, n: usize) -> usize {
        debug_assert!(n > 0);
        ((u128::from(self.next_u64()) * n as u128) >> 64) as usize
    }

    /// Uniform integer in `lo..=hi`.
    pub fn range(&mut self, lo: u32, hi: u32) -> u32 {
        lo + self.below((hi - lo) as usize + 1) as u32
    }

    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.uniform() < p
    }

    /// Standard normal (Box–Muller; one of the pair is discarded for simplicity).
    pub fn normal(&mut self) -> f64 {
        let u = 1.0 - self.uniform(); // (0, 1]
        let v = self.uniform();
        (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
    }

    /// One of `A C G T` as ASCII.
    pub fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4)]
    }

    /// `len` random bases as ASCII.
    pub fn bases(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.base()).collect()
    }

    /// `len` random base codes (`0..4`, see [`crate::alphabet`]).
    pub fn codes(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.below(4) as u8).collect()
    }

    /// A different canonical base than `base` (uniform over the other three).
    pub fn substitute(&mut self, base: u8) -> u8 {
        let others: [u8; 3] = match base {
            b'A' => *b"CGT",
            b'C' => *b"AGT",
            b'G' => *b"ACT",
            _ => *b"ACG",
        };
        others[self.below(3)]
    }
}

/// A non-negative integer-valued distribution (lengths, trims).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Count {
    Fixed(u32),
    /// `DiscreteUniform(lo, hi)`, both inclusive.
    Uniform {
        lo: u32,
        hi: u32,
    },
    /// Failures before the first success with probability `p`.
    Geometric {
        p: f64,
    },
    /// Failures before the `r`-th success with probability `p`.
    NegativeBinomial {
        r: u32,
        p: f64,
    },
    /// IgSim's `DomainFlankMix`: `short` with probability `p_short`, else `long`.
    Mixture {
        short: (u32, u32),
        long: (u32, u32),
        p_short: f64,
    },
}

impl Count {
    pub fn sample(&self, rng: &mut Rng) -> u32 {
        match *self {
            Count::Fixed(n) => n,
            Count::Uniform { lo, hi } => rng.range(lo, hi),
            Count::Geometric { p } => geometric(rng, p),
            Count::NegativeBinomial { r, p } => (0..r).map(|_| geometric(rng, p)).sum(),
            Count::Mixture {
                short,
                long,
                p_short,
            } => {
                let (lo, hi) = if rng.bernoulli(p_short) { short } else { long };
                rng.range(lo, hi)
            }
        }
    }
}

/// Inversion: `⌊ln U / ln(1 − p)⌋` with `U ∈ (0, 1]`.
fn geometric(rng: &mut Rng, p: f64) -> u32 {
    if p >= 1.0 {
        return 0;
    }
    let u = 1.0 - rng.uniform();
    (u.ln() / (1.0 - p).ln()).floor().min(f64::from(u32::MAX)) as u32
}

/// A per-base rate in `[0, 1]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Rate {
    Fixed(f64),
    Uniform { lo: f64, hi: f64 },
    LogNormal { mu: f64, sigma: f64 },
}

impl Rate {
    pub fn sample(&self, rng: &mut Rng) -> f64 {
        let x = match *self {
            Rate::Fixed(x) => x,
            Rate::Uniform { lo, hi } => lo + (hi - lo) * rng.uniform(),
            Rate::LogNormal { mu, sigma } => (mu + sigma * rng.normal()).exp(),
        };
        x.clamp(0.0, 1.0)
    }
}

/// IgSim's `GatedRate`: with probability `p` draw a rate from `rate`, else 0.
///
/// A real library is a mixture of unmutated reads and a mutated memory tail;
/// the gate expresses that zero inflation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GatedRate {
    pub p: f64,
    pub rate: Rate,
}

impl GatedRate {
    pub const NONE: Self = Self {
        p: 0.0,
        rate: Rate::Fixed(0.0),
    };

    pub fn sample(&self, rng: &mut Rng) -> f64 {
        if rng.bernoulli(self.p) {
            self.rate.sample(rng)
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mean(samples: impl Iterator<Item = f64>, n: usize) -> f64 {
        samples.take(n).sum::<f64>() / n as f64
    }

    #[test]
    fn distributions_have_their_documented_means() {
        let mut rng = Rng::new(1);
        let n = 200_000;
        let geometric = Count::Geometric { p: 0.2 };
        let m = mean(
            std::iter::repeat_with(|| f64::from(geometric.sample(&mut rng))),
            n,
        );
        assert!((m - 4.0).abs() < 0.05, "Geometric(0.2) mean {m}");
        let nb = Count::NegativeBinomial { r: 3, p: 0.35 };
        let m = mean(std::iter::repeat_with(|| f64::from(nb.sample(&mut rng))), n);
        assert!(
            (m - 3.0 * 0.65 / 0.35).abs() < 0.05,
            "NegativeBinomial mean {m}"
        );
        let uniform = Count::Uniform { lo: 3, hi: 7 };
        let m = mean(
            std::iter::repeat_with(|| f64::from(uniform.sample(&mut rng))),
            n,
        );
        assert!((m - 5.0).abs() < 0.02);
        let lognormal = Rate::LogNormal {
            mu: -2.65,
            sigma: 0.5,
        };
        let m = mean(std::iter::repeat_with(|| lognormal.sample(&mut rng)), n);
        let expected = (-2.65f64 + 0.125).exp();
        assert!(
            (m - expected).abs() < 0.002,
            "LogNormal mean {m} vs {expected}"
        );
        assert!((0..1000).all(|_| GatedRate::NONE.sample(&mut rng) == 0.0));
    }

    #[test]
    fn seeds_are_reproducible() {
        let (mut a, mut b) = (Rng::new(42), Rng::new(42));
        assert!((0..100).all(|_| a.next_u64() == b.next_u64()));
        assert_ne!(Rng::new(1).next_u64(), Rng::new(2).next_u64());
    }
}
