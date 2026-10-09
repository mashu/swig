//! Simulation parameters and IgSim.jl's presets.
//!
//! Every field is a distribution, sampled once per read (trims, lengths,
//! rates) or once per base (substitutions and indels at the sampled rate).
//! The presets reproduce IgSim.jl's `train_params`, `easy_params`,
//! `mid_params` and `hard_params` exactly, parameter for parameter; IgSim
//! has no constant region, so `include_c` and `c_whole_flank` are this
//! crate's, the same in every preset.

use super::random::{Count, GatedRate, Rate};

/// Instrument substitutions on the V(D)J body, ramping linearly from `p5` at
/// the first body base to `p3` at the last (IgSim's `IlluminaError`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IlluminaError {
    pub p5: f64,
    pub p3: f64,
}

impl IlluminaError {
    pub const NONE: Self = Self { p5: 0.0, p3: 0.0 };
    /// Merged paired-end MiSeq, uniform ~Q30.
    pub const MISEQ: Self = Self {
        p5: 0.001,
        p3: 0.001,
    };
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimParams {
    pub v_trim_5p: Count,
    pub v_trim_3p: Count,
    pub d_trim_5p: Count,
    pub d_trim_3p: Count,
    pub j_trim_5p: Count,
    /// N1 (V–D, or V–J without D).
    pub n1_length: Count,
    /// N2 (D–J); also the N between two tandem Ds.
    pub n2_length: Count,
    /// Probability that the read carries a D.
    pub include_d: f64,
    /// Probability of a second D joining the DJ (DJ → DDJ → VDDJ).
    pub include_dd: f64,
    pub flank_5p: Count,
    pub flank_3p: Count,
    /// Per-read SHM substitution rate on body bases.
    pub body_error_rate: GatedRate,
    /// Per-read single-base indel rate on body bases.
    pub indel_rate: GatedRate,
    pub illumina_error: IlluminaError,
    pub min_length: usize,
    pub max_length: usize,
    pub max_retries: usize,
    /// Probability that the emitted read is reverse-complemented. IgSim emits
    /// germline-sense reads only (0.0); the assigner orients reads itself, so
    /// evaluations usually set this to 0.5.
    pub reverse_complement: f64,
    /// Probability that the read reaches the constant region, when the
    /// simulator has a C database ([`super::Simulator::append_constant`]).
    pub include_c: f64,
    /// Probability that C fills the whole 3′ flank; otherwise its length is
    /// uniform in `0..=flank`.
    pub c_whole_flank: f64,
}

/// IgSim `shm_igm()`: 40 % of reads mutated at a rate uniform in 0.5–7 %.
pub const SHM_IGM: GatedRate = GatedRate {
    p: 0.40,
    rate: Rate::Uniform {
        lo: 0.005,
        hi: 0.07,
    },
};

/// IgSim `shm_igg()`: 97 % of reads mutated, log-normal rate (median ≈ 7 %).
pub const SHM_IGG: GatedRate = GatedRate {
    p: 0.97,
    rate: Rate::LogNormal {
        mu: -2.65,
        sigma: 0.50,
    },
};

impl SimParams {
    /// IgSim `train_params()` (and `mid_params()`): IgM-like library.
    pub fn train() -> Self {
        Self {
            v_trim_5p: Count::Fixed(0),
            v_trim_3p: Count::Uniform { lo: 0, hi: 10 },
            d_trim_5p: Count::Geometric { p: 0.20 },
            d_trim_3p: Count::Geometric { p: 0.21 },
            j_trim_5p: Count::Uniform { lo: 0, hi: 14 },
            n1_length: Count::Uniform { lo: 0, hi: 18 },
            n2_length: Count::Uniform { lo: 0, hi: 18 },
            include_d: 0.99,
            include_dd: 0.001,
            flank_5p: Count::Mixture {
                short: (10, 60),
                long: (70, 140),
                p_short: 0.55,
            },
            flank_3p: Count::Mixture {
                short: (40, 110),
                long: (80, 160),
                p_short: 0.55,
            },
            body_error_rate: SHM_IGM,
            indel_rate: GatedRate {
                p: 0.08,
                rate: Rate::Uniform { lo: 0.0, hi: 0.002 },
            },
            illumina_error: IlluminaError::MISEQ,
            min_length: 80,
            max_length: 900,
            max_retries: 32,
            reverse_complement: 0.0,
            // A fifth of the reads stop before C: they test that C-free
            // reads get no C call.
            include_c: 0.8,
            c_whole_flank: 0.6,
        }
    }

    /// IgSim `mid_params()` (identical to `train_params()`).
    pub fn mid() -> Self {
        Self::train()
    }

    /// IgSim `easy_params()`: light trims, short N, no SHM, no indels.
    pub fn easy() -> Self {
        Self {
            v_trim_3p: Count::Uniform { lo: 0, hi: 4 },
            d_trim_5p: Count::Uniform { lo: 0, hi: 2 },
            d_trim_3p: Count::Uniform { lo: 0, hi: 2 },
            j_trim_5p: Count::Uniform { lo: 0, hi: 3 },
            n1_length: Count::Uniform { lo: 0, hi: 4 },
            n2_length: Count::Uniform { lo: 0, hi: 4 },
            body_error_rate: GatedRate::NONE,
            indel_rate: GatedRate::NONE,
            flank_5p: Count::Mixture {
                short: (10, 55),
                long: (70, 130),
                p_short: 0.5,
            },
            flank_3p: Count::Mixture {
                short: (40, 100),
                long: (80, 140),
                p_short: 0.5,
            },
            ..Self::train()
        }
    }

    /// IgSim `hard_params()`: heavier trims, negative-binomial N, IgG SHM.
    pub fn hard() -> Self {
        Self {
            v_trim_5p: Count::Uniform { lo: 0, hi: 1 },
            v_trim_3p: Count::Uniform { lo: 0, hi: 18 },
            j_trim_5p: Count::Uniform { lo: 0, hi: 14 },
            n1_length: Count::NegativeBinomial { r: 3, p: 0.35 },
            n2_length: Count::NegativeBinomial { r: 3, p: 0.35 },
            body_error_rate: SHM_IGG,
            indel_rate: GatedRate {
                p: 0.15,
                rate: Rate::Uniform { lo: 0.0, hi: 0.004 },
            },
            flank_5p: Count::Mixture {
                short: (10, 60),
                long: (70, 150),
                p_short: 0.6,
            },
            flank_3p: Count::Mixture {
                short: (35, 110),
                long: (80, 160),
                p_short: 0.6,
            },
            ..Self::train()
        }
    }

    /// Preset by IgSim name: `easy`, `mid`, `train`, `hard`.
    pub fn preset(name: &str) -> Option<Self> {
        match name {
            "easy" => Some(Self::easy()),
            "mid" => Some(Self::mid()),
            "train" => Some(Self::train()),
            "hard" => Some(Self::hard()),
            _ => None,
        }
    }

    /// Same parameters with a fixed SHM rate on every read (for stratified sweeps).
    pub fn with_shm(self, rate: f64) -> Self {
        Self {
            body_error_rate: GatedRate {
                p: 1.0,
                rate: Rate::Fixed(rate),
            },
            ..self
        }
    }

    /// Same parameters with each read reverse-complemented with probability `p`.
    pub fn with_reverse_complement(self, p: f64) -> Self {
        Self {
            reverse_complement: p,
            ..self
        }
    }
}

impl Default for SimParams {
    fn default() -> Self {
        Self::train()
    }
}
