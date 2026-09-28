//! Swig calling profiles as configuration values.
//!
//! [`EngineConfig::default`](crate::config::EngineConfig::default) is R-optimized,
//! which is Swig's default (`aer_robust` + `r_optimized`).

use crate::adjudicate::DPresencePenalty;
use crate::config::EngineConfig;
use crate::germline::SegmentKind;
use crate::scoring::AffineScoring;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CallingProfile {
    /// Legacy human-IGH tuples. V gap-open −5, `min_d_match` 6, no D-presence cost.
    TruthOptimized,
    /// Swig default. V gap-open −9, `min_d_match` 4, distributed-D, `P_D` 12 relaxed to 10.
    ROptimized,
    /// IgBLAST-agreement: D `(2,−4,−11,−1)`, J `(2,−4,−13,−1)`, `min_d_match` 5.
    IgBlastAgreement,
    /// R-optimized scoring with a fixed D-presence cost of 10.
    SensitiveD,
}

impl CallingProfile {
    pub const ALL: [CallingProfile; 4] = [
        Self::TruthOptimized,
        Self::ROptimized,
        Self::IgBlastAgreement,
        Self::SensitiveD,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Self::TruthOptimized => "truth",
            Self::ROptimized => "r",
            Self::IgBlastAgreement => "igblast",
            Self::SensitiveD => "sensitive-d",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "truth" | "truth_optimized" => Some(Self::TruthOptimized),
            "r" | "r_optimized" | "r-optimized" => Some(Self::ROptimized),
            "igblast" | "igblast_compatible" => Some(Self::IgBlastAgreement),
            "sensitive-d" | "sensitive_d" => Some(Self::SensitiveD),
            _ => None,
        }
    }

    pub fn config(self) -> EngineConfig {
        let mut config = EngineConfig::default();
        match self {
            Self::ROptimized => {}
            Self::TruthOptimized => {
                config.v.scoring = AffineScoring::new(2, -3, -5, -1);
                config.d.min_exact_run = 6;
                config.d.distributed = None;
                config.min_retained[SegmentKind::D.index()] = 6;
                config.objective.d_presence = DPresencePenalty::NONE;
                config.options.junction_tie = false;
            }
            Self::IgBlastAgreement => {
                config.v.scoring = AffineScoring::new(2, -3, -5, -1);
                config.d.scoring = AffineScoring::new(2, -4, -11, -1);
                config.j.scoring = AffineScoring::new(2, -4, -13, -1);
                config.d.min_exact_run = 5;
                config.d.distributed = None;
                config.min_retained[SegmentKind::D.index()] = 5;
                config.objective.d_presence = DPresencePenalty::NONE;
                config.options.junction_tie = false;
            }
            Self::SensitiveD => {
                config.objective.d_presence = DPresencePenalty::SENSITIVE_D;
            }
        }
        config
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_round_trip_and_profiles_differ() {
        for profile in CallingProfile::ALL {
            assert_eq!(CallingProfile::from_name(profile.name()), Some(profile));
        }
        let configs: Vec<EngineConfig> = CallingProfile::ALL.iter().map(|p| p.config()).collect();
        for (i, a) in configs.iter().enumerate() {
            for b in &configs[i + 1..] {
                assert_ne!(a, b);
            }
        }
        assert!(CallingProfile::ROptimized.config().d.distributed.is_some());
        assert!(CallingProfile::TruthOptimized
            .config()
            .d
            .distributed
            .is_none());
        let r = CallingProfile::ROptimized.config();
        assert_eq!(r.v.scoring, AffineScoring::new(2, -3, -9, -1));
        assert_eq!(r.d.min_exact_run, 4);
        assert_eq!(r.objective.d_presence, DPresencePenalty::R_OPTIMIZED);
        assert!(r.options.junction_tie);
        assert_eq!(EngineConfig::default().d.min_exact_run, 4);
        let sensitive = CallingProfile::SensitiveD.config();
        assert_eq!(sensitive.v.scoring, r.v.scoring);
        assert_eq!(
            sensitive.objective.d_presence,
            DPresencePenalty::SENSITIVE_D
        );
        let igblast = CallingProfile::IgBlastAgreement.config();
        assert_eq!(igblast.j.scoring, AffineScoring::new(2, -4, -13, -1));
        assert_eq!(igblast.d.min_exact_run, 5);
        assert_eq!(
            CallingProfile::from_name("r_optimized"),
            Some(CallingProfile::ROptimized)
        );
        assert_eq!(
            CallingProfile::from_name("sensitive_d"),
            Some(CallingProfile::SensitiveD)
        );
    }
}
