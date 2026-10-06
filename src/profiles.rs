//! Calling profiles as configuration values.
//!
//! The profiles keep the scoring tuples and gates of the method note's
//! profiles; the D-presence cost of `r` and `sensitive-d` is derived from the
//! D database ([`DPrior`]) instead of the note's constants.
//! [`EngineConfig::default`](crate::config::EngineConfig::default) is `r`.

use crate::config::{DPrior, EngineConfig};
use crate::germline::SegmentKind;
use crate::scoring::AffineScoring;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CallingProfile {
    /// Legacy human-IGH tuples. V gap-open −5, `min_d_match` 6, no D-presence cost.
    TruthOptimized,
    /// Default. V gap-open −9, `min_d_match` 4, `P_D` derived from the D
    /// database at even prior odds (11 for human IGH).
    ROptimized,
    /// IgBLAST-agreement: D `(2,−4,−11,−1)`, J `(2,−4,−13,−1)`, `min_d_match` 5.
    IgBlastAgreement,
    /// R-optimized scoring with D twice as likely a priori (`P_D` 9 for human IGH).
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
                config.min_retained[SegmentKind::D.index()] = 6;
                config.objective.d_prior = DPrior::Fixed(0);
            }
            Self::IgBlastAgreement => {
                config.v.scoring = AffineScoring::new(2, -3, -5, -1);
                config.d.scoring = AffineScoring::new(2, -4, -11, -1);
                config.j.scoring = AffineScoring::new(2, -4, -13, -1);
                config.d.min_exact_run = 5;
                config.min_retained[SegmentKind::D.index()] = 5;
                config.objective.d_prior = DPrior::Fixed(0);
            }
            Self::SensitiveD => {
                config.objective.d_prior = DPrior::Derived {
                    window: 25,
                    prior_odds: 2.0,
                };
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
        assert!(CallingProfile::ALL
            .iter()
            .all(|p| p.config().d.distributed.is_none()));
        let r = CallingProfile::ROptimized.config();
        assert_eq!(r.v.scoring, AffineScoring::new(2, -3, -9, -1));
        assert_eq!(r.d.min_exact_run, 4);
        assert_eq!(r.objective.d_prior, DPrior::DERIVED);
        assert_eq!(EngineConfig::default().d.min_exact_run, 4);
        let sensitive = CallingProfile::SensitiveD.config();
        assert_eq!(sensitive.v.scoring, r.v.scoring);
        assert_ne!(sensitive.objective.d_prior, r.objective.d_prior);
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
