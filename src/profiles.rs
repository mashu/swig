//! The method note's calling profiles as configuration values.

use crate::adjudicate::DPresencePenalty;
use crate::admission::DistributedEvidence;
use crate::config::EngineConfig;
use crate::scoring::AffineScoring;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CallingProfile {
    /// Default production tuples: V (2,−3,−5,−1), D (2,−3,−13,−1), J (2,−3,−17,−2);
    /// exact-run D admission only (distributed-D is an AER-R feature).
    TruthOptimized,
    /// AER-R: V gap-open −9, `min_d_match` 5, distributed-D gate on.
    ROptimized,
    /// D (2,−4,−11,−1).
    IgBlastAgreement,
    /// Fixed `P_D` 10, `min_d_match` 4.
    SensitiveD,
}

impl CallingProfile {
    pub const ALL: [CallingProfile; 4] =
        [Self::TruthOptimized, Self::ROptimized, Self::IgBlastAgreement, Self::SensitiveD];

    pub fn name(self) -> &'static str {
        match self {
            Self::TruthOptimized => "truth",
            Self::ROptimized => "r",
            Self::IgBlastAgreement => "igblast",
            Self::SensitiveD => "sensitive-d",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|p| p.name() == name)
    }

    pub fn config(self) -> EngineConfig {
        let mut config = EngineConfig::default();
        match self {
            Self::TruthOptimized => {}
            Self::ROptimized => {
                config.v.scoring = AffineScoring::new(2, -3, -9, -1);
                config.d.min_exact_run = 5;
                config.d.distributed = Some(DistributedEvidence::default());
            }
            Self::IgBlastAgreement => {
                config.d.scoring = AffineScoring::new(2, -4, -11, -1);
            }
            Self::SensitiveD => {
                config.d.min_exact_run = 4;
                config.d.distributed = Some(DistributedEvidence::default());
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
        assert!(CallingProfile::TruthOptimized.config().d.distributed.is_none());
    }
}
