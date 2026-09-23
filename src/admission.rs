//! D admission gates: which D placements may enter adjudication.
//!
//! The method note's hard gates are expressed as small composable values:
//! `AnyOf<ExactRunFloor, Option<DistributedEvidence>>` is "exact-run floor, or
//! distributed-D evidence when enabled". Turning a gate into a soft prior later
//! means moving it into [`crate::adjudicate::PathModel`] and admitting all.

use crate::evidence::rescore;
use crate::placement::Placement;
use crate::scoring::AffineScoring;

/// Exact-run floor, or distributed-D evidence when enabled (`Some`).
pub type DefaultAdmission = AnyOf<ExactRunFloor, Option<DistributedEvidence>>;

pub trait Admission: Send + Sync {
    /// `query` is the oriented query, `reference` the primary allele's codes.
    fn admit(&self, placement: &Placement, query: &[u8], reference: &[u8]) -> bool;
}

/// Admit every placement.
#[derive(Clone, Copy, Debug, Default)]
pub struct AdmitAll;

impl Admission for AdmitAll {
    fn admit(&self, _: &Placement, _: &[u8], _: &[u8]) -> bool {
        true
    }
}

/// `min_d_match`: longest run of consecutive exact matches must reach `min_run`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExactRunFloor {
    pub min_run: u32,
}

impl Admission for ExactRunFloor {
    fn admit(&self, placement: &Placement, _: &[u8], _: &[u8]) -> bool {
        placement.evidence.longest_exact_run >= self.min_run
    }
}

/// AER-R distributed-D gate: many aligned columns at high identity, rescored
/// under a fixed tuple.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistributedEvidence {
    pub min_columns: u32,
    pub min_matches: u32,
    pub min_identity: f64,
    pub min_rescore: i32,
    pub scoring: AffineScoring,
}

impl Default for DistributedEvidence {
    fn default() -> Self {
        Self {
            min_columns: 16,
            min_matches: 14,
            min_identity: 0.80,
            min_rescore: 19,
            scoring: AffineScoring::new(2, -3, -13, -1),
        }
    }
}

impl Admission for DistributedEvidence {
    fn admit(&self, placement: &Placement, query: &[u8], reference: &[u8]) -> bool {
        let e = &placement.evidence;
        e.match_columns >= self.min_columns
            && e.matches >= self.min_matches
            && e.identity() >= self.min_identity
            && rescore(placement.alignment.columns(), query, reference, &self.scoring) >= self.min_rescore
    }
}

/// A disabled gate (`None`) admits nothing, so `AnyOf<G, Option<H>>` reduces to `G`.
impl<A: Admission> Admission for Option<A> {
    fn admit(&self, placement: &Placement, query: &[u8], reference: &[u8]) -> bool {
        self.as_ref().is_some_and(|gate| gate.admit(placement, query, reference))
    }
}

/// Admit if either gate admits.
#[derive(Clone, Copy, Debug, Default)]
pub struct AnyOf<A, B>(pub A, pub B);

impl<A: Admission, B: Admission> Admission for AnyOf<A, B> {
    fn admit(&self, placement: &Placement, query: &[u8], reference: &[u8]) -> bool {
        self.0.admit(placement, query, reference) || self.1.admit(placement, query, reference)
    }
}

/// Admit only if both gates admit.
#[derive(Clone, Copy, Debug, Default)]
pub struct AllOf<A, B>(pub A, pub B);

impl<A: Admission, B: Admission> Admission for AllOf<A, B> {
    fn admit(&self, placement: &Placement, query: &[u8], reference: &[u8]) -> bool {
        self.0.admit(placement, query, reference) && self.1.admit(placement, query, reference)
    }
}
