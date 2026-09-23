//! Post-adjudication refinement of the winning path only.
//!
//! The note's target architecture refines boundaries on the winning path
//! alone (clip-rescore or a short banded re-alignment). Refiners receive the
//! oriented query, may append new placements (e.g. a re-aligned segment) and
//! rewrite the path. Natural implementations: terminal-substitution extension
//! of V 3′ / J 5′ edges, the R-optimised high-SHM score-tie extension, or a
//! short full-matrix re-alignment of D inside its final cut interval.

use crate::adjudicate::CallPath;
use crate::placement::Placement;

pub trait Refiner: Send + Sync {
    fn refine(&self, query: &[u8], placements: &mut Vec<Placement>, path: &mut CallPath);
}

/// Leave the adjudicated path unchanged.
#[derive(Clone, Copy, Debug, Default)]
pub struct NoRefine;

impl Refiner for NoRefine {
    fn refine(&self, _: &[u8], _: &mut Vec<Placement>, _: &mut CallPath) {}
}
