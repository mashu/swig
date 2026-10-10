//! Helpers shared by the unit tests.

use crate::align::{Alignment, Op};
use crate::clip::ClipIndex;
use crate::evidence::Evidence;
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::{AffineScoring, SubstitutionProfile};
use crate::sim::Rng;

/// A deterministic stream of random words.
pub fn rng(seed: u64) -> impl FnMut() -> u64 {
    let mut rng = Rng::new(seed);
    move || rng.next_u64()
}

/// `len` random base codes from `next`.
pub fn random_codes(next: &mut impl FnMut() -> u64, len: usize) -> Vec<u8> {
    (0..len).map(|_| (next() % 4) as u8).collect()
}

/// An ungapped placement of `kind` aligning `reference` to `query` from
/// `q_start`, scored with `(2, −3, −9, −1)`.
pub fn placement(kind: SegmentKind, query: &[u8], reference: &[u8], q_start: u32) -> Placement {
    let scoring = AffineScoring::new(2, -3, -9, -1);
    let profile = SubstitutionProfile::build(&scoring, reference);
    let len = reference.len() as u32;
    let score = (0..reference.len())
        .map(|k| profile.at(query[q_start as usize + k], k))
        .sum();
    let alignment = Alignment {
        q_start,
        q_end: q_start + len,
        r_start: 0,
        r_end: len,
        score,
        ops: vec![Op::Match; reference.len()],
    };
    Placement {
        kind,
        alleles: vec![0],
        clip: ClipIndex::build(&alignment, query, &profile, scoring.gaps),
        evidence: Evidence::default(),
        anchor: true,
        alignment,
    }
}
