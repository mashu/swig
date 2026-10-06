//! Post-adjudication refinement of the winning path only.
//!
//! The note's target architecture refines boundaries on the winning path
//! alone (clip-rescore or a short banded re-alignment). Refiners receive the
//! oriented query, may append new placements (e.g. a re-aligned segment) and
//! rewrite the path.
//!
//! [`extend_outer`] also lives here: it stitches the co-linear V prefix or J
//! suffix back onto a local alignment and adds its substitution score, so a
//! V (J) alignment reaches the germline's 5′ (3′) end or the read's end,
//! whichever comes first. The flank carries allele evidence that a local
//! alignment would drop (a mismatch-rich 5′ end of V), and ranking on the
//! extended score uses it.

use crate::adjudicate::CallPath;
use crate::align::{Alignment, Op};
use crate::placement::Placement;
use crate::scoring::SubstitutionProfile;

pub trait Refiner: Send + Sync {
    fn refine(&self, query: &[u8], placements: &mut Vec<Placement>, path: &mut CallPath);
}

/// Leave the adjudicated path unchanged.
#[derive(Clone, Copy, Debug, Default)]
pub struct NoRefine;

impl Refiner for NoRefine {
    fn refine(&self, _: &[u8], _: &mut Vec<Placement>, _: &mut CallPath) {}
}

/// Which unaligned flank of a local alignment to stitch back on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OuterExtend {
    None,
    /// V: co-linear bases before the alignment.
    Prefix,
    /// J: co-linear bases after the alignment.
    Suffix,
}

/// Stitch co-linear outer bases onto `alignment` and add their substitution scores.
///
/// The extension length is the shorter of the unaligned query flank and the
/// unaligned reference flank, so the added columns stay on one diagonal.
/// Returns whether any base was added.
pub fn extend_outer(
    alignment: &mut Alignment,
    query: &[u8],
    reference: &[u8],
    profile: &SubstitutionProfile,
    side: OuterExtend,
) -> bool {
    match side {
        OuterExtend::None => false,
        OuterExtend::Prefix => extend_flank(alignment, query, reference, profile, true),
        OuterExtend::Suffix => extend_flank(alignment, query, reference, profile, false),
    }
}

fn extend_flank(
    alignment: &mut Alignment,
    query: &[u8],
    reference: &[u8],
    profile: &SubstitutionProfile,
    left: bool,
) -> bool {
    let length = if left {
        (alignment.q_start as usize).min(alignment.r_start as usize)
    } else {
        (query.len().saturating_sub(alignment.q_end as usize))
            .min(reference.len().saturating_sub(alignment.r_end as usize))
    };
    if length == 0 {
        return false;
    }
    let (q0, r0) = if left {
        (
            alignment.q_start as usize - length,
            alignment.r_start as usize - length,
        )
    } else {
        (alignment.q_end as usize, alignment.r_end as usize)
    };
    let mut added = Vec::with_capacity(length);
    for offset in 0..length {
        alignment.score += profile.at(query[q0 + offset], r0 + offset);
        added.push(Op::Match);
    }
    if left {
        added.append(&mut alignment.ops);
        alignment.ops = added;
        alignment.q_start -= length as u32;
        alignment.r_start -= length as u32;
    } else {
        alignment.ops.append(&mut added);
        alignment.q_end += length as u32;
        alignment.r_end += length as u32;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::scoring::AffineScoring;

    fn profile(reference: &[u8]) -> SubstitutionProfile {
        SubstitutionProfile::build(&AffineScoring::new(2, -3, -9, -1), reference)
    }

    #[test]
    fn outer_prefix_adds_substitution_score() {
        let reference = encode(b"ACGTACGT");
        let query = encode(b"ACGAACGT");
        let mut alignment = Alignment {
            q_start: 4,
            q_end: 8,
            r_start: 4,
            r_end: 8,
            score: 8,
            ops: vec![Op::Match; 4],
        };
        assert!(extend_outer(
            &mut alignment,
            &query,
            &reference,
            &profile(&reference),
            OuterExtend::Prefix
        ));
        assert_eq!(
            (alignment.q_start, alignment.r_start, alignment.ops.len()),
            (0, 0, 8)
        );
        // Three matches (+2) and one mismatch (−3) on top of the original 8.
        assert_eq!(alignment.score, 8 + 3 * 2 - 3);
    }
}
