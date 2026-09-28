//! Post-adjudication refinement of the winning path only.
//!
//! The note's target architecture refines boundaries on the winning path
//! alone (clip-rescore or a short banded re-alignment). Refiners receive the
//! oriented query, may append new placements (e.g. a re-aligned segment) and
//! rewrite the path.
//!
//! Two Swig endpoint rules live here as plain functions, because they change
//! coordinates (and, for the outer extension, the score) of an alignment the
//! dynamic program already produced:
//!
//! * [`extend_outer`] stitches the co-linear V prefix or J suffix back onto a
//!   local alignment and adds its substitution score. Swig does this for every
//!   profile, before ranking.
//! * [`extend_junction_tie`] writes a zero-score five-base block onto the V 3′
//!   end or the J 5′ start when V identity is at most 0.88. Score does not
//!   change. R-optimized and Sensitive-D only.

use crate::adjudicate::{CallPath, CalledSegment};
use crate::align::{Alignment, Op};
use crate::alphabet::N;
use crate::clip::ClipIndex;
use crate::evidence::Evidence;
use crate::placement::Placement;
use crate::scoring::{GapCosts, Score, SubstitutionProfile};

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

/// V identity at or below this selects the high-SHM junction tie.
pub const JUNCTION_TIE_IDENTITY: f64 = 0.88;
/// Length of the zero-score block written at a junction-facing end.
pub const JUNCTION_TIE_LENGTH: usize = 5;

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

/// Add exactly `length` co-linear Match columns. `left` prepends them at `(q, r)`.
fn append_matches(
    alignment: &mut Alignment,
    query: &[u8],
    profile: &SubstitutionProfile,
    q: usize,
    r: usize,
    length: usize,
    left: bool,
) {
    let mut added = Vec::with_capacity(length);
    for offset in 0..length {
        alignment.score += profile.at(query[q + offset], r + offset);
        added.push(Op::Match);
    }
    if left {
        added.append(&mut alignment.ops);
        alignment.ops = added;
        alignment.q_start = q as u32;
        alignment.r_start = r as u32;
    } else {
        alignment.ops.append(&mut added);
        alignment.q_end = (q + length) as u32;
        alignment.r_end = (r + length) as u32;
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

/// Materialize a junction-facing five-base block whose substitution score is exactly zero.
///
/// The block must lie entirely inside sequence not occupied by the next selected
/// segment (`query_boundary` is that segment's facing edge), every base must be
/// canonical, and no shorter prefix of the block may have a positive score.
/// The called span must already end (or start) on the alignment edge: this
/// extends the selected alignment, it does not un-clip one.
pub fn extend_junction_tie(
    placement: &mut Placement,
    called: &mut CalledSegment,
    query: &[u8],
    reference: &[u8],
    profile: &SubstitutionProfile,
    gaps: GapCosts,
    extend_right: bool,
    query_boundary: u32,
) -> bool {
    let span = called.span;
    let alignment = &placement.alignment;
    let length = JUNCTION_TIE_LENGTH;
    let (q, r) = if extend_right {
        if span.q_hi != alignment.q_end || span.r_hi != alignment.r_end {
            return false;
        }
        if span.q_hi as usize + length > query.len()
            || span.q_hi + length as u32 > query_boundary
            || span.r_hi as usize + length > reference.len()
        {
            return false;
        }
        (span.q_hi as usize, span.r_hi as usize)
    } else {
        if span.q_lo != alignment.q_start || span.r_lo != alignment.r_start {
            return false;
        }
        if (alignment.q_start as usize) < length
            || (alignment.r_start as usize) < length
            || alignment.q_start - (length as u32) < query_boundary
        {
            return false;
        }
        (
            alignment.q_start as usize - length,
            alignment.r_start as usize - length,
        )
    };
    if !zero_score_block(query, reference, q, r, length, profile, extend_right) {
        return false;
    }
    append_matches(
        &mut placement.alignment,
        query,
        profile,
        q,
        r,
        length,
        !extend_right,
    );
    placement.clip = ClipIndex::build(&placement.alignment, query, profile, gaps);
    placement.evidence = Evidence::from_columns(placement.alignment.columns(), query, reference);
    if extend_right {
        called.exit = placement.alignment.q_end;
    } else {
        called.entry = placement.alignment.q_start;
    }
    called.span = placement
        .clip
        .span(called.entry, called.exit)
        .expect("a zero-score extension is a run of Match columns");
    true
}

/// True when the `length` bases at `(q, r)` sum to 0 and no edge-adjacent prefix is positive.
fn zero_score_block(
    query: &[u8],
    reference: &[u8],
    q: usize,
    r: usize,
    length: usize,
    profile: &SubstitutionProfile,
    extend_right: bool,
) -> bool {
    let mut total = 0;
    let mut cumulative = 0;
    let mut maximum_prefix = Score::MIN;
    for step in 0..length {
        let offset = if extend_right {
            step
        } else {
            length - step - 1
        };
        let (query_base, reference_base) = (query[q + offset], reference[r + offset]);
        if query_base >= N || reference_base >= N {
            return false;
        }
        let contribution = profile.at(query_base, r + offset);
        total += contribution;
        cumulative += contribution;
        maximum_prefix = maximum_prefix.max(cumulative);
    }
    total == 0 && maximum_prefix <= 0
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

    #[test]
    fn junction_block_rejects_a_positive_prefix() {
        let reference = encode(b"AAAAA");
        let query = encode(b"AACAA");
        let built = profile(&reference);
        // Walking from the right edge: A, A, C, A, A. The first base matches (+2).
        assert!(!zero_score_block(
            &query, &reference, 0, 0, 5, &built, false
        ));
        // Extend right, reference AAAAA, query GGAAA: two mismatches then three matches.
        // Prefix sums stay negative and the total is zero.
        let query = encode(b"GGAAA");
        assert!(zero_score_block(&query, &reference, 0, 0, 5, &built, true));
    }
}
