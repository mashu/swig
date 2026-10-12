//! Outer-flank extension of local alignments.
//!
//! [`extend_outer`] stitches the co-linear V prefix or J suffix back onto a
//! local alignment and adds its substitution score, so a V (J) alignment
//! reaches the germline's 5′ (3′) end or the read's end, whichever comes
//! first, unless the flank stops looking like germline (below). The flank carries allele evidence that a local alignment would drop
//! (a mismatch-rich 5′ end of V), and ranking on the extended score uses it.
//!
//! The stitch stops where the flank stops looking like germline: once its
//! running score falls more than `X_DROP` below the alignment edge. A
//! 5′-truncated V behind an adapter, a UMI or random flank would otherwise
//! take those bases as V: they would move the V start, inflate the read's
//! divergence, and let junk bases decide between alleles that differ only
//! there. The local alignment is optimal, so every flank prefix scores at
//! most 0. Germline mutated at 30 % still gains 0.5 a column on average at
//! (2, −3) and rarely falls that far, while non-V bases lose 1.75 a column
//! and are cut after about 17. Junk within that reach of the edge is still
//! stitched on.

use crate::align::{Alignment, Op};
use crate::scoring::{Score, SubstitutionProfile};

/// How far the flank's running score may fall below the alignment edge
/// while it is stitched on: 10 net mismatches at (2, −3). A drop of 10
/// (the audit's first proposal) cut real flanks at 30 % SHM and cost 1.3
/// points of V and 1.1 of J allele accuracy there; 30 costs none.
const X_DROP: Score = 30;

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
/// The extension walks outward over the shorter of the unaligned query flank
/// and the unaligned reference flank, so the added columns stay on one
/// diagonal, and stops before the column at which the running score falls
/// below `−X_DROP`. Returns whether any base was added.
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
    let room = if left {
        (alignment.q_start as usize).min(alignment.r_start as usize)
    } else {
        (query.len().saturating_sub(alignment.q_end as usize))
            .min(reference.len().saturating_sub(alignment.r_end as usize))
    };
    // The outermost column of the flank at distance `k` from the edge.
    let column = |k: usize| {
        if left {
            (
                alignment.q_start as usize - 1 - k,
                alignment.r_start as usize - 1 - k,
            )
        } else {
            (alignment.q_end as usize + k, alignment.r_end as usize + k)
        }
    };
    let (mut running, mut length, mut kept) = (0, 0, 0);
    for k in 0..room {
        let (q, r) = column(k);
        running += profile.at(query[q], r);
        if running < -X_DROP {
            break;
        }
        (length, kept) = (k + 1, running);
    }
    if length == 0 {
        return false;
    }
    alignment.score += kept;
    let added = vec![Op::Match; length];
    if left {
        alignment.ops.splice(0..0, added);
        alignment.q_start -= length as u32;
        alignment.r_start -= length as u32;
    } else {
        alignment.ops.extend(added);
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

    /// A 40-base flank in front of a 5′-truncated V: stitched in full when
    /// it is germline mutated at 25 %, cut within `X_DROP / 3` columns when
    /// every column mismatches.
    #[test]
    fn outer_prefix_stops_where_the_flank_is_not_v() {
        let mut next = crate::testing::rng(2);
        let reference = crate::testing::random_codes(&mut next, 60);
        let local = |flank: &[u8]| {
            let mut alignment = Alignment {
                q_start: 40,
                q_end: 60,
                r_start: 40,
                r_end: 60,
                score: 40,
                ops: vec![Op::Match; 20],
            };
            let query = [flank, &reference[40..]].concat();
            extend_outer(
                &mut alignment,
                &query,
                &reference,
                &profile(&reference),
                OuterExtend::Prefix,
            );
            assert_eq!(
                alignment.ops.len() as u32,
                alignment.q_end - alignment.q_start
            );
            alignment
        };
        let junk: Vec<u8> = reference[..40].iter().map(|&c| 3 - c).collect();
        let cut = local(&junk);
        assert_eq!(
            cut.q_start as Score,
            40 - X_DROP / 3,
            "stitched from {}",
            cut.q_start
        );
        assert!(cut.score >= 40 - X_DROP, "score {}", cut.score);
        // Every fourth base mutated: +3 every four columns.
        let mutated: Vec<u8> = (0..40)
            .map(|i| {
                if i % 4 == 3 {
                    3 - reference[i]
                } else {
                    reference[i]
                }
            })
            .collect();
        let whole = local(&mutated);
        assert_eq!(
            (whole.q_start, whole.r_start, whole.score),
            (0, 0, 40 + 10 * 3)
        );
    }
}
