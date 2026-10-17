//! Co-optimal placements of a called segment: the ties a report names, and
//! the placements its posterior counts.
//!
//! Every placement of the called kind whose clip to the called interval
//! `[entry, exit)` is worth exactly the called clip under the call's model
//! is tied with it: on the fixed scale an equal score, under a per-read
//! call ([`super::CallModel::PerRead`]) also equal matches and mismatches,
//! which makes the value equal at every divergence. That is "equal value
//! under this caller", not equal biological probability: the probability
//! is the posterior mass of the counted placements, those tied with the
//! call and those whose germline bases between the call's neighbours its
//! primary allele contains. The report ([`crate::report`]) and the C stage
//! ([`crate::constant`]) both read these rules.

use std::ops::Range;

use super::CalledSegment;
use crate::clip::ClipSpan;
use crate::evidence::Evidence;
use crate::germline::{AlleleId, GermlineSet, SegmentKind};
use crate::placement::Placement;

/// The placements tied with a called segment on its interval, and those its
/// posterior counts (see the module documentation).
pub(crate) struct CoOptimal<'a> {
    /// Tied placements, the primary — reported first — at index 0.
    pub tied: Vec<&'a Placement>,
    /// The allele reported first: the call's own ([`CalledSegment::allele`]),
    /// or the primary placement's first.
    pub first: AlleleId,
    /// The alleles the report names, `first` first: every allele of a tied
    /// placement, and of a J only those at the first allele's germline
    /// coordinates ([`positional`]).
    pub named: Vec<AlleleId>,
    /// Indices of the placements of the called kind whose probability the
    /// posterior sums: tied with the call, or carried by its primary allele.
    pub counted: Vec<usize>,
}

impl CoOptimal<'_> {
    /// The named alleles' names, the first allele's first, each once.
    pub fn names<'s>(&self, set: &'s GermlineSet) -> Vec<&'s str> {
        let mut names: Vec<&str> = Vec::new();
        for &allele in &self.named {
            for name in set.allele(allele).names() {
                if !names.contains(&name) {
                    names.push(name);
                }
            }
        }
        names
    }
}

/// Whether a called kind's ties are positional: a J allele tied with the
/// call is named only when it explains the called bases at the first
/// allele's own germline coordinates. J4*02 and J5*02 trimmed to their
/// shared 3′ part align the same read bases from different germline
/// offsets: each is a different trim of a different gene, and the call's
/// support (the DAG's tie-break) chose one. On GenAIRR reads, whose truth
/// compares alleles at the true allele's germline index, this raises IGH J
/// by 0.26–0.38 points; positional D or V ties gain nothing (D −0.05 to
/// −0.08, V ±0.00).
fn positional(kind: SegmentKind) -> bool {
    kind == SegmentKind::J
}

/// Co-optimal alleles of `called` among `placements` (of every kind; only
/// the called one's count), `set` holding their alleles.
///
/// Tied: worth exactly the call on its interval `[entry, exit)` — an equal
/// score, and under a per-read call (`per_read`) equal matches and
/// mismatches. The primary, whose span and alignment a report gives,
/// explains the called read bases: a tied clip that starts later or ends
/// earlier (on the fixed scale the bases it leaves can net 0) would drop
/// read bases the path's N/P regions do not hold. Among those, it is the
/// one with the most germline support (unclipped score, the DAG's own
/// tie-break), then the longest exact run and most matches on the
/// interval, then the lowest allele id; a call that names its allele stays
/// first, and that allele is reported first.
///
/// Counted by the posterior: the tied placements, and every placement that
/// explains most of the read bases the call explains (more than half of the
/// called span, within `region`, the read the call's neighbours leave to
/// it) with germline bases there that the first allele contains. If the
/// read came from such a placement, the call still carries its germline
/// bases. A placement that only reaches into the region by a few bases (a
/// D inside the called J, a tandem D's partner) is another hypothesis, so
/// it does not count.
pub(crate) fn co_optimal<'a>(
    placements: &'a [Placement],
    query: &[u8],
    set: &GermlineSet,
    called: &CalledSegment,
    per_read: bool,
    region: Range<u32>,
) -> CoOptimal<'a> {
    let called_placement = &placements[called.placement];
    let kind = called_placement.kind;
    let called_score = called.span.score;
    let columns = |p: &Placement, span: &ClipSpan| {
        let reference = &set.allele(p.primary_allele()).codes;
        let e = Evidence::from_columns(p.clip.columns(span), query, reference);
        (e.matches, e.mismatches)
    };
    let called_columns = per_read.then(|| columns(called_placement, &called.span));
    let tied_indices: Vec<usize> = (0..placements.len())
        .filter(|&i| {
            let p = &placements[i];
            p.kind == kind
                && p.clip.span(called.entry, called.exit).is_some_and(|s| {
                    s.score == called_score && called_columns.is_none_or(|c| columns(p, &s) == c)
                })
        })
        .collect();
    // The called placement first, so it wins ties of the sort key; one
    // whose allele the call names stays first.
    let mut tied: Vec<&Placement> = std::iter::once(called_placement)
        .chain(
            tied_indices
                .iter()
                .filter(|&&i| i != called.placement)
                .map(|&i| &placements[i]),
        )
        .collect();
    tied[usize::from(called.allele.is_some())..].sort_by_cached_key(|p| {
        let own = p
            .clip
            .span(called.entry, called.exit)
            .expect("tied placement clips to the called interval");
        let e = Evidence::from_columns(
            p.clip.columns(&own),
            query,
            &set.allele(p.primary_allele()).codes,
        );
        (
            (own.q_lo, own.q_hi) != (called.span.q_lo, called.span.q_hi),
            std::cmp::Reverse(p.score()),
            std::cmp::Reverse(e.longest_exact_run),
            std::cmp::Reverse(e.matches),
            p.primary_allele(),
        )
    });
    let first = called.allele.unwrap_or(tied[0].primary_allele());
    // A tied allele's germline coordinates on the called interval.
    let frame = |p: &Placement, allele: AlleleId| {
        let own = p
            .clip
            .span(called.entry, called.exit)
            .expect("tied placement clips to the called interval");
        let shift = i64::from(p.germline_shift(allele));
        (i64::from(own.r_lo) + shift, i64::from(own.r_hi) + shift)
    };
    let at = frame(tied[0], first);
    let mut named = vec![first];
    for p in &tied {
        for &allele in &p.alleles {
            if !named.contains(&allele) && (!positional(kind) || frame(p, allele) == at) {
                named.push(allele);
            }
        }
    }
    let reference = &set.allele(first).codes;
    let called_bases = called.span.q_lo..called.span.q_hi;
    let carried = |p: &Placement| {
        let Some(own) = p.clip.span(region.start, region.end) else {
            return false;
        };
        let shared = own
            .q_hi
            .min(called_bases.end)
            .saturating_sub(own.q_lo.max(called_bases.start));
        if 2 * shared <= called_bases.len() as u32 {
            return false;
        }
        let codes = &set.allele(p.primary_allele()).codes;
        let bases = &codes[own.r_lo as usize..own.r_hi as usize];
        bases.is_empty() || reference.windows(bases.len()).any(|w| w == bases)
    };
    let counted = (0..placements.len())
        .filter(|&i| {
            let p = &placements[i];
            p.kind == kind && (tied_indices.contains(&i) || carried(p))
        })
        .collect();
    CoOptimal {
        tied,
        first,
        named,
        counted,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Alignment, Op};
    use crate::alphabet::encode;
    use crate::clip::ClipIndex;
    use crate::evidence::Evidence;
    use crate::scoring::{AffineScoring, Score, SubstitutionProfile};
    use crate::testing::placement;

    /// IGHJ5*02 holds IGHJ4*02's 3′ part three bases further from its
    /// start. Trimmed to that part, both align the same read bases equally
    /// well: tied, but from different germline coordinates. A J call names
    /// only the first allele (the other stays a tied placement, and a
    /// near-tied alternative of the report); a D call names both.
    #[test]
    fn j_ties_are_named_at_the_first_alleles_germline_coordinates() {
        let core = "ACTACTTTGACTACTGGGGCCAGG";
        let set = GermlineSet::from_fasta(
            SegmentKind::J,
            &format!(">IGHJ4*02\n{core}\n>IGHJ5*02\nGGT{core}\n"),
        )
        .unwrap();
        let query = encode(format!("TTTTT{core}").as_bytes());
        let first = placement(SegmentKind::J, &query, &encode(core.as_bytes()), 5);
        // IGHJ5*02 from its fourth base.
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let reference = &set.allele(1).codes;
        let profile = SubstitutionProfile::build(&scoring, reference);
        let len = core.len() as u32;
        let score: Score = (0..len)
            .map(|k| profile.at(query[5 + k as usize], 3 + k as usize))
            .sum();
        let alignment = Alignment {
            q_start: 5,
            q_end: 5 + len,
            r_start: 3,
            r_end: 3 + len,
            score,
            ops: vec![Op::Match; core.len()],
        };
        let shifted = Placement {
            kind: SegmentKind::J,
            alleles: vec![1],
            shifted: Vec::new(),
            ends: (0, reference.len() as i32),
            clip: ClipIndex::build(&alignment, &query, &profile, scoring.gaps),
            local_score: score,
            evidence: Evidence::default(),
            local_evidence: Evidence::default(),
            anchor: Some(Score::MIN),
            alignment,
        };
        assert_eq!(first.score(), shifted.score());
        for (kind, named) in [
            (SegmentKind::J, vec!["IGHJ4*02"]),
            (SegmentKind::D, vec!["IGHJ4*02", "IGHJ5*02"]),
        ] {
            let placements = [first.clone(), shifted.clone()].map(|mut p| {
                p.kind = kind;
                p
            });
            let called = CalledSegment {
                placement: 0,
                entry: 5,
                exit: 5 + len,
                span: placements[0].clip.full(),
                allele: None,
            };
            let co = co_optimal(
                &placements,
                &query,
                &set,
                &called,
                true,
                0..query.len() as u32,
            );
            assert_eq!(co.tied.len(), 2, "{kind:?}");
            assert_eq!(co.names(&set), named, "{kind:?}");
        }
    }
}
