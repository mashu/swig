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
use crate::germline::GermlineSet;
use crate::placement::Placement;

/// The placements tied with a called segment on its interval, and those its
/// posterior counts (see the module documentation).
pub(crate) struct CoOptimal<'a> {
    /// Tied placements, the primary — reported first — at index 0.
    pub tied: Vec<&'a Placement>,
    /// Indices of the placements of the called kind whose probability the
    /// posterior sums: tied with the call, or carried by its primary allele.
    pub counted: Vec<usize>,
}

impl CoOptimal<'_> {
    /// The tied placements' allele names, primary first, each once.
    pub fn names<'s>(&self, set: &'s GermlineSet) -> Vec<&'s str> {
        let mut names: Vec<&str> = Vec::new();
        for p in &self.tied {
            for &allele in &p.alleles {
                for name in set.allele(allele).names() {
                    if !names.contains(&name) {
                        names.push(name);
                    }
                }
            }
        }
        names
    }
}

/// Co-optimal alleles of `called` among `placements` (of every kind; only
/// the called one's count), `set` holding their alleles.
///
/// Tied: worth exactly the call on its interval `[entry, exit)` — an equal
/// score, and under a per-read call (`per_read`) equal matches and
/// mismatches. The primary is the one with the most germline support
/// (unclipped score, the DAG's own tie-break), then the longest exact run
/// and most matches on the interval, then the lowest allele id; a call
/// chosen by its allele's posterior marginal stays first.
///
/// Counted by the posterior: the tied placements, and every placement that
/// explains most of the read bases the call explains (more than half of the
/// called span, within `region`, the read the call's neighbours leave to
/// it) with germline bases there that the primary allele contains. If the
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
    // chosen by its allele's posterior marginal stays first.
    let mut tied: Vec<&Placement> = std::iter::once(called_placement)
        .chain(
            tied_indices
                .iter()
                .filter(|&&i| i != called.placement)
                .map(|&i| &placements[i]),
        )
        .collect();
    tied[usize::from(called.by_marginal)..].sort_by_cached_key(|p| {
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
            std::cmp::Reverse(p.score()),
            std::cmp::Reverse(e.longest_exact_run),
            std::cmp::Reverse(e.matches),
            p.primary_allele(),
        )
    });
    let reference = &set.allele(tied[0].primary_allele()).codes;
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
    CoOptimal { tied, counted }
}
