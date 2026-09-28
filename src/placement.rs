//! Placements: aligned segment hits, the nodes of the adjudication DAG.
//!
//! A placement is `(kind, alleles, alignment, clip index, evidence)` — the
//! method note's node `(allele, [a, b), S)` plus everything needed to clip and
//! rescore it in O(1). Alleles whose alignments are column-for-column identical
//! (same query columns, same per-column scores) are merged into one placement,
//! so co-optimal alleles do not multiply DAG work.

use crate::align::Alignment;
use crate::clip::ClipIndex;
use crate::evidence::Evidence;
use crate::germline::{AlleleId, SegmentKind};
use crate::scoring::Score;
use crate::segment::PreparedSegment;

#[derive(Clone, Debug)]
pub struct Placement {
    pub kind: SegmentKind,
    /// Alleles sharing this exact alignment; the first is the one that was aligned.
    pub alleles: Vec<AlleleId>,
    /// Query coordinates are absolute (whole oriented query).
    pub alignment: Alignment,
    pub clip: ClipIndex,
    /// Evidence of the unclipped alignment.
    pub evidence: Evidence,
    /// The exact ungapped D tract occurs in at least two distinct D sequences.
    /// R-optimized uses this to relax the D-presence penalty. Other segments leave it false.
    pub relaxes_d_cost: bool,
}

impl Placement {
    pub fn new(
        kind: SegmentKind,
        allele: AlleleId,
        alignment: Alignment,
        query: &[u8],
        segment: &PreparedSegment,
    ) -> Self {
        let profile = segment.profile(allele);
        let clip = ClipIndex::build(&alignment, query, profile, segment.gaps());
        let reference = &segment.set().allele(allele).codes;
        let evidence = Evidence::from_columns(alignment.columns(), query, reference);
        Self {
            kind,
            alleles: vec![allele],
            alignment,
            clip,
            evidence,
            relaxes_d_cost: false,
        }
    }

    #[inline]
    pub fn q_start(&self) -> u32 {
        self.alignment.q_start
    }

    #[inline]
    pub fn q_end(&self) -> u32 {
        self.alignment.q_end
    }

    #[inline]
    pub fn score(&self) -> Score {
        self.alignment.score
    }

    #[inline]
    pub fn primary_allele(&self) -> AlleleId {
        self.alleles[0]
    }
}

/// Merge placements in `placements[from..]` whose clip indices are equivalent,
/// keeping the first occurrence (order is preserved).
pub fn merge_identical(placements: &mut Vec<Placement>, from: usize) {
    let mut kept = from;
    for i in from..placements.len() {
        let target = (from..kept).find(|&k| {
            placements[k].kind == placements[i].kind
                && placements[k].clip.equivalent(&placements[i].clip)
        });
        match target {
            Some(k) => {
                let extra = std::mem::take(&mut placements[i].alleles);
                for allele in extra {
                    if !placements[k].alleles.contains(&allele) {
                        placements[k].alleles.push(allele);
                    }
                }
            }
            None => {
                placements.swap(kept, i);
                kept += 1;
            }
        }
    }
    placements.truncate(kept);
}

/// Mark D placements whose exact ungapped tract is shared by two distinct D sequences.
///
/// `references[allele]` is that allele's codes. Identical sequences count once,
/// matching Swig's template-class relaxation.
pub fn mark_repeated_d_tracts(placements: &mut [Placement], query: &[u8], references: &[&[u8]]) {
    for placement in placements.iter_mut().filter(|p| p.kind == SegmentKind::D) {
        let Some(reference) = references.get(placement.primary_allele() as usize).copied() else {
            continue;
        };
        if !exact_canonical_tract(placement, query, reference) {
            continue;
        }
        let tract = &query[placement.q_start() as usize..placement.q_end() as usize];
        if tract.is_empty() {
            continue;
        }
        let mut classes = 0u32;
        let mut seen: Vec<&[u8]> = Vec::new();
        for candidate in references {
            if candidate.windows(tract.len()).any(|window| window == tract)
                && !seen.iter().any(|prior| *prior == *candidate)
            {
                seen.push(candidate);
                classes += 1;
                if classes >= 2 {
                    placement.relaxes_d_cost = true;
                    break;
                }
            }
        }
    }
}

fn exact_canonical_tract(placement: &Placement, query: &[u8], reference: &[u8]) -> bool {
    placement
        .alignment
        .ops
        .iter()
        .all(|op| *op == crate::align::Op::Match)
        && placement.alignment.columns().all(|(_, q, r)| {
            let (a, b) = (query[q as usize], reference[r as usize]);
            a < crate::alphabet::N && a == b
        })
}
