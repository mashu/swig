//! Placements: aligned segment hits, the nodes of the adjudication DAG.
//!
//! A placement is `(kind, alleles, alignment, clip index, evidence)` — the
//! method note's node `(allele, [a, b), S)` plus everything needed to clip and
//! rescore it in O(1). Alleles whose alignments are column-for-column identical
//! (same query columns, same per-column scores) are merged into one placement,
//! so co-optimal alleles do not multiply DAG work.

use crate::align::Alignment;
use crate::clip::{ClipIndex, ClipSpan};
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
    /// The score a clipped span of this placement must keep to anchor a
    /// call, the anchor gate of its window ([`crate::segment::EvalueGates`]);
    /// `None` when its local alignment never reached it.
    pub anchor: Option<Score>,
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
            anchor: None,
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

    /// Whether `span`, a clip of this placement, anchors a call: the bases
    /// a path keeps are significant on their own. A significant local
    /// alignment is not enough, since the outer flank stitched onto it and
    /// the cuts of a path can leave few of its bases, or junk.
    pub fn anchors(&self, span: &ClipSpan) -> bool {
        self.anchor.is_some_and(|floor| span.score >= floor)
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
