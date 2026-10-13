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
    /// The merged alleles whose germline coordinates differ from the
    /// first's, with the difference ([`Placement::germline_shift`]).
    pub shifted: Vec<(AlleleId, i32)>,
    /// Query coordinates are absolute (whole oriented query).
    pub alignment: Alignment,
    pub clip: ClipIndex,
    /// Evidence of the unclipped alignment.
    pub evidence: Evidence,
    /// Evidence of the local alignment, before the stage stitched its outer
    /// flank on ([`crate::flank`]): what the identity floor judges
    /// ([`crate::pipeline::AssignerOptions::min_identity`]).
    pub local_evidence: Evidence,
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
            shifted: Vec::new(),
            alignment,
            clip,
            evidence,
            local_evidence: evidence,
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

    /// Germline coordinate of `allele`, one of [`Placement::alleles`], minus
    /// the first's. Merged alignments share their query columns, operations
    /// and per-column scores, so an allele of another length before the
    /// aligned bases has the first's coordinates shifted by a constant.
    pub fn germline_shift(&self, allele: AlleleId) -> i32 {
        self.shifted
            .iter()
            .find(|&&(merged, _)| merged == allele)
            .map_or(0, |&(_, shift)| shift)
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
/// keeping the first occurrence (order is preserved) and the germline shift
/// of every allele merged into it.
pub fn merge_identical(placements: &mut Vec<Placement>, from: usize) {
    let mut kept = from;
    for i in from..placements.len() {
        let target = (from..kept).find(|&k| {
            placements[k].kind == placements[i].kind
                && placements[k].clip.equivalent(&placements[i].clip)
        });
        match target {
            Some(k) => {
                let start = |p: &Placement| i64::from(p.alignment.r_start);
                let shift = (start(&placements[i]) - start(&placements[k])) as i32;
                let merged = std::mem::take(&mut placements[i].alleles);
                for allele in merged {
                    if !placements[k].alleles.contains(&allele) {
                        let own = shift + placements[i].germline_shift(allele);
                        placements[k].alleles.push(allele);
                        if own != 0 {
                            placements[k].shifted.push((allele, own));
                        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::testing::placement;

    /// An allele whose germline starts 3 bases earlier aligns the same read
    /// bases column for column: merged, its coordinates are the first's
    /// shifted by 3, also when it reaches the first through another merge.
    #[test]
    fn merged_alleles_keep_their_germline_shift() {
        let query = encode(b"ACGTTGCAAGTCGGAT");
        let at = |allele: AlleleId, r_start: u32| {
            let mut p = placement(SegmentKind::D, &query, &query[2..12], 2);
            p.alleles = vec![allele];
            p.alignment.r_start += r_start;
            p.alignment.r_end += r_start;
            p
        };
        let mut merged = vec![at(0, 0), at(1, 3), at(2, 0)];
        merge_identical(&mut merged, 0);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].alleles, [0, 1, 2]);
        assert_eq!([0, 1, 2].map(|a| merged[0].germline_shift(a)), [0, 3, 0]);
        // A placement merged earlier brings its own alleles' shifts along.
        let mut later = vec![at(5, 1), merged.remove(0)];
        merge_identical(&mut later, 0);
        assert_eq!(later[0].alleles, [5, 0, 1, 2]);
        assert_eq!(
            [5, 0, 1, 2].map(|a| later[0].germline_shift(a)),
            [0, -1, 2, -1]
        );
    }
}
