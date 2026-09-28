//! Where to search for D: junction windows derived from V/J geometry.
//!
//! D is never searched on the whole query (the note: no three-way product
//! search). A [`DWindowPolicy`] proposes one or more query intervals from the
//! current V and J placements; the pipeline merges overlapping windows so a D
//! placement is never found twice.

use std::ops::Range;

use crate::germline::GermlineSet;
use crate::placement::Placement;

pub trait DWindowPolicy: Send + Sync {
    /// `v` and `j` are sorted best-first. Append windows to `out`.
    fn windows(
        &self,
        query_len: usize,
        v: &[Placement],
        j: &[Placement],
        out: &mut Vec<Range<usize>>,
    );
}

fn residual(query_len: usize, v: &[Placement], j: &[Placement], margin: u32) -> Range<usize> {
    let v_best = v.first();
    let j_best = j.first();
    let floor = v_best.map_or(0, |p| p.q_start() as usize);
    let start = v_best
        .map_or(0, |p| p.q_end().saturating_sub(margin) as usize)
        .max(floor);
    let end = j_best
        .map_or(query_len, |p| (p.q_start() + margin) as usize)
        .min(query_len);
    start..end
}

/// The ordinary residual `[b_V − margin, a_J + margin)` of the best V and J.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResidualWindow {
    pub overlap_margin: u32,
}

impl DWindowPolicy for ResidualWindow {
    fn windows(
        &self,
        query_len: usize,
        v: &[Placement],
        j: &[Placement],
        out: &mut Vec<Range<usize>>,
    ) {
        let window = residual(query_len, v, j, self.overlap_margin);
        if !window.is_empty() {
            out.push(window);
        }
    }
}

/// Residual widened by `reclaim` bases into V and J.
///
/// This replaces the note's gated *boundary rescue*: D hits that overlap an
/// over-extended V or J are found up front, and the adjudicator's exact cut
/// optimisation decides — with the D-presence prior and the strict "with-D must
/// beat no-D" tie rule — whether reclaiming those bases improves the joint
/// score. No gap detector, no separate re-alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExpandedResidual {
    pub reclaim: u32,
    pub min_width: usize,
}

impl Default for ExpandedResidual {
    fn default() -> Self {
        Self {
            reclaim: 24,
            min_width: 4,
        }
    }
}

impl DWindowPolicy for ExpandedResidual {
    fn windows(
        &self,
        query_len: usize,
        v: &[Placement],
        j: &[Placement],
        out: &mut Vec<Range<usize>>,
    ) {
        let window = residual(query_len, v, j, self.reclaim);
        if window.len() >= self.min_width {
            out.push(window);
        }
    }
}

/// Adds the conserved-anchor window `[q(V Cys) + 3, q(J Trp/Phe))` — the
/// CDR3 interior — as a second location hypothesis (AER-R). Anchors are
/// projected through the best V and J alignments (extrapolated along their
/// diagonals when the anchor lies outside the aligned part).
pub struct WithAnchorWindow<P> {
    pub inner: P,
    v_anchor: Vec<Option<u32>>,
    j_anchor: Vec<Option<u32>>,
}

impl<P> WithAnchorWindow<P> {
    pub fn new(inner: P, v: &GermlineSet, j: &GermlineSet) -> Self {
        Self {
            inner,
            v_anchor: v.alleles().iter().map(|a| a.anchor).collect(),
            j_anchor: j.alleles().iter().map(|a| a.anchor).collect(),
        }
    }

    /// Same policy with no anchor window (anchors ignored).
    pub fn disabled(inner: P) -> Self {
        Self {
            inner,
            v_anchor: Vec::new(),
            j_anchor: Vec::new(),
        }
    }
}

impl<P: DWindowPolicy> DWindowPolicy for WithAnchorWindow<P> {
    fn windows(
        &self,
        query_len: usize,
        v: &[Placement],
        j: &[Placement],
        out: &mut Vec<Range<usize>>,
    ) {
        self.inner.windows(query_len, v, j, out);
        let (Some(vb), Some(jb)) = (v.first(), j.first()) else {
            return;
        };
        let v_anchor = self
            .v_anchor
            .get(vb.primary_allele() as usize)
            .copied()
            .flatten();
        let j_anchor = self
            .j_anchor
            .get(jb.primary_allele() as usize)
            .copied()
            .flatten();
        let (Some(va), Some(ja)) = (v_anchor, j_anchor) else {
            return;
        };
        let start = (vb.clip.project_reference(va) + 3).clamp(0, query_len as i64) as usize;
        let end = jb.clip.project_reference(ja).clamp(0, query_len as i64) as usize;
        if start < end {
            out.push(start..end);
        }
    }
}

/// Sort and merge overlapping or touching windows in place.
pub fn merge_windows(windows: &mut Vec<Range<usize>>) {
    windows.sort_by_key(|w| (w.start, w.end));
    let mut kept = 0;
    for i in 0..windows.len() {
        if kept > 0 && windows[i].start <= windows[kept - 1].end {
            windows[kept - 1].end = windows[kept - 1].end.max(windows[i].end);
        } else {
            windows[kept] = windows[i].clone();
            kept += 1;
        }
    }
    windows.truncate(kept);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_overlapping_windows() {
        let mut w = vec![10..20, 0..5, 18..30, 5..7];
        merge_windows(&mut w);
        assert_eq!(w, vec![0..7, 10..30]);
    }
}
