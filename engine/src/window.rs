//! Where to search for D: junction windows derived from V/J geometry.
//!
//! D is never searched on the whole query (the note: no three-way product
//! search). A [`DWindowPolicy`] proposes one or more query intervals from the
//! V and J placements that bound the D ([`bounding_v`]); the pipeline merges
//! overlapping windows so a D placement is never found twice.

use std::ops::Range;

use crate::germline::GermlineSet;
use crate::placement::Placement;

pub trait DWindowPolicy: Send + Sync {
    /// Append to `out` the windows of a D between the V placement `v` and
    /// the J placement `j` that bound it; `None`: no V (J) bounds it, and
    /// the D may lie from the read's start (up to its end).
    fn windows(
        &self,
        query_len: usize,
        v: Option<&Placement>,
        j: Option<&Placement>,
        out: &mut Vec<Range<usize>>,
    );
}

/// The V placement that bounds a read's D windows, among its V placements
/// `v` (best first): the best one, or, when V hits can anchor a call
/// (`anchoring`: the V stage has E-value gates), the best anchor
/// ([`Placement::anchor`]). A V hit that is not significant on its own may
/// be a chance hit downstream of the real D, in a read that starts inside
/// CDR3, as for the J window ([`DownstreamWindow`]); without an anchor, D
/// is searched as on a read without V.
pub fn bounding_v(v: &[Placement], anchoring: bool) -> Option<&Placement> {
    if anchoring {
        v.iter().find(|p| p.anchor.is_some())
    } else {
        v.first()
    }
}

/// `[b_V − margin, a_J + margin)` of the bounding V and J, clamped to the
/// read and to V's start; empty (never inverted) when V reaches past the J's
/// start by more than the margin.
fn residual(
    query_len: usize,
    v: Option<&Placement>,
    j: Option<&Placement>,
    margin: u32,
) -> Range<usize> {
    let floor = v.map_or(0, |p| p.q_start() as usize);
    let start = v
        .map_or(0, |p| p.q_end().saturating_sub(margin) as usize)
        .max(floor);
    let end = j
        .map_or(query_len, |p| (p.q_start() + margin) as usize)
        .min(query_len);
    start..end.max(start)
}

/// The ordinary residual `[b_V − margin, a_J + margin)` of the bounding V and J.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResidualWindow {
    pub overlap_margin: u32,
}

impl DWindowPolicy for ResidualWindow {
    fn windows(
        &self,
        query_len: usize,
        v: Option<&Placement>,
        j: Option<&Placement>,
        out: &mut Vec<Range<usize>>,
    ) {
        let window = residual(query_len, v, j, self.overlap_margin);
        if !window.is_empty() {
            out.push(window);
        }
    }
}

/// Residual widened by `reclaim` bases into V and J, searched when it holds
/// at least `min_width` bases (and at least one).
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
    /// A V or J alignment reaches into the junction only through chance
    /// matches, a few bases; 12 covers that with margin. A wider window
    /// hurts as well as costs: a D alignment free to stretch into a V or J
    /// tail can trade its junction core for chance matches that the cut
    /// then removes.
    fn default() -> Self {
        Self {
            reclaim: 12,
            min_width: 4,
        }
    }
}

impl DWindowPolicy for ExpandedResidual {
    fn windows(
        &self,
        query_len: usize,
        v: Option<&Placement>,
        j: Option<&Placement>,
        out: &mut Vec<Range<usize>>,
    ) {
        let window = residual(query_len, v, j, self.reclaim);
        if window.len() >= self.min_width.max(1) {
            out.push(window);
        }
    }
}

/// Adds the conserved-anchor window `[q(V Cys) + 3, q(J Trp/Phe))` — the
/// CDR3 interior — as a second location hypothesis for D. Anchors are
/// projected through the bounding V and J alignments (extrapolated along
/// their diagonals when the anchor lies outside the aligned part).
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
        v: Option<&Placement>,
        j: Option<&Placement>,
        out: &mut Vec<Range<usize>>,
    ) {
        self.inner.windows(query_len, v, j, out);
        let (Some(vb), Some(jb)) = (v, j) else {
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

/// Where J is searched: right of V, reaching `reclaim` bases back into it.
///
/// J follows V on every path through a V, and a J alignment can only reach
/// into V by chance matches at the junction edge, which [`ExpandedResidual`]
/// bounds by the same margin. When some V placement is significant on its
/// own ([`Placement::anchor`]), the window starts `reclaim` bases before the
/// earliest V end, so it contains every J that a path through any V
/// placement can use up to that overlap. Otherwise the read may hold no V at
/// all (it starts inside CDR3, or is a J–C fragment), and its V placements
/// are chance hits that can lie downstream of the real J, in its 3′ flank or
/// C region; a path may skip V, so the whole read is searched. So it is
/// when no hit can anchor (a V stage without E-value gates).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DownstreamWindow {
    pub reclaim: u32,
}

impl Default for DownstreamWindow {
    fn default() -> Self {
        Self {
            reclaim: ExpandedResidual::default().reclaim,
        }
    }
}

impl DownstreamWindow {
    pub fn window(&self, query_len: usize, upstream: &[Placement]) -> Range<usize> {
        if !upstream.iter().any(|p| p.anchor.is_some()) {
            return 0..query_len;
        }
        let start = upstream
            .iter()
            .map(|p| p.q_end().saturating_sub(self.reclaim) as usize)
            .min()
            .unwrap_or(0);
        start.min(query_len)..query_len
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
    use crate::germline::SegmentKind;

    #[test]
    fn merges_overlapping_windows() {
        let mut w = vec![10..20, 0..5, 18..30, 5..7];
        merge_windows(&mut w);
        assert_eq!(w, vec![0..7, 10..30]);
        let mut empty: Vec<Range<usize>> = Vec::new();
        merge_windows(&mut empty);
        assert!(empty.is_empty());
        let mut nested = vec![0..100, 10..20, 30..40];
        merge_windows(&mut nested);
        assert_eq!(nested, vec![0..100]);
    }

    #[test]
    fn j_window_without_v_or_past_the_read() {
        let window = DownstreamWindow { reclaim: 12 };
        assert_eq!(window.window(80, &[]), 0..80);
        assert_eq!(window.window(0, &[]), 0..0);
    }

    /// A V hit that is not significant on its own may be a chance hit
    /// downstream of the real J (a read that starts inside CDR3): only a
    /// significant V bounds where J can start.
    #[test]
    fn j_window_starts_after_significant_v_only() {
        let mut next = crate::testing::rng(4);
        let query = crate::testing::random_codes(&mut next, 218);
        let mut v = crate::testing::placement(SegmentKind::V, &query, &query[..150], 0);
        let window = DownstreamWindow { reclaim: 12 };
        v.anchor = None;
        assert_eq!(window.window(218, std::slice::from_ref(&v)), 0..218);
        let mut anchor = v.clone();
        anchor.anchor = Some(0);
        assert_eq!(window.window(218, &[anchor, v]), 138..218);
    }

    /// A read that starts inside CDR3 (a D at 5..29, then its J at 36..96)
    /// may hold a chance V hit downstream of its J. That hit is not
    /// significant on its own, so it does not bound the D, which is searched
    /// from the read's start; a V stage that cannot anchor bounds D by its
    /// best hit.
    #[test]
    fn d_window_starts_after_significant_v_only() {
        let mut next = crate::testing::rng(8);
        let query = crate::testing::random_codes(&mut next, 192);
        let mut chance = crate::testing::placement(SegmentKind::V, &query, &query[81..117], 81);
        chance.anchor = None;
        let j = crate::testing::placement(SegmentKind::J, &query, &query[36..96], 36);
        let policy = WithAnchorWindow::disabled(ExpandedResidual::default());
        let windows = |v: Option<&Placement>| {
            let mut out = Vec::new();
            policy.windows(192, v, Some(&j), &mut out);
            out
        };
        let v = std::slice::from_ref(&chance);
        assert!(bounding_v(v, true).is_none());
        let d = 5..29;
        let found = windows(bounding_v(v, true));
        assert!(
            found.iter().any(|w| w.start <= d.start && d.end <= w.end),
            "{found:?}"
        );
        assert_eq!(windows(bounding_v(v, false)), Vec::<Range<usize>>::new());
        let mut anchor = crate::testing::placement(SegmentKind::V, &query, &query[..30], 0);
        anchor.anchor = Some(0);
        let v = [chance.clone(), anchor];
        assert_eq!(bounding_v(&v, true).map(Placement::q_end), Some(30));
        assert_eq!(bounding_v(&v, false).map(Placement::q_end), Some(117));
    }

    /// A V that reaches past the J's start by more than the margin makes the
    /// residual empty; no policy may propose an inverted window.
    #[test]
    fn crossing_v_and_j_propose_no_window() {
        let mut next = crate::testing::rng(6);
        let query = crate::testing::random_codes(&mut next, 100);
        let v = crate::testing::placement(SegmentKind::V, &query, &query[..80], 0);
        let j = crate::testing::placement(SegmentKind::J, &query, &query[40..60], 40);
        for min_width in [0, 1, 4] {
            let mut out = Vec::new();
            ExpandedResidual {
                reclaim: 12,
                min_width,
            }
            .windows(100, Some(&v), Some(&j), &mut out);
            assert!(
                out.iter().all(|w| w.start < w.end),
                "min_width {min_width}: {out:?}"
            );
        }
    }
}
