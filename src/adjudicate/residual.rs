//! The ordinary path of the method note: commit the best V×J pair, then place D
//! inside the residual `[b_V, a_J)` and add `S_D / 2`.
//!
//! Kept as a baseline: it is fast, and it exhibits the documented "V/J eat D"
//! failure mode, which makes it a useful regression reference for the DAG.
//! The V×J step uses the exact optimal cut when V and J overlap.

use super::{Adjudicator, CallPath, CalledSegment, PathScore};
use crate::clip::ClipSpan;
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResidualAdjudicator {
    pub pair_bonus: Score,
    /// V/J pairs overlapping by more query bases are rejected.
    pub max_vj_overlap: u32,
    /// Minimum query bases per clipped segment, indexed by [`SegmentKind::index`].
    pub min_retained: [u32; 3],
}

impl Default for ResidualAdjudicator {
    fn default() -> Self {
        Self {
            pair_bonus: 24,
            max_vj_overlap: 16,
            min_retained: [10, 3, 6],
        }
    }
}

#[derive(Debug, Default)]
pub struct ResidualScratch {
    v: Vec<usize>,
    d: Vec<usize>,
    j: Vec<usize>,
}

#[derive(Clone, Copy)]
struct Pick {
    placement: usize,
    entry: u32,
    exit: u32,
    span: ClipSpan,
}

impl Pick {
    fn segment(self) -> CalledSegment {
        CalledSegment {
            placement: self.placement,
            entry: self.entry,
            exit: self.exit,
            span: self.span,
        }
    }
}

impl ResidualAdjudicator {
    fn keeps(&self, kind: SegmentKind, span: &ClipSpan) -> bool {
        span.q_len() >= self.min_retained[kind.index()].max(1)
    }

    fn clipped(&self, placements: &[Placement], i: usize, entry: u32, exit: u32) -> Option<Pick> {
        let p = &placements[i];
        let span = p.clip.span(entry, exit)?;
        self.keeps(p.kind, &span).then_some(Pick {
            placement: i,
            entry,
            exit,
            span,
        })
    }

    fn full(&self, placements: &[Placement], i: usize) -> Option<Pick> {
        let p = &placements[i];
        self.clipped(placements, i, p.q_start(), p.q_end())
    }

    /// Best V/J pair with the optimal cut through their overlap.
    fn best_pair(
        &self,
        placements: &[Placement],
        v: usize,
        j: usize,
    ) -> Option<(Score, Pick, Pick)> {
        let (pv, pj) = (&placements[v], &placements[j]);
        if pj.q_start() <= pv.q_start() {
            return None;
        }
        if pv.q_end() <= pj.q_start() {
            let (a, b) = (self.full(placements, v)?, self.full(placements, j)?);
            return Some((a.span.score + b.span.score, a, b));
        }
        if pv.q_end() - pj.q_start() > self.max_vj_overlap {
            return None;
        }
        (pj.q_start()..=pv.q_end())
            .filter_map(|c| {
                let a = self.clipped(placements, v, pv.q_start(), c)?;
                let b = self.clipped(placements, j, c, pj.q_end())?;
                Some((a.span.score + b.span.score, a, b))
            })
            .fold(None, |best: Option<(Score, Pick, Pick)>, x| match best {
                Some(b) if b.0 >= x.0 => Some(b),
                _ => Some(x),
            })
    }
}

impl Adjudicator for ResidualAdjudicator {
    type Scratch = ResidualScratch;

    fn adjudicate(
        &self,
        query_len: usize,
        placements: &[Placement],
        scratch: &mut ResidualScratch,
    ) -> Option<CallPath> {
        let ResidualScratch { v, d, j } = scratch;
        for list in [&mut *v, &mut *d, &mut *j] {
            list.clear();
        }
        for (i, p) in placements.iter().enumerate() {
            match p.kind {
                SegmentKind::V => v.push(i),
                SegmentKind::D => d.push(i),
                SegmentKind::J => j.push(i),
            }
        }

        // 1. Commit the flanks: best V×J (+ pair bonus), or a lone V / J.
        let mut best: Option<(Score, Option<Pick>, Option<Pick>)> = None;
        let mut consider = |score: Score, pv: Option<Pick>, pj: Option<Pick>| {
            if best.map_or(true, |(s, _, _)| score > s) {
                best = Some((score, pv, pj));
            }
        };
        for &i in v.iter() {
            if let Some(p) = self.full(placements, i) {
                consider(p.span.score, Some(p), None);
            }
        }
        for &i in j.iter() {
            if let Some(p) = self.full(placements, i) {
                consider(p.span.score, None, Some(p));
            }
        }
        for &vi in v.iter() {
            for &ji in j.iter() {
                if let Some((score, pv, pj)) = self.best_pair(placements, vi, ji) {
                    consider(score + self.pair_bonus, Some(pv), Some(pj));
                }
            }
        }
        let (flank_score, pv, pj) = best?;

        // 2. D inside the residual left by the committed flanks, at half weight.
        let lo = pv.map_or(0, |p| p.span.q_hi);
        let hi = pj.map_or(query_len as u32, |p| p.span.q_lo);
        let pd = d
            .iter()
            .filter_map(|&i| self.clipped(placements, i, lo, hi))
            .fold(None, |best: Option<Pick>, x| match best {
                Some(b) if b.span.score >= x.span.score => Some(b),
                _ => Some(x),
            });

        let d_score = pd.map_or(0, |p| p.span.score.max(0) / 2);
        let segments: Vec<CalledSegment> = [pv, pd, pj]
            .into_iter()
            .flatten()
            .map(Pick::segment)
            .collect();
        Some(CallPath {
            segments,
            score: PathScore {
                score: flank_score + d_score,
                d_count: u8::from(pd.is_some()),
                support: [pv, pd, pj]
                    .into_iter()
                    .flatten()
                    .map(|p| placements[p.placement].score())
                    .sum(),
            },
        })
    }
}
