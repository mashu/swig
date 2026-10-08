//! Sibling bounds: skip DP for alleles of a gene that cannot enter the
//! top-N local scores.
//!
//! The alleles of a gene share a *frame* (see [`GeneFrames`]): frame
//! position `f = j + o_X` for position `j` of allele `X` at offset `o_X`, so
//! `X` covers `span(X) = [o_X, o_X + |X|)`. With `s_X(q, j)` the profile
//! score of query code `q` at position `j`, the bound of a sibling `A` given
//! a representative `R` of the same frame is
//!
//! ```text
//! Δ(R, A) = Σ_{f ∈ span(A) ∩ span(R)} max(0, max_q [s_A(q, f − o_A) − s_R(q, f − o_R)])
//!         + Σ_{f ∈ span(A) \ span(R)} max(0, max_q s_A(q, f − o_A))
//! ```
//!
//! Positions only `R` covers add nothing; for two alleles of equal span it
//! is the per-position maximum profile gain.
//!
//! **Theorem U.** When the two bands are the same frame diagonals
//! ([`FamilyBoundTable::frame_band`]), `S*_{B_A}(A) ≤ S*_{B_R}(R) + Δ(R, A)`.
//!
//! *Proof.* Split an optimal path of `A` into the columns on positions only
//! `A` covers (5′ and 3′ of `span(R)`) and the middle part over the overlap.
//! Strip the middle part's leading and trailing gap columns (each scores
//! ≤ 0). What remains starts and ends with a match column, so every gap
//! column in it scores as before, and shifted into `R`'s coordinates it lies
//! inside `R`'s matrix and band: a local path of `R`, which scores at most
//! `S*(R)` with `R`'s profile and at most that plus the overlap gains with
//! `A`'s. In the stripped and overhanging columns a substitution scores at
//! most `max_q s_A` and a gap ≤ 0. An empty middle part needs only
//! `S*(R) ≥ 0`. ∎
//!
//! The proof uses nothing about the offsets but that the bands and `Δ` use
//! the same placement: the V table takes it from the retriever's own frames,
//! whose bands are shifted by each allele's offset. A sibling in another
//! frame, or with an inconsistent band, has no bound and is always scored.
//!
//! **Theorem L.** Let `π` be a traced local path of `R` in `B_R` (any path,
//! not only an optimal one), `π_A` its columns from the first to the last
//! match column on a position of `span(A)`, and `D(R, A)` the positions of
//! `span(A) ∩ span(R)` whose profile columns differ. With the bands the
//! same frame diagonals,
//!
//! ```text
//! L(π, A) = max(0, S_R(π_A) + Σ_{f ∈ D(R, A) matched to q in π_A} [s_A(q, f − o_A) − s_R(q, f − o_R)])
//!         ≤ S*_{B_A}(A)
//! ```
//!
//! *Proof.* Every position between the two end columns lies in `span(A)`
//! (an interval) and in `span(R)` (`π` is a path of `R`), so `π_A` shifted
//! into `A`'s coordinates lies inside `A`'s matrix, and inside `B_A` by
//! the band condition. It starts with a match column, so every gap column
//! follows the same column as in `π` and scores alike. Its substitution
//! columns score the same with either profile off `D`. So `L` before the
//! `max` is the score of a local path of `A`, and a local optimum is at
//! least that and at least 0. ∎
//!
//! [`FramePath`] stores the prefix sums of `π`'s column scores, so `L`
//! costs `O(|D|)`: the clipped sum plus a correction at the differing
//! positions. The column scores are the DP's when `open ≤ extend` (a gap
//! column extends a gap of its own kind and opens one otherwise); with
//! `open > extend` the DP can score a path higher than this sum, and `L`
//! stays a lower bound.
//!
//! `L ≤ S*(A) ≤ U` together pin the exact score when they meet, and the
//! `top`-th largest of the lower bounds of a request's candidates — one
//! value per candidate — never exceeds the `top`-th largest exact score:
//! the cutoff that prunes siblings in [`crate::segment`].
//!
//! `Δ` and `D` depend only on the profiles, so they are computed once per
//! same-frame pair when the stage is built, both orders in one pass.

use crate::align::{Alignment, Op};
use crate::alphabet::ALPHABET_SIZE;
use crate::band::Band;
use crate::germline::AlleleId;
use crate::retrieve::{Candidate, GeneFrames};
use crate::scoring::{GapCosts, Score, SubstitutionProfile};
use crate::segment::PreparedSegment;

/// `Δ(R, A)` and `D(R, A)` of every ordered pair of alleles that share a frame.
#[derive(Clone, Debug, Default)]
pub struct FamilyBoundTable {
    /// `(frame, offset)` of each allele.
    placement: Vec<(u32, u32)>,
    /// Position of each allele among its frame's members.
    rank: Vec<u32>,
    /// Per frame: its first entry in `delta` and `diffs` and its member
    /// count `n`; the pair of member ranks `(r, a)` is entry `first + r·n + a`.
    frames: Vec<(usize, usize)>,
    delta: Vec<Score>,
    /// `D(R, A)` as a range of `diff_positions` (the same for both orders).
    diffs: Vec<(u32, u32)>,
    /// Frame positions, ascending within each pair's range.
    diff_positions: Vec<u32>,
}

impl FamilyBoundTable {
    /// Alleles placed in the gene frames whose diagonals their bands come
    /// from (the V retriever's [`GeneFrames`]).
    pub fn with_frames(segment: &PreparedSegment, frames: &GeneFrames) -> Self {
        let placement = (0..segment.set().len())
            .map(|allele| {
                let allele = AlleleId::try_from(allele).expect("allele count fits u32");
                (frames.frame(allele), frames.offset(allele))
            })
            .collect();
        Self::from_placement(segment, placement)
    }

    /// One frame per gene, every allele at offset 0. Any placement is
    /// consistent with unbanded candidates, so this serves a stage whose
    /// bands are all [`Band::Full`] (J).
    pub fn build(segment: &PreparedSegment) -> Self {
        let placement = segment
            .set()
            .alleles()
            .iter()
            .map(|a| (a.gene, 0))
            .collect();
        Self::from_placement(segment, placement)
    }

    /// The table of `(frame, offset)` per allele; frame ids are dense.
    fn from_placement(segment: &PreparedSegment, placement: Vec<(u32, u32)>) -> Self {
        let frame_count = placement.iter().map(|&(f, _)| f as usize + 1).max();
        let mut members: Vec<Vec<AlleleId>> = vec![Vec::new(); frame_count.unwrap_or(0)];
        let mut rank = vec![0; placement.len()];
        for (allele, &(frame, _)) in (0..).zip(&placement) {
            let list = &mut members[frame as usize];
            rank[allele as usize] = u32::try_from(list.len()).expect("allele count fits u32");
            list.push(allele);
        }
        let covers: Vec<Vec<Score>> = (0..placement.len())
            .map(|allele| cover(segment.profile(allele as AlleleId)))
            .collect();
        let placed = |allele: AlleleId| Placed {
            profile: segment.profile(allele),
            offset: placement[allele as usize].1 as usize,
            cover: &covers[allele as usize],
        };
        let mut frames = Vec::with_capacity(members.len());
        let (mut delta, mut diffs, mut diff_positions) = (Vec::new(), Vec::new(), Vec::new());
        let mut gains = Gains::default();
        for list in &members {
            let (first, n) = (delta.len(), list.len());
            frames.push((first, n));
            delta.resize(first + n * n, 0);
            diffs.resize(first + n * n, (0, 0));
            for (i, &r) in list.iter().enumerate() {
                for (j, &a) in list.iter().enumerate().skip(i + 1) {
                    let start = diff_positions.len();
                    let (r_to_a, a_to_r) = gains.deltas(placed(r), placed(a), &mut diff_positions);
                    delta[first + i * n + j] = r_to_a;
                    delta[first + j * n + i] = a_to_r;
                    let [start, end] = [start, diff_positions.len()]
                        .map(|p| u32::try_from(p).expect("diff positions fit u32"));
                    diffs[first + i * n + j] = (start, end);
                    diffs[first + j * n + i] = (start, end);
                }
            }
        }
        Self {
            placement,
            rank,
            frames,
            delta,
            diffs,
            diff_positions,
        }
    }

    /// Entry of the ordered pair in `delta` and `diffs`; `None` across frames.
    #[inline]
    fn pair(&self, representative: AlleleId, sibling: AlleleId) -> Option<usize> {
        let (frame, _) = self.placement[representative as usize];
        (frame == self.placement[sibling as usize].0).then(|| {
            let (first, n) = self.frames[frame as usize];
            let (r, a) = (
                self.rank[representative as usize] as usize,
                self.rank[sibling as usize] as usize,
            );
            first + r * n + a
        })
    }

    /// `Δ(R, A)` for two alleles of one frame; `None` across frames.
    pub fn delta(&self, representative: AlleleId, sibling: AlleleId) -> Option<Score> {
        self.pair(representative, sibling)
            .map(|pair| self.delta[pair])
    }

    /// Whether the two candidates' bands are the same set of frame
    /// diagonals: the alleles share a frame, and both bands are unbanded, or
    /// both are diagonal ranges shifted by the offset difference
    /// (`lo_A − lo_R = hi_A − hi_R = o_A − o_R`). A frame diagonal `i − f`
    /// is the allele diagonal `i − j` minus the allele's offset.
    pub fn frame_band(&self, representative: &Candidate, sibling: &Candidate) -> bool {
        let (frame_r, offset_r) = self.placement[representative.allele as usize];
        let (frame_a, offset_a) = self.placement[sibling.allele as usize];
        if frame_r != frame_a {
            return false;
        }
        match (representative.band, sibling.band) {
            (Band::Full, Band::Full) => true,
            (Band::Diagonals { lo: lo_r, hi: hi_r }, Band::Diagonals { lo: lo_a, hi: hi_a }) => {
                let shift = i64::from(offset_a) - i64::from(offset_r);
                i64::from(lo_a) - i64::from(lo_r) == shift
                    && i64::from(hi_a) - i64::from(hi_r) == shift
            }
            _ => false,
        }
    }

    /// Theorem U: `S*(R) + Δ(R, A)` bounds the sibling's local score in its
    /// band, given `score`, the representative's in its band, when
    /// [`Self::frame_band`] holds; otherwise `None`.
    pub fn upper_bound(
        &self,
        representative: &Candidate,
        sibling: &Candidate,
        score: Score,
    ) -> Option<Score> {
        if !self.frame_band(representative, sibling) {
            return None;
        }
        self.delta(representative.allele, sibling.allele)
            .map(|delta| score.saturating_add(delta))
    }

    /// Index `alignment`, a traced path of `representative` against `query`
    /// with `profile` and `gaps`, by frame position into `out`.
    pub fn frame_path(
        &self,
        representative: AlleleId,
        alignment: &Alignment,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
        out: &mut FramePath,
    ) {
        out.first = alignment.r_start + self.placement[representative as usize].1;
        out.matched.clear();
        out.column.clear();
        out.prefix.clear();
        out.prefix.push(0);
        let (mut q, mut r) = (alignment.q_start as usize, alignment.r_start as usize);
        let (mut total, mut previous) = (0, None);
        for (column, &op) in (0..).zip(&alignment.ops) {
            let gap = if previous == Some(op) {
                gaps.extend
            } else {
                gaps.open
            };
            total += match op {
                Op::Match => {
                    let code = query[q];
                    out.matched.push(code);
                    out.column.push(column);
                    (q, r) = (q + 1, r + 1);
                    profile.at(code, r - 1)
                }
                Op::Del => {
                    out.matched.push(DELETED);
                    out.column.push(column);
                    r += 1;
                    gap
                }
                Op::Ins => {
                    q += 1;
                    gap
                }
            };
            out.prefix.push(total);
            previous = Some(op);
        }
        debug_assert!(
            gaps.open > gaps.extend || total == alignment.score,
            "path columns sum to {total}, the alignment scores {}",
            alignment.score
        );
    }

    /// Theorem L: a lower bound on the sibling's local score in its band from
    /// `path`, a traced path of the representative, when [`Self::frame_band`]
    /// holds; otherwise `None`. Costs `O(|D(R, A)|)` past the clip ends.
    pub fn path_lower_bound(
        &self,
        path: &FramePath,
        representative: &Candidate,
        sibling: &Candidate,
        segment: &PreparedSegment,
    ) -> Option<Score> {
        if !self.frame_band(representative, sibling) {
            return None;
        }
        let (r, a) = (representative.allele, sibling.allele);
        let (offset_r, offset_a) = (
            self.placement[r as usize].1 as usize,
            self.placement[a as usize].1 as usize,
        );
        let profile_a = segment.profile(a);
        // Path positions (frame position − `first`) inside span(A).
        let first = path.first as usize;
        let lo = offset_a.saturating_sub(first);
        let hi = (offset_a + profile_a.len())
            .saturating_sub(first)
            .min(path.matched.len());
        let matched = |&t: &usize| path.matched[t] != DELETED;
        let Some(start) = (lo..hi).find(matched) else {
            return Some(0);
        };
        let end = (start..hi).rfind(matched).expect("start is matched");
        let mut score =
            path.prefix[path.column[end] as usize + 1] - path.prefix[path.column[start] as usize];
        let profile_r = segment.profile(r);
        let (from, to) = self.diffs[self.pair(r, a).expect("frame_band holds")];
        for &f in &self.diff_positions[from as usize..to as usize] {
            let f = f as usize;
            let Some(t) = f.checked_sub(first).filter(|&t| t >= start) else {
                continue;
            };
            if t > end {
                break;
            }
            let code = path.matched[t];
            if code != DELETED {
                score += profile_a.at(code, f - offset_a) - profile_r.at(code, f - offset_r);
            }
        }
        Some(score.max(0))
    }
}

/// [`FramePath::matched`] at a position the path deletes.
const DELETED: u8 = u8::MAX;

/// A traced path of a representative indexed by frame position, for
/// [`FamilyBoundTable::path_lower_bound`]; reused across paths.
#[derive(Clone, Debug, Default)]
pub struct FramePath {
    /// Frame position of the path's first reference base.
    first: u32,
    /// Per frame position from `first`: the query code the path matches
    /// there, or [`DELETED`].
    matched: Vec<u8>,
    /// Per frame position from `first`: the path column that consumes it.
    column: Vec<u32>,
    /// `prefix[k]`: the representative's score of the path's columns `..k`.
    prefix: Vec<Score>,
}

/// `cover[j] = Σ_{t < j} max(0, max_q s(q, t))`: the most a stretch of
/// substitution columns over positions `[a, b)` can score is
/// `cover[b] − cover[a]`.
fn cover(profile: &SubstitutionProfile) -> Vec<Score> {
    let mut best = vec![0; profile.len()];
    for code in 0..ALPHABET_SIZE as u8 {
        for (slot, &s) in best.iter_mut().zip(profile.row(code)) {
            *slot = (*slot).max(s);
        }
    }
    let mut prefix = Vec::with_capacity(best.len() + 1);
    prefix.push(0);
    let mut total = 0;
    for b in best {
        total += b;
        prefix.push(total);
    }
    prefix
}

/// An allele in its frame.
#[derive(Clone, Copy)]
struct Placed<'a> {
    profile: &'a SubstitutionProfile,
    offset: usize,
    /// [`cover`] of the profile.
    cover: &'a [Score],
}

impl Placed<'_> {
    fn end(&self) -> usize {
        self.offset + self.profile.len()
    }

    /// The most substitution columns can score on the allele's positions
    /// outside the frame interval `[lo, hi)`.
    fn outside(&self, lo: usize, hi: usize) -> Score {
        let len = self.profile.len();
        let lo = lo.clamp(self.offset, self.end()) - self.offset;
        let hi = hi.clamp(self.offset, self.end()) - self.offset;
        self.cover[lo] + self.cover[len] - self.cover[hi.max(lo)]
    }
}

/// Per-position gains of a pair over its overlap, reused across pairs.
#[derive(Default)]
struct Gains {
    forward: Vec<Score>,
    backward: Vec<Score>,
}

impl Gains {
    /// `(Δ(R, A), Δ(A, R))`, appending `D(R, A)` to `diffs`. One pass over
    /// the five profile rows of the overlap computes `max(0, max_q d)` and
    /// `max(0, max_q −d)` for `d = s_A − s_R` together; a column differs
    /// where either is positive. The overhangs come from [`cover`].
    fn deltas(&mut self, r: Placed<'_>, a: Placed<'_>, diffs: &mut Vec<u32>) -> (Score, Score) {
        let lo = r.offset.max(a.offset);
        let hi = r.end().min(a.end()).max(lo);
        let width = hi - lo;
        self.forward.clear();
        self.forward.resize(width, 0);
        self.backward.clear();
        self.backward.resize(width, 0);
        for code in 0..ALPHABET_SIZE as u8 {
            let from = &r.profile.row(code)[lo - r.offset..hi - r.offset];
            let to = &a.profile.row(code)[lo - a.offset..hi - a.offset];
            for (((f, b), &s_r), &s_a) in self
                .forward
                .iter_mut()
                .zip(self.backward.iter_mut())
                .zip(from)
                .zip(to)
            {
                *f = (*f).max(s_a - s_r);
                *b = (*b).max(s_r - s_a);
            }
        }
        diffs.extend(
            (lo..hi)
                .zip(self.forward.iter().zip(&self.backward))
                .filter(|&(_, (&f, &b))| f > 0 || b > 0)
                .map(|(position, _)| u32::try_from(position).expect("frame positions fit u32")),
        );
        (
            self.forward.iter().sum::<Score>() + a.outside(lo, hi),
            self.backward.iter().sum::<Score>() + r.outside(lo, hi),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{Aligner, BandedGotoh};
    use crate::alphabet::{decode, encode};
    use crate::germline::{GermlineSet, SegmentKind};
    use crate::scoring::AffineScoring;

    const SCORING: AffineScoring = AffineScoring::new(2, -3, -9, -1);

    /// A frame of `alleles` (`(sequence, offset)`), all of gene IGHV1-1.
    fn frame_segment(alleles: &[(&[u8], u32)], scoring: &AffineScoring) -> PreparedSegment {
        let mut set = GermlineSet::new(SegmentKind::V);
        for (k, (sequence, _)) in alleles.iter().enumerate() {
            set.push(&format!("IGHV1-1*{:02}", k + 1), sequence, None)
                .unwrap();
        }
        PreparedSegment::new(set, scoring, scoring.gaps)
    }

    fn frame_table(segment: &PreparedSegment, alleles: &[(&[u8], u32)]) -> FamilyBoundTable {
        FamilyBoundTable::from_placement(segment, alleles.iter().map(|&(_, o)| (0, o)).collect())
    }

    /// Δ by its definition, position by position.
    fn direct_delta(
        r: &SubstitutionProfile,
        offset_r: u32,
        a: &SubstitutionProfile,
        offset_a: u32,
    ) -> Score {
        (0..a.len())
            .map(|j| {
                let f = j as i64 + i64::from(offset_a);
                let at_r = f - i64::from(offset_r);
                let covered = (0..r.len() as i64).contains(&at_r);
                (0..ALPHABET_SIZE as u8)
                    .map(|q| a.at(q, j) - if covered { r.at(q, at_r as usize) } else { 0 })
                    .max()
                    .unwrap_or(0)
                    .max(0)
            })
            .sum()
    }

    fn candidate(allele: AlleleId, band: Band) -> Candidate {
        Candidate {
            allele,
            band,
            votes: 0,
        }
    }

    #[test]
    fn deltas_count_gains_and_overhangs() {
        let alleles: [(&[u8], u32); 4] = [
            (b"ACGTACGT", 2),
            (b"AGGTACGT", 2),   // one SNP
            (b"TTACGTACGT", 0), // 2-nt 5′ extension
            (b"ACGTAC", 2),     // 3′ truncation
        ];
        let segment = frame_segment(&alleles, &SCORING);
        let table = frame_table(&segment, &alleles);
        assert_eq!(table.delta(0, 0), Some(0));
        assert_eq!(table.delta(0, 1), Some(5)); // 2 − (−3)
        assert_eq!(table.delta(1, 0), Some(5));
        assert_eq!(table.delta(0, 2), Some(4)); // two overhanging matches
        assert_eq!(table.delta(2, 0), Some(0)); // positions only R covers add 0
        assert_eq!(table.delta(0, 3), Some(0));
        assert_eq!(table.delta(3, 0), Some(4));
        // Another frame: no bound.
        let mut placement: Vec<_> = alleles.iter().map(|&(_, o)| (0, o)).collect();
        placement[3].0 = 1;
        let split = FamilyBoundTable::from_placement(&segment, placement);
        assert_eq!(split.delta(0, 3), None);
        assert_eq!(split.delta(3, 3), Some(0));
    }

    /// The one-pass table equals the definition on random pairs with
    /// offsets, overhangs on either side, N and hotspot profiles; at equal
    /// span it is the per-position maximum gain.
    #[test]
    fn table_deltas_match_their_definition() {
        let mut next = crate::testing::rng(13);
        let hotspot = crate::scoring::HotspotScoring::for_shm(SCORING);
        for case in 0..200 {
            let base = crate::testing::random_codes(&mut next, 40);
            let alleles: Vec<(Vec<u8>, u32)> = (0..4)
                .map(|_| {
                    let start = (next() % 8) as usize;
                    let end = base.len() - (next() % 8) as usize;
                    let mut codes = base[start..end].to_vec();
                    for _ in 0..next() % 4 {
                        let at = (next() as usize) % codes.len();
                        codes[at] = (next() % 5) as u8;
                    }
                    (codes, start as u32)
                })
                .collect();
            let mut set = GermlineSet::new(SegmentKind::V);
            for (k, (codes, _)) in alleles.iter().enumerate() {
                set.push(&format!("IGHV1-1*{k:02}"), decode(codes).as_bytes(), None)
                    .unwrap();
            }
            let segment = if case % 2 == 0 {
                PreparedSegment::new(set, &SCORING, SCORING.gaps)
            } else {
                PreparedSegment::new(set, &hotspot, SCORING.gaps)
            };
            let table = FamilyBoundTable::from_placement(
                &segment,
                alleles.iter().map(|&(_, o)| (0, o)).collect(),
            );
            for r in 0..4 {
                for a in 0..4 {
                    let (pr, pa) = (segment.profile(r), segment.profile(a));
                    let expected =
                        direct_delta(pr, alleles[r as usize].1, pa, alleles[a as usize].1);
                    assert_eq!(table.delta(r, a), Some(expected), "case {case} ({r}, {a})");
                }
            }
        }
    }

    /// The same-length bound the table generalises:
    /// `Σ_j max(0, max_q (s_A(q, j) − s_R(q, j)))`.
    fn sibling_delta(r: &SubstitutionProfile, a: &SubstitutionProfile) -> Score {
        assert_eq!(r.len(), a.len());
        (0..a.len())
            .map(|j| {
                (0..ALPHABET_SIZE as u8)
                    .map(|q| a.at(q, j) - r.at(q, j))
                    .max()
                    .unwrap_or(0)
                    .max(0)
            })
            .sum()
    }

    /// Alleles of one length at one offset: Δ is the same-length sibling
    /// bound, whether the frame comes from the gene ([`FamilyBoundTable::build`])
    /// or from a placement at any common offset.
    #[test]
    fn frame_delta_equals_sibling_delta_at_equal_length() {
        let mut next = crate::testing::rng(29);
        let hotspot = crate::scoring::HotspotScoring::for_shm(SCORING);
        for case in 0..100 {
            let base = crate::testing::random_codes(&mut next, 30);
            let alleles: Vec<Vec<u8>> = (0..4)
                .map(|_| {
                    let mut codes = base.clone();
                    for _ in 0..next() % 5 {
                        let at = (next() as usize) % codes.len();
                        codes[at] = (next() % 5) as u8;
                    }
                    codes
                })
                .collect();
            let mut set = GermlineSet::new(SegmentKind::V);
            for (k, codes) in alleles.iter().enumerate() {
                set.push(&format!("IGHV1-1*{k:02}"), decode(codes).as_bytes(), None)
                    .unwrap();
            }
            let segment = if case % 2 == 0 {
                PreparedSegment::new(set, &SCORING, SCORING.gaps)
            } else {
                PreparedSegment::new(set, &hotspot, SCORING.gaps)
            };
            let offset = (next() % 4) as u32;
            let tables = [
                FamilyBoundTable::build(&segment),
                FamilyBoundTable::from_placement(&segment, vec![(0, offset); 4]),
            ];
            for r in 0..4 {
                for a in 0..4 {
                    let expected = sibling_delta(segment.profile(r), segment.profile(a));
                    for table in &tables {
                        assert_eq!(table.delta(r, a), Some(expected), "case {case} ({r}, {a})");
                    }
                }
            }
        }
    }

    #[test]
    fn frame_bands_are_the_same_frame_diagonals() {
        let alleles: [(&[u8], u32); 3] = [(b"ACGTACGT", 2), (b"TTACGTACGT", 0), (b"ACGTAC", 5)];
        let segment = frame_segment(&alleles, &SCORING);
        let table = frame_table(&segment, &alleles);
        let banded = |allele, lo, hi| candidate(allele, Band::Diagonals { lo, hi });
        let rep = banded(0, 10, 20);
        assert!(table.frame_band(&rep, &banded(1, 8, 18)));
        assert!(table.frame_band(&rep, &banded(2, 13, 23)));
        assert!(!table.frame_band(&rep, &banded(1, 10, 20)));
        assert!(!table.frame_band(&rep, &banded(1, 8, 19)));
        assert!(!table.frame_band(&rep, &candidate(1, Band::Full)));
        assert!(table.frame_band(&candidate(0, Band::Full), &candidate(2, Band::Full)));
        assert_eq!(table.upper_bound(&rep, &banded(1, 8, 18), 30), Some(34));
        assert_eq!(table.upper_bound(&rep, &banded(1, 9, 18), 30), None);
    }

    /// Exhaustive oracle on tiny alleles: every query of length ≤ 5 over
    /// ACGTN and of length 6 over ACGT, every small frame band shifted into
    /// each allele's coordinates (and the full matrix), every ordered pair of
    /// six alleles — a SNP and an N, a 5′ extension, a 3′ truncation, an
    /// inner start with a 3′ overhang, a 2-nt allele — under three scorings,
    /// one with `open = extend`: `L(π_R, A) ≤ S*(A) ≤ S*(R) + Δ(R, A)` for
    /// the traced path `π_R`, whose frame path sums to its score.
    #[test]
    fn frame_bounds_hold_exhaustively() {
        let alleles: [(&[u8], u32); 6] = [
            (b"ACGTTGCA", 2),
            (b"ACGATGNA", 2),
            (b"GTACGTTGCA", 0),
            (b"ACGTTG", 2),
            (b"GTTCCAGG", 4),
            (b"CC", 7),
        ];
        let mut queries: Vec<Vec<u8>> = vec![Vec::new()];
        for len in 1..=6u32 {
            let alphabet: u32 = if len <= 5 { 5 } else { 4 };
            for code in 0..alphabet.pow(len) {
                queries.push(
                    (0..len)
                        .map(|i| (code / alphabet.pow(i) % alphabet) as u8)
                        .collect(),
                );
            }
        }
        let mut frame_bands = vec![None];
        for lo in -9..=4 {
            for width in 0..=3 {
                frame_bands.push(Some((lo, lo + width)));
            }
        }
        let mut workspace = Default::default();
        let mut path = FramePath::default();
        let mut checked = 0u64;
        for scoring in [
            SCORING,
            AffineScoring::new(2, -3, -2, -1),
            AffineScoring::new(1, -1, -1, -1),
        ] {
            let segment = frame_segment(&alleles, &scoring);
            let table = frame_table(&segment, &alleles);
            for frame_band in &frame_bands {
                let candidates: Vec<Candidate> = (0..alleles.len() as AlleleId)
                    .map(|allele| {
                        let offset = alleles[allele as usize].1 as i32;
                        let band = frame_band.map_or(Band::Full, |(lo, hi)| Band::Diagonals {
                            lo: lo + offset,
                            hi: hi + offset,
                        });
                        candidate(allele, band)
                    })
                    .collect();
                for query in &queries {
                    let paths: Vec<Option<Alignment>> = candidates
                        .iter()
                        .map(|c| {
                            let profile = segment.profile(c.allele);
                            BandedGotoh.align(query, profile, scoring.gaps, c.band, &mut workspace)
                        })
                        .collect();
                    let scores: Vec<Score> = paths
                        .iter()
                        .map(|path| path.as_ref().map_or(0, |path| path.score))
                        .collect();
                    for r in &candidates {
                        let traced = paths[r.allele as usize].as_ref().map(|alignment| {
                            let profile = segment.profile(r.allele);
                            table.frame_path(
                                r.allele,
                                alignment,
                                query,
                                profile,
                                scoring.gaps,
                                &mut path,
                            );
                            assert_eq!(path.prefix.last(), Some(&alignment.score));
                        });
                        for a in &candidates {
                            let upper = table.upper_bound(r, a, scores[r.allele as usize]).unwrap();
                            let lower = traced.map_or(0, |()| {
                                table.path_lower_bound(&path, r, a, &segment).unwrap()
                            });
                            let exact = scores[a.allele as usize];
                            assert!(
                                lower <= exact && exact <= upper,
                                "{scoring:?} {frame_band:?} {query:?}: S({}) = {exact} ∉ [{lower}, {upper}]",
                                a.allele
                            );
                            checked += 1;
                        }
                    }
                }
            }
        }
        assert!(checked > 30_000_000, "{checked} cases");
    }

    /// `L` by its definition: the path's columns from the first to the last
    /// match on a position of `span(A)`, each scored with `A`'s profile.
    fn direct_lower_bound(
        path: &Alignment,
        query: &[u8],
        offsets: (u32, u32),
        sibling: &SubstitutionProfile,
        gaps: GapCosts,
    ) -> Score {
        let shift = i64::from(offsets.0) - i64::from(offsets.1);
        let at = |r: u32| {
            let j = i64::from(r) + shift;
            (0..sibling.len() as i64).contains(&j).then_some(j as usize)
        };
        let columns: Vec<_> = path.columns().collect();
        let inside = |&(op, _, r): &(Op, u32, u32)| op == Op::Match && at(r).is_some();
        let (Some(start), Some(end)) = (
            columns.iter().position(inside),
            columns.iter().rposition(inside),
        ) else {
            return 0;
        };
        let mut score = 0;
        for k in start..=end {
            let (op, q, r) = columns[k];
            score += match op {
                Op::Match => sibling.at(query[q as usize], at(r).expect("inside the clip")),
                gap if columns[k - 1].0 == gap => gaps.extend,
                _ => gaps.open,
            };
        }
        score.max(0)
    }

    /// Longer alleles with realistic bands: mutated, shifted, indel-carrying
    /// copies of either allele of a pair. Both bounds hold, and the
    /// `O(|D|)` lower bound equals its definition.
    #[test]
    fn bounds_bracket_the_local_score() {
        let rep = encode(b"CAGGTCACCTTGAAGGAGTCTGGTCCTGCCTGGTGAAACCCACACAGACC");
        let mut sib = rep[3..].to_vec();
        sib[20] = 1;
        sib[23] = 1;
        sib[41] = 3;
        sib.extend(encode(b"GT"));
        let alleles: [(&[u8], u32); 2] = [
            (&decode(&rep).into_bytes(), 0),
            (&decode(&sib).into_bytes(), 3),
        ];
        let segment = frame_segment(&alleles, &SCORING);
        let table = frame_table(&segment, &alleles);
        let mut words = crate::testing::rng(3);
        let mut next = || words() as usize;
        let mut workspace = Default::default();
        let mut path = FramePath::default();
        let mut checked = 0;
        for _ in 0..300 {
            let source = if next() % 2 == 0 { &rep } else { &sib };
            let mut query: Vec<u8> = (0..next() % 12).map(|_| (next() % 4) as u8).collect();
            let lead = query.len() as i32;
            query.extend(source.iter().map(|&c| {
                if next() % 8 == 0 {
                    (next() % 5) as u8
                } else {
                    c
                }
            }));
            if next() % 2 == 0 {
                let at = next() % query.len();
                query.insert(at, (next() % 4) as u8);
            }
            let (lo, hi) = match next() % 3 {
                0 => (lead - 2, lead + 2),
                _ => {
                    let centre = lead + (next() % 7) as i32 - 3;
                    (centre - 8, centre + 8)
                }
            };
            let full = next() % 4 == 0;
            let pair: Vec<Candidate> = (0..2)
                .map(|allele| {
                    let offset = alleles[allele as usize].1 as i32;
                    let band = if full {
                        Band::Full
                    } else {
                        Band::Diagonals {
                            lo: lo + offset,
                            hi: hi + offset,
                        }
                    };
                    candidate(allele, band)
                })
                .collect();
            let paths: Vec<Option<Alignment>> = pair
                .iter()
                .map(|c| {
                    let profile = segment.profile(c.allele);
                    BandedGotoh.align(&query, profile, SCORING.gaps, c.band, &mut workspace)
                })
                .collect();
            let score = |c: &Candidate| paths[c.allele as usize].as_ref().map_or(0, |p| p.score);
            for (r, a) in [(&pair[0], &pair[1]), (&pair[1], &pair[0])] {
                let upper = table.upper_bound(r, a, score(r)).unwrap();
                assert!(score(a) <= upper, "{:?}: {} > {upper}", a.band, score(a));
                let Some(alignment) = &paths[r.allele as usize] else {
                    continue;
                };
                let profile = segment.profile(r.allele);
                table.frame_path(
                    r.allele,
                    alignment,
                    &query,
                    profile,
                    SCORING.gaps,
                    &mut path,
                );
                let lower = table.path_lower_bound(&path, r, a, &segment).unwrap();
                let offsets = (alleles[r.allele as usize].1, alleles[a.allele as usize].1);
                let profile = segment.profile(a.allele);
                let direct = direct_lower_bound(alignment, &query, offsets, profile, SCORING.gaps);
                assert_eq!(lower, direct, "{:?}", a.band);
                assert!(lower <= score(a), "{:?}: {} < {lower}", a.band, score(a));
                checked += u32::from(lower > 0);
            }
        }
        assert!(checked > 300, "{checked} positive lower bounds");
    }
}
