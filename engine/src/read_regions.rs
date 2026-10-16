//! Per-read FR/CDR/CDR3/junction coordinates and productivity.
//!
//! Germline region starts (from headers or [`crate::annotate`]) are projected
//! onto the read through the *called* V and J alignments. Nothing here runs
//! during assignment: it is evaluated lazily when a call is reported, only for
//! the two called segments, and each boundary is one binary search in the
//! alignment's clip index ([`crate::clip::ClipIndex::project_reference`]) —
//! O(log L) per boundary, about a microsecond per read in total.
//!
//! Projection rules: a germline position inside the aligned part maps to its
//! aligned read base (the next read base across a deletion); outside it,
//! along the alignment's end diagonal, however far. A projected anchor
//! counts only on its own side of the other segment's alignment: Cys104
//! starts before the J's, Trp/Phe118 ends after the V's; and neither counts
//! when their codons overlap, since one of them was then extrapolated
//! through bases the other segment explains. FR1–FR3 never start before the
//! read's V alignment, FR4 never ends after the J alignment, FR3 ends where
//! Cys104's codon does and FR4 starts at Trp/Phe118, so the two never
//! overlap. The junction (Cys104 … Trp/Phe118, inclusive) is reported only
//! when both anchors count and it lies inside the read; under the pair
//! guards ([`crate::pipeline::PairGuards`]), as in SwiftIG, only when the
//! called V aligns Cys104's first base (or its clipped Cys104 is
//! recovered, below) and the called J Trp/Phe118's, and at most
//! `max_junction_span + 6` bases separate them (`complete_vdj` too then
//! needs the junction). Coordinates are 0-based half-open on the oriented
//! query.
//!
//! **Clipped Cys104** (as in SwiftIG 0.38.12). A local V alignment can stop
//! short of Cys104 because clustered substitutions before it lower its
//! score. Under the pair guards Cys104 then still counts, at its projection
//! along the called V's end diagonal, when all of these hold
//! ([`clipped_cys104`]): the allele's anchor codon is a cysteine and ends
//! at most 18 germline bases past the called V; the read has an intact
//! cysteine codon (`TGT`/`TGC`) at the projection, ending at or before the
//! first base of the called D or J; every read and germline base of the
//! co-linear stretch, from up to 12 gapless aligned columns before the V's
//! end through the codon, is unambiguous; and no alignment of that germline
//! stretch holding a gap, anchored at its start and free to end anywhere
//! up to 18 read bases past the codon, scores as much as the co-linear one.
//! No other cysteine is searched for, and an in-frame junction plays no
//! part in the decision. Calls, scores, alignments and segment boundaries
//! do not change, only the junction and what follows from it (CDR3,
//! `vj_in_frame`, `productive`, the reading frame, `d_frame`,
//! `complete_vdj`); the stop-codon and frameshift vetoes still apply. The
//! codon's coordinate is structural: junctional additions can rebuild a
//! deleted codon, so it does not prove that V contributed every base.
//!
//! The reading frame is the junction's (Cys104's codon); without a junction
//! it is the called V's, from its anchor or given coding frame
//! (`frame=`/`SWIGMETA=`). A segment has a frameshift when its called
//! alignment's net indel length is not a multiple of 3, and a productive
//! rearrangement is in frame, free of stop codons and of frameshifts.

use std::ops::Range;

use crate::adjudicate::CalledSegment;
use crate::align::Op;
use crate::alphabet::N;
use crate::germline::{Allele, SegmentKind};
use crate::pipeline::{AlleleCatalog, Call};
use crate::placement::Placement;
use crate::scoring::{GapCosts, Score};
use crate::segment::PreparedSegment;
use crate::translate::{STOP, codon};

/// Germline bases at most between the called V's end and the end of a
/// recovered clipped Cys104 codon ([`clipped_cys104`]), and read bases at
/// most past the codon where a rival gapped alignment may end.
const CLIPPED_TAIL: u32 = 18;
/// Gapless aligned columns at most before the called V's end that the
/// recovery compares again, so an indel across the end cannot pass for a
/// co-linear tail.
const CLIPPED_CONTEXT: u32 = 12;
/// Germline bases the recovery compares at most.
const CLIPPED_WINDOW: usize = (CLIPPED_CONTEXT + CLIPPED_TAIL) as usize;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadRegions {
    pub fwr1: Option<Range<u32>>,
    pub cdr1: Option<Range<u32>>,
    pub fwr2: Option<Range<u32>>,
    pub cdr2: Option<Range<u32>>,
    pub fwr3: Option<Range<u32>>,
    pub cdr3: Option<Range<u32>>,
    pub fwr4: Option<Range<u32>>,
    pub junction: Option<Range<u32>>,
    pub vj_in_frame: Option<bool>,
    pub stop_codon: Option<bool>,
    pub productive: Option<bool>,
    /// The reading frame: offset (0–2) on the oriented query of the first
    /// base of a codon.
    pub frame: Option<u32>,
    pub v_frameshift: Option<bool>,
    pub j_frameshift: Option<bool>,
    /// The called D's reading frame relative to its own coding frame, 0–2
    /// (AIRR `d_frame` − 1), when both are known.
    pub d_frame: Option<u32>,
    /// The V starts within 3 bases of its germline start and the J ends
    /// within 3 of its germline end.
    pub complete_vdj: Option<bool>,
}

struct Projected<'a> {
    called: &'a CalledSegment,
    placement: &'a Placement,
    allele: &'a Allele,
}

impl Projected<'_> {
    fn at(&self, germline: Option<u32>) -> Option<i64> {
        germline.map(|r| self.placement.clip.project_reference(r))
    }

    /// The read base of conserved anchor `germline`: with `aligned`, only
    /// one the called alignment aligns to it (SwiftIG's rule), else
    /// projected however far ([`Self::at`]).
    fn anchor(&self, germline: Option<u32>, aligned: bool) -> Option<i64> {
        if aligned {
            let r = germline?;
            let q = self.placement.clip.aligned_query_in(&self.called.span, r)?;
            Some(i64::from(q))
        } else {
            self.at(germline)
        }
    }

    /// Whether the called alignment's net indel length is not a multiple of 3.
    fn frameshift(&self) -> bool {
        let net: i64 = self
            .placement
            .clip
            .columns(&self.called.span)
            .map(|(op, _, _)| match op {
                Op::Match => 0,
                Op::Ins => 1,
                Op::Del => -1,
            })
            .sum();
        net.rem_euclid(3) != 0
    }

    /// The query offset (0–2) of codon starts by this germline's coding
    /// frame, `germline_frame`, at the called alignment's first column.
    fn query_frame(&self, germline_frame: u32) -> u32 {
        let span = &self.called.span;
        let reference_phase = (span.r_lo + 3 - germline_frame) % 3;
        (span.q_lo + 3 - reference_phase) % 3
    }
}

fn called<'a, C: AlleleCatalog>(
    call: &'a Call,
    catalog: &'a C,
    kind: SegmentKind,
) -> Option<Projected<'a>> {
    let path = call.path.as_ref()?;
    let mut segments = path
        .segments
        .iter()
        .filter(|s| call.placements[s.placement].kind == kind);
    let called = if kind == SegmentKind::J {
        segments.next_back()?
    } else {
        segments.next()?
    };
    let placement = &call.placements[called.placement];
    let allele = catalog
        .segment(kind)?
        .set()
        .allele(placement.primary_allele());
    Some(Projected {
        called,
        placement,
        allele,
    })
}

/// `[start, end)` clamped to `[floor, ceil)`, if non-empty.
fn interval(start: Option<i64>, end: Option<i64>, floor: i64, ceil: i64) -> Option<Range<u32>> {
    let (start, end) = (start?.max(floor), end?.min(ceil));
    (start < end).then_some(start as u32..end as u32)
}

/// Whether `codes` hold a cysteine codon (`TGT`/`TGC`) at `at`.
fn is_cys(codes: &[u8], at: u32) -> bool {
    let at = at as usize;
    codes
        .get(at..at + 3)
        .is_some_and(|c| codon(c[0], c[1], c[2]) == b'C')
}

/// The read position of Cys104 of a called V that ends before it, when
/// the clipped-Cys104 checks of the [module documentation](self) hold:
/// the projection along the called V's end diagonal. `boundary` is the
/// first read base of the called D or J; `segment` is the V database.
fn clipped_cys104(
    query: &[u8],
    v: &Projected<'_>,
    segment: &PreparedSegment,
    boundary: u32,
) -> Option<u32> {
    let (span, germline) = (&v.called.span, v.allele.codes.as_slice());
    let anchor = v.allele.anchor?;
    if anchor < span.r_hi || anchor + 3 - span.r_hi > CLIPPED_TAIL || !is_cys(germline, anchor) {
        return None;
    }
    let tail = anchor + 3 - span.r_hi;
    let boundary = boundary.min(query.len() as u32);
    if span.q_hi + tail > boundary {
        return None;
    }
    let cys = span.q_hi + anchor - span.r_hi;
    if !is_cys(query, cys) {
        return None;
    }
    let context = v
        .placement
        .clip
        .columns(span)
        .rev()
        .take(CLIPPED_CONTEXT as usize)
        .take_while(|&(op, ..)| op == Op::Match)
        .count() as u32;
    let (q0, r0) = (
        (span.q_hi - context) as usize,
        (span.r_hi - context) as usize,
    );
    let len = (context + tail) as usize;
    let read = &query[q0..(boundary as usize).min(q0 + len + CLIPPED_TAIL as usize)];
    if read[..len]
        .iter()
        .chain(&germline[r0..r0 + len])
        .any(|&c| c >= N)
    {
        return None;
    }
    let profile = segment.profile(v.placement.primary_allele());
    let score = |i: usize, j: usize| profile.at(read[i], r0 + j);
    (!gapped_rival(read.len(), len, score, segment.gaps())).then_some(cys)
}

/// Whether an alignment of germline bases `0..len` against a read of
/// `read_len` bases that holds a gap scores at least as much as the
/// co-linear one (read base `i` against germline base `i`, all `len`).
/// Both start at base 0; the read may end anywhere. `score(i, j)` scores
/// read base `i` against germline base `j`.
///
/// Affine DP over (read, germline) cells keeping only paths that hold a
/// gap: such a path leaves the co-linear diagonal once, by a gap, so the
/// diagonal's prefix scores enter each cell beside it as a gap's source.
/// There is no local reset.
fn gapped_rival(
    read_len: usize,
    len: usize,
    score: impl Fn(usize, usize) -> Score,
    gaps: GapCosts,
) -> bool {
    const NONE: Score = Score::MIN / 4;
    let mut diagonal = [0; CLIPPED_WINDOW + 1];
    for i in 1..=len {
        diagonal[i] = diagonal[i - 1] + score(i - 1, i - 1);
    }
    let co_linear = diagonal[len];
    // Row 0: the germline deleted from the start.
    let mut previous = [NONE; CLIPPED_WINDOW + 1];
    for (j, cell) in previous.iter_mut().enumerate().take(len + 1).skip(1) {
        *cell = gaps.cost(j as u32);
    }
    if previous[len] >= co_linear {
        return true;
    }
    let mut current = [NONE; CLIPPED_WINDOW + 1];
    let mut insertion = [NONE; CLIPPED_WINDOW + 1];
    for i in 1..=read_len {
        current[0] = gaps.cost(i as u32);
        let mut deletion = NONE;
        for j in 1..=len {
            let above = previous[j].max(if i - 1 == j { diagonal[j] } else { NONE });
            insertion[j] = (above + gaps.open).max(insertion[j] + gaps.extend);
            let left = current[j - 1].max(if i + 1 == j { diagonal[i] } else { NONE });
            deletion = (left + gaps.open).max(deletion + gaps.extend);
            current[j] = (previous[j - 1] + score(i - 1, j - 1))
                .max(insertion[j])
                .max(deletion);
        }
        if current[len] >= co_linear {
            return true;
        }
        std::mem::swap(&mut previous, &mut current);
    }
    false
}

impl ReadRegions {
    pub fn from_call<C: AlleleCatalog>(call: &Call, catalog: &C) -> Self {
        let n = call.query.len() as i64;
        let v = called(call, catalog, SegmentKind::V);
        let j = called(call, catalog, SegmentKind::J);
        let mut out = Self::default();

        let v_end = v.as_ref().map(|v| i64::from(v.called.span.q_hi));
        let j_start = j.as_ref().map(|j| i64::from(j.called.span.q_lo));
        // Under the pair guards, a Cys104 the called V does not align may
        // still be recovered just past its end, before the called D or J.
        let clipped = match (
            catalog.pair_guards(),
            &v,
            &j,
            catalog.segment(SegmentKind::V),
        ) {
            (Some(_), Some(v), Some(j), Some(segment))
                if v.anchor(v.allele.anchor, true).is_none() =>
            {
                let d = called(call, catalog, SegmentKind::D);
                let boundary = d.map_or(j.called.span.q_lo, |d| {
                    d.called.span.q_lo.min(j.called.span.q_lo)
                });
                clipped_cys104(&call.query, v, segment, boundary).map(i64::from)
            }
            _ => None,
        };
        let cys104 = v
            .as_ref()
            .and_then(|v| clipped.or_else(|| v.at(v.allele.anchor)))
            .filter(|&c| j_start.is_none_or(|s| c < s));
        let trp118 = j
            .as_ref()
            .and_then(|j| j.at(j.allele.anchor))
            .filter(|&w| v_end.is_none_or(|e| w + 3 > e));
        let (cys104, trp118) = match (cys104, trp118) {
            (Some(c), Some(w)) if c + 3 > w => (None, None),
            anchors => anchors,
        };
        // Under the pair guards the junction needs both anchors aligned in
        // their called segments (or Cys104 recovered past a clipped V), and
        // it spans at most the junction span and the codons between them,
        // as SwiftIG's did: an anchor projected past its alignment's end
        // (a V cut well before Cys104, a J aligned only after Trp/Phe118)
        // bounds regions but reports no junction.
        let junction_holds = |c: i64, w: i64| match catalog.pair_guards() {
            None => true,
            Some(guards) => {
                let aligned = |s: &Option<Projected<'_>>, at: i64| {
                    s.as_ref().and_then(|s| s.anchor(s.allele.anchor, true)) == Some(at)
                };
                (aligned(&v, c) || clipped == Some(c))
                    && aligned(&j, w)
                    && w - c <= i64::from(guards.max_junction_span) + 6
            }
        };
        if let Some(v) = &v {
            let r = &v.allele.regions;
            let floor = i64::from(v.called.span.q_lo);
            let [fr1, cdr1, fr2, cdr2, fr3] =
                [r.fr1, r.cdr1, r.fr2, r.cdr2, r.fr3].map(|x| v.at(x));
            out.fwr1 = interval(fr1, cdr1, floor, n);
            out.cdr1 = interval(cdr1, fr2, floor, n);
            out.fwr2 = interval(fr2, cdr2, floor, n);
            out.cdr2 = interval(cdr2, fr3, floor, n);
            // Without a placed Cys104, FR3 ends with the V alignment, as in
            // SwiftIG: a V trimmed past its anchor still has its FR3 bases.
            // It stops where FR4 starts, at Trp/Phe118, which may lie up to
            // two bases inside the V's end.
            let fr3_end = cys104.map_or_else(
                || {
                    let v_end = i64::from(v.called.span.q_hi);
                    trp118.map_or(v_end, |w| v_end.min(w))
                },
                |c| c + 3,
            );
            out.fwr3 = interval(fr3, Some(fr3_end), floor, n);
        }
        if let Some(j) = &j {
            out.fwr4 = interval(trp118, Some(i64::from(j.called.span.q_hi)), 0, n);
        }
        if let (Some(c), Some(w)) = (cys104, trp118)
            && c >= 0
            && w + 3 <= n
            && junction_holds(c, w)
        {
            out.junction = Some(c as u32..(w + 3) as u32);
            out.cdr3 = interval(Some(c + 3), Some(w), 0, n);
            out.vj_in_frame = Some((w + 3 - c) % 3 == 0);
        }
        // Stop codons in the read, as a ribosome reads it (no indel
        // correction), in the frame of the Cys104 codon (sequence_frame
        // with a junction) from the V alignment start to the J end.
        if let (Some(v), Some(c)) = (&v, cys104) {
            let end = j
                .as_ref()
                .map_or(i64::from(v.called.span.q_hi), |j| {
                    i64::from(j.called.span.q_hi)
                })
                .min(n);
            let start = i64::from(v.called.span.q_lo);
            let first = start + (c - start).rem_euclid(3);
            let q = &call.query;
            let stop = (first..end - 2)
                .step_by(3)
                .any(|p| codon(q[p as usize], q[p as usize + 1], q[p as usize + 2]) == STOP);
            out.stop_codon = Some(stop);
        }
        out.v_frameshift = v.as_ref().map(Projected::frameshift);
        out.j_frameshift = j.as_ref().map(Projected::frameshift);
        out.productive = match (out.vj_in_frame, out.stop_codon) {
            (Some(in_frame), Some(stop)) => Some(
                in_frame
                    && !stop
                    && out.v_frameshift != Some(true)
                    && out.j_frameshift != Some(true),
            ),
            _ => None,
        };
        out.frame = match (cys104, &v) {
            (Some(c), _) if c >= 0 && out.junction.is_some() => Some((c % 3) as u32),
            (_, Some(v)) => {
                let germline = v.allele.anchor.map(|a| a % 3);
                Some(v.query_frame(germline.or(v.allele.frame.map(u32::from)).unwrap_or(0)))
            }
            _ => None,
        };
        if let (Some(frame), Some(d)) = (out.frame, called(call, catalog, SegmentKind::D))
            && let Some(own) = d.allele.frame
        {
            out.d_frame = Some((frame + 3 - d.query_frame(u32::from(own))) % 3);
        }
        // Under the pair guards, as in SwiftIG, only with a junction.
        if let (Some(v), Some(j)) = (&v, &j)
            && (catalog.pair_guards().is_none() || out.junction.is_some())
        {
            let j_len = j.allele.codes.len() as u32;
            out.complete_vdj = Some(v.called.span.r_lo <= 2 && j.called.span.r_hi + 3 >= j_len);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;
    use crate::annotate::AnnotateOptions;
    use crate::config::DefaultAssigner;
    use crate::pipeline::VdjAssigner;
    use crate::prepare::prepare;
    use crate::profiles::CallingProfile;
    use crate::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

    fn database() -> SyntheticDatabase {
        SyntheticDatabase::generate(&SyntheticSpec::default(), 5)
    }

    fn assigner(db: &SyntheticDatabase) -> DefaultAssigner {
        assigner_with(db, true)
    }

    /// The R-optimized assigner of `db`, with its D database or without.
    fn assigner_with(db: &SyntheticDatabase, with_d: bool) -> DefaultAssigner {
        let database = prepare(
            db.v.clone(),
            with_d.then(|| db.d.clone()),
            db.j.clone(),
            AnnotateOptions::default(),
        )
        .unwrap();
        let config = CallingProfile::ROptimized.config();
        DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap()
    }

    /// The called V's and J's read spans, `[q_lo, q_hi)`.
    fn spans(call: &Call, assigner: &DefaultAssigner) -> [Option<Range<u32>>; 2] {
        [SegmentKind::V, SegmentKind::J].map(|kind| {
            called(call, assigner, kind).map(|s| s.called.span.q_lo..s.called.span.q_hi)
        })
    }

    /// The rules every read's regions keep: FR3 ends where FR4 may start,
    /// and a junction holds two whole conserved codons, Cys104 starting
    /// before the J's alignment and Trp/Phe118 ending after the V's.
    fn check(regions: &ReadRegions, [v, j]: &[Option<Range<u32>>; 2], read: &str) {
        if let (Some(fr3), Some(fr4)) = (&regions.fwr3, &regions.fwr4) {
            assert!(
                fr3.end <= fr4.start,
                "FR3 {fr3:?} overlaps FR4 {fr4:?}: {read}"
            );
        }
        if let Some(junction) = &regions.junction {
            assert!(junction.len() >= 6, "junction {junction:?}: {read}");
            let (v, j) = (v.as_ref().unwrap(), j.as_ref().unwrap());
            assert!(junction.start < j.start && junction.end > v.end, "{read}");
        }
    }

    /// A V cut before its Cys104 codon followed by a J from 4–5 bases
    /// before its Trp118 codon: Cys104 projects into the J alignment, onto
    /// Trp118's codon, so neither the junction nor FR3 reaches past it.
    #[test]
    fn crossed_anchor_codons_give_no_junction() {
        let db = database();
        let assigner = assigner(&db);
        let (v, j) = (&db.v.alleles()[0], &db.j.alleles()[0]);
        let (cys, trp) = (v.anchor.unwrap() as usize, j.anchor.unwrap() as usize);
        for before in [4, 5] {
            let read =
                crate::alphabet::decode(&[&v.codes[..cys - 3], &j.codes[trp - before..]].concat());
            let call = assigner.assign(&encode(read.as_bytes()), &mut assigner.workspace());
            let regions = ReadRegions::from_call(&call, &assigner);
            check(&regions, &spans(&call, &assigner), &read);
            assert_eq!(regions.junction, None, "{read}");
            assert_eq!(regions.productive, None, "{read}");
        }
    }

    /// A light-chain-like read whose V is trimmed past Cys104 and whose J
    /// is trimmed past Phe118, abutting: each anchor would project into
    /// the other segment's alignment, so neither counts.
    #[test]
    fn anchors_projected_into_the_other_segment_do_not_count() {
        let db = database();
        let assigner = assigner(&db);
        let (v, j) = (&db.v.alleles()[0], &db.j.alleles()[0]);
        let (cys, trp) = (v.anchor.unwrap() as usize, j.anchor.unwrap() as usize);
        let read = crate::alphabet::decode(&[&v.codes[..cys - 6], &j.codes[trp + 6..]].concat());
        let call = assigner.assign(&encode(read.as_bytes()), &mut assigner.workspace());
        let regions = ReadRegions::from_call(&call, &assigner);
        let spans = spans(&call, &assigner);
        check(&regions, &spans, &read);
        assert_eq!(regions.junction, None);
        assert_eq!(regions.productive, None);
        let (v_span, j_span) = (spans[0].clone().unwrap(), spans[1].clone().unwrap());
        assert!(
            regions
                .fwr3
                .as_ref()
                .is_none_or(|fr3| fr3.end <= j_span.start)
        );
        assert!(
            regions
                .fwr4
                .as_ref()
                .is_none_or(|fr4| fr4.start + 3 > v_span.end)
        );
    }

    /// A read of V allele 0 whose 4 bases before Cys104 are substituted,
    /// so the V alignment stops 7 germline bases before the codon's end,
    /// then the codon (`codon`, else the germline's), 10 junction bases and
    /// J allele 0. Returns the read and Cys104's position on it.
    fn clipped_v_read(db: &SyntheticDatabase, codon: Option<[u8; 3]>) -> (Vec<u8>, usize) {
        let (v, j) = (&db.v.alleles()[0], &db.j.alleles()[0]);
        let anchor = v.anchor.unwrap() as usize;
        let substituted: Vec<u8> = v.codes[anchor - 4..anchor].iter().map(|&c| 3 - c).collect();
        let codon = codon.map_or_else(|| v.codes[anchor..anchor + 3].to_vec(), Vec::from);
        let read = [
            &v.codes[..anchor - 4],
            &substituted,
            &codon,
            &encode(b"CCCCAAAGGG"),
            &j.codes,
        ]
        .concat();
        (read, anchor)
    }

    /// The called V's span and the read's regions.
    fn regions_of(read: &[u8]) -> (Range<u32>, ReadRegions) {
        regions_with(read, true)
    }

    /// [`regions_of`] with the D database or without.
    fn regions_with(read: &[u8], with_d: bool) -> (Range<u32>, ReadRegions) {
        let db = database();
        let assigner = assigner_with(&db, with_d);
        let call = assigner.assign(read, &mut assigner.workspace());
        let regions = ReadRegions::from_call(&call, &assigner);
        let spans = spans(&call, &assigner);
        check(&regions, &spans, &crate::alphabet::decode(read));
        (spans[0].clone().unwrap(), regions)
    }

    /// Clustered substitutions before Cys104 stop the V's alignment short
    /// of it; the intact codon there is still the junction's start.
    #[test]
    fn a_clipped_v_keeps_its_cys104() {
        let db = database();
        let (read, cys) = clipped_v_read(&db, None);
        assert!(is_cys(&read, cys as u32));
        let (v, regions) = regions_of(&read);
        assert!(
            v.end as usize <= cys - 4,
            "V {v:?} reaches into the substitutions"
        );
        let junction = regions.junction.expect("Cys104 recovered");
        assert_eq!(junction.start as usize, cys);
        assert_eq!(regions.cdr3.unwrap().start as usize, cys + 3);
        assert_eq!(regions.frame, Some((cys % 3) as u32));
        assert!(regions.productive.is_some());
    }

    /// No other codon stands in for a substituted cysteine, and an unknown
    /// base in the clipped tail leaves Cys104 unresolved.
    #[test]
    fn a_clipped_v_needs_an_intact_unambiguous_tail() {
        let db = database();
        let (read, cys) = clipped_v_read(&db, Some([3, 2, 2]));
        let (_, regions) = regions_of(&read);
        assert_eq!(regions.junction, None, "TGG is no cysteine");
        assert_eq!(regions.productive, None);
        let (mut read, cys_at) = clipped_v_read(&db, None);
        assert_eq!(cys, cys_at);
        read[cys - 2] = N;
        let (_, regions) = regions_of(&read);
        assert_eq!(regions.junction, None, "N in the tail");
    }

    /// A V clipped more than 18 germline bases before Cys104's codon ends
    /// keeps no junction. Without a D database no chance D in the tail
    /// stops the recovery first.
    #[test]
    fn a_long_clipped_tail_is_not_recovered() {
        let db = database();
        let (v, j) = (&db.v.alleles()[0], &db.j.alleles()[0]);
        let anchor = v.anchor.unwrap() as usize;
        // Every other base of the 18 before the codon substituted.
        let mut tail = v.codes[anchor - 18..anchor].to_vec();
        for c in tail.iter_mut().step_by(2) {
            *c = 3 - *c;
        }
        let read = [
            &v.codes[..anchor - 18],
            &tail,
            &v.codes[anchor..anchor + 3],
            &encode(b"CCCCAAAGGG"),
            &j.codes,
        ]
        .concat();
        let (span, regions) = regions_with(&read, false);
        assert!(
            anchor + 3 - span.end as usize > CLIPPED_TAIL as usize,
            "V {span:?} ends within the tail limit"
        );
        assert_eq!(regions.junction, None);
    }

    /// The rival search: a gapped alignment scoring as much as the
    /// co-linear one blocks the recovery, one scoring less does not.
    #[test]
    fn gapped_rivals_block_the_recovery() {
        use crate::scoring::AffineScoring;
        let rival = |read: &[u8], germline: &[u8], scoring: AffineScoring| {
            let (read, germline) = (encode(read), encode(germline));
            gapped_rival(
                read.len(),
                germline.len(),
                |i, j| scoring.substitution(read[i], germline[j]),
                scoring.gaps,
            )
        };
        let scoring = AffineScoring::new(2, -3, -5, -1);
        // Co-linear: two mismatches (−6) among matches; nothing gapped
        // comes close.
        assert!(!rival(b"ACGTTCATGCAGTGT", b"ACGTACATGGAGTGT", scoring));
        // One read base inserted: co-linear is mostly mismatches.
        assert!(rival(b"ACGTAACATGGAGTGT", b"ACGTACATGGAGTGT", scoring));
        // A mismatched germline end costs −9 co-linear, −7 deleted.
        assert!(rival(b"ACGTACATGGAGAAA", b"ACGTACATGGAGTTT", scoring));
        // Exact tie: the mismatch (−3) equals an insertion and a match.
        assert!(rival(b"CA", b"A", scoring));
        assert!(!rival(b"CA", b"A", AffineScoring::new(2, -3, -6, -1)));
    }

    /// Every simulated read of the hard preset keeps the region rules.
    #[test]
    fn simulated_reads_keep_the_region_rules() {
        let db = database();
        let assigner = assigner(&db);
        let params = SimParams::preset("hard").unwrap();
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params).unwrap();
        let mut workspace = assigner.workspace();
        for read in sim.reads(600, 11).unwrap() {
            let call = assigner.assign(&encode(&read.sequence), &mut workspace);
            let regions = ReadRegions::from_call(&call, &assigner);
            check(
                &regions,
                &spans(&call, &assigner),
                &String::from_utf8_lossy(&read.sequence),
            );
        }
    }
}
