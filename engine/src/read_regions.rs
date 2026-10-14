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
//! called V aligns Cys104's first base and the called J Trp/Phe118's, and
//! at most `max_junction_span + 6` bases separate them (`complete_vdj` too
//! then needs the junction). Coordinates are 0-based half-open on the
//! oriented query.
//!
//! The reading frame is the junction's (Cys104's codon); without a junction
//! it is the called V's, from its anchor or given coding frame
//! (`frame=`/`SWIGMETA=`). A segment has a frameshift when its called
//! alignment's net indel length is not a multiple of 3, and a productive
//! rearrangement is in frame, free of stop codons and of frameshifts.

use std::ops::Range;

use crate::adjudicate::CalledSegment;
use crate::align::Op;
use crate::germline::{Allele, SegmentKind};
use crate::pipeline::{AlleleCatalog, Call};
use crate::placement::Placement;
use crate::translate::{STOP, codon};

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

impl ReadRegions {
    pub fn from_call<C: AlleleCatalog>(call: &Call, catalog: &C) -> Self {
        let n = call.query.len() as i64;
        let v = called(call, catalog, SegmentKind::V);
        let j = called(call, catalog, SegmentKind::J);
        let mut out = Self::default();

        let v_end = v.as_ref().map(|v| i64::from(v.called.span.q_hi));
        let j_start = j.as_ref().map(|j| i64::from(j.called.span.q_lo));
        let cys104 = v
            .as_ref()
            .and_then(|v| v.at(v.allele.anchor))
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
        // their called segments, and it spans at most the junction span and
        // the codons between them, as SwiftIG's did: an anchor projected
        // past its alignment's end (a V cut before Cys104, a J aligned only
        // after Trp/Phe118) bounds regions but reports no junction.
        let junction_holds = |c: i64, w: i64| match catalog.pair_guards() {
            None => true,
            Some(guards) => {
                let aligned = |s: &Option<Projected<'_>>, at: i64| {
                    s.as_ref().and_then(|s| s.anchor(s.allele.anchor, true)) == Some(at)
                };
                aligned(&v, c) && aligned(&j, w) && w - c <= i64::from(guards.max_junction_span) + 6
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
        let database = prepare(
            db.v.clone(),
            Some(db.d.clone()),
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
