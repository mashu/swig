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
//! along the alignment's end diagonal. FR1–FR3 never start before the read's
//! V alignment, FR4 never ends after the J alignment, and the junction
//! (Cys104 … Trp/Phe118, inclusive) is reported only when it lies inside the
//! read. Coordinates are 0-based half-open on the oriented query.

use std::ops::Range;

use crate::adjudicate::CalledSegment;
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
    /// AIRR names of the FR/CDR regions, in output column order.
    pub const NAMES: [&'static str; 7] = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4", "cdr3"];

    /// The FR/CDR intervals in [`ReadRegions::NAMES`] order.
    pub fn intervals(&self) -> [Option<&Range<u32>>; 7] {
        [
            &self.fwr1, &self.cdr1, &self.fwr2, &self.cdr2, &self.fwr3, &self.fwr4, &self.cdr3,
        ]
        .map(Option::as_ref)
    }

    pub fn from_call<C: AlleleCatalog>(call: &Call, catalog: &C) -> Self {
        let n = call.query.len() as i64;
        let v = called(call, catalog, SegmentKind::V);
        let j = called(call, catalog, SegmentKind::J);
        let mut out = Self::default();

        let cys104 = v.as_ref().and_then(|v| v.at(v.allele.anchor));
        if let Some(v) = &v {
            let r = &v.allele.regions;
            let floor = i64::from(v.called.span.q_lo);
            let [fr1, cdr1, fr2, cdr2, fr3] =
                [r.fr1, r.cdr1, r.fr2, r.cdr2, r.fr3].map(|x| v.at(x));
            out.fwr1 = interval(fr1, cdr1, floor, n);
            out.cdr1 = interval(cdr1, fr2, floor, n);
            out.fwr2 = interval(fr2, cdr2, floor, n);
            out.cdr2 = interval(cdr2, fr3, floor, n);
            out.fwr3 = interval(fr3, cys104.map(|c| c + 3), floor, n);
        }
        let trp118 = j.as_ref().and_then(|j| j.at(j.allele.anchor));
        if let Some(j) = &j {
            out.fwr4 = interval(trp118, Some(i64::from(j.called.span.q_hi)), 0, n);
        }
        if let (Some(c), Some(w)) = (cys104, trp118)
            && c >= 0
            && w + 3 <= n
            && c < w
        {
            out.junction = Some(c as u32..(w + 3) as u32);
            out.cdr3 = interval(Some(c + 3), Some(w), 0, n);
            out.vj_in_frame = Some((w + 3 - c) % 3 == 0);
        }
        // Stop codons in the V reading frame from the V alignment start to the J end.
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
        out.productive = out
            .vj_in_frame
            .zip(out.stop_codon)
            .map(|(in_frame, stop)| in_frame && !stop);
        out
    }
}
