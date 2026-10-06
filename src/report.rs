//! AIRR Rearrangement output.
//!
//! Coordinates follow the AIRR convention: 1-based closed intervals on the
//! *oriented* query (the reverse complement of `sequence` when `rev_comp`).
//! CIGAR strings use `S` for query bases before/after the segment, `N` for
//! germline bases before the alignment, and `=`/`X`/`I`/`D` for columns.
//!
//! Co-optimal labels are recomputed on the called cut interval: every
//! placement of the same kind whose clip to `[entry, exit)` scores exactly the
//! called score contributes its alleles. As the note says, this means "equal
//! score under this caller", not equal biological probability.

use std::fmt::Write as _;

use crate::adjudicate::CalledSegment;
use crate::align::Op;
use crate::alphabet::{decode, decode_base, N};
use crate::evidence::Evidence;
use crate::germline::SegmentKind;
use crate::pipeline::{AlleleCatalog, Call};
use crate::read_regions::ReadRegions;
use crate::scoring::Score;
use crate::stats::KarlinAltschul;
use crate::translate::translate;

#[derive(Clone, Debug, PartialEq)]
pub struct SegmentReport {
    /// Comma-separated co-optimal allele names, called allele first.
    pub call: String,
    pub score: Score,
    pub identity: f64,
    /// Karlin–Altschul E-value, when a model is supplied for the segment kind.
    pub support: Option<f64>,
    pub cigar: String,
    pub sequence_start: u32,
    pub sequence_end: u32,
    pub germline_start: u32,
    pub germline_end: u32,
    pub sequence_alignment: String,
    pub germline_alignment: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AirrRecord {
    pub sequence_id: String,
    pub sequence: String,
    pub rev_comp: bool,
    pub v: Option<SegmentReport>,
    pub d: Option<SegmentReport>,
    pub d2: Option<SegmentReport>,
    pub j: Option<SegmentReport>,
    pub np1: Option<String>,
    pub np2: Option<String>,
    pub np3: Option<String>,
    pub path_score: Option<Score>,
    /// FR/CDR/CDR3/junction coordinates and productivity.
    pub regions: ReadRegions,
    /// Oriented query as ASCII (region sequences are slices of it).
    oriented: String,
}

/// Region columns in output order, with their accessor.
const REGION_COLUMNS: [(&str, fn(&ReadRegions) -> &Option<std::ops::Range<u32>>); 7] = [
    ("fwr1", |r| &r.fwr1),
    ("cdr1", |r| &r.cdr1),
    ("fwr2", |r| &r.fwr2),
    ("cdr2", |r| &r.cdr2),
    ("fwr3", |r| &r.fwr3),
    ("fwr4", |r| &r.fwr4),
    ("cdr3", |r| &r.cdr3),
];

/// Per-kind statistics for support values (indexed by [`SegmentKind::index`]).
pub type SupportModels = [Option<KarlinAltschul>; 3];

const SEGMENT_FIELDS: [&str; 11] = [
    "call",
    "score",
    "identity",
    "support",
    "cigar",
    "sequence_start",
    "sequence_end",
    "germline_start",
    "germline_end",
    "sequence_alignment",
    "germline_alignment",
];

fn push_run(cigar: &mut String, symbol: char, length: u32) {
    if length > 0 {
        let _ = write!(cigar, "{length}{symbol}");
    }
}

fn segment_report<C: AlleleCatalog>(
    call: &Call,
    called: &CalledSegment,
    catalog: &C,
    support: &SupportModels,
) -> Option<SegmentReport> {
    let dag_placement = &call.placements[called.placement];
    let kind = dag_placement.kind;
    let segment = catalog.segment(kind)?;
    let set = segment.set();
    let called_score = called.span.score;

    // Co-optimal alleles on the called interval. The primary name is the one
    // with the most germline support (unclipped score, the DAG's own tie-break),
    // then the longest exact run and most matches on the interval.
    let mut tied: Vec<_> = std::iter::once(dag_placement)
        .chain(
            call.placements
                .iter()
                .filter(|p| p.kind == kind)
                .filter(|p| {
                    p.clip
                        .span(called.entry, called.exit)
                        .is_some_and(|s| s.score == called_score)
                }),
        )
        .collect();
    tied.sort_by_key(|p| {
        let own = p
            .clip
            .span(called.entry, called.exit)
            .expect("tied placement clips to the called interval");
        let e = Evidence::from_columns(
            p.clip.columns(&own),
            &call.query,
            &set.allele(p.primary_allele()).codes,
        );
        (
            std::cmp::Reverse(p.score()),
            std::cmp::Reverse(e.longest_exact_run),
            std::cmp::Reverse(e.matches),
            p.primary_allele(),
        )
    });
    let placement = tied[0];
    let span = placement
        .clip
        .span(called.entry, called.exit)
        .expect("primary tied placement clips to the called interval");
    let reference = &set.allele(placement.primary_allele()).codes;
    let mut names: Vec<&str> = Vec::new();
    for p in &tied {
        for &allele in &p.alleles {
            for name in set.allele(allele).names() {
                if !names.contains(&name) {
                    names.push(name);
                }
            }
        }
    }

    // CIGAR and alignment strings.
    let mut cigar = String::new();
    push_run(&mut cigar, 'S', span.q_lo);
    push_run(&mut cigar, 'N', span.r_lo);
    let mut sequence_alignment = String::with_capacity(span.q_len() as usize + 8);
    let mut germline_alignment = String::with_capacity(span.r_len() as usize + 8);
    let mut run: Option<(char, u32)> = None;
    for (op, q, r) in placement.clip.columns(&span) {
        let (symbol, qs, rs) = match op {
            Op::Match => {
                let (a, b) = (call.query[q as usize], reference[r as usize]);
                let symbol = if a < N && a == b { '=' } else { 'X' };
                (
                    symbol,
                    char::from(decode_base(a)),
                    char::from(decode_base(b)),
                )
            }
            Op::Ins => ('I', char::from(decode_base(call.query[q as usize])), '-'),
            Op::Del => ('D', '-', char::from(decode_base(reference[r as usize]))),
        };
        sequence_alignment.push(qs);
        germline_alignment.push(rs);
        run = match run {
            Some((s, len)) if s == symbol => Some((s, len + 1)),
            Some((s, len)) => {
                push_run(&mut cigar, s, len);
                Some((symbol, 1))
            }
            None => Some((symbol, 1)),
        };
    }
    if let Some((s, len)) = run {
        push_run(&mut cigar, s, len);
    }
    push_run(&mut cigar, 'S', call.query.len() as u32 - span.q_hi);

    let evidence = Evidence::from_columns(placement.clip.columns(&span), &call.query, reference);
    let db_len: usize = set.alleles().iter().map(|a| a.codes.len()).sum();
    let evalue =
        support[kind.index()].map(|ka| ka.evalue(span.score, call.query.len(), db_len, set.len()));

    Some(SegmentReport {
        call: names.join(","),
        score: span.score,
        identity: evidence.identity(),
        support: evalue,
        cigar,
        sequence_start: span.q_lo + 1,
        sequence_end: span.q_hi,
        germline_start: span.r_lo + 1,
        germline_end: span.r_hi,
        sequence_alignment,
        germline_alignment,
    })
}

impl AirrRecord {
    pub fn from_call<C: AlleleCatalog>(
        sequence_id: &str,
        sequence: &str,
        call: &Call,
        catalog: &C,
        support: &SupportModels,
    ) -> Self {
        let mut record = AirrRecord {
            sequence_id: sequence_id.to_string(),
            sequence: sequence.to_string(),
            rev_comp: call.orientation.is_reverse(),
            v: None,
            d: None,
            d2: None,
            j: None,
            np1: None,
            np2: None,
            np3: None,
            path_score: None,
            regions: ReadRegions::from_call(call, catalog),
            oriented: decode(&call.query),
        };
        let Some(path) = &call.path else {
            return record;
        };
        record.path_score = Some(path.score.score);

        let mut junctions = Vec::with_capacity(3);
        for (i, called) in path.segments.iter().enumerate() {
            if i > 0 {
                let (from, to) = (path.segments[i - 1].span.q_hi, called.span.q_lo);
                junctions.push(if to > from {
                    decode(&call.query[from as usize..to as usize])
                } else {
                    String::new()
                });
            }
            let report = segment_report(call, called, catalog, support);
            let slot = match call.placements[called.placement].kind {
                SegmentKind::V => &mut record.v,
                SegmentKind::D if record.d.is_none() => &mut record.d,
                SegmentKind::D => &mut record.d2,
                SegmentKind::J => &mut record.j,
            };
            *slot = report;
        }
        let mut junctions = junctions.into_iter();
        record.np1 = junctions.next();
        record.np2 = junctions.next();
        record.np3 = junctions.next();
        record
    }

    pub fn tsv_header() -> String {
        let mut columns = vec![
            "sequence_id".to_string(),
            "sequence".to_string(),
            "rev_comp".to_string(),
        ];
        for prefix in ["v", "d", "d2", "j"] {
            columns.extend(SEGMENT_FIELDS.iter().map(|f| format!("{prefix}_{f}")));
        }
        for np in ["np1", "np2", "np3"] {
            columns.push(np.to_string());
            columns.push(format!("{np}_length"));
        }
        columns.push("path_score".to_string());
        for (name, _) in REGION_COLUMNS {
            columns.extend([
                name.to_string(),
                format!("{name}_start"),
                format!("{name}_end"),
            ]);
        }
        columns.extend(
            [
                "cdr3_aa",
                "junction",
                "junction_length",
                "junction_aa",
                "vj_in_frame",
                "stop_codon",
                "productive",
            ]
            .map(str::to_string),
        );
        columns.join("\t")
    }

    pub fn to_tsv(&self) -> String {
        let mut fields = vec![
            self.sequence_id.clone(),
            self.sequence.clone(),
            (if self.rev_comp { "T" } else { "F" }).to_string(),
        ];
        for segment in [&self.v, &self.d, &self.d2, &self.j] {
            match segment {
                Some(s) => fields.extend([
                    s.call.clone(),
                    s.score.to_string(),
                    format!("{:.4}", s.identity),
                    s.support.map_or_else(String::new, |e| format!("{e:.3e}")),
                    s.cigar.clone(),
                    s.sequence_start.to_string(),
                    s.sequence_end.to_string(),
                    s.germline_start.to_string(),
                    s.germline_end.to_string(),
                    s.sequence_alignment.clone(),
                    s.germline_alignment.clone(),
                ]),
                None => fields.extend(std::iter::repeat(String::new()).take(SEGMENT_FIELDS.len())),
            }
        }
        for np in [&self.np1, &self.np2, &self.np3] {
            match np {
                Some(s) => fields.extend([s.clone(), s.len().to_string()]),
                None => fields.extend([String::new(), String::new()]),
            }
        }
        fields.push(self.path_score.map_or_else(String::new, |s| s.to_string()));
        let slice =
            |range: &std::ops::Range<u32>| &self.oriented[range.start as usize..range.end as usize];
        for (_, get) in REGION_COLUMNS {
            match get(&self.regions) {
                Some(range) => fields.extend([
                    slice(range).to_string(),
                    (range.start + 1).to_string(),
                    range.end.to_string(),
                ]),
                None => fields.extend([String::new(), String::new(), String::new()]),
            }
        }
        let aa = |range: &std::ops::Range<u32>| {
            String::from_utf8(translate(
                &crate::alphabet::encode(slice(range).as_bytes()),
                0,
            ))
            .expect("amino-acid letters are ASCII")
        };
        fields.push(self.regions.cdr3.as_ref().map_or_else(String::new, aa));
        match &self.regions.junction {
            Some(range) => {
                fields.extend([slice(range).to_string(), range.len().to_string(), aa(range)])
            }
            None => fields.extend([String::new(), String::new(), String::new()]),
        }
        let flag = |value: Option<bool>| {
            value.map_or_else(String::new, |v| (if v { "T" } else { "F" }).to_string())
        };
        fields.extend([
            flag(self.regions.vj_in_frame),
            flag(self.regions.stop_codon),
            flag(self.regions.productive),
        ]);
        fields.join("\t")
    }
}
