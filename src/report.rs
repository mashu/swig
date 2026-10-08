//! AIRR Rearrangement output.
//!
//! Coordinates follow the AIRR convention: 1-based closed intervals on the
//! *oriented* query (the reverse complement of `sequence` when `rev_comp`).
//! CIGAR strings use `S` for query bases before/after the segment, `N` for
//! germline bases before the alignment, and `=`/`X`/`I`/`D` for columns.
//!
//! Co-optimal labels are recomputed on the called cut interval: every
//! placement of the same kind whose clip to `[entry, exit)` is worth exactly
//! the called clip under the call's model contributes its alleles. On the
//! fixed scale that is an equal score; under a per-read call
//! ([`crate::adjudicate::CallModel::PerRead`]) also equal matches and
//! mismatches, which makes the value `λ_r·S_sub + μ_r·n + λ0·S_gap` equal at
//! every divergence. As the note says, this means "equal value under this
//! caller", not equal biological probability. The
//! probability is the `posterior` column ([`crate::adjudicate::Marginals`]):
//! the posterior mass of every placement of that kind that is tied on the
//! interval, or whose germline bases between the call's neighbours the
//! primary allele contains: the probability that the reported names include
//! the read's allele up to germline equivalence. `d_presence_posterior` is the probability that the read has a D.

use std::fmt::Write as _;
use std::ops::Range;

use crate::adjudicate::{CallModel, CalledSegment};
use crate::align::Op;
use crate::alphabet::{N, decode, decode_base};
use crate::clip::ClipSpan;
use crate::evidence::Evidence;
use crate::germline::{GermlineSet, SegmentKind};
use crate::pipeline::{AlleleCatalog, Call};
use crate::placement::Placement;
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
    /// Posterior probability that the reported names include the read's
    /// allele up to germline equivalence (see the module documentation).
    pub posterior: Option<f64>,
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
    /// Posterior probability that the read contains a D segment.
    pub d_presence_posterior: Option<f64>,
    /// FR/CDR/CDR3/junction coordinates and productivity.
    pub regions: ReadRegions,
    /// Oriented query as ASCII (region sequences are slices of it).
    oriented: String,
}

/// Statistics for support values per rearranged kind ([`SegmentKind::REARRANGED`]),
/// indexed by [`SegmentKind::index`].
pub type SupportModels = [Option<KarlinAltschul>; 3];

const SEGMENT_FIELDS: [&str; 12] = [
    "call",
    "score",
    "identity",
    "support",
    "posterior",
    "cigar",
    "sequence_start",
    "sequence_end",
    "germline_start",
    "germline_end",
    "sequence_alignment",
    "germline_alignment",
];

/// A probability with six significant digits of its distance from 0 or 1,
/// so `0.999999` and `0.9999999` stay distinguishable.
fn probability(p: f64) -> String {
    let tail = p.min(1.0 - p);
    if tail > 0.0 && tail < 1e-4 {
        format!("{p:.10}")
    } else {
        format!("{p:.6}")
    }
}

fn push_run(cigar: &mut String, symbol: char, length: u32) {
    if length > 0 {
        let _ = write!(cigar, "{length}{symbol}");
    }
}

/// The placements tied with a called segment on its interval, and those its
/// posterior counts (see the module documentation).
pub(crate) struct CoOptimal<'a> {
    /// Tied placements, the primary — reported first — at index 0.
    pub tied: Vec<&'a Placement>,
    /// Indices of the placements of the called kind whose probability the
    /// posterior sums: tied with the call, or carried by its primary allele.
    pub counted: Vec<usize>,
}

impl CoOptimal<'_> {
    /// The tied placements' allele names, primary first, each once.
    pub fn names<'s>(&self, set: &'s GermlineSet) -> Vec<&'s str> {
        let mut names: Vec<&str> = Vec::new();
        for p in &self.tied {
            for &allele in &p.alleles {
                for name in set.allele(allele).names() {
                    if !names.contains(&name) {
                        names.push(name);
                    }
                }
            }
        }
        names
    }
}

/// Co-optimal alleles of `called` among `placements` (of every kind; only
/// the called one's count), `set` holding their alleles.
///
/// Tied: worth exactly the call on its interval `[entry, exit)` — an equal
/// score, and under a per-read call (`per_read`) equal matches and
/// mismatches. The primary is the one with the most germline support
/// (unclipped score, the DAG's own tie-break), then the longest exact run
/// and most matches on the interval, then the lowest allele id; a call
/// chosen by its allele's posterior marginal stays first.
///
/// Counted by the posterior: the tied placements, and every placement that
/// explains most of the read bases the call explains (more than half of the
/// called span, within `region`, the read the call's neighbours leave to
/// it) with germline bases there that the primary allele contains. If the
/// read came from such a placement, the call still carries its germline
/// bases. A placement that only reaches into the region by a few bases (a
/// D inside the called J, a tandem D's partner) is another hypothesis, so
/// it does not count.
pub(crate) fn co_optimal<'a>(
    placements: &'a [Placement],
    query: &[u8],
    set: &GermlineSet,
    called: &CalledSegment,
    per_read: bool,
    region: Range<u32>,
) -> CoOptimal<'a> {
    let called_placement = &placements[called.placement];
    let kind = called_placement.kind;
    let called_score = called.span.score;
    let columns = |p: &Placement, span: &ClipSpan| {
        let reference = &set.allele(p.primary_allele()).codes;
        let e = Evidence::from_columns(p.clip.columns(span), query, reference);
        (e.matches, e.mismatches)
    };
    let called_columns = per_read.then(|| columns(called_placement, &called.span));
    let tied_indices: Vec<usize> = (0..placements.len())
        .filter(|&i| {
            let p = &placements[i];
            p.kind == kind
                && p.clip.span(called.entry, called.exit).is_some_and(|s| {
                    s.score == called_score && called_columns.is_none_or(|c| columns(p, &s) == c)
                })
        })
        .collect();
    // The called placement first, so it wins ties of the sort key; one
    // chosen by its allele's posterior marginal stays first.
    let mut tied: Vec<&Placement> = std::iter::once(called_placement)
        .chain(
            tied_indices
                .iter()
                .filter(|&&i| i != called.placement)
                .map(|&i| &placements[i]),
        )
        .collect();
    tied[usize::from(called.by_marginal)..].sort_by_cached_key(|p| {
        let own = p
            .clip
            .span(called.entry, called.exit)
            .expect("tied placement clips to the called interval");
        let e = Evidence::from_columns(
            p.clip.columns(&own),
            query,
            &set.allele(p.primary_allele()).codes,
        );
        (
            std::cmp::Reverse(p.score()),
            std::cmp::Reverse(e.longest_exact_run),
            std::cmp::Reverse(e.matches),
            p.primary_allele(),
        )
    });
    let reference = &set.allele(tied[0].primary_allele()).codes;
    let called_bases = called.span.q_lo..called.span.q_hi;
    let carried = |p: &Placement| {
        let Some(own) = p.clip.span(region.start, region.end) else {
            return false;
        };
        let shared = own
            .q_hi
            .min(called_bases.end)
            .saturating_sub(own.q_lo.max(called_bases.start));
        if 2 * shared <= called_bases.len() as u32 {
            return false;
        }
        let codes = &set.allele(p.primary_allele()).codes;
        let bases = &codes[own.r_lo as usize..own.r_hi as usize];
        bases.is_empty() || reference.windows(bases.len()).any(|w| w == bases)
    };
    let counted = (0..placements.len())
        .filter(|&i| {
            let p = &placements[i];
            p.kind == kind && (tied_indices.contains(&i) || carried(p))
        })
        .collect();
    CoOptimal { tied, counted }
}

fn segment_report<C: AlleleCatalog>(
    call: &Call,
    called: &CalledSegment,
    catalog: &C,
    support: &SupportModels,
) -> Option<SegmentReport> {
    let kind = call.placements[called.placement].kind;
    debug_assert!(
        kind.is_rearranged(),
        "a C placement reached the V(D)J report"
    );
    let segment = catalog.segment(kind)?;
    let set = segment.set();
    let per_read = call
        .path
        .as_ref()
        .is_some_and(|path| path.model == CallModel::PerRead);
    let region = call
        .path
        .as_ref()
        .map_or(0..call.query.len() as u32, |path| {
            let at = path.segments.iter().position(|s| s == called);
            let before = at
                .and_then(|i| i.checked_sub(1))
                .map(|i| path.segments[i].span.q_hi);
            let after = at
                .and_then(|i| path.segments.get(i + 1))
                .map(|s| s.span.q_lo);
            before.unwrap_or(0)..after.unwrap_or(call.query.len() as u32)
        });
    let co = co_optimal(&call.placements, &call.query, set, called, per_read, region);
    let placement = co.tied[0];
    let span = placement
        .clip
        .span(called.entry, called.exit)
        .expect("primary tied placement clips to the called interval");
    let reference = &set.allele(placement.primary_allele()).codes;
    let posterior = call.marginals.as_ref().map(|m| {
        co.counted
            .iter()
            .filter_map(|&i| m.placement.get(i))
            .sum::<f64>()
            // Rounding; and under a tandem grammar, a path holding two such
            // D placements, so the sum is then an upper bound.
            .min(1.0)
    });
    let names = co.names(set);

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
    let evalue = support[kind.index()]
        .map(|ka| ka.evalue(span.score, call.query.len(), set.total_len(), set.len()));

    Some(SegmentReport {
        call: names.join(","),
        score: span.score,
        identity: evidence.identity(),
        support: evalue,
        posterior,
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
            d_presence_posterior: None,
            regions: ReadRegions::from_call(call, catalog),
            oriented: decode(&call.query),
        };
        let Some(path) = &call.path else {
            return record;
        };
        record.path_score = Some(path.score.score);
        record.d_presence_posterior = call.marginals.as_ref().map(|m| m.d_present());

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
                SegmentKind::C => unreachable!("C placements never enter a V(D)J path"),
            };
            *slot = report;
        }
        let mut junctions = junctions.into_iter();
        record.np1 = junctions.next();
        record.np2 = junctions.next();
        record.np3 = junctions.next();
        record
    }

    /// The junction (V Cys104 … J Trp/Phe118) as nucleotides of the oriented read.
    pub fn junction_sequence(&self) -> Option<&str> {
        let range = self.regions.junction.as_ref()?;
        self.oriented.get(range.start as usize..range.end as usize)
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
        columns.push("d_presence_posterior".to_string());
        for name in ReadRegions::NAMES {
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
                    s.posterior.map_or_else(String::new, probability),
                    s.cigar.clone(),
                    s.sequence_start.to_string(),
                    s.sequence_end.to_string(),
                    s.germline_start.to_string(),
                    s.germline_end.to_string(),
                    s.sequence_alignment.clone(),
                    s.germline_alignment.clone(),
                ]),
                None => fields.extend(std::iter::repeat_n(String::new(), SEGMENT_FIELDS.len())),
            }
        }
        for np in [&self.np1, &self.np2, &self.np3] {
            match np {
                Some(s) => fields.extend([s.clone(), s.len().to_string()]),
                None => fields.extend([String::new(), String::new()]),
            }
        }
        fields.push(self.path_score.map_or_else(String::new, |s| s.to_string()));
        fields.push(
            self.d_presence_posterior
                .map_or_else(String::new, probability),
        );
        let slice =
            |range: &std::ops::Range<u32>| &self.oriented[range.start as usize..range.end as usize];
        for interval in self.regions.intervals() {
            match interval {
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
