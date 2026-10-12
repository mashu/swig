//! AIRR Rearrangement output in Swig's schema.
//!
//! Every row has the same columns, in the order Swig's interface and
//! command line have always read ([`AirrRecord::COLUMNS`]): the AIRR
//! Rearrangement fields, Swig's `c_prefix_identity`, `region_definition`,
//! `v_annotation_source`, `j_annotation_source` and `*_alternatives`, then
//! the engine's posterior probabilities. Without a C database the `c_`
//! columns are empty.
//!
//! Coordinates follow the AIRR convention: 1-based closed intervals on the
//! *oriented* query (the reverse complement of `sequence` when `rev_comp`).
//! `sequence` is the read as given, upper-cased, `U` read as `T`. CIGAR
//! strings use `S` for query bases before/after the segment, `N` for
//! germline bases before the alignment, and `=`/`X`/`I`/`D` for columns.
//! `*_identity` is matches over alignment columns (a fraction), `*_support`
//! an E-value.
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
//! the read's allele up to germline equivalence. `d_presence_posterior` is
//! the probability that the read has a D.
//!
//! `*_alternatives` lists the alleles the call does not name whose
//! placement, clipped to the called interval, scores within
//! `max(2, |S|/100)` of the call's `S`, as `name|score|identity|start|end|
//! germline_start|germline_end`, `;`-separated, best first.
//!
//! **Constant region.** With a C stage ([`crate::constant`]) the C call is
//! reported like the others, with three differences.
//! * `c_support` is the anchored E-value `N_C·(2s + 1)·e^(−λS)` — a union
//!   bound over C genes and start diagonals — not the Karlin–Altschul
//!   `K·m·n·e^(−λS)` of `v_`, `d_` and `j_support`; the columns are not
//!   comparable.
//! * Ties are on the fixed scale, and `c_posterior` weighs C's placements
//!   at C's own divergence, given that the read reaches C.
//! * Coordinates are those of the C alignment clipped to the J exit:
//!   `c_germline_start` 2 flags a reference whose base 1 repeats the J's
//!   last base (codon-aligned records).
//!
//! `c_prefix_identity` is the identity of the C germline from its first
//! base: the uncovered germline prefix compared base for base with the read
//! bases before the alignment, then the alignment's columns.
//!
//! A short C remnant matches the shared 5′ end of several genes (IGHG1,
//! IGHG2 and IGHG4 share `GCCTCCACCAAG`), so `c_call` then lists several
//! genes and its *first* name is merely first by the tie-break (longest
//! exact run, most matches, lowest allele id); read the isotype from the
//! whole list.

use std::fmt::Write as _;
use std::ops::Range;

use crate::adjudicate::{CallModel, CallPath, CalledSegment, Marginals};
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
use crate::translate::{UNKNOWN, codon};

/// An allele near-tied with a call (see the module documentation).
#[derive(Clone, Debug, PartialEq)]
pub struct Alternative {
    pub name: String,
    pub score: Score,
    pub identity: f64,
    /// 1-based closed on the oriented query, and on the germline.
    pub sequence_start: u32,
    pub sequence_end: u32,
    pub germline_start: u32,
    pub germline_end: u32,
}

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
    /// Near-tied alleles the call does not list, best first.
    pub alternatives: Vec<Alternative>,
    /// C only: identity of the germline from its first base.
    pub prefix_identity: Option<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AirrRecord {
    pub sequence_id: String,
    /// The read as given, upper-cased, `U` as `T`.
    pub sequence: String,
    /// FASTQ qualities, in the read's given orientation.
    pub quality: String,
    pub rev_comp: bool,
    /// The locus the path's alleles share; `None` when unknown.
    pub locus: Option<&'static str>,
    pub v: Option<SegmentReport>,
    pub d: Option<SegmentReport>,
    pub j: Option<SegmentReport>,
    /// The constant region ([`crate::constant`]).
    pub c: Option<SegmentReport>,
    pub np1: Option<String>,
    pub np2: Option<String>,
    /// V, D and J stitched in read order, N/P bases against `N`.
    pub sequence_alignment: String,
    pub germline_alignment: String,
    /// Posterior probability that the read contains a D.
    pub d_presence_posterior: Option<f64>,
    /// FR/CDR/CDR3/junction coordinates, frames and productivity.
    pub regions: ReadRegions,
    /// `IMGT` when the called V has region coordinates.
    pub region_definition: Option<&'static str>,
    /// Where the called V's and J's coordinates come from.
    pub v_annotation_source: Option<String>,
    pub j_annotation_source: Option<String>,
    /// Oriented query as ASCII (region sequences are slices of it).
    oriented: String,
}

/// Statistics for support values per rearranged kind ([`SegmentKind::REARRANGED`]),
/// indexed by [`SegmentKind::index`].
pub type SupportModels = [Option<KarlinAltschul>; 3];

/// A probability with six decimals, or ten when it lies within 10⁻⁴ of 0
/// or 1, so `0.999999` and `0.9999999` stay distinguishable. A distance
/// below 5·10⁻¹¹ rounds away at ten decimals, so it prints as an exact 0 or
/// 1 does: a last-bit difference between platforms' `exp` cannot change
/// the text.
pub(crate) fn probability(p: f64) -> String {
    let tail = p.min(1.0 - p);
    if (5e-11..1e-4).contains(&tail) {
        format!("{p:.10}")
    } else {
        format!("{p:.6}")
    }
}

pub(crate) fn push_run(cigar: &mut String, symbol: char, length: u32) {
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

/// The CIGAR, the sequence and germline alignment strings, and the evidence
/// of `placement`'s columns on `span` (`reference` its primary allele).
fn aligned(
    query: &[u8],
    placement: &Placement,
    span: &ClipSpan,
    reference: &[u8],
) -> (String, String, String, Evidence) {
    let mut cigar = String::new();
    push_run(&mut cigar, 'S', span.q_lo);
    push_run(&mut cigar, 'N', span.r_lo);
    let mut sequence_alignment = String::with_capacity(span.q_len() as usize + 8);
    let mut germline_alignment = String::with_capacity(span.r_len() as usize + 8);
    let mut run: Option<(char, u32)> = None;
    for (op, q, r) in placement.clip.columns(span) {
        let (symbol, qs, rs) = match op {
            Op::Match => {
                let (a, b) = (query[q as usize], reference[r as usize]);
                let symbol = if a < N && a == b { '=' } else { 'X' };
                (
                    symbol,
                    char::from(decode_base(a)),
                    char::from(decode_base(b)),
                )
            }
            Op::Ins => ('I', char::from(decode_base(query[q as usize])), '-'),
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
    push_run(&mut cigar, 'S', query.len() as u32 - span.q_hi);
    let evidence = Evidence::from_columns(placement.clip.columns(span), query, reference);
    (cigar, sequence_alignment, germline_alignment, evidence)
}

/// An identity, a fraction with up to six decimals, trailing zeros dropped.
pub(crate) fn fraction(value: f64) -> String {
    let mut text = format!("{value:.6}");
    while text.ends_with('0') {
        text.pop();
    }
    if text.ends_with('.') {
        text.pop();
    }
    text
}

/// An E-value as C's `%.6e` prints it (`1.234560e-05`).
fn evalue(value: f64) -> String {
    let text = format!("{value:.6e}");
    match text.split_once('e') {
        Some((mantissa, exponent)) => {
            let (sign, digits) = match exponent.strip_prefix('-') {
                Some(digits) => ('-', digits),
                None => ('+', exponent),
            };
            format!("{mantissa}e{sign}{digits:0>2}")
        }
        None => text,
    }
}

/// The read as given, upper-cased, `U` as `T`, and every symbol that is no
/// IUPAC nucleotide as `N`.
fn normalized(sequence: &str) -> String {
    sequence
        .bytes()
        .map(|b| match b.to_ascii_uppercase() {
            b'U' => 'T',
            c @ (b'A' | b'C' | b'G' | b'T' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B'
            | b'D' | b'H' | b'V' | b'N') => char::from(c),
            _ => 'N',
        })
        .collect()
}

/// Translate the ASCII `bases` from offset `frame`; with `gaps`, a codon
/// holding `-` reads `-`. Unknown codons read `X`.
fn translate_text(bases: &str, frame: usize, gaps: bool) -> String {
    let bytes = bases.as_bytes();
    let code = |b: u8| crate::alphabet::encode_base(b);
    bytes
        .get(frame..)
        .unwrap_or_default()
        .chunks_exact(3)
        .map(|c| {
            if gaps && c.contains(&b'-') {
                '-'
            } else {
                let aa = codon(code(c[0]), code(c[1]), code(c[2]));
                char::from(if aa == 0 { UNKNOWN } else { aa })
            }
        })
        .collect()
}

/// `S_alt ≥ S − max(2, |S|/100)`: the alternatives' score floor.
fn near(score: Score) -> Score {
    score - (score.abs() / 100).max(2)
}

/// The alleles of `kind` near-tied with `called` that `names` does not
/// list, each placement clipped to `[entry, exit)`, best first.
fn alternatives(
    placements: &[Placement],
    query: &[u8],
    set: &GermlineSet,
    called: &CalledSegment,
    names: &[&str],
) -> Vec<Alternative> {
    let kind = placements[called.placement].kind;
    let floor = near(called.span.score);
    let mut found: Vec<Alternative> = Vec::new();
    for p in placements.iter().filter(|p| p.kind == kind) {
        let Some(span) = p.clip.span(called.entry, called.exit) else {
            continue;
        };
        if span.score < floor {
            continue;
        }
        let reference = &set.allele(p.primary_allele()).codes;
        let identity = Evidence::from_columns(p.clip.columns(&span), query, reference).identity();
        for &allele in &p.alleles {
            for name in set.allele(allele).names() {
                if names.contains(&name) || found.iter().any(|a| a.name == name) {
                    continue;
                }
                found.push(Alternative {
                    name: name.to_string(),
                    score: span.score,
                    identity,
                    sequence_start: span.q_lo + 1,
                    sequence_end: span.q_hi,
                    germline_start: span.r_lo + 1,
                    germline_end: span.r_hi,
                });
            }
        }
    }
    found.sort_by_key(|a| std::cmp::Reverse(a.score));
    found
}

/// The report of `called`: the co-optimal names, the primary's clipped
/// alignment and coordinates, and the near-tied alternatives.
fn called_report(
    call: &Call,
    placements: &[Placement],
    set: &GermlineSet,
    co: &CoOptimal<'_>,
    called: &CalledSegment,
    support: Option<f64>,
    posterior: Option<f64>,
) -> SegmentReport {
    let query = &call.query;
    let placement = co.tied[0];
    let span = placement
        .clip
        .span(called.entry, called.exit)
        .expect("primary tied placement clips to the called interval");
    let reference = &set.allele(placement.primary_allele()).codes;
    let (cigar, sequence_alignment, germline_alignment, evidence) =
        aligned(query, placement, &span, reference);
    let names = co.names(set);
    SegmentReport {
        alternatives: alternatives(placements, query, set, called, &names),
        call: names.join(","),
        score: span.score,
        identity: evidence.identity(),
        support,
        posterior,
        cigar,
        sequence_start: span.q_lo + 1,
        sequence_end: span.q_hi,
        germline_start: span.r_lo + 1,
        germline_end: span.r_hi,
        sequence_alignment,
        germline_alignment,
        prefix_identity: None,
    }
}

/// The report of `called`, a segment of `path` (the call's, or another
/// reading of the same placements) with posterior `marginals`.
pub(crate) fn segment_report<C: AlleleCatalog>(
    call: &Call,
    path: &CallPath,
    marginals: Option<&Marginals>,
    called: &CalledSegment,
    catalog: &C,
    support: &SupportModels,
) -> Option<SegmentReport> {
    let kind = call.placements[called.placement].kind;
    assert!(
        kind.is_rearranged(),
        "a C placement reached the V(D)J report"
    );
    let segment = catalog.segment(kind)?;
    let set = segment.set();
    let per_read = path.model == CallModel::PerRead;
    let region = {
        let at = path.segments.iter().position(|s| s == called);
        let before = at
            .and_then(|i| i.checked_sub(1))
            .map(|i| path.segments[i].span.q_hi);
        let after = at
            .and_then(|i| path.segments.get(i + 1))
            .map(|s| s.span.q_lo);
        before.unwrap_or(0)..after.unwrap_or(call.query.len() as u32)
    };
    let co = co_optimal(&call.placements, &call.query, set, called, per_read, region);
    let posterior = marginals.map(|m| {
        co.counted
            .iter()
            .filter_map(|&i| m.placement.get(i))
            .sum::<f64>()
            // Rounding; and under a tandem grammar, a path holding two such
            // D placements, so the sum is then an upper bound.
            .min(1.0)
    });
    let evalue = support[kind.index()]
        .map(|ka| ka.evalue(called.span.score, call.query.len(), set.total_len()));
    Some(called_report(
        call,
        &call.placements,
        set,
        &co,
        called,
        evalue,
        posterior,
    ))
}

/// Identity of `allele` from its first germline base on `span`: the
/// germline bases before the alignment against the read bases before it,
/// then the alignment's columns. `None` when the read does not cover them.
fn prefix_identity(
    call: &Call,
    set: &GermlineSet,
    placement: &Placement,
    span: &ClipSpan,
) -> Option<f64> {
    let codes = &set.allele(placement.primary_allele()).codes;
    let origin = span.q_lo.checked_sub(span.r_lo)? as usize;
    let before = &call.query[origin..span.q_lo as usize];
    let leading = before
        .iter()
        .zip(codes)
        .filter(|&(&q, &r)| q < N && q == r)
        .count() as u32;
    let evidence = Evidence::from_columns(placement.clip.columns(span), &call.query, codes);
    let columns = span.r_lo + evidence.columns;
    (columns > 0).then(|| f64::from(leading + evidence.matches) / f64::from(columns))
}

/// The report of the read's C call: support is its anchored E-value
/// ([`crate::constant`]), not a K·m·n E-value.
fn constant_report<C: AlleleCatalog>(call: &Call, catalog: &C) -> Option<SegmentReport> {
    let c = call.c.as_ref()?;
    let set = catalog.segment(SegmentKind::C)?.set();
    let called = &c.segment;
    let co = co_optimal(
        &call.constant,
        &call.query,
        set,
        called,
        false,
        called.entry..called.exit,
    );
    let mut report = called_report(
        call,
        &call.constant,
        set,
        &co,
        called,
        Some(c.evalue),
        Some(c.posterior),
    );
    let span = co.tied[0]
        .clip
        .span(called.entry, called.exit)
        .expect("primary tied placement clips to the called interval");
    report.prefix_identity = prefix_identity(call, set, co.tied[0], &span);
    Some(report)
}

/// One TSV row under construction.
pub(crate) struct Row {
    pub(crate) line: String,
    /// Fields written so far.
    fields: usize,
}

impl Row {
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self {
            line: String::with_capacity(capacity),
            fields: 0,
        }
    }

    fn separate(&mut self) {
        if self.fields > 0 {
            self.line.push('\t');
        }
        self.fields += 1;
    }

    pub(crate) fn field(&mut self, value: &str) {
        self.separate();
        self.line.push_str(value);
    }

    /// A field of free text: tabs and line breaks become spaces.
    pub(crate) fn text(&mut self, value: &str) {
        self.separate();
        self.line.extend(value.chars().map(|c| {
            if matches!(c, '\t' | '\r' | '\n') {
                ' '
            } else {
                c
            }
        }));
    }

    pub(crate) fn some<T>(&mut self, value: Option<T>, show: impl FnOnce(T) -> String) {
        match value {
            Some(v) => self.field(&show(v)),
            None => self.field(""),
        }
    }

    pub(crate) fn flag(&mut self, value: Option<bool>) {
        self.field(match value {
            Some(true) => "T",
            Some(false) => "F",
            None => "",
        });
    }

    pub(crate) fn number(&mut self, value: Option<u32>) {
        self.some(value, |v| v.to_string());
    }
}

impl AirrRecord {
    /// The columns of every row, in order.
    pub const COLUMNS: [&'static str; 111] = [
        "sequence_id",
        "sequence",
        "quality",
        "sequence_aa",
        "rev_comp",
        "productive",
        "vj_in_frame",
        "stop_codon",
        "complete_vdj",
        "locus",
        "v_call",
        "d_call",
        "j_call",
        "c_call",
        "sequence_alignment",
        "sequence_alignment_aa",
        "germline_alignment",
        "germline_alignment_aa",
        "junction",
        "junction_aa",
        "np1",
        "np2",
        "cdr1",
        "cdr1_aa",
        "cdr2",
        "cdr2_aa",
        "cdr3",
        "cdr3_aa",
        "fwr1",
        "fwr1_aa",
        "fwr2",
        "fwr2_aa",
        "fwr3",
        "fwr3_aa",
        "fwr4",
        "fwr4_aa",
        "v_score",
        "v_identity",
        "v_support",
        "v_cigar",
        "d_score",
        "d_identity",
        "d_support",
        "d_cigar",
        "j_score",
        "j_identity",
        "j_support",
        "j_cigar",
        "c_score",
        "c_identity",
        "c_prefix_identity",
        "c_support",
        "c_cigar",
        "v_sequence_start",
        "v_sequence_end",
        "v_germline_start",
        "v_germline_end",
        "d_sequence_start",
        "d_sequence_end",
        "d_germline_start",
        "d_germline_end",
        "j_sequence_start",
        "j_sequence_end",
        "j_germline_start",
        "j_germline_end",
        "c_sequence_start",
        "c_sequence_end",
        "c_germline_start",
        "c_germline_end",
        "cdr1_start",
        "cdr1_end",
        "cdr2_start",
        "cdr2_end",
        "cdr3_start",
        "cdr3_end",
        "fwr1_start",
        "fwr1_end",
        "fwr2_start",
        "fwr2_end",
        "fwr3_start",
        "fwr3_end",
        "fwr4_start",
        "fwr4_end",
        "v_sequence_alignment",
        "d_sequence_alignment",
        "j_sequence_alignment",
        "c_sequence_alignment",
        "v_germline_alignment",
        "d_germline_alignment",
        "j_germline_alignment",
        "c_germline_alignment",
        "junction_length",
        "junction_aa_length",
        "np1_length",
        "np2_length",
        "v_frameshift",
        "j_frameshift",
        "d_frame",
        "sequence_frame",
        "region_definition",
        "v_annotation_source",
        "j_annotation_source",
        "v_alternatives",
        "d_alternatives",
        "j_alternatives",
        "c_alternatives",
        "v_posterior",
        "d_posterior",
        "j_posterior",
        "c_posterior",
        "d_presence_posterior",
    ];

    /// The record of `call` for the read `sequence_id` with symbols
    /// `sequence` (as given); `quality` is empty until set.
    pub fn from_call<C: AlleleCatalog>(
        sequence_id: &str,
        sequence: &str,
        call: &Call,
        catalog: &C,
        support: &SupportModels,
    ) -> Self {
        let oriented = decode(&call.query);
        let mut record = AirrRecord {
            sequence_id: sequence_id.to_string(),
            sequence: normalized(sequence),
            quality: String::new(),
            rev_comp: call.orientation.is_reverse(),
            locus: None,
            v: None,
            d: None,
            j: None,
            c: constant_report(call, catalog),
            np1: None,
            np2: None,
            sequence_alignment: String::new(),
            germline_alignment: String::new(),
            d_presence_posterior: None,
            regions: ReadRegions::from_call(call, catalog),
            region_definition: None,
            v_annotation_source: None,
            j_annotation_source: None,
            oriented,
        };
        let Some(path) = &call.path else {
            return record;
        };
        record.d_presence_posterior = call.marginals.as_ref().map(|m| m.d_present());

        // The kind and exit of the previous segment.
        let mut previous: Option<(SegmentKind, u32)> = None;
        let mut loci = crate::locus::Loci::ANY;
        for called in path.segments.iter() {
            let placement = &call.placements[called.placement];
            let kind = placement.kind;
            if let Some(set) = catalog.segment(kind).map(|s| s.set()) {
                let own = placement
                    .alleles
                    .iter()
                    .fold(crate::locus::Loci::NONE, |l, &a| {
                        l.union(set.allele(a).loci)
                    });
                loci = loci.intersection(own);
            }
            let from = previous.map_or(called.span.q_lo, |(_, at)| at);
            let to = called.span.q_lo;
            if to > from {
                record.sequence_alignment += &record.oriented[from as usize..to as usize];
                record
                    .germline_alignment
                    .extend(std::iter::repeat_n('N', (to - from) as usize));
            }
            if let Some((before, _)) = previous {
                // AIRR names an N/P region by the segment it follows.
                let np = match before {
                    SegmentKind::V => &mut record.np1,
                    _ => &mut record.np2,
                };
                *np = Some(record.oriented[from.min(to) as usize..to as usize].to_string());
            }
            let report = segment_report(
                call,
                path,
                call.marginals.as_ref(),
                called,
                catalog,
                support,
            );
            if let Some(report) = &report {
                record.sequence_alignment += &report.sequence_alignment;
                record.germline_alignment += &report.germline_alignment;
            }
            let slot = match kind {
                SegmentKind::V => &mut record.v,
                SegmentKind::D => &mut record.d,
                SegmentKind::J => &mut record.j,
                SegmentKind::C => unreachable!("segment_report refuses C"),
            };
            if slot.is_none() {
                *slot = report;
            }
            previous = Some((kind, called.span.q_hi));
        }
        record.locus = (loci != crate::locus::Loci::ANY)
            .then(|| loci.first())
            .flatten()
            .map(|l| l.name());
        let source = |kind: SegmentKind| {
            let segment = path
                .segments
                .iter()
                .find(|s| call.placements[s.placement].kind == kind)?;
            let allele = call.placements[segment.placement].primary_allele();
            let set = catalog.segment(kind)?.set();
            Some(set.allele(allele))
        };
        if let Some(v) = source(SegmentKind::V) {
            record.v_annotation_source = v.source.clone();
            record.region_definition =
                (v.anchor.is_some() || !v.regions.is_empty()).then_some("IMGT");
        }
        record.j_annotation_source = source(SegmentKind::J).and_then(|j| j.source.clone());
        record
    }

    /// `sequence` as rows report it: upper-cased, `U` as `T`, every symbol
    /// that is no IUPAC nucleotide as `N`.
    pub fn normalized(sequence: &str) -> String {
        normalized(sequence)
    }

    /// The junction (V Cys104 … J Trp/Phe118) as nucleotides of the oriented read.
    pub fn junction_sequence(&self) -> Option<&str> {
        let range = self.regions.junction.as_ref()?;
        self.oriented.get(range.start as usize..range.end as usize)
    }

    /// The header line of [`Self::to_tsv`], without a line break.
    pub fn tsv_header() -> String {
        Self::COLUMNS.join("\t")
    }

    /// The `[start, end)` slice of the oriented read.
    fn slice(&self, range: &Range<u32>) -> &str {
        &self.oriented[range.start as usize..range.end as usize]
    }

    /// The row: one value per column of [`Self::COLUMNS`], without a line
    /// break.
    pub fn to_tsv(&self) -> String {
        let r = &self.regions;
        let frame = r.frame.map(|f| f as usize);
        let mut row = Row::with_capacity(4096);
        row.text(&self.sequence_id);
        row.field(&self.sequence);
        row.text(&self.quality);
        row.some(frame, |f| translate_text(&self.oriented, f, false));
        row.flag(Some(self.rev_comp));
        row.flag(r.productive);
        row.flag(r.vj_in_frame);
        row.flag(r.stop_codon);
        row.flag(r.complete_vdj);
        row.field(self.locus.unwrap_or(""));
        let segments = [&self.v, &self.d, &self.j, &self.c];
        for s in segments {
            row.field(s.as_ref().map_or("", |s| &s.call));
        }
        // The stitched alignment translates from the Cys104 column (or the
        // first in-frame column without a junction).
        let alignment_frame = frame.and_then(|f| {
            let start = self.v.as_ref().map(|v| v.sequence_start - 1)?;
            let target = r.junction.as_ref().map(|j| j.start);
            let mut position = start;
            self.sequence_alignment
                .bytes()
                .enumerate()
                .find_map(|(column, b)| {
                    if b == b'-' {
                        return None;
                    }
                    let hit = match target {
                        Some(cys) => position == cys,
                        None => position as usize % 3 == f,
                    };
                    position += 1;
                    hit.then_some(column % 3)
                })
        });
        row.field(&self.sequence_alignment);
        row.some(alignment_frame, |f| {
            translate_text(&self.sequence_alignment, f, true)
        });
        row.field(&self.germline_alignment);
        row.some(alignment_frame, |f| {
            translate_text(&self.germline_alignment, f, true)
        });
        row.some(r.junction.as_ref(), |j| self.slice(j).to_string());
        row.some(r.junction.as_ref(), |j| {
            translate_text(self.slice(j), 0, false)
        });
        row.field(self.np1.as_deref().unwrap_or(""));
        row.field(self.np2.as_deref().unwrap_or(""));
        // A V region translates in the read's frame; CDR3 and FWR4 from
        // their first base, which the junction puts in frame.
        let region = |row: &mut Row, range: Option<&Range<u32>>, in_frame: bool| {
            row.some(range, |x| self.slice(x).to_string());
            row.some(range, |x| {
                let offset = match (in_frame, frame) {
                    (true, Some(f)) => (f + 3 - x.start as usize % 3) % 3,
                    _ => 0,
                };
                translate_text(self.slice(x), offset, false)
            });
        };
        region(&mut row, r.cdr1.as_ref(), true);
        region(&mut row, r.cdr2.as_ref(), true);
        region(&mut row, r.cdr3.as_ref(), false);
        region(&mut row, r.fwr1.as_ref(), true);
        region(&mut row, r.fwr2.as_ref(), true);
        region(&mut row, r.fwr3.as_ref(), true);
        region(&mut row, r.fwr4.as_ref(), false);
        for (s, kind) in segments.into_iter().zip(SegmentKind::ALL) {
            row.some(s.as_ref(), |s| s.score.to_string());
            row.some(s.as_ref(), |s| fraction(s.identity));
            if kind == SegmentKind::C {
                row.some(s.as_ref().and_then(|s| s.prefix_identity), fraction);
            }
            row.some(s.as_ref().and_then(|s| s.support), evalue);
            row.field(s.as_ref().map_or("", |s| &s.cigar));
        }
        for s in segments {
            row.number(s.as_ref().map(|s| s.sequence_start));
            row.number(s.as_ref().map(|s| s.sequence_end));
            row.number(s.as_ref().map(|s| s.germline_start));
            row.number(s.as_ref().map(|s| s.germline_end));
        }
        for range in [
            &r.cdr1, &r.cdr2, &r.cdr3, &r.fwr1, &r.fwr2, &r.fwr3, &r.fwr4,
        ] {
            row.number(range.as_ref().map(|x| x.start + 1));
            row.number(range.as_ref().map(|x| x.end));
        }
        for s in segments {
            row.field(s.as_ref().map_or("", |s| &s.sequence_alignment));
        }
        for s in segments {
            row.field(s.as_ref().map_or("", |s| &s.germline_alignment));
        }
        row.number(r.junction.as_ref().map(|j| j.end - j.start));
        row.number(r.junction.as_ref().map(|j| (j.end - j.start) / 3));
        row.some(self.np1.as_ref(), |np| np.len().to_string());
        row.some(self.np2.as_ref(), |np| np.len().to_string());
        row.flag(r.v_frameshift);
        row.flag(r.j_frameshift);
        row.number(r.d_frame.map(|f| f + 1));
        row.number(r.frame.map(|f| f + 1));
        row.field(self.region_definition.unwrap_or(""));
        row.field(self.v_annotation_source.as_deref().unwrap_or(""));
        row.field(self.j_annotation_source.as_deref().unwrap_or(""));
        for s in segments {
            let mut text = String::new();
            for (i, a) in s.iter().flat_map(|s| &s.alternatives).enumerate() {
                if i > 0 {
                    text.push(';');
                }
                let _ = write!(
                    text,
                    "{}|{}|{}|{}|{}|{}|{}",
                    a.name,
                    a.score,
                    fraction(a.identity),
                    a.sequence_start,
                    a.sequence_end,
                    a.germline_start,
                    a.germline_end
                );
            }
            row.field(&text);
        }
        for s in segments {
            row.some(s.as_ref().and_then(|s| s.posterior), probability);
        }
        row.some(self.d_presence_posterior, probability);
        row.line
    }
}

#[cfg(test)]
mod tests {
    use super::probability;

    #[test]
    fn probabilities_show_their_tail_but_not_rounding_noise() {
        assert_eq!(probability(0.5), "0.500000");
        assert_eq!(probability(0.99999), "0.9999900000");
        assert_eq!(probability(1e-7), "0.0000001000");
        // Within rounding of 0 or 1, the text is that of 0 or 1.
        assert_eq!(probability(1.0), "1.000000");
        assert_eq!(probability(1.0 - f64::EPSILON), "1.000000");
        assert_eq!(probability(0.0), "0.000000");
        assert_eq!(probability(1e-17), "0.000000");
    }
}
