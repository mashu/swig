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
//! `sequence_aa` translates the oriented read in `sequence_frame`, and
//! `stop_codon` is a stop codon in that translation (the Cys104 codon's
//! frame) between `v_sequence_start` and `j_sequence_end`.
//! `sequence_alignment_aa` and `germline_alignment_aa` translate each
//! alignment row in its own frame, the read's and the called V germline's
//! coding frame, a run of three gap columns as `-`; the shorter is padded
//! with `-` so the residues at the Cys104 column line up. Their stops can
//! differ from `stop_codon`, which reads the read without indel correction.
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

use std::fmt::{self, Write as _};
use std::ops::Range;

use crate::adjudicate::{CallModel, CallPath, CalledSegment, CoOptimal, Marginals, co_optimal};
use crate::align::Op;
use crate::alphabet::{N, decode_base};
use crate::clip::ClipSpan;
use crate::evidence::Evidence;
use crate::germline::{GermlineSet, SegmentKind};
use crate::locus::Locus;
use crate::pipeline::{AlleleCatalog, Call};
use crate::placement::Placement;
use crate::read_regions::ReadRegions;
use crate::scoring::Score;
use crate::stats::KarlinAltschul;
use crate::translate::codon;

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
    /// The primary allele's coding frame: the offset (0–2) of a codon start
    /// on the germline, from its anchor or given frame, when it has one.
    pub germline_frame: Option<u32>,
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
    pub locus: Option<Locus>,
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
    /// The oriented read's symbols, as `sequence` gives them (region
    /// sequences are slices of it).
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

/// The CIGAR, the sequence and germline alignment strings, and the evidence
/// of `placement`'s columns on `span` (`reference` its primary allele); the
/// sequence alignment quotes the oriented read's symbols `text`.
fn aligned(
    query: &[u8],
    text: &[u8],
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
                    char::from(text[q as usize]),
                    char::from(decode_base(b)),
                )
            }
            Op::Ins => ('I', char::from(text[q as usize]), '-'),
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
/// IUPAC nucleotide as `N`: one symbol per byte, as
/// [`crate::alphabet::encode`] reads it.
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

/// The reverse complement of [`normalized`] `bases`, IUPAC codes included.
fn reverse_complement_text(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|b| {
            char::from(match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                b'R' => b'Y',
                b'Y' => b'R',
                b'K' => b'M',
                b'M' => b'K',
                b'B' => b'V',
                b'V' => b'B',
                b'D' => b'H',
                b'H' => b'D',
                // S, W and N are their own complements.
                other => other,
            })
        })
        .collect()
}

/// Translate the ASCII `bases` from offset `frame`; unknown codons read `X`.
fn translate_text(bases: &str, frame: usize) -> String {
    let code = crate::alphabet::encode_base;
    bases
        .as_bytes()
        .get(frame..)
        .unwrap_or_default()
        .chunks_exact(3)
        .map(|c| char::from(codon(code(c[0]), code(c[1]), code(c[2]))))
        .collect()
}

/// Translate the gapped alignment row `text` in the frame of its own bases:
/// skip `skip` bases, then every three bases make a codon, whatever gap
/// columns lie between them, so an indel that is not a multiple of 3 never
/// shifts the row out of its register. Every run of three gap columns
/// reads `-`, so an in-frame indel keeps both rows aligned. Also returns
/// the number of residues before column `anchor`.
fn translate_row(text: &str, skip: usize, anchor: usize) -> (String, usize) {
    let code = crate::alphabet::encode_base;
    let mut protein = String::with_capacity(text.len() / 3 + 1);
    let mut at_anchor = None;
    let (mut skip, mut codon_bases, mut filled, mut gaps) = (skip, [0u8; 3], 0, 0);
    for (column, b) in text.bytes().enumerate() {
        if column == anchor {
            at_anchor = Some(protein.len());
        }
        if b == b'-' {
            gaps += 1;
            if gaps == 3 {
                protein.push('-');
                gaps = 0;
            }
            continue;
        }
        gaps = 0;
        if skip > 0 {
            skip -= 1;
            continue;
        }
        codon_bases[filled] = code(b);
        filled += 1;
        if filled == 3 {
            let [a, b, c] = codon_bases;
            protein.push(char::from(codon(a, b, c)));
            filled = 0;
        }
    }
    let residues = at_anchor.unwrap_or(protein.len());
    (protein, residues)
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
/// alignment and coordinates, and the near-tied alternatives; no support or
/// posterior yet.
fn called_report(
    call: &Call,
    text: &[u8],
    placements: &[Placement],
    set: &GermlineSet,
    co: &CoOptimal<'_>,
    called: &CalledSegment,
) -> SegmentReport {
    let query = &call.query;
    let placement = co.tied[0];
    let span = placement
        .clip
        .span(called.entry, called.exit)
        .expect("primary tied placement clips to the called interval");
    let allele = set.allele(placement.primary_allele());
    let (cigar, sequence_alignment, germline_alignment, evidence) =
        aligned(query, text, placement, &span, &allele.codes);
    let names = co.names(set);
    SegmentReport {
        alternatives: alternatives(placements, query, set, called, &names),
        call: names.join(","),
        score: span.score,
        identity: evidence.identity(),
        support: None,
        posterior: None,
        cigar,
        sequence_start: span.q_lo + 1,
        sequence_end: span.q_hi,
        germline_start: span.r_lo + 1,
        germline_end: span.r_hi,
        sequence_alignment,
        germline_alignment,
        germline_frame: allele.anchor.map(|a| a % 3).or(allele.frame.map(u32::from)),
        prefix_identity: None,
    }
}

/// The report of `called`, a segment of `path` (the call's, or another
/// reading of the same placements) with posterior `marginals`; `text` is
/// the oriented read.
fn segment_report<C: AlleleCatalog>(
    call: &Call,
    text: &[u8],
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
    let evalue = support[kind.index()].map(|ka| {
        ka.evalue(
            called.span.score,
            call.query.len(),
            set.total_len(),
            set.len(),
        )
    });
    Some(SegmentReport {
        support: evalue,
        posterior,
        ..called_report(call, text, &call.placements, set, &co, called)
    })
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

/// The report of the read's C call (`text` the oriented read): support is
/// its anchored E-value ([`crate::constant`]), not a K·m·n E-value.
fn constant_report<C: AlleleCatalog>(
    call: &Call,
    text: &[u8],
    catalog: &C,
) -> Option<SegmentReport> {
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
    let mut report = SegmentReport {
        support: Some(c.evalue),
        posterior: Some(c.posterior),
        ..called_report(call, text, &call.constant, set, &co, called)
    };
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

    /// A field `write` writes into the line.
    pub(crate) fn write(&mut self, write: impl FnOnce(&mut String)) {
        self.separate();
        write(&mut self.line);
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

    /// A field holding `value`, written in place; empty for `None`.
    pub(crate) fn number(&mut self, value: Option<impl fmt::Display>) {
        self.separate();
        if let Some(value) = value {
            let _ = write!(self.line, "{value}");
        }
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
    /// `sequence` (as given, the text `call.query` was encoded from: the
    /// alignments, regions and junction quote its symbols, IUPAC codes
    /// included); `quality` is empty until set.
    ///
    /// # Panics
    ///
    /// If `sequence` is not as long as `call.query`.
    pub fn from_call<C: AlleleCatalog>(
        sequence_id: &str,
        sequence: &str,
        call: &Call,
        catalog: &C,
        support: &SupportModels,
    ) -> Self {
        assert_eq!(
            sequence.len(),
            call.query.len(),
            "the read's symbols are those the call's query was encoded from"
        );
        let sequence = normalized(sequence);
        let oriented = if call.orientation.is_reverse() {
            reverse_complement_text(&sequence)
        } else {
            sequence.clone()
        };
        let mut record = AirrRecord {
            sequence_id: sequence_id.to_string(),
            sequence,
            quality: String::new(),
            rev_comp: call.orientation.is_reverse(),
            locus: None,
            v: None,
            d: None,
            j: None,
            c: constant_report(call, oriented.as_bytes(), catalog),
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
                record.oriented.as_bytes(),
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
            .flatten();
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

    /// Empty the C columns unless `c_prefix_identity` reaches `min`
    /// (`--c-prefix-identity`).
    pub fn gate_c_prefix(&mut self, min: f64) {
        if self
            .c
            .as_ref()
            .is_some_and(|c| c.prefix_identity.is_none_or(|p| p < min))
        {
            self.c = None;
        }
    }

    /// The oriented read's symbols, as `sequence` gives them: the text the
    /// coordinates refer to.
    pub(crate) fn oriented(&self) -> &str {
        &self.oriented
    }

    /// The column of `sequence_alignment`, which starts at read base
    /// `start`, that holds read base `position`.
    fn alignment_column(&self, start: u32, position: u32) -> Option<usize> {
        let mut at = start;
        self.sequence_alignment
            .bytes()
            .enumerate()
            .filter(|&(_, b)| b != b'-')
            .find_map(|(column, _)| {
                at += 1;
                (at - 1 == position).then_some(column)
            })
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
        row.some(frame, |f| translate_text(&self.oriented, f));
        row.flag(Some(self.rev_comp));
        row.flag(r.productive);
        row.flag(r.vj_in_frame);
        row.flag(r.stop_codon);
        row.flag(r.complete_vdj);
        row.field(self.locus.map_or("", Locus::name));
        let segments = [&self.v, &self.d, &self.j, &self.c];
        for s in segments {
            row.field(s.as_ref().map_or("", |s| &s.call));
        }
        // Each stitched row translates in its own frame, the read's
        // (`sequence_frame`) and the called V germline's coding frame, and
        // the shorter is padded with `-` so the residues at the Cys104
        // column (the first column without a junction) line up.
        let translated = self.v.as_ref().zip(frame).map(|(v, f)| {
            let phase = |frame: usize, start: u32| (frame + 3 - start as usize % 3) % 3;
            let read_skip = phase(f, v.sequence_start - 1);
            let germline_skip = phase(v.germline_frame.unwrap_or(0) as usize, v.germline_start - 1);
            let anchor = r
                .junction
                .as_ref()
                .and_then(|j| self.alignment_column(v.sequence_start - 1, j.start))
                .unwrap_or(0);
            let (read, read_at) = translate_row(&self.sequence_alignment, read_skip, anchor);
            let (germline, germline_at) =
                translate_row(&self.germline_alignment, germline_skip, anchor);
            let at = read_at.max(germline_at);
            (
                "-".repeat(at - read_at) + &read,
                "-".repeat(at - germline_at) + &germline,
            )
        });
        row.field(&self.sequence_alignment);
        row.field(translated.as_ref().map_or("", |t| &t.0));
        row.field(&self.germline_alignment);
        row.field(translated.as_ref().map_or("", |t| &t.1));
        row.field(r.junction.as_ref().map_or("", |j| self.slice(j)));
        row.some(r.junction.as_ref(), |j| translate_text(self.slice(j), 0));
        row.field(self.np1.as_deref().unwrap_or(""));
        row.field(self.np2.as_deref().unwrap_or(""));
        // A V region translates in the read's frame; CDR3 and FWR4 from
        // their first base, which the junction puts in frame.
        let region = |row: &mut Row, range: Option<&Range<u32>>, in_frame: bool| {
            row.field(range.map_or("", |x| self.slice(x)));
            row.some(range, |x| {
                let offset = match (in_frame, frame) {
                    (true, Some(f)) => (f + 3 - x.start as usize % 3) % 3,
                    _ => 0,
                };
                translate_text(self.slice(x), offset)
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
            row.number(s.as_ref().map(|s| s.score));
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
        row.number(self.np1.as_ref().map(String::len));
        row.number(self.np2.as_ref().map(String::len));
        row.flag(r.v_frameshift);
        row.flag(r.j_frameshift);
        row.number(r.d_frame.map(|f| f + 1));
        row.number(r.frame.map(|f| f + 1));
        row.field(self.region_definition.unwrap_or(""));
        row.field(self.v_annotation_source.as_deref().unwrap_or(""));
        row.field(self.j_annotation_source.as_deref().unwrap_or(""));
        for s in segments {
            row.write(|text| {
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
            });
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
    use std::collections::HashMap;

    use super::*;
    use crate::alphabet::{decode, encode};
    use crate::annotate::AnnotateOptions;
    use crate::config::DefaultAssigner;
    use crate::pipeline::VdjAssigner;
    use crate::prepare::prepare;
    use crate::profiles::CallingProfile;
    use crate::sim::{SyntheticDatabase, SyntheticSpec};
    use crate::translate::translate;

    /// The default assigner over the synthetic database.
    fn synthetic_assigner() -> DefaultAssigner {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 5);
        let database = prepare(db.v, Some(db.d), db.j, AnnotateOptions::default()).unwrap();
        let config = CallingProfile::ROptimized.config();
        DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap()
    }

    /// The bases of the first allele of `kind`, and its anchor.
    fn germline(assigner: &DefaultAssigner, kind: SegmentKind) -> (String, Option<u32>) {
        let allele = assigner.segment(kind).unwrap().set().allele(0);
        (decode(&allele.codes), allele.anchor)
    }

    /// The row of `read`, by column name.
    fn row(assigner: &DefaultAssigner, read: &str) -> HashMap<&'static str, String> {
        let call = assigner.assign(&encode(read.as_bytes()), &mut assigner.workspace());
        let record = AirrRecord::from_call("q", read, &call, assigner, &[None, None, None]);
        AirrRecord::COLUMNS
            .into_iter()
            .zip(record.to_tsv().split('\t').map(str::to_string))
            .collect()
    }

    /// An indel that is not a multiple of 3 shifts neither translated
    /// alignment row out of its own frame: `germline_alignment_aa` reads the
    /// germline's protein and `sequence_alignment_aa` the read's, as
    /// `sequence_aa` does.
    #[test]
    fn translated_alignment_rows_keep_their_own_frames() {
        let assigner = synthetic_assigner();
        let (v, anchor) = germline(&assigner, SegmentKind::V);
        let (d, _) = germline(&assigner, SegmentKind::D);
        let (j, _) = germline(&assigner, SegmentKind::J);
        let anchor = anchor.expect("an annotated V") as usize;
        let protein = |bases: &str, frame: usize| {
            String::from_utf8(translate(&encode(bases.as_bytes()), frame)).unwrap()
        };
        let germline_v = protein(&v[..anchor], anchor % 3);
        let inserted = format!("{}T{}", &v[..120], &v[120..]);
        let deleted_and_inserted = format!("{}{}A{}", &v[..100], &v[101..200], &v[200..]);
        for edited in [inserted, deleted_and_inserted] {
            let read = format!("{edited}GGTCA{d}CCATG{j}");
            let row = row(&assigner, &read);
            assert!(row["v_cigar"].contains('I'), "{}", row["v_cigar"]);
            let germline_aa = row["germline_alignment_aa"].trim_start_matches('-');
            assert!(
                germline_aa.starts_with(&germline_v),
                "germline_alignment_aa {germline_aa} against the germline {germline_v}"
            );
            let read_aa = row["sequence_alignment_aa"].replace('-', "");
            assert!(
                row["sequence_aa"].contains(&read_aa),
                "sequence_alignment_aa {read_aa} not in sequence_aa {}",
                row["sequence_aa"]
            );
        }
    }

    /// IUPAC codes in a read reach every column that quotes its bases, on
    /// either strand: each such column, gaps removed, is a slice of the
    /// oriented read.
    #[test]
    fn iupac_codes_reach_alignments_regions_and_junction() {
        let assigner = synthetic_assigner();
        let (v, anchor) = germline(&assigner, SegmentKind::V);
        let (d, _) = germline(&assigner, SegmentKind::D);
        let (j, _) = germline(&assigner, SegmentKind::J);
        let anchor = anchor.expect("an annotated V") as usize;
        let mut read = format!("{v}GGTCA{d}CCATG{j}").into_bytes();
        let n = read.len();
        for (at, code) in [(20, b'Y'), (anchor - 15, b'K'), (n - 6, b'M')] {
            read[at] = code;
        }
        let forward = String::from_utf8(read).unwrap();
        for (given, reverse) in [
            (forward.clone(), false),
            (reverse_complement_text(&forward), true),
        ] {
            let row = row(&assigner, &given);
            assert_eq!(row["rev_comp"], if reverse { "T" } else { "F" });
            assert_eq!(row["sequence"], given);
            for column in [
                "sequence_alignment",
                "v_sequence_alignment",
                "j_sequence_alignment",
                "junction",
                "cdr3",
                "fwr1",
                "fwr3",
                "fwr4",
            ] {
                let bases = row[column].replace('-', "");
                assert!(!bases.is_empty(), "{column} is empty");
                assert!(
                    forward.contains(&bases),
                    "{column} {bases} is not read from {forward}"
                );
            }
            for (column, code) in [("fwr1", 'Y'), ("fwr3", 'K'), ("sequence_alignment", 'M')] {
                assert!(row[column].contains(code), "{column} {}", row[column]);
            }
        }
    }

    /// The C gate keeps C where `c_prefix_identity` reaches it, and only
    /// there.
    #[test]
    fn the_c_gate_keeps_c_from_its_prefix_identity_up() {
        let assigner = synthetic_assigner();
        let (v, _) = germline(&assigner, SegmentKind::V);
        let (j, _) = germline(&assigner, SegmentKind::J);
        let read = format!("{v}GGTCA{j}");
        let call = assigner.assign(&encode(read.as_bytes()), &mut assigner.workspace());
        let record = AirrRecord::from_call("q", &read, &call, &assigner, &[None, None, None]);
        let with_c = |prefix_identity: Option<f64>| {
            let c = record.v.clone().map(|v| SegmentReport {
                prefix_identity,
                ..v
            });
            AirrRecord {
                c,
                ..record.clone()
            }
        };
        for (identity, kept) in [(Some(0.97), true), (Some(0.969), false), (None, false)] {
            let mut gated = with_c(identity);
            gated.gate_c_prefix(0.97);
            assert_eq!(gated.c.is_some(), kept, "{identity:?}");
        }
        let mut without = with_c(None);
        without.c = None;
        without.gate_c_prefix(0.97);
        assert_eq!(without, AirrRecord { c: None, ..record });
    }

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
