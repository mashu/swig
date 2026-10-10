//! Swig's double-D screen: VDDJ evidence beside the ordinary call.
//!
//! Screening leaves every AIRR row unchanged. Each read whose call places V
//! and J is read a second time over the same placements under the tandem-D
//! grammar ([`crate::adjudicate::VdjGrammar::TANDEM_D`]): VJ, VDJ and VDDJ
//! compete in one placement DAG, a second D paying its own `P_D` and the
//! tandem penalty on top of its alignment score. A read whose best tandem
//! path holds two D segments is reported in a sidecar table when the pair
//! also passes the screen's evidence rules:
//!
//! - **span**: in `long-span` mode, the call's V exit to J entry is at least
//!   `min_vj_span` bases;
//! - **seed**: each D holds an exact run of at least `seed_length` bases;
//! - **gain**: the two D scores together beat the best single D placement
//!   on the V–J span by at least `min_score_gain`;
//! - **pseudo-tandem**: the read from D1's start to D2's end, `pseudo_trim`
//!   bases shorter, is more than `pseudo_mismatches` Hamming mismatches from
//!   every substring of a D allele of the locus (the IgScout Δ-distance
//!   test), so no single D explains it.
//!
//! The sidecar row of read `i` of a batch starts with
//! `swig_batch_record_index` `i`, and `swig_double_d_posterior` is the
//! posterior probability of two D under the tandem grammar.

use crate::alphabet::decode;
use crate::evidence::Evidence;
use crate::germline::SegmentKind;
use crate::locus::{Loci, Locus};
use crate::pipeline::{AlleleCatalog, Call};
use crate::report::{Row, SegmentReport, SupportModels, fraction, probability, segment_report};
use crate::scoring::Score;

/// When the screen runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DoubleDMode {
    Off,
    /// Every read whose call places V and J.
    All,
    /// Only reads whose V–J span is at least `min_vj_span`.
    LongSpan,
}

impl DoubleDMode {
    pub fn name(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::All => "all",
            Self::LongSpan => "long_span",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "off" => Some(Self::Off),
            "all" => Some(Self::All),
            "long-span" | "long_span" => Some(Self::LongSpan),
            _ => None,
        }
    }
}

/// The screen's settings (see the module documentation).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DoubleDOptions {
    pub mode: DoubleDMode,
    pub min_vj_span: u32,
    pub seed_length: u32,
    pub pseudo_trim: u32,
    pub pseudo_mismatches: u32,
    pub min_score_gain: Score,
}

impl Default for DoubleDOptions {
    fn default() -> Self {
        Self {
            mode: DoubleDMode::Off,
            min_vj_span: 40,
            seed_length: 11,
            pseudo_trim: 5,
            pseudo_mismatches: 3,
            min_score_gain: 8,
        }
    }
}

/// The sidecar's columns.
pub const COLUMNS: [&str; 38] = [
    "swig_batch_record_index",
    "sequence_id",
    "locus",
    "v_call",
    "standard_d_call",
    "d_call",
    "d2_call",
    "j_call",
    "d_score",
    "d_identity",
    "d_cigar",
    "d2_score",
    "d2_identity",
    "d2_cigar",
    "d_sequence_start",
    "d_sequence_end",
    "d_germline_start",
    "d_germline_end",
    "d2_sequence_start",
    "d2_sequence_end",
    "d2_germline_start",
    "d2_germline_end",
    "d_sequence_alignment",
    "d_germline_alignment",
    "d2_sequence_alignment",
    "d2_germline_alignment",
    "np2",
    "np3",
    "np2_length",
    "np3_length",
    "swig_double_d_vj_span",
    "swig_double_d_seed_length",
    "swig_double_d_pair_score",
    "swig_double_d_best_single_score",
    "swig_double_d_score_gain",
    "swig_double_d_pseudo_distance",
    "swig_double_d_mode",
    "swig_double_d_posterior",
];

/// The sidecar's header line, without a line break.
pub fn header() -> String {
    COLUMNS.join("\t")
}

/// Fewest mismatches between `span` shortened by `trim` bases, at any
/// offset, and any equally long substring of an allele in `alleles`;
/// `None` when the shortened span is shorter than `min_len` or longer than
/// every allele.
fn delta_distance<'a>(
    span: &[u8],
    alleles: impl Iterator<Item = &'a [u8]> + Clone,
    trim: usize,
    min_len: usize,
) -> Option<u32> {
    let len = span
        .len()
        .checked_sub(trim)
        .filter(|&l| l > 0 && l >= min_len)?;
    let mut best: Option<u32> = None;
    for start in 0..=trim {
        let window = &span[start..start + len];
        for allele in alleles.clone() {
            for germline in allele.windows(len) {
                let bound = best.unwrap_or(u32::MAX);
                let mut distance = 0;
                for (a, b) in window.iter().zip(germline) {
                    distance += u32::from(a != b);
                    if distance >= bound {
                        break;
                    }
                }
                best = Some(distance.min(bound));
            }
        }
    }
    best
}

/// The sidecar row of `call` (read `index` of its batch, `sequence_id`),
/// whose main row reports `v`, `d` and `j`, if the tandem reading passes
/// the screen. The line has no line break.
#[allow(clippy::too_many_arguments)]
pub fn screen<C: AlleleCatalog>(
    index: usize,
    sequence_id: &str,
    locus: Option<&str>,
    main: [Option<&SegmentReport>; 3],
    call: &Call,
    catalog: &C,
    support: &SupportModels,
    options: &DoubleDOptions,
) -> Option<String> {
    if options.mode == DoubleDMode::Off {
        return None;
    }
    let tandem = call.tandem.as_ref()?;
    let d_set = catalog.segment(SegmentKind::D)?.set();
    let [v, d, j] = main;
    let (v, j) = (v?, j?);
    // The call's V exit and J entry bound the span the screen reads.
    let (exit, entry) = (v.sequence_end, j.sequence_start - 1);
    let vj_span = entry.saturating_sub(exit);
    if options.mode == DoubleDMode::LongSpan && vj_span < options.min_vj_span {
        return None;
    }
    let ds: Vec<_> = tandem
        .path
        .segments
        .iter()
        .filter(|s| call.placements[s.placement].kind == SegmentKind::D)
        .collect();
    let [d1, d2] = ds[..] else { return None };
    let run = |s: &crate::adjudicate::CalledSegment| {
        let p = &call.placements[s.placement];
        let codes = &d_set.allele(p.primary_allele()).codes;
        Evidence::from_columns(p.clip.columns(&s.span), &call.query, codes).longest_exact_run
    };
    if run(d1) < options.seed_length || run(d2) < options.seed_length {
        return None;
    }
    let pair_score = d1.span.score + d2.span.score;
    let best_single = call
        .placements
        .iter()
        .filter(|p| p.kind == SegmentKind::D)
        .filter_map(|p| p.clip.span(exit, entry.max(exit)))
        .map(|s| s.score)
        .max()
        .unwrap_or(0)
        .max(0);
    let gain = pair_score - best_single;
    if gain < options.min_score_gain {
        return None;
    }
    let loci = locus
        .and_then(|name| Locus::ALL.into_iter().find(|l| l.name() == name))
        .map_or(Loci::ANY, Loci::of);
    let alleles = d_set
        .alleles()
        .iter()
        .filter(move |a| a.loci.pairs_with(loci))
        .map(|a| a.codes.as_slice());
    let span = &call.query[d1.span.q_lo as usize..d2.span.q_hi as usize];
    let pseudo = delta_distance(
        span,
        alleles,
        options.pseudo_trim as usize,
        options.seed_length as usize,
    );
    if pseudo.is_some_and(|p| p <= options.pseudo_mismatches) {
        return None;
    }
    let marginals = tandem.marginals.as_ref();
    let report = |s| segment_report(call, &tandem.path, marginals, s, catalog, support);
    let (r1, r2) = (report(d1)?, report(d2)?);
    let between = |from: u32, to: u32| decode(&call.query[from as usize..to.max(from) as usize]);
    let tandem_entry = tandem.path.segments.last().map_or(entry, |s| s.span.q_lo);
    let np2 = between(d1.span.q_hi, d2.span.q_lo);
    let np3 = between(d2.span.q_hi, tandem_entry);

    let mut row = Row::with_capacity(1024);
    row.field(&index.to_string());
    row.text(sequence_id);
    row.field(locus.unwrap_or(""));
    row.field(&v.call);
    row.field(d.map_or("", |d| &d.call));
    row.field(&r1.call);
    row.field(&r2.call);
    row.field(&j.call);
    for r in [&r1, &r2] {
        row.field(&r.score.to_string());
        row.field(&fraction(r.identity));
        row.field(&r.cigar);
    }
    for r in [&r1, &r2] {
        for value in [
            r.sequence_start,
            r.sequence_end,
            r.germline_start,
            r.germline_end,
        ] {
            row.field(&value.to_string());
        }
    }
    for r in [&r1, &r2] {
        row.field(&r.sequence_alignment);
        row.field(&r.germline_alignment);
    }
    row.field(&np2);
    row.field(&np3);
    row.field(&np2.len().to_string());
    row.field(&np3.len().to_string());
    for value in [
        Score::try_from(vj_span).unwrap_or(Score::MAX),
        Score::try_from(options.seed_length).unwrap_or(Score::MAX),
        pair_score,
        best_single,
        gain,
    ] {
        row.field(&value.to_string());
    }
    row.some(pseudo, |p| p.to_string());
    row.field(options.mode.name());
    row.some(
        marginals.and_then(|m| m.d_count.get(2).copied()),
        probability,
    );
    Some(row.line)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delta_distance_is_the_fewest_mismatches_of_a_trimmed_span() {
        let d1 = b"ACGTACGTACGTAAAA".as_slice();
        let d2 = b"TTTTGGGGCCCC".as_slice();
        // One allele's substring, one base changed, two foreign bases at
        // the end: trimming 2 aligns the rest at one mismatch.
        let span = b"ACGTACCTACGTGG";
        assert_eq!(delta_distance(span, [d1, d2].into_iter(), 2, 8), Some(1));
        assert_eq!(delta_distance(span, [d1, d2].into_iter(), 0, 8), Some(3));
        // Too short once trimmed, or longer than every allele.
        assert_eq!(delta_distance(b"ACGTAC", [d1].into_iter(), 2, 8), None);
        assert_eq!(delta_distance(&[0; 20], [d2].into_iter(), 2, 8), None);
    }

    #[test]
    fn modes_read_their_names() {
        for mode in [DoubleDMode::Off, DoubleDMode::All, DoubleDMode::LongSpan] {
            assert_eq!(DoubleDMode::from_name(mode.name()), Some(mode));
        }
        assert_eq!(
            DoubleDMode::from_name("long-span"),
            Some(DoubleDMode::LongSpan)
        );
        assert_eq!(header().split('\t').count(), COLUMNS.len());
    }
}
