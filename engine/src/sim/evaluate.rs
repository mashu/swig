//! Scoring calls against simulated truth.
//!
//! The metrics are defined so they can be reported in a paper and compared
//! across tools: a call is an [`Assignment`], built from this crate's
//! [`AirrRecord`] or from any AIRR Rearrangement TSV row.
//!
//! Allele identity, D presence and the junction do not depend on how the
//! caller oriented the read; read positions do.
//!
//! * **Exact** — the primary (first) called allele is the simulated allele.
//! * **Equivalent** — the primary allele is *germline-equivalent* on this
//!   read: its sequence contains the trimmed, pre-mutation germline body that
//!   was simulated. No caller can tell such alleles apart from the read, so
//!   this is the fair allele accuracy; *exact* additionally rewards luck.
//! * **Truth listed** — the simulated allele is among the reported names.
//! * **Gene** — the primary allele's gene is the gene of an equivalent allele.
//! * **Boundaries** — absolute error of the called segment ends against the
//!   labelled read positions (only for reads whose reported orientation,
//!   `rev_comp`, is the true one).
//! * **Tandem D** — on a read with two simulated Ds, a single D call (most
//!   callers report one D) is scored against the true D it covers, and the
//!   other true D is the missed second D.
//! * **D presence**, stratified by the length of the simulated D remnant
//!   (the longer one of a tandem DD): a 2-nt remnant is not identifiable, a
//!   20-nt one must be found.
//! * **Junction** — the called junction sequence equals the true one
//!   (Cys104 … Trp/Phe118 inclusive), when both anchor codons are templated
//!   in the read.
//! * **Calibration** of reported posteriors (`*_posterior`,
//!   `d_presence_posterior`): reads binned by posterior, with the observed
//!   frequency of the event per bin, the expected calibration error and the
//!   Brier score. A segment's event is "some reported name is
//!   germline-equivalent to the truth"; D presence's event is "the read has
//!   a simulated D", however short its remnant.
//! * **Constant region**, when the reads were simulated with C, stratified
//!   by the length of the C part ([`C_LENGTH_BINS`]) and scored against
//!   *anchored* equivalence ([`super::truth`]). Short C remnants are shared
//!   by several genes (`GCCTCCACCAAG` opens IGHG1, IGHG2 and IGHG4), so
//!   allele equivalence alone says little about the isotype; the table adds
//!   the gene of the first name (what a user reading the isotype off
//!   `c_call` gets), whether the reported genes are exactly the equivalent
//!   genes, over-calls (a reported name outside the equivalent set), and the
//!   oracle's resolution (genes in the equivalent set). On reads without C,
//!   a C call is a false call.

use std::collections::BTreeSet;
use std::fmt;
use std::ops::Range;

use super::simulate::{SimulatedRead, Simulator};
use super::truth::{Slot, TrueSegment, Truth};
use crate::germline::GermlineSet;
use crate::pipeline::{AlleleCatalog, Call};
use crate::report::{AirrRecord, SegmentReport};

fn gene(name: &str) -> &str {
    GermlineSet::gene_of_name(name)
}

/// The distinct genes of `names`.
fn genes(names: &[String]) -> BTreeSet<&str> {
    names.iter().map(|n| gene(n)).collect()
}

/// One called segment: reported allele names (primary first) and its read
/// interval, 0-based half-open on the oriented read.
#[derive(Clone, Debug, PartialEq)]
pub struct SegmentCall {
    pub names: Vec<String>,
    pub span: Option<Range<usize>>,
    /// Reported probability that the segment is one of `names`.
    pub posterior: Option<f64>,
}

/// A caller's answer for one read, from any tool.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Assignment {
    /// The reported orientation; `None` when the caller does not report one.
    pub rev_comp: Option<bool>,
    /// Indexed by [`Slot::index`].
    pub segments: [Option<SegmentCall>; 5],
    pub junction: Option<String>,
    /// Reported probability that the read contains a D.
    pub d_presence: Option<f64>,
}

impl Assignment {
    pub fn from_record(record: &AirrRecord) -> Self {
        let segment = |report: &Option<SegmentReport>| {
            report.as_ref().map(|r| SegmentCall {
                names: r.call.split(',').map(str::to_string).collect(),
                span: Some(r.sequence_start.saturating_sub(1) as usize..r.sequence_end as usize),
                posterior: r.posterior,
            })
        };
        Self {
            rev_comp: Some(record.rev_comp),
            segments: [
                segment(&record.v),
                segment(&record.d),
                None,
                segment(&record.j),
                segment(&record.c),
            ],
            junction: record.junction_sequence().map(str::to_string),
            d_presence: record.d_presence_posterior,
        }
    }

    /// Parse one AIRR Rearrangement TSV row; returns `(sequence_id, assignment)`.
    /// Only `sequence_id` is required; without a `rev_comp` value the
    /// orientation is unknown. `{v,d,d2,j,c}_call` may list several
    /// comma-separated names (primary first); coordinates are AIRR 1-based closed.
    pub fn from_airr(header: &[&str], line: &str) -> Result<(String, Self), String> {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| {
            header
                .iter()
                .position(|h| *h == name)
                .and_then(|i| fields.get(i).copied())
                .filter(|v| !v.is_empty())
        };
        let id = get("sequence_id").ok_or("AIRR row without sequence_id")?;
        let probability = |name: &str| get(name).and_then(|v| v.parse::<f64>().ok());
        let mut assignment = Self {
            rev_comp: get("rev_comp").and_then(|v| match v {
                "T" | "TRUE" | "true" | "True" => Some(true),
                "F" | "FALSE" | "false" | "False" => Some(false),
                _ => None,
            }),
            junction: get("junction").map(|j| j.to_ascii_uppercase()),
            d_presence: probability("d_presence_posterior"),
            ..Self::default()
        };
        for slot in Slot::ALL {
            let p = slot.label();
            let Some(call) = get(&format!("{p}_call")) else {
                continue;
            };
            let coordinate = |key: &str| get(&format!("{p}_sequence_{key}"))?.parse::<usize>().ok();
            assignment.segments[slot.index()] = Some(SegmentCall {
                names: call.split(',').map(|n| n.trim().to_string()).collect(),
                span: coordinate("start")
                    .zip(coordinate("end"))
                    .map(|(start, end)| start.saturating_sub(1)..end),
                posterior: probability(&format!("{p}_posterior")),
            });
        }
        Ok((id.to_string(), assignment))
    }

    pub fn segment(&self, slot: Slot) -> Option<&SegmentCall> {
        self.segments[slot.index()].as_ref()
    }

    /// Any V or J call (a read without one is unassigned).
    pub fn assigned(&self) -> bool {
        self.segment(Slot::V).is_some() || self.segment(Slot::J).is_some()
    }

    /// Whether the reported orientation is `truth`'s, so that read
    /// positions are comparable.
    pub fn oriented(&self, truth: &Truth) -> bool {
        self.rev_comp == Some(truth.reverse_complemented)
    }
}

/// The true segment each slot of `call` is scored against, by
/// [`Slot::index`]: the read's own, except on a read with two true Ds and a
/// single D call. That call answers the true D it covers: the one its span
/// overlaps more (when the read positions are comparable), else the one its
/// primary name is equivalent to, the 5′ D first. The other true D is then
/// the missed second D.
fn scored_against<'a>(truth: &'a Truth, call: &Assignment) -> [Option<&'a TrueSegment>; 5] {
    let mut segments = Slot::ALL.map(|slot| truth.segment(slot));
    let (Some(first), Some(second), Some(d), None) = (
        truth.segment(Slot::D),
        truth.segment(Slot::D2),
        call.segment(Slot::D),
        call.segment(Slot::D2),
    ) else {
        return segments;
    };
    let overlap = |t: &TrueSegment| {
        let span = d.span.as_ref().filter(|_| call.oriented(truth))?;
        Some(
            span.end
                .min(t.span.end)
                .saturating_sub(span.start.max(t.span.start)),
        )
    };
    let covers_second = match (overlap(first), overlap(second)) {
        (Some(a), Some(b)) if a != b => b > a,
        _ => !primary_equivalent(first, d) && primary_equivalent(second, d),
    };
    if covers_second {
        segments.swap(Slot::D.index(), Slot::D2.index());
    }
    segments
}

/// The length of the longest simulated D remnant (0 without D).
fn d_remnant(truth: &Truth) -> usize {
    let d2 = truth.segment(Slot::D2).map_or(0, |d| d.span.len());
    truth.d_length.max(d2)
}

/// Whether the primary allele of `call` is germline-equivalent to `truth`.
fn primary_equivalent(truth: &TrueSegment, call: &SegmentCall) -> bool {
    call.names
        .first()
        .is_some_and(|p| truth.equivalent.contains(p))
}

/// Per-segment tallies.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SegmentStats {
    /// Reads whose truth contains this segment.
    pub reads: usize,
    /// ...of which the caller reported the segment.
    pub called: usize,
    pub exact: usize,
    pub equivalent: usize,
    pub truth_listed: usize,
    pub gene: usize,
    /// Reported names summed over called reads (mean ambiguity = names / called).
    pub names: usize,
    /// Calls on reads whose truth lacks the segment.
    pub spurious: usize,
    /// Boundary errors over called reads whose reported orientation is the
    /// true one.
    pub start_error: usize,
    pub end_error: usize,
    pub start_exact: usize,
    pub end_exact: usize,
    pub boundary_reads: usize,
}

impl SegmentStats {
    fn rate(part: usize, whole: usize) -> f64 {
        if whole == 0 {
            f64::NAN
        } else {
            100.0 * part as f64 / whole as f64
        }
    }

    /// Percentage of reads with the segment whose primary call is equivalent.
    pub fn equivalent_rate(&self) -> f64 {
        Self::rate(self.equivalent, self.reads)
    }

    pub fn exact_rate(&self) -> f64 {
        Self::rate(self.exact, self.reads)
    }

    pub fn gene_rate(&self) -> f64 {
        Self::rate(self.gene, self.reads)
    }

    pub fn called_rate(&self) -> f64 {
        Self::rate(self.called, self.reads)
    }

    fn add(&mut self, truth: Option<&TrueSegment>, call: Option<&SegmentCall>, oriented: bool) {
        let Some(truth) = truth else {
            self.spurious += usize::from(call.is_some());
            return;
        };
        self.reads += 1;
        let Some(call) = call else { return };
        self.called += 1;
        self.names += call.names.len();
        let Some(primary) = call.names.first() else {
            return;
        };
        self.exact += usize::from(*primary == truth.allele);
        self.equivalent += usize::from(primary_equivalent(truth, call));
        self.truth_listed += usize::from(call.names.contains(&truth.allele));
        self.gene += usize::from(truth.equivalent.iter().any(|e| gene(e) == gene(primary)));
        if let Some(span) = call.span.as_ref().filter(|_| oriented) {
            self.boundary_reads += 1;
            self.start_error += span.start.abs_diff(truth.span.start);
            self.end_error += span.end.abs_diff(truth.span.end);
            self.start_exact += usize::from(span.start == truth.span.start);
            self.end_exact += usize::from(span.end == truth.span.end);
        }
    }
}

/// C length strata: none, below the gate's reach, then doubling lengths.
/// The gate needs a clipped score of about 18 (`S ≥ 18` for the 10 human
/// IGHC genes), nine exact bases.
pub const C_LENGTH_BINS: [(usize, usize); 7] = [
    (0, 0),
    (1, 7),
    (8, 15),
    (16, 31),
    (32, 63),
    (64, 127),
    (128, usize::MAX),
];

/// Constant-region tallies of one C length stratum, against anchored
/// equivalence. In the stratum without C, `called` counts false calls.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConstantStats {
    pub reads: usize,
    /// Reads with a C call.
    pub called: usize,
    /// The primary name is equivalent.
    pub equivalent: usize,
    /// The primary name's gene is the simulated allele's gene.
    pub first_gene: usize,
    /// The reported genes are exactly the genes of the equivalent alleles.
    pub gene_set: usize,
    /// Calls reporting a name outside the equivalent alleles.
    pub over_call: usize,
    /// Distinct genes reported, summed over calls.
    pub genes_listed: usize,
    /// Distinct genes of the equivalent alleles, summed over reads with C.
    pub genes_equivalent: usize,
}

impl ConstantStats {
    fn add(&mut self, truth: Option<&TrueSegment>, call: Option<&SegmentCall>) {
        self.reads += 1;
        if let Some(truth) = truth {
            self.genes_equivalent += genes(&truth.equivalent).len();
        }
        let Some(call) = call else { return };
        self.called += 1;
        self.genes_listed += genes(&call.names).len();
        let (Some(truth), Some(primary)) = (truth, call.names.first()) else {
            return;
        };
        let listed_equivalent = call.names.iter().all(|n| truth.equivalent.contains(n));
        self.over_call += usize::from(!listed_equivalent);
        self.equivalent += usize::from(primary_equivalent(truth, call));
        self.first_gene += usize::from(gene(primary) == gene(&truth.allele));
        self.gene_set += usize::from(genes(&call.names) == genes(&truth.equivalent));
    }
}

/// D remnant length strata.
pub const D_LENGTH_BINS: [(usize, usize); 6] =
    [(0, 0), (1, 3), (4, 7), (8, 11), (12, 15), (16, usize::MAX)];

/// Upper bounds of the posterior bins of a [`Calibration`] (the last bin
/// is `[0.999, 1]`).
pub const CALIBRATION_BINS: [f64; 7] = [0.5, 0.8, 0.9, 0.95, 0.99, 0.999, f64::INFINITY];

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CalibrationBin {
    pub count: usize,
    /// Sum of the reported posteriors in the bin.
    pub posterior: f64,
    /// Reads whose event occurred.
    pub correct: usize,
}

/// Reported posteriors against observed outcomes.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Calibration {
    pub bins: [CalibrationBin; CALIBRATION_BINS.len()],
    /// Sum of squared errors `(p − outcome)²`.
    pub squared_error: f64,
}

impl Calibration {
    pub fn add(&mut self, posterior: f64, correct: bool) {
        let bin = CALIBRATION_BINS
            .iter()
            .position(|&hi| posterior < hi)
            .unwrap_or(CALIBRATION_BINS.len() - 1);
        let b = &mut self.bins[bin];
        b.count += 1;
        b.posterior += posterior;
        b.correct += usize::from(correct);
        self.squared_error += (posterior - f64::from(u8::from(correct))).powi(2);
    }

    pub fn count(&self) -> usize {
        self.bins.iter().map(|b| b.count).sum()
    }

    /// Expected calibration error: `Σ_bins n_b · |mean posterior − frequency| / n`.
    pub fn ece(&self) -> f64 {
        let n = self.count();
        if n == 0 {
            return f64::NAN;
        }
        self.bins
            .iter()
            .filter(|b| b.count > 0)
            .map(|b| (b.posterior - b.correct as f64).abs())
            .sum::<f64>()
            / n as f64
    }

    /// Mean squared error of the posteriors.
    pub fn brier(&self) -> f64 {
        self.squared_error / self.count() as f64
    }
}

/// Calibration of every reported posterior.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PosteriorCalibration {
    pub v: Calibration,
    pub d: Calibration,
    pub j: Calibration,
    pub d_presence: Calibration,
    pub c: Calibration,
}

impl PosteriorCalibration {
    /// The calibrations to show: C only when some read reported it.
    fn targets(&self) -> Vec<(&'static str, &Calibration)> {
        let mut targets = vec![
            ("V", &self.v),
            ("D", &self.d),
            ("J", &self.j),
            ("D presence", &self.d_presence),
        ];
        if self.c.count() > 0 {
            targets.push(("C", &self.c));
        }
        targets
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DPresence {
    pub reads: usize,
    /// Reads with any D call.
    pub called: usize,
    /// Reads whose primary D call is germline-equivalent to the true D it
    /// is scored against (on a tandem DD, the one it covers).
    pub correct: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Evaluation {
    pub reads: usize,
    pub oriented: usize,
    /// Reads without any V(D)J path.
    pub unassigned: usize,
    pub v: SegmentStats,
    pub d: SegmentStats,
    pub d2: SegmentStats,
    pub j: SegmentStats,
    pub d_presence: [DPresence; D_LENGTH_BINS.len()],
    /// Reads whose junction is templated at both anchors.
    pub junction_reads: usize,
    pub junction_exact: usize,
    /// By [`C_LENGTH_BINS`]; empty unless the reads were simulated with C.
    pub c: [ConstantStats; C_LENGTH_BINS.len()],
    pub calibration: PosteriorCalibration,
}

impl Evaluation {
    /// The evaluation of `calls`, one per read of `reads` (simulated by
    /// `sim`), as `catalog` reports them in its AIRR records.
    pub fn of_calls<C: AlleleCatalog>(
        sim: &Simulator<'_>,
        reads: &[SimulatedRead],
        calls: &[Call],
        catalog: &C,
    ) -> Self {
        let mut evaluation = Self::default();
        for (read, call) in reads.iter().zip(calls) {
            let sequence = String::from_utf8_lossy(&read.sequence);
            let record = AirrRecord::from_call("r", &sequence, call, catalog, &[None, None, None]);
            let truth = Truth::from_read(sim, read, "r");
            evaluation.add(&truth, &Assignment::from_record(&record));
        }
        evaluation
    }

    pub fn add(&mut self, truth: &Truth, call: &Assignment) {
        self.reads += 1;
        let oriented = call.oriented(truth);
        self.oriented += usize::from(oriented);
        self.unassigned += usize::from(!call.assigned());
        let scored = scored_against(truth, call);
        for (slot, stats) in [
            (Slot::V, &mut self.v),
            (Slot::D, &mut self.d),
            (Slot::D2, &mut self.d2),
            (Slot::J, &mut self.j),
        ] {
            stats.add(scored[slot.index()], call.segment(slot), oriented);
        }

        let bin = D_LENGTH_BINS
            .iter()
            .position(|&(lo, hi)| (lo..=hi).contains(&d_remnant(truth)))
            .expect("bins cover every length");
        let presence = &mut self.d_presence[bin];
        presence.reads += 1;
        let called_d = call.segment(Slot::D);
        presence.called += usize::from(called_d.is_some());
        if let (Some(called), Some(true_d)) = (called_d, scored[Slot::D.index()]) {
            presence.correct += usize::from(primary_equivalent(true_d, called));
        }

        if let Some(junction) = &truth.junction {
            self.junction_reads += 1;
            self.junction_exact += usize::from(call.junction.as_deref() == Some(junction.as_str()));
        }

        if truth.constant {
            let length = truth.segment(Slot::C).map_or(0, |c| c.span.len());
            let bin = C_LENGTH_BINS
                .iter()
                .position(|&(lo, hi)| (lo..=hi).contains(&length))
                .expect("bins cover every length");
            self.c[bin].add(truth.segment(Slot::C), call.segment(Slot::C));
        }

        for (slot, calibration) in [
            (Slot::V, &mut self.calibration.v),
            (Slot::D, &mut self.calibration.d),
            (Slot::J, &mut self.calibration.j),
            (Slot::C, &mut self.calibration.c),
        ] {
            let Some(segment) = call.segment(slot) else {
                continue;
            };
            if let Some(posterior) = segment.posterior {
                let correct = scored[slot.index()]
                    .is_some_and(|t| segment.names.iter().any(|n| t.equivalent.contains(n)));
                calibration.add(posterior, correct);
            }
        }
        if let Some(posterior) = call.d_presence {
            self.calibration
                .d_presence
                .add(posterior, truth.segment(Slot::D).is_some());
        }
    }

    pub fn junction_rate(&self) -> f64 {
        SegmentStats::rate(self.junction_exact, self.junction_reads)
    }
}

impl fmt::Display for Evaluation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let pct = SegmentStats::rate;
        writeln!(
            f,
            "reads {}  oriented {:.2}%  unassigned {}  junction exact {:.2}% (of {})",
            self.reads,
            pct(self.oriented, self.reads),
            self.unassigned,
            self.junction_rate(),
            self.junction_reads
        )?;
        writeln!(
            f,
            "seg   reads  called%  equiv%  exact%  listed%   gene%  names  spurious  start|err|  end|err|  start=%  end=%"
        )?;
        for (label, s) in [
            ("V", &self.v),
            ("D", &self.d),
            ("D2", &self.d2),
            ("J", &self.j),
        ] {
            let mean = |x: usize| {
                if s.boundary_reads == 0 {
                    f64::NAN
                } else {
                    x as f64 / s.boundary_reads as f64
                }
            };
            writeln!(
                f,
                "{label:<4}{:>7}  {:>7.2} {:>7.2} {:>7.2} {:>8.2} {:>7.2} {:>6.2} {:>9} {:>11.2} {:>9.2} {:>8.2} {:>6.2}",
                s.reads,
                s.called_rate(),
                s.equivalent_rate(),
                s.exact_rate(),
                pct(s.truth_listed, s.reads),
                s.gene_rate(),
                if s.called == 0 {
                    f64::NAN
                } else {
                    s.names as f64 / s.called as f64
                },
                s.spurious,
                mean(s.start_error),
                mean(s.end_error),
                pct(s.start_exact, s.boundary_reads),
                pct(s.end_exact, s.boundary_reads),
            )?;
        }
        write!(f, "D remnant  reads  D called%  D correct%")?;
        for (&(lo, hi), p) in D_LENGTH_BINS.iter().zip(&self.d_presence) {
            let label = match (lo, hi) {
                (0, 0) => "none".to_string(),
                (lo, usize::MAX) => format!("{lo}+"),
                (lo, hi) => format!("{lo}-{hi}"),
            };
            write!(
                f,
                "\n  {label:<8}{:>6} {:>10.2} {:>11.2}",
                p.reads,
                pct(p.called, p.reads),
                pct(p.correct, p.reads)
            )?;
        }
        if self.c.iter().any(|c| c.reads > 0) {
            write!(
                f,
                "\nC length   reads  called%  equiv%  1st gene%  gene set%  over-call%  genes  oracle"
            )?;
            for (&(lo, hi), c) in C_LENGTH_BINS.iter().zip(&self.c) {
                let label = match (lo, hi) {
                    (0, 0) => "none".to_string(),
                    (lo, usize::MAX) => format!("{lo}+"),
                    (lo, hi) => format!("{lo}-{hi}"),
                };
                write!(
                    f,
                    "\n  {label:<8}{:>6} {:>8.2}",
                    c.reads,
                    pct(c.called, c.reads)
                )?;
                if lo == 0 {
                    continue;
                }
                let mean = |sum: usize, n: usize| sum as f64 / n as f64;
                write!(
                    f,
                    " {:>7.2} {:>10.2} {:>10.2} {:>11.2} {:>6.2} {:>7.2}",
                    pct(c.equivalent, c.reads),
                    pct(c.first_gene, c.reads),
                    pct(c.gene_set, c.reads),
                    pct(c.over_call, c.called),
                    mean(c.genes_listed, c.called),
                    mean(c.genes_equivalent, c.reads),
                )?;
            }
        }
        let targets = self.calibration.targets();
        if targets.iter().all(|(_, c)| c.count() == 0) {
            return Ok(());
        }
        write!(f, "\nposterior   ")?;
        for (label, _) in &targets {
            write!(f, " {label:>22}")?;
        }
        write!(f, "\n  bin       ")?;
        for _ in &targets {
            write!(f, " {:>8} {:>6} {:>6}", "n", "mean%", "obs%")?;
        }
        let mut lo = 0.0;
        for (i, &hi) in CALIBRATION_BINS.iter().enumerate() {
            let label = if hi.is_finite() {
                format!("{lo:.3}-{hi:.3}")
            } else {
                format!("{lo:.3}-1    ")
            };
            write!(f, "\n  {label:<10}")?;
            for (_, c) in &targets {
                let b = &c.bins[i];
                let mean = if b.count == 0 {
                    f64::NAN
                } else {
                    100.0 * b.posterior / b.count as f64
                };
                write!(
                    f,
                    " {:>8} {:>6.2} {:>6.2}",
                    b.count,
                    mean,
                    pct(b.correct, b.count)
                )?;
            }
            lo = hi;
        }
        write!(f, "\n  ECE %     ")?;
        for (_, c) in &targets {
            write!(f, " {:>22.3}", 100.0 * c.ece())?;
        }
        write!(f, "\n  Brier     ")?;
        for (_, c) in &targets {
            write!(f, " {:>22.5}", c.brier())?;
        }
        Ok(())
    }
}

/// One yes/no metric of two callers on the same reads.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PairedOutcome {
    pub both: usize,
    /// Reads only the first caller gets right.
    pub first_only: usize,
    /// Reads only the second caller gets right.
    pub second_only: usize,
    pub neither: usize,
}

impl PairedOutcome {
    pub fn add(&mut self, first: bool, second: bool) {
        *match (first, second) {
            (true, true) => &mut self.both,
            (true, false) => &mut self.first_only,
            (false, true) => &mut self.second_only,
            (false, false) => &mut self.neither,
        } += 1;
    }

    pub fn reads(&self) -> usize {
        self.both + self.first_only + self.second_only + self.neither
    }

    /// Accuracy of the first and the second caller, in percent.
    pub fn rates(&self) -> (f64, f64) {
        let n = self.reads();
        (
            SegmentStats::rate(self.both + self.first_only, n),
            SegmentStats::rate(self.both + self.second_only, n),
        )
    }

    /// Accuracy difference, second − first, in percentage points.
    pub fn difference(&self) -> f64 {
        let (first, second) = self.rates();
        second - first
    }

    /// Half-width of the 95 % confidence interval of [`Self::difference`]
    /// (Wald interval for paired proportions), in percentage points. Only
    /// the discordant reads carry information about the difference.
    pub fn half_width(&self) -> f64 {
        let n = self.reads() as f64;
        let (b, c) = (self.first_only as f64, self.second_only as f64);
        let d = (c - b) / n;
        100.0 * 1.96 * (((b + c) / n - d * d).max(0.0) / n).sqrt()
    }

    /// Two-sided exact McNemar test: the probability of discordant counts at
    /// least this unequal if both callers were equally accurate.
    pub fn p_value(&self) -> f64 {
        let n = self.first_only + self.second_only;
        let k = self.first_only.min(self.second_only);
        if n == 0 {
            return 1.0;
        }
        // P(X ≤ k) for X ~ Binomial(n, 1/2), summed in log space.
        let ln_half = -(n as f64) * std::f64::consts::LN_2;
        let (mut ln_choose, mut tail) = (0.0f64, f64::NEG_INFINITY);
        for i in 0..=k {
            if i > 0 {
                ln_choose += ((n - i + 1) as f64).ln() - (i as f64).ln();
            }
            let term = ln_choose + ln_half;
            let (hi, lo) = (tail.max(term), tail.min(term));
            tail = hi + (lo - hi).exp().ln_1p();
        }
        (2.0 * tail.exp()).min(1.0)
    }
}

/// The Brier score of two callers' posteriors on the same reads, compared
/// read by read: the mean of `d = (p₂ − e₂)² − (p₁ − e₁)²`, with a paired
/// t statistic `t = d̄ / (s_d / √n)`. Each caller's event `e` is its own
/// (its reported names may differ).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PairedBrier {
    pub reads: usize,
    /// `Σ (p − e)²` of the first and the second caller.
    pub first: f64,
    pub second: f64,
    /// `Σ d` and `Σ d²`.
    pub difference: f64,
    pub squared: f64,
}

impl PairedBrier {
    /// One read: each caller's posterior and whether its event occurred.
    pub fn add(&mut self, first: (f64, bool), second: (f64, bool)) {
        let error = |(p, event): (f64, bool)| (p - f64::from(u8::from(event))).powi(2);
        let (a, b) = (error(first), error(second));
        self.reads += 1;
        self.first += a;
        self.second += b;
        self.difference += b - a;
        self.squared += (b - a).powi(2);
    }

    /// Brier score of the first and the second caller.
    pub fn scores(&self) -> (f64, f64) {
        let n = self.reads.max(1) as f64;
        (self.first / n, self.second / n)
    }

    /// Mean difference `d̄`, second − first (negative: the second is better).
    pub fn mean_difference(&self) -> f64 {
        self.difference / self.reads.max(1) as f64
    }

    /// Standard error of [`Self::mean_difference`].
    fn standard_error(&self) -> f64 {
        let n = self.reads as f64;
        if self.reads < 2 {
            return f64::NAN;
        }
        let mean = self.difference / n;
        ((self.squared - n * mean * mean).max(0.0) / (n - 1.0) / n).sqrt()
    }

    /// Half-width of the 95 % confidence interval of the mean difference.
    pub fn half_width(&self) -> f64 {
        1.96 * self.standard_error()
    }

    /// Paired t statistic; `|t| > 2.58` is a difference at p < 0.01.
    pub fn t(&self) -> f64 {
        self.mean_difference() / self.standard_error()
    }
}

/// D remnant strata of the paired D allele table: below `min_d_match` (4 nt
/// by default) a remnant cannot be found, so a call counted correct there
/// rests on chance matches; 4–7 nt is short but detectable.
pub const D_REMNANT_STRATA: [(usize, usize); 3] = [(1, 3), (4, 7), (8, usize::MAX)];

/// Two callers compared read by read on the same simulated reads: allele
/// equivalence of V, D and J (reads with that segment), D by remnant length
/// ([`D_REMNANT_STRATA`]), no D call on reads without one (D allele accuracy
/// alone rewards calling D everywhere), junction exactness (reads with a
/// templated junction), and the Brier scores of the D allele and D presence
/// posteriors (reads where both callers report one). Reads simulated with C
/// add C allele equivalence (reads with C) and no C call (reads without).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Paired {
    pub v: PairedOutcome,
    pub d: PairedOutcome,
    pub d_by_remnant: [PairedOutcome; D_REMNANT_STRATA.len()],
    pub no_d: PairedOutcome,
    pub j: PairedOutcome,
    pub junction: PairedOutcome,
    pub c: PairedOutcome,
    pub no_c: PairedOutcome,
    pub d_brier: PairedBrier,
    pub d_presence_brier: PairedBrier,
}

impl Paired {
    pub fn add(&mut self, truth: &Truth, first: &Assignment, second: &Assignment) {
        // Each caller's D answers the true D it covers ([`scored_against`]).
        let callers = [first, second];
        let scored = callers.map(|call| scored_against(truth, call));
        let correct = |caller: usize, slot: Slot| {
            let call = callers[caller].segment(slot);
            scored[caller][slot.index()]
                .zip(call)
                .is_some_and(|(t, c)| primary_equivalent(t, c))
        };
        for (slot, outcome) in [
            (Slot::V, &mut self.v),
            (Slot::D, &mut self.d),
            (Slot::J, &mut self.j),
            (Slot::C, &mut self.c),
        ] {
            if truth.segment(slot).is_some() {
                outcome.add(correct(0, slot), correct(1, slot));
            }
        }
        if truth.segment(Slot::D).is_some()
            && let Some(stratum) = D_REMNANT_STRATA
                .iter()
                .position(|&(lo, hi)| (lo..=hi).contains(&d_remnant(truth)))
        {
            self.d_by_remnant[stratum].add(correct(0, Slot::D), correct(1, Slot::D));
        }
        let silent = |call: &Assignment, slot: Slot| call.segment(slot).is_none();
        if truth.segment(Slot::D).is_none() {
            self.no_d
                .add(silent(first, Slot::D), silent(second, Slot::D));
        }
        if truth.constant && truth.segment(Slot::C).is_none() {
            self.no_c
                .add(silent(first, Slot::C), silent(second, Slot::C));
        }
        // The calibration events of [`Evaluation`]: some reported D name is
        // equivalent to the truth; the read has a D.
        let d_posterior = |caller: usize| {
            let segment = callers[caller].segment(Slot::D)?;
            let listed = scored[caller][Slot::D.index()]
                .is_some_and(|t| segment.names.iter().any(|n| t.equivalent.contains(n)));
            Some((segment.posterior?, listed))
        };
        if let (Some(a), Some(b)) = (d_posterior(0), d_posterior(1)) {
            self.d_brier.add(a, b);
        }
        let has_d = truth.segment(Slot::D).is_some();
        if let (Some(a), Some(b)) = (first.d_presence, second.d_presence) {
            self.d_presence_brier.add((a, has_d), (b, has_d));
        }
        if let Some(junction) = &truth.junction {
            let exact = |call: &Assignment| call.junction.as_deref() == Some(junction.as_str());
            self.junction.add(exact(first), exact(second));
        }
    }
}

impl fmt::Display for Paired {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "metric       reads   first%  second%   Δ pp   ±95% CI   first only  second only  McNemar p"
        )?;
        let [short, detectable, long] = &self.d_by_remnant;
        let mut rows = vec![
            ("V allele", &self.v),
            ("D allele", &self.d),
            ("D 1-3 nt", short),
            ("D 4-7 nt", detectable),
            ("D 8+ nt", long),
            ("no D call", &self.no_d),
            ("J allele", &self.j),
            ("junction", &self.junction),
        ];
        if self.c.reads() + self.no_c.reads() > 0 {
            rows.extend([("C allele", &self.c), ("no C call", &self.no_c)]);
        }
        for (label, o) in rows {
            let (first, second) = o.rates();
            write!(
                f,
                "\n{label:<10}{:>8} {:>8.2} {:>8.2} {:>+7.2} {:>9.2} {:>12} {:>12} {:>10.2e}",
                o.reads(),
                first,
                second,
                o.difference(),
                o.half_width(),
                o.first_only,
                o.second_only,
                o.p_value()
            )?;
        }
        write!(
            f,
            "\nBrier             reads    first   second         Δ   ±95% CI   paired t"
        )?;
        for (label, b) in [
            ("D allele", &self.d_brier),
            ("D presence", &self.d_presence_brier),
        ] {
            let (first, second) = b.scores();
            write!(
                f,
                "\n{label:<14}{:>9} {:>8.5} {:>8.5} {:>+9.5} {:>9.5} {:>10.2}",
                b.reads,
                first,
                second,
                b.mean_difference(),
                b.half_width(),
                b.t()
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn outcome(first_only: usize, second_only: usize, both: usize) -> PairedOutcome {
        PairedOutcome {
            both,
            first_only,
            second_only,
            neither: 0,
        }
    }

    #[test]
    fn mcnemar_is_the_exact_binomial_tail() {
        // 2 vs 8 discordant: P(X ≤ 2) = (1 + 10 + 45) / 1024, doubled.
        let p = outcome(2, 8, 90).p_value();
        assert!((p - 2.0 * 56.0 / 1024.0).abs() < 1e-12, "{p}");
        assert_eq!(outcome(0, 0, 100).p_value(), 1.0);
        assert_eq!(outcome(5, 5, 0).p_value(), 1.0);
        // Large and lopsided stays finite and tiny.
        let p = outcome(10, 400, 5000).p_value();
        assert!(p > 0.0 && p < 1e-90, "{p}");
    }

    #[test]
    fn paired_brier_is_a_paired_t_on_squared_errors() {
        let mut b = PairedBrier::default();
        // Errors (first, second): (0.04, 0.01), (0.25, 0.09), (0.01, 0.01).
        b.add((0.8, true), (0.9, true));
        b.add((0.5, false), (0.3, false));
        b.add((0.1, false), (0.9, true));
        let (first, second) = b.scores();
        assert!((first - 0.10).abs() < 1e-12 && (second - 0.11 / 3.0).abs() < 1e-12);
        // d = −0.03, −0.16, 0: mean −0.0633, sd 0.0850, se 0.0491.
        assert!((b.mean_difference() + 0.19 / 3.0).abs() < 1e-12);
        assert!((b.t() + 1.2898).abs() < 1e-3, "{}", b.t());
        assert!(PairedBrier::default().t().is_nan());
    }

    #[test]
    fn paired_difference_and_interval() {
        let o = outcome(10, 30, 960);
        assert_eq!(o.rates(), (97.0, 99.0));
        assert!((o.difference() - 2.0).abs() < 1e-12);
        // sqrt((40/1000 − 0.02²) / 1000) · 1.96 = 1.233 pp.
        assert!((o.half_width() - 1.2334).abs() < 1e-3, "{}", o.half_width());
    }

    fn segment(allele: &str, span: Range<usize>) -> TrueSegment {
        TrueSegment {
            allele: allele.into(),
            equivalent: vec![allele.into()],
            span,
        }
    }

    /// A read with V, J and a templated junction, emitted reverse
    /// complemented.
    fn reversed_truth() -> Truth {
        Truth {
            id: "r".into(),
            reverse_complemented: true,
            segments: [
                Some(segment("V1*01", 0..9)),
                None,
                None,
                Some(segment("J1*01", 40..60)),
                None,
            ],
            constant: false,
            junction: Some("TGTGCG".into()),
            d_length: 0,
            shm_rate: 0.0,
        }
    }

    /// Allele identity and the junction do not depend on the reported
    /// orientation; boundaries need it, and a calls file without `rev_comp`
    /// reports none: a perfect caller without the column is perfect.
    #[test]
    fn orientation_gates_boundaries_only() {
        let truth = reversed_truth();
        let header = [
            "sequence_id",
            "v_call",
            "v_sequence_start",
            "v_sequence_end",
            "j_call",
            "junction",
        ];
        let row = "r\tV1*01\t1\t9\tJ1*01\ttgtgcg";
        let (_, unknown) = Assignment::from_airr(&header, row).unwrap();
        assert_eq!(unknown.rev_comp, None);
        let with_column = [&header[..], &["rev_comp"]].concat();
        let (_, reported) = Assignment::from_airr(&with_column, &format!("{row}\tT")).unwrap();
        assert_eq!(reported.rev_comp, Some(true));
        let (_, wrong) = Assignment::from_airr(&with_column, &format!("{row}\tF")).unwrap();
        assert_eq!(wrong.rev_comp, Some(false));
        for (call, oriented) in [(unknown, false), (wrong, false), (reported, true)] {
            let mut e = Evaluation::default();
            e.add(&truth, &call);
            assert_eq!(e.oriented, usize::from(oriented));
            assert_eq!(
                (e.v.equivalent, e.j.equivalent, e.junction_exact),
                (1, 1, 1)
            );
            assert_eq!(e.v.boundary_reads, usize::from(oriented));
            let mut paired = Paired::default();
            paired.add(&truth, &call, &call);
            assert_eq!((paired.v.both, paired.junction.both), (1, 1));
        }
    }

    /// A read with V, a 10-nt 5′ D, a 12-nt 3′ D and J, in germline sense.
    fn tandem_truth() -> Truth {
        Truth {
            id: "r".into(),
            reverse_complemented: false,
            segments: [
                Some(segment("V1*01", 0..9)),
                Some(segment("D1*01", 10..20)),
                Some(segment("D2*01", 25..37)),
                Some(segment("J1*01", 40..60)),
                None,
            ],
            constant: false,
            junction: None,
            d_length: 10,
            shm_rate: 0.0,
        }
    }

    /// A forward call of `calls` (slot, primary name, span), posterior 0.9.
    fn call(calls: &[(Slot, &str, Option<Range<usize>>)]) -> Assignment {
        let mut call = Assignment {
            rev_comp: Some(false),
            ..Assignment::default()
        };
        for (slot, name, span) in calls {
            call.segments[slot.index()] = Some(SegmentCall {
                names: vec![name.to_string()],
                span: span.clone(),
                posterior: Some(0.9),
            });
        }
        call
    }

    /// A single D call on a tandem-DD read answers the true D it covers,
    /// whichever of the two that is; the other is the missed second D.
    #[test]
    fn a_single_d_call_is_scored_against_the_true_d_it_covers() {
        let truth = tandem_truth();
        let d = |name, span: Option<Range<usize>>| call(&[(Slot::D, name, span)]);
        let evaluate = |call: &Assignment| {
            let mut e = Evaluation::default();
            e.add(&truth, call);
            e
        };
        // The 3′ D, by its span and without one (by its name).
        for call in [d("D2*01", Some(25..37)), d("D2*01", None)] {
            let e = evaluate(&call);
            assert_eq!((e.d.reads, e.d.equivalent, e.d.exact), (1, 1, 1));
            assert_eq!((e.d2.reads, e.d2.called), (1, 0));
            // Binned by the longer remnant (12 nt).
            let presence = DPresence {
                reads: 1,
                called: 1,
                correct: 1,
            };
            assert_eq!(e.d_presence[4], presence);
            assert_eq!(e.calibration.d.bins[3].correct, 1, "0.9 is in [0.9, 0.95)");
        }
        let e = evaluate(&d("D2*01", Some(25..37)));
        let boundaries = (e.d.boundary_reads, e.d.start_error, e.d.end_error);
        assert_eq!(boundaries, (1, 0, 0));
        // The 5′ D is still right.
        let e = evaluate(&d("D1*01", Some(10..20)));
        assert_eq!((e.d.equivalent, e.d.start_error, e.d.end_error), (1, 0, 0));
        // A wrong allele over the 3′ D is wrong, with exact boundaries.
        let e = evaluate(&d("D3*01", Some(25..37)));
        assert_eq!((e.d.equivalent, e.d.start_error, e.d.end_error), (0, 0, 0));
        // Two D calls are scored in read order.
        let both = call(&[
            (Slot::D, "D2*01", Some(25..37)),
            (Slot::D2, "D1*01", Some(10..20)),
        ]);
        let e = evaluate(&both);
        assert_eq!((e.d.equivalent, e.d2.equivalent, e.d2.called), (0, 0, 1));
        // Paired: each caller's D is matched on its own.
        let mut paired = Paired::default();
        paired.add(&truth, &d("D2*01", Some(25..37)), &d("D1*01", Some(10..20)));
        assert_eq!((paired.d.both, paired.d_by_remnant[2].both), (1, 1));
        assert_eq!(paired.d_brier.reads, 1);
        assert!(paired.d_brier.difference.abs() < 1e-12 && paired.d_brier.first < 0.011);
    }

    #[test]
    fn constant_calls_are_scored_by_gene_set_and_first_name() {
        let names = |list: &[&str]| list.iter().map(|n| n.to_string()).collect::<Vec<_>>();
        let truth = TrueSegment {
            allele: "G1*01".into(),
            equivalent: names(&["G1*01", "G1*02", "G2*01"]),
            span: 100..110,
        };
        let call = |list: &[&str]| SegmentCall {
            names: names(list),
            span: None,
            posterior: None,
        };
        let mut stats = ConstantStats::default();
        // Equivalent, but the first name is the other gene.
        stats.add(Some(&truth), Some(&call(&["G2*01", "G1*01"])));
        // Right gene first; a gene missing from the list.
        stats.add(Some(&truth), Some(&call(&["G1*02"])));
        // An over-call: G3 is not equivalent.
        stats.add(Some(&truth), Some(&call(&["G1*01", "G2*01", "G3*01"])));
        // No call is not an over-call.
        stats.add(Some(&truth), None);
        let expected = ConstantStats {
            reads: 4,
            called: 3,
            equivalent: 3,
            first_gene: 2,
            gene_set: 1,
            over_call: 1,
            genes_listed: 2 + 1 + 3,
            genes_equivalent: 4 * 2,
        };
        assert_eq!(stats, expected);
    }
}
