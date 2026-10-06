//! Scoring calls against simulated truth.
//!
//! The metrics are defined so they can be reported in a paper and compared
//! across tools: a call is an [`Assignment`], built from this crate's
//! [`AirrRecord`] or from any AIRR Rearrangement TSV row.
//!
//! * **Exact** — the primary (first) called allele is the simulated allele.
//! * **Equivalent** — the primary allele is *germline-equivalent* on this
//!   read: its sequence contains the trimmed, pre-mutation germline body that
//!   was simulated. No caller can tell such alleles apart from the read, so
//!   this is the fair allele accuracy; *exact* additionally rewards luck.
//! * **Truth listed** — the simulated allele is among the reported names.
//! * **Gene** — the primary allele's gene is the gene of an equivalent allele.
//! * **Boundaries** — absolute error of the called segment ends against the
//!   labelled read positions (only for correctly oriented reads).
//! * **D presence**, stratified by the length of the simulated D remnant:
//!   a 2-nt remnant is not identifiable, a 20-nt one must be found.
//! * **Junction** — the called junction sequence equals the true one
//!   (Cys104 … Trp/Phe118 inclusive), when both anchor codons are templated
//!   in the read.

use std::fmt;
use std::ops::Range;

use super::truth::{Slot, Truth};
use crate::germline::GermlineSet;
use crate::report::{AirrRecord, SegmentReport};

fn gene(name: &str) -> &str {
    GermlineSet::gene_of_name(name)
}

/// One called segment: reported allele names (primary first) and its read
/// interval, 0-based half-open on the oriented read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SegmentCall {
    pub names: Vec<String>,
    pub span: Option<Range<usize>>,
}

/// A caller's answer for one read, from any tool.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Assignment {
    pub rev_comp: bool,
    /// Indexed by [`Slot::index`].
    pub segments: [Option<SegmentCall>; 4],
    pub junction: Option<String>,
}

impl Assignment {
    pub fn from_record(record: &AirrRecord) -> Self {
        let segment = |report: &Option<SegmentReport>| {
            report.as_ref().map(|r| SegmentCall {
                names: r.call.split(',').map(str::to_string).collect(),
                span: Some(r.sequence_start.saturating_sub(1) as usize..r.sequence_end as usize),
            })
        };
        Self {
            rev_comp: record.rev_comp,
            segments: [
                segment(&record.v),
                segment(&record.d),
                segment(&record.d2),
                segment(&record.j),
            ],
            junction: record.junction_sequence().map(str::to_string),
        }
    }

    /// Parse one AIRR Rearrangement TSV row; returns `(sequence_id, assignment)`.
    /// Only `sequence_id` is required. `{v,d,d2,j}_call` may list several
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
        let mut assignment = Self {
            rev_comp: matches!(get("rev_comp"), Some("T" | "TRUE" | "true" | "True")),
            junction: get("junction").map(|j| j.to_ascii_uppercase()),
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
    /// Boundary errors over called, correctly oriented reads.
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

    fn add(
        &mut self,
        truth: Option<&super::truth::TrueSegment>,
        call: Option<&SegmentCall>,
        oriented: bool,
    ) {
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
        if !oriented {
            return;
        }
        self.exact += usize::from(*primary == truth.allele);
        self.equivalent += usize::from(truth.equivalent.contains(primary));
        self.truth_listed += usize::from(call.names.contains(&truth.allele));
        self.gene += usize::from(truth.equivalent.iter().any(|e| gene(e) == gene(primary)));
        if let Some(span) = &call.span {
            self.boundary_reads += 1;
            self.start_error += span.start.abs_diff(truth.span.start);
            self.end_error += span.end.abs_diff(truth.span.end);
            self.start_exact += usize::from(span.start == truth.span.start);
            self.end_exact += usize::from(span.end == truth.span.end);
        }
    }
}

/// D remnant length strata.
pub const D_LENGTH_BINS: [(usize, usize); 6] =
    [(0, 0), (1, 3), (4, 7), (8, 11), (12, 15), (16, usize::MAX)];

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DPresence {
    pub reads: usize,
    /// Reads with any D call.
    pub called: usize,
    /// Reads whose primary D call is germline-equivalent to the first true D.
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
}

impl Evaluation {
    pub fn add(&mut self, truth: &Truth, call: &Assignment) {
        self.reads += 1;
        let oriented = call.rev_comp == truth.reverse_complemented;
        self.oriented += usize::from(oriented);
        self.unassigned += usize::from(!call.assigned());
        for (slot, stats) in [
            (Slot::V, &mut self.v),
            (Slot::D, &mut self.d),
            (Slot::D2, &mut self.d2),
            (Slot::J, &mut self.j),
        ] {
            stats.add(truth.segment(slot), call.segment(slot), oriented);
        }

        let bin = D_LENGTH_BINS
            .iter()
            .position(|&(lo, hi)| (lo..=hi).contains(&truth.d_length))
            .expect("bins cover every length");
        let presence = &mut self.d_presence[bin];
        presence.reads += 1;
        let called_d = call.segment(Slot::D);
        presence.called += usize::from(called_d.is_some());
        if let (Some(called), Some(true_d)) = (called_d, truth.segment(Slot::D)) {
            let primary = called.names.first();
            presence.correct +=
                usize::from(oriented && primary.is_some_and(|p| true_d.equivalent.contains(p)));
        }

        if let Some(junction) = &truth.junction {
            self.junction_reads += 1;
            self.junction_exact +=
                usize::from(oriented && call.junction.as_deref() == Some(junction.as_str()));
        }
    }

    /// Merge another evaluation (e.g. from another thread or stratum).
    pub fn merge(&mut self, other: &Evaluation) {
        fn segment(a: &mut SegmentStats, b: &SegmentStats) {
            a.reads += b.reads;
            a.called += b.called;
            a.exact += b.exact;
            a.equivalent += b.equivalent;
            a.truth_listed += b.truth_listed;
            a.gene += b.gene;
            a.names += b.names;
            a.spurious += b.spurious;
            a.start_error += b.start_error;
            a.end_error += b.end_error;
            a.start_exact += b.start_exact;
            a.end_exact += b.end_exact;
            a.boundary_reads += b.boundary_reads;
        }
        self.reads += other.reads;
        self.oriented += other.oriented;
        self.unassigned += other.unassigned;
        segment(&mut self.v, &other.v);
        segment(&mut self.d, &other.d);
        segment(&mut self.d2, &other.d2);
        segment(&mut self.j, &other.j);
        for (a, b) in self.d_presence.iter_mut().zip(&other.d_presence) {
            a.reads += b.reads;
            a.called += b.called;
            a.correct += b.correct;
        }
        self.junction_reads += other.junction_reads;
        self.junction_exact += other.junction_exact;
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
                if s.called == 0 { f64::NAN } else { s.names as f64 / s.called as f64 },
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
        Ok(())
    }
}
