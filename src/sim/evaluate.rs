//! Scoring AIRR calls against simulated truth.
//!
//! The metrics are defined so they can be reported in a paper and compared
//! across tools (anything that yields an [`AirrRecord`]-like call):
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
//! * **Junction** — the called junction interval equals the true one
//!   (Cys104 … Trp/Phe118 inclusive), when both anchor codons are templated
//!   in the read.

use std::fmt;

use super::simulate::{Region, SimulatedRead, Simulator};
use crate::germline::{GermlineSet, SegmentKind};
use crate::report::{AirrRecord, SegmentReport};

/// Names of alleles whose sequence contains `body` (germline-equivalent).
fn equivalent_names<'s>(set: &'s GermlineSet, body: &[u8]) -> Vec<&'s str> {
    let codes = crate::alphabet::encode(body);
    set.alleles()
        .iter()
        .filter(|a| codes.is_empty() || a.codes.windows(codes.len()).any(|w| w == codes.as_slice()))
        .flat_map(|a| a.names())
        .collect()
}

fn gene(name: &str) -> &str {
    GermlineSet::gene_of_name(name)
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
        truth: Option<(&str, Vec<&str>, std::ops::Range<usize>)>,
        call: Option<&SegmentReport>,
        oriented: bool,
    ) {
        let Some((truth_name, equivalent, span)) = truth else {
            self.spurious += usize::from(call.is_some());
            return;
        };
        self.reads += 1;
        let Some(call) = call else { return };
        self.called += 1;
        let names: Vec<&str> = call.call.split(',').collect();
        self.names += names.len();
        let primary = names[0];
        if !oriented {
            return;
        }
        self.exact += usize::from(primary == truth_name);
        self.equivalent += usize::from(equivalent.contains(&primary));
        self.truth_listed += usize::from(names.contains(&truth_name));
        self.gene += usize::from(equivalent.iter().any(|&e| gene(e) == gene(primary)));
        let start = (call.sequence_start as usize).saturating_sub(1);
        let end = call.sequence_end as usize;
        self.boundary_reads += 1;
        self.start_error += start.abs_diff(span.start);
        self.end_error += end.abs_diff(span.end);
        self.start_exact += usize::from(start == span.start);
        self.end_exact += usize::from(end == span.end);
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
    pub fn add(&mut self, sim: &Simulator<'_>, read: &SimulatedRead, record: &AirrRecord) {
        self.reads += 1;
        let oriented = record.rev_comp == read.reverse_complemented;
        self.oriented += usize::from(oriented);
        self.unassigned += usize::from(record.path_score.is_none());

        let truth = |kind: SegmentKind, region: Region, allele: Option<usize>| {
            let set = sim.set(kind)?;
            let allele = allele?;
            let body = sim.germline_body(read, region)?;
            let span = read.span(region)?;
            Some((
                set.alleles()[allele].name.as_str(),
                equivalent_names(set, body),
                span,
            ))
        };
        self.v.add(
            truth(SegmentKind::V, Region::V, Some(read.v)),
            record.v.as_ref(),
            oriented,
        );
        self.d.add(
            truth(SegmentKind::D, Region::D, read.d),
            record.d.as_ref(),
            oriented,
        );
        self.d2.add(
            truth(SegmentKind::D, Region::D2, read.d2),
            record.d2.as_ref(),
            oriented,
        );
        self.j.add(
            truth(SegmentKind::J, Region::J, Some(read.j)),
            record.j.as_ref(),
            oriented,
        );

        let d_length = read.span(Region::D).map_or(0, |s| s.len());
        let bin = D_LENGTH_BINS
            .iter()
            .position(|&(lo, hi)| (lo..=hi).contains(&d_length))
            .expect("bins cover every length");
        let presence = &mut self.d_presence[bin];
        presence.reads += 1;
        presence.called += usize::from(record.d.is_some());
        if let (Some(call), Some((_, equivalent, _))) =
            (&record.d, truth(SegmentKind::D, Region::D, read.d))
        {
            let primary = call.call.split(',').next().unwrap_or_default();
            presence.correct += usize::from(oriented && equivalent.contains(&primary));
        }

        let v_anchor = sim
            .set(SegmentKind::V)
            .and_then(|s| s.alleles()[read.v].anchor);
        let j_anchor = sim
            .set(SegmentKind::J)
            .and_then(|s| s.alleles()[read.j].anchor);
        if let (Some(va), Some(ja)) = (v_anchor, j_anchor) {
            let start = read.locate(Region::V, va);
            let end = read.locate(Region::J, ja + 2);
            let intact = start.is_some_and(|s| read.locate(Region::V, va + 2) == Some(s + 2))
                && end.is_some_and(|e| read.locate(Region::J, ja) == Some(e - 2));
            if let (true, Some(start), Some(end)) = (intact, start, end) {
                self.junction_reads += 1;
                self.junction_exact += usize::from(
                    oriented && record.regions.junction == Some(start as u32..end as u32 + 1),
                );
            }
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
