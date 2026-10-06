//! Ground truth of simulated reads, independent of the simulator object.
//!
//! [`Truth`] is everything evaluation needs about one read: the simulated
//! alleles with their germline-equivalent alternatives, the labelled read
//! intervals, the junction, and the D remnant length. It is built from a
//! [`SimulatedRead`] and written to and read from a TSV, so reads simulated
//! once can score any tool's AIRR output (`vdj-assign simulate` / `evaluate`).
//!
//! Coordinates are 0-based half-open on the read *in germline sense*, which
//! is the oriented query of a caller that orients the read correctly.

use std::ops::Range;

use super::simulate::{Region, SimulatedRead, Simulator};
use crate::germline::{GermlineSet, SegmentKind};

/// Segment slots of a call, in read order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Slot {
    V,
    D,
    /// Second D of a tandem DD.
    D2,
    J,
}

impl Slot {
    pub const ALL: [Slot; 4] = [Slot::V, Slot::D, Slot::D2, Slot::J];

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn label(self) -> &'static str {
        match self {
            Slot::V => "v",
            Slot::D => "d",
            Slot::D2 => "d2",
            Slot::J => "j",
        }
    }

    fn region(self) -> Region {
        match self {
            Slot::V => Region::V,
            Slot::D => Region::D,
            Slot::D2 => Region::D2,
            Slot::J => Region::J,
        }
    }

    fn kind(self) -> SegmentKind {
        match self {
            Slot::V => SegmentKind::V,
            Slot::D | Slot::D2 => SegmentKind::D,
            Slot::J => SegmentKind::J,
        }
    }
}

/// One simulated segment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TrueSegment {
    /// The simulated allele.
    pub allele: String,
    /// Every allele name (aliases included) whose sequence contains the
    /// trimmed, pre-mutation germline body: no read can separate them.
    pub equivalent: Vec<String>,
    /// Labelled read interval.
    pub span: Range<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Truth {
    pub id: String,
    /// The emitted read is the reverse complement of the germline-sense read.
    pub reverse_complemented: bool,
    /// Indexed by [`Slot::index`].
    pub segments: [Option<TrueSegment>; 4],
    /// Junction (V Cys104 … J Trp/Phe118, inclusive) when both anchor codons
    /// are templated in the read, as germline-sense nucleotides.
    pub junction: Option<String>,
    /// Length of the first D remnant in the read (0 without D).
    pub d_length: usize,
    pub shm_rate: f64,
}

/// Names of alleles whose sequence contains `body`.
fn equivalent_names(set: &GermlineSet, body: &[u8]) -> Vec<String> {
    let codes = crate::alphabet::encode(body);
    set.alleles()
        .iter()
        .filter(|a| {
            a.codes
                .windows(codes.len().max(1))
                .any(|w| w == codes.as_slice())
        })
        .flat_map(|a| a.names().map(str::to_string))
        .collect()
}

impl Truth {
    pub fn from_read(sim: &Simulator<'_>, read: &SimulatedRead, id: &str) -> Self {
        let allele = |slot: Slot| match slot {
            Slot::V => Some(read.v),
            Slot::D => read.d,
            Slot::D2 => read.d2,
            Slot::J => Some(read.j),
        };
        let segments = Slot::ALL.map(|slot| {
            let set = sim.set(slot.kind())?;
            let index = allele(slot)?;
            let body = sim.germline_body(read, slot.region())?;
            Some(TrueSegment {
                allele: set.alleles()[index].name.clone(),
                equivalent: equivalent_names(set, body),
                span: read.span(slot.region())?,
            })
        });
        let anchor = |kind: SegmentKind, index: usize| sim.set(kind)?.alleles()[index].anchor;
        let junction = (|| {
            let (va, ja) = (
                anchor(SegmentKind::V, read.v)?,
                anchor(SegmentKind::J, read.j)?,
            );
            let start = read.locate(Region::V, va)?;
            let end = read.locate(Region::J, ja + 2)? + 1;
            let intact = read.locate(Region::V, va + 2) == Some(start + 2)
                && read.locate(Region::J, ja) == Some(end - 3);
            intact.then(|| String::from_utf8_lossy(&read.forward[start..end]).into_owned())
        })();
        Self {
            id: id.to_string(),
            reverse_complemented: read.reverse_complemented,
            segments,
            junction,
            d_length: read.span(Region::D).map_or(0, |s| s.len()),
            shm_rate: read.shm_rate,
        }
    }

    pub fn segment(&self, slot: Slot) -> Option<&TrueSegment> {
        self.segments[slot.index()].as_ref()
    }

    /// Column names of [`Truth::to_tsv`].
    pub fn tsv_header() -> String {
        let mut columns = vec!["sequence_id".to_string(), "rev_comp".to_string()];
        for slot in Slot::ALL {
            let p = slot.label();
            columns.extend([
                format!("{p}_call"),
                format!("{p}_equivalent"),
                format!("{p}_start"),
                format!("{p}_end"),
            ]);
        }
        columns.extend(["junction", "d_length", "shm_rate"].map(str::to_string));
        columns.join("\t")
    }

    /// One TSV row; spans are 0-based half-open on the germline-sense read.
    pub fn to_tsv(&self) -> String {
        let mut fields = vec![
            self.id.clone(),
            (if self.reverse_complemented { "T" } else { "F" }).to_string(),
        ];
        for segment in &self.segments {
            match segment {
                Some(s) => fields.extend([
                    s.allele.clone(),
                    s.equivalent.join(","),
                    s.span.start.to_string(),
                    s.span.end.to_string(),
                ]),
                None => fields.extend(std::iter::repeat(String::new()).take(4)),
            }
        }
        fields.push(self.junction.clone().unwrap_or_default());
        fields.push(self.d_length.to_string());
        fields.push(self.shm_rate.to_string());
        fields.join("\t")
    }

    /// Parse a row written by [`Truth::to_tsv`] (columns located by `header`).
    pub fn from_tsv(header: &[&str], line: &str) -> Result<Self, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| -> Result<&str, String> {
            let index = header
                .iter()
                .position(|h| *h == name)
                .ok_or_else(|| format!("truth TSV lacks column '{name}'"))?;
            Ok(fields.get(index).copied().unwrap_or(""))
        };
        let number = |name: &str| -> Result<usize, String> {
            get(name)?
                .parse()
                .map_err(|_| format!("bad {name} '{}'", get(name).unwrap_or("")))
        };
        let mut segments: [Option<TrueSegment>; 4] = Default::default();
        for slot in Slot::ALL {
            let p = slot.label();
            let allele = get(&format!("{p}_call"))?;
            if allele.is_empty() {
                continue;
            }
            segments[slot.index()] = Some(TrueSegment {
                allele: allele.to_string(),
                equivalent: get(&format!("{p}_equivalent"))?
                    .split(',')
                    .filter(|n| !n.is_empty())
                    .map(str::to_string)
                    .collect(),
                span: number(&format!("{p}_start"))?..number(&format!("{p}_end"))?,
            });
        }
        let junction = get("junction")?;
        Ok(Self {
            id: get("sequence_id")?.to_string(),
            reverse_complemented: get("rev_comp")? == "T",
            segments,
            junction: (!junction.is_empty()).then(|| junction.to_string()),
            d_length: number("d_length")?,
            shm_rate: get("shm_rate")?
                .parse()
                .map_err(|_| "bad shm_rate".to_string())?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::{SimParams, SyntheticDatabase, SyntheticSpec};

    #[test]
    fn truth_round_trips_through_tsv() {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 2);
        let params = SimParams::mid().with_reverse_complement(0.5);
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params).unwrap();
        let header = Truth::tsv_header();
        let columns: Vec<&str> = header.split('\t').collect();
        for (i, read) in sim.reads(50, 3).unwrap().iter().enumerate() {
            let truth = Truth::from_read(&sim, read, &format!("r{i}"));
            let v = truth.segment(Slot::V).unwrap();
            assert!(v.equivalent.contains(&v.allele));
            let back = Truth::from_tsv(&columns, &truth.to_tsv()).unwrap();
            assert_eq!(back.segments, truth.segments);
            assert_eq!(
                (back.junction, back.d_length),
                (truth.junction, truth.d_length)
            );
        }
    }
}
