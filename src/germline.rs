//! Germline allele databases.

use std::collections::HashMap;
use std::fmt;
use std::io::{self, Write};

use crate::alphabet;
use crate::fasta::{FastaError, parse_fasta};

/// Gene segment type. The derive order (V < D < J) is genomic order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SegmentKind {
    V,
    D,
    J,
}

impl SegmentKind {
    pub const ALL: [SegmentKind; 3] = [SegmentKind::V, SegmentKind::D, SegmentKind::J];

    /// Dense index for per-kind arrays.
    #[inline]
    pub const fn index(self) -> usize {
        match self {
            SegmentKind::V => 0,
            SegmentKind::D => 1,
            SegmentKind::J => 2,
        }
    }

    pub const fn label(self) -> &'static str {
        match self {
            SegmentKind::V => "V",
            SegmentKind::D => "D",
            SegmentKind::J => "J",
        }
    }
}

/// Dense allele id within one [`GermlineSet`].
pub type AlleleId = u32;
/// Dense gene id within one [`GermlineSet`].
pub type GeneId = u32;

/// 0-based germline start of each V region (IMGT definitions: FR1, CDR1,
/// FR2, CDR2, FR3). FR3 ends after the Cys104 codon, i.e. at `anchor + 3`;
/// CDR3 starts there. `None` where unknown.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RegionStarts {
    pub fr1: Option<u32>,
    pub cdr1: Option<u32>,
    pub fr2: Option<u32>,
    pub cdr2: Option<u32>,
    pub fr3: Option<u32>,
}

impl RegionStarts {
    /// Header keys, in germline order.
    pub const KEYS: [&'static str; 5] = ["fr1", "cdr1", "fr2", "cdr2", "fr3"];

    pub fn values(&self) -> [Option<u32>; 5] {
        [self.fr1, self.cdr1, self.fr2, self.cdr2, self.fr3]
    }

    pub fn from_values(values: [Option<u32>; 5]) -> Self {
        let [fr1, cdr1, fr2, cdr2, fr3] = values;
        Self {
            fr1,
            cdr1,
            fr2,
            cdr2,
            fr3,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.values().iter().all(Option::is_none)
    }

    /// Known starts must be non-decreasing, end no later than `anchor` and
    /// lie inside the sequence.
    pub fn check(&self, anchor: Option<u32>, len: usize) -> Result<(), String> {
        let mut previous: Option<(&str, u32)> = None;
        for (key, value) in Self::KEYS.iter().zip(self.values()) {
            let Some(value) = value else { continue };
            if value as usize >= len {
                return Err(format!("{key}={value} is outside the {len}-nt sequence"));
            }
            if let Some((before, at)) = previous.filter(|&(_, at)| value < at) {
                return Err(format!("{key}={value} precedes {before}={at}"));
            }
            previous = Some((key, value));
        }
        match (previous, anchor) {
            (Some((key, at)), Some(anchor)) if at > anchor => {
                Err(format!("{key}={at} follows anchor={anchor}"))
            }
            _ => Ok(()),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Allele {
    /// Full label, e.g. `IGHV1-2*02`.
    pub name: String,
    pub gene: GeneId,
    /// Encoded germline sequence (see [`crate::alphabet`]).
    pub codes: Vec<u8>,
    /// Conserved-anchor position in germline coordinates: 0-based first base
    /// of the V Cys104 codon or the J Trp/Phe118 codon. Fixes the reading frame.
    pub anchor: Option<u32>,
    /// V region boundaries (unused for D and J; J's FR4 starts at `anchor`).
    pub regions: RegionStarts,
    /// Other names of an identical sequence merged into this allele.
    pub aliases: Vec<String>,
}

impl Allele {
    /// The name followed by its aliases.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        std::iter::once(self.name.as_str()).chain(self.aliases.iter().map(String::as_str))
    }
}

/// An identical sequence merged into an earlier allele.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MergedDuplicate {
    pub kept: String,
    pub merged: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GermlineError {
    Fasta(FastaError),
    EmptySequence {
        name: String,
    },
    BadAttribute {
        name: String,
        key: String,
        value: String,
    },
    BadRegions {
        name: String,
        reason: String,
    },
    ConflictingDuplicate {
        name: String,
    },
}

impl fmt::Display for GermlineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fasta(e) => write!(f, "{e}"),
            Self::EmptySequence { name } => write!(f, "allele '{name}' has an empty sequence"),
            Self::BadAttribute { name, key, value } => {
                write!(f, "allele '{name}' has an invalid {key}='{value}'")
            }
            Self::BadRegions { name, reason } => write!(f, "allele '{name}': {reason}"),
            Self::ConflictingDuplicate { name } => {
                write!(
                    f,
                    "allele name '{name}' occurs twice with different sequences"
                )
            }
        }
    }
}

impl std::error::Error for GermlineError {}

impl From<FastaError> for GermlineError {
    fn from(e: FastaError) -> Self {
        Self::Fasta(e)
    }
}

/// All alleles of one segment kind for one locus.
#[derive(Clone, Debug)]
pub struct GermlineSet {
    kind: SegmentKind,
    alleles: Vec<Allele>,
    genes: Vec<String>,
    gene_ids: HashMap<String, GeneId>,
}

fn parse_position(name: &str, key: &str, value: &str) -> Result<u32, GermlineError> {
    value
        .parse::<u32>()
        .map_err(|_| GermlineError::BadAttribute {
            name: name.to_string(),
            key: key.to_string(),
            value: value.to_string(),
        })
}

impl GermlineSet {
    pub fn new(kind: SegmentKind) -> Self {
        Self {
            kind,
            alleles: Vec::new(),
            genes: Vec::new(),
            gene_ids: HashMap::new(),
        }
    }

    /// Gene part of an allele label (`IGHV1-2*02` → `IGHV1-2`).
    pub fn gene_of_name(name: &str) -> &str {
        name.split('*').next().unwrap_or(name)
    }

    /// Add an allele from ASCII sequence.
    pub fn push(
        &mut self,
        name: &str,
        sequence: &[u8],
        anchor: Option<u32>,
    ) -> Result<AlleleId, GermlineError> {
        if sequence.is_empty() {
            return Err(GermlineError::EmptySequence {
                name: name.to_string(),
            });
        }
        self.push_codes(name, alphabet::encode(sequence), anchor)
    }

    fn push_codes(
        &mut self,
        name: &str,
        codes: Vec<u8>,
        anchor: Option<u32>,
    ) -> Result<AlleleId, GermlineError> {
        let gene_name = Self::gene_of_name(name);
        let next_gene = GeneId::try_from(self.genes.len()).expect("gene count fits u32");
        let gene = *self
            .gene_ids
            .entry(gene_name.to_string())
            .or_insert_with(|| {
                self.genes.push(gene_name.to_string());
                next_gene
            });
        let id = AlleleId::try_from(self.alleles.len()).expect("allele count fits u32");
        self.alleles.push(Allele {
            name: name.to_string(),
            gene,
            codes,
            anchor,
            regions: RegionStarts::default(),
            aliases: Vec::new(),
        });
        Ok(id)
    }

    /// Load from FASTA. Header: `>NAME [anchor=POS] [fr1= cdr1= fr2= cdr2= fr3=] [aliases=A,B]`
    /// (0-based positions; unknown keys are ignored).
    pub fn from_fasta(kind: SegmentKind, text: &str) -> Result<Self, GermlineError> {
        let mut set = Self::new(kind);
        for record in parse_fasta(text)? {
            let name = record.id();
            let position = |key: &str| {
                record
                    .attribute(key)
                    .map(|v| parse_position(name, key, v))
                    .transpose()
            };
            let anchor = position("anchor")?;
            let mut values = [None; 5];
            for (slot, key) in values.iter_mut().zip(RegionStarts::KEYS) {
                *slot = position(key)?;
            }
            let regions = RegionStarts::from_values(values);
            regions
                .check(anchor, record.sequence.len())
                .map_err(|reason| GermlineError::BadRegions {
                    name: name.to_string(),
                    reason,
                })?;
            if anchor.is_some_and(|a| a as usize + 3 > record.sequence.len()) {
                return Err(GermlineError::BadRegions {
                    name: name.to_string(),
                    reason: "the anchor codon runs past the sequence end".to_string(),
                });
            }
            let id = set.push(name, &record.sequence, anchor)?;
            let allele = &mut set.alleles[id as usize];
            allele.regions = regions;
            if let Some(list) = record.attribute("aliases") {
                allele.aliases = list
                    .split(',')
                    .filter(|a| !a.is_empty())
                    .map(str::to_string)
                    .collect();
            }
        }
        Ok(set)
    }

    /// Merge identical sequences into the first allele carrying them (later
    /// names become aliases). Lossless for assignment: the merged alleles would
    /// always score identically, and reports list every name. A name used for
    /// two different sequences is an error.
    pub fn deduplicate(&mut self) -> Result<Vec<MergedDuplicate>, GermlineError> {
        let mut by_sequence: HashMap<Vec<u8>, usize> = HashMap::with_capacity(self.alleles.len());
        let mut by_name: HashMap<String, usize> = HashMap::with_capacity(self.alleles.len());
        let mut kept: Vec<Allele> = Vec::with_capacity(self.alleles.len());
        let mut merged = Vec::new();
        for allele in std::mem::take(&mut self.alleles) {
            match by_sequence.get(&allele.codes) {
                Some(&index) => {
                    let target = &mut kept[index];
                    for name in allele.names() {
                        if name != target.name && !target.aliases.iter().any(|a| a == name) {
                            target.aliases.push(name.to_string());
                            merged.push(MergedDuplicate {
                                kept: target.name.clone(),
                                merged: name.to_string(),
                            });
                        }
                    }
                    by_name.entry(allele.name.clone()).or_insert(index);
                }
                None => {
                    if by_name.contains_key(&allele.name) {
                        return Err(GermlineError::ConflictingDuplicate { name: allele.name });
                    }
                    by_name.insert(allele.name.clone(), kept.len());
                    by_sequence.insert(allele.codes.clone(), kept.len());
                    kept.push(allele);
                }
            }
        }
        let mut rebuilt = Self::new(self.kind);
        for allele in kept {
            let id = rebuilt.push_codes(&allele.name, allele.codes, allele.anchor)?;
            let target = &mut rebuilt.alleles[id as usize];
            target.regions = allele.regions;
            target.aliases = allele.aliases;
        }
        *self = rebuilt;
        Ok(merged)
    }

    /// Write as FASTA with every known annotation in the header; `extra`
    /// appends further `key=value` tokens per allele.
    pub fn write_fasta<W: Write>(
        &self,
        out: &mut W,
        extra: impl Fn(AlleleId) -> Option<String>,
    ) -> io::Result<()> {
        for (id, allele) in self.alleles.iter().enumerate() {
            write!(out, ">{}", allele.name)?;
            if let Some(anchor) = allele.anchor {
                write!(out, " anchor={anchor}")?;
            }
            for (key, value) in RegionStarts::KEYS.iter().zip(allele.regions.values()) {
                if let Some(value) = value {
                    write!(out, " {key}={value}")?;
                }
            }
            if !allele.aliases.is_empty() {
                write!(out, " aliases={}", allele.aliases.join(","))?;
            }
            if let Some(tokens) = extra(id as AlleleId) {
                write!(out, " {tokens}")?;
            }
            writeln!(out, "\n{}", alphabet::decode(&allele.codes))?;
        }
        Ok(())
    }

    /// Mutable access for annotation passes (sequence and name stay fixed).
    pub fn annotate_mut(
        &mut self,
    ) -> impl Iterator<Item = (&[u8], &mut Option<u32>, &mut RegionStarts)> {
        self.alleles
            .iter_mut()
            .map(|a| (a.codes.as_slice(), &mut a.anchor, &mut a.regions))
    }

    pub fn kind(&self) -> SegmentKind {
        self.kind
    }

    pub fn len(&self) -> usize {
        self.alleles.len()
    }

    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }

    pub fn alleles(&self) -> &[Allele] {
        &self.alleles
    }

    #[inline]
    pub fn allele(&self, id: AlleleId) -> &Allele {
        &self.alleles[id as usize]
    }

    pub fn gene_count(&self) -> usize {
        self.genes.len()
    }

    pub fn gene_name(&self, gene: GeneId) -> &str {
        &self.genes[gene as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_alleles_by_gene() {
        let text = ">IGHV1-2*01\nACGT\n>IGHV1-2*02 anchor=1\nACGA\n>IGHV3-23*01\nTTTT\n";
        let set = GermlineSet::from_fasta(SegmentKind::V, text).unwrap();
        assert_eq!(set.len(), 3);
        assert_eq!(set.gene_count(), 2);
        assert_eq!(set.allele(0).gene, set.allele(1).gene);
        assert_eq!(set.allele(1).anchor, Some(1));
        assert_eq!(set.gene_name(set.allele(2).gene), "IGHV3-23");
    }

    #[test]
    fn parses_regions_and_round_trips_annotated_fasta() {
        let text =
            ">V1*01 anchor=12 fr1=0 cdr1=3 fr2=6 cdr2=9 fr3=9 aliases=V1D*01\nACGTACGTACGTTGT\n";
        let set = GermlineSet::from_fasta(SegmentKind::V, text).unwrap();
        let allele = set.allele(0);
        assert_eq!(
            allele.regions.values(),
            [Some(0), Some(3), Some(6), Some(9), Some(9)]
        );
        assert_eq!(allele.names().collect::<Vec<_>>(), ["V1*01", "V1D*01"]);
        let mut written = Vec::new();
        set.write_fasta(&mut written, |_| Some("annotation=header".into()))
            .unwrap();
        let again = GermlineSet::from_fasta(SegmentKind::V, std::str::from_utf8(&written).unwrap())
            .unwrap();
        assert_eq!(again.allele(0).regions, allele.regions);
        assert_eq!(again.allele(0).anchor, Some(12));
        assert!(GermlineSet::from_fasta(SegmentKind::V, ">V anchor=2 fr2=5\nACGTACGT\n").is_err());
        assert!(GermlineSet::from_fasta(SegmentKind::V, ">V cdr1=x\nACGT\n").is_err());
    }

    #[test]
    fn deduplicates_identical_sequences_into_aliases() {
        let text = ">A*01\nACGTAC\n>B*01\nACGTAC\n>C*01\nTTTT\n>A*01\nACGTAC\n";
        let mut set = GermlineSet::from_fasta(SegmentKind::V, text).unwrap();
        let merged = set.deduplicate().unwrap();
        assert_eq!(set.len(), 2);
        assert_eq!(set.allele(0).aliases, ["B*01"]);
        assert_eq!(
            merged,
            [MergedDuplicate {
                kept: "A*01".into(),
                merged: "B*01".into()
            }]
        );
        let mut conflict =
            GermlineSet::from_fasta(SegmentKind::V, ">A*01\nACGT\n>A*01\nTTTT\n").unwrap();
        assert!(matches!(
            conflict.deduplicate(),
            Err(GermlineError::ConflictingDuplicate { .. })
        ));
    }
}
