//! Germline allele databases.

use std::collections::HashMap;
use std::fmt;
use std::io::{self, Write};

use crate::alphabet;
use crate::fasta::{FastaError, parse_fasta};
use crate::translate::codon;

/// Gene segment type. The derive order (V < D < J < C) is genomic order.
///
/// V, D and J are the [`REARRANGED`](Self::REARRANGED) kinds: the placement
/// lattice, its posteriors and every per-kind array of the adjudicator are
/// built over them alone. C (the constant region) is joined to J by an RNA
/// splice at a fixed phase, with no N/P additions and no hypermutation, so it
/// is assigned after adjudication by [`crate::constant`] and never enters a
/// lattice. Non-exhaustive: further kinds may follow.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[non_exhaustive]
pub enum SegmentKind {
    V,
    D,
    J,
    C,
}

impl SegmentKind {
    pub const ALL: [SegmentKind; 4] = [
        SegmentKind::V,
        SegmentKind::D,
        SegmentKind::J,
        SegmentKind::C,
    ];

    /// The kinds a recombination lattice is built over; per-rearranged-kind
    /// arrays are indexed by [`index`](Self::index) below
    /// `REARRANGED.len()`.
    pub const REARRANGED: [SegmentKind; 3] = [SegmentKind::V, SegmentKind::D, SegmentKind::J];

    /// Whether the kind takes part in V(D)J recombination.
    #[inline]
    pub const fn is_rearranged(self) -> bool {
        !matches!(self, SegmentKind::C)
    }

    /// Dense index for per-kind arrays.
    #[inline]
    pub const fn index(self) -> usize {
        match self {
            SegmentKind::V => 0,
            SegmentKind::D => 1,
            SegmentKind::J => 2,
            SegmentKind::C => 3,
        }
    }

    pub const fn label(self) -> &'static str {
        match self {
            SegmentKind::V => "V",
            SegmentKind::D => "D",
            SegmentKind::J => "J",
            SegmentKind::C => "C",
        }
    }
}

/// IMGT unique numbering: 0-based gapped nucleotide columns of the FR1,
/// CDR1, FR2, CDR2 and FR3 starts (IMGT codons 1, 27, 39, 56 and 66).
pub const IMGT_REGION_COLUMNS: [u32; 5] = [0, 78, 114, 165, 195];
/// IMGT codons of the conserved Cys23, Trp41 and Cys104 (column `3(k − 1)`).
pub const IMGT_CONSERVED_CODONS: [u32; 3] = [23, 41, 104];

/// What one source knows about where a V region starts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RegionSlot {
    /// 0-based position of the region's first base.
    Known(u32),
    /// The region starts before the sequence does: a 5′-partial sequence
    /// that begins after the region's first base. Nothing may fill it.
    BeforeStart,
    /// Not known; a value derived from motifs or relatives may fill it.
    Unknown,
}

impl RegionSlot {
    pub fn known(self) -> Option<u32> {
        match self {
            Self::Known(position) => Some(position),
            Self::BeforeStart | Self::Unknown => None,
        }
    }
}

/// The gap runs of one record placed on its gapped columns.
struct GapColumns {
    /// `(first column, length, bases before)` per run.
    runs: Vec<(u32, u32, u32)>,
    len: u32,
}

impl GapColumns {
    fn new(gaps: &[(u32, u32)], len: u32) -> Self {
        let mut dots = 0;
        let runs = gaps
            .iter()
            .map(|&(before, run)| {
                let start = before + dots;
                dots += run;
                (start, run, before)
            })
            .collect();
        Self { runs, len }
    }

    /// The base at or after `column` (`None` past the last base), and
    /// whether `column` itself holds it.
    fn locate(&self, column: u32) -> (Option<u32>, bool) {
        let mut dots = 0;
        for &(start, run, before) in &self.runs {
            if column < start {
                break;
            }
            if column < start + run {
                return ((before < self.len).then_some(before), false);
            }
            dots += run;
        }
        let base = column - dots;
        ((base < self.len).then_some(base), true)
    }

    /// The base in `column`, if the column holds one.
    fn base(&self, column: u32) -> Option<u32> {
        let (base, exact) = self.locate(column);
        base.filter(|_| exact)
    }

    /// First base of the codon in columns `column..column + 3`, if all
    /// three hold bases.
    fn codon(&self, column: u32) -> Option<u32> {
        self.base(column)
            .filter(|&first| self.base(column + 2) == Some(first + 2))
    }
}

/// The IMGT unique numbering of a V sequence, read from the gap columns of
/// an IMGT-gapped FASTA record: one `.` for each IMGT position the sequence
/// lacks, so gapped column `c` (0-based) is IMGT nucleotide `c + 1`.
///
/// A region whose first column lies inside the leading gap run starts
/// before the sequence ([`RegionSlot::BeforeStart`]). FR1 does so only once
/// the leading run reaches CDR1; before that it starts at the first
/// complete IMGT codon, base `(3 − lead mod 3) mod 3`. A first column
/// inside an interior run maps to the next base, and one past the last base
/// is [`RegionSlot::Unknown`]. The mapping stays exact across frameshift
/// markers, since every column, gap or base, is one IMGT position.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImgtNumbering {
    /// FR1, CDR1, FR2, CDR2 and FR3 starts.
    pub slots: [RegionSlot; 5],
    /// First base of IMGT codon 104 (Cys104) when its three columns hold bases.
    pub anchor: Option<u32>,
    /// The last base lies before the last column of codon 104: a 3′-partial
    /// sequence that holds no Cys104.
    pub ends_before_cys104: bool,
    /// Columns where an interior run of one or two `.` starts: IMGT's
    /// frameshift markers in pseudogenes, which keep later codons in frame.
    pub frameshift_markers: Vec<u32>,
    /// Every other interior run spans whole codons from a codon boundary, as
    /// IMGT gaps do. A test of the database's format, not of one record.
    pub codon_aligned: bool,
    /// Residues of IMGT codons 23, 41 and 104 where their columns hold bases.
    pub conserved: [Option<u8>; 3],
}

impl ImgtNumbering {
    /// Numbering of a record with gap runs `gaps` (as
    /// [`crate::fasta::FastaRecord::gaps`]) and encoded bases `codes`;
    /// `None` for an ungapped record, which carries no numbering.
    pub fn read(gaps: &[(u32, u32)], codes: &[u8]) -> Option<Self> {
        if gaps.is_empty() {
            return None;
        }
        let len = u32::try_from(codes.len()).expect("allele length fits u32");
        let columns = GapColumns::new(gaps, len);
        let lead = gaps.first().filter(|g| g.0 == 0).map_or(0, |g| g.1);
        let slots = IMGT_REGION_COLUMNS.map(|column| {
            let base = if column == 0 && lead < IMGT_REGION_COLUMNS[1] {
                Some((3 - lead % 3) % 3).filter(|&b| b < len)
            } else if column < lead {
                return RegionSlot::BeforeStart;
            } else {
                columns.locate(column).0
            };
            base.map_or(RegionSlot::Unknown, RegionSlot::Known)
        });
        let mut frameshift_markers = Vec::new();
        let mut codon_aligned = true;
        for &(start, run, before) in &columns.runs {
            if before == 0 || before >= len {
                continue;
            }
            if run < 3 {
                frameshift_markers.push(start);
            } else if run % 3 != 0 || start % 3 != 0 {
                codon_aligned = false;
            }
        }
        let cys104 = 3 * (IMGT_CONSERVED_CODONS[2] - 1);
        let dots_before_last: u32 = gaps.iter().filter(|g| g.0 < len).map(|g| g.1).sum();
        Some(Self {
            slots,
            anchor: columns.codon(cys104),
            ends_before_cys104: len - 1 + dots_before_last < cys104 + 2,
            frameshift_markers,
            codon_aligned,
            conserved: IMGT_CONSERVED_CODONS.map(|k| {
                columns.codon(3 * (k - 1)).map(|b| {
                    let b = b as usize;
                    codon(codes[b], codes[b + 1], codes[b + 2])
                })
            }),
        })
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
    /// IMGT numbering read from the input's gap columns (V only), as read:
    /// annotation decides whether the database's gaps follow IMGT and
    /// merges them under the header keys.
    pub imgt: Option<Box<ImgtNumbering>>,
}

impl Allele {
    /// The name followed by its aliases.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        std::iter::once(self.name.as_str()).chain(self.aliases.iter().map(String::as_str))
    }

    /// Take what this allele lacks from `copy`, an identical sequence: its
    /// anchor when this allele's region starts fit it, its region starts
    /// (only as a whole, and only when they fit this allele's anchor) and
    /// its IMGT numbering. Values this allele has are kept, so the allele
    /// stays consistent ([`RegionStarts::check`]).
    fn absorb(&mut self, copy: Allele) {
        let len = self.codes.len();
        if self.anchor.is_none() && self.regions.check(copy.anchor, len).is_ok() {
            self.anchor = copy.anchor;
        }
        if self.regions.is_empty() && copy.regions.check(self.anchor, len).is_ok() {
            self.regions = copy.regions;
        }
        self.imgt = self.imgt.take().or(copy.imgt);
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
    /// Summed allele lengths (the search space of E-values).
    total_len: usize,
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
            total_len: 0,
        }
    }

    /// Allele label of a FASTA id. IMGT reference headers read
    /// `accession|allele|species|functionality|…`, so a `|`-separated id
    /// gives its first field that names an allele (`GENE*NN`), else its
    /// first field.
    pub fn allele_label(id: &str) -> &str {
        let mut fields = id.split('|');
        match id.split('|').find(|field| field.contains('*')) {
            Some(allele) => allele,
            None => fields.next().unwrap_or(id),
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
        // In u64: on 32-bit targets `a as usize + 3` can wrap.
        if anchor.is_some_and(|a| u64::from(a) + 3 > codes.len() as u64) {
            return Err(GermlineError::BadRegions {
                name: name.to_string(),
                reason: "the anchor codon runs past the sequence end".to_string(),
            });
        }
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
        self.total_len += codes.len();
        self.alleles.push(Allele {
            name: name.to_string(),
            gene,
            codes,
            anchor,
            regions: RegionStarts::default(),
            aliases: Vec::new(),
            imgt: None,
        });
        Ok(id)
    }

    /// Load from FASTA. Header: `>NAME [anchor=POS] [fr1= cdr1= fr2= cdr2= fr3=] [aliases=A,B]`
    /// (0-based positions; unknown keys are ignored). A gapped V record
    /// keeps its IMGT numbering in [`Allele::imgt`].
    pub fn from_fasta(kind: SegmentKind, text: &str) -> Result<Self, GermlineError> {
        let mut set = Self::new(kind);
        for record in parse_fasta(text)? {
            let name = Self::allele_label(record.id());
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
            let id = set.push(name, &record.sequence, anchor)?;
            let allele = &mut set.alleles[id as usize];
            allele.regions = regions;
            if kind == SegmentKind::V {
                allele.imgt = ImgtNumbering::read(&record.gaps, &allele.codes).map(Box::new);
            }
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

    /// The same alleles — ids, names, aliases and genes — with an `N` anchor
    /// column prepended to every sequence (`C′ = [N] + C`), the reference of
    /// an anchored alignment ([`crate::scoring::AnchorColumn`]). Reference
    /// position `r` of an allele is position `r + 1` of its copy; annotations
    /// are not carried.
    pub fn with_anchor_column(&self) -> Self {
        let mut set = Self::new(self.kind);
        for allele in &self.alleles {
            let codes = std::iter::once(alphabet::N)
                .chain(allele.codes.iter().copied())
                .collect();
            set.push_codes(&allele.name, codes, None)
                .expect("a copied allele is valid");
            set.alleles.last_mut().expect("just pushed").aliases = allele.aliases.clone();
        }
        set
    }

    /// Merge identical sequences into the first allele carrying them (later
    /// names become aliases). Lossless for assignment: the merged alleles would
    /// always score identically, and reports list every name. The kept
    /// allele takes from later copies what it lacks, as long as it stays
    /// consistent: an anchor, region starts, an IMGT numbering; values it
    /// has are kept. A name or alias used for two different sequences is
    /// an error.
    pub fn deduplicate(&mut self) -> Result<Vec<MergedDuplicate>, GermlineError> {
        let mut by_sequence: HashMap<Vec<u8>, usize> = HashMap::with_capacity(self.alleles.len());
        let mut by_name: HashMap<String, usize> = HashMap::with_capacity(self.alleles.len());
        let mut kept: Vec<Allele> = Vec::with_capacity(self.alleles.len());
        let mut merged = Vec::new();
        for allele in std::mem::take(&mut self.alleles) {
            let index = by_sequence
                .get(&allele.codes)
                .copied()
                .unwrap_or(kept.len());
            if let Some(name) = allele
                .names()
                .find(|name| by_name.get(*name).is_some_and(|&at| at != index))
            {
                let name = name.to_string();
                return Err(GermlineError::ConflictingDuplicate { name });
            }
            by_name.extend(allele.names().map(|name| (name.to_string(), index)));
            if index == kept.len() {
                by_sequence.insert(allele.codes.clone(), index);
                kept.push(allele);
                continue;
            }
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
            target.absorb(allele);
        }
        let mut rebuilt = Self::new(self.kind);
        for allele in kept {
            let id = rebuilt.push_codes(&allele.name, allele.codes, allele.anchor)?;
            let target = &mut rebuilt.alleles[id as usize];
            target.regions = allele.regions;
            target.aliases = allele.aliases;
            target.imgt = allele.imgt;
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

    /// Summed length of all alleles.
    pub fn total_len(&self) -> usize {
        self.total_len
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

    /// Per-rearranged-kind arrays (`[_; 3]` indexed by `index()`) cover V, D
    /// and J exactly; C indexes past them and is the one non-rearranged kind.
    #[test]
    fn rearranged_kinds_index_the_per_kind_arrays() {
        for (i, kind) in SegmentKind::REARRANGED.into_iter().enumerate() {
            assert_eq!(kind.index(), i);
            assert!(kind.is_rearranged());
        }
        assert_eq!(SegmentKind::C.index(), SegmentKind::REARRANGED.len());
        assert!(!SegmentKind::C.is_rearranged());
        assert_eq!(SegmentKind::ALL[..3], SegmentKind::REARRANGED);
        assert_eq!(SegmentKind::C.label(), "C");
    }

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
    fn reads_allele_names_from_imgt_headers() {
        let text = ">M99641|IGHV1-18*01|Homo_sapiens|F|V-REGION|188..483|296 nt|1| | |\nACGT\n\
                    >X|IGHV1-18*02|Homo_sapiens|F|\nACGA\n>plain*01 anchor=1\nACGG\n>A|B\nTTTT\n";
        let set = GermlineSet::from_fasta(SegmentKind::V, text).unwrap();
        let names: Vec<&str> = set.alleles().iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, ["IGHV1-18*01", "IGHV1-18*02", "plain*01", "A"]);
        assert_eq!(set.gene_count(), 3);
        assert_eq!(set.allele(2).anchor, Some(1));
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
        let conflicts = [
            ">A*01\nACGT\n>A*01\nTTTT\n",
            // A*01 names X, then joins Y's allele as an alias.
            ">A*01\nACGT\n>B*01\nTTTT\n>A*01\nTTTT\n",
            ">A*01 aliases=C*01\nACGT\n>C*01\nTTTT\n",
        ];
        for text in conflicts {
            let mut set = GermlineSet::from_fasta(SegmentKind::V, text).unwrap();
            assert!(
                matches!(
                    set.deduplicate(),
                    Err(GermlineError::ConflictingDuplicate { .. })
                ),
                "{text}"
            );
        }
    }

    /// A header value of a merged copy is kept: the copy is the same
    /// sequence, so its anchor and region starts hold for the kept allele.
    #[test]
    fn deduplication_keeps_the_annotation_of_merged_copies() {
        let sequence = "ACGTACGTACGTACGTACGTACGT";
        let text = format!(">V1*01\n{sequence}\n>V1*02 anchor=12 fr1=0 cdr1=3\n{sequence}\n");
        let mut set = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
        set.deduplicate().unwrap();
        let kept = set.allele(0);
        assert_eq!((kept.name.as_str(), kept.anchor), ("V1*01", Some(12)));
        assert_eq!(kept.regions.values(), [Some(0), Some(3), None, None, None]);
        // A value the kept allele has wins.
        let text = format!(">V1*01 anchor=9\n{sequence}\n>V1*02 anchor=12\n{sequence}\n");
        let mut set = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
        set.deduplicate().unwrap();
        assert_eq!(set.allele(0).anchor, Some(9));
        // A copy's anchor before the kept allele's region starts is left
        // out, so the allele still passes the checks a reload applies.
        let text = format!(">V1*01 fr1=0 cdr1=20\n{sequence}\n>V1*02 anchor=12\n{sequence}\n");
        let mut set = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
        set.deduplicate().unwrap();
        let kept = set.allele(0);
        assert_eq!(kept.anchor, None);
        assert!(kept.regions.check(kept.anchor, kept.codes.len()).is_ok());
        let mut written = Vec::new();
        set.write_fasta(&mut written, |_| None).unwrap();
        let reloaded = String::from_utf8(written).unwrap();
        assert!(GermlineSet::from_fasta(SegmentKind::V, &reloaded).is_ok());
        // Anchors past the end are refused, also where `usize` is 32 bits.
        assert!(set.push("V*01", b"ACGTAC", Some(4)).is_err());
        assert!(set.push("V*02", b"ACGTAC", Some(u32::MAX)).is_err());
    }

    /// A random IMGT-like gapped V record: 106 IMGT codons with Cys23,
    /// Trp41 and Cys104 most of the time, whole-codon gaps where IMGT
    /// leaves positions empty in some loci, a leading run of 0–90 columns
    /// (5′-partial), and optionally a one-column frameshift marker and a 3′
    /// truncation.
    fn random_gapped(rng: &mut crate::sim::Rng) -> Vec<u8> {
        const GAPPABLE: [usize; 13] = [10, 31, 32, 33, 34, 35, 57, 58, 59, 60, 73, 81, 82];
        let mut columns = Vec::with_capacity(318);
        for codon in 1..=106 {
            match codon {
                _ if GAPPABLE.contains(&codon) && rng.below(2) == 0 => columns.extend(b"..."),
                23 | 104 if rng.below(8) > 0 => columns.extend(b"TGT"),
                41 if rng.below(8) > 0 => columns.extend(b"TGG"),
                _ => columns.extend(rng.bases(3)),
            }
        }
        let lead = rng.below(91);
        columns[..lead].fill(b'.');
        if rng.below(2) == 0 {
            let at = lead + 1 + rng.below(columns.len() - lead - 1);
            columns[at] = b'.';
        }
        if rng.below(3) == 0 {
            let first_base = columns.iter().position(|&c| c != b'.').unwrap();
            let end = first_base + 1 + rng.below(columns.len() - first_base);
            columns.truncate(end);
        }
        columns
    }

    /// The numbering by walking the gapped columns one by one.
    fn column_walk(gapped: &[u8]) -> ImgtNumbering {
        let column_of: Vec<u32> = (0..gapped.len() as u32)
            .filter(|&c| gapped[c as usize] != b'.')
            .collect();
        let codes = alphabet::encode(
            &gapped
                .iter()
                .copied()
                .filter(|&c| c != b'.')
                .collect::<Vec<_>>(),
        );
        let len = column_of.len() as u32;
        let lead = column_of[0];
        let in_column = |c: u32| column_of.iter().position(|&x| x == c).map(|b| b as u32);
        let codon_at = |c: u32| {
            let first = in_column(c)?;
            (in_column(c + 1)? == first + 1 && in_column(c + 2)? == first + 2).then_some(first)
        };
        let slots = IMGT_REGION_COLUMNS.map(|column| {
            if column == 0 {
                return match (3 - lead % 3) % 3 {
                    _ if lead >= 78 => RegionSlot::BeforeStart,
                    b if b < len => RegionSlot::Known(b),
                    _ => RegionSlot::Unknown,
                };
            }
            if column < lead {
                return RegionSlot::BeforeStart;
            }
            column_of
                .iter()
                .position(|&c| c >= column)
                .map_or(RegionSlot::Unknown, |b| RegionSlot::Known(b as u32))
        });
        let mut frameshift_markers = Vec::new();
        let mut codon_aligned = true;
        let mut column = 0;
        while column < gapped.len() {
            let run = gapped[column..].iter().take_while(|&&c| c == b'.').count();
            let interior = column > 0 && gapped[column + run..].iter().any(|&c| c != b'.');
            if run > 0 && interior {
                if run < 3 {
                    frameshift_markers.push(column as u32);
                } else if run % 3 != 0 || column % 3 != 0 {
                    codon_aligned = false;
                }
            }
            column += run.max(1);
        }
        ImgtNumbering {
            slots,
            anchor: codon_at(309),
            ends_before_cys104: *column_of.last().unwrap() < 311,
            frameshift_markers,
            codon_aligned,
            conserved: [23, 41, 104].map(|k| {
                codon_at(3 * (k - 1)).map(|b| {
                    let b = b as usize;
                    codon(codes[b], codes[b + 1], codes[b + 2])
                })
            }),
        }
    }

    #[test]
    fn imgt_numbering_matches_a_column_walk() {
        let mut rng = crate::sim::Rng::new(23);
        let mut seen = [0usize; 5];
        for _ in 0..3000 {
            let gapped = random_gapped(&mut rng);
            let text = format!(">V\n{}\n", String::from_utf8_lossy(&gapped));
            let set = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
            let Some(numbering) = set.allele(0).imgt.as_deref() else {
                assert!(!gapped.contains(&b'.'), "{text}");
                continue;
            };
            assert_eq!(*numbering, column_walk(&gapped), "{text}");
            seen[0] += usize::from(numbering.slots[1] == RegionSlot::BeforeStart);
            seen[1] += usize::from(numbering.ends_before_cys104);
            seen[2] += usize::from(!numbering.frameshift_markers.is_empty());
            seen[3] += usize::from(numbering.anchor.is_some());
            seen[4] += usize::from(!numbering.codon_aligned);
        }
        assert!(
            seen.iter().all(|&n| n > 20),
            "every case exercised: {seen:?}"
        );
    }

    /// The gap layouts `(bases before, run)` of three IMGT human IGHV
    /// records, on random bases (the numbering depends on the layout only).
    #[test]
    fn imgt_layouts_of_real_records() {
        let mut rng = crate::sim::Rng::new(24);
        let mut read = |len: usize, gaps: &[(u32, u32)]| {
            ImgtNumbering::read(gaps, &rng.codes(len)).expect("gapped")
        };
        let known = |starts: [u32; 5]| starts.map(RegionSlot::Known);
        // IGHV1-2*02: full length.
        let full = read(296, &[(27, 3), (87, 12), (162, 6), (195, 3)]);
        assert_eq!(full.slots, known([0, 75, 99, 150, 174]));
        assert_eq!(full.anchor, Some(285));
        assert!(full.codon_aligned && full.frameshift_markers.is_empty());
        assert!(!full.ends_before_cys104);
        // IGHV5-10-1*02: a one-column frameshift marker in FR3 (codon 96).
        let shifted = read(295, &[(27, 3), (87, 12), (162, 6), (195, 3), (263, 1)]);
        assert_eq!(shifted.slots, known([0, 75, 99, 150, 174]));
        assert_eq!(shifted.anchor, Some(284));
        assert_eq!(shifted.frameshift_markers, [287]);
        assert!(shifted.codon_aligned);
        // IGHV3-72*02: partial in 5′ (starts inside CDR1) and in 3′.
        let partial = read(165, &[(0, 84), (6, 12), (120, 3)]);
        assert_eq!(
            partial.slots,
            [
                RegionSlot::BeforeStart,
                RegionSlot::BeforeStart,
                RegionSlot::Known(18),
                RegionSlot::Known(69),
                RegionSlot::Known(99)
            ]
        );
        assert_eq!(partial.anchor, None);
        assert!(partial.ends_before_cys104);
        assert!(ImgtNumbering::read(&[], &rng.codes(296)).is_none());
    }

    #[test]
    fn numbering_is_read_for_v_only_and_survives_deduplication() {
        let gapped = "CAGGTG...CAGTGT";
        let plain = gapped.replace('.', "");
        let text = format!(">V1*01\n{plain}\n>V1*02\n{gapped}\n>V2*01 fr1=0\n{gapped}A\n");
        let mut v = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
        assert!(v.allele(0).imgt.is_none());
        let numbered = v.allele(1).imgt.clone().expect("gapped record");
        // Header keys and gap numbering are stored side by side, unmerged.
        assert_eq!(v.allele(2).regions.fr1, Some(0));
        assert!(v.allele(2).imgt.is_some());
        v.deduplicate().unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v.allele(0).aliases, ["V1*02"]);
        assert_eq!(v.allele(0).imgt, Some(numbered));
        let j = GermlineSet::from_fasta(SegmentKind::J, &format!(">J1*01\n{gapped}\n")).unwrap();
        assert!(j.allele(0).imgt.is_none());
    }
}
