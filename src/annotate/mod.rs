//! Germline database annotation: anchors, reading frames and V region boundaries.
//!
//! V and J germline genes are in frame with their conserved anchors (V Cys104,
//! J Trp/Phe118), so an anchor fixes the reading frame and, for V, the region
//! boundaries follow from three conserved residues (IMGT definitions):
//!
//! ```text
//! FR1  = [frame start, Cys23 + 4 codons)      CDR1 = [Cys23 + 4, Trp41 − 2)
//! FR2  = [Trp41 − 2, Trp41 + 15)              CDR2 = [Trp41 + 15, FR3 start)
//! FR3  = [Cys104 − (L − 1), Cys104 + 1)       L = 38 codons heavy, 36 light
//! CDR3 = [Cys104 + 1, J Trp/Phe118)           FR4  = [J Trp/Phe118, J end)
//! ```
//!
//! Sources, in order of trust (a given value is never overridden):
//!
//! 1. `anchor=` / `fr1=` … `fr3=` in the FASTA header, then the IMGT
//!    numbering of an IMGT-gapped V set ([`ImgtGaps`]). A region the
//!    numbering places before a 5′-partial sequence stays empty, and an
//!    allele whose numbering ends before Cys104 gets no anchor.
//! 2. Motif search on a one-pass three-frame translation ([`find_v`], [`find_j`]).
//! 3. The same search with a database-consensus prior on the anchor's
//!    distance from the 3′ end, for alleles whose motifs were weak.
//! 4. Transfer from the most similar confidently annotated allele by
//!    alignment ([`RelativeIndex`], [`transfer()`]) — rescues pseudogenes and
//!    mutated motifs. When the relative's anchor codon maps at or past the
//!    allele's 3′ end (an IMGT 3′-partial sequence), the allele holds no
//!    anchor and none is reported: a weak motif hit there is another codon.
//!
//! The chain type (which fixes the FR3 length) is voted by J residues
//! (`W` heavy, `F` light) and Trp41 contexts (`W.R` heavy, `W.Q` light).
//! Everything runs once per database: one translation per allele plus an
//! alignment for the rare transferred allele — milliseconds for thousands of
//! alleles, and nothing is added to per-read work.

mod given;
mod j_motif;
mod transfer;
mod v_motif;

use std::collections::HashMap;
use std::io::{self, Write};

pub use given::ImgtGaps;
pub use j_motif::{JHit, find_j};
pub use transfer::{
    MIN_SHARED_FRACTION, MIN_SHARED_KMERS, RelativeIndex, kmers, shared, transfer, transfer_aligned,
};
pub use v_motif::{VHit, VSearch, find_v};

use crate::align::GotohWorkspace;
use crate::germline::{GermlineSet, RegionSlot, RegionStarts};
use crate::translate::codon;
use given::Given;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChainType {
    Heavy,
    Light,
}

impl ChainType {
    /// FR3 length in codons including Cys104: IMGT 66–104 minus the constant
    /// gaps (heavy: position 73; kappa/lambda: 73, 81, 82).
    pub fn fr3_codons(self) -> u32 {
        match self {
            Self::Heavy => 38,
            Self::Light => 36,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Heavy => "heavy",
            Self::Light => "light",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "heavy" => Some(Self::Heavy),
            "light" => Some(Self::Light),
            _ => None,
        }
    }
}

/// How an allele's anchor was obtained.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AnchorSource {
    Header,
    /// The IMGT numbering of a gapped input record.
    ImgtGaps,
    Motif,
    Consensus,
    Transfer {
        from: String,
    },
    LowConfidence,
    Missing,
}

impl AnchorSource {
    pub fn label(&self) -> String {
        match self {
            Self::Header => "header".into(),
            Self::ImgtGaps => "imgt-gaps".into(),
            Self::Motif => "motif".into(),
            Self::Consensus => "consensus".into(),
            Self::Transfer { from } => format!("transfer:{from}"),
            Self::LowConfidence => "low-confidence".into(),
            Self::Missing => "missing".into(),
        }
    }

    /// Reliable enough to serve as a transfer donor.
    fn trusted(&self) -> bool {
        matches!(
            self,
            Self::Header | Self::ImgtGaps | Self::Motif | Self::Consensus
        )
    }
}

#[derive(Clone, Debug)]
pub struct AlleleAnnotation {
    pub source: AnchorSource,
    pub notes: Vec<String>,
    /// Region starts (FR1, CDR1, FR2, CDR2, FR3) that may be wrong: derived
    /// for an allele whose anchor is low-confidence or missing, a CDR1
    /// start from an ambiguous Cys23, or a derived bound of a CDR of
    /// unusual length. Given values are never uncertain. An error in a
    /// start not flagged here is a silent error.
    pub uncertain: [bool; 5],
}

impl AlleleAnnotation {
    fn new(source: AnchorSource, notes: Vec<String>) -> Self {
        Self {
            source,
            notes,
            uncertain: [false; 5],
        }
    }

    /// The notes, then `uncertain=` and the keys of the uncertain starts
    /// (`uncertain=cdr1,fr3`), as one `; `-separated field.
    pub fn notes_field(&self) -> String {
        let keys: Vec<&str> = RegionStarts::KEYS
            .iter()
            .zip(self.uncertain)
            .filter_map(|(&key, uncertain)| uncertain.then_some(key))
            .collect();
        let mut fields = self.notes.clone();
        if !keys.is_empty() {
            fields.push(format!("uncertain={}", keys.join(",")));
        }
        fields.join("; ")
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct AnnotateOptions {
    /// Chain type; `None` votes it from the database.
    pub chain: Option<ChainType>,
}

/// Provenance of every allele's annotation, aligned with allele ids.
#[derive(Clone, Debug, Default)]
pub struct DatabaseAnnotation {
    pub chain: Option<ChainType>,
    /// Consensus distance (nt) from the anchor codon to the 3′ end.
    pub v_tail: Option<u32>,
    pub j_tail: Option<u32>,
    pub v: Vec<AlleleAnnotation>,
    pub j: Vec<AlleleAnnotation>,
    /// Whether the V set's IMGT gap columns were read.
    pub imgt_gaps: ImgtGaps,
}

fn mode(values: impl Iterator<Item = u32>) -> Option<u32> {
    mode_with_share(values).0
}

/// Most frequent value (smallest on ties) and its share of all values.
fn mode_with_share(values: impl Iterator<Item = u32>) -> (Option<u32>, f64) {
    let mut counts: HashMap<u32, usize> = HashMap::new();
    let mut total = 0usize;
    values.for_each(|v| {
        *counts.entry(v).or_default() += 1;
        total += 1;
    });
    counts
        .into_iter()
        .max_by_key(|&(value, count)| (count, std::cmp::Reverse(value)))
        .map_or((None, 0.0), |(value, count)| {
            (Some(value), count as f64 / total as f64)
        })
}

fn residue_at(codes: &[u8], anchor: u32) -> u8 {
    let a = anchor as usize;
    codes
        .get(a..a + 3)
        .map_or(b'?', |c| codon(c[0], c[1], c[2]))
}

/// Region starts implied by a V hit (IMGT definitions, see module docs).
fn regions_from_hit(hit: &VHit, chain: Option<ChainType>) -> RegionStarts {
    RegionStarts {
        fr1: Some(hit.frame as u32),
        cdr1: hit.cys23_nt().map(|c| c + 12),
        fr2: hit.trp41_nt().and_then(|w| w.checked_sub(6)),
        cdr2: hit.trp41_nt().map(|w| w + 45),
        fr3: chain.and_then(|c| hit.anchor().checked_sub(3 * (c.fr3_codons() - 1))),
    }
}

/// Given values (header, then IMGT gaps) win; derived values fill the
/// unknown slots, never one the IMGT numbering places before the sequence
/// start. Derived values that break the ordering are dropped, latest first,
/// with a note. A CDR of unusual length is noted, and its derived bounds
/// are marked uncertain.
fn merge_regions(
    given: &[RegionSlot; 5],
    derived: RegionStarts,
    anchor: Option<u32>,
    len: usize,
    annotation: &mut AlleleAnnotation,
) -> RegionStarts {
    let notes = &mut annotation.notes;
    let mut values = given.map(RegionSlot::known);
    for ((value, slot), derived) in values.iter_mut().zip(given).zip(derived.values()) {
        if *slot == RegionSlot::Unknown {
            *value = derived;
        }
    }
    let mut merged = RegionStarts::from_values(values);
    while let Err(reason) = merged.check(anchor, len) {
        let Some(index) = (0..5)
            .rev()
            .find(|&i| merged.values()[i].is_some() && given[i].known().is_none())
        else {
            notes.push(format!("header regions inconsistent: {reason}"));
            break;
        };
        notes.push(format!(
            "dropped derived {}: {reason}",
            RegionStarts::KEYS[index]
        ));
        let mut values = merged.values();
        values[index] = None;
        merged = RegionStarts::from_values(values);
    }
    let values = merged.values();
    for (region, [start, end], usual) in [("CDR1", [1, 2], 4..=14), ("CDR2", [3, 4], 0..=12)] {
        let Some(codons) = values[start]
            .zip(values[end])
            .map(|(a, b)| (i64::from(b) - i64::from(a)) / 3)
        else {
            continue;
        };
        if !usual.contains(&codons) {
            notes.push(format!("unusual {region} length {codons} codons"));
            for bound in [start, end] {
                annotation.uncertain[bound] |= given[bound].known().is_none();
            }
        }
    }
    merged
}

/// Positions carried over from the nearest eligible relative.
struct Relative {
    donor: usize,
    mapped: Vec<Option<u32>>,
    /// The target's and the donor's bases inside their local alignment.
    aligned: std::ops::Range<u32>,
    donor_aligned: std::ops::Range<u32>,
}

impl Relative {
    /// Whether the target ends, on the side of the donor's position `at`,
    /// while still aligned to the donor (to within a codon): a partial
    /// sequence, not an alignment that stopped short.
    fn runs_out_before(&self, at: u32, len: usize) -> bool {
        if at < self.donor_aligned.start {
            self.aligned.start <= 3
        } else {
            self.aligned.end as usize + 3 >= len
        }
    }
}

/// Transfer `positions` from the nearest eligible relative.
fn from_relative(
    set: &GermlineSet,
    target: usize,
    eligible: &[bool],
    index: &mut RelativeIndex,
    positions_of: impl Fn(usize) -> Vec<Option<u32>>,
    workspace: &mut GotohWorkspace,
) -> Option<Relative> {
    let donor = index.nearest(target, eligible)?;
    let alleles = set.alleles();
    let (mapped, alignment) = transfer_aligned(
        &alleles[target].codes,
        &alleles[donor].codes,
        &positions_of(donor),
        workspace,
    )?;
    Some(Relative {
        donor,
        mapped,
        aligned: alignment.q_start..alignment.q_end,
        donor_aligned: alignment.r_start..alignment.r_end,
    })
}

/// The relative index of `set`, built on first use.
fn relatives<'a>(index: &'a mut Option<RelativeIndex>, set: &GermlineSet) -> &'a mut RelativeIndex {
    index.get_or_insert_with(|| {
        RelativeIndex::build(set.alleles().iter().map(|a| a.codes.as_slice()))
    })
}

fn annotate_j(
    j: &mut GermlineSet,
    frames: &mut [Vec<u8>; 3],
    workspace: &mut GotohWorkspace,
) -> (Vec<AlleleAnnotation>, Option<u32>, i32) {
    let mut result: Vec<AlleleAnnotation> = Vec::with_capacity(j.len());
    let mut anchors: Vec<Option<u32>> = Vec::with_capacity(j.len());
    let mut weak: Vec<Option<JHit>> = Vec::with_capacity(j.len());
    let mut vote = 0i32;
    for allele in j.alleles() {
        let mut notes = Vec::new();
        if let Some(anchor) = allele.anchor {
            let residue = residue_at(&allele.codes, anchor);
            if !matches!(residue, b'W' | b'F') {
                notes.push(format!("header anchor codon is {}", char::from(residue)));
            }
            result.push(AlleleAnnotation::new(AnchorSource::Header, notes));
            anchors.push(Some(anchor));
            weak.push(None);
            continue;
        }
        let hit = find_j(&allele.codes, None, frames);
        match hit.filter(|h| h.confident) {
            Some(h) => {
                vote += 2 * i32::from(h.chain_vote());
                result.push(AlleleAnnotation::new(AnchorSource::Motif, notes));
                anchors.push(Some(h.anchor()));
                weak.push(None);
            }
            None => {
                result.push(AlleleAnnotation::new(AnchorSource::Missing, notes));
                anchors.push(None);
                weak.push(hit);
            }
        }
    }
    let tail = mode(
        j.alleles()
            .iter()
            .zip(&anchors)
            .filter_map(|(a, &x)| x.map(|x| a.codes.len() as u32 - x)),
    );
    // Consensus prior for weak alleles.
    for (index, allele) in j.alleles().iter().enumerate() {
        if anchors[index].is_some() || tail.is_none() {
            continue;
        }
        if let Some(h) = find_j(&allele.codes, tail, frames).filter(|h| h.confident) {
            anchors[index] = Some(h.anchor());
            result[index].source = AnchorSource::Consensus;
            weak[index] = None;
        }
    }
    // Sibling transfer, then low-confidence motif hits. K-mer profiles are
    // built only if some allele needs a relative.
    let trusted: Vec<bool> = result.iter().map(|r| r.source.trusted()).collect();
    let mut relative_index = None;
    let any_donor = trusted.iter().any(|&t| t);
    for index in 0..j.len() {
        if anchors[index].is_some() {
            continue;
        }
        let donor = if any_donor {
            let relatives = relatives(&mut relative_index, j);
            from_relative(
                j,
                index,
                &trusted,
                relatives,
                |d| vec![anchors[d]],
                workspace,
            )
        } else {
            None
        };
        let codes_len = j.alleles()[index].codes.len();
        let fits = |r: &Relative| r.mapped[0].is_some_and(|a| a as usize + 3 <= codes_len);
        let partial =
            |r: &Relative| anchors[r.donor].is_some_and(|a| r.runs_out_before(a, codes_len));
        match donor {
            Some(Relative { donor, mapped, .. })
                if mapped[0].is_some_and(|a| a as usize + 3 <= codes_len) =>
            {
                let anchor = mapped[0].expect("checked above");
                let residue = residue_at(&j.alleles()[index].codes, anchor);
                if !matches!(residue, b'W' | b'F') {
                    result[index]
                        .notes
                        .push(format!("conserved residue is {}", char::from(residue)));
                }
                result[index].source = AnchorSource::Transfer {
                    from: j.alleles()[donor].name.clone(),
                };
                anchors[index] = Some(anchor);
            }
            // A partial allele: the relative's anchor codon maps outside it,
            // and the allele ends while still aligned to the relative.
            Some(r) if !fits(&r) && partial(&r) => {
                result[index].notes.push(format!(
                    "partial: the anchor codon of {} maps outside the sequence",
                    j.alleles()[r.donor].name
                ));
            }
            _ => {
                if let Some(h) = weak[index] {
                    result[index].source = AnchorSource::LowConfidence;
                    result[index]
                        .notes
                        .push(format!("weak motif {}G.G", char::from(h.residue)));
                    anchors[index] = Some(h.anchor());
                }
            }
        }
    }
    for ((_, anchor, _), found) in j.annotate_mut().zip(anchors) {
        *anchor = anchor.or(found);
    }
    (result, tail, vote)
}

/// Annotate `v` (anchors and region starts) and `j` (anchors) in place.
pub fn annotate(
    v: &mut GermlineSet,
    j: Option<&mut GermlineSet>,
    options: AnnotateOptions,
) -> DatabaseAnnotation {
    let mut frames = [Vec::new(), Vec::new(), Vec::new()];
    let mut workspace = GotohWorkspace::default();
    let (j_annotations, j_tail, mut vote) = match j {
        Some(j) => annotate_j(j, &mut frames, &mut workspace),
        None => (Vec::new(), None, 0),
    };

    // What the input states: header keys, then IMGT gap columns if the
    // database's gaps follow the IMGT numbering.
    let imgt_gaps = ImgtGaps::check(v);
    let mut result: Vec<AlleleAnnotation> = Vec::with_capacity(v.len());
    let given: Vec<Given> = v
        .alleles()
        .iter()
        .map(|allele| {
            let mut notes = Vec::new();
            let given = Given::read(allele, imgt_gaps.used(), &mut notes);
            result.push(AlleleAnnotation::new(AnchorSource::Missing, notes));
            given
        })
        .collect();

    // V pass 1: given anchors fix the frame; otherwise free motif search. An
    // allele whose IMGT numbering ends before Cys104 is not searched.
    let mut hits: Vec<Option<VHit>> = Vec::with_capacity(v.len());
    for ((allele, given), annotation) in v.alleles().iter().zip(&given).zip(&mut result) {
        let hit = if given.ends_before_cys104 {
            None
        } else {
            let search = VSearch {
                fixed_anchor: given.anchor,
                ..VSearch::default()
            };
            find_v(&allele.codes, search, &mut frames)
        };
        annotation.source = match (&given.source, hit) {
            (Some(source), _) => source.clone(),
            (None, Some(h)) if h.confident => AnchorSource::Motif,
            _ => AnchorSource::Missing,
        };
        if let (Some(anchor), Some(source)) = (given.anchor, &given.source) {
            let residue = char::from(residue_at(&allele.codes, anchor));
            if residue != 'C' {
                annotation.notes.push(match source {
                    AnchorSource::ImgtGaps => format!("IMGT 104 codon is {residue}"),
                    _ => format!("header anchor codon is {residue}"),
                });
            }
        }
        if annotation.source.trusted() {
            vote += hit.map_or(0, |h| i32::from(h.chain_vote));
        }
        hits.push(hit);
    }
    let chain = options.chain.or(match vote.signum() {
        1 => Some(ChainType::Heavy),
        -1 => Some(ChainType::Light),
        _ => None,
    });

    // Database consensus from pass 1: the 3′ tail after Cys104, the frame of
    // full-length alleles and the codon index of their Cys23.
    let trusted: Vec<VHit> = (0..v.len())
        .filter(|&i| result[i].source.trusted())
        .filter_map(|i| hits[i])
        .collect();
    let (tail, tail_share) = mode_with_share(
        (0..v.len())
            .filter(|&i| result[i].source.trusted())
            .filter_map(|i| {
                given[i]
                    .anchor
                    .or(hits[i].map(|h| h.anchor()))
                    .map(|x| v.alleles()[i].codes.len() as u32 - x)
            }),
    );
    let (frame_mode, _) = mode_with_share(trusted.iter().map(|h| h.frame as u32));
    let (cys23_mode, _) = mode_with_share(
        trusted
            .iter()
            .filter(|h| Some(h.frame as u32) == frame_mode)
            .filter_map(|h| h.cys23)
            .map(|c| c as u32),
    );
    let cys23_prior = cys23_mode
        .zip(frame_mode)
        .map(|(c, f)| (f as usize, c as usize));
    let v_tail = tail;

    // V pass 2: every allele again with the database priors. The tail prior
    // applies to confident alleles only when one tail dominates the database.
    for (index, allele) in v.alleles().iter().enumerate() {
        if given[index].ends_before_cys104 {
            continue;
        }
        let was_confident = result[index].source.trusted();
        let search = VSearch {
            fixed_anchor: given[index].anchor,
            tail_prior: tail.filter(|_| !was_confident || tail_share >= 0.5),
            cys23_prior,
            fr3_codons: chain.map(ChainType::fr3_codons),
        };
        let refined = find_v(&allele.codes, search, &mut frames);
        match (&result[index].source, refined) {
            (AnchorSource::Header | AnchorSource::ImgtGaps, Some(h)) => hits[index] = Some(h),
            (AnchorSource::Motif, Some(h)) if h.confident => hits[index] = Some(h),
            (AnchorSource::Missing, Some(h)) if h.confident => {
                hits[index] = Some(h);
                result[index].source = AnchorSource::Consensus;
            }
            (AnchorSource::Motif, _) => {
                result[index].source = AnchorSource::Missing;
                result[index]
                    .notes
                    .push("motifs inconsistent with the database consensus".into());
            }
            _ => {}
        }
    }
    // Final anchors and regions for trusted alleles.
    let len_of = |set: &GermlineSet, i: usize| set.alleles()[i].codes.len();
    let mut anchors: Vec<Option<u32>> = Vec::with_capacity(v.len());
    let mut regions: Vec<RegionStarts> = Vec::with_capacity(v.len());
    for (index, given) in given.iter().enumerate() {
        let trusted_hit = hits[index].filter(|_| result[index].source.trusted());
        let anchor = given.anchor.or(trusted_hit.map(|h| h.anchor()));
        let merged = match (anchor, trusted_hit) {
            (Some(anchor), Some(hit)) => {
                let annotation = &mut result[index];
                if hit.cys23_alternatives > 0 && given.regions[1] == RegionSlot::Unknown {
                    annotation.notes.push(format!(
                        "Cys23 ambiguous ({} other C in window); CDR1 start from the best-scoring one",
                        hit.cys23_alternatives
                    ));
                    annotation.uncertain[1] = true;
                }
                merge_regions(
                    &given.regions,
                    regions_from_hit(&hit, chain),
                    Some(anchor),
                    len_of(v, index),
                    annotation,
                )
            }
            _ => given.known(),
        };
        anchors.push(anchor);
        regions.push(merged);
    }

    // V pass 3: an ambiguous Cys23 in an allele the positional prior could
    // not reach (5′-truncated, other frame) is taken from the most similar
    // unambiguous allele, if it lands on a Cys codon in this allele's frame.
    let mut relative_index = None;
    let unambiguous: Vec<bool> = (0..v.len())
        .map(|i| {
            result[i].source.trusted()
                && regions[i].cdr1.is_some()
                && hits[i].is_some_and(|h| h.cys23_alternatives == 0)
        })
        .collect();
    let any_unambiguous = unambiguous.iter().any(|&u| u);
    for index in 0..v.len() {
        if !any_unambiguous {
            break;
        }
        let Some(hit) =
            hits[index].filter(|h| h.cys23_alternatives > 0 && result[index].source.trusted())
        else {
            continue;
        };
        if Some(hit.frame as u32) == frame_mode || given[index].regions[1] != RegionSlot::Unknown {
            continue;
        }
        let relatives = relatives(&mut relative_index, v);
        let Some(Relative { donor, mapped, .. }) = from_relative(
            v,
            index,
            &unambiguous,
            relatives,
            |d| vec![regions[d].cdr1],
            &mut workspace,
        ) else {
            continue;
        };
        let codes = &v.alleles()[index].codes;
        let lands_on_cys = mapped[0]
            .filter(|&p| p >= 12 && p % 3 == hit.frame as u32)
            .filter(|&p| residue_at(codes, p - 12) == b'C');
        if let Some(cdr1) = lands_on_cys {
            let mut values = regions[index].values();
            values[1] = Some(cdr1);
            if RegionStarts::from_values(values)
                .check(anchors[index], codes.len())
                .is_ok()
            {
                regions[index] = RegionStarts::from_values(values);
                let annotation = &mut result[index];
                annotation
                    .notes
                    .retain(|n| !n.starts_with("Cys23 ambiguous"));
                annotation.notes.push(format!(
                    "ambiguous Cys23 resolved from {}",
                    v.alleles()[donor].name
                ));
                annotation.uncertain[1] = false;
            }
        }
    }

    // V pass 4: sibling transfer of anchor and all region starts.
    let trusted: Vec<bool> = result
        .iter()
        .zip(&anchors)
        .map(|(r, a)| r.source.trusted() && a.is_some())
        .collect();
    let any_donor = trusted.iter().any(|&t| t);
    for index in 0..v.len() {
        if anchors[index].is_some() || given[index].ends_before_cys104 {
            continue;
        }
        let positions_of = |d: usize| {
            let mut p = vec![anchors[d]];
            p.extend(regions[d].values());
            p
        };
        let donor = if any_donor {
            let relatives = relatives(&mut relative_index, v);
            from_relative(v, index, &trusted, relatives, positions_of, &mut workspace)
        } else {
            None
        };
        let len = len_of(v, index);
        let partial = |r: &Relative| anchors[r.donor].is_some_and(|a| r.runs_out_before(a, len));
        match donor {
            Some(Relative { donor, mapped, .. })
                if mapped[0].is_some_and(|a| a as usize + 3 <= len) =>
            {
                let anchor = mapped[0].expect("checked above");
                let derived = RegionStarts::from_values([
                    mapped[1], mapped[2], mapped[3], mapped[4], mapped[5],
                ]);
                let annotation = &mut result[index];
                let residue = residue_at(&v.alleles()[index].codes, anchor);
                if residue != b'C' {
                    annotation
                        .notes
                        .push(format!("conserved residue is {}", char::from(residue)));
                }
                regions[index] = merge_regions(
                    &given[index].regions,
                    derived,
                    Some(anchor),
                    len,
                    annotation,
                );
                result[index].source = AnchorSource::Transfer {
                    from: v.alleles()[donor].name.clone(),
                };
                anchors[index] = Some(anchor);
            }
            // The closest relative's anchor codon maps at or past the 3′ end
            // and the allele ends while still aligned to it: a 3′-partial
            // allele, which holds no Cys104. A weak motif hit would be some
            // other codon, so the anchor stays missing; the relative's region
            // starts inside the allele are kept.
            Some(r) if partial(&r) => {
                let derived = RegionStarts::from_values([
                    r.mapped[1],
                    r.mapped[2],
                    r.mapped[3],
                    r.mapped[4],
                    r.mapped[5],
                ]);
                let annotation = &mut result[index];
                annotation.notes.push(format!(
                    "3′-partial: the anchor codon of {} maps past the 3′ end",
                    v.alleles()[r.donor].name
                ));
                regions[index] =
                    merge_regions(&given[index].regions, derived, None, len, annotation);
            }
            relative => {
                if let Some(hit) = hits[index] {
                    let annotation = &mut result[index];
                    annotation.notes.push(match relative {
                        Some(r) => format!(
                            "weak motifs; the anchor of {} does not fit and their alignment stops short",
                            v.alleles()[r.donor].name
                        ),
                        None => "weak motifs; no similar annotated allele".into(),
                    });
                    regions[index] = merge_regions(
                        &given[index].regions,
                        regions_from_hit(&hit, chain),
                        Some(hit.anchor()),
                        len,
                        annotation,
                    );
                    result[index].source = AnchorSource::LowConfidence;
                    anchors[index] = Some(hit.anchor());
                }
            }
        }
    }
    // Every derived start of an allele without a trusted anchor is uncertain.
    for ((annotation, given), regions) in result.iter_mut().zip(&given).zip(&regions) {
        if annotation.source.trusted() || matches!(annotation.source, AnchorSource::Transfer { .. })
        {
            continue;
        }
        let derived = given.regions.iter().zip(regions.values());
        for (flag, (slot, value)) in annotation.uncertain.iter_mut().zip(derived) {
            *flag |= value.is_some() && slot.known().is_none();
        }
    }
    for ((_, anchor, slot), (found, merged)) in
        v.annotate_mut().zip(anchors.into_iter().zip(regions))
    {
        *anchor = anchor.or(found);
        *slot = merged;
    }
    DatabaseAnnotation {
        chain,
        v_tail,
        j_tail,
        v: result,
        j: j_annotations,
        imgt_gaps,
    }
}

impl DatabaseAnnotation {
    fn counts(list: &[AlleleAnnotation]) -> String {
        let mut counts: Vec<(String, usize)> = Vec::new();
        for annotation in list {
            let label = match &annotation.source {
                AnchorSource::Transfer { .. } => "transfer".to_string(),
                other => other.label(),
            };
            match counts.iter_mut().find(|(l, _)| *l == label) {
                Some((_, count)) => *count += 1,
                None => counts.push((label, 1)),
            }
        }
        counts
            .iter()
            .map(|(l, c)| format!("{c} {l}"))
            .collect::<Vec<_>>()
            .join(", ")
    }

    /// One line per segment for logs.
    pub fn summary(&self) -> String {
        let chain = self.chain.map_or("unknown", ChainType::name);
        let mut text = format!("chain {chain}; V anchors: {}", Self::counts(&self.v));
        if !self.j.is_empty() {
            text += &format!("; J anchors: {}", Self::counts(&self.j));
        }
        match &self.imgt_gaps {
            ImgtGaps::Absent => {}
            ImgtGaps::Used { gapped } => {
                text += &format!("; IMGT numbering read from {gapped} gapped V alleles");
            }
            ImgtGaps::Ignored { gapped, reason } => {
                text +=
                    &format!("; gaps of {gapped} V alleles ignored, not IMGT numbering ({reason})");
            }
        }
        text
    }

    /// Alleles whose anchor is not from a header, IMGT gaps, motif or consensus.
    pub fn uncertain(&self) -> usize {
        self.v
            .iter()
            .chain(&self.j)
            .filter(|a| !a.source.trusted())
            .count()
    }

    /// Per-allele report (TSV).
    pub fn write_tsv<W: Write>(
        &self,
        out: &mut W,
        v: &GermlineSet,
        j: Option<&GermlineSet>,
    ) -> io::Result<()> {
        writeln!(
            out,
            "segment\tallele\taliases\tlength\tframe\tanchor\tresidue\tfr1\tcdr1\tfr2\tcdr2\tfr3\tcdr1_codons\tcdr2_codons\tsource\tnotes"
        )?;
        let fmt = |x: Option<u32>| x.map_or(String::new(), |x| x.to_string());
        let rows = std::iter::once((v, &self.v, "V")).chain(j.map(|j| (j, &self.j, "J")));
        for (set, annotations, segment) in rows {
            for (allele, annotation) in set.alleles().iter().zip(annotations) {
                let r = &allele.regions;
                let span = |a: Option<u32>, b: Option<u32>| {
                    a.zip(b)
                        .map_or(String::new(), |(a, b)| ((b - a) / 3).to_string())
                };
                let residue = allele.anchor.map_or(String::new(), |a| {
                    char::from(residue_at(&allele.codes, a)).to_string()
                });
                writeln!(
                    out,
                    "{segment}\t{}\t{}\t{}\t{}\t{}\t{residue}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    allele.name,
                    allele.aliases.join(","),
                    allele.codes.len(),
                    allele.anchor.map_or(String::new(), |a| (a % 3).to_string()),
                    fmt(allele.anchor),
                    fmt(r.fr1),
                    fmt(r.cdr1),
                    fmt(r.fr2),
                    fmt(r.cdr2),
                    fmt(r.fr3),
                    span(r.cdr1, r.fr2),
                    span(r.cdr2, r.fr3),
                    annotation.source.label(),
                    annotation.notes_field(),
                )?;
            }
        }
        Ok(())
    }
}
