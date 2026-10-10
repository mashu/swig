//! Germline database annotation: anchors, reading frames and V region boundaries.
//!
//! V and J germline genes are in frame with their conserved anchors (V Cys104,
//! J Trp/Phe118), so an anchor fixes the reading frame and, for V, the region
//! boundaries follow from three conserved residues (IMGT definitions):
//!
//! ```text
//! FR1  = [frame start, Cys23 + 4 codons)      CDR1 = [Cys23 + 4, Trp41 − 2)
//! FR2  = [Trp41 − 2, Trp41 + 15)              CDR2 = [Trp41 + 15, FR3 start)
//! FR3  = [Cys104 − (L − 1), Cys104 + 1)       L = FR3 codons (below)
//! CDR3 = [Cys104 + 1, J Trp/Phe118)           FR4  = [J Trp/Phe118, J end)
//! ```
//!
//! The CDR2/FR3 boundary has no conserved residue. An FR3 start that is not
//! given comes from the nearest relative whose FR3 start is given, else from
//! the FR3 length the database's given starts imply for the allele's
//! Trp41 → Cys104 spacing ([`Fr3Model`]), else from the chain default
//! (38 codons heavy, 36 light; wrong for some TR and IGL families).
//! Likewise, when the Cys23 window holds several `C` codons or none, CDR1
//! starts where the nearest relative with a settled Cys23 puts it, else
//! where the database's settled alleles say Cys23 sits ([`Cys23Model`]).
//! Every start that may be wrong is marked in
//! [`AlleleAnnotation::uncertain`].
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
//! The chain type (which fixes the default FR3 length) is voted by J residues
//! (`W` heavy, `F` light) and Trp41 contexts (`W.R` heavy, `W.Q` light).
//! Everything runs once per database: one translation per allele plus an
//! alignment for the rare transferred allele — milliseconds for thousands of
//! alleles, and nothing is added to per-read work.

mod cys23;
mod fr3;
mod given;
mod j_motif;
mod transfer;
mod v_motif;

use std::collections::HashMap;
use std::io::{self, Write};

pub use cys23::Cys23Model;
pub use fr3::Fr3Model;
pub use given::ImgtGaps;
pub use j_motif::{JHit, find_j};
pub use transfer::{
    MIN_SHARED_FRACTION, MIN_SHARED_KMERS, RelativeIndex, kmers, shared, transfer, transfer_aligned,
};
pub use v_motif::{VHit, VSearch, find_v};

use crate::align::GotohWorkspace;
use crate::germline::{GermlineSet, RegionSlot, RegionStarts};
use crate::translate::{codon, translate};
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
    /// for an allele whose anchor is low-confidence or missing; a CDR1
    /// start placed without a Cys23 codon, or from an ambiguous Cys23 that
    /// no relative settled; a derived bound of a CDR of unusual length; an
    /// FR3 length of model posterior below 0.95 or outside the model, or
    /// the chain default under [`DatabaseAnnotation::fr3_warning`]; a start
    /// transferred from a relative whose own start is uncertain. Given
    /// values are never uncertain. An error in a start not flagged here is
    /// a silent error.
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
    /// FR3 lengths learned from the alleles whose FR3 start is given.
    pub fr3_model: Option<Fr3Model>,
    /// Most V allele names are TR or IGL ones: loci whose FR3 length IMGT
    /// varies by gene family, so one chain default is wrong for some.
    pub fr3_varies: bool,
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

/// A probability for a note, rounded down to two decimals so that only a
/// certain choice reads `1.00`.
fn probability(p: f64) -> String {
    format!("{:.2}", (p * 100.0).floor() / 100.0)
}

fn residue_at(codes: &[u8], anchor: u32) -> u8 {
    let a = anchor as usize;
    codes
        .get(a..a + 3)
        .map_or(b'?', |c| codon(c[0], c[1], c[2]))
}

/// Region starts implied by a V hit (IMGT definitions, see module docs),
/// with the FR3 start `fr3`.
fn regions_from_hit(hit: &VHit, fr3: Fr3Start) -> RegionStarts {
    let fr3 = match fr3 {
        Fr3Start::Codons(codons) => hit.anchor().checked_sub(3 * (codons - 1)),
        Fr3Start::At(position) => Some(position),
        Fr3Start::Unknown => None,
    };
    RegionStarts {
        fr1: Some(hit.frame as u32),
        cdr1: hit.cys23_nt().map(|c| c + 12),
        fr2: hit.trp41_nt().and_then(|w| w.checked_sub(6)),
        cdr2: hit.trp41_nt().map(|w| w + 45),
        fr3,
    }
}

/// A derived FR3 start: an FR3 length before Cys104, or a position.
#[derive(Clone, Copy, Debug)]
enum Fr3Start {
    Codons(u32),
    At(u32),
    Unknown,
}

/// Where the FR3 start of an allele without a given one came from.
#[derive(Clone, Copy, Debug)]
enum Fr3Choice {
    /// The nearest relative with a given FR3 start (allele index).
    Relative(usize),
    /// The database model's most probable length and its posterior.
    Model { codons: u32, p: f64 },
    /// The model admits no FR3 length for the allele's spacing.
    OutsideModel,
    /// No model: the chain default.
    ChainDefault,
}

/// The FR3 length for `hit` when its FR3 start is not given: the database
/// model's most probable length for the Trp41 → Cys104 spacing, else the
/// chain default.
fn fr3_length(
    hit: &VHit,
    chain: Option<ChainType>,
    model: Option<&Fr3Model>,
) -> (Fr3Start, Fr3Choice) {
    let default = chain.map_or(Fr3Start::Unknown, |c| Fr3Start::Codons(c.fr3_codons()));
    let Some(model) = model else {
        return (default, Fr3Choice::ChainDefault);
    };
    let spacing = hit.trp41.and_then(|w| hit.cys104.checked_sub(w));
    match model.length_for(spacing.map(|d| d as u32)) {
        Some((codons, p)) => (Fr3Start::Codons(codons), Fr3Choice::Model { codons, p }),
        None => (default, Fr3Choice::OutsideModel),
    }
}

/// [`merge_regions`] with the regions of `hit` and the first FR3 start the
/// merge keeps of `relative` and [`fr3_length`]. Returns the merged starts
/// and, when the FR3 start was not given but derived and kept, its source.
#[allow(clippy::too_many_arguments)]
fn merge_hit(
    given: &[RegionSlot; 5],
    hit: &VHit,
    relative: Option<(u32, usize)>,
    chain: Option<ChainType>,
    model: Option<&Fr3Model>,
    anchor: Option<u32>,
    len: usize,
    annotation: &mut AlleleAnnotation,
) -> (RegionStarts, Option<Fr3Choice>) {
    let open = given[4] == RegionSlot::Unknown;
    if let Some((fr3, donor)) = relative.filter(|_| open) {
        let mut trial = annotation.clone();
        let merged = merge_regions(
            given,
            regions_from_hit(hit, Fr3Start::At(fr3)),
            anchor,
            len,
            &mut trial,
        );
        if merged.fr3 == Some(fr3) {
            *annotation = trial;
            return (merged, Some(Fr3Choice::Relative(donor)));
        }
    }
    let (fr3, choice) = fr3_length(hit, chain, model);
    let merged = merge_regions(given, regions_from_hit(hit, fr3), anchor, len, annotation);
    (merged, (open && merged.fr3.is_some()).then_some(choice))
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

/// A start transferred from a relative is as uncertain as the relative's
/// own (a CDR1 from an ambiguous Cys23, an FR3 from a doubtful length);
/// given starts never are.
fn inherit_uncertain(annotation: &mut AlleleAnnotation, given: &[RegionSlot; 5], donor: [bool; 5]) {
    for ((flag, slot), inherited) in annotation.uncertain.iter_mut().zip(given).zip(donor) {
        *flag |= inherited && slot.known().is_none();
    }
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

/// Least posterior of a model FR3 length that is not marked uncertain. The
/// model assumes CDR2 and FR3 lengths independent across the labelled
/// alleles, which fails for a family with no label: human TRGV1–8 get
/// posteriors of 0.58–0.92 for a wrong length when only TRGV9–11 are
/// annotated. On 13 loci with 10 % of alleles annotated, one model value in
/// nine with a posterior of 0.8–0.95 was wrong, and none above.
const FR3_CONFIDENT: f64 = 0.95;

/// Writes the provenance note and uncertainty of a derived FR3 start, after
/// the merge and only for the value kept.
struct Fr3Notes<'a> {
    set: &'a GermlineSet,
    model: Option<&'a Fr3Model>,
    /// The chain default is used in a locus whose FR3 length varies by family.
    default_doubtful: bool,
}

impl Fr3Notes<'_> {
    fn record(&self, annotation: &mut AlleleAnnotation, choice: Option<Fr3Choice>) {
        let Some(choice) = choice else { return };
        let (note, uncertain) = match choice {
            Fr3Choice::Relative(donor) => (
                Some(format!("FR3 start from {}", self.set.alleles()[donor].name)),
                false,
            ),
            Fr3Choice::Model { codons, p } => (
                Some(format!(
                    "FR3 length {codons} codons from the database model (n={}, p={})",
                    self.model.map_or(0, Fr3Model::labels),
                    probability(p)
                )),
                p < FR3_CONFIDENT,
            ),
            Fr3Choice::OutsideModel => (
                Some("FR3 length outside the database model; chain default".into()),
                true,
            ),
            Fr3Choice::ChainDefault => (None, self.default_doubtful),
        };
        annotation.notes.extend(note);
        annotation.uncertain[4] |= uncertain;
    }
}

/// For every trusted allele (`hits`) whose FR3 start is open, the FR3 start
/// of its nearest relative with a given FR3 start and an anchor, and that
/// relative: the relative's Cys104 → FR3 distance when the two share a
/// Trp41 → Cys104 spacing, else its FR3 start mapped by alignment. The FR3
/// start is a homologous position, so an indel between the two falls in
/// the region their alignment puts it in. Kept when it lies in frame before
/// the anchor and makes an FR3 of IMGT length.
fn fr3_from_relatives(
    set: &GermlineSet,
    given: &[Given],
    anchors: &[Option<u32>],
    hits: &[Option<VHit>],
    index: &mut Option<RelativeIndex>,
    workspace: &mut GotohWorkspace,
) -> Vec<Option<(u32, usize)>> {
    let fr3_of = |i: usize| given[i].regions[4].known().filter(|_| anchors[i].is_some());
    let eligible: Vec<bool> = (0..set.len()).map(|i| fr3_of(i).is_some()).collect();
    let mut found = vec![None; set.len()];
    if !eligible.contains(&true) {
        return found;
    }
    let spacing = |i: usize| {
        let trp41 = hits[i]?.trp41_nt()?;
        anchors[i]?.checked_sub(trp41)
    };
    for (target, slot) in found.iter_mut().enumerate() {
        let Some(anchor) = anchors[target].filter(|_| hits[target].is_some()) else {
            continue;
        };
        if given[target].regions[4] != RegionSlot::Unknown {
            continue;
        }
        let Some(donor) = relatives(index, set).nearest(target, &eligible) else {
            continue;
        };
        let donor_fr3 = fr3_of(donor).expect("eligible donor");
        let donor_anchor = anchors[donor].expect("eligible donor");
        let fr3 = if spacing(target).is_some() && spacing(target) == spacing(donor) {
            donor_anchor
                .checked_sub(donor_fr3)
                .and_then(|d| anchor.checked_sub(d))
        } else {
            let alleles = set.alleles();
            transfer(
                &alleles[target].codes,
                &alleles[donor].codes,
                &[Some(donor_fr3)],
                workspace,
            )
            .and_then(|mapped| mapped[0])
        };
        *slot = fr3
            .filter(|&f| {
                anchor
                    .checked_sub(f)
                    .is_some_and(|d| d.is_multiple_of(3) && Fr3Model::admits(d / 3 + 1))
            })
            .map(|f| (f, donor));
    }
    found
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

    // FR3 and CDR2 lengths of the alleles whose FR3 start is given.
    let fr3_varies = 2 * v
        .alleles()
        .iter()
        .filter(|a| a.name.starts_with("TR") || a.name.starts_with("IGL"))
        .count()
        > v.len();
    let fr3_model = Fr3Model::learn((0..v.len()).filter_map(|i| {
        let fr3 = given[i].regions[4].known()?;
        let hit = hits[i].filter(|_| result[i].source.trusted());
        let anchor = given[i].anchor.or(hit.map(|h| h.anchor()))?;
        let codons = anchor.checked_sub(fr3).filter(|d| d.is_multiple_of(3))? / 3 + 1;
        let cdr2 = given[i].regions[3]
            .known()
            .or(hit.and_then(|h| h.trp41_nt()).map(|w| w + 45));
        let cdr2_codons = cdr2
            .and_then(|c| fr3.checked_sub(c))
            .filter(|d| d.is_multiple_of(3))
            .map(|d| d / 3);
        Some((codons, cdr2_codons))
    }));
    // The Trp41 search covers every FR3 length the chain or the model admits.
    let fr3_range = match (chain.map(ChainType::fr3_codons), fr3_model.as_ref()) {
        (Some(l), Some(model)) => {
            let (lo, hi) = model.observed();
            Some((l.min(lo), l.max(hi)))
        }
        (Some(l), None) => Some((l, l)),
        (None, model) => model.map(Fr3Model::observed),
    };

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
            fr3_codons: fr3_range,
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
    let trusted_hits: Vec<Option<VHit>> = (0..v.len())
        .map(|i| hits[i].filter(|_| result[i].source.trusted()))
        .collect();
    let mut anchors: Vec<Option<u32>> = given
        .iter()
        .zip(&trusted_hits)
        .map(|(given, hit)| given.anchor.or(hit.map(|h| h.anchor())))
        .collect();
    let mut relative_index = None;
    let fr3_relatives = fr3_from_relatives(
        v,
        &given,
        &anchors,
        &trusted_hits,
        &mut relative_index,
        &mut workspace,
    );
    let fr3_notes = Fr3Notes {
        set: v,
        model: fr3_model.as_ref(),
        default_doubtful: fr3_model.is_none() && fr3_varies,
    };
    let mut regions: Vec<RegionStarts> = Vec::with_capacity(v.len());
    for (index, given) in given.iter().enumerate() {
        let trusted_hit = trusted_hits[index];
        let anchor = anchors[index];
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
                let (merged, choice) = merge_hit(
                    &given.regions,
                    &hit,
                    fr3_relatives[index],
                    chain,
                    fr3_model.as_ref(),
                    Some(anchor),
                    len_of(v, index),
                    annotation,
                );
                fr3_notes.record(annotation, choice);
                merged
            }
            _ => given.known(),
        };
        regions.push(merged);
    }

    // V pass 3: Cys23 where the motif search could not settle it (several
    // `C` codons in the window, or none) and no header or gap column gives
    // CDR1: the nearest relative with a settled Cys23, else the database's
    // Cys23 model.
    let open_cys23: Vec<Option<VHit>> = (0..v.len())
        .map(|i| {
            trusted_hits[i]
                .filter(|h| h.trp41.is_some() && (h.cys23_alternatives > 0 || h.cys23.is_none()))
                .filter(|_| given[i].regions[1] == RegionSlot::Unknown)
        })
        .collect();
    if open_cys23.iter().any(Option::is_some) {
        let settled: Vec<bool> = (0..v.len())
            .map(|i| {
                regions[i].cdr1.is_some()
                    && trusted_hits[i].is_some()
                    && (given[i].regions[1].known().is_some() || open_cys23[i].is_none())
            })
            .collect();
        let translations: Vec<Option<Vec<u8>>> = trusted_hits
            .iter()
            .enumerate()
            .map(|(i, h)| h.map(|h| translate(&v.alleles()[i].codes, h.frame)))
            .collect();
        let model = Cys23Model::learn((0..v.len()).filter(|&i| settled[i]).filter_map(|i| {
            let h = trusted_hits[i]?;
            let cys23 = regions[i].cdr1?.checked_sub(12 + h.frame as u32)?;
            Some(cys23::Label {
                aa: translations[i].as_deref()?,
                frame: h.frame,
                cys23: cys23.is_multiple_of(3).then_some(cys23 as usize / 3)?,
                trp41: h.trp41?,
            })
        }));
        let any_settled = settled.contains(&true);
        for index in 0..v.len() {
            let Some(hit) = open_cys23[index] else {
                continue;
            };
            let missing = hit.cys23.is_none();
            let codes = &v.alleles()[index].codes;
            let set_cdr1 = |regions: &mut RegionStarts, cdr1: u32| {
                let mut values = regions.values();
                values[1] = Some(cdr1);
                let candidate = RegionStarts::from_values(values);
                let fits = candidate.check(anchors[index], codes.len()).is_ok();
                if fits {
                    *regions = candidate;
                }
                fits
            };
            // 1. The nearest relative with a settled Cys23, in any frame. With
            // no Cys in the window, its CDR1 start is taken as mapped (it may
            // precede position 12: a 5′-partial allele that starts after Cys23).
            let relative = if any_settled {
                from_relative(
                    v,
                    index,
                    &settled,
                    relatives(&mut relative_index, v),
                    |d| vec![regions[d].cdr1],
                    &mut workspace,
                )
            } else {
                None
            };
            if let Some(Relative { donor, mapped, .. }) = relative {
                let lands = mapped[0]
                    .filter(|&p| p % 3 == hit.frame as u32)
                    .filter(|&p| missing || (p >= 12 && residue_at(codes, p - 12) == b'C'));
                if let Some(cdr1) = lands
                    && set_cdr1(&mut regions[index], cdr1)
                {
                    let name = &v.alleles()[donor].name;
                    let annotation = &mut result[index];
                    if missing {
                        annotation
                            .notes
                            .push(format!("no Cys23; CDR1 start from {name}"));
                        annotation.uncertain[1] = true;
                    } else {
                        annotation
                            .notes
                            .retain(|n| !n.starts_with("Cys23 ambiguous"));
                        annotation
                            .notes
                            .push(format!("ambiguous Cys23 resolved from {name}"));
                        annotation.uncertain[1] = false;
                    }
                    continue;
                }
            }
            // 2. The database's Cys23 model, over the `C` codons of the
            // window; a mutated Cys23 is placed only at an offset where a
            // labelled allele has its Cys23 (a full-length 5′ end), and only
            // when it takes half the posterior.
            let (Some(model), Some(aa), Some(trp41)) =
                (&model, translations[index].as_deref(), hit.trp41)
            else {
                continue;
            };
            let window = trp41.saturating_sub(20)..=(trp41.saturating_sub(10)).min(aa.len() - 1);
            let candidates = window
                .filter(|&k| aa[k] == b'C' || (missing && model.seen_offset(3 * k + hit.frame)));
            let Some((k, p)) = model.choose(aa, hit.frame, trp41, candidates) else {
                continue;
            };
            if (missing && p < 0.5)
                || !set_cdr1(&mut regions[index], (3 * k + hit.frame) as u32 + 12)
            {
                continue;
            }
            let annotation = &mut result[index];
            annotation.uncertain[1] = true;
            if missing {
                annotation.notes.push(format!(
                    "no Cys23 ({} at its position); CDR1 start from the database model (n={}, p={})",
                    char::from(aa[k]),
                    model.labels(),
                    probability(p)
                ));
            } else {
                for note in annotation
                    .notes
                    .iter_mut()
                    .filter(|n| n.starts_with("Cys23 ambiguous"))
                {
                    *note = format!(
                        "Cys23 ambiguous ({} other C in window); CDR1 start from the database model (n={}, p={})",
                        hit.cys23_alternatives,
                        model.labels(),
                        probability(p)
                    );
                }
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
                let donor_uncertain = result[donor].uncertain;
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
                inherit_uncertain(annotation, &given[index].regions, donor_uncertain);
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
                let donor_uncertain = result[r.donor].uncertain;
                let annotation = &mut result[index];
                annotation.notes.push(format!(
                    "3′-partial: the anchor codon of {} maps past the 3′ end",
                    v.alleles()[r.donor].name
                ));
                regions[index] =
                    merge_regions(&given[index].regions, derived, None, len, annotation);
                inherit_uncertain(annotation, &given[index].regions, donor_uncertain);
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
                    let (merged, choice) = merge_hit(
                        &given[index].regions,
                        &hit,
                        None,
                        chain,
                        fr3_model.as_ref(),
                        Some(hit.anchor()),
                        len,
                        annotation,
                    );
                    regions[index] = merged;
                    fr3_notes.record(annotation, choice);
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
        fr3_model,
        fr3_varies,
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
        if let Some(model) = &self.fr3_model {
            text += &format!("; FR3 lengths learned from {} alleles", model.labels());
        }
        text
    }

    /// No V allele has a given FR3 start, and the names are those of a
    /// locus (TR, IGL) whose FR3 length varies by gene family: one chain
    /// default places some FR3 starts wrong, and every such start is
    /// marked uncertain.
    pub fn fr3_warning(&self) -> bool {
        self.fr3_model.is_none() && self.fr3_varies && self.chain.is_some()
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
