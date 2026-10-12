//! What the input states about each V allele before any search: header
//! keys, then the IMGT numbering of gapped records.
//!
//! An IMGT-gapped V set carries the IMGT unique numbering in its gap
//! columns, so every region start and Cys104 can be read off them
//! ([`ImgtNumbering`]). Any codon alignment written with `.` would parse the
//! same way, so the numbering is used only when the database looks like
//! IMGT's: a majority of gapped alleles is codon-aligned, and at each of
//! IMGT codons 23, 41 and 104 a majority of the gapped alleles holding the
//! codon have the conserved Cys, Trp, Cys there (IMGT human IGH: 97.8, 95.3
//! and 98.5 %). Another alignment's columns put other residues there.
//!
//! Header keys win over gap columns; if their union fails
//! [`RegionStarts::check`], the gaps of that allele are ignored with a note.

use crate::germline::{
    Allele, GermlineSet, IMGT_CONSERVED_CODONS, ImgtNumbering, RegionSlot, RegionStarts,
};

use super::{AnchorSource, Note};

/// The conserved residues at [`IMGT_CONSERVED_CODONS`].
const CONSERVED: [u8; 3] = *b"CWC";

/// How the V set's IMGT gap columns were used.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum ImgtGaps {
    /// No V record is gapped.
    #[default]
    Absent,
    /// The numbering of every gapped allele was read.
    Used { gapped: usize },
    /// The gaps do not follow the IMGT numbering and were ignored.
    Ignored { gapped: usize, reason: String },
}

impl ImgtGaps {
    /// Check the gapped alleles of `v` against the IMGT format (module docs).
    pub fn check(v: &GermlineSet) -> Self {
        let numbered: Vec<&ImgtNumbering> = v
            .alleles()
            .iter()
            .filter_map(|a| a.imgt.as_deref())
            .collect();
        let gapped = numbered.len();
        if gapped == 0 {
            return Self::Absent;
        }
        let aligned = numbered.iter().filter(|n| n.codon_aligned).count();
        if 2 * aligned <= gapped {
            return Self::Ignored {
                gapped,
                reason: format!("{aligned} of {gapped} gapped V alleles are codon-aligned"),
            };
        }
        for (index, (codon, residue)) in IMGT_CONSERVED_CODONS.iter().zip(CONSERVED).enumerate() {
            let held = numbered.iter().filter_map(|n| n.conserved[index]);
            let (total, conserved) =
                held.fold((0, 0), |(t, c), r| (t + 1, c + usize::from(r == residue)));
            if 2 * conserved <= total {
                return Self::Ignored {
                    gapped,
                    reason: format!(
                        "{conserved} of {total} gapped V alleles have {} at IMGT {codon}",
                        char::from(residue)
                    ),
                };
            }
        }
        Self::Used { gapped }
    }

    pub fn used(&self) -> bool {
        matches!(self, Self::Used { .. })
    }
}

/// The given annotation of one V allele.
#[derive(Clone, Debug)]
pub struct Given {
    pub anchor: Option<u32>,
    /// [`AnchorSource::Header`] or [`AnchorSource::ImgtGaps`] when the
    /// anchor is given.
    pub source: Option<AnchorSource>,
    /// The header anchor is `SWIGMETA`'s ([`Allele::anchor_from_meta`]).
    pub from_meta: bool,
    pub regions: [RegionSlot; 5],
    /// The IMGT numbering ends before Cys104: the allele holds no anchor,
    /// and none may be searched for or transferred.
    pub ends_before_cys104: bool,
}

impl Given {
    /// Header keys, then (with `use_gaps`) the allele's IMGT numbering in
    /// the slots and anchor the header leaves open. Notes go to `notes`.
    pub fn read(allele: &Allele, use_gaps: bool, notes: &mut Vec<Note>) -> Self {
        let header = Self {
            anchor: allele.anchor,
            source: allele.anchor.map(|_| AnchorSource::Header),
            from_meta: allele.anchor_from_meta,
            regions: allele
                .regions
                .values()
                .map(|v| v.map_or(RegionSlot::Unknown, RegionSlot::Known)),
            ends_before_cys104: false,
        };
        let Some(imgt) = allele.imgt.as_deref().filter(|_| use_gaps) else {
            return header;
        };
        let mut merged = header.clone();
        for (slot, numbered) in merged.regions.iter_mut().zip(imgt.slots) {
            if *slot == RegionSlot::Unknown {
                *slot = numbered;
            }
        }
        if merged.anchor.is_none() {
            merged.anchor = imgt.anchor;
            merged.source = imgt.anchor.map(|_| AnchorSource::ImgtGaps);
            merged.ends_before_cys104 = imgt.ends_before_cys104;
        }
        if let Err(reason) = merged.known().check(merged.anchor, allele.codes.len()) {
            notes.push(Note::ImgtGapsIgnored { reason });
            return header;
        }
        notes.extend(
            imgt.frameshift_markers
                .iter()
                .map(|column| Note::ImgtFrameshiftMarker {
                    codon: column / 3 + 1,
                }),
        );
        if merged.ends_before_cys104 {
            notes.push(Note::ImgtEndsBeforeCys104);
        }
        merged
    }

    /// The known region starts.
    pub fn known(&self) -> RegionStarts {
        RegionStarts::from_values(self.regions.map(RegionSlot::known))
    }
}
