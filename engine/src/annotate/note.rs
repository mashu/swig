//! What annotation reports about one allele: the notes column of
//! `annotation.tsv` ([`super::AlleleAnnotation::notes_field`]).
//!
//! Notes are values, so the passes that revise an earlier finding (a
//! relative that settles an ambiguous Cys23) and the tests match them by
//! variant, never by their wording.

use std::fmt;

use crate::germline::RegionStarts;

/// A posterior for a note, rounded down to two decimals so that only a
/// certain choice reads `1.00`.
fn probability(p: f64) -> String {
    format!("{:.2}", (p * 100.0).floor() / 100.0)
}

/// Where the CDR1 start of an allele with several `C` codons in its Cys23
/// window came from.
#[derive(Clone, Debug, PartialEq)]
pub enum Cys23Choice {
    /// The best-scoring `C` of the window.
    BestScoring,
    /// The database's Cys23 model ([`super::Cys23Model`]) of `labels`
    /// alleles, with posterior `p`.
    Model { labels: u32, p: f64 },
}

#[derive(Clone, Debug, PartialEq)]
pub enum Note {
    /// A header (key or `SWIGMETA`) anchor on a codon that is not the
    /// conserved residue (V `C`, J `W`/`F`); the codon's residue.
    HeaderResidue(u8),
    /// IMGT codon 104 of a gapped record is not `C`.
    ImgtResidue(u8),
    /// A transferred anchor on a codon that is not the conserved residue.
    TransferredResidue(u8),
    /// J: the relative's anchor codon maps outside this partial allele.
    PartialJ { relative: String },
    /// J: the anchor of a `[WF]G.G` motif that is not confident.
    WeakJMotif { residue: u8 },
    /// V: motifs confident on their own but not under the database's
    /// priors (second pass); the anchor is dropped.
    InconsistentMotifs,
    /// V: several `C` codons in the Cys23 window; CDR1 from `choice`.
    Cys23Ambiguous { others: usize, choice: Cys23Choice },
    /// V: an ambiguous Cys23 settled by the nearest relative whose Cys23 is.
    Cys23Resolved { relative: String },
    /// V: no `C` in the Cys23 window; CDR1 from the nearest relative.
    NoCys23 { relative: String },
    /// V: no `C` in the Cys23 window; CDR1 from the database's Cys23 model
    /// at a codon of residue `residue`.
    NoCys23Model { residue: u8, labels: u32, p: f64 },
    /// V: the relative's anchor codon maps past the 3′ end of this allele,
    /// which holds no Cys104.
    ThreePrimePartial { relative: String },
    /// V: the IMGT numbering ends before Cys104.
    ImgtEndsBeforeCys104,
    /// V: weak motifs used as a low-confidence anchor; `relative` is the
    /// nearest relative when its anchor did not fit.
    WeakMotifs { relative: Option<String> },
    /// V: the given region starts break [`RegionStarts::check`].
    HeaderRegionsInconsistent { reason: String },
    /// V: a derived start (index into [`RegionStarts::KEYS`]) dropped to keep
    /// the starts ordered.
    DroppedDerived { region: usize, reason: String },
    /// V: CDR1 (`cdr` 1) or CDR2 (2) of an unusual length.
    UnusualCdrLength { cdr: u8, codons: i64 },
    /// V: FR3 start from the nearest relative with a given one.
    Fr3FromRelative { relative: String },
    /// V: FR3 length from the database's FR3 model of `labels` alleles.
    Fr3FromModel { codons: u32, labels: u32, p: f64 },
    /// V: the FR3 model admits no length for the allele; chain default.
    Fr3OutsideModel,
    /// V: the record's IMGT gap numbering contradicts its header keys.
    ImgtGapsIgnored { reason: String },
    /// V: an IMGT frameshift marker (a one- or two-column gap) in `codon`.
    ImgtFrameshiftMarker { codon: u32 },
}

impl fmt::Display for Note {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let residue = |r: &u8| char::from(*r);
        match self {
            Self::HeaderResidue(r) => write!(f, "header anchor codon is {}", residue(r)),
            Self::ImgtResidue(r) => write!(f, "IMGT 104 codon is {}", residue(r)),
            Self::TransferredResidue(r) => write!(f, "conserved residue is {}", residue(r)),
            Self::PartialJ { relative } => write!(
                f,
                "partial: the anchor codon of {relative} maps outside the sequence"
            ),
            Self::WeakJMotif { residue: r } => write!(f, "weak motif {}G.G", residue(r)),
            Self::InconsistentMotifs => {
                write!(f, "motifs inconsistent with the database consensus")
            }
            Self::Cys23Ambiguous { others, choice } => {
                write!(f, "Cys23 ambiguous ({others} other C in window); ")?;
                match choice {
                    Cys23Choice::BestScoring => write!(f, "CDR1 start from the best-scoring one"),
                    Cys23Choice::Model { labels, p } => write!(
                        f,
                        "CDR1 start from the database model (n={labels}, p={})",
                        probability(*p)
                    ),
                }
            }
            Self::Cys23Resolved { relative } => {
                write!(f, "ambiguous Cys23 resolved from {relative}")
            }
            Self::NoCys23 { relative } => write!(f, "no Cys23; CDR1 start from {relative}"),
            Self::NoCys23Model {
                residue: r,
                labels,
                p,
            } => write!(
                f,
                "no Cys23 ({} at its position); CDR1 start from the database model (n={labels}, p={})",
                residue(r),
                probability(*p)
            ),
            Self::ThreePrimePartial { relative } => write!(
                f,
                "3′-partial: the anchor codon of {relative} maps past the 3′ end"
            ),
            Self::ImgtEndsBeforeCys104 => {
                write!(f, "3′-partial: the IMGT numbering ends before Cys104")
            }
            Self::WeakMotifs { relative } => match relative {
                Some(relative) => write!(
                    f,
                    "weak motifs; the anchor of {relative} does not fit and their alignment stops short"
                ),
                None => write!(f, "weak motifs; no similar annotated allele"),
            },
            Self::HeaderRegionsInconsistent { reason } => {
                write!(f, "header regions inconsistent: {reason}")
            }
            Self::DroppedDerived { region, reason } => {
                write!(
                    f,
                    "dropped derived {}: {reason}",
                    RegionStarts::KEYS[*region]
                )
            }
            Self::UnusualCdrLength { cdr, codons } => {
                write!(f, "unusual CDR{cdr} length {codons} codons")
            }
            Self::Fr3FromRelative { relative } => write!(f, "FR3 start from {relative}"),
            Self::Fr3FromModel { codons, labels, p } => write!(
                f,
                "FR3 length {codons} codons from the database model (n={labels}, p={})",
                probability(*p)
            ),
            Self::Fr3OutsideModel => {
                write!(f, "FR3 length outside the database model; chain default")
            }
            Self::ImgtGapsIgnored { reason } => write!(f, "IMGT gap numbering ignored: {reason}"),
            Self::ImgtFrameshiftMarker { codon } => {
                write!(f, "IMGT frameshift marker in codon {codon}")
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Probabilities round down, so only a certain choice reads `1.00`.
    #[test]
    fn notes_read_as_in_annotation_tsv() {
        let model = Note::Fr3FromModel {
            codons: 38,
            labels: 12,
            p: 0.999,
        };
        assert_eq!(
            model.to_string(),
            "FR3 length 38 codons from the database model (n=12, p=0.99)"
        );
        let dropped = Note::DroppedDerived {
            region: 4,
            reason: "fr3=9 follows anchor=6".into(),
        };
        assert_eq!(
            dropped.to_string(),
            "dropped derived fr3: fr3=9 follows anchor=6"
        );
        assert_eq!(
            Note::UnusualCdrLength { cdr: 2, codons: 13 }.to_string(),
            "unusual CDR2 length 13 codons"
        );
    }
}
