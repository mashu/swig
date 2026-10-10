//! Receptor loci: which chain an allele belongs to, and which alleles may
//! form one rearrangement.
//!
//! One database may hold several loci, as Swig composes them: every heavy
//! and light chain of a BCR study (IGH, IGK, IGL), or every TR locus. A
//! rearrangement joins segments of one locus only, so each allele carries
//! the [`Loci`] its name places it in, and the pipeline keeps a read's J, D
//! and C candidates and every adjudicated path within one locus
//! ([`crate::pipeline`]).
//!
//! A locus is read from the allele label (`IGHV1-2*02`, `TRBJ2-7*01`), or,
//! when the label names none, from the whole FASTA header. Shared TRA/TRD
//! V genes (`TRAV14/DV4`) belong to both loci and pair with either. An
//! allele whose locus cannot be read belongs to every locus, so a database
//! of unrecognised names is assigned as a single locus.

use std::fmt;

/// An immunoglobulin or T-cell receptor locus.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Locus {
    Igh,
    Igk,
    Igl,
    Tra,
    Trb,
    Trd,
    Trg,
}

impl Locus {
    /// Every locus, in the order names are searched.
    pub const ALL: [Locus; 7] = [
        Locus::Igh,
        Locus::Igk,
        Locus::Igl,
        Locus::Tra,
        Locus::Trb,
        Locus::Trd,
        Locus::Trg,
    ];

    /// The AIRR `locus` value.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Igh => "IGH",
            Self::Igk => "IGK",
            Self::Igl => "IGL",
            Self::Tra => "TRA",
            Self::Trb => "TRB",
            Self::Trd => "TRD",
            Self::Trg => "TRG",
        }
    }

    /// Whether the locus has D genes (IGH, TRB, TRD).
    pub const fn has_d(self) -> bool {
        matches!(self, Self::Igh | Self::Trb | Self::Trd)
    }

    const fn bit(self) -> u8 {
        1 << self as u8
    }
}

impl fmt::Display for Locus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// A set of loci. [`Loci::ANY`] is the set of an allele whose locus is
/// unknown: it pairs with every allele.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Loci(u8);

impl Loci {
    pub const NONE: Self = Self(0);
    pub const ANY: Self = Self((1 << Locus::ALL.len()) - 1);

    pub const fn of(locus: Locus) -> Self {
        Self(locus.bit())
    }

    pub const fn contains(self, locus: Locus) -> bool {
        self.0 & locus.bit() != 0
    }

    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    pub const fn intersection(self, other: Self) -> Self {
        Self(self.0 & other.0)
    }

    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Whether an allele of these loci may join one of `other`'s.
    pub const fn pairs_with(self, other: Self) -> bool {
        !self.intersection(other).is_empty()
    }

    /// The loci in the set, in [`Locus::ALL`] order.
    pub fn iter(self) -> impl Iterator<Item = Locus> {
        Locus::ALL.into_iter().filter(move |&l| self.contains(l))
    }

    /// The one locus of the set, if it has exactly one.
    pub fn single(self) -> Option<Locus> {
        (self.0.count_ones() == 1).then(|| self.first().expect("one locus"))
    }

    /// The first locus of the set in [`Locus::ALL`] order.
    pub fn first(self) -> Option<Locus> {
        self.iter().next()
    }

    /// The loci named in `text`: the locus whose code (`IGH`, …, `TRG`)
    /// occurs first, case-insensitively, and TRD as well for a shared
    /// `TRAV…/DV…` gene. Empty when no code occurs.
    pub fn named_in(text: &str) -> Self {
        let upper = text.to_ascii_uppercase();
        let first = Locus::ALL
            .iter()
            .filter_map(|&locus| upper.find(locus.name()).map(|at| (at, locus)))
            .min();
        match first {
            None => Self::NONE,
            Some((_, Locus::Tra)) if upper.contains("/DV") => {
                Self::of(Locus::Tra).union(Self::of(Locus::Trd))
            }
            Some((_, locus)) => Self::of(locus),
        }
    }

    /// The loci of an allele labelled `label` in a record with FASTA header
    /// `header`: those its label names, else those the header names, else
    /// [`Loci::ANY`].
    pub fn infer(label: &str, header: &str) -> Self {
        [label, header]
            .into_iter()
            .map(Self::named_in)
            .find(|loci| !loci.is_empty())
            .unwrap_or(Self::ANY)
    }
}

impl Default for Loci {
    fn default() -> Self {
        Self::ANY
    }
}

impl fmt::Debug for Loci {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if *self == Self::ANY {
            return f.write_str("Loci(any)");
        }
        let names: Vec<&str> = self.iter().map(Locus::name).collect();
        write!(f, "Loci({})", names.join("|"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loci_are_read_from_labels_then_headers() {
        let of = Loci::of;
        assert_eq!(Loci::infer("IGHV1-2*02", ""), of(Locus::Igh));
        assert_eq!(Loci::infer("igkj1*01", ""), of(Locus::Igk));
        assert_eq!(Loci::infer("TRBD2*01", ""), of(Locus::Trb));
        assert_eq!(
            Loci::infer("TRAV14/DV4*01", ""),
            of(Locus::Tra).union(of(Locus::Trd))
        );
        assert_eq!(Loci::infer("TRDV2*01", ""), of(Locus::Trd));
        // The label wins over the header; the header is read when the
        // label names no locus.
        assert_eq!(
            Loci::infer("IGLV3-1*01", "X1 IGLV3-1*01 IGH-like"),
            of(Locus::Igl)
        );
        assert_eq!(Loci::infer("V1", "V1 Homo sapiens IGKV1-5"), of(Locus::Igk));
        assert_eq!(Loci::infer("V1", "V1 synthetic"), Loci::ANY);
    }

    #[test]
    fn sets_pair_through_a_shared_locus() {
        let tra_trd = Loci::of(Locus::Tra).union(Loci::of(Locus::Trd));
        assert!(tra_trd.pairs_with(Loci::of(Locus::Trd)));
        assert!(!Loci::of(Locus::Igh).pairs_with(Loci::of(Locus::Igk)));
        assert!(Loci::ANY.pairs_with(Loci::of(Locus::Trg)));
        assert_eq!(tra_trd.single(), None);
        assert_eq!(tra_trd.first(), Some(Locus::Tra));
        assert_eq!(Loci::of(Locus::Igl).single(), Some(Locus::Igl));
        assert_eq!(Loci::ANY.iter().count(), 7);
        assert!(Locus::Trd.has_d() && !Locus::Igk.has_d());
    }
}
