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
//! A locus is read from a `locus=` header key (as `swig prepare` writes it),
//! else from the gene name in the allele label (`IGHV1-2*02`, `TRBJ2-7*01`),
//! or, when the label names none, from a gene name in the FASTA header's
//! free text ([`Loci::infer`]). A locus code counts only as the start of a
//! gene name, so words such as `light` or `strain` name no locus. Shared
//! TRA/TRD V genes (`TRAV14/DV4`) belong to both loci and pair with either.
//! An allele whose locus cannot be read belongs to every locus, so a
//! database of unrecognised names is assigned as a single locus.

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

    /// The locus of an AIRR `locus` value, case-insensitively.
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|locus| locus.name().eq_ignore_ascii_case(name))
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

    /// The loci of a `,`-separated list of locus names (`IGH`, `TRA,TRD`),
    /// as [`Self::list`] writes it; `None` if a name is no locus or the
    /// list is empty.
    pub fn parse_list(text: &str) -> Option<Self> {
        let loci = text.split(',').try_fold(Self::NONE, |loci, name| {
            Locus::from_name(name).map(|locus| loci.union(Self::of(locus)))
        })?;
        (!loci.is_empty()).then_some(loci)
    }

    /// The names of the loci in the set, `,`-separated.
    pub fn list(self) -> String {
        self.iter().map(Locus::name).collect::<Vec<_>>().join(",")
    }

    /// The loci named in `text`: the locus of its first gene name, and TRD
    /// as well for a shared `TRAV…/DV…` gene. A gene name is a locus code
    /// (`IGH`, …, `TRG`, case-insensitively) that starts the text or follows
    /// a character other than a letter or digit, followed by a segment
    /// letter (`V`, `D`, `J`, `C`) or, for IGH, an isotype letter (`A`,
    /// `E`, `G`, `M`: `IGHG1`, `IGHEP1`). Empty when no gene name occurs: a
    /// code inside a word (`LIGHT`, `STRAIN`) or before another letter
    /// (`TRANSFER`) names none.
    pub fn named_in(text: &str) -> Self {
        let upper = text.to_ascii_uppercase();
        let bytes = upper.as_bytes();
        let gene_at = |at: usize| {
            if at > 0 && bytes[at - 1].is_ascii_alphanumeric() {
                return None;
            }
            let locus = Locus::ALL
                .into_iter()
                .find(|locus| bytes[at..].starts_with(locus.name().as_bytes()))?;
            let next = *bytes.get(at + locus.name().len())?;
            let isotype = locus == Locus::Igh && matches!(next, b'A' | b'E' | b'G' | b'M');
            (matches!(next, b'V' | b'D' | b'J' | b'C') || isotype).then_some(locus)
        };
        match (0..bytes.len()).find_map(gene_at) {
            None => Self::NONE,
            Some(Locus::Tra) if upper.contains("/DV") => {
                Self::of(Locus::Tra).union(Self::of(Locus::Trd))
            }
            Some(locus) => Self::of(locus),
        }
    }

    /// The loci of an allele labelled `label` in a record with FASTA header
    /// `header`: those its label names ([`Self::named_in`]), else those
    /// the header's free text names, else [`Loci::ANY`]. Free text is the
    /// whitespace-separated tokens that hold a digit, as gene names do
    /// (`IGKV1-5`, `M99641|IGKV1-5*01|Homo_sapiens`), and are no `key=value`
    /// attribute (`annotation=transfer:IGHV1-2*01` names a relative's gene).
    pub fn infer(label: &str, header: &str) -> Self {
        let named = Self::named_in(label);
        if !named.is_empty() {
            return named;
        }
        header
            .split_whitespace()
            .filter(|token| !token.contains('=') && token.bytes().any(|b| b.is_ascii_digit()))
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
        assert_eq!(
            Loci::infer("V1", "V1|IGLVI-20*01|Homo sapiens"),
            of(Locus::Igl)
        );
        assert_eq!(
            Loci::infer("V1", "V1 (TRAV14/DV4*01)"),
            of(Locus::Tra).union(of(Locus::Trd))
        );
        assert_eq!(
            Loci::infer("KV1-5*01", "KV1-5*01 Mus musculus strain BALB/c IGKV1-5"),
            of(Locus::Igk)
        );
        assert_eq!(Loci::infer("Homo_sapiens_IGHV1-2*02", ""), of(Locus::Igh));
    }

    /// Gene names of every segment and IGH isotype name their locus; a code
    /// inside a word, before another letter, in a `key=value` attribute or
    /// in a header word without a digit names none.
    #[test]
    fn only_gene_names_name_a_locus() {
        let of = Loci::of;
        for (label, locus) in [
            ("IGHG1*01", Locus::Igh),
            ("IGHM*01", Locus::Igh),
            ("IGHEP1*01", Locus::Igh),
            ("IGHA2*01", Locus::Igh),
            ("IGKC*01", Locus::Igk),
            ("IGLV(I)-20*01", Locus::Igl),
            ("TRBC1*01", Locus::Trb),
            ("TRDD1*01", Locus::Trd),
            ("TRGJP1*01", Locus::Trg),
        ] {
            assert_eq!(Loci::infer(label, ""), of(locus), "{label}");
        }
        for header in [
            "VH1-18*01 Homo sapiens germline, donor strain unknown",
            "VK1 human kappa light chain",
            "V1 full-length transcript, Australia",
            "ACC1.1 Homo sapiens kappa light chain, extra high weight",
            "V1*01 anchor=285 annotation=transfer:IGHV1-2*01",
            "V1*01 trace amounts; track 2",
            "V1 strain=BALB_c",
        ] {
            let label = header.split_whitespace().next().unwrap();
            assert_eq!(Loci::infer(label, header), Loci::ANY, "{header}");
        }
        assert_eq!(Loci::named_in("LIGHT STRAIN TRANSFER"), Loci::NONE);
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

    #[test]
    fn lists_of_loci_round_trip() {
        let tra_trd = Loci::of(Locus::Tra).union(Loci::of(Locus::Trd));
        assert_eq!(tra_trd.list(), "TRA,TRD");
        assert_eq!(Loci::parse_list("TRA,TRD"), Some(tra_trd));
        assert_eq!(Loci::parse_list("igk"), Some(Loci::of(Locus::Igk)));
        assert_eq!(Loci::parse_list(&Loci::ANY.list()), Some(Loci::ANY));
        for bad in ["", "IGH,", "XYZ", "IGH,TRX"] {
            assert_eq!(Loci::parse_list(bad), None, "{bad}");
        }
    }
}
