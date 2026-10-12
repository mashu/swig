//! Swig's reference metadata: the `SWIGMETA=` header field.
//!
//! Swig prepares every reference it ships, and every FASTA a user loads,
//! with region coordinates transferred from IMGT-annotated relatives
//! (`src/germline-preprocess.ts` of the web app). It writes them into the
//! header as thirteen comma-separated integers:
//!
//! | Index | Value |
//! |---|---|
//! | 0 | coding-frame start: 0-based offset (0–2) of the first complete codon, −1 if unknown |
//! | 1 | J only: 0-based last base of CDR3; the Trp/Phe118 codon starts one base later; −1 otherwise |
//! | 2–11 | V only: FWR1, CDR1, FWR2, CDR2, FWR3 as 0-based half-open `[start, end)` pairs, −1 if unknown |
//! | 12 | provenance: 1 IMGT-gapped, 2 AIRR-C, 3 IMGT-boundary-transfer, 4 validated-J-motif, 5 provided, 6 J-anchor-transfer |
//!
//! FWR3 ends after the Cys104 codon, so a V anchor is FWR3's end − 3. The
//! values are *given*, like header keys: they win over every value the
//! database annotation ([`crate::annotate`]) would derive, except an anchor
//! the allele's own sequence confidently contradicts, and fill only what an
//! explicit `anchor=`, `fr1=`, … key leaves open.

/// The parsed `SWIGMETA=` field of one record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SwigMeta {
    /// Offset of the first complete codon (0–2).
    pub frame: Option<u8>,
    /// J: 0-based last base of CDR3.
    pub cdr3_stop: Option<u32>,
    /// FWR1, CDR1, FWR2, CDR2, FWR3 as `[start, end)`.
    pub regions: [Option<(u32, u32)>; 5],
    /// Provenance of the coordinates (index 12).
    pub source: Option<&'static str>,
}

impl SwigMeta {
    /// The field after `SWIGMETA=`: exactly thirteen integers, the frame in
    /// −1…2 and the CDR3 stop at least −1. `None` when malformed.
    pub fn parse(text: &str) -> Option<Self> {
        let values: Vec<i64> = text
            .split(',')
            .map(|v| v.parse().ok())
            .collect::<Option<_>>()?;
        let values: [i64; 13] = values.try_into().ok()?;
        if !(-1..=2).contains(&values[0]) || values[1] < -1 {
            return None;
        }
        let position = |v: i64| u32::try_from(v).ok();
        let region = |k: usize| {
            let (start, end) = (values[2 + 2 * k], values[3 + 2 * k]);
            if start >= 0 && end > start {
                Some((position(start)?, position(end)?))
            } else {
                None
            }
        };
        Some(Self {
            frame: u8::try_from(values[0]).ok(),
            cdr3_stop: position(values[1]),
            regions: [region(0), region(1), region(2), region(3), region(4)],
            source: match values[12] {
                1 => Some("IMGT-gapped"),
                2 => Some("AIRR-C"),
                3 => Some("IMGT-boundary-transfer"),
                4 => Some("validated-J-motif"),
                5 => Some("provided"),
                6 => Some("J-anchor-transfer"),
                _ => Some("unclassified"),
            },
        })
    }

    /// The V region starts (FR1, CDR1, FR2, CDR2, FR3).
    pub fn region_starts(&self) -> [Option<u32>; 5] {
        self.regions.map(|r| r.map(|(start, _)| start))
    }

    /// The V anchor (first base of Cys104): FWR3's end − 3.
    pub fn v_anchor(&self) -> Option<u32> {
        self.regions[4].and_then(|(_, end)| end.checked_sub(3))
    }

    /// The J anchor (first base of Trp/Phe118): one past CDR3's last base.
    pub fn j_anchor(&self) -> Option<u32> {
        self.cdr3_stop.map(|stop| stop + 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_v_and_j_fields() {
        let v = SwigMeta::parse("0,-1,0,75,75,99,99,150,150,171,171,285,1").unwrap();
        assert_eq!(v.frame, Some(0));
        assert_eq!(
            v.region_starts(),
            [Some(0), Some(75), Some(99), Some(150), Some(171)]
        );
        assert_eq!(v.v_anchor(), Some(282));
        assert_eq!(v.j_anchor(), None);
        assert_eq!(v.source, Some("IMGT-gapped"));
        let j = SwigMeta::parse("1,18,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4").unwrap();
        assert_eq!(
            (j.frame, j.j_anchor(), j.v_anchor()),
            (Some(1), Some(19), None)
        );
        assert_eq!(j.source, Some("validated-J-motif"));
        let unknown = SwigMeta::parse("-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0").unwrap();
        assert_eq!(unknown.frame, None);
        assert_eq!(unknown.source, Some("unclassified"));
    }

    #[test]
    fn rejects_malformed_fields() {
        for text in [
            "",
            "0,-1",
            "3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1",
            "0,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1",
            "0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,7",
            "0,x,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1",
        ] {
            assert_eq!(SwigMeta::parse(text), None, "{text}");
        }
        // An empty or reversed interval is unknown, not an error.
        let m = SwigMeta::parse("0,-1,5,5,9,3,-1,-1,-1,-1,-1,-1,5").unwrap();
        assert_eq!(m.regions, [None; 5]);
    }
}
