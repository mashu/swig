//! FR3 length learned from the database's annotated alleles.
//!
//! The CDR2/FR3 boundary (IMGT 66) has no conserved residue. It is fixed by
//! the IMGT gaps of FR3, which are constant within IG heavy (38 codons) and
//! kappa (36), but vary by gene family elsewhere: lambda 36/38, TRA 32–34,
//! TRB 37/38, TRG and TRD 37–39 (tcrdist3's IMGT tables). The motif search
//! measures the Trp41 → Cys104 spacing `d` in codons, and IMGT numbering
//! gives
//!
//! ```text
//! d = 14 + c + L
//! ```
//!
//! for a CDR2 of `c` codons and an FR3 of `L` codons (Cys104 included): CDR2
//! starts 15 codons after Trp41, and Cys104 is the last of FR3's `L`. With
//! the two lengths independent across the database's annotated alleles,
//!
//! ```text
//! P(L | d) ∝ P(L) · P_c(d − 14 − L)
//! P(L)     = (n_L + ½) / (N + ½·|S|)    for L ∈ S = [L_min − 1, min(L_max + 1, 39)]
//! P_c(c)   = (m_c + ½) / (M + ½·13)     for c ∈ 0..=12
//! ```
//!
//! where `n_L` of the `N` labels have FR3 length `L`, `m_c` of the `M`
//! labels with a CDR2 in range have CDR2 length `c`, and `L_min`, `L_max`
//! are the shortest and longest labels. Both are Jeffreys (½) estimates:
//! the support of `P(L)` reaches one codon beyond the observed lengths (but
//! not past 39, IMGT positions 66–104), so three labels of one length do
//! not make that length certain. The decision and its reported posterior
//! come from the same distribution.
//!
//! Labels are the alleles whose FR3 start is given (header or IMGT gaps)
//! and whose anchor is known. The model serves an allele with no annotated
//! relative (a relative's FR3 start is preferred: it is a homologous
//! position), and the Trp41 search of the second motif pass covers
//! `[14 + L_min, 26 + L_max]` codons together with the chain default.
//! Independence fails for a family with no label (human TRGV1–8 when only
//! TRGV9–11 are annotated, each family with its own CDR2 and FR3), whose
//! errors carry posteriors up to about 0.95; the annotation marks lengths
//! below 0.95 uncertain.

/// CDR2 lengths (codons) the model admits; the motif search uses the same range.
const CDR2: std::ops::RangeInclusive<u32> = 0..=12;
/// Longest FR3 (codons): IMGT positions 66–104.
const MAX_FR3: u32 = 39;
/// Codons from Trp41 to Cys104 not in CDR2 or FR3: `d = 14 + c + L`.
const SPACING_OFFSET: u32 = 14;

/// FR3 and CDR2 lengths of the annotated alleles (module docs).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fr3Model {
    /// Labels per FR3 length, indexed by codons (`1..=39`).
    lengths: [u32; MAX_FR3 as usize + 1],
    /// Labels per CDR2 length, indexed by codons.
    cdr2: [u32; 13],
    labels: u32,
}

impl Fr3Model {
    /// Learn from `(FR3 codons, CDR2 codons)` of annotated alleles. A length
    /// outside 1–39 is no IMGT FR3 and is skipped, as is a CDR2 length
    /// outside 0–12; `None` without a usable label.
    pub fn learn(samples: impl IntoIterator<Item = (u32, Option<u32>)>) -> Option<Self> {
        let mut model = Self {
            lengths: [0; MAX_FR3 as usize + 1],
            cdr2: [0; 13],
            labels: 0,
        };
        for (fr3, cdr2) in samples {
            if !Self::admits(fr3) {
                continue;
            }
            model.lengths[fr3 as usize] += 1;
            model.labels += 1;
            if let Some(c) = cdr2.filter(|c| CDR2.contains(c)) {
                model.cdr2[c as usize] += 1;
            }
        }
        (model.labels > 0).then_some(model)
    }

    /// Whether `codons` can be an IMGT FR3 length (1–39).
    pub fn admits(codons: u32) -> bool {
        (1..=MAX_FR3).contains(&codons)
    }

    /// Number of labelled alleles.
    pub fn labels(&self) -> u32 {
        self.labels
    }

    /// Shortest and longest labelled FR3 length.
    pub fn observed(&self) -> (u32, u32) {
        let seen = |l: &u32| self.lengths[*l as usize] > 0;
        let lo = (1..=MAX_FR3).find(seen).expect("at least one label");
        let hi = (1..=MAX_FR3).rev().find(seen).expect("at least one label");
        (lo, hi)
    }

    /// `S`: the observed lengths widened by one codon on each side, within 1–39.
    fn support(&self) -> std::ops::RangeInclusive<u32> {
        let (lo, hi) = self.observed();
        lo.saturating_sub(1).max(1)..=(hi + 1).min(MAX_FR3)
    }

    fn p_length(&self, l: u32) -> f64 {
        let bins = self.support().count() as f64;
        (f64::from(self.lengths[l as usize]) + 0.5) / (f64::from(self.labels) + 0.5 * bins)
    }

    fn p_cdr2(&self, c: i64) -> f64 {
        let Some(c) = u32::try_from(c).ok().filter(|c| CDR2.contains(c)) else {
            return 0.0;
        };
        let total: u32 = self.cdr2.iter().sum();
        let bins = CDR2.count() as f64;
        (f64::from(self.cdr2[c as usize]) + 0.5) / (f64::from(total) + 0.5 * bins)
    }

    /// The MAP FR3 length for a Trp41 → Cys104 spacing of `spacing` codons
    /// (`None`: no Trp41, so the prior `P(L)` alone) and its posterior
    /// probability. Ties go to the more often labelled length, then the
    /// shorter. `None` when no length in `S` leaves a CDR2 of 0–12 codons.
    pub fn length_for(&self, spacing: Option<u32>) -> Option<(u32, f64)> {
        let weights: Vec<(u32, f64)> = self
            .support()
            .map(|l| {
                let cdr2 = spacing.map_or(1.0, |d| {
                    self.p_cdr2(i64::from(d) - i64::from(SPACING_OFFSET) - i64::from(l))
                });
                (l, self.p_length(l) * cdr2)
            })
            .collect();
        let total: f64 = weights.iter().map(|w| w.1).sum();
        if total <= 0.0 {
            return None;
        }
        weights
            .iter()
            .max_by(|a, b| {
                let count = |l: u32| self.lengths[l as usize];
                a.1.total_cmp(&b.1)
                    .then(count(a.0).cmp(&count(b.0)))
                    .then(b.0.cmp(&a.0))
            })
            .map(|&(l, w)| (l, w / total))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute force from raw counts: every `(L, c)` with `14 + c + L = d`,
    /// `L` one codon around the labels (within 1–39) and `c` in 0–12,
    /// weighted by the Jeffreys estimates written out again.
    fn oracle(labels: &[(u32, Option<u32>)], spacing: u32) -> Option<(u32, f64)> {
        let lo = labels.iter().map(|l| l.0).min().unwrap();
        let hi = labels.iter().map(|l| l.0).max().unwrap();
        let support: Vec<u32> = (lo.max(2) - 1..=(hi + 1).min(39)).collect();
        let n = |l: u32| labels.iter().filter(|x| x.0 == l).count() as f64;
        let m = |c: u32| labels.iter().filter(|x| x.1 == Some(c)).count() as f64;
        let with_cdr2 = labels
            .iter()
            .filter(|x| x.1.is_some_and(|c| c <= 12))
            .count() as f64;
        let mut weights = Vec::new();
        for &l in &support {
            for c in 0..=12 {
                if 14 + c + l == spacing {
                    let p_l = (n(l) + 0.5) / (labels.len() as f64 + 0.5 * support.len() as f64);
                    let p_c = (m(c) + 0.5) / (with_cdr2 + 6.5);
                    weights.push((l, p_l * p_c));
                }
            }
        }
        let total: f64 = weights.iter().map(|w| w.1).sum();
        let best = weights.iter().copied().max_by(|a, b| {
            a.1.total_cmp(&b.1)
                .then(n(a.0).total_cmp(&n(b.0)))
                .then(b.0.cmp(&a.0))
        })?;
        Some((best.0, best.1 / total))
    }

    #[test]
    fn map_length_matches_brute_force() {
        // TRB-like: two FR3 lengths with family-specific CDR2 lengths, and
        // one rare length.
        let labels: Vec<(u32, Option<u32>)> = (0..40)
            .map(|i| {
                if i % 2 == 0 {
                    (37, Some(6))
                } else {
                    (38, Some(6))
                }
            })
            .chain([(37, Some(5)), (38, Some(7)), (36, Some(6)), (36, None)])
            .collect();
        let model = Fr3Model::learn(labels.iter().copied()).unwrap();
        for spacing in 40..80 {
            let got = model.length_for(Some(spacing));
            let want = oracle(&labels, spacing);
            assert_eq!(got.map(|g| g.0), want.map(|w| w.0), "{spacing}");
            if let (Some(g), Some(w)) = (got, want) {
                assert!((g.1 - w.1).abs() < 1e-12, "{spacing}: {g:?} {w:?}");
            }
        }
        assert_eq!(model.length_for(Some(14 + 6 + 37)).unwrap().0, 37);
        assert_eq!(model.length_for(Some(14 + 6 + 38)).unwrap().0, 38);
        assert_eq!(model.observed(), (36, 38));
    }

    /// Three labels of one length leave room for its neighbours: the
    /// posterior never reaches certainty where another length is feasible.
    #[test]
    fn sparse_labels_are_not_certain() {
        let model = Fr3Model::learn([(38, Some(8)), (38, Some(7)), (38, Some(8))]).unwrap();
        // S = 37..=39, so d = 14 + c + L with c in 0..=12 spans 51..=65, and
        // two or more lengths are feasible for d in 52..=64.
        for spacing in 52..=64 {
            let (_, p) = model.length_for(Some(spacing)).unwrap();
            assert!(p < 0.999, "{spacing}: {p}");
        }
        let (l, p) = model.length_for(Some(14 + 8 + 38)).unwrap();
        assert_eq!(l, 38);
        assert!(p > 0.5, "{p}");
        assert_eq!(model.length_for(Some(50)), None);
        assert_eq!(model.length_for(Some(66)), None);
        // Without a spacing the prior decides.
        assert_eq!(model.length_for(None).unwrap().0, 38);
        // IMGT FR3 holds at most 39 codons; longer labels are not used.
        let capped = Fr3Model::learn([(39, None), (45, None)]).unwrap();
        assert_eq!((capped.labels(), capped.observed()), (1, (39, 39)));
        assert!(Fr3Model::learn([(0, None), (40, Some(3))]).is_none());
    }
}
