//! The per-read model: what every term of a path is worth for one read, in
//! integer units that the Viterbi call ([`super::PlacementDag`]) and the
//! forward–backward pass ([`super::posterior`]) both read.
//!
//! **Log-odds at the read's divergence.** `(2, −3)` scores are log-odds of
//! "templated by this germline" against "random" only at one divergence,
//! `δ0 = 1 − e^{2λ0}/4 = 11.2 %` ([`EvidenceModel::log_odds_divergence`]).
//! At the read's own divergence `p̂` ([`divergence`]) a match is worth
//! `a(p̂) = ln 4(1 − p̂)` nats and a mismatch `b(p̂) = ln(4p̂/3)`. For a
//! match/mismatch scoring that is an affine map of each substitution
//! column's score `s`, `λ_r·s + μ_r` ([`EvidenceScale::at_divergence`]), so
//! clip values stay separable in their cuts. An `N` column (`s = 0`) is worth
//! 0. Gaps and the prior terms (pair bonus, tandem, `ι`) stay at `λ0`: `p̂`
//! measures substitutions, and log prior odds are scale-free.
//!
//! **Integer units.** One unit is `1/F` nats, `F = 2¹⁶`. Each kind's affine
//! map is quantised as a whole, not column value by column value:
//!
//! ```text
//! column(kind, s) = s·u_r + v_r      u_r = round(F·λ_r),  v_r = round(F·μ_r)   (0 when s = 0)
//! prior(k)        = k·u0             u0  = round(F·λ0)
//! ```
//!
//! Every term of a path is then an integer, so the DAG's maximum and its ties
//! are exact, and the posterior reads the same integers as nats,
//! `P(path, cuts) ∝ e^{units/F}` ([`ReadModel::nats`]): the call is a mode of
//! the reported distribution. At `p̂ = δ0`, `u_r = u0` and `v_r = 0`, so the
//! per-read model is the fixed model times `u0`, and the fixed model
//! ([`ReadModel::fixed`], one unit per score unit) is its exact special case.
//! `F` is a power of two, so `units/F` is exact in `f64`; a mismatch at
//! `p̂ = 10⁻⁴` is about `−6·10⁵` units, so `i64` holds any read.
//!
//! **The D cost.** `P*_D` ([`crate::config::DPrior`]) is derived as prior odds,
//! but the call uses it as a significance threshold against chance D hits in
//! the random N/P sequence of a junction. Their null distribution does not
//! depend on `p̂` and is dominated by exact runs. A cost that is constant in
//! nats would make specificity depend on divergence: an exact 5-mer is worth
//! `5·a(p̂) = 6.93` nats at low `p̂`, above `λ0·P*_D = 6.62` for
//! `P*_D = 10.45`, while the fixed scale rejects it (10 < 10.45). The cost
//! is therefore read at the read's scale, as `P*_D/m` matched D columns
//! (`m` the D match score, 2):
//!
//! ```text
//! c_D(p̂) = (P*_D / m) · a(p̂)       units: ⌊P*_D · (m·u_D + v_D) / m⌋
//! ```
//!
//! An exact `k`-mer is worth `k·a(p̂)`, so exact runs of `⌊P*_D/m⌋` bases lose
//! and runs of `⌊P*_D/m⌋ + 1` win at every `p̂`, by `min(f, 1 − f)·a(p̂)` with
//! `f` the fraction of `P*_D/m` (0.17 and 0.22 for the synthetic and the
//! human IGH database): the fixed model's admitted exact-run classes,
//! without relying on ties. At `δ0` it is `λ0·P*_D`. The
//! floor is the faithful rounding: an integer value beats a real cost `c`
//! exactly when it beats `⌊c⌋`. Conditioning the cost on `p̂` does not use
//! the junction twice, because `p̂` comes from columns outside it.

use super::CalledSegment;
use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// An integer value in a [`ReadModel`]'s units.
pub type Units = i64;

/// `F`: units per nat of the per-read model.
pub const UNITS_PER_NAT: f64 = 65536.0;

/// Which model the Viterbi call maximises.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CallModel {
    /// The read's own model ([`ReadModel::at`]): the same integers as the
    /// posterior, so the call is a mode of the reported distribution.
    PerRead,
    /// The fixed scores ([`ReadModel::fixed`]). The posterior still reads
    /// evidence at the read's divergence ([`EvidenceModel::nats`]), so the
    /// call need not be its mode.
    #[default]
    FixedScale,
}

/// Nats per unit of alignment evidence: a substitution column scoring
/// `s ≠ 0` is worth `lambda · s + mu`, an `N` column (`s = 0`) nothing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvidenceScale {
    pub lambda: f64,
    pub mu: f64,
}

impl EvidenceScale {
    /// Scores read as `lambda` nats per unit, with no per-column term.
    pub fn fixed(lambda: f64) -> Self {
        Self { lambda, mu: 0.0 }
    }

    /// Log-odds of "templated at divergence `p`" against "random" (uniform
    /// background) for a match/mismatch scoring: a match is worth
    /// `ln 4(1 − p)` nats and a mismatch `ln(4p/3)`.
    pub fn at_divergence(p: f64, matched: Score, mismatch: Score) -> Self {
        let (a, b) = ((4.0 * (1.0 - p)).ln(), (4.0 * p / 3.0).ln());
        let lambda = (a - b) / f64::from(matched - mismatch);
        Self {
            lambda,
            mu: a - lambda * f64::from(matched),
        }
    }
}

/// The scorings a read's model derives from: `lambda`, the nats per unit of
/// the prior terms and gaps (the Karlin–Altschul λ0 of the V scoring), and
/// each kind's match and mismatch score, to read its alignments as log-odds
/// at the read's divergence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvidenceModel {
    pub lambda: f64,
    /// Match and mismatch score per [`SegmentKind::index`].
    pub substitution: [(Score, Score); 3],
}

impl EvidenceModel {
    /// Evidence scales of the three segment kinds at divergence `p`.
    pub fn evidence(&self, p: f64) -> [EvidenceScale; 3] {
        self.substitution
            .map(|(matched, mismatch)| EvidenceScale::at_divergence(p, matched, mismatch))
    }

    /// `δ0` of `kind`: the divergence at which its scores are log-odds at
    /// `lambda` (`λ_r = λ0`, `μ_r = 0`), where a match is worth `λ0·m`.
    pub fn log_odds_divergence(&self, kind: SegmentKind) -> f64 {
        let matched = f64::from(self.substitution[kind.index()].0);
        1.0 - (self.lambda * matched).exp() / 4.0
    }

    /// The posterior of a [`CallModel::FixedScale`] call: evidence at `p`,
    /// unquantised, and the prior terms, the D cost `P*_D` included, at
    /// `lambda`.
    pub fn nats(&self, p: f64, d_presence: f64) -> PathNats {
        PathNats {
            lambda: self.lambda,
            evidence: self.evidence(p),
            d_cost: self.lambda * d_presence,
        }
    }
}

/// A read's path model in nats: the posterior is `P(path, cuts) ∝ e^{nats}`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathNats {
    /// Nats per score unit of the prior terms and of gaps.
    pub lambda: f64,
    /// Alignment evidence per [`SegmentKind::index`].
    pub evidence: [EvidenceScale; 3],
    /// What including a D costs.
    pub d_cost: f64,
}

/// One read's model in integer units (see the module documentation).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReadModel {
    /// Units per nat.
    per_nat: f64,
    /// Units per score unit of a prior term or gap.
    u0: Units,
    /// `(u_r, v_r)` per [`SegmentKind::index`].
    affine: [(Units, Units); 3],
    /// What including a D costs.
    d_cost: Units,
}

impl ReadModel {
    /// The fixed scale: one unit per score unit (`u0 = u_r = 1`, `v_r = 0`),
    /// `λ` nats each, and the D cost `⌊P*_D⌋`.
    pub fn fixed(lambda: f64, d_presence: f64) -> Self {
        Self {
            per_nat: 1.0 / lambda,
            u0: 1,
            affine: [(1, 0); 3],
            d_cost: d_presence.floor() as Units,
        }
    }

    /// The read's own model at divergence `p`, in `1/F` nats, with the D cost
    /// `c_D(p)` derived from `P*_D` (`d_presence`).
    pub fn at(model: &EvidenceModel, p: f64, d_presence: f64) -> Self {
        let quantise = |nats: f64| (UNITS_PER_NAT * nats).round() as Units;
        let affine = model
            .evidence(p)
            .map(|e| (quantise(e.lambda), quantise(e.mu)));
        let matched = model.substitution[SegmentKind::D.index()].0;
        let (u, v) = affine[SegmentKind::D.index()];
        let column = (Units::from(matched) * u + v) as f64;
        Self {
            per_nat: UNITS_PER_NAT,
            u0: quantise(model.lambda),
            affine,
            d_cost: (d_presence * column / f64::from(matched)).floor() as Units,
        }
    }

    /// A substitution column of `kind` scoring `s`.
    #[inline]
    pub fn column(&self, kind: SegmentKind, s: Score) -> Units {
        if s == 0 {
            return 0;
        }
        let (u, v) = self.affine[kind.index()];
        Units::from(s) * u + v
    }

    /// A prior term or gap run scoring `k`.
    #[inline]
    pub fn prior(&self, k: Score) -> Units {
        Units::from(k) * self.u0
    }

    /// What including a D costs.
    pub fn d_cost(&self) -> Units {
        self.d_cost
    }

    /// The same model in nats, `units / F` exactly (up to `f64` rounding of
    /// sums), for the posterior.
    pub fn nats(&self) -> PathNats {
        let nats = |units: Units| units as f64 / self.per_nat;
        PathNats {
            lambda: nats(self.u0),
            evidence: self.affine.map(|(u, v)| EvidenceScale {
                lambda: nats(u),
                mu: nats(v),
            }),
            d_cost: nats(self.d_cost),
        }
    }
}

/// Divergence of a read from its called germlines: the fraction of
/// mismatched `Match` columns in the called V and J segments (a D is too
/// short to tell), with half a pseudo-mismatch so that an unmutated read gets
/// a small positive rate. Clamped to `[10⁻⁴, 0.5]`.
pub fn divergence(placements: &[Placement], segments: &[CalledSegment]) -> f64 {
    let (mut columns, mut mismatches) = (0usize, 0usize);
    for segment in segments {
        let p = &placements[segment.placement];
        if p.kind == SegmentKind::D {
            continue;
        }
        let span = segment.span.q_lo..segment.span.q_hi;
        for (_, entry, exit) in p.clip.matches().filter(|(q, ..)| span.contains(q)) {
            // N columns (score 0) carry no evidence either way.
            columns += usize::from(exit != entry);
            mismatches += usize::from(exit < entry);
        }
    }
    ((mismatches as f64 + 0.5) / (columns as f64 + 1.0)).clamp(1e-4, 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::AffineScoring;
    use crate::stats::KarlinAltschul;

    fn model() -> EvidenceModel {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        EvidenceModel {
            lambda: KarlinAltschul::ungapped(&scoring).expect("λ").lambda,
            substitution: [(2, -3); 3],
        }
    }

    fn assert_close(a: f64, b: f64, what: &str) {
        assert!(
            (a - b).abs() <= 1e-9 * (1.0 + a.abs()),
            "{what}: {a} vs {b}"
        );
    }

    #[test]
    fn evidence_scale_is_log_odds_at_the_divergence() {
        let at = |p: f64| EvidenceScale::at_divergence(p, 2, -3);
        let nats = |e: EvidenceScale, s: Score| e.lambda * f64::from(s) + e.mu;
        for p in [0.001, 0.05, 0.2] {
            assert_close(nats(at(p), 2), (4.0 * (1.0 - p)).ln(), "match");
            assert_close(nats(at(p), -3), (4.0 * p / 3.0).ln(), "mismatch");
        }
        // (2, −3) is log-odds at 11.2 % divergence: the fixed-λ model.
        let delta0 = model().log_odds_divergence(SegmentKind::V);
        assert!((delta0 - 0.112).abs() < 1e-3, "{delta0}");
        let base = at(delta0);
        assert!(
            (base.lambda - model().lambda).abs() < 1e-12 && base.mu.abs() < 1e-12,
            "{base:?}"
        );
    }

    /// At δ0 every term is the fixed model's times `u0`: the fixed scale is
    /// the per-read model's exact special case (for an integer `P*_D`).
    #[test]
    fn at_the_log_odds_divergence_units_are_fixed_scores_times_u0() {
        let m = model();
        let delta0 = m.log_odds_divergence(SegmentKind::D);
        for d_presence in [10.0, 11.0] {
            let read = ReadModel::at(&m, delta0, d_presence);
            let fixed = ReadModel::fixed(m.lambda, d_presence);
            let u0 = read.prior(1);
            assert_eq!(u0, (UNITS_PER_NAT * m.lambda).round() as Units);
            for kind in SegmentKind::ALL {
                for s in [-3, 0, 2, 5] {
                    assert_eq!(read.column(kind, s), fixed.column(kind, s) * u0);
                }
            }
            assert_eq!(read.d_cost(), fixed.d_cost() * u0);
            let lambda = (read.nats().lambda, fixed.nats().lambda);
            assert!(
                (lambda.0 - lambda.1).abs() <= 0.5 / UNITS_PER_NAT,
                "{lambda:?}"
            );
        }
    }

    /// Exact runs of `⌊P*_D/m⌋` bases lose and one base more wins at every
    /// divergence, by `min(f, 1 − f)·a(p)`, `f` the fraction of `P*_D/m`.
    #[test]
    fn d_cost_keeps_the_fixed_scale_exact_run_classes() {
        let m = model();
        for d_presence in [10.34f64, 10.45, 10.55, 11.2] {
            let runs = (d_presence / 2.0).floor() as Units;
            let f = d_presence / 2.0 - runs as f64;
            let share = f.min(1.0 - f) - 1e-3;
            let fixed = ReadModel::fixed(m.lambda, d_presence);
            // The fixed scale: 2k > ⌊P*⌋ exactly for k > runs.
            assert!(2 * runs <= fixed.d_cost() && 2 * (runs + 1) > fixed.d_cost());
            for i in 0..64 {
                let p = 1e-4 * (0.5e4f64).powf(f64::from(i) / 63.0);
                let read = ReadModel::at(&m, p, d_presence);
                let a = read.column(SegmentKind::D, 2);
                let margin = (a as f64 * share) as Units;
                assert!(runs * a + margin < read.d_cost(), "{d_presence} at {p}");
                assert!(
                    (runs + 1) * a > read.d_cost() + margin,
                    "{d_presence} at {p}"
                );
            }
        }
    }

    #[test]
    fn nats_are_units_over_f() {
        let m = model();
        let read = ReadModel::at(&m, 0.02, 10.45);
        let nats = read.nats();
        let exact = m.evidence(0.02)[0];
        assert!((nats.evidence[0].lambda - exact.lambda).abs() <= 0.5 / UNITS_PER_NAT);
        assert!((nats.evidence[0].mu - exact.mu).abs() <= 0.5 / UNITS_PER_NAT);
        assert_eq!(
            nats.d_cost * UNITS_PER_NAT,
            read.d_cost() as f64,
            "a power-of-two F divides exactly"
        );
        // c_D(p) = (P*_D/m)·a(p) to within the quantisation.
        let a = (4.0f64 * 0.98).ln();
        assert!(
            (nats.d_cost - 10.45 / 2.0 * a).abs() < 1e-3,
            "{}",
            nats.d_cost
        );
    }
}
