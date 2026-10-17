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
//! the reported distribution. At `p̂ = δ0`, `u_r = u0` and `v_r = 0` for
//! every kind whose match and mismatch scores have `λ0` as their own
//! Karlin–Altschul λ (all three by default), so the per-read model is then
//! the fixed model times `u0`, and the fixed model ([`ReadModel::fixed`], one
//! unit per score unit) is its exact special case.
//! `F` is a power of two, so `units/F` is exact in `f64`; a mismatch at
//! `p̂ = 10⁻⁴` is about `−6·10⁵` units, so `i64` holds any read.
//!
//! **The presence costs.** `P*_D` and `P*_J` ([`crate::config::DPrior`]) are
//! derived as prior odds, but the call uses them as significance thresholds
//! against chance hits in random sequence: the N/P bases of a junction, or
//! whatever a read holds where its J would be. Their null distribution does
//! not depend on `p̂` and is dominated by exact runs. A cost that is
//! constant in nats would make specificity depend on divergence: an exact
//! 5-mer is worth `5·a(p̂) = 6.93` nats at low `p̂`, above `λ0·P*_D = 6.62`
//! for `P*_D = 10.45`, while the fixed scale rejects it (10 < 10.45). Each
//! kind's cost is therefore read at the read's scale, as `P*_k/m` matched
//! columns of that kind (`m` its match score, 2):
//!
//! ```text
//! c_k(p̂) = (P*_k / m) · a(p̂)       units: ⌊P*_k · (m·u_k + v_k) / m⌋
//! ```
//!
//! An exact `n`-mer is worth `n·a(p̂)`, so exact runs of `⌊P*_k/m⌋` bases lose
//! and runs of `⌊P*_k/m⌋ + 1` win at every `p̂`, by `min(f, 1 − f)·a(p̂)` with
//! `f` the fraction of `P*_k/m` (for D, 0.17 and 0.22 for the synthetic and
//! the human IGH database): the fixed model's admitted exact-run classes,
//! without relying on ties. At `δ0` it is `λ0·P*_k`. The
//! floor is the faithful rounding: an integer value beats a real cost `c`
//! exactly when it beats `⌊c⌋`. Conditioning the cost on `p̂` does not use
//! the junction twice, because `p̂` comes from columns outside it.

use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::{AffineScoring, Score};
use crate::stats::KarlinAltschul;

/// An integer value in a read model's units (see the module documentation).
pub type Units = i64;

/// `F`: units per nat of the per-read model.
pub const UNITS_PER_NAT: f64 = 65536.0;

/// Which model the Viterbi call maximises.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CallModel {
    /// The read's own model (`ReadModel::at`): the same integers as the
    /// posterior, so the call is a mode of the reported distribution.
    #[default]
    PerRead,
    /// The fixed scores (`ReadModel::fixed`). The posterior still reads
    /// evidence at the read's divergence ([`EvidenceModel::nats`]), so the
    /// call need not be its mode.
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
    /// Match and mismatch score per rearranged kind ([`SegmentKind::REARRANGED`]), indexed by
    /// [`SegmentKind::index`].
    pub substitution: [(Score, Score); 3],
    /// Nats per score unit of each kind's presence cost `P*_k`: the
    /// Karlin–Altschul λ of the kind's scoring, the one
    /// [`crate::config::DPrior::cost`] derives `P*_k` with (1 without one).
    pub presence_lambda: [f64; 3],
}

impl EvidenceModel {
    /// The model of the V, D and J scorings, `None` when the V scoring has
    /// no Karlin–Altschul λ (its scores are then not log-odds).
    pub fn new(scorings: [&AffineScoring; 3]) -> Option<Self> {
        let lambda = |s: &AffineScoring| KarlinAltschul::ungapped(s).map(|ka| ka.lambda);
        Some(Self {
            lambda: lambda(scorings[0])?,
            substitution: scorings.map(|s| (s.matched, s.mismatch)),
            presence_lambda: scorings.map(|s| lambda(s).unwrap_or(1.0)),
        })
    }

    /// Evidence scales of the three rearranged kinds at divergence `p`.
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
    /// unquantised, the prior terms at `lambda`, and the presence costs at
    /// the λ they were derived with, `λ_k·P*_k = ln(N·L) − ln odds`.
    pub fn nats(&self, p: f64, presence: [f64; 3]) -> PathNats {
        self.nats_with(p, presence, DTerms::default())
    }

    /// [`EvidenceModel::nats`] with the D terms `d`; the D cost rises by
    /// what its `P*_D/m` columns gain from `d.base` and falls by `d.credit`
    /// ([`ReadModel::at_with`]).
    pub fn nats_with(&self, p: f64, presence: [f64; 3], d: DTerms) -> PathNats {
        let mut costs: [f64; 3] = std::array::from_fn(|k| self.presence_lambda[k] * presence[k]);
        let (matched, _) = self.substitution[SegmentKind::D.index()];
        costs[SegmentKind::D.index()] +=
            presence[SegmentKind::D.index()] / f64::from(matched) * d.base - d.credit;
        PathNats {
            lambda: self.lambda,
            evidence: self.evidence(p),
            presence: costs,
            d,
        }
    }
}

/// The D terms of a path model, in nats: geometric priors on how many
/// germline bases a D clip trims, and on the N/P lengths around it.
///
/// * `trim` — what a D clip pays per germline base it trims from its
///   allele's 5′ and 3′ ends ([`Placement::trims`]): `ln((μ + 1)/μ)` for a
///   geometric trim of mean `μ`, relative to an untrimmed D.
/// * `base` — what a D clip gains per read base it keeps: with geometric
///   N/P lengths (mean `μ`) the junction around a D pays `ln((μ + 1)/μ)`
///   per N/P base, and every base the D keeps is one fewer. Only D clips
///   gain it, so a V or J cut is valued as before.
/// * `credit` — what the D cost is lowered by: the trim cost of a D with the
///   prior's mean trims, `Σ μ·ln((μ + 1)/μ)` (1.85 nats for means 5 and 7).
///   The trims then re-rank Ds without making a D with typical trims harder
///   to call, while a chance hit deep inside a long allele pays more than
///   before.
///
/// `trim` and `base` are linear in a cut's germline or read position, so a
/// clip stays worth `exit(x) − entry(y)` and the DAG and the posterior keep
/// their recursions ([`super::PlacementDag`]).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DTerms {
    pub trim: [f64; 2],
    pub base: f64,
    pub credit: f64,
}

impl DTerms {
    /// The terms of geometric 5′ and 3′ trims of means `trim` and N/P
    /// lengths of mean `np` (`None`: that prior off).
    pub fn geometric(trim: Option<(f64, f64)>, np: Option<f64>) -> Self {
        let per_base = |mean: f64| ((mean + 1.0) / mean).ln();
        Self {
            trim: trim.map_or([0.0; 2], |(t5, t3)| [per_base(t5), per_base(t3)]),
            base: np.map_or(0.0, per_base),
            credit: trim.map_or(0.0, |(t5, t3)| t5 * per_base(t5) + t3 * per_base(t3)),
        }
    }

    pub fn is_off(&self) -> bool {
        *self == Self::default()
    }

    /// [`ReadModel::d_span`] in nats.
    #[cfg(test)]
    pub(crate) fn span(&self, placement: &Placement, span: &crate::clip::ClipSpan) -> f64 {
        if placement.kind != SegmentKind::D {
            return 0.0;
        }
        let (t5, t3) = placement.trims(span.r_lo, span.r_hi);
        f64::from(span.q_len()) * self.base
            - f64::from(t5) * self.trim[0]
            - f64::from(t3) * self.trim[1]
    }

    /// [`ReadModel::d_cut`] in nats, for D `placement` (the trims and the
    /// read offset `o` of the column, as integers, then the nats).
    #[inline]
    pub(crate) fn cut(&self, placement: &Placement, o: u32, r: u32) -> DCut {
        let (t5, t3) = placement.trims(r, r + 1);
        DCut {
            t5,
            t3,
            o,
            entry: f64::from(t5) * self.trim[0] + f64::from(o) * self.base,
            exit: f64::from(o + 1) * self.base - f64::from(t3) * self.trim[1],
        }
    }
}

/// The D terms of the cuts around one `Match` column ([`DTerms::cut`]).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct DCut {
    pub t5: u32,
    pub t3: u32,
    pub o: u32,
    /// Nats an entry cut before the column adds to its prefix.
    pub entry: f64,
    /// Nats an exit cut after the column adds to its prefix.
    pub exit: f64,
}

/// A read's path model in nats: the posterior is `P(path, cuts) ∝ e^{nats}`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathNats {
    /// Nats per score unit of the prior terms and of gaps.
    pub lambda: f64,
    /// Alignment evidence per rearranged kind ([`SegmentKind::REARRANGED`]), indexed by
    /// [`SegmentKind::index`].
    pub evidence: [EvidenceScale; 3],
    /// What including a segment costs, per rearranged kind.
    pub presence: [f64; 3],
    /// The D trim and N/P length terms.
    pub d: DTerms,
}

/// One read's model in integer units (see the module documentation).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ReadModel {
    /// Units per nat.
    per_nat: f64,
    /// Units per score unit of a prior term or gap.
    u0: Units,
    /// `(u_r, v_r)` per [`SegmentKind::index`].
    affine: [(Units, Units); 3],
    /// What including a segment costs, per [`SegmentKind::index`].
    presence: [Units; 3],
    /// [`DTerms`]: per trimmed 5′ and 3′ germline base, per kept read base
    /// (the credit is in the D cost).
    d_trim: [Units; 2],
    d_base: Units,
}

impl ReadModel {
    /// The fixed scale: one unit per score unit (`u0 = u_r = 1`, `v_r = 0`),
    /// `λ` nats each, and the presence costs `⌊P*_k⌋`.
    pub fn fixed(lambda: f64, presence: [f64; 3]) -> Self {
        Self {
            per_nat: 1.0 / lambda,
            u0: 1,
            affine: [(1, 0); 3],
            presence: presence.map(|cost| cost.floor() as Units),
            d_trim: [0; 2],
            d_base: 0,
        }
    }

    /// The read's own model at divergence `p`, in `1/F` nats, with the
    /// presence costs `c_k(p)` derived from `P*_k` (`presence`).
    #[cfg(test)]
    pub fn at(model: &EvidenceModel, p: f64, presence: [f64; 3]) -> Self {
        Self::at_with(model, p, presence, DTerms::default())
    }

    /// [`ReadModel::at`] with the D terms `d`. A matched D column is then
    /// worth `a(p̂) + d.base`, and the D cost is read as `P*_D/m` such
    /// columns, so untrimmed exact runs keep their classes
    /// ([`super::units`]), less the credit `d.credit`.
    pub fn at_with(model: &EvidenceModel, p: f64, presence: [f64; 3], d: DTerms) -> Self {
        let quantise = |nats: f64| (UNITS_PER_NAT * nats).round() as Units;
        let affine = model
            .evidence(p)
            .map(|e| (quantise(e.lambda), quantise(e.mu)));
        let d_base = quantise(d.base);
        let cost = |kind: SegmentKind| {
            let matched = model.substitution[kind.index()].0;
            let (u, v) = affine[kind.index()];
            let (base, credit) = match kind {
                SegmentKind::D => (d_base, quantise(d.credit)),
                _ => (0, 0),
            };
            let column = (Units::from(matched) * u + v + base) as f64;
            (presence[kind.index()] * column / f64::from(matched)).floor() as Units - credit
        };
        Self {
            per_nat: UNITS_PER_NAT,
            u0: quantise(model.lambda),
            affine,
            presence: SegmentKind::REARRANGED.map(cost),
            d_trim: d.trim.map(quantise),
            d_base,
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

    /// What including a segment of `kind` costs.
    #[inline]
    pub fn presence(&self, kind: SegmentKind) -> Units {
        self.presence[kind.index()]
    }

    /// Whether the D terms are on ([`ReadModel::d_cut`]).
    #[inline]
    pub fn has_d_terms(&self) -> bool {
        self.d_trim != [0; 2] || self.d_base != 0
    }

    /// The D terms of a cut of D `placement` on its `Match` column at read
    /// offset `o` (from the placement's start) and reference position `r`:
    /// what an entry cut before the column and an exit cut after it add to
    /// their prefix values ([`DTerms`]).
    #[inline]
    pub fn d_cut(&self, placement: &Placement, o: u32, r: u32) -> (Units, Units) {
        debug_assert_eq!(placement.kind, SegmentKind::D);
        let (t5, t3) = placement.trims(r, r + 1);
        (
            Units::from(t5) * self.d_trim[0] + Units::from(o) * self.d_base,
            Units::from(o + 1) * self.d_base - Units::from(t3) * self.d_trim[1],
        )
    }

    /// The D terms of a whole clip of `placement` retaining `span`: what
    /// [`ReadModel::d_cut`] adds at its two cuts, from the span's bounds.
    #[cfg(test)]
    pub fn d_span(&self, placement: &Placement, span: &crate::clip::ClipSpan) -> Units {
        if placement.kind != SegmentKind::D {
            return 0;
        }
        let (t5, t3) = placement.trims(span.r_lo, span.r_hi);
        Units::from(span.q_len()) * self.d_base
            - Units::from(t5) * self.d_trim[0]
            - Units::from(t3) * self.d_trim[1]
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
            presence: self.presence.map(nats),
            d: DTerms {
                trim: self.d_trim.map(nats),
                base: nats(self.d_base),
                // Already in the D cost.
                credit: 0.0,
            },
        }
    }
}

/// Informative columns and mismatches counted for the read's divergence
/// ([`divergence_counts`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Columns {
    pub columns: usize,
    pub mismatches: usize,
}

impl Columns {
    /// `p̂ = (k + ½)/(n + 1)` for `k` mismatches among `n` columns, clamped
    /// to `[10⁻⁴, 0.5]` ([`divergence`]).
    pub fn divergence(self) -> f64 {
        ((self.mismatches as f64 + 0.5) / (self.columns as f64 + 1.0)).clamp(1e-4, 0.5)
    }
}

/// The columns behind a read's divergence, counted over all its placements
/// and over its anchored ones ([`divergence_counts`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DivergenceCounts {
    /// What [`divergence`] counts.
    pub all: Columns,
    /// The same count with the V and J restricted to placements whose local
    /// alignment reached the anchor gate ([`Placement::anchor`]); the count
    /// of all placements when no V or J is anchored.
    pub anchored: Columns,
}

/// The read's divergence `p̂`, estimated before adjudication from the columns
/// no path contests: the best-scoring V placement's `Match` columns before
/// every D and J placement starts, and the best J's after every V and D
/// placement ends. Whichever path is called explains them by that V or J,
/// and they exclude the junction, so a D cost conditioned on `p̂` does not
/// use the junction's bases twice.
///
/// `p̂ = (k + ½)/(n + 1)` for `k` mismatches among `n` informative columns
/// (`N` columns, scoring 0, excluded): the posterior mean under a Jeffreys
/// `Beta(½, ½)` prior, i.e. the predictive mismatch probability of a new
/// column, such as a D's. With fewer than 50 such columns (a spurious D or J
/// placement far out), every `Match` column of the two placements counts.
/// Clamped to `[10⁻⁴, 0.5]`.
pub fn divergence(placements: &[Placement]) -> f64 {
    divergence_counts(placements).all.divergence()
}

/// The columns [`divergence`] counts, and the same count over the anchored
/// V and J placements only, in one pass over the placements. Junk V or J
/// hits in random sequence (which never anchor) would otherwise set the
/// mismatch rate that the open-flank cost reads ([`FlankTable`]).
pub fn divergence_counts(placements: &[Placement]) -> DivergenceCounts {
    const MIN_COLUMNS: usize = 50;
    type Best<'a> = Option<&'a Placement>;
    fn keep<'a>(best: &mut Best<'a>, p: &'a Placement) {
        if best.is_none_or(|b| p.score() > b.score()) {
            *best = Some(p);
        }
    }
    // Best V and J, over all placements and over the anchored ones.
    let (mut v, mut j, mut anchored_v, mut anchored_j): (Best, Best, Best, Best) =
        (None, None, None, None);
    let (mut lo, mut hi) = (u32::MAX, 0u32);
    for p in placements {
        match p.kind {
            SegmentKind::V => {
                keep(&mut v, p);
                if p.anchor.is_some() {
                    keep(&mut anchored_v, p);
                }
            }
            SegmentKind::J => {
                keep(&mut j, p);
                if p.anchor.is_some() {
                    keep(&mut anchored_j, p);
                }
            }
            _ => {}
        }
        if p.kind != SegmentKind::V {
            lo = lo.min(p.q_start());
        }
        if p.kind != SegmentKind::J {
            hi = hi.max(p.q_end());
        }
    }
    // (informative columns, mismatches) of V before `lo` and J from `hi` on,
    // or of the whole V and J when that leaves fewer than MIN_COLUMNS.
    let count = |v: Best, j: Best| {
        let within = |lo: u32, hi: u32| {
            let mut counted = Columns::default();
            for (placement, kept) in [(v, 0..lo), (j, hi..u32::MAX)] {
                let matches = placement.into_iter().flat_map(|p| p.clip.matches());
                for (_, entry, exit) in matches.filter(|(q, ..)| kept.contains(q)) {
                    counted.columns += usize::from(exit != entry);
                    counted.mismatches += usize::from(exit < entry);
                }
            }
            counted
        };
        let counted = within(lo, hi);
        if counted.columns < MIN_COLUMNS {
            within(u32::MAX, 0)
        } else {
            counted
        }
    };
    let all = count(v, j);
    let same = |a: Best, b: Best| match (a, b) {
        (Some(a), Some(b)) => std::ptr::eq(a, b),
        (None, None) => true,
        _ => false,
    };
    // Nothing anchored, or the anchored bests are the bests: the same count.
    let unchanged = same(anchored_v, v) && same(anchored_j, j);
    let anchored = if (anchored_v.is_none() && anchored_j.is_none()) || unchanged {
        all
    } else {
        count(anchored_v, anchored_j)
    };
    DivergenceCounts { all, anchored }
}

/// What a read's adjudication reads besides its placements: its divergence
/// `p̂` ([`divergence`]), its length and its anchored columns, which the
/// open-flank cost reads ([`FlankTable`]). One read's placement subsets
/// (the loci of a multi-locus database) share it, so their paths are
/// hypotheses of one model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReadContext {
    pub divergence: f64,
    /// Bases of the read (as oriented).
    pub query_len: u32,
    pub anchored: Columns,
}

impl ReadContext {
    /// The context of a read of `query_len` bases with `placements`.
    pub fn of(placements: &[Placement], query_len: u32) -> Self {
        let counts = divergence_counts(placements);
        Self {
            divergence: counts.all.divergence(),
            query_len,
            anchored: counts.anchored,
        }
    }

    /// The context of `placements` without their read: the read is taken to
    /// end where the last placement ends.
    pub fn of_placements(placements: &[Placement]) -> Self {
        let end = placements.iter().map(Placement::q_end).max().unwrap_or(0);
        Self::of(placements, end)
    }
}

/// The open end of a path whose first or last segment is a D: the missing
/// V before it, or the missing J after it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flank {
    /// A D that starts its path: no V before it.
    Start = 0,
    /// A D that ends its path: no J after it.
    End = 1,
}

/// The open-flank model with the V and J stages' gates
/// ([`crate::config::OpenFlank`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedOpenFlank {
    /// `N_open`: read bases beyond a D explained without a missing segment
    /// (N/P bases and the D's unaligned germline).
    pub n_open: u32,
    /// `π`: prior that the bases beyond an anchor are no junction.
    pub junk_prior: f64,
    /// Most visible bases of the missing segment the miss model reads.
    pub max_visible: u32,
    /// `ε_min`: the miss probability never drops below it (0 = no floor).
    pub miss_floor: f64,
    /// `S_acc` of the missing segment per [`Flank`] (V for the start, J for
    /// the end): the smallest score the stage keeps in a 64-base window.
    /// `None` (a stage without significance) turns that side off.
    pub accept: [Option<Score>; 2],
    /// The missing segment's match and mismatch scores per [`Flank`].
    pub substitution: [(Score, Score); 2],
    /// λ of the D scoring, the one `P*_D` was derived with: the fixed scale
    /// charges the cost as `c/λ_D` score units.
    pub d_lambda: f64,
}

/// `P(X ≥ k)` for `X ~ BetaBinomial(n, a, b)`, summed in linear space from
/// `P(X = 0) = Π_{i<n} (b + i)/(a + b + i)` by the pmf's ratio steps.
pub fn beta_binomial_tail(n: u32, a: f64, b: f64, k: u32) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n {
        return 0.0;
    }
    let pmf0 = (0..n).fold(1.0, |p, i| p * (b + f64::from(i)) / (a + b + f64::from(i)));
    tail_from(n, a, b, k, pmf0)
}

/// `P(X ≥ k)` given `P(X = 0)`, `1 ≤ k ≤ n`.
#[inline]
fn tail_from(n: u32, a: f64, b: f64, k: u32, pmf0: f64) -> f64 {
    let mut pmf = pmf0;
    let mut tail = 0.0;
    for x in 0..n {
        if x >= k {
            tail += pmf;
        }
        let (x, rest) = (f64::from(x), f64::from(n - x));
        pmf *= rest / (x + 1.0) * (x + a) / (rest - 1.0 + b);
    }
    (tail + pmf).min(1.0)
}

/// `c = ln(1 + π·(1 − ε)/ε)` nats at miss probability `eps`.
#[inline]
fn flank_cost(eps: f64, junk_prior: f64) -> f64 {
    if eps >= 1.0 {
        return 0.0;
    }
    let eps = eps.max(f64::MIN_POSITIVE);
    (junk_prior * (1.0 - eps) / eps).ln_1p()
}

/// `v ≥ 0` rounded to the nearest integer, halves up (`f64::round` for
/// non-negative values, without a library call).
#[inline]
fn round_non_negative(v: f64) -> i64 {
    let whole = v as i64;
    whole + i64::from(v - whole as f64 >= 0.5)
}

/// One cost of a [`FlankTable`] in each scale.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct FlankEntry {
    /// `c` in nats.
    cost: f64,
    /// In the call's units.
    units: Units,
    /// In the posterior's nats.
    nats: f64,
    /// On the fixed scale, `Δ`.
    score: Score,
}

/// One read's open-flank costs ([`ResolvedOpenFlank`]), per [`Flank`] and
/// per number `x` of read bases beyond the D. A D that ends its path with
/// `x` bases after it, or starts it with `x` before it, pays
///
/// ```text
/// c(x) = ln(1 + π·(1 − ε_eff(x))/ε_eff(x))
/// ```
///
/// the log-ratio of "no junction there, or one whose V or J was missed"
/// against "a junction whose V or J was missed". `ε(ℓ)` is the probability
/// that `ℓ` visible bases of the missing segment stay undetected
/// ([`FlankTable::miss_probabilities`]), at least `ε_min`, and the first
/// `N_open` bases are free: `ε_eff(x) = ε(min(x − N_open, ℓ_max))`
/// (`ε(ℓ ≤ 0) = 1`). Every cost is non-negative, non-decreasing in `x`, 0
/// within the allowance, constant from `x = N_open + ℓ_max` on and at most
/// `ln(1 + π/ε_min)`.
///
/// Each cost is kept in three scales: the call's units (per-read:
/// `round(c·F)`; fixed: `Δ`), the posterior's nats (per-read: those units
/// over `F`; fixed: `c`) and the fixed scores `Δ = ⌊P*_D + c/λ_D⌋ −
/// ⌊P*_D⌋`, so that an integer D score `S` beats `⌊P*_D⌋ + Δ` exactly when
/// it beats `P*_D + c/λ_D`.
#[derive(Clone, Debug, Default)]
pub(crate) struct FlankTable {
    query_len: u32,
    /// The read and model the costs were computed for; `None`: every cost 0.
    key: Option<(ReadContext, ResolvedOpenFlank)>,
    /// The D presence cost `P*_D` and call model the scales were derived for.
    scale: Option<(f64, CallModel)>,
    /// Per side: `ε` by visible length `ℓ`, the last one standing for every
    /// longer one.
    missed: [Vec<f64>; 2],
    /// Per side and index: the cost, computed for the indices a D uses
    /// (`cost` NaN elsewhere).
    entries: [Vec<FlankEntry>; 2],
    /// Per side: the visible lengths the entries cover.
    covers: [usize; 2],
    /// Per placement of the last [`FlankTable::fill`]: whether it can open
    /// (a V) or close (a J) a D's flank ([`FlankTable::closes`]).
    closes: Vec<bool>,
    /// `N_open` and `ℓ_max` of the model ([`index`]).
    n_open: u32,
    max_visible: u32,
    /// What the last [`FlankTable::fill`] read ([`FlankTable::refill`]).
    filled: Option<Filled>,
}

/// The arguments of a [`FlankTable::fill`], the placements by address and
/// length.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Filled {
    flank: Option<ResolvedOpenFlank>,
    presence: [u64; 3],
    read: ReadModel,
    calls: CallModel,
    context: ReadContext,
    placements: (usize, usize),
}

impl Filled {
    fn of(
        flank: Option<&ResolvedOpenFlank>,
        presence: [f64; 3],
        read: &ReadModel,
        calls: CallModel,
        context: &ReadContext,
        placements: &[Placement],
    ) -> Self {
        Self {
            flank: flank.copied(),
            presence: presence.map(f64::to_bits),
            read: *read,
            calls,
            context: *context,
            placements: (placements.as_ptr() as usize, placements.len()),
        }
    }
}

impl FlankTable {
    /// The table for a read with context `context` and D placements among
    /// `placements`, under `flank` (`None`: every cost 0), a D presence
    /// cost `P*_D` of `d_presence`, and a call maximising `calls`: the
    /// costs of every D in `placements` at either end of a path. They are
    /// kept while the read and model stay the same.
    pub fn fill(
        &mut self,
        flank: Option<&ResolvedOpenFlank>,
        presence: [f64; 3],
        read: &ReadModel,
        calls: CallModel,
        context: &ReadContext,
        placements: &[Placement],
    ) {
        self.fill_anew(flank, presence, read, calls, context, placements);
        self.filled = Some(Filled::of(
            flank, presence, read, calls, context, placements,
        ));
    }

    /// [`FlankTable::fill`], unless the last fill had the same arguments
    /// and the same placement slice: the posterior and the decoder of a
    /// call reuse the table its Viterbi pass filled, with the read's
    /// placements unchanged in between. A fill for another read or slice
    /// (by address and length, or context) fills anew; so does every
    /// [`FlankTable::fill`], which a call's adjudication starts with, so a
    /// buffer reused by the next read never finds the last read's table.
    pub fn refill(
        &mut self,
        flank: Option<&ResolvedOpenFlank>,
        presence: [f64; 3],
        read: &ReadModel,
        calls: CallModel,
        context: &ReadContext,
        placements: &[Placement],
    ) {
        if self.filled
            != Some(Filled::of(
                flank, presence, read, calls, context, placements,
            ))
        {
            self.fill(flank, presence, read, calls, context, placements);
        } else if cfg!(debug_assertions) {
            let mut fresh = FlankTable::default();
            fresh.fill(flank, presence, read, calls, context, placements);
            debug_assert!(
                fresh.closes == self.closes && fresh.same_costs(self, placements),
                "a flank table reused for other placements"
            );
        }
    }

    /// Whether `self` and `other` charge every D of `placements` alike.
    fn same_costs(&self, other: &Self, placements: &[Placement]) -> bool {
        placements.iter().all(|p| {
            [Flank::Start, Flank::End].into_iter().all(|side| {
                let (lo, hi) = (p.q_start(), p.q_end());
                self.units(side, p.kind, lo, hi) == other.units(side, p.kind, lo, hi)
                    && self.nats(side, p.kind, lo, hi).to_bits()
                        == other.nats(side, p.kind, lo, hi).to_bits()
                    && self.score(side, p.kind, lo, hi) == other.score(side, p.kind, lo, hi)
            })
        })
    }

    /// [`FlankTable::fill`] without recording its arguments.
    fn fill_anew(
        &mut self,
        flank: Option<&ResolvedOpenFlank>,
        presence: [f64; 3],
        read: &ReadModel,
        calls: CallModel,
        context: &ReadContext,
        placements: &[Placement],
    ) {
        let Some(flank) = flank else {
            self.query_len = context.query_len;
            self.clear();
            return;
        };
        let d_presence = presence[SegmentKind::D.index()];
        // A V or J that cannot pay its own presence cost under the call's
        // model is no evidence of the junction's V or J.
        self.closes.clear();
        self.closes.extend(placements.iter().map(|p| match p.kind {
            SegmentKind::V => pays(p, read, flank.substitution[Flank::Start as usize]),
            SegmentKind::J => pays(p, read, flank.substitution[Flank::End as usize]),
            _ => true,
        }));
        let ds = || placements.iter().filter(|p| p.kind == SegmentKind::D);
        // The farthest D from each read end.
        let beyond = |p: &Placement| [p.q_start(), context.query_len.saturating_sub(p.q_end())];
        let mut needed = [None::<u32>; 2];
        for x in ds().map(beyond) {
            for (needed, x) in needed.iter_mut().zip(x) {
                *needed = Some(needed.map_or(x, |n| n.max(x)));
            }
        }
        self.prepare(flank, d_presence, calls, context, needed);
        for x in ds().map(beyond) {
            for (side, x) in x.into_iter().enumerate() {
                self.cost(flank, side, index(flank, x));
            }
        }
    }

    /// [`FlankTable::fill`] with every cost up to `needed[side]` bases
    /// beyond a D on each side (`None`: none).
    #[cfg(test)]
    pub fn fill_to(
        &mut self,
        flank: &ResolvedOpenFlank,
        d_presence: f64,
        calls: CallModel,
        context: &ReadContext,
        needed: [Option<u32>; 2],
    ) {
        self.filled = None;
        self.prepare(flank, d_presence, calls, context, needed);
        for (side, needed) in needed.into_iter().enumerate() {
            for i in needed.map_or(0..0, |x| 0..index(flank, x) + 1) {
                self.cost(flank, side, i);
            }
        }
    }

    /// The miss probabilities up to `needed[side]` bases beyond a D, unless
    /// the table holds them for this read and model, and the scales.
    fn prepare(
        &mut self,
        flank: &ResolvedOpenFlank,
        d_presence: f64,
        calls: CallModel,
        context: &ReadContext,
        needed: [Option<u32>; 2],
    ) {
        self.query_len = context.query_len;
        (self.n_open, self.max_visible) = (flank.n_open, flank.max_visible);
        let key = (*context, *flank);
        let covered = (0..2).all(|side| {
            flank.accept[side].is_none()
                || needed[side].is_none_or(|x| index(flank, x) < self.covers[side])
        });
        if self.key != Some(key) || !covered {
            self.key = Some(key);
            for (side, needed) in needed.into_iter().enumerate() {
                self.misses(flank, context.anchored, side, needed);
            }
        }
        if self.scale != Some((d_presence, calls)) {
            self.scale = Some((d_presence, calls));
            for entry in self.entries.iter_mut().flatten() {
                if !entry.cost.is_nan() {
                    *entry = scaled(entry.cost, flank, d_presence, calls);
                }
            }
        }
    }

    /// The cost at index `i` of `side`, unless computed.
    #[inline]
    fn cost(&mut self, flank: &ResolvedOpenFlank, side: usize, i: usize) {
        let (missed, entries) = (&self.missed[side], &mut self.entries[side]);
        let (Some(&last), Some((d_presence, calls))) = (missed.last(), self.scale) else {
            return;
        };
        if let Some(entry) = entries.get_mut(i)
            && entry.cost.is_nan()
        {
            let missed = missed.get(i).copied().unwrap_or(last);
            let cost = flank_cost(missed, flank.junk_prior);
            *entry = scaled(cost, flank, d_presence, calls);
        }
    }

    /// `ε` of `side` for up to `needed` bases beyond a D, from the anchored
    /// columns `counts`, with its entries uncomputed; none when the side is
    /// off or no D needs it.
    fn misses(
        &mut self,
        flank: &ResolvedOpenFlank,
        counts: Columns,
        side: usize,
        needed: Option<u32>,
    ) {
        self.missed[side].clear();
        self.covers[side] = 0;
        self.entries[side].clear();
        let (m, mu) = flank.substitution[side];
        if let (Some(accept), Some(needed), true) = (flank.accept[side], needed, m > 0 && m > mu) {
            let covers = index(flank, needed) + 1;
            let missed = &mut self.missed[side];
            miss_probabilities(
                missed,
                accept,
                (m, mu),
                counts,
                covers as u32 - 1,
                flank.miss_floor,
            );
            self.covers[side] = covers;
            let uncomputed = FlankEntry {
                cost: f64::NAN,
                ..FlankEntry::default()
            };
            self.entries[side].resize(covers, uncomputed);
        }
    }

    /// Every cost 0.
    pub fn clear(&mut self) {
        self.filled = None;
        self.key = None;
        self.scale = None;
        self.entries.iter_mut().for_each(Vec::clear);
        self.missed.iter_mut().for_each(Vec::clear);
        self.covers = [0; 2];
        self.closes.clear();
    }

    /// Per placement of the last [`FlankTable::fill`]: whether a V can open
    /// a D's flank (precede a D) and a J close one (follow a D); empty when
    /// the cost is off. A V or J none of whose clips (keeping its outer end)
    /// is worth more than its presence cost under the call's model never
    /// pays for itself, so it is no evidence of the junction's V or J:
    /// without this, a junk J of a base or two after the D would close its
    /// flank for less than the cost.
    pub fn closes(&self) -> &[bool] {
        &self.closes
    }

    /// The entry of a segment of `kind` on `[q_start, q_end)` at the
    /// `flank` end of its path; `None` when it costs nothing.
    #[inline]
    fn entry(
        &self,
        flank: Flank,
        kind: SegmentKind,
        q_start: u32,
        q_end: u32,
    ) -> Option<&FlankEntry> {
        if kind != SegmentKind::D {
            return None;
        }
        let entries = &self.entries[flank as usize];
        let x = match flank {
            Flank::Start => q_start,
            Flank::End => self.query_len.saturating_sub(q_end),
        };
        let entry = entries.get(x.saturating_sub(self.n_open).min(self.max_visible) as usize);
        debug_assert!(
            entries.is_empty() || entry.is_some_and(|e| !e.cost.is_nan()),
            "a D whose open-flank cost was not computed"
        );
        entry
    }

    /// The cost, in the call's units, of a segment of `kind` on
    /// `[q_start, q_end)` at the `flank` end of its path.
    #[inline]
    pub fn units(&self, flank: Flank, kind: SegmentKind, q_start: u32, q_end: u32) -> Units {
        self.entry(flank, kind, q_start, q_end)
            .map_or(0, |e| e.units)
    }

    /// The same cost in the posterior's nats.
    #[inline]
    pub fn nats(&self, flank: Flank, kind: SegmentKind, q_start: u32, q_end: u32) -> f64 {
        self.entry(flank, kind, q_start, q_end)
            .map_or(0.0, |e| e.nats)
    }

    /// The same cost on the fixed scale.
    #[inline]
    pub fn score(&self, flank: Flank, kind: SegmentKind, q_start: u32, q_end: u32) -> Score {
        self.entry(flank, kind, q_start, q_end)
            .map_or(0, |e| e.score)
    }
}

/// `ε(ℓ)` for `ℓ = 0..=upto` into `eps`: the probability that `ℓ`
/// visible bases of a missing segment stay below its accept score
/// `S_acc` = `accept`, scored ungapped `(m, μ)`, when each base
/// mismatches with a probability drawn from `Beta(k + ½, n − k + ½)`
/// of the read's anchored columns (`k` of `n` mismatch). The span then
/// needs `k(ℓ) = ⌊(m·ℓ − S_acc)/(m − μ)⌋ + 1` mismatches (none can hide
/// it while `m·ℓ < S_acc`), which `X_ℓ ~ BetaBinomial(ℓ, a, b)` reaches
/// with probability `P(X_ℓ ≥ k(ℓ))`. A segment shows in `ℓ` bases when
/// any prefix of them does, so `ε(ℓ)` is the smallest of these over
/// `ℓ' ≤ ℓ` (non-increasing, still an upper bound of the true miss
/// probability), and at least `floor`.
///
/// The tail is carried along `ℓ` in O(1) per step, by the Pólya-urn
/// identities of the beta-binomial (`s = a + b + ℓ − 1`):
///
/// ```text
/// P(X_ℓ ≥ k)     = P(X_{ℓ−1} ≥ k) + P(X_{ℓ−1} = k−1)·(a + k − 1)/s
/// P(X_ℓ = k)     = (P(X_{ℓ−1} = k)·(b + ℓ − 1 − k) + P(X_{ℓ−1} = k−1)·(a + k − 1))/s
/// P(X_ℓ ≥ k + 1) = P(X_ℓ ≥ k) − P(X_ℓ = k)
/// ```
///
/// starting from `P(X_ℓ ≥ k) = 1 − Σ_{j<k} P(X_ℓ = j)` where `m·ℓ` first
/// reaches `S_acc`. The subtraction keeps an absolute error of a few ulps
/// of 1, so a tail below [`SUMMED_BELOW`] is summed instead, from
/// `P(X_ℓ = k)` up until its terms fall below `10⁻⁷` of the sum or the sum
/// reaches the running minimum (with `k(ℓ)` unchanged the tail only grows,
/// so `ε` then changes only where `k` steps). Once `ε` reaches `floor` it
/// stays there, and the list ends.
fn miss_probabilities(
    eps: &mut Vec<f64>,
    accept: Score,
    (m, mu): (Score, Score),
    counts: Columns,
    upto: u32,
    floor: f64,
) {
    let a = counts.mismatches as f64 + 0.5;
    let b = counts.columns.saturating_sub(counts.mismatches) as f64 + 0.5;
    // Both below 1 (no anchored columns): the pmf is U-shaped, so a
    // summed tail is summed whole.
    let unimodal = a >= 1.0 || b >= 1.0;
    eps.clear();
    // P(X_ℓ = 0).
    let mut zero = 1.0;
    // From the first ℓ with m·ℓ ≥ S_acc on: k(ℓ), P(X_ℓ = k − 1),
    // P(X_ℓ = k), and P(X_ℓ ≥ k) while it is carried (`None`: summed).
    let mut walk: Option<(u32, f64, f64, Option<f64>)> = None;
    let mut best = 1.0f64;
    for l in 0..=upto {
        let lf = f64::from(l);
        if l > 0 {
            let n = lf - 1.0;
            let step = 1.0 / (a + b + n);
            zero *= (b + n) * step;
            walk = walk.map(|(k, below, at, tail)| {
                let k_f = f64::from(k);
                let up = below * (k_f - 1.0 + a) * step;
                (
                    k,
                    // P(X_ℓ = k − 1) = P(X_{ℓ−1} = k − 1)·ℓ/(ℓ − k + 1)·(ℓ − k + b)/s.
                    below * lf / (lf - k_f + 1.0) * (lf - k_f + b) * step,
                    at * (n - k_f + b) * step + up,
                    tail.map(|t| t + up),
                )
            });
        }
        let score = m * l as Score;
        if score < accept {
            eps.push(best.max(floor));
            continue;
        }
        let k = ((score - accept) / (m - mu)) as u32 + 1;
        if k > l {
            // Not even every base mismatching hides the span.
            best = 0.0;
        } else {
            let (mut at_k, mut below, mut at, mut tail) = match walk {
                Some(w) => w,
                None => {
                    // The first span long enough to show: P(X_ℓ = j) for
                    // j < k by the pmf's ratio steps, and 1 minus their sum.
                    let (mut p, mut sum) = (zero, 0.0);
                    let mut below = 0.0;
                    for j in 0..k {
                        let x = f64::from(j);
                        sum += p;
                        below = p;
                        p *= (lf - x) / (x + 1.0) * (x + a) / (lf - x - 1.0 + b);
                    }
                    (k, below, p, Some(1.0 - sum))
                }
            };
            let stepped = at_k != k;
            while at_k < k {
                // P(X_ℓ ≥ k + 1) = P(X_ℓ ≥ k) − P(X_ℓ = k), and the pmf one up.
                tail = tail.map(|t| t - at);
                let x = f64::from(at_k);
                below = at;
                at *= (lf - x) / (x + 1.0) * (x + a) / (lf - x - 1.0 + b);
                at_k += 1;
            }
            tail = tail.filter(|&t| t >= SUMMED_BELOW);
            walk = Some((k, below, at, tail));
            let tail = match tail {
                Some(t) => t,
                // With k unchanged a summed tail only grows.
                None if !stepped && best < 1.0 => best,
                None => {
                    // Only a tail below the running minimum changes ε.
                    let (mut tail, mut term) = (0.0, at);
                    for x in k..=l {
                        tail += term;
                        if x == l || tail >= best {
                            break;
                        }
                        let xf = f64::from(x);
                        let next =
                            term * ((lf - xf) * (xf + a)) / ((xf + 1.0) * (lf - xf - 1.0 + b));
                        if unimodal && next < term && next < 1e-7 * tail {
                            break;
                        }
                        term = next;
                    }
                    tail
                }
            };
            best = best.min(tail);
        }
        let missed = best.max(floor);
        eps.push(missed);
        if floor > 0.0 && missed <= floor {
            break;
        }
    }
}

/// The smallest tail [`miss_probabilities`] carries by differences (whose
/// absolute error is a few ulps of 1); a smaller one is summed.
const SUMMED_BELOW: f64 = 1e-6;

/// The entry of cost `c` nats in each scale ([`FlankTable`]).
#[inline]
fn scaled(c: f64, flank: &ResolvedOpenFlank, d_presence: f64, calls: CallModel) -> FlankEntry {
    // P*_D ≥ 0 and c ≥ 0: truncation is the floor.
    let score = (d_presence + c / flank.d_lambda) as Score - d_presence as Score;
    let (units, nats) = match calls {
        CallModel::PerRead => {
            let units = round_non_negative(UNITS_PER_NAT * c);
            (units, units as f64 / UNITS_PER_NAT)
        }
        CallModel::FixedScale => (Units::from(score), c),
    };
    FlankEntry {
        cost: c,
        units,
        nats,
        score,
    }
}

/// Whether a clip of V or J placement `p` that keeps its outer end (a V
/// keeps its start, a J its end, as the first and last segment of a path)
/// is worth more than the kind's presence cost under `read`, its match and
/// mismatch scores being `(m, μ)`.
fn pays(p: &Placement, read: &ReadModel, (m, mu): (Score, Score)) -> bool {
    let cost = read.presence(p.kind);
    // The whole alignment, one such clip: its substitution columns and its
    // gaps (the rest of its score; `N` columns score 0 and are worth 0).
    let evidence = &p.evidence;
    let (matches, mismatches) = (evidence.matches as Score, evidence.mismatches as Score);
    let gaps = p.score() - m * matches - mu * mismatches;
    let whole = Units::from(matches) * read.column(p.kind, m)
        + Units::from(mismatches) * read.column(p.kind, mu)
        + read.prior(gaps);
    if whole > cost {
        return true;
    }
    // Otherwise every clip, along the columns.
    let mut steps = super::lattice::column_steps(p);
    match p.kind {
        // The best prefix: from the start (after any leading gap) to an exit.
        SegmentKind::V => {
            let Some((_, _, _, first)) = steps.next() else {
                return false;
            };
            let mut value = read.column(p.kind, first);
            if value > cost {
                return true;
            }
            for (_, _, gap, column) in steps {
                value += read.prior(gap) + read.column(p.kind, column);
                if value > cost {
                    return true;
                }
            }
            false
        }
        // The best suffix: from its best entry to the end.
        _ => {
            let (mut value, mut lowest) = (0, Units::MAX);
            for (_, _, gap, column) in steps {
                value += read.prior(gap);
                lowest = lowest.min(value);
                value += read.column(p.kind, column);
            }
            lowest != Units::MAX && value - lowest > cost
        }
    }
}

/// Where the cost of `x` bases beyond a D sits in a side's entries: its
/// visible length, `min(x − N_open, ℓ_max)`.
#[inline]
fn index(flank: &ResolvedOpenFlank, x: u32) -> usize {
    x.saturating_sub(flank.n_open).min(flank.max_visible) as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::N;
    use crate::scoring::AffineScoring;

    /// An ungapped placement of `kind` on `query[q_start..q_start + len]`
    /// with mismatches at the `mismatched` offsets.
    fn placement(
        kind: SegmentKind,
        query: &[u8],
        q_start: u32,
        len: u32,
        mismatched: &[u32],
    ) -> Placement {
        let mut reference = query[q_start as usize..(q_start + len) as usize].to_vec();
        for &k in mismatched {
            reference[k as usize] = (reference[k as usize] + 1) & 3;
        }
        crate::testing::placement(kind, query, &reference, q_start)
    }

    #[test]
    fn divergence_counts_uncontested_v_and_j_columns() {
        let query = crate::sim::Rng::new(3).codes(400);
        let v = |mismatched: &[u32]| placement(SegmentKind::V, &query, 0, 300, mismatched);
        let j = placement(SegmentKind::J, &query, 330, 60, &[59]);
        let d = placement(SegmentKind::D, &query, 290, 20, &[]);
        // An unmutated V alone: half a pseudo-mismatch in 301.
        assert_eq!(divergence(&[v(&[])]), 0.5 / 301.0);
        // V before the D (0..290) and J after V and D (330..390): the V
        // mismatch at 295 and the J one at 389 are contested or not.
        let both = [v(&[10, 295]), d, j.clone()];
        assert_eq!(divergence(&both), 2.5 / 351.0);
        // The best-scoring V counts, not the first.
        let worse = v(&[1, 2, 3, 4]);
        assert_eq!(divergence(&[worse, v(&[10]), j]), 2.5 / 361.0);
        // Nothing to count: the 0.5 clamp.
        assert_eq!(divergence(&[]), 0.5);
    }

    /// A spurious J far to the left leaves fewer than 50 uncontested V
    /// columns: every column of the best V and J counts instead.
    #[test]
    fn divergence_falls_back_to_whole_placements() {
        let query = crate::sim::Rng::new(4).codes(400);
        let v = placement(SegmentKind::V, &query, 0, 300, &[100, 200]);
        let j = placement(SegmentKind::J, &query, 40, 30, &[]);
        assert_eq!(divergence(&[v, j]), 2.5 / 331.0);
        // All-N query: no informative column.
        let blank = vec![N; 100];
        let n = placement(SegmentKind::V, &blank, 0, 100, &[]);
        assert_eq!(divergence(&[n]), 0.5);
    }

    fn model() -> EvidenceModel {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        EvidenceModel::new([&scoring; 3]).expect("λ")
    }

    /// Under `--calls fixed` the posterior charges a presence cost at the λ
    /// it was derived with: `P*_D = ln(N·L)/λ_D` costs `ln(N·L)` nats even
    /// when the D scoring's λ is not the V scoring's.
    #[test]
    fn fixed_scale_posterior_charges_the_derived_presence_prior_in_nats() {
        use crate::config::{DPrior, DSearchSpace};
        let (v, d) = (
            AffineScoring::new(2, -3, -13, -1),
            AffineScoring::new(2, -4, -13, -1),
        );
        let space = DSearchSpace {
            genes: 30,
            offsets: 1,
        };
        let presence = [
            0.0,
            DPrior::DERIVED.cost(space, &d),
            DPrior::DERIVED.cost(space, &v),
        ];
        let model = EvidenceModel::new([&v, &d, &v]).expect("λ");
        assert_ne!(model.presence_lambda[1], model.lambda);
        for p in [0.01, 0.112, 0.3] {
            let nats = model.nats(p, presence);
            assert_eq!(nats.presence[0], 0.0);
            for k in [1, 2] {
                assert!((nats.presence[k] - 750f64.ln()).abs() < 1e-9, "{nats:?}");
            }
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
            let presence = [0.0, d_presence, d_presence - 3.0];
            let read = ReadModel::at(&m, delta0, presence);
            let fixed = ReadModel::fixed(m.lambda, presence);
            let u0 = read.prior(1);
            assert_eq!(u0, (UNITS_PER_NAT * m.lambda).round() as Units);
            for kind in SegmentKind::REARRANGED {
                for s in [-3, 0, 2, 5] {
                    assert_eq!(read.column(kind, s), fixed.column(kind, s) * u0);
                }
            }
            for kind in SegmentKind::REARRANGED {
                assert_eq!(read.presence(kind), fixed.presence(kind) * u0);
            }
            let lambda = (read.nats().lambda, fixed.nats().lambda);
            assert!(
                (lambda.0 - lambda.1).abs() <= 0.5 / UNITS_PER_NAT,
                "{lambda:?}"
            );
        }
    }

    /// Exact runs of `⌊P*_k/m⌋` bases lose and one base more wins at every
    /// divergence, by `min(f, 1 − f)·a(p)`, `f` the fraction of `P*_k/m`:
    /// for a D and for a J.
    #[test]
    fn presence_costs_keep_the_fixed_scale_exact_run_classes() {
        let m = model();
        for cost in [7.27f64, 7.91, 10.34, 10.45, 10.55, 11.2] {
            let runs = (cost / 2.0).floor() as Units;
            let f = cost / 2.0 - runs as f64;
            let share = f.min(1.0 - f) - 1e-3;
            let presence = [0.0, cost, cost];
            let fixed = ReadModel::fixed(m.lambda, presence);
            for kind in [SegmentKind::D, SegmentKind::J] {
                // The fixed scale: 2k > ⌊P*⌋ exactly for k > runs.
                let charged = fixed.presence(kind);
                assert!(2 * runs <= charged && 2 * (runs + 1) > charged);
                for i in 0..64 {
                    let p = 1e-4 * (0.5e4f64).powf(f64::from(i) / 63.0);
                    let read = ReadModel::at(&m, p, presence);
                    let a = read.column(kind, 2);
                    let margin = (a as f64 * share) as Units;
                    let charged = read.presence(kind);
                    assert!(runs * a + margin < charged, "{cost} at {p}");
                    assert!((runs + 1) * a > charged + margin, "{cost} at {p}");
                }
            }
            assert_eq!(fixed.presence(SegmentKind::V), 0);
        }
    }

    #[test]
    fn nats_are_units_over_f() {
        let m = model();
        let read = ReadModel::at(&m, 0.02, [0.0, 10.45, 7.91]);
        let nats = read.nats();
        let exact = m.evidence(0.02)[0];
        assert!((nats.evidence[0].lambda - exact.lambda).abs() <= 0.5 / UNITS_PER_NAT);
        assert!((nats.evidence[0].mu - exact.mu).abs() <= 0.5 / UNITS_PER_NAT);
        for kind in SegmentKind::REARRANGED {
            assert_eq!(
                nats.presence[kind.index()] * UNITS_PER_NAT,
                read.presence(kind) as f64,
                "a power-of-two F divides exactly"
            );
        }
        // c_k(p) = (P*_k/m)·a(p) to within the quantisation.
        let a = (4.0f64 * 0.98).ln();
        let d = nats.presence[SegmentKind::D.index()];
        assert!((d - 10.45 / 2.0 * a).abs() < 1e-3, "{d}");
        let j = nats.presence[SegmentKind::J.index()];
        assert!((j - 7.91 / 2.0 * a).abs() < 1e-3, "{j}");
    }

    /// `ln Γ(x)` for `x > 0` (Lanczos, g = 7): the reference of the
    /// beta-binomial tests.
    fn ln_gamma(x: f64) -> f64 {
        const G: [f64; 9] = [
            0.999_999_999_999_809_9,
            676.520_368_121_885_1,
            -1_259.139_216_722_402_8,
            771.323_428_777_653_1,
            -176.615_029_162_140_6,
            12.507_343_278_686_905,
            -0.138_571_095_265_720_12,
            9.984_369_578_019_572e-6,
            1.505_632_735_149_311_6e-7,
        ];
        if x < 0.5 {
            let pi = std::f64::consts::PI;
            return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
        }
        let x = x - 1.0;
        let t = x + 7.5;
        let sum = G[1..]
            .iter()
            .enumerate()
            .fold(G[0], |s, (i, g)| s + g / (x + i as f64 + 1.0));
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + sum.ln()
    }

    /// `P(X ≥ k)` of `BetaBinomial(n, a, b)` summed from its pmf
    /// `C(n, x)·B(x + a, n − x + b)/B(a, b)` in log-gamma form.
    fn direct_tail(n: u32, a: f64, b: f64, k: u32) -> f64 {
        let ln_beta = |x: f64, y: f64| ln_gamma(x) + ln_gamma(y) - ln_gamma(x + y);
        let ln_choose = |n: u32, x: u32| {
            ln_gamma(f64::from(n) + 1.0)
                - ln_gamma(f64::from(x) + 1.0)
                - ln_gamma(f64::from(n - x) + 1.0)
        };
        (k..=n)
            .map(|x| {
                let (xf, rest) = (f64::from(x), f64::from(n - x));
                (ln_choose(n, x) + ln_beta(xf + a, rest + b) - ln_beta(a, b)).exp()
            })
            .sum()
    }

    fn close_to(got: f64, want: f64, relative: f64) -> bool {
        (got - want).abs() <= relative * want.abs() + 1e-300
    }

    const SHAPES: [(f64, f64); 6] = [
        (0.5, 0.5),
        (0.5, 300.5),
        (13.5, 32.5),
        (2.5, 47.5),
        (40.0, 3.0),
        (1.5, 120.5),
    ];

    #[test]
    fn beta_binomial_tail_matches_direct_sum() {
        for (a, b) in SHAPES {
            for n in 0..=20u32 {
                for k in 0..=n + 1 {
                    let (got, want) = (beta_binomial_tail(n, a, b, k), direct_tail(n, a, b, k));
                    assert!(
                        close_to(got, want, 1e-9),
                        "n {n} k {k} a {a} b {b}: {got} {want}"
                    );
                }
            }
        }
        // The table's walk along ℓ gives the same tails up to 64 bases,
        // as their running minimum over prefixes, and with a floor.
        let mut eps = Vec::new();
        for (mismatches, columns) in [(0, 0), (0, 45), (13, 45), (3, 300), (60, 120), (2, 4000)] {
            let counts = Columns {
                columns,
                mismatches,
            };
            let (a, b) = (mismatches as f64 + 0.5, (columns - mismatches) as f64 + 0.5);
            for floor in [0.0, 1e-3] {
                miss_probabilities(&mut eps, 16, (2, -3), counts, 64, floor);
                // The list ends where it reaches the floor.
                let len = eps.len();
                assert!(len == 65 || eps[len - 1] == floor);
                let mut best = 1.0f64;
                for l in 0..=64u32 {
                    if 2 * l >= 16 {
                        best = best.min(beta_binomial_tail(l, a, b, (2 * l - 16) / 5 + 1));
                    }
                    let want = best.max(floor);
                    let got = eps[(l as usize).min(len - 1)];
                    assert!(
                        close_to(got, want, 1e-6),
                        "{columns} {mismatches} {l}: {got} {want}"
                    );
                }
            }
        }
    }

    /// The open-flank model of the tests: IGH-like accept scores (V 24,
    /// J 16), `(2, −3)`, `N_open = 24`, `π = ½`, 64 visible bases,
    /// `ε_min = 10⁻³`.
    fn open_flank() -> ResolvedOpenFlank {
        ResolvedOpenFlank {
            n_open: 24,
            junk_prior: 0.5,
            max_visible: 64,
            miss_floor: 1e-3,
            accept: [Some(24), Some(16)],
            substitution: [(2, -3); 2],
            d_lambda: 0.6337,
        }
    }

    /// `X_max = N_open + ℓ_max`, from which the cost is constant.
    const X_MAX: u32 = 24 + 64;

    /// `c(x)` per side for a read of anchored `columns` with `mismatches`,
    /// up to 40 bases past `X_max`, read through the table's lookups.
    fn costs(flank: &ResolvedOpenFlank, columns: usize, mismatches: usize) -> [Vec<f64>; 2] {
        let context = ReadContext {
            divergence: 0.01,
            query_len: 1000,
            anchored: Columns {
                columns,
                mismatches,
            },
        };
        let mut table = FlankTable::default();
        table.fill_to(
            flank,
            10.77,
            CallModel::FixedScale,
            &context,
            [Some(X_MAX); 2],
        );
        let at = |side: Flank, x: u32| match side {
            Flank::Start => table.nats(side, SegmentKind::D, x, x + 10),
            Flank::End => table.nats(side, SegmentKind::D, 990 - x, 1000 - x),
        };
        [Flank::Start, Flank::End].map(|side| (0..=X_MAX + 40).map(|x| at(side, x)).collect())
    }

    #[test]
    fn open_flank_cost_shape() {
        let ceiling = (1.0f64 + 0.5 * (1.0 - 1e-3) / 1e-3).ln();
        let flank = open_flank();
        let x_max = X_MAX as usize;
        let clean = costs(&flank, 300, 1);
        for (side, accept) in [(0, 24), (1, 16)] {
            let c = &clean[side];
            // 0 while fewer than S_acc/m bases of the segment can show.
            let free = 24 + accept / 2;
            assert!(c[..free].iter().all(|&c| c == 0.0), "{c:?}");
            assert!(c[free] > 0.0);
            // Non-decreasing, constant past X_max, bounded by ln(1 + π/ε_min),
            // which a clean read far out reaches.
            assert!(c.windows(2).all(|w| w[0] <= w[1]), "{c:?}");
            assert!(c[x_max..].iter().all(|&v| v == c[x_max]));
            assert!(c.iter().all(|&v| v <= ceiling));
            assert!(ceiling - c[x_max] <= 1e-12, "{}", c[x_max]);
        }
        // Non-increasing in the mismatch rate.
        let by_rate: Vec<[Vec<f64>; 2]> = (0..=30)
            .map(|mismatches| costs(&flank, 60, mismatches))
            .collect();
        for pair in by_rate.windows(2) {
            for (lo, hi) in pair[0].iter().zip(&pair[1]) {
                assert!(lo.iter().zip(hi).all(|(a, b)| b <= a));
            }
        }
        // A J-only anchor at 30 % SHM (45 columns, 13 mismatches): the J
        // could have been missed, so 64 visible bases cost under a nat,
        // against the 6.2-nat ceiling of a clean read.
        let mutated = costs(&flank, 45, 13);
        let far = mutated[1][x_max];
        assert!(far < 1.0, "{far}");
        assert!(costs(&flank, 45, 0)[1][x_max] > 6.0);
        // A stage without significance: no cost on its side.
        let blind = ResolvedOpenFlank {
            accept: [None, Some(16)],
            ..flank
        };
        let c = costs(&blind, 300, 1);
        assert!(c[0].iter().all(|&v| v == 0.0) && c[1][x_max] > 0.0);
        // No model: an empty table, every cost 0.
        let mut table = FlankTable::default();
        let context = ReadContext::of(&[], 100);
        let read = ReadModel::fixed(1.0, [0.0, 10.0, 7.0]);
        table.fill(
            None,
            [0.0, 10.0, 7.0],
            &read,
            CallModel::PerRead,
            &context,
            &[],
        );
        assert_eq!(table.units(Flank::End, SegmentKind::D, 0, 10), 0);
    }

    /// The fixed scale charges `Δ = ⌊P*_D + c/λ_D⌋ − ⌊P*_D⌋`: an integer
    /// path term `S` beats `⌊P*_D⌋ + Δ` exactly when it beats the real
    /// cost `P*_D + c/λ_D`. The per-read model charges `round(c·F)` units
    /// and its posterior those units in nats; the fixed posterior `c`.
    #[test]
    fn fixed_scale_flank_is_exact() {
        let flank = open_flank();
        let x_max = X_MAX;
        let context = ReadContext {
            divergence: 0.01,
            query_len: 1000,
            anchored: Columns {
                columns: 200,
                mismatches: 4,
            },
        };
        let (mut fixed, mut per_read) = (FlankTable::default(), FlankTable::default());
        for d_presence in [0.0, 6.17, 10.55, 10.77, 11.0, 14.2] {
            fixed.fill_to(
                &flank,
                d_presence,
                CallModel::FixedScale,
                &context,
                [Some(x_max); 2],
            );
            per_read.fill_to(
                &flank,
                d_presence,
                CallModel::PerRead,
                &context,
                [Some(x_max); 2],
            );
            for side in [Flank::Start, Flank::End] {
                for x in 0..=x_max {
                    let (q_start, q_end) = match side {
                        Flank::Start => (x, x + 10),
                        Flank::End => (990 - x, 1000 - x),
                    };
                    let c = fixed.nats(side, SegmentKind::D, q_start, q_end);
                    let delta = fixed.score(side, SegmentKind::D, q_start, q_end);
                    let real = d_presence + c / flank.d_lambda;
                    let floor = d_presence.floor() as Score;
                    assert_eq!(floor + delta, real.floor() as Score);
                    for s in floor - 2..floor + delta + 3 {
                        assert_eq!(
                            s > floor + delta,
                            f64::from(s) > real,
                            "{d_presence} {c} {s}"
                        );
                    }
                    // The call's units are the scale's.
                    assert_eq!(
                        fixed.units(side, SegmentKind::D, q_start, q_end),
                        Units::from(delta)
                    );
                    let units = per_read.units(side, SegmentKind::D, q_start, q_end);
                    assert_eq!(units, (c * UNITS_PER_NAT).round() as Units);
                    let nats = per_read.nats(side, SegmentKind::D, q_start, q_end);
                    assert_eq!(nats * UNITS_PER_NAT, units as f64);
                    assert_eq!(per_read.score(side, SegmentKind::D, q_start, q_end), delta);
                    // Only a D pays.
                    for kind in [SegmentKind::V, SegmentKind::J] {
                        assert_eq!(per_read.units(side, kind, q_start, q_end), 0);
                    }
                }
            }
        }
    }

    /// The posterior and the decoder of a call reuse the table its Viterbi
    /// pass filled ([`FlankTable::refill`]), but only for the same
    /// arguments and placement slice: other placements (the next read in
    /// the same buffer is filled by its own adjudication first) or another
    /// context are filled anew, with their own costs.
    #[test]
    fn refill_reuses_only_the_same_reads_table() {
        let flank = open_flank();
        let query = crate::sim::Rng::new(5).codes(400);
        let mut v = placement(SegmentKind::V, &query, 0, 250, &[10]);
        v.anchor = Some(40);
        let near = placement(SegmentKind::D, &query, 260, 15, &[]);
        let far = placement(SegmentKind::D, &query, 300, 15, &[3]);
        let (presence, read) = ([0.0, 10.0, 7.0], ReadModel::fixed(1.0, [0.0, 10.0, 7.0]));
        let fill = |table: &mut FlankTable, placements: &[Placement], context: &ReadContext| {
            table.refill(
                Some(&flank),
                presence,
                &read,
                CallModel::PerRead,
                context,
                placements,
            )
        };
        let a = [v.clone(), near.clone()];
        let context = ReadContext::of(&a, 400);
        let mut table = FlankTable::default();
        fill(&mut table, &a, &context);
        let filled = table.filled;
        assert!(filled.is_some());
        fill(&mut table, &a, &context);
        assert_eq!(table.filled, filled);
        let end =
            |t: &FlankTable, p: &Placement| t.units(Flank::End, p.kind, p.q_start(), p.q_end());
        assert!(end(&table, &near) > 0);
        for (b, context) in [
            (vec![v.clone(), far.clone()], context),
            (a.to_vec(), ReadContext::of(&a, 380)),
        ] {
            fill(&mut table, &b, &context);
            assert_ne!(table.filled, filled);
            let mut fresh = FlankTable::default();
            fill(&mut fresh, &b, &context);
            for p in &b {
                assert_eq!(end(&table, p), end(&fresh, p));
            }
        }
        // A fill records its own arguments.
        table.fill(
            Some(&flank),
            presence,
            &read,
            CallModel::PerRead,
            &context,
            &a,
        );
        assert_eq!(table.filled, filled);
    }

    /// The anchored count reads only the V and J that reached the anchor
    /// gate: a junk J in random sequence does not lower it, and without an
    /// anchored V or J it is the count of all placements.
    #[test]
    fn anchored_columns_skip_unanchored_placements() {
        let query = crate::sim::Rng::new(3).codes(400);
        let mut v = placement(SegmentKind::V, &query, 0, 300, &[10, 20]);
        let mut junk = placement(
            SegmentKind::J,
            &query,
            340,
            60,
            &(0..30).map(|k| 2 * k).collect::<Vec<_>>(),
        );
        (v.anchor, junk.anchor) = (None, None);
        let counts = divergence_counts(&[v.clone(), junk.clone()]);
        assert_eq!(counts.anchored, counts.all);
        v.anchor = Some(40);
        let counts = divergence_counts(&[v.clone(), junk.clone()]);
        assert_eq!(counts.all.mismatches, 2 + 30);
        assert_eq!(
            counts.anchored,
            Columns {
                columns: 300,
                mismatches: 2
            }
        );
        junk.anchor = Some(40);
        let counts = divergence_counts(&[v, junk]);
        assert_eq!(counts.anchored, counts.all);
    }
}
