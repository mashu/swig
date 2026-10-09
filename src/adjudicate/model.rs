//! Path grammar (which topologies are legal) and path model (how they score).
//!
//! Separating the two keeps "VDDJ allowed?" (a grammar question) apart from
//! "what does a second D cost?" (a scoring question). The DAG needs one
//! structural property from a [`PathModel`]: a segment contributes
//! *clip score + a bonus that does not depend on the cut positions*, and a
//! junction contributes a score linear in its inserted length. That
//! separability is what lets running maxima (or sums) evaluate every pair of
//! cuts in O(1) amortised per cut (see [`super::PlacementDag`]).

use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

/// Which segment sequences form a path. Transitions between kinds must go
/// forward in [`SegmentKind::REARRANGED`] order (V, D, J): the DAG visits kinds in
/// that order ([`super::PlacementDag`]).
pub trait Grammar: Send + Sync {
    fn can_start(&self, kind: SegmentKind) -> bool;
    /// `d_used` counts D segments up to and including `prev`.
    fn can_follow(&self, prev: SegmentKind, next: SegmentKind, d_used: u8) -> bool;
    /// Whether a path from `first` to `last` with `d_count` D segments is complete.
    fn can_end(&self, first: SegmentKind, last: SegmentKind, d_count: u8) -> bool;
}

/// V → D{0..max_d} → J with optional truncation at either end.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VdjGrammar {
    pub max_d: u8,
    /// Paths must start with V.
    pub require_v: bool,
    /// Paths must end with J.
    pub require_j: bool,
}

impl VdjGrammar {
    /// Heavy-chain default: VJ or VDJ, partial reads allowed.
    pub const VDJ: Self = Self {
        max_d: 1,
        require_v: false,
        require_j: false,
    };
    /// Tandem D as a first-class topology: VJ, VDJ, VDDJ.
    pub const TANDEM_D: Self = Self {
        max_d: 2,
        require_v: false,
        require_j: false,
    };
    /// Light chains: no D.
    pub const VJ: Self = Self {
        max_d: 0,
        require_v: false,
        require_j: false,
    };
}

impl Default for VdjGrammar {
    fn default() -> Self {
        Self::VDJ
    }
}

impl Grammar for VdjGrammar {
    fn can_start(&self, kind: SegmentKind) -> bool {
        match kind {
            SegmentKind::V => true,
            SegmentKind::D => !self.require_v && self.max_d > 0,
            SegmentKind::J => !self.require_v,
            SegmentKind::C => false,
        }
    }

    fn can_follow(&self, prev: SegmentKind, next: SegmentKind, d_used: u8) -> bool {
        use SegmentKind::{D, J, V};
        match (prev, next) {
            (V, D) | (D, D) => d_used < self.max_d,
            (V, J) | (D, J) => true,
            _ => false,
        }
    }

    fn can_end(&self, first: SegmentKind, last: SegmentKind, _d_count: u8) -> bool {
        let has_flank = first != SegmentKind::D || last != SegmentKind::D;
        has_flank && (!self.require_j || last == SegmentKind::J)
    }
}

pub trait PathModel: Send + Sync {
    /// Prior cost of including a D (`P*_D`), in score units and unrounded.
    /// Every D on a path pays it. On the fixed scale an integer path score
    /// beats it exactly when it beats `⌊P*_D⌋` ([`PathModel::d_presence_score`]).
    fn d_presence(&self) -> f64;

    /// `⌊P*_D⌋`: for integer scores `S`, `S > P*_D ⟺ S > ⌊P*_D⌋`, so this is
    /// the fixed-scale cost of a D, not an approximation of it.
    fn d_presence_score(&self) -> Score {
        self.d_presence().floor() as Score
    }

    /// Added to the clip score of `placement` besides the D presence cost;
    /// `d_ordinal` counts D segments up to and including this one (so the
    /// first D sees 1).
    fn segment_bonus(&self, placement: &Placement, d_ordinal: u8) -> Score;

    /// Score per unassigned query base (N/P region) between two segments.
    /// Junction scores are linear in the inserted length, which keeps both
    /// cuts of a transition separable (see [`super::PlacementDag`]).
    fn junction_per_base(&self, _prev: SegmentKind, _next: SegmentKind) -> Score {
        0
    }

    /// Edge score for `inserted` unassigned query bases between two segments.
    fn junction(&self, prev: SegmentKind, next: SegmentKind, inserted: u32) -> Score {
        self.junction_per_base(prev, next) * inserted as Score
    }

    /// Whole-path term (e.g. the V+J pair bonus).
    fn end(&self, _first: SegmentKind, _last: SegmentKind, _d_count: u8) -> Score {
        0
    }
}

impl<M: PathModel> PathModel for &M {
    fn d_presence(&self) -> f64 {
        (**self).d_presence()
    }

    fn segment_bonus(&self, placement: &Placement, d_ordinal: u8) -> Score {
        (**self).segment_bonus(placement, d_ordinal)
    }

    fn junction_per_base(&self, prev: SegmentKind, next: SegmentKind) -> Score {
        (**self).junction_per_base(prev, next)
    }

    fn end(&self, first: SegmentKind, last: SegmentKind, d_count: u8) -> Score {
        (**self).end(first, last, d_count)
    }
}

/// The path objective, `Σ S_clip − P_D·#D − tandem·[2nd D] + pair_bonus·[V…J]
/// + per-base N/P cost`, read as a log-likelihood ratio.
///
/// Substitution scores are log-odds (in units of `1/λ` nats) of "templated by
/// this germline" against "random sequence"; N/P bases score 0 under both.
/// The remaining terms are log prior odds of the path's hypotheses:
///
/// * `d_presence` (`P*_D`) — a D is one of `N` genes at one of `L` positions
///   of the junction, so its prior is spread over `N·L` hypotheses: including
///   the best one is worth it when `λ·S_D > ln(N·L) − ln(prior odds of D)`.
///   [`crate::config::DPrior`] derives it from the database at build time.
///   It is kept unrounded: the posterior uses it as it is, and the fixed-scale
///   call as `⌊P*_D⌋` ([`PathModel::d_presence_score`]).
/// * `tandem_d_penalty` — prior odds of a second D (DJ → DDJ), on top of its `P_D`.
/// * `pair_bonus` — prior odds of a complete V…J rearrangement over a partial
///   read. 0 by default: on full-length reads the J's own clip score already
///   completes the path, and 0–48 changed no call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EngineeringObjective {
    pub pair_bonus: Score,
    pub d_presence: f64,
    /// Extra prior cost for every D beyond the first.
    pub tandem_d_penalty: Score,
    /// Score per unassigned junction base (0 = no N/P length prior).
    pub insertion_per_base: Score,
}

impl Default for EngineeringObjective {
    fn default() -> Self {
        Self {
            pair_bonus: 0,
            d_presence: 11.0,
            tandem_d_penalty: 12,
            insertion_per_base: 0,
        }
    }
}

impl EngineeringObjective {
    /// Smallest D local score that can improve a fixed-scale path: `P_D + 1`,
    /// with `P_D = ⌊P*_D⌋`. It bounds nothing for a per-read call: there a D
    /// scoring 10 ≤ `P_D` (8 matches, 2 mismatches) is worth 6.40 nats at
    /// `p̂ = 0.3`, above its cost `c_D` of 5.38 for `P*_D = 10.45`
    /// ([`super::ReadModel`]).
    ///
    /// A clipped span of an optimal local alignment never outscores the
    /// whole alignment (its prefixes and suffixes score ≥ 0), so a D scoring
    /// at most `P_D` adds at most 0. Dropping it leaves its neighbours' cuts
    /// valid (the exit before it is ≤ the entry after it), so the no-D path
    /// scores at least as much and "with-D must strictly beat no-D" rejects
    /// it. With a negative per-base junction score a D could pay for itself
    /// by replacing N bases, so no floor is derived then.
    pub fn d_floor(&self) -> Score {
        if self.insertion_per_base < 0 {
            Score::MIN
        } else {
            self.d_presence_score() + 1
        }
    }
}

impl PathModel for EngineeringObjective {
    fn d_presence(&self) -> f64 {
        self.d_presence
    }

    fn segment_bonus(&self, placement: &Placement, d_ordinal: u8) -> Score {
        if placement.kind == SegmentKind::D && d_ordinal >= 2 {
            -self.tandem_d_penalty
        } else {
            0
        }
    }

    fn junction_per_base(&self, _prev: SegmentKind, _next: SegmentKind) -> Score {
        self.insertion_per_base
    }

    fn end(&self, first: SegmentKind, last: SegmentKind, _d_count: u8) -> Score {
        if first == SegmentKind::V && last == SegmentKind::J {
            self.pair_bonus
        } else {
            0
        }
    }
}
