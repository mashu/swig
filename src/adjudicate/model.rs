//! Path grammar (which topologies are legal) and path model (how they score).
//!
//! Separating the two keeps "VDDJ allowed?" (a grammar question) apart from
//! "what does a second D cost?" (a scoring question). The DAG needs one
//! structural property from a [`PathModel`]: a segment contributes
//! *clip score + a bonus that does not depend on the cut positions*, and a
//! junction contributes a function of its inserted length. That separability
//! is what makes every cut evaluable in O(log n) (see [`super::PlacementDag`]).

use crate::germline::SegmentKind;
use crate::placement::Placement;
use crate::scoring::Score;

pub trait Grammar: Send + Sync {
    fn can_start(&self, kind: SegmentKind) -> bool;
    /// `d_used` counts D segments up to and including `prev`.
    fn can_follow(&self, prev: SegmentKind, next: SegmentKind, d_used: u8) -> bool;
    fn can_end(&self, first: SegmentKind, last: SegmentKind) -> bool;
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

    fn can_end(&self, first: SegmentKind, last: SegmentKind) -> bool {
        let has_flank = first != SegmentKind::D || last != SegmentKind::D;
        has_flank && (!self.require_j || last == SegmentKind::J)
    }
}

pub trait PathModel: Send + Sync {
    /// Added to the clip score of `placement`; `d_ordinal` counts D segments up
    /// to and including this one (so the first D sees 1).
    fn segment_bonus(&self, placement: &Placement, d_ordinal: u8) -> Score;

    /// Edge score for `inserted` unassigned query bases between two segments (N/P region).
    fn junction(&self, _prev: SegmentKind, _next: SegmentKind, _inserted: u32) -> Score {
        0
    }

    /// Whole-path term (e.g. the V+J pair bonus).
    fn end(&self, _first: SegmentKind, _last: SegmentKind, _d_count: u8) -> Score {
        0
    }
}

/// The path objective, `Σ S_clip − P_D·#D − tandem·[2nd D] + pair_bonus·[V…J]
/// + per-base N/P cost`, read as a log-likelihood ratio.
///
/// Substitution scores are log-odds (in units of `1/λ` nats) of "templated by
/// this germline" against "random sequence"; N/P bases score 0 under both.
/// The remaining terms are log prior odds of the path's hypotheses:
///
/// * `d_presence` (`P_D`) — a D is one of `N` genes at one of `L` positions of
///   the junction, so its prior is spread over `N·L` hypotheses: including
///   the best one is worth it when `λ·S_D > ln(N·L) − ln(prior odds of D)`.
///   [`crate::config::DPrior`] derives it from the database at build time.
/// * `tandem_d_penalty` — prior odds of a second D (DJ → DDJ), on top of its `P_D`.
/// * `pair_bonus` — prior odds of a complete V…J rearrangement over a partial read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EngineeringObjective {
    pub pair_bonus: Score,
    pub d_presence: Score,
    /// Extra prior cost for every D beyond the first.
    pub tandem_d_penalty: Score,
    /// Score per unassigned junction base (0 = no N/P length prior).
    pub insertion_per_base: Score,
}

impl Default for EngineeringObjective {
    fn default() -> Self {
        Self {
            pair_bonus: 24,
            d_presence: 11,
            tandem_d_penalty: 12,
            insertion_per_base: 0,
        }
    }
}

impl PathModel for EngineeringObjective {
    fn segment_bonus(&self, placement: &Placement, d_ordinal: u8) -> Score {
        match placement.kind {
            SegmentKind::D => {
                let tandem = if d_ordinal >= 2 {
                    self.tandem_d_penalty
                } else {
                    0
                };
                -self.d_presence - tandem
            }
            SegmentKind::V | SegmentKind::J => 0,
        }
    }

    fn junction(&self, _prev: SegmentKind, _next: SegmentKind, inserted: u32) -> Score {
        self.insertion_per_base * inserted as Score
    }

    fn end(&self, first: SegmentKind, last: SegmentKind, _d_count: u8) -> Score {
        if first == SegmentKind::V && last == SegmentKind::J {
            self.pair_bonus
        } else {
            0
        }
    }
}
