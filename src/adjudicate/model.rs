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

/// `P_D`: base penalty for introducing D, relaxed for a strong hypothesis.
///
/// Swig relaxes R-optimized's 12 down to 10 when the D scores at least 20, or
/// when its exact ungapped tract occurs in two different D sequences
/// ([`Placement::relaxes_d_cost`]). Sensitive-D uses a fixed 10. Legacy
/// profiles charge nothing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DPresencePenalty {
    pub base: Score,
    pub relaxed: Score,
    /// Unclipped D score that earns `relaxed`.
    pub strong_score: Score,
}

impl DPresencePenalty {
    /// R-optimised: 12, relaxed to 10 when the D scores at least 20 or is a repeated exact tract.
    pub const R_OPTIMIZED: Self = Self {
        base: 12,
        relaxed: 10,
        strong_score: 20,
    };
    /// Sensitive-D: fixed 10.
    pub const SENSITIVE_D: Self = Self {
        base: 10,
        relaxed: 10,
        strong_score: 0,
    };
    /// Truth and IgBLAST-agreement: D is free apart from having to beat the no-D path.
    pub const NONE: Self = Self {
        base: 0,
        relaxed: 0,
        strong_score: 0,
    };

    pub fn penalty(&self, placement: &Placement) -> Score {
        if placement.score() >= self.strong_score || placement.relaxes_d_cost {
            self.relaxed
        } else {
            self.base
        }
    }
}

/// The note's engineering objective, as a path model:
/// `Σ S_clip − P_D·[D] − tandem·[2nd D] + pair_bonus·[V…J] + per-base N/P cost`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EngineeringObjective {
    pub pair_bonus: Score,
    pub d_presence: DPresencePenalty,
    /// Extra prior cost for every D beyond the first.
    pub tandem_d_penalty: Score,
    /// Score per unassigned junction base (0 reproduces the note; negative
    /// values express an N/P length prior).
    pub insertion_per_base: Score,
    /// Extra presence cost charged when a D's longest exact run is below
    /// [`Self::d_short_run_floor`]. Suppresses distributed-D false positives
    /// that admit on identity but lack a solid ungapped core.
    pub d_short_run_penalty: Score,
    /// Exact-run length at or above which no short-run surcharge applies.
    pub d_short_run_floor: u32,
    /// Per-base penalty for aligned reference span above
    /// [`Self::d_long_span_soft_max`]. Prefers compact D alleles over ones that
    /// eat into N-addition; zeroed when soft-max is 0.
    pub d_long_span_soft_max: u32,
    pub d_long_span_penalty: Score,
}

impl Default for EngineeringObjective {
    fn default() -> Self {
        Self {
            pair_bonus: 24,
            d_presence: DPresencePenalty::R_OPTIMIZED,
            tandem_d_penalty: 12,
            insertion_per_base: 0,
            // Mild surcharge: enough to cut short-core FPs, light enough that
            // true short Ds on IgG still beat the no-D path.
            d_short_run_penalty: 1,
            d_short_run_floor: 6,
            // Off by default: even mild over-long penalties cost D recall.
            d_long_span_soft_max: 0,
            d_long_span_penalty: 0,
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
                let run = placement.evidence.longest_exact_run;
                let short = if run < self.d_short_run_floor {
                    self.d_short_run_penalty * (self.d_short_run_floor - run) as Score
                } else {
                    0
                };
                let span = placement
                    .alignment
                    .r_end
                    .saturating_sub(placement.alignment.r_start);
                let long = if self.d_long_span_soft_max > 0 && span > self.d_long_span_soft_max
                {
                    self.d_long_span_penalty * (span - self.d_long_span_soft_max) as Score
                } else {
                    0
                };
                -self.d_presence.penalty(placement) - tandem - short - long
            }
            SegmentKind::V => 0,
            SegmentKind::J => placement.path_adjust,
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
