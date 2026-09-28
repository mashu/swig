//! Query orientation from retrieval strength.
//!
//! The query is oriented before any DP: V seed votes are collected on the
//! forward strand and, unless the forward evidence is already overwhelming,
//! on the reverse complement. Retrieval on the reverse strand is cheap (one
//! sort of packed votes) compared with aligning both strands.

/// Strand of the query that the call refers to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Forward,
    Reverse,
}

impl Orientation {
    pub fn is_reverse(self) -> bool {
        self == Orientation::Reverse
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrientationPolicy {
    /// Consider the reverse complement at all.
    pub try_reverse: bool,
    /// Skip the reverse strand when the forward strand's best V candidate
    /// already has at least this many weighted votes.
    pub prune_reverse_at: Option<u32>,
}

impl Default for OrientationPolicy {
    fn default() -> Self {
        Self {
            try_reverse: true,
            prune_reverse_at: Some(160),
        }
    }
}

impl OrientationPolicy {
    /// Forward only (pre-oriented input).
    pub const FORWARD_ONLY: Self = Self {
        try_reverse: false,
        prune_reverse_at: None,
    };

    pub fn wants_reverse(&self, forward_strength: u32) -> bool {
        self.try_reverse
            && self
                .prune_reverse_at
                .map_or(true, |limit| forward_strength < limit)
    }

    /// Strictly stronger reverse evidence flips the strand; ties stay forward.
    pub fn choose(forward_strength: u32, reverse_strength: u32) -> Orientation {
        if reverse_strength > forward_strength {
            Orientation::Reverse
        } else {
            Orientation::Forward
        }
    }
}
