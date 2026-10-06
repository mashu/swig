//! Declarative engine configuration and construction of the default assigner.
//!
//! Configuration is plain data; [`crate::profiles::CallingProfile`] produces
//! the note's calling profiles as values of [`EngineConfig`]. Settings that
//! can be *derived* losslessly from the admission gates (the D L-mer length
//! and the D traceback score floor) are derived here, so a profile change can
//! never make the prefilters silently drop admissible D hits.

use std::fmt;
use std::sync::Arc;

use crate::adjudicate::{EngineeringObjective, PlacementDag, VdjGrammar};
use crate::admission::{AnyOf, DistributedEvidence, ExactRunFloor};
use crate::align::InterSequenceGotoh;
use crate::germline::{GermlineSet, SegmentKind};
use crate::index::IndexError;
use crate::orient::StrandIndex;
use crate::pipeline::{Assigner, AssignerOptions};
use crate::refine::{NoRefine, OuterExtend};
use crate::retrieve::{
    AmbiguityBudget, ExhaustiveRetriever, LmerMaskRetriever, SeedVoteParams, SeedVoteRetriever,
    UngappedParams, UngappedPrefilter, MAX_LMER,
};
use crate::scoring::{AffineScoring, Score};
use crate::seed::SeedShape;
use crate::segment::{HitAcceptance, KeepPolicy, PreparedSegment, SegmentStage};
use crate::window::{DownstreamWindow, ExpandedResidual, WithAnchorWindow};

/// V: seeded, banded, vote-ranked, optionally ungapped-pre-ranked.
#[derive(Clone, Debug, PartialEq)]
pub struct VSegmentConfig {
    pub scoring: AffineScoring,
    pub primary: SeedShape,
    pub fallback: Option<SeedShape>,
    pub seeds: SeedVoteParams,
    pub budget: AmbiguityBudget,
    /// Ungapped diagonal pre-ranking before banded DP (`None` = off).
    /// Off by default: it is lossy.
    pub ungapped: Option<UngappedParams>,
    pub keep: KeepPolicy,
    /// Minimum substitution columns and identity of an accepted local hit.
    pub accept: HitAcceptance,
    /// Skip score-only DP for same-gene siblings that cannot enter the top-N,
    /// using a sound upper bound from the gene's representative.
    pub family_delta: bool,
}

/// J: every allele, unbanded, in the read right of V.
#[derive(Clone, Debug, PartialEq)]
pub struct JSegmentConfig {
    pub scoring: AffineScoring,
    pub keep: KeepPolicy,
    pub accept: HitAcceptance,
    pub window: DownstreamWindow,
}

/// D: L-mer prefilter, unbanded DP inside junction windows, admission gates.
#[derive(Clone, Debug, PartialEq)]
pub struct DSegmentConfig {
    pub scoring: AffineScoring,
    /// L-mer length of the prefilter; `None` derives the largest L that is
    /// lossless for the configured gates (see [`DSegmentConfig::lossless_lmer`]).
    pub lmer: Option<usize>,
    pub keep: KeepPolicy,
    /// `min_d_match`.
    pub min_exact_run: u32,
    /// Distributed-D gate (`None` = disabled; the note makes it AER-R only).
    pub distributed: Option<DistributedEvidence>,
}

impl DSegmentConfig {
    /// Largest L for which "shares ≥ 1 L-mer" is implied by every enabled gate:
    /// an exact run of `min_exact_run` contains an L-mer for `L ≤ min_exact_run`;
    /// ≥ 80 % identity forces an exact run ≥ 4 (pigeonhole), so distributed-D caps L at 4.
    pub fn lossless_lmer(&self) -> usize {
        let run = self.min_exact_run.max(1) as usize;
        let distributed_cap = if self.distributed.is_some() {
            4
        } else {
            usize::MAX
        };
        run.min(distributed_cap).min(MAX_LMER)
    }

    /// Traceback score floor implied by the gates: an alignment containing an
    /// exact run of `k` has local score `≥ k·m`; a distributed-D hit rescored
    /// under the D tuple itself scores `≥ min_rescore`. The configured floor is
    /// lowered to that bound, never raised above it.
    pub fn lossless_min_score(&self) -> Score {
        let run_floor = self.scoring.matched * self.min_exact_run as Score;
        let floor = match &self.distributed {
            None => run_floor,
            Some(gate) if gate.scoring == self.scoring => run_floor.min(gate.min_rescore),
            Some(_) => self.keep.min_score.min(run_floor),
        };
        floor.max(self.keep.min_score.min(floor))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EngineConfig {
    pub v: VSegmentConfig,
    pub d: DSegmentConfig,
    pub j: JSegmentConfig,
    pub d_window: ExpandedResidual,
    /// Add the conserved-anchor D window when alleles carry `anchor=` positions.
    pub anchor_window: bool,
    pub grammar: VdjGrammar,
    pub objective: EngineeringObjective,
    /// Minimum retained query bases after clipping, per [`SegmentKind::index`].
    pub min_retained: [u32; 3],
    pub options: AssignerOptions,
}

impl Default for VSegmentConfig {
    /// Gap-open −9, at least 24 substitution columns, identity ≥ 0.60. The
    /// vote budget keeps 3 genes and walks up to 16 while votes stay within
    /// the ambiguity margin; the 3 best accepted alignments (and ties) are kept.
    fn default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -9, -1),
            primary: SeedShape::contiguous(9),
            fallback: Some(SeedShape::contiguous(7)),
            seeds: SeedVoteParams::default(),
            budget: AmbiguityBudget::default(),
            ungapped: None,
            keep: KeepPolicy {
                top: 3,
                min_score: 1,
                cap: 64,
            },
            accept: HitAcceptance {
                min_subst: 24,
                min_identity: 0.60,
            },
            family_delta: true,
        }
    }
}

impl Default for JSegmentConfig {
    /// `(2, −3, −17, −2)`, at least 10 substitution columns, identity ≥ 0.60.
    fn default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -17, -2),
            keep: KeepPolicy {
                top: 3,
                min_score: 1,
                cap: 64,
            },
            accept: HitAcceptance {
                min_subst: 10,
                min_identity: 0.60,
            },
            window: DownstreamWindow::default(),
        }
    }
}

impl Default for DSegmentConfig {
    fn default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -13, -1),
            lmer: None,
            // Every admissible D hit (up to 32): a distributed-D hit below the
            // best few must still reach adjudication.
            keep: KeepPolicy {
                top: 32,
                min_score: 0,
                cap: 32,
            },
            min_exact_run: 4,
            distributed: Some(crate::admission::DistributedEvidence::default()),
        }
    }
}

impl Default for EngineConfig {
    fn default() -> Self {
        Self {
            v: VSegmentConfig::default(),
            d: DSegmentConfig::default(),
            j: JSegmentConfig::default(),
            d_window: ExpandedResidual::default(),
            anchor_window: true,
            grammar: VdjGrammar::VDJ,
            objective: EngineeringObjective::default(),
            // V ≥ 24 aligned bases, J ≥ 10, D ≥ min_d_match. Query bases
            // retained after a clip use the same floors.
            min_retained: [24, 4, 10],
            options: AssignerOptions::default(),
        }
    }
}

#[derive(Debug)]
pub enum BuildError {
    Index {
        kind: SegmentKind,
        source: IndexError,
    },
    WrongKind {
        expected: SegmentKind,
        found: SegmentKind,
    },
    Empty(SegmentKind),
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Index { kind, source } => write!(f, "{} index: {source}", kind.label()),
            Self::WrongKind { expected, found } => {
                write!(
                    f,
                    "expected a {} germline set, got {}",
                    expected.label(),
                    found.label()
                )
            }
            Self::Empty(kind) => write!(f, "the {} germline set is empty", kind.label()),
        }
    }
}

impl std::error::Error for BuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Index { source, .. } => Some(source),
            _ => None,
        }
    }
}

fn check(set: &GermlineSet, expected: SegmentKind) -> Result<(), BuildError> {
    if set.kind() != expected {
        return Err(BuildError::WrongKind {
            expected,
            found: set.kind(),
        });
    }
    if set.is_empty() {
        return Err(BuildError::Empty(expected));
    }
    Ok(())
}

/// Retriever of the default V stage.
pub type DefaultVRetriever = UngappedPrefilter<SeedVoteRetriever>;

fn v_stage(
    set: GermlineSet,
    config: &VSegmentConfig,
) -> Result<SegmentStage<DefaultVRetriever>, BuildError> {
    let kind = set.kind();
    let seeds = SeedVoteRetriever::build(
        &set,
        config.primary.clone(),
        config.fallback.clone(),
        config.seeds,
        config.budget,
    )
    .map_err(|source| BuildError::Index { kind, source })?;
    let segment = Arc::new(PreparedSegment::new(
        set,
        &config.scoring,
        config.scoring.gaps,
    ));
    let retriever = UngappedPrefilter::new(seeds, Arc::clone(&segment), config.ungapped);
    let family_delta = config
        .family_delta
        .then(|| crate::family_bound::FamilyBoundTable::build(segment.set()));
    let mut stage = SegmentStage::new(segment, retriever, config.keep);
    stage.accept = config.accept;
    stage.outer = OuterExtend::Prefix;
    stage.family_delta = family_delta;
    Ok(stage)
}

fn j_stage(set: GermlineSet, config: &JSegmentConfig) -> SegmentStage<ExhaustiveRetriever> {
    let retriever = ExhaustiveRetriever { alleles: set.len() };
    let segment = PreparedSegment::new(set, &config.scoring, config.scoring.gaps);
    let mut stage = SegmentStage::new(Arc::new(segment), retriever, config.keep);
    stage.accept = config.accept;
    stage.outer = OuterExtend::Suffix;
    stage
}

fn d_stage(set: GermlineSet, config: &DSegmentConfig) -> SegmentStage<LmerMaskRetriever> {
    let lmer = config.lmer.unwrap_or_else(|| config.lossless_lmer());
    let retriever = LmerMaskRetriever::build(&set, lmer, 1);
    let segment = PreparedSegment::new(set, &config.scoring, config.scoring.gaps);
    let keep = KeepPolicy {
        min_score: config.lossless_min_score(),
        ..config.keep
    };
    let mut stage = SegmentStage::new(Arc::new(segment), retriever, keep);
    stage.accept = HitAcceptance {
        min_subst: config.min_exact_run,
        min_identity: 0.60,
    };
    stage
}

/// The default engine: batched banded Gotoh (AVX2 inter-sequence kernel when
/// available), k-mer orientation, seed-vote V retrieval, exhaustive J right
/// of V, L-mer D prefilter, expanded residual + anchor D windows,
/// exact-run-or-distributed admission, placement DAG.
pub type DefaultAssigner = Assigner;

impl Assigner {
    pub fn from_config(
        config: &EngineConfig,
        v: GermlineSet,
        d: Option<GermlineSet>,
        j: GermlineSet,
    ) -> Result<Self, BuildError> {
        check(&v, SegmentKind::V)?;
        check(&j, SegmentKind::J)?;
        if let Some(d) = &d {
            check(d, SegmentKind::D)?;
        }
        let d_window = if config.anchor_window {
            WithAnchorWindow::new(config.d_window, &v, &j)
        } else {
            WithAnchorWindow::disabled(config.d_window)
        };
        let grammar = if d.is_some() {
            config.grammar
        } else {
            VdjGrammar {
                max_d: 0,
                ..config.grammar
            }
        };
        let strand = StrandIndex::build([&v, &j], StrandIndex::K);
        Ok(Assigner {
            v: v_stage(v, &config.v)?,
            j: j_stage(j, &config.j),
            d: d.map(|set| d_stage(set, &config.d)),
            strand,
            j_window: config.j.window,
            aligner: InterSequenceGotoh::new(),
            d_window,
            d_admission: AnyOf(
                ExactRunFloor {
                    min_run: config.d.min_exact_run,
                },
                config.d.distributed,
            ),
            adjudicator: PlacementDag::new(grammar, config.objective, config.min_retained),
            refiner: NoRefine,
            options: config.options,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derived_d_settings_are_lossless_bounds() {
        let mut d = DSegmentConfig::default();
        assert_eq!((d.lossless_lmer(), d.lossless_min_score()), (4, 8));
        d.distributed = None;
        d.min_exact_run = 6;
        assert_eq!((d.lossless_lmer(), d.lossless_min_score()), (6, 12));
        d.distributed = Some(DistributedEvidence::default());
        assert_eq!((d.lossless_lmer(), d.lossless_min_score()), (4, 12));
    }
}
