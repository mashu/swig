//! Declarative engine configuration and construction of the default assigner.
//!
//! Configuration is plain data; [`crate::profiles::CallingProfile`] produces
//! the note's calling profiles as values of [`EngineConfig`]. Settings that
//! can be *derived* losslessly from the admission gates (the D L-mer length
//! and the D traceback score floor) are derived here, so a profile change can
//! never make the prefilters silently drop admissible D hits.

use std::fmt;
use std::sync::Arc;

use crate::adjudicate::{DPresencePenalty, EngineeringObjective, PlacementDag, VdjGrammar};
use crate::admission::{AnyOf, DistributedEvidence, ExactRunFloor};
use crate::align::InterSequenceGotoh;
use crate::germline::{GermlineSet, SegmentKind};
use crate::index::IndexError;
use crate::pipeline::{Assigner, AssignerOptions};
use crate::refine::{NoRefine, OuterExtend};
use crate::retrieve::{
    AmbiguityBudget, LmerMaskRetriever, SeedVoteParams, SeedVoteRetriever, UngappedParams,
    UngappedPrefilter, MAX_LMER,
};
use crate::scoring::{AffineScoring, Score};
use crate::seed::SeedShape;
use crate::segment::{HitAcceptance, KeepPolicy, PreparedSegment, SegmentStage};
use crate::window::{ExpandedResidual, WithAnchorWindow};

/// V or J: seeded, banded, vote-ranked, optionally ungapped-pre-ranked.
#[derive(Clone, Debug, PartialEq)]
pub struct VjSegmentConfig {
    pub scoring: AffineScoring,
    pub primary: SeedShape,
    pub fallback: Option<SeedShape>,
    pub seeds: SeedVoteParams,
    pub budget: AmbiguityBudget,
    /// Ungapped diagonal pre-ranking before banded DP (`None` = off).
    /// Off by default: it is lossy, and Swig aligns the seed candidates directly.
    pub ungapped: Option<UngappedParams>,
    pub keep: KeepPolicy,
    /// Co-linear flank stitched back onto the local alignment before ranking.
    pub outer: OuterExtend,
    /// Minimum substitution columns and identity of an accepted local hit.
    pub accept: HitAcceptance,
    /// Softer acceptance used while placing, then discarded when the best V
    /// identity is above [`crate::refine::JUNCTION_TIE_IDENTITY`]. `None` =
    /// always use [`Self::accept`]. For J under SHM: retrieve/score soft, keep
    /// strict gates on clean (high-identity) V reads.
    pub accept_shm: Option<HitAcceptance>,
    /// Skip score-only DP for same-gene siblings that cannot enter the top-N
    /// local scores, using a sound upper bound from the gene's representative.
    pub family_delta: bool,
    /// Build substitution profiles with [`crate::scoring::HotspotScoring`]
    /// (AID hotspot / coldspot mismatch costs). Intended for J under SHM.
    pub hotspot: bool,
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
    pub v: VjSegmentConfig,
    pub d: DSegmentConfig,
    pub j: VjSegmentConfig,
    pub d_window: ExpandedResidual,
    /// Add the conserved-anchor D window when alleles carry `anchor=` positions.
    pub anchor_window: bool,
    pub grammar: VdjGrammar,
    pub objective: EngineeringObjective,
    /// Minimum retained query bases after clipping, per [`SegmentKind::index`].
    pub min_retained: [u32; 3],
    pub options: AssignerOptions,
}

impl VjSegmentConfig {
    /// Swig R-optimized V: gap-open −9, at least 24 substitution columns, identity ≥ 0.60.
    /// The vote budget keeps 3 genes and walks up to 16 while votes stay within the AER margin.
    pub fn v_default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -9, -1),
            primary: SeedShape::contiguous(9),
            fallback: Some(SeedShape::contiguous(7)),
            seeds: SeedVoteParams::default(),
            budget: AmbiguityBudget::default(),
            ungapped: None,
            // Score-screen to Swig's top_v of 3. Traceback stops once three
            // alignments pass acceptance, including the tie group at that score.
            keep: KeepPolicy {
                top: 3,
                slack: 0,
                min_score: 1,
                cap: 64,
                prescore_above: 4,
                screen: true,
            },
            outer: OuterExtend::Prefix,
            accept: HitAcceptance {
                min_subst: 24,
                min_identity: 0.60,
            },
            accept_shm: None,
            family_delta: true,
            hotspot: false,
        }
    }

    /// J: `(2, −3, −17, −2)`. Primary is a 5-mer (not 7): the locus is ~50 bp
    /// with few alleles, and under SHM a 7-mer primary can lock onto decoys and
    /// skip the short-seed fallback. No shorter fallback — 4-mers are too shared
    /// across the tiny J set and drown the true allele. V keeps longer seeds.
    /// Profiles use SHM hotspot mismatch costs so the true J can win under
    /// heavy mutation when decoys disagree on cold/neutral positions.
    /// Clean reads keep Swig-like acceptance; high-SHM V reads keep the soft gate.
    pub fn j_default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -17, -2),
            primary: SeedShape::contiguous(5),
            fallback: None,
            seeds: SeedVoteParams::default(),
            budget: AmbiguityBudget {
                base: 3,
                max: 5,
                ..AmbiguityBudget::default()
            },
            ungapped: None,
            // Score-screen a few more J candidates under SHM ambiguity.
            keep: KeepPolicy {
                top: 3,
                slack: 0,
                min_score: 1,
                cap: 64,
                prescore_above: 3,
                screen: true,
            },
            outer: OuterExtend::Suffix,
            accept: HitAcceptance {
                min_subst: 10,
                min_identity: 0.60,
            },
            accept_shm: Some(HitAcceptance {
                min_subst: 8,
                min_identity: 0.50,
            }),
            // J genes rarely have enough alleles for the bound to pay off.
            family_delta: false,
            hotspot: true,
        }
    }
}

impl Default for DSegmentConfig {
    fn default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -13, -1),
            lmer: None,
            keep: KeepPolicy {
                top: 16,
                slack: Score::MAX / 4,
                min_score: 0,
                cap: 32,
                prescore_above: 0,
                screen: false,
            },
            min_exact_run: 4,
            distributed: Some(crate::admission::DistributedEvidence::default()),
        }
    }
}

impl Default for EngineConfig {
    fn default() -> Self {
        Self {
            v: VjSegmentConfig::v_default(),
            d: DSegmentConfig::default(),
            j: VjSegmentConfig::j_default(),
            d_window: ExpandedResidual::default(),
            anchor_window: true,
            grammar: VdjGrammar::VDJ,
            objective: EngineeringObjective {
                pair_bonus: 24,
                d_presence: DPresencePenalty::R_OPTIMIZED,
                tandem_d_penalty: 12,
                insertion_per_base: 0,
                d_short_run_penalty: 1,
                d_short_run_floor: 6,
                d_long_span_soft_max: 0,
                d_long_span_penalty: 0,
            },
            // Swig: V ≥ 24 aligned bases, J ≥ 10, D ≥ min_d_match. Query bases
            // retained after a clip use the same floors.
            min_retained: [24, 4, 10],
            options: AssignerOptions {
                junction_tie: true,
                ..AssignerOptions::default()
            },
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

/// Retriever of the default V and J stages.
pub type DefaultVjRetriever = UngappedPrefilter<SeedVoteRetriever>;

fn vj_stage(
    set: GermlineSet,
    config: &VjSegmentConfig,
) -> Result<SegmentStage<DefaultVjRetriever>, BuildError> {
    let kind = set.kind();
    let seeds = SeedVoteRetriever::build(
        &set,
        config.primary.clone(),
        config.fallback.clone(),
        config.seeds,
        config.budget,
    )
    .map_err(|source| BuildError::Index { kind, source })?;
    let segment = Arc::new(if config.hotspot {
        PreparedSegment::new(
            set,
            &crate::scoring::HotspotScoring::for_shm(config.scoring),
            config.scoring.gaps,
        )
    } else {
        PreparedSegment::new(set, &config.scoring, config.scoring.gaps)
    });
    let retriever = UngappedPrefilter::new(seeds, Arc::clone(&segment), config.ungapped);
    let family_delta = config
        .family_delta
        .then(|| crate::family_bound::FamilyBoundTable::build(segment.set()));
    let mut stage = SegmentStage::new(segment, retriever, config.keep);
    // Place with the soft SHM gate when configured; clean-V filtering restores
    // `accept` after placement (see pipeline).
    stage.accept = config.accept_shm.unwrap_or(config.accept);
    stage.accept_clean = config.accept_shm.map(|_| config.accept);
    stage.outer = config.outer;
    stage.family_delta = family_delta;
    Ok(stage)
}

fn d_stage(set: GermlineSet, config: &DSegmentConfig) -> SegmentStage<LmerMaskRetriever> {
    let lmer = config.lmer.unwrap_or_else(|| config.lossless_lmer());
    let retriever = LmerMaskRetriever::build(&set, lmer, 1);
    let segment = PreparedSegment::new(set, &config.scoring, config.scoring.gaps);
    let floor = config.lossless_min_score();
    // A shared L-mer is an exact run of L, so every retrieved candidate already
    // scores ≥ L·m. When that meets the floor, the score-only pre-pass can only
    // prune through `top`, so it runs only above `top` candidates.
    let floor_implied = lmer as Score * config.scoring.matched >= floor;
    let prescore_above = if floor_implied {
        config.keep.prescore_above.max(config.keep.top)
    } else {
        config.keep.prescore_above
    };
    let keep = KeepPolicy {
        min_score: floor,
        prescore_above,
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
/// available), seed-vote + ungapped V/J retrieval,
/// L-mer D prefilter, expanded residual + anchor D windows,
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
        Ok(Assigner {
            v: vj_stage(v, &config.v)?,
            j: vj_stage(j, &config.j)?,
            d: d.map(|set| d_stage(set, &config.d)),
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
