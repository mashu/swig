//! Declarative engine configuration and construction of the default assigner.
//!
//! Configuration is plain data; [`crate::profiles::CallingProfile`] produces
//! the note's calling profiles as values of [`EngineConfig`]. Settings that
//! can be *derived* losslessly from the admission gates (the D L-mer length
//! and the D traceback score floor) are derived here, so a profile change can
//! never make the prefilters silently drop admissible D hits.

use std::fmt;
use std::sync::Arc;

use crate::adjudicate::{
    CallModel, DAlleleDecoding, EngineeringObjective, EvidenceModel, PlacementDag, VdjGrammar,
};
use crate::admission::{AnyOf, DefaultAdmission, DistributedEvidence, ExactRunFloor};
use crate::align::{Backend, InterSequenceGotoh};
use crate::constant::ConstantStage;
use crate::family_bound::FamilyBoundTable;
use crate::flank::OuterExtend;
use crate::germline::{GermlineSet, SegmentKind};
use crate::index::IndexError;
use crate::locus::{Loci, Locus};
use crate::orient::StrandIndex;
use crate::pipeline::{Assigner, AssignerOptions};
use crate::retrieve::{
    AmbiguityBudget, ExhaustiveRetriever, LmerMaskRetriever, MAX_LMER, SeedVoteParams,
    SeedVoteRetriever, UngappedParams, UngappedPrefilter,
};
use crate::scoring::{AffineScoring, Score};
use crate::seed::SeedShape;
use crate::segment::{
    EvalueGates, HitAcceptance, KeepPolicy, PreparedSegment, SegmentStage, Significance, TieOrder,
};
use crate::stats::KarlinAltschul;
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
    /// E-value gates (`None` = no significance floor and no anchoring).
    pub evalue: Option<EvalueGates>,
    /// Skip DP for same-frame siblings that cannot enter the top-N or whose
    /// score their bounds pin (see [`crate::family_bound`]); lossless.
    pub family_delta: bool,
}

/// J: every allele, unbanded, in the read right of V.
#[derive(Clone, Debug, PartialEq)]
pub struct JSegmentConfig {
    pub scoring: AffineScoring,
    pub keep: KeepPolicy,
    pub accept: HitAcceptance,
    /// E-value gates (`None` = no significance floor and no anchoring).
    pub evalue: Option<EvalueGates>,
    pub window: DownstreamWindow,
    /// Sibling bounds as for V, each gene one frame; lossless.
    pub family_delta: bool,
}

/// C: splice-anchored alignment right of the called J ([`crate::constant`]).
#[derive(Clone, Debug, PartialEq)]
pub struct CSegmentConfig {
    pub scoring: AffineScoring,
    /// `s`: C's germline base 0 lies within `s` bases of the J exit.
    pub spread: u32,
    /// Call C iff its anchored E-value is at most this.
    pub evalue: f64,
    /// Anchored alignments kept per read. `min_score` is an anchored score;
    /// the stage adds the anchor bonus.
    pub keep: KeepPolicy,
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

    /// Traceback score floor of the D stage: the configured keep floor,
    /// raised to the lowest local score an admissible hit can have. An
    /// alignment holding an exact run of `k` scores at least `k·m` (the run
    /// alone is a local alignment); a distributed-D hit rescored under the D
    /// tuple itself scores at least `min_rescore`; under another tuple it
    /// only has a positive score. Raising the floor to that bound therefore
    /// drops nothing the gates would admit.
    pub fn lossless_min_score(&self) -> Score {
        let run = self.scoring.matched * self.min_exact_run as Score;
        let admissible = match &self.distributed {
            None => run,
            Some(gate) if gate.scoring == self.scoring => run.min(gate.min_rescore),
            Some(_) => 1,
        };
        self.keep.min_score.max(admissible)
    }
}

/// What the D search spans, for the D prior ([`DPrior::cost`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DSearchSpace {
    /// Genes of the D database.
    pub genes: usize,
    /// `n_D = Σ_a (L_a − k + 1)⁺`: the offsets of the D alleles the search
    /// aligns (each distinct sequence once) at which an exact run of
    /// `k = min_d_match` bases, the admission gate, can start.
    pub offsets: usize,
}

impl DSearchSpace {
    /// The search space of `d` under the exact-run gate `min_exact_run`.
    pub fn of(d: &GermlineSet, min_exact_run: u32) -> Self {
        let run = min_exact_run.max(1) as usize;
        Self {
            genes: d.gene_count(),
            offsets: d
                .alleles()
                .iter()
                .map(|a| (a.codes.len() + 1).saturating_sub(run))
                .sum(),
        }
    }
}

/// Prior cost of including a D, in score units (see [`EngineeringObjective`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DPrior {
    /// A fixed cost (0 = D only has to beat the no-D path).
    Fixed(Score),
    /// `P*_D = [ln(N · window) − ln(prior_odds)] / λ`: `N` D genes in the
    /// database, `window` junction positions, `λ` of the D substitution
    /// scores, and the prior odds of a D against none. For human IGH (32
    /// genes, 25 positions, even odds) this is 10.55; for TRB (2 D genes)
    /// 6.17. It is not rounded: the fixed-scale call compares integer scores
    /// with `⌊P*_D⌋`, which decides `S > P*_D` exactly.
    Derived { window: u32, prior_odds: f64 },
    /// The D search's own chance-hit count against a D-presence rate `π`:
    ///
    /// ```text
    /// P*_D = [ln(K · m · n_D) − ln(π / (1 − π))] / λ
    /// ```
    ///
    /// `K` and `λ` are the (lattice) Karlin–Altschul parameters of the D
    /// scoring, `m = window` read positions and `n_D` the database offsets
    /// at which an admissible hit can start ([`DSearchSpace`]). The search
    /// expects `E(S) = K·m·n_D·e^{−λS}` chance hits scoring `S` or more. A
    /// present D scores above the threshold almost surely, so the best hit
    /// is more likely a D than chance when `π·1 > (1 − π)·E(S)`, i.e. when
    /// `S > P*_D`. `Derived` is the same rule with `K·m·n_D` replaced by
    /// `N·L`, genes × positions, which undercounts a search that tries every
    /// allele and allele offset, and with `π = ½`. The two errors offset:
    /// at even odds the search's count
    /// puts `P*_D` at 14.2 for human IGH, and a measured `π` (0.9–0.99)
    /// lowers it by 3.5–7.3. Against a brute-force count of chance hits
    /// (random windows, one best hit per allele) `K·m·n_D·e^{−λS}` is high
    /// by a factor 1.3 at even `S` and 1.85 at odd `S`: exact runs put most
    /// chance hits on even scores. Read per read like every D cost, as
    /// `P*_D/2` matched D columns ([`crate::adjudicate`]).
    Search { window: u32, presence: f64 },
}

impl DPrior {
    pub const DERIVED: Self = Self::Derived {
        window: 25,
        prior_odds: 1.0,
    };

    /// [`DPrior::Search`] over the default 25 read positions at D-presence
    /// rate `presence`, in (0, 1).
    pub const fn search(presence: f64) -> Self {
        Self::Search {
            window: 25,
            presence,
        }
    }

    /// The cost for a D search over `space` scored with `scoring`. Falls back
    /// to `λ = 1` (nats per score unit) and `K = 1` when the tuple has no
    /// Karlin–Altschul statistics (it is then not a log-odds score).
    pub fn cost(&self, space: DSearchSpace, scoring: &AffineScoring) -> f64 {
        let ka = KarlinAltschul::ungapped(scoring);
        let lambda = ka.map_or(1.0, |ka| ka.lambda);
        let nats = match *self {
            Self::Fixed(cost) => return f64::from(cost),
            Self::Derived { window, prior_odds } => {
                (space.genes.max(1) as f64 * f64::from(window.max(1))).ln() - prior_odds.ln()
            }
            Self::Search { window, presence } => {
                let k = ka.map_or(1.0, |ka| ka.k);
                (k * f64::from(window.max(1)) * space.offsets.max(1) as f64).ln()
                    - (presence / (1.0 - presence)).ln()
            }
        };
        (nats / lambda).max(0.0)
    }
}

/// Declarative form of [`EngineeringObjective`]; `d_prior` is resolved
/// against the D database when the assigner is built.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObjectiveConfig {
    pub pair_bonus: Score,
    pub d_prior: DPrior,
    pub tandem_d_penalty: Score,
    pub insertion_per_base: Score,
}

impl Default for ObjectiveConfig {
    fn default() -> Self {
        Self {
            pair_bonus: 0,
            d_prior: DPrior::DERIVED,
            tandem_d_penalty: 12,
            insertion_per_base: 0,
        }
    }
}

impl ObjectiveConfig {
    pub fn resolve(
        &self,
        d: Option<&GermlineSet>,
        d_config: &DSegmentConfig,
    ) -> EngineeringObjective {
        let space = d.map_or(
            DSearchSpace {
                genes: 1,
                offsets: 1,
            },
            |d| DSearchSpace::of(d, d_config.min_exact_run),
        );
        EngineeringObjective {
            pair_bonus: self.pair_bonus,
            d_presence: self.d_prior.cost(space, &d_config.scoring),
            tandem_d_penalty: self.tandem_d_penalty,
            insertion_per_base: self.insertion_per_base,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EngineConfig {
    pub v: VSegmentConfig,
    pub d: DSegmentConfig,
    pub j: JSegmentConfig,
    /// Used by [`Assigner::with_constant`] when a C database is given.
    pub c: CSegmentConfig,
    pub d_window: ExpandedResidual,
    /// Add the conserved-anchor D window when alleles carry `anchor=` positions.
    pub anchor_window: bool,
    pub grammar: VdjGrammar,
    pub objective: ObjectiveConfig,
    /// Minimum retained query bases after clipping, per rearranged kind
    /// ([`SegmentKind::REARRANGED`]), indexed by [`SegmentKind::index`].
    pub min_retained: [u32; 3],
    /// The model the DAG's call maximises: each read's own (the model its
    /// posteriors come from) or the fixed scores.
    pub calls: CallModel,
    /// How the reported D allele is chosen: the call's, or the one of
    /// largest posterior marginal given the call's topology. The latter
    /// needs `posterior`; without it the D allele is the call's.
    pub d_allele: DAlleleDecoding,
    /// Report posterior probabilities (forward–backward over the placement
    /// DAG, priors at the V scoring's λ, alignment evidence at each read's
    /// divergence): per called segment and for D presence. With posteriors,
    /// or with a per-read call, D placements below the objective's floor
    /// (`P_D + 1`) are kept, since they carry probability and can win a
    /// per-read call; otherwise they are pruned, which never changes a
    /// fixed-scale call.
    pub posterior: bool,
    /// Instruction set of the alignment kernel; `None` picks the fastest
    /// this CPU runs ([`Backend::detect`]).
    pub simd: Option<Backend>,
    pub options: AssignerOptions,
}

impl Default for VSegmentConfig {
    /// Gap-open −9. Every 9-mer votes (7-mers too while the leader holds
    /// fewer than 25 seeds). The vote budget keeps 3 genes and walks up to 16
    /// while votes stay within two standard deviations of the leader's
    /// (`AmbiguityBudget`); the 3 best alignments that pass
    /// the E-value gate (and ties) are kept. There is no column or identity
    /// floor: a positive `(2, −3)` local alignment already has identity above
    /// 3/5, the E-value gate implies a minimum length, and the note's floors
    /// (24 columns, 60 %) changed no call (docs/ASSESSMENT.md §5).
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
            accept: HitAcceptance::PERMISSIVE,
            evalue: Some(EvalueGates::default()),
            family_delta: true,
        }
    }
}

impl Default for JSegmentConfig {
    /// `(2, −3, −17, −2)`; hits gated by E-value only, as for V.
    fn default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -17, -2),
            keep: KeepPolicy {
                top: 3,
                min_score: 1,
                cap: 64,
            },
            accept: HitAcceptance::PERMISSIVE,
            evalue: Some(EvalueGates::default()),
            window: DownstreamWindow::default(),
            family_delta: true,
        }
    }
}

impl Default for CSegmentConfig {
    /// The V scoring, so the gate and the posterior read V's λ. `s = 3`:
    /// reference conventions put C's base 0 within one base of the J exit
    /// (IMGT CH1 records after the J's last base, codon-aligned records on
    /// it, 3′-partial J records one base early), and SHM indels near the J's
    /// 3′ end misplace the exit by up to two more. The gate is the anchor
    /// gate of [`EvalueGates`]: C is called on its own evidence.
    fn default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -9, -1),
            spread: 3,
            evalue: EvalueGates::default().anchor,
            keep: KeepPolicy {
                top: 3,
                min_score: 1,
                cap: 64,
            },
        }
    }
}

impl Default for DSegmentConfig {
    fn default() -> Self {
        Self {
            scoring: AffineScoring::new(2, -3, -13, -1),
            lmer: None,
            // The 32 best admissible D hits reach adjudication; inadmissible
            // ones are traced (up to 64 tracebacks) but take no slot.
            keep: KeepPolicy {
                top: 32,
                min_score: 0,
                cap: 64,
            },
            min_exact_run: 4,
            // The note's distributed-D gate admits nothing the exact run and
            // the objective's D floor do not (docs/ASSESSMENT.md §5).
            distributed: None,
        }
    }
}

impl Default for EngineConfig {
    fn default() -> Self {
        Self {
            v: VSegmentConfig::default(),
            d: DSegmentConfig::default(),
            j: JSegmentConfig::default(),
            c: CSegmentConfig::default(),
            d_window: ExpandedResidual::default(),
            anchor_window: true,
            grammar: VdjGrammar::VDJ,
            objective: ObjectiveConfig::default(),
            // A clipped V or J keeps at least one base, a D `min_d_match`.
            min_retained: [1, 4, 1],
            calls: CallModel::default(),
            d_allele: DAlleleDecoding::default(),
            posterior: true,
            simd: None,
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
    /// The segment's scoring has no Karlin–Altschul `λ`.
    NoStatistics(SegmentKind),
    /// The requested SIMD backend does not run on this CPU.
    Unsupported(Backend),
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
            Self::NoStatistics(kind) => write!(
                f,
                "the {} scoring has no Karlin–Altschul statistics (its expected score must be negative)",
                kind.label()
            ),
            Self::Unsupported(backend) => {
                write!(f, "this CPU does not support the {} kernel", backend.name())
            }
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
    let significance = config
        .evalue
        .and_then(|gates| Significance::new(&set, &config.scoring, gates));
    let segment = Arc::new(PreparedSegment::new(
        set,
        &config.scoring,
        config.scoring.gaps,
    ));
    // Bounds are sound for bands shifted through the frames the bands come from.
    let family_delta = config
        .family_delta
        .then(|| FamilyBoundTable::with_frames(&segment, seeds.frames()));
    let retriever = UngappedPrefilter::new(seeds, Arc::clone(&segment), config.ungapped);
    let mut stage = SegmentStage::new(segment, retriever, config.keep);
    stage.accept = config.accept;
    stage.significance = significance;
    stage.outer = OuterExtend::Prefix;
    // The order the V stage has always produced: gene leaders were scored,
    // and so ranked, before their siblings.
    stage.ties = TieOrder::GeneLeadersFirst;
    stage.family_delta = family_delta;
    Ok(stage)
}

fn j_stage(set: GermlineSet, config: &JSegmentConfig) -> SegmentStage<ExhaustiveRetriever> {
    let retriever = ExhaustiveRetriever { alleles: set.len() };
    let significance = config
        .evalue
        .and_then(|gates| Significance::new(&set, &config.scoring, gates));
    let segment = PreparedSegment::new(set, &config.scoring, config.scoring.gaps);
    // J bands are unbanded, which any placement of the alleles bounds.
    let family_delta = config
        .family_delta
        .then(|| FamilyBoundTable::build(&segment));
    let mut stage = SegmentStage::new(Arc::new(segment), retriever, config.keep);
    stage.accept = config.accept;
    stage.significance = significance;
    stage.outer = OuterExtend::Suffix;
    // Ties in retrieval (allele id) order, as J has always ranked them.
    stage.ties = TieOrder::Retrieval;
    stage.family_delta = family_delta;
    stage
}

fn d_stage(
    set: GermlineSet,
    config: &DSegmentConfig,
    objective: &EngineeringObjective,
    keep_below_floor: bool,
) -> SegmentStage<LmerMaskRetriever, DefaultAdmission> {
    let lmer = config.lmer.unwrap_or_else(|| config.lossless_lmer());
    let retriever = LmerMaskRetriever::build(&set, lmer, 1);
    let segment = PreparedSegment::new(set, &config.scoring, config.scoring.gaps);
    let keep = KeepPolicy {
        // The objective's D floor prunes only what cannot win a fixed-scale
        // call (see `EngineeringObjective::d_floor`).
        min_score: if keep_below_floor {
            config.lossless_min_score()
        } else {
            config.lossless_min_score().max(objective.d_floor())
        },
        ..config.keep
    };
    let admission = AnyOf(
        ExactRunFloor {
            min_run: config.min_exact_run,
        },
        config.distributed,
    );
    SegmentStage::new(Arc::new(segment), retriever, keep).with_accept(admission)
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
        let objective = config.objective.resolve(d.as_ref(), &config.d);
        let substitution = |s: &AffineScoring| (s.matched, s.mismatch);
        let evidence = KarlinAltschul::ungapped(&config.v.scoring).map(|ka| EvidenceModel {
            lambda: ka.lambda,
            substitution: [
                substitution(&config.v.scoring),
                substitution(&config.d.scoring),
                substitution(&config.j.scoring),
            ],
        });
        let adjudicator = PlacementDag {
            evidence,
            calls: config.calls,
            posterior: config.posterior && evidence.is_some(),
            d_allele: config.d_allele,
            ..PlacementDag::new(grammar, objective, config.min_retained)
        };
        // A D database of several loci: a read of one locus prices a D by
        // that locus' D genes.
        let d_loci: Vec<Locus> = d.as_ref().map_or(Vec::new(), |d| {
            Locus::ALL
                .into_iter()
                .filter(|&l| d.alleles().iter().any(|a| a.loci == Loci::of(l)))
                .collect()
        });
        let locus_adjudicators: Vec<(Locus, PlacementDag)> = match &d {
            Some(d) if d_loci.len() > 1 => d_loci
                .iter()
                .map(|&locus| {
                    let (part, _) = d.subset(|a| a.loci.contains(locus));
                    let model = config.objective.resolve(Some(&part), &config.d);
                    (
                        locus,
                        PlacementDag {
                            model,
                            ..adjudicator
                        },
                    )
                })
                .collect(),
            _ => Vec::new(),
        };
        // The D floor must not prune a D any locus would call.
        let floor = locus_adjudicators
            .iter()
            .map(|(_, a)| a.model)
            .fold(objective, |low, m| {
                if m.d_presence < low.d_presence {
                    m
                } else {
                    low
                }
            });
        Ok(Assigner {
            v: v_stage(v, &config.v)?,
            j: j_stage(j, &config.j),
            // A D below the objective's floor carries posterior mass, and can
            // win a per-read call.
            d: d.map(|set| {
                let keep = adjudicator.posterior || !adjudicator.calls_on_fixed_scale();
                d_stage(set, &config.d, &floor, keep)
            }),
            strand,
            j_window: config.j.window,
            aligner: match config.simd {
                None => InterSequenceGotoh::new(),
                Some(backend) => InterSequenceGotoh::with_backend(backend)
                    .ok_or(BuildError::Unsupported(backend))?,
            },
            d_window,
            adjudicator,
            locus_adjudicators,
            c: None,
            options: config.options,
        })
    }
}

impl<A, RV, RJ, RD, W, Adm, Adj> Assigner<A, RV, RJ, RD, W, Adm, Adj> {
    /// Add constant-region assignment against C database `c`
    /// ([`crate::constant`]); V, D and J calls are unchanged.
    pub fn with_constant(
        mut self,
        c: GermlineSet,
        config: &CSegmentConfig,
    ) -> Result<Self, BuildError> {
        check(&c, SegmentKind::C)?;
        let stage =
            ConstantStage::new(c, config).ok_or(BuildError::NoStatistics(SegmentKind::C))?;
        self.c = Some(stage);
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adjudicate::PathModel;

    /// A search space of `genes` genes and as many offsets.
    fn genes(genes: usize) -> DSearchSpace {
        DSearchSpace {
            genes,
            offsets: genes,
        }
    }

    #[test]
    fn d_prior_scales_with_the_d_database() {
        let scoring = AffineScoring::new(2, -3, -13, -1);
        let cost = |prior: DPrior, space| (prior.cost(space, &scoring) * 100.0).round() / 100.0;
        assert_eq!(cost(DPrior::DERIVED, genes(32)), 10.55); // human IGH
        assert_eq!(cost(DPrior::DERIVED, genes(2)), 6.17); // TRB
        assert_eq!(cost(DPrior::Fixed(0), genes(32)), 0.0);
        let more_likely = DPrior::Derived {
            window: 25,
            prior_odds: 4.0,
        };
        assert_eq!(cost(more_likely, genes(32)), 8.36);
        // The fixed-scale call charges ⌊P*_D⌋: an integer D score must exceed it.
        let objective = EngineeringObjective {
            d_presence: DPrior::DERIVED.cost(genes(32), &scoring),
            ..EngineeringObjective::default()
        };
        assert_eq!(
            (objective.d_presence_score(), objective.d_floor()),
            (10, 11)
        );
    }

    /// The 20,000-read human IGH database: 30 genes, 37 distinct alleles
    /// with 802 offsets for an exact run of 4.
    const HUMAN: DSearchSpace = DSearchSpace {
        genes: 30,
        offsets: 802,
    };

    #[test]
    fn search_prior_counts_the_d_search() {
        let scoring = AffineScoring::new(2, -3, -13, -1);
        let ka = KarlinAltschul::ungapped(&scoring).unwrap();
        let cost = |presence| DPrior::search(presence).cost(HUMAN, &scoring);
        // At even odds the threshold is where one chance hit is expected.
        let even = cost(0.5);
        assert!((ka.k * 25.0 * 802.0 * (-ka.lambda * even).exp() - 1.0).abs() < 1e-9);
        assert_eq!((even * 100.0).round(), 1422.0);
        // At presence π, π/(1 − π) chance hits.
        assert_eq!((cost(0.9) * 100.0).round(), 1075.0);
        assert!((cost(0.9) - (even - 9f64.ln() / ka.lambda)).abs() < 1e-9);
        // Derived counts genes × positions: the same rule with K·n_D = N.
        let derived = DPrior::DERIVED.cost(HUMAN, &scoring);
        let as_genes = DSearchSpace {
            offsets: (30.0 / ka.k).round() as usize,
            ..HUMAN
        };
        assert!((DPrior::search(0.5).cost(as_genes, &scoring) - derived).abs() < 0.01);
    }

    #[test]
    fn search_space_counts_gate_offsets_of_distinct_alleles() {
        let text = ">IGHD1-1*01\nACGTACGTAC\n>IGHD1-1*02\nACGTACGTAA\n>IGHD2-2*01\nACG\n";
        let d = GermlineSet::from_fasta(SegmentKind::D, text).unwrap();
        // 10 − 4 + 1 offsets for each 10-mer; the 3-mer holds no run of 4.
        assert_eq!(
            DSearchSpace::of(&d, 4),
            DSearchSpace {
                genes: 2,
                offsets: 14
            }
        );
    }

    /// Alleles (each of the D database) with a chance hit scoring at
    /// least `score` against `window`: the best ungapped local alignment, per
    /// allele as the search keeps it, among those holding an exact run of
    /// `run` (the admission gate).
    fn chance_hits(
        alleles: &[Vec<u8>],
        window: &[u8],
        scoring: &AffineScoring,
        run: u32,
        score: Score,
    ) -> usize {
        let best = |allele: &[u8]| {
            let mut best = 0;
            for shift in -(allele.len() as isize - 1)..window.len() as isize {
                let (mut sum, mut exact, mut longest) = (0, 0, 0);
                for (i, &base) in allele.iter().enumerate() {
                    let Some(&read) = usize::try_from(i as isize + shift)
                        .ok()
                        .and_then(|q| window.get(q))
                    else {
                        continue;
                    };
                    (sum, exact) = if base == read {
                        (sum + scoring.matched, exact + 1)
                    } else {
                        (sum + scoring.mismatch, 0)
                    };
                    if sum <= 0 {
                        (sum, longest) = (0, 0);
                        continue;
                    }
                    longest = longest.max(exact);
                    if longest >= run {
                        best = best.max(sum);
                    }
                }
            }
            best
        };
        alleles.iter().filter(|a| best(a) >= score).count()
    }

    /// The closed form against a brute-force count of chance hits of random
    /// D-like alleles (lengths 12–31) in random 25-base windows: it brackets
    /// the count from above within a factor of two around the threshold.
    #[test]
    fn search_count_matches_brute_force_chance_hits() {
        let scoring = AffineScoring::new(2, -3, -13, -1);
        let ka = KarlinAltschul::ungapped(&scoring).unwrap();
        let mut rng = crate::sim::Rng::new(7);
        let mut random = |len: usize| (0..len).map(|_| rng.below(4) as u8).collect::<Vec<_>>();
        let alleles: Vec<Vec<u8>> = (0..20).map(|i| random(12 + i)).collect();
        let offsets: usize = alleles.iter().map(|a| a.len() - 3).sum();
        let (window, trials) = (25, 4000);
        let windows: Vec<Vec<u8>> = (0..trials).map(|_| random(window)).collect();
        for score in [10, 11, 12] {
            let brute = windows
                .iter()
                .map(|w| chance_hits(&alleles, w, &scoring, 4, score))
                .sum::<usize>() as f64
                / trials as f64;
            let closed =
                ka.k * window as f64 * offsets as f64 * (-ka.lambda * f64::from(score)).exp();
            assert!(
                brute < closed && closed < 2.0 * brute,
                "S = {score}: brute force {brute:.3}, closed form {closed:.3}"
            );
        }
    }

    #[test]
    fn derived_d_settings_are_lossless_bounds() {
        let mut d = DSegmentConfig::default();
        assert_eq!((d.lossless_lmer(), d.lossless_min_score()), (4, 8));
        d.min_exact_run = 6;
        assert_eq!((d.lossless_lmer(), d.lossless_min_score()), (6, 12));
        d.distributed = Some(DistributedEvidence::default());
        assert_eq!((d.lossless_lmer(), d.lossless_min_score()), (4, 12));
        // A distributed gate under another tuple bounds nothing.
        d.distributed = Some(DistributedEvidence {
            scoring: AffineScoring::new(1, -1, -5, -1),
            ..DistributedEvidence::default()
        });
        assert_eq!(d.lossless_min_score(), 1);
        // A configured floor above the bound is kept.
        d.keep.min_score = 15;
        assert_eq!(d.lossless_min_score(), 15);
    }
}
