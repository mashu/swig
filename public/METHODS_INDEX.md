# Swig methods index

This index is the implementation-facing methods map for Swig 0.39.3. It covers every runnable analysis block from file intake through phylogenetic UCA inference. The UI may summarize a method in one sentence; the linked document is the specification.

## Evidence labels used in these documents

| Label | Meaning |
|---|---|
| **Direct implementation** | Swig implements the stated equations or algorithm itself. It may still use different storage or parallelization. |
| **Browser port** | The scientific core is deliberately matched to named source code or a paper, with listed engineering differences. |
| **Compatible** | Swig reproduces a defined decision rule or output convention, but is not the named software. |
| **Custom** | The method is Swig-specific. Papers are context, validation targets, or sources of parameter values—not claims of identity. |

## Interface-to-method audit

| Interface block | Complete method document |
|---|---|
| Datasets, concatenated gzip choice, study names, FASTA/FASTQ/AIRR parsing, FASTQ QC, random subsampling | [Input, QC, and study structure](methods/01_INPUT_QC_AND_STUDY.md) |
| Biological search space, reference matrix, online/local reference preparation, allele exclusion | [Reference preparation](methods/02_REFERENCE_PREPARATION.md) |
| Swig aligner, calling profile, strand, identity floor, C gate, segment E-values and posteriors, workers, AIRR fields, record detail | [V(D)J assignment](methods/03_VDJ_ASSIGNMENT.md) |
| Optional two-D/VDDJ screen and Double-D explorer | [Double-D evidence screen](methods/04_DOUBLE_D_SCREEN.md) |
| AIRR result index, repertoire dashboard, sequence filters, metadata re-indexing, cumulative working set, repertoire selection | [Storage, summaries, and selection](methods/05_STORAGE_DASHBOARD_SELECTION.md) |
| Resolve ambiguous germline calls by pooling repertoire evidence | [Repertoire allele refinement](REPERTOIRE_ALLELE_REFINEMENT.md) |
| Exact collapse, FAD-compatible denoising, conservative error model, indel-aware error model | [Collapse and denoising](methods/06_COLLAPSE_AND_DENOISING.md) |
| CHMMAIRRa fitting, posterior filter, detailed Viterbi reconstruction | [Chimera inference](methods/07_CHIMERA_INFERENCE.md) |
| CDR3-based lineage assignment, lineage table, explicit merges | [Lineage assignment](methods/08_LINEAGE_ASSIGNMENT.md) |
| SHM distributions, two-pass possible-missing-V screen, and personalized expressed V-set inference | [Post-lineage SHM and germline evidence](methods/09_SHM_AND_GERMLINE_DIAGNOSTICS.md) |
| Targeted sequence search, query V/J inference, fixed-point expansion, lineage neighbours | [Retrieval and neighbour search](methods/10_RETRIEVAL_AND_NEIGHBOURS.md) |
| Closest-member/consensus lineage root, quick/Kalign alignment, Alivibe round trip, FastTree and rooting | [Alignment and phylogeny](methods/11_ALIGNMENT_AND_PHYLOGENY.md) |
| Fixed-tree UCA model, HMM, ML/grid/Gibbs-MH placement, posterior logos/tracks/exports | [Phylogenetic UCA inference](PHYLO_UCA_INFERENCE.md) |
| Interactive versus pipeline ordering, sessions/projects, exports, cancellation and rollback | [Execution, persistence, and cancellation](methods/12_EXECUTION_PERSISTENCE_CANCELLATION.md) |
| Portable CLI configuration, standalone execution, and lazy lineage-study sessions | [CLI and lineage-study interchange](methods/13_CLI_AND_LINEAGE_STUDY.md) |

## Important scope distinctions

- V(D)J assignment is per rearranged sequence. Repertoire allele refinement is the optional first cross-record model.
- Applying a stage changes a typed downstream mask or a call overlay; it never deletes the immutable assigned AIRR rows.
- Default allele-model fits are within donor and may combine samples from that donor. Cross-donor scopes are explicit overrides.
- Collapse and lineage scopes are independent hard candidate-generation boundaries.
- The ordinary rooted tree and the phylogenetic UCA calculation are distinct. UCA inference removes the guide, fixes a tree of observed lineage sequences, then places an unobserved ancestor against that tree.

## Reproducibility

Saved sessions retain all settings needed to reproduce applied downstream state and store compact masks/vectors rather than a second AIRR table. Machine-readable outputs preserve original calls alongside any policy-selected calls. Defaults are starting values, not universal biological constants; each method document identifies the settings that require assay-, organism-, or dataset-specific sensitivity checks.


## Separate experimental joint germline method (0.38.9)

[Joint inherited / SHM model](methods/10_JOINT_INHERITED_SHM_MODEL.md) defines its single likelihood, split-data evidence, numerical bounds and limitations. It runs alongside the existing personalized-reference box.
