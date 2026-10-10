# Repertoire-level germline-call refinement

## Purpose

V(D)J assignment is normally performed one sequence at a time. With a deliberately overcomplete reference database, a read may not distinguish several related alleles even when the repertoire collectively supports only a subset. Swig's optional allele-pooling step combines:

1. local alignment evidence retained for each read; and
2. the allele-usage distribution inferred from other active records in the same repertoire pool.

The method refines assignments. It does **not** certify a diploid genotype, infer genomic copy number, discover an allele absent from the reference set, or overwrite the original AIRR result.

## Independent pools

Models are fitted independently by user-selected study boundary (donor/subject by default), receptor locus, and segment. Under the default, every sample, timepoint, and compartment carrying the same donor ID contributes to the same fit, but evidence never crosses donor IDs. V and J are enabled by default. D is available as an experimental option because short templated spans, trimming, and N addition make D emissions much less identifiable.

No posterior information crosses these boundaries. Donor pooling is usually appropriate for longitudinal or compartmental data; dataset and sample scopes are available when fits must remain narrower. Cohort and entire-study scopes are explicitly labeled **cross-donor overrides** because they deliberately allow information to cross participant IDs.

## Reference graph and local evidence

Reference records with identical nucleotide sequences are collapsed into one unresolved node carrying all equivalent labels. The remaining nodes are connected within a small bounded nucleotide edit radius. This graph constructs sparse candidate support; Swig never assigns every read nonzero evidence against the complete database.

Let \(E_{ra}\) denote local evidence for read \(r\) and reference node \(a\).

- Every literal co-optimal AIRR call starts at weight 1. A two-way explicit tie is therefore 50:50 before repertoire pooling.
- A retained alternative with alignment score \(S_a\) receives the weight below. In rows written by the Swig aligner, the retained alternatives are the entries of `*_alternatives`: alleles whose alignment, clipped to the called interval, scores within \(\max(2, \lfloor|S|/100\rfloor)\) of the call's score \(S\).

  \[
  \min\left(1,\exp\left[\frac{S_a-S_{\max}}{T}\right]\right),
  \]

  where \(T\) is an adjustable score temperature.
- An alternative without a numeric score receives a separate fallback weight.

These values form an evidence kernel. The Swig aligner's raw alignment scores are not claimed to be calibrated read likelihoods, and the kernel does not read the aligner's per-read `*_posterior` columns.

### SHM-adaptive neighbour leakage

An unreported nearby allele should not receive an artificial exact zero merely because it fell outside the caller's retained hypothesis set. For a same-length one-substitution neighbour, Swig uses per-edit evidence odds

\[
o_1(\mu)=\min\left(o_{\max},\;o_0+s\frac{\mu}{3(1-\mu)}\right).
\]

- \(o_0\): zero-SHM neighbour odds (default 0.01). This is irreducible assignment-model leakage, **not** a sequencing-error estimate.
- \(\mu\): read-level SHM estimate, normally one minus the best retained V-reference identity, clamped at a configurable maximum.
- \(1/3\): equal-substitution probability of landing on the particular nucleotide distinguishing the neighbouring allele.
- \(s\): SHM-sensitivity multiplier (default 1; zero disables SHM adaptation).
- \(o_{\max}\): conservative cap for highly mutated or poorly matched reads.

The term \(\mu/[3(1-\mu)]\) is the likelihood ratio of a particular diagnostic substitution to retention under a context-averaged equal-substitution approximation. It is intentionally simple and inspectable, not a five-mer or other context-dependent model.

A substitution-only neighbour at distance \(d\) receives relative leakage \(o_1(\mu)^d\). Reference neighbours involving an indel receive the non-SHM baseline only because the substitution derivation does not apply. Direct, alternative, and neighbour weights are combined, capped to a maximum sparse row width, and normalized.

### Interactive evidence-kernel inspection

The advanced settings include a reference-sequence inspector so parameter changes can be interpreted before fitting the repertoire model. The user chooses a V, D, or J database allele and an assumed read SHM fraction. Swig then evaluates the exact sparse evidence-row constructor used by the worker with that allele as the sole literal primary call.

- The selected reference sequence remains visible as the alignment anchor, but its dominant probability is deliberately not drawn as a bar.
- Every non-primary candidate receives a horizontal bar and its absolute probability after normalization over the complete row, including the primary.
- Alternative bar lengths are scaled relative to the largest non-primary candidate. The printed percentage, rather than bar length, is the absolute probability.
- Candidate sequences are globally aligned to the selected reference through a common anchor coordinate. Columns that are gaps in every displayed sequence are removed.
- Nucleotide coloring shows ordinary bases. Difference-highlighter coloring makes matches neutral, substitutions red, and insertion/deletion columns amber.
- Sequence-identical database labels remain one unresolved reference node and are displayed together.

The diagnostic responds to the zero-SHM floor, SHM sensitivity and cap, assumed SHM, edit radius, maximum neighbour odds, and candidate cap. Dirichlet alpha, hurdle-model presence/tail controls, and fitting stopping controls do not alter this local pre-repertoire kernel. Numeric settings commit on Enter or when focus leaves the field, matching the rest of Swig's large-data controls; the inexpensive diagnostic then updates immediately.

## Selectable repertoire models

The model selector leaves the original continuous Dirichlet implementation intact and adds a fast hurdle active-set alternative. Both models consume exactly the same sparse local-evidence matrix, pool at the same donor/study boundary, and produce the same per-record MAP probability interface. Reassignment policies and every hard-count/export diagnostic therefore work identically after either fit.

## Dirichlet mixture

For one independent pool, let \(\boldsymbol\theta\) be expressed reference usage:

\[
\boldsymbol\theta \sim \operatorname{Dirichlet}(\boldsymbol\alpha),
\qquad
p(z_r=a\mid\boldsymbol\theta,E_r)\propto\theta_aE_{ra}.
\]

Mean-field coordinate ascent uses

\[
q(z_r=a) \propto E_{ra}\exp\left(\psi(\gamma_a)-\psi\left(\sum_b\gamma_b\right)\right)
\]

and

\[
\gamma_a=\alpha_a+\sum_rw_rq(z_r=a).
\]

The default \(w_r=1\) gives every active record one vote. Optional abundance weighting uses `duplicate_count`, but allows PCR abundance and clonal expansion to influence the inferred mixture.

The symmetric prior covers every locus-matched reference node in the selected segment, including an allele absent from every read's sparse candidate neighbourhood. Such an inactive node has zero row likelihood, receives no assignment responsibility, and retains \(\gamma_a=\alpha_a\). Swig accounts for the combined prior mass of these inactive nodes in the Dirichlet normalization without materializing them once per donor or once per read. Nodes with unknown locus metadata are included only when they occur in an observed candidate row, because there is no defensible locus pool to which an otherwise inactive custom record can be assigned.

This resembles sparse variational machinery used for LDA, but the biological model is a finite mixture with one latent germline node per read, not a document containing multiple topics.

## Fast hurdle active set

The alternative model separates repertoire activity from positive usage. An excluded allele has exactly zero mixture weight; a retained allele has positive frequency under a long-tailed slab. For candidate active set (A),

\[
p(r\mid A,\boldsymbol\theta)=\sum_{a\in A}\theta_aE_{ra}.
\]

The configurable active prior is the per-database-allele probability of contributing measurable rearrangements in this donor/locus/segment pool. It is not a claim of genomic presence. Conditional positive frequency uses a discretized gamma slab whose default shape is below one, preserving substantial density in the low-frequency tail. The quadrature floor is numerical integration support, not a biological minimum frequency or read-count threshold.

Inference avoids subset enumeration:

1. sparse EM fits frequencies over the current active set;
2. each allele receives a one-dimensional conditional inclusion test obtained by removing its current weight and integrating a possible positive frequency over a small log-spaced quadrature;
3. alleles below the configurable inclusion-posterior threshold are pruned in parallel, except when removal would leave any observed sparse evidence pattern with no active explanation; and
4. frequencies are refitted over the smaller exact non-zero set.

For a candidate frequency (f), records outside that allele's sparse neighbourhood contribute one aggregated \(\log(1-f)\) term. Only evidence patterns containing the candidate are visited explicitly. Consequently a complete inclusion sweep remains proportional to the sparse matrix size rather than the number of reads times the complete database. The saved model and sidecar report the approximate inclusion probability and exact active/inactive state separately from the per-record MAP probability used by reassignment.

### Before/after hard-assignment figure

After fitting, Swig displays paired horizontal count bars for every allele receiving at least one hard assignment in a selected donor/study-boundary, locus, and segment pool. Each modeled record contributes its complete configured weight to exactly one reference node in each series:

- **Local best** is the argmax of that record's normalized sparse local-evidence row before repertoire pooling.
- **After policy** is the fitted posterior MAP node. Under the confidence-gated policy, a MAP node below the selected confidence threshold is held at the local-best node in this hard-count projection.

With unique-record weighting, a bar value is a record count. With abundance weighting, it is a `duplicate_count`-weighted assignment count. The bars are therefore neither mixture-frequency posterior means nor sums of soft responsibilities. Sequence-identical reference labels remain one unresolved node and are displayed together.

The figure can show all hard-assigned alleles, only alleles whose counts change, or only alleles whose local-best count falls to zero after the selected policy. The latter is the exact “vanishes” set for the hard projection. CSV export includes every row and explicit `vanishes`/`appears` flags even when the on-screen row limit is smaller; the visible figure exports as SVG. The compact fitted-model summary remains available as a download rather than a truncated table in the page. No short, thresholded allele table is shown below the chart.

For confidence-gated application, the actual AIRR overlay retains the immutable original call string below threshold. An original string may contain multiple co-optimal alleles, so it is deliberately not represented as a single hard-count bar. The chart uses the local-evidence argmax there to keep the projection one-record/one-allele; the note below the chart states this distinction.

## Sparse browser implementation

- Candidate evidence uses compressed sparse row storage.
- Repeated sparse evidence patterns are collapsed during coordinate updates and multiplied by their accumulated weight.
- Donor/locus/segment pools are fitted separately.
- Dirichlet expected-log weights use a stable digamma approximation; hurdle frequencies use sparse EM and small one-dimensional inclusion quadratures. Responsibilities are normalized in log space in both models.
- The model runs in a dedicated worker. V, D, and J matrices are built and released sequentially to bound peak memory.

For \(N\) modeled records and average sparse width \(K\), iterative work is proportional to \(NK\), rather than \(N\) times the complete reference count.

## Applying and exporting results

Fitting never changes downstream calls. The user explicitly chooses a reassignment policy before applying:

- **Best posterior for every modeled record** replaces each modeled V or J call with its posterior MAP node, regardless of its maximum probability.
- **Best posterior if confidence passes** replaces a modeled V or J call only when its maximum posterior reaches the configurable threshold (80% by default). A call below threshold remains exactly as originally reported.
- A record not modeled for that segment always retains its original call.
- Reset restores immutable original calls.
- Existing downstream collapse, chimera, selection, lineage, SHM, and query results are reset when the upstream overlay changes, preventing partitions derived from one call policy from being reused under another.

Exports include a complete model summary, a long-form sparse per-record posterior sidecar, a refined AIRR table, paired hard-count chart data as CSV, the paired hard-count figure as SVG, and a surviving-allele FASTA for the selected fitted pool. The FASTA uses the current reassignment policy and confidence threshold. Its minimum post-reassignment count defaults to zero and excludes nodes whose hard after-policy count is below the chosen value; sequence-identical labels are emitted as separate FASTA names over their shared reference sequence.

The sidecar reconstructs every sparse candidate responsibility from the saved mixture parameters while streaming the AIRR input, so the complete reads-by-candidate posterior—including zero responsibility for hurdle-excluded candidates—is available without retaining a second large matrix in interactive memory. It records model activity/inclusion probability, the selected reassignment policy, confidence threshold when applicable, and whether each MAP row is selected by that policy. The refined table places policy-selected calls in the ordinary call columns while `swig_original_*` and `swig_repertoire_*` columns retain provenance. Options, compact MAP/entropy and hard-projection vectors, mixture summaries, policy, threshold, and apply/reset state are included in Swig sessions.

## Pipeline position

Allele pooling is an optional first assignment stage. Its model is fitted on the complete assigned input before any collapse, chimera exclusion, or repertoire selection mask. Once the user applies a reassignment policy, collapse and denoising use those selected V/D/J calls when their key or partition depends on germline assignments. The automatic pipeline therefore runs:

`assignment → allele pooling/reassignment → collapse/denoising → chimera exclusion → selection → lineage and diagnostics`.

Later filters do not clear the fitted allele model or revert its calls. Conversely, applying, restoring, or changing an already-applied allele policy resets downstream stages whose partitions may have become stale.

## Interpretation and limitations

- Posterior usage is expressed rearrangement usage, not genomic genotype frequency or copy number.
- Expanded lineages violate independent-read assumptions. Unique weighting is the conservative default.
- A missing true allele cannot be recovered unless represented by a modeled reference node. Use Swig's missing-allele screen and a dedicated workflow such as IgDiscover, TIgGER, partis, or related methods.
- D pooling is experimental.
- The SHM term is a context-averaged approximation. Its sensitivity control exists for explicit sensitivity analysis, not to imply universal calibration.
- A symmetric Dirichlet prior does not force exact zero usage: every known-locus database node retains its prior mass. Local assignment probability is nevertheless exactly zero outside the configured sparse evidence neighbourhood; broadening that neighbourhood is an evidence-model decision, not something the repertoire prior should silently do.
- The hurdle model's “active” state means detectably used in the expressed repertoire. Expressed-read evidence alone cannot distinguish genomic absence from a genomically present but silent allele.
- Hurdle inclusion probabilities use conditional active-set empirical-Bayes quadrature rather than an exhaustive posterior over all allele subsets. Correlated near-identical alleles should be checked for threshold sensitivity.

## Related methods

This narrower assignment-refinement model is informed by the broader subject-specific germline/genotype problem:

- [TIgGER genotype inference](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.00129/full)
- [partis per-sample germline inference](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1007133)
- [PIgLET allele sequence clusters](https://academic.oup.com/nar/article/51/16/e86/7238142)
