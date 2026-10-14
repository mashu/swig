# Post-lineage analysis: SHM, missing-V evidence, and personalized V inference

## SHM metrics

Swig compares `v_sequence_alignment` with `v_germline_alignment` (falling back to complete aligned fields). Only A/C/G/T aligned pairs enter nucleotide denominators. Query coordinates are advanced through query bases and AIRR FWR1/2/3 and CDR1/2 intervals classify positions.

Reported choices are:

- V nucleotide mismatch count and rate;
- translated codon replacement count and rate;
- synonymous changed-codon count;
- CDR1/2 nucleotide mismatch rate;
- framework 1/2/3 nucleotide mismatch rate.

Codons use the AIRR `sequence_frame`/`v_frame`, counted from the read's first base, so the V alignment's first codon may be partial; rows without a frame (IgBLAST's) take the V germline's, whose codons start at germline base 1 (`v_germline_start`); an indel clears the partial aligned codon and resumes at the next frame boundary. A changed codon is counted once as synonymous or replacement, not once per changed nucleotide. Missing/ambiguous codons are excluded.

`duplicate_count` weights abundance summaries, means, medians, and 95th percentiles. Scalar analyzed/skipped counts cover every active row. Plot/session records are memory bounded: at most the configured number per lineage × sample and a deterministic global reservoir of 100,000. Thus a plotted point sample can be bounded even while aggregate record/abundance totals remain complete.

These are **direct descriptive mismatch summaries**, not a context-dependent SHM targeting model, selection test, or phylogenetic mutation reconstruction. Yaari et al.'s S5F work is relevant biological context but is not implemented here: [Models of somatic hypermutation targeting and substitution](https://pmc.ncbi.nlm.nih.gov/articles/PMC3828525/).

## Possible missing V allele screen

This diagnostic requires lineage assignments and treats each donor × lineage as one independent unit. It is a two-pass warning system, not genotype inference.

### Pass 1: discovery representatives

For each eligible low-SHM V alignment, Swig maps observed mismatches to germline coordinates. Per reported V allele it identifies recurrent linked substitution patterns subject to minimum aligned bases, unit/coverage support, allele fraction, binomial-tail probability, maximum SNP count, and diversity in J calls, junction lengths, and CDR3 fingerprints. AID WRCY/RGYW context is annotated but does not waive statistical thresholds.

Default discovery requirements include at least 6 independent supporting lineages, 20 covered units, allele fraction 0.2, SHM rate at most 0.08, at least 180 aligned bases, at most 6 candidate SNPs, and binomial survival probability at most \(10^{-6}\).

### Pass 2: all-member validation

For each proposed linked pattern, Swig scans every retained member of the relevant lineages. It requires joint coverage and linked alternate support, counts units retaining the reference or a conflicting base, requires at least 3 near-germline units, minimum linked fraction 0.9, and maximum other-alternate fraction 0.02. It also requires default diversity of 3 J calls, 3 junction lengths, and 6 distinct CDR3 fingerprints. This second pass prevents one convenient lineage representative from hiding contradictory descendants.

Candidates can be exported as proposed V FASTA after explicit selection. No reference is added to the running assignment database automatically.

## Personalized expressed V-set inference

This is a separate model from the high-specificity missing-V warning. It uses every assigned lineage with a usable V alignment, with exactly one vote per lineage.

### Representative rule

Swig selects the member with the lowest nucleotide mismatch rate in its **current** `v_sequence_alignment` versus `v_germline_alignment`. Ties prefer more aligned V bases and then the earlier AIRR ordinal. There is no random draw, temporal “early” requirement, allele-invariant mask, or circularity correction in this selection step. Lower-SHM observations are used simply because they carry more information about the rearranged germline allele.

### Candidate graph and homologous coordinates

Each subject/locus is processed separately. Same-length V references may compete across gene labels only along explicitly substitution-compatible paths. The affine global alignment scores are match +2, mismatch -3, gap open -5, extension -1. An ungapped path may tie the optimum; a strictly better gapped path excludes that pair. This prevents compensating indels from authorizing raw-coordinate swaps. Every observation scores only directly compatible candidates, even inside a connected component.

SHM calibration uses the published 1,024-context human silent five-mer targeting and substitution tables (HS5F; provenance in `src/shm-model/README.md`). Observed reference sequences are globally aligned to a common longest reference within each subject/locus/V family. Only mapped homologous positions and the same reference-to-alternate event share calibration. The tested physical gene label is excluded. Gapped columns without an anchor mapping do not borrow evidence.

### Novel-haplotype validity test

A proposal is a linked substitution haplotype. Mutation burden excludes precisely its complete proposed site set and the uncertain last 12 V positions; it never masks the union of other proposals. Recurrent patterns are proposed separately as interior cores, observable boundary extensions and complete patterns. The default fraction cutoff is zero, so rare alleles are not excluded merely because their gene is abundant. Four independent lineage representatives, the candidate cap and the evidence test still apply.

The exposure model distinguishes an unmutated component from an SHM-active component. A zero-inflated Gamma-Poisson model is estimated by moments from other genes in the same subject/locus. The active intensity has an exponential prior. The zero proportion is corrected for the configured sequencing-error rate; under each candidate, the posterior unmutated probability and active mean are calculated from noncandidate positions, including the Poisson-error/geometric-SHM count convolution. This is an empirical approximation, not a claim that all repertoire mutation burdens follow one exact distribution. With fewer than 20 external representatives, the unmutated component is not assumed and a diffuse active exposure is used.

HS5F supplies the baseline alternate hazards. Homologous-site, leave-gene-out counts calibrate their location and between-gene dispersion, with equal gene votes. For proposal testing, multiplicative rate deviations use a Student-t prior on log rate (four degrees of freedom), rather than an exponentially decaying prior on rate itself. This admits exceptional gene-specific SHM/selection effects. An SHM event cannot mutate the explicitly unmutated component; sequencing errors can still occur there.

The null also allows a coordinated substitution event inside each connected group of overlapping 5-mer contexts (successive changed positions at most four bases apart). Its hazard rises with SHM exposure and vanishes in the unmutated component. This prevents tandem/nearby SHM from receiving independent-site linkage evidence as though it were an allele. These event-rate nuisance terms occur under both hypotheses. Longer-range dependence is represented only by the shared unmutated/active state, so this remains a restricted SHM model.

The alternative adds one lineage-level mixture component carrying the whole haplotype. Covered site states, active exposure (rounded to 0.0005), and posterior unmutated probability (rounded to 0.001) form sufficient-statistic groups. Null and alternative use identical noncandidate evidence. Scalar rate optimization brackets multiple coarse local maxima; EM and coordinate updates fit the linked mixture.

The test charges half log N for the added mixture frequency plus `log C(L,s) + s log 3 + log P` for discovering s changes among L positions and P expressed parent sequences. Nested proposals are tested against the most specific accepted subset; a terminal extension cannot rescue an unsupported interior core. A second conditional test checks weaker hypotheses on different reference backbones against stronger accepted sequences. These are penalized proposal-validity statistics, not calibrated posterior probabilities or false-discovery rates.

### Common final emission model and set selection

All admitted known and novel sequences compete using the same categorical mutation model: HS5F alternate-specific rates multiplied by the mean of the leave-gene-out Gamma calibration. A lineage's exposure is anchored to its actual current parent alignment, excluding the uncertain terminal window, and cannot change with the size of the proposed candidate catalogue. The three alternate hazards define a total mutation hazard and a normalized substitution distribution; the sequencing-error floor completes the nucleotide probabilities.

The final mixture does **not** multiply a ratio from the haplotype null into an uncorrected parent emission. That would make scores from different parent alleles incomparable. The joint validity test and the final common emission model are distinct documented stages.

Greedy additions optimize a new mixture frequency, followed by an EM refit. Backward removal refits the remaining components. The set objective charges mixture-frequency BIC and the full sequence-search code length for each novel sequence. The active-allele setting is a generous computational guard per gene label, scaled to the number of observed genes in a component; it is not a diploid genotype or copy-number prior.

### Terminal V uncertainty

AIRR `sequence`, `rev_comp`, and `v_sequence_start` reconstruct raw downstream bases from an aligned anchor 12 nt before the V reference end. Query insertions advance raw query coordinates but never become a second observation at a germline position. Both matching and nonmatching downstream bases enter symmetrically. Without raw sequence/anchor coverage, terminal evidence is missing.

Raw boundary evidence stops before the aligned D start (or J start when D is unavailable), so templated downstream-segment bases are not treated as candidate V bases.

For terminal proposal validity, the full raw window is scored under both hypotheses. Each hypothesis profiles the tested gene's distribution over V deletions 0–12 and its initial inserted dinucleotide distribution. Later inserted bases use a leave-gene-out Markov transition model. The alternative adds a shared allele frequency; the same trimming and junction nuisance parameters are fitted under both hypotheses. This permits gene-specific end processing and insertion bias instead of forcing discrepancies with a borrowed prior into an allele call.

Changes confined to the final two V bases are observationally equivalent to a different initial junction word under this conservative null. They are not admitted as novel allele identities, even if a numerical SHM/context approximation would distinguish them. More internal terminal changes can be inferred using downstream templated evidence. Supported interior candidates retain the parent record's terminal bases as placeholders; those inherited bases have not been independently discovered. This explicitly trades some exact terminal recovery for specificity.

The final whole-sequence frequency fit remains a separate, simpler stage: it marginalizes the terminal window with pooled trimming/composition estimates. Its frequencies are model-dependent; the stricter terminal proposal test controls admission of novel sequences. Neither stage is a complete V(D)J recombination model, and a full reassignment remains necessary.

Biological context: gene-dependent trimming and structured insertions are modeled in [IGoR](https://www.nature.com/articles/s41467-018-02832-w); [tandem SHM substitutions](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.807015/full) and [short-range mutational associations](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.618409/full) motivate the correlated-event null. Those studies motivate these model classes; their quantitative parameters have not been transferred as a validated macaque generative model.

### CLI

`swig-cli personalized-germline --airr sample.processed.airr.tsv.gz --v-reference V.fasta --out personalized`

Use the complete processed AIRR table with numeric `clone_id` or `lineage_id`, not a lineage-viewer export limited to the top displayed lineages. The command writes the dashboard JSON (including proposal evidence), evidence TSV and per-pool V FASTA. Web and CLI call the same TypeScript engine. New inference remains in a lazy Web Worker in the browser; V(D)J assignment remains WASM.

The downloadable V FASTA is conservative: inferred active alleles replace only genes actually tested for the selected subject/locus. Untested genes and other loci remain unchanged. Novel records inherit the exact parent's `SWIGMETA` annotation. The download is a candidate reference and must be used in a complete new assignment run; the current analysis is never silently changed.

SHM, possible-missing-V evidence, and personalized inference are separate full-width result sections. They share a post-lineage location because all consume lineage assignments; none modifies another result.

## Prominent incomplete-reference escalation

Swig displays a blocking-style visual warning when either condition holds:

- at least one retained candidate is supported by more than 50 independent lineages; or
- more than five candidates with support at or below 50 independent lineages pass the two evidence stages.

The warning tells the user to rerun with a more complete germline database or use personalized germline discovery with **IgDiscover**. It does not present the candidate FASTA as a validated genotype and does not automatically rerun assignment.

## Literature relationship

The screen is **custom and deliberately conservative**. It shares the problem of subject-specific allele discovery with TIgGER, IgDiscover, and partis germline inference but does not implement their genotype/novel-allele models:

- [TIgGER subject-specific allele identification](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.00129/full)
- [IgDiscover individualized V databases](https://doi.org/10.1038/ncomms13642)
- [partis per-sample germline inference](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1007133)

## Limitations

- AIRR alignment or lineage errors can create recurrent artifacts.
- Clonal dependence is reduced, not eliminated, by one donor-lineage unit.
- Expressed repertoires cannot prove genomic presence/absence or copy number.
- The binomial background is simplified and substitutions are not modeled with S5F context probabilities.
- A warning must be validated by dedicated germline inference and, where important, genomic evidence.
- Personalized inference covers expressed substitution-compatible V alleles. Silent genes, unsupported terminal bases and indel-containing alternatives remain unidentifiable.
- The HS5F calibration, log-rate MAP test, fixed burden estimates, approximate boundary background and two-stage search are computational approximations. Paired dropout controls do not establish the truth of stable background candidates.

### 0.38.8 implementation details

Proposal seeds require linked support, not repeated pristine whole sequences. Their bounded queue prioritizes expected unmutated support using the same exposure model; final-two-base variants do not consume it. Partial haplotypes must compete against admitted strict supersets. Terminal changes must pass both profiled rearrangement and flexible SHM null comparisons. Stable direct log-survival calculations prevent high mutation hazards from producing zero-probability/NaN mixture rows. The CLI exposes `--max-candidates N` (default 64), matching the advanced Web control. See the 0.38.8 benchmark report for remaining biological and computational limitations.
