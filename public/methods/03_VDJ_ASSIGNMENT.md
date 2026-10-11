# V(D)J assignment and AIRR annotation

**Evidence label:** custom. The aligner is Swig's own method; Karlin–Altschul statistics and the AIRR Rearrangement schema supply its significance form and output conventions.

## Core execution

Since 0.39.0 the aligner is the Rust crate in `engine/` (crate and native binary `swig`). It is compiled with WebAssembly SIMD128 to `public/swig.wasm`, which Swig Web runs in a bounded pool of compute workers and `swig-cli` runs in worker threads when it does not run the native binary (since 0.39.1 the release executables run the native `swig` as a child process and keep the module as the fallback; see [the CLI methods](13_CLI_AND_LINEAGE_STUDY.md)). Both use one wrapper, `src/aligner-engine.ts`, which passes run settings by the option names of the native `swig` command line, so the browser, the CLI and the native binary share one option vocabulary, one validation and one mapping to the engine configuration. The module has no imports from its host.

Each worker instantiates its own module and loads the same composed V/D/J/C references once. Loading prepares the database: records with identical sequences are merged into one allele that reports every name, the germlines are annotated (below), and the seed indexes are built. The worker then assigns independent record batches, and Swig commits the resulting AIRR batches in input order. Apart from its identifier and qualities, a row depends only on the read's sequence. Each worker therefore assigns identical sequences once within a batch and keeps the rows of the last 4,096 distinct sequences for repeats in later batches. This cache changes no output.

The run settings are explicit:

- **calling profile**: scores, D exact-run floor and D prior (below);
- **strand**: `both` (decided per read), `forward` (the read as given) or `reverse` (its reverse complement only);
- **minimum identity**: V, D and J hits whose matches cover less than this fraction of their alignment columns are dropped before they shape search windows or enter adjudication. Swig Web and `swig-cli` default to 0.6; the engine alone has no floor. It is an acceptance floor for candidate alignments, not a read-trimming heuristic, and it does not apply to C;
- **covered C-prefix identity**: an optional C gate, off by default (below); and
- the optional [double-D screen](04_DOUBLE_D_SCREEN.md), which writes a sidecar and never changes these rows.

The `sequence` field holds the input read, upper-cased, with `U` read as `T` and any non-IUPAC symbol as `N`. When `rev_comp` is `T`, every coordinate refers to the read's reverse complement. V/J alignment, region and junction fields are derived outputs.

## Search

The complete model and search are specified in [`engine/docs/ALGORITHM.md`](../../engine/docs/ALGORITHM.md). That document also states the guarantee of each search step: exact, lossless (provably drops nothing the model could choose), bounded, or heuristic. In outline:

1. **Orientation.** With strand `both`, the aligner counts, on each strand, the read positions whose 11-mer occurs anywhere in the V or J germline. The larger count wins when it exceeds its Poisson null \(\lambda_0 + 3\sqrt{\lambda_0}\), where \(\lambda_0\) is the number of 11-mers times the germline's 11-mer density. Otherwise, as for heavily mutated reads, V seed votes on both strands decide.
2. **V retrieval (heuristic).** The alleles of a gene share one coordinate frame (a *gene frame*) when their seeds agree on one offset. SNP and 5′/3′-truncated alleles join the frame; an allele with an internal indel, or a divergent one, gets its own. Every read 9-mer votes for gene frames with weight 4. While the leading frame holds fewer than 25 nine-mer seeds, 7-mers also vote, with weight 1. Retrieval keeps 3 genes, and up to 16 while a gene's votes stay within \(z\sqrt{2\cdot4\cdot L}\) of the leader's \(L\) votes (\(z=2\)). Every allele of a kept gene becomes a candidate, aligned in a band of the seed cluster's diagonals ±8.
3. **V and J scoring (lossless).** Score-only dynamic programming runs on one representative per gene frame. Exact upper and lower bounds then skip, settle or align each sibling allele without changing the ranking that full scoring would give. V and J keep the 3 best alignments that pass the E-value gate, plus ties. Each V alignment is extended co-linearly without gaps to the germline's 5′ end or the read's start, and each J alignment to the germline's 3′ end or the read's end. Mismatches in these extensions count.
4. **J (bounded).** Every J allele of the V candidates' locus is aligned, unbanded, from 12 bases before the earliest V end to the end of the read.
5. **D (bounded, then lossless).** D is searched only in junction windows. The first is \([\text{best V end} - 12, \text{best J start} + 12)\); when germline anchors are known, the conserved-anchor window \([\text{V Cys104} + 3, \text{J Trp/Phe118})\) is added. A prefilter keeps D alleles that share an \(L\)-mer with a window, with \(L\) equal to the profile's exact-run floor, so it cannot drop an admissible hit. Each D allele contributes its optimal local alignment per window. Up to 32 alignments with at least 60% identity that hold an exact run at least as long as the profile's floor are admitted.
6. **Alignment kernel (exact).** Affine-gap (Gotoh) local alignment is vectorised across sequences: each 16-bit SIMD lane computes the recurrence for a different read/allele job. Each stage runs over a chunk of 64 reads while every decision stays per read, so a read's result does not depend on its chunk. Swig Web and `swig-cli` run the WebAssembly SIMD128 kernel; the native binary selects AVX-512BW, AVX2 or NEON at run time. Every backend reproduces the scalar Gotoh kernel's scores and tracebacks exactly.
7. **Significance gates.** A V or J hit is kept only when its Karlin–Altschul E-value for the searched window and the whole segment database is at most 1. A read is called only when its path contains an *anchor*: a V or J hit with E ≤ 10⁻³. Random sequence therefore receives no call. A J downstream of a significant V needs only E ≤ 1. D has no E-value gate; the exact-run floor admits it and the D prior prices it.
8. **Adjudication (exact over the retained hits).** Every retained V, J and admitted D alignment becomes a placement in one directed acyclic graph. A path is V → D (zero or one) → J; partial paths (V only, J only) are allowed. Between consecutive placements there are two cut points: the first segment keeps the read bases before its exit, the next keeps those from its entry on, and the bases in between are N/P. Overlapping placements can therefore split the overlap at any base or leave part of it to N/P. The clipped score of a placement separates into one term per cut. This lets an exact semi-Markov Viterbi recursion choose the topology (VJ or VDJ), the V, D and J alleles, and both cut points of every junction together, as an exact optimum over every grammatical path and canonical cut pair. No V or J endpoint is fixed before D is considered. Among equal values a path with fewer D wins, so a D must strictly beat the best path without one. Remaining ties go to the path with more germline support (summed unclipped scores). The main call allows at most one D; tandem D is examined only by the separate double-D screen.
9. **Per-read divergence.** By default the call maximises each path's log-odds at the read's own divergence \(\hat p\). \(\hat p\) is estimated before adjudication from V and J columns that no path contests, as \((k+\tfrac12)/(n+1)\), clamped to \([10^{-4}, 0.5]\). A matching column is then worth \(\ln 4(1-\hat p)\) nats and a mismatch \(\ln(4\hat p/3)\). An `N` column is worth 0, and gaps stay at the prior \(\lambda\). The D prior cost is converted to \(P^*_D/m\) matched D columns (\(m\) the D match score), so that its specificity against chance exact runs does not depend on divergence. For +2/−3 scoring, at \(\hat p = 11.2\%\) the model coincides with the fixed scores. The reported `*_score` values are the fixed-scale integer scores of the called, clipped alignments.
10. **Posteriors.** The same objective defines a distribution over paths and cut points. Forward–backward over the same graph gives the marginal probability of every placement and of the D count, and the call is a mode of that distribution. `v_posterior`, `d_posterior` and `j_posterior` give the probability that the reported names include the read's allele, up to germline equivalence. The sum covers the placements tied with the call and those that explain most of the called bases with germline bases the primary allele contains. `d_presence_posterior` is \(1 - P(\text{no D})\). If a read's evidence range would leave floating-point range, its posterior fields are left empty rather than reported wrongly.
11. **Co-optimal names and alternatives.** A `*_call` lists, comma-separated, every allele whose placement clipped to the called interval is worth exactly the call: an equal score, and under the per-read model also equal matches and mismatches. The best-supported allele comes first. `*_alternatives` lists the remaining alleles whose clipped placement scores at least \(S - \max(2, \lfloor |S|/100 \rfloor)\), where \(S\) is the call's score. Entries are written as `name|score|identity|start|end|germline_start|germline_end`, `;`-separated, best first. Equal names mean equal value under this caller, not equal biological probability; the `*_posterior` columns carry the probability.

## Calling profiles

A profile changes only the scores, the D exact-run floor and the D prior. Scoring tuples are written match/mismatch/gap open/gap extend; a gap of length \(L\) costs open \(+ (L-1)\cdot\) extend.

| Profile (config name) | V | D | J | D exact-run floor | D prior cost \(P^*_D\) |
|---|---|---|---|---:|---|
| **R-optimized** (`r_optimized`, default) | +2/−3/−9/−1 | +2/−3/−13/−1 | +2/−3/−17/−2 | 4 nt | derived from the D genes of the read's locus, even odds: 10.55 for human IGH |
| **Sensitive-D** (`sensitive_d`) | +2/−3/−9/−1 | +2/−3/−13/−1 | +2/−3/−17/−2 | 4 nt | the same with a D twice as likely a priori: 9.45 for human IGH |
| **Truth-optimized** (`truth_optimized`) | +2/−3/−5/−1 | +2/−3/−13/−1 | +2/−3/−17/−2 | 6 nt | none: a D must only beat the no-D path |
| **IgBLAST-agreement** (`igblast_compatible`) | +2/−3/−5/−1 | +2/−4/−11/−1 | +2/−4/−13/−1 | 5 nt | none: a D must only beat the no-D path |

A clipped D must also keep at least the floor's number of read bases.

**The derived D prior.** A D hypothesis is one of \(N\) D genes placed at one of \(L = 25\) junction positions. A MAP decision includes the best D when \(\lambda S_D > \ln(NL) - \ln \pi_D\), where \(\pi_D\) is the prior odds of a D against none. The cost in score units is therefore

\[
P^*_D = \frac{\ln(N L) - \ln \pi_D}{\lambda},
\]

with \(\lambda = 0.634\) for the +2/−3 D substitution scores. \(N\) counts the D genes of the read's locus only, so in a multi-locus database each locus is priced by its own D set. R-optimized uses even odds (\(\pi_D = 1\)): 10.55 score units for the 32 human IGH D genes, and 6.17 for a locus with 2 D genes such as human TRB. Sensitive-D uses \(\pi_D = 2\), which lowers the cost by \(\ln 2\) nats (1.09 score units): 9.45 for human IGH. On the fixed scale a D must therefore score more than \(\lfloor P^*_D \rfloor\). Under the per-read model it must be worth more than \(P^*_D/2\) matched D columns. Without mismatches, that is a D tract of at least 6 nt under R-optimized and 5 nt under Sensitive-D for human IGH. For a 2-gene locus the 4-nt floor is the binding limit. The cost is a prior, not a fitted threshold; the profile tuples were selected on simulated reads in earlier Swig releases.

**Removed in 0.39.0.** The SwiftIG C++20 core, its four assignment strategies (standard, RIAT-MP, AER and the experimental AER-R) and the IgBLAST-balanced profile were removed. There is one search. Portable configs that still name a strategy (`annotation.assignerStrategy`) load with it silently dropped, and `swig-cli --vdj` rejects `--assigner` as an unsupported option. A saved browser session or lineage-study manifest whose profile the aligner does not have restores as R-optimized.

The IgBLAST-agreement profile is not an IgBLAST implementation. IgBLAST uses its own BLAST-based germline search and annotation logic; see Ye et al., [IgBLAST](https://doi.org/10.1093/nar/gkt382). The profile name records the objective its tuples were chosen for (agreement with IgBLAST calls on a supplied simulation), not equivalence of calls.

## Segment support values

`v_support`, `d_support` and `j_support` are Karlin–Altschul expectation values: the expected number of chance local alignments with at least the reported segment score. Smaller is stronger. They are neither posterior probabilities nor repertoire-level evidence that an allele is present.

For the reported (clipped) score \(S\) of a called segment, the aligner evaluates

\[
E = K\,m'\,n'\,e^{-\lambda S}, \qquad
\ell = \frac{\ln(K m n)}{H}, \quad m' = \max(1, m-\ell), \quad n' = \max(1, n-N\ell).
\]

- \(m\) is the length of the whole oriented read, for V, D and J alike.
- \(n\) is the total nucleotide length of that segment's database after identical records are merged; it includes every locus in the database.
- \(N\) is the number of distinct allele sequences in that database.
- \(\lambda\), \(K\) and the relative entropy \(H\) are the exact ungapped Karlin–Altschul parameters of the segment's match/mismatch scores under uniform base composition. They are computed when the run starts for each scoring tuple in use, including D/J mismatch overrides: \(\lambda\) by bisection and \(K\) from the Karlin–Altschul series. The computation reproduces NCBI's ungapped blastn values (reward 1/penalty −3: \(\lambda = 1.374\), \(K = 0.711\)). For +2/−3, \(\lambda = 0.634\) and \(K = 0.408\).
- \(\ell\) is the BLAST finite-length correction.

The search gates above use the same \(\lambda\) and \(K\), with the length of the searched window and without the length correction. The reported value is computed afterwards for the reported score and the whole read, so it need not equal the gating quantity exactly.

Gapped alignments score slightly higher than ungapped ones, so these values are mildly anti-conservative. Near-identical alleles enlarge \(n\), which is conservative. D has no significance gate, so a short D that the joint model calls can have `d_support` well above 1. `c_support` is a different, anchored quantity (below), and the columns are not comparable.

The form follows Karlin and Altschul's [local-alignment statistics](https://doi.org/10.1073/pnas.87.6.2264) and BLAST's treatment of search spaces ([Altschul et al. 1997](https://doi.org/10.1093/nar/25.17.3389)). The values are **not IgBLAST E-values** and are not expected to equal them, because IgBLAST uses its own scoring, gapped statistics, databases, query ranges and hit construction.

## Multi-locus databases

Swig composes one V, D, J and C database across the selected loci, for example IGH, IGK and IGL of a BCR study, or every TR locus. Each allele's locus is read from its label (`IGHV1-2*02`, `TRBJ2-7*01`) or, when the label names none, from its FASTA header. Shared TRA/TRD V genes (`TRAV14/DV4`) belong to both loci. An allele whose locus cannot be read pairs with every locus.

A rearrangement joins segments of one locus. J candidates must pair with the V candidates' locus, D candidates with the V's and J's, and C with the called J's. When a read's placements span several loci, each locus is adjudicated over its own placements. The call then comes from the locus with the largest total likelihood (or, without posteriors, the best path), and the posteriors are the likelihood-weighted mixture over loci. Germline annotation runs per locus. The D prior and the C gate count only the genes of the read's locus, so adding loci to a database does not change the D or C threshold of an IGH read. The V/J significance gates and the reported V/D/J E-values, however, use each segment database's full length. The `locus` column reports the locus the path's alleles share.

## Constant-region calling and isotype labels

With a C database, C is assigned after adjudication, from finished calls only. C never enters the placement graph, so it cannot change V, D or J calls, coordinates, strand or locus. A run without C does no C work.

J and C are joined by an RNA splice at a fixed phase: there are no N/P bases and no hypermutation. Let \(x\) be the called J's exit; C is considered only when J ends the path and read bases follow it. Each hypothesis places base 0 of a C allele from the J's locus at read position \(x+\delta\), \(|\delta| \le 3\). The tolerance absorbs reference conventions, which differ by a base (IMGT CH1 records start after the J's last base, codon-aligned records repeat it, 3′-partial J records end a base early), and indels near the J 3′ end, which move \(x\).

Each alignment is anchored at C base 0 on one of the 7 start diagonals and has a free end. The ordinary local kernels compute it exactly through an anchor column, and it is then clipped to the read bases after \(x\), so C explains only bases the J does not. Alleles that are identical over the reachable length form one class, aligned once. The free end lets the alignment stop before a primer, adapter or other downstream sequence, so downstream disagreement does not invalidate an otherwise contiguous C tract.

C is called when its anchored E-value is at most \(10^{-3}\):

\[
E_C = N_C \cdot 7 \cdot e^{-\lambda S},
\]

where \(S\) is the clipped C score under +2/−3/−9/−1 (\(\lambda = 0.634\)) and \(N_C\) is the number of C genes of the J's locus in the C database. This is a union bound over genes and start diagonals, reported as `c_support`. For the 10 human IGHC genes it requires \(S \ge 18\), nine exact bases. Because \(N_C\) counts genes, a smaller C database supports a shorter tract; no fixed biological minimum length is asserted.

`c_call` lists the co-optimal names. Ranking is by clipped score, then the longest exact run, matches and allele order. A short remnant often matches a 5′ end that several genes share (IGHG1, IGHG2 and IGHG4 open with `GCCTCCACCAAG`); `c_call` then lists all of them, and the first name is first only by that tie-break. `c_posterior` is the probability that the reported names include the read's C allele, given that the read reaches C. It weighs each kept C placement by its log-odds at C's own divergence, \(p_C = (\text{mismatches} + \tfrac12)/(\text{informative columns} + 1)\) clamped to \([10^{-4}, \tfrac12]\), because C is not hypermutated. No C-presence probability is reported.

**Covered C-prefix identity.** `c_prefix_identity` is the identity of the called C germline from its first base. Any germline prefix that the alignment does not cover is compared base for base with the read bases before the alignment, followed by the alignment's own columns. The optional gate (0, the default, turns it off) discards a C call whose value is below the chosen fraction, or whose prefix the read does not cover. Failing the gate clears only the C fields.

**References and isotype labels.** The bundled reference builder assembles one secreted/coding path per IMGT C allele and excludes membrane-only and UTR exons. The caller can therefore distinguish only the constant genes and alleles in the selected C FASTA, and only when the sequenced tract covers discriminating bases.

The browser derives an `isotype` label from `c_call` when `c_identity` is at least 0.90 and the C alignment covers either at least 30 read bases, or at least 12 bases with `c_support` ≤ 10⁻⁴. A short called C must therefore be ten times stronger than the calling gate before it receives an isotype label. Ambiguity across IgG or IgA subclasses is collapsed to `IgG` or `IgA` rather than choosing the first subclass arbitrarily; other mixed labels are joined with `/`.

Secreted-versus-membrane isoform assignment requires a reference FASTA that explicitly contains the alternative exon paths, and reads that cover them; it is not inferred from an unobserved tail.

## Reference metadata versus read assignment

FWR/CDR coordinate transfer for a custom germline is a reference-preparation operation, not part of read-to-allele mapping. Strict, permissive and best-guess metadata matching select an annotated template from which V-region boundaries or J frame/CDR3-anchor information are projected into `SWIGMETA`. They do not alter the allele's nucleotide sequence or add or remove reference alleles.

The aligner treats `SWIGMETA` (and explicit `anchor=`/`fr1=`… header keys) as given values. Its own germline annotation fills only what they leave open. In order of trust, that annotation draws on:

1. the IMGT numbering of an IMGT-gapped V set;
2. a motif search placing Cys23, Trp41 and Cys104, and the J Trp/Phe118;
3. the same search with a database-consensus prior; and
4. transfer from the most similar confidently annotated allele.

An FR3 start that no relative supplies comes from a database model of FR3 length, or from the chain default. Starts that may be wrong are flagged. All of this runs once per database, before any read. `v_annotation_source` and `j_annotation_source` name each called allele's source: the `SWIGMETA` provenance label (for example `IMGT-gapped` or `IMGT-boundary-transfer`), or the aligner's own label (`header`, `imgt-gaps`, `motif`, `consensus`, `transfer:<allele>`, `low-confidence`, `missing`).

Germline metadata enters assignment in one place only: the anchors define the conserved-anchor D window, which is added to the V/J-derived window. A different anchor can therefore add or remove D candidates. Metadata never changes seeds, alignment scores, V or J candidates, or the adjudication objective.

After assignment, the chosen alleles' region coordinates are projected through their already-selected read alignments. This yields the FWR/CDR sequences, junction/CDR3 boundaries, reading frame, translation and productivity. A permissive metadata choice can thus change these derived fields, and analyses that filter or group on them, without changing the V/D/J alignment. `swig-cli prepare-reference` persists this schema once; its diagnostics report the transfer template, identity, structural checks and unresolved reason for every allele.

## Output and record detail

The aligner writes Swig's AIRR Rearrangement schema of 111 columns, identical in every row, and Swig adds `swig_*` study/provenance fields:

- **AIRR fields.** The AIRR Rearrangement fields come first.
- **Swig additions.** The Swig-specific columns are `c_prefix_identity`, `region_definition`, `v_annotation_source`, `j_annotation_source` and `v/d/j/c_alternatives`.
- **Posterior columns (new in 0.39.0).** The schema ends with `v_posterior`, `d_posterior`, `j_posterior`, `c_posterior` and `d_presence_posterior`. Probabilities print with six decimals, or ten within 10⁻⁴ of 0 or 1.
- **Empty C fields.** Without a C database, every `c_` column is empty.

Coordinates are 1-based closed intervals on the oriented read. CIGAR strings use `S` for read bases outside the segment, `N` for germline bases before the alignment, and `=`/`X`/`I`/`D` for columns. `*_identity` is matches over alignment columns, as a fraction.

Region fields are projected through the called V and J alignments from `SWIGMETA` or, where it is absent, from the aligner's own annotation. The fields are defined as follows:

- **Junction.** The junction runs from V Cys104 through J Trp/Phe118, inclusive, and is reported only when it lies inside the read.
- **Frameshift.** A segment has a frameshift when its called alignment's net indel length is not a multiple of 3.
- **`productive`.** True only for a rearrangement that is in frame, has no stop codon in the V reading frame from the V start to the J end, and has no V or J frameshift.
- **`complete_vdj`.** True when V starts within 3 bases of its germline start and J ends within 3 bases of its germline end.
- **`region_definition`.** `IMGT` when the called V has region coordinates.

The interactive record view renders these committed fields: query and germline coordinates, FWR/CDR intervals, the stitched V(D)J alignment, individual segment alignments, segment E-values, junction decomposition, co-optimal names and near-tied alternatives. It does not rerun assignment when opened.

The schema follows the AIRR Rearrangement representation so downstream tools can consume the table: [AIRR Community standardized representations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6173121/). Swig adds namespaced columns where the AIRR core schema has no equivalent. AIRR compatibility does not make the Swig aligner another AIRR-producing caller such as IgBLAST or IMGT/V-QUEST.

## Accuracy and speed

Measured comparisons of the 0.39.0 aligner with SwiftIG 0.38.11 are in [`BENCHMARK_SWIG_ALIGNER_0.39.0.md`](../../BENCHMARK_SWIG_ALIGNER_0.39.0.md), covering call accuracy on simulated reads and throughput. The aligner's evaluation protocol is in [`ALGORITHM.md` §7](../../engine/docs/ALGORITHM.md): IgSim simulation with base-level germline provenance, germline-equivalent allele accuracy, paired comparisons with exact McNemar tests, D accuracy reported with D-free specificity, and posterior calibration.

## Performance and cancellation

Each worker prepares the references and builds its indexes once. Records are processed in bounded batches, and results are appended to browser storage or a user-selected writable stream. Cancellation aborts input reading, terminates the worker pool, aborts unfinalized output and clears the temporary result store; no partial run is opened as a completed analysis.

The batched pipeline equals the per-read scalar pipeline, and its output is byte-identical across SIMD backends and across x86-64 and AArch64. The WebAssembly module gives the same rows, with one exception: a posterior whose value lies within a last-bit difference of `exp` of a printed rounding tie can differ in that digit (2 of 260,000 rows in the engine's tests).

## Limitations

- The profile tuples and gates were selected on simulated reads, not by external biological validation. No simulator matches every repertoire's SHM, indel, trimming or junction distributions.
- Short D segments and highly mutated reads can remain intrinsically ambiguous. A D remnant shorter than the profile's exact-run floor cannot be found, and very short D tracts are hard to distinguish from chance junction matches.
- R-optimized can call a short chance D on reads whose D was fully trimmed. This is most frequent in loci with few, short D genes such as TRB: with 2 genes the derived cost (6.17 score units) is paid by any exact 4-nt run, which scores 8. Sensitive-D lowers the cost further. Read `d_presence_posterior` and `d_posterior` alongside such calls.
- Posteriors are model-based and conditional on the retained candidates: an allele that retrieval or the keep policy dropped receives no mass. They assume one divergence per read, uniform across V, D and J, and inherit the D prior. They are not externally calibrated probabilities, and D presence is the least well calibrated of them.
- E-values use ungapped \(\lambda\) and \(K\) and the total database length. They gate hits and describe chance expectation; they are not calibrated p-values.
- C calls depend on sufficient post-J evidence and on whether the selected C FASTA contains the relevant gene or isoform path. Net indel drift of more than 3 bases inside C truncates the C alignment. Very short reads may be insufficient, especially with a large C database.
- Important study-specific calls should be benchmarked against an independent workflow and suitable truth/control material.

## Assignment progress

Multi-file progress uses read counts only when all inputs have known counts. Otherwise every file uses the same expanded-file-size workload proxy; an unknown count is never assigned a one-read weight. Gzip footer sizes are read without a counting pass, with a rough compressed-size expansion fallback for large/wrapped sizes.

Within a file, parsed character progress is discounted by the fraction of parsed records actually acknowledged by the AIRR store. Compressed read-ahead and dispatched batches are not completed assignments. End-of-input switches to the exact committed/parsed count. Subsampling separates the scan from selected-read commits. Index finalization reserves the last percentage; 100% is shown only when results are ready. These are approximate work fractions, not elapsed-time predictions.
