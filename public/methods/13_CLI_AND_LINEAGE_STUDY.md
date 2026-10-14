# Standalone CLI, portable pipeline configuration, and lazy lineage studies

**Evidence label:** direct reuse of Swig analysis kernels plus custom execution/interchange engineering.

## Runtime and distribution

`swig-cli` is compiled into one operating-system executable. The release binary embeds the Bun runtime, the native aligner built for its platform (the `swig` binary of the Rust crate in `engine/`), the same `swig.wasm` aligner module used by Swig Web (that crate compiled with WebAssembly SIMD128; it needs no WASI or other host imports), the built-in reference pack, the CLI worker, and the shared TypeScript post-analysis kernels.

Since 0.39.1 the CLI runs the native aligner wherever the executable carries one. A file embedded in a Bun executable cannot be executed in place, so the CLI copies it once per version and content hash to the user cache directory (`SWIG_CLI_CACHE_DIR`, else the platform cache directory, else a directory in the temporary directory that is private to the user and checked to be so before anything in it runs), checks that `swig --version` names the CLI's own version, and starts it as a child process with the prepared germline FASTA, the same options the WebAssembly engine receives, and `--threads` equal to the worker count. Records stream to its standard input and AIRR rows return on its standard output in input order; the engine writes each 16,384-record batch as soon as it is assigned, so the CLI keeps more than that many records in flight past the batch it waits for, and it waits only when the next batch is ready to send (the short last batch reaches the engine with the end of input). A record without an identifier is named `sequence_N` by its position in the input before it is batched, so names do not depend on the route or the worker count. On identical prepared references the native and WebAssembly routes write the same AIRR (byte-identical in the tests on 17,500 reads and on FASTA, FASTQ, and AIRR pipeline datasets; a posterior within a last-bit difference of `exp` of a printed rounding tie can differ in that digit, as described in [V(D)J assignment](03_VDJ_ASSIGNMENT.md)); per core the native route is about 1.7 times faster, mostly from 256-bit AVX2 vectors and native code generation, and it avoids the per-process JIT compilation. `--aligner auto` (the default; `SWIG_CLI_ALIGNER` sets it, and `annotation.aligner` in JSON) falls back to WebAssembly, saying why, when no usable binary exists or the child fails to start; `--aligner native` fails instead; `--aligner wasm` always uses WebAssembly. `SWIG_CLI_ALIGNER_PATH` names another native binary of the same version, which takes precedence over the embedded one. The Linux binaries link glibc statically and the Windows binary its C runtime, so the native route needs nothing newer than the executable itself (glibc 2.17 on Linux). Double-D screening runs only on WebAssembly, because the native route does not carry the sidecar. The pipeline starts one child per dataset that needs annotation, since the engine detects FASTA, FASTQ, or AIRR once per stream. The person running an analysis does not install Node, npm, Bun, Python, a container engine, or a separate WASM runtime.

Release automation builds Linux x86-64, Linux ARM64, macOS Intel, macOS Apple silicon, and Windows x86-64 binaries. It executes the end-to-end smoke analysis against the actual Linux and macOS artifacts on matching runners; the Windows artifact is cross-compiled on Linux and runs the native-route tests on a Windows runner before the release is published. Baseline x86-64 targets avoid requiring AVX2. Bun is a build-time packaging choice, not a distinct scientific implementation: annotation and post-analysis call the repository's shared code. Releases are built with Bun 1.3.14; Bun 1.4.2 closes `node:readline` while a paused consumer still holds unread lines, which truncates streamed input.

The repository also emits an unpackaged JavaScript CLI during developer tests. That artifact exists to make local regression tests cheap; it is not the documented end-user installation path.

## Configuration and defaults

The configuration is JSON with `schema: 1`. Every analysis block is optional when a default is acceptable; `output.directory` is deliberately required unless `--out` supplies it. `swig-cli init` writes a complete example. Swig Web exposes **Export CLI config** both before execution in **Analyze → Execution mode → Pipeline** and after execution in Results. It writes:

- exact V/D/J/C FASTA used by the browser run;
- species and receptor/locus scope;
- FASTQ expected-error filtering, 3′ trimming, exact per-dataset subsample size, and base random seed;
- calling profile, AIRR preserve/reannotate policy, strand, V/D/J identity floor, optional covered C-prefix identity floor, worker count, and Double-D settings;
- enabled state and parameters for allele refinement/reassignment, collapse or denoising, CHMMAIRRa, repertoire selection, lineage assignment, SHM, and possible-missing-V evidence;
- input paths and formats, concatenated-gzip member ranges when applicable, plus dataset, biological-sample, donor/subject, cohort, timepoint, and compartment labels; and
- output naming and ordinary AIRR compression (`none` or `gzip`).

Input paths may be edited or replaced for batch use. `subjectId` is the donor grouping key: different files/samples with the same nonempty value are in the same donor. Collapse, allele-pooling, and lineage scopes remain separate explicit settings, exactly as in the browser. `annotation.workers: 0` (the hand-written/default value) selects up to eight host threads; a positive integer is an exact worker count with no CLI-side maximum. `--workers N` is applied after reading JSON and therefore always overrides `annotation.workers`. The resolved count controls assignment and the independently parallelizable collapse/denoising and CHMMAIRRa work as well. Assignment-only `--vdj` uses every core (`availableParallelism()`) when neither thread option is present or with `--workers 0`; an explicit `--workers 32` starts 32 native threads or WebAssembly workers. If no config is supplied, one input may be annotated with `--sample`, `--donor`, `--dataset`, and `--out` flags. Every pipeline run must name an output directory explicitly through `output.directory` or `--out`; normalization no longer turns an omitted destination into a silent working-directory write.

Local replacement FASTA in `references.files` is prepared before the aligner indexes it, by default. The CLI calls the same germline preprocessor and species-first taxonomic template tiers as Swig Web, using the embedded fixed reference pack to transfer and validate V FWR/CDR boundaries and J frame/CDR3-anchor metadata. D/C files receive the same parsing, locus/segment checks, normalization, and metadata validation. No reference or read is sent over the network. Set `references.prepareMetadata` to `false` only when the file must be passed through without metadata transfer:

```json
{
  "references": {
    "species": "Homo sapiens",
    "scope": "IGH",
    "prepareMetadata": true,
    "files": {"V": "custom-V.fasta", "J": "custom-J.fasta"}
  }
}
```

Swig Web exports already-prepared exact FASTA through `references.inline`; those strings are not reprocessed. If both an inline and file value exist for one segment, the file remains the replacement and follows the `prepareMetadata` setting. Preparation counts and warnings are written to stderr and the run summary, while the resolved config records the effective flag.

### Reusable prepared references

`swig-cli prepare-reference` runs metadata preparation without processing reads. It accepts the same source-FASTA, organism, receptor-scope, and optional IgBLAST metadata inputs as the assignment-only route. It writes normalized per-segment FASTAs whose headers contain validated `SWIGMETA`, a portable `.swig-reference.json` manifest, and an `.annotation-diagnostics.tsv` table. The manifest stores input/output SHA-256 values, fixed-pack release and hash, effective matching settings, annotation counts, per-record provenance, optional IgBLAST FWR4 offsets, and relative paths to the prepared FASTAs. `swig-cli --vdj --prepared-reference MANIFEST` resolves those paths relative to the manifest and verifies every output hash before starting workers.

```bash
swig-cli prepare-reference \
  -germline_db_V V.fasta -germline_db_D D.fasta -germline_db_J J.fasta \
  -organism human -ig_seqtype Ig --out-prefix refs/custom
swig-cli --vdj -query reads.fasta.gz \
  --prepared-reference refs/custom.swig-reference.json -out calls.airr.tsv
```

The matching presets are:

| Mode | V same-gene / nearest identity | J same-gene / nearest identity | Ranked fallback candidates |
|---|---:|---:|---:|
| `strict` | 0.80 / 0.72 | 0.75 / 0.68 | 12 |
| `permissive` | 0.65 / 0.55 | 0.60 / 0.50 | 32 |
| `best-guess` | 0 / 0 | 0 / 0 | 64 |

Each threshold and fallback count can be overridden explicitly. Best-guess removes identity as a rejection rule; it does not fabricate a V region map or J frame when projected coordinates are contradictory or the required frame-consistent conserved anchor is absent. Diagnostics distinguish identity rejection, unmappable/nonmonotonic boundaries, empty projected regions, absent V cysteine anchors, incomplete J-anchor projection, target J-motif mismatch, and invalid projected metadata. `--require-complete` writes the manifest and diagnostic artifacts but exits nonzero if any V/J record remains unresolved.

These options belong strictly to germline **metadata transfer**. They never change the supplied reference sequences or allele collection and are not consulted by seed retrieval, affine alignment, scoring, or the adjudication objective. The aligner maps a read against the submitted nucleotide sequences first; the selected allele's metadata is then projected through that alignment to populate FWR/CDR, junction, translation, frame, and productivity fields. Consequently, metadata choices can affect downstream analyses that consume those derived fields. The one place where germline anchors enter assignment is the conserved-anchor D search window, so a different anchor can add or remove D candidates; metadata never changes seeds, alignment scores, or V/J candidates (see [V(D)J assignment](03_VDJ_ASSIGNMENT.md)).

The fixed order is annotation → optional allele fit and policy reassignment → collapse/denoising → optional chimera filter → optional explicit selection → lineage assignment → SHM → possible-missing-V evidence. The CLI deliberately does not infer lineage phylogenies or UCAs because those analyses require choosing and inspecting a lineage MSA.

FASTA, wrapped FASTA, FASTQ, wrapped FASTQ, AIRR TSV/CSV, and gzip-compressed versions are streamed in bounded batches. Annotation uses a bounded worker pool. The annotated AIRR header is written after the first completed input batch and each subsequent batch body is appended in input order; preserved AIRR inputs therefore may not introduce a new column after that header has been committed. The retained table is serialized in bounded row chunks after global policy/filter decisions instead of constructing one repertoire-sized string. `output.airrCompression: "gzip"` or the later `--airr-compression gzip` override wraps both ordinary AIRR writers in bounded gzip streams and appends `.gz`; `none` is the default. The separately linked lineage-study AIRR stays uncompressed because its manifest indexes exact byte ranges. FASTQ filtering and seeded reservoir sampling use the same shared functions as the browser. A `gzipRange` identifies one complete gzip member by end-exclusive compressed byte offsets, allowing two config inputs to address separate members of one source file. Small data pasted into Analyze are stored in that exported input's `inline` field, making the config self-contained instead of referring to a nonexistent `pasted-sequences.txt`; file inputs remain path based. AIRR outputs include collision-safe dataset/sample/donor fields and retain the source sequence ID; call reassignment occurs before call-dependent collapse and lineage grouping. The resolved config and a machine-readable run summary are always emitted. Repertoire-wide enabled methods still retain the compact/full row state those methods actually consume; streaming the file writer does not turn a global allele, deduplication, chimera, or lineage calculation into an online algorithm.

`annotation.airrMode` is deliberately explicit. Browser exports use `reannotate`, matching browser upload: the AIRR `sequence` is sent through the aligner again and old assignments are not trusted. The default for a hand-written/legacy CLI config is `preserve`, which treats existing AIRR calls as the annotation result. Double-D screening with AIRR requires reannotation because the sparse evidence calculation needs the aligner pass.

## Assignment-only `--vdj` route

`swig-cli --vdj` is a separate execution route for high-volume assignment. It does not call `annotateAirrBatch`, construct study manifests, parse results into post-analysis records, retain repertoire rows, or run allele pooling/collapse/chimera/selection/lineage/SHM/missing-allele code. The query parser supplies bounded FASTA batches to the native aligner's input stream or to the WebAssembly worker pool, completed batches are consumed in submission order, and each AIRR body is written immediately after its single header. `-out`/`--out` is mandatory; `-out -` is an explicit streamed stdout destination. An output path ending in `.gz` is gzip-compressed while plain `.tsv` and stdout remain uncompressed. AIRR writers replace tab/CR/LF characters inside user-controlled identifier or quality values so every record retains the header's tab-delimited field count. `-query -` reads FASTA from stdin. `-outfmt` currently accepts only `19` because other IgBLAST renderers are not Swig output formats.

The route accepts the common IgBLAST names `-query`, `-germline_db_V`, `-germline_db_D`, `-germline_db_J`, `-c_region_db`, `-out`, `-outfmt`, `-organism`, `-ig_seqtype`, `-domain_system`, `-strand`, `-num_threads`, `-custom_internal_data`, `-auxiliary_data`, and `-d_frame_data`. Source germline arguments are nucleotide FASTA (plain or gzip), not BLAST database prefixes. `--workers` has higher precedence than `-num_threads`. Without either, or with `--workers 0`, the route uses every core: native threads, or one WebAssembly worker per core.

Swig-specific `--batch-records`, `--minimum-identity` (the V/D/J identity floor, 0.6 by default), `--minimum-c-prefix-identity`, and `--calling-profile` remain available. `--calling-profile` accepts `r_optimized` (the default), `sensitive_d`, `truth_optimized`, and `igblast_compatible`. Pipeline JSON validates the same names in `annotation.callingProfile`, and an unknown name fails the run. The CLI help states the R-optimized and Sensitive-D D priors: 10.55 and 9.45 score units for human IGH. Since 0.39.0 the aligner has one search and no assignment strategy. `--assigner` therefore fails as an unsupported option, and JSON configs drop a legacy `annotation.assignerStrategy` silently.

Four IgBLAST options override the selected profile's D/J settings:

- `-min_D_match` (at least 4) sets the D exact-run floor;
- `-D_penalty` and `-J_penalty` (integers from −1000 to −1) set the D and J mismatch scores; and
- `-min_J_length` sets the fewest aligned J columns.

The profile's other settings stay in force, including the rule from which its D prior is derived. `-num_alignments_D` and `-num_alignments_J` are validated and accepted for IgBLAST compatibility, but they do not limit the search: every admissible D and the best J alignments enter adjudication, and near-tied alleles are reported in `*_alternatives`. `-strand` accepts `both`, `plus`, or `minus`, and `-domain_system` accepts only `imgt`. Unsupported BLAST scoring, search, or report options fail rather than being silently ignored.

A supplied C database uses the same splice-anchored rule as Swig Web. C is aligned only after the called J's end, with the first C base within 7 bases of it (SwiftIG's tolerance; 3 bases in 0.39.0–0.39.3). It is called at an anchored E-value of at most 10⁻³ over the C genes of the J's locus and 15 start offsets, on at least 12 aligned columns at 90% identity or more (SwiftIG's floors), and it never changes V, D, or J. The optional 0..1 `--minimum-c-prefix-identity` / `annotation.minimumConstantPrefixIdentity` control examines the complete covered prefix and, on failure, clears only C. The general `--minimum-identity` option remains the V/D/J floor; under the plausibility guards it also applies to the part of a V or J alignment that a call keeps.

The pipeline route and `--vdj`, native or WebAssembly, run the aligner with SwiftIG's plausibility guards, which have no `swig-cli` switch: V hits of at least 24 aligned columns and J hits of at least 10; a V–J pair only when it faces its junction, holds a long split or aligns both conserved codons in order, with at most 225 bases between V and J; a D only beside a V or a J, and at the read's open end only when it pays for the J or V the read would show beyond it (the open-flank cost); and a junction (with CDR3, `vj_in_frame`, `productive` and `complete_vdj`) only when the called V and J align the first bases of Cys104 and Trp/Phe118, at most 231 bases apart. A chimeric or truncated read then keeps its better segment instead of becoming a rearrangement, and a read that stops at its D or starts after its V keeps its D, while random bases beyond a V or before a J rarely get one. [V(D)J assignment](03_VDJ_ASSIGNMENT.md) gives the rules and their measured effect.

Both the pipeline route and `--vdj` write the aligner's Karlin–Altschul expectation values in AIRR `v_support`, `d_support`, and `j_support`, the anchored C E-value in `c_support`, and the posterior probabilities `v_posterior`, `d_posterior`, `j_posterior`, `c_posterior`, and `d_presence_posterior`. The same WebAssembly module computes them in Swig Web; IgBLAST is not invoked. See [V(D)J assignment](03_VDJ_ASSIGNMENT.md) for the definitions and why the support values are not IgBLAST E-values.

Reference/annotation behavior is intentionally four-way:

| Invocation | Metadata used by the aligner |
|---|---|
| no annotation file and no `--swigannots` | Headers are canonicalized, IMGT `.` gaps (and any `-` gaps) are removed, valid IUPAC bases are retained, and existing `SWIGMETA` is removed. The aligner then annotates the ungapped records itself (conserved motifs, a database consensus, and transfer from annotated relatives). FWR/CDR/junction fields are therefore filled wherever that annotation finds anchors, and `v_annotation_source`/`j_annotation_source` name its source. Because the gaps are removed first, IMGT gap numbering is not used on this route. |
| `-custom_internal_data` and/or `-auxiliary_data` | V `.ndm.imgt` entries are exact-identifier matches with 1-based inclusive FWR/CDR intervals converted to 0-based half-open `SWIGMETA` bounds. J `.aux` frame/CDR3-stop fields remain 0-based. Its optional fifth field is used to map the FWR4 end before each batch is written. A custom V file must cover every selected V identifier, matching IgBLAST's documented custom-data requirement. If J auxiliary data are given without custom V data, the embedded fixed species pack supplies validated V metadata in place of IgBLAST's installed organism directory. |
| `--swigannots` | The same species-first, progressively broadened and anchor-validated preprocessor used for a Swig Web upload. This is mutually exclusive with the two V/J IgBLAST annotation inputs. |
| `--prepared-reference MANIFEST` | Uses hash-verified, already-prepared V/D/J/C FASTAs and saved J FWR4 offsets from `prepare-reference`; no homology metadata search is repeated. |

`-d_frame_data` can additionally write exact 0-based D frame-one starts. Auxiliary/internal chain labels are checked against the locus inferred from each germline identifier. Duplicate `.aux` identifiers follow IgBLAST's map behavior (the last record wins); the official human file itself contains a duplicate. The route does not claim to reproduce IgBLAST's BLAST search, result ranking, or non-AIRR formatting. It reuses IgBLAST's public file/interface conventions to initialize the aligner's germline metadata.

The read kernel has no cross-read dependency, but total throughput includes single-threaded input/decompression, batch construction, worker messaging, ordered commit, and one output stream. On the WebAssembly route each worker also owns a complete WASM memory and germline index, and the pool compiles the module once and shares it with every worker (per-worker compilation competed with the workers for the CPU; under Bun 1.3.14 four workers then ran no faster than two). The native route shares one germline index across its threads. The pipeline's automatic worker count therefore stays at most eight, while `--vdj` uses every core. Workers take the next available batch rather than receiving round-robin queues, while the parent commits completed batches strictly in input order. Hot requests and responses are flat primitive/string messages so the standalone Bun runtime can use its optimized worker-message path. Bounded 1 MiB read/decompression chunks and an 8 MiB output high-water mark reduce callback/backpressure stalls without changing the AIRR bytes.

An omitted or zero `annotation.batchRecords`/`--batch-records` selects 2,000 records for one or two workers, 1,000 for three or four, and 500 for five or more. A positive value remains an exact override. This avoids starving a large worker pool on finite datasets without imposing small-batch overhead on a serial run. The measured 1–32-worker audit, memory curve, Amdahl estimate, and identified execution bottlenecks are recorded in [`BENCHMARK_CLI_SCALING_0.36.0.md`](../../BENCHMARK_CLI_SCALING_0.36.0.md); the byte-identical 0.36.1 before/after benchmarks are in [`BENCHMARK_CLI_THROUGHPUT_0.36.1.md`](../../BENCHMARK_CLI_THROUGHPUT_0.36.1.md). Both were measured with the SwiftIG core of those releases; measured comparisons of the 0.39.0 aligner with SwiftIG 0.38.11 are in [`BENCHMARK_SWIG_ALIGNER_0.39.0.md`](../../BENCHMARK_SWIG_ALIGNER_0.39.0.md).

Standalone releases use `--minify` without `--bytecode` (Bun 1.3.14 does not embed the worker entry correctly with bytecode in a multi-entry build). Linux x64 retains the existing baseline/Nehalem-compatible artifact and also publishes a Haswell/AVX2 modern artifact; the modern target affects only the JavaScript runtime. Neither changes the aligner's scoring or AIRR serialization: the embedded native aligner selects AVX-512BW, AVX2, NEON, or portable lanes at run time with byte-identical output, and the WebAssembly module is the same on every target.

## Audited browser-to-CLI option coverage

| Browser pipeline area | Portable config and CLI behavior |
|---|---|
| Input and study structure | Format, source/member name, path or pasted inline content, optional gzip-member range, dataset/sample/donor/cohort/timepoint/compartment, and study design are retained. |
| FASTQ preprocessing | Expected-error threshold, Phred offset, all 3′ trim settings, per-dataset reservoir size, and base seed run through shared browser code. |
| References | Exact browser-composed V/D/J/C FASTA is embedded unchanged. Hand-written `references.files` replacements default to the same local metadata preparation used by the browser; `references.prepareMetadata: false` is the explicit opt-out. |
| Assignment | Workers, batch size/default, calling profile, V/D/J identity floor, optional covered C-prefix identity floor, strand, AIRR mode, and every Double-D threshold are executed. |
| Allele refinement | Model, pooling boundary, segments, weighting, every advanced prior/kernel/inference option, and best-versus-confidence reassignment policy are executed before collapse. |
| Collapse/denoising | Exact key and every FAD/conservative/indel parameter, unresolved policy, scope, candidate cap, and constant-call partition are executed. Deterministic exact-key shards and independent denoising partitions use the configured worker pool without changing representatives or tie breaks. |
| CHMMAIRRa | Segment, BW/DB/auto choice, exact MSA, prior/mutation states/switch, DFR floor, filter threshold, unevaluated policy, and detailed Viterbi labels are executed. Row batches use the configured worker pool and are committed in ordinal order. Selected unaligned browser references are aligned with Kalign and embedded during export. A row whose call allele refinement reassigned is threaded against the allele the aligner aligned (`swig_original_v_call` / `swig_original_j_call`), whose alignment columns the row carries; before 0.39.4 threading failed on such rows, and with `retainUnevaluated` false they were dropped. |
| Repertoire selection | The complete interactive selection schema—not only the compact pre-run fields—is retained by Results export and applied by the CLI. |
| Lineage, SHM, missing-V evidence | Lineage scope/call policy/identity/productivity/candidate cap, SHM metric, and every two-pass missing-V threshold are executed. |
| Outputs | An explicit directory/prefix, batch-streamed annotated AIRR, chunk-written retained AIRR, method summaries, and optional lazy lineage-study bundle are controlled by config; `--workers` overrides JSON. |

The audit intentionally excludes display-only state, repertoire plots, ad hoc interactive queries, and lineage phylogenetics/UCA. Those are not unattended pipeline stages. The browser still exports the config from Results after they have been used, but the CLI stops at non-phylogenetic repertoire and lineage outputs.

## Lazy lineage-study bundle

When lineage assignment is enabled, the CLI writes two linked files:

1. an uncompressed AIRR TSV sorted by numeric lineage ID; and
2. a compressed lineage-study manifest containing exact references/options, dataset metadata, the bounded lineage summary table, optional lineage SHM mean/q95, and byte ranges into that AIRR file.

Every range stores its row count and SHA-256. The manifest also stores the complete AIRR size and SHA-256. Swig Web's **Load lineage study** first reads only the small manifest and checks the selected linked AIRR size. It can then search/filter the summary table without importing the AIRR rows. When a lineage is selected, the browser reads exactly that byte slice, verifies its row count and range SHA-256, and constructs a temporary result store containing only those rows. The normal lineage alignment, FastTree, and phylogenetic-UCA workbench then operates on that small store.

This route deliberately avoids the full-session metadata-index rebuild. It does not make arbitrary whole-repertoire filtering available. Opening one lineage reads one range; selecting several checkboxes and choosing **Open together** reads and verifies only those ranges, then appends their rows to one temporary store for the ordinary combined alignment/tree/UCA workbench. A later selection replaces that temporary store. The manifest contains at most the same 10,000 largest lineage summaries retained by the normal interactive explorer; assignment itself still covered every eligible record when the CLI ran. The manifest also carries each lineage's representative CDR3 and, when SHM was enabled, the exact lowest-SHM member's CDR3.

The linked AIRR must be uncompressed because byte offsets into gzip streams are not independently seekable without an additional block index. It is an analysis interchange file, not a replacement for the ordinary final AIRR export.

## Failure and reproducibility behavior

- A malformed or unsupported config is normalized only where a documented default exists; missing input or unusable references fail the run.
- With `references.prepareMetadata: true`, malformed identifiers, duplicate names, locus/segment mismatches, invalid IUPAC symbols, and invalid embedded `SWIGMETA` fail before the aligner indexes the references. A V/J record for which no transfer validates is retained without invented `SWIGMETA` and reported as a warning. The aligner's own germline annotation may then derive its anchor and regions, and `v_annotation_source`/`j_annotation_source` report that source.
- CHMMAIRRa requires an aligned reference MSA. Browser export builds and embeds it before download when selected references are unaligned; an already aligned explicit MSA is retained byte-for-byte.
- Pipeline output requires a dedicated explicit directory; assignment-only output requires an explicit file or `-` stdout destination. The resolved pipeline config records all inherited defaults.
- A lineage range is rejected if its byte count, row count, or SHA-256 does not agree with the manifest.
- The standalone binary has no network dependency during analysis unless a future explicit input option says otherwise.

## Engineering references and relationship

- [Bun standalone executables](https://bun.sh/docs/bundler/executables) provide the runtime embedding, asset embedding, worker bundling, and per-platform compilation used for distribution. Bun does not define any biological method.
- [AIRR Rearrangement schema](https://docs.airr-community.org/en/latest/datarep/rearrangements.html) supplies interoperable field conventions. The JSON config, lineage-sort order, byte-range manifest, range hashes, and lazy browser store are Swig-specific.
- [NCBI's IgBLAST setup documentation](https://ncbi.github.io/igblast/cook/How-to-set-up.html) defines the roles of `.ndm.imgt`, `.aux`, `-custom_internal_data`, `-auxiliary_data`, and `-d_frame_data`; the [official 1.22.0 release source/package](https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/1.22.0/) supplies their exact coordinate comments and parser behavior. Swig consumes those metadata conventions but does not implement the IgBLAST search engine.
