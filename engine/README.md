# swig: the V(D)J aligner of Swig

This crate is Swig's V(D)J aligner. Swig Web runs it as WebAssembly
(`public/swig.wasm`) in its browser workers, swig-cli runs the same module
in its workers, and it is also available natively as the `swig` binary. It
assigns V, D, J and, optionally, C germline alleles to reads by seeded,
banded, affine alignment with placement-DAG adjudication, and writes Swig's
AIRR Rearrangement TSV.

It originated as vdj_assign, a Rust implementation of the SwiftIG / Swig
0.38 method note. Its decision rules were then re-derived from one
log-likelihood-ratio objective (see [docs/ALGORITHM.md](docs/ALGORITHM.md))
instead of tuned patches. Every component is a replaceable, statically
dispatched trait, and every search shortcut is either lossless or measured
on simulated repertoires.

The crate and binary are `swig` 0.39.3. The workspace member `wasm/` (crate
`swig_wasm`) is the WebAssembly module. There are no dependencies; Rust ≥
1.89, edition 2024. The alignment kernel runs on AVX-512BW, AVX2, NEON,
WebAssembly SIMD or portable lanes, chosen at run time (on WebAssembly, at
build time), and calls do not depend on the kernel.

- [docs/ALGORITHM.md](docs/ALGORITHM.md): the model, the search and its
  guarantees, structured for a methods section.
- [docs/ASSESSMENT.md](docs/ASSESSMENT.md): the critical assessment behind
  the current design. It covers what was ad hoc, what replaced it, the
  measurements, and the remaining plan.
- [docs/REVIEW.md](docs/REVIEW.md): a code review of soundness and design,
  what it fixed, and the open findings ranked.

## Install

The engine has no releases of its own. Each
[Swig release](https://github.com/mashu/swig/releases) attaches the native
aligner as one executable per platform, built by
[`.github/workflows/release-cli.yml`](../.github/workflows/release-cli.yml)
(`cargo build --release --locked --bin swig`) and run with `--version` on
its own platform:

| Asset | Platform |
|---|---|
| `swig-aligner-linux-x64` | Linux x86-64 |
| `swig-aligner-linux-arm64` | Linux AArch64 |
| `swig-aligner-macos-x64` | macOS Intel |
| `swig-aligner-macos-arm64` | macOS Apple silicon |
| `swig-aligner-windows-x64.exe` | Windows x86-64 |

The binaries are not signed. Make a downloaded file executable
(`chmod +x swig-aligner-linux-x64`); on macOS, also clear the quarantine
flag before the first run (`xattr -d com.apple.quarantine
swig-aligner-macos-arm64`). The same releases attach `swig-cli-*`, Swig's
command line, which carries the WebAssembly build of this engine and does
not need the native binary.

## Build and test

From the root of the Swig repository:

```sh
npm run test:engine        # = cargo test --profile checked --workspace --manifest-path engine/Cargo.toml
                           #   unit, annotation, loci, robustness, synthetic, simulated-accuracy
                           #   and WebAssembly-ABI tests, on release code with debug assertions
                           #   and overflow checks (the `checked` profile); then builds
                           #   engine/target/release/swig for the CLI tests
cargo test --release --workspace --manifest-path engine/Cargo.toml -- --ignored
                           # only the ignored test: accuracy at scale, 20,000 reads per preset
cargo build --release --manifest-path engine/Cargo.toml    # engine/target/release/swig

rustup target add wasm32-unknown-unknown
bash scripts/build-swig-wasm.sh                            # public/swig.wasm (npm run build:swig-wasm)
```

[`scripts/build-swig-wasm.sh`](../scripts/build-swig-wasm.sh) builds the
`swig_wasm` crate for `wasm32-unknown-unknown` with `+simd128` and the
`wasm` profile (fat LTO, stripped) into `engine/target/wasm-simd`, and
copies the module to `public/swig.wasm` (or `$SWIG_WASM_OUTPUT`).

## Command line

```sh
swig --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta reads.fastq > calls.tsv

swig --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --c IGHC.fasta \
    --profile sensitive-d --strand forward --threads 8 reads.fastq > calls.tsv

zcat reads.fastq.gz | swig --v V.fasta --d D.fasta --j J.fasta - > calls.tsv

# one database of several loci: every rearrangement stays within one locus
cat IGHV.fasta IGKV.fasta IGLV.fasta > V.fasta
cat IGHJ.fasta IGKJ.fasta IGLJ.fasta > J.fasta
swig --v V.fasta --d IGHD.fasta --j J.fasta reads.fasta > calls.tsv

# Swig's double-D screen: AIRR rows unchanged, VDDJ evidence in a sidecar
swig --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta \
    --double-d all --double-d-out double_d.tsv reads.fasta > calls.tsv

# optional: write the annotated database for review or curation
swig prepare --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta [--c IGHC.fasta] [--chain heavy|light] \
    --out prepared/

# simulate reads with known truth (IgSim model), then score any tool's AIRR TSV against it
swig simulate --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta [--c IGHC.fasta] \
    --preset hard --reads 10000 --seed 1 [--rc 0.5] [--shm RATE] [--dd P] --out sim
swig --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta sim.fasta > calls.tsv
swig evaluate --truth sim.truth.tsv calls.tsv              # works for IgBLAST etc. too
swig evaluate --truth sim.truth.tsv calls.tsv other.tsv    # paired: 95 % CIs, exact McNemar tests
```

`--v` and `--j` are required; `--d` and `--c` are optional. Assignment
options (`swig --help` prints them):

| Option | Meaning |
|---|---|
| `--profile NAME` | `r` (default), `truth`, `igblast` or `sensitive-d`; also `r_optimized`, `truth_optimized`, `igblast_compatible`, `sensitive_d` (see *Calling profiles*) |
| `--threads N` | worker threads (default: all cores) |
| `--chain heavy\|light` | chain type for annotation (default: voted per locus from the database) |
| `--strand both\|forward\|reverse` | strands a read is assigned on (default `both`, chosen per read) |
| `--min-identity F` | drop V, D and J hits whose local alignment, before the stitched outer flank, is below identity F (0–1; default none) |
| `--no-anchor-window` | do not add the conserved-anchor D window |
| `--no-posterior` | skip the V, D, J and D-presence posteriors (~30 % faster) |
| `--calls per-read\|fixed` | the model the call maximises: each read's own (default; the model its posteriors come from) or the fixed `(2, −3)` scores |
| `--d-allele map\|max-marginal` | report the call's D allele (default) or the one of largest posterior marginal given the call's D presence and topology (needs posteriors) |
| `--d-presence-rate R` | D prior from the D search's own chance-hit count at D-presence rate R, 0 < R < 1 (see *The model in brief*) |
| `--d-open-flank on\|off` | charge a D that ends (starts) its path for the J (V) the read would show in the bases after (before) it (default `on`; see *The model in brief*) |
| `--min-d-match N` | shortest exact run of an admissible D (default: the profile's) |
| `--d-mismatch S`, `--j-mismatch S` | D or J mismatch score, negative (default: the profile's) |
| `--min-j-length N` | fewest aligned J columns of a J hit (default none) |
| `--c-evalue E` | call C when its anchored E-value is at most E (default 10⁻³; needs `--c`) |
| `--c-prefix-identity F` | report C only where the germline from its first base reaches identity F (default 0, off; above 0 needs `--c`) |
| `--double-d off\|all\|long-span` | Swig's double-D screen (default off); needs `--double-d-out FILE`, and `all` and `long-span` need `--d` |
| `--double-d-min-span N`, `--double-d-seed-length N`, `--double-d-pseudo-trim N`, `--double-d-pseudo-mismatches N`, `--double-d-min-gain S` | the screen's evidence rules (defaults 40, 11, 5, 3, 8) |
| `--simd NAME` | alignment kernel: `avx512`, `avx2`, `neon`, `portable` or `scalar` (default: the fastest this CPU runs, AVX-512 excepted) |

`simulate` presets are those of IgSim.jl: `easy`, `mid` (alias `train`) and
`hard`. `--rc` reverse-complements that fraction of reads (default 0.5),
`--shm` fixes every read's SHM rate, and `--dd` sets the tandem-DD
fraction. It writes `PREFIX.fasta` and `PREFIX.truth.tsv`.

**Reads** may be FASTA, FASTQ or an AIRR Rearrangement table. The format is
detected from the first non-blank character (`>`, `@`, anything else a
table). A table is tab-separated, or comma-separated without quoting when
its header holds no tab, and is read from its `sequence`, `sequence_id` and
`quality` columns: AIRR input is re-annotated from its sequences. A header
that starts with a tab starts with an empty column name, as pandas writes
its index column.
Multi-line FASTQ is supported; qualities are validated and kept for the
`quality` column. FASTQ sequence lines may not hold spaces, tabs or control
characters (FASTA drops them, but a FASTQ quality counts each base); no-call
symbols such as `.` read as N. Reads are streamed in batches of 16,384, so
memory stays bounded on large files, and `-` reads stdin. Gzip input is
recognised and rejected with a `zcat` hint, since the crate has no
dependencies.

**Repeated sequences** are assigned once: within a batch, and across
batches through a cache of the last 4,096 distinct sequences' rows. A row
depends only on the read's sequence, besides its identifier and qualities,
so the cache changes no output.

The command line rejects an option that takes a value given twice (a
repeated flag changes nothing), a FASTQ file passed as a germline database,
and `--double-d` without `--double-d-out` and the reverse. An option that
filters or adds to the output is refused without the database it needs, as
a usage error before any database is read: `--c-evalue` and
`--c-prefix-identity` above 0 without `--c`, and `--double-d all` or
`long-span` without `--d`. An option that only tunes a stage the run does
not have is run and named in a warning, since it changes nothing: the D
options (`--min-d-match`, `--d-mismatch`, `--d-presence-rate`, `--d-allele`,
`--no-anchor-window`) without `--d`, the screen's rules without
`--double-d`, and `--d-allele max-marginal` with `--no-posterior`. A
profile is a bundle and is never reported. With the screen on, D alleles
shorter than `--double-d-seed-length` are named too, since the screen
cannot see them. The WebAssembly module refuses the same options (an empty
D or C text is no database) and adds the warnings to its summary
(`RunOptions::check_databases`, `Run::warnings`). A debug build prints a
warning, because it is 10–50× slower. The database summary, warnings and
throughput go to stderr.

**Germline databases** are FASTA and may hold several loci. IMGT reference
headers (`accession|allele|species|…`) are read as the allele name, and
identical sequences are merged and reported under all their names. Anchors
and V region boundaries come from the headers where present and are
detected otherwise, on every start (see *Germline annotation and regions*).
Names should read `GENE*ALLELE`, as IMGT's do: the part before `*` is the
gene. V retrieval keeps candidates by gene (3 genes, up to 16 within the
leader's vote noise), so a V name without `*` is given a gene by sequence:
it joins the gene of its closest `*`-named relative when their local score
is at least 0.95 of the smaller self-score (about 2 % divergence), else the
names without `*` group among themselves at the same level, a group taking
its first name as gene (`gene=` in a prepared header). Reported names do
not change. On human IGH at 30 % SHM, renaming every `*` to `_` cost 8.2
points of V allele accuracy before this grouping and +0.10 (not
significant) with it. `prepare` and assignment warn when more than 16 V
names hold no `*` and have no close relative.
As with the reads, a leading UTF-8 byte-order mark is skipped and a byte
that is no UTF-8 (a Latin-1 description) reads as U+FFFD.

## Output: Swig's AIRR schema

Every row has the same 111 tab-separated columns (`AirrRecord::COLUMNS`),
header first: Swig's 106 columns in the order Swig's interface has always
read, then the engine's five posterior columns.

| Columns | Content |
|---|---|
| 1–10 | `sequence_id`, `sequence`, `quality`, `sequence_aa`, `rev_comp`, `productive`, `vj_in_frame`, `stop_codon`, `complete_vdj`, `locus` |
| 11–22 | `v_call`, `d_call`, `j_call`, `c_call`, `sequence_alignment`, `sequence_alignment_aa`, `germline_alignment`, `germline_alignment_aa`, `junction`, `junction_aa`, `np1`, `np2` |
| 23–36 | `cdr1`, `cdr2`, `cdr3`, `fwr1`, `fwr2`, `fwr3`, `fwr4`, each followed by its `_aa` |
| 37–53 | `v_`, `d_` and `j_` `score`, `identity`, `support`, `cigar`; `c_score`, `c_identity`, `c_prefix_identity`, `c_support`, `c_cigar` |
| 54–69 | `v_`, `d_`, `j_` and `c_` `sequence_start`, `sequence_end`, `germline_start`, `germline_end` |
| 70–83 | `cdr1`, `cdr2`, `cdr3`, `fwr1`–`fwr4` `_start` and `_end` |
| 84–91 | `v_`, `d_`, `j_`, `c_sequence_alignment`, then `v_`, `d_`, `j_`, `c_germline_alignment` |
| 92–99 | `junction_length`, `junction_aa_length`, `np1_length`, `np2_length`, `v_frameshift`, `j_frameshift`, `d_frame`, `sequence_frame` |
| 100–106 | `region_definition`, `v_annotation_source`, `j_annotation_source`, `v_alternatives`, `d_alternatives`, `j_alternatives`, `c_alternatives` |
| 107–111 | `v_posterior`, `d_posterior`, `j_posterior`, `c_posterior`, `d_presence_posterior` |

- **Reads.** `sequence` is the read as given, upper-cased, `U` as `T`, any
  symbol that is no IUPAC nucleotide as `N`; `quality` holds its FASTQ (or
  AIRR) qualities as given. A read without an identifier is `sequence_N`,
  N its 1-based position among the reads the run has read: in the input for
  `swig`, but in a WebAssembly instance among those of all its
  `swig_annotate` calls. Swig Web and swig-cli, which spread an input over
  several instances, therefore name such reads by input position before
  batching. A read without a significant V or J (see *Significance* below)
  keeps these columns and `rev_comp`, with empty calls.
- **Coordinates** are 1-based closed intervals on the *oriented* query (the
  reverse complement of `sequence` when `rev_comp` is `T`) and on the
  germline. CIGARs use `S`/`N` for unaligned query and germline bases and
  `=`/`X`/`I`/`D` for columns. Flags are `T`/`F`, empty when unknown.
- **Calls.** `v_call`, `d_call`, `j_call` and `c_call` list every
  co-optimal allele, comma-separated: each placement of that kind whose clip
  to the called interval is worth exactly the call under the call's model.
  The best-supported allele that spans the called read bases is first, and
  the segment's coordinates and alignment are its. `*_score` is the clipped score,
  `*_identity` matches over alignment columns (a fraction), and `*_support`
  the Karlin–Altschul E-value `K·m·n·e^{−λS}` of the reported span over the
  read length `m` and the database length `n`, the formula of the V/J gates
  (`c_support` is C's anchored E-value instead).
- **N/P and stitched alignments.** `np1` follows V and `np2` follows D, as
  AIRR names an N/P region by the segment it follows.
  `sequence_alignment` and `germline_alignment` stitch V, D and J in read
  order, N/P bases against `N`.
- **Regions.** `fwr1`–`fwr4`, `cdr1`–`cdr3` and `junction` are projected
  from the called V's and J's germline coordinates. The coordinates come
  from `SWIGMETA=` (or header keys) when the reference gives them, and from
  the engine's own annotation otherwise; `v_annotation_source` and
  `j_annotation_source` name the source (a `SWIGMETA` provenance such as
  `IMGT-gapped`, or `header`, `imgt-gaps`, `motif`, `consensus`,
  `transfer:NAME`, `low-confidence`, `missing`). `region_definition` is
  `IMGT` when the called V has an anchor or region starts. `productive`
  means in frame, free of stop codons and of frameshifts; `d_frame` and
  `sequence_frame` are 1-based; `complete_vdj` means the V starts within 3
  bases of its germline start and the J ends within 3 of its germline end.
- **Alternatives.** `*_alternatives` lists the alleles the call does not
  name whose placement of the same kind, clipped to the called interval,
  scores within `max(2, |S|/100)` of the call's score `S`, best first, as
  `name|score|identity|sequence_start|sequence_end|germline_start|germline_end`
  entries separated by `;` (coordinates as above).
- **Posteriors.** `*_posterior` is the posterior probability that the
  reported names include the read's allele up to germline equivalence, and
  `d_presence_posterior` the probability that the read has a D (see
  *Posterior probabilities* below). They print with six decimals, or ten
  within 10⁻⁴ of 0 or 1. `--no-posterior` leaves the V, D, J and D-presence
  columns empty; `c_posterior` comes with every C call.

**Constant region (isotype).** The `c_` columns are filled only with `--c`
(model: *Constant region* below). `c_support` is the anchored E-value
`N_C·7·e^(−λS)` over the C genes of the J's locus and the 7 start
diagonals, and C is called at `c_support ≤ 10⁻³` (`--c-evalue`).
`c_prefix_identity` is the identity of the C germline from its first base:
the germline bases the alignment leaves uncovered, compared base for base
with the read bases before it, then the alignment's columns.
`--c-prefix-identity F` empties the C columns of a row whose
`c_prefix_identity` is below F. Give each locus' spliced coding C regions,
CH1 onward (IMGT or codon-aligned records such as abstar's both work), since
every gene of the locus counts. For a short C part shared by several genes
(IGHG1/2/4 open alike), `c_call`'s *first* name is any of them: read the
isotype from the genes listed.

## Double-D screen

`--double-d all` or `--double-d long-span` runs Swig's double-D screen
(`double_d.rs`; Swig's method note:
[04_DOUBLE_D_SCREEN.md](../public/methods/04_DOUBLE_D_SCREEN.md)). It leaves
every AIRR row unchanged: the DAG calls VJ or VDJ, and the screen reports
evidence for a second D in a sparse sidecar TSV (`--double-d-out FILE`).

A read is screened when its call places V and J, its locus has D genes
(IGH, TRB, TRD) or is unknown, and the V–J span (the bases between the
called V's end and J's start) is at least two seed lengths; in `long-span`
mode the span must also reach `--double-d-min-span` (40). The screen's
evidence rules:

- **Seeds.** Each D hit holds an exact `--double-d-seed-length`-word (11;
  6–24) of a D allele of the read's locus, or a `max(6, L − 3)` rescue word
  whose ungapped extension reaches `L` aligned bases. Seeds extend without
  gaps (match +2, mismatch −2) to their best-scoring prefix on each side;
  hits on one allele diagonal count once. A D allele shorter than the seed
  length is not indexed, and the run names it in a warning: at the default
  11 the screen cannot see TRDD1 (8 nt) or TRDD2 (9 nt), so it finds no
  TRD tandem.
- **Gain.** Two ordered, non-overlapping hits score at least
  `--double-d-min-gain` (8) above the best single hit.
- **Pseudo-tandem veto.** The read from D1's start to D2's end, shortened
  by `--double-d-pseudo-trim` (5) bases, must be more than
  `--double-d-pseudo-mismatches` (3) Hamming mismatches from every equally
  long substring of a D allele of the locus (the IgScout Δ-distance test),
  so that no single D explains it.

Pairs are tried in evidence order (pair score, gain, matched bases, shorter
NP2, names). The first that passes is reported, and other passing pairs
within two score units are listed as alternatives (up to 12). A hit's call
lists every allele with the same score and coordinates.

The sidecar has 38 columns: `swig_batch_record_index` (the read's 0-based
index in the input; in the WebAssembly ABI, in the batch), `sequence_id`,
`locus`, `v_call`, `standard_d_call` (the AIRR row's `d_call`), `d_call`,
`d2_call`, `j_call`; for D1 and D2 the score, identity and CIGAR, the
sequence and germline start and end, and the sequence and germline
alignments; `np2` and `np3` (the bases between D1 and D2, and between D2
and J) and their lengths; then `swig_double_d_vj_span`, `_seed_length`,
`_pair_score`, `_best_single_score`, `_score_gain`, `_pseudo_distance`
(empty when the trimmed span is shorter than a seed or longer than every D
allele, where the veto does not apply), `_mode` (`all` or `long_span`) and
`swig_double_d_alternatives`, `;`-separated `D1+D2|pair score|gain|pseudo
distance` entries (−1 for no distance).

## As a library

```rust
use swig::prelude::*;

let database = prepare(v_set, Some(d_set), j_set, AnnotateOptions::default())?; // dedup + anchors
let config = CallingProfile::ROptimized.config();
let assigner = DefaultAssigner::from_config(&config, database.v, database.d, database.j)?;
let mut workspace = assigner.workspace();                  // one per thread
let call = assigner.assign(&encode(b"CAGGTG..."), &mut workspace);
let calls = assign_all(&assigner, &queries, 0);            // parallel, input order
```

or a whole run as the command line and the WebAssembly module perform it,
from FASTA text to AIRR TSV (`swig::run`; options by their command-line
names, without `--`):

```rust
let mut options = RunOptions::default();
options.set("profile", Some("sensitive-d"))?;
options.set("double-d", Some("all"))?;
options.set("no-posterior", None)?;                        // a flag takes no value
let texts = GermlineTexts { v: &v_fasta, d: Some(&d_fasta), j: &j_fasta, c: None };
let database = PreparedDatabase::from_fasta(texts, options.annotate())?;
let mut run = Run::new(&options, database)?;
let mut sidecar = Vec::new();
let tally = run.write_airr(reads.as_bytes(), &mut out, Some(&mut sidecar))?;
// or, batch by batch without headers: let rows = run.assign(&records);
//   rows.airr, rows.double_d, rows.tally; run.airr_header(), run.double_d_header()
```

`RunOptions::NAMES` lists the 24 option names and whether each takes a
value; `RunOptions::config` maps them onto the profile's `EngineConfig`.

## In Swig (WebAssembly)

The `swig_wasm` crate (`wasm/src/lib.rs`) exposes the engine as a small C
ABI over linear memory, without a bindings generator. It is built for
`wasm32-unknown-unknown` with `simd128` and needs no imports: Swig
instantiates it with `WebAssembly.instantiate(module, {})`. The module,
`public/swig.wasm`, is about 1 MB (0.33 MB gzipped). Swig Web's aligner
workers (`src/aligner-worker.ts`, `src/aligner-compute-worker.ts`) and
swig-cli's workers (`cli-src/swig-worker.mjs`) each hold one instance
through the TypeScript wrapper
[`src/aligner-engine.ts`](../src/aligner-engine.ts).

Each instance holds one engine. Text crosses the boundary as UTF-8 in
buffers from `swig_alloc`; a function borrows the buffers it is given, and
the host frees them with `swig_free`. A function returns −1 on failure, with
the message in the error buffer. Results stay readable until the next call.

| Export | Does |
|---|---|
| `swig_alloc(len) → ptr`, `swig_free(ptr, len)` | buffers for the host to fill and free |
| `swig_set_option(name, name_len, value, value_len) → 0 \| −1` | set an option by its command-line name without `--` (a flag with an empty value); the run is rebuilt at the next use; a host that wants the defaults back starts a new instance |
| `swig_init_database(v, v_len, d, d_len, j, j_len, c, c_len) → records \| −1` | load the germline FASTA texts (V and J required, D and C empty when absent) and prepare them; the result buffer holds `{"summary": …, "warnings": […]}` as JSON |
| `swig_annotate(reads, len) → reads \| −1` | assign FASTA, FASTQ or AIRR-table text; the result buffer holds the AIRR TSV, header first, then one row per input record in input order |
| `swig_result_ptr()`, `swig_result_len()` | the result buffer |
| `swig_double_d_result_ptr()`, `swig_double_d_result_len()`, `swig_double_d_count()` | the double-D sidecar of the last annotation (header first, a row per screened read) and its row count |
| `swig_error_ptr()`, `swig_error_len()` | the last error message |
| `swig_panic_message()` | after a trap, leave the panic's message in the error buffer; the host then discards the instance |
| `swig_version()` | leave `{"version": …, "backend": …}` in the result buffer |

- **Same engine, same rows.** The module runs `run::Run`, as the native
  binary does: the same option names and validation, the same database
  preparation and warnings, the same rows. Every batch carries its own
  header line, which the wrapper splits off. The crate's test
  `the_abi_runs_the_library_run` checks that the ABI's rows equal
  `Run::write_airr` for FASTA and AIRR-table input. In the standalone
  project, the module's TSV equalled the native command line's byte for
  byte on 45,000 rows; a posterior can still differ in its last printed
  digit where it lands within a last-bit difference of `exp` (WebAssembly's
  portable libm against the platform's) of a rounding tie (2 of 260,000
  rows in a test, every call identical).
- **The wrapper.** `AlignerEngine.create(module, settings, germlines)` maps
  Swig's settings to options (`profile`, `strand`, `min-identity`,
  `c-prefix-identity`, the `double-d` options, and the `min-d-match`,
  `d-mismatch`, `j-mismatch` and `min-j-length` overrides that
  `swig-cli --vdj` uses), loads the germlines and reads the summary.
  `annotate(reads, count, doubleD)` checks that the module returned one row
  per record of the batch and returns the header, the rows and, with the
  screen on, the sidecar.
- **SIMD.** WebAssembly has no run-time feature detection, so the kernel is
  chosen at build time: the `simd128` backend, 16 lanes in two `v128`
  registers. Swig ships only this build, so a host without WebAssembly SIMD
  cannot instantiate it; a build without `simd128` would use the scalar
  kernel, since emulated lanes there run several times slower.
- **swig-cli runs the native binary.** Since Swig 0.39.1 each standalone
  `swig-cli` executable embeds the native `swig` for its platform, extracts it
  to the user cache and streams records through it (stdin FASTA/FASTQ/AIRR,
  stdout AIRR TSV, `--threads N`, the same option names). `Run::write_airr`
  flushes the rows of every `BATCH` (16,384 records) as soon as they are
  assigned, so the CLI keeps more than that many records in flight. The
  WebAssembly module stays the fallback and the route for the double-D
  sidecar; the rows are the same, up to the last-digit posterior exception
  above. The Linux release binaries link glibc statically
  (`-C target-feature=+crt-static`) and the Windows one its C runtime
  (`.cargo/config.toml`), so neither needs a newer system library than
  swig-cli itself.
- **Threads.** `wasm32-unknown-unknown` has none: `assign_map` runs every
  chunk on the calling thread, and Swig parallelises by running one
  instance per worker. In the standalone project's measurements, one
  worker under V8 assigned 1,100–1,450 simulated reads/s, about 60 % of the
  native AVX2 kernel on one thread.

## Databases of several loci

One database may hold several loci, as Swig composes them: IGH, IGK and IGL
for a BCR study, or the TR loci (TRA, TRB, TRD, TRG). Each allele's locus is
found in this order (`locus.rs`, `locus_infer.rs`, `prepare.rs`):

1. a `locus=` header key (`locus=IGK`, `locus=TRA,TRD`), as `swig prepare`
   writes it;
2. its name (`IGHV1-2*02`, `TRBJ2-7*01`: the first locus code that starts a
   word and is followed by a segment letter, case-insensitively; for IGH
   also an isotype letter, `IGHG1`, `IGHT2D`, `IGHZ`, `IGHO`, `IGHY`), else
   a gene name in the header's free text. Shared TRA/TRD V genes
   (`TRAV14/DV4`) belong to both loci;
3. for a V, J or C allele whose name places it nowhere (`novel_1`,
   `sample7_J1*01`), the locus of its closest relatives: the database's own
   named alleles of the kind, else an embedded panel of 269 human and mouse
   V, J and C family representatives. A relative must align clearly (a local
   score of at least 30, and 0.3 of the allele's score against itself for
   the database's alleles; for the panel's 0.1 for V and C, 0.5 for the
   short J); every locus within 15 of the best counts. As the panel is of
   another species for most databases, a J or C it places keeps only the
   loci of the database's V alleles, and a V it places weakly in TRA alone
   or TRD alone takes both (the two loci share V genes);
4. otherwise the loci the database names: a D allele (too short to place by
   sequence) the loci with D genes (IGH in a BCR database, TRB and TRD in a
   TCR one), a V, J or C allele every locus the database names. In a
   database of several loci this is warned about once per kind, listing the
   alleles. Only an allele of a database that names no locus at all keeps
   every locus, which is warned about too.

`swig prepare` writes the result as `locus=` (with `locus_source=relative`
or `fallback` for the inferred ones), so a reload does not infer again. Steps
3 and 4 cost nothing when every name places its locus, as in every shipped
reference. **Give custom-named alleles a `locus=` key** when the warnings
appear: the fallback cannot tell a custom TRB D from a TRD D, so TRDD alleles
then compete on TRB reads. Swig Web and swig-cli write `locus=` for the
custom records of a reference cell or of `prepare-reference` themselves, and
`swig-cli --vdj` keeps the locus a header gives (`locus=TRB`, or a locus word
such as `TRB`).

- **One locus per rearrangement.** J candidates must share a locus with the
  read's anchoring V placements (every locus when no V placement anchors),
  and D candidates with the anchoring V and J. When the anchors hold
  several loci, as on an scFv read, the J alleles of each locus are searched
  and kept on their own, so one chain's J cannot crowd out the other's. When
  every placement of a read shares a locus, the read is adjudicated once.
  Otherwise each locus is adjudicated over its own placements, all at the
  divergence of the whole read, and the call is the best of their paths in
  the model the call maximises (`--calls`), with or without posteriors; an
  exact tie goes to the path whose last segment comes first in the read (by
  kind, then position), as within one locus. The posteriors are
  the mixture over loci weighted by their likelihoods `Z`, so under
  `--calls fixed` the called locus can hold less than half of the posterior.
  C joins the locus the path's V, D and J share (an allele of unknown locus
  constrains nothing), else its J's. `locus` reports the locus the path's
  alleles share, empty when it is unknown (a fallback to several loci
  counts as unknown).
- **Per-locus annotation.** The V and J alleles of each locus are annotated
  on their own: chain type, FR3 length model and Cys23 priors differ between
  heavy and light chains and between TR loci.
- **Per-locus priors and gates.** A D database of several loci (TRB and TRD)
  prices a D by its own locus' D genes (`P_D`, below), and the C gate's
  `N_C` counts the C genes of the read's locus only. Alleles placed by the
  fallback (step 4) can be called in every locus they pair with but count
  towards no locus that has alleles of a definite locus.

`tests/loci.rs` assigns IGH and IGK reads against one IGH + IGK database:
every called allele belongs to the read's locus, and no IGK read gets a D,
also when the D (or J) alleles have custom names.

## Germline annotation and regions

V and J germline genes are in frame with their conserved anchors — V Cys104
and J Trp/Phe118 (IMGT numbering) — so an anchor fixes the reading frame, and
the V region boundaries follow from three conserved residues:

```text
FR1  = [frame start, Cys23 + 4 codons)      CDR1 = [Cys23 + 4, Trp41 − 2)
FR2  = [Trp41 − 2, Trp41 + 15)              CDR2 = [Trp41 + 15, FR3 start)
FR3  = [Cys104 − (L − 1), Cys104 + 1)       L = FR3 length in codons
CDR3 = [Cys104 + 1, J Trp/Phe118)           FR4  = [J Trp/Phe118, J end)
junction = [Cys104, Trp/Phe118 + 1)         (both conserved codons included)
```

These are the IMGT/AIRR definitions. FR2 has no IMGT gaps. FR3 (IMGT 66–104)
has no conserved residue at its start, and its length varies: 38 codons in
IG heavy and 36 in kappa, but by gene family in lambda (36/38), TRA (32–34),
TRB (37/38), TRG and TRD (37–39). The CDR2/FR3 boundary therefore comes from
annotated relatives or a length learned from the database (below).

### Header keys and `SWIGMETA`

All positions are 0-based in the germline sequence; unknown keys are ignored.

| Key | Meaning |
|---|---|
| `anchor=` | first base of V Cys104 / J Trp-Phe118 |
| `fr1=` `cdr1=` `fr2=` `cdr2=` `fr3=` | V region starts |
| `frame=` | offset (0–2) of the first complete codon (a V or J anchor also fixes it; a D needs it for `d_frame`) |
| `aliases=A,B` | other names of this sequence |
| `SWIGMETA=` | Swig's prepared coordinates (below) |

Swig prepares every reference it ships, and every FASTA a user loads, with
region coordinates transferred from IMGT-annotated relatives
(`src/germline-preprocess.ts`), and writes them into the header as
`SWIGMETA=` followed by 13 comma-separated integers (`metadata.rs`):

| Index | Value |
|---|---|
| 0 | coding-frame start: offset (0–2) of the first complete codon, −1 if unknown |
| 1 | J only: last base of CDR3 (the Trp/Phe118 codon starts one base later), else −1 |
| 2–11 | V only: FWR1, CDR1, FWR2, CDR2 and FWR3 as `[start, end)` pairs, −1 if unknown |
| 12 | provenance: 1 `IMGT-gapped`, 2 `AIRR-C`, 3 `IMGT-boundary-transfer`, 4 `validated-J-motif`, 5 `provided`, 6 `J-anchor-transfer` (other values `unclassified`) |

FWR3 ends after the Cys104 codon, so a V anchor is FWR3's end − 3, and a J
anchor is one past CDR3's last base. `SWIGMETA` values are given values,
like header keys: they fill what the explicit keys leave open (when the
result stays consistent) and win over everything the engine would derive,
with one exception. A given anchor is checked against the allele's own
motifs: a J anchor off `[WF]G.G` while a confident motif lies elsewhere in
its frame, or a V anchor with weak motifs in another frame than all of its
given region starts, is *contradicted*. A `SWIGMETA` anchor then gives way
to the motif hit, and the allele reports the annotation's provenance; an
explicit `anchor=` (or a `SWIGMETA` anchor no motif can replace) is kept but
flagged and counted in the uncertain-anchor warning. Swig's shipped
references hold such anchors (TRAJ24*01–03 and TRAJ48*01, among 50 J
allele records, sit 15–21 bases early; IGLV2-11*03, IGLV2-8*03,
IGHV4-30-4*05, TRBV5-8*02 and TRBV11-3*03 were given another frame), so
reads of those alleles used to get a wrong junction.
A malformed field is an error. The provenance becomes the allele's
`v_annotation_source` or `j_annotation_source`; an allele without one
reports the annotation's own source.

### IMGT-gapped input

An IMGT-gapped V set (the IMGT reference format: `.` for the gap columns) is
numbered from its gap columns. The FR1, CDR1, FR2, CDR2 and FR3 starts are
IMGT columns 1, 79, 115, 166 and 196, and Cys104 is column 310.

- The gaps are used only when most gapped alleles are codon-aligned and
  most of those that hold IMGT codon 23, 41 or 104 carry Cys, Trp or Cys
  there. A codon that no gapped allele holds (Cys104 when all are
  3′-partial) is not tested, but one of the three must be. Otherwise the
  gaps are ignored (a codon alignment in other columns) and the summary
  line says why.
- A region that the numbering places before a 5′-partial sequence stays
  empty; a motif hit cannot fill it.
- An allele whose numbering ends before Cys104 gets no anchor (`3′-partial`).
- `prepare` writes both facts into the ungapped headers (`before=fr1,cdr1`,
  `partial=3`), so a reloaded database keeps them.
- A one- or two-column interior gap (an IMGT frameshift marker) is noted, and
  the columns after it still map exactly.
- Header keys and `SWIGMETA` win over the gap columns.
- Ungapped alleles appended to a gapped set are annotated as below, and take
  their FR3 start from their gapped relatives.

### Sources, in order of trust

A given value is never overridden. With several loci, each locus is
annotated on its own.

1. **Header keys and `SWIGMETA`**, then the **IMGT gap columns** (source
   `imgt-gaps`).
2. **Motif search** on a one-pass three-frame translation. For V, Cys23, Trp41
   and Cys104 are placed jointly under spacing constraints, scored by
   conserved context: `[ST]C`; `W` with `Q`+3, `P`+5, `G`+6 (`WVRQAPG` /
   `WYQQKPG`); and `YYC`, with `Y` weighted above `F`/`H`. Upstream stop codons
   are penalised, which also selects the frame. For J, `[WF]G.G` with an open
   FR4.
3. **Database consensus.** Once pass 1 has run, the chain type is voted:
   J `W`/`F` and Trp41 context `W.R`/`W.Q`. That turns the Trp41 → Cys104 window
   into a CDR2-length constraint, `[14 + L, 26 + L]`, over every FR3 length `L`
   the chain default or the given FR3 starts admit. The dominant V 3′ tail
   (applied to confident alleles only when it holds ≥ 50 % of the database) and
   the Cys23 codon index of full-length alleles (only in their frame) become
   priors for a second pass.
4. **Sibling alignment.**
   - An FR3 start that is not given is taken from the most similar allele
     whose FR3 start is given. If the two share a Trp41 → Cys104 spacing, the
     Cys104 → FR3 distance is copied; otherwise the FR3 start is mapped by
     alignment, so a CDR2 indel stays in CDR2.
   - A Cys23 window with several Cys codons, or none (a mutated Cys23, or a
     5′-partial allele that starts after it), takes CDR1 from the most
     similar allele whose Cys23 is settled, in any frame.
   - An allele whose motifs stay weak (pseudogene, mutated Cys) gets its anchor
     and every region start transferred by alignment from the most similar
     trusted allele.
   - A relative must share ≥ 50 % of the target's distinct 8-mers, so a chance
     match can never donate.
   - When the relative's anchor codon maps at or past the allele's 3′ end
     (an IMGT 3′-partial sequence), the allele holds no Cys104 and gets no
     anchor. Its weak motif hit would be another codon: on human IGH, 22
     partial V alleles (IGHV4-31\*06–\*09, IGHV2-70\*02, …) had been given an
     out-of-frame TGT ~30 nt upstream, which shifted every junction they
     produced.
5. **Database models**, for alleles that no relative settles. Both are
   learned once per database (per locus), with Jeffreys (½) pseudo-counts:
   - *FR3 length.* Trp41 → Cys104 spans `d = 14 + c + L` codons for a CDR2 of
     `c` and an FR3 of `L` codons. The model picks the `L` that maximises
     `P(L) · P(c = d − 14 − L)`, both estimated from the given FR3 starts. `L`
     ranges over the observed lengths ±1 codon, up to 39. Without any given
     FR3 start, the chain default (38 heavy, 36 light) is used.
   - *Cys23.* Each candidate codon is scored by the database's Cys23
     nucleotide offsets, its CDR1 lengths and the two residues before
     Cys23, against background frequencies. TRB Cys23 follows `R`, not the
     IG `[ST]`. A mutated Cys23 is placed only at an offset the database
     has, with posterior ≥ ½.
6. **Low-confidence motif hit (no relative), or missing.** Both are reported.

Every derived start that may be wrong is flagged in `annotation.tsv` as
`uncertain=…`; given values never are. The flagged cases are:
- every start of a low-confidence or missing anchor;
- an FR3 length with model posterior below 0.95, or outside the model;
- the chain default in a TR or IGL database without any given FR3 start
  (`prepare` and assignment also warn);
- a CDR1 from an ambiguous Cys23 that no relative settled, or from no Cys23;
- a CDR of unusual length;
- a start transferred from a relative whose own start is flagged.

### Speed

- **Codon table:** translation is one lookup per codon in a 64-entry table on
  the 2-bit codon index.
- **Candidate scoring:** a V allele costs a few dozen candidate scores per pass.
- **Relative search:** an inverted 8-mer index (CSR over 65,536 keys, built
  lazily) is searched only for an unresolved allele or an FR3 start that is
  not given while others are, touching only alleles that share k-mers.

| Database | Preparation time |
|---|---|
| Human IGH, IMGT-gapped (415 V) | 1.1 ms |
| Human IGH, ungapped | 29 ms |
| IGH-sized, anchors in headers | 0.8 ms |
| IGH-sized, worst case (no recognisable motifs) | 1.3 ms |
| 3,400 alleles, worst case | 11.7 ms |
| 3,000 synthetic alleles with motifs, including transfers | ~40 ms |

Per read, regions are projected lazily when a call is reported, only through
the called V and J alignments. Each boundary is one binary search in the
clip index. This costs under a microsecond per read, against ~450 µs for
assignment and 11.8 µs for the whole AIRR record, and nothing is added to
assignment itself.

### Reliability (tests/annotation.rs)

Synthetic genes are built from heavy and light framework templates:
- random CDR1/CDR2 of realistic lengths over all 20 residues;
- 12 % mutation of non-conserved framework residues, which creates decoy Cys
  and Trp;
- random synonymous codons;
- 5′ truncation of 1–2 nt, which shifts the frame.

| Test | Result |
|---|---|
| 300 heavy + 300 light alleles | every boundary exact |
| 12,000 alleles over 40 databases | every anchor and region start exact (3 flagged CDR1 errors before the Cys23 model) |
| IMGT-gapped heavy and light sets with a Cys104 → Ser pseudogene, a 3′-partial, a 5′-partial starting in CDR1, a frameshift marker and a header override | every anchor and start exact; no FR1 or CDR1 before the sequence |
| TRB (FR3 37/38) and TRA (34/32) sets, 3 headers per family | every boundary exact; without headers `fr3_warning` and every FR3 start flagged |
| TRB and TRA sets with one header per family, 5′-truncated alleles | 0 errors (6 flagged CDR1 errors before the Cys23 model) |
| Ungapped alleles appended to an IMGT-gapped 36/38 light set; a CDR2 deletion | FR3 from the gapped relative, exact |
| Truncated TRB allele with an upstream `SC` decoy; Cys23 → Tyr; 5′-partial after Cys23 | CDR1 exact |
| Pseudogene with Cys104 and Trp41 destroyed | anchor and all regions transferred exactly from its sibling |
| 3′-partial allele ending inside Cys104, with an out-of-frame TGT decoy | no anchor, flagged `3′-partial` |
| 200 simulated reads, half reverse-complemented | FR1–FR4, CDR1–CDR3, junction, frame, stop codon and productive flag exact |

### `prepare`

`swig prepare` writes `V.fasta`, `D.fasta` and `J.fasta` with every
anchor, frame and region start in the headers as keys (plus
`annotation=SOURCE`, and for a partial IMGT-gapped record `before=` the
starts before its first base and `partial=3` an end before Cys104),
`C.fasta` with `--c`, and `annotation.tsv` with one row per allele:

- frame, anchor and its residue;
- region starts, CDR1 and CDR2 lengths in codons;
- source (`header`, `imgt-gaps`, `motif`, `consensus`, `transfer:NAME`,
  `low-confidence`, `missing`) and notes. The notes end with an
  `uncertain=cdr1,fr3`-style token naming the starts that may be wrong.

With `--c`, `prepare` and assignment warn about database conventions that
keep C from being called: a C allele missing more than 3 bases at its 5′ end
against its gene, a C record over 2 kb or without a stop-free reading frame
(a genomic record with introns, not the spliced coding region), and a J
allele whose bases past its anchor differ from most J alleles' by more than 3
(bases past the splice donor move the J exit off the splice).

Reloading the written FASTA uses the headers as-is, so manual corrections
stick (a `SWIGMETA` provenance is not written, so reloaded alleles report
`header`). Assignment runs log a one-line summary. They warn when any allele
is transferred, low-confidence or missing, and when a TR or IGL locus gives
no FR3 start, so that every FR3 start rests on one chain default.

## Architecture

```text
 chunk of reads ─► per read: strand from germline k-mer support; V candidates (gene frames)
                ─► V stage, whole chunk: score batch, sibling bounds, traceback waves
                ─► per read: J window right of V, J alleles of each anchoring V's locus ─► J stage, whole chunk
                ─► per read: D windows from V/J geometry, D candidates of the locus (L-mer bitsets)
                ─► D stage, whole chunk ─► per read: admission gates
                ─► per read: placement DAG over VJ / VDJ, per locus ─► anchor check ─► Call
                ─► with --c, whole chunk: C anchored at the J exit (one job per reach class)
                ─► AIRR row; with --double-d, the double-D screen ─► sidecar row
```

Every per-read decision is taken per read. Batch aligners must answer
each job exactly as a single call would, so a read's result never depends on
its chunk: `assign(q)` is literally a chunk of one, and this is tested. The
stages are the methods `stage_orient`, `stage_v`, `stage_j`, `stage_d` and
`stage_adjudicate` of `Assigner`.

| Module | Responsibility |
|---|---|
| `alphabet`, `fasta`, `germline`, `reads` | 2-bit codes (N = 4), germline FASTA and sets with anchors, streaming FASTA/FASTQ/AIRR-table reads |
| `locus`, `metadata` | Loci from allele names and which alleles may join; Swig's `SWIGMETA=` field |
| `translate`, `annotate`, `prepare`, `read_regions` | Codon table, anchors and V regions per locus (header/`SWIGMETA` > IMGT gaps > motif > consensus > relatives > database models), per-read FR/CDR/junction projection |
| `scoring`, `stats` | Log-odds substitution profiles (uniform or position-specific), affine gaps, Karlin–Altschul λ/K/E-values |
| `orient` | Strand from germline k-mer presence on both strands, with a significance test |
| `seed`, `index`, `retrieve/frames` | Rolling k-mers, CSR k-mer index over gene frames (alleles of a gene in one coordinate system) |
| `retrieve/seed_vote`, `retrieve/budget` | Radix-sorted diagonal votes per frame, gene budget, candidate bands |
| `retrieve/exhaustive`, `retrieve/lmer_mask`, `retrieve/ungapped` | J: every allele; D: lossless L-mer bitsets; optional ungapped pre-ranking |
| `band`, `align/gotoh`, `align/inter_seq`, `align/lanes` | Diagonal bands; scalar banded Gotoh; inter-sequence Gotoh written once over a SIMD lane trait (AVX-512BW 32 lanes; AVX2, NEON, `simd128`, portable 16) |
| `segment`, `family_bound` | One search stage: score → sibling bounds along the gene frame → traceback waves → acceptance and E-value gates → outer extension |
| `clip`, `evidence`, `placement` | O(1) clip-and-rescore at any cut; exact-run / identity evidence |
| `window`, `admission` | Where J and D are searched; D evidence gates |
| `constant`, `retrieve/reach` | Constant region after adjudication: splice window, anchor-column DP, anchored E-value and posterior; one candidate per class of alleles the read cannot reach apart |
| `adjudicate/{lattice,dag,posterior,units,model,decode,ties,oracle}` | The placement lattice; Viterbi and forward–backward over it; the per-read model in integer units both read; grammar and log-odds path model; max-marginal D decoding; the co-optimal alleles a call reports; the exhaustive oracle the tests check the DAG against |
| `flank` | Outer-flank extension of V and J alignments |
| `pipeline`, `config`, `profiles` | Stage-wise composition with per-locus adjudication, declarative configuration, calling profiles |
| `run` | A whole run, shared by the command line and the WebAssembly module: options by their command-line names, the assigner they build, the repeated-sequence cache, reads streamed to AIRR TSV |
| `report`, `double_d`, `batch` | Swig's AIRR schema; the double-D screen and its sidecar; parallel batches |
| `sim/{random,params,simulate,synthetic,truth,evaluate}` | IgSim port, synthetic germlines, truth records, accuracy evaluation |
| `wasm/` (crate `swig_wasm`) | The WebAssembly module's C ABI |

## The model in brief

The full description, with the guarantee of each search step, is in
[docs/ALGORITHM.md](docs/ALGORITHM.md).

- **Scores are log-odds.** `(2, −3)` substitution scores are the
  log-likelihood ratio, in units of 1/λ nats with λ = 0.634, of "templated by
  this germline at 11.2 % divergence" against "random sequence". N and P bases
  score 0 under both hypotheses. Affine gaps price indels.
- **The path objective** is `Σ clipped segment scores − P_V·#V − P_D·#D − P_J·#J`,
  maximised exactly over the topologies VJ and VDJ (and the partial V…, …J
  and D…J), over the placements, and over both cuts between neighbours (the
  earlier one exits at `x`, the later enters at `y ≥ x`, and the bases
  between are N/P). `P_V` and `P_J` price a V and a J as `P_D` prices a D,
  over their own genes, so a read that starts after its V or stops before
  its J does not get one from a few chance bases. Ties go to fewer D, then
  to more germline support (the summed unclipped scores). An AIRR row holds
  one D; Swig's double-D screen reports a second one beside the call and
  leaves the call alone.
- **Calls come from the read's own model.** Each path is valued at the read's
  divergence `p̂`, estimated before adjudication from the V and J columns no
  path contests: a match is worth `ln 4(1−p̂)` nats, a mismatch `ln(4p̂/3)`, a
  V, D or J costs `P_k/2` matched columns of its kind. Values are integers (1/2¹⁶ nat), so the
  maximum is exact, and the posterior reads the same integers: the call is a
  mode of the reported distribution. `--calls fixed` keeps the fixed scores,
  of which the per-read model at 11.2 % divergence is an exact multiple.
  Paired against fixed-scale calls on 20,000 reads per preset, it left
  4.0–4.3 pp more D-free reads without a D on mid sets and improved the
  D-presence Brier score on every set (ASSESSMENT §9.8).
- **`P_D` is derived, not tuned.** A D is one of `N` genes at one of `L`
  junction positions. Including it pays off when `λ·S_D > ln(N·L) − ln(prior
  odds)`, which gives `P*_D` = 10.55 for human IGH (32 genes, L = 25, even
  odds) and 6.17 for TRB, kept unrounded; on integer scores the fixed scale
  charges `⌊P*_D⌋`, which decides `S > P*_D` exactly. A D database of several
  loci gives each locus its own `P*_D`. `--d-presence-rate π` counts what the
  D search tries instead, `P*_D = [ln(K·m·n_D) − ln(π/(1 − π))]/λ` over
  alleles × offsets × positions. At π = 0.9 it improved the D-presence Brier
  score on all six sets of 20,000 reads, but that rate was found by
  measuring; at the simulator's own rate (0.94) D-free specificity fell
  14–26 pp. The rate has to come from the user, so it is not the default
  (ASSESSMENT §9.9).
- **A D at a path's open end pays for the J (V) it implies.** A path that
  ends on a D claims a junction whose J the read does not show in the bases
  after it (a path that starts on a D, a V missed before it). Past 24 free
  bases (N/P and the D's unaligned tail), the D pays
  `ln(1 + ½·(1 − ε)/ε)` nats, where `ε` is the probability that the missing
  segment's visible bases stay below its stage's accept score at the read's
  mismatch rate (beta-binomial over the anchored V and J columns). `ε` never
  drops below the floor `ε_min = 10⁻³`, so the cost is at most 6.2 nats even
  when the read's J (V) is present but missing from the database; a locus
  whose database part has no J (V) at all pays nothing on that side, and a
  V or J too weak to pay for itself cannot open or close a D's flank. With
  a derived D prior (the default profile `r`, and `sensitive-d`), a chance
  D far from the anchor of a D-free read goes (a V followed by 100 random
  bases: 93.5 % → 6.8 % D calls, `sensitive-d` 99.3 % → 9.5 %; 100 random
  bases then a J: 97.7 % → 10.2 %, `sensitive-d` 99.7 % → 18.2 %);
  full-length reads keep their calls. The profiles `truth` and `igblast`
  charge no D prior (`P*_D = 0`), so a chance D scoring a dozen points
  outweighs the 6.2-nat ceiling and most stay (`truth`: 93.7 % → 92.8 % and
  97.8 % → 96.2 %). A chance D within about 30 bases of the anchor sits
  where a real D sits and stays. Real N/P regions longer than about 24
  bases at a read that ends in its J lose some D calls (1.1 pp at 26–30
  bases). `--d-open-flank off` restores 0.39.3's calls (docs/ALGORITHM.md
  §2.3).
- **Significance.** A V or J hit is kept at E ≤ 1 against its database and
  window. A read is called only when its path contains an *anchor*, a V or J
  hit at E ≤ 10⁻³ whose clip on the path still scores that much, so random
  sequence receives no call.
- **Constant region.** C is spliced to J at a fixed phase, so C's base 0
  lies within `s = 3` bases of the J's exit. A bonus column in front of C
  turns every local aligner into an exact anchored one, and C is called at
  an anchored E-value `N_C·(2s + 1)·e^(−λS) ≤ 10⁻³`, after adjudication, so
  it never moves a V, D or J call. On 20,000 simulated reads per preset
  (human IGH with abstar's IGHC; synthetic): ≤ 2 false calls on ~4,100
  C-free reads, the primary C allele anchored-equivalent ≥ 99.6 % from 16 nt,
  `c_posterior` ECE ≤ 0.18 %, 9 % more instructions with `--c`
  (ASSESSMENT §11).
- **Posterior probabilities.** Read as a distribution, `P(path, cuts) ∝
  exp(nats)`, the same model gives calibrated probabilities by
  forward–backward over the same DAG; the other priors stay at their
  derivation λ. `--d-allele max-marginal` reports the D allele of largest
  posterior marginal instead of the call's (off by default: it gains only on
  heavily mutated reads, +0.25 pp on hard reads).

## Extension points

Every component is a trait. `Assigner<A, RV, RJ, RD, W, Adm, Adj>` is generic
over all of them, and the defaults form `DefaultAssigner`. Swapping one changes
the type, and the whole pipeline monomorphises. There are no trait objects and
no runtime type checks.

| Trait | Default | Alternatives shipped |
|---|---|---|
| `Aligner` (+ `score_batch`, `align_batch`) | `InterSequenceGotoh` (best `Backend` of the CPU) | `BandedGotoh`; batch methods default to single-job loops |
| `Retriever` | V: `UngappedPrefilter<SeedVoteRetriever>` (prefilter off); J: `ExhaustiveRetriever`; D: `LmerMaskRetriever` | `SeedVoteRetriever` alone |
| `CandidateBudget` | `AmbiguityBudget{3, 16, z = 2}` | `FixedBudget` |
| `SubstitutionModel` | `AffineScoring` | `HotspotScoring` (position-specific; the kernel's profile path) |
| `DWindowPolicy` | `WithAnchorWindow<ExpandedResidual>` | `ResidualWindow`, `ExpandedResidual` |
| `Admission` | `AnyOf<ExactRunFloor, Option<DistributedEvidence>>` (distributed-D off by default) | `AdmitAll`, `AllOf`, any `Option<A>` |
| `Adjudicator` | `PlacementDag<VdjGrammar, EngineeringObjective>` | your own |
| `Grammar`, `PathModel` | `VdjGrammar::VDJ`, `EngineeringObjective` | `VdjGrammar::VJ` (or any `max_d`); your own path grammar / objective |
| `VdjAssigner` (+ `assign_many`) | chunk-staged `Assigner` | default `assign_many` loops `assign` |

```rust
let experiment = DefaultAssigner::from_config(&config, v, Some(d), j)?
    .with_d_window(ResidualWindow { overlap_margin: 2 })
    .with_d_admission(AllOf(ExactRunFloor { min_run: 5 }, MyNewGate));
```

## Acceleration tricks

Measured single-threaded on simulated reads (IgSim `mid` preset, half
reverse-complemented, 450–650 nt) against a human IGH germline set of
541 V, 90 D and 16 J alleles, in a shared 4-core container (±5 %):

| Stage of optimisation | reads/s/thread |
|---|---|
| Starting point (method-note pipeline, AVX2, chunk batching) | 849 |
| + k-mer orientation, V retrieved on one strand, J searched right of V only | ~1000 |
| + band = seed cluster ± 8 diagonals (was ≥ 57 diagonals wide) | ~1230 |
| + gene frames and radix-sorted votes (V retrieval 267 → 73 µs/read) | ~1620 |
| + D floor `P_D + 1`, junction windows reclaim 12 instead of 24 | ~1830 |
| + uniform-profile fast path in the SIMD kernel | ~2020 |
| + clip index in three allocations instead of seven | **~2200** |
| + exact DAG over every entry cut (+3 % instructions) | ~1900 |
| + posterior decoding, on by default (`--no-posterior` skips it) | ~1300 |

Rows after the first seven were measured later in a busier container
(±10 %); in instructions the exact DAG costs 3 % and the posterior pass 30 %
of a hard-preset run (it also keeps every admissible D in the lattice).
Four threads gave ~7300 reads/s without posteriors and ~5400 with them. The
noise-based V gene margin with every k-mer as a seed (more accurate at high
SHM, see *V seeds and gene margin*) then added 4 % of instructions per read
on the hard preset and 6 % at 30 % SHM; wall-clock differences were within
the container's noise. On the IMGT human IGH reference, 20,000 hard reads,
one thread, by kernel backend (identical calls): AVX-512BW ~2340, AVX2
~1870, portable ~1500, scalar ~730 reads/s. Calls are identical across the
pure-speed steps; the steps that changed calls are listed in
[docs/ASSESSMENT.md](docs/ASSESSMENT.md) with their accuracy effect. These
timings were taken on the standalone engine.

### Data and retrieval

**Nucleotides and k-mers as integers.** Rolling 2-bit k-mer keys in `u64`, with
N resetting the window. Spaced seeds are data (a mask), not code.

**V seeds and gene margin.** Every 9-mer of the read votes for the gene
frames that hold it (weight 4), and 7-mers too (weight 1) while the leader
holds fewer than 25 nine-mer seeds, where a count's Poisson relative error
`1/√n` exceeds 20 %. A gene is kept while its votes are within
`z·√(2·4·L)` of the leader's `L`, two standard deviations (`z = 2`) of the
difference between two seed counts near the leader's: relative to the leader,
but proportionally wider where counts are small (mutated reads). Measured
jointly on human IGH (hard and mid presets, fixed 20 and 30 % SHM, 200-nt
3′ fragments) and the synthetic database, against the previous release
(fixed margin of 8 votes or 5 %, 16-vote trigger, closed syncmers): mean V
gene over the seven conditions 85.2 → 88.1 %, 30 % SHM 75.7 → 84.7 %, hard
preset 98.36 → 98.50 %. A proportional margin (10–30 % of `L`) loses 4–9
points at 30 % SHM; z above 2 adds at most 0.1 point (the 16-gene cap binds).

Sampled seeds (`SeedSampling`: closed and open syncmers, minimisers) remain
available and correct. Closed syncmers are context-free and have a window
guarantee (every unmutated 12-base stretch yields a shared 9-mer seed at
`w = 4`), but under the noise-based margin every k-mer beats them: 1.1
points of mean V gene (4.4 at 30 % SHM) for 0.4 % fewer instructions on the
hard preset (ASSESSMENT §8).

**Strand from a k-mer bitset.** A presence bitset over all 4¹¹ k-mers of the V
and J germline (512 KiB). One pass computes forward and reverse-complement
keys together, so each strand costs one bit test per base. Each read position
counts at most once, so repeats cannot inflate the evidence. The count is
trusted only when it exceeds its Poisson null (`λ₀ + 3√λ₀`); otherwise the
read falls back to seed votes on both strands. Only the chosen strand is
searched for V. `--strand forward` or `reverse` skips the test.

**Gene frames.** The alleles of a gene share almost all seeds. A frame maps
each allele into the gene's longest allele by the offset at which its shared
seeds agree; an allele with an internal indel, or any divergent allele, gets
its own frame. The k-mer index is built over frames with identical postings
merged. A seed shared by 20 alleles is one posting and one vote, while
allele-specific seeds still vote for their gene.

**CSR k-mer index.** Offsets and packed `u32` postings (frame, position).
Keys use a direct-address table up to 24 key bits, and a sorted key array with
binary search beyond. Seeds with more than 256 postings are skipped.

**Radix-sorted packed votes.** Each vote is one `u64`, `(frame · 2^b + diagonal
+ bias) << 8 | weight`, with `b` sized to the query. The sort key is ~21 bits,
so two 11-bit LSD radix passes order the votes in linear time. One scan then
yields the best diagonal cluster per frame. A seed repeated within a frame
does not vote (it cannot locate the diagonal).

**Bands from the seeds.** A candidate's band is its frame cluster's diagonal
range, shifted by the allele's frame offset and widened by 8 diagonals.
Narrowing from ≥ 57 to ~17–25 diagonals changed no call on reads with indels.

**Exhaustive J, small window.** The J set is tiny (10–20 alleles per locus)
and under heavy SHM keeps few intact seeds, so every J allele of the V's
locus is aligned, unbanded, in the read right of V (from 12 bases before the
earliest V end).

**Ungapped pre-ranking** (`retrieve/ungapped`): Kadane's maximum segment sum
along diagonals around the band centre ranks candidates before DP. It is
lossy, so it is off by default.

### Dynamic programming

**Banded Gotoh with band-only traceback.** O(n·w) time and O(n·w) trace bytes.
Row offsets are allocated once per call from the band's cell count.

**Branch-free inner loop.**
- LLVM lowers `i32::max` in this loop to branches that mispredict on real
  alignments.
- A micro-benchmark replaying one alignment hides this, because the branch
  predictor memorises it: 1.85 ns/cell there versus 6.2 ns/cell in situ.
- Maxima use bit arithmetic. Traceback codes are computed arithmetically,
  because they are ordered like the tie-break preference.
- Rows are pre-sliced, so the zipped loop has no bounds checks.

**Shortened loop-carried chain.** With `t = max(0, diag, I)`, the horizontal
gap state satisfies `D' = max(t + g_o, D + max(g_o, g_e))`.
- This is exact for any gap costs and is checked against a textbook
  three-matrix Gotoh oracle.
- The serial path per cell drops from three maxima to one add and one max.

**Score first, traceback in waves.** Score-only DP ranks candidates. Traceback
then runs in waves down the ranking until `top` alignments pass acceptance
(plus ties), so a high score that fails the length/identity floor does not
take a slot. With at most `min(top + 1, cap)` candidates the scoring pass is
skipped (tracing all of them is cheaper).

**Sibling bounds along the gene frame (`family_bound`).** The alleles of a
gene share a frame (`frame position = allele position + offset`), and V and J
score one representative per frame first. A sibling `A` aligned in the same
frame diagonals as its representative `R` is bracketed exactly:

- `S*(A) ≤ S*(R) + Δ(R, A)`, where Δ sums the best profile gain at each frame
  position both cover and the best substitution score where only `A` does;
- `S*(A) ≥ L`, the traceback of the request's best representative clipped to
  `A`'s span and rescored with `A`'s profile, in `O(|D|)` from prefix sums
  (`D`: the positions where the two profiles differ). Any traced path serves,
  so the strongest representative is aligned instead of scored.

The `top`-th largest lower bound (one per candidate) never exceeds the true
`top`-th score, so a sibling whose upper bound falls below it is skipped; one
whose bounds meet needs no DP; one whose lower bound reaches it is aligned
outright and its traceback reused. The rest is scored. The ranking from the
`top`-th score up is that of full scoring, so calls are byte-identical, at
9.0 % fewer instructions on 2000 hard reads and 8.3 % on mid. Under a length
or identity acceptance floor, siblings are pruned only below the score floor.

**Inter-sequence SIMD (`align/inter_seq`).**
- Each 16-bit lane holds the same cell of a *different* job. Every lane
  runs the ordinary recurrence, so traceback codes are the scalar ones and the
  scalar traceback walk is reused verbatim.
- Diagonal bands use band-relative coordinates, so all lanes walk the same
  grid although their bands start on different diagonals. Unbanded jobs (J, D)
  use a rectangular layout.
- **Uniform fast path.** When every lane's profile depends on base identity
  alone (every `AffineScoring` profile), the group stores reference codes,
  and a cell's score is one load, one comparison and one blend. Otherwise the
  five interleaved profile rows are selected by query-code masks, which is
  exact for any `SubstitutionModel`.
- Out-of-band and out-of-matrix cells read as H = 0, I = D = −∞, exactly like
  the scalar kernel. Saturating i16 arithmetic, with a per-job admission check;
  jobs that could approach the sentinel fall back to the scalar kernel.
- **One kernel, every CPU.** The recurrence is written once over
  `align::lanes::Lanes<L>`, a register of `L` lanes with an associated mask
  type (vector masks on AVX2/NEON/portable, `k` registers on AVX-512), and
  compiled per backend inside a `#[target_feature]` entry point:
  AVX-512BW (32 lanes, opt-in with `--simd avx512` because 512-bit
  instructions lower the clock on many CPUs), AVX2 (16), NEON (16 in two
  registers), WebAssembly `simd128` (16 in two registers, in builds with
  that feature) and portable Rust arrays (16) that the compiler vectorises
  for the target's baseline SIMD. `Backend::detect` chooses at run time;
  on WebAssembly without `simd128`, which has no SIMD to vectorise portable
  lanes into and where emulated lanes run about 4× slower than one scalar
  job, it chooses the scalar kernel.

**Stage-wise batching across reads.** A read has only a few V, J and D
candidates, too few to fill 16 lanes. So every stage runs over a chunk of 64
reads and hands the aligner one batch per wave. Keep policies, windows,
admission and adjudication stay per read.

### D search

**Lossless D prefilter.** Per-L-mer allele bitsets (one OR per query L-mer).
L is derived from the admission gates so it can never drop an admissible hit:
`L = min_d_match` with only the exact-run gate (the default), `L = 4` with the
optional distributed-D gate (≥ 80 % identity over ≥ 16 columns forces an
exact run ≥ 4, by pigeonhole).

**Score floor from the objective.** A clipped span of an optimal local
alignment never outscores the whole alignment. So a D scoring at most `P_D`
can only add ≤ 0 to a path, and the strict with-D-beats-no-D rule rejects it;
dropping a D leaves its neighbours' cuts valid. The floor `P_D + 1` is
therefore exact for a fixed-scale call. A D below it still carries posterior
mass and can win a per-read call, so the floor applies only to fixed-scale
calls without posteriors (`--no-posterior --calls fixed`): posteriors keep
every admissible D (+1 % of instructions), and so do per-read calls. With
several loci the floor is the lowest `P_D` of a locus that has D genes.

### Adjudication

**O(1) clip-and-rescore.** Prefix sums over alignment columns give the score
of any sub-alignment between two query cuts in O(1), with exact affine gap
accounting at the new edges.

**One lattice, two semirings, one model.** Viterbi (max-plus) and
forward–backward (sum-product) run the same recursion over the same lattice
and read the same per-read integer values (each kind's affine
log-odds map `s·u_r + v_r` quantised whole, priors as multiples of `u0`), so
the call is exactly a mode of the posterior: one slot per
placement, D count and first segment, with a value per canonical entry cut.
Placements are visited by segment kind, then read start, so a D can precede
a J whose alignment starts before it; only the cuts order them.
Clip scores are separable in (entry, exit), the N/P score is linear in the
gap, and the support tie-break is separable too, so the inner maximum (or
sum) over entries is a running one along the placement and the one over
exits `x ≤ y` another, read in O(1) per entry. On an edge, every entry
beyond the predecessor's end sees the same best exit (or exit mass), which
is added once to a running maximum (or difference array): an edge costs
O(overlap). Every pair of cuts is evaluated, and both passes equal
exhaustive enumeration (tested). The posterior keeps linear values with one
log-scale per slot, carries N/P weights relative to each cut so they never
exceed 1 (a gap between non-overlapping neighbours joins the log-scale),
walks cut weights as running products along the alignment, and evaluates no
transcendental function per cut. A read whose evidence would leave the
`f64` range (a placement spanning > 1 200 nats on one side, or an overlap
> 600 nats, far beyond any junction) gets no posteriors rather than wrong
ones.

### Parallelism and allocation

One workspace per thread means no DP allocation after warm-up. Work is
self-scheduled in chunks of 64 from an atomic counter over scoped threads;
the calling thread is one of the workers, and a target without threads
(WebAssembly) runs every chunk on it. Results come back in input order, and
reporting is formatted on the workers (`assign_map`).

## Simulation and evaluation

`swig::sim` is a dependency-free port of
[IgSim.jl](https://github.com/mashu/IgSim.jl). It has the same generative
model and the `easy`/`mid`/`hard` presets:
- uniform allele sampling, exonuclease trims, N additions, random flanks;
- gated SHM with single-base indels, and Illumina substitutions;
- optional tandem DD.

It extends IgSim in four ways:
- **Per-base germline origins.** Boundaries, junctions and germline-equivalent
  alleles are exact truth even across indels.
- **Optional reverse complementation.**
- **Optional constant region** (`simulate --c`). 80 % of reads continue
  through the J–C splice into a C allele, from the base that follows the J
  in the reading frame (so IMGT and codon-aligned C records both work), for
  the whole 3′ flank or a uniform part of it. It draws from a random stream
  of its own: without `--c`, reads and truth are byte-identical, and with it
  their V(D)J part is.
- **`SyntheticDatabase`.** An IGH-like germline with a family → gene → allele
  similarity hierarchy and conserved anchors, plus constant regions whose
  subclasses share CH1's opening, for tests that cannot ship IMGT-derived
  data.

`sim::Truth` holds one read's ground truth and reads and writes it as TSV.
`sim::Assignment` is a call, built from an `AirrRecord` or from any tool's
AIRR Rearrangement TSV row. `swig simulate` and `swig evaluate`
expose both, so other callers can be scored on identical reads.
`sim::Evaluation` scores calls against the truth:
- an allele call is *equivalent* when the primary allele contains the
  trimmed, pre-mutation germline body (no read can separate such alleles);
- *gene*, *listed* (truth among the reported names) and *exact*;
- boundary errors, D presence by remnant length, junction exactness;
- with C, by C length: calls on reads without C, the primary name's
  *anchored* equivalence (the allele reads the simulated bases from its own
  splice offset; containment anywhere would make almost any allele
  equivalent to a short C part), the first name's gene, exact gene sets,
  over-calls and the oracle's resolution.

`sim::Paired` compares two callers read by read on the same truth: for V, D
and J allele equivalence, D-free reads left without a D call, and junction
exactness, it reports both accuracies, their difference with a 95 %
interval, and an exact McNemar p-value from the reads only one caller gets
right. D allele accuracy alone rewards calling D everywhere, so the D-free
column belongs next to it. D is also split by the true remnant's length
(1–3, 4–7, 8+ nt; below `min_d_match` = 4 a remnant cannot be found), and the
D allele and D-presence posteriors are compared by Brier score with a paired
t statistic. On 20,000 reads per preset, for example, the residual
adjudicator (removed since) gained +0.4–0.6 pp of D accuracy (p < 10⁻⁴) only
by giving 21–25 pp more D-free reads a D; see
[docs/ASSESSMENT.md §9.5](docs/ASSESSMENT.md) for the full table.

The measurements below were made on the standalone engine, on single-locus
human IGH databases; [docs/ASSESSMENT.md](docs/ASSESSMENT.md) gives the
protocols, and its §4 the accuracy before and after the redesign (on human
IGH, random reads called fell from ~30 % to 0, and hard-preset junctions
rose from 98.39 to 99.55 %).

At 20,000 reads per preset on the public IMGT human IGH reference (the one
shipped in the `airrship` PyPI package: 415 V / 37 D / 13 J after merging
identical sequences, ORFs and pseudogenes included; 95 % intervals ≤ ±0.5 pp):

| Preset | V | V gene | D | J | Junction | D-free reads given a D | Random reads called |
|---|---|---|---|---|---|---|---|
| easy | 99.88 | 100.00 | 99.69 | 99.94 | 100.00 | 0.0 % | 0.04 % |
| mid | 98.61 | 99.97 | 93.11 | 99.11 | 99.97 | 24.5 % | 0.04 % |
| hard | 93.23 | 99.61 | 85.91 | 94.25 | 99.65 | 16.0 % | 0.04 % |

**Against IgBLAST 1.22** (Bioconda build, defaults, `-organism human` with
the same germline as BLAST databases; the same simulated reads, all forward
strand, 5000 per preset, scored by `swig evaluate`):

| Preset | Tool | V | V gene | D | J | Junction | D-free reads given a D | Throughput, 4 threads |
|---|---|---|---|---|---|---|---|---|
| easy | swig | **99.88** | **100.00** | **99.29** | 99.84 | **100.00** | 0 % | **5200 reads/s** |
| easy | IgBLAST | 99.12 | 99.78 | 98.86 | **99.86** | 99.77 | 0 % | 67 reads/s |
| mid | swig | **97.88** | **99.72** | 93.35 | **98.74** | **99.97** | **23 %** | **5100 reads/s** |
| mid | IgBLAST | 97.18 | 99.62 | **93.73** | 98.60 | 99.72 | 43 % | 66 reads/s |
| hard | swig | **89.44** | 98.84 | **83.00** | **92.80** | **99.69** | **16 %** | **5200 reads/s** |
| hard | IgBLAST | 89.16 | **98.90** | 81.75 | 91.84 | 98.07 | 30 % | 71 reads/s |

IgBLAST's higher mid-preset D accuracy comes with twice the D calls on D-free
junctions. Throughput is wall clock for the whole command, posteriors on,
measured with the previous V retrieval defaults (the current ones use 4 % more
instructions).
IgBLAST takes J anchors from its auxiliary file, which lacks three of this
database's J alleles; the comparison adds them (frame and CDR3 end from
`swig prepare`, which agrees with the file on all 13 shared alleles).
Without them IgBLAST reports no junction for 18 % of reads.

The examples measure and diagnose; without `--v`/`--j` they use the
synthetic database:

```sh
# accuracy and throughput on simulated reads (also --shm, --seed, --rc, --profile,
# --simd, --tail and the ablation switches listed in examples/benchmark.rs);
# --random N adds a negative control: N random reads with the simulated length
# distribution, reported as the share with a path, a V call and a J call
cargo run --release --manifest-path engine/Cargo.toml --example benchmark -- \
    [--v V.fasta --d D.fasta --j J.fasta] [--preset easy|mid|hard] [--reads N] [--threads N] [--random N]
# search errors vs model errors against an exhaustive oracle
cargo run --release --manifest-path engine/Cargo.toml --example diagnose -- \
    --v V.fasta --d D.fasta --j J.fasta --preset hard [--reads N] [--show K]
# where D errors arise: adjudication, not found, outside the D windows, eaten by J or V
cargo run --release --manifest-path engine/Cargo.toml --example d_errors -- \
    --v V.fasta --d D.fasta --j J.fasta [--preset hard] [--reads N] [--no-anchor-window] [--j-reclaim N]
```

`examples/diagnose` separates errors by aligning every allele unbanded under
the same objective. If that exhaustive optimum is right but the pipeline is
wrong, it is a *search* error; if the optimum itself is wrong or tied, it is
a *model* error or tie.

| Human IGH, 2000 reads | correct | ties | search errors | model errors |
|---|---|---|---|---|
| V, mid | 97.4 % | 1.3 % | **0** | 1.3 % |
| J, mid | 98.8 % | 0.6 % | **0** | 0.7 % |
| V, hard | 88.9 % | 5.5 % | 0.4 % | 5.3 % |
| J, hard | 93.4 % | 2.8 % | 0.1 % | 3.7 % |

Nearly all remaining errors are where SHM mimics an allele's distinguishing
SNPs, so no score-based caller can do better.

**Posterior calibration** (same reads; expected calibration error in %, then
Brier score). The event is "some reported name is germline-equivalent to the
truth" for V, D and J (the call lists all co-optimal alleles), and "the read
has a simulated D" for D presence:

| Preset | V | D | J | D presence |
|---|---|---|---|---|
| easy | 0.04 / 0.0008 | 0.49 / 0.0015 | 0.07 / 0.0000 | 0.20 / 0.0029 |
| mid | 0.29 / 0.0088 | 1.29 / 0.040 | 0.30 / 0.0059 | 1.98 / 0.026 |
| hard | 1.13 / 0.043 | 1.16 / 0.072 | 0.69 / 0.031 | 1.24 / 0.024 |

With one fixed λ instead of each read's divergence, and only co-optimal
placements pooled, mid-preset ECEs were 4.6 / 7.7 / 3.9 / 2.1.

**Where D errors arise** (`examples/d_errors`, D remnants ≥ 8 nt): on the hard
preset 87.5 % of reads get the right D; 10.2 % get another allele although
the right one was a candidate (similar D genes under SHM), 1.8 % never had
the right candidate, and none lose the D to a J alignment that covers it
(0.00 %). No true D lay outside the D windows, with or without the anchor
window; searching J over the whole read lets 0.05 % fall outside.

## Validation

- **Posterior exactness:** `marginals_equal_exhaustive_summation` compares
  ln Z, every placement marginal and the D-count distribution with exhaustive
  summation over paths and canonical cuts (VJ, VDJ, VDDJ, a per-base N/P cost,
  a per-column evidence term, an N column, an insertion), and checks that the
  Viterbi score is the best path score; `strong_junction_costs_stay_in_range`
  does the same with N/P costs of −4 and −8 per base across a 300-base V,
  `long_costly_gaps_keep_their_mass` across a 300-base gap that every path
  crosses, and `overlaps_beyond_the_f64_range_are_refused` checks that a
  540-base overlap gets no posteriors; `per_read_call_is_a_mode_of_the_posterior`
  checks that with one integer model the marginals are those of `e^{units/F}`
  and the Viterbi key their exhaustive maximum; `posteriors_are_calibrated`
  bounds V/J calibration error on IgM- and IgG-like reads.
- **DAG optimality:** `dag_equals_brute_force_over_all_paths_and_cuts` and
  `dag_equals_brute_force_on_hypermutated_igsim_reads` enumerate every
  VJ/VDJ/VDDJ path and every pair of canonical exit and entry cuts, with and
  without a per-base insertion cost and with the lexicographic support
  tie-break, and require the same optimum; the IgSim case fails on a DAG that
  tries entries only inside the overlap. Constructed reads need two cuts
  (`two_cuts_leave_a_shared_mismatch_stretch_to_n`), an entry beyond the
  overlap (`entries_beyond_the_overlap_are_tried`), and a D before a J whose
  alignment starts earlier (`a_d_can_precede_a_j_that_starts_earlier`).
  `per_read_dag_equals_exhaustive_enumeration` repeats the enumeration in the
  per-read model's integer units at five divergences (800 optima, 135 of 640
  not the fixed-scale path; at 11.2 % the fixed path exactly), and
  `conditioned_dag_equals_exhaustive_enumeration` checks the masked,
  topology-restricted pass behind max-marginal D decoding.
- **Kernel exactness:**
  - The scalar kernel equals a textbook Gotoh across three gap regimes.
  - Every traceback rescores to its reported score.
  - Every backend this CPU runs (AVX-512, AVX2, NEON, portable) equals the
    scalar kernel on 160 random mixed batches of up to five registers, and on
    degenerate jobs (empty sides, all-N queries, one-row queries, bands that
    miss the matrix). Its
    scores, end cells and complete tracebacks are compared over full and
    banded layouts, N bases, indels, four scoring schemes, and hotspot profiles
    mixed into half of the batches. Both the uniform and the profile path are
    covered.
  - The whole chunk-batched SIMD pipeline equals the per-read scalar pipeline
    on 400 simulated reads (paths, orientations and every placement).
- **Simulated accuracy (`tests/simulated.rs`):** IgSim presets on the
  synthetic database, half reverse-complemented. Accuracy floors sit ~1.5 points
  below the measured values for orientation, V/D/J allele and gene, junctions
  and long D remnants. Random reads must receive no call. `accuracy_at_scale`
  (ignored by default; `--ignored`) runs 20,000 reads per preset with floors
  ~3 standard errors under the measured rates: 500-read samples have 95 %
  intervals of ±1–4 pp and hide most regressions.
- **Robustness (`tests/robustness.rs`):** empty, 1-base, all-N, IUPAC,
  low-complexity, 20 kb, chimeric reads and every prefix and suffix of a read
  through C, through every profile with the C stage, with AIRR-record
  invariants; C without a J, V–J reads ending in nothing, one base or N, and
  adapters, homopolymers and splice donors after every J are never called C.
- **Constant region:** the anchor-column alignment equals a textbook anchored
  Gotoh on 500 random cases through every backend; reach-class members align
  identically; V, D and J calls, placements and posteriors are identical with
  and without the C stage, and batched SIMD C calls equal per-read scalar
  ones; anchored truth, C-length strata and the simulator's reading frame are
  tested, and `tests/simulated.rs` bounds C accuracy, false calls and
  `c_posterior` calibration.
- **Loci (`tests/loci.rs`):** one IGH + IGK database; every rearrangement
  stays within the read's locus, each locus is annotated on its own, and no
  IGK read gets a D.
- **The run and the module:** `run.rs` checks that options are read by
  their command-line names, that a run writes the header and one row per
  read equal to `AirrRecord::from_call`, that repeated sequences within and
  across batches reuse their rows under their own identifiers and
  qualities, that the double-D screen leaves the AIRR rows unchanged and
  writes well-formed sidecar rows for tandem-D reads, and that options are
  refused or reported against the databases a run has (a table of every
  refusal and warning, and the D alleles the screen cannot see).
  `wasm/src/lib.rs` checks that the ABI's rows equal `Run::write_airr` for
  FASTA and AIRR-table input, that options rebuild the run and bad ones are
  refused, and that an option without its database is refused at load or
  named in the summary's warnings; `tests/cli.rs` that the command line
  refuses one as a usage error before reading any database. `metadata.rs` checks that
  `SWIGMETA` fields are read and malformed ones rejected.
- **Platforms:** in the standalone project's CI the suite passed on x86-64
  and AArch64 Linux, with AArch64 output equal to x86-64 output byte for
  byte, and on `wasm32-wasip1` with `simd128` under Wasmtime, where the
  kernel tests compare the `simd128`, portable and scalar kernels. Swig's CI
  runs the suite natively (`npm run test:engine`) before building
  `public/swig.wasm`, then tests the built module in Node.js
  (`tests/swig-wasm.test.mjs`: profiles, FASTA/FASTQ/AIRR input, heavy,
  light and TR loci, strands, constant regions, independence from batching,
  the repeated-sequence cache and the double-D screen).
- **Synthetic calling (`tests/synthetic.rs`):** D genes, D-free junctions,
  hypermutated V, reverse complements, reads that stop before their J.
- **Batch consistency:** chunk-batched parallel output equals per-read
  sequential output.
- **Input parsing:** multi-line FASTA/FASTQ, `@`/`+` at the start of quality
  lines, empty reads, CRLF, AIRR tables, batch resumption, and every error
  case.
- **Karlin–Altschul:** the series reproduces K = 1/3 for the ±1 walk and
  NCBI's ungapped blastn table (reward 1 / penalty −3 → λ 1.374, K 0.711, H 1.31).
- **Germline annotation and read regions:** see *Reliability* under
  *Germline annotation and regions*.

## Calling profiles

The default is `r`. Swig's interface names the profiles `r_optimized`,
`truth_optimized`, `igblast_compatible` and `sensitive_d`; both spellings
are accepted.

| Profile | V tuple | D tuple | J tuple | `min_d_match` | `P_D` (human IGH) |
|---|---|---|---|---|---|
| `r` (default) | (2,−3,−9,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 4 | derived, even odds: 10.55 |
| `sensitive-d` | (2,−3,−9,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 4 | derived, odds 2: 9.45 |
| `truth` | (2,−3,−5,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 6 | 0 |
| `igblast` | (2,−3,−5,−1) | (2,−4,−11,−1) | (2,−4,−13,−1) | 5 | 0 |

Every profile gates V and J hits by E-value (accept 1, anchor 10⁻³) and
stitches the co-linear V prefix and J suffix back onto the local alignment
before ranking, and calls each read with its own model. A with-D path must
strictly beat the no-D path; ties go to fewer D segments, then to more
germline support. `--min-d-match`, `--d-mismatch`, `--j-mismatch` and
`--min-j-length` override a profile's D and J settings.

## Assumptions and deviations from the method note

Each is a named field you can change.

1. **D presence.** `P_D` is derived from the D database (`DPrior::Derived`)
   instead of the note's 12-relaxed-to-10 rule. A sweep of constant costs put
   the relaxed rule on the same precision/recall curve as a single constant.
2. **Tandem D.** Not part of the call: an AIRR row holds one D. Swig's
   double-D screen reports a supported second D in a sidecar. The DAG itself
   handles any `max_d`, each D paying its own `P_D`.
3. **Retrieval.** Every 9-mer votes, and 7-mers too while the leader holds
   fewer than 25 seeds. V keeps 3 genes and walks up to 16 while gene votes
   stay within two standard deviations of the leader's (`z·√(2·4·L)` votes
   for a leader with `L`; the note used a fixed 8 or 5 %); every allele of a
   kept gene is a candidate. J is exhaustive. D uses the lossless L-mer
   prefilter instead of the note's k = 5/3 seeds.
4. **Keep.** V and J trace until the top 3 alignments that pass the E-value
   gate (plus ties), J per locus of the read's anchors; D keeps up to 32
   admissible alignments.
5. **Windows.** D windows reclaim 12 bases into V and J (the note's rescue of
   24/10 is realised by this window and the DAG's exact cut choice). The
   conserved-anchor window is always added when anchors are present. J is
   searched from 12 bases before the earliest V end when some V is an
   anchor, and over the whole read when none is; D windows start from the
   best anchoring V, and from the read's start when no V anchors; only
   anchors choose the locus of the J and D candidates. Without a J (or V)
   the D window reaches 60 bases past the V (before the J), the D and its N.
   When the anchors share no locus (a read holding two chains, as an scFv
   does), each of their loci with D genes bounds a D window by its own V
   and J.
6. **Significance.** V and J hits need E ≤ 1, and a call needs a V or J hit
   at E ≤ 10⁻³ (`EvalueGates`). The note has no reporting threshold.
7. **Tie-breaks.** Equal-score paths go to fewer D, then to more germline
   support; co-optimal alleles are listed with the best-supported first.
8. **Removed note mechanisms.**
   - The high-SHM junction "tie" (a zero-score block written onto V/J).
   - The ungapped J path bonus, the V-identity-switched J acceptance gates and
     the J hotspot profile.
   - The D short-run / long-span surcharges.
   - Inert by sensitivity analysis, so neutral by default (the fields
     remain): the V/J acceptance floors (24/10 columns, 60 % identity, which
     is the break-even identity of `(2, −3)` and so implied by a positive
     local score), the V…J pair bonus (24), the distributed-D gate, and the
     24/10-base retention floors.

   Ablation on simulated reads showed no accuracy benefit for any of them
   (see the assessment).
9. **Ungapped pre-ranking.** Available, and off by default, because it is
   lossy. Set `ungapped: Some(UngappedParams::default())` to enable it.
10. **Gap cost constraint.** `GapCosts::new` requires `open ≤ extend ≤ 0`.

## Known limitations

- **One alignment per D allele per window.** A D's optimal local alignment in
  its window is not always its optimal *clipped* alignment at the chosen cut.
  The 12-base reclaim keeps this rare; see the assessment.
- **Posteriors are conditional on the candidates.** An allele that retrieval
  or the keep policy dropped gets no mass. The call is a mode of the
  posterior over paths and cuts, not the allele of largest marginal, so a
  reported posterior can still fall below 0.5. D-presence calibration is the
  weakest (ECE up to 1.4 % on human IGH, 3.3 % on the synthetic database):
  one D cost per database cannot follow what the residual tracks, junction
  length (see ALGORITHM.md §8).
- **E-values** use exact *ungapped* λ and K and the total database length.
  Gapped scores are slightly higher (anti-conservative); near-identical alleles
  inflate the database length (conservative). They gate hits; they are not
  calibrated p-values.
- **Segment type is not checked from FASTA.** FASTA carries no segment type, so
  the command line cannot detect a D file passed as `--j`.
- **Compressed input.** Gzip must be piped through `zcat`; there is no built-in
  decompression.
- **FR3 starts need some annotated alleles.** The CDR2/FR3 boundary has no
  conserved residue. Across 13 tcrdist3 loci, FR3 starts not given were
  scored against IMGT:
  - with no given FR3 start, the chain default is wrong for 539 of 1,695
    alleles (32 %), all flagged but IGHV4-31\*05, whose FR3 lacks a codon;
  - with `fr3=` on the \*01 alleles, 1 of 872 others is wrong, flagged;
  - with `fr3=` on 10 %, 73 of 1,519 are wrong, all flagged, most of them
    behind a low-confidence anchor.
  Supply an IMGT-gapped V set, `SWIGMETA`, or `fr3=` on a few alleles per
  gene family.
- **Loci come from names.** Chain type is voted per locus, and a locus is
  read from allele names or headers. A database whose names carry no locus
  code is annotated and assigned as one locus, so a mixed heavy + light
  database with such names needs `--chain` or separate runs.
- **FR1 start.** FR1 starts at the first complete codon. A V sequence that
  includes a leader needs `fr1=` in the header.
- **A lone decoy Cys23.** When the Cys23 window holds exactly one Cys
  codon, it is taken as Cys23. An allele whose Cys23 is mutated or cut
  away, with a mutation-made Cys in the window, would get a wrong CDR1
  start without a flag. None occurred in the tests or the 13 loci.
- **IUPAC codes.** Germline IUPAC codes are stored as N, and a written
  database spells them N.
- **The constant region.**
  - The first `c_call` name is arbitrary among genes that share a short C
    part (IGHG1/2/4 open alike); the listed genes are right.
  - Bases after the C part (adapter, primer) can extend some alleles by a
    chance match and narrow the list: 2.2 % of simulated C parts from 16 nt
    followed by other bases.
  - A C part needs about 9 exact bases to pass the gate (10 human IGHC
    genes), and net indel drift beyond 3 bases inside C truncates its
    alignment.
  - 3′-partial J records (IMGT IGHJ6\*02/\*03) lose J calls to their longer
    siblings on reads that continue into C, whose first base completes them
    (−0.5 pp J on hard simulated reads, with or without `--c`).
- **The double-D screen** is evidence beside the call, not a VDDJ call:
  exact seeds lose sensitivity under mutation, extension is ungapped, and the
  pseudo-tandem veto is a bounded alternative-explanation test, not a
  probability over recombination histories.

## Roadmap

The earlier roadmap is done: calls from the per-read model, a D prior from
the D search's statistics (as `--d-presence-rate`), the AVX-512BW, NEON and
portable backends of one kernel, allele rescoring along the gene frame, the
annotation priors, C-segment assignment and the WebAssembly build
(ASSESSMENT §9.7–§11). Open, in order of expected value; reasons and
measurements are in [docs/ASSESSMENT.md](docs/ASSESSMENT.md):

1. **A prior on N/P lengths**, which the remaining D-presence
   miscalibration follows (ASSESSMENT §9.9).
2. **8-bit lanes for D jobs** (32 per AVX2 register).
3. **The open review findings**: [docs/REVIEW.md](docs/REVIEW.md) §3 ranks
   them and §4 gives the design recommendations.
