# vdj_assign

Seeded, banded, affine V(D)J germline assignment with placement-DAG
adjudication. It started as a Rust implementation of the SwiftIG / Swig 0.38.x
method note; its decision rules are now derived from one log-likelihood-ratio
objective (see [docs/ALGORITHM.md](docs/ALGORITHM.md)) instead of tuned patches.
Every component is a replaceable, statically dispatched trait, and every
search shortcut is either lossless or measured on simulated repertoires.

Zero dependencies, Rust ≥ 1.89, edition 2024. The alignment kernel runs on
AVX-512BW, AVX2, NEON, WebAssembly SIMD or portable lanes, chosen at run
time; output is byte-identical across them and across x86-64, AArch64 and
WebAssembly, so the same engine also runs in a browser (see *In the
browser* below).

- [docs/ALGORITHM.md](docs/ALGORITHM.md): the model, the search and its
  guarantees, structured for a methods section.
- [docs/ASSESSMENT.md](docs/ASSESSMENT.md): the critical assessment behind
  the current design. It covers what was ad hoc, what replaced it, the measurements, and
  the remaining plan.

## Build and run

```sh
cargo build --release                          # always run the release build
cargo test --release                           # unit, annotation, robustness, synthetic and simulated-accuracy tests
cargo test --release -- --ignored              # + accuracy at scale: 20,000 reads per preset

target/release/vdj-assign --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta \
    [--profile r|truth|igblast|sensitive-d] [--threads N] [--tandem-d] \
    [--chain heavy|light] [--forward-only] [--no-anchor-window] [--no-posterior] [--residual] \
    [--calls per-read|fixed] [--d-allele map|max-marginal] [--d-presence-rate R] \
    [--simd avx512|avx2|neon|portable|scalar] [--c IGHC.fasta [--c-evalue E]] \
    reads.fastq > calls.tsv

zcat reads.fastq.gz | target/release/vdj-assign --v V.fasta --d D.fasta --j J.fasta - > calls.tsv

# optional: write the annotated database for review or curation
target/release/vdj-assign prepare --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta [--c IGHC.fasta] \
    --out prepared/

# simulate reads with known truth (IgSim model), then score any tool's AIRR TSV against it
target/release/vdj-assign simulate --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta [--c IGHC.fasta] \
    --preset hard --reads 10000 --seed 1 --out sim          # sim.fasta + sim.truth.tsv
target/release/vdj-assign --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta sim.fasta > calls.tsv
target/release/vdj-assign evaluate --truth sim.truth.tsv calls.tsv   # works for IgBLAST etc. too
# two callers (or settings) on the same reads: paired differences, 95 % CIs, exact McNemar tests
target/release/vdj-assign evaluate --truth sim.truth.tsv calls.tsv other.tsv

# accuracy and throughput on simulated reads (synthetic database unless --v/--d/--j are given)
cargo run --release --example benchmark -- [--v V.fasta --d D.fasta --j J.fasta] \
    [--preset easy|mid|hard] [--reads N] [--threads N] [--random N]
# search errors vs model errors against an exhaustive oracle
cargo run --release --example diagnose -- --v V.fasta --d D.fasta --j J.fasta --preset hard
# where D errors arise: adjudication, not found, outside the D windows, eaten by J or V
cargo run --release --example d_errors -- --v V.fasta --d D.fasta --j J.fasta --preset hard
```

Reads may be FASTA or FASTQ. The format is detected from the first
character; multi-line FASTQ is supported and qualities are validated, then
ignored. FASTQ sequence lines must hold letters only (FASTA drops anything
else, but a FASTQ quality counts each base). Reads are streamed in batches of 16,384, so memory stays bounded on
large files, and `-` reads stdin. Gzip input is recognised and rejected with
a `zcat` hint, since the crate has no dependencies.

The CLI rejects:
- a repeated option, such as `--v` given twice;
- a FASTQ file passed as a germline database.

A debug build prints a warning, because it is 10–50× slower.

Germline databases are FASTA; IMGT reference headers
(`accession|allele|species|…`) are read as the allele name. No IMGT
annotation is needed: anchors and V
region boundaries are taken from the headers where present and detected
otherwise, on every start, in milliseconds (see *Germline annotation* below).
Identical sequences are merged and reported under all their names. Anchors
also enable the anchor-bounded D window. Output is an AIRR Rearrangement TSV.
Coordinates are 1-based closed intervals on the oriented query. CIGARs use
`S`/`N` for unaligned query and germline bases and `=`/`X`/`I`/`D` for
columns. `np1`–`np3` hold the junctional N/P regions, and `*_support` holds
Karlin–Altschul E-values. `*_posterior` is the posterior probability that
the call is germline-equivalent to the read's allele, and
`d_presence_posterior` the probability that the read has a D at all (see
*Posterior probabilities* below; `--no-posterior` skips them, ~30 % faster).
Region columns follow: `fwr1`–`fwr4`, `cdr1`–`cdr3`
(sequence, `_start`, `_end`), `cdr3_aa`, `junction`, `junction_length`,
`junction_aa`, `vj_in_frame`, `stop_codon` and `productive`. A read without a
significant V or J (see *Significance* below) is written with empty calls.

**Constant region (isotype).** With `--c C.fasta`, twelve `c_` columns
follow `productive` (`c_call`, `c_score`, `c_identity`, `c_support`,
`c_posterior`, `c_cigar`, sequence and germline coordinates and alignments);
without it the output is unchanged. C's germline base 0 is anchored within 3
bases of the called J's end and aligned exactly from there to the read's
end; V, D and J calls never change. `c_support` is the anchored E-value
`N_C·7·e^(−λS)` over the file's C genes and start diagonals, and C is called
at `c_support ≤ 10⁻³` (`--c-evalue`). Give one locus' spliced coding C
regions, CH1 onward (IMGT or codon-aligned records such as abstar's both
work), since every gene in the file counts. `c_call` lists every co-optimal
allele: for a short C part shared by several genes (IGHG1/2/4 open alike)
the *first* name is any of them, so read the isotype from the genes listed.

As a library:

```rust
use vdj_assign::prelude::*;

let database = prepare(v_set, Some(d_set), j_set, AnnotateOptions::default())?; // dedup + anchors
let config = CallingProfile::ROptimized.config();
let assigner = DefaultAssigner::from_config(&config, database.v, database.d, database.j)?;
let mut workspace = assigner.workspace();                  // one per thread
let call = assigner.assign(&encode(b"CAGGTG..."), &mut workspace);
let calls = assign_all(&assigner, &queries, 0);            // parallel, input order
```

or a whole run as the command line performs it, from FASTA text to AIRR TSV
(`vdj_assign::run`; options by their command-line names):

```rust
let mut options = RunOptions::default();
options.set("tandem-d", None)?;
let texts = GermlineTexts { v: &v_fasta, d: Some(&d_fasta), j: &j_fasta, c: None };
let database = PreparedDatabase::from_fasta(texts, options.annotate())?;
let tally = Run::new(&options, database)?.write_airr(reads.as_bytes(), &mut out)?;
```

## In the browser (WebAssembly)

`wasm/` compiles the engine to a WebAssembly module, and `web/` holds its
JavaScript API, a Web Worker and a page that assigns reads in the browser.
The page reads the files and runs the module locally; nothing is uploaded.

```sh
rustup target add wasm32-unknown-unknown
wasm/build.sh                                 # web/vdj_assign_simd.wasm, web/vdj_assign.wasm
python3 -m http.server -d web 8000            # then open http://localhost:8000
node web/test.mjs target/release/vdj-assign   # the modules against the command line
```

```js
import { load, assignInWorkers } from './vdj-assign.js';

const vdj = await load();                     // the SIMD build where the engine has SIMD
const run = vdj.createRun({ v, d, j, c, options: { profile: 'r', 'tandem-d': true } });
const tsv = run.assign(fastqText);            // AIRR TSV, header included
run.free();

// In a browser, one run on Web Workers, rows in input order:
const result = await assignInWorkers({ germlines: { v, d, j, options }, reads, workers: 4 });
```

- **Same engine, same output.** The module runs `run::Run`, as the command
  line does: the same options by the same names (`profile`, `tandem-d`,
  `c-evalue`, …), the same database preparation and warnings, the same
  AIRR rows. On 45,000 rows (three presets, three option sets) and in
  headless Chromium with four workers, its TSV equals the command line's
  byte for byte; CI checks it on every push (`web/test.mjs`).
- **Two modules.** WebAssembly has no run-time feature detection, so
  `vdj_assign_simd.wasm` uses the `simd128` kernel (every current browser,
  Node.js) and `vdj_assign.wasm` the scalar kernel; `load()` validates a
  SIMD probe and fetches the right one. 0.9 MB each, 0.3 MB gzipped, no
  imports.
- **Speed.** With SIMD under V8, 1,100–1,450 reads/s per worker on
  simulated reads, about 60 % of the native AVX2 kernel (1,850–1,970 on the
  same hard-preset file, one thread). `assignInWorkers` gives each worker
  every `n`-th block of 256 reads; workers share nothing, and each holds a
  copy of the reads.
- **Zero dependencies there too.** The module exposes a small C ABI over
  linear memory (`wasm/src/lib.rs`), so it builds with `cargo` and the
  `wasm32-unknown-unknown` target alone; `web/vdj-assign.js` wraps it, and
  runs in browsers, Web Workers and Node.js 18 or later.
- **The page.** Germline FASTA files (or demo data: a synthetic database
  and simulated reads), reads as FASTA or FASTQ (gzip decompressed by the
  browser), every assignment option and a worker count; it gives the AIRR
  TSV, the prepared database files and a preview. CI uploads `web/` with
  both modules as the `vdj-assign-web` artifact, ready for static hosting.

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

### Header keys

All positions are 0-based in the germline sequence; unknown keys are ignored.

| Key | Meaning |
|---|---|
| `anchor=` | first base of V Cys104 / J Trp-Phe118 |
| `fr1=` `cdr1=` `fr2=` `cdr2=` `fr3=` | V region starts |
| `aliases=A,B` | other names of this sequence |

### IMGT-gapped input

An IMGT-gapped V set (the IMGT reference format: `.` for the gap columns) is
numbered from its gap columns. The FR1, CDR1, FR2, CDR2 and FR3 starts are
IMGT columns 1, 79, 115, 166 and 196, and Cys104 is column 310.

- The gaps are used only when most gapped alleles are codon-aligned and
  carry Cys, Trp and Cys at IMGT 23, 41 and 104. Otherwise they are ignored
  (a codon alignment in other columns) and the summary line says why.
- A region that the numbering places before a 5′-partial sequence stays
  empty; a motif hit cannot fill it.
- An allele whose numbering ends before Cys104 gets no anchor (`3′-partial`).
- A one- or two-column interior gap (an IMGT frameshift marker) is noted, and
  the columns after it still map exactly.
- Header keys win over the gap columns.
- Ungapped alleles appended to a gapped set are annotated as below, and take
  their FR3 start from their gapped relatives.

### Sources, in order of trust

A header value is never overridden.

1. **Header keys**, then the **IMGT gap columns** (source `imgt-gaps`).
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
   learned once per database, with Jeffreys (½) pseudo-counts:
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
  (`prepare` and `assign` also warn);
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
clip index. This costs under a microsecond per read, against ~450 µs for assignment and
11.8 µs for the whole AIRR record, and nothing is added to assignment itself.

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

`vdj-assign prepare` writes `V.fasta`, `D.fasta` and `J.fasta` with every
anchor and region start in the headers (plus `annotation=SOURCE`), `C.fasta`
with `--c`, and `annotation.tsv` with one row per allele:

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
stick. Assignment runs log a one-line summary. They warn when any allele is
transferred, low-confidence or missing, and when a TR or IGL database gives
no FR3 start, so that every FR3 start rests on one chain default.

## Architecture

```text
 chunk of reads ─► per read: strand from germline k-mer support; V candidates (gene frames)
                ─► V stage, whole chunk: score batch, sibling bounds, traceback waves
                ─► per read: J window right of V ─► J stage, whole chunk (every allele)
                ─► per read: D windows from V/J geometry, D candidates (L-mer bitsets)
                ─► D stage, whole chunk ─► per read: admission gates
                ─► per read: placement DAG over VJ / VDJ / VDDJ ─► anchor check ─► Call
                ─► with --c, whole chunk: C anchored at the J exit (one job per reach class) ─► AIRR
```

Every per-read decision is taken per read. Batch aligners must answer
each job exactly as a single call would, so a read's result never depends on
its chunk: `assign(q)` is literally a chunk of one, and this is tested. The
stages are the methods `stage_orient`, `stage_v`, `stage_j`, `stage_d` and
`stage_adjudicate` of `Assigner`.

| Module | Responsibility |
|---|---|
| `alphabet`, `fasta`, `germline`, `reads` | 2-bit codes (N = 4), germline FASTA and sets with anchors, streaming FASTA/FASTQ |
| `translate`, `annotate`, `prepare`, `read_regions` | Codon table, anchors and V regions (header > IMGT gaps > motif > consensus > relatives > database models), per-read FR/CDR/junction projection |
| `scoring`, `stats` | Log-odds substitution profiles (uniform or position-specific), affine gaps, Karlin–Altschul λ/K/E-values |
| `orient` | Strand from germline k-mer presence on both strands, with a significance test |
| `seed`, `index`, `retrieve/frames` | Rolling k-mers, CSR k-mer index over gene frames (alleles of a gene in one coordinate system) |
| `retrieve/seed_vote`, `retrieve/budget` | Radix-sorted diagonal votes per frame, gene budget, candidate bands |
| `retrieve/exhaustive`, `retrieve/lmer_mask`, `retrieve/ungapped` | J: every allele; D: lossless L-mer bitsets; optional ungapped pre-ranking |
| `band`, `align/gotoh`, `align/inter_seq`, `align/lanes` | Diagonal bands; scalar banded Gotoh; inter-sequence Gotoh written once over a SIMD lane trait (AVX-512BW 32 lanes; AVX2, NEON, portable 16) |
| `segment`, `family_bound` | One search stage: score → sibling bounds along the gene frame → traceback waves → acceptance and E-value gates → outer extension |
| `clip`, `evidence`, `placement` | O(1) clip-and-rescore at any cut; exact-run / identity evidence |
| `window`, `admission` | Where J and D are searched; D evidence gates |
| `constant`, `retrieve/reach` | Constant region after adjudication: splice window, anchor-column DP, anchored E-value and posterior; one candidate per class of alleles the read cannot reach apart |
| `adjudicate/{lattice,dag,posterior,units,model,decode,residual}` | The placement lattice; Viterbi and forward–backward over it; the per-read model in integer units both read; grammar and log-odds path model; max-marginal D decoding; residual baseline |
| `flank` | Outer-flank extension of V and J alignments |
| `pipeline`, `config`, `profiles` | Stage-wise composition, declarative configuration, calling profiles |
| `run` | A whole run, shared by the command line and the WebAssembly module: options by their command-line names, the assigner they build, reads streamed to AIRR TSV |
| `report`, `batch` | AIRR output, parallel batches |
| `sim/{random,params,simulate,synthetic,evaluate}` | IgSim port, synthetic germlines, accuracy evaluation |

## The model in brief

The full description, with the guarantee of each search step, is in
[docs/ALGORITHM.md](docs/ALGORITHM.md).

- **Scores are log-odds.** `(2, −3)` substitution scores are the
  log-likelihood ratio, in units of 1/λ nats with λ = 0.634, of "templated by
  this germline at 11.2 % divergence" against "random sequence". N and P bases
  score 0 under both hypotheses. Affine gaps price indels.
- **The path objective** is `Σ clipped segment scores − P_D·#D − tandem·[2nd D]`,
  maximised exactly over the topologies VJ / VDJ / VDDJ,
  over the placements, and over both cuts between neighbours (the earlier one
  exits at `x`, the later enters at `y ≥ x`, and the bases between are N/P).
  Ties go to fewer D, then to more germline support (the summed unclipped
  scores).
- **Calls come from the read's own model.** Each path is valued at the read's
  divergence `p̂`, estimated before adjudication from the V and J columns no
  path contests: a match is worth `ln 4(1−p̂)` nats, a mismatch `ln(4p̂/3)`, a
  D costs `P_D/2` matched D columns. Values are integers (1/2¹⁶ nat), so the
  maximum is exact, and the posterior reads the same integers: the call is a
  mode of the reported distribution. `--calls fixed` keeps the fixed scores,
  of which the per-read model at 11.2 % divergence is an exact multiple.
- **`P_D` is derived, not tuned.** A D is one of `N` genes at one of `L`
  junction positions. Including it pays off when `λ·S_D > ln(N·L) − ln(prior
  odds)`, which gives `P*_D` = 10.55 for human IGH (32 genes, L = 25, even
  odds) and 6.17 for TRB, kept unrounded; on integer scores the fixed scale
  charges `⌊P*_D⌋`, which decides `S > P*_D` exactly. `--d-presence-rate π`
  counts what the D search tries instead (Karlin–Altschul `K·m·n_D` over
  alleles × offsets × positions) at the D-presence rate `π`; the rate has to
  come from the user, so it is not the default.
- **Significance.** A V or J hit is kept at E ≤ 1 against its database and
  window. A read is called only when its path contains an *anchor*, a V or J
  hit at E ≤ 10⁻³, so random sequence receives no call.
- **Constant region.** C is spliced to J at a fixed phase, so C's base 0
  lies within `s = 3` bases of the J's exit. A bonus column in front of C
  turns every local aligner into an exact anchored one, and C is called at
  an anchored E-value `N_C·(2s + 1)·e^(−λS) ≤ 10⁻³`, after adjudication, so
  it never moves a V, D or J call.
- **Posterior probabilities.** Read as a distribution, `P(path, cuts) ∝
  exp(nats)`, the same model gives calibrated probabilities by
  forward–backward over the same DAG; the other priors stay at their
  derivation λ. `--d-allele max-marginal` reports the D allele of largest
  posterior marginal instead of the call's (off by default: it gains only on
  heavily mutated reads).

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
| `Adjudicator` | `PlacementDag<VdjGrammar, EngineeringObjective>` | `ResidualAdjudicator` |
| `Grammar`, `PathModel` | `VdjGrammar::{VDJ, TANDEM_D, VJ}`, `EngineeringObjective` | your own path grammar / objective |
| `VdjAssigner` (+ `assign_many`) | chunk-staged `Assigner` | default `assign_many` loops `assign` |

```rust
let experiment = DefaultAssigner::from_config(&config, v, Some(d), j)?
    .with_adjudicator(ResidualAdjudicator::default())
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
pure-speed steps; the steps
that changed calls are listed in [docs/ASSESSMENT.md](docs/ASSESSMENT.md)
with their accuracy effect.

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

Sampled seeds (`SeedSampling`) remain available and correct: closed syncmers
(Edgar 2021; the smallest `s`-mer of the seed's key, `s = k − w`, first or
last) are context-free, so a seed shared by read and germline is kept in both
or neither, and have a *window guarantee*: any `w` consecutive k-mers hold
one on a contiguous shape, so every unmutated stretch of `k + w − 1` bases
(12 for 9-mers, `w = 4`) yields a shared seed (0 seedless 12-base windows on
IGHV). Open syncmers (the smallest `s`-mer in the middle, Shaw & Yu 2022) are
context-free without that guarantee (5 % of 12-base IGHV windows seedless at
`w = 2`); minimisers have the guarantee but are not context-free. Under the
noise-based margin every k-mer beats closed and open syncmers (minimisers
were not re-measured): closed syncmers (`w = 4`, 40 % of k-mers) are 1.1
points of mean V gene behind (4.4 at 30 % SHM) for 0.4 % fewer instructions
on the hard preset (and slightly more at 30 % SHM). Their earlier parity came
from the fixed margin and was fragile: a different hash order alone moved
them 0.3 points on the hard preset.

**Strand from a k-mer bitset.** A presence bitset over all 4¹¹ k-mers of the V
and J germline (512 KiB). One pass computes forward and reverse-complement
keys together, so each strand costs one bit test per base. Each read position
counts at most once, so repeats cannot inflate the evidence. The count is
trusted only when it exceeds its Poisson null (`λ₀ + 3√λ₀`); otherwise the
read falls back to seed votes on both strands. Only the chosen strand is
searched for V.

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

**Exhaustive J, small window.** The J set is tiny (10–20 alleles) and under
heavy SHM keeps few intact seeds, so every J allele is aligned, unbanded, in
the read right of V (from 12 bases before the earliest V end).

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
take a slot. With at most `top + 1` candidates the scoring pass is skipped.

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
  on a target without baseline SIMD (WebAssembly without `simd128`), where
  emulated lanes run about 4× slower than one scalar job, it chooses the
  scalar kernel.

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
therefore exact for the call. A D below it still carries posterior mass, so
the floor applies only with `--no-posterior`; posteriors keep every
admissible D (+1 % of instructions).

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
self-scheduled in chunks of 64 from an atomic counter over scoped threads.
Results come back in input order, and reporting is formatted on the workers
(`assign_map`).

## Simulation and evaluation

`vdj_assign::sim` is a dependency-free port of
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
AIRR Rearrangement TSV row. `vdj-assign simulate` and `vdj-assign evaluate`
expose both, so other callers can be scored on identical reads.
`sim::Evaluation` scores calls against the truth:
- an allele call is *equivalent* when the primary allele contains the
  trimmed, pre-mutation germline body (no read can separate such alleles);
- *gene*, *listed* (truth among the reported names) and *exact*;
- boundary errors, D presence by remnant length, junction exactness;
- a random-read negative control;
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
t statistic. On 20,000 reads per preset, for example, the
residual adjudicator gains +0.4–0.6 pp of D accuracy (p < 10⁻⁴) only by
giving 21–25 pp more D-free reads a D; see
[docs/ASSESSMENT.md §9.5](docs/ASSESSMENT.md) for the full table.

Human IGH (541 V / 90 D / 16 J; a local copy, not shipped), 5000 reads per
preset, germline-equivalent accuracy in %, before → after the redesign:

| Preset (IgSim) | Orientation | V | V gene | D | J | Junction | Random reads called |
|---|---|---|---|---|---|---|---|
| easy (no SHM) | 99.94 → **100** | 99.54 → **99.84** | 99.92 → **99.98** | 98.78 → **99.29** | 99.90 → **99.94** | 99.94 → **100** | 30.5 → **0** |
| mid (IgM-like) | 100 → 100 | 97.76 → **97.94** | 99.84 → 99.84 | 91.96 → **93.00** | 98.84 → 98.80 | 99.93 → **99.97** | 30.2 → **0** |
| hard (IgG-like) | 99.98 → 99.98 | 89.74 → **89.96** | 98.44 → **98.50** | 81.72 → **83.58** | 92.70 → **93.32** | 98.39 → **99.55** | 31.2 → **0** |

Junctions "after" use the corrected anchors of 3′-partial V alleles, which
leave those alleles' reads without a templated junction (hard: 2416 → 2214
scored reads); with the old annotation the same calls scored 98.47 %.

The synthetic database (`SyntheticSpec::default()`, seed 1) gives, after:
easy V/D/J 99.98/99.82/99.98, mid 99.54/87.90/98.92, hard 97.18/79.69/94.06,
and junctions ≥ 99.8 %.

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
strand, 5000 per preset, scored by `vdj-assign evaluate`):

| Preset | Tool | V | V gene | D | J | Junction | D-free reads given a D | Throughput, 4 threads |
|---|---|---|---|---|---|---|---|---|
| easy | vdj_assign | **99.88** | **100.00** | **99.29** | 99.84 | **100.00** | 0 % | **5200 reads/s** |
| easy | IgBLAST | 99.12 | 99.78 | 98.86 | **99.86** | 99.77 | 0 % | 67 reads/s |
| mid | vdj_assign | **97.88** | **99.72** | 93.35 | **98.74** | **99.97** | **23 %** | **5100 reads/s** |
| mid | IgBLAST | 97.18 | 99.62 | **93.73** | 98.60 | 99.72 | 43 % | 66 reads/s |
| hard | vdj_assign | **89.44** | 98.84 | **83.00** | **92.80** | **99.69** | **16 %** | **5200 reads/s** |
| hard | IgBLAST | 89.16 | **98.90** | 81.75 | 91.84 | 98.07 | 30 % | 71 reads/s |

IgBLAST's higher mid-preset D accuracy comes with twice the D calls on D-free
junctions. Throughput is wall clock for the whole command, posteriors on,
measured with the previous V retrieval defaults (the current ones use 4 % more
instructions).
IgBLAST takes J anchors from its auxiliary file, which lacks three of this
database's J alleles; the comparison adds them (frame and CDR3 end from
`vdj-assign prepare`, which agrees with the file on all 13 shared alleles).
Without them IgBLAST reports no junction for 18 % of reads.

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
  ones; anchored
  truth, C-length strata and the simulator's reading frame are tested, and
  `tests/simulated.rs` bounds C accuracy, false calls and `c_posterior`
  calibration.
- **Cross-architecture:** the suite passes on x86-64 and on AArch64 (CI runs
  both), and AArch64 output equals x86-64 output byte for byte. It also
  passes on `wasm32-wasip1` with `simd128` under Wasmtime (CI), where the
  kernel tests compare the `simd128`, portable and scalar kernels, and the
  browser modules' output equals the native command line's
  (`web/test.mjs`).
- **Synthetic calling (`tests/synthetic.rs`):** D genes, D-free junctions,
  tandem VDDJ, hypermutated V, reverse complements, residual adjudicator.
- **Batch consistency:** chunk-batched parallel output equals per-read
  sequential output.
- **Input parsing:** multi-line FASTA/FASTQ, `@`/`+` at the start of quality
  lines, empty reads, CRLF, batch resumption, and every error case.
- **Karlin–Altschul:** the series reproduces K = 1/3 for the ±1 walk and
  NCBI's ungapped blastn table (reward 1 / penalty −3 → λ 1.374, K 0.711, H 1.31).
- **Germline annotation and read regions:** see *Reliability* under
  *Germline annotation and regions*.

## Calling profiles

The default is `r`. `r_optimized` is accepted as an alias.

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
germline support.

## Assumptions and deviations from the method note

Each is a named field you can change.

1. **D presence.** `P_D` is derived from the D database (`DPrior::Derived`)
   instead of the note's 12-relaxed-to-10 rule. A sweep of constant costs put
   the relaxed rule on the same precision/recall curve as a single constant.
2. **Tandem D.** A second D pays its own `P_D` plus 12 (`tandem_d_penalty`).
3. **Retrieval.** Every 9-mer votes, and 7-mers too while the leader holds
   fewer than 25 seeds. V keeps 3 genes and walks up to 16 while gene votes
   stay within two standard deviations of the leader's (`z·√(2·4·L)` votes
   for a leader with `L`; the note used a fixed 8 or 5 %); every allele of a
   kept gene is a candidate. J is exhaustive. D uses the lossless L-mer
   prefilter instead of the note's k = 5/3 seeds.
4. **Keep.** V and J trace until the top 3 alignments that pass the E-value
   gate (plus ties); D keeps up to 32 admissible alignments.
5. **Windows.** D windows reclaim 12 bases into V and J (the note's rescue of
   24/10 is realised by this window and the DAG's exact cut choice). The
   conserved-anchor window is always added when anchors are present. J is
   searched from 12 bases before the earliest V end.
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
  the CLI cannot detect a D file passed as `--j`.
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
  Supply an IMGT-gapped V set, or `fr3=` on a few alleles per gene family.
- **Chain type is one per database.** It is voted per database, so a mixed
  IGK + IGL database is fine, but a mixed heavy + light database needs
  `--chain` or separate runs.
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

## Roadmap

In order of expected value; reasons and measurements are in
[docs/ASSESSMENT.md](docs/ASSESSMENT.md).

1. ~~**Calls from the per-read model.**~~ Done: one integer model at the
   read's divergence drives the Viterbi call and forward–backward, so the
   call is a mode of its posterior. Paired on 20,000 reads per preset:
   4.0–4.3 pp more D-free reads left without a D on mid sets, D-presence
   Brier better on every set (t −5 to −24), D allele within noise for
   remnants of 4 nt or more (the overall −0.3–0.4 pp on the synthetic
   database lies in 1–3 nt remnants, below `min_d_match`). Max-marginal D
   alleles are a switch (+0.25 pp on hard reads only). See ASSESSMENT §9.8.
2. ~~**D prior from the D search's statistics.**~~ Done as a switch,
   `--d-presence-rate π`: `P*_D = [ln(K·m·n_D) − ln(π/(1 − π))]/λ` counts
   every distinct allele, offset and read position the search tries, and is
   read per read like the default. At π = 0.9 it meets every criterion on
   20,000 reads per preset (D-presence Brier better on all six sets, no D
   call +0.77 pp on human hard, D on remnants of 4 nt or more within
   noise), but that rate was found by measuring; at the simulator's own
   rate (0.94) exact 5-mers win and D-free specificity falls 14–26 pp. The
   default stays `ln(N·L)` at even odds. Next on that line: a prior on N/P
   lengths, which the remaining miscalibration follows. See ASSESSMENT §9.9.
3. ~~**Wider SIMD.**~~ Done: AVX-512BW, NEON and portable backends of one
   kernel. Next on that line: 8-bit lanes for D jobs (32 per AVX2 register).
4. ~~**Allele rescoring along the gene frame.**~~ Done: siblings are bounded
   from both sides along their representative's traceback, and DP runs only
   where the bounds straddle the cutoff (V and J); calls unchanged, 9 % fewer
   instructions on hard reads.
5. ~~**Annotation.**~~ Done: IMGT gap numbering, FR3 from annotated relatives
   or a learned length, a learned Cys23 model, and structured `uncertain=`
   flags (see the assessment, §10).
6. ~~**C-segment assignment.**~~ Done: `--c` anchors C's base 0 at the
   called J's exit and aligns it exactly with any kernel, after
   adjudication, so V/D/J output is byte-identical with and without it and
   costs nothing without it. On 20,000 simulated reads per preset (human IGH
   with abstar's IGHC; synthetic): ≤ 2 false calls on ~4,100 C-free reads,
   primary C allele anchored-equivalent ≥ 99.6 % from 16 nt, `c_posterior`
   ECE ≤ 0.18 %, 9 % of instructions with `--c`. See ASSESSMENT §11.
7. ~~**WebAssembly.**~~ Done: a `simd128` kernel backend, a zero-dependency
   module (`wasm/`) and a browser page (`web/`) whose output equals the
   command line's byte for byte, at about 60 % of native AVX2 speed per
   worker. Next on that line: shared-memory threads (one copy of the reads
   for all workers), which need cross-origin isolation from the host.
