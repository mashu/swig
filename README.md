# vdj_assign

Seeded, banded, affine V(D)J germline assignment with placement-DAG
adjudication. It started as a Rust implementation of the SwiftIG / Swig 0.38.x
method note; its decision rules are now derived from one log-likelihood-ratio
objective (see [docs/ALGORITHM.md](docs/ALGORITHM.md)) instead of tuned patches.
Every component is a replaceable, statically dispatched trait, and every
search shortcut is either lossless or measured on simulated repertoires.

Zero dependencies, Rust ≥ 1.75, edition 2021.

- [docs/ALGORITHM.md](docs/ALGORITHM.md): the model, the search and its
  guarantees, structured for a methods section.
- [docs/ASSESSMENT.md](docs/ASSESSMENT.md): the critical assessment behind
  the current design. It covers what was ad hoc, what replaced it, the measurements, and
  the remaining plan.

## Build and run

```sh
cargo build --release                          # always run the release build
cargo test                                     # 67 unit, 7 annotation, 10 synthetic, 7 simulated-accuracy tests

target/release/vdj-assign --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta \
    [--profile r|truth|igblast|sensitive-d] [--threads N] [--tandem-d] \
    [--chain heavy|light] [--forward-only] [--no-anchor-window] [--no-posterior] [--residual] \
    reads.fastq > calls.tsv

zcat reads.fastq.gz | target/release/vdj-assign --v V.fasta --d D.fasta --j J.fasta - > calls.tsv

# optional: write the annotated database for review or curation
target/release/vdj-assign prepare --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --out prepared/

# simulate reads with known truth (IgSim model), then score any tool's AIRR TSV against it
target/release/vdj-assign simulate --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta \
    --preset hard --reads 10000 --seed 1 --out sim          # sim.fasta + sim.truth.tsv
target/release/vdj-assign --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta sim.fasta > calls.tsv
target/release/vdj-assign evaluate --truth sim.truth.tsv calls.tsv   # works for IgBLAST etc. too

# accuracy and throughput on simulated reads (synthetic database unless --v/--d/--j are given)
cargo run --release --example benchmark -- [--v V.fasta --d D.fasta --j J.fasta] \
    [--preset easy|mid|hard] [--reads N] [--threads N] [--random N]
# search errors vs model errors against an exhaustive oracle
cargo run --release --example diagnose -- --v V.fasta --d D.fasta --j J.fasta --preset hard
```

Reads may be FASTA or FASTQ. The format is detected from the first
character; multi-line FASTQ is supported and qualities are validated, then
ignored. Reads are streamed in batches of 16,384, so memory stays bounded on
large files, and `-` reads stdin. Gzip input is recognised and rejected with
a `zcat` hint, since the crate has no dependencies.

The CLI rejects:
- a repeated option, such as `--v` given twice;
- a FASTQ file passed as a germline database.

A debug build prints a warning, because it is 10–50× slower.

Germline databases are FASTA. No IMGT reference is needed: anchors and V
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
*Posterior probabilities* below; `--no-posterior` skips them, ~15 % faster).
Region columns follow: `fwr1`–`fwr4`, `cdr1`–`cdr3`
(sequence, `_start`, `_end`), `cdr3_aa`, `junction`, `junction_length`,
`junction_aa`, `vj_in_frame`, `stop_codon` and `productive`. A read without a
significant V or J (see *Significance* below) is written with empty calls.

As a library:

```rust
use vdj_assign::prelude::*;

let config = CallingProfile::ROptimized.config();
let assigner = DefaultAssigner::from_config(&config, v_set, Some(d_set), j_set)?;
let mut workspace = assigner.workspace();                  // one per thread
let call = assigner.assign(&encode(b"CAGGTG..."), &mut workspace);
let calls = assign_all(&assigner, &queries, 0);            // parallel, input order
```

## Germline annotation and regions

V and J germline genes are in frame with their conserved anchors — V Cys104
and J Trp/Phe118 (IMGT numbering) — so an anchor fixes the reading frame, and
the V region boundaries follow from three conserved residues:

```text
FR1  = [frame start, Cys23 + 4 codons)      CDR1 = [Cys23 + 4, Trp41 − 2)
FR2  = [Trp41 − 2, Trp41 + 15)              CDR2 = [Trp41 + 15, FR3 start)
FR3  = [Cys104 − (L − 1), Cys104 + 1)       L = 38 codons heavy, 36 light
CDR3 = [Cys104 + 1, J Trp/Phe118)           FR4  = [J Trp/Phe118, J end)
junction = [Cys104, Trp/Phe118 + 1)         (both conserved codons included)
```

These are the IMGT/AIRR definitions. FR2 has no IMGT gaps, and FR3 has
constant gaps per locus (heavy: 73; kappa and lambda: 73, 81, 82).

### Header keys

All positions are 0-based in the germline sequence; unknown keys are ignored.

| Key | Meaning |
|---|---|
| `anchor=` | first base of V Cys104 / J Trp-Phe118 |
| `fr1=` `cdr1=` `fr2=` `cdr2=` `fr3=` | V region starts |
| `aliases=A,B` | other names of this sequence |

### Sources, in order of trust

A header value is never overridden.

1. **Header keys.**
2. **Motif search** on a one-pass three-frame translation. For V, Cys23, Trp41
   and Cys104 are placed jointly under spacing constraints, scored by
   conserved context: `[ST]C`; `W` with `Q`+3, `P`+5, `G`+6 (`WVRQAPG` /
   `WYQQKPG`); and `YYC`, with `Y` weighted above `F`/`H`. Upstream stop codons
   are penalised, which also selects the frame. For J, `[WF]G.G` with an open
   FR4.
3. **Database consensus.** Once pass 1 has run, the chain type is voted:
   J `W`/`F` and Trp41 context `W.R`/`W.Q`. That turns the Trp41 → Cys104 window
   into a CDR2-length constraint, `[14 + L, 26 + L]`. The dominant V 3′ tail
   (applied to confident alleles only when it holds ≥ 50 % of the database) and
   the Cys23 codon index of full-length alleles (only in their frame) become
   priors for a second pass.
4. **Sibling alignment.**
   - An ambiguous Cys23 in a 5′-truncated allele is taken from the most similar
     unambiguous allele.
   - An allele whose motifs stay weak (pseudogene, mutated Cys) gets its anchor
     and every region start transferred by alignment from the most similar
     trusted allele.
   - A relative must share ≥ 50 % of the target's distinct 8-mers, so a chance
     match can never donate.
5. **Low-confidence motif hit, or missing.** Both are reported.

### Speed

- **Codon table:** translation is one lookup per codon in a 64-entry table on
  the 2-bit codon index.
- **Candidate scoring:** a V allele costs a few dozen candidate scores per pass.
- **Relative search:** an inverted 8-mer index (CSR over 65,536 keys, built
  lazily) is searched only for the rare unresolved allele, touching only
  alleles that share k-mers.

| Database | Preparation time |
|---|---|
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
| 12,000 alleles over 40 databases | 0 anchor, FR2, CDR2 or FR3 errors |
| CDR1 starts in the same 12,000 | 5 wrong, all flagged `Cys23 ambiguous` (truncated alleles without a sibling); 0 silent errors |
| Pseudogene with Cys104 and Trp41 destroyed | anchor and all regions transferred exactly from its sibling |
| 200 simulated reads, half reverse-complemented | FR1–FR4, CDR1–CDR3, junction, frame, stop codon and productive flag exact |

### `prepare`

`vdj-assign prepare` writes `V.fasta`, `D.fasta` and `J.fasta` with every
anchor and region start in the headers (plus `annotation=SOURCE`), and
`annotation.tsv` with one row per allele:

- frame, anchor and its residue;
- region starts, CDR1 and CDR2 lengths in codons;
- source (`header`, `motif`, `consensus`, `transfer:NAME`, `low-confidence`,
  `missing`) and notes.

Reloading the written FASTA uses the headers as-is, so manual corrections
stick. Assignment runs log a one-line summary and warn when any allele is
transferred, low-confidence or missing.

## Architecture

```text
 chunk of reads ─► per read: strand from germline k-mer support; V candidates (gene frames)
                ─► V stage, whole chunk: score batch, sibling bounds, traceback waves
                ─► per read: J window right of V ─► J stage, whole chunk (every allele)
                ─► per read: D windows from V/J geometry, D candidates (L-mer bitsets)
                ─► D stage, whole chunk ─► per read: admission gates
                ─► per read: placement DAG over VJ / VDJ / VDDJ ─► anchor check ─► Call ─► AIRR
```

Every per-read decision is taken per read. Batch aligners must answer
each job exactly as a single call would, so a read's result never depends on
its chunk: `assign(q)` is literally a chunk of one, and this is tested. The
stages are the methods `stage_orient`, `stage_v`, `stage_j`, `stage_d` and
`stage_adjudicate` of `Assigner`.

| Module | Responsibility |
|---|---|
| `alphabet`, `fasta`, `germline`, `reads` | 2-bit codes (N = 4), germline FASTA and sets with anchors, streaming FASTA/FASTQ |
| `translate`, `annotate`, `prepare`, `read_regions` | Codon table, anchors and V regions (header > motif > consensus > sibling transfer), per-read FR/CDR/junction projection |
| `scoring`, `stats` | Log-odds substitution profiles (uniform or position-specific), affine gaps, Karlin–Altschul λ/K/E-values |
| `orient` | Strand from germline k-mer presence on both strands, with a significance test |
| `seed`, `index`, `retrieve/frames` | Rolling k-mers, CSR k-mer index over gene frames (alleles of a gene in one coordinate system) |
| `retrieve/seed_vote`, `retrieve/budget` | Radix-sorted diagonal votes per frame, gene budget, candidate bands |
| `retrieve/exhaustive`, `retrieve/lmer_mask`, `retrieve/ungapped` | J: every allele; D: lossless L-mer bitsets; optional ungapped pre-ranking |
| `band`, `align/gotoh`, `align/inter_seq` | Diagonal bands; scalar banded Gotoh; AVX2 inter-sequence Gotoh (16 jobs per register) |
| `segment`, `family_bound` | One search stage: score → sibling upper bounds → traceback waves → acceptance and E-value gates → outer extension |
| `clip`, `evidence`, `placement` | O(1) clip-and-rescore at any cut; exact-run / identity evidence |
| `window`, `admission` | Where J and D are searched; D evidence gates |
| `adjudicate/{dag,model,residual}` | Semi-Markov Viterbi over placements; grammar and log-odds path model; residual baseline |
| `refine` | Outer-flank extension; `Refiner` extension point |
| `pipeline`, `config`, `profiles` | Stage-wise composition, declarative configuration, calling profiles |
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
- **`P_D` is derived, not tuned.** A D is one of `N` genes at one of `L`
  junction positions. Including it pays off when `λ·S_D > ln(N·L) − ln(prior
  odds)`, which gives 11 for human IGH (32 genes, L = 25, even odds) and 6
  for TRB.
- **Significance.** A V or J hit is kept at E ≤ 1 against its database and
  window. A read is called only when its path contains an *anchor*, a V or J
  hit at E ≤ 10⁻³, so random sequence receives no call.
- **Posterior probabilities.** Read as a distribution, `P(path, cuts) ∝
  exp(nats)`, the objective gives calibrated probabilities by forward–backward
  over the same DAG. Alignment evidence is read as log-odds at the read's own
  divergence (a match is worth `ln 4(1−p̂)` nats, a mismatch `ln(4p̂/3)`); the
  priors stay at their derivation λ.

## Extension points

Every component is a trait. `Assigner<A, RV, RJ, RD, W, Adm, Adj, Rf>` is generic
over all of them, and the defaults form `DefaultAssigner`. Swapping one changes
the type, and the whole pipeline monomorphises. There are no trait objects and
no runtime type checks.

| Trait | Default | Alternatives shipped |
|---|---|---|
| `Aligner` (+ `score_batch`, `align_batch`) | `InterSequenceGotoh` (AVX2, scalar fallback) | `BandedGotoh`; batch methods default to single-job loops |
| `Retriever` | V: `UngappedPrefilter<SeedVoteRetriever>` (prefilter off); J: `ExhaustiveRetriever`; D: `LmerMaskRetriever` | `SeedVoteRetriever` alone |
| `CandidateBudget` | `AmbiguityBudget{3, 16, 5 %, 8}` | `FixedBudget` |
| `SubstitutionModel` | `AffineScoring` | `HotspotScoring` (position-specific; the kernel's profile path) |
| `DWindowPolicy` | `WithAnchorWindow<ExpandedResidual>` | `ResidualWindow`, `ExpandedResidual` |
| `Admission` | `AnyOf<ExactRunFloor, Option<DistributedEvidence>>` (distributed-D off by default) | `AdmitAll`, `AllOf`, any `Option<A>` |
| `Adjudicator` | `PlacementDag<VdjGrammar, EngineeringObjective>` | `ResidualAdjudicator` |
| `Grammar`, `PathModel` | `VdjGrammar::{VDJ, TANDEM_D, VJ}`, `EngineeringObjective` | your own path grammar / objective |
| `Refiner` | `NoRefine` | extension point |
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

Four threads give ~9000 reads/s. Calls are identical across the
pure-speed steps; the steps that changed calls are listed in
[docs/ASSESSMENT.md](docs/ASSESSMENT.md) with their accuracy effect.

### Data and retrieval

**Nucleotides and k-mers as integers.** Rolling 2-bit k-mer keys in `u64`, with
N resetting the window. Spaced seeds are data (a mask), not code.

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

**Sibling upper bounds (`family_bound`).** After a gene's representative is
scored, a same-length sibling satisfies `S*(A) ≤ S*(R) + Δ(R, A)`, where Δ
sums the per-position maximum profile gain. Siblings whose bound cannot reach
the current top-N cutoff are not scored.

**Inter-sequence SIMD (`align/inter_seq`).**
- Each 16-bit AVX2 lane holds the same cell of a *different* job. Every lane
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
- Runtime AVX2 detection, with a scalar fallback on other CPUs.

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
can only add ≤ 0 to a path, and the strict with-D-beats-no-D rule rejects it.
The D stage therefore keeps only scores ≥ `P_D + 1` (and ≥ the gates' own
floor). Dropping a D leaves its neighbours' cuts valid, so this is exact.

### Adjudication

**O(1) clip-and-rescore.** Prefix sums over alignment columns give the score
of any sub-alignment between two query cuts in O(1), with exact affine gap
accounting at the new edges.

**DAG with running maxima.** Clip scores are separable in (entry, exit), the
N/P score is linear in the gap, and the support tie-break is separable too.
For every exit cut, the best predecessor state is a lexicographic prefix
maximum found by one binary search; for every entry cut of the next
placement, the best exit before it is a running maximum read in O(1). Every
pair of cuts is evaluated, and the result provably equals exhaustive
enumeration (tested).

**Forward–backward with the same sums.** The posterior pass replaces every
running maximum by a running sum. On an edge, all entries beyond the
predecessor's end see the same exit mass, which forward adds once to a
difference array and backward reads from a suffix sum, so an edge costs
O(overlap). Values are linear with one log-scale per slot; cut weights are
running products along the alignment and prior weights a table lookup, so
no transcendental function is evaluated per cut.

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

It extends IgSim in three ways:
- **Per-base germline origins.** Boundaries, junctions and germline-equivalent
  alleles are exact truth even across indels.
- **Optional reverse complementation.**
- **`SyntheticDatabase`.** An IGH-like germline with a family → gene → allele
  similarity hierarchy and conserved anchors, for tests that cannot ship
  IMGT-derived data.

`sim::Truth` holds one read's ground truth and reads and writes it as TSV.
`sim::Assignment` is a call, built from an `AirrRecord` or from any tool's
AIRR Rearrangement TSV row. `vdj-assign simulate` and `vdj-assign evaluate`
expose both, so other callers can be scored on identical reads.
`sim::Evaluation` scores calls against the truth:
- an allele call is *equivalent* when the primary allele contains the
  trimmed, pre-mutation germline body (no read can separate such alleles);
- *gene*, *listed* (truth among the reported names) and *exact*;
- boundary errors, D presence by remnant length, junction exactness;
- a random-read negative control.

Human IGH (541 V / 90 D / 16 J; a local copy, not shipped), 5000 reads per
preset, germline-equivalent accuracy in %, before → after the redesign:

| Preset (IgSim) | Orientation | V | V gene | D | J | Junction | Random reads called |
|---|---|---|---|---|---|---|---|
| easy (no SHM) | 99.94 → **100** | 99.54 → **99.84** | 99.92 → **99.98** | 98.78 → **99.27** | 99.90 → **99.94** | 99.94 → **100** | 30.5 → **0** |
| mid (IgM-like) | 100 → 100 | 97.76 → **97.94** | 99.84 → 99.84 | 91.96 → **93.00** | 98.84 → 98.80 | 99.93 → 99.93 | 30.2 → **0** |
| hard (IgG-like) | 99.98 → 99.98 | 89.74 → **89.82** | 98.44 → 98.24 | 81.72 → **83.53** | 92.70 → **93.32** | 98.39 → **98.51** | 31.2 → **0** |

The synthetic database (`SyntheticSpec::default()`, seed 1) gives, after:
easy V/D/J 99.98/99.82/99.98, mid 99.52/87.81/98.92, hard 97.10/79.36/94.00,
and junctions ≥ 99.8 %.

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
Brier score). The event is "the call is germline-equivalent to the truth"
for V, D and J, and "the read has a simulated D" for D presence:

| Preset | V | D | J | D presence |
|---|---|---|---|---|
| easy | 0.04 / 0.0008 | 0.46 / 0.0015 | 0.07 / 0.0000 | 0.25 / 0.0032 |
| mid | 0.28 / 0.0088 | 1.68 / 0.041 | 0.25 / 0.0060 | 2.17 / 0.028 |
| hard | 1.34 / 0.045 | 1.80 / 0.073 | 0.73 / 0.031 | 1.46 / 0.027 |

With one fixed λ instead of each read's divergence, and only co-optimal
placements pooled, mid-preset ECEs were 4.6 / 7.7 / 3.9 / 2.1.

## Validation

- **Posterior exactness:** `marginals_equal_exhaustive_summation` compares
  ln Z, every placement marginal and the D-count distribution with exhaustive
  summation over paths and canonical cuts (VJ, VDJ, VDDJ, a per-base N/P cost,
  a per-column evidence term); `posteriors_are_calibrated` bounds V/J
  calibration error on IgM- and IgG-like reads.
- **DAG optimality:** `dag_equals_brute_force_over_all_paths_and_cuts`
  enumerates every VJ/VDJ/VDDJ path and every pair of exit and entry cuts,
  with a per-base insertion cost and the lexicographic support tie-break, and
  requires the same optimum. `two_cuts_leave_a_shared_mismatch_stretch_to_n`
  checks a constructed read where no single cut is optimal.
- **Kernel exactness:**
  - The scalar kernel equals a textbook Gotoh across three gap regimes.
  - Every traceback rescores to its reported score.
  - The AVX2 kernel equals the scalar kernel on 160 random mixed batches. Its
    scores, end cells and complete tracebacks are compared over full and
    banded layouts, N bases, indels, four scoring schemes, and hotspot profiles
    mixed into half of the batches. Both the uniform and the profile path are
    covered.
  - The whole chunk-batched SIMD pipeline equals the per-read scalar pipeline
    on 400 simulated reads (paths, orientations and every placement).
- **Simulated accuracy (`tests/simulated.rs`):** IgSim presets on the
  synthetic database, half reverse-complemented. Accuracy floors sit ~1.5 points
  below the measured values for orientation, V/D/J allele and gene, junctions
  and long D remnants. Random reads must receive no call.
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
| `r` (default) | (2,−3,−9,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 4 | derived, even odds: 11 |
| `sensitive-d` | (2,−3,−9,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 4 | derived, odds 2: 9 |
| `truth` | (2,−3,−5,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 6 | 0 |
| `igblast` | (2,−3,−5,−1) | (2,−4,−11,−1) | (2,−4,−13,−1) | 5 | 0 |

Every profile gates V and J hits by E-value (accept 1, anchor 10⁻³) and
stitches the co-linear V prefix and J suffix back onto the local alignment
before ranking. A with-D path must strictly beat the no-D path; ties go to
fewer D segments, then to more germline support.

## Assumptions and deviations from the method note

Each is a named field you can change.

1. **D presence.** `P_D` is derived from the D database (`DPrior::Derived`)
   instead of the note's 12-relaxed-to-10 rule. A sweep of constant costs put
   the relaxed rule on the same precision/recall curve as a single constant.
2. **Tandem D.** A second D pays its own `P_D` plus 12 (`tandem_d_penalty`).
3. **Retrieval.** V keeps 3 genes and walks up to 16 while gene votes stay
   within `max(8, 5 %)` of the leader; every allele of a kept gene is a
   candidate. J is exhaustive. D uses the lossless L-mer prefilter instead of
   the note's k = 5/3 seeds.
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
  or the keep policy dropped gets no mass, and the call and the posterior use
  different evidence scales (fixed vs per-read), so a reported posterior can
  fall below 0.5. D-presence calibration is the weakest (ECE up to 2.2 %),
  because the D prior counts genes × positions (see ALGORITHM.md §8).
- **E-values** use exact *ungapped* λ and K and the total database length.
  Gapped scores are slightly higher (anti-conservative); near-identical alleles
  inflate the database length (conservative). They gate hits; they are not
  calibrated p-values.
- **Segment type is not checked from FASTA.** FASTA carries no segment type, so
  the CLI cannot detect a D file passed as `--j`.
- **SIMD targets.** The kernel is AVX2 only. Other CPUs, including ARM/NEON,
  use the exact scalar path.
- **Compressed input.** Gzip must be piped through `zcat`; there is no built-in
  decompression.
- **FR3 length rule.** The CDR2/FR3 boundary uses the FR3 length (38 codons
  heavy, 36 light). That length holds across human IGHV/IGKV/IGLV families
  but is unverified for TR loci. Give `fr3=` in the header where it differs;
  CDR2 lengths outside 0–12 codons are noted.
- **Chain type is one per database.** It is voted per database, so a mixed
  IGK + IGL database is fine, but a mixed heavy + light database needs
  `--chain` or separate runs.
- **FR1 start.** FR1 starts at the first complete codon. A V sequence that
  includes a leader needs `fr1=` in the header.
- **Ambiguous Cys23.** A 5′-truncated allele with no sibling and an adjacent
  decoy Cys may get a wrong CDR1 start. This is always flagged in
  `annotation.tsv`.
- **IUPAC codes.** Germline IUPAC codes are stored as N, and a written
  database spells them N.

## Roadmap

In order of expected value; reasons and measurements are in
[docs/ASSESSMENT.md](docs/ASSESSMENT.md).

1. **Calls from the per-read model.** Viterbi with per-read evidence weights,
   or maximum-posterior-marginal calls, so the call and its posterior agree.
2. **D prior from the D search's statistics** (Karlin–Altschul `K·m·n` and a
   D-inclusion rate), checked by D-presence calibration.
3. **Wider SIMD.** An AVX-512BW kernel (32 lanes) and a NEON port (8 lanes) of
   the same design.
4. **Allele rescoring along the gene frame.** Score siblings along their
   representative's traceback and run DP only where bounds straddle the cutoff.
5. **Annotation.** A per-database FR3 length learned from header-annotated
   alleles (TR loci). A consensus position prior for Cys23.
6. **C-segment assignment.**
