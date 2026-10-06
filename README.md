# vdj_assign

Seeded, banded, affine V(D)J germline assignment with placement-DAG adjudication:
a Rust implementation of the SwiftIG / Swig 0.38.x method note, built so that
every heuristic is a replaceable, statically dispatched component.

Zero dependencies, Rust ≥ 1.75, edition 2021.

## Build and run

```sh
cargo build --release                             # always benchmark/run the release build
cargo test                                        # 45 unit + 7 annotation + 10 integration tests
cargo run --release --example throughput 20000    # synthetic IGH-sized benchmark

target/release/vdj-assign --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta \
    [--profile r|truth|igblast|sensitive-d] [--threads N] [--tandem-d] \
    [--chain heavy|light] [--forward-only] [--no-anchor-window] [--residual] reads.fastq > calls.tsv

zcat reads.fastq.gz | target/release/vdj-assign --v V.fasta --d D.fasta --j J.fasta - > calls.tsv

# optional: write the annotated database for review or curation
target/release/vdj-assign prepare --v IGHV.fasta --d IGHD.fasta --j IGHJ.fasta --out prepared/
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
also enable the anchor-bounded D window. Output is an AIRR Rearrangement TSV. Coordinates are 1-based closed
intervals on the oriented query. CIGARs use `S`/`N` for unaligned query and
germline bases and `=`/`X`/`I`/`D` for columns. `np1`–`np3` hold the
junctional N/P regions, and `*_support` holds Karlin–Altschul E-values.
Region columns follow: `fwr1`–`fwr4`, `cdr1`–`cdr3` (sequence, `_start`,
`_end`), `cdr3_aa`, `junction`, `junction_length`, `junction_aa`,
`vj_in_frame`, `stop_codon` and `productive`.

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
clip index. This costs 0.72 µs per read, against ~326 µs for assignment and
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
 chunk of reads ─► per read: orient (V votes, fwd vs rc; J breaks ties), retrieve V and J candidates
                ─► V stage, whole chunk: one score batch + one traceback batch (16 jobs per AVX2 register)
                ─► J stage, whole chunk
                ─► per read: D windows from V/J geometry (+ anchor window), D candidates
                ─► D stage, whole chunk ─► per read: admission gates
                ─► per read: placement DAG over VJ / VDJ / VDDJ ─► refiner ─► Call ─► AIRR
```

Every per-read decision is still taken per read. Batch aligners must answer
each job exactly as a single call would, so a read's result never depends on
its chunk: `assign(q)` is literally a chunk of one, and this is tested.

| Module | Responsibility | Method-note section |
|---|---|---|
| `alphabet`, `fasta`, `germline` | 2-bit codes (N = 4), germline FASTA, germline sets with anchors | Input |
| `reads` | Streaming FASTA/FASTQ reads with format detection | Input |
| `translate` | 64-entry codon table on 2-bit codon indices; one-pass three-frame translation | Annotation |
| `annotate`, `annotate/{v_motif,j_motif,transfer}` | Anchors, frames, V region starts: header > motif > consensus > sibling transfer | Annotation |
| `prepare` | Deduplicate + annotate; write annotated FASTA and `annotation.tsv` | Annotation |
| `read_regions` | Lazy per-read FR/CDR/junction projection, frame, stop codon, productive | Reporting |
| `scoring` | Affine tuples, substitution profiles; N scores 0 | Scoring tuples |
| `seed`, `index` | Rolling k-mer integers, spaced seeds, CSR k-mer index | Seeding |
| `retrieve/seed_vote` | Packed-vote diagonal clustering, per-gene budget | Candidate retrieval, votes |
| `retrieve/budget` | Fixed / ambiguity-triggered candidate budget (4 → 16) | Budget |
| `retrieve/ungapped` | Ungapped diagonal pre-ranking decorator | Acceleration (new) |
| `retrieve/lmer_mask` | D prefilter: per-L-mer allele bitsets | D retrieval |
| `band`, `align/gotoh` | Diagonal bands (bins of 4, w = clamp(span+8, 28, 256)); scalar banded Gotoh | Aligner, banding |
| `align/inter_seq` | AVX2 inter-sequence Gotoh, 16 jobs per register, exact scalar semantics | Aligner (acceleration) |
| `clip`, `evidence`, `placement` | O(1) clip-and-rescore at any cut; exact-run / identity evidence | Boundary scoring |
| `segment` | Retrieve → pre-score → keep → traceback → merge identical, batched over many reads | Per-segment search |
| `window` | Where D is searched: expanded residual, anchor window | Residual D search, rescue |
| `admission` | Exact-run floor (`min_d_match`), distributed-D gate, combinators | D admission, AER-R |
| `adjudicate/dag` | Semi-Markov Viterbi over placements with exact cut choice | Target architecture |
| `adjudicate/model` | Grammar (VJ/VDJ/VDDJ), path model, D-presence penalty `P_D` | AER-R objective |
| `adjudicate/residual` | Ordinary commit-V×J-then-D path | Ordinary path |
| `refine` | Post-adjudication polishing of the winning path | Boundary refinement |
| `pipeline`, `config`, `profiles` | Composition, configuration, calling profiles | Profiles |
| `stats`, `report`, `batch` | E-values, AIRR output, parallel batches | Reporting |

## Extension points

Every component is a trait. `Assigner<A, RV, RJ, RD, W, Adm, Adj, Rf>` is generic
over all of them, and the defaults form `DefaultAssigner`. Swapping one changes
the type, and the whole pipeline monomorphises. There are no trait objects and
no runtime type checks.

| Trait | Default | Alternatives shipped |
|---|---|---|
| `Aligner` (+ `score_batch`, `align_batch`) | `InterSequenceGotoh` (AVX2, scalar fallback) | `BandedGotoh`; batch methods default to single-job loops |
| `Retriever` | `UngappedPrefilter<SeedVoteRetriever>` (V/J), `LmerMaskRetriever` (D) | `SeedVoteRetriever` alone, `params: None` |
| `CandidateBudget` | `AmbiguityBudget{4, 16, 5 %, 8}` | `FixedBudget` |
| `DWindowPolicy` | `WithAnchorWindow<ExpandedResidual>` | `ResidualWindow`, `ExpandedResidual` |
| `Admission` | `AnyOf<ExactRunFloor, Option<DistributedEvidence>>` | `AdmitAll`, `AllOf`, any `Option<A>` |
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

Measured single-threaded in a 1-core container, on a synthetic IGH-sized
database. The database has 7 V families × 25 genes (~15 % within-family
divergence) × 2 alleles = 350 V, 27 D and 12 J alleles. Reads carry 5 % SHM
and half are reverse-complemented (`examples/throughput.rs`).

| Stage of optimisation | reads/s/thread |
|---|---|
| First complete pipeline | 508 |
| + ungapped V/J pre-ranking, lossless D L-mer and score floor | 1162 |
| + branch-free scalar kernel with shortened dependency chain | 1393 |
| + pre-sized traceback, arithmetic trace codes, conditional pre-score pass | 1724 |
| + AVX2 inter-sequence kernel, batched per read (~3 V jobs per register) | ~2000 |
| + stage-wise batching across 64-read chunks (~190 V jobs per batch) | **3345** |

In-situ DP time per read, scalar versus chunk-batched SIMD:

| Segment | Scalar | SIMD, chunk-batched | Speedup |
|---|---|---|---|
| V | 188 µs | 38.5 µs | 4.9× |
| J | 36 µs | 10.0 µs | 3.6× |
| D | 132 µs | 16.9 µs | 7.8× |

DP is now ~65 of ~280 µs per read. The rest is retrieval (V seed votes on
both strands, ~80 µs), placement construction and adjudication; that is
where the next speedups are.

Kernel-level numbers are from 400 *distinct* V-like jobs; repeating one job
lets the branch predictor memorise it and flatters scalar code:

| Jobs per register | Score-only | Traceback |
|---|---|---|
| Scalar | 2.40 ns/cell | 3.91 ns/cell |
| 16 (full) | 0.53 ns/cell (4.6×) | 0.84 ns/cell (4.7×) |
| 8 | — | 1.39 ns/cell |
| 4 | — | 2.72 ns/cell |

### Data and retrieval

**Nucleotides and k-mers as integers.** Rolling 2-bit k-mer keys in `u64`, with
N resetting the window. Spaced seeds are data (a mask), not code.

**CSR k-mer index.** Offsets and packed `u32` postings (allele, position).
Keys use a direct-address table up to 24 key bits, and a sorted key array with
binary search beyond. Frequent seeds (> 256 occurrences) are skipped, as in the
note.

**Packed votes with sort-and-scan clustering.** Each seed hit is one `u64`
(allele | diagonal bin | weight). A single `sort_unstable` followed by a linear
scan yields diagonal clusters per allele: O(h log h), no hash maps.

**Ungapped before gapped (BLAST staging).**
- Each V/J candidate is scored by Kadane's maximum segment sum along the
  diagonals around its band centre, at about 1/50 of the banded DP cost.
- Only candidates within `slack` of the best, and always the top `keep_top`,
  reach DP.
- This is a lossy speed knob (defaults: top 3, slack 30); `None` disables it.

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
- Together with the branch-free loop: 6.2 → 2.05 ns/cell on real D jobs.

**Score first, traceback later — only when it can prune.**
- Score-only DP ranks candidates, and only kept ones are traced.
- When a stage has no more candidates than it keeps (`prescore_above`), it
  aligns directly: the pre-pass could not prune.

**Inter-sequence SIMD (`align/inter_seq`).**
- Each 16-bit AVX2 lane holds the same cell of a *different* job. Every lane
  runs the ordinary recurrence, so traceback codes are the scalar ones and the
  scalar traceback walk is reused verbatim. There is no striped lazy-F loop and
  none of its traceback problems.
- Diagonal bands use band-relative coordinates (row i, offset k = i − j − lo),
  so all lanes walk the same grid even though their bands start on different
  diagonals. The reference position of cell (i, k) depends only on i − k.
  Profile columns are laid out once per batch, interleaved, and each cell's
  substitution scores are then one contiguous load per base code.
- Unbanded jobs (D windows) use a rectangular layout.
- Scores come from the per-allele profiles, selected by per-row query-code
  masks. This makes the kernel exact for any `SubstitutionModel`, not only
  match/mismatch.
- Out-of-band and out-of-matrix cells read as H = 0, I = D = −∞, exactly like
  the scalar kernel. Narrower bands are masked per lane. Cells past a lane's
  query or reference get a −∞ substitution, so they can never feed back or
  become a strict new maximum.
- Saturating i16 arithmetic, with a per-job admission check on the attainable
  score; jobs that could approach the sentinel fall back to the scalar kernel.
- Runtime AVX2 detection, with a scalar fallback on other CPUs.

**Stage-wise batching across reads.** After prefiltering, a read has ~3 V, ~3 J
and ~9 D candidates, far too few to fill 16 lanes. So the pipeline runs every
stage over a chunk of 64 reads (the batch scheduler's chunk) and hands the
aligner one batch per stage. Jobs are sorted by shape into registers. Keep
policies, D windows, admission and adjudication stay per read.

### D search

**Lossless D prefilter.** The prefilter uses per-L-mer allele bitsets (one OR
per query L-mer). L and the traceback score floor are derived from the
admission gates, so they can never drop an admissible hit:
- L = `min_d_match` when only the exact-run gate is on.
- L = 4 when distributed-D is on: ≥ 80 % identity forces an exact run ≥ 4, by
  pigeonhole.
- Score floor = `m · min_d_match`, because an exact run of k scores ≥ k·m.
- When L·m already meets that floor, every retrieved candidate passes it, so
  the score-only pre-pass is skipped unless there are more candidates than
  the keep policy's `top`.

### Adjudication

**O(1) clip-and-rescore.** Prefix sums over alignment columns give the score
of any sub-alignment between two query cuts in O(1), with exact affine gap
accounting at the new edges.

**DAG with prefix maxima.** Clip scores are separable in (entry, exit). So for
every overlap cut, the best predecessor state is a prefix maximum found by one
binary search. Every cut through every overlap is evaluated, and the result
provably equals brute-force enumeration.

### Parallelism and allocation

**Workspaces and batching.** One workspace per thread means no DP allocation
after warm-up. Work is self-scheduled in chunks of 64 from an atomic counter
over scoped threads. Results come back in input order, and reporting is formatted
on the workers (`assign_map`).

## Validation

- **DAG optimality:** `dag_equals_brute_force_over_all_paths_and_cuts`
  enumerates every VJ/VDJ/VDDJ path and every cut, with a per-base insertion
  cost, and requires the same optimum.
- **Kernel exactness:**
  - The scalar kernel equals a textbook Gotoh across three gap regimes.
  - Every traceback rescores to its reported score.
  - Wide bands reproduce the full matrix.
  - The AVX2 kernel equals the scalar kernel — scores, end cells and complete
    tracebacks — on 120 random mixed batches of 1–40 jobs. These cover full
    and banded layouts, band widths 7–66, N bases, indels and four scoring
    schemes.
  - The whole chunk-batched SIMD pipeline equals the per-read scalar pipeline
    on 400 simulated reads (paths, orientations and every placement).
- **Input parsing:** multi-line FASTA/FASTQ, `@`/`+` at the start of quality
  lines, empty reads, CRLF, batch resumption, and every error case.
- **Karlin–Altschul:** the series reproduces K = 1/3 for the ±1 walk and
  NCBI's ungapped blastn table (reward 1 / penalty −3 → λ 1.374, K 0.711, H 1.31).
- **Synthetic calling accuracy:**

| Test | Result |
|---|---|
| D gene | 60/60 correct |
| D-free junctions | false D in 0/40 |
| Tandem VDDJ, both D genes | 30/30 correct |
| V gene at 8 % SHM | always correct |
| Reverse-complemented reads | identical calls and CIGARs |

- **Batch consistency:** chunk-batched parallel output equals per-read
  sequential output.
- **Germline annotation and read regions:** see *Reliability* under
  *Germline annotation and regions*.

## Calling profiles

The default is Swig's default: AER-R adjudication with the R-optimized profile.
`r_optimized` is accepted as an alias of `r`.

| Profile | V tuple | D tuple | J tuple | `min_d_match` | Distributed-D | `P_D` | Junction tie |
|---|---|---|---|---|---|---|---|
| `r` (default) | (2,−3,−9,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 4 | on | 12, or 10 if score ≥ 20 or the exact tract is in two D sequences | on |
| `sensitive-d` | (2,−3,−9,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 4 | on | 10 fixed | on |
| `truth` | (2,−3,−5,−1) | (2,−3,−13,−1) | (2,−3,−17,−2) | 6 | off | 0 | off |
| `igblast` | (2,−3,−5,−1) | (2,−4,−11,−1) | (2,−4,−13,−1) | 5 | off | 0 | off |

Every profile requires 60% identity, at least 24 V substitution columns and 10 J
substitution columns, and stitches the co-linear V prefix and J suffix back
onto the local alignment before ranking. A with-D path must strictly beat the
no-D path; ties go to fewer D segments.

## Assumptions and deviations — please review

The note leaves several of these unspecified; each is a named field you can change.

1. **Tandem D.** A second D costs an extra 12. Swig's double-D screen is a
   separate opt-in and is not this penalty.
2. **Retrieval constants.**
   - V keeps 3 genes and walks up to 16 while votes stay within `max(8, 5%)`.
   - J keeps 3 genes.
   - Traceback keeps the top 3 V and top 2 J local scores, plus whatever is tied with the last of those. A hit that fails the length or identity floor does not occupy a slot.
   - V family-delta: score one allele per gene first; skip same-length siblings whose profile upper bound cannot enter that top-N. J leaves this off.
   - `weak_votes` 16 triggers the fallback seed.
   - `max_occurrences` 256.
   - Minimum substitution columns: V 24, J 10, D `min_d_match`. Identity floor 0.60.
3. **D keep.** D keeps a wide lossless set (top 16, cap 32) so a distributed
   hit below the ordinary top-2 is not dropped. Swig reaches the same hits
   with a second D pass.
4. **IgBLAST J tuple** is `(2, −4, −13, −1)`, matching Swig's agreement profile.
5. **Boundary rescue.** The 24/10 rescue is realised as `ExpandedResidual`
   (reclaim 24 into V/J) plus the DAG's exact cut choice, not as a separate
   post-hoc step.
6. **Anchor window.** It is always added when anchors are present, not only
   when the residual is weak.
7. **D retrieval.** The note's k = 5/3 D seeds are replaced by the lossless
   L-mer bitset prefilter.
8. **Distributed-D.** On for `r` and `sensitive-d`. Off for `truth` and `igblast`,
   which are the legacy non-AER-R profiles.
9. **Ungapped pre-ranking.** Available, and off by default, because it can drop
   an allele Swig would align. Set `ungapped: Some(UngappedParams::default())` to enable it.
10. **Gap cost constraint.** `GapCosts::new` requires `open ≤ extend ≤ 0`. With
    cheaper opening, Gotoh prices gap runs as repeated opens and the "consecutive
    gap columns extend" convention (used by clipping and reporting) no longer
    describes the optimum. DP scores themselves stay exact.

## Known limitations

- **No reporting threshold.** A random 120-nt read receives weak V and J calls.
  Its V E-value (~0.2 versus ~1e-144 for real reads) exposes this, but nothing
  suppresses it yet. An E-value or minimum-score gate per segment is the
  obvious next component, either as a `Refiner` or in `report`.
- **Approximate E-values.** They use exact *ungapped* parameters and are
  mildly anti-conservative for gapped alignments. Calibrate λ and K by
  simulation for real use. They are never used for calling.
- **Segment type is not checked from FASTA.** FASTA carries no segment type, so
  the CLI cannot detect a D file passed as `--j`.
- **Benchmark conditions.** Measurements come from a 1-core container on
  synthetic data; real IGHV families and read lengths will shift them.
- **SIMD targets.** The kernel is AVX2 only; other CPUs, including ARM/NEON,
  use the exact scalar path. A NEON port is 8 lanes of the same design.
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

## Not implemented (roadmap)

**Faster non-DP work.** It is now ~75 % of the time.
- V retrieval on both strands (~80 µs per read).
- Placement construction.

**Retrieval.**
- Minimizers / syncmers.
- BMI2 `pext` for spaced seeds.
- Suffix automaton / matching statistics for D.

**Kernels.**
- Cheaper substitution selection (byte-shuffle lookup instead of five masked loads).
- AVX-512 (32 lanes) and a NEON port.
- WFA for near-identical V ends.

**Annotation.**
- A per-database FR3 length learned from header-annotated alleles (TR loci).
- A consensus position prior for Cys23 relative to other conserved residues,
  so that truncated alleles benefit without a sibling.

**Engineering.**
- Arena allocation for placements.
- A reporting threshold on E-value / minimum score.
- C-segment assignment.
