# Critical assessment of vdj_assign

This assessment covers what the code did at commit `a0ddd70`, what was
principled and what was patchwork, how each finding was measured, what
replaced it, and what remains. All numbers come from the simulation harness
added for this purpose (`src/sim`, `examples/benchmark.rs`,
`examples/diagnose.rs`). The real-database numbers use a human IGH germline
set (541 V, 90 D and 16 J alleles, including inferred novel alleles) that is
kept locally and not shipped. Reads follow the IgSim presets, with half of
them reverse-complemented.

## 1. How the tool was assessed

1. **Ground truth.** IgSim.jl was ported to Rust with its generative model and
   presets intact, and two additions:
   - each read base records its germline origin, so boundaries and junctions
     are exact even across indels;
   - reads can be reverse-complemented.

   A synthetic IGH-like database generator (family → gene → allele) makes
   the same tests shippable.
2. **Fair metrics.** An allele call counts as correct when it is
   *germline-equivalent* to the truth on the read: the called allele contains
   the trimmed, pre-mutation body that was simulated. D presence is
   stratified by remnant length, because a 2-nt remnant is not identifiable.
3. **Search vs model errors.** An oracle aligns every allele, unbanded, under
   the same objective. Errors where the oracle is right are search errors,
   which better engineering can fix. Errors where the oracle is wrong or tied
   are model errors, which no caller using this objective can fix.
4. **Ablation.** Each patch was switched off on its own and in combination.
5. **Profiling.** Callgrind (instruction counts) and per-stage wall clock.

## 2. Verdict at the start

**Principled and kept.**
- The placement DAG: semi-Markov Viterbi with exact clip-rescore through
  prefix sums, and separability turned into prefix maxima. It is provably
  equal to brute force.
- The lossless D prefilter, by pigeonhole on the admission gates.
- The sound sibling upper bound (`family_bound`).
- The inter-sequence SIMD kernel, exact with respect to the scalar kernel.
- Stage-wise batching, whose results do not depend on the chunk.
- Trait-based components with static dispatch.

**Patchwork, measured and resolved.** Accuracy deltas are on the real-database
hard preset, 3000 reads, unless stated.

| Mechanism | What it was | Measured effect | Resolution |
|---|---|---|---|
| Ungapped J path bonus | `(ungapped − min) / 4` added to J on reads with V identity ≤ 0.88 | none | removed |
| Switched J acceptance | soft gate during placement, strict gate re-applied when V identity > 0.88 | none | removed |
| J hotspot profile | mismatch −2 at AID hotspots, J only | +0.07 pp when removed (the simulator has no hotspots) | removed from defaults; `HotspotScoring` kept as an optional model |
| Junction tie | a zero-score 5-base block written onto V/J at ≤ 0.88 identity | none on boundaries | removed |
| D short-run / long-span surcharges | `1 × (6 − run)` presence surcharge | +0.3 pp D when removed | removed |
| `P_D` 12, relaxed to 10 | relaxed when the score ≥ 20 or the exact tract occurs in two D sequences | on the same precision/recall curve as a constant | derived `P_D` (§3.4) |
| `supplement_all_alleles` for J | J seed retrieval, then every allele added unbanded anyway | redundant | exhaustive J by design (§3.2) |
| Two keep code paths | "screen" for V/J and "select + slack" for D | — | one keep algorithm |
| `band_min` 28 | at least 57 diagonals per V band | the same calls with 17–25 | band = seed cluster ± 8 |

**Defects found.**
1. **Mis-orientation of unmutated reads (0.06–0.10 %).** The reverse strand
   was skipped once forward V votes reached 160. Votes counted postings, so
   the antisense of IGHJ6's `TACTACTAC` repeat collected ~200 votes on
   IGHV3-48's `AGTAGTAGT` CDR2, and the true strand (1084 votes) was never
   examined.
2. **Arbitrary tie-breaking picks the wrong allele.** When two V alleles
   differ only in 3′ length, a junction base can belong to the V or to the D
   with equal path score. The report then chose by exact-run length and allele
   id, so the 2-nt-shorter allele won. Unmutated reads lost 0.2–0.3 % of V
   calls this way, and the same mechanism affected D.
3. **Random reads were assigned.** 30 % of random reads received a path.
   Nothing gated significance, and any long random sequence holds a
   10-column, 60 %-identity J alignment. Once J became exhaustive, 93 % of
   random reads had a J call and 1.6 % a V call, before the gates of §3.8.
4. **Wide D windows hurt.** A D's optimal local alignment in a 24-base-reclaim
   window could trade its junction core for chance matches in a V tail, which
   the cut then removed.

**Speed profile (real database, mid preset).** 849 reads/s per thread. Per
read:
- orientation + V retrieval: 279 µs; the `u64` sort of ~50 000 votes,
  inflated by alleles sharing every seed, dominated.
- V DP: 371 µs; every allele of every retrieved gene was scored in bands of
  ≥ 57 diagonals. The family bound pruned nothing, because its cutoff is the
  3rd-best score, which every sibling of the top genes clears.
- J: 104 µs, every allele across the whole read.
- D: 192 µs.

## 3. What replaced it

Each change was measured on the same reads before it was kept.

1. **Strand by germline k-mer presence** (`orient.rs`). Count, for each strand,
   the read positions whose 11-mer occurs in the V/J germline. One bit test
   per base; each position counts once. The count is trusted when it is
   significant against its Poisson null, and otherwise falls back to seed
   votes on both strands. Easy-preset orientation is now 100 %.
2. **J searched exhaustively right of V** (`retrieve/exhaustive.rs`,
   `window::DownstreamWindow`). This replaced J seed retrieval and all three J
   patches. J on the hard preset rose 92.53 → 93.20 %.
3. **Germline support as the last tie-break** (`adjudicate::PathScore`,
   `report.rs`). It is a cut-marginal argument: the better-supported allele
   admits more equal-score cuts. The DAG stays exact, because the key is
   lexicographic and separable. V on easy rose 99.63 → 99.83 %, D +0.3 to
   +0.9 pp.
4. **Derived `P_D`** (`config::DPrior`). Scores are log-odds (λ = 0.634, i.e.
   11.2 % divergence). A D hypothesis is spread over `N` genes × `L`
   positions, so `P_D = (ln(N·L) − ln odds)/λ`: 11 for human IGH, which is the
   note's constant recovered from first principles, and 6 for TRB. A per-read
   `L` was tried and traced the same curve, so the simpler constant per
   database was kept.
5. **Gene frames + radix-sorted votes** (`retrieve/frames.rs`,
   `retrieve/seed_vote.rs`). Postings and votes now scale with genes. Joining a
   frame requires position-consistent seeds: a first version let an allele
   with an internal codon deletion join its frame, which voided half the gene's
   seeds, and the consistency rule fixes that. V retrieval went 267 → 73 µs.
6. **Seed-span bands** (± 8 diagonals) halve V DP time with identical calls.
7. **D floor from the objective** (Lemma 1 in ALGORITHM.md) and a **12-base
   reclaim**. D on the hard preset rose 83.1 → 83.5 % (mid: 92.5 → 93.2 %),
   and D time fell ~20 %.
8. **E-value gates with anchoring** (`segment::EvalueGates`). Hits are kept at
   E ≤ 1, and a call needs a V or J hit at E ≤ 10⁻³. Random reads called fell
   from 30 % to 0 %. A flat 10⁻³ floor had cost 0.3 pp of junctions on the
   hard preset; conditioning on the anchor recovers them, and junctions now
   beat the original.
9. **Uniform-profile SIMD fast path** and a **three-allocation clip index**:
   +10–20 % each, identical calls.
10. **Structure.** `assign_many` became five named stages. `segment.rs` has
    one keep algorithm. `config` separates V (seeded) and J (exhaustive)
    settings and declarative objective settings (`ObjectiveConfig`) from
    resolved numbers (`EngineeringObjective`).

## 4. Results

Human IGH, 5000 reads per preset, germline-equivalent accuracy in %:

| Preset | Orientation | V | V gene | D | J | Junction | Random reads called |
|---|---|---|---|---|---|---|---|
| easy | 99.94 → 100 | 99.54 → 99.84 | 99.92 → 99.98 | 98.78 → 99.27 | 99.90 → 99.94 | 99.94 → 100 | 30.5 → 0 |
| mid | 100 → 100 | 97.76 → 97.94 | 99.84 → 99.84 | 91.96 → 92.93 | 98.84 → 98.80 | 99.93 → 99.93 | 30.2 → 0 |
| hard | 99.98 → 99.98 | 89.74 → 89.80 | 98.44 → 98.24 | 81.72 → 83.55 | 92.70 → 93.32 | 98.39 → 98.47 | 31.2 → 0 |

Throughput is 849 → ~2200 reads/s per thread, and ~3600 → ~9000 reads/s on 4
threads.

Error decomposition after the changes (2000 reads):

| | search errors | ties | model errors |
|---|---|---|---|
| V, mid | 0 | 1.3 % | 1.3 % |
| V, hard | 0.4 % | 5.5 % | 5.3 % |
| J, hard | 0.1 % (was 0.8 %) | 2.8 % | 3.7 % |

V and J are at the limit of the scoring model. Under i.i.d. SHM the
`(2, −3)` score ranks alleles by mismatches over the same region, which is the
maximum-likelihood rule. The remaining errors are mutations that land on an
allele's distinguishing positions.

## 5. What is still weak

These are listed honestly, for the paper's limitations section and for future
work.

1. **V gene at extreme SHM: −0.2 pp on the hard preset** against the start.
   Gene-frame votes are unions over alleles, and at 15–25 % SHM the gene budget
   sometimes ranks a sibling gene above the true one. Search errors went from
   0.1 % to 0.4 % of reads. A base budget of 5 genes instead of 3 recovers
   it (98.44 %) at about 10 % of throughput, so the default stays at 3 and the
   budget is the knob. Votes normalised by frame size are the alternative to
   try.
2. **D presence on D-free junctions.** 12–24 % of reads with no D receive one,
   because a 9-nt chance match is indistinguishable from a 9-nt D remnant.
   `P_D` sets the trade-off, and the prior odds are the single knob.
   Precision/recall curves, not one number, belong in the paper.
3. **One alignment per D allele and window**, and **a single cut per
   overlap**. Both are approximations of the segmentation model (§8 of
   ALGORITHM.md).
4. **Constants still inherited from the note.** These have not been derived
   and should get a sensitivity analysis:
   - `pair_bonus` 24 and `tandem` 12;
   - acceptance at 60 % identity with 24/10 columns;
   - the budget (3–16 genes, 5 %/8 votes) and seed lengths 9/7;
   - `max_occurrences` 256;
   - distributed-D thresholds.
5. **E-values are approximate.** They use ungapped λ/K (anti-conservative for
   gapped scores) and the total database length (conservative for
   near-identical alleles). They gate hits; they are not p-values.
6. **Simulator realism.** IgSim has no AID hotspot targeting, no
   codon-length indels and no allele-frequency skew, so hotspot-aware scoring
   cannot be judged on it. Validation on real data with an independent caller
   (IgBLAST, partis) is still needed.
7. **SIMD portability.** AVX2 only; AVX-512 hardware runs at half width, and
   ARM falls back to scalar.

## 6. Plan

In order of expected value:

1. **Posterior decoding** over the DAG (log-sum-exp prefix sums replace prefix
   maxima, plus a backward pass). It gives calibrated per-allele and D-presence
   probabilities and maximum-posterior-marginal calls. It turns limitation 2
   into a reported confidence rather than a hard call, and generalises the
   support tie-break.
2. **Two cuts per overlap with N/P between**, so the DAG is exact for the
   segmentation model and the D floor is exact unconditionally.
3. **Sensitivity analysis** of the inherited constants on the simulator,
   across presets and both databases, with each constant either derived or
   shown to be flat.
4. **Cross-tool benchmark.** Run IgBLAST, partis and Swig through
   `vdj-assign simulate` / `evaluate` (implemented: any AIRR TSV is scored with
   the same metrics on identical reads) at several SHM levels, and report
   precision/recall curves for D.
5. **AVX-512BW (32 lanes) and NEON (8 lanes)** kernels of the same design.
6. **Allele rescoring along the gene frame.** Trace one representative per
   frame, rescore siblings along its path, and run DP only where bounds
   straddle the top-N cutoff. Most of the remaining V time is siblings of the
   top gene.
7. **Calibrated statistics.** Fit gapped λ/K and an effective database size
   by simulation.
