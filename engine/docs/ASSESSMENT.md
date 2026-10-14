# Critical assessment of swig

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
   database was kept. It is now kept unrounded (`P*_D` = 10.55 for human
   IGH): the fixed scale charges `⌊P*_D⌋`, which decides `S > P*_D` exactly
   on integer scores, the posterior `P*_D` itself, and the per-read call
   reads it at the read's scale (§9.8).
5. **Gene frames + radix-sorted votes** (`retrieve/frames.rs`,
   `retrieve/seed_vote.rs`). Postings and votes now scale with genes. Joining a
   frame requires position-consistent seeds: a first version let an allele
   with an internal codon deletion join its frame, which voided half the gene's
   seeds, and the consistency rule fixes that. V retrieval went 267 → 73 µs.
6. **Seed-span bands** (± 8 diagonals) halve V DP time with identical calls.
7. **D floor from the objective** (Lemma 1 in ALGORITHM.md) and a **12-base
   reclaim**. D on the hard preset rose 83.1 → 83.5 % (mid: 92.5 → 93.2 %),
   and D time fell ~20 %. The floor is now applied only without posteriors
   (item 14): it never changes a call (verified: identical D calls with and
   without it on 5000 mid and hard reads), but a D below it still carries
   probability.
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
11. **Two cuts between neighbours** (`adjudicate/dag.rs`). The DAG used one
    cut per overlap with no N/P between overlapping placements, an
    approximation of the segmentation model that also made the D floor exact
    only when D's neighbours did not overlap. Each placement now has its own
    exit and entry cut (`x ≤ y`, N/P between). A running maximum over exits
    keeps the cost per entry O(1). On 3000-read benchmarks, accuracy moved by
    at most ±0.2 pp (real hard junction 98.28 → 98.35 %, synthetic hard D
    79.02 → 79.24 %) and throughput was unchanged within noise.

    *This version was not yet exact.* It tried entries only inside the
    overlap, arguing that a later entry "clips a prefix of an optimal local
    alignment", which scores ≥ 0. Adversarial review showed the argument
    wrong: the clipped bases are a middle stretch of the next alignment, not
    a prefix, and can score below zero (mismatches right after V's end).
    The brute-force oracle of the time searched 12 bases past the overlap,
    wide enough to see it, but its synthetic reads never put mismatches
    right after a neighbour's end, so it never fired. Item 13 fixes the DAG
    and adds an oracle test on hypermutated IgSim reads, which fails on the
    restricted DAG.
12. **Posterior decoding** (`adjudicate/posterior.rs`). Forward–backward
    over the same DAG gives the marginal probability of every placement and
    D count under `P(path, cuts) ∝ exp(nats)`; it equals exhaustive
    summation (tested). Two modelling decisions made the probabilities
    calibrated, each found by inspecting miscalibrated reads:
    - *Evidence at the read's divergence.* With one λ, `(2, −3)` scores are
      log-odds at 11.2 % divergence only, so V/J posteriors on IgM-like reads
      were underconfident (mid-preset ECE V 4.6 %, J 3.9 %). Rescoring each
      `Match` column as `λ_r·s + μ_r`, the log-odds at the read's own
      mismatch rate, brought them to 0.4 % and 0.8 %. Priors keep their λ.
    - *The event users are scored on.* Pooling only co-optimal placements
      left D underconfident (77 % predicted, 99.7 % observed on easy reads):
      a runner-up allele that lacks one chance-matching flank base is
      germline-equivalent to the call on the true remnant. The reported
      posterior is therefore the mass of every placement whose germline
      bases, between the call's neighbours, the called allele contains.
      Clipping the competitor to the called interval instead was too lenient
      (four bases occur in any D) and made D overconfident.
    The pass costs 18 % of instructions and about 30 % of single-thread wall
    clock (`--no-posterior` skips it). Final calibration after item 14 is in
    §4.
13. **Exact Viterbi over every entry cut** (`adjudicate/dag.rs`, rewritten).
    A node's slot (D count, first kind) now holds one value per canonical
    entry cut, and the recursion is the max-plus twin of the posterior's:
    a running maximum over entries for each exit, a running maximum over
    exits per successor kind for each entry, and, for entries beyond the
    predecessor's end, one value added at the first such entry and carried
    by a prefix maximum. An edge still costs `O(overlap + 1)`. The oracle in
    `tests/synthetic.rs` now enumerates every canonical exit and every
    canonical entry with no window. A new oracle test runs it on the
    pipeline's own placements for 40 hypermutated IgSim reads; it asserts
    that some optimum enters beyond the overlap, and it fails on the
    restricted DAG. Viterbi is now exactly the mode of the posterior at the
    fixed λ.

    Effect on identical 5000-read samples: V and J calls unchanged; D
    accuracy +0.00 to +0.07 pp; D calls on D-free junctions +0.4 to
    +0.8 pp (real mid 24.7 → 25.4 %, hard 12.2 → 12.9 %). The exact
    optimum finds a few with-D paths the restricted search could not reach:
    the objective, through `P_D`, decides them now. Cost: +3 % instructions
    over a whole run.
14. **Posterior numerics and evidence, from adversarial review.** Each
    finding was confirmed by an independent verifier before it was fixed.
    The numerical ones are now covered by exhaustive-summation tests (a
    junction with an `N` column and an insertion; a strong N/P cost across a
    300-base V); the reporting ones only by the calibration tests:
    - *`N` columns* scored 0 under the substitution scores but were given
      `μ_r` under per-read evidence; they now carry 0 nats. *Gaps* were
      rescaled by `λ_r`, although `p̂` measures substitutions; they stay at
      the prior λ.
    - *Range.* One centre per placement put a long V's far cuts outside the
      `f64` range, and one such placement voided the whole read's
      posteriors. Entry and exit weights are now centred separately over the
      cuts a path can use, and junction recursions keep only factors ≤ 1
      (a strong per-base N/P cost across a 300-base V had overflowed).
    - *D floor.* Dropping Ds below `P_D + 1` is lossless for the maximum but
      not for the sum; posteriors now keep every admissible D.
    - *Reported event.* A tandem D could be counted twice, once for each
      region; a placement now counts only if it overlaps the call's region
      (tightened in item 17). The event is "the reported names include the
      read's allele, up to germline equivalence", which is what §4 scores.
15. **Closed-syncmer seeds** (`seed.rs`), adopted here and reverted in
    item 18. See §8.
16. **3′-partial V alleles got a wrong anchor** (`annotate/mod.rs`). Found
    through the IgBLAST comparison (§8), whose junctions disagreed with the
    truth on 4.5 % of reads. IMGT 3′-partial alleles end inside or before
    their Cys104 codon (IGHV4-31\*09 ends `…AVYY`). The sibling transfer
    rejected the anchor because the codon did not fit, and the annotation
    fell back to the weak motif hit, an out-of-frame TGT about 30 nt
    upstream. Because the simulator's truth uses the same annotation, the
    harness had scored those junctions as correct. Now, when the closest
    trusted relative's anchor maps at or past the 3′ end, the allele gets no
    anchor and a `3′-partial` note (22 V alleles of the human IGH set). The
    J pass follows the same rule. A test builds a partial allele with a TGT
    decoy; it fails on the old code. Item 17 refines the rule.
17. **Second adversarial review** (of items 13–16; each finding reproduced
    by an independent verifier, or by me where the verifier failed):
    - *A D could not precede a J that starts earlier.* Nodes were sorted by
      read start, and an edge `b → e` needed `e` to start after `b`. A J
      alignment reaching back over the D, the "J eats D" configuration, so
      excluded V → D → J with that J. Both brute-force oracles had the same
      restriction. Nodes are now sorted by kind, then start: transitions
      between kinds depend only on the cuts. A constructed test fails on the
      old order. Measured: missing Ds covered by the called J fall to 0.00 %
      on mid and hard (`d_errors`), and reads whose called J covers ≥ 3
      bases of the true D from 0.44 / 0.95 % when first measured (before
      item 13) to 0.21 / 0.42 %; D calls move by ≤ 0.03 pp.
    - *Ties across predecessors depended on evaluation order* (a predecessor
      overlapping the entry beat an earlier one ending before it). They are
      now compared explicitly: path start, then node order. The documented
      rule was also wrong about junction bases that score the same on either
      side: the later segment gets them (ALGORITHM §4 now says so).
    - *Posterior range.* A gap between non-overlapping neighbours was a
      linear factor and underflowed for long costly gaps; it now joins the
      log-scale. An overlap of 540+ matched bases gave marginals above 1;
      such reads are now refused (`OVERLAP_SPREAD`). Neither occurs on
      immunoglobulin reads; both have tests that fail on the old code.
    - *D floor keyed on the setting, not the adjudicator.* With `--residual`
      (posteriors on by default) the floor was dropped although the residual
      adjudicator has no floor of its own, so it could call sub-floor Ds. The
      floor now follows whether the DAG computes posteriors, and
      `--residual` turns them off.
    - *Reported event.* A placement reaching 1–3 bases into the call's
      region counted in full, since a few bases occur in almost any allele.
      A placement now counts only if it covers more than half of the called
      span. D-allele ECE fell on every SHM preset (real mid 1.45 → 1.29 %,
      hard 1.47 → 1.16 %; synthetic 2.63 → 2.44 %, 2.94 → 2.43 %).
    - *3′-partial rule.* An allele was declared partial whenever its
      relative's anchor mapped past its end, which also happens when a
      deletion stops the sibling alignment short. It is now partial only if
      it ends while still aligned to the relative (within a codon); otherwise
      the weak motif hit stands, flagged. Partial alleles also keep the
      relative's FR1–FR3 starts, so their reads still get FR1–CDR2 columns
      (26 human IGH alleles). Anchors on the human set are unchanged.
    - *Docs.* D-error numbers had been measured on the previous commit;
      "whole-read J changes no call" no longer held; the item 11 account of
      the old oracle was wrong. All are corrected below with measurements of
      this code.
18. **Gene margin from sampling noise; every k-mer; short seeds by count**
    (`retrieve/budget.rs`, `retrieve/seed_vote.rs`). The V gene budget kept
    genes within `max(8 votes, 5 % of the leader)`; 8 votes are two 9-mer
    seeds whatever the read. A seed count `n` fluctuates by about `√n`, so
    two genes near the leader's `L` votes (4 per seed) differ by chance by
    about `√(2·4·L)`. The margin is now `z` such standard deviations:
    relative to the leader, but proportionally wider for weakly supported,
    mutated reads. The 7-mer trigger (16 votes, i.e. 4 seeds) became a seed
    count: 7-mers vote while the leader holds fewer than 25 nine-mer seeds,
    where the count's relative error `1/√n` exceeds 20 %. Measured jointly
    with the seed sampling on human IGH (hard and mid presets, fixed 20 and
    30 % SHM with two seeds each, 200-nt 3′ fragments) and the synthetic
    database, V gene in %:

    | Setting | Hard | 30 % SHM | Mean of 7 conditions |
    |---|---|---|---|
    | previous release: closed syncmers, fixed margin, 4-seed trigger | 98.36 | 75.70 | 85.16 |
    | the same after the seed audit's new order (444d966) | 98.06 | 75.37 | 85.24 |
    | closed syncmers, proportional margin 10 / 30 % of `L`, 7-mers always | 98.06 / 98.40 | 70.88 / 75.98 | 83.44 / 85.25 |
    | every k-mer, `z = 1` / `1.5` / `2` / `3`, 7-mers always* | 98.38–98.40 | 82.44–84.72 | 87.29–88.19 |
    | **every k-mer, `z = 2`, 25-seed trigger (default)** | **98.50** | **84.65** | **88.12** |
    | closed syncmers, `z = 2`, 25-seed trigger | 98.48 | 80.25 | 87.05 |
    | open syncmers (`w = 2`), `z = 2`, 7-mers always | 98.38 | 80.27 | 86.95 |
    | every k-mer, `z = 2`, 7-mers always | 98.40 | 84.70 | 88.17 |

    \*At 30 % SHM, 7-mers always equals the 25-seed trigger (84.70 vs
    84.65); on the hard preset it is 0.1 point lower. Against the previous
    release the default gains 0.14 points of V gene on the hard preset, 9.0
    at 30 % SHM and 3.0 on the mean of seven conditions.

    With closed syncmers a proportional margin must reach 50 % of `L` to
    match the fixed one at 30 % SHM (4-seed trigger: 72.47 at 30 %, 76.19
    at 50 %, against 75.37); with every k-mer 15–30 % already matches it,
    but even 50 % stays below `z = 2` (83.98 against 84.70, 7-mers always).
    Above `z = 2` the 16-gene cap binds: `z = 3` adds at most 0.1 point
    (synthetic, 30 % SHM), and a cap of 32 adds one point there only.
    Accuracy saturates between 25 and 50 seeds; always adding 7-mers is no
    better and costs 13 % more instructions per read than the 25-seed
    trigger on the hard preset (0.3 % at 30 % SHM, where the trigger fires
    on most reads anyway).

    *Sampling reversed.* Under the noise-based margin, closed syncmers
    trail every k-mer by 1.1 points of mean V gene (4.4 at 30 % SHM) for
    0.4 % fewer instructions on the hard preset, and open syncmers' edge
    (item 15, §8) disappears: the syncmers' parity with every k-mer, and
    open syncmers' lead, came from how sparser seeds interacted with the
    fixed margin and trigger. Every k-mer is the default; `SeedSampling`
    stays as an option. Minimisers were not re-measured under the new
    margin.

    Cost against the previous defaults: +4.1 % instructions per read on the
    hard preset, +6.1 % at 30 % SHM. Position accuracy: mean V start/end
    error on the hard preset is 0.55 / 0.93 bases (0.55–0.59 / 0.88–1.02
    across every setting swept). On the hard preset D, J and D-free false
    positives did not change; at 30 % SHM D and J rose by 0.7–1.3 points
    and D-free false positives by 1.3–2.3 (more reads reach a V call). V
    calibration improved (hard ECE 1.31 → 1.13 % against the previous
    release).

## 4. Results

Human IGH, 5000 reads per preset, germline-equivalent accuracy in %:

| Preset | Orientation | V | V gene | D | J | Junction | Random reads called |
|---|---|---|---|---|---|---|---|
| easy | 99.94 → 100 | 99.54 → 99.84 | 99.92 → 99.98 | 98.78 → 99.29 | 99.90 → 99.94 | 99.94 → 100 | 30.5 → 0 |
| mid | 100 → 100 | 97.76 → 97.94 | 99.84 → 99.84 | 91.96 → 93.00 | 98.84 → 98.80 | 99.93 → 99.97 | 30.2 → 0 |
| hard | 99.98 → 99.98 | 89.74 → 89.96 | 98.44 → 98.50 | 81.72 → 83.58 | 92.70 → 93.32 | 98.39 → 99.55 | 31.2 → 0 |

Junctions "after" are scored with the corrected anchors (item 16), which
remove the reads of 3′-partial alleles from the denominator (hard: 2416 →
2214 reads); with the old annotation the same calls scored 98.47 %.

Throughput is 849 → ~1900 reads/s per thread without posteriors and ~1300
with them, and ~3600 → ~7300 / ~5400 reads/s on 4 threads (a busier
container than the earlier ~2200 / ~9000; in instructions the exact DAG
costs 3 % and posteriors 18 %). Item 18 adds 4 % of instructions on the
hard preset; wall clock moved within the container's noise.

Calibration on real IGH (ECE %, easy / mid / hard), with the event of items
14 and 17:

| | V | D | J | D presence |
|---|---|---|---|---|
| ECE | 0.04 / 0.29 / 1.13 | 0.49 / 1.29 / 1.16 | 0.07 / 0.30 / 0.69 | 0.20 / 1.98 / 1.24 |
| Brier | 0.001 / 0.009 / 0.043 | 0.002 / 0.040 / 0.072 | 0.000 / 0.006 / 0.031 | 0.003 / 0.026 / 0.024 |

On the synthetic database, V and J ECE stay ≤ 1 % and D ≤ 2.5 %; D presence
is the weakest (5.0 % mid, 3.2 % hard), because the synthetic D genes are
shorter and more alike than human ones, which the derived `P_D` does not
model (§5.2). Items 13, 14 and 17 together improved D-allele ECE on every
preset with SHM (real 1.68 → 1.29 and 1.80 → 1.16 %, synthetic 3.45 → 2.44
and 4.21 → 2.43 %).

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

1. **V at extreme SHM.** The hard-preset V gene shortfall against the
   starting code (98.44 %; the previous release had 98.36 %) is gone: the
   noise-based margin (item 18) gives 98.50 %. At a fixed 30 % SHM, V gene
   is 85.2 % and V allele 53.8 %; 3.3 % of reads are misoriented and 3.9 %
   get no V or J call. The V posterior is overconfident
   there (ECE 16 %): when retrieval misses the true gene, no candidate can
   carry its mass. The 16-gene cap binds at that divergence; raising it
   helped only the synthetic database.
2. **D presence on D-free junctions.** 13–26 % of reads with no D receive one
   (19 % on human mid reads, 4 pp fewer since calls use the per-read model,
   §9.8), because a 9-nt chance match is indistinguishable from a 9-nt D
   remnant.
   `P_D` sets the trade-off, and the prior odds are the single knob.
   Precision/recall curves, not one number, belong in the paper; the
   reported `d_presence_posterior` now carries this uncertainty per read.
   The derivation undercounts the search: `ln(N·L)` counts genes × read
   positions, but a local alignment tries every allele, allele offset and
   read position (`K·m·n_D = e^{9.0}` against `e^{6.6}` for human IGH), and
   even prior odds understate IgSim's D inclusion; the two errors cancel.
   The search's own count is now a switch, `--d-presence-rate` (§9.9). It
   cannot replace the default: no presence rate follows from the search,
   at the simulator's rate exact 5-mers win (no-D −14 to −26 pp), and every
   constant cost moves along the same curve of D accuracy against D-free
   specificity. What the residual miscalibration tracks is junction length,
   a matter for an N/P length prior.
3. ~~**Calls and posteriors use different evidence scales.**~~ Resolved
   (§9.8): one integer model at the read's divergence drives both, so the
   call is a mode of its posterior. What remains is the gap between a mode
   and a marginal. Taking D presence from the posterior gives a D to a
   third of D-free reads (−35 pp on 20,000 reads), so presence stays the
   call's. Taking only the D allele from it, with a conditioned Viterbi pass
   for the reported alignment, is now a switch (`--d-allele max-marginal`):
   on top of per-read calls it gains +0.25 pp on human hard reads and is
   within noise elsewhere, so it is off by default.
4. **One alignment per D allele and window.** A D's optimal local alignment
   in its window is not always its optimal clipped alignment at the chosen
   cuts; the 12-base reclaim keeps this rare. It is the one remaining
   approximation of the segmentation model.
5. **Constants inherited from the note** — now measured (§7): five are
   inert on full-length simulated reads and default to neutral values; the
   rest are speed knobs or have a measured optimum at the inherited value.
   The V/J column and retention floors are not inert on truncated,
   chimeric or cross-locus reads, and since 0.39.4 they are back as part
   of the plausibility guards (§7.1). `tandem` (12) was not swept:
   VDDJ reads are 0.1 % of IgSim reads.
6. **E-values are approximate.** They use ungapped λ/K (anti-conservative for
   gapped scores) and the total database length (conservative for
   near-identical alleles). They gate hits; they are not p-values.
7. **Simulator realism.** IgSim has no AID hotspot targeting, no
   codon-length indels and no allele-frequency skew, so hotspot-aware scoring
   cannot be judged on it. Validation on real data with an independent caller
   (IgBLAST, partis) is still needed.
8. ~~**SIMD portability.**~~ Resolved (§9.4): one kernel over a lane
   abstraction runs AVX-512BW, AVX2, NEON and portable lanes.
9. **D errors are mostly the wrong allele, not a missing D.** With a remnant
   of ≥ 8 nt (§8), 2.7 % (mid) and 10.2 % (hard) of reads get a wrong D
   allele although a placement of the right one was found; the search loses
   0.2 % and 1.8 %. These are model errors: short, mutated D remnants that
   another allele explains as well or better.
10. **Posterior cost.** Measured again in §9.6: posteriors cost 30 % of all
    instructions (10.1 G vs 7.0 G on 2000 hard reads), not 18 %, because
    they also keep every admissible D in the lattice. No slot is negligible
    (every admissible D scores ≥ 8 against `P_D` 11), so pruning would be
    an approximation; a single-precision pass is the remaining option.
11. **Anchors without a relative.** Read without its gap columns, the
    human IGH set leaves seven V alleles (pseudogenes and ORFs) with a
    low-confidence motif anchor. One of them, IGHV1-38-4\*01 (Cys104 → Tyr),
    is placed on a TGT in another frame. Without a relative there is no
    evidence to place it better, so it is reported as low-confidence rather
    than silently trusted. IMGT-gapped input numbers all of them exactly
    (§10).

## 6. Plan

In order of expected value:

1. ~~**Posterior decoding**~~ — done (§3, item 12), and calls from the
   per-read model (§9.8). The D prior from the D search's own statistics is
   a switch (§9.9). Next on that line: a prior on N/P lengths, which the
   D-presence residual points to (short junctions are D-free more often than
   the model allows).
2. ~~**Two cuts per overlap with N/P between**~~ — done (§3, item 11).
3. ~~**Sensitivity analysis**~~ — done (§7).
4. **Cross-tool benchmark.** IgBLAST done (§8). partis and Swig remain, at
   several SHM levels, with precision/recall curves for D.
5. ~~**AVX-512BW (32 lanes) and NEON**~~ — done (§9.4), with a portable
   backend for every other target.
6. ~~**Allele rescoring along the gene frame.**~~ — done (§9.6): siblings
   are bracketed by an upper bound from their frame representative's score
   and a lower bound from its traceback, and DP runs only where the bounds
   straddle the top-N cutoff, for V and J. Calls are byte-identical.
7. **Calibrated statistics.** Fit gapped λ/K and an effective database size
   by simulation.

## 7. Sensitivity analysis of the inherited constants

One constant at a time around the defaults, on real human IGH at the mid and
hard presets (3000 reads each, plus 1000 random reads); `examples/benchmark`
takes each constant as a flag. "Flat" means every accuracy metric (V, V
gene, D, J, junction, D-free false positives, random-read calls) and every
calibration error was unchanged.

| Constant | Inherited value | Range tested | Effect | Now |
|---|---|---|---|---|
| V…J pair bonus | 24 | 0, 12, 36, 48 | flat | 0 |
| V acceptance | 24 columns, 60 % | 12–48 columns; 50–70 % | flat | none in the stage; 24 columns under the pair guards since 0.39.4 (§7.1) |
| J acceptance | 10 columns, 60 % | 6–18 columns; 50–70 % | flat | none in the stage; 10 columns under the pair guards since 0.39.4 (§7.1) |
| Distributed-D gate | on (16/14/80 %/19) | off; looser; tighter | flat | off |
| Retained bases V / J | 24 / 10 | V 1–36, J 1–14 | flat | 1 / 1; 24 / 10 under the pair guards since 0.39.4 (§7.1) |
| Retained bases D | 4 | 3, 5 | D-presence ECE ±0.1 pp | 4 (= `min_d_match`) |
| Seed lengths | 9 / 7 | 11/7, 11/9, 9/–, 8/6, 7/– | V ≤ ±0.07 pp | 9 / 7 |
| `max_occurrences` | 256 | 64–1024 | flat | 256 |
| Band padding | 8 | 4, 12 | flat (speed only) | 8 |
| Budget base genes | 3 | 2, 4, 5, 8 | hard V −0.20 / +0.10 / +0.20 / +0.20 pp; slower | 3 (knob) |
| Budget margin | 5 % / 8 votes | 2–10 %; 4–16 | ≤ 0.07 pp on the presets; 4–9 pp of V gene at 30 % SHM (item 18) | `z·√(2·4·L)`, `z = 2` |
| 7-mer trigger | 16 votes (4 seeds) | 4–50 seeds, always | presets flat; +4.8 pp V gene at 30 % SHM from 4 to 25 seeds (item 18) | 25 seeds |
| `min_d_match` | 4 | 3, 5, 6 | D −0.03 to −0.92 pp | 4 (optimum) |
| `P_D` | 11 (derived) | 5–17 | D accuracy vs D-free false positives (mid: 94.2 % / 61 % at 5, 93.2 % / 23 % at 11, 87.2 % / 0 % at 17) | 11 |

Two results matter for the paper:

- **Five inherited mechanisms were redundant** given the E-value gates and the
  D floor, and now default to neutral values (the fields remain for profiles
  and experiments). The 60 % identity floor is no accident: it is the
  break-even identity of `(2, −3)` (`2m − 3x > 0 ⇔ m/(m + x) > 3/5`), so a
  positive local score already implies it. All five together changed no call
  on real or synthetic databases at any preset, left calibration unchanged
  (within 0.03 pp), and kept random-read calls at 0.04 % of 20,000. These
  were full-length simulated reads only; §7.1 qualifies the conclusion.
- **The derived `P_D` is the calibration optimum.** D-presence calibration
  error is lowest at `P_D = 11` on both presets (mid 2.25 %, hard 1.66 %; 2.71 %
  and 2.59 % at 9, 3.36 % and 5.37 % at 13), so the posterior independently
  supports the value the derivation gives. Accuracy alone cannot choose it:
  D accuracy and D-free false positives trade monotonically.

### 7.1 Note (0.39.4): the floors on reads the simulator does not make

The sweeps above, and the ablation that called the V/J floors redundant,
used IgSim reads only: full-length rearrangements between random flanks.
On those, a positive local score implies 60 % identity and the E-value
gates imply a length, and the floors changed no call. They did not test reads whose V or J is truncated,
chimeric reads (two templates joined) or reads that cross loci, and there
the floors and SwiftIG's pairing rule decide whether a read becomes a
rearrangement. Swig 0.39.0 dropped them with the rest of SwiftIG; 0.39.4
restores them as the plausibility guards (ALGORITHM §4c, on by default,
`--no-pair-guards`): V and J hits of at least 24 and 10 local columns,
clipped V and J keeping at least 24 and 10 read bases, SwiftIG's
`supports_rearranged_pair` with its 225-base junction span, a D only
between a V and a J (since 0.39.5 also at a read's open end where it pays
the open-flank cost of ALGORITHM §2.3), and junctions only from aligned
anchors at most 231 bases apart. The 60 % identity floor is Swig Web's and swig-cli's
(`--min-identity 0.6`), and under the guards it also judges the part of a
V or J that a call keeps.

Against 0.39.3, on the human IMGT references, 20,000 simulated BCR reads
(IGH, IGK, IGL) per preset:

| | mid 0.39.3 | mid 0.39.4 | hard 0.39.3 | hard 0.39.4 |
|---|---|---|---|---|
| V / D / J allele | 98.70 / 93.46 / 99.11 % | same | 93.77 / 85.88 / 94.92 % | 93.77 / 85.83 / 94.91 % |
| reads whose V, D or J call or coordinates changed | | 0 | | 9 |
| templated junctions reported exactly | 99.97 % | 99.13 % | 99.63 % | 95.67 % |
| reads without a templated junction given one | 90.9 % | 24.1 % | 89.2 % | 14.0 % |
| `productive = T` without a templated junction | 1,122 | 215 | 673 | 90 |

On constructed reads, 0.39.3 paired the first 70 bases of a V with a J 20
random bases later, made a V and a J 500 random bases apart into a
rearrangement with a 543-nt junction, called 8 J bases after 12 random ones
a J on 18 of 20 reads, and gave a V followed by 100 random bases a D on 19
of 20 reads (D-presence posterior 0.91–0.99); `tests/guards.rs` checks each
case. On 1,000 mid IGH reads cut 30 or 100 random bases after their V, a D
was called on 70.0 and 92.5 % of reads and is on 0.3 and 0.5 %; behind 30
or 100 random bases before their J, a D on 81.0 and 95.9 % and a V on 7 of
2,000, and now on none. CPU time is unchanged within ±2 % (10,000 reads,
one thread). The costs: the templated junctions whose anchor codon is so
mutated that the alignment stops before it (0.84 points on mid reads,
3.95 on hard ones), and the D of reads that hold one but start after their
V, which SwiftIG did not call either (on 952 mid IGH reads with a D, cut
where their V ends, 878 right D calls in 0.39.3 and none now; with the
last 30, 40 or 60 V bases kept, 805, 866 and 880 against 881, 881 and
880).

So the conclusion above holds for full-length reads, and it does not
license dropping a check whose purpose is reads the simulator does not
make. Two gaps remain:

- **Short-V cross-locus chimeras.** J is searched only in the locus of the
  read's anchoring V, so the 30–60 bases of an IGHV that a chimera carries
  before an IGK J find a chance IGHJ, which SwiftIG, searching every locus,
  never paired. The guards add a rule SwiftIG did not need (a V under 80
  columns pairs only with a J significant on its own), but an IGHJ that
  aligns to the IGK J well enough to be significant passes it: of 168
  such constructed reads, 10 are called productive (0.39.3: 30; SwiftIG:
  none).
- **J-mode CHMMAIRRa.** Swig's PCR-chimera filter now extends the V
  alignment gaplessly toward its germline 3′ end, and the J toward its 5′
  start (a J by a run scoring at least 2), before it threads a read
  ([07_CHIMERA_INFERENCE.md](../../public/methods/07_CHIMERA_INFERENCE.md)). That
  restores V mode (420 of 500 and 1,698 of 2,000 simulated IGH V chimeras
  flagged; SwiftIG 420 and 1,696; 0.39.3 416 and 1,676), but J mode flags
  35 of 60 and 130 of 200 J chimeras against SwiftIG's 38 and 138 (0.39.3:
  35 and 125): the engine prefers the 3′ parent's J with a hard 5′ clip,
  so the 5′ parent's bases are not in the J alignment the filter reads.

## 8. Seed sampling, the J/D interplay, and IgBLAST

**Would minimisers or syncmers lose sensitivity?** Seeds are sampled alike in
the gene-frame index and in the read (`seed::SeedSampling`). The literature's
recommendations, read from the primary sources:

- Edgar (2021, PeerJ 9:e10805): closed syncmers (a k-mer is kept when the
  first smallest of its `s`-mers is its first or last one) beat minimisers
  at any density they can reach; open syncmers with an offset beat
  minimisers at lower densities. Closed syncmers are context-free (a seed
  shared by read and germline is kept in both or neither) and have a
  *window guarantee*: any `w = k − s` consecutive k-mers hold one, so every
  unmutated stretch of `k + w − 1` bases yields a shared seed (12 for
  9-mers at `w = 4`; on IGHV, 0 seedless 12-base windows).
- Shaw & Yu (2022, Bioinformatics 38:4659, Theorem 8): among open
  syncmers, the smallest `s`-mer in the middle maximises conservation. Open
  syncmers keep context-freedom but lose the window guarantee (5 % of
  12-base IGHV windows seedless at `w = 2`, 24 % at `w = 4`).
- Dutta et al. (2022) and strobealign (Sahlin 2022) use other shapes for
  long or short DNA reads; minimisers keep the window guarantee but are not
  context-free (9–11 % of unmutated k-mers selected on one side only).

Our closed syncmers matched Edgar's definition and tie rule; an audit found
the order fixed poly-A s-mers as smallest, spaced shapes ranked the raw
sequence instead of the key, `w` was clamped silently, and the docs stated a
loose window bound (13 bases). All are fixed and tested (`seed.rs`). On a
spaced shape, selection by key keeps read and germline consistent but has no
window guarantee (neighbouring keys share no s-mers); the default shapes are
contiguous.

With the old fixed gene margin, closed syncmers (`w = 4`, 40 % of k-mers)
matched every k-mer on the presets and edged it at high SHM (fixed 30 % SHM,
V gene 76.35 vs 74.90 %), for 2.3 % fewer instructions, so they became the
default (item 15); open syncmers (`w = 2`) measured better still. That
parity was fragile: the audit's new order alone, which changes which 40 % of
k-mers are kept, moved closed syncmers on the hard preset from 98.36 to
98.06 % V gene (V ECE 1.31 → 1.58 %) while every k-mer stayed at 98.24 %,
and it came from how sparser seeds interacted with the fixed margin and vote
trigger. Under the
noise-based margin (item 18) every k-mer is better: mean V gene over seven
conditions 88.12 vs 87.05 % (closed `w = 4`) and 86.95 % (open `w = 2`,
7-mers always), 84.65 vs 80.25 % at 30 % SHM, equal on the presets, for
0.4 % more instructions on the hard preset. Every k-mer is the default
again; the sampled schemes remain options. The window guarantee does not
affect boundaries: seeds only choose candidates and the alignment band
(cluster ± 8 diagonals), and mean V boundary errors were 0.55–0.59 /
0.88–1.02 bases on the hard preset for every scheme. It could matter only
for an indel longer than 8 bases with no shared seed beyond it, which the
simulator does not produce.

**Does searching J downstream of V let a misplaced J eat the D?** The J is
searched in `[min V end − 12, read end)`, every allele unbanded, before the D.
The D is not placed jointly with the J, but three mechanisms keep the J from
taking its bases:

1. the D windows reach 12 bases back into the J (`reclaim`), and the
   conserved-anchor window `[V Cys104 + 3, J Trp/Phe118)` is added, which
   depends on where the J's conserved codon lands, not on how far its
   alignment reaches to the left;
2. adjudication is joint: the DAG gives each placement its own exit and
   entry cut and orders segments by kind, so a D can precede a J whose
   alignment reaches back over it (item 17), any overlap between them goes
   to whichever explains it better, and bases neither explains go to N/P;
3. the J is still an alignment, not a call: the top three J alignments
   (E ≤ 1, plus ties) are kept, and the path chooses among them, with or
   without a D.

`examples/d_errors.rs` attributes every D error on simulated reads. Human
IGH, 5000 reads, D remnant ≥ 8 nt:

| | mid | hard |
|---|---|---|
| D correct | 97.08 % | 87.52 % |
| wrong D, a correct placement existed (adjudication) | 2.66 % | 10.22 % |
| wrong D, not found | 0.23 % | 1.77 % |
| true D outside the D windows | 0.00 % | 0.00 % |
| missing D, the called J covers its bases | 0.00 % | 0.00 % |
| missing D, the called V covers its bases | 0.00 % | 0.02 % |
| called J covers ≥ 3 bases of the true D (any call) | 0.21 % | 0.42 % |

Counterfactuals on the hard preset: without the anchor window no true D
falls outside the windows (D correct 87.43 %); with J searched over the
whole read, 0.05 % do, with or without the anchor window. Whole-read J
search is 18 % slower; it changes no call on mid and a few on hard (D called
−0.04 pp, J −0.02 pp, one more D-free read given a D). The failure mode
seen before is not reproduced. Since the starting commit, J moved from seed
retrieval plus every allele over the whole read to an exhaustive search
right of V, the single shared cut became two independent cuts, and a D can
now precede a J that starts earlier (item 17); the anchor window was already
there. Whole-read J is the only configuration measured that lets a true D
fall outside the windows. So the joint search is kept where it matters, in
adjudication, and J retrieval right of V costs nothing measurable. The D
errors that remain are wrong alleles (§5.9).

**Against IgBLAST 1.22.** Same reads (5000 per preset, forward strand), same
germline as BLAST databases, `-organism human -ig_seqtype Ig -outfmt 19`,
scored by `swig evaluate`:

| Preset | Tool | V | V gene | D | J | Junction | D called on D-free reads | Reads/s, 4 threads (previous V retrieval) |
|---|---|---|---|---|---|---|---|---|
| easy | swig | 99.88 | 100.00 | 99.29 | 99.84 | 100.00 | 0 % | 5200 |
| easy | IgBLAST | 99.12 | 99.78 | 98.86 | 99.86 | 99.77 | 0 % | 67 |
| mid | swig | 97.88 | 99.72 | 93.35 | 98.74 | 99.97 | 23 % | 5100 |
| mid | IgBLAST | 97.18 | 99.62 | 93.73 | 98.60 | 99.72 | 43 % | 66 |
| hard | swig | 89.44 | 98.84 | 83.00 | 92.80 | 99.69 | 16 % | 5200 |
| hard | IgBLAST | 89.16 | 98.90 | 81.75 | 91.84 | 98.07 | 30 % | 71 |

IgBLAST's higher mid D accuracy comes with twice the D calls on D-free
junctions, and its V gene is 0.06 pp better on the hard preset. Two setup
issues had to be fixed for a fair junction comparison. IgBLAST takes J
anchors from its auxiliary file, which lacks three of this database's J
alleles, so it reported no junction for 18 % of reads; the entries were
added, with frames and CDR3 ends from `swig prepare`, which agrees
with the file on all 13 shared alleles. The remaining disagreements exposed
item 16 on our side.

## 9. Engineering review (edition 2024, Rust 1.97)

A line-by-line review of the whole crate for soundness, duplication, dead
abstractions and portability, with every claim re-measured at a scale where
the comparisons are statistically meaningful.

### 9.1 How it was checked

- **20,000 reads per preset**, not 2000–5000. At 5000 reads a 93 % rate has
  a 95 % interval of ±0.7 pp, so most sub-point differences in §§3–8 were
  within noise. Every comparison below is **paired** (both configurations
  on the same reads) and tested with an exact McNemar test on the reads only
  one of them gets right (`swig evaluate --truth T a.tsv b.tsv`,
  `sim::Paired`).
- **Real germline**: the IMGT human IGH reference shipped in the `airrship`
  PyPI package (424 V including ORFs and pseudogenes, 44 D, 13 J; 415 / 37 /
  13 after merging identical sequences), read directly from IMGT headers.
- **Regression oracle**: after every refactoring step, the AIRR output of
  60,000 reads (20,000 per preset) was compared byte for byte with the
  output of the starting commit. Every refactoring below left it identical.
- **Debug build** (overflow checks and the clip-index rescoring assertion)
  over 40,000 real-germline reads with the tandem-D grammar; **callgrind**
  instruction counts for speed; **qemu** for AArch64.

### 9.2 Defects found and fixed

1. **Sibling bounds were applied across bands.** `S*(A) ≤ S*(R) + Δ(R, A)`
   maps A's optimal path onto R cell for cell, so it holds only when both
   were aligned in the same band. The V stage applied it whenever lengths
   matched, including siblings in another gene frame or at another frame
   offset, whose bands differ. It now bounds same-band siblings only; Δ is
   precomputed per same-length allele pair (it was recomputed per read), and
   each sibling finds its representative's score in O(1) (it scanned every
   representative of the chunk). A property test checks the bound for random
   queries and every band both alleles share.
2. **The D keep floor ignored its configuration.** `lossless_min_score`
   ended in `floor.max(keep.min_score.min(floor))`, which is always `floor`,
   and its doc described the opposite rule. It is now the configured floor
   raised to the lowest score an admissible hit can have.
3. **IMGT headers broke gene grouping.** The allele name was the first
   whitespace token, i.e. the whole `accession|allele|species|…` string, so
   every IMGT allele became its own gene (gene frames, the gene budget and
   sibling bounds silently disabled) and calls could not be evaluated
   (0 % "accuracy" by name). The name is now the field that names an allele;
   raw IMGT files give byte-identical calls to cleaned ones.
4. **An unsound extension point.** `Refiner` ran after the posterior
   marginals were computed, so any refiner that rewrote the path would have
   reported probabilities of a different call; its only implementation was a
   no-op. Removed (one fewer type parameter on `Assigner`); the outer-flank
   extension it shared a module with is now `flank`.
5. **Wrong documentation**: Karlin–Altschul statistics were described as
   reporting-only (they gate every call); a BMI2 `pext` seed path that does
   not exist; `HotspotScoring` described as part of the J profile; the cut
   recursion as O(log n) (it is O(1) amortised); the posterior cost (§9.6).

No crash, overflow or assertion failure was found: degenerate reads (empty,
1 base, all-N, IUPAC, low-complexity, 20 kb, chimeric, every prefix and
suffix of a read) pass every profile (`tests/robustness.rs`), and the 25
functional IMGT alleles the annotator leaves without an anchor are exactly
the ones IMGT flags "partial in 3′".

### 9.3 Design

| Before | After |
|---|---|
| `dag.rs` and `posterior.rs` each with their own node order, layout, slot groups and incoming/outgoing iterators | one `adjudicate/lattice.rs`; each pass keeps only its arithmetic (−262 lines) |
| AVX2 intrinsics in one function; every other CPU scalar | one kernel over `align::lanes::Lanes<L>` (§9.4) |
| five hand-copied `with_*` builders | one field-moving macro |
| three examples with their own option parsing and database setup | `examples/common` |
| six private random generators in tests | `sim::Rng` everywhere |
| stage 1 retrieved V again after a non-decisive orientation | the chosen strand's candidates are kept |
| database length summed per reported segment | stored once in `GermlineSet` |
| traceback buffers zero-filled per group | grow-only (−2.2 % instructions) |

### 9.4 Portability

The inter-sequence kernel is written once, generic over `Lanes<L>` — a
register of `L` 16-bit lanes with an associated mask type, so AVX-512's `k`
registers and AVX2/NEON vector masks both map one-to-one onto
instructions — and compiled per backend inside a `#[target_feature]` entry
point. `Backend::detect` picks one at run time; `--simd` overrides it.

| Backend | Lanes | Single thread, 20,000 hard reads | Calls |
|---|---|---|---|
| AVX-512BW (`--simd avx512`) | 32 | 2403 reads/s (2339 after PR #2) | identical |
| AVX2 (default on x86-64) | 16 | 2103 reads/s (1869 after PR #2) | identical |
| NEON (default on AArch64) | 16 (two `q` registers) | — (tested under qemu) | identical to x86-64 |
| portable (any target; compiler-vectorised) | 16 | 1621 reads/s (1498 after PR #2) | identical |
| scalar | 1 | 810 reads/s (728 after PR #2) | identical |

Every backend equals the scalar kernel cell for cell on random batches up to
five registers wide and on degenerate jobs; the full test suite passes on
x86-64 and on AArch64; AArch64 output equals x86-64 output byte for byte.
AVX-512 is opt-in: +14 % on one thread but ~+4 % on four (within noise in
this container), and 512-bit instructions lower the clock on many CPUs. The
generic kernel costs the same instructions as the hand-written AVX2 one.

### 9.5 Results at 20,000 reads, and what the ablations really show

Human IGH (IMGT, `airrship`), 20,000 reads per preset, half
reverse-complemented, germline-equivalent accuracy in %:

| Preset | Oriented | V | V gene | D | J | Junction | D-free reads given a D | Random reads called |
|---|---|---|---|---|---|---|---|---|
| easy | 100.00 | 99.88 | 100.00 | 99.69 | 99.94 | 100.00 | 0.0 % | 0.04 % |
| mid | 100.00 | 98.61 | 99.97 | 93.11 | 99.11 | 99.97 | 24.5 % | 0.04 % |
| hard | 99.98 | 93.23 | 99.61 | 85.91 | 94.25 | 99.65 | 16.0 % | 0.04 % |

These rows include PR #2 (every k-mer as a V seed, noise-based gene margin),
merged into this branch; paired against the code before it on the same
reads, it gains +0.17 pp hard V (p 1e−9) and changes nothing else
significantly. Calibration (ECE %, V / D / J / D presence): easy 0.06 /
0.30 / 0.04 / 0.14, mid 0.09 / 1.19 / 0.18 / 1.96, hard 0.32 / 1.41 / 0.18 /
1.10. On the
synthetic database D-presence ECE is 7.3 % (mid) and 4.1 % (hard) at this
sample size, the clearest sign that the D prior needs §5.2's statistics.

Paired against the default on the same 20,000 reads (Δ in pp, exact
McNemar p; "no D" is D-free reads left without a D call):

| Change | Preset | D allele | no D | V | J | Junction |
|---|---|---|---|---|---|---|
| residual adjudicator | mid | +0.37 (p 7e−5) | −24.8 (p 5e−79) | −0.02 | −0.02 | 0 |
| | hard | +0.62 (p 3e−8) | −21.2 (p 1e−56) | 0 | −0.05 | −0.01 |
| `sensitive-d` | mid | +0.77 (p 4e−44) | −18.4 (p 2e−62) | 0 | 0 | 0 |
| | hard | +1.01 (p 3e−58) | −14.1 (p 1e−44) | 0 | 0 | 0 |
| `igblast` profile | hard | +0.26 (p 0.10) | −75.4 | −0.16 (p 9e−6) | −0.64 (p 7e−25) | −0.40 (p 6e−11) |
| `truth` profile | hard | −0.60 (p 1e−4) | — | −0.17 (p 2e−6) | −0.07 (p 1e−4) | −0.30 (p 1e−8) |
| no anchor D window | mid | +0.05 (p 0.013) | −0.36 (p 0.13) | 0 | 0 | 0 |
| | hard | −0.01 (p 1.0) | −0.19 (p 0.5) | 0 | 0 | 0 |
| tandem-D grammar | hard | +0.01 (p 0.69) | — | 0 | 0 | 0 |

Read with both columns, every gain in D accuracy is bought with D calls on
D-free reads: the residual baseline's +0.4–0.6 pp is not a better model but
a more liberal D threshold, and `sensitive-d` sits on the same curve. The
default is the only configuration that does not trade one for the other
blindly; where on the curve to sit is the prior odds of a D (§5.2). The
anchor D window and the tandem grammar change nothing measurable on IgSim
reads (and the anchor window costs nothing measurable either); both stay as
cheap insurance for reads the simulator does not produce.

### 9.6 Where the time goes

Callgrind, 2000 hard reads, one thread (10.1 G instructions, including
database preparation; 10.5 G after PR #2's every-k-mer seeds, of which the
posterior is 3.1 G):

| Stage | Share |
|---|---|
| V (retrieval, scoring, traceback) | 21 % |
| J (every allele right of V) | 15 % |
| D (prefilter, alignment) | 14 % |
| adjudication: Viterbi + posterior | 35 % (posterior alone 30 % of the run) |
| database preparation, I/O, reporting | the rest |

Remaining opportunities, in order of expected gain: (1) 8-bit lanes for D
jobs — D scores fit in `i8`, so AVX2 would run 32 D jobs per register;
(2) a single-precision posterior pass; (3) banded J from J seeds when the V
is unmutated (exhaustive J is the safe default under SHM); (4) overlapping
input parsing with assignment.

*Allele rescoring along the gene frame (plan item 6), done.* Every sibling
of a frame representative is bracketed, `L ≤ S*(A) ≤ S*(R) + Δ(R, A)`
(`family_bound`): the upper bound now holds across a gene frame, not only
at equal length, and the lower bound rescores the representative's
traceback along the sibling. The top-N cutoff is taken over the lower
bounds, so it prunes more, siblings whose bounds meet need no DP, and those
that will be traced are aligned once instead of scored and then aligned.
J gets the same bounds. Instructions on 2000 reads, one thread, each step
byte-identical to the last on the 20,000-read sets:

| Step | hard | mid |
|---|---|---|
| main | 10,549 M | 10,967 M |
| cutoff only under permissive acceptance, explicit tie key | 10,556 M | 10,976 M |
| Δ across the gene frame | 10,488 M | 10,966 M |
| traced lower bounds, cutoff over them, traceback reuse | 10,064 M | 10,550 M |
| known and early-traced siblings | 9,953 M | 10,436 M |
| J bounds | 9,708 M | 10,165 M |
| bookkeeping | 9,705 M | 10,162 M |
| strongest representative aligned, not scored | 9,655 M | 10,112 M |
| one representative per frame | 9,627 M | 10,082 M |
| representative ties toward the frame's medoid | 9,606 M | 10,057 M |

In all −9.0 % (hard) and −8.3 % (mid): V and J DP fell from 4.04 G to
3.21 G instructions on hard, while the stage's own bookkeeping rose from
0.02 G to 0.05 G (frame paths 12 M, the table 6 M, lower bounds 3 M). Two
experiments
did not lower the count and were dropped: gene-frame placement for J (J's
unbanded bounds accept any placement), and raising the cutoff with the
early-traced siblings' exact scores before scoring the rest. The cutoff
applies only when acceptance has no length or identity floor: a rejected
alignment makes traceback refill from below the N-th score, where a pruned
sibling may lie (with V acceptance of 60 columns and 93 % identity, main's
bounds changed the allele set of 1 of 20,000 hard reads).

### 9.7 Roadmap status

Of the roadmap in the README, the wider-SIMD item is done (AVX-512BW, NEON
and a portable backend), and so are allele rescoring along the gene frame
(§9.6), calls from the per-read model (§9.8) and the annotation priors
(§10). The D prior from the D search's own statistics is implemented and
measured, and stays a switch (§9.9): the search's count is derived, the
presence rate is not. C-segment assignment is done too (§11).

### 9.8 Calls from the per-read model

The Viterbi call used the fixed `(2, −3)` scores and the posterior the
read's own log-odds, so the reported probabilities belonged to a model that
had not made the call. Now one integer model at the read's divergence `p̂`
drives both (ALGORITHM §4, *Values*): substitution columns `s·u_r + v_r`
(each kind's affine log-odds map quantised whole, 1/2¹⁶ nat), gaps and priors
as multiples of `u0 = round(2¹⁶·λ)`, and a D cost of `P*_D/m` matched D
columns. The fixed scale is the exact special case at 11.2 %, and stays as
`--calls fixed`. `p̂` is estimated before adjudication from the V and J
columns no path contests, so the call is one Viterbi pass. Exhaustive oracles
check the integer DAG (800 optima on mid and hard IgSim placements, 135 of
640 off the fixed path) and that its key is the maximum of the posterior's
`e^{units/F}`.

The fixed scale had hidden a second D penalty: it undervalues an exact
match at low divergence (2λ = 1.27 nats against `ln 4(1 − p̂)` = 1.38), so
an exact 5-mer tied `P_D` = 10 and lost on the "fewer D" rule. Read at the
read's scale, a constant D cost in nats lets those 5-mers in (a prototype
gave 12–22 pp more D-free reads a D); the cost `c_D(p̂)` keeps the fixed
scale's exact-run classes at every `p̂` instead.

Default (per-read call) against main, paired on 20,000 reads per preset (Δ
pp, exact McNemar p; D split by the true remnant's length, of which 1–3 nt
is below `min_d_match` = 4 and cannot be found):

| Set | V | D allele | D 1–3 nt | D 4–7 nt | D 8+ nt | no D | J | Junction |
|---|---|---|---|---|---|---|---|---|
| igh.easy | 0 | −0.01 (1.0) | — | +1.45 (1.0) | −0.01 (0.69) | −0.96 (0.5) | 0 | 0 |
| igh.mid | 0 (1.0) | −0.12 (0.10) | −6.51 (2e−8) | +0.48 (0.46) | +0.05 (0.31) | **+4.29** (2e−10) | +0.01 (0.61) | −0.01 (1.0) |
| igh.hard | +0.01 (0.75) | −0.01 (1.0) | −0.48 (0.63) | +0.07 (1.0) | +0.01 (1.0) | 0 (1.0) | +0.02 (0.45) | +0.01 (1.0) |
| syn.easy | +0.01 (1.0) | −0.01 (1.0) | — | — | −0.01 (1.0) | 0 | −0.01 (1.0) | 0 |
| syn.mid | 0 | −0.39 (4e−5) | −7.34 (1e−18) | +1.02 (0.024) | −0.09 (0.13) | **+3.99** (2e−15) | **+0.07** (0.003) | 0 |
| syn.hard | +0.02 (0.45) | −0.33 (1e−5) | −2.52 (1e−6) | −0.47 (0.12) | −0.12 (0.09) | +0.31 (0.44) | 0 (1.0) | +0.03 (0.5) |

(igh.easy's no-D change is 2 of 209 reads.) Calibration, main → default:

| Set | D presence Brier (paired t) | D presence ECE % | D allele ECE % | D allele Brier t |
|---|---|---|---|---|
| igh.easy | 0.00065 → 0.00035 (−6.2) | 0.15 → 0.10 | 0.31 → 0.32 | −1.1 |
| igh.mid | 0.02810 → 0.02462 (−13.3) | 2.23 → 1.36 | 1.02 → 0.66 | −3.6 |
| igh.hard | 0.02261 → 0.02168 (−5.3) | 0.99 → 0.47 | 0.98 → 0.87 | +1.3 |
| syn.easy | 0.00062 → 0.00029 (−6.4) | 0.16 → 0.12 | 0.58 → 0.59 | +0.8 |
| syn.mid | 0.05266 → 0.04547 (−23.7) | 4.89 → 3.27 | 2.43 → 1.67 | −4.5 |
| syn.hard | 0.04353 → 0.04079 (−12.8) | 2.99 → 1.82 | 1.68 → 1.50 | −1.9 |

V and J ECE move by at most 0.05 (human hard V 0.33 → 0.37, J 0.35 → 0.37).

The steps separately, each paired against the one before:

- *Unrounded `P*_D` in the posterior* (the call uses `⌊P*_D⌋`, which equals
  the old rounding on both databases): calls identical; D-presence Brier t
  −6.9 / −15.2 / −7.6 (human easy / mid / hard), −6.9 / −26.3 / −18.0
  (synthetic).
- *`p̂` before adjudication* (fixed-scale call): calls identical; D-presence
  Brier +0.00003 to +0.00011 (t +1.1 to +4.7), V / D / J ECE within 0.04.
- *Per-read call*: the tables above.

Why it is the default: no D-free read set gets more D calls (two sets get
4 pp fewer), D-presence Brier improves on every set, V, J and junction lose
nothing, and D allele accuracy on remnants of 4 nt or more is within noise
on every set. The overall D allele drop on the synthetic database (−0.39
and −0.33 pp) is significant only in the 1–3 nt stratum, remnants that
`min_d_match` = 4 makes undetectable, so a call counted correct there rests
on chance matches; the per-read call simply calls fewer of them. The one
doubtful number: on syn.hard the strata of 4 nt and more lose 29 reads
together (p 0.016), not significant at 0.01 but not clearly zero.

**Max-marginal D allele** (`--d-allele max-marginal`), paired against the
per-read call: D allele +0.05 (p 0.12), **+0.25** (6e−5), +0.04 (0.24),
+0.14 (0.017) on human mid / hard and synthetic mid / hard, all of it on
remnants of 8 nt or more; no D-free read changes; D allele Brier +0.00004 to
+0.00035 (t ≤ 2.2). On the fixed-scale call it gains +0.20 (p 7e−4) and
+0.27 (p 2e−4) on human mid and hard: the per-read call is already the mode
of the same posterior, which leaves less to gain. Significant on one set
only, it stays off by default.

Instructions on 2000 hard reads: main 10.550 G; per-read call (default)
10.659 G (+1.0 %); fixed-scale call 10.605 G; per-read call with
max-marginal D 10.713 G (+1.5 %).

### 9.9 A D prior from the D search's statistics

`P*_D = [ln(N·L) − ln π_odds]/λ` counts genes × positions at even odds.
`DPrior::Search` (`--d-presence-rate π`) counts what the search tries
instead (ALGORITHM §2.3): `P*_D = [ln(K·m·n_D) − ln(π/(1 − π))]/λ`, the
score at which `π/(1 − π)` chance hits are expected, with the lattice
`K = 0.408` and `λ = 0.634` of `(2, −3)`, `m = 25` positions and
`n_D = Σ_a (L_a − 3)⁺`, the offsets of the distinct D alleles at which an
exact run of 4 can start: 802 for the 37 human alleles and 659 for the 36
synthetic ones. `ln(K·m·n_D)` is 9.01 and 8.81 nats, against `ln(N·L)` =
6.62 and 6.55. The per-read reading (`P*_D/2` matched D columns, §9.8) is
unchanged.

**The count against brute force** (random 25-base windows against the
human alleles, best ungapped hit per allele, the exact-run gate; 1500
windows):

| S | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|
| chance hits ≥ S, brute force | 10.85 | 4.17 | 3.06 | 1.13 | 0.79 |
| `K·m·n_D·e^{−λS}` | 14.49 | 7.69 | 4.08 | 2.17 | 1.15 |

Exact runs carry 85–95 % of the chance hits at even scores, so the
closed form is high by 1.35 there and by 1.85 at odd scores; BLAST's
length adjustment would make it low by 0.55. A unit test checks the bracket
on random alleles.

**What the posterior wants**, measured offline on the default's output by
moving the D cost with the calls held: P*_D + 0.25 to + 0.5 improves
D-presence ECE and Brier on all six sets; the Brier optimum is + 1.1 (human
mid), + 0.25 (human hard), + 1.5 (synthetic mid), + 0.75 (synthetic hard).
By junction length, junctions of up to 8 nt want 1.1–1.7 nats more (up to
12 nt on the mid sets) and those of 17 nt or more want less: IgSim's D-free
reads have short junctions. By `p̂`, the mid sets want a flat shift and the hard sets more
cost at high `p̂`.

**What the calls allow**, paired against the default on 20,000 reads per
preset:

- A constant + 0.5 (P*_D 10.95 and 10.84): D-presence Brier better on all
  six sets, but D on 4–7 nt remnants −0.56 pp on human hard (p 0.008) and
  −0.39 pp on synthetic hard (p 0.002). Every discordant read moves one
  way: the same D accuracy against D-free specificity curve as every prior
  change before.
- The count with the read's own D window, `m = W` (π = 0.94): Brier worse
  on all six sets (t + 6 to + 23), and 8.9–9.8 pp fewer D-free synthetic
  reads without a D. Long junctions carry a D more often than the window's
  larger chance-hit count allows for.
- The count at the simulator's rate, π = 0.94 (P*_D 9.88 and 9.57, below
  10): exact 5-mers win. No D call −24.0 / −14.5 / −26.3 / −18.4 pp (human
  mid / hard, synthetic mid / hard), Brier worse on all six sets.
- The count at π = 0.9 (P*_D 10.75 and 10.44):

| Set | V | D allele | D 1–3 nt | D 4–7 nt | D 8+ nt | no D | J | Junction |
|---|---|---|---|---|---|---|---|---|
| igh.easy | 0 | 0 | — | 0 | 0 | 0 | 0 | 0 |
| igh.mid | 0 | −0.05 (0.002) | −1.17 (0.016) | −0.14 (0.5) | −0.01 (1.0) | +0.18 (0.5) | 0 | 0 |
| igh.hard | 0 | −0.09 (2e−5) | −1.77 (0.001) | −0.21 (0.25) | −0.02 (0.25) | **+0.77** (0.008) | 0 | 0 |
| syn.easy | 0 | 0 | — | — | 0 | 0 | 0 | 0 |
| syn.mid | 0 | −0.02 (0.12) | −0.17 (0.5) | −0.08 (0.5) | 0 | +0.15 (0.25) | 0 | 0 |
| syn.hard | 0 | −0.03 (0.06) | −0.34 (0.12) | −0.04 (1.0) | 0 | +0.20 (0.12) | 0 | 0 |

| Set | D presence Brier (paired t) | D presence ECE % | D allele ECE % | D allele Brier t |
|---|---|---|---|---|
| igh.easy | 0.00035 → 0.00029 (−5.2) | 0.10 → 0.09 | 0.32 → 0.32 | −1.4 |
| igh.mid | 0.02462 → 0.02398 (−9.1) | 1.36 → 1.13 | 0.66 → 0.63 | −5.9 |
| igh.hard | 0.02168 → 0.02158 (−1.5) | 0.47 → 0.39 | 0.87 → 0.78 | −5.8 |
| syn.easy | 0.00029 → 0.00027 (−5.5) | 0.12 → 0.12 | 0.59 → 0.59 | +6.7 |
| syn.mid | 0.04547 → 0.04485 (−20.0) | 3.27 → 3.08 | 1.67 → 1.63 | −13.9 |
| syn.hard | 0.04079 → 0.04041 (−11.9) | 1.82 → 1.60 | 1.50 → 1.44 | −11.0 |

(The synthetic easy D allele Brier moves by less than 10⁻⁵.) V, J and the
junction do not change; V and J ECE move by at most 0.01.

At π = 0.9 the switch meets every adoption criterion: D-presence
calibration better on all six sets, no D call and D on remnants of 4 nt or
more nowhere worse, the overall D allele loss on human mid and hard all in
1–3 nt remnants (undetectable under `min_d_match` = 4) or within noise.

**Why it is not the default.** The count is derived from the search, the
rate is not, and the rate decides. 0.9 was found by measuring on the
simulator; the rate a user would know or estimate breaks specificity. An
EM estimate over the batch, offline, recovers IgSim's rate (0.942 against
0.944 on human mid, 0.899 against 0.901 on synthetic mid, 0.934 / 0.948
and 0.872 / 0.902 on the hard sets) and is batch-dependent; with it the
count admits exact 5-mers as above. What passes is a band, P*_D from 10.45
to about 10.8, while the count is uncertain by more than that: the mean D
window is 52–58 positions, not 25, and Σ L_a instead of `n_D` adds 0.13
nats, together 1.4 score units. So the default stays the per-read
model's `ln(N·L)` at even odds; `--d-presence-rate` serves a user who
knows the repertoire's rate. The posterior's residual tracks junction
length, which a prior on N/P lengths would model.

Instructions on 2000 hard reads: 9.715 G with the default and with
`--d-presence-rate 0.9` (the prior is one constant per database).

## 10. Annotation priors (roadmap item 5)

### 10.1 What was wrong

- **The harness's own germline set was misread.** `IGHV.fasta` is
  IMGT-gapped, and the parser dropped the dots. With them went the
  numbering that fixes every region start and Cys104. Twelve V alleles got a
  wrong anchor:
  - IGHV1-68\*01–03, IGHV3-30-2/-33/-52, IGHV3-32, IGHV3-33-2 and IGHV1-38-4
    were pseudogenes placed on a TGT two codons upstream;
  - three 3′-partial IGHV1/OR16 alleles got an anchor they do not hold.
- **The truth inherited those errors.** The simulator's junction truth uses
  the same annotation, so the harness scored those reads' junctions as right
  (ALGORITHM.md §7).
- **FR3 starts in TR and lambda loci came from one chain default.** A whole
  family whose FR3 length differs (TRA 32–34, TRB 37/38, TRG/TRD 37–39,
  lambda 36/38 codons) was wrong without a flag.
- **"Flagged" had no definition.** The tests and the scoring scripts used
  different notes.

### 10.2 What replaced it

- **IMGT gap numbering.** The parser keeps the gap runs, and `ImgtNumbering`
  maps the IMGT columns onto bases.
  - Each start is known, before the sequence, or unknown.
  - Frameshift markers are accepted.
  - The policy lives in `annotate()`: the gaps are used only when most
    gapped alleles are codon-aligned with C/W/C at IMGT 23/41/104.
  - Header keys win over the gaps.
- **FR3 starts that are not given** come from:
  1. the nearest relative with a given FR3 start;
  2. else `Fr3Model`, a Jeffreys posterior over FR3 and CDR2 lengths
     (ALGORITHM.md §3a);
  3. else the chain default.

  The provenance note is written after the merge, only for the value kept.
- **Cys23 when the window holds several `C` or none** comes from the nearest
  relative with a settled Cys23, in any frame, else from `Cys23Model`.
- **A structured `uncertain` bit per region start**, written as
  `uncertain=…`. Tests and scripts count an error as silent only when the
  bit is clear.

### 10.3 Measurements

Annotation is scored against IMGT's gap columns: the human IGH set, and 13
tcrdist3 loci (human and mouse TRA/TRB/TRG/TRD where present, IGH, IGK and
IGL). "Main" is the code before this item.

| Check | Main | Now |
|---|---|---|
| Gapped human IGH (415 V): anchors | 370 right, 9 wrong, 3 spurious | 379 right, 0 wrong; the 36 without one are 3′-partial |
| — CDR1 / CDR2 / FR3 starts | 373 / 389 / 399 right of 413 / 412 / 412 | all right |
| — preparation | 28 ms | 1.0 ms |
| Same set, dots stripped from every non-\*01 / half / 90 % of alleles: silent anchor, CDR1 or FR3 errors among them | — | 0 / 0 / 0 |
| — CDR1 missing among them | 22 / 18 / 38 | 1 / 4 / 16 |
| IGH without dots, 13 loci without headers | — | anchor, frame, FR1, FR2, CDR2, FR3 and source identical; 29 CDR1 starts missing → right, none wrong |
| Headers on \*01 alleles: FR3 of the others | 190 silent, 5 flagged | 0 silent, 1 flagged |
| Headers on 10 %: FR3 of the others | 341 silent, 74 flagged, 1082 right | 0 silent, 73 flagged, 1432 right |
| — CDR1 missing (\*01 / 10 %) | 17 / 64 | 5 / 50 |
| A gapped set with dots moved one codon in FR1 (not IMGT) | — | falls back, `annotation.tsv` identical to the ungapped run |
| `every_error_is_flagged` (12,000 alleles) | 3 flagged CDR1 errors | 0 errors |

Wherever headers are given, anchor, frame, FR1, FR2 and CDR2 are
byte-identical to main.

Without any given FR3 start, 539 of 1,695 FR3 starts are wrong (32 %). All
but one are flagged; the exception is IGHV4-31\*05, whose FR3 lacks a codon
and which the gapped set numbers correctly.

**Why the FR3 flag sits at p < 0.95, not ½.** With 10 % of alleles
annotated, model values were right as follows:

| Posterior | Right |
|---|---|
| `p ≥ 0.95` | 309 of 309 |
| `0.8 ≤ p < 0.95` | 102 of 115 |
| `p < 0.8` | 11 of 22 |

Relatives were right 721 of 721 times. At ½, 14 errors would stay silent:
TRGV1–8 at `p` = 0.58–0.92, when only TRGV9–11 carry labels, and
TRAV2/40 at 0.84. A transferred start now inherits its donor's flag. Before, 12 mouse
TRGV starts copied from an unflagged chain-default donor were silent.

**Reads** (20,000 per preset, `measure.sh`, paired with exact McNemar):

- **Allele calls.** V, D and J have 0 discordant reads. "No D call" differs
  on 1 hard read (87.22 → 87.13 %, p = 1.0). Every read whose call text
  changed comes from one of the 12 re-anchored alleles: mostly reordered
  D ties.
- **Junctions.** The truth was regenerated with the new build (seed 7,
  reads byte-identical), and both builds were scored on it:

  | Preset | Main | Now | Reads that moved | McNemar p |
  |---|---|---|---|---|
  | easy | 97.76 % | 100.00 % | 409 | 1.5e−123 |
  | mid | 97.50 % | 99.99 % | 368 | 3.3e−111 |
  | hard | 96.94 % | 99.84 % | 246 | 1.8e−74 |

  Every moved read went from wrong to right. On the stored truth, the same
  outputs read 100.00 → 97.00, 99.99 → 96.26 and 99.82 → 92.96 %. That gap
  is the relabelling, not a loss: §9.5's junction column was scored on the
  old truth.
- **Synthetic database.** Outputs are byte-identical.
- **Calibration** is unchanged (ECE %, V / D / J / D presence): mid
  0.09 / 1.02 / 0.16 / 2.23, hard 0.33 / 0.98 / 0.35 / 0.99.
- **Speed.** Per read: 1933 / 1967 vs 1932 / 1994 reads/s/thread, hard, two
  runs each.

**Database preparation** (best of 7–15 runs):

| Database | Main | Now |
|---|---|---|
| Gapped IGH | 28 ms | 1.0 ms |
| Ungapped IGH | 26 ms | 29 ms |
| Mouse TRA, 10 % headers | 16–17 ms | 18–19 ms |
| Human TRB, 10 % headers | 4.8 ms | 7.0 ms |
| 3,000 synthetic alleles (test; bound 200 ms) | ~40 ms | 43–79 ms |

### 10.4 What is still weak

- **An FR3 family with no label is extrapolated.** It is flagged, but at
  `p < 0.95` about one model value in ten is wrong.
- **One `C` in the Cys23 window is taken as Cys23 without question.** A
  mutated Cys23 next to a mutation-made decoy would be a silent CDR1 error.
  None was seen.
- **Pre-existing silent errors on stripped IGH alleles:**
  - IGHV3-41\*02's CDR2 start, transferred from IGHV3-7\*03 across an indel;
  - IGHV1/OR15-4\*01's FR2 start (stripped half only).

## 11. The constant region (roadmap item 6)

### 11.1 What was decided

- **C is a segment kind, not a lattice kind.** `SegmentKind::C` exists for
  databases, reports and placements, but adjudication iterates
  `SegmentKind::REARRANGED` only, and debug assertions keep C placements out
  of the lattice and the V(D)J report. A C kind in the DAG would let C
  evidence move the J exit or the J allele (the prototype saw 10 J flips in
  6,000 reads when a C's repeated `G` completed IGHJ6\*02/\*03), every read
  would pay for it, and the read's SHM divergence does not apply to C.
- **A stage after adjudication** (`constant.rs`) reads finished calls only:
  V, D and J are unchanged by construction, and a run without `--c` does no
  C work. The model, the anchor column that makes the anchored DP exact on
  every kernel, the reach classes and the gate are in ALGORITHM.md §4b.
- **The gate is derived, not tuned:** `E_C = N_C·(2s + 1)·e^(−λS) ≤ 10⁻³`,
  the anchor gate of V and J.
- **The prototype's accuracy metric was near-tautological.** Its truth
  counted an allele equivalent when it *contained* the C part anywhere, and
  its gene tally accepted any equivalent allele's gene. On IGHC a 16–32 nt C
  part then spans 2.6 genes on average, so "100 % at 16 nt" mostly said that
  the hit lay in the right prefix class. The truth is now *anchored* (the
  allele's germline from its own splice offset equals the simulated bases),
  and the table adds the first name's gene, exact gene sets, over-calls and
  the oracle's resolution.
- **Other corrections of the design review.** No sibling bound table for C
  (its Δ sums over CH2/CH3, beyond any read's reach; reach classes already
  remove the redundancy). The 16-bit lane limit depends on the longest C
  record, not on read length, so `prepare --c` warns about records over
  2 kb or without a stop-free frame (genomic C-GENE records). `s = 3` rested
  on the measured histogram below; since 0.39.4 `s = 7`, with SwiftIG's
  floors (§11.4). The truth gains `c_` columns only with
  `--c`, so existing datasets stay byte-identical.

### 11.2 Measurements

**Datasets.** `simulate --c`, seed 7, 20,000 reads per preset, on the same
V(D)J model and seeds as the six existing sets:
- IMGT human IGH with abstar's human IGHC (88 records, 81 distinct, 10 genes,
  codon-aligned: each record repeats the J's last base);
- the synthetic database (seed 1) with its synthetic C (25 alleles, 14 genes
  in 5 isotypes whose subclasses share CH1's opening; IMGT-style records).

80 % of reads read into C, the whole 3′ flank for 60 % of those, else a
uniform part of it. About 4,100 reads per set carry no C base.

**Invariance (must hold).**
- Without `--c`, AIRR output is byte-identical to the previous release on
  all six existing sets, for each of the AVX-512, AVX2, portable and scalar
  kernels; `simulate` without `--c` writes byte-identical reads and truth.
- With `--c`, the V/D/J columns are byte-identical to the run without it on
  the six existing sets (AVX2 and scalar) and on the six C sets.

**C accuracy** (per C-length stratum, the range over easy/mid/hard; equivalence
is anchored; gene set = the reported genes are exactly the equivalent ones;
over-call = a reported name outside the equivalent set, per call):

| Database | C part | Reads | Called | Primary equivalent | First name's gene | Gene set exact | Over-call | Genes listed / equivalent |
|---|---|---|---|---|---|---|---|---|
| IGH | none | 4109–4141 | 0.00–0.05 % | | | | | |
| IGH | 1–7 nt | 508–550 | 1.0–2.6 % | 1.0–2.6 % | | | 0 | |
| IGH | 8–15 nt | 592–643 | 90.3–91.9 % | 90.0–91.6 % | 60.3–61.5 % | 89.5–91.1 % | ≤ 0.37 % | 2.6 / 2.6 |
| IGH | 16–31 nt | 1168–1246 | 99.8–100 % | 99.60–100 % | 66.5–68.7 % | 97.9–99.0 % | ≤ 0.40 % | 2.6 / 2.6 |
| IGH | 32–63 nt | 3867–4330 | ≥ 99.95 % | 99.86–99.90 % | 90.3–90.5 % | 97.9–98.6 % | ≤ 0.28 % | 1.73 / 1.75 |
| IGH | 64–127 nt | 7455–8347 | ≥ 99.88 % | 99.79–99.96 % | 94.2–95.1 % | 99.8–99.9 % | ≤ 0.12 % | 1.4 / 1.4 |
| IGH | 128+ nt | 1051–1896 | ≥ 99.88 % | 99.82–99.95 % | 96.7–97.4 % | 99.8–100 % | ≤ 0.11 % | 1.1 / 1.1 |
| synthetic | none | 4001–4138 | 0.02–0.05 % | | | | | |
| synthetic | 8–15 nt | 575–649 | 78.4–81.6 % | 78.4–81.6 % | 28.7–30.7 % | 74.9–77.4 % | 0 | 2.9 / 3.0 |
| synthetic | 16–31 nt | 1164–1272 | ≥ 99.84 % | 99.84–100 % | 50.5–52.2 % | 91.7–93.6 % | ≤ 0.26 % | 2.4 / 2.6 |
| synthetic | 32–63 nt | 3846–4284 | ≥ 99.91 % | 99.81–99.93 % | 95.1–96.1 % | 98.9–99.0 % | ≤ 0.26 % | 1.14 / 1.15 |
| synthetic | 64+ nt | 9231–9755 | ≥ 99.80 % | 99.75–100 % | 99.75–100 % | 99.7–100 % | ≤ 0.12 % | 1.0 / 1.0 |

- **No C call on reads without C:** 0, 0, 2, 2, 2 and 1 false calls in about
  4,100 C-free reads per set (≥ 99.95 %). On the six existing sets, whose
  3′ flanks are random, 4–8 reads in 20,000 get a C call.
- **Primary equivalent ≥ 99.6 % in every stratum from 16 nt.** The misses
  are instrument errors that turn the true allele's base into another
  allele's (8–14 per set) and, on hard reads, reads without a J call, a
  correct orientation or a J exit next to C (16–23 per hard set).
- **8–15 nt** is limited by the gate: `S ≥ 18` needs 9 exact bases on IGH
  (10 genes), `S ≥ 19` needs 10 on the synthetic database (14 genes); every
  call there is equivalent. (0.39.4: the 12-column floor leaves 8–11 nt
  uncalled, §11.4.)
- **Exact gene sets trail the anchored oracle** at 16–63 nt (IGH 97.9–99.0 %,
  synthetic 91.7–99.0 %). Of the misses from 16 nt, 679 are on reads whose
  C part is followed by other bases (2.2 % of 30,220 such reads), 34 on reads
  ending in C (0.06 % of 58,239): a chance match after the C part extends some
  alleles, which no read-level caller can tell from a C base. The oracle
  knows where C ends; the read does not.
- **The first name is a guess among genes for short C parts:** 61 % on IGH at
  8–15 nt and 67–69 % at 16–31 nt, where the listed genes are right 98–99 %
  of the time. Read the isotype from the gene set of `c_call`, not its
  first name.

**Against the local K·m·n gate** (same calls, kept when
`K·m·n·e^(−λS) ≤ 10⁻³` with `m` the bases after the J exit and `n` the C
database length; paired, exact McNemar):

| Set | C allele, K·m·n → anchored | p | No C call | p |
|---|---|---|---|---|
| IGH easy / mid / hard | +3.93 / +4.06 / +4.25 pp | < 10⁻¹⁸⁷ | 0 / 0 / −0.05 pp | 1 / 1 / 0.5 |
| synthetic easy / mid / hard | +3.01 / +3.00 / +3.18 pp | < 10⁻¹⁴² | −0.05 / −0.05 / −0.02 pp | 0.5 / 0.5 / 1 |

**Calibration of `c_posterior`** (event: some reported name is
anchored-equivalent): ECE 0.13–0.18 %, Brier 0.0002–0.0006 on every set.
Only 7–14 reads per set have `c_posterior < 0.95`, too few to judge
calibration below it: with C's own divergence near zero, ties are the only
real ambiguity, and ties are listed.

**Read content without `--c`** (C reads against the same reads without C,
paired on the same V(D)J truth): V, D, no D call and junction are within
noise everywhere. J allele is −0.01 / −0.10 / −0.52 pp on IGH easy / mid /
hard (p = 0.5, 1·10⁻⁶, 5·10⁻³⁰) and within noise on the synthetic database. Every
lost read is an IGHJ6\*02 or \*03 read called IGHJ6\*01 or \*04: those IMGT
records are 3′-partial (one base short), and the base that follows them in a
read through C (the C's `G`, which the record lacks) completes the longer
alleles. It is a property of the reads and the J records, not of the C
stage, and `--c` changes no J call.

**Cost** (callgrind, 2,000 IGH C reads, one thread):

| Reads | Without `--c` | With `--c` | C share | of which the stage / report / database |
|---|---|---|---|---|
| mid | 9,684,669,685 | 10,661,256,792 | 9.2 % | 742 M / 108 M / ≈ 100 M |
| hard | 9,230,962,862 | 10,176,966,565 | 9.3 % | 710 M / 112 M / ≈ 100 M |

The database part is a fixed cost per run: the gene frames behind `prepare`'s
truncation warning. Without `--c` the C stage costs nothing: 9,230,582,300
instructions on the standard 2,000 hard reads against 9,233,069,238 before
the stage existed (and 9,330,252,064 before this item's kernel
prerequisite, −1.1 %).

**Where the J exit lies** (`x` minus the read position of C's germline
base 0; IGH reads have `o = 1`, except after the 3′-partial IGHJ6\*02/\*03):

| Set | −1 | 0 | +1 | +2 | other |
|---|---|---|---|---|---|
| IGH mid | 1 | 2,426 | 13,442 | 1 | 0 |
| IGH hard | 5 | 2,334 | 13,520 | 16 | 7 (+91 to +100) |
| synthetic mid | 4 | 15,850 | 8 | 0 | 0 |
| synthetic hard | 31 | 15,922 | 22 | 1 | 3 (−23, +4, +10) |

The IGH conventions take one base of the spread; SHM indels add one more,
rarely two.

**Sensitivity** (20,000 hard reads; false calls on about 4,000 C-free reads,
equivalence at 8–15 nt and from 16 nt, exact gene set from 16 nt):

| Set | s | E = 1 | 10⁻¹ | 10⁻² | 10⁻³ | 10⁻⁴ | 8–15 nt at 10⁻³ | ≥ 16 nt | gene set ≥ 16 nt |
|---|---|---|---|---|---|---|---|---|---|
| IGH | 1 | 928 | 104 | 19 | 1 | 0 | 91.1 % | 99.69 % | 99.21 % |
| IGH | 3 | 812 | 101 | 16 | 2 | 0 | 90.8 % | 99.80 % | 99.30 % |
| IGH | 6 | 967 | 107 | 14 | 3 | 0 | 78.5 % | 99.80 % | 99.29 % |
| synthetic | 1 | 983 | 99 | 5 | 1 | 1 | 91.2 % | 99.78 % | 98.95 % |
| synthetic | 3 | 639 | 92 | 7 | 1 | 0 | 78.4 % | 99.79 % | 98.94 % |
| synthetic | 6 | 458 | 99 | 12 | 3 | 0 | 78.4 % | 99.80 % | 98.93 % |

The equivalence from 16 nt does not depend on the gate between 1 and 10⁻⁴.
False calls stay 2–6× below the bound `E` per read at every gate from 1 to
10⁻³. `s = 1` misses the IGH reads two bases off (−0.11 pp from 16 nt on
hard); `s = 6` adds start hypotheses, which raise the gate by a base at
8–15 nt. `s = 3` costs the synthetic database's 10-base threshold at
8–15 nt (78 % against 91 % at `s = 1`).

**Deterministic null.** After every one of the 13 human IGHJ alleles, eight
fixed tails (TruSeq read 1 and 2 adapters, Nextera, poly-A, poly-G, poly-T,
and two splice donors with intron) score at most 7 against abstar's IGHC
(clipped, any allele and start); the gate needs 18. `tests/robustness.rs`
runs the same check on the synthetic database (best 4–6; a gate of E ≤ 10
would call them).

### 11.3 What is still weak

- **The first `c_call` name** among genes sharing a short C part is
  arbitrary. The gene set is right; the first name is not an isotype call.
- **Bases after C extend some alleles by chance**, so the reported genes can
  be fewer than the read's true ambiguity (2.2 % of C parts followed by
  other bases, from 16 nt). `c_posterior` does not see it: the extension is
  real evidence to the model.
- **Short C parts** need 9–10 exact bases to pass the gate (0.39.4: 12
  aligned columns at 90 % identity, §11.4).
- **Net indel drift** inside C beyond `s ∓ δ` truncates the C alignment.
- **3′-partial J records** (IGHJ6\*02/\*03) lose J calls on reads that
  continue into C, with or without `--c` (−0.52 pp J on hard IGH reads).
- **`c_posterior` below 0.95** is too rare in simulation (≤ 14 reads per
  20,000) to check its calibration there.

### 11.4 Note (0.39.4): SwiftIG's C anchoring and floors

`s = 3` above rests on simulated reads, whose J ends at its germline end
and whose C starts at its gene's base 0 (the simulator draws only C
alleles that lack at most 3 of their gene's first bases). SwiftIG, the
aligner of Swig ≤ 0.38, anchored C more loosely, taking C's first aligned
base up to 8 bases after the J's end and up to 3 bases into C, and called C
on at least 12 aligned columns at 90 % identity. On 2,552 constructed reads
(every human IG and TR C gene, 15–300 nt of C, J 3′ ends cut by 1–6 bases,
1–5 bases between J and C, C's first bases deleted, 2–20 % C divergence,
both strands):

| | SwiftIG | 0.39.3 (`s = 3`) | 0.39.4 (`s = 7`, floors) |
|---|---|---|---|
| C right | 2,035 | 2,083 | 2,074 |
| wrong C gene | 4 | 41 | 1 |
| right in SwiftIG, not here | — | 163 | 33 |

0.39.3 lost the C of reads whose J end lay more than 3 bases from C's
start, called IGLC3 for IGLC2 and TRBC2 for TRBC1 there, could not name
IGLC1\*01 or IGLC3\*01 (6 bases shorter than their genes), and called
IGHG2\*16 at 79 % identity on a read whose C part was IGKC. 0.39.4 sets
`s = 7` (15 start diagonals, so the gate needs `S ≥ 19` on human IGHC
instead of 18) and adds SwiftIG's floors, applied to each C placement
before ranking (`--c-min-identity` sets the identity floor). Of the 33
reads still behind SwiftIG, 29 have no J, and SwiftIG called their C only
after a chance J; 4 have C at 12–20 % divergence, where the anchored
alignment counts leading mismatches that SwiftIG's local one trimmed.
`prepare` now warns about C records missing more than 7 leading bases and
J records whose tail strays by more than 7.

The price is short C parts. On 20,000 simulated human BCR reads per preset
(primary C allele anchored-equivalent, by C length; false calls on reads
without C):

| C part | mid 0.39.3 | mid 0.39.4 | hard 0.39.3 | hard 0.39.4 |
|---|---|---|---|---|
| none: calls | 1 of 4,097 | 0 | 1 of 4,006 | 0 |
| 1–7 nt: called | 5.2 % | 0 | 8.1 % | 0 |
| 8–15 nt | 94.1 % | 55.6 % | 93.2 % | 55.5 % |
| 16–31 nt | 99.83 % | 99.30 % | 99.69 % | 99.30 % |
| 32 nt and more | 99.84–100 % | within 0.02 pp | 99.53–99.71 % | within 0.02 pp |
| C allele, all C reads | 96.71 % | 94.94 % | 96.04 % | 94.15 % |

Of the 584 reads that lost a right C call (281 mid, 303 hard), 580 have a
C part shorter than 32 nt (the other four: 33–102 nt; two hard reads gained
one); 8–11 nt fall under the 12-column floor, as they did in SwiftIG, and
Swig Web's isotype column never labelled them. `c_posterior` ECE is 0.13–0.14 %. `tests/simulated.rs`
pins the synthetic database's 8–15 nt floors at 47.83 % (mid) and
63.89 % (hard) accordingly.
