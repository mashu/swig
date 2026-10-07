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
2. **D presence on D-free junctions.** 12–24 % of reads with no D receive one,
   because a 9-nt chance match is indistinguishable from a 9-nt D remnant.
   `P_D` sets the trade-off, and the prior odds are the single knob.
   Precision/recall curves, not one number, belong in the paper; the
   reported `d_presence_posterior` now carries this uncertainty per read.
   The derivation itself undercounts the search: `ln(N·L)` counts genes ×
   read positions, but a local alignment tries every allele, allele offset
   and read position (Karlin–Altschul `K·m·n ≈ e^{10}` against `e^{6.7}` for
   human IGH; ~0.26 expected chance 9-mer matches per junction, in line with
   the 22 % false-positive rate on mid). Even prior odds understate IgSim's
   99 % D inclusion, and the two errors roughly cancel there. Both terms
   should be estimated, not assumed.
3. **Calls and posteriors use different evidence scales.** The Viterbi call
   reads scores at the fixed λ, the posterior at the read's divergence, so
   the posterior's mode can disagree with the call (a short exact D against
   a longer one with a mismatch on an unmutated read). Measured on 4000
   reads per preset: taking D presence *and* allele from the posterior
   (maximum marginal) gains 0.2–0.8 pp of D accuracy but raises the D-free
   false-positive rate by 13–16 pp, so total D errors do not fall (mid 287
   vs 271, hard 681 vs 685). Keeping the Viterbi presence decision and
   taking only the allele from the posterior gains 0.1–0.3 pp at an
   unchanged false-positive rate; it is not adopted because the reported
   alignment would then need a second, conditioned Viterbi pass.
4. **One alignment per D allele and window.** A D's optimal local alignment
   in its window is not always its optimal clipped alignment at the chosen
   cuts; the 12-base reclaim keeps this rare. It is the one remaining
   approximation of the segmentation model.
5. **Constants inherited from the note** — now measured (§7): five are
   inert and default to neutral values; the rest are speed knobs or have a
   measured optimum at the inherited value. `tandem` (12) was not swept:
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
11. **Anchors without a relative.** Seven V alleles of the human IGH set
    (pseudogenes and ORFs) keep a low-confidence motif anchor. One of them,
    IGHV1-38-4\*01 (Cys104 → Tyr), is placed on a TGT in another frame;
    without a relative there is no evidence to place it better, so it is
    reported as low-confidence rather than silently trusted.

## 6. Plan

In order of expected value:

1. ~~**Posterior decoding**~~ — done (§3, item 12). Next on that line:
   calls from the per-read model (limitation 3) and a D prior from the D
   search's own statistics (limitation 2).
2. ~~**Two cuts per overlap with N/P between**~~ — done (§3, item 11).
3. ~~**Sensitivity analysis**~~ — done (§7).
4. **Cross-tool benchmark.** IgBLAST done (§8). partis and Swig remain, at
   several SHM levels, with precision/recall curves for D.
5. ~~**AVX-512BW (32 lanes) and NEON**~~ — done (§9.4), with a portable
   backend for every other target.
6. **Allele rescoring along the gene frame.** Trace one representative per
   frame, rescore siblings along its path, and run DP only where bounds
   straddle the top-N cutoff. Most of the remaining V time is siblings of the
   top gene.
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
| V acceptance | 24 columns, 60 % | 12–48 columns; 50–70 % | flat | none |
| J acceptance | 10 columns, 60 % | 6–18 columns; 50–70 % | flat | none |
| Distributed-D gate | on (16/14/80 %/19) | off; looser; tighter | flat | off |
| Retained bases V / J | 24 / 10 | V 1–36, J 1–14 | flat | 1 / 1 |
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
  (within 0.03 pp), and kept random-read calls at 0.04 % of 20,000.
- **The derived `P_D` is the calibration optimum.** D-presence calibration
  error is lowest at `P_D = 11` on both presets (mid 2.25 %, hard 1.66 %; 2.71 %
  and 2.59 % at 9, 3.36 % and 5.37 % at 13), so the posterior independently
  supports the value the derivation gives. Accuracy alone cannot choose it:
  D accuracy and D-free false positives trade monotonically.

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
scored by `vdj-assign evaluate`:

| Preset | Tool | V | V gene | D | J | Junction | D called on D-free reads | Reads/s, 4 threads (previous V retrieval) |
|---|---|---|---|---|---|---|---|---|
| easy | vdj_assign | 99.88 | 100.00 | 99.29 | 99.84 | 100.00 | 0 % | 5200 |
| easy | IgBLAST | 99.12 | 99.78 | 98.86 | 99.86 | 99.77 | 0 % | 67 |
| mid | vdj_assign | 97.88 | 99.72 | 93.35 | 98.74 | 99.97 | 23 % | 5100 |
| mid | IgBLAST | 97.18 | 99.62 | 93.73 | 98.60 | 99.72 | 43 % | 66 |
| hard | vdj_assign | 89.44 | 98.84 | 83.00 | 92.80 | 99.69 | 16 % | 5200 |
| hard | IgBLAST | 89.16 | 98.90 | 81.75 | 91.84 | 98.07 | 30 % | 71 |

IgBLAST's higher mid D accuracy comes with twice the D calls on D-free
junctions, and its V gene is 0.06 pp better on the hard preset. Two setup
issues had to be fixed for a fair junction comparison. IgBLAST takes J
anchors from its auxiliary file, which lacks three of this database's J
alleles, so it reported no junction for 18 % of reads; the entries were
added, with frames and CDR3 ends from `vdj-assign prepare`, which agrees
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
  one of them gets right (`vdj-assign evaluate --truth T a.tsv b.tsv`,
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
input parsing with assignment; (5) allele rescoring along the gene frame
(plan item 6).

### 9.7 Roadmap status

Of the roadmap in the README, the wider-SIMD item is done (AVX-512BW, NEON
and a portable backend). Calls from the per-read model, a D prior from the D
search's own statistics, allele rescoring along the gene frame, annotation
priors and C-segment assignment remain open; §9.5 shows the D prior is the
one with the clearest measured need.
