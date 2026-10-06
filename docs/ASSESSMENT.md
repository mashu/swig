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
15. **Closed-syncmer seeds** (`seed.rs`). See §8.
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

## 4. Results

Human IGH, 5000 reads per preset, germline-equivalent accuracy in %:

| Preset | Orientation | V | V gene | D | J | Junction | Random reads called |
|---|---|---|---|---|---|---|---|
| easy | 99.94 → 100 | 99.54 → 99.84 | 99.92 → 99.98 | 98.78 → 99.29 | 99.90 → 99.94 | 99.94 → 100 | 30.5 → 0 |
| mid | 100 → 100 | 97.76 → 97.92 | 99.84 → 99.82 | 91.96 → 93.00 | 98.84 → 98.80 | 99.93 → 99.97 | 30.2 → 0 |
| hard | 99.98 → 99.98 | 89.74 → 89.86 | 98.44 → 98.36 | 81.72 → 83.58 | 92.70 → 93.32 | 98.39 → 99.50 | 31.2 → 0 |

Junctions "after" are scored with the corrected anchors (item 16), which
remove the reads of 3′-partial alleles from the denominator (hard: 2416 →
2214 reads); with the old annotation the same calls scored 98.47 %.

Throughput is 849 → ~1900 reads/s per thread without posteriors and ~1300
with them, and ~3600 → ~7300 / ~5400 reads/s on 4 threads (a busier
container than the earlier ~2200 / ~9000; in instructions the exact DAG
costs 3 % and posteriors 18 %).

Calibration on real IGH (ECE %, easy / mid / hard), with the event of items
14 and 17:

| | V | D | J | D presence |
|---|---|---|---|---|
| ECE | 0.04 / 0.31 / 1.31 | 0.49 / 1.29 / 1.16 | 0.07 / 0.30 / 0.69 | 0.20 / 1.98 / 1.24 |
| Brier | 0.001 / 0.009 / 0.045 | 0.002 / 0.040 / 0.072 | 0.000 / 0.006 / 0.031 | 0.003 / 0.026 / 0.024 |

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
8. **SIMD portability.** AVX2 only; AVX-512 hardware runs at half width, and
   ARM falls back to scalar.
9. **D errors are mostly the wrong allele, not a missing D.** With a remnant
   of ≥ 8 nt (§8), 2.7 % (mid) and 10.2 % (hard) of reads get a wrong D
   allele although a placement of the right one was found; the search loses
   0.2 % and 1.8 %. These are model errors: short, mutated D remnants that
   another allele explains as well or better.
10. **Posterior cost.** The forward–backward pass is 18 % of instructions
    but ~30 % of single-thread wall clock (dependent floating-point
    multiply–adds and per-slot normalisation). A single-precision pass
    with the same scaling, or skipping slots whose mass is negligible,
    would cut it.
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
5. **AVX-512BW (32 lanes) and NEON (8 lanes)** kernels of the same design.
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
| Budget margin | 5 % / 8 votes | 2–10 %; 4–16 | ≤ 0.07 pp | 5 % / 8 |
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
the gene-frame index and in the read (`seed::SeedSampling`). Closed syncmers
(a k-mer is kept when its smallest `s`-mer, `s = k − 4`, is first or last)
depend on the k-mer alone. A seed shared by read and germline is therefore
kept in both or neither, and any `k + 4` unmutated bases contain one (13 for
the 9-mer index, 11 for the 7-mer fallback). Minimisers lack this property:
whether a k-mer is kept depends on its neighbours, so a mutation can unselect
a shared seed on one side only. Measured against every k-mer on human IGH:

| Condition | every k-mer | closed syncmers (40 %) |
|---|---|---|
| presets, three 5000-read hard samples, synthetic database | — | V and V gene within ±0.08 pp |
| fixed SHM 15 % (V / gene) | 100.00 / 96.25 | 100.00 / 96.60 |
| fixed SHM 20 % | 99.95 / 94.50 | 99.95 / 95.10 |
| fixed SHM 25 % | 98.70 / 87.00 | 98.85 / 87.40 |
| fixed SHM 30 % | 93.05 / 74.90 | 93.45 / 76.35 |
| 3′ fragments covering 30–80 V bases | — | within ±0.4 pp, both signs |
| instructions, 1000 hard reads | 6.09 × 10⁹ | 5.95 × 10⁹ (−2.4 %) |

So no sensitivity is lost; at high SHM the sparser seeds are marginally
better, because fewer chance seeds dilute the true gene's votes. Minimisers
(w = 4) and open syncmers gave the same accuracy within noise, so closed
syncmers were chosen for their guarantee, not for a measured gain. The gain
is small because V retrieval is 5.5 % of a run (voting, sorting and
clustering became about 60 % cheaper); alignment and adjudication dominate.

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
| D correct | 97.08 % | 87.50 % |
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

| Preset | Tool | V | V gene | D | J | Junction | D called on D-free reads | Reads/s, 4 threads |
|---|---|---|---|---|---|---|---|---|
| easy | vdj_assign | 99.88 | 100.00 | 99.29 | 99.84 | 100.00 | 0 % | 5200 |
| easy | IgBLAST | 99.12 | 99.78 | 98.86 | 99.86 | 99.77 | 0 % | 67 |
| mid | vdj_assign | 97.88 | 99.72 | 93.35 | 98.74 | 99.97 | 23 % | 5100 |
| mid | IgBLAST | 97.18 | 99.62 | 93.73 | 98.60 | 99.72 | 43 % | 66 |
| hard | vdj_assign | 89.26 | 98.58 | 83.00 | 92.80 | 99.69 | 16 % | 5200 |
| hard | IgBLAST | 89.16 | 98.90 | 81.75 | 91.84 | 98.07 | 30 % | 71 |

IgBLAST's higher mid D accuracy comes with twice the D calls on D-free
junctions, and its V gene is 0.3 pp better on the hard preset. Two setup
issues had to be fixed for a fair junction comparison. IgBLAST takes J
anchors from its auxiliary file, which lacks three of this database's J
alleles, so it reported no junction for 18 % of reads; the entries were
added, with frames and CDR3 ends from `vdj-assign prepare`, which agrees
with the file on all 13 shared alleles. The remaining disagreements exposed
item 16 on our side.
