# The vdj_assign algorithm

This document is a methods-section description of what vdj_assign computes
(the model) and how it finds it (the search). It gives the guarantee of each
search step: *exact*, *lossless* (provably drops nothing the model could
choose), *bounded* (drops only what lies outside a stated geometric bound), or
*heuristic* (measured, not proven). Code references are to `src/`.

## 1. Problem

Given a read `q` and germline databases `V`, `D` and `J` (alleles, grouped
into genes), report:
- the strand of `q`;
- the alleles of its V, D (zero, one or two) and J segments;
- the read coordinates of each segment;
- hence the N/P regions and the junction.

The read may carry somatic hypermutation (substitutions and short indels),
sequencing errors, trimmed segment ends, non-templated junction bases, and
unrelated flanking sequence.

## 2. Model

### 2.1 Segments are local alignments scored as log-odds

A segment hypothesis (allele `a`, read interval `I`) is scored by affine-gap
local alignment of `q[I]` to `a`. With match score `m` and mismatch `μ`, a
substitution column scores

```text
s(x, y) = m  if x = y,   μ  if x ≠ y,   0  if x or y is ambiguous (N).
```

These are log-likelihood ratios. Let `λ` solve `¼e^{λm} + ¾e^{λμ} = 1` (the
Karlin–Altschul λ for uniform base composition). Then `e^{λ s(x,y)}` is the
likelihood ratio of "`x` was templated by `y` at divergence δ" against "`x` is
random", with `1 − δ = ¼e^{λm}`. For the default `(m, μ) = (2, −3)`,
`λ = 0.634` nats per score unit and `δ = 11.2 %` (`stats.rs`). Gap runs cost
`g_o + (L−1)·g_e`. N/P and flank bases score 0, since they are random under
both hypotheses.

V and J alignments are extended co-linearly, without gaps, to the germline's
outer end (V 5′, J 3′) or the read's end, whichever comes first
(`refine::extend_outer`). This makes the alignment end-to-end on its outer
side, where the germline must begin (V) or end (J). The substitution scores of
the extension, mismatches included, count towards the ranking.

### 2.2 Paths

A *placement* is an aligned segment hypothesis. A *path* is a sequence of
placements in read order that the grammar allows (`adjudicate/model.rs`):

```text
V → D{0..max_d} → J,     partial paths allowed (V…, …J, V…J);  max_d = 1 (VDJ) or 2 (VDDJ)
```

Between two consecutive placements there are two *cuts*: the first keeps
read bases before its *exit* `x`, the second keeps bases from its *entry*
`y ≥ x` on, and the `y − x` bases between them are N/P. Each placement is
then *clipped* to the alignment columns it keeps. Overlapping placements may
therefore split the overlap at any point, or leave part of it to neither.

### 2.3 Objective

```text
score(path) = Σ_segments S_clip(segment)
            − P_D · #D − tandem · [second D]
            + pair_bonus · [path starts with V and ends with J]  (0 by default)
            + insertion_per_base · #N/P bases                    (0 by default)
```

The non-alignment terms are log prior odds in score units:

- **`P_D`, D presence.** A D hypothesis is one of `N` D genes placed at one of
  `L` junction positions. Its prior is spread over `N·L` alternatives. A MAP
  decision includes the best D when
  `λ·S_D > ln(N·L) − ln π_D`, where `π_D` is the prior odds of a D against
  none. Hence `P_D = ⌈(ln(N·L) − ln π_D) / λ⌉` (`config::DPrior::Derived`).
  This is 11 for human IGH (N = 32, L = 25, π_D = 1) and 6 for TRB (N = 2).
  The legacy profiles set `P_D = 0`.
- **`tandem`**: the prior odds of a second D (DJ → DDJ), on top of its own `P_D`.
- **`pair_bonus`**: the prior odds of a complete V…J rearrangement against a
  partial read. It is 0 by default: on full-length reads the J's own clip
  score completes the path, and values 0–48 changed no call.

The adjudicator returns the path maximising, lexicographically:

1. `score`;
2. fewer D segments ("with-D must strictly beat no-D");
3. *support*, the sum of the placements' unclipped scores.

Key 3 separates paths that explain the read equally well. A typical case is
two alleles of different 3′ length at a junction where one base can belong to
either side. The allele with more support admits more cuts that reach the same
score, so it has the larger likelihood summed over cuts. Co-optimal alleles
(equal clipped score on the called interval) are all reported, the
best-supported first.

### 2.4 Significance

Each V and J hit has a Karlin–Altschul E-value against its database for the
query window that was searched:
`E = K·m·n·e^{−λS}`, with `m` the window length and `n` the database length
(`segment::Significance`).

- A hit is kept only at **E ≤ 1**: at most one chance hit of that score is
  expected.
- A read is *called* only when its path contains an **anchor**, a V or J hit
  at **E ≤ 10⁻³** (`segment::EvalueGates`, `Placement::anchor`).

Random sequence therefore receives no call. A J downstream of a significant V
needs only E ≤ 1, because given the V a J is expected there.

## 3. Placements and the clip index

**Lemma 1 (clipping never helps).** Let `A` be an optimal local alignment and
`A[a..b]` a contiguous sub-alignment beginning and ending in match columns.
Then `S(A[a..b]) ≤ S(A)`.

*Proof.* Every prefix and every suffix of an optimal local alignment has a
non-negative score; otherwise removing it would give a higher-scoring local
alignment. `S(A[a..b]) = S(A) − S(prefix) − S(suffix) ≤ S(A)`. ∎

**Clip index** (`clip.rs`). Flatten `A` into columns with prefix scores `P[k]`
and two query-offset tables: the first match column at or after a position,
and the last match column before it. A clip to `[lo, hi)` is two lookups and
one subtraction: `S = P[exit(hi)] − P[entry(lo)]`. It is exact, because a clip
starts and ends on match columns and so never splits a gap run.

**Separability.** The clipped score is `exit(hi).prefix − entry(lo).prefix`:
one term per cut. The support key is cut-independent.

## 4. Adjudication: exact semi-Markov Viterbi (`adjudicate/dag.rs`)

Placements are nodes, sorted by read start. A state is
`(node, entry cut, D count, first segment kind)`. Its value is the best path
score before the node's own segment, with support as a second key. With a
linear N/P score `ι·(y − x)`, a transition `b → e` that exits `b` at `x`
and enters `e` at `y ≥ x` has the value

```text
max over x ≤ y, and states s of b with q_lo(s) ≤ exit_b(x).q − min_retained(b):
    [score(s) − entry_prefix(s), support(s)]
  + exit_b(x).prefix − ι·exit_b(x).q
  + bonus(b) + ι·entry_e(y).q
```

Every term is a function of `s`, of `x`, or of `y` alone. So:

1. *Best state for an exit.* The maximand depends on `s` only through
   `(score − entry_prefix, support)`, and the constraint is a threshold on
   `q_lo(s)`. Sorting a node's states by `q_lo` once and keeping running
   lexicographic maxima answers each exit with one binary search.
2. *Best exit for an entry.* The best `x ≤ y` is a running maximum over the
   node's exits, built once per node and successor kind and read in O(1)
   per entry.

Entries are tried only inside the overlap, or at `e`'s first base when the
placements do not overlap. Moving an entry further right clips a prefix of an
optimal local alignment, which scores ≥ 0 (Lemma 1), and lengthens the N/P
run, which scores ≤ 0 (`ι ≤ 0`: a per-base N/P score is a length cost, never
a reward), so it never helps.

Complexity is `O(S log S + L log S + T·C)` per node, with `S` states, `L`
exit cuts, `T` successors and `C` entry cuts, instead of `O(S·T·L·C)`.

**Theorem 1.** The DAG returns the optimum of §2.3 over all grammatical
paths, placements and pairs of cuts. This is checked by an exhaustive,
memoised enumeration in `tests/synthetic.rs` that tries every exit and every
entry `y ≥ x` within 12 bases of the meeting point (wider than the DAG's own
entry range), and by a constructed case in `adjudicate/dag.rs` where leaving
a shared mismatch stretch to N beats every single cut.

## 4a. Posterior decoding (`adjudicate/posterior.rs`)

The objective is a log-likelihood ratio, so it also defines a distribution:

```text
P(path, cuts | read) = exp(nats(path, cuts)) / Z,     Z = Σ_{path, cuts} exp(nats)
```

over every grammatical path and every pair of cuts between neighbours. The
Viterbi path of §4 is its mode when the evidence is read at the fixed λ.
Forward–backward over the same DAG gives the marginal probability of every
placement and of every D count.

**Evidence at the read's divergence.** `(2, −3)` is log-odds only at 11.2 %
divergence: a fixed λ makes unmutated reads underconfident (a SNP between
alleles is worth `ln(3(1−p)/p) = 8.0` nats at `p = 0.1 %`, not 3.2) and
heavily mutated ones overconfident. Let `p̂` be the fraction of mismatched
`Match` columns in the called V and J (with half a pseudo-mismatch, clamped
to `[10⁻⁴, 0.5]`). A `Match` column scoring `s` is worth `λ_r·s + μ_r` nats,
with `(λ_r, μ_r)` chosen so that a match is worth `ln 4(1 − p̂)` and a
mismatch `ln(4p̂/3)`, the log-odds of "templated at divergence `p̂`" against
"uniform random"; gaps are worth `λ_r·g`. Each segment kind uses its own
match/mismatch scores. The priors (`P_D`, pair bonus, tandem, `ι`) stay at the
λ they were derived with. At `p̂ = 11.2 %`, `λ_r = λ` and `μ_r = 0`. Because
the map is affine per column, clip scores stay separable in their cuts.

**Domain.** The first segment keeps its start and the last its end, as in §4
(V and J are extended to the read's ends on purpose). Between neighbours every
pair of cuts `x ≤ y` is summed, including entries beyond the overlap that the
Viterbi search skips because they cannot win. Each segmentation is counted
once by taking *canonical* cuts: an entry on the first retained base and an
exit right after the last (one pair per distinct retained sub-alignment).

**Recursions.** For node `b` in group `g` with forward mass `F(b, y)` per
entry,

```text
X(b, x)  = e^{exit(x) + bonus} · Σ_{y ≤ x − min} F(b, y) · e^{−entry(y)}
F(e, y') = Σ_b Σ_{x ≤ y'} X(b, x) · e^{ι(y' − x)}          (+ 1 at a path start)
Z        = Σ_b X(b, end of b) · e^{end(g, b)}
```

(all terms in nats), and the backward pass mirrors them. Every inner sum is a
prefix or suffix sum, as every maximum of §4 is a running maximum. On an edge
`b → e`, every entry of `e` beyond `b`'s end sees the same exit mass, so
forward adds it once to a difference array and backward reads it from a
suffix sum: an edge costs `O(overlap + 1)`.

**Arithmetic.** Values are linear with one log-scale per `(node, group)`
slot, as in scaled HMM recursions. Cut weights are running products along the
alignment (steps `e^{±(λ_r d + μ_r)}` tabulated for small `d`), and prior
weights come from an `e^{λk}` table, so no transcendental function is
evaluated per cut.

**Reported probabilities.** A segment's `posterior` is the mass of every
placement of its kind that is tied with the call on the called interval, or
whose germline bases between the call's neighbours the primary allele
contains: the probability that the call is germline-equivalent to the read's
allele, the event §7 scores. `d_presence_posterior = 1 − P(no D)`.

**Theorem 2.** The marginals equal exhaustive summation over all grammatical
paths and canonical cut pairs, for VJ, VDJ and VDDJ grammars, a per-base N/P
cost and a per-column evidence term `μ ≠ 0` (unit test in
`adjudicate/posterior.rs`, which also checks that the Viterbi score is the
best path score).

## 5. Search

Each stage narrows the hypotheses the DAG sees. The table states what each
step can lose.

| Step | Where | Guarantee | Parameter and justification |
|---|---|---|---|
| Strand | `orient.rs` | heuristic + fallback | Count read positions whose 11-mer occurs in the V/J germline, per strand. Trusted when the winner exceeds its Poisson null `λ₀ + 3√λ₀` (`λ₀` = k-mers × germline k-mer density); otherwise compare V seed votes on both strands. |
| V retrieval | `retrieve/seed_vote.rs`, `retrieve/frames.rs` | heuristic | Seed votes per gene frame, best diagonal cluster per frame. Keep 3 genes, and up to 16 while votes are within `max(8, 5 %)` of the leader. Every allele of a kept gene becomes a candidate. |
| V band | `retrieve/seed_vote.rs` | bounded | Cluster's diagonal range, shifted by the allele's frame offset, ± 8 diagonals: holds unless an indel longer than 8 escapes the seeds. |
| V scoring | `segment.rs`, `family_bound.rs` | lossless | Score-only DP; a same-length sibling is skipped when `S(rep) + Δ(rep, sibling)` cannot reach the current top-3 cutoff (sound upper bound). |
| V/J keep | `segment.rs` | defined | Trace back until the top 3 alignments with E ≤ 1, plus ties. No length or identity floor: a positive `(2, −3)` local alignment has identity > 3/5 and the E-value implies a length. |
| J | `retrieve/exhaustive.rs`, `window.rs` | bounded | Every J allele, unbanded, in `[min V end − 12, read end)`: a J can reach back into V only through chance edge matches. |
| D windows | `window.rs` | bounded | `[best V end − 12, best J start + 12)` plus the conserved-anchor window `[V Cys104 + 3, J Trp/Phe118)`. |
| D prefilter | `retrieve/lmer_mask.rs` | lossless w.r.t. gates | Shares ≥ 1 L-mer, with `L` derived from the gates (`min_d_match`; 4 with the optional distributed-D gate, by pigeonhole). |
| D floor | `config.rs`, `adjudicate/model.rs` | lossless | Score ≥ max(gate floor, `P_D + 1`): by Lemma 1, a D at or below `P_D` adds ≤ 0, and dropping it leaves its neighbours' cuts valid (exit ≤ entry), so the no-D path is at least as good (requires `ι ≥ 0`). |
| D admission | `admission.rs` | defined | Exact run ≥ `min_d_match` = 4 (3, 5 and 6 are all worse). The note's distributed-evidence alternative is optional and off: it admitted no D that changed a call. |

*Gene frames.* The alleles of a gene are placed in one coordinate system,
`frame position = allele position + offset`. An allele joins its gene's frame
only when ≥ 30 % of its seeds occur in the gene's longest allele, and ≥ 90 %
of those agree on one offset. SNP alleles and 5′/3′-truncated alleles join;
an allele with an internal indel, or a divergent one, gets its own frame. The
seed index is built over frames with identical postings merged. Retrieval
work therefore scales with genes, not alleles, and allele-specific seeds
still vote.

*One alignment per allele and window.* Each D allele contributes its optimal
local alignment within the window. In a wide window that alignment can trade
the junction core for chance matches in a V or J tail, which the cut then
removes. The 12-base reclaim keeps this rare. With 24 bases, D accuracy was
0.4–0.7 points lower.

## 6. Implementation

- **Inter-sequence SIMD.** Each AVX2 16-bit lane runs the scalar Gotoh
  recurrence for a different job. Traceback codes are bit-identical to the
  scalar kernel's, so the same traceback walk serves both.
  - Banded jobs use band-relative coordinates.
  - With uniform profiles, a cell's substitution is a code comparison;
    otherwise the five profile rows are selected by mask, which is exact for
    any substitution model.
- **Stage-wise batching.** Each stage runs over a chunk of 64 reads, giving
  one score batch and a few traceback waves per stage, while every decision
  stays per read. A read's result does not depend on its chunk (tested).
- **Exactness.** The scalar kernel equals a textbook Gotoh; the SIMD kernel
  equals the scalar kernel (scores, end cells, tracebacks); the batched
  pipeline equals the per-read scalar pipeline.

## 7. Evaluation protocol

Reads come from the IgSim generative model (`sim/`): uniform allele choice,
exonuclease trims, N additions, random flanks, gated SHM with single-base
indels, and instrument substitutions. The presets are `easy` (no SHM), `mid`
(IgM-like) and `hard` (IgG-like); half of the reads are reverse-complemented.
Every read base carries its germline origin.

Metrics (`sim/evaluate.rs`):
- **germline-equivalent** allele accuracy: the primary call contains the
  trimmed, pre-mutation germline body, so no read could separate it from the
  truth;
- gene accuracy;
- boundary error;
- D presence by remnant length;
- junction exactness;
- a random-read negative control;
- **calibration** of any reported posterior: reads binned by posterior, the
  observed frequency of the event per bin, the expected calibration error
  (ECE, `Σ_bins n_b·|mean posterior − frequency| / n`) and the Brier score.
  A segment's event is "the primary call is germline-equivalent to the
  truth"; D presence's is "the read has a simulated D", however short.

**Search error vs model error** (`examples/diagnose.rs`). Each read is also
aligned to every allele, unbanded, under the same objective. If that
exhaustive optimum is equivalent to the truth but the call is not, the search
lost it. If the optimum itself is wrong, or ties the call, no caller using
this objective can do better.

## 8. Open problems

1. **Calls from the per-read model.** The Viterbi call uses fixed scores and
   the posterior uses evidence at the read's divergence, so the posterior's
   mode can differ from the call (at low SHM a short exact D can beat a longer
   one with a mismatch). Viterbi with real-valued per-read weights, or
   maximum-posterior-marginal calls, would make them agree.
2. **D prior from the D search's own statistics.** `P_D = ln(N·L) − ln odds`
   counts genes × read positions. A local search actually tries every allele,
   allele offset and read position, so the chance-match burden is closer to
   the Karlin–Altschul `K·m·n` of the D window against the D database (about
   `e^{10}` for human IGH, against `e^{6.7}`); even prior odds understate how
   often IgSim includes a D (99 %). On IgSim the two errors roughly cancel
   (D-presence ECE 0.2–2.2 %). Both should come from data.
3. **Calibrated statistics.** Fit gapped λ and K, and an effective database
   size for near-identical alleles, by simulation.
