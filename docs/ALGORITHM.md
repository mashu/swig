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
(`flank::extend_outer`). This makes the alignment end-to-end on its outer
side, where the germline must begin (V) or end (J). The substitution scores of
the extension, mismatches included, count towards the ranking.

### 2.2 Paths

A *placement* is an aligned segment hypothesis. A *path* is a sequence of
placements that the grammar allows (`adjudicate/model.rs`), with their
retained bases in read order:

```text
V → D{0..max_d} → J,     partial paths allowed (V…, …J, V…J);  max_d = 1 (VDJ) or 2 (VDDJ)
```

Between two consecutive placements there are two *cuts*: the first keeps
read bases before its *exit* `x`, the second keeps bases from its *entry*
`y ≥ x` on, and the `y − x` bases between them are N/P. Each placement is
then *clipped* to the alignment columns it keeps. Overlapping placements may
therefore split the overlap at any point, or leave part of it to neither.
Only the cuts order two segments of different kinds: a D may precede a J
whose alignment starts before the D's, as long as the J is entered after the
D's exit. Two segments of one kind (tandem D) follow the order of their
alignment starts.

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
  none. Hence `P*_D = (ln(N·L) − ln π_D) / λ` (`config::DPrior::Derived`),
  kept unrounded: 10.55 for human IGH (N = 32, L = 25, π_D = 1), 10.45 and
  10.34 for the human and synthetic databases of the 20,000-read
  measurements (30 and 28 genes), 6.17 for TRB (N = 2). For integer scores
  `S > P*_D ⟺ S > ⌊P*_D⌋`, so the fixed scale charges `P_D = ⌊P*_D⌋` (10 for
  all three heavy-chain sets); rounding to nearest is not faithful when the
  fraction is ≥ ½. The posterior uses `P*_D` itself, and the per-read call
  reads it at the read's scale (§4). The legacy profiles set `P_D = 0`.
- **`P*_D` from the D search** (`DPrior::Search`, `--d-presence-rate π`).
  `N·L` undercounts the hypotheses a local D search tries: every distinct
  allele, every allele offset and every read position. The search expects
  `E(S) = K·m·n_D·e^{−λS}` chance hits scoring `S` or more, with `K` and `λ`
  the lattice Karlin–Altschul parameters of the D scoring (0.408 and 0.634
  for `(2, −3)`), `m = L` read positions and
  `n_D = Σ_a (L_a − k + 1)⁺` the allele offsets at which an exact run of
  `k = min_d_match` bases, the admission gate, can start (802 for the human
  set, 659 for the synthetic one). A present D clears the threshold almost
  surely, so the best hit is more likely a D than chance when
  `π > (1 − π)·E(S)`:

  ```text
  P*_D = [ln(K·m·n_D) − ln(π / (1 − π))] / λ
  ```

  `ln(K·m·n_D)` exceeds `ln(N·L)` by 2.4 nats on human IGH; the even odds
  of `Derived` offset it. Against a brute-force count of chance hits
  (random 25-base windows, one best hit per allele, the run gate) the closed
  form is high by 1.30–1.35 at even scores and 1.85 at odd ones: exact runs
  put most chance hits on even scores (`config` tests, and the real D sets).
  The cost is read per read like the default (§4). It is a switch, not the
  default, because the rate is not a property of the search: at the
  simulator's own rate (0.94 on mid reads) `P*_D` falls below 10, so exact
  5-mers win, and the rate at which every measured criterion holds (0.9)
  was found by measuring (ASSESSMENT §9.9). The per-read window, `m = W`,
  measured worse on every set: long junctions carry a D more often, which
  outweighs their larger chance-hit count.
- **`tandem`**: the prior odds of a second D (DJ → DDJ), on top of its own `P_D`.
- **`pair_bonus`**: the prior odds of a complete V…J rearrangement against a
  partial read. It is 0 by default: on full-length reads the J's own clip
  score completes the path, and values 0–48 changed no call.

The adjudicator returns the path maximising, lexicographically:

1. its value: by default the read's own log-odds (§4, *Values*), with
   `--calls fixed` the score above;
2. fewer D segments ("with-D must strictly beat no-D");
3. *support*, the sum of the placements' unclipped scores.

Key 3 separates paths that explain the read equally well. A typical case is
two alleles of different 3′ length at a junction where one base can belong to
either side. The allele with more support admits more cuts that reach the same
score, so it has the larger likelihood summed over cuts. Co-optimal alleles
(equal value on the called interval: an equal clipped score, and under a
per-read call equal matches and mismatches too) are all reported, the
best-supported first. The reported `path_score` is the fixed-scale score
above, whichever model chose the path.

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

Placements are nodes, sorted by segment kind (V, D, J), then read start.
Every transition between kinds goes forward in this order wherever the two
alignments start, and the cuts alone decide whether it is feasible
(`x ≤ y`); two segments of one kind (tandem D) follow the order of their
starts. The paths through a node with a given D count and first segment kind
form a **slot**. A slot holds one value
per canonical entry cut `y` of the node: the best key (score, then germline
support) of a path prefix that enters there. With a linear N/P score
`ι·(y − x)` the recursion is

```text
X(b, x)  = exit_b(x) + bonus(b) + max_{y ≤ x − min(b)} [F(b, y) − entry_b(y)]
F(e, y') = max_b max_{x ≤ y'} [X(b, x) + ι·(y' − x)]          (0 at a path start)
best     = max_b X(b, end of b) + end(g, b)
```

over every exit cut `x` of `b` and every entry cut `y' ≥ x` of `e`. Each term
depends on one cut only (separability, §3), so:

1. *Best entry for an exit.* `max_{y ≤ x − min} [F(b, y) − entry_b(y)]` is a
   running maximum along the node, read at `x − min`.
2. *Best exit for an entry.* `max_{x ≤ y'} [X(b, x) − ι·x]` is a second running
   maximum along `b`, one per successor kind, read at `y'`.
3. *Entries beyond the overlap.* Every entry `y' ≥ end(b)` sees the same best
   exit of `b`, so the edge `b → e` adds it once at `e`'s first entry at or
   after `end(b)`, and a prefix maximum over `e`'s entries (shifted by `ι`
   per base) carries it to the rest.

An edge costs `O(overlap + 1)` and a read
`O(Σ_slots len + Σ_edges (overlap + 1))`. The search runs over **every**
canonical entry of `e`, not only those inside the overlap. An earlier version
tried entries only inside the overlap, arguing that moving an entry further
right "clips a prefix of an optimal local alignment, which scores ≥ 0"
(Lemma 1). That is wrong: Lemma 1 bounds prefixes that start at the
alignment's own start, but the bases between `b`'s end and a later entry are
a middle stretch of `e`'s alignment and can score below zero (two mismatches
right after V's end, as a mutated read has). Entering after them and leaving
them to N/P then wins. Adversarial review found the gap, and a hypermutated
IgSim oracle (below) confirms that real reads hit it. A second review found
that nodes sorted by read start alone let a D precede only a J that starts
after it, so a J alignment reaching back over the D (the "J eats D"
configuration) excluded the path V → D → J with that J. The brute-force
oracles had the same restriction; sorting by kind removes it.

**Values: the read's own model** (`adjudicate/units.rs`). The default call
(`CallModel::PerRead`) maximises each path at the read's own divergence `p̂`,
in the model §4a's posterior reads, so the call is a mode of the reported
distribution. Every term is an integer in units of `1/F` nats, `F = 2¹⁶`:

```text
column(kind, s)  = s·u_r + v_r     u_r = round(F·λ_r(p̂)),  v_r = round(F·μ_r(p̂))   (0 for s = 0)
gap, ι, tandem, pair bonus scoring k:   k·u0,   u0 = round(F·λ)
D cost  c_D(p̂)   = (P*_D / m)·ln 4(1 − p̂)  =  ⌊P*_D·(m·u_D + v_D)/m⌋ units
```

with `(λ_r, μ_r)` the per-column log-odds of §4a and `m` the D match score.
Each kind's affine map is quantised as a whole, so at `p̂ = δ` (`λ_r = λ`,
`μ_r = 0`) every value is the fixed score times `u0`: the fixed scale
(`CallModel::FixedScale`, `--calls fixed`: one unit per score unit, D cost
`⌊P*_D⌋`) is an exact special case. A clip is still worth
`exit(x) − entry(y)`, its prefix values accumulated along the columns
(canonical cuts never split a gap run), so the recursion, its ties and
Theorem 1 carry over to integers unchanged, and every tie is exact.

- *Divergence.* `p̂` is estimated before adjudication from the columns no
  path contests: the best V's match columns before every D and J placement
  and the best J's after every V and D placement, `N` excluded, as
  `(k + ½)/(n + 1)` (the Jeffreys predictive rate), from all columns of the
  two placements when fewer than 50 remain, clamped to `[10⁻⁴, 0.5]`. One
  Viterbi pass follows. The junction's bases are not among those columns,
  so a D cost conditioned on `p̂` does not use them twice.
- *D cost.* `P*_D` is derived as prior odds, but the call uses it as a
  threshold against chance D hits in random N/P sequence, whose null
  distribution does not depend on `p̂` and is dominated by exact runs. A cost
  constant in nats would make specificity depend on divergence: an exact
  5-mer is worth `5·ln 4(1 − p̂) = 6.93` nats at low `p̂`, above
  `λ·P*_D = 6.62`, while the fixed scale rejects it (10 < 10.45). Read as
  `P*_D/m` matched D columns, the cost keeps the fixed scale's admitted
  exact-run classes at every `p̂`: runs of `⌊P*_D/m⌋ = 5` bases lose and runs
  of 6 win, by `min(f, 1 − f)·ln 4(1 − p̂)` with `f` the fraction of `P*_D/m`
  (0.22 and 0.17 here). The design prototypes measured the alternative: with
  the constant cost `λ·P_D` (6.34 nats) the per-read call gave 12–22 pp more
  D-free reads a D, most of them exact 5-mers.
- *D floor.* For a per-read call the floor `P_D + 1` bounds nothing: 8
  matches and 2 mismatches score 10 but are worth 6.40 nats at `p̂ = 0.3`,
  above `c_D = 5.38`. Every admissible D is kept, as with posteriors.
- *Co-optimal names* are the placements of equal clipped score *and* equal
  matches and mismatches on the called interval, which makes
  `λ_r·S_sub + μ_r·n + λ·S_gap` equal at every divergence.

Ties: complete paths compare by `PathScore` (value, then fewer D, then
germline support), so a path with a D must strictly beat the best path
without one. Among equal keys: a node keeps its earliest entry, so a junction
base that scores the same on either side goes to the later segment; for that
entry, the predecessor's latest exit; across predecessors, a path start, then
the first predecessor in node order (compared explicitly, so the result does
not depend on the order candidates are folded in).

**Theorem 1.** The DAG returns the optimum of §2.3 over all grammatical
paths, placements and pairs of canonical cuts. Checks:

- an exhaustive, memoised enumeration in `tests/synthetic.rs` over every
  grammatical path, every canonical exit `x` and every canonical entry
  `y ≥ x` (no window), on 40 reads with 0–2 D segments under the tandem-D
  grammar at `ι = −1` (VJ, VDJ and VDDJ optima all occur);
- the same oracle on the pipeline's own placements for 40 hypermutated IgSim
  reads (`dag_equals_brute_force_on_hypermutated_igsim_reads`, top 3 V,
  6 D, 3 J, `ι ∈ {0, −1}`). The test asserts that some optimum enters a
  placement beyond the overlap, and it fails on the earlier restricted search;
- three constructed cases in `adjudicate/dag.rs`: leaving a shared mismatch
  stretch to N beats every single cut; an entry in the middle of J's
  alignment, past V's end, beats every entry inside the overlap; and
  V → D → J beats V → J when the J's alignment starts one base before the D
  (it fails when nodes are sorted by read start alone);
- the per-read model (`per_read_dag_equals_exhaustive_enumeration`): the
  same memoised oracle valued in integer units, on the placements of 80 mid
  and hard IgSim reads under the tandem grammar at `ι ∈ {0, −1}`, at the
  read's own `p̂` and at 10⁻³, 0.02 and 0.3 (640 optima, 135 of them not the
  fixed-scale path), and at `δ` with an integer `P_D`, where the path is the
  fixed call's and its key the fixed score times `u0`; every traced path
  re-scores to its key;
- two constructed per-read cases: an exact 7-mer D (score 14) against a
  10-mer with a mismatch (15), where the fixed scale takes the 10-mer and the
  read's model at `p̂ = 1.2 %` the 7-mer (9.6 against 8.3 nats); and a D at
  the fixed floor that wins at `p̂ = 0.3`.

The exhaustive test on synthetic reads alone separates neither restricted
search from the exact one: its reads never put mismatches right after a
neighbour's end, nor a J alignment over a D. The hypermutated-read test and
the constructed cases are the ones that do.

## 4a. Posterior decoding (`adjudicate/posterior.rs`)

The objective is a log-likelihood ratio, so it also defines a distribution:

```text
P(path, cuts | read) = exp(nats(path, cuts)) / Z,     Z = Σ_{path, cuts} exp(nats)
```

over every grammatical path and every pair of canonical cuts between
neighbours: the domain §4 maximises over. Under the default per-read call
both passes read one integer model (§4, *Values*), `nats = units/F`, so the
Viterbi path of §4 is its mode. Forward–backward over the same DAG gives the
marginal probability of every placement and of every D count.

**Evidence at the read's divergence.** `(2, −3)` is log-odds only at 11.2 %
divergence: a fixed λ makes unmutated reads underconfident (a SNP between
alleles is worth `ln(3(1−p)/p) = 8.0` nats at `p = 0.1 %`, not 3.2) and
heavily mutated ones overconfident. Let `p̂` be the read's divergence,
estimated before adjudication from the uncontested V and J columns (§4,
*Values*). A substitution column scoring
`s ≠ 0` is worth `λ_r·s + μ_r` nats, with `(λ_r, μ_r)` chosen so that a match
is worth `ln 4(1 − p̂)` and a mismatch `ln(4p̂/3)`: the log-odds of
"templated at divergence `p̂`" against "uniform random". An `N` column
(`s = 0`) is worth 0, as it is to the scoring: it supports neither
hypothesis. Gaps stay at the prior λ, because `p̂` measures substitutions, not
indels. Each segment kind uses its own match/mismatch scores. The prior
terms (pair bonus, tandem, `ι`) stay at the λ they were derived with, and a D
costs `c_D(p̂)` (§4). At `p̂ = 11.2 %`, `λ_r = λ` and `μ_r = 0`. The map is
affine per column, so clip scores stay separable in their cuts. With
`--calls fixed` the call keeps the fixed scores while the posterior reads
this evidence unquantised with the D cost `λ·P*_D`; its mode can then differ
from the call.

**Domain.** The first segment keeps its start and the last its end, as in §4
(V and J are extended to the read's ends on purpose). Between neighbours every
pair of canonical cuts `x ≤ y` is summed. A canonical entry lies on the first
retained base and a canonical exit right after the last, so each distinct
retained sub-alignment is one pair and each segmentation is counted once.

**Recursions.** For node `b` in group `g` with forward mass `F(b, y)` per
entry,

```text
X(b, x)  = e^{exit(x) + bonus} · Σ_{y ≤ x − min} F(b, y) · e^{−entry(y)}
F(e, y') = Σ_b Σ_{x ≤ y'} X(b, x) · e^{ι(y' − x)}          (+ 1 at a path start)
Z        = Σ_b X(b, end of b) · e^{end(g, b)}
```

(all terms in nats), the sum-product twin of §4, and the backward pass mirrors
them. Every inner sum is a prefix or suffix sum where §4 has a running
maximum, and entries beyond `b`'s end share one exit mass, so an edge costs
`O(overlap + 1)` here too.

**Arithmetic.** Values are linear, with one log-scale per `(node, group)`
slot, as in scaled HMM recursions.

- *Cut weights* are running products along the alignment (steps
  `e^{±(λ_r d + μ_r)}` tabulated for small `d`; prior weights from an
  `e^{λk}` table), so no transcendental function is evaluated per cut. A
  product that leaves the `f64` range is reseeded exactly at the next cut.
- *Centring.* Entry and exit weights are centred separately, each on the
  midpoint of the evidence range over the cuts a path can actually use on
  that side (entries past the start only if some kind can precede the
  placement, exits before the end only if some kind can follow it). The
  difference of the two centres is carried as a per-node shift. A placement
  whose range on one side exceeds 1 200 nats (`MAX_SPREAD`; weights within
  `e^{±600}`) would leave the `f64` range, and the read then reports no
  posteriors.
- *Junctions.* The N/P factor `e^{ι(y' − x)}` is accumulated as a relative
  recursion (`acc ← acc·e^{ι} + v`), whose factors are ≤ 1, and every prefix
  sum is normalised before it multiplies a cut weight. A strong per-base cost
  across a long V therefore cannot overflow. The cost of a gap between
  neighbours that do not overlap joins the slot's log-scale instead of being
  a factor, so a long, costly gap that every path crosses cannot underflow.
- *Overlaps.* Inside an overlap, the forward mass rises with the
  predecessor's exit evidence and the backward mass falls with the
  successor's entry evidence, so their product is about `e^{−g}` relative to
  its factors' maxima, with `g` the smaller of the two gains across the
  overlap. A read with `g > 600` nats (`OVERLAP_SPREAD`; about 430 matched
  bases, against a few dozen nats for a junction) gets no posteriors rather
  than wrong ones.

**Reported probabilities.** A segment's `posterior` is the mass of every
placement of its kind that is tied with the call on the called interval, or
that explains most of the read bases the call explains (more than half of
the called span, within the region between the call's neighbours) with
germline bases there that the primary allele contains. A placement that only
reaches into the region by a few bases (a D inside the called J) is another
hypothesis: a few bases occur in almost any allele. The event is "the
reported names include the read's allele, up to germline equivalence"; §7
scores exactly this event. Under the tandem grammar a path can hold two such
D placements, so the sum is clamped to 1 and is then an upper bound.
`d_presence_posterior = 1 − P(no D)`.

**Max-marginal D allele** (`adjudicate/decode.rs`, `--d-allele
max-marginal`, off by default). The call is the single best explanation;
its D allele need not be the most probable one. Under 0/1 loss on the D
allele the Bayes act is `a* = argmax_a Σ_{D placements p ∋ a} P(p on the
path)`. D presence stays the Viterbi decision (the presence marginal sums
the look-elsewhere mass of every chance hit), so only a call with exactly one
D is decoded, and a tie within a relative 10⁻⁹ keeps the call's allele.
When `a*` differs, the reported path is re-traced by the same DAG under the
call's model, on the lattice masked to the V and J placements and the D
placements of `a*`, with the grammar restricted to the call's first and last
kind and D count (`Within`). An exhaustive oracle checks this conditioned
pass for every D placement of 40 hypermutated reads. On top of the per-read
call it changes 0.3–1.1 % of rows and gains D allele accuracy significantly
only on human hard reads (+0.25 pp, p 6·10⁻⁵; syn.hard +0.14, p 0.017; mid
sets within noise), with D-free reads untouched; on the fixed-scale call it
gains +0.20 and +0.27 pp (mid, hard). It needs posteriors, so
`--no-posterior` and `--residual` report the call's allele.

**Theorem 2.** The marginals equal exhaustive summation over all grammatical
paths and canonical cut pairs (unit tests in `adjudicate/posterior.rs`):

- VJ, VDJ and VDDJ grammars, with and without a per-base N/P cost;
- evidence at a read's divergence (`μ ≠ 0`);
- a junction with an `N` column and a J with an insertion, at fixed λ and at
  `p̂ = 0.5 %` (the oracle reads `N` columns as 0 and gaps at the prior λ;
  the test fails if `N` columns are scored as `μ`);
- a 300-base V with every 7th base mutated and a strong N/P cost
  (`ι ∈ {−4, −8}`), which overflowed before the junction recursion was made
  relative;
- a 300-base N/P gap at `ι = −4` that every path crosses
  (`require_v`, `require_j`), which underflowed to no posteriors when the gap
  was a linear factor;
- a D before a J whose alignment starts one base earlier, under VDJ and
  VDDJ grammars.

The same tests check that the Viterbi key is the brute-force maximum and that
the traced path re-scores to it, and a 540-base V/J overlap is refused
(`None`) instead of giving marginals above 1. With one integer model for both
passes (`per_read_call_is_a_mode_of_the_posterior`), the marginals equal the
exhaustive sums of `e^{units/F}` and the Viterbi key their maximum, on the
three junction fixtures under VJ, VDJ and VDDJ grammars, with and without an
N/P cost, at `p̂` = 0.5 %, `δ` and 30 %: the call is a mode of the reported
posterior.

## 5. Search

Each stage narrows the hypotheses the DAG sees. The table states what each
step can lose.

| Step | Where | Guarantee | Parameter and justification |
|---|---|---|---|
| Strand | `orient.rs` | heuristic + fallback | Count read positions whose 11-mer occurs in the V/J germline, per strand. Trusted when the winner exceeds its Poisson null `λ₀ + 3√λ₀` (`λ₀` = k-mers × germline k-mer density); otherwise compare V seed votes on both strands. |
| V seeds | `seed.rs` | bounded | Every 9-mer (7-mers as below). `SeedSampling` also offers closed syncmers (a seed is kept when the first smallest `s`-mer of its key, `s = k − w`, in a fixed random order, is first or last; context-free; on a contiguous shape any `w` consecutive k-mers hold one, so every unmutated `k + w − 1` bases yield a shared seed), open syncmers (smallest `s`-mer in the middle, the conservation-optimal offset; context-free, no window guarantee) and minimisers (window guarantee, not context-free). Under the gene margin below, closed syncmers (`w = 4`) are 1.1 points of mean V gene behind every k-mer (4.4 at 30 % SHM) for 0.4 % fewer instructions, so every k-mer is the default. |
| V retrieval | `retrieve/seed_vote.rs`, `retrieve/budget.rs`, `retrieve/frames.rs` | heuristic | Seed votes per gene frame (4 per 9-mer, 1 per 7-mer), best diagonal cluster per frame. 7-mers vote while the leader holds fewer than 25 nine-mer seeds, where a count's Poisson relative error `1/√n` exceeds 20 %. Keep 3 genes, and up to 16 while votes are within `z·√(2·4·L)` of the leader's `L` (z = 2): two standard deviations of the difference between two seed counts near the leader's, so relative to the leader but proportionally wider for weakly supported (mutated) reads. A fixed margin (8 votes, 5 %) or a proportional one (10–30 % of `L`) loses 4–9 points of V gene at 30 % SHM (same seeds and short-seed rule); accuracy saturates at z ≈ 2 (the 16-gene cap binds). Every allele of a kept gene becomes a candidate. |
| V band | `retrieve/seed_vote.rs` | bounded | Cluster's diagonal range, shifted by the allele's frame offset, ± 8 diagonals: holds unless an indel longer than 8 escapes the seeds. |
| V scoring | `segment.rs`, `family_bound.rs` | lossless | Score-only DP on one representative `R` per gene frame (most votes, then the frame's medoid, then lowest id). A sibling `A` in `R`'s frame diagonals satisfies `L ≤ S*(A) ≤ S*(R) + Δ(R, A)`: `Δ` sums the best profile gain per frame position both cover plus the best score where only `A` does; `L` is a traced representative path clipped to `A`'s span and rescored with `A`'s profile. The cutoff `c` is the 3rd largest of one lower bound per candidate, so `c ≤` the exact 3rd score. `A` is skipped when `U < c`, known when `U = L`, aligned outright when `L ≥ c` (its traceback reused), scored otherwise. The ranking from the 3rd score up is that of full scoring; ties follow an explicit key (gene leaders first, then retrieval order), never the schedule. Under a length or identity floor only `U <` floor prunes. |
| J scoring | `segment.rs`, `family_bound.rs` | lossless | The same bounds with each gene one frame, its alleles at offset 0 (sound because J is unbanded); ties in allele order. |
| V/J keep | `segment.rs` | defined | Trace back until the top 3 alignments with E ≤ 1, plus ties. No length or identity floor: a positive `(2, −3)` local alignment has identity > 3/5 and the E-value implies a length. |
| J | `retrieve/exhaustive.rs`, `window.rs` | bounded | Every J allele, unbanded, in `[min V end − 12, read end)`: a J can reach back into V only through chance edge matches. Whole-read J search is 18 % slower and changes a few hard-preset calls (D −0.04 pp, J −0.02 pp), none on mid. The J is not placed jointly with the D, yet it does not "eat" it: the D windows reach 12 bases into J and include the conserved-anchor window, and the DAG orders segments by kind, so a D can precede a J whose alignment reaches back over it and the cuts give any J/D overlap to whichever explains it better. With a D remnant ≥ 8 nt (`examples/d_errors.rs`, human IGH), a missing D whose bases the called J covers is 0.00 % of reads on mid and hard, and a true D outside the windows 0 %; whole-read J lets 0.05 % fall outside. |
| D windows | `window.rs` | bounded | `[best V end − 12, best J start + 12)` plus the conserved-anchor window `[V Cys104 + 3, J Trp/Phe118)`. |
| D prefilter | `retrieve/lmer_mask.rs` | lossless w.r.t. gates | Shares ≥ 1 L-mer, with `L` derived from the gates (`min_d_match`; 4 with the optional distributed-D gate, by pigeonhole). |
| D floor | `config.rs`, `adjudicate/model.rs` | lossless for the call | Only when the DAG computes no posteriors (`--no-posterior`; `--residual` turns posteriors off): score ≥ max(gate floor, `P_D + 1`). By Lemma 1 a D at or below `P_D` adds ≤ 0, and dropping it leaves its neighbours' cuts valid (exit ≤ entry), so the no-D path is at least as good (requires `ι = 0`, the default). With posteriors every admissible D is kept: a D below the floor cannot be called but carries mass `e^{s − P_D}`, and dropping it would bias `d_presence_posterior` down. |
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

- **Inter-sequence SIMD.** Each 16-bit lane runs the scalar Gotoh
  recurrence for a different job. Traceback codes are bit-identical to the
  scalar kernel's, so the same traceback walk serves both. The kernel is
  written once over a lane abstraction and compiled for AVX-512BW (32
  lanes), AVX2, NEON and portable lanes (16), chosen at run time.
  - Banded jobs use band-relative coordinates.
  - With uniform profiles, a cell's substitution is a code comparison;
    otherwise the five profile rows are selected by mask, which is exact for
    any substitution model.
- **Stage-wise batching.** Each stage runs over a chunk of 64 reads, giving
  one score batch and a few traceback waves per stage, while every decision
  stays per read. A read's result does not depend on its chunk (tested).
- **Exactness.** The scalar kernel equals a textbook Gotoh; every SIMD
  backend equals the scalar kernel (scores, end cells, tracebacks); the
  batched pipeline equals the per-read scalar pipeline, and output is
  byte-identical across backends and across x86-64 and AArch64.

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
  A segment's event is "the reported names include the read's allele, up to
  germline equivalence" (§4a); D presence's is "the read has a simulated D",
  however short.

**Comparisons** (`sim::Paired`). Two callers or settings are compared on the
same reads, at least 20,000 per preset: for each metric, the accuracy
difference with its 95 % interval (paired Wald) and an exact two-sided
McNemar test on the discordant reads. D allele accuracy is always reported
with the fraction of D-free reads left without a D call, since either can
be raised at the other's expense, and split by the true remnant's length
(1–3, 4–7, 8+ nt): a remnant shorter than `min_d_match` (4) cannot be found,
so a call counted correct there rests on chance matches. Posteriors are
compared by the paired Brier difference `(p₂ − e₂)² − (p₁ − e₁)²` per read,
with its mean, 95 % interval and t statistic (D allele and D presence).

**D error attribution** (`examples/d_errors.rs`). Every wrong or missing D
call is attributed to where it arose: the true D lies outside the D windows
(recomputed from the call's own V and J placements); inside them, but no
placement of an equivalent allele overlaps it; or such a placement exists and
the path took another D or none. Missing D calls are further split by which
called segment covers the true D's bases (J, V or neither). A J-covered miss
is the failure mode of a J alignment "eating" the D. Counterfactual flags
search J over the whole read (`--j-reclaim 100000`) or drop the anchor
window (`--no-anchor-window`).

**Cross-tool comparison.** IgBLAST 1.22 (Bioconda) on the same IgSim reads
and germline database (`-organism human -ig_seqtype Ig -outfmt 19`, custom
V/D/J databases from the same FASTA), scored by the same evaluator after
mapping IgBLAST's names back. IgBLAST takes J frames and CDR3 ends from its
auxiliary file, which must list every J allele of the database: three were
missing (IGHJ4\*02_g23-, IGHJ5\*02_c33t, IGHJ5\*04) and were added from
`vdj-assign prepare`'s anchors (CDR3 end = anchor − 1, frame = anchor mod 3;
these agree with the file on all 13 shared alleles). Without them IgBLAST
reports no junction for 18 % of reads.

**Search error vs model error** (`examples/diagnose.rs`). Each read is also
aligned to every allele, unbanded, under the same objective. If that
exhaustive optimum is equivalent to the truth but the call is not, the search
lost it. If the optimum itself is wrong, or ties the call, no caller using
this objective can do better.

## 8. Open problems

1. **The D threshold on the read's scale.** Calls come from the per-read
   model (§4), but its D cost `c_D(p̂)` converts `P*_D` so that the fixed
   scale's exact-run classes survive; it is a level-preserving threshold.
   One `p̂` per read also assumes uniform mutation: a D-region rate differing
   from V and J is untested, and the plug-in `p̂` under-weights D remnants
   with two mismatches on unmutated reads against the Beta-binomial
   predictive (by up to `ln 3`).
2. **D prior from data.** The D search's chance-hit count replaces `N·L`
   behind `--d-presence-rate` (§2.3), but the presence rate has to be
   supplied, and the count and the rate are not separately identifiable from
   the calls: every constant cost moves along one curve of D accuracy
   against D-free specificity (ASSESSMENT §9.9). What the posterior's
   residual miscalibration tracks is junction length (short junctions want
   1–2 nats more cost, long ones less), which points to a prior on N/P
   lengths, not to the D search.
3. **Calibrated statistics.** Fit gapped λ and K, and an effective database
   size for near-identical alleles, by simulation.
