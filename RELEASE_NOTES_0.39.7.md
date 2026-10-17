# Swig 0.39.7

Swig 0.39.7 makes D and J calls more accurate at no measurable cost in
speed. It changes two things:
- **D priors.** The placement DAG, which chooses V, D and J and their cut
  points jointly, now weighs how many germline bases a D loses at each end
  and how long the N/P regions around it are
  (`engine/src/adjudicate/`).
- **Positional J names.** A J allele tied with the call names it only when
  it explains the read at the called allele's own germline coordinates
  (`engine/src/adjudicate/ties.rs`).

On reads simulated by [GenAIRR](https://github.com/MuteJester/GenAIRR) and
scored against its ambiguity-aware truth, human IGH D calls rise by
0.30–3.12 points, IGH J by 0.24–0.36 and TRB D by 0.79–0.87. Swig now leads
IgBLAST 1.22 and MiXCR 3.0.13 on D and J on every IG set.
`BENCHMARK_GENAIRR_0.39.6.md` gives the benchmark, the scoring and the
analysis that led to both changes.

The aligner runs in Swig Web, in `swig-cli` (native and WebAssembly) and in
the `swig` binary. The changes apply to all three under the R-optimized
(default) and Sensitive-D profiles.

## What changes for users

- **D trim prior.** A D clip now pays for every germline base it trims
  from its allele's ends:
  - `ln(6/5)` = 0.18 nats per base at the 5′ end, a geometric trim of
    mean 5;
  - `ln(8/7)` = 0.13 nats per base at the 3′ end, a geometric trim of
    mean 7.

  Before, a D placement's cuts cost only their alignment evidence. A short
  exact match deep inside a long allele was then worth as much as an
  untrimmed short allele: a mutated true D (`ACTACAATAA`) lost to another
  gene's exact 7-mer (`CTACAAT`) inside it. Choosing the wrong gene inside
  the right D tract was most of swig's D loss: 4.7 points at 5 % SHM and
  10.0 at 10 %. With the priors it is 3.0 and 6.8.

- **N/P length prior around a D.** With geometric N/P lengths of mean 8,
  every base a D keeps is one N/P base fewer, worth `ln(9/8)` = 0.12 nats.
  The D cost is read as `P*_D/m` matched columns, each now worth that much
  more. An untrimmed exact run therefore wins or loses as before at every
  divergence (ALGORITHM §4, *Values*).

- **The D cost is credited with a typical D's trim cost.** The trim prior
  alone made a D harder to call: swig declined a D on 420 of 10,000
  GenAIRR reads at 10 % SHM, against 100. So the D cost is lowered by what
  a D with the prior's mean trims pays, `5·ln(6/5) + 7·ln(8/7)` = 1.85
  nats. The trims then re-rank Ds without changing how often a typical D
  is called, and a chance hit deep inside a long allele pays more than
  before. Several credits were measured on IgSim (see *Validation*). This
  one keeps the D-presence Brier score of 0.39.6 on IGH, and gives the
  most right D calls.

- **Both priors stay inside the joint placement.** Each term depends on a
  single cut, through its germline or read position. A clip therefore
  stays worth `exit(x) − entry(y)`, and the DAG and the posterior keep
  their recursions, which are exact. The V and J cuts are still chosen
  jointly with the D's, and the posteriors read the same terms. The
  terms are added to each D placement's cut values once, when the
  placement is laid out. The recursions do no extra work, and V and J
  columns do no extra work at all.

- **Positional J names.** IGHJ4\*02 and IGHJ5\*02 trimmed to the 3′ part
  they share align the same read bases equally well, from different
  germline offsets. Swig named both; it now names the one the DAG chose by
  germline support. The other stays a tied placement and appears in
  `j_alternatives` at the call's score with its own coordinates. Positional
  D or V names were measured too: D lost 0.05–0.08 points and V did not
  change, so only J is positional.

- **Options.** The `swig` binary, and the WebAssembly module's option
  vocabulary, take `--d-trim M5,M3|off` (default `5,7`) and
  `--np-length M|off` (default `8`). The Truth-optimized and
  IgBLAST-agreement profiles keep both off: they have no D-presence cost
  for the credit to lower. Swig Web and `swig-cli` use the profile's
  setting.

## What does not change

- V calls change by at most 0.01 points on GenAIRR. The IGK and IGL
  scores on GenAIRR do not change.
- The fixed-scale call (`--calls fixed`) has no D priors: their terms are
  fractions of a score unit. Its posterior reads them, as it reads the
  evidence at the read's divergence.
- `swig --d-trim off --np-length off` gives 0.39.6's rows except J names:
  on 10,000 GenAIRR IGH reads at 5 % SHM, 108 rows differ, in `j_call` and
  `j_alternatives` only.
  With `--d-open-flank off` as well, the rows on the CLI tests' fixtures
  are 0.39.4's, byte for byte (`engine/tests/cli.rs`).
- Speed. On 1,000 IGH reads at 10 % SHM, one thread, the binary executes
  0.6 % more instructions than 0.39.6 (callgrind, database preparation
  included). With the priors off it executes 0.1 % more. Wall time is
  within run-to-run noise: 6.12 s (0.39.6) against 6.22 s (0.39.7) per
  10,000 reads, mean of four interleaved single-thread runs, ±8 % between
  runs.

## Validation

**GenAIRR.** GenAIRR 2.3.1 simulates the reads, 10,000 per set, with the
same seed and reference for every tool. The table gives the mean Jaccard
score of each segment's call against GenAIRR's tie set, in percent; both
empty scores 1. Bold marks the best tool on each set.

| Set | Tool | V | D | J |
|---|---|---:|---:|---:|
| IGH naive | IgBLAST 1.22 | 99.86 | 97.80 | 99.48 |
| | MiXCR 3.0.13, all top ties | 99.76 | 94.56 | 99.53 |
| | swig 0.39.6 | **99.98** | 98.17 | 99.52 |
| | swig 0.39.7 | **99.98** | **98.47** | **99.76** |
| IGH 2 % SHM | IgBLAST 1.22 | 99.80 | 96.55 | 99.36 |
| | swig 0.39.6 | **99.85** | 96.38 | 99.39 |
| | swig 0.39.7 | **99.85** | **97.28** | **99.70** |
| IGH 5 % SHM | IgBLAST 1.22 | 99.64 | 93.24 | 99.14 |
| | swig 0.39.6 | **99.68** | 93.16 | 99.15 |
| | swig 0.39.7 | 99.67 | **95.42** | **99.47** |
| IGH 10 % SHM | IgBLAST 1.22 | 99.36 | 84.28 | 98.48 |
| | swig 0.39.6 | 99.42 | 87.12 | 98.70 |
| | swig 0.39.7 | **99.43** | **90.24** | **99.06** |
| IGH realistic | IgBLAST 1.22 | 99.53 | 92.56 | 99.05 |
| | swig 0.39.6 | **99.61** | 92.78 | 99.08 |
| | swig 0.39.7 | 99.60 | **94.96** | **99.44** |
| TRB naive | IgBLAST 1.22 | **99.97** | 98.17 | 99.93 |
| | MiXCR 3.0.13, all top ties | 99.94 | 94.10 | **100.00** |
| | swig 0.39.7 | **99.97** | **99.91** | 99.89 |
| TRB realistic | IgBLAST 1.22 | **99.94** | 97.65 | **99.92** |
| | swig 0.39.7 | **99.94** | **99.50** | 99.87 |

The light-chain sets (IGK and IGL, 5 % SHM and realistic) are unchanged:
swig leads on V and J there, as in 0.39.6. On TRB J swig still trails
IgBLAST and MiXCR by 0.04–0.11 points, because of extreme J trims below
its 10-column J floor; 0.39.7 does not change that.

**IgSim.** The same 5,000-read sets per preset as 0.39.x's assessment,
scored by `swig evaluate`:
- *D* is the right D up to germline equivalence;
- *D-free given a D* is the share of reads without a D that get one;
- *Brier* is the D-presence posterior's Brier score;
- *J listed* is the share of reads whose true J allele is among the names.

| Set | D 0.39.6 → 0.39.7 | D-free given a D | D-presence Brier | J listed |
|---|---|---|---|---|
| human IGH mid | 93.01 → 93.01 | 24.05 → 19.41 | 0.0229 → 0.0234 | 99.06 → 98.96 |
| human IGH hard | 85.79 → 86.74 | 16.40 → 16.40 | 0.0218 → 0.0236 | 95.60 → 95.06 |
| macaque (KIMDB) IGH mid | 92.37 → 93.17 | 19.38 → 19.38 | 0.0221 → 0.0220 | 98.58 → 98.62 |
| macaque (KIMDB) IGH hard | 83.40 → 84.35 | 14.03 → 12.23 | 0.0245 → 0.0254 | 95.78 → 95.68 |
| mouse IGH mid | 88.07 → 88.75 | 23.86 → 23.52 | 0.0540 → 0.0551 | 100.00 → 100.00 |
| Macaca IGH mid | 91.56 → 92.05 | 17.42 → 14.33 | 0.0291 → 0.0273 | 99.92 → 99.92 |
| alpaca IGH mid | 93.19 → 93.89 | 19.31 → 22.28 | 0.0354 → 0.0354 | 99.88 → 99.86 |
| human TRB mid | 91.89 → 92.64 | 22.82 → 26.91 | 0.0638 → 0.0658 | 99.88 → 99.88 |
| human TRB hard | 83.24 → 84.32 | 12.90 → 16.62 | 0.0818 → 0.0824 | 99.28 → 99.28 |
| human IGK mid | – | – | – | 96.38 → 88.96 |

V and J accuracy up to germline equivalence do not change on any set.

**The credit.** The D cost credit sets the operating point between D
sensitivity and D-free specificity. Less credit trades right D calls for
fewer D calls on D-free reads (human IGH mid, D accuracy at
D-free-given-a-D):
- no credit: 90.57 % at 5.49 %;
- half the credit: 92.08 % at 14.77 %;
- the full credit (0.39.7): 93.01 % at 19.41 %;
- 0.39.6: 93.01 % at 24.05 %.

With the full credit, 0.39.7 calls as many right Ds as 0.39.6 on this set
with fewer D-free reads given a D, and more right Ds on every other IGH
set.

A credit per allele length, the prior's expected trim cost for the
allele, was also tried. It left TRB's D-free rate unchanged, but cut IGH
D accuracy by up to 1.2 points and worsened the IGH D-presence Brier
score, so 0.39.7 uses the single credit.

**Exactness.** The DAG and the posterior with the D terms were checked
against exhaustive enumeration over every path and cut pair:
- `per_read_dag_equals_exhaustive_enumeration` covers pipeline placements
  of mid and hypermutated reads, with and without a per-base N/P cost;
- `d_terms_keep_the_call_a_posterior_mode` covers the marginals and the
  Viterbi key on fixtures whose Ds trim germline at both ends, one with a
  deletion that steps the 5′ trim by two;
- `j_ties_are_named_at_the_first_alleles_germline_coordinates` covers the
  J naming rule: an IGHJ4\*02/IGHJ5\*02 pair names one J, and the same
  pair as D names both.

The test changes:
- `the_open_flank_cost_can_be_turned_off` now pins 0.39.4's bytes with the
  priors off.
- `max_marginal_d_keeps_the_topology_and_raises_the_d_mass` runs with the
  priors on and off.
- The WebAssembly test `every J allele is aligned...` accepts a true J
  that is a tied alternative from other germline coordinates (sim_0000292:
  IGHJ5-4 behind the called IGHJ4-3).

## Known trade-offs

- **TRB and alpaca D-free reads get a D more often.** On TRB this rises by
  3.7–4.1 points (26.91 % mid, 16.62 % hard), and on the alpaca set by
  3.0. TRB's two D genes are 12 and 16 bases long, so a true TRB D trims
  less than the prior expects, and the credit overpays it. Per-allele
  credits remove this on TRB but cost IGH D accuracy (above).
- **Fewer J names.** The true J allele is listed less often on IgSim:
  −7.4 points on human IGK, −0.1 to −0.5 on IGH. IgSim's truth counts
  these alleles as germline-equivalent: the first name's accuracy does not
  change, and the dropped alleles are in `j_alternatives`.
- **The trim means are generic.** The means 5 and 7 suit human IGH. They
  were not fitted to GenAIRR or IgSim; GenAIRR's own trim tables did no
  better in the prototype.

## Reproduce

```sh
cd engine && cargo test --release           # exactness tests, oracles
cd .. && npm test                           # WebAssembly and CLI suites
# GenAIRR benchmark: benchmarks/genairr/README.md
swig --v V.fasta --d D.fasta --j J.fasta reads.fasta                          # 0.39.7
swig --v V.fasta --d D.fasta --j J.fasta --d-trim off --np-length off reads.fasta   # 0.39.6's calls
```
