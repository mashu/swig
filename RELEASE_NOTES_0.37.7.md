# Swig 0.37.7

## AER-R distributed D evidence

- Experimental AER-R can now retain a long, high-identity D alignment when substitutions interrupt every calling-profile exact-run seed. The alternative evidence rule requires at least 16 aligned substitution columns, 14 matches, 80% identity, and score 19 under the fixed `+2/−3/−13/−1` tuple.
- The supplied complete-KIMDB macaque regression now calls `IGHD6-39*01` over query positions 386–402 with the 17-nt alignment `GGTCTAGAAGCATCTAC` versus `GGTATAGCAGCAGCTAC` under both AER-R calling profiles.
- Rescue is based on the quality of the surviving D hypotheses, not on `d_call` being empty. A short D call therefore cannot block evaluation of a stronger distributed match.
- AER-R first checks its bounded deferred pool with an exact score-only floor. If that remains unresolved, it forces the 3-mer candidate tier across the complete locus-matched D set, scores only new gene/diagonal hypotheses, and traces only score-capable candidates before the unchanged joint V–D–J comparison.
- A deterministic adversarial regression with fifty perfect 7-nt decoys verifies under both calling profiles that short calls cannot exhaust the candidate/traceback budget.

Ordinary AER, RIAT-MP, standard assignment, and their output are unchanged. AER-R remains opt-in as `aer_robust` in browser configuration and `--assigner aer_robust` in `swig-cli`; RIAT-MP remains the default.

## Validation and cost

- Across 14,000 deterministic KIMDB simulations, 26 previously wrong/absent D calls became retained-tract compatible, no previously compatible D call became incompatible, three wrong calls changed to another wrong/absent result, and one compatible call gained a better boundary. These are development/tuning simulations, not held-out validation.
- On 5,000 fully trimmed-D controls, the D-call rate changed from 53.50% in pinned 0.37.6 to 53.52%, one additional chance 14/17 match. The high residual rate is a pre-existing consequence of short exact matches and remains a separate D-presence model-selection limitation.
- Six alternating warmed one-worker runs measured +2.0% on 5,000 clean reads and +1.3% on a difficult 5,000-read SHM/indel cohort versus pinned 0.37.6. The rejected unconditional safety-pool scan was about 44% slower on clean reads and is not shipped.
- The all-reference candidate oracle now has precise naming in code and documentation: it removes candidate-count pruning but retains the production biological evidence gates.

See [`BENCHMARK_AER_ROBUST_0.37.7.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_AER_ROBUST_0.37.7.md) for the complete audit, controls, and reproduction commands.
