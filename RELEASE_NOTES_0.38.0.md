# Swig 0.38.0

## R-optimized AER-R profile

- Adds an opt-in `R-optimized` calling profile for AER-R in Swig Web, exported pipeline JSON, the general CLI, and streaming `swig-cli --vdj`.
- The profile corrects the audited junction-facing V over-extension mechanism, adds an evidence-conditioned short-D model-selection cost, and reports a narrowly calibrated same-span, one-point D uncertainty set. It uses the existing candidate search and SIMD affine kernels.
- The profile is rejected unless the selected assigner is `aer_robust`. Switching away from AER-R in the browser safely returns the profile to `truth_optimized`.
- Legacy AER-R, AER, RIAT-MP, standard assignment, and all three previous calling profiles remain available and unchanged. RIAT-MP plus truth-optimized calling remains the default.

## Validation

- Uses a sequence-equivalence-aware multiclass Brier score for ambiguous calls and separate endpoint losses; tandem-D truths are not incorrectly scored as one-class D calls.
- On a 3,844-record clean test subset excluding every previously inspected target case, V/D/J Brier improves by 1.23%, 10.54%, and 4.98%. V-end MAE improves `0.539 → 0.427 nt`; J-start MAE improves `0.595 → 0.468 nt`.
- D-presence sensitivity improves by +0.574 percentage points and no-D specificity by +12.155 points. The previous AER-R profile remains available.
- On 760 dramatic prior V-stop failures, exact endpoints improve `4.21% → 75.66%`, within-one accuracy `10.39% → 93.09%`, and MAE `4.518 → 0.376 nt`.
- All 43 previous J misses remain recovered, with 100% specificity among their 42 no-D truths. On the 93 prior D misses, truth-in-set improves `18.28% → 22.58%` and Brier loss improves `0.800 → 0.678`. The strong perfect-D and distributed-D regression fixtures remain called.
- Four alternating same-process one-worker validation runs per profile show no slowdown: median 7.133 s legacy versus 6.369 s R-optimized. Run-to-run noise is large, so the apparent speedup is not treated as a performance claim.

See [`BENCHMARK_R_OPTIMIZED_0.38.0.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_R_OPTIMIZED_0.38.0.md) for the scoring rule, frozen settings, split/exclusion details, tradeoffs, regression cohorts, and reproduction command.
