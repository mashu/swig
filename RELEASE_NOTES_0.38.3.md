# Swig 0.38.3

## Promoted R-optimized D candidate recovery

- Promotes the former Sensitive-D 4-nt exact D-candidate floor into the AER-R-only R-optimized profile because it improves fair D accuracy on both supplied complete simulations.
- Retains R-optimized's calibrated 12-point joint D-presence cost, relaxed to 10 only for raw D score at least 20 or exact support from two distinct D template sequences.
- The supplied mutation-interrupted `IGHD5-32` regression is now recovered by main R-optimized as well as Sensitive-D.

## Redefined Sensitive-D

- Redefines `sensitive_d` as R-optimized with a fixed 10-point final joint D-presence cost. Candidate search, the 4-nt exact-run gate, distributed-D fallback, V/D/J scoring, candidate depth, and D ambiguity reporting are unchanged.
- Sensitive-D still requires admitted D alignment evidence and a D-bearing complete V–D–J partition that beats no-D; it does not force a D call on every read.
- On the supplied 20k IgM-like simulation it adds 378 D-positive reports and recovers 279 of R-optimized's 575 remaining truth-D misses. On the 20k IgG-like simulation it adds 489 reports and recovers 356 of 778 remaining misses.
- This is an intentional sensitivity/specificity trade-off. Fair D accuracy falls by 0.1959 and 0.1993 percentage points, while no-D false calls increase by 99 and 133 on the IgM-like and IgG-like data, respectively.

## Accuracy and boundary summary

Using `1/K` credit when one of `K` unique predicted sequence classes matches truth, 0.38.3 R-optimized V/D/J accuracy is 94.7154% / 76.8977% / 91.4748% on IgM-like and 87.4435% / 68.8658% / 86.2089% on IgG-like. Sensitive-D is 94.7154% / 76.7019% / 91.4698% and 87.4426% / 68.6664% / 86.1914%, respectively. Tandem-D truths are excluded only from the primary single-D score.

The D decision change does not reinstate the old trimming bug. Rejected D hypotheses remain unable to alter V/J boundaries; only the winning complete partition can allocate overlapping query bases. Relative to 0.38.3 R-optimized, Sensitive-D moves V-end or J-start MAE by at most 0.00065 nt across either complete simulation, and the full ≥5-nt error tails are essentially unchanged.

## Interfaces and ABI

- Browser and CLI version strings are 0.38.3. Both surfaces expose the revised `R-optimized` and `Sensitive-D` meanings, enforce their AER-R requirement, retain the profile in saved/exported configuration, and reuse it for post-analysis query inference.
- Adds `swig_set_aer_r_profile_decision_tuning_v2(penalty, relaxation)` so explicit CLI tuning can preserve either the R-optimized 12→10 decision or Sensitive-D's fixed 10 decision. The original exported decision setters remain available.
- Adds deterministic WASM regressions for profile promotion, the lower post-signal Sensitive-D decision, a no-D decoy, and V/J partition isolation.

## Speed

Three strictly serial one-worker repetitions found no meaningful regression versus 0.38.1 R-optimized. Median CPU time changed from 29.874 seconds to 29.440 seconds for 0.38.3 R-optimized and 29.804 seconds for Sensitive-D on IgM-like (−1.45% and −0.23%). On IgG-like it changed from 52.449 seconds to 52.869 and 52.810 seconds (+0.80% and +0.69%). R-optimized and Sensitive-D are effectively speed-equivalent because the threshold change adds no search or alignment pass.

See [`BENCHMARK_SENSITIVE_D_0.38.3.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_SENSITIVE_D_0.38.3.md) for the scoring rule, threshold sweep, full call/detection tables, endpoint distributions, supplied-fixture results, and repeated timing.
