# Swig 0.37.6

## Experimental AER-R assignment

- Adds **AER-R**, a new separately selectable experimental assignment strategy (`aer_robust`). RIAT-MP remains the Web/CLI default.
- Keeps ordinary AER, RIAT-MP, and standard SwiftIG intact; pinned-output tests confirm ordinary AER AIRR output is unchanged.
- Jointly scores non-overlapping V–D–J partitions, including the exact affine-score cost of clipping independently optimized V/J endpoints.
- Retains both nominal and selectively expanded conserved-anchor D-location hypotheses, rather than letting one local optimum hide the other.
- Adds a gated complete-J retry for strongly V-supported orientations lacking any valid V/J partition, exact cutoff-tie retention, retained-D-substring ambiguity recovery, and scored J/C overlap splitting.
- Adds a mutation-tolerant shorter-seed rescue to the optional double-D screen only when AER-R is selected; the extension must still reach the configured ordinary seed length.
- Exposes AER-R in the Web selector, portable JSON, generic CLI, and `swig-cli --vdj --assigner aer_robust`.

## Validation and performance

- Adds a deterministic development-only V(D)J simulator with realistic trimming, junction formation, SHM, indels/errors, partial reads, reverse complements, tandem D, and exact provenance.
- Adds a compile-time exhaustive all-reference oracle for diagnosing candidate-pruning versus joint-scoring failures. Neither simulator nor oracle is part of runtime WebAssembly.
- Adds direct regression coverage for the supplied KIMDB 17-nt D case, simulation-level accuracy, double-D rescue, strategy plumbing, and optimized/reference-kernel equivalence.
- On the audited difficult 5,000-record simulation, AER-R improved retained-D-tract compatibility from 88.77% to 89.91% and strong-D detection from 98.17% to 98.88%, while taking 2.7% longer because fallbacks were frequently exercised.
- On an ordinary clean 5,000-read benchmark, AER-R and AER had the same AIRR SHA-256 and their two-run means differed by 0.13%, within timing noise.

## Important limitation

AER-R increased D calls on simulated records retaining fewer than six true D bases from 2.36% to 4.32%. The simulator is a development stress test, not held-out biological validation. AER-R therefore remains opt-in; use ordinary RIAT-MP/AER when preserving released behavior is more important than testing the new robustness tradeoff. Full methods and measurements are in [`BENCHMARK_AER_ROBUST_0.37.6.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_AER_ROBUST_0.37.6.md).
