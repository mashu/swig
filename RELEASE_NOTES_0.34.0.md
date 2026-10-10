# Swig 0.34.0 release notes

## Calibrated AIRR segment support

- Added standard AIRR `v_support`, `d_support`, `j_support`, and `c_support` columns to assignments produced by the shared SwiftIG core. This covers Swig Web, full `swig-cli` pipelines, and streaming `swig-cli --vdj` output.
- Values are BLAST-form expectation values calibrated to SwiftIG's exact affine scores, actual searched query span, and supplied segment-reference search space. They are not posterior probabilities and do not claim numerical identity with IgBLAST.
- Shipped truth-optimized and IgBLAST-agreement/balanced scoring tuples have separate deterministic offline calibrations. Unknown tuning-only score tuples leave support empty instead of reusing incompatible constants.
- The browser record panel shows each populated E-value beside identity; all four fields remain available in raw AIRR detail and downloads.

## Calibration and performance

- Added reproducible fixed-seed null simulation and fitting tools plus held-out calibration checks in [`BENCHMARK_AIRR_SUPPORT_0.34.0.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_AIRR_SUPPORT_0.34.0.md).
- Database totals are cached once per worker. Runtime performs no extra alignment, candidate search, or reference scan.
- Repeated 50,000-read, four-worker human-IGH benchmarks measured a 0.21% mean runtime difference, within run noise. Four new columns increased the benchmark AIRR output by 1.27%.

## Verification

- Added WebAssembly tests for field population, finite values, absent-segment behavior, and reference-search-space scaling.
- Added streaming CLI assertions for support values and absence of a C value when no constant reference was supplied.
