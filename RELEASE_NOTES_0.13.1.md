# Swig 0.13.1

## D/J accuracy tuning

SwiftIG's browser-core D and J defaults were tuned against the supplied 100,000 simulated human IGH rearrangements and their exact simulation references. D now requires a seven-nucleotide exact run and uses a stronger mismatch penalty. J uses stronger mismatch and gap-open penalties. V scoring, candidate depths, alignment algorithms, AIRR fields, ambiguity merging, and the opt-in double-D path are otherwise unchanged.

On the complete supplied tuning set, first-call accuracy changed from 63.750% to 70.825% for D and from 96.021% to 98.053% for J. The ambiguity-penalizing fair scores changed from 63.844% to 70.882% and from 95.838% to 97.712%, respectively. D multi-call frequency fell rather than rose. See [`BENCHMARK_DJ_0.13.1.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_DJ_0.13.1.md) for scoring definitions, D-positive/D-negative stratification, IgBLAST comparison, hashes, parameter search, limitations, and reproduction scripts.

## Verification

- The tuned production binary matches the independently source-built core on the deterministic 20,017-read confirmation subset.
- The production core shows no throughput regression in the single-runtime checks reported with the benchmark.
- All 85 automated tests pass, including heavy/light BCR, all supported TCR loci, uploaded segment swaps, macaque KIMDB CDR annotation, streaming/storage, post-analysis, denoising, phylogeny, and double-D isolation.
- `npm run build` completes the static GitHub Pages artifact in `dist/`.
