# Swig 0.38.10

This is a result-preserving performance release for AER-R with the R-optimized calling profile.

## Faster assignment

- Uses a shorter but algebraically identical affine-DP recurrence wherever all active SIMD lanes share the same band interior; the original boundary recurrence remains in force at every band edge.
- Hoists reference pointers and lengths out of DP cells.
- Removes redundant seed-hit packing/unpacking and integer division during dense vote reduction.
- Uses exact partial ranking when only the first bounded candidate set can be returned.
- Reuses one annotation engine and its allocation workspaces across worker batches.
- Reuses a previously serialized annotation for an identical nucleotide sequence through a bounded 4,096-entry FIFO cache. Sequence IDs and FASTQ qualities are emitted from the current record, and all cache state is discarded whenever the database or any calling option changes.

No score, threshold, candidate limit, search fallback, alignment tie, boundary decision, AIRR field, worker setting, or default has changed. AER-R and R-optimized remain the assignment defaults.

## Validation

- Both supplied macaque repertoires were benchmarked separately against the complete rhesus KIMDB reference using one worker, 1,000-record batches, AER-R, and R-optimized.
- Paired 5,000-read wall-time throughput improved by approximately 1.40x on both inputs.
- The resulting AIRR files were byte-for-byte identical to 0.38.9.
- The permanent optimized/reference WASM oracle and a new cross-batch FASTQ cache regression require exact output equality.

See [`BENCHMARK_AER_R_SPEED_0.38.10.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_AER_R_SPEED_0.38.10.md) for timings, hashes, and scope.
