# Swig 0.37.0

## Hyper-optimized AER and RIAT-MP kernels

- AER and RIAT-MP now evaluate independent affine alignments four at a time with exact WebAssembly SIMD128 recurrences and deterministic, reference-equivalent tracebacks.
- Exact score-only SIMD screens wide safety pools; full traceback is restricted to candidates that can still reach the requested output cutoff, including ties and invalid-hit fallthrough.
- Production seed search replaces repeated hash lookups with direct-address CSR k-mer hits, generation-stamped vote/reduction workspaces, and a linear two-pass per-gene reduction.
- Orientation candidate rankings are reused by the selected strand rather than recomputed.
- RIAT-MP batches root alignments, prelinks tree topology once, and uses parent-linked trace arenas for near-optimal multipath enumeration.
- AIRR rows append directly to the output buffer without constructing roughly 100 temporary strings per record.

The original scalar/hash/ostream implementations remain available as reference paths. A permanent end-to-end WASM test requires byte-identical reference and optimized AIRR for both engines.

On the documented one-worker 10,000-read release-WASM control, AER throughput increased 2.20x and RIAT-MP throughput increased 1.77x; peak RSS also fell by roughly 13–14%. See [`BENCHMARK_SWIFTIG_KERNELS_0.37.0.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_SWIFTIG_KERNELS_0.37.0.md) for the call-level profiles, protocol, timings, exact output hashes, and limitations.
