# Swig 0.13.3

## IgBLAST-balanced calling profile

Analysis parameters now include an optional **IgBLAST-balanced** profile. It maximizes agreement with the supplied IgBLAST calls subject to two constraints on the same 100,000-record simulation: mean V/D/J first-call truth accuracy and mean V/D/J ambiguity-aware truth accuracy must each exceed IgBLAST's corresponding score.

The profile starts from the 0.13.2 IgBLAST-agreement D/J settings. It removes a D call only when the selected D alignment has exactly a five-nucleotide longest exact run and `j_sequence_start - v_sequence_end <= 11`. Longer exact runs and five-base hits in longer V-J spans are unchanged. The rule affected 2,168 of 100,000 tuning records.

When the rule removes a D call, Swig also clears D scores, coordinates, CIGAR, frame and alternatives; combines the intervening sequence into `np1`; clears `np2`; and rebuilds the stitched nucleotide and amino-acid V-J alignments. The double-D evidence table retains its two-D result while its `standard_d_call` field is reconciled with the AIRR row.

The optimized 0.13.2 WebAssembly core is unchanged. The balanced evidence rule is applied batchwise in the compute worker, preserving the existing multi-worker and streaming path. On the deterministic 5,004-record benchmark, the complete consistency transformation reduced single-worker benchmark throughput by approximately 7% relative to the agreement-only profile.

The truth-optimized profile remains the default and is unchanged. The agreement-only profile also remains independently selectable and unchanged.

See [`BENCHMARK_IGBLAST_BALANCED_0.13.3.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_IGBLAST_BALANCED_0.13.3.md) for scoring definitions, complete results and limitations.
