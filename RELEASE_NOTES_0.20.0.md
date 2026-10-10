# Swig 0.20.0

## Assignment strategies

- Added an explicit **Assignment strategy** control under advanced analysis
  parameters. It is independent of the existing calling-profile control.
- Added **AER (Adaptive Exact Refinement)** and made it the default. AER keeps
  exact affine V alignment and expands candidate depth only when the leading
  k-mer ranking is ambiguous.
- Added **RIAT-MP (Root-Indexed Allele Tree, Multipath)**. It indexes close
  V-allele cluster roots, propagates sparse allele changes, and performs no
  descendant V alignments. A bounded second root traceback is considered only
  for selected indel-bearing cases.
- Retained **Standard SwiftIG** as a selectable fixed-depth comparator.
- Retained calibrated exact D and J assignment for every strategy. A structured
  J prototype was tested and is documented, but was not selected for release.

The selected strategy is shown in the run manifest and Results header, saved in
portable sessions and project checkpoints, and reused when post-analysis query
sequences auto-infer V/J constraints. Sessions made before 0.20.0 restore as
Standard SwiftIG so their original assignment semantics are represented
correctly.

## Truth-profile D/J defaults

The compiled truth-profile defaults now match the jointly calibrated supplied
low-SHM and IgG simulations:

- D scoring `+2/-3/-13/-1`, two reported hits, minimum six aligned bases.
- J scoring `+2/-3/-17/-2`, two reported hits, minimum ten aligned bases.

IgBLAST-oriented profiles remain separate explicit options and retain their own
D/J settings.

## Verification

See [`BENCHMARK_ASSIGNERS_0.20.0.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_ASSIGNERS_0.20.0.md) for the
algorithm definitions, deterministic benchmark sample, first-pick and
ambiguity-aware V/D/J results, IgBLAST agreement, throughput, structured-J
experiment, input hashes, and limitations.
