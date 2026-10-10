# Swig 0.13.2

## Optional IgBLAST-agreement calling profile

Analysis parameters now include an explicit **Calling profile** selector. The existing ground-truth-optimized settings remain the default. The optional IgBLAST-agreement profile changes only D/J scoring, exact-run floors, and candidate depths; it was selected without consulting the supplied truth columns.

On the complete supplied 100,000-record human IGH tuning set, D first-call/fair agreement with IgBLAST rises from 83.339%/83.328% to 98.825%/98.902%. J rises from 99.298%/99.320% to 99.361%/99.382%. The profile, scoring rules, ambiguity accounting, hashes, and limitations are documented in [`BENCHMARK_IGBLAST_PROFILE_0.13.2.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_IGBLAST_PROFILE_0.13.2.md).

The selected profile is passed to every parallel WASM worker, displayed in the run manifest and results header, stored in saved sessions, restored with older sessions defaulting safely to the truth-optimized profile, and reused when post-analysis queries auto-infer V/J constraints.

## Regression protection

- The dual-profile production WASM resets deterministically between profiles and rejects unknown profile identifiers.
- Default calls are byte-for-byte identical to the prior production core on the deterministic 5,004-record regression subset.
- Profile switching is covered by the WASM tests alongside FASTA/FASTQ/AIRR, BCR/TCR, KIMDB annotations, and double-D isolation.
