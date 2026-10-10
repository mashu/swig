# Swig 0.39.0

Swig's V(D)J aligner is now the Rust Swig aligner in `engine/`, which
replaces the SwiftIG C++ core. Swig Web and `swig-cli` run it as
`public/swig.wasm`, built for `wasm32-unknown-unknown` with SIMD128 and
needing no imports. Releases also attach it as a native command-line program
(`swig-aligner-<platform>`). Every other Swig tool is unchanged: collapse,
chimeras, lineages, SHM, germline inference, alignment and phylogeny.

## What changes for users

- **One search.** The aligner retrieves V by seed votes over gene frames,
  aligns J and D exhaustively in their windows, and chooses three things
  together, exactly, in one placement DAG: the topology (VJ or VDJ), the
  alleles, and both cut points of every junction. Each read is called
  under its own divergence model. The assignment strategies (standard,
  RIAT-MP, AER, AER-R) and their selector are gone. A saved session or
  config that names one loads with it dropped, and `swig-cli --vdj` rejects
  `--assigner`.
- **Calling profiles.** These are R-optimized (default), Sensitive-D,
  Truth-optimized and IgBLAST-agreement. The IgBLAST-balanced profile is
  removed: a saved session that names it opens with R-optimized, and a CLI
  config that names it is rejected. Under
  R-optimized, a D is included only when it beats the no-D path by a prior
  cost derived from the D genes of the read's own locus: 10.55 score units
  for human IGH. Sensitive-D makes a D twice as likely a priori (9.45).
  The profiles are no longer tied to a strategy.
- **Speed.** On simulated human reads, one worker assigns 2.7× (IGH),
  3.2× (mixed BCR) and 1.7× (mixed TCR) more reads per second.
  `swig-cli --vdj` with 4 workers finishes 10,000 IGH reads in 3.3 s
  instead of 8.6 s.
- **Accuracy**, paired against SwiftIG 0.38.11 on the same simulated reads:
  - Junctions are exact on 99.7–100% of reads (SwiftIG: 89.9–97.4%).
  - D alleles are right more often on every set: +0.8 to +2.1 pp for IGH,
    +15.9 pp for TRB.
  - C alleles improve by 3.3–8.5 pp.
  - V and J are unchanged within noise.
  - The price is more D calls on reads whose D was fully trimmed (−2.8 to
    −10.8 pp on those reads).

  See [`BENCHMARK_SWIG_ALIGNER_0.39.0.md`](BENCHMARK_SWIG_ALIGNER_0.39.0.md).
- **Posterior probabilities.** New AIRR columns `v_posterior`,
  `d_posterior`, `j_posterior` and `c_posterior` give the probability that
  the reported names include the read's allele, up to germline equivalence.
  `d_presence_posterior` gives the probability that the read has a D at all.
- **Annotation without metadata.** When a germline record has no `SWIGMETA`,
  the aligner derives FWR/CDR boundaries and anchors from IMGT gaps (where
  the FASTA keeps them), conserved motifs and annotated relatives. As a result,
  `swig-cli --vdj` without annotation data now fills FWR/CDR, junction and
  productivity fields. `v_annotation_source` and `j_annotation_source` name
  where each allele's coordinates came from.
- **Constant region.** C base 1 is anchored within 3 bases of the called J's
  end, in the J's locus. C is called when its anchored E-value over the
  locus' C genes is at most 10⁻³, instead of the earlier ≥ 90% identity and
  ≥ 30 nt (or short-route) rules. Short exact C tracts after J are called,
  and C never changes the V(D)J call. The optional covered-prefix identity
  gate is unchanged.
- **Support values.** `*_support` are Karlin–Altschul E-values with λ and K
  computed exactly for each profile's ungapped scores, not fitted offline.
  The numbers therefore differ from 0.38.x; they remain E-values over the
  supplied segment database.
- **Double-D screen.** The same evidence rules (seeds, pair gain,
  pseudo-tandem veto) and the same 38-column sidecar. It now runs between
  the Swig aligner's V end and J start, and the shorter rescue seed is
  always on (previously AER-R only). On three simulated cohorts it finds as
  many tandem-D reads as SwiftIG's AER-R screen, with the same number of
  false positives or fewer.
- **Output.** The AIRR table keeps Swig's 106 columns in their order and
  adds the five posterior columns. CIGARs use `=`/`X` for matches and
  mismatches. `productive` also requires that neither V nor J has a
  frameshift.

## Gap analysis: the Rust engine (vdj_assign) versus SwiftIG as Swig used it

Before the switch, the Rust engine lacked the following. Each item was added
in this release:

| Swig needed | Added |
| --- | --- |
| One database holding several loci (IGH+IGK+IGL, or the TR loci), with J, D and C kept within the V's locus | `engine/src/locus.rs`; per-locus annotation, D prior and C significance |
| `SWIGMETA` reference metadata (frame, J CDR3 stop, FWR1–FWR3 intervals, provenance) | `engine/src/metadata.rs`; headers win over detection |
| Swig's AIRR schema (106 columns in order, `j_call` at column 13), `*_alternatives`, `quality`, `locus`, `sequence_alignment`/`germline_alignment`, `region_definition`, `v/j_annotation_source` | `engine/src/report.rs` |
| `c_prefix_identity` and its optional gate | `c-prefix-identity` option (0 turns it off) |
| Strand both/forward/reverse and the V/D/J identity floor | `strand`, `min-identity` options |
| FASTQ quality passthrough and AIRR-table input (`sequence` column, TSV or CSV) | `engine/src/reads.rs` |
| A bounded repeated-sequence cache | 4,096 sequences per worker |
| The double-D sidecar | `engine/src/double_d.rs` |
| IgBLAST-style tuning (`-min_D_match`, `-D_penalty`, `-J_penalty`, `-min_J_length`) | `min-d-match`, `d-mismatch`, `j-mismatch`, `min-j-length` options |
| A browser ABI with one row per input record, the header in every batch, and the database summary | `engine/wasm/src/lib.rs`, `src/aligner-engine.ts` |

These SwiftIG features were intentionally not carried over:

- **The four assignment strategies and IgBLAST-balanced.** They are
  superseded by one exact adjudication.
- **SwiftIG's offline-fitted gapped E-value constants.** They are replaced
  by exact ungapped λ and K.
- **The 90%-identity C rule.** It is replaced by the anchored E-value.
- **SwiftIG's AER-R boundary rescues** (joint-boundary D refinement,
  distributed-D rule, J retry, junction-facing V/J repairs). The placement
  DAG evaluates every cut point jointly and needs none of them.

## Removed

- The SwiftIG C++ sources (`wasm/src`, `wasm/include/swiftig`,
  `wasm/UPSTREAM.md`), `scripts/build-swiftig-wasm.sh`,
  `public/swiftig.wasm`, and the E-value calibration tooling.
- The strategy benchmarks and simulators that exercised SwiftIG strategies
  (`tests/benchmark-aer-robustness.mjs`, `benchmark-r-optimized.mjs`,
  `benchmark-simulated-accuracy.mjs` and their helpers).
- The reports that measured SwiftIG strategies and profiles. Older release
  notes link to their last version at commit `a08a2c3`.
- `src/balanced-calling-profile.ts` and its tests.

The Alivibe MSA module keeps its C++/WASI build.

## Build and CI

- `npm run build:swig-wasm` (Rust 1.89+, `wasm32-unknown-unknown`) builds the
  aligner module, and `npm run test:engine` runs the engine's tests.
- The Pages workflow tests the engine before building.
- The release workflow builds `swig.wasm` for the standalone CLI and adds
  native `swig` aligner binaries for Linux x86-64/ARM64, macOS Intel/Apple
  silicon and Windows x86-64.
- `npm run benchmark:50k` uses distinct reads, so the repeated-sequence cache
  does not answer them.
