# Swig aligner 0.39.0 versus SwiftIG 0.38.11

Swig 0.39.0 replaces the SwiftIG C++ core (`public/swiftig.wasm`, AER-R
strategy with the R-optimized profile, the 0.38.11 default) with the Rust
Swig aligner (`engine/`, `public/swig.wasm`, R-optimized profile, the 0.39.0
default). This report compares the two WebAssembly modules on identical
simulated reads, references and settings. All reads are simulated; none of
these numbers is external biological validation.

## Setup

- **References.** Human IG and TR V, D, J and C alleles from the bundled
  IMGT pack (`public/references/imgt-202632-7-swig-0.7.json.gz`), written as
  FASTA with their `SWIGMETA` metadata exactly as Swig Web composes them.
  IGH: 480 V, 48 D, 15 J, 111 C. BCR database: IGH + IGK + IGL (741 V,
  48 D, 39 J, 135 C). TCR database: TRA + TRB + TRD + TRG (367 V, 6 D,
  101 J, 27 C).
- **Reads.** `swig simulate` (the engine's IgSim model: allele usage,
  deletions, P/N additions, hotspot-weighted SHM, indels, sequencing
  errors, partial reads, flanks, half reverse-complemented, 80% of reads
  through the J–C splice into a C allele), with the truth of every read:

  | Set | Reads | Command |
  | --- | ---: | --- |
  | IGH mid | 5,000 | `swig simulate --v IGH/V.fasta --d IGH/D.fasta --j IGH/J.fasta --c IGH/C.fasta --preset mid --reads 5000 --seed 11` |
  | IGH hard | 5,000 | as above, `--preset hard --seed 12` (higher SHM and indel rates) |
  | BCR mixed | 7,500 | 2,500 each of IGH mid (seed 11), IGK mid and IGL mid (seed 13), against the BCR database |
  | TCR mixed | 5,000 | 2,500 each of TRB and TRA, `--preset mid --shm 0` (seeds 21, 22), against the TCR database |

- **Settings.** Both strands, minimum identity 0.6, C database included,
  default calling profile, double-D screen off. One WebAssembly instance
  (one worker) in Node.js 22.22, reads passed in batches of 1,000 as Swig's
  workers do; initialization (database indexing) is timed separately. Every
  read is distinct, so neither module's repeated-sequence cache applies.
  Host: 4-core x86-64 Linux container.
- **Modules.** SwiftIG 0.38.11: `public/swiftig.wasm` at commit
  `a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f` (SHA-256 prefix
  `98191e3a7a46b7a0`). Swig aligner 0.39.0: `scripts/build-swig-wasm.sh`
  (SHA-256 prefix `85ecc859d703f2d0`).
- **Scoring.** `swig evaluate --truth SET.truth.tsv swiftig.tsv swig.tsv`
  scores any AIRR TSV against the simulated truth and compares the two
  callers read by read. An allele is correct when the first reported name is
  germline-equivalent to the simulated allele on that read (its sequence
  contains the trimmed, pre-mutation germline that was simulated), so luck
  among identical alleles is not rewarded. The junction is correct when the
  called junction equals the true Cys104…Trp/Phe118 sequence (reads with
  both anchor codons templated). C is scored against anchored equivalence.
  Δ is the second caller minus the first in percentage points with its 95%
  confidence interval; *p* is an exact McNemar test on the reads only one
  caller gets right.

## Speed

Reads per second for one worker, median of three runs:

| Set | SwiftIG 0.38.11 | Swig aligner 0.39.0 | Speed-up |
| --- | ---: | ---: | ---: |
| IGH mid | 405 | 1,080 | 2.7× |
| IGH hard | 327 | 1,024 | 3.1× |
| BCR mixed | 417 | 1,350 | 3.2× |
| TCR mixed | 757 | 1,280 | 1.7× |

Database initialization per worker takes 40–70 ms for SwiftIG and
80–160 ms for the Swig aligner, which annotates every allele from IMGT
gaps, motifs and relatives on load; it is paid once per worker.

`swig-cli --vdj` with 4 workers on 10,000 IGH reads (the IGH mid and hard
sets; V, D, J and 107 C alleles), wall time including start-up and
reference preparation: 8.6 s with 0.38.11, 3.3 s with 0.39.0 (2.6×). The
native `swig` binary (`engine/`, attached to releases as
`swig-aligner-*`) assigns the same 10,000 reads with all 111 C alleles in
1.8 s on 4 threads (5,900 reads/s) and 5.9 s on one.

## Accuracy

Paired comparisons (first: SwiftIG 0.38.11; second: Swig aligner 0.39.0).

### IGH mid (5,000 reads)

| Metric | Reads | SwiftIG % | Swig % | Δ pp | ±95% CI | *p* |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| V allele | 5,000 | 98.54 | 98.60 | +0.06 | 0.12 | 0.51 |
| D allele | 4,755 | 92.74 | 93.52 | +0.78 | 0.42 | 3.9e-4 |
| D, 1–3 nt remnant | 139 | 41.01 | 47.48 | +6.47 | 5.71 | 0.049 |
| D, 4–7 nt remnant | 303 | 48.18 | 54.46 | +6.27 | 4.28 | 6.6e-3 |
| D, 8+ nt remnant | 4,313 | 97.54 | 97.75 | +0.21 | 0.30 | 0.22 |
| No D call (reads without D) | 245 | 82.86 | 75.51 | −7.35 | 3.81 | 2.8e-4 |
| J allele | 5,000 | 98.62 | 98.70 | +0.08 | 0.18 | 0.52 |
| Junction | 3,874 | 97.37 | 99.97 | +2.61 | 0.50 | 7.9e-31 |
| C allele | 4,002 | 88.31 | 96.15 | +7.85 | 0.84 | 1.9e-91 |
| No C call (reads without C) | 998 | 100.00 | 100.00 | 0.00 | 0.00 | 1 |

### IGH hard (5,000 reads)

| Metric | Reads | SwiftIG % | Swig % | Δ pp | ±95% CI | *p* |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| V allele | 5,000 | 92.14 | 92.60 | +0.46 | 0.31 | 4.4e-3 |
| D allele | 4,741 | 83.13 | 85.26 | +2.13 | 0.60 | 2.2e-12 |
| D, 1–3 nt remnant | 142 | 20.42 | 36.62 | +16.20 | 6.94 | 1.5e-5 |
| D, 4–7 nt remnant | 340 | 35.88 | 42.35 | +6.47 | 3.93 | 2.1e-3 |
| D, 8+ nt remnant | 4,259 | 88.99 | 90.30 | +1.31 | 0.53 | 1.4e-6 |
| No D call (reads without D) | 259 | 91.12 | 80.31 | −10.81 | 3.78 | 7.5e-9 |
| J allele | 5,000 | 92.84 | 92.94 | +0.10 | 0.29 | 0.59 |
| Junction | 2,165 | 94.18 | 99.72 | +5.54 | 0.96 | 1.5e-36 |
| C allele | 3,994 | 88.33 | 96.82 | +8.49 | 0.86 | 1.8e-102 |
| No C call (reads without C) | 1,006 | 100.00 | 100.00 | 0.00 | 0.00 | 1 |

### BCR mixed: IGH, IGK and IGL against one database (7,500 reads)

| Metric | Reads | SwiftIG % | Swig % | Δ pp | ±95% CI | *p* |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| V allele | 7,500 | 98.53 | 98.59 | +0.05 | 0.07 | 0.29 |
| D allele | 2,374 | 92.42 | 93.30 | +0.88 | 0.56 | 3.1e-3 |
| No D call (reads without D) | 5,126 | 99.67 | 99.51 | −0.16 | 0.12 | 0.02 |
| J allele | 7,500 | 99.20 | 99.20 | 0.00 | 0.08 | 1 |
| Junction | 4,606 | 95.09 | 99.93 | +4.84 | 0.62 | 1.5e-67 |
| C allele | 5,968 | 90.63 | 97.03 | +6.40 | 0.62 | 2.0e-113 |
| No C call (reads without C) | 1,532 | 100.00 | 99.93 | −0.07 | 0.13 | 1 |

Neither caller placed a D on an IGK or IGL read or reported a locus other
than the simulated one.

### TCR mixed: TRB and TRA against one database (5,000 reads)

| Metric | Reads | SwiftIG % | Swig % | Δ pp | ±95% CI | *p* |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| V allele | 5,000 | 99.90 | 99.92 | +0.02 | 0.09 | 1 |
| D allele | 2,095 | 77.04 | 92.94 | +15.89 | 1.57 | 9.6e-99 |
| D, 1–3 nt remnant | 260 | 6.15 | 56.92 | +50.77 | 6.08 | 3.7e-40 |
| D, 4–7 nt remnant | 516 | 54.46 | 93.02 | +38.57 | 4.23 | 1.3e-58 |
| D, 8+ nt remnant | 1,319 | 99.85 | 100.00 | +0.15 | 0.21 | 0.5 |
| No D call (reads without D) | 2,905 | 99.79 | 96.97 | −2.82 | 0.60 | 4.1e-25 |
| J allele | 5,000 | 100.00 | 100.00 | 0.00 | 0.00 | 1 |
| Junction | 4,534 | 89.92 | 100.00 | +10.08 | 0.88 | 5.4e-138 |
| C allele | 3,983 | 93.07 | 96.33 | +3.26 | 0.55 | 1.5e-39 |

No TRA read received a D from either caller; the 82 additional D calls on
reads without D are TRB reads whose D was fully trimmed (88 of 405 for the
Swig aligner, 6 for SwiftIG).

### Reading the results

- **V and J** are unchanged within noise on every set, and V improves on the
  hard set (+0.46 pp).
- **Junctions** are exact on 99.7–100% of reads (SwiftIG: 89.9–97.4%): the
  placement DAG chooses both cut points of each junction together with the
  alleles, instead of fixing V and J first.
- **D** improves on every set, most on short remnants and on TRB, whose
  few short D genes the per-locus D prior prices correctly. The price is
  more D calls on reads whose D was fully trimmed (−2.8 to −10.8 pp on that
  minority of reads); the reported `d_presence_posterior` quantifies that
  uncertainty per read.
- **C** improves by 3.3–8.5 pp: C is anchored at the called J's end and
  called by an anchored E-value instead of a 90% identity and length rule,
  so short exact C tracts after J are called.

## Double-D screen

The screen (`public/methods/04_DOUBLE_D_SCREEN.md`) was compared with
SwiftIG's AER-R screen, the 0.38.11 default, on the development simulator
used by `tests/swig-wasm.test.mjs` (`tests/vdj-simulator.mjs`, rhesus macaque
KIMDB 1.1 IGH, 20% tandem-D reads, forward reads, default evidence options,
mode `all`):

| Seed | Reads (tandem) | SwiftIG found | Swig found | Both D right (SwiftIG / Swig) | Single-D false positives (SwiftIG / Swig) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 404 | 1,000 (196) | 112 | 113 | 92 / 92 | 2 / 2 |
| 7 | 1,000 (205) | 117 | 118 | 94 / 95 | 5 / 4 |
| 99 | 1,000 (202) | 135 | 135 | 112 / 112 | 5 / 5 |

The screen's evidence rules are SwiftIG's; it runs between the Swig
aligner's V end and J start, and with the shorter rescue seed always on. The
screened run took 1.0–1.1 s against SwiftIG's 1.6–1.9 s.

## Limitations

- All reads are simulated by one model (`swig simulate`), which also comes
  from the Swig aligner's code base; its SHM, indel and trimming
  distributions are not those of every repertoire. The double-D cohorts use
  the separate development simulator in `tests/`.
- The no-D trade-off is real: on reads whose D was deleted entirely, the
  R-optimized prior more often reports a short chance D match, especially
  on TRB.
- Timings are from one 4-core container and vary by about ±5% between
  runs; ratios are more stable than absolute rates.
