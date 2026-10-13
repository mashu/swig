# Swig aligner 0.39.2 versus upstream Swig 0.38.11

Upstream [MurrellGroup/swig](https://github.com/MurrellGroup/swig) is at
commit `a08a2c3` (Swig 0.38.11). Its aligner is the SwiftIG C++ core
(`public/swiftig.wasm`) with the AER-R strategy and the R-optimized profile.
This fork replaced that core with the Rust Swig aligner (`engine/`) in 0.39.0,
added the native route in 0.39.1, and fixed the defects an audit found in
0.39.2.

This report covers three things:
- upstream's aligner against 0.39.1 and 0.39.2, on identical reads,
  references and settings;
- the algorithmic differences between the two aligners;
- a third code review of 0.39.2.

All reads are simulated; none of these numbers is external biological
validation.

## Setup

- **Reads.** The frozen suite: the four sets of
  [`BENCHMARK_SWIG_ALIGNER_0.39.0.md`](BENCHMARK_SWIG_ALIGNER_0.39.0.md)
  (IGH mid and hard, BCR mixed, TCR mixed, same seeds) plus three more:
  - IGH easy (3,000 reads);
  - IGH at 30 % SHM (3,000 reads);
  - 1,000 IGK reads run against the IGH database.

  All are made by `swig simulate` with their truth.
- **References.** Human IG and TR V, D, J and C from the bundled IMGT pack,
  written with their `SWIGMETA` as Swig composes them:
  - IGH: 480 V, 48 D, 15 J, 111 C;
  - BCR: 741 V, 48 D, 39 J, 135 C;
  - TCR: 367 V, 6 D, 101 J, 27 C.
- **Settings.** Each version's own Swig Web defaults:
  - SwiftIG: AER-R, R-optimized;
  - Swig: R-optimized;
  - both: minimum identity 0.6, both strands, C database included,
    double-D off.

  One WebAssembly instance (one worker) in Node.js 22, reads passed in
  batches of 1,000, as Swig's workers send them.
- **Modules** (SHA-256 prefixes):
  - SwiftIG 0.38.11: upstream `public/swiftig.wasm`, `98191e3a7a46b7a0`
    (the module of the 0.39.0 report);
  - Swig 0.39.1: `b91371cf56898728`;
  - Swig 0.39.2: `48ba508f7db738cf`.
- **Scoring.** `swig evaluate` (0.39.2) scores each AIRR TSV against the
  truth and compares two callers read by read.
  - An allele is correct when the first reported name is
    germline-equivalent to the simulated one.
  - The junction is correct when it equals the true Cys104…Trp/Phe118
    sequence.
  - Δ is 0.39.2 minus SwiftIG, in percentage points, with its 95 % CI;
    *p* is an exact McNemar test.

  SwiftIG's figures reproduce the 0.39.0 report's.

## Speed

Timing runs are in progress; the next commit fills in this section.

## Accuracy

Paired comparisons on the same reads. Δ and *p* compare 0.39.2 with SwiftIG.

### IGH easy (3,000 reads)

| Metric | Reads | SwiftIG 0.38.11 % | Swig 0.39.1 % | Swig 0.39.2 % | Δ 0.39.2 − SwiftIG (pp) | ±95% CI | *p* |
| −-- | −--: | −--: | −--: | −--: | −--: | −--: | −--: |
| V allele | 3,000 | 99.73 | 99.93 | 99.93 | +0.20 | 0.16 | 0.031 |
| D allele | 2,962 | 99.70 | 99.76 | 99.76 | +0.07 | 0.16 | 0.69 |
| D, 4–7 nt remnant | 5 | 60.00 | 80.00 | 80.00 | +20.00 | 35.06 | 1 |
| D, 8+ nt remnant | 2,957 | 99.76 | 99.80 | 99.80 | +0.03 | 0.15 | 1 |
| No D call (reads without D) | 38 | 100.00 | 100.00 | 100.00 | +0.00 | 0.00 | 1 |
| J allele | 3,000 | 99.97 | 99.90 | 99.90 | −0.07 | 0.09 | 0.5 |
| Junction | 2,761 | 97.36 | 100.00 | 99.75 | +2.39 | 0.62 | 1.9e-15 |
| C allele | 2,370 | 88.19 | 96.24 | 96.24 | +8.06 | 1.10 | 6.4e-58 |
| No C call (reads without C) | 630 | 100.00 | 99.68 | 99.68 | −0.32 | 0.44 | 0.5 |

### IGH mid (5,000 reads)

| Metric | Reads | SwiftIG 0.38.11 % | Swig 0.39.1 % | Swig 0.39.2 % | Δ 0.39.2 − SwiftIG (pp) | ±95% CI | *p* |
| −-- | −--: | −--: | −--: | −--: | −--: | −--: | −--: |
| V allele | 5,000 | 98.54 | 98.60 | 98.62 | +0.08 | 0.12 | 0.34 |
| D allele | 4,755 | 92.74 | 93.52 | 93.52 | +0.78 | 0.42 | 3.9e-4 |
| D, 1–3 nt remnant | 139 | 41.01 | 47.48 | 47.48 | +6.47 | 5.71 | 0.049 |
| D, 4–7 nt remnant | 303 | 48.18 | 54.46 | 54.46 | +6.27 | 4.28 | 6.6e-3 |
| D, 8+ nt remnant | 4,313 | 97.54 | 97.75 | 97.75 | +0.21 | 0.30 | 0.22 |
| No D call (reads without D) | 245 | 82.86 | 75.51 | 75.51 | −7.35 | 3.81 | 2.8e-4 |
| J allele | 5,000 | 98.62 | 98.70 | 98.70 | +0.08 | 0.18 | 0.52 |
| Junction | 3,874 | 97.37 | 99.97 | 99.85 | +2.48 | 0.52 | 2.6e-24 |
| C allele | 4,002 | 88.31 | 96.15 | 96.15 | +7.85 | 0.84 | 1.9e-91 |
| No C call (reads without C) | 998 | 100.00 | 100.00 | 100.00 | +0.00 | 0.00 | 1 |

### IGH hard (5,000 reads)

| Metric | Reads | SwiftIG 0.38.11 % | Swig 0.39.1 % | Swig 0.39.2 % | Δ 0.39.2 − SwiftIG (pp) | ±95% CI | *p* |
| −-- | −--: | −--: | −--: | −--: | −--: | −--: | −--: |
| V allele | 5,000 | 92.14 | 92.60 | 92.66 | +0.52 | 0.30 | 1.1e-3 |
| D allele | 4,741 | 83.25 | 85.38 | 85.40 | +2.15 | 0.60 | 1.1e-12 |
| D, 1–3 nt remnant | 142 | 20.42 | 36.62 | 36.62 | +16.20 | 6.94 | 1.5e-5 |
| D, 4–7 nt remnant | 339 | 35.99 | 42.48 | 42.48 | +6.49 | 3.95 | 2.1e-3 |
| D, 8+ nt remnant | 4,260 | 89.11 | 90.42 | 90.45 | +1.34 | 0.53 | 8.3e-7 |
| No D call (reads without D) | 259 | 91.12 | 80.31 | 80.31 | −10.81 | 3.78 | 7.4e-9 |
| J allele | 5,000 | 92.84 | 92.94 | 92.98 | +0.14 | 0.29 | 0.41 |
| Junction | 2,165 | 94.18 | 99.72 | 99.54 | +5.36 | 0.97 | 1.1e-32 |
| C allele | 3,994 | 88.33 | 96.82 | 96.82 | +8.49 | 0.86 | 1.8e-102 |
| No C call (reads without C) | 1,006 | 100.00 | 100.00 | 100.00 | +0.00 | 0.00 | 1 |

### IGH, 30 % SHM (3,000 reads)

| Metric | Reads | SwiftIG 0.38.11 % | Swig 0.39.1 % | Swig 0.39.2 % | Δ 0.39.2 − SwiftIG (pp) | ±95% CI | *p* |
| −-- | −--: | −--: | −--: | −--: | −--: | −--: | −--: |
| V allele | 3,000 | 58.13 | 62.03 | 63.90 | +5.77 | 1.22 | 1.5e-20 |
| D allele | 2,871 | 33.02 | 41.55 | 41.55 | +8.53 | 1.52 | 6.9e-28 |
| D, 1–3 nt remnant | 90 | 55.56 | 67.78 | 67.78 | +12.22 | 9.15 | 0.019 |
| D, 4–7 nt remnant | 211 | 9.95 | 15.64 | 15.64 | +5.69 | 5.20 | 0.05 |
| D, 8+ nt remnant | 2,570 | 34.12 | 42.76 | 42.76 | +8.64 | 1.61 | 1.5e-25 |
| No D call (reads without D) | 129 | 55.81 | 33.33 | 33.33 | −22.48 | 9.43 | 1.5e-5 |
| J allele | 3,000 | 55.80 | 59.47 | 59.20 | +3.40 | 1.28 | 2.5e-7 |
| Junction | 2,274 | 41.42 | 84.56 | 85.05 | +43.62 | 2.18 | 5.6e-251 |
| C allele | 2,412 | 77.36 | 83.54 | 83.42 | +6.05 | 1.66 | 1.5e-12 |
| No C call (reads without C) | 588 | 100.00 | 99.83 | 99.83 | −0.17 | 0.33 | 1 |

### BCR: IGH, IGK and IGL against one database (7,500 reads)

| Metric | Reads | SwiftIG 0.38.11 % | Swig 0.39.1 % | Swig 0.39.2 % | Δ 0.39.2 − SwiftIG (pp) | ±95% CI | *p* |
| −-- | −--: | −--: | −--: | −--: | −--: | −--: | −--: |
| V allele | 7,500 | 98.53 | 98.59 | 98.59 | +0.05 | 0.07 | 0.29 |
| D allele | 2,374 | 92.42 | 93.30 | 93.30 | +0.88 | 0.56 | 3.1e-3 |
| D, 1–3 nt remnant | 78 | 41.03 | 50.00 | 50.00 | +8.97 | 8.09 | 0.065 |
| D, 4–7 nt remnant | 162 | 46.30 | 52.47 | 52.47 | +6.17 | 5.04 | 0.031 |
| D, 8+ nt remnant | 2,134 | 97.80 | 97.99 | 97.99 | +0.19 | 0.39 | 0.48 |
| No D call (reads without D) | 5,126 | 99.67 | 99.51 | 99.51 | −0.16 | 0.12 | 0.021 |
| J allele | 7,500 | 99.20 | 99.20 | 99.20 | +0.00 | 0.08 | 1 |
| Junction | 4,606 | 95.09 | 99.93 | 99.35 | +4.26 | 0.66 | 5.8e-40 |
| C allele | 5,968 | 90.63 | 97.03 | 97.03 | +6.40 | 0.62 | 1.9e-113 |
| No C call (reads without C) | 1,532 | 100.00 | 99.93 | 99.93 | −0.07 | 0.13 | 1 |

### TCR: TRA and TRB against one database (5,000 reads)

| Metric | Reads | SwiftIG 0.38.11 % | Swig 0.39.1 % | Swig 0.39.2 % | Δ 0.39.2 − SwiftIG (pp) | ±95% CI | *p* |
| −-- | −--: | −--: | −--: | −--: | −--: | −--: | −--: |
| V allele | 5,000 | 99.90 | 99.92 | 99.92 | +0.02 | 0.09 | 1 |
| D allele | 2,095 | 77.14 | 92.98 | 92.98 | +15.85 | 1.57 | 1.9e-98 |
| D, 1–3 nt remnant | 260 | 6.15 | 56.92 | 56.92 | +50.77 | 6.08 | 3.7e-40 |
| D, 4–7 nt remnant | 515 | 54.56 | 93.20 | 93.20 | +38.64 | 4.24 | 1.3e-58 |
| D, 8+ nt remnant | 1,320 | 99.92 | 100.00 | 100.00 | +0.08 | 0.15 | 1 |
| No D call (reads without D) | 2,905 | 99.79 | 96.97 | 96.97 | −2.82 | 0.60 | 4.1e-25 |
| J allele | 5,000 | 100.00 | 100.00 | 100.00 | +0.00 | 0.00 | 1 |
| Junction | 4,534 | 89.92 | 100.00 | 97.31 | +7.39 | 1.01 | 7.3e-47 |
| C allele | 3,983 | 93.07 | 96.33 | 96.33 | +3.26 | 0.55 | 1.5e-39 |
| No C call (reads without C) | 1,017 | 100.00 | 100.00 | 100.00 | +0.00 | 0.00 | 1 |

### Junctions: coverage, not precision

SwiftIG writes no junction unless the first base of each conserved anchor
lies inside its V and J alignments. Swig projects the anchors along the
alignments. The table splits the reads with a true junction by outcome,
excluding reads of the corrected-anchor alleles (next section):

| Set | Reads | SwiftIG correct / wrong / none | 0.39.2 correct / wrong / none | Exact among reported: SwiftIG / 0.39.2 |
| --- | ---: | --- | --- | --- |
| IGH easy | 2,754 | 97.39 / 0.00 / 2.61 | 100.00 / 0.00 / 0.00 | 100.00 / 100.00 |
| IGH mid | 3,869 | 97.36 / 0.00 / 2.64 | 99.97 / 0.00 / 0.03 | 100.00 / 100.00 |
| IGH hard | 2,161 | 94.26 / 0.00 / 5.74 | 99.72 / 0.09 / 0.19 | 100.00 / 99.91 |
| IGH 30 % SHM | 2,268 | 41.49 / 0.35 / 58.16 | 85.27 / 3.84 / 10.89 | 99.16 / 95.70 |
| BCR | 4,579 | 95.09 / 0.02 / 4.89 | 99.93 / 0.04 / 0.02 | 99.98 / 99.96 |
| TCR | 4,412 | 89.69 / 0.00 / 10.31 | 100.00 / 0.00 / 0.00 | 100.00 / 100.00 |

- **The gain is coverage.** Almost every junction SwiftIG gets "wrong" is
  one it does not report. When it reports one, it is as exact as Swig's, and
  more exact at 30 % SHM (99.2 % against 95.7 %).
- **Where Swig differs.** It reports a junction on almost every read that
  has one: 89 % of them at 30 % SHM, against 42 % for SwiftIG.

### Corrected anchors and the truth

Swig 0.39.2 replaces misplaced `SWIGMETA` anchors with the allele's own
motif hit. That covers 50 J allele records (human TRAJ24\*01–03 and
TRAJ48\*01 among them) and 8 V alleles. SwiftIG accepts the same misplaced
anchors.

The frozen truths were simulated from the misplaced anchors. So on reads of
those alleles, 0.39.2's corrected junctions score as wrong:
- TCR: 122 reads;
- BCR: 27 reads;
- IGH: 4–7 reads per set.

With those reads excluded, 0.39.1 and 0.39.2 junctions are identical on
every set except 30 % SHM, where 0.39.2 is +0.66 pp.

### Reads of the wrong locus

The test ran 1,000 IGK reads against the IGH database. SwiftIG calls an IGH
J on 919 of them; Swig 0.39.2 makes no call.

| Version | V | D | J | Junction |
| --- | ---: | ---: | ---: | ---: |
| SwiftIG 0.38.11 | 19 | 0 | 919 | 0 |
| Swig 0.39.1 | 0 | 5 | 5 | 0 |
| Swig 0.39.2 | 0 | 0 | 0 | 0 |

With the `swig` binary's own default (no identity floor), 0.39.1 gave 72 V
calls and 59 junctions on these reads; 0.39.2 gives none.

### Reading the results

- **V** is the same within noise on ordinary sets. It is +0.5 pp on hard
  reads and +5.8 pp at 30 % SHM, where 0.39.2 adds 1.9 pp over 0.39.1.
- **D** improves on every set: +0.1 to +2.2 pp on IGH and BCR, +8.5 pp at
  30 % SHM, and +15.9 pp on TCR.
  - Most of the gain is on short remnants (1–7 nt), where a correct call
    partly rests on chance matches.
  - On 8+ nt remnants the gain is significant only on hard reads and at
    30 % SHM.
- **The price is more D calls on reads whose D was fully deleted.**
  - IGH: −7.4 and −10.8 pp.
  - 30 % SHM: −22.5 pp.
  - TRB: 88 of 405 D-free reads get a D, against 6 for SwiftIG. The TRB D
    prior is about 6, so one exact 4-mer (score 8) clears it.
- **J** is the same, except +3.4 pp at 30 % SHM.
- **C** is +3.3 to +8.5 pp. C is anchored at the J exit and called by an
  anchored E-value, not a 90 %-identity rule, so short exact C tracts are
  called. False C calls on C-free reads: at most 0.32 pp.
- **Junctions:** see above. The gain is coverage.

## What changed algorithmically

| Stage | SwiftIG 0.38.11 (upstream) | Swig aligner 0.39.2 (this fork) |
| --- | --- | --- |
| Germline annotation | Anchors and regions only from `SWIGMETA`, prepared by Swig's TypeScript. No junction when an anchor is missing; no motif fallback. Locus from a substring anywhere in the header. | `SWIGMETA` and header keys first, then IMGT gaps, joint Cys23/Trp41/Cys104 and `[WF]G.G` motifs, relatives and database models, per locus. Contradicted `SWIGMETA` anchors are corrected. Locus from `locus=` or the gene name. |
| Orientation | Seed votes on both strands; when not lopsided, annotates both strands fully and keeps the higher joint score. | Decided once: 11-mer germline presence, then the strongest V seed cluster. Cheaper, but 3.3 % misoriented at 30 % SHM (open item). |
| V retrieval | 9-mer seeds (7-mers on weak reads), one candidate per allele, top 3 alleles extended to 16 within max(8, 5 %) votes. | Gene frames: alleles of a gene share coordinates. Genes are ranked: 3, up to 16 within z·√(2·4·L) votes. Every allele of a kept gene is a candidate. Lossless sibling bounds. |
| Alignment | Gotoh local, R-optimized tuples, band ≥ 28 around the best seed bin, 4-lane SIMD screen. V 5′ and J 3′ extended to the end without a stop. | Same tuples (C uses V's). Band = seed cluster ± 8, spanning a seeded indel's far side. J: every allele, unbanded. D: every allele sharing a 4-mer. Anchored DP for C. 16/32-lane inter-sequence SIMD, bit-identical to scalar. X-drop flank extension. |
| D admission and prior | Exact run ≥ 4 or distributed evidence; a D costs 12 (10 when strong), for every locus. | Exact run ≥ 4 as the stage's acceptance rule (distributed gate off). P*_D = ln(N·25)/λ per locus (10.55 on IGH, about 6 on TRB), read at the read's own divergence. A J pays P*_J too. |
| Boundaries and junction | V and J fixed independently. For each V×J pair and D inside it, two partitions: D takes the overlap, or V/J keep it. Plus: a +24 pair bonus, 24/10-base floors, a geometry check and a 225-nt junction cap. | One exact placement DAG (semi-Markov Viterbi) over topology, alleles and every pair of cuts, with N/P free. A call needs a V or J clip with E ≤ 10⁻³. |
| Ties and alternatives | No D first, then a V restore rule, longer D run, alphabetical names. Wide D lists of exact-substring alleles. | Value, then fewer D, then support. Co-optimal names need equal matches and mismatches; best-supported first. |
| Statistics | Offline-fitted gapped E-values for five tuples. `v_support` is empty under the default R-optimized profile. No posteriors, no gate. | Exact ungapped λ, K for any tuple; E-values gate hits (E ≤ 1) and calls (anchor E ≤ 10⁻³). Posteriors for V, D, J, C and D presence by forward–backward over the same DAG. |
| C | Post-J seeded local alignment; identity ≥ 0.90 and an E-value route. | Anchored at the J exit (±3), anchored E-value per locus, `c_posterior`. |
| Multi-locus | V and J across all loci (top 3 V, top 2 J); a pair needs one locus; one D cost for all. | J and D restricted to the anchoring placements' loci; per-locus adjudication, D/J priors and C gates. |
| Double-D screen | Seeded screen between V end and J start. | The same rules and sidecar. |

SwiftIG safeguards that this fork does not carry over:
- **Removed:** the pair geometry check, the 225-nt junction cap, and the
  24/10-base V/J floors.
- **Replaced:** two-strand annotation on ambiguous reads, the minimum band
  of 28, the distributed-D rule (off), and the ≥ 30-base ordinary C route.

The fork measured several of these as having no effect on its simulated
reads, but they guard against chimeras, long adapters and unseeded long
indels, which the simulator does not produce.

## What this fork changed

- **0.39.0** replaced the SwiftIG core with the Rust Swig aligner. This
  brought the placement DAG, posteriors, gene-frame retrieval, inter-sequence
  SIMD, motif and IMGT annotation, multi-locus databases, the anchored C
  call, and one calling profile family instead of four strategies.
- **0.39.1** made `swig-cli` run the native aligner (statically linked) with
  WebAssembly as the fallback. WebAssembly workers share one compiled module.
- **0.39.2** fixed what an audit found (see
  [`RELEASE_NOTES_0.39.2.md`](RELEASE_NOTES_0.39.2.md)):
  - misplaced reference anchors;
  - the anchor gate judged on the kept clip;
  - a J presence cost;
  - J search after chance V hits;
  - the 64-allele V cap;
  - gene-frame indel alleles and long V indels;
  - outer-flank X-drop;
  - extended-score keep;
  - D admission as acceptance;
  - one E-value formula;
  - input handling;
  - removal of the residual adjudicator and the tandem-D grammar.

## Third review of 0.39.2

Five reviewers covered the input layer, retrieval, alignment, adjudication
and the pipeline/interfaces. A separate verifier tried to refute each
finding; 22 of 23 were confirmed.

| Severity | Finding | Effect | Fix |
| --- | --- | --- | --- |
| Major | An AIRR table whose header starts with a tab (pandas `to_csv` with its default index, or a CSV with an index column converted by `swig-cli`) is read one column off: `reads.rs` `detect()` eats the header's leading tab. | Every read gets no call, with exit 0, on all routes. | Keep the first line's indentation in `detect()`. |
| Major | D windows are still bounded by the best V, anchoring or not (`window.rs` `residual()` and `WithAnchorWindow`). | On reads that start after V, a chance V empties the D window. On such cut reads, D allele +10.9 pp (209/0) and D-presence Brier 0.109 → 0.005 with the fix. | Take the D geometry from an anchoring V, as the J window and locus filter already do. |
| Minor | `--min-identity` (0.6 in Swig Web and `swig-cli`) is measured on the flank-extended alignment, after the keep. | Correct high-SHM J calls are dropped for the flank the aligner stitched on itself. Testing the local alignment instead recovers J +0.2 pp at 30 % SHM (floor 0.6) and +0.7 pp on hard reads (floor 0.8). | Test the floor on the local alignment. |
| Minor | The extended-score keep counts identical allele copies toward `top`, and hits at the read end pay no flank. | 3 J calls lost at 30 % SHM. | Keep the union of the top by local and by extended score. |
| Minor | J-less (or V-less) reads: the D window runs to the read's end while P*_D is priced for 25 positions. | 64–97 % of reads cut after V with a 25–100-base random tail get a D. | Bound the open side; price D for the window searched. |
| Minor | V pays no presence cost. | About 2 % of J-anchored reads with a random 5′ prefix get a junk V (E up to 10⁸). | Price V like J. |
| Minor | A read holding two loci (scFv) gets one J request for all loci. | One chain's J alleles crowd out the other's: 20 of 800 scFv reads lose the light J. | One J request per anchoring locus. |
| Minor | The cross-locus choice uses log Z under per-subset models. | The called chain can flip with `--no-posterior` (17 of 800 scFv reads). | Choose on the model the call maximises. |
| Minor | Gene identity comes only from `*` in allele names. | Without `*` the 16-gene budget becomes a 16-allele cap: V −8.2 pp at 30 % SHM, with no warning. | Warn; document GENE\*ALLELE. |
| Minor | `--calls fixed`: a tied placement with a shorter span becomes the report's primary. | The row drops read bases and names another gene. | Require ties to cover the called span. |
| Minor | Germline FASTA ignores the reads' BOM skip and UTF-8 tolerance. | Files Swig Web accepts fail in the binary and `swig-cli`. | Decode lossily; strip the BOM. |
| Minor | `swig prepare` does not persist IMGT before-start slots or 3′-partial status. | A reloaded 5′-partial allele gains `fr1=0`. | Write and read markers for both. |
| Minor | `ImgtGaps::check` fails when no gapped allele holds a codon (0 of 0). | The part loses its exact numbering. | Skip codons with no data. |
| Minor | Swig Web names unnamed reads per worker (audit B2, untracked). | Duplicate `sequence_N` ids. | Name records before batching, as `swig-cli` does. |
| Minor | Debug invariants never run in CI, and one property test breaks the stage's precondition. | Assurance only. | Run a debug-assertion test pass in CI; fix the test. |
| Info | `write_batches` drops the valid reads before a malformed record in the same batch. | The run still fails, naming the line. | Write them before failing, or document it. |
| Info | IGHT/IGHZ/IGHO constant alleles get no locus. | No shipped scope is affected. | Add the isotype letters. |
| Info | `join_indel` thresholds use mixed 9-/7-mer units after the fallback. | Output unchanged on the suite. | Count primary seeds. |
| Info | ALGORITHM.md still says C records over about 2 kb leave the SIMD lanes; the bound is now about 15 kb. | Documentation only. | Update the text. |
| Info | `--c-prefix-identity` and `--double-d` are silently ignored without the database they need. | No effect. | Warn. |
| Info | `swig_version` / `AlignerEngine.version()` have no caller and mishandle a stopped instance. | No effect. | Remove them or return the status. |

None of these findings changes the comparison above on the frozen suite,
which has no cut reads, index columns, scFv reads or unnamed records.

## Limitations

- All reads come from this fork's simulator (`swig simulate`), on which its
  priors were measured. There is no real-data comparison.
- The junction gain is coverage of reported junctions, not more precise
  boundaries.
- The D gain is mostly on short remnants. It comes with more D calls on
  D-free reads, notably on TRB.
- Timings are from one 4-core container and vary between runs; ratios are
  more stable than absolute rates.
