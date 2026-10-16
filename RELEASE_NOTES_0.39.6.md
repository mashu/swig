# Swig 0.39.6

Swig 0.39.6 brings Swig 0.38.12 of the upstream repository
([MurrellGroup/swig](https://github.com/MurrellGroup/swig), commits
`eb95362` and `6efdc45`) to the Rust aligner. 0.38.12 changed SwiftIG's C++
core, which 0.39.0 replaced, so its fix is ported rather than merged:
- a V whose local alignment stops a few bases short of an intact Cys104
  codon gets its junction, CDR3 and productivity back
  (`engine/src/read_regions.rs`);
- Swig Web shows an empty `productive` field as **Productivity
  unresolved** instead of calling the record non-productive;
- `swig-cli --version` is tested against `package.json`.

No call, score, alignment, segment boundary, posterior or NP length
changes. On 80,000 simulated reads, 167 reads gain a junction, and no
other read changes in any column.

The aligner runs in Swig Web, in `swig-cli` (native and WebAssembly) and in
the `swig` binary, and the recovery applies to all three, given the same
reference and calling profile.

## What changes for users

- **A clipped V keeps its Cys104.** Clustered substitutions just before
  Cys104 can end a V's local alignment short of it. Since 0.39.4 the pair
  guards then withheld the junction, though the read holds an intact
  cysteine where the V puts it. Swig's own regression read is an alpaca
  IGHV3S66\*01 read named `Dyl` (`tests/fixtures/clipped-v-anchor.fasta`).
  With the bundled Vicugna pacos IGH reference it now gets the results
  below, on both strands:
  - `productive=T`;
  - `cdr3_start=289` and `cdr3_end=339`;
  - CDR3 `AARFAGIVAGPWDEYNY`.

  Its V call (IGHV3S66\*01), V score (452) and V end (276) do not change.

  The anchor is recovered, at its projection along the called V's end
  diagonal, only when all of these hold:
  - the allele's anchor codon is a cysteine and ends at most 18 germline
    bases past the V's end;
  - the read has an intact `TGT` or `TGC` at that coordinate, ending before
    the first base of the called D or J;
  - every compared base is unambiguous, in the read and in the germline;
  - no gapped alignment of the stretch scores as well or better than the
    co-linear one. The stretch runs from up to 12 gapless aligned V columns
    before the V's end through the codon. The gapped alignments compared
    are anchored at its start and may end anywhere up to 18 read bases past
    the projection. A tie leaves the anchor unresolved.

  No other cysteine is searched for, and the junction need not be in
  frame. The junction, CDR3, `vj_in_frame`, `productive`, `complete_vdj`,
  the reading frame, `d_frame` and the translated columns follow the
  recovered codon. The stop-codon and V/J frameshift vetoes still decide
  `productive`. ALGORITHM §4c (*Clipped Cys104*) states the rule; REVIEW
  item 20 records what it leaves open.

- **The rival check uses the calling profile's V scores**, as SwiftIG's
  check used its profile's. R-optimized (the default) and Sensitive-D open
  a V gap at −9, as SwiftIG's R-optimized and Sensitive-D did.
  Truth-optimized and IgBLAST-agreement open one at −5. On `Dyl` the best
  gapped explanation then scores 18, as the co-linear one does, against 10
  at −9. So `Dyl` is recovered under R-optimized and Sensitive-D and stays
  unresolved under the other two profiles.

- **Productivity unresolved.** A record whose junction cannot be placed has
  an empty `productive` field in AIRR output, as before. Swig Web's record
  detail now labels it **Productivity unresolved** in a neutral tag,
  instead of "Non-productive". The results table gives it a neutral grey
  dot instead of the red one. Explicit `F` calls still read "Non-productive".
  Upstream kept the warning colour on the detail tag; Swig makes it neutral
  to match the table.

- **CLI version check.** `tests/cli-smoke.test.mjs` asserts that
  `swig-cli --version` prints `package.json`'s version. The check covers the
  repository bundle and, through `SWIG_CLI_EXECUTABLE`, the standalone
  executables. Upstream added it when its 0.38.12 branch was found to
  report 0.38.11. `swig-cli` 0.39.5 reported its version correctly; the
  test keeps it so.

## What does not change

- **Calls and alignments.** V, D, J and C calls, scores, CIGARs,
  coordinates, support and posteriors, alternatives and NP lengths are
  unchanged. The recovery runs after the call, only for a read whose
  called V does not align Cys104's first base.
- **Reads that already had a junction** are unchanged, column for column.
- **Without the pair guards** (`--no-pair-guards`), anchors are projected
  along the alignment however far, as before, so the recovery does not
  apply.
- **J side.** A J aligned only after Trp/Phe118 still gives no junction.
- **Mutated codons.** A mutated cysteine codon is not recovered, nor is a V
  that stops more than 18 germline bases before the codon's end.
- **Repertoire selection.** Its "Productive: no" option still keeps rows
  whose productivity is unresolved, as in 0.39.5. Only the results-table
  filter, which matches `F` exactly, separates them.

## Validation

**Regression read and vetoes** (`tests/swig-wasm.test.mjs`, the WebAssembly
module Swig Web runs). These are upstream's cases on `Dyl`:

| Case | Expected | Result |
|---|---|---|
| `Dyl` and its reverse complement, R-optimized and Sensitive-D | the CDR3 above, `productive=T`; V alignment unchanged | ✓ |
| The same reads, Truth-optimized and IgBLAST-agreement | V alignment unchanged, CDR3 and `productive` empty (the −5 gap tie) | ✓ |
| Stop codon in FR2 | CDR3 recovered, `stop_codon=T`, `productive=F` | ✓ |
| One base inserted in the junction | CDR3 recovered, `vj_in_frame=F`, `productive=F` | ✓ |
| One base inserted in the V | CDR3 recovered, `v_frameshift=T`, `productive=F` | ✓ |
| Cys104 deleted, `TGT` placed later in the junction | V end 276, no CDR3, `productive` empty | ✓ |
| `N` in the clipped tail | no CDR3, `productive` empty | ✓ |
| One base inserted before the codon | no CDR3, `productive` empty | ✓ |
| A cysteine at the projection with a competing gapped mapping | no CDR3, `productive` empty | ✓ |

Against the 0.39.5 module, the first two tests fail and the last two pass,
as they should. The new tests in `read_regions.rs` build clipped reads from
a synthetic database and check the DP's tie rule directly. Five mutants
were each caught by a failing test:
- the recovery disabled;
- no cysteine check on the read;
- no ambiguity check;
- no 18-base tail limit;
- a strict instead of a non-strict rival comparison.

**Paired before/after** (`node tests/benchmark-junction-anchor.mjs
OLD_SWIG NEW_SWIG`). This runs native 0.39.5 and 0.39.6 release builds on
5,000 reads per set, simulated by IgSim (`swig simulate`, seed 765127,
`mid` and `hard` presets). Both builds use R-optimized with Swig Web's 0.6
identity floor. The germline sets are human IGH, IGK, IGL and TRB, alpaca
IGH, macaque (AG07107) IGH, KIMDB 1.1 rhesus IGH and mouse (BALB/c) IGH.
Across the 80,000 reads, 167 rows changed:
- every changed row is a read that had no junction and now has one;
- only junction-derived columns changed;
- no read that had a junction changed, and none lost one.

| Set | Preset | Recovered | Templated, equal | Templated, different | Rebuilt codon | Productive |
|---|---|---:|---:|---:|---:|---:|
| human IGH | mid | 0 | 0 | 0 | 0 | 0 |
| human IGH | hard | 13 | 2 | 0 | 11 | 0 |
| human IGK | mid | 0 | 0 | 0 | 0 | 0 |
| human IGK | hard | 1 | 1 | 0 | 0 | 0 |
| human IGL | mid | 0 | 0 | 0 | 0 | 0 |
| human IGL | hard | 0 | 0 | 0 | 0 | 0 |
| human TRB | mid | 0 | 0 | 0 | 0 | 0 |
| human TRB | hard | 14 | 5 | 0 | 9 | 2 |
| alpaca IGH | mid | 51 | 0 | 0 | 51 | 3 |
| alpaca IGH | hard | 19 | 1 | 0 | 18 | 0 |
| macaque IGH | mid | 4 | 2 | 0 | 2 | 0 |
| macaque IGH | hard | 14 | 4 | 0 | 10 | 0 |
| KIMDB rhesus IGH | mid | 6 | 2 | 0 | 4 | 1 |
| KIMDB rhesus IGH | hard | 15 | 2 | 0 | 13 | 0 |
| mouse IGH | mid | 9 | 1 | 0 | 8 | 0 |
| mouse IGH | hard | 21 | 3 | 0 | 18 | 2 |
| **all** | | **167** | **23** | **0** | **144** | **8** |

"Templated" reads have a junction whose two anchor codons the simulator
drew from germline. All 23 such recoveries equal the simulated junction.
The other 144 hold, at the called V's projected Cys104 coordinate, a
cysteine codon rebuilt by junction bases that are not V's:
- 116 wholly, the V having been trimmed before the codon;
- 28 partly, with 1 or 2 of its bases from V.

Their V call is the simulated allele, or one equivalent to it, on 141 of
the 144. The coordinate is structural, as upstream says: an observed
cysteine does not prove that V supplied its bases. Of the 167, 50 are in
frame and 8 are productive, all 8 with a rebuilt codon. The other 159 keep
`productive=F`, through a stop codon (139) or a frame or frameshift veto.

`swig evaluate` measures the effect on junction accuracy:
- **exact** is the share of reads with a templated junction whose called
  junction is right;
- **spurious** is the share of reads without a templated junction that get
  one anyway.

| Set | Exact, mid | Exact, hard | Spurious, mid | Spurious, hard |
|---|---|---|---|---|
| human IGH | 99.29 → 99.29 % | 96.27 → 96.36 % | 82.57 → 82.57 % | 21.63 → 22.05 % |
| human IGK | 99.05 → 99.05 % | 95.00 → 95.04 % | 6.10 → 6.10 % | 7.62 → 7.62 % |
| human IGL | 99.25 → 99.25 % | 94.53 → 94.53 % | 5.91 → 5.91 % | 5.66 → 5.66 % |
| human TRB | 99.88 → 99.88 % | 97.66 → 97.80 % | 71.43 → 71.43 % | 36.22 → 36.85 % |
| alpaca IGH | 98.63 → 98.63 % | 92.91 → 93.14 % | 23.76 → 24.95 % | 13.02 → 13.41 % |
| macaque IGH | 99.44 → 99.49 % | 96.71 → 96.88 % | 87.27 → 87.50 % | 21.93 → 22.31 % |
| KIMDB rhesus IGH | 99.44 → 99.49 % | 97.03 → 97.11 % | 93.51 → 93.96 % | 19.96 → 20.46 % |
| mouse IGH | 99.54 → 99.57 % | 96.90 → 97.05 % | 60.08 → 60.55 % | 17.44 → 18.03 % |

Exact templated junctions gain up to 0.24 points. Junctions on a codon that
is not wholly templated gain up to 1.2 points. The high spurious rates on
mid reads were there before 0.39.6. On the mid IGH and TRB sets, 89–99.9 %
of those reads have a V that keeps 1 or 2 bases of its trimmed codon and
still aligns Cys104's first base. These are structural checks on simulated
reads, not proof of germline provenance or biological function.

**Speed.** CPU time is unchanged within noise. On 5,000 hard human IGH
reads, one thread, six interleaved runs each, 0.39.6 took 2.54–3.00 s and
0.39.5 took 2.65–2.95 s. The rival DP covers at most 30 × 48 cells, and
runs only when the read has a cysteine codon at the projection.

**Test suites.** All of these pass:
- `npm run test:engine`: the engine's unit and integration tests, with the
  4 new `read_regions` tests;
- `npm test`: the 4 new WebAssembly tests and the CLI version test;
- `npm run build`.

## Differences from SwiftIG 0.38.12

- The comparison scores with the called V allele's substitution profile
  and the profile's gap costs. SwiftIG used the profile's uniform V
  scores. For the shipped profiles the two agree.
- The V alignment is the called V's span, after the placement DAG has
  clipped it. SwiftIG used the V's local alignment. The boundary is the
  first read base of the called D or J span. SwiftIG used the start of
  the D or J alignment.
- Upstream's `benchmark-junction-anchor.mjs` drove two SwiftIG WebAssembly
  builds. Swig's version of it compares two native `swig` binaries on the
  sets above.

## Reproduce

```sh
npm run test:engine          # builds engine/target/release/swig (0.39.6)
git worktree add ../swig-0.39.5 72489f6 && \
  cargo build --release --manifest-path ../swig-0.39.5/engine/Cargo.toml --bin swig
node tests/benchmark-junction-anchor.mjs ../swig-0.39.5/engine/target/release/swig \
  engine/target/release/swig 5000
```
