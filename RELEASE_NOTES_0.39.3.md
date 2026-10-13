# Swig 0.39.3

Swig 0.39.3 fixes the defects that a third review of the Swig aligner, made
after 0.39.2, found: 2 major and 13 minor. Thirteen are fixed outright and
two in part. Each fix has a regression test that fails on 0.39.2. The
aligner runs in Swig Web, in `swig-cli` (native and WebAssembly) and in the
`swig` binary. Every change to calls and AIRR columns below applies to all
three. Changes to reading files, options and exports name the route they
apply to.

On ordinary reads almost nothing changes. On 28,500 simulated IGH, BCR and
TCR reads, with default settings:
- 1,000 IGK reads run against the IGH database give byte-identical
  output;
- the TCR set keeps every call. On 378 reads whose candidates span more
  than one locus, the posteriors move by at most 0.006;
- on the IGH easy, mid and hard and BCR sets, 22 rows change `v_call`.
  21 keep their first name and add full-length alleles tied with it. One
  names another allele first;
- on the IGH set at 30 % SHM, 79 rows change a call, a segment boundary or
  the junction. J accuracy rises by 0.10 points; the rest moves within
  noise.

The changes matter for:
- reads that start after their V, or stop before their J;
- reads that hold two chains, such as scFv reads;
- heavily mutated reads under the hosts' default identity floor (0.6);
- `--calls fixed`;
- AIRR tables with an index column, and reference files with a byte-order
  mark or non-UTF-8 bytes;
- IMGT-gapped and `swig prepare`d V databases;
- reads without a name in Swig Web.

## What changes for users

- **D on reads that start after their V.** The D search window was
  bounded by the best V placement, even when that placement was a chance
  hit downstream of the real D. On a read that starts in CDR3, the window
  then began after the D, and the read lost a D it clearly holds, reporting
  a D presence of 0. Only a V that anchors the read (E ≤ 10⁻³) now bounds
  the window; without one, the window starts at the read's start. On 1,913
  simulated IGH reads that start after their V:

  | | 0.39.2 | 0.39.3 |
  |---|---|---|
  | D allele right | 81.86 % | 93.10 % (215 gained, 0 lost) |
  | D presence Brier score | 0.109 | 0.003 |

  The wider window also lets more chance Ds in on reads that have none. On
  600 D-free IGH reads that start in CDR3 behind 30 random bases, a D is
  called on 79.8 % of reads instead of 68.7 %; behind 100 random bases, on
  97.7 % instead of 88.7 %. Real IGH reads almost always hold a D, so the
  trade favours them. The chance-D rate is still open; see "Still open".
- **D on reads without a J or without a V.** The D window of a read that
  stops before its J ran to the read's end, and that of a read that starts
  after its V ran from the read's start. It now ends 60 bases after the V,
  or starts 60 bases before the J: the reach of a D and its N region. This
  bounds where a chance D can lie, not how often one is called. On 600
  D-free reads cut 100 random bases after their V, a D is called on 93.5 %
  of reads instead of 97.2 %; with 60 bases or fewer, nothing changes.
- **Reads that hold two chains.** On a read with two loci, such as an scFv,
  the J alleles of all loci were searched in one request and kept as one
  top 3, so one chain's alleles could crowd out the other chain's J. Each
  locus is now searched and kept on its own. The locus of such a read was
  also chosen differently with and without posteriors: by the likelihood of
  each locus' placements, each under its own estimate of the read's
  divergence, or by the best fixed-scale path. Each locus is now
  adjudicated at the divergence of the whole read, and the call takes the
  best path in the model that `--calls` names. On 800 simulated scFv reads:

  | | 0.39.2 | 0.39.3 |
  |---|---|---|
  | light-chain calls without a J | 20 | 0 |
  | calls with a junction | 774 | 794 |
  | locus differs with `--no-posterior` (default calls) | 17 | 0 |
  | locus differs with `--no-posterior` (`--calls fixed`) | 17 | 0 |

  `--calls fixed` now takes effect on these reads: the default and fixed
  calls differ in locus, V or J on 32 of them instead of 13. Under `--calls fixed` with
  posteriors, 20 calls hold a V posterior below 0.5, because the
  fixed-scale call need not be the posterior's mode. Default calls hold
  none.
- **Identity floor.** `--min-identity` (default 0.6 in Swig Web and
  `swig-cli`; off in the `swig` binary) was measured on the alignment
  including its stitched outer flank. A heavily mutated J whose local
  alignment was 20 of 25 bases, followed by a divergent flank, was dropped.
  The floor now judges a hit's local alignment. At 0.6, against 0.39.2, the
  IGH set at 30 % SHM gains:
  - J +0.30 points (9 reads gained, 0 lost; McNemar p 0.004);
  - junction +0.31 points (9 gained, 2 lost);
  - C +0.21 points (9 gained, 4 lost).

  With 0.39.3, output at 0.6 is byte-identical to output without a floor
  on all seven suite sets.
- **Best local alignments kept.** V and J alignments were kept by their
  score with the outer flank, and identical allele copies each took a
  place in the top 3. A hit that ends at the read's end pays no flank. So
  three copies of a short chance J at the read's end could push out the
  real J, whose divergent flank lowered its score. Each stage now keeps
  the top 3 by extended score and the top 3 by local score. At 30 % SHM, 3
  reads gain their right J. On the other sets, 21 reads called as a
  5′-partial V allele now also list the full-length alleles tied with it on
  the called bases; their first name is unchanged.
- **`--calls fixed` rows.** A placement tied with the call but spanning
  fewer read bases could become the row's primary. The AIRR row then
  dropped read bases from `sequence_alignment` and named another gene than
  the call. The primary is now a tied placement that spans the called
  bases.

  | Rows that do not stitch (`--calls fixed`) | 0.39.2 | 0.39.3 |
  |---|---|---|
  | IGH mid (5,000 reads) | 7 | 0 |
  | IGH hard (5,000 reads) | 19 | 0 |
  | IGH 30 % SHM (3,000 reads) | 40 | 0 |

  In default output one suite row changes for the same reason: IGH mid
  read sim3553 now starts V where its path does and lists IGHV4-30-4\*05
  before \*06. The truth is \*06, so V allele accuracy on that set falls by
  0.02 points.
- **A V pays a presence cost.** A D and a J pay one; a V did not, so a read
  that starts after its V got a V from a few chance bases. A V now pays
  `P*_V`, derived as `P*_D` and `P*_J` are (12.69 for the 125 human IGHV
  genes of the suite). On 600 IGH reads that start at their J behind 30
  random bases, a V is called on 1 read instead of 6; behind 100 bases, on
  4 instead of 17. In the fix's own test, on 27,000 reads that keep only
  the last 8–90 bases of their V, it removed 129 V calls, all wrong, and
  lost no right call.
- **AIRR tables with an index column.** A table whose first column name
  is empty, as pandas writes its index, was read one column off: each read
  took the column left of `sequence` as its sequence and got no calls.
  Such tables are now read correctly on every route.
- **Reference files with a byte-order mark or non-UTF-8 bytes.** The
  `swig` binary and `swig-cli` refused a germline FASTA with a UTF-8
  byte-order mark, or with a non-UTF-8 byte such as a Latin-1 letter in a
  description. Swig Web read them. Now the byte-order mark is skipped and
  invalid bytes are replaced, as for reads. The binary's `evaluate` reads
  truth and calls files the same way.
- **IMGT-gapped V alleles.** When no gapped V allele held one of the
  checked codons, as when the only gapped alleles are partial, the check of
  the IMGT gaps failed and IMGT numbering was dropped for every allele. A
  codon that no gapped allele holds is now skipped; at least one must be
  held.
- **`swig prepare` round trip.** A prepared database lost two facts on
  reload: that a 5′-partial allele starts after FR1 or CDR1, and that an
  allele is 3′-partial. A reloaded 5′-partial allele gained a wrong FR1,
  and a second `prepare` changed the database. They are now written to the
  header as `before=` and `partial=3`, and a database prepares to itself.
- **Allele names without `*`.** A gene is the allele name's part before
  `*`. Without it, every allele is a gene of its own, so the V search keeps
  at most 16 alleles instead of 16 genes. On the suite's IGH set with every
  `*` replaced, V allele accuracy at 30 % SHM falls by 8.2 points. A V
  database where more than 16 names lack `*` now gets a warning: the `swig`
  binary prints it, and the WebAssembly module returns it with its other
  database warnings. No shipped reference triggers it. The accuracy cost
  itself is still open.
- **Reads without a name in Swig Web.** Each Swig Web worker numbered
  unnamed reads from 1, so a run gave several reads the same
  `sequence_N`. Reads are now named by their position in the input before
  they are batched, for FASTA, FASTQ and AIRR tables, as `swig-cli` and the
  `swig` binary already did. A table without a `sequence_id` column gets
  one. `swig-cli` output is unchanged.
- **`swig-cli` help.** `--minimum-identity` is described as the local
  alignment's identity.

## For developers

- `npm run test:engine` runs the engine tests with debug assertions and
  overflow checks on (`--profile checked`, optimized), then builds the
  release binary. Debug invariants did not run in CI before, and one
  property test broke its stage's precondition; it is fixed.
- ALGORITHM.md now states the SIMD bound for C alignments correctly: a C
  record stays in the SIMD lanes up to about 15 kb, not 2 kb, and one
  longer record sends every C job to the scalar kernel.

## Output

The native and WebAssembly routes still write byte-identical AIRR (release
smoke tests). The tables compare 0.39.3 with 0.39.2 on the same reads, both
scored by the 0.39.1 evaluator, in percentage points.

Default settings of the `swig` binary (no identity floor):

| Set | Reads | V allele | D allele | J allele | C allele | Junction exact |
|---|---|---|---|---|---|---|
| IGH easy | 3,000 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 |
| IGH mid | 5,000 | −0.02 | +0.00 | +0.00 | +0.00 | +0.00 |
| IGH hard | 5,000 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 |
| IGH 30 % SHM | 3,000 | +0.00 | +0.03 | +0.10 | −0.04 | +0.00 |
| BCR (IGH+IGK+IGL) | 7,500 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 |
| TCR (TRA+TRB+TRD+TRG database) | 5,000 | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 |

At 30 % SHM, V gains 2 reads and loses 2, D gains 7 and loses 6, J gains 3
and loses none, and C gains 3 and loses 4. None of these differences is
significant.

With `--min-identity 0.6`, the default of Swig Web and `swig-cli`, only the
IGH 30 % SHM row differs from the table above: V +0.00, D +0.00, J +0.30,
C +0.21, junction +0.31. One IGH hard row changes its J, without changing
accuracy.

## Timing

Speed is unchanged. These runs use the `swig` native binary on the IGH mid
and hard sets: 5,000 reads each, with human IGH V, D, J and C. Each figure
is the median of 5 interleaved runs on a 4-core cloud VM. Wall time
includes database preparation.

| Set | Threads | 0.39.2 | 0.39.3 |
|---|---|---|---|
| IGH mid | 1 | 3.26 s | 3.17 s |
| IGH hard | 1 | 3.10 s | 3.05 s |
| IGH mid | 4 | 0.98 s | 0.97 s |
| IGH hard | 4 | 0.89 s | 0.90 s |

Single runs varied by up to 12 % from the median. This machine was slower
than the one that timed 0.39.2, so compare the columns here, not with the
0.39.2 notes.

## Still open

`engine/docs/REVIEW.md` §0 lists what this release does not fix. From this
review:
- the chance-D rate on reads without a D: `P*_D` is priced for a
  junction-sized window, so most D-free reads with 30 or more random bases
  next to their D window still get a D. It needs an N/P length prior;
- the accuracy cost of V databases whose allele names lack `*`;
- 7 of 1,200 V-less test reads that still get a V, because the read's
  divergence is estimated partly from that V's own columns;
- smaller items: valid reads before a malformed record in the same batch
  are dropped; IGHT, IGHZ and IGHO constant alleles get no locus;
  `--c-prefix-identity` and `--double-d` are ignored without the database
  they need;
- `swig-cli` with `prepareMetadata` rejects a V database written by
  `swig prepare`.
