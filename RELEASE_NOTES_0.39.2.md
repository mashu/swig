# Swig 0.39.2

Swig 0.39.2 fixes most of the defects that an audit of the Swig aligner,
made after 0.39.1, found. The aligner runs in Swig Web, in `swig-cli`
(native and WebAssembly) and in the `swig` binary. Every change to calls and
AIRR columns below applies to all three. Changes to reading files, options
and exports name the route they apply to.

On ordinary reads, the J and C calls do not change and the V and D calls
barely do. On 25,500 simulated IGH, BCR and TCR reads at ordinary SHM:
- 44 V calls change, mostly the order of tied IGHV4 alleles; the first V
  gene changes on 4 reads;
- 3 D calls change.

The changes hit:
- reads of alleles whose reference anchors are misplaced;
- reads that hold little of their V or J;
- reads run against a database of the wrong locus;
- heavily mutated reads;
- reads with a long V indel or a 5′ adapter;
- malformed input.

Single-threaded, the aligner takes 15–16 % less time. Every other Swig tool
is unchanged.

## What changes for users

- **Junctions of alleles with misplaced anchors.** Swig's shipped references
  misplace the `SWIGMETA` anchors of some alleles:
  - 50 J allele records, across species, have an anchor on an aromatic
    codon before their real Trp/Phe118. Human TRAJ24\*01–03 are 15 bases
    early and TRAJ48\*01 is 21 bases early.
  - Some 5′-partial V alleles have the wrong frame, among them
    IGHV4-30-4\*05, IGLV2-11\*03, IGLV2-8\*03, TRBV5-8\*02 and TRBV11-3\*03.

  Reads of those J alleles got a junction 3–9 codons short. Reads of those V
  alleles got a junction 2–44 bases off and out of frame, so never
  productive.

  The aligner now checks a given anchor against the allele's own motifs:
  - A J anchor is contradicted when it lies off `[WF]G.G` while a confident
    motif lies elsewhere in its frame.
  - A V anchor is contradicted when its motifs are weak and it lies in
    another frame than its given CDR1–FR3 starts, which share one frame.
  - A contradicted `SWIGMETA` anchor gives way to a confident motif hit, and
    the row's `v_annotation_source` or `j_annotation_source` names the
    annotation's source.
  - Without a confident motif hit, the anchor is kept and flagged, and it
    counts in the uncertain-anchor warning. So is an explicit `anchor=`.

  This corrects all 50 J allele records and 8 V alleles. Five V alleles
  still get the wrong junction:
  - IGHV1/OR15-5\*01 and rhesus IGHV1-197\*02 have no confident motif to
    take the anchor's place, so the anchor is kept and flagged.
  - Mouse TRAV7-3\*02 and TRAV9D-4\*03 and sheep TRAV45-3\*02 are not
    caught.

  The references themselves are unchanged in this release.
- **No calls without evidence.**
  - A read is called only when its path keeps a V or J part that is
    significant on its own (E ≤ 10⁻³). Before, a significant local alignment
    was enough even when the path kept one base of it.
    - Light-chain reads run against a heavy-chain database were called with
      junctions and posteriors near 1. 1,000 IGK reads against the IGH
      database got 72 V calls and 59 junctions, 12 of them productive. They
      now get no call.
    - Reverse reads under a forced strand were also called, often with
      posteriors near 1: 67 of 5,000 reversed IGH reads under
      `--strand forward`. They now get no call.
  - A J now pays a presence cost, as a D does. Before, any J alignment of
    positive score completed a path for free, so a read that stops before
    its J got a J from a few chance bases, often with a junction. The test
    used 6,000 IGH reads cut where their J starts: 1,000 mid and 1,000 hard
    reads, each with 0, 30 or 100 random bases after the cut. A spurious J
    went from 7.9 % to 1.0 % of the mid reads and from 10.4 % to 2.4 % of
    the hard reads.
  - A projected Cys104 now counts only when it starts before the J's
    alignment, and a projected Trp/Phe118 only when it ends after the V's.
    Neither counts when their codons overlap. Before, when V or J was
    trimmed into its anchor codon, the other segment's anchor could land
    inside it. That gave a junction of 4–5 bases, FR3 running into FR4, or a
    junction shifted into the D, often with `productive` T. The benchmark
    sets lose 357 such junctions, none of them in the truth.
- **Reads that start inside CDR3 or J.** A chance V hit no longer bounds the
  J search or chooses the J's locus. The test cut 1,197 IGH reads to start
  6 bases before their J and ran them against the composed IGH+IGK+IGL
  database. The right J gene is now called on 1,179 reads instead of 641.
- **V search.**
  - Every allele of a kept V gene is a candidate. Before, a cap of 64
    alleles made calls depend on the order of the V FASTA file. On reads
    that cover only the V 3′ end, it dropped the best allele in up to two
    thirds of reads.
  - An allele with an indel near one end no longer silences the seeds its
    gene holds past that indel. Exact 3′ fragments of IGHV4-31 were called
    IGHV4-34 at posterior 0.998.
  - When an allele's seeds tie between two offsets, the offset nearest 0
    wins. Before, the tie could resolve differently from one run to the
    next.
  - V indels of 9 to about 60 bases are aligned across. Before, V was
    aligned on one side of the indel only:
    - when the 5′ side was kept, the rest of V became N region and the
      junction moved by the indel's length;
    - otherwise the V start moved, or a paralog was called.
  - V and J alignments are kept by their score with the outer flank, the
    score adjudication uses. Before, they were kept by local score and
    extended afterwards. A 5′-partial allele could then be dropped behind
    its full-length siblings even when its extended score was higher.
- **Outer flanks.** The bases stitched before a V and after a J now stop
  where the running score falls more than 30 below the alignment's edge.
  Before, every co-linear base up to the germline's or the read's end was
  taken.
  - This keeps most of an adapter, UMI or random flank off a 5′-truncated V.
    An 8-base flank is still taken whole, and about 15 bases of a longer
    one.
  - It also cuts some real, divergent V starts. At 30 % SHM, exact V starts
    fall from 98.26 % to 97.04 % and the mean V start error rises from 0.42
    to 1.34 bases.
  - Other sets are unchanged or nearly so.
- **D candidates.** The D stage keeps only admissible D alignments among its
  32: those with an exact run of `--min-d-match` bases or distributed-D
  evidence. Before, inadmissible hits took slots first and were dropped
  afterwards, so admissible Ds behind them could be lost.
- **`*_support` columns.**
  - V, D and J support is now `K·m·n·e^(−λS)` over the read length `m` and
    the database length `n`, the formula the significance gates use.
  - Before, support applied BLAST's length adjustment. That takes 10–13
    bases off every allele, nearly the whole of a D. For example, it shrank
    the TRB+TRD D set from 74 to about 12 bases.
  - Values change on every row: by about 1.1× for V, 1.3–1.5× for J and
    2.2× for D (7× on the TCR set). `c_support` is unchanged.
  - No call changes because of this.
- **Posteriors.** Most change only in their last digits. The median change
  is 3·10⁻⁵, and 93 % of changed values move by less than 0.01 (86 % at
  30 % SHM). This comes from the J presence cost and the other model fixes.
  - `c_posterior` is unchanged where C is still called.
  - `d_presence_posterior` now prints `0.000000` instead of `-0.000000` for
    reads without a D hypothesis, such as light chains.
- **Other columns.**
  - `sequence_alignment_aa` and `germline_alignment_aa` are each translated
    in their own frame. Before, an indel that is not a multiple of 3 put one
    of them out of register.
  - IUPAC ambiguity codes of the read (R, Y, K, …) now appear as such in
    every derived column: alignments, regions, junction, N/P and the
    double-D sidecar. Before, they became N there while `sequence` kept
    them. Any other non-base symbol reads as N everywhere, `sequence`
    included.
  - The `*_alternatives` columns list an allele merged into another's
    alignment at its own germline coordinates.
  - `v_annotation_source` and `j_annotation_source` report a `SWIGMETA`
    provenance only for coordinates that `SWIGMETA` supplied. On the
    benchmark sets about 500 rows change for this reason besides the
    corrected anchors. For example, TRBJ2-2P\*01 and TRAJ16\*01–02 now report
    `low-confidence` instead of `unclassified`.
- **Reads.**
  - An AIRR table's `sequence` keeps every symbol, as FASTQ does: a no-call
    `.` reads as N instead of being dropped, which used to shift every later
    coordinate. An AIRR `quality`, when given, must match its sequence in
    length, as a FASTQ quality already had to.
  - FASTQ and AIRR sequences and qualities must be printable ASCII (`!` to
    `~`).
    - A quality holding a tab or a non-ASCII byte is now an error naming
      its line. Before, a tab became a space.
    - Before, a non-ASCII character in an AIRR sequence was silently
      dropped.
  - In the `swig` binary:
    - A file with classic Mac (CR-only) line endings is now an error naming
      the line and suggesting a conversion. Before, a FASTA file read as one
      empty read with exit status 0. LF and CRLF files are unaffected.
    - A UTF-8 byte-order mark at the start is skipped.
    - Bytes that are not UTF-8 read as U+FFFD instead of stopping the run
      without naming the line.
    - So the binary now reads such files as Swig Web does. `swig-cli`
      reads FASTA and FASTQ itself, and that reader is unchanged.
- **Germline databases.**
  - These records are now refused with a message:
    - an allele name or alias holding `,`, `;` or `|` (or `+` in a D name),
      since these separate names in the report;
    - a record without a name;
    - a line holding a carriage return. A germline file with CR-only line
      endings now fails naming the line, with a conversion hint. Before,
      its message repeated the whole file.
  - The locus comes from a new `locus=` header key (for example `locus=IGH`
    or `locus=TRA,TRD`). Otherwise it comes from a gene name such as `IGHV…`
    or `TRBJ…` in the allele name or the header. An unknown code in
    `locus=` is an error naming the allele.
  - Before, when a name gave no locus, any header word holding a locus code
    counted: "strain" made an allele TRA, and "light" made it IGH.
  - An allele whose name and header give no locus is annotated in the
    context of its nearest relative of known locus, when it has one. It
    still pairs with every locus.
  - `swig prepare` writes `locus=` for every allele of known locus, an
    `annotation=` that names the annotation's source, and `uncertain=` for
    the starts the annotation doubts. A prepared database reads back with
    the loci, trust and provenance its annotation gave.
- **Options.**
  - `--d-mismatch` and `--j-mismatch` take −1000 to −1; values below −1000
    are refused with a message. This applies to the `swig` binary and the
    WebAssembly module. `swig-cli --vdj` checks `-D_penalty` and
    `-J_penalty` against the same range. Before, any negative value was
    accepted. In 0.39.1, `--j-mismatch -2147483648` aborted the run. With
    `--j-mismatch -1000000`, two reads took 4 s and one reported a J score of
    −2,999,908.
  - `--d-allele max-marginal` names the decoded D allele first, with its own
    germline coordinates. A decoding is kept only when its reported D
    posterior is at least the call's.
  - `--calls fixed` reads each presence prior at its own scoring's λ. With
    `--d-mismatch`, `--j-mismatch` or the `igblast` profile, posteriors
    change: at `--d-mismatch -4` the D prior was 6.30 nats instead of 6.62.
    Calls do not change.
- **The WebAssembly module** (Swig Web and the `swig-cli` WebAssembly route).
  - Batches still work once the module's memory passes 2 GiB. Pointers and
    lengths are now read as unsigned.
  - `swig_alloc` returns 0 when it cannot allocate instead of trapping, so
    the host reports that the aligner ran out of WebAssembly memory.
  - An aligner panic is reported with its message instead of "unreachable".
    Later calls on that instance refuse.
  - `swig_set_option` refuses a flag given a value, such as `no-posterior`
    with `false`, and an option given none. Before, the value was dropped
    and the flag switched on.
  - `swig_version` reports the kernel the run uses.
  - `swig_reset_options` is gone; no host called it.
- **Removed.**
  - The `swig` binary's `--residual` option, the commit-V/J-then-D
    adjudicator, is gone.
  - The library's tandem-D grammar is gone. Swig's double-D screen, which
    reports a second D beside the call, stays.
- **The `swig` binary's other modes.**
  - A command line that lacks an input its mode needs now exits with status
    2 and prints the usage, as other bad command lines do. Before, it exited
    with status 1. This covers:
    - assignment without `--v`, `--j` or the reads;
    - `prepare` or `simulate` without `--out`;
    - `evaluate` without `--truth`;
    - `--double-d` without `--double-d-out`.
  - `simulate`:
    - keeps each read within one locus;
    - draws no V or J allele without a trusted anchor, and no C allele that
      `prepare` warns about (5′-truncated or genomic);
    - checks `--rc`, `--shm` and `--dd` against their domains;
    - warns when `--dd` is given for a database with fewer than two D
      alleles.
  - `evaluate`:
    - reports a failed write;
    - no longer lets allele identity depend on the reported orientation;
    - scores a single D call on a tandem-DD read against the D it covers.

## Output

The native and WebAssembly routes still write byte-identical AIRR (release
smoke tests).

Against 0.39.1, the comparison is paired on the same 28,500 simulated reads,
and both versions are scored by the 0.39.1 evaluator. Changes are in
percentage points.

| Set | Reads | V allele | D allele | J allele | C allele | Junction exact |
|---|---|---|---|---|---|---|
| IGH easy | 3,000 | +0.00 | +0.00 | +0.00 | +0.00 | −0.25 |
| IGH mid | 5,000 | +0.02 | +0.00 | +0.00 | +0.00 | −0.13 |
| IGH hard | 5,000 | +0.06 | +0.02 | +0.00 | +0.00 | −0.18 |
| IGH 30 % SHM | 3,000 | **+1.80** | −0.03 | −0.43 | −0.46 | +0.04 |
| BCR (IGH+IGK+IGL) | 7,500 | +0.00 | +0.00 | +0.00 | +0.00 | −0.59 |
| TCR (TRA and TRB reads, TRA+TRB+TRD+TRG database) | 5,000 | +0.00 | +0.00 | +0.00 | +0.00 | −2.69 |

- **V at 30 % SHM.** The +1.80 points are 69 reads gained and 15 lost
  (McNemar p 2·10⁻⁹), mostly from the removed 64-allele cap.
- **J and C at 30 % SHM.** 20 calls are no longer made: their paths hold no
  significant segment. 15 of them had a wrong or missing V gene, but 16 had
  the right J gene. With the J presence cost and the shorter stitched
  flanks, 23 reads lose a right J and 10 gain one, and 25 lose a right C and
  14 gain one.
- **Junction losses.** On the IGH easy, mid and hard, BCR and TCR sets,
  every exact junction lost belongs to a read of a corrected-anchor allele.
  The frozen truths were simulated from the misplaced anchors, so they still
  hold the old junctions. At 30 % SHM, 14 exact junctions are lost and 15
  gained; 4 of the lost are reads of IGHV4-30-4\*05.
- **IGK reads against the IGH database.** The 1,000 off-locus reads are not
  scored here; see "No calls without evidence".

## Timing

These runs use the `swig` native binary on the IGH mid and hard sets: 5,000
reads each, with human IGH V, D, J and C. Each figure is the median of 5
interleaved runs on a 4-core Intel Xeon at 2.10 GHz (cloud VM). Wall time
includes database preparation.

| Set | Threads | 0.39.1 | 0.39.2 |
|---|---|---|---|
| IGH mid | 1 | 3.05 s | 2.57 s |
| IGH hard | 1 | 2.92 s | 2.49 s |
| IGH mid | 4 | 0.87 s | 0.78 s |
| IGH hard | 4 | 0.84 s | 0.83 s |

Single runs varied by ±10 % on this machine.

The D stage now aligns each candidate once instead of scoring all of them
and then aligning most of them again. Other changes save work too:
- the best cell is tracked only on traceback passes;
- anchored C profiles stay on the uniform SIMD path;
- the adjudication recursion no longer rescans its slots;
- rows are written in place.

SIMD lanes now take scores up to 30,000 instead of 8,000. The C anchor bonus
grows with the longest C record, so one record of about 4 kb sent every C
alignment to the scalar kernel. The shipped Salmo salar IGH set has two such
records, IGHD\*01 and IGHDD\*01. Its C alignments run on SIMD again; output is
unchanged.

## Still open

`engine/docs/REVIEW.md` §0 lists what this release does not fix. Among
them:
- the strand test's null model;
- posterior overconfidence at high SHM;
- J and C calls dropped with a read whose V fails the anchor gate;
- frame checks when annotation is transferred between siblings, and the
  five V alleles above that keep a wrong anchor;
- FR/CDR regions taken from another allele than `v_call`'s first name;
- public fields that bypass validation;
- serial reading and writing at 4 threads;
- regenerating Swig's references with correct anchors.
