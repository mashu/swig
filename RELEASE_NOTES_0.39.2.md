# Swig 0.39.2

Swig 0.39.2 fixes the defects found by an audit of the Swig aligner made
after 0.39.1. The aligner runs in Swig Web, in `swig-cli` (native and
WebAssembly) and in the `swig` binary, so every change below applies to all
three. On ordinary reads, the V, D, J and C calls do not change. The
changes hit:

- reads of alleles whose reference anchors are misplaced;
- reads that hold little of their V or J;
- reads run against a database of the wrong locus;
- malformed input.

The aligner is also about 13 % faster single-threaded. Every other Swig
tool is unchanged.

## What changes for users

- **Junctions of alleles with misplaced anchors.** Swig's shipped references
  misplace the `SWIGMETA` anchors of some alleles:
  - 50 J alleles have an anchor on an aromatic codon before their real
    Trp/Phe118. Human TRAJ24\*01–03 are 15 bases early and TRAJ48\*01 is 21
    bases early.
  - 13 5′-partial V alleles have the wrong frame, among them IGLV2-11\*03,
    IGLV2-8\*03, IGHV4-30-4\*05, TRBV5-8\*02 and TRBV11-3\*03.

  Every read of those alleles got a junction several codons short or a wrong
  `productive` flag. The aligner now checks a given anchor against the
  allele's own motifs:
  - A J anchor is contradicted when it lies off `[WF]G.G` while a confident
    motif lies elsewhere in its frame.
  - A V anchor is contradicted when its motifs are weak and it lies in a
    different frame from all of its given region starts.
  - A contradicted `SWIGMETA` anchor gives way to the motif hit, and the row's
    `v_annotation_source` or `j_annotation_source` names the annotation's
    source.
  - An explicit `anchor=` key is kept, but flagged and counted in the
    uncertain-anchor warning.

  The references themselves are unchanged in this release.
- **No calls without evidence.**
  - A read is called only when its path keeps a V or J part that is
    significant on its own (E ≤ 10⁻³). Before, a significant local alignment
    was enough even when the path kept one base of it. Light-chain reads run
    against a heavy-chain database were called with junctions and posteriors
    near 1:
    - 1,000 IGK reads against the IGH database got 72 V calls and
      59 junctions, 12 of them productive.
    - They now get no call.
    - Reverse reads under a forced strand were called the same way.
  - A J now pays a presence cost, as a D does. Before, any J alignment of
    positive score completed a path for free, so a read that stops before its
    J got a J from a few chance bases, often with a junction. On 6,000 IGH
    reads cut before their J, that happened to 7.9–9.6 % of reads; it now
    happens to 0.9–1.8 %.
  - A junction is no longer reported when the projected Cys104 and Trp/Phe118
    codons overlap. Such a junction is 4–5 bases long, and its FR3 runs into
    FR4. The benchmark sets had 358 of them, none in the truth.
- **Reads that start inside CDR3 or J.** A chance V hit no longer bounds the
  J search or chooses the J's locus. On 1,197 IGH reads cut to start 6 bases
  before their J, against the composed IGH+IGK+IGL database, the right J gene
  is called on 1,179 reads instead of 641.
- **V search.**
  - Every allele of a kept V gene is a candidate. A cap of 64 alleles made
    calls depend on the order of the V FASTA file and dropped the best allele
    of reads that cover only the V 3′ end.
  - An allele with an internal indel no longer silences its gene's own seeds.
    Exact 3′ fragments of IGHV4-31 were called IGHV4-34 at posterior 0.998.
  - V indels of 9 bases or more are aligned across. Before, the alignment
    stopped at the indel and the junction moved by the indel's length.
  - A 5′ adapter, UMI or random flank is no longer stitched onto V.
- **`*_support` columns.** V, D and J support is now `K·m·n·e^(−λS)` over the
  read length `m` and the database length `n`, the formula the significance
  gates use. Before, support applied BLAST's length adjustment. For D sets
  that adjustment collapsed to a constant: a TRB+TRD D call scoring 12
  reported 1.01 instead of 6.62. Values change on every row, by about 1.1×
  for V, 1.3–1.5× for J and 2.2× for D (7× on TCR). `c_support` is
  unchanged. No call changes because of this.
- **Posteriors.** They change in their last digits: the median change is
  10⁻⁵, and 95 % of changed values move by less than 0.01. This comes from
  the J presence cost and the other model fixes. `d_presence_posterior` now
  prints `0.000000` instead of `-0.000000` for reads without a D hypothesis,
  such as light chains.
- **Other columns.**
  - `sequence_alignment_aa` and `germline_alignment_aa` are each translated
    in their own frame. Before, an indel that is not a multiple of 3 put one
    of them out of register.
  - IUPAC codes and other symbols of the read now appear as such in every
    derived column (alignments, regions, junction, N/P). Before, they became
    N there while `sequence` kept them.
  - The `*_alternatives` columns list an allele merged into another's
    alignment at its own germline coordinates.
- **Reads.**
  - An AIRR table's `sequence` now keeps every symbol, as FASTQ does: a
    no-call `.` reads as N instead of being dropped, which used to shift
    every later coordinate.
  - In FASTQ and AIRR input, sequences and qualities must be printable ASCII
    (`!` to `~`). A quality must match its sequence in length. A quality
    holding a tab or a non-ASCII byte is now an error naming its line;
    before, a tab was rewritten to a space, and a non-ASCII byte in an AIRR
    sequence silently became `NNN`.
  - A file with classic Mac (CR-only) line endings is now an error naming
    the line and suggesting a conversion. Before, a FASTA file read as one
    empty read with exit status 0. LF and CRLF files are unaffected.
  - A UTF-8 byte-order mark at the start is skipped. Header bytes that are
    not UTF-8 read as U+FFFD instead of stopping the run. The command line
    now reads such files as Swig Web does.
- **Germline databases.**
  - These records are now refused with a message:
    - an allele name or alias holding `,`, `;` or `|` (or `+` in a D name),
      since these separate names in the report;
    - a record without a name;
    - a sequence line holding a carriage return.
  - The locus comes only from gene names such as `IGHV…` or `TRBJ…`. Before,
    header words counted too: "strain" made a database TRA, and "light"
    made it IGH.
  - An allele whose name has no locus code is annotated with its nearest
    relative's locus.
  - A database written by `swig prepare` reads back with the loci, trust and
    provenance its annotation gave.
- **Removed.**
  - The `swig` binary's `--residual` option, which used the
    commit-V/J-then-D adjudicator, is gone.
  - The library's tandem-D grammar is gone. Swig's double-D screen is
    unchanged.
  - The WebAssembly export `swig_reset_options` is gone; no host called it.
- **`swig simulate` and `swig evaluate`** (the engine binary).
  - `simulate`:
    - keeps each read within one locus;
    - draws only alleles a rearranged transcript carries;
    - checks `--rc`, `--shm` and `--dd` against their domains;
    - warns about `--dd` only when it is given.
  - `evaluate`:
    - reports a failed write;
    - treats missing inputs as usage errors;
    - no longer lets allele identity depend on the reported orientation;
    - scores a single D call on a tandem-DD read against the D it covers.

## Output

The native and WebAssembly routes still write byte-identical AIRR (release
smoke tests).

Against 0.39.1, the comparison is paired on the same reads: 28,500 simulated
reads plus 1,000 IGK reads run against the IGH database. Both versions are
scored by the 0.39.1 evaluator.

| Set | Reads | V allele | D allele | J allele | C allele | Junction exact |
|---|---|---|---|---|---|---|
| IGH easy | 3,000 | +0.00 | +0.00 | +0.00 | +0.00 | −0.25 |
| IGH mid | 5,000 | +0.02 | +0.00 | +0.00 | +0.00 | −0.13 |
| IGH hard | 5,000 | +0.06 | +0.02 | +0.00 | +0.00 | −0.18 |
| IGH 30 % SHM | 3,000 | **+1.80** | −0.03 | −0.43 | −0.46 | +0.04 |
| BCR (IGH+IGK+IGL) | 7,500 | +0.00 | +0.00 | +0.00 | +0.00 | −0.59 |
| TCR (TRA+TRB+TRD+TRG) | 5,000 | +0.00 | +0.00 | +0.00 | +0.00 | −2.69 |

Changes are in percentage points.

- **V at 30 % SHM.** The +1.80 points are 69 reads gained and 15 lost
  (McNemar p 2·10⁻⁹), mostly from the removed 64-allele cap.
- **J and C at 30 % SHM.** Two kinds of call are no longer made: calls
  whose path holds no significant segment, and J calls too weak to pay the
  J presence cost. Most of the affected reads had the wrong V or the wrong
  J gene.
- **Junction losses on the in-locus sets.** Every exact junction lost there
  belongs to a read of a corrected-anchor allele. The frozen truths were
  simulated from the misplaced anchors, so they still hold the short
  junctions. No exact junction was lost on any other read.

## Timing

The `swig` native binary, single-threaded, on the IGH mid and hard sets
(5,000 reads each, human IGH V/D/J/C): median of 3 interleaved runs, 4-core
Intel Xeon at 2.10 GHz.

| Set | 0.39.1 | 0.39.2 |
|---|---|---|
| IGH mid | 3.21 s | 2.78 s |
| IGH hard | 3.10 s | 2.76 s |

The D stage now aligns each candidate once instead of scoring all of them
and then aligning most of them again. Other changes save work too:
- the best cell is tracked only on traceback passes;
- anchored C profiles stay on the uniform SIMD path;
- the adjudication recursion no longer rescans its slots;
- rows are written in place.

## Still open

`engine/docs/REVIEW.md` §0 lists what this release does not fix:
- the strand test's null model;
- posterior overconfidence at high SHM;
- frame checks when annotation is transferred between siblings;
- public fields that bypass validation;
- serial reading and writing at 4 threads;
- regenerating Swig's references with correct anchors.
