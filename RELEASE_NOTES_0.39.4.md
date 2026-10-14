# Swig 0.39.4

Swig 0.39.4 restores the plausibility guards and the constant-region (C)
rules of SwiftIG, the aligner of Swig ≤ 0.38, which the switch to the Rust
aligner in 0.39.0 dropped. Each guard is a comparison per hit or per V–J
pair: CPU time is unchanged (within ±2 % on 10,000 simulated reads, one
thread). It also restores the sensitivity of the CHMMAIRRa PCR-chimera
filter, and fixes defects in Swig Web and `swig-cli` that predate the
switch. The aligner runs in Swig Web, in `swig-cli` (native and
WebAssembly) and in the `swig` binary; every change to calls below applies
to all three.

On ordinary reads, calls do not change. On 20,000 simulated human BCR reads
(IGH, IGK, IGL) of the IgM-like preset (mid, mean SHM 1.5 %) and 20,000 of
the IgG-like preset with SHM indels (hard, mean SHM 7.7 %), against 0.39.3:

| | mid 0.39.3 | mid 0.39.4 | hard 0.39.3 | hard 0.39.4 |
|---|---|---|---|---|
| V, D, J allele right | 98.70 / 93.46 / 99.11 % | same | 93.77 / 85.88 / 94.92 % | 93.77 / 85.83 / 94.91 % |
| reads whose V, D or J call or coordinates changed | | 0 | | 9 |
| templated junctions reported exactly | 99.97 % | 99.13 % | 99.63 % | 95.67 % |
| reads without a templated junction that get one | 90.9 % | 24.1 % | 89.2 % | 14.0 % |
| `productive=T` without a templated junction | 1,122 | 215 | 673 | 90 |
| C allele right | 96.71 % | 94.94 % | 96.04 % | 94.15 % |

Of the 584 C calls lost, 580 are on C tracts shorter than 32 nt, where
SwiftIG's floors (12 aligned columns at 90 % identity) now apply: 8–11 nt
of C are no longer called, nor short tracts below 90 % identity. Swig Web's
isotype column never labelled those either.

## What changes for users

- **V–J pairs must be plausible (SwiftIG's `supports_rearranged_pair`).**
  A V and a J form a rearrangement only when the V starts before the J,
  at most 225 read bases lie between them, and their local alignments
  hold at least 80 aligned columns together and either face the junction
  (the V reaches within 45 bases of its germline 3′ end, the J starts
  within 30 of its 5′ end) or hold a long split (V ≥ 80, J ≥ 30 columns).
  A pair failing that still forms one when V's Cys104 and J's Trp/Phe118
  codons are both aligned, in that order. When the best call breaks a
  guard, the read gets the best call that keeps them: its V alone, its J
  alone, or a V–J pair that passes, as in SwiftIG. In 0.39.3 the first 70
  bases of a V paired with a J 20 random bases later, and V and J 500
  random bases apart formed a rearrangement with a 543-nt junction.
- **A short V does not vouch for a chance J.** Unless both conserved
  codons are aligned, a V of fewer than 80 aligned columns pairs only with
  a J that is significant on its own (E ≤ 10⁻³). The aligner searches J
  only in the locus of the read's V, so the 30–60 bases of an IGHV that a
  chimera carries before an IGK J found a chance IGHJ there; SwiftIG,
  which searched every locus, did not pair them. On 168 such reads, 30
  were called productive in 0.39.3 and 10 now; SwiftIG called none. The
  other 10 pair with an IGHJ that aligns to the IGK J well enough to be
  significant on its own (see below).
- **Hit floors.** A V hit needs 24 aligned columns and a J hit 10 (local
  alignment, before the stitched outer flank), and a V or J clipped by the
  call keeps a positive score and the identity floor (`--min-identity`,
  0.6 in Swig Web and `swig-cli`). In 0.39.3 eight bases of J after 12
  random ones were called J on 18 of 20 reads.
- **A D only between a V and a J.** As in SwiftIG, which paired V and J
  before it looked for a D, a call with a D also holds a V and a J that
  pair. On 1,000 mid IGH reads cut 30 or 100 random bases after their V,
  0.39.3 called a D on 70.0 and 92.5 % of reads, and 0.39.4 on 0.3 and
  0.5 % (each with a chance J). On reads that start at their J behind 30
  or 100 random bases, and so hold no D, 0.39.3 called one on 81.0 and
  95.9 %, and 0.39.4 on none.
- **Reads that start after their V lose their D again.** 0.39.3 had
  restored the D of such reads (right D on 878 of 952 mid IGH reads cut
  where their V ends). Since a D on a read without a V could not be told
  from chance (81–96 % of the D-free reads above got one), 0.39.4 returns
  to SwiftIG's rule: those reads get no D (0 of 952). A read that keeps
  the last 30, 40 or 60 bases of its V keeps its D (805, 866 and 880 right
  calls; 0.39.3: 881, 881, 880).
- **Junctions only from aligned anchors (SwiftIG's `annotate_junction`).**
  The junction, CDR3, `vj_in_frame`, `productive` and `complete_vdj` are
  reported only when the called V aligns the first base of Cys104 and the
  called J that of Trp/Phe118, at most 231 bases apart. 0.39.3 projected an
  anchor along the alignment's end diagonal however far, so a V trimmed
  past Cys104, or a J aligned only after Trp/Phe118, got a junction made of
  bases neither segment explains (see the table above). The cost is the
  templated junctions whose anchor codon is mutated, so that the alignment
  stops before it: 0.84 points on the mid preset, 3.95 points on the hard
  preset.
- **FWR3 without Cys104.** When Cys104 is not placed, FWR3 ends where the
  V alignment does, as in SwiftIG. 0.39.3 dropped FWR3, and Swig Web's
  alignment viewer then reported no region map, on 40 of 400 simulated TRA
  and TRD reads.
- **Constant region as SwiftIG anchored it.** C's base 0 may lie within 7
  bases of the J's end (SwiftIG took C's first aligned base up to 8 bases
  after J and 3 into C; 0.39.3 allowed 3), and C is called on at least 12
  aligned columns at 90 % identity, SwiftIG's floors, besides the anchored
  E-value gate (now `N_C · 15 · e^(−λS) ≤ 10⁻³`, 15 start diagonals). On
  2,552 constructed reads (every human IG and TR C gene, 15–300 nt of C,
  J 3′ ends cut by 1–6 bases, 1–5 bases between J and C, C's first bases
  deleted, 2–20 % C divergence, both strands):

  | | SwiftIG | 0.39.3 | 0.39.4 |
  |---|---|---|---|
  | C right | 2,035 | 2,083 | 2,074 |
  | wrong C gene | 4 | 41 | 1 |
  | right in SwiftIG, not here | — | 163 | 33 |

  0.39.3 lost the C of reads whose J end lay more than 3 bases from C's
  start, called IGLC3 for IGLC2 and TRBC2 for TRBC1 there, could not name
  IGLC1\*01 or IGLC3\*01 (6 bases shorter than their genes), and called
  IGHG2\*16 at 79 % identity on an IGKC sequence. Of the 33 reads left, 29
  have no J, and SwiftIG called their C only after a chance J; 4 have C at
  12–20 % divergence, where the anchored alignment counts leading
  mismatches that SwiftIG's local one trimmed. `--c-min-identity F`
  (default 0.9; 0 turns the floor off, the 12-column floor remains) sets
  the identity floor.
- **`--no-pair-guards`.** Turns the pair guards off (`swig`, and the
  WebAssembly module's `no-pair-guards` option); calls then are those of
  0.39.3 except for C and FWR3.

## PCR-chimera filter (CHMMAIRRa)

The Rust aligner ends a V alignment before trailing mismatches it gives to
the junction, and starts a J after leading ones; a late chimera breakpoint
puts the second parent's informative bases exactly there. CHMMAIRRa now
extends the V alignment gaplessly toward its germline 3′ end, and the J
toward its 5′ start, by the best-scoring run of the read bases beside it
(+1 match, −1 mismatch), before it measures the distance from reference and
threads the read onto the MSA. A J extension must score at least 2: on a J
of about 50 bases, one junction base that matched by chance flagged
ordinary reads; one more base of a V did not. On 2,500 simulated IGH V
chimeras (two V alleles of different genes joined at a homologous
position), 260 J chimeras and 1,400 ordinary reads:

| | SwiftIG | 0.39.3 | 0.39.4 |
|---|---|---|---|
| V chimeras flagged, set 1 | 420 / 500 | 416 / 500 | 420 / 500 |
| V chimeras flagged, set 2 | 1,696 / 2,000 | 1,676 / 2,000 | 1,698 / 2,000 |
| J chimeras flagged (J mode) | 38 / 60, 138 / 200 | 35 / 60, 125 / 200 | 35 / 60, 130 / 200 |
| ordinary reads flagged | 0 | 0 | 0 |

On 11,000 further ordinary simulated reads (4,000 BCR at IgM-like and
4,000 at IgG-like SHM, 3,000 TCR), the extension flags no read that the
filter did not flag without it (V mode 4, J mode 6; SwiftIG's alignments:
9 and 11).

`swig-cli` also threaded every read whose allele the pipeline's allele
refinement had reassigned against the new allele's MSA row, with the old
allele's alignment: threading failed on all of them, and with
`retainUnevaluated` false they were dropped. It now threads against the
allele the aligner aligned (`swig_original_v_call`).

## Fixes in Swig Web and `swig-cli` that predate 0.39

- **SHM replacement and silent counts** phased codons from the read's
  first base instead of the V alignment's, so about two reads in three
  counted them in the wrong frame. Rows without `sequence_frame` (IgBLAST
  output) now take the V germline's frame from `v_germline_start`.
- **Denoising** cut the V–J slice of a reverse-complemented read with
  reverse-strand coordinates from the read as given. The reverse
  complement now keeps IUPAC codes (R↔Y, K↔M, B↔V, D↔H), both strands are
  upper-cased alike, and `rev_comp` is read as T, TRUE or 1 everywhere
  (denoising, CHMMAIRRa, personalized germline inference).
- **Isotype.** The isotype column required `c_support ≤ 10⁻⁴` for C
  stretches under 30 bases, SwiftIG's E-value; the aligner's gate is
  10⁻³ on its own E-value, so some C calls had no isotype.
- **Alignment viewer.** A long `=`/`X` CIGAR now wraps inside its segment
  card instead of running under the Query/Germline/Score block.

## Behaviour that stays different from SwiftIG

- `stop_codon` and `productive` look for stop codons through the end of
  the J alignment (as IgBLAST does), not only to Trp/Phe118; TRAJ8\*02 and
  TRAJ51\*01, which carry an in-frame stop in FR4, give non-productive
  rearrangements.
- Random sequence gets no call (SwiftIG called a J on 65 of 100 random
  reads), and the J of a full V is not paired with a junk J of the V's
  locus when the read's J belongs to another locus.
- A V of 30–60 bases before a J of another locus can still pair with a J
  of its own locus that is significant on its own (10 of 168 such reads
  productive; SwiftIG 0): the aligner searches J only in the V's locus.
- J-mode CHMMAIRRa still misses a few J chimeras that SwiftIG flagged: the
  aligner prefers the 3′ parent J with a hard 5′ clip.

## Tests

- `engine/tests/guards.rs`: each guard on constructed reads, with and
  without the guards: a V that does not reach its junction, V and J beyond
  the junction span, a J under 10 columns, a D without a V or a J, a
  junction without aligned anchors, a short V before a chance J; C 7 bases
  off the splice; C floors.
- `tests/simulated.rs`: the junction floors reflect the anchor rule, and a
  new check requires the guards to report less than half the spurious
  junctions of the unguarded engine on the hard preset.
- `swig simulate` keeps drawing only C alleles that start within 3 bases of
  their gene, so simulations are byte-identical to 0.39.3's;
  `swig evaluate` also prints the spurious-junction rate.
- Web: SHM codon phase with and without `sequence_frame`, denoising of
  reverse-complemented and IUPAC reads, `rev_comp` spellings, isotype gate,
  CHMMAIRRa extension and its J threshold, the `c-min-identity` option of
  the WebAssembly module.
