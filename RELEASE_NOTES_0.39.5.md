# Swig 0.39.5

Swig 0.39.5 gives back the D of reads that stop at their D or start after
their V, which 0.39.4's plausibility guards removed, while D-free reads with
long random tails still rarely get one. Every other rule of 0.39.4 stays:
reads that show both a V and a J keep its V–J pair guards, hit floors,
junction and constant-region rules. It also fixes four problems that
`engine/docs/REVIEW.md` listed as open after 0.39.3, one of them only in
part:
- chance D calls on D-free reads that stop before their J or start after
  their V (in part: a chance D within about 30 bases of the V or J stays);
- IgT, IgZ and IgO constant alleles that got no locus;
- V names without `*`, which cost V accuracy on mutated reads;
- options that were silently ignored without the database they need.

It also gives a locus to alleles with custom names, which used to pair with
every locus, and fixes five gaps in how Swig Web and `swig-cli` prepare
references and report the aligner's warnings. Every fix has a regression
test that fails without it.

The aligner runs in Swig Web, in `swig-cli` (native and WebAssembly) and in
the `swig` binary. Every change the aligner itself makes to calls and AIRR
columns below applies to all three, given the same reference. Changes to
reading files, preparing references, options and messages name the route
they apply to.

With the shipped references, calls change against 0.39.4 only on reads
where a D has no J after it or no V before it: reads that stop at their D or
start after their V, and full-length reads whose J (or V) is not found or
fails the guards.
On the 29,500 simulated reads of the frozen suite:
- the IGH easy and mid, BCR and TCR sets, and 1,000 IGK reads run against
  the IGH database, give 0.39.4's AIRR byte for byte;
- on IGH hard, 3 reads gain a right D call and some posteriors move;
- on the IGH set at 30 % SHM, 250 full-length reads gain a D call: 240
  called with a V and no J (on 122 of them the guards drop the J that
  0.39.3 called; on the other 118, 0.39.3 found no J either) and 10 called
  with a J and no V (see "Full-length reads whose J the guards drop").

## What changes for users

- **A D at a read's open end is called again.** 0.39.4 called a D only
  between a V and a J that pair, as SwiftIG did. That removed every chance
  D of a read that stops before its J or starts after its V, and with them
  every real D of such reads. 0.39.5 calls a D that ends its path (no J
  after it) or starts it (no V before it), but makes it pay for the J (or
  V) the read would show in the bases beyond it if one were there. Two
  things set the cost:
  - how many such bases the read has, beyond a free stretch of 24 bases
    that N and P regions can fill;
  - the read's own measured divergence, so a mutated J that the read does
    show is not mistaken for a missing one.

  The cost has a ceiling (6.2 nats), so a database that lacks the read's J
  or V cannot erase every D. Such a D must also lie where a read without
  that J (V) would have searched it: it ends at most 60 bases after the V
  (starts at most 60 bases before the J). The D is still one of the
  alignments found in the window up to that J (back to that V), one per D
  allele. So a D that ends more than 60 bases after the V, or whose allele
  aligns better farther on, is not called even in part, where the same
  read without the J calls its part within the 60 bases. A V and a J are
  still paired only under 0.39.4's guards.

  With the frozen suite's IGH database (per read; `--calls fixed` in
  brackets where measured):

  | | 0.39.3 | 0.39.4 | 0.39.5 |
  |---|---|---|---|
  | right D, reads cut at the end of their D (4,309 IgM-like) | 90.2 % (90.1) | 0.0 % (0.0) | 92.2 % (92.2) |
  | the same, IgG-like (4,251) | 84.2 % (84.1) | 0.0 % (0.0) | 86.3 % (86.3) |
  | right D, reads ending 0–20 bases into their J (1,500 IgM-like / 1,500 IgG-like) | 91.2 / 83.1 % | 49.7 / 41.5 % | 91.7 / 83.8 % |
  | right D, reads starting after their V (1,913) | 93.1 % (93.5) | 0.0 % (0.0) | 93.1 % (93.5) |
  | D called on D-free reads, a V then 100 random bases (600) | 93.5 % (98.0) | 0.2 % (0.5) | 3.2 % (2.2) |
  | 100 random bases then a J | 97.7 % (99.5) | 0.0 % (0.0) | 9.5 % (10.7) |
  | a V then 60 random bases | 89.5 % (94.0) | 0.7 % (2.0) | 77.3 % (86.7) |
  | a V then 30 random bases | 72.0 % (81.8) | 0.0 % (0.3) | 72.8 % (82.8) |
  | a V then 10 or 25 random bases | 27.2, 63.7 % | 0.0, 0.3 % | 27.5, 64.8 % |
  | 30 random bases then a J | 79.8 % (88.0) | 0.0 % (0.0) | 79.3 % (88.2) |

  What this does not fix:
  - A chance D within about 30 bases of the V or J stays, as in 0.39.3. It
    sits exactly where the real D of a truncated read sits, so no length
    term can separate them; only a higher D price for every read could, at
    the cost of real D calls. 0.39.4 removed these chance Ds, and the real
    ones beside them.
  - Under `--profile truth` and `--profile igblast`, which charge no D
    price, most chance Ds stay.
  - The D-presence posterior is a marginal over paths: on 100 of 600 D-free
    reads of 100 random bases then a J, the call has no D but a D presence
    above 0.5 (0.39.4: 13). Under `--calls fixed`, a V then 100 random
    bases gets a mean D presence of 0.22 on the engine's synthetic test
    fixture (0.17 without the pair guards): a read that falls back to its
    V–D path has no chance-J paths in its posterior.

  What it costs, against 0.39.3 (0.39.4 called none of these Ds):
  - **Real D behind a long N region.** The cost assumes non-templated
    regions of at most about 24 bases. On 1,500 simulated reads per bin
    that end 0–20 bases into their J, right D calls fall by 1.20 points
    with an N2 of 26–30 bases (23 lost, 5 gained) and by 4.60 points with
    31–40 bases; with 20–25 bases they rise by 0.73 points. Reads that
    start inside a long N1 lose 0.07 points at 26–30 bases and 1.13 points
    at 31–40.
  - **A database with no J of the read's locus.** On IGH reads run against
    a database whose J alleles are all from another locus under custom
    names, right D calls fall from 4,408 to 4,003 (9.2 %; 0.39.4: 174).
    Dropping all J genes but IGHJ6 loses none, and nor does replacing the
    IGH J alleles with TRB J alleles under their own names, or removing the
    TRB J alleles from the TCR database.
  - **J posteriors at high SHM.** Paths without a J lose posterior mass, so
    a called J gains some, right or not: at 30 % SHM the J posterior rises
    on 802 of 2,597 reads with the same J call as 0.39.4 (falls on 8), and
    its Brier score on them goes from 0.2797 to 0.2814.

  `--d-open-flank off` turns the cost off and with it the D at a read's open
  end: on the suite's databases and the shipped references, that gives
  0.39.4's calls byte for byte. Names that 0.39.5 reads differently can
  still change calls: IGH constant genes of isotype F, O, T, W, X, Y or Z,
  which 0.39.4 left without a locus; alleles whose name, `locus=` key and
  header give no locus; and V names without `*`, which can now join a
  relative's gene. With `--no-pair-guards` as well, the calls are 0.39.3's
  apart from C and FWR3.
- **Full-length reads whose J the guards drop.** When a read's J fails the
  pair guards (a heavily mutated J, a chance J, or one more than 225 bases
  from the V), 0.39.4 called the V alone, 0.39.3 the V, a D and the J.
  0.39.5 may also call the V with a D within 60 bases of it, which pays the
  open-flank cost over the bases after it; a D further on is not called
  beside the V, as on the same read without the J. Unlike that read, it
  also loses a D that runs past the 60 bases, or whose allele aligns better
  past them; the read without the J calls the part within them.

  On the 3,000 IGH reads at 30 % SHM against 0.39.4, 250 reads gain a D
  call. 240 of them are called with a V and no J. On 122 of these 0.39.3
  had called the J; on the other 118 it found no J either and called the V
  and a D. The other 10 are called with a J and no V (5 had a V in 0.39.3).
  82 gain a right D (23 of the 122, 54 of the 118 and 5 of the 10; none
  lose one), and D allele accuracy rises from 37.76 to 40.61 %, the
  D-presence Brier score falls from 0.123 to 0.042. But most of the new
  calls are wrong: 155 name a wrong allele (96 of them away from the true
  D) and 13 fall on D-free reads (right "no D" calls 46.51 → 36.43 %),
  mostly because a heavily mutated V aligns short of its end and the D
  found within 60 bases of it is often a chance hit.

  Against 0.39.3, the 3,000 reads lose 28 right D calls and gain 2 (−0.91
  points). On all 28, 0.39.3 called a J that 0.39.5 does not. 21 of them,
  and both gains, are among the 250. On the other 7, neither 0.39.4 nor
  0.39.5 calls a D: on 6 the D lies more than 60 bases past the V, and the
  seventh gets no call at all.
- **IgT, IgZ and IgO constant genes.** An IGH constant gene is now
  recognised by any isotype letter: A, D, E, F, G, M, O, T, W, X, Y, Z
  (`IGHT2D`, `IGHZ`, `IGHO`).
  - After any of these letters except D, which also starts D gene names
    (`IGHD1-1`), at most one more letter may come before the next digit,
    `*` or other non-letter (`IGHMD`, `IGHGP`). A word such as `IGHTEST1`
    therefore names no locus.
  - That limit is new for A, E, G and M, which 0.39.4 accepted with any
    number of letters after them, so `IGHGPX1` no longer names IGH. No name
    in Swig's IMGT reference pack loses its locus.
  - In that pack, 8 constant alleles in 6 of its 64 species and strain sets
    (5 species) had no locus: IGHT (6), IGHZ (1) and IGHO (1). They are now
    IGH.
  - Calls on the shipped references do not change.
  - In a database of several loci, an allele of no locus could join any
    locus' segments and counted in every locus' C gate.
- **Alleles with custom names.** An allele whose name, `locus=` key and
  header name no locus (`novel_1`, `sample7_V3`) used to pair with every
  locus. In a database of several loci, that let light-chain and TRA reads
  get a D, and light-chain reads get a heavy C.

  The locus is now taken from, in order:
  1. the `locus=` key;
  2. the name;
  3. the header's gene names;
  4. the allele's closest named relative in the database, or in a small
     built-in panel of human and mouse V, J and C genes;
  5. failing all of these, the loci the database names; a D takes only
     those of them with D genes (IGH, TRB, TRD). In a database of several
     loci, a warning names these alleles and asks for a `locus=` key. When
     steps 1–4 place no allele of the database in a locus, the alleles
     still pair with every locus, now with a warning.

  `swig prepare` writes the result as `locus=`. On the suite's databases
  with custom names:

  | Custom-named alleles | 0.39.4 (as 0.39.3) | 0.39.5 |
  |---|---|---|
  | every D of the BCR database: light-chain reads with a D | 1,077 of 5,000 | 0, as with IMGT names |
  | every D of the TCR database: TRA reads with a D | 839 | 0; TRB reads still get TRD D alleles, as in 0.39.4: D allele accuracy is 3.63 points below IMGT names (76 reads), and 53 TRB reads with no D get one (see Limits) |
  | every J of the BCR database: calls against IMGT names | 4 junctions lost (0.39.3: C allele −0.74 points) | equal |
  | every V of the BCR database: calls against IMGT names | 4 junctions lost, as with custom J names; 7 of 7,500 reads change: 6 lose their junction, 4 list more tied alleles (0.39.3: same accuracy; 4 reads list more tied alleles) | equal on all 7,500 reads |
  | every allele of an IGH database: reads with no `locus` | all 5,000 | none |

  Limits:
  - D genes are too short to place by relatives. Custom-named TRB and TRD D
    alleles both fall back to TRB and TRD, so they need a `locus=` key. With
    the key, calls on the suite's 5,000 TCR reads equal those with IMGT names.
  - Loading 20 custom-named alleles takes about 12 ms more. Loading 705
    custom-named human BCR V alleles takes about 60 ms more, but a run on
    such a database is now faster: 7,500 BCR reads take 3.2 s instead of
    5.8 s on one thread.
- **V names without `*`.** A gene is the allele name's part before `*`.
  With no `*`, each allele counted as a gene of its own, which costs V
  accuracy on mutated reads.
  - A V name without `*` now joins the gene of its closest `*`-named
    relative within about 2 % sequence divergence. Otherwise it forms a
    gene with the other V alleles without `*` that it reaches through steps
    of about 2 % or less.
  - On the suite's IGH reads at 30 % SHM, with every `*` removed, V allele
    accuracy is now 0.10 points above IMGT names; in 0.39.4, as in 0.39.3,
    it was 8.20 points below.
  - With 20 % of names hidden, it is 0.03 points above; it was 1.77 points
    below in both.
  - Names with `*`, as in every shipped reference, keep their genes.
- **Custom references in Swig Web and `swig-cli`.**
  - A Swig Web reference cell is one locus, so a custom record in it
    (a bare `novel_1`, say) now gets that locus. Before, it was refused
    with "does not identify one of the supported IG/TR loci".
  - `swig-cli prepare-reference` and the config route place a custom
    record by its closest template and the locus words of its header. A
    record nothing places is refused with a message that asks for a
    `locus=` key.
  - The `swig-cli --vdj` route rewrote every germline header to the bare
    name, which dropped a `locus=` key. Custom TRB D alleles then gave TRD
    reads a TRB D. It now keeps the key.
  - The later template tiers of reference preparation read each record's
    original header, so a record placed in the first tier is not refused in
    a later one.
  - Reference preparation in Swig Web and `swig-cli` accepts IMGT's IgD
    constant gene (`IGHD*01`) and J–C orphon records (`IGLJ-C/OR18*01`) in
    a C database. Before, they were refused as a D or J gene.
- **Options without their database.** Two options were ignored when the
  run had no database for them; they are now refused, as `--c-evalue`
  already was:
  - `--c-prefix-identity` above 0 needs a C database;
  - `--double-d all` or `long-span` needs a D database.

  Where the refusal happens on each route:
  - `swig` binary: exit 2 with the usage, before any database is read.
  - WebAssembly module: `swig_init_database` returns −1 with the message in
    the error buffer. When the option is set after the database is loaded,
    the next `swig_annotate` does this instead.
  - Swig Web: refuses the C gate before a run without a C reference.
  - `swig-cli`: on the config route, `annotation.minimumConstantPrefixIdentity`
    above 0 needs C germlines and `annotation.doubleD` (`all` or
    `long-span`) needs D germlines. On the `--vdj` route,
    `--minimum-c-prefix-identity` above 0 needs a C database
    (`-c_region_db`); that route has no double-D screen. `swig-cli` checks
    these before it starts the aligner, so its message is the same with
    `--aligner native` and `--aligner wasm`.

  Options that only tune a stage the run lacks are run with a warning. These
  are the D options without a D database (`--min-d-match`, `--d-mismatch`,
  `--d-presence-rate`, `--d-allele`, `--no-anchor-window`,
  `--d-open-flank off`), 0.39.4's `--c-min-identity` without a C database,
  and the double-D rules without `--double-d`. So is
  `--d-allele max-marginal` with `--no-posterior`, and a D allele too short
  for the double-D screen's seeds (TRDD1 and TRDD2 at the default seed
  length). These warnings change no calls: a warned option is ignored, so
  the run writes the same AIRR as without it.
- **Aligner warnings are shown.** Swig Web and `swig-cli` never showed the
  aligner's database and run warnings: unknown or fallback loci, names
  without `*`, J and C convention problems, and the warnings above.
  - `swig-cli` now prints each once, as `Aligner warning: …`, on both
    routes. The pipeline also lists them in `<prefix>.summary.json` under
    `references.alignerWarnings`.
  - Swig Web shows them in the run summary.


## Output

The native and WebAssembly routes still write byte-identical AIRR (release
smoke tests). Against 0.39.4, on the frozen suite with default settings,
paired on the same reads and scored by the 0.39.1 evaluator:

| Set | Reads | Rows that change a call | V, D, J, C allele and junction |
|---|---|---|---|
| IGH easy | 3,000 | 0 (same bytes) | unchanged |
| IGH mid | 5,000 | 0 (same bytes) | unchanged |
| IGH hard | 5,000 | 3 (posteriors also move) | D allele +0.06 points (3 gained) |
| IGH 30 % SHM | 3,000 | 250 | D allele 37.76 → 40.61 % (82 gained, none lost); right "no D" 46.51 → 36.43 % |
| BCR (IGH+IGK+IGL) | 7,500 | 0 (same bytes) | unchanged |
| TCR (TRA+TRB+TRD+TRG database) | 5,000 | 0 (same bytes) | unchanged |
| IGK reads against the IGH database | 1,000 | 0 (same bytes) | — |

Against 0.39.3 the losses are 0.39.4's: junctions from aligned anchors
(IGH mid −0.49 points, hard −3.46, 30 % SHM −40.94, BCR −0.61, TCR 1
read), its C rules (C allele on IGH easy −1.27 points, mid −1.50, hard
−1.90, 30 % SHM −5.80, BCR −2.11, TCR −1.86, and 1 C called on a C-free
IGH hard read) and its pair guards (at 30 % SHM, J allele −1.80 points with
57 lost and 3 gained, and 1 V lost; on IGH hard, 1 J lost). At 30 % SHM the
reads whose J the guards drop (above) also lose 28 right D calls and gain 2
(−0.91 points; 0.39.4 lost 108). Their D allele and D-presence Brier
scores rise from 0.2085 to 0.2109 and from 0.0367 to 0.0424 (paired t 3.21
and 3.52), almost all of it on those reads.

Other real and simulated data:
- **Macaque.** On two macaque repertoires of 12,229 and 13,698 D calls, no D
  call changes against 0.39.4.
- **Shipped references.** On 11 sets of 2,000 reads simulated from the IMGT
  pack's references (human and mouse BCR and TCR, six sets of four fish
  species and one of the platypus), 1 call of 22,000 changes against
  0.39.4: a platypus read with no J call gains a D, the true allele tied
  with one other. Elsewhere only J posteriors move, on 17 reads.

## Timing

With databases whose alleles name their locus, as every shipped reference
does, speed is unchanged. Instruction counts (callgrind, one thread)
against 0.39.4:
- 1,000 IGH hard reads: +0.38 %;
- the 600 reads of a V followed by 100 random bases: −13.5 %.

Loading a database does a little more work even when every allele names
its locus: each header is checked for a `gene=` key, and the open-flank
cost is set up for each locus. On the suite's human references that adds
about 0.2 % to the load's instructions for IGH, 0.3 % for BCR and 0.9 % for
TCR (callgrind, one thread, no reads; each figure moves by up to 0.1 point
with memory layout), under 1 ms. Alleles that name no locus add more (see
"Alleles with custom names").

## For developers

- New `swig` option: `--d-open-flank on|off` (default `on`). The WebAssembly
  module takes it through the generic option call.
- The posterior's forward–backward, the exhaustive oracle and the DAG all
  include the new cost. Property tests check that they agree in both call
  models, with the cost on and off, and under the pair guards:
  `guarded_calls_with_open_flanks_equal_the_exhaustive_subset_optima`
  checks every guarded call against the exhaustive optimum of the subset
  that gave it, and against the rule written out on its own
  (`CombineRule`). `PlacementDag::path_score` and `PlacementDag::posterior`
  are now crate-internal. `Adjudicator::adjudicate` charges the open end
  only up to the last placement; callers that know the read's length use
  `adjudicate_at`. `DWindowPolicy::open_reach` gives the reach of an
  open-end D (`ExpandedResidual`: its `reach`, 60).
- When the best path breaks the guards, the candidate subsets are 0.39.4's
  and, under the cost, V placements with the D placements within their
  reach and D placements with the J placements they lie within the reach
  of. Of 435 calls of the oracle test's fixtures, a path the rule admits
  beats the call on none under the per-read model and on 1 under the fixed
  scale (0.39.4's subsets: the VDJ optimum's V clip stretches the span past
  225 bases and the whole subset is passed over).
- New module `locus_infer` and the embedded panel `engine/data/locus_panel.fa`
  (269 human and mouse V, J and C records, built from the shipped IMGT pack
  by `scripts/build-locus-data.mjs`). `public/swig.wasm` grows by about
  0.16 MB.
- `RunOptions::check_databases` returns an option's refusal as an error,
  or else its warnings, and `Run::warnings` returns the run's warnings. The
  WebAssembly summary adds the run's warnings to the database's in
  `warnings`. A refusal fails `swig_init_database` (or, after an option
  change, the next `swig_annotate`) with its message in the error buffer.

## Still open

`engine/docs/REVIEW.md` §0 lists what this release does not fix. Among them:
- chance Ds within about 30 bases of the V or J of a truncated read;
- wrong D calls beside a heavily mutated V whose J the guards drop (above);
- a chance J between the V and a D at the read's end, which ends the D
  window before that D (on simulated IGH reads cut after their D, about 5.5
  points of right D calls);
- real Ds behind N regions longer than about 24 bases at a read's open end,
  and the J-posterior calibration at high SHM described above;
- custom-named D alleles in a database holding TRB and TRD, which need a
  `locus=` key;
- the double-D screen's seed length, which cannot see TRDD1 and TRDD2 (it is
  now reported, not changed);
- the strand test's null model, posterior overconfidence at high SHM, and
  the other items carried over from 0.39.3 and 0.39.4.
