# Swig 0.39.5

Swig 0.39.5 fixes four problems that `engine/docs/REVIEW.md` listed as open
after 0.39.3, one of them only in part:
- chance D calls on D-free reads that stop before their J or start after
  their V (in part: a chance D within about 30 bases of the V or J stays);
- IgT, IgZ and IgO constant alleles that got no locus;
- V names without `*`, which cost V accuracy on mutated reads;
- options that were silently ignored without the database they need.

It also gives a locus to alleles with custom names, which used to pair with
every locus.

Checking these fixes end to end also turned up five gaps in how Swig Web and
`swig-cli` prepare references and report the aligner's warnings. Those are
fixed too. Every fix has a regression test that fails on 0.39.3.

The aligner runs in Swig Web, in `swig-cli` (native and WebAssembly) and in
the `swig` binary. Every change the aligner itself makes to calls and AIRR
columns below applies to all three, given the same reference. Changes to
reading files, preparing references, options and messages name the route
they apply to; the fixes to how Swig Web and `swig-cli` prepare custom
references also change calls, but only on the routes they name.

With the shipped references, calls change only on reads where 0.39.3
called a D with no J after it or no V before it. On full-length reads that
happens only when the aligner misses the J or V. On the 29,500 simulated
reads of the frozen suite:
- the IGH easy, mid and hard, BCR and TCR sets, and 1,000 IGK reads run
  against the IGH database, keep every call. On IGH hard only posteriors
  move;
- on the IGH set at 30 % SHM, 14 rows change a call (13 D, 1 J), all reads
  in which 0.39.3 found no J. No right D call is lost, and right "no D"
  calls rise from 33.33 % to 34.11 %.

With databases whose alleles name their locus, as every shipped reference
does, speed is unchanged.

## What changes for users

- **Chance D on reads without a D.** A read that stops before its J, or
  starts after its V, used to get a D from a chance match in its random
  bases on most reads. The D was priced for a junction-sized window, and the
  bases beyond such a read's last segment are not one.

  A D that ends its path (no J after it), or starts it (no V before it), now
  pays for the J (or V) the read would show in the bases beyond it if one
  were there. Two things set the cost:
  - how many such bases the read has, beyond a free stretch of 24 bases
    that N and P regions can fill;
  - the read's own measured divergence, so a mutated J that the read does
    show is not mistaken for a missing one.

  The cost has a ceiling, so a database that lacks the read's J or V cannot
  erase every D. On 600 D-free IGH reads per row:

  | D-free reads | 0.39.3 | 0.39.4 |
  |---|---|---|
  | a V followed by 100 random bases | 93.5 % get a D | 6.8 % |
  | 100 random bases followed by a J | 97.7 % | 10.2 % |
  | a V followed by 60 random bases | 89.5 % | 73.5 % |
  | a V followed by 10, 25 or 30 random bases | 27.2, 63.7, 72.0 % | unchanged |
  | 30 random bases followed by a J | 79.8 % | unchanged |

  Under `--calls fixed` the first two rows read 98.0 → 6.5 % and
  99.5 → 10.8 %.

  What this does not fix:
  - A chance D within about 30 bases of the V or J stays. It sits exactly
    where the real D of a truncated read sits, so no length term can
    separate them; only a higher D price for every read could, at the cost
    of real D calls.
  - Under `--profile truth` and `--profile igblast`, which charge no D
    price, most chance Ds stay. Behind a V: `truth` 93.7 → 92.8 %,
    `igblast` 100 → 99.5 %.

  What it costs:
  - **Real D behind a long N region.** The cost assumes non-templated
    regions of at most about 24 bases. On 1,500 simulated reads per bin
    that end 0–20 bases into their J, right D calls fall by 1.07 points
    with an N2 of 26–30 bases (17 lost, 1 gained) and by 3.80 points with
    31–40 bases. With 20–25 bases the share of right D calls does not
    change (2 lost, 2 gained), but 5 of the 1,500 reads lose their D call
    and the D-presence Brier score rises from 0.0056 to 0.0126. Reads that
    start inside a long N1 lose 0.07 points at 26–30 bases and 1.13 points
    at 31–40.
  - **D-presence calibration on truncated reads.** On 1,500 simulated
    reads per set that end inside their J, 3 D calls change per set and
    right D calls stay at 91.20 % (mid) and 83.13 % (hard, 1 lost and 1
    gained), but the D-presence Brier score rises by 0.0045 (mid) and
    0.0014 (hard).
  - **J posteriors at high SHM.** Paths without a J lose posterior mass, so
    a called J gains some, right or not. At 30 % SHM the J posterior rises on
    908 of 2,742 reads with the same J call, and on those reads its Brier
    score goes from 0.1870 to 0.1887 (paired t = 3.8). No right J call
    changes; one read that had no J now gets a wrong one.
  - **A database with no J of the read's locus.** On IGH reads run against a
    database whose J alleles are all from another locus under custom names,
    7.3 % of right D calls are lost.
    - Dropping all J genes but IGHJ6 loses none.
    - So does replacing the IGH J alleles with TRB J alleles under their
      own names, or removing the TRB J alleles from the TCR database.

  `--d-open-flank off` turns the cost off. On the suite's databases and the
  shipped references, that gives 0.39.3's calls. Names that 0.39.4 reads
  differently can still change calls: IGH constant genes of isotype F, O,
  T, W, X, Y or Z, which 0.39.3 left without a locus; alleles whose name,
  `locus=` key and header give no locus; and V names without `*`, which can
  now join a relative's gene.
- **IgT, IgZ and IgO constant genes.** An IGH constant gene is now
  recognised by any isotype letter: A, D, E, F, G, M, O, T, W, X, Y, Z
  (`IGHT2D`, `IGHZ`, `IGHO`).
  - After any of these letters except D, which also starts D gene names
    (`IGHD1-1`), at most one more letter may come before the next digit,
    `*` or other non-letter (`IGHMD`, `IGHGP`). A word such as `IGHTEST1`
    therefore names no locus.
  - That limit is new for A, E, G and M, which 0.39.3 accepted with any
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

  | Custom-named alleles | 0.39.3 | 0.39.4 |
  |---|---|---|
  | every D of the BCR database: light-chain reads with a D | 1,077 of 5,000 | 0, as with IMGT names |
  | every D of the TCR database: TRA reads with a D | 839 | 0 |
  | every J of the BCR database: C allele against IMGT names | −0.74 points | equal |
  | every V of the BCR database: calls against IMGT names | same accuracy; 4 of 7,500 reads list 1 or 2 extra tied alleles | equal on all 7,500 reads |
  | every allele of an IGH database: reads with no `locus` | all 5,000 | none |

  Limits:
  - D genes are too short to place by relatives. Custom-named TRB and TRD D
    alleles both fall back to TRB and TRD, so they need a `locus=` key.
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
    accuracy is now 0.10 points above IMGT names; it was 8.20 points below.
  - With 20 % of names hidden, it is 0.03 points above; it was 1.77 points
    below.
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
  - WebAssembly module: refuses to start.
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
  `--d-open-flank off`) and the double-D rules without `--double-d`. So is
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
smoke tests). Against 0.39.3, on the frozen suite with default settings,
paired on the same reads and scored by the 0.39.1 evaluator:

| Set | Reads | Rows that change a call | V, D, J, C allele and junction |
|---|---|---|---|
| IGH easy | 3,000 | 0 | unchanged |
| IGH mid | 5,000 | 0 | unchanged |
| IGH hard | 5,000 | 0 (posteriors only) | unchanged |
| IGH 30 % SHM | 3,000 | 14 (13 D, 1 J) | unchanged; right "no D" 33.33 → 34.11 % |
| BCR (IGH+IGK+IGL) | 7,500 | 0 | unchanged |
| TCR (TRA+TRB+TRD+TRG database) | 5,000 | 0 | unchanged |
| IGK reads against the IGH database | 1,000 | 0 | — |

Other real and simulated data:
- **Macaque.** On two macaque repertoires of 12,229 and 13,698 D calls, no D
  call changes.
- **Reads that start after their V.** On 1,913 simulated IGH reads that hold
  a D and a J but no V, no call changes; only some D and J posteriors
  move.
- **Shipped references of other species.** Calls on the IMGT pack's sets do
  not change. Some posteriors change in their last digits, from the
  chance-D cost.

## Timing

With databases whose alleles name their locus, as every shipped reference
does, speed is unchanged. Instruction counts (callgrind, one thread)
against 0.39.3:
- 1,000 IGH hard reads: +0.39 %;
- the 600 reads of a V followed by 100 random bases: −0.46 %.

Wall time of the `swig` native binary on the IGH mid and hard sets (5,000
reads each, human IGH V, D, J and C): each figure is the median of 5
interleaved runs on a 4-core cloud VM, database preparation included.

| Set | Threads | 0.39.3 | 0.39.4 |
|---|---|---|---|
| IGH mid | 1 | 3.25 s | 3.19 s |
| IGH hard | 1 | 3.13 s | 3.05 s |
| IGH mid | 4 | 0.97 s | 1.00 s |
| IGH hard | 4 | 0.91 s | 0.91 s |

Single runs varied by up to 20 % from the median (one 0.39.3 run by 49 %),
so the differences are noise. Loading a database does a little more work
even when every allele names its locus: each header is checked for a
`gene=` key, and the open-flank cost is set up for each locus. On the
suite's human references that adds about 0.2 % to the load's instructions
for IGH, 0.4 % for BCR and 0.9 % for TCR, under 1 ms. Alleles that name no
locus add more (see "Alleles with custom names").

## For developers

- New `swig` option: `--d-open-flank on|off` (default `on`). The WebAssembly
  module takes it through the generic option call.
- The posterior's forward–backward, the exhaustive oracle and the DAG all
  include the new cost. Property tests check that they agree in both call
  models, with the cost on and off. `PlacementDag::path_score` and
  `PlacementDag::posterior` are now crate-internal. `Adjudicator::adjudicate`
  charges the open end only up to the last placement; callers that know the
  read's length use `adjudicate_at`.
- New module `locus_infer` and the embedded panel `engine/data/locus_panel.fa`
  (269 human and mouse V, J and C records, built from the shipped IMGT pack
  by `scripts/build-locus-data.mjs`). `public/swig.wasm` grows from 1.0 to
  1.16 MB.
- `RunOptions::check_databases` returns an option's refusal as an error,
  or else its warnings, and `Run::warnings` returns the run's warnings. The
  WebAssembly summary adds the run's warnings to the database's in
  `warnings`. A refusal fails `swig_init_database` (or, after an option
  change, the next `swig_annotate`) with its message in the error buffer.

## Still open

`engine/docs/REVIEW.md` §0 lists what this release does not fix. Among them:
- chance Ds within about 30 bases of the V or J of a truncated read;
- real Ds behind N regions longer than about 24 bases at a read's open end,
  and the J-posterior calibration at high SHM described above;
- custom-named D alleles in a database holding TRB and TRD, which need a
  `locus=` key;
- the double-D screen's seed length, which cannot see TRDD1 and TRDD2 (it is
  now reported, not changed);
- the strand test's null model, posterior overconfidence at high SHM, and
  the other items carried over from 0.39.3.
