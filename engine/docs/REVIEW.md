# Code review: soundness, design and the WebAssembly port

A review of the whole crate (≈30,600 lines at `31c0c0f`, when it was the
standalone vdj_assign project), done before adding the WebAssembly build.
Five subsystems were reviewed line by line in parallel: adjudication,
alignment kernels, search, data and IO, and
simulation and tests. The composition layer (pipeline, configuration,
batching and the binary) was read alongside them. Every finding below
marked *confirmed* was reproduced, or follows directly from the code
quoted. Findings marked *suspected* were reasoned about but not
reproduced. File and function names are given rather than line numbers,
since the fixes moved lines.

## 0. Status after the 0.39.1 audit and the 0.39.2 review

A second audit, made after Swig 0.39.1, covered the input layer, retrieval,
alignment, adjudication, the pipeline, the interfaces and the simulator. It
re-checked this review's open items as well. Its fixes are one commit each,
with a regression test that fails without them. Of this review's open
items, these are now fixed:

| Item | Fix |
|---|---|
| 3.1.1 D traces taken before admission | admission is the D stage's acceptance predicate |
| 3.1.2 residual adjudicator | removed, with `--residual` |
| 3.1.3 `prepare` round trip | a prepared database reads back with the trust, loci and provenance its annotation gave |
| 3.1.4 pseudogene-heavy synthetic reads | `simulate` draws only alleles a rearranged transcript carries, within one locus |
| 3.1.5 `evaluate` and `rev_comp` | allele identity no longer depends on the reported orientation |
| 3.1.7 inverted D window | never proposed |
| 3.1.9 silent 64-allele cap | every allele of a kept V gene is a candidate |
| 3.1.10 max-marginal D decoding | the decoded D carries its allele and is kept only if its posterior holds |
| 3.1.12 frame ties, σ series, trace memory, BOM and non-UTF-8 input, `--tail`, timing test, floors | fixed (σ as two binomial tails; scores up to 30,000 in lanes with a bound per traceback group) |
| 4.1–4.3, 4.7, 4.8 | one acceptance mechanism; the adjudication mode as a validated type; the lattice recursion written once; typed notes; one test kit |
| 4.5 (part) | the tandem-D grammar and the residual path are gone |
| §0 13 (part) IGHT, IGHZ and IGHO without a locus | every IGH isotype letter names IGH; a custom-named V, J or C takes its relatives' locus (the database's named alleles, else an embedded panel, which places a J or C only in loci of the database's V), a custom D the loci with D genes, with a warning asking for `locus=`; C follows the locus the path shares; Swig Web and swig-cli write `locus=` for custom records |

The audit's own fixes are listed in the commit log; the README and
ALGORITHM describe the behaviour they changed. Still open, by priority:

1. **Strand test** (`orient`). Its null is random sequence, which the
   antisense strand of an IG read is not, and its fallback ignores J:
   reads with very high SHM or a short V can lose their orientation.
   Needs a model and a measurement.
2. **Posterior calibration under SHM.** At 30 % SHM the top bin
   (posterior ≥ 0.999) is 50–80× overconfident: the posterior is exact
   over its candidates, but the V/J keep policy has no score margin and the
   report credits shorter co-optimal placements. D presence is also
   overconfident because every D trim gets full prior weight. A
   posterior-only `ln T(p̂)` term was measured and not adopted: it would
   break "the call is a mode of the posterior", and a constant shift fits as
   well. An N/P and trim-length prior should replace it.
3. **Anchor gate and wrong V calls.** A call now needs a V or J whose clip
   on the path is significant on its own. At 30 % SHM, all 16 calls this
   drops had a wrong or missing V gene, but 12 had the right J gene. With the J presence
   cost, the set loses 0.43 pp of J and 0.46 pp of C accuracy. Keeping J and
   C when only V fails the gate is worth measuring.
4. **Annotation transfer.** Sibling transfer copies region starts through a
   local alignment with no frame or motif check. The 3′-partial
   IGHV1/OR16-3/-4/-22 alleles get a low-confidence Cys104 in a non-coding
   frame. Five V alleles keep a misplaced given anchor:
   - IGHV1/OR15-5\*01 and rhesus IGHV1-197\*02 have no confident motif to
     replace it;
   - mouse TRAV7-3\*02 and TRAV9D-4\*03 and sheep TRAV45-3\*02 escape the
     frame test, since their truncated region starts share the anchor's
     frame.
5. **Public fields bypass validation** (3.1.8): literal `GapCosts` with
   positive costs, `Allele` regions, `DistributedEvidence` (and with it
   3.1.6, `lossless_lmer` fixed at 4).
6. **FR/CDR from another allele than `v_call`'s first name** (3.1.11):
   regions project through the called placement's primary allele.
7. **Throughput.** The calling thread reads, builds rows and writes alone,
   leaving three cores idle for 10–13 % of wall time at 4 threads.
   Adjudication scratch and about 360 allocations per read remain.
8. **Structure.** `segment.rs` and `annotate/mod.rs` should be split
   (4.4); dead API remains (`InterSequenceGotoh::scalar_only`,
   `KmerIndex::build_from`, `FixedBudget`, `GeneFrames::is_empty`); the
   public surface is wider than any caller needs (4.6); ALGORITHM.md still
   carries review narrative (4.9).
9. **Swig's shipped references** (outside the aligner). `SWIGMETA`
   misplaces the anchors of 50 J allele records (TRAJ24*01–03 and
   TRAJ48*01 among them) and gives several 5′-partial V alleles the wrong
   frame. At load the engine replaces all 50 J anchors and 8 V anchors with
   the motif hit; five V alleles keep theirs (item 4). The references
   should be regenerated.

A third review, after Swig 0.39.2, found 2 major and 13 minor defects. Swig
0.39.3 fixes 13 of them outright and 2 in part, each with a regression test
that fails on 0.39.2; `RELEASE_NOTES_0.39.3.md` lists them. Still open from
that review:

10. **Chance D on reads without a J or a V.** Resolved for chance Ds more
    than about 30 bases from the anchor: a D that ends its path (no J after
    it) or starts it (no V before it) now pays for the J (V) its read would
    show in the bases beyond it (ALGORITHM §2.3, `--d-open-flank`). This
    holds with a derived D prior (the default profile and `sensitive-d`);
    `--profile truth` and `igblast` charge no D prior, so a chance D of a
    dozen points outweighs the cost's 6.2-nat ceiling and most stay (92.8
    and 99.5 % behind a V). A V followed by 100 random bases gets a D in
    6.8 % of D-free reads instead of 93.5 %, 100 random bases followed by
    a J in 10.2 % instead of 97.7 % (`--calls fixed`: 6.5 % and 10.8 %
    instead of 98.0 % and 99.5 %), 60 random bases after the V in 73.5 %
    instead of 89.5 %. The
    full-length suite sets keep every call; igh_hard, igh_shm30 and the
    cut reads lose no D allele, and the macaque lineage reads change no D
    call. Open by design for the near ones: 10, 25 and 30 random bases
    after the V still give a D in 27.2, 63.7 and 72.0 % of reads, and 30
    before the J in 79.8 %. Such a chance D sits where the real D of a
    truncated read sits, so no length or window term separates them; only
    `P*_D` itself moves them, along the curve of D accuracy against D-free
    specificity (ASSESSMENT §9.9). The cost assumes N/P regions of at most
    about 24 bases: real D calls behind a 26–30-base N2 at a read that
    ends in its J drop by 1.1 pp, behind 31–40 bases by 3.8 pp (N1 of
    31–40 bases at a read that starts there: 1.1 pp); reads that end
    inside their J keep their D calls, at a D-presence Brier 0.0045 (mid)
    and 0.0015 (hard) higher. The cost also moves the J posterior: V-D
    paths without a J lose mass, so a called J gains it whether right or
    not. On igh_shm30 the J posterior rises on 908 of 2,742 rows with the
    same J call (falls on 1), and its Brier score goes from 0.1870 to
    0.1890 (paired t = 3.8); 18 of the 21 rows that move by more than 0.05
    are wrong J calls (igh_hard: t = 0.8). With the frozen suite's IGH
    database the V side prices less than the J side: a V stays hidden in
    12 visible bases (accept score 24), a J in 8 (16), as measured on reads
    cut inside them, and a read whose only anchor is its J reads the V's
    divergence from a few dozen columns. That is why 100 random bases
    before a J keep more chance Ds (10.2 %) than 100 after a V (6.8 %). The
    N1/N2 gap hinge first proposed here was measured and rejected (−8 to
    −12 pp D allele at 30 % SHM); a geometric N-length allowance traced the
    same trade-off as the fixed one.
11. **Gene identity from allele names.** *Resolved:* a V name without `*`
    joins the gene of a `*`-named relative within about 2 %, or groups with
    other such names at that level (`gene=`; ALGORITHM §3a). The 8.2 pp V
    allele loss at 30 % SHM with every `*` replaced is now +0.10 pp (not
    significant); names with `*` keep their genes.
12. **Junk V by circular divergence.** With the V presence cost, 7 of 1,200
    V-less test reads still get a V: the read's divergence is estimated
    partly from the junk V's own columns, which raises its per-read value.
13. **Smaller items.** `write_batches` drops the valid reads that precede a
    malformed record in the same batch;
    `join_indel` compares weighted votes after the 7-mer fallback, not the
    primary seed counts its guard describes; `swig_version` and
    `AlignerEngine.version()` have no caller.
14. **V keep on reads holding two chains.** The V stage keeps the top 3
    alignments over the whole read, so a lightly mutated light V can crowd
    every heavy V out (seen on synthetic scFv reads with a mutated heavy
    chain); the read is then called as its light chain. The J stage keeps a top 3 per locus;
    the V stage could too.
15. **swig-cli and prepared databases** (outside the aligner). With
    `prepareMetadata` on, swig-cli rejects a V database written by
    `swig prepare` ("contains invalid SWIGMETA coordinates").
16. **The double-D screen cannot see TRDD1 or TRDD2.** It indexes only D
    alleles at least `--double-d-seed-length` long (11 by default), and its
    rescue words come from those alleles too, so the 8- and 9-nt TRDD1 and
    TRDD2 are never seeded; every human TRD tandem holds one of them. On
    4,000 simulated TCR reads it found 0 of 199 TRD tandem-D reads at seed
    11 and 38 (19 %) at seed 8 (TRB: 52 % and 83 %). A run now names the
    alleles the screen cannot see. A per-locus default seed length (or
    indexing short alleles at their own length) needs a false-positive
    measurement on single-D TRD reads first.

## 1. Verdict

**Is the algorithm sound?** Yes. The model and its search hold up under
independent checks that go beyond the existing test oracles:

- **The placement-DAG Viterbi.** It was fuzzed against the exhaustive
  oracle on 60 random reads with gapped alignments. The checks covered
  every grammar, ι ∈ {0, −1, −3}, three divergences and the fixed scale,
  about 3,600 in all. Mismatches: 0.
- **The forward–backward pass.** It was compared with a log-space oracle
  up to ln Z ≈ 1,250 nats. Largest marginal error: 10⁻¹⁰.
- **The sibling bounds** (`family_bound`). 3,600 extra randomised cases
  were run, with hotspot profiles, E-value floors, both tie orders, caps
  of 2, 5 and 64, and truncated alleles. The bounded ranking equalled
  full scoring every time.
- **The clip index.** It is exact.
- **The integer per-read model.** It is the same for the call and the
  posterior, so the call is a mode of the reported distribution.
- **The rest of the search.** The D′ shortcut in the Gotoh recurrence,
  band-relative SIMD coordinates, i16 lane admission, the anchor-column
  bound for C, reach-class exactness, radix packing and the CSR index all
  check out.

The defects are at the edges, not in the core. A traceback shortcut
dropped a D candidate, AIRR N/P regions were mislabelled, and malformed
input reached the TSV unchecked. These and the other edge defects are
fixed in this branch (§2) or listed in §3.

**Is the code properly designed?** Largely, yes.

- Each heuristic is a trait with a default and alternatives, and the
  pipeline monomorphises with no trait objects.
- The SIMD kernel is written once over a lane trait. A new instruction set
  is one ~90-line `Lanes` implementation; the `simd128` backend added here
  was exactly that.
- Configuration is plain data, and its lossless settings are derived, not
  tuned.
- The documentation states each search step's guarantee: exact, lossless,
  bounded or heuristic.

The weaknesses are listed in §4. The main ones:

- two acceptance mechanisms for one stage;
- a cluster of interacting mode flags on the DAG that leak into the
  construction of the D stage;
- the Viterbi, forward and backward recursions written three times;
- public fields that bypass validated constructors.

**Legacy or repeated code?** Some of each.

- *Legacy:* knobs kept for compatibility with earlier output, such as
  `TieOrder::GeneLeadersFirst` ("the order the V stage has always
  produced"). Measured-and-rejected alternatives are still compiled in: the
  residual adjudicator, open syncmers and minimisers, and the ungapped
  prefilter, which sits in the default V retriever type while switched
  off.
- *Repeated:* the lattice edge walk (three copies); job construction in
  `segment.rs` (four copies); the J rescue in `annotate`, which repeats V
  pass 4; synthetic-assigner setup across tests (about ten copies).

**Monolithic files?** Few. Most files are cohesive and below 600 lines.
The exceptions are:

| Where | What |
|---|---|
| `annotate::annotate` | one 481-line function running four passes over shared vectors |
| `segment.rs` | 1,600 lines in layers: preparation, policy, ranking, bounding, waves, finishing |
| `sim::simulate::attempt` | 159 lines |
| `Display for sim::Evaluation` | 144 lines |
| the examples' `main` functions | 131–219 lines each |

**Ground-up work or patchwork?** Ground-up in structure, patchwork in
places.

- **Ground-up.** The architecture is principled and was clearly designed
  as a whole: one objective, one lattice, one kernel, guarantees per step.
  Several earlier ad hoc mechanisms were replaced rather than layered
  over; the README's "removed note mechanisms" lists them.
- **Patchwork.** The iteration history shows in the places listed under
  *legacy* above. Rejected options stay in the code instead of being
  removed. A flag set on one component changes another
  (`posterior`/`calls` decide the D stage's floor). ALGORITHM.md still
  narrates review history ("An earlier version…", "A second review
  found…").
- **Before this branch, the application layer.** Mapping options to the
  configuration, phrasing warnings, choosing the adjudicator and streaming
  AIRR lived in the binary, where a second front end would have had to
  copy them. This branch moves it into the library (`run.rs`).

## 2. Fixed in this branch

Each fix has its own commit and a regression test that fails without it.

| Finding | Effect | Fix |
|---|---|---|
| **D candidate dropped by the score-skipping shortcut** (`segment.rs`, `needs_scores`) | With at most `top + 1` candidates, scoring was skipped and every candidate traced in retrieval order. The wave stops at `cap`, so with `cap < top + 1` the last candidate was never aligned. D runs with `top = cap = 32`, so a D window with exactly 33 candidates lost one allele whatever its score. | Score whenever `cap` binds below `top + 1`. No call changed on 10,000 synthetic reads; only some D posterior digits did. |
| **AIRR `np1`/`np2`/`np3` filled by position** (`report.rs`) | A read without a V reported its D–J region as `np1`; AIRR defines `np1` after V, `np2` after the first D and `np3` after the second. | Name each region by the segment it follows. |
| **FASTQ sequence bytes unchecked** (`reads.rs`) | A tab inside a FASTQ sequence line reached the AIRR `sequence` column and shifted every later field of the row. | Reject sequence lines holding whitespace or a control character, naming the line and the symbol; no-call symbols such as `.` still read as N. |
| **Deduplication dropped header annotation** (`germline.rs`) | A merged copy's `anchor=` and region starts were discarded, against "a header value is never overridden". A name could also be attached to two sequences (`>A X`, `>B Y`, `>A Y` was accepted), and aliases were never checked. `GermlineSet::push` accepted an anchor past the sequence end, which later underflowed. | The kept allele takes what it lacks from its copies. Every name and alias is checked. `push` validates the anchor. |
| **Slow default kernel on WebAssembly without SIMD** (`Backend::detect`) | Detection fell back to portable lanes, which there run about 4× slower than the scalar kernel because nothing vectorises them. A catch-all dispatch arm would also have dropped jobs silently for a backend not built for the target. | Fall back to the scalar kernel on WebAssembly without `simd128`; other targets keep portable lanes. Name every backend in dispatch and fail loudly on one not built. Remove two dead `#[allow]`s. |
| **Probability text depended on the platform's `exp`** (`report.rs`, `probability`) | A distance from 0 or 1 below 10⁻⁴ printed with ten decimals even when it rounds away: `1 − 10⁻¹⁶` printed as `1.0000000000` and an exact 1 as `1.000000`. Platforms whose `exp` differs in the last bit (glibc against WebAssembly's portable libm) printed different text for the same call. | Distances below 5·10⁻¹¹ print as an exact 0 or 1 does. Native and WebAssembly output are now identical on 45,000 rows. A value that lands within a last-bit difference of a rounding tie can still print differently: 2 of 260,000 rows differed in the sixth decimal of `c_posterior`. |
| **The quick start skipped `prepare`** (`lib.rs`) | `from_config` on raw sets has no anchors and no regions: no junction, no FR/CDR, silently. | The example prepares the database. The crate map listed a nonexistent `refine` module and missed ten real ones; it is rewritten. |
| **Application logic in the binary** (then `bin/assign.rs`, now `src/bin/swig.rs`) | Another front end (WebAssembly) would have duplicated it. | `run.rs`: `RunOptions` reads options by their command-line names; `PreparedDatabase::from_fasta`, `warnings` and `files` handle the database; `Run` writes AIRR rows over any `BufRead`/`Write`. The binary's output is byte-identical before and after on all four subcommands and eight option sets. |
| **`assign_map` always spawned threads** (`batch.rs`) | It panics where threads do not exist (WebAssembly), and it spawned a pointless thread for `threads = 1`. | The caller is one of the workers. A target that cannot spawn runs every chunk on the caller. |

## 3. Open findings, by priority

### 3.1 Behaviour: needs a decision or a measurement

1. **The D stage keeps 32 traces before admission is applied** (confirmed
   by reading; impact measured on a random proxy).
   - `SegmentStage` keeps the top 32 D alignments by local score, and the
     exact-run gate (`Admission`) runs only afterwards
     (`pipeline::stage_d`, `retain_from`). Inadmissible hits therefore take
     traceback slots.
   - On random D-like windows at the posterior-mode floor, the cap bound
     in 85 % of windows and dropped an admissible allele in 74 %.
   - This contradicts "posteriors keep every admissible D" (README,
     ALGORITHM §5) and biases `d_presence_posterior` downward.
   - Fix: make admission the stage's acceptance predicate (§4.1), then
     measure on 20,000 reads per preset.
2. **The residual adjudicator path** (confirmed).
   - `--residual` turns posteriors off, but the D stage keeps sub-floor Ds
     unless `--calls fixed` is also given, because the floor is decided
     from the DAG's flags (`config::from_config`).
   - The residual adjudicator then calls a D whose clipped score is
     negative (reproduced: score −8 on a constructed read).
   - Its `path_score` follows a different objective (pair bonus 24,
     `S_D/2`, no D cost) from the one documented for `CallPath.score`.
   - Fix: demote it to a test baseline, or fix all three.
3. **A `prepare` → reload round trip launders uncertain annotation**
   (confirmed).
   - The written FASTA carries `annotation=SOURCE` but no `uncertain=`,
     and neither is parsed back. On reload every value is a trusted
     header: never flagged, never warned about, eligible as a transfer
     donor.
   - Reproduced: 136 uncertain values became 0.
   - Fix: write `uncertain=` into the headers, parse both keys, and treat
     derived values as given but untrusted.
4. **The synthetic database is mostly pseudogenes** (confirmed).
   - `sim::synthetic` mutates V and J without codon checks: 88 % of V
     alleles and 9 % of J alleles carry in-frame stops (50 seeds).
   - Consequences:
     - `productive` is never exercised meaningfully;
     - `simulate::tests::constant_regions_continue_the_reading_frame…`
       passes only for its seed;
     - the standalone project's browser demo, which used it (since
       removed), showed about 99 % unproductive reads (it said so).
   - Fix: substitute synonymous-or-sense codons, as C already does. This
     shifts every synthetic accuracy floor, so re-measure them with it.
5. **`evaluate` scores a calls file without `rev_comp` as misoriented**
   (confirmed).
   - A missing column parses as forward, and allele metrics are gated on
     orientation. A perfect caller without the column scores 44 % V
     equivalence at `--rc 0.5`.
   - Fix: require the column, or score allele identity regardless of
     orientation.
6. **`lossless_lmer` hard-codes 4 for any distributed-D gate** (confirmed
   counter-example at 75 % identity).
   - It is lossless only for the default gate parameters.
   - Fix: derive L from the gate's fields.
7. **An inverted D window panics** when `ExpandedResidual { min_width: 0 }`
   is configured (confirmed; the default 4 is safe).
   - Fix: clamp `start..end.max(start)` and drop empty windows.
8. **Public fields bypass validation** (confirmed).
   - `GapCosts { open: 10, extend: 10 }` built literally saturates every
     SIMD lane (32767 against the scalar 77,990).
   - `Allele`'s public fields admit inconsistent regions.
   - `DistributedEvidence` fields break finding 6.
   - Fix: private fields with validated constructors; lane admission
     should also reject positive gap costs.
9. **Retrieval caps candidates at 64 alleles silently**
   (`seed_vote::retrieve`), against "every allele of a kept gene becomes a
   candidate". Either document it or cap per gene.
10. **Max-marginal D decoding** (suspected).
    - `CalledSegment::by_marginal` carries no allele id, so a merged
      non-primary allele can win yet be listed second.
    - Its conditioning is on presence only, not on the topology its doc
      states.
11. **FR/CDR regions may come from a different allele than `v_call`'s
    first name** (suspected).
    - `read_regions` projects through the called placement, while the
      report re-sorts co-optimal alleles.
    - They differ only for alleles with indels relative to each other.
12. **Smaller issues**:
    - `frames.rs` breaks a `HashMap` tie by `RandomState` order, so frames
      can differ between processes;
    - the orientation test checks the winner against its null but not the
      margin between strands (reads under ~17 nt);
    - the Karlin–Altschul σ series stops at 500 terms (5 % error in K for
      `5/−2`);
    - traceback memory per lane group is bounded only by
      `16,000² × lanes` (latent; 4 GiB is all of wasm32);
    - required AIRR columns are missing: `sequence_alignment`,
      `germline_alignment`, `locus` (fixed in Swig 0.39.0: Swig's schema);
    - `AirrRecord::tsv_header()` omits `c_` columns that rows can carry
      (fixed in Swig 0.39.0: one 111-column schema for every row);
    - a UTF-8 BOM, non-UTF-8 headers and a reader used again after an
      error are mishandled;
    - `examples/benchmark --tail` scores against the untruncated truth;
    - `tests/annotation.rs` has a wall-clock bound (< 200 ms) that failed
      once under load;
    - several accuracy floors sit within 1–3 reads of the measured value.

### 3.2 Platform notes from the port

- The library's IO was already almost platform-free: germlines from
  `&str`, reads from any `BufRead`, AIRR to `String`. Only `batch` (threads)
  and `PreparedDatabase::write` (files) touched the OS.
- Node's WASI host segfaults intermittently on one test binary
  (`tests/simulated.rs`, `family_bounds_agree…`). The same binary passes
  under Wasmtime every time, and V8 itself ran 45,000 rows of the browser
  build cleanly, so the fault is Node's WASI layer. The standalone
  project's CI therefore ran the WebAssembly suite under Wasmtime.
- `sim::random`'s doc claims "a seed fixes every simulated read on every
  platform", but `ln`, `exp` and `cos` come from each platform's libm. On
  `wasm32-wasip1`, 105 of 20,000 hard-preset truth rows differed in the
  last bit of `shm_rate` (read sequences were identical in that sample).

## 4. Design recommendations, by value

1. **One acceptance mechanism.**
   - Make `SegmentStage` generic over an acceptance predicate and give the
     D stage its `Admission` directly.
   - This fixes §3.1-1, removes `retain_from` and the magic
     `min_identity: 0.60`, and lets inadmissible hits free their slots.
2. **Adjudication modes as a type.**
   - `PlacementDag`'s `evidence`, `calls`, `posterior` and `d_allele`
     interact, and invalid combinations degrade silently.
   - `PlacementDag::new` defaults to the fixed scale while
     `CallModel::default()` is per-read.
   - Replace them with a validated mode type and private fields, and build
     the D stage's floor from the adjudicator actually used
     (`Adjudicator::d_score_floor`). This fixes §3.1-2 at its root.
3. **Write the lattice recursion once.**
   - Viterbi, forward and backward repeat the start group, the
     inside-overlap entry walk, the beyond-end tail and the per-column cut
     valuation.
   - A `Lattice::edge(b, e)` yielding entry and predecessor offsets once,
     plus a column walker parametrised by valuation, would remove the
     triplicated index arithmetic. A fully generic semiring is not worth
     it, because of the back-pointers, tie rules and scaling.
4. **Split the two monoliths.**
   - `annotate()`: one function per pass over a state struct. This also
     removes its `too_many_arguments` allowance and the duplicated J/V
     rescue.
   - `segment.rs`:

     | Module | Contents |
     |---|---|
     | `prepared.rs` | `PreparedSegment`; this breaks the segment ↔ family_bound and segment → retrieve → segment cycles |
     | `segment/policy.rs` | keep, acceptance and tie order |
     | `segment/rank.rs` | ranking |
     | `segment/bounded.rs` | ranking with sibling bounds, with `family_bound` |
     | `segment/waves.rs` | traceback waves |

     `finish` and the outer extension move into `placement`/`flank`, and
     one `run_jobs` helper replaces the four job-building copies.
5. **Remove what was measured and rejected, or move it to tests.**
   - `TieOrder::GeneLeadersFirst`: one tie key for every stage and request
     size.
   - The residual adjudicator.
   - `SeedSampling::{Open, Minimizer}`.
   - `UngappedPrefilter` in the default type.
   - `ResidualWindow`, which equals `ExpandedResidual { min_width: 1 }`.
   - Unused API: `Evaluation::merge` and its siblings, `Backend::lanes`,
     `KarlinAltschul::bit_score`, `alphabet::is_canonical`.
6. **Narrow the public surface.**
   - `Within`, `ReadModel`, `Topology` and `PlacementDag::{viterbi,
     conditioned}` are public only for `tests/synthetic.rs`'s brute-force
     oracle. Move the oracle into unit tests.
   - Make the `annotate` internals and `family_bound` helpers
     `pub(crate)`.
   - `HotspotScoring`, test-only by its own doc, should leave the prelude.
7. **Structured notes.** Annotation notes are matched by text
   (`starts_with("Cys23 ambiguous")`), and tests assert on their wording.
   Use an enum.
8. **A shared test kit.**
   - Setting up the synthetic assigner is repeated about ten times.
   - The simulate → assign → evaluate loop is in a test and in an example.
   - Accuracy floors are fixed numbers tied to one RNG stream; computing
     them as three standard errors below the measured rate would make them
     robust.
9. **Documentation hygiene.** Move ALGORITHM.md's review narrative to
   ASSESSMENT.md. Restate `posterior.rs`'s header in terms of the
   per-read model (it still says `exp(λ·score)`).

## 5. What this branch adds for WebAssembly

The standalone project shipped its own JavaScript API, Web Worker and
browser page with the module. Since the engine moved into Swig, those are
removed: Swig's wrapper and workers use the module instead. The table says
where each piece lives now.

| Piece | Where | Notes |
|---|---|---|
| `simd128` kernel backend | `align/lanes.rs`, `align/inter_seq.rs` | `Lanes<16>` over two `v128`. It is compiled only into `simd128` builds, since WebAssembly has no run-time feature detection. The whole suite passed on `wasm32-wasip1` + `simd128` under Wasmtime; the kernel tests compare `simd128`, portable and scalar. |
| Thread-optional batches | `batch.rs` | The caller works; helpers are spawned only where possible. |
| Shared run | `run.rs`, `prepare.rs` | `RunOptions`, `PreparedDatabase::{from_fasta, warnings, files}`, `Run`. Used by the binary and the module alike. |
| Module | `wasm/` (crate `swig_wasm`) | Zero-dependency `cdylib` with a C ABI over linear memory and no imports. The review built a SIMD and a scalar module (0.9 MB each, 0.3 MB gzipped); Swig builds only the SIMD one, `public/swig.wasm`, with `scripts/build-swig-wasm.sh`. |
| JavaScript API | Swig's `src/aligner-engine.ts` | One engine per worker over the ABI, used by Swig Web's aligner workers and swig-cli's workers. It replaces the standalone project's API. |
| Page | Swig Web | Replaces the standalone project's self-contained page (files or demo data in; AIRR TSV, prepared database and a preview out). |
| Checks | `wasm/src/lib.rs` (`the_abi_runs_the_library_run`), Swig's `tests/swig-wasm.test.mjs` | At review time, a Node.js check and the standalone CI showed that both modules agree, sessions reassemble, errors surface, and the module output equals the native command line: 45,000 rows across presets and options were identical, and so was headless Chromium with four workers. |

Throughput under V8 with SIMD was 1,100–1,450 reads/s per worker on
simulated reads. That is about 60 % of the native AVX2 kernel on the same
file, single-threaded.
