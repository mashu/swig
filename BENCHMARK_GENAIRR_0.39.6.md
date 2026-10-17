# GenAIRR benchmark: swig 0.39.6, IgBLAST 1.22 and MiXCR 3.0.13

Reads simulated by GenAIRR 2.3.1, the three tools run against the same
GenAIRR reference, and every V, D and J call scored against GenAIRR's
ambiguity-aware truth. The report also measures where swig loses accuracy
and how much two changes, prototyped as post-processing of swig's output,
would gain. All reads are simulated; none of this is external biological
validation.

## Summary

- **swig 0.39.6 already leads or ties on V everywhere, and beats IgBLAST
  on J for every IG locus.** MiXCR's top hit is ahead on IGH J at low SHM,
  because it does not list J4\*02 and J5\*02 together. On IGH D swig
  leads under no SHM, 10 % SHM and the realistic set. It is level with
  IgBLAST at 2 % and 5 % SHM (−0.17 and −0.08 points).
- **Most of the loss is D**: 2–13 points of D against at most 0.6 of V or
  J. About 90 % of swig's wrong-gene D calls overlap the true D. Swig finds
  the D but picks the wrong allele or frame there, by local alignment
  score rather than by what the read supports over the whole D. Swig
  also declines to call a D on 0.3–1 % of reads, though the truth always
  has one.
- **J loses over-reported ties.** When J4\*02 and J5\*02 are trimmed to
  their shared 3′ part, they align identically. GenAIRR's truth
  nonetheless lists only the sampled one, because it compares alleles at
  the same germline index (see Scoring).
- **Two prototypes beat both tools everywhere on IGH D and IG J:**
  - a D call from a posterior over allele, offset and both trims,
    reporting the allele set of the highest expected score;
  - J tie sets recomputed position by position at swig's placement.

  With generic geometric trim and N-length priors, IGH D rises by 0.3–4.9
  points over swig 0.39.6 and 0.7–7.7 over IgBLAST; IGH J rises by
  0.2–0.4. GenAIRR's own trimming tables do no better.
- **Reporting genuine ties pays, guessing does not.** MiXCR's top hit
  alone scores 0.2–8.4 points lower on V, and 7.6–8.2 lower on IGH D,
  than all of its top-tied hits.
  Swig's current co-optimal names are the right idea; the gains come from
  getting the tie set right, not from listing more alleles.

## Setup

- **Simulator.** GenAIRR 2.3.1 (`pip install GenAIRR`), Rust kernel, with
  its built-in references: human IGH (OGRDB, 198 V, 33 D, 7 J), IGK (168
  V, 8 J), IGL (181 V, 10 J) and TRB (98 V, 3 D, 16 J; with
  `allow_curatable_refdata()`). 10,000 reads per set, seed 7, forward
  strand. Scenarios (`benchmarks/genairr/genairr_bench.py`):

  | Scenario | Steps after `recombine()` |
  |---|---|
  | naive | none |
  | shm2, shm5, shm10 | `mutate(model="s5f", rate=0.02 / 0.05 / 0.10)` |
  | realistic | S5F 5 %, 5′ primer trim 0–15, 3′ trim 0–5, 0–1 polymerase indel, PCR errors 0.2 %, sequencing errors 0.2 %, 0–2 N calls |
  | tcr_realistic | the realistic artefacts without SHM (TRB) |

- **Reference for the tools.** The GenAIRR reference exported to FASTA:
  - **swig** gets each allele's anchor and region starts in its header;
    an anchor is kept only on its conserved codon, since some TRB ORFs
    carry placeholders.
  - **IgBLAST** gets `makeblastdb` databases and a J auxiliary file.
  - **MiXCR** gets a repseqio library with per-allele anchor points.
- **Tools** (4 threads each, defaults otherwise):
  - **swig 0.39.6**: native binary, R-optimized profile.
  - **IgBLAST 1.22.0** (bioconda): `-organism human -outfmt 19`, `-ig_seqtype Ig` or `TCR`.
  - **MiXCR 3.0.13**: built from source, the last release that runs
    without a license key. `align`, then `exportAlignments` of all hits
    with scores, reported two ways: the top hit, and all hits tied with
    the top score. MiXCR 4.x needs a license and was not run.
- **Run time** for 10,000 reads on 4 threads: swig 0.9–1.8 s, IgBLAST
  84–294 s, MiXCR 6.6–24 s (JVM start included).

## Scoring

GenAIRR's `v_call`, `d_call` and `j_call` are tie sets. Each allele earns a
point for every read base that matches it at the *same germline index* as
the sampled allele, over the segment's true read region. Adjacent N bases
count only while they narrow the tie. Every allele at the maximum is in the
set. Mutations can make another allele win outright: 2 % of reads at 5 %
SHM do not have the sampled V in their truth set, 3 % not the sampled D.

A call C against truth T scores `|T ∩ C| / |T ∪ C|` per read and segment.
Both empty scores 1, one empty scores 0. With a single shared allele this
is `1 / |union|`: one right allele against a two-allele tie scores ½, and so
does a two-allele call against a one-allele truth. The truth always holds a
D for IGH and TRB, so a missing D call scores 0. The tables give the mean
over reads, in %.

Index-based scoring has a consequence. Alleles of different lengths are
compared off by their length difference, so two alleles whose trimmed
remnants are identical are not tied if their germline coordinates differ.
J4\*02 and J5\*02 (48 and 51 bases) are the common case: from the read
alone they are one allele, and the truth names only the sampled one.

## Results

Bold: best of the row group. "+ D posterior" is swig 0.39.6 with both
prototypes applied: positional J tie sets and the D posterior. V is swig's
own call.

| Run | Tool | V | D | J |
|---|---|---:|---:|---:|
| IGH naive | IgBLAST 1.22 | 99.86 | 97.80 | 99.48 |
| | MiXCR 3.0.13 | 99.58 | 86.42 | 99.74 |
| | MiXCR, all top ties | 99.76 | 94.56 | 99.53 |
| | swig 0.39.6 | **99.98** | 98.17 | 99.52 |
| | + D posterior, geometric priors | **99.98** | **98.51** | **99.76** |
| IGH 2 % SHM | IgBLAST 1.22 | 99.80 | 96.55 | 99.36 |
| | MiXCR 3.0.13 | 98.55 | 85.14 | 99.59 |
| | MiXCR, all top ties | 99.65 | 93.34 | 99.37 |
| | swig 0.39.6 | **99.85** | 96.38 | 99.39 |
| | + D posterior, geometric priors | **99.85** | **97.42** | **99.70** |
| IGH 5 % SHM | IgBLAST 1.22 | 99.64 | 93.24 | 99.14 |
| | MiXCR 3.0.13 | 97.23 | 82.03 | 99.05 |
| | MiXCR, all top ties | 99.30 | 90.22 | 98.81 |
| | swig 0.39.6 | **99.68** | 93.16 | 99.15 |
| | + D posterior, geometric priors | **99.68** | **95.67** | **99.47** |
| | + D posterior, GenAIRR's priors | **99.68** | 95.27 | **99.47** |
| IGH 10 % SHM | IgBLAST 1.22 | 99.36 | 84.28 | 98.48 |
| | MiXCR 3.0.13 | 94.63 | 73.49 | 96.77 |
| | MiXCR, all top ties | 98.45 | 81.07 | 96.60 |
| | swig 0.39.6 | **99.42** | 87.12 | 98.70 |
| | + D posterior, geometric priors | **99.42** | **92.02** | **99.06** |
| | + D posterior, GenAIRR's priors | **99.42** | 90.72 | **99.06** |
| IGH realistic | IgBLAST 1.22 | 99.53 | 92.56 | 99.05 |
| | MiXCR 3.0.13 | 95.02 | 81.36 | 98.98 |
| | MiXCR, all top ties | 99.01 | 89.55 | 98.73 |
| | swig 0.39.6 | **99.61** | 92.78 | 99.08 |
| | + D posterior, geometric priors | **99.61** | **95.24** | **99.45** |
| | + D posterior, GenAIRR's priors | **99.61** | 94.82 | **99.45** |
| IGK 5 % SHM | IgBLAST 1.22 | 99.63 | – | 98.88 |
| | MiXCR, all top ties | 99.14 | – | 98.35 |
| | swig 0.39.6 (+ positional J) | **99.74** | – | **99.43** |
| IGK realistic | IgBLAST 1.22 | 99.48 | – | 98.11 |
| | MiXCR, all top ties | 97.99 | – | 97.09 |
| | swig 0.39.6 / + positional J | **99.63** | – | 99.02 / **99.03** |
| IGL 5 % SHM | IgBLAST 1.22 | 99.53 | – | 98.32 |
| | MiXCR, all top ties | 98.48 | – | 98.16 |
| | swig 0.39.6 (+ positional J) | **99.63** | – | **99.10** |
| IGL realistic | IgBLAST 1.22 | 99.34 | – | 98.32 |
| | MiXCR, all top ties | 97.16 | – | 96.49 |
| | swig 0.39.6 / + positional J | **99.51** | – | 98.49 / **98.50** |
| TRB naive | IgBLAST 1.22 | **99.97** | 98.17 | 99.93 |
| | MiXCR, all top ties | 99.94 | 94.10 | **100.00** |
| | swig 0.39.6 | **99.97** | 99.12 | 99.89 |
| | + D posterior, geometric priors | **99.97** | **99.15** | 99.89 |
| TRB realistic | IgBLAST 1.22 | **99.94** | 97.65 | **99.92** |
| | MiXCR, all top ties | 99.80 | 93.66 | **99.92** |
| | swig 0.39.6 | **99.94** | 98.63 | 99.87 |
| | + D posterior, geometric priors | **99.94** | **98.67** | 99.88 |

Light chains: MiXCR's top hit alone scores V 89.6–93.6, J 87.9–97.2.
`benchmarks/genairr/summary.py` prints every row.

### Other metrics, IGH 5 % SHM

The same calls, scored three ways:
- **Jaccard**: the rule above.
- **first**: the first reported allele is in the truth set.
- **any**: the call shares an allele with the truth set.

| Tool | V Jaccard / first / any | D Jaccard / first / any | J Jaccard / first / any |
|---|---|---|---|
| swig 0.39.6 | 99.68 / 99.96 / 99.96 | 93.16 / 94.08 / 94.53 | 99.15 / 99.47 / 99.68 |
| IgBLAST 1.22 | 99.64 / 99.89 / 99.89 | 93.24 / 93.36 / 95.03 | 99.14 / 99.44 / 99.69 |
| MiXCR, top hit | 97.23 / 99.66 / 99.66 | 82.03 / 90.62 / 90.62 | 99.05 / 99.06 / 99.06 |
| MiXCR, all top ties | 99.30 / 99.66 / 99.75 | 90.22 / 90.62 / 92.18 | 98.81 / 99.06 / 99.25 |

## Where swig loses accuracy

`benchmarks/genairr/loss.py` splits each tool's loss (`1 − score`, in
points of the segment's accuracy) by how the call relates to the truth
set.

| IGH 5 % SHM, class | swig 0.39.6 | IgBLAST | swig + prototypes |
|---|---:|---:|---:|
| **D** wrong gene | 4.73 | 4.86 | 3.05 |
| D over-reported tie (superset) | 0.90 | 1.32 | 0.85 |
| D missing call | 0.73 | 0.09 | 0 |
| D under-reported tie (subset) | 0.42 | 0.37 | 0.40 |
| **J** over-reported tie | 0.53 | 0.55 | 0 |
| J wrong gene | 0.32 | 0.30 | 0.53 |
| **V** under-reported tie | 0.21 | 0.24 | 0.21 |
| V over-reported tie | 0.07 | 0.01 | 0.07 |

The realistic and 10 % sets split the same way (D wrong gene 5.0 and 10.0
points for swig 0.39.6).

1. **D allele choice inside the right tract.** 95 of 103 wrong-gene D
   calls in a 2,000-read sample overlap the true D. Swig places the D
   where it is, but chooses among alleles by local alignment score. That
   score clips at the first mismatch and charges the D's mismatch
   penalty. So a mutated true D (`ACTACAATAA`, one substitution) loses to
   another gene's exact 7-mer (`CTACAAT`) within it. The truth counts every
   base of the true D.
2. **Missing D.** Swig's D-presence model rejects 5–10-base Ds and mutated
   short ones. GenAIRR always has a D, so these score 0. IgBLAST nearly
   always calls one, but on these reads its call scores 0 in 82–85 % of
   them.
3. **J ties across germline offsets.** For J4\*02 against J5\*02 on a J
   trimmed by 14–17 bases, swig's two co-optimal placements are equally
   good. GenAIRR's index-based truth names only the sampled one (see
   Scoring), so the pair scores ½. Choosing one by its germline frame
   scores 1 or 0.
4. **V and J tie boundaries.** The bases that separate tied V alleles sit
   at the V's 3′ end, often in the 1–9 bases after a mutated stretch where
   swig's local alignment stops. GenAIRR counts them; swig does not. These
   cost 0.2–0.6 points of V.

## Prototypes

Both run on swig's AIRR output (`benchmarks/genairr/recall.py`,
`d_posterior.py`).

**Positional J tie set.** The tie set is recomputed GenAIRR's way at
swig's placement: every J allele scored base by base at the called
allele's germline index, through the call's CIGAR (insertions skipped, as
GenAIRR does), over swig's aligned span. Ignoring the CIGAR breaks reads
with indels: IGK realistic J then fell from 99.02 to 98.29. The same
re-call for V changes nothing beyond ±0.01. Extending V or J past the
aligned span gave no consistent gain in our tests, whether by GenAIRR's
narrowing walk, counting every base, or a posterior over the segment's
end.

**D posterior.** The read between swig's V end and J start is explained
by hypotheses h = (allele, offset, 5′ trim, 3′ trim) with N bases around.
- **Likelihood:** `(1 − μ)/¼` per matched and `(μ/3)/¼` per mismatched D
  base, where μ is `1 − v_identity` of the read, clipped to [0.005, 0.3].
- **Priors:** geometric trims (means 5 and 7 bases), negative-binomial N
  lengths (mean 8), and uniform genes with each gene's mass split among
  its alleles.
- **Report:** each hypothesis implies GenAIRR's tie set for that
  placement. The reported set maximises the expected score
  `Σ_h P(h | read) · |T_h ∩ R| / |T_h ∪ R|` over the top tie sets and
  their pairwise unions.

What the priors are worth, on 500 reads at 5 % SHM:

| Priors | D |
|---|---:|
| flat (uniform trims and N lengths) | 92.85 |
| negative binomial, shape 4 | 94.67 |
| shape 2 | 95.17 |
| **geometric** (shape 1) | **95.50** |
| GenAIRR's trim, N and usage tables | 95.30 |

The trim prior matters most, the N-length prior next. A peaked trim prior
fails on short D genes: an untrimmed IGHD7-27 (`CTAACTGGGGA`) loses to a
trimmed IGHD1-20 (`GGTATAACTGGAACGAC`) that contains the same `TAACTGG`. No simulator parameter is needed.

**Junction boundaries.** Swig's V end is exact in 75 % of reads and J
start in 62 % (5 % SHM), with J often one base early or 2–5 late. The D
posterior takes both as given. Modelling V tail, N1, D, N2 and J head
together should add D accuracy at high SHM.

## Implementing it in swig

1. **Positional tie sets** in `report.rs` (`called_report`). Replace or
   filter `CoOptimal::names` with the alleles of maximal positional
   matches at the primary allele's germline index along the called
   alignment. That is about 200 alleles × 300 positions per V call, tens
   of microseconds.
2. **D refinement** after the call, for paths with a V and a J: the D
   posterior above with the D database, the read's divergence and the
   geometric priors. The best (allele, offset) over all trims takes
   O(alleles × offsets × length): the 5′ and 3′ trims are independent
   given the D's ends, so running maxima and prefix sums replace the trim
   double loop, about 0.1 ms per read. The top hypotheses' tie sets give
   the reported names. The best hypothesis gives `d_sequence_start/end`,
   the D CIGAR and the NP lengths. `d_presence_posterior` can keep gating
   reads that may have no D.
3. **Next**: a joint junction model (V tail, N1, D, N2, J head); and
   allele-usage priors learned from the repertoire itself, which matter
   when one gene is common (GenAIRR draws genes uniformly).

## Reproduce

`benchmarks/genairr/README.md` gives the setup: GenAIRR and numpy in a
venv, IgBLAST from bioconda, MiXCR 3.0.13 built from source. Then:

```sh
cd benchmarks/genairr
PYTHON=venv/bin/python ./run_matrix.sh 10000 7
PYTHON=venv/bin/python ./improve.sh work/ref_igh work/runs/igh_shm5 parametric:1.0 HUMAN_IGH_OGRDB parametric
python3 merge_calls.py work/runs/igh_shm5.recall.tsv work/runs/igh_shm5.dpost-parametric.tsv work/runs/igh_shm5.final-parametric.tsv
python3 summary.py work/runs
```

`improve.sh` re-calls J positionally and leaves V as swig's (`v=keep`).
For the TRB set, pass `HUMAN_TCRB_IMGT` as the config; it is used only by
the `genairr` priors.
