# swig paper

A research-paper draft about the Rust V(D)J aligner in `engine/` (crate
`swig`, version 0.39.7). The manuscript covers only the aligner: its model,
exact adjudication, search, implementation and evaluation.

| File | Content |
|---|---|
| `main.tex` | The manuscript, self-contained: TikZ figures and plot data are inline |
| `references.bib` | Bibliography (natbib, `unsrtnat`) |
| `main.pdf` | The compiled draft |
| `scripts/` | Everything that produces the paper's numbers (below) |

## Build

```sh
cd paper
latexmk -pdf main.tex
```

Needs a TeX Live with `pgfplots`, `siunitx`, `booktabs`, `cleveref`,
`natbib` and `tabularx` (on Ubuntu: `texlive-latex-extra`,
`texlive-science`, `texlive-pictures`, `latexmk`).

## Status

A starting point, not a submission. Every claim was checked against the
code, the regenerated outputs and the references in three adversarial review
rounds (seven, four and two reviewers, each finding re-derived by an
independent verifier; 109 corrections applied). Before submitting:

- authors, affiliations, acknowledgements and funding are placeholders;
- every result is on simulated reads; a real-data evaluation is the main
  missing piece (see the Limitations section);
- partis, IMGT/V-QUEST and MiXCR 4 were not compared;
- re-check the bibliography against the publishers' records (the GenAIRR,
  partis and IMGT/GENE-DB entries were checked; the others are standard
  references entered from their DOIs);
- pick a journal template (the draft uses `article` so that it moves easily).

## Where the numbers come from

All results except the MiXCR columns were produced for this draft with swig
0.39.7 (`engine/`, commit of this branch), IgBLAST 1.22.0 (Bioconda) and
GenAIRR 2.3.1, on a 4-vCPU Intel Xeon VM at 2.1 GHz (AVX-512BW), 15 GB RAM,
Ubuntu 24.04, Rust 1.97.0.

| Paper item | Produced by | Output |
|---|---|---|
| GenAIRR table and figure (`tab:genairr`, `fig:genairr`) | `scripts/genairr.sh` | `work/genairr/scores.txt` |
| MiXCR columns of the GenAIRR table and figure | `BENCHMARK_GENAIRR_0.39.6.md` (MiXCR 3.0.13 built from source; not re-run) | same reads: GenAIRR seed 7; the regenerated reads reproduce that report's IgBLAST and swig scores to every digit |
| Out-of-sample GenAIRR table (`tab:oos`; seed 50000, reads disjoint from seed 7) | `scripts/genairr.sh` | `work/genairr/scores.txt`, section `seed50000` |
| IgSim table, swig vs IgBLAST (`tab:igsim`) | `scripts/igsim.sh` | `work/igsim/results.json`, `igh_*.swig.tsv` / `igh_*.igblast.tsv` |
| Specificity (random and IGK reads) | `scripts/igsim.sh` | `work/igsim/controls.txt` |
| Truncated reads (`tab:truncated`) | `scripts/igsim.sh`, `scripts/truncated.py` | `work/igsim/trunc/summary.txt` |
| Ablations (`tab:ablation`, `fig:operating`) | `scripts/igsim.sh` | `work/igsim/results.json` (`paired`) |
| D errors (`tab:derrors`) and search vs model errors (separate sample, seed 1) | `engine/examples/d_errors.rs`, `diagnose.rs` via `scripts/igsim.sh` | `work/igsim/d_errors_*.txt`, `diagnose_*.txt` |
| Calibration (`tab:calibration`, `fig:calibration`) | `scripts/igsim.sh` | `work/igsim/results.json` (`ece`, `brier`, `bins`) |
| Constant region paragraph | `scripts/igsim.sh` (`ighc_*` sets) | `work/igsim/results.json` (`c_length`) |
| Speed (`tab:speed`) | `scripts/timing.sh` on an idle machine | `work/igsim/timing/summary.txt` |

Design measurements quoted in the methods to justify a constant (for
example the open-flank allowance, the 435 checked guarded calls, the
prototype's 93.12 against 93.16, the 37–47 % cost of a separate D pass, the
0.6 % instruction cost of the D priors) are taken from `engine/docs/ALGORITHM.md`, `BENCHMARK_GENAIRR_0.39.6.md` and
`RELEASE_NOTES_0.39.7.md`, which state them for the current design.

## Reproduce

```sh
bash paper/scripts/setup.sh      # build swig, IgBLAST 1.22 from Bioconda, GenAIRR venv, IMGT FASTA
bash paper/scripts/igsim.sh      # IgSim experiments (about 20 min on 4 cores, mostly IgBLAST)
bash paper/scripts/genairr.sh    # GenAIRR experiments (about 60 min on 4 cores, mostly IgBLAST)
bash paper/scripts/timing.sh     # throughput; run alone on an idle machine
```

`WORK` (default `paper/work`, ignored by git) holds every intermediate file.
The germline for IgSim is the human set of the IMGT/GENE-DB release 202632-7
pack shipped in `public/references/`, exported as plain FASTA; swig annotates
its anchors itself. IgBLAST's auxiliary file lacks IGHJ5\*03 and IGHJ5\*04;
`igsim.sh` adds them with IGHJ5\*01's frame and CDR3 end, which swig's
anchors reproduce.
