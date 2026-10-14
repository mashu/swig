# Swig lab-meeting course

`swig-lab-meeting.html` is a self-contained, interactive slide deck that teaches
how the Swig aligner works, for an audience of biologists: alignments and the DP
matrix, log-odds scoring (BLOSUM and Swig's +2/−3 at the read's own divergence),
affine gaps, Needleman–Wunsch / Smith–Waterman / Gotoh, banded alignment, k-mer
seeding and indexing, SIMD, the placement DAG that places V, D and J jointly,
posteriors and E-values, the IgSim-based simulator, and the benchmark results
against IgBLAST and the previous Swig core.

Open the file in any modern browser. No build step and no server; it works
offline (fonts fall back to system faces without network).

Keys: ← → / space to move, `N` speaker notes, `O` overview, `F` fullscreen.
The interactive demos (DP matrix, banding, k-mers, dot plot, D prior, placement
DAG, simulator) run small JavaScript re-implementations of the ideas on toy
germline fragments; every benchmark number is quoted from
`engine/docs/ASSESSMENT.md`, `engine/README.md` and
`BENCHMARK_SWIG_ALIGNER_0.39.2.md`.
