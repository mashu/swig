# CHMMAIRRa PCR-chimera inference

## Input projection

CHMMAIRRa runs on V (default) or J for the current cumulative working set. Its reference input must be an equal-length aligned FASTA with at least two records. Swig can build this MSA from the exact assignment references using Kalign WASM or validate a user-supplied MSA.

For each AIRR row, Swig locates the assigned allele (or gene-level fallback) in the MSA. The local AIRR germline alignment is located in the degapped MSA reference and first extended gaplessly toward the junction — a V to its germline 3′ end, a J to its germline 5′ start — by the best-scoring prefix (+1 match, −1 mismatch) of the read bases beside it (`sequence`, oriented by `rev_comp`, after `v_sequence_end` / before `j_sequence_start`). The aligner may end a V (start a J) before trailing mismatches it attributes to the junction, and a chimera's last informative columns can lie there; the extension keeps them regardless of where the aligner put the boundary and never takes a flank that does not score positive. Query insertions relative to that germline are then removed, and the observed query is threaded onto the allele's MSA columns. Outside observed coverage, gaps and `N` are non-informative emissions. A missing reference/alignment is reported as unevaluated rather than assigned an arbitrary posterior.

Distance from reference (DFR) exactly counts unequal characters in the (extended) local sequence/germline alignment up to their common length, including gap and `N` differences, following CHMMAIRRa's `add_DFR_column!`. Rows below the chosen DFR are not evaluated; the default minimum is 1.

## HMM

At each MSA column the hidden reference is one aligned allele. Under the discretized model, state also contains one of the configured mutation-rate bins. An informative observed nucleotide equals the state reference with probability \(1-\mu\) and each alternative nucleotide with probability \(\mu/3\). Gap/`N` observations contribute emission 1.

The transition probability of switching reference is the configured chimera prior divided by alignment length; transition to a different reference enters the `chimeric` state, which remains chimeric thereafter. Mutation-rate switching is modeled separately. A scaled forward pass returns posterior probability that the final path has switched reference at least once.

- **Baum–Welch (`BW`, IG default):** estimates one mutation probability per reference with smoothed expected match/mismatch counts, then evaluates the reference-switch posterior.
- **Discretized Bayesian (`DB`, TCR default):** marginalizes over the displayed mutation-rate grid. TCR starts with a single 0.005 rate; IG's default grid spans 0–0.25.

The run stores probabilities for all evaluated ordinals and only a bounded heap of the top 500 flagged records for interactive display. The threshold (0.95 default) is a decision layer; changing the later exclusion threshold does not rerun the HMM.

## Detailed reconstruction

Opening one flagged record reruns a full Viterbi traceback on demand. It reports the starting reference, each reference-switch boundary, complete reference path, and displayed parent sequences. The identity-highlighter colors query agreement with the displayed parents; color is not itself the hidden state. Triple-gap columns may be hidden from the SVG only, never from the HMM.

## Literature and code relationship

This is a **browser port of the CHMMAIRRa/CHMMera structured HMM**, including reference/rate states, absorbing chimera indicator, DFR convention, BW/DB modes, and posterior definition. Swig differs in its AIRR-to-MSA threading, browser worker orchestration, bounded top-record retention, optional J mode, and explicit cumulative-mask/filter UI. It is not claimed to be byte-for-byte identical to the Julia program. See Chernyshev et al., [Detection of PCR chimeras in adaptive immune receptor repertoire data](https://academic.oup.com/bioinformatics/article/41/11/btaf576/8297098), and the [CHMMAIRRa.jl source](https://github.com/MurrellGroup/CHMMAIRRa.jl).

## Limitations

- The model explains template switching among references present in the MSA; missing parental alleles can distort inference.
- A high posterior is model evidence, not experimental confirmation of a PCR chimera.
- Mutation emissions are context independent.
- Filtering may retain unevaluated rows or exclude them, by explicit user choice.
- D segments are not modeled.
