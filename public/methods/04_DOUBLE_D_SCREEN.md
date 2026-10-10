# Optional Double-D / VDDJ evidence screen

## Scope

The screen runs inside the aligner (`engine/src/double_d.rs`) after the ordinary V(D)J call of a read. It applies only when that call places both V and J, and only on loci with D segments (IGH, TRB and TRD); a read whose locus is unknown is screened against every D allele. The screen is off by default. A positive result is written to a sparse sidecar table and shown as an interactive overlay; the ordinary AIRR row is identical with the screen on or off. This separation is intentional, because tandem-D inference from short, mutable junctions is substantially less certain than ordinary V/J assignment. The ordinary call itself allows at most one D.

## Candidate construction

The search window is the oriented read between the called V's end and the called J's start, as reported in that read's AIRR row (`v_sequence_end` to `j_sequence_start`, both excluded). Its length is the V–J span. In `long span` mode the read is skipped unless the span reaches the configured threshold (40 nt by default); `all` mode screens every eligible read. A window shorter than twice the seed length \(k\) cannot hold two seeds and is skipped.

Every \(k\)-mer without `N` of each D allele at least \(k\) long is indexed; \(k = 11\) by default and can be set from 6 to 24. Only alleles whose locus pairs with the read's locus are used. An exact seed hit is extended left and right without gaps, with match `+2` and mismatch `−2`, keeping the highest-scoring extension in each direction. A hit's score is that of the extended ungapped alignment.

Shorter rescue words of length \(\max(6, k-3)\) (8 for the default \(k = 11\)) are always indexed as well. A hit found only through a rescue word is discarded unless its ungapped extension reaches at least \(k\) aligned read bases. The shorter word therefore recovers seeds disrupted by a mutation or sequencing error without lowering the minimum span of evaluated evidence. When \(k = 6\) the rescue word is not shorter, and only exact \(k\)-mers are used. Exact seeds are taken first, and each allele diagonal contributes one hit, whichever word found it.

Two hits form a candidate pair \((D_1, D_2)\) only when they are ordered and non-overlapping on the read. The pair score is the sum of their two extension scores. It must exceed the best single hit's score by at least the configured minimum score gain (8 by default).

## Pseudo-tandem veto

Two seeds from one D allele with a mutation or deletion between them can mimic two D segments. To reject that explanation, Swig takes the complete read span from the start of \(D_1\) through the end of \(D_2\). It removes a configurable delta length (5 nt by default), at every split between the two ends, and finds the smallest Hamming distance between the remaining substring and any equally long substring of a locus-compatible D allele. If that distance is at most the configured maximum (3 mismatches by default), the pair is vetoed. When the shortened span is shorter than \(k\), or longer than every D allele, no single D can be compared; the pair is then not vetoed and its pseudo distance is reported as empty.

Pairs are tested in evidence order: pair score, score gain, matched bases, shorter NP2, then allele names. The first pair that survives the veto is reported. Surviving pairs within two score units of it, up to 12, are listed as alternatives; lower-scoring combinations are not tested. Each D call lists every allele aligned with the same score and coordinates.

## Output

The sidecar has 38 columns, from `swig_batch_record_index` through `swig_double_d_alternatives`. `swig_batch_record_index` is the read's 0-based index within its assignment batch; Swig uses it to attach the row to its AIRR record. The other columns report:

- the read's `sequence_id`, the locus, and the V, J and ordinary (`standard_d_call`) calls;
- D1 (`d_*`) and D2 (`d2_*`) calls, scores, identities, CIGARs, read and germline coordinates, and aligned sequences;
- the NP2 (D1–D2) and NP3 (D2–J) sequences and lengths; and
- the V–J span, seed length, pair score, best single score, score gain, pseudo distance and mode.

Alternatives are written as `D1+D2|score|gain|pseudo distance` (−1 when no distance applies), `;`-separated.

## Literature relationship

The 11-nt exact seed starting value and the single-D mimic concern are informed by IgScout/tandem-D work, but this is a **custom Swig evidence screen**, not an implementation of IgScout. See Safonova and Pevzner, [V(DD)J recombination is an important and evolutionarily conserved mechanism](https://pmc.ncbi.nlm.nih.gov/articles/PMC7605257/), and Briney et al., [Frequency and genetic characterization of V(DD)J recombinants](https://pmc.ncbi.nlm.nih.gov/articles/PMC3449247/).

## Explorer and downstream use

The Double-D explorer filters the sparse positive sidecar; opening a row overlays the tandem interpretation on the original AIRR record. Repertoire selection may explicitly require or exclude positives. Where both D alignments can be safely projected, lineage-root construction can use V–D1–D2–J rather than allowing the unchanged baseline single-D composite to erase D2.

## Limitations

- A positive is supported VDDJ evidence, not a definitive genomic recombination call.
- The window depends on the ordinary call's V and J ends: a V or J alignment that reaches into a D removes those bases from the screen.
- Exact seeds lose sensitivity under mutation or sequencing error. The rescue words recover some of it, at the cost of more candidate hits; shorter configured seeds reduce specificity.
- The ungapped seed extension is deliberately fast and is not a full pair-HMM.
- The pseudo-tandem veto is a bounded alternative-explanation test, not a full probability ratio over recombination histories.
- Defaults are starting values and should be checked against negative controls and organism-specific D references.
