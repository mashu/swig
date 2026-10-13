# Input, FASTQ quality control, subsampling, and study structure

## Accepted inputs

Swig accepts FASTA, FASTQ, AIRR rearrangement tables containing a `sequence` column, and gzip-compressed versions of those formats. Parsing is streaming and record aware: wrapped FASTA sequences and wrapped FASTQ sequence/quality blocks are reconstructed and validated before batching. A malformed record stops the run with the record name where available. A record without an identifier (a blank FASTA or FASTQ header, or a blank or missing AIRR `sequence_id`) is named `sequence_N` by its 1-based position among the input's records before batching, so its name does not depend on the number of workers.

Gzip inspection validates member boundaries before analysis. If one uploaded `.gz` contains multiple independently compressed members, Swig presents those members as multiple recoverable input files and asks whether they are technical pieces of one merged sample or separate samples with independently editable metadata. It does not silently discard bytes after the first member or accept trailing junk as sequence data.

Every accepted input becomes a dataset/library with a stable dataset ID and editable sample, donor/subject, cohort, timepoint, and compartment. Files may share a biological sample ID. A directory tree initializes donor labels from its first nested directory; this is only an initializer and remains editable.

## FASTQ filtering

FASTQ QC is disabled by default. FASTA and AIRR rows have no Phred scores and pass this stage unchanged.

For a retained FASTQ prefix of length \(L\), Swig computes the expected number of base-call errors directly as

\[
E=\sum_{i=1}^{L}10^{-Q_i/10}.
\]

The read is retained when \(E\) is no greater than the configured maximum. Phred+33 and Phred+64 are explicit choices; a character outside the chosen encoding range is an error rather than being silently clipped. The shipped threshold is 0.01 but the entire filter is off until enabled.

Optional 3′ trimming uses a terminal running window. While the current terminal window's mean Phred is below the configured threshold, one terminal base is removed and both the window sum and expected-error sum are updated. Trimming stops at the first passing suffix boundary. A trimmed read below the configured minimum retained length is rejected. Defaults are a 4-base window, mean Phred 20, and minimum length 50; trimming is also off by default.

The operation order is:

1. optional 3′ quality trimming;
2. minimum retained-length check;
3. full retained-prefix expected-error test;
4. optional reservoir sampling;
5. V(D)J assignment.

This is a **direct expected-error implementation**, not a USEARCH/VSEARCH wrapper. The expected-error concept follows Edgar and Flyvbjerg, but Swig's streaming parser, user controls, and trimming implementation are its own: [Error filtering, pair assembly and error correction for next-generation sequencing reads](https://doi.org/10.1093/bioinformatics/btv401).

## Exact random subsampling

When enabled, each dataset is sampled independently using seeded one-pass reservoir sampling. If a dataset has \(N\) eligible records and the requested size is \(k<N\), each record has probability \(k/N\) of inclusion. The selected records are restored to original input order before batching. The complete input must still be scanned, but memory is \(O(k)\) records per dataset.

The standalone CLI imports the same parser, quality-filter function, seeded generator, and per-dataset seed derivation used by the browser. Its config records the complete `preprocessing.fastqFilter` and `preprocessing.subsample` blocks. Quality filtering therefore precedes the reservoir in both execution routes; the requested sample size applies to QC-eligible reads, not raw FASTQ records.

This is a standard reservoir-sampling application; the seeded generator and per-dataset policy are Swig engineering choices. See Vitter, [Random sampling with a reservoir](https://doi.org/10.1145/3147.3165).

## Multi-dataset boundaries

Metadata are written on every combined AIRR row. `sequence_id` is prefixed by the stable dataset ID to prevent cross-file collisions, while the original identifier is retained separately. Downstream donor, sample, cohort, and whole-study scopes use these values as hard partition boundaries—not merely labels.

Changing metadata after assignment rewrites the browser-local AIRR metadata overlay and compact indexes in one IndexedDB transaction. Assignment fields do not change, but downstream state that depended on the old boundaries is invalidated. If the transaction is cancelled or fails, it is aborted and the old metadata/index state remains active.

## Complexity and limitations

- Parsing, QC, fingerprinting, and unsampled batching are \(O(B)\) in input bytes.
- Reservoir sampling is \(O(N)\) time and \(O(k)\) memory per dataset.
- Quality trimming is 3′-only. It is not adapter removal, paired-read merging, UMI consensus, or a learned error model.
- The 0.01 expected-error default is deliberately strict and assay dependent.
- A concatenated gzip member is a file container boundary, not biological evidence that members are separate samples; Swig therefore asks.
- A browser-exported CLI config stores each separated member as an end-exclusive compressed-byte `gzipRange` against the original file and retains its editable display name separately. This reproduces the choice without materializing temporary member files.
