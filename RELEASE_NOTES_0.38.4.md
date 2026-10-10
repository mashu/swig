# Swig 0.38.4

## Constant-region calling

- C is now a strictly post-J annotation: SwiftIG searches from the documented one-base J/C overlap instead of aligning long C references against the complete V(D)J read.
- A C hypothesis can no longer clip or otherwise move the selected J boundary.
- Reported C calls independently require at least 90% local identity and calibrated chance evidence. Ordinary/unanchored hits retain the ≥30-nt and `c_support <= 1e-5` route; short gapless 5′-anchored hits instead use a database-size/offset-adjusted `c_support <= 1e-4` gate, so primer-limited C prefixes are not rejected by a universal 30-nt rule.
- A separate **covered C-prefix identity** threshold is available in Web advanced settings, pipeline JSON (`annotation.minimumConstantPrefixIdentity`), and assignment-only CLI (`--minimum-c-prefix-identity`). It is off by default. When enabled it includes leading mismatches and indels that local alignment could omit, reports `c_prefix_identity`, and removes only the C assignment on failure.
- C evidence no longer contributes to V/D/J boundary, call, strand, or orientation selection; even a rejected C hypothesis cannot change the rearrangement assignment.
- Local C alignment may end before a mismatching primer, adapter, or downstream tail; a good contiguous match is not discarded merely because later sequence disagrees.
- Ambiguous IgG- or IgA-subclass labels are shown conservatively as `IgG` or `IgA` in the web interface rather than taking the first comma-separated subclass.

On the first 1,000 constant-negative records of each supplied simulation, false C calls changed from 29 to 0 (IgM-like) and 40 to 0 (IgG-like). The final archive's adaptive short-prefix implementation also produced zero false C calls among 5,000 records of either dataset. Those runs used 24.712 and 31.282 CPU-seconds, respectively. C-enabled CPU time improved 6.60-fold and 4.77-fold on the original 1,000-record comparisons because upstream V(D)J sequence is no longer included in C alignment.

Without a C database, complete 20,000-record IgM-like and IgG-like fair V/D/J accuracies remain unchanged; final-archive CPU changed by -0.72% and +0.18% versus the paired 0.38.3 runs. See [`BENCHMARK_CONSTANT_CALLING_0.38.4.md`](https://github.com/mashu/swig/blob/a08a2c3fb2d42dba9dd1ac97e6fda783089b5b3f/BENCHMARK_CONSTANT_CALLING_0.38.4.md) for the protocol, accuracy table, false-positive controls, and timings.

## AIRR TSV and gzip output

- `swig-cli --vdj -out name.tsv.gz` now writes actual gzip-compressed AIRR TSV instead of plain TSV with a misleading suffix.
- Pipeline runs accept `output.airrCompression: "gzip"` or `--airr-compression gzip`; annotated and processed AIRR tables are compressed incrementally and their resolved names are reported. The lazy lineage-study AIRR remains plain because the manifest uses seekable byte offsets.
- Web Results and post-analysis exports offer None/Gzip for TSV, CSV, and JSONL. Save-As output uses a bounded native compression stream; the IndexedDB download worker can stream `.tsv.gz` directly without assembling a repertoire-sized Blob.
- User-controlled AIRR identifier/quality fields are delimiter-safe: embedded tab, carriage-return, or newline characters are replaced with spaces by both retained output implementations.
- CLI regression coverage decompresses the result, checks every row against the header field count, and verifies that `cut -f13` selects `j_call`.

Plain `.tsv` output and stdout remain uncompressed by default. Comma-separated values inside a call field are deliberate co-optimal calls, not delimiters; `grep -v ','` continues to exclude those ambiguous calls as expected.

## Interface and compatibility

The browser and CLI version strings are 0.38.4. Both use the same rebuilt WebAssembly engine. The web advanced controls, saved sessions, and exported CLI configs retain the optional C-prefix gate; CLI help documents the adaptive C evidence and compression options; and the full method document describes post-J placement, output behavior, and isoform limitations.
