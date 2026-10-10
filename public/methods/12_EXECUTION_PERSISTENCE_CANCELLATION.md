# Execution, persistence, exports, and cancellation

**Evidence label:** direct implementation of Swig's browser-local state machine.

This document describes orchestration rather than a biological estimator. It defines stage order, which state is authoritative, what sessions contain, and what “Cancel” guarantees.

## Immutable input and derived state

The assigned AIRR rows are the immutable base table. Downstream operations do not delete or rewrite that table. They create one or more of:

- a call overlay containing policy-selected V/D/J calls;
- a cumulative active-record bit mask;
- representative ordinals and abundance vectors from collapse/denoising;
- chimera, selection, and lineage masks or assignments;
- compact dashboards and diagnostics;
- selected-lineage alignments, trees, and UCA results; and
- provenance tying every derived object to input ordinals, lineage membership, alignment fingerprints, and parameters.

Resetting a stage restores or recomputes derived state; it does not reconstruct the original AIRR file because that table was never destructively changed.

## Interactive and pipeline order

The automatic pipeline uses the same operations and runs them in this fixed order:

1. per-sequence V(D)J assignment;
2. optional repertoire allele fit and policy-selected call reassignment;
3. exact collapse or denoising;
4. optional CHMMAIRRa exclusion;
5. repertoire selection;
6. lineage assignment; and
7. SHM and lineage-aware missing-allele diagnostics.

Allele refinement is deliberately first among cross-record stages because a selected V/D/J call can change collapse partitions, chimera/reference grouping, and lineage call compatibility. Targeted query, individual lineage alignment, ordinary phylogeny, and phylogenetic UCA inference remain on demand because they require a selected target.

Interactive stages consume the current cumulative mask. Applying or changing an upstream stage clears dependent downstream objects whose partitions could be stale. Merely fitting or previewing a model does not apply it.

Each post-analysis card may be skipped independently. A skipped card is omitted from guided next-step navigation and makes no state transition. Step 04, repertoire selection, is skipped by default; this avoids requiring a no-op selection scan in the ordinary path. Include/skip choices are session state and can be changed later.

## Browser storage

Large AIRR text is held in chunked IndexedDB storage, with a compact per-record index for filtering and lookup. Sequence payloads needed by repertoire-scale algorithms are streamed in batches instead of constructing a second full in-memory AIRR table. Workers receive typed arrays, compact records, or bounded sequence arenas appropriate to each operation.

Exact collapse uses deterministic key shards for large inputs. FAD/conservative/indel denoising dynamically schedules complete independent study/C/locus/V/J partitions across nested workers, while CHMMAIRRa schedules bounded row batches. A single denoising partition is intentionally indivisible because its parent/centroid decisions are global within that partition. The CLI uses equivalent worker-thread scheduling. Assignment, post-analysis, and CHMMAIRRa pools all honor the configured worker count, but they run at different pipeline stages rather than concurrently.

Study-metadata re-indexing is one IndexedDB read/write transaction. Cancellation aborts that transaction, so the previous index remains authoritative.

## Live assignment telemetry

The browser assignment card reports three scheduler-level quantities without changing assignment or storage behavior:

- **Queries / s** is a time-aware exponential moving average with a two-second time constant, sampled no more than twice per second. Its numerator is the number of AIRR records acknowledged after the browser has appended their batch to the result store. It therefore measures end-to-end committed throughput, including output/indexing backpressure, rather than only the aligner's kernel speed.
- **Active workers** is the exact number of compute-worker slots currently executing an annotation batch, shown against the effective worker-pool size. During per-dataset initialization it instead indicates that the concurrently launched workers are instantiating the aligner module and preparing their own germline databases (annotation and seed indexes). It is not an operating-system CPU-utilization estimate.
- **Assignment state** distinguishes per-dataset initialization, concurrent input streaming, draining already-dispatched final batches, and result-index finalization. A falling active-worker count during draining is expected; a sustained low count while streaming identifies scheduler or result-store backpressure for follow-up profiling.

The browser creates a fresh compute pool for each dataset. Compiling the shared WebAssembly module and concurrently constructing each worker's germline index can therefore produce a brief all-core burst at every dataset boundary before steady-state query processing begins. The telemetry exposes that burst as initialization rather than counting it as assignment throughput.

The aligner module (`swig.wasm`, compiled from the Rust crate in `engine/` with WebAssembly SIMD128) is instantiated with no imports: it needs no WASI or other host interface. WASI, through a browser shim, is used only by the separate Alivibe MSA module of the lineage-alignment workbench.

## Portable sessions and linked AIRR

A portable session stores references, settings, study metadata, call overlays, masks, representative/lineage vectors, dashboards, manual alignments, trees, UCA results, and other compact derived state. It deliberately does **not** duplicate the main AIRR table. The session instead records the linked AIRR filename, header, record count, byte information, and content fingerprint.

Restoration therefore verifies and indexes the user-selected linked AIRR file before installing the saved analysis state. It is not rerunning V(D)J assignment: it is rebuilding browser-local random-access indexes and then applying saved vectors to their exact ordinal positions. Session restoration uses a temporary result store; cancellation or a mismatch clears that temporary store and leaves the current app session unchanged.

Project-directory checkpoints use the same logical separation while keeping the linked artifacts in the chosen directory. Generated alignments are reproducible and need not be duplicated unless a downstream result depends on their exact bytes; manually edited alignments are retained.

## Cancellation contract

An actively running compute-intensive operation exposes **Cancel**. Cancellation is cooperative for streamed JavaScript loops and immediate for dedicated-worker computations by terminating the worker. A replacement worker is created for the next run.

The user-visible guarantee is **commit on success**: after cancellation completes, committed analysis state is the state that existed before the cancelled action began. The Cancel label can temporarily become **Restoring…** while the prior index and overlays are reapplied.

| Operation class | Cancellation mechanism | State after cancellation |
|---|---|---|
| Input/gzip-member inspection | abort signal between file/member reads; no dataset is committed until inspection completes | previous dataset list and next dataset ID |
| Reference preparation | abort propagated through catalog fetch/preprocessing workers | previous compiled reference set |
| Initial V(D)J assignment | analysis workers are terminated; the in-progress result store is cleared | no partial analysis session is installed |
| Metadata re-indexing | aborts the single IndexedDB transaction | previous metadata index and downstream state |
| Allele fitting, collapse/denoising, chimera, selection, lineages, diagnostics, query, neighbours | active worker/loop abort plus a runtime checkpoint of mask, call overlay, representatives, lineage assignments, and dashboards | complete pre-run checkpoint; any partially mutated index is rebuilt before controls unlock |
| Automatic pipeline | one outer transaction spans all enabled post-assignment stages | no partial pipeline result or working-set change |
| Kalign, FastTree, and phylogenetic UCA | dedicated outer worker is terminated; results are assigned only after success | previous alignment/tree/UCA result remains installed |
| Session restoration | temporary store plus abort signal | previous app state; temporary indexed rows are cleared |
| Portable-session save | abortable snapshot checks and streamed gzip compression; download starts only after complete encoding | no session download is initiated |
| Large AIRR/analysis exports | abort signal between chunks; partial writable is aborted where the File System Access API supports it; fallback Blob is not downloaded | no completed download is initiated by Swig |

Some rollbacks require scanning the AIRR index to reapply the checkpoint. Controls remain locked during that restoration so a second operation cannot observe a half-restored state. Cancellation does not pretend that CPU work already performed never occurred, and the browser or operating system can retain transient cache pages, but no partial scientific result becomes authoritative.

On-demand local-index record searches expose Cancel and retain the last completed page when aborted. Small synchronous UI transforms that complete atomically do not expose a running state. Once a browser download has already been handed to the browser's download manager, Swig cannot retract that completed hand-off.

## Failure and navigation behavior

An ordinary computational error follows the same transaction boundary as cancellation: uncommitted derived state is discarded or rolled back, while the error is reported. A tab close, reload, browser Back/Forward gesture, or link that would navigate away while loaded data/results or a run are present triggers Swig's unsaved-work guard. The user can stay on the page or explicitly leave. Browser security permits only the browser's standard confirmation text for a true tab close/reload.

This guard protects accidental navigation; it is not durable persistence. Save a portable session or project checkpoint before deliberately closing the page.

## Exports

AIRR TSV, CSV, and JSON Lines exports stream the selected population from the immutable table plus the active overlays. Refined exports preserve original calls and add policy/provenance fields. Collapse and lineage exports add the corresponding representative, abundance, or lineage annotations. JSON/model/TSV/SVG artifacts are generated from the committed result object shown in the interface.

Pipeline mode also exposes a pre-run CLI-config export. It serializes current preprocessing, reference, annotation, study, and stage settings without starting assignment. When selected-reference CHMMAIRRa is enabled, the export performs only the required abortable Kalign preparation and embeds that exact MSA so the standalone CLI has no external alignment dependency. The Results-page exporter remains separate and translates committed interactive settings—including advanced denoising, chimera, lineage, and evidence thresholds—back into the same config fields.

## Reproducibility limits

- A saved session is valid only with the exact linked AIRR content it fingerprints.
- Browser and filesystem quotas can still prevent storage or export; a failed write is not treated as success.
- Cancellation checks occur between bounded batches for main-thread loops, so response time depends on the current batch.
- External browser actions—closing the process, force-killing a tab, revoking a directory permission, or deleting a file outside Swig—cannot be made transactional by the app.

## Reference and exact relationship

- Rubelt F et al. [Adaptive Immune Receptor Repertoire Community recommendations for sharing immune-repertoire sequencing data](https://doi.org/10.1038/ni.3873). *Nature Immunology* (2017). Swig uses AIRR tabular conventions for interoperable assignment output; its IndexedDB representation, linked-session design, rollback checkpoints, and export streaming are Swig engineering and are not specified by AIRR.
