# Swig 0.39.1

`swig-cli` now runs the native Swig aligner. Each standalone executable
carries the `swig` binary built for its platform and streams reads through
it. The WebAssembly module stays as the fallback. On 10,000 human IGH reads
on four cores, `swig-cli --vdj` takes 1.5 s instead of 7.3 s. Swig Web and
every other tool are unchanged.

## What changes for users

- **Native aligner in `swig-cli`.** On first use the executable copies its
  embedded aligner to the user cache, once per version and content hash:
  `~/.cache/swig-cli`, `~/Library/Caches/swig-cli` or
  `%LOCALAPPDATA%\swig-cli`. If none of these is usable, it uses a directory
  in the temporary directory that is private to the user, and checks that
  it is private before running anything from it. `SWIG_CLI_CACHE_DIR`
  overrides the location. The copy must report the CLI's own version
  (`swig --version`). It then runs as a child process:
  - It gets the prepared germline FASTA and the same options the
    WebAssembly engine gets.
  - `--threads` equals the worker count.
  - Records go to its standard input, and AIRR rows come back in input
    order.
  - The pipeline starts one child per dataset.
- **Choosing the aligner.** `--aligner auto|native|wasm` selects the aligner;
  `annotation.aligner` in a config file and the `SWIG_CLI_ALIGNER`
  environment variable do the same.
  - `auto` (the default) falls back to WebAssembly and says why when no
    usable binary exists or the child fails to start.
  - `native` fails instead of falling back.
  - `wasm` always uses WebAssembly.
  - Double-D screening runs only on WebAssembly: `auto` switches to it, and
    `native` refuses.
  - `SWIG_CLI_ALIGNER_PATH` names another `swig` binary of the same version
    and takes precedence over the embedded one. The repository's `cli/`
    bundle uses the native aligner this way.
- **Static native binaries.** The Linux aligners now link glibc statically.
  Built dynamically on the release runners they needed glibc 2.34, while
  `swig-cli` runs from glibc 2.17. Without this, RHEL 8, Ubuntu 20.04,
  Debian 11 and many HPC hosts would silently have stayed on WebAssembly.
  The Windows aligner links its C runtime statically, so it runs without the
  Visual C++ Redistributable. Both changes apply to the `swig-aligner-*`
  release assets too.
- **One WebAssembly compilation.** The CLI compiles `swig.wasm` once and
  shares the module with every worker. Before, each worker compiled its own
  copy, and under the release runtime four workers ran no faster than two.
  The WebAssembly route now scales with workers.
- **`--vdj` uses every core.** It uses every core when neither `--workers`
  nor `-num_threads` is given, and with `--workers 0`. An explicit count is
  still exact. The pipeline's `annotation.workers: 0` still selects up to
  eight, because that count also sizes the collapse, denoising and
  CHMMAIRRa pools.
- **Records without an identifier.** A FASTA or FASTQ record with an empty
  header, or an AIRR row with an empty or missing `sequence_id`, is named
  `sequence_N` by its position N in the input, before batching. Before, each
  WebAssembly worker numbered its own records. That repeated names, and the
  names depended on the worker count.
- **Blank AIRR sequences.** An AIRR row whose `sequence` is only whitespace
  is skipped, in Swig Web and in the CLI, as the aligner always did. Before,
  it counted as a record that the aligner then dropped.
- **Worker threads on Windows.** The Windows executable could not start any
  worker thread: the WebAssembly aligner pool, or the pools for post-analysis
  and CHMMAIRRa. Bun 1.3 resolves an embedded worker given as a file URL, but
  then reads its backslash path (`B:\~BUN\root\…`) from disk and fails with
  ENOENT. The 0.39.0 executable has the same code and Bun version. The CLI
  now names embedded workers with a relative specifier, which Bun maps to
  the embedded entry. The new Windows release job runs these workers.
- **Cleaner failures.** When the aligner fails, the CLI prints one message
  and exits with status 1. It no longer adds an unhandled-rejection trace
  or a dump of bundle source, or loses the real error behind "worker pool
  closed". A failed WebAssembly fallback no longer leaves the process
  running. Native start-up errors report the engine's own message, not its
  usage text, and a binary that fails `--version` reports the loader's
  reason.

## Output

The native and WebAssembly routes write the same AIRR. In the release tests,
on 17,500 synthetic reads and on FASTA, FASTQ and AIRR pipeline datasets,
the output is byte-identical. On the timing set below, both 0.39.0 and
0.39.1 on either route are byte-identical. The one known difference is the
one the methods already describe: a posterior that lies within a last-bit
difference of `exp` of a printed rounding tie can differ in that digit,
because WebAssembly uses a portable libm. On simulated IGH reads with
`-J_penalty -2`, 1 row in 10,000 differed, in the tenth decimal of
`j_posterior`.

## Timing

| Executable | Aligner | Wall time (median of 5) |
|---|---|---|
| `swig-cli-linux-x64` 0.39.0 | WebAssembly, 4 workers | 7.32 s |
| `swig-cli-linux-x64-modern` 0.39.0 | WebAssembly, 4 workers | 7.11 s |
| `swig-cli` 0.39.1, baseline target | native, 4 threads (default) | 1.53 s |
| `swig-cli` 0.39.1, baseline target | `--aligner wasm`, 4 workers | 3.74 s |
| `swig` 0.39.1 native binary alone | 4 threads | 1.35 s |

Method:
- **Input.** 10,000 simulated human IGH reads: the IGH mid and hard sets of
  [`BENCHMARK_SWIG_ALIGNER_0.39.0.md`](BENCHMARK_SWIG_ALIGNER_0.39.0.md),
  with human IGH V, D and J and no C.
- **Command.** `swig-cli --vdj -query reads.fasta -germline_db_V V.fasta
  -germline_db_D D.fasta -germline_db_J J.fasta -outfmt 19 -out out.tsv`.
- **Machine.** A 4-core Intel Xeon at 2.10 GHz (cloud VM), Linux x86-64.
- **Runs.** Wall time includes start-up and germline preparation. Each
  configuration had one warm-up run (which also fills the extraction cache),
  then five rounds of all configurations, interleaved.
- **Executables.** The 0.39.0 rows are the published release assets. The
  0.39.1 executable was built from this release with Bun 1.3.14
  (`--target=bun-linux-x64-baseline`), embedding the statically linked
  aligner.
- **Native binary alone.** This row reads the same FASTA through
  `swig --v --d --j --threads 4`. Its references are not the CLI's prepared
  ones, so its rows are not compared.

## Release build

Releases stay on Bun 1.3.14. Bun 1.4.2 closes `node:readline` while a
paused consumer still holds unread lines ("readline was closed"), which
truncates streamed input. 0.39.0's own code fails the same way when built
with 1.4.2. The Windows executable is still cross-compiled on Linux, but
the release now runs its native-route tests on a Windows runner before it
publishes.
