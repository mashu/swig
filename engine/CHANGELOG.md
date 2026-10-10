# Changelog

Each release's section is its GitHub release notes; the tag `vX.Y.Z` must
match the crate version (`.github/workflows/release.yml`). Merging a version
bump and its section to `main` publishes the release and creates the tag; a
pushed tag or a hand run of the Release workflow on `main` does the same.

## [0.1.1] - 2026-10-10

### The browser page needs no server

- The page is one self-contained HTML file, `vdj-assign-0.1.1.html`
  (2.5 MB): its scripts and both WebAssembly modules are inlined and its
  workers start from blobs, so it fetches nothing and runs opened straight
  from disk (double-click it) or from any static host. The 0.1.0 page
  needed an HTTP server, because browsers do not let a page opened from disk
  load module scripts, fetch `.wasm` files or start workers from files.
- The demo is deployed to GitHub Pages from `main` whenever CI passes:
  https://mashu.github.io/vdj_assign/
- `assignInWorkers` takes a `createWorker` option; its default, a module
  worker on `worker.js`, is unchanged.
- The browser test opens the page from disk and from a server in Chromium,
  Firefox and WebKit, and checks in Chromium, which reports every load,
  that the page opened from disk requests nothing but itself and its blobs.

### Releases

- Merging a version bump and its changelog section to `main` publishes the
  release; a pushed tag or a hand run of the Release workflow still does.

### Release assets

| File | Contents |
|---|---|
| `vdj-assign-0.1.1-x86_64-unknown-linux-gnu.tar.gz` | Linux x86-64 (glibc 2.35 or later) |
| `vdj-assign-0.1.1-aarch64-unknown-linux-gnu.tar.gz` | Linux AArch64 (glibc 2.35 or later) |
| `vdj-assign-0.1.1-x86_64-apple-darwin.tar.gz` | macOS on Intel |
| `vdj-assign-0.1.1-aarch64-apple-darwin.tar.gz` | macOS on Apple silicon |
| `vdj-assign-0.1.1-x86_64-pc-windows-msvc.zip` | Windows x86-64 (statically linked C runtime) |
| `vdj-assign-0.1.1.html` | The browser page as one file: open it, no server needed |
| `vdj-assign-web-0.1.1.zip` | The same page (`index.html`), the JavaScript API and both WebAssembly modules |
| `SHA256SUMS` | SHA-256 of every file above |

## [0.1.0] - 2026-10-09

The first release.

### Command line (`vdj-assign`)

- **Assignment**: V(D)J (and optionally constant-region) germline alleles for
  FASTA or FASTQ reads, from files or stdin, as AIRR Rearrangement TSV with
  per-segment posterior probabilities, FWR1–FWR4 and CDR1–CDR3 regions and
  the junction.
- **`prepare`**: the deduplicated, annotated germline database (anchors and
  IMGT region starts, from headers, IMGT gaps, motifs or sibling alleles)
  with an `annotation.tsv` for review.
- **`simulate`** and **`evaluate`**: reads from the IgSim model with their
  truth, and germline-equivalence-aware scoring of any tool's AIRR TSV,
  with a paired McNemar comparison of two tools on the same reads.
- **Kernels**: AVX2, NEON or portable SIMD lanes, chosen at run time;
  AVX-512BW and the scalar kernel with `--simd`. Output is byte-identical
  across them and across x86-64 and AArch64 Linux; on macOS and Windows,
  as in the browser, a posterior can rarely differ in its last printed
  digit where the platform's `exp` differs in the last bit.
- `-V` / `--version` prints the version.

### WebAssembly

- Two modules, `vdj_assign_simd.wasm` (WebAssembly SIMD) and
  `vdj_assign.wasm` (scalar), with no imports. `vdj-assign.js` fetches the
  right one in browsers and Web Workers; in Node.js 18 or later (the
  bundle's `package.json` marks it an ES module) give `load()` the bytes.
- A page (`index.html`) that assigns reads in the browser, on Web Workers,
  with demo data; nothing is uploaded. Its output equals the command line's
  up to a rare last-digit difference in a posterior.

### Release assets

| File | Contents |
|---|---|
| `vdj-assign-0.1.0-x86_64-unknown-linux-gnu.tar.gz` | Linux x86-64 (glibc 2.35 or later) |
| `vdj-assign-0.1.0-aarch64-unknown-linux-gnu.tar.gz` | Linux AArch64 (glibc 2.35 or later) |
| `vdj-assign-0.1.0-x86_64-apple-darwin.tar.gz` | macOS on Intel |
| `vdj-assign-0.1.0-aarch64-apple-darwin.tar.gz` | macOS on Apple silicon |
| `vdj-assign-0.1.0-x86_64-pc-windows-msvc.zip` | Windows x86-64 (statically linked C runtime) |
| `vdj-assign-web-0.1.0.zip` | The page, the JavaScript API and both WebAssembly modules |
| `SHA256SUMS` | SHA-256 of every archive above |
