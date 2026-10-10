# Changelog

Each release's section is its GitHub release notes; the tag `vX.Y.Z` must
match the crate version (`.github/workflows/release.yml`). Push the tag, or
run the Release workflow by hand on `main`, which creates it.

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
- A page that assigns reads in the browser, on Web Workers, with demo
  data; nothing is uploaded. It is one self-contained HTML file, both
  modules included, that runs opened from disk (no server) or from any
  static host. Its output equals the command line's up to a rare
  last-digit difference in a posterior.

### Release assets

| File | Contents |
|---|---|
| `vdj-assign-0.1.0-x86_64-unknown-linux-gnu.tar.gz` | Linux x86-64 (glibc 2.35 or later) |
| `vdj-assign-0.1.0-aarch64-unknown-linux-gnu.tar.gz` | Linux AArch64 (glibc 2.35 or later) |
| `vdj-assign-0.1.0-x86_64-apple-darwin.tar.gz` | macOS on Intel |
| `vdj-assign-0.1.0-aarch64-apple-darwin.tar.gz` | macOS on Apple silicon |
| `vdj-assign-0.1.0-x86_64-pc-windows-msvc.zip` | Windows x86-64 (statically linked C runtime) |
| `vdj-assign-0.1.0.html` | The browser page as one file: open it, no server needed |
| `vdj-assign-web-0.1.0.zip` | The same page (`index.html`), the JavaScript API and both WebAssembly modules |
| `SHA256SUMS` | SHA-256 of every file above |
