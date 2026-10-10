#!/usr/bin/env bash
# Build Swig's aligner (the Rust crate in engine/) as a WebAssembly module
# with SIMD: public/swig.wasm, or $SWIG_WASM_OUTPUT. Needs a Rust toolchain
# with the wasm32-unknown-unknown target (rustup target add wasm32-unknown-unknown).
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output="${SWIG_WASM_OUTPUT:-$project_dir/public/swig.wasm}"
target=wasm32-unknown-unknown
target_dir="$project_dir/engine/target/wasm-simd"

if ! command -v cargo >/dev/null 2>&1; then
  echo "Install a Rust toolchain (https://rustup.rs) to build the aligner." >&2
  exit 2
fi

RUSTFLAGS="${RUSTFLAGS:-} -C target-feature=+simd128" cargo build \
  --manifest-path "$project_dir/engine/Cargo.toml" \
  -p swig_wasm --profile wasm --target "$target" --target-dir "$target_dir"

mkdir -p "$(dirname "$output")"
cp "$target_dir/$target/wasm/swig_wasm.wasm" "$output"
ls -l "$output"
