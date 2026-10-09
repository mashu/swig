#!/bin/sh
# Build the WebAssembly modules into web/: vdj_assign_simd.wasm for engines
# with WebAssembly SIMD (every current browser, Node.js 16.4+) and
# vdj_assign.wasm for those without. web/vdj-assign.js picks one at load.
#
#   rustup target add wasm32-unknown-unknown
#   wasm/build.sh
#   python3 -m http.server -d web 8000     # then open http://localhost:8000
set -eu
cd "$(dirname "$0")/.."
target=wasm32-unknown-unknown
build() {
    # Each feature set gets its own target directory, so neither rebuilds
    # the other's dependencies.
    cargo build -p vdj_assign_wasm --profile wasm --target "$target" --target-dir "target/$1"
    cp "target/$1/$target/wasm/vdj_assign_wasm.wasm" "web/$2"
}
build wasm vdj_assign.wasm
RUSTFLAGS="${RUSTFLAGS:-} -C target-feature=+simd128" build wasm-simd vdj_assign_simd.wasm
ls -l web/vdj_assign.wasm web/vdj_assign_simd.wasm
