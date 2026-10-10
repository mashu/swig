#!/bin/sh
# Build the WebAssembly modules into web/: vdj_assign_simd.wasm for engines
# with WebAssembly SIMD (every current browser, Node.js 16.4+) and
# vdj_assign.wasm for those without. web/vdj-assign.js picks one at load.
# Then write web/index.html: web/page.html with app.js, vdj-assign.js,
# worker.js and both modules inlined (web/standalone.js), one file that runs
# opened from disk or from any static host, with no server.
#
#   rustup target add wasm32-unknown-unknown
#   wasm/build.sh
#   then open web/index.html in a browser
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

fail() {
    echo "wasm/build.sh: $*" >&2
    exit 1
}

# web/standalone.js with each `//@ name` line replaced by that file: a module
# as a base64 string, vdj-assign.js without `export` (and with import.meta,
# a module-only token, replaced: the page never takes the paths it names),
# worker.js and app.js without their import of vdj-assign.js.
script() {
    while IFS= read -r line || [ -n "$line" ]; do
        case $line in
        '//@ vdj_assign_simd.wasm' | '//@ vdj_assign.wasm')
            file=${line#'//@ '}
            key=baseline
            [ "$file" = vdj_assign_simd.wasm ] && key=simd
            printf "  %s: '" "$key"
            base64 < "web/$file" | tr -d '\n'
            printf "',\n"
            ;;
        '//@ vdj-assign.js')
            sed -e 's/^export //' -e 's/import\.meta\.url/self.location.href/g' web/vdj-assign.js
            ;;
        '//@ worker.js' | '//@ app.js')
            sed "/^import .* from '.\/vdj-assign.js';\$/d" "web/${line#'//@ '}"
            ;;
        '//@ '*)
            fail "web/standalone.js: unknown inclusion '$line'"
            ;;
        *)
            printf '%s\n' "$line"
            ;;
        esac
    done < web/standalone.js
}

marker='<script type="module" src="app.js"></script>'
[ "$(grep -cxF "$marker" web/page.html)" = 1 ] ||
    fail "web/page.html must hold the line '$marker' once"
bundle=web/index.html.tmp
script > "$bundle.js"
# The script is inline: nothing in it may end or escape the script element,
# and no module syntax may remain.
if grep -qi -e '</script' -e '<!--' "$bundle.js"; then
    fail "the inlined script holds '</script' or '<!--'"
fi
if grep -q -e 'import\.meta' -e '^import ' -e '^export ' "$bundle.js"; then
    fail "the inlined script still holds module syntax"
fi
{
    awk -v m="$marker" '$0 == m { exit } { print }' web/page.html
    echo '<script>'
    cat "$bundle.js"
    echo '</script>'
    awk -v m="$marker" 'after { print } $0 == m { after = 1 }' web/page.html
} > "$bundle"
rm "$bundle.js"
mv "$bundle" web/index.html
ls -l web/index.html
