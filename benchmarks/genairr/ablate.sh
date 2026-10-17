#!/usr/bin/env bash
# ablate.sh REFDIR RUNPREFIX
#
# What the D posterior's gain is made of: runs dposterior/ (build it with
# `cargo build --release --manifest-path dposterior/Cargo.toml`) on
# RUNPREFIX.swig.tsv under each variant below, keeps swig's V and J calls,
# and scores D. With RUNPREFIX.fasta present it also runs swig with
# `--d-allele max-marginal`, the engine's own posterior decoding.
#
# Variants: which priors are on (trim, N length; the gene-split allele prior
# is on unless all are flat) and how the call is decided (mapset: the tie
# set of the best hypothesis; jaccard: the set of best expected score).
# `oneoff` keeps each allele at its offset of best local score only, as the
# engine's DAG holds one D placement per allele.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REF=$1; P=$2
BIN=${DPOSTERIOR:-$HERE/dposterior/target/release/dposterior}
SWIG=${SWIG:-$HERE/../../engine/target/release/swig}
export NATS=${NATS:-8} TOPK=${TOPK:-50}
VARIANTS=(
  "flat:TRIM=flat NP=flat ALLELE=flat DECIDE=mapset"
  "trim:NP=flat DECIDE=mapset"
  "trim+np:DECIDE=mapset"
  "trim+np-geom:NP=geom DECIDE=mapset"
  "trim+np,jaccard:"
  "oneoff,flat:ONEOFF=1 TRIM=flat NP=flat ALLELE=flat DECIDE=mapset"
  "oneoff,trim:ONEOFF=1 NP=flat DECIDE=mapset"
  "oneoff,trim+np-geom:ONEOFF=1 NP=geom DECIDE=mapset"
)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
args=("swig 0.39.6=$P.swig.tsv")
if [[ -f $P.fasta ]]; then
  "$SWIG" --v "$REF/V.swig.fasta" --d "$REF/D.swig.fasta" --j "$REF/J.swig.fasta" --threads "${THREADS:-4}" \
    --d-allele max-marginal "$P.fasta" > "$tmp/swig-mm.tsv" 2> /dev/null
  args+=("swig --d-allele max-marginal=$tmp/swig-mm.tsv")
fi
i=0
for v in "${VARIANTS[@]}"; do
  name=${v%%:*}; read -ra envs <<< "${v#*:}"
  env "${envs[@]}" "$BIN" "$REF/D.plain.fasta" "$REF/J.plain.fasta" "$P.swig.tsv" "$tmp/$i.d.tsv" 2> /dev/null
  python3 "$HERE/merge_calls.py" "$P.swig.tsv" "$tmp/$i.d.tsv" "$tmp/$i.tsv"
  args+=("$name=$tmp/$i.tsv"); i=$((i + 1))
done
python3 "$HERE/genairr_bench.py" score "$P.truth.tsv" "${args[@]}"
