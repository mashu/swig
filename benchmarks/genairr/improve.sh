#!/usr/bin/env bash
# improve.sh REFDIR RUNPREFIX PRIOR CONFIG LABEL
#
# Prototype post-processing of RUNPREFIX.swig.tsv:
#   RUNPREFIX.recall.tsv        positional V and J tie sets at swig's placement
#   RUNPREFIX.dpost-LABEL.tsv   plus the D posterior call (d_posterior.py),
#                               PRIOR = parametric:SHAPE | genairr, CONFIG =
#                               GenAIRR DataConfig name for `genairr`
# The D posterior runs as $JOBS processes (default 4).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REF=$1; P=$2; PRIOR=$3; CONFIG=$4; LABEL=$5
PYTHON=${PYTHON:-python3}; JOBS=${JOBS:-4}
python3 "$HERE/recall.py" "$REF" "$P.swig.tsv" "$P.recall.tsv" v=keep,d=keep,j=core
if [[ ! -f $REF/D.plain.fasta ]]; then cp "$P.recall.tsv" "$P.dpost-$LABEL.tsv"; exit 0; fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
head -1 "$P.swig.tsv" > "$tmp/header"
tail -n +2 "$P.swig.tsv" | split -n l/"$JOBS" - "$tmp/part."
for part in "$tmp"/part.??; do
  cat "$tmp/header" "$part" > "$part.tsv"
  ( cd "$HERE" && "$PYTHON" d_posterior.py "$REF" "$part.tsv" "$P.recall.tsv" "$part.out" "$PRIOR" "$CONFIG" ) &
done
wait
{ head -1 "$tmp/part.aa.out"; for f in "$tmp"/part.??.out; do tail -n +2 "$f"; done; } > "$P.dpost-$LABEL.tsv"
