#!/usr/bin/env bash
# run_matrix.sh [N] [SEED]
#
# Exports GenAIRR's references, builds the IgBLAST databases and the MiXCR
# library, simulates every locus:scenario in $MATRIX (N reads, SEED) and
# runs the tools on each (run_tools.sh). Output goes to $WORK (default
# ./work). Needs: $PYTHON with GenAIRR and numpy, IGBLAST_ROOT, MIXCR_DIR.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
N=${1:-10000}; SEED=${2:-7}
PYTHON=${PYTHON:-python3}
WORK=${WORK:-$HERE/work}
MATRIX=${MATRIX:-"igh:naive igh:shm2 igh:shm5 igh:shm10 igh:realistic igk:shm5 igk:realistic igl:shm5 igl:realistic tcrb:naive tcrb:tcr_realistic"}
mkdir -p "$WORK/runs"
cd "$WORK"
for item in $MATRIX; do
  locus=${item%%:*}; ref=ref_$locus
  [[ -f $ref/reference.json ]] || "$PYTHON" "$HERE/genairr_bench.py" export human_$locus $ref
  if [[ ! -d $ref/igdb && -n ${IGBLAST_ROOT:-} ]]; then
    mkdir -p $ref/igdb
    for k in V D J; do
      [[ -f $ref/$k.plain.fasta ]] && LD_LIBRARY_PATH=$IGBLAST_ROOT/lib $IGBLAST_ROOT/bin/makeblastdb \
        -in $ref/$k.plain.fasta -dbtype nucl -parse_seqids -out $ref/igdb/$k > /dev/null
    done
  fi
  if [[ ! -f $ref/genairr.json && -n ${MIXCR_DIR:-} ]]; then
    ( cd $ref && java -cp "$MIXCR_DIR/target/mixcr-3.0.13.jar:$MIXCR_DIR/target/deps/*" io.repseq.cli.Main \
        compile -f genairr.raw.json genairr.json > /dev/null 2>&1 )
  fi
done
for item in $MATRIX; do
  locus=${item%%:*}; scen=${item##*:}; p=runs/${locus}_${scen}
  [[ -f $p.truth.tsv ]] || "$PYTHON" "$HERE/genairr_bench.py" simulate human_$locus $scen $N $SEED $p || { echo "simulate $item failed"; continue; }
  rm -f $p.times
  "$HERE/run_tools.sh" ref_$locus $p.fasta $p
  echo "== $item $(tr '\n' ' ' < $p.times)"
done
