#!/usr/bin/env bash
# genairr.sh: the GenAIRR benchmark of the paper.
#   1. simulate the 11 sets (seed 7) and run swig and IgBLAST on them
#      (benchmarks/genairr/run_matrix.sh; MiXCR is not run here, see README.md)
#   2. swig with the D priors off, and with the max-marginal D decoder
#   3. a fresh seed (8) for the IGH sets: swig, swig with the D priors off, IgBLAST
#   4. score everything against GenAIRR's tie sets -> $WORK/genairr/scores.txt
# Run setup.sh first.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
WORK=${WORK:-$ROOT/paper/work}
B=$ROOT/benchmarks/genairr
SW=$ROOT/engine/target/release/swig
W=$WORK/genairr
export PYTHON=$WORK/venv/bin/python IGBLAST_ROOT=$WORK/igblast/root TOOLS="swig igblast" WORK=$W
mkdir -p $W

"$B/run_matrix.sh" 10000 7

cd $W
swig_run() { # ref reads out [args]
  local ref=$1 reads=$2 out=$3; shift 3
  local d=(); [[ -f $ref/D.swig.fasta ]] && d=(--d $ref/D.swig.fasta)
  [[ -s $out ]] || "$SW" --v $ref/V.swig.fasta "${d[@]}" --j $ref/J.swig.fasta --threads 4 "$@" "$reads" > "$out" 2> /dev/null
}
for run in igh_naive igh_shm2 igh_shm5 igh_shm10 igh_realistic tcrb_naive tcrb_tcr_realistic; do
  swig_run ref_${run%%_*} runs/$run.fasta runs/$run.swig-priors-off.tsv --d-trim off --np-length off
  swig_run ref_${run%%_*} runs/$run.fasta runs/$run.swig-maxmarg.tsv --d-allele max-marginal
done

mkdir -p seed8
for scen in naive shm5 shm10 realistic; do
  p=seed8/igh_$scen
  [[ -f $p.truth.tsv ]] || $PYTHON $B/genairr_bench.py simulate human_igh $scen 10000 8 $p
  [[ -s $p.igblast.tsv ]] || "$B/run_tools.sh" ref_igh $p.fasta $p
  swig_run ref_igh $p.fasta $p.swig-priors-off.tsv --d-trim off --np-length off
done

cd $B
{
for run in igh_naive igh_shm2 igh_shm5 igh_shm10 igh_realistic igk_shm5 igk_realistic igl_shm5 igl_realistic tcrb_naive tcrb_tcr_realistic; do
  args=()
  for t in swig swig-priors-off swig-maxmarg igblast; do
    [[ -f $W/runs/$run.$t.tsv ]] && args+=("$t=$W/runs/$run.$t.tsv")
  done
  echo "## $run"
  python3 genairr_bench.py score $W/runs/$run.truth.tsv "${args[@]}"
done
for scen in naive shm5 shm10 realistic; do
  p=$W/seed8/igh_$scen
  echo "## seed8 igh_$scen"
  python3 genairr_bench.py score $p.truth.tsv swig=$p.swig.tsv swig-priors-off=$p.swig-priors-off.tsv igblast=$p.igblast.tsv
done
} > $W/scores.txt 2>&1
cat $W/scores.txt
