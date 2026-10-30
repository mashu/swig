#!/usr/bin/env bash
# timing.sh: throughput of swig and IgBLAST on an otherwise idle machine,
# on the IgSim sets of igsim.sh (run it first).
# Median of three runs for swig; IgBLAST once. Wall clock of the whole command.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
WORK=${WORK:-$ROOT/paper/work}
G=$WORK/germ
SW=$ROOT/engine/target/release/swig
R=$WORK/igblast/root
W=$WORK/igsim
DB=(--v $G/IGHV.fasta --d $G/IGHD.fasta --j $G/IGHJ.fasta)
mkdir -p $W/timing && cd $W/timing

wall() { local t0; t0=$(date +%s.%N); "$@" > /dev/null 2>&1; echo "$(date +%s.%N) - $t0" | bc; }
median3() { printf '%s\n' "$@" | sort -g | sed -n 2p; }

{
for set in igh_mid igh_hard; do
  n=$(grep -c '^>' $W/$set.fasta)
  for threads in 1 4; do
    for post in on off; do
      args=(); [[ $post == off ]] && args=(--no-posterior)
      t=$(median3 $(for i in 1 2 3; do wall "$SW" "${DB[@]}" --threads $threads "${args[@]}" $W/$set.fasta; done))
      echo "swig $set threads=$threads posterior=$post reads=$n seconds=$t reads_per_s=$(echo "$n / $t" | bc)"
    done
  done
  for simd in avx512 avx2 portable scalar; do
    t=$(median3 $(for i in 1 2 3; do wall "$SW" "${DB[@]}" --threads 1 --simd $simd $W/$set.fasta; done))
    "$SW" "${DB[@]}" --threads 1 --simd $simd $W/$set.fasta > $set.$simd.tsv 2> /dev/null
    echo "swig $set simd=$simd threads=1 reads=$n seconds=$t reads_per_s=$(echo "$n / $t" | bc)"
  done
  for simd in avx512 portable scalar; do
    if cmp -s $set.avx2.tsv $set.$simd.tsv; then echo "identical $set avx2 $simd"; else echo "DIFFER $set avx2 $simd"; fi
  done
  if cmp -s $set.avx2.tsv $W/$set.swig.tsv; then echo "identical $set 1-thread 4-thread"; else echo "DIFFER $set 1-thread 4-thread"; fi
  t=$(wall env IGDATA=$R/share/igblast LD_LIBRARY_PATH=$R/lib $R/share/igblast/bin/igblastn \
    -germline_db_V $W/db/IGHV -germline_db_D $W/db/IGHD -germline_db_J $W/db/IGHJ -auxiliary_data $W/db/human_gl.aux \
    -organism human -ig_seqtype Ig -outfmt 19 -num_threads 4 -query $W/$set.fasta -out igblast.$set.tsv)
  echo "igblast $set threads=4 reads=$n seconds=$t reads_per_s=$(echo "$n / $t" | bc -l | cut -c1-6)"
done
} > $W/timing/summary.txt
echo done
