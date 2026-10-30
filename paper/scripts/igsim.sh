#!/usr/bin/env bash
# igsim.sh: every IgSim experiment of the paper on the human IMGT IGH set.
#   1. swig 0.39.7 against IgBLAST 1.22 (easy, mid, hard, 30 % SHM), forward strand
#   2. negative controls: IGK reads against the IGH database, random reads
#   3. swig on both strands with C (orientation, calibration, C)
#   4. ablations paired against the default (mid, hard)
#   5. truncated reads (open-flank cost), D error attribution, search vs model errors
# Run setup.sh first. Output: $WORK/igsim. Accuracy does not depend on load;
# run timing.sh separately on an idle machine.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
WORK=${WORK:-$ROOT/paper/work}
G=$WORK/germ
SW=$ROOT/engine/target/release/swig
EX=$ROOT/engine/target/release/examples
R=$WORK/igblast/root
W=$WORK/igsim
DB=(--v $G/IGHV.fasta --d $G/IGHD.fasta --j $G/IGHJ.fasta)
mkdir -p $W/db && cd $W

# IgBLAST databases; its auxiliary file lacks IGHJ5*03 and *04, which share
# IGHJ5*01's frame and CDR3 end (swig prepare's anchors agree)
for k in V D J; do
  [[ -f db/IGH$k.nsq || -f db/IGH$k.00.nsq ]] || LD_LIBRARY_PATH=$R/lib $R/bin/makeblastdb -in $G/IGH$k.fasta -dbtype nucl -parse_seqids -out db/IGH$k > /dev/null
done
AUX=db/human_gl.aux
if [[ ! -f $AUX ]]; then
  cp $R/share/igblast/optional_file/human_gl.aux $AUX
  row=$(grep -P '^IGHJ5\*01\t' $AUX | cut -f2-)
  printf 'IGHJ5*03\t%s\nIGHJ5*04\t%s\n' "$row" "$row" >> $AUX
fi
igblast() { # reads out
  env IGDATA=$R/share/igblast LD_LIBRARY_PATH=$R/lib $R/share/igblast/bin/igblastn \
    -germline_db_V db/IGHV -germline_db_D db/IGHD -germline_db_J db/IGHJ -auxiliary_data $AUX \
    -organism human -ig_seqtype Ig -outfmt 19 -num_threads 4 -query "$1" -out "$2" 2> "$2.log"
}
swig() { # reads out [args]
  local reads=$1 out=$2; shift 2
  [[ -s $out ]] || "$SW" "${DB[@]}" --threads 4 "$@" "$reads" > "$out" 2> "$out.log"
}

# 1-2. Cross-tool sets and negative controls
for p in easy mid hard; do
  [[ -f igh_$p.fasta ]] || "$SW" simulate "${DB[@]}" --preset $p --reads 10000 --seed 11 --rc 0 --out igh_$p 2> /dev/null
done
[[ -f igh_shm30.fasta ]] || "$SW" simulate "${DB[@]}" --preset hard --shm 0.3 --reads 3000 --seed 11 --rc 0 --out igh_shm30 2> /dev/null
[[ -f igk_mid.fasta ]] || "$SW" simulate --v $G/IGKV.fasta --j $G/IGKJ.fasta --preset mid --reads 1000 --seed 11 --rc 0 --out igk_mid 2> /dev/null
[[ -f random.fasta ]] || python3 - igh_mid.fasta random.fasta <<'EOF'
import random, sys
lens = [len(l.strip()) for l in open(sys.argv[1]) if not l.startswith('>')]
rng = random.Random(11)
with open(sys.argv[2], 'w') as f:
    for i in range(2000):
        f.write(f">rand{i}\n{''.join(rng.choice('ACGT') for _ in range(rng.choice(lens)))}\n")
EOF
for set in igh_easy igh_mid igh_hard igh_shm30 igk_mid random; do
  swig $set.fasta $set.swig.tsv
  [[ -s $set.igblast.tsv ]] || igblast $set.fasta $set.igblast.tsv
done

# 3. Both strands, with C
for p in easy mid hard; do
  [[ -f ighc_$p.fasta ]] || "$SW" simulate "${DB[@]}" --c $G/IGHC.fasta --preset $p --reads 10000 --seed 12 --rc 0.5 --out ighc_$p 2> /dev/null
  swig ighc_$p.fasta ighc_$p.swig.tsv --c $G/IGHC.fasta
done

# 4. Ablations
declare -A ABL=(
  [priors_off]="--d-trim off --np-length off"
  [open_flank_off]="--d-open-flank off"
  [no_guards]="--no-pair-guards"
  [fixed]="--calls fixed"
  [max_marginal]="--d-allele max-marginal"
  [sensitive_d]="--profile sensitive-d"
  [truth_profile]="--profile truth"
  [igblast_profile]="--profile igblast"
  [no_posterior]="--no-posterior"
  [no_anchor_window]="--no-anchor-window"
)
for p in mid hard; do
  for name in "${!ABL[@]}"; do
    swig igh_$p.fasta igh_$p.$name.tsv ${ABL[$name]}
  done
done

# 5. Truncated reads, D errors, search vs model errors
mkdir -p trunc
for p in mid hard; do
  mkdir -p trunc/$p
  [[ -f trunc/$p/vcut100.fasta ]] || python3 $HERE/truncated.py make igh_$p.truth.tsv igh_$p.fasta trunc/$p
  for set in dcut cdr3start vcut30 vcut100 jcut100; do
    for cfg in default open_flank_off no_guards no_guards_open_flank_off; do
      case $cfg in
        default) args=();;
        open_flank_off) args=(--d-open-flank off);;
        no_guards) args=(--no-pair-guards);;
        no_guards_open_flank_off) args=(--no-pair-guards --d-open-flank off);;
      esac
      swig trunc/$p/$set.fasta trunc/$p/$set.$cfg.tsv "${args[@]}"
      echo "$p $set $cfg $(python3 $HERE/truncated.py score trunc/$p/$set.truth.tsv trunc/$p/$set.$cfg.tsv)"
    done
  done
done > trunc/summary.txt
for p in mid hard; do
  [[ -s d_errors_$p.txt ]] || "$EX/d_errors" "${DB[@]}" --preset $p --reads 10000 > d_errors_$p.txt 2>&1
  [[ -s diagnose_$p.txt ]] || "$EX/diagnose" "${DB[@]}" --preset $p --reads 2000 > diagnose_$p.txt 2>&1
done

python3 $HERE/collect.py $W $W/results.json "$SW"
python3 $HERE/controls.py $W > $W/controls.txt
echo "igsim done: $W"
