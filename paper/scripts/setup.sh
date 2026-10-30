#!/usr/bin/env bash
# setup.sh: build swig, fetch IgBLAST 1.22 (Bioconda), create a Python venv
# with GenAIRR 2.3.1, and export the human IMGT germline sets as FASTA.
#
# Everything goes under $WORK (default: paper/work, ignored by git).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
WORK=${WORK:-$ROOT/paper/work}
mkdir -p "$WORK"

cargo build --release --manifest-path "$ROOT/engine/Cargo.toml" --bin swig --examples

if [[ ! -x $WORK/igblast/root/share/igblast/bin/igblastn ]]; then
  mkdir -p "$WORK/igblast/root"
  ( cd "$WORK/igblast"
    curl -sSO https://conda.anaconda.org/bioconda/linux-64/igblast-1.22.0-h6a68c12_1.tar.bz2
    curl -sSO https://conda.anaconda.org/bioconda/linux-64/ncbi-vdb-3.1.1-h9948957_3.tar.bz2
    for f in *.tar.bz2; do tar -xjf "$f" -C root; done )
fi

if [[ ! -x $WORK/venv/bin/python ]]; then
  python3 -m venv "$WORK/venv"
  "$WORK/venv/bin/pip" install -q GenAIRR==2.3.1 numpy
fi

# Human IG and TR germline sets from the IMGT/GENE-DB release 202632-7 pack
mkdir -p "$WORK/germ"
python3 - "$ROOT/public/references/imgt-202632-7-swig-0.7.json.gz" "$WORK/germ" <<'EOF'
import gzip, json, os, sys
pack, out = sys.argv[1:3]
d = json.load(gzip.open(pack))
human = next(s for s in d["species"] if s["name"] == "Homo sapiens")
for locus, segments in human["loci"].items():
    for segment, records in segments.items():
        with open(os.path.join(out, f"{locus}{segment}.fasta"), "w") as f:
            for name, sequence, *_ in records:
                f.write(f">{name}\n{sequence}\n")
EOF
echo "setup done: $WORK"
