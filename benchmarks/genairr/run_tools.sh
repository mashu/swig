#!/usr/bin/env bash
# run_tools.sh REFDIR READS.fasta OUTPREFIX [extra swig args...]
#
# Runs the tools named in $TOOLS (default: swig igblast mixcr) on READS
# against the reference REFDIR (from `genairr_bench.py export`) and writes
#   OUTPREFIX.swig.tsv        swig AIRR output
#   OUTPREFIX.igblast.tsv     IgBLAST AIRR output (-outfmt 19)
#   OUTPREFIX.mixcr.tsv       MiXCR, top hit per segment
#   OUTPREFIX.mixcr-ties.tsv  MiXCR, every hit tied with the top score
#   OUTPREFIX.times           wall-clock seconds per tool
#
# Environment:
#   SWIG          swig binary (default: engine/target/release/swig)
#   IGBLAST_ROOT  unpacked bioconda igblast + ncbi-vdb (bin/, lib/, share/igblast/)
#   MIXCR_DIR     MiXCR 3.0.13 source tree built with `mvn package` and
#                 `mvn dependency:copy-dependencies -DoutputDirectory=target/deps`
#   THREADS       threads per tool (default 4)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REF=$(cd "$1" && pwd); READS=$(realpath "$2"); OUT=$(realpath -m "$3"); shift 3
SWIG=${SWIG:-$HERE/../../engine/target/release/swig}
THREADS=${THREADS:-4}
TOOLS=${TOOLS:-swig igblast mixcr}
HAS_D=$(python3 -c "import json;print(json.load(open('$REF/reference.json'))['has_d'])")
LOCUS=$(python3 -c "import json;print(json.load(open('$REF/reference.json'))['locus'])")
timed() { local label=$1; shift; local t0; t0=$(date +%s.%N); "$@"; local rc=$?; echo "$label $(echo "$(date +%s.%N) - $t0" | bc) s" >> "$OUT.times"; return $rc; }

if [[ " $TOOLS " == *" swig "* ]]; then
  D_ARGS=(); [[ $HAS_D == True ]] && D_ARGS=(--d "$REF/D.swig.fasta")
  timed swig "$SWIG" --v "$REF/V.swig.fasta" "${D_ARGS[@]}" --j "$REF/J.swig.fasta" --threads "$THREADS" "$@" "$READS" \
    > "$OUT.swig.tsv" 2> "$OUT.swig.log" || { echo "swig failed"; tail "$OUT.swig.log"; }
fi

if [[ " $TOOLS " == *" igblast "* ]]; then
  : "${IGBLAST_ROOT:?set IGBLAST_ROOT}"
  SEQTYPE=Ig; [[ $LOCUS == TR* ]] && SEQTYPE=TCR
  D_DB=(-germline_db_D "$REF/igdb/V"); [[ $HAS_D == True ]] && D_DB=(-germline_db_D "$REF/igdb/D")
  timed igblast env IGDATA="$IGBLAST_ROOT/share/igblast" LD_LIBRARY_PATH="$IGBLAST_ROOT/lib" "$IGBLAST_ROOT/share/igblast/bin/igblastn" \
    -germline_db_V "$REF/igdb/V" "${D_DB[@]}" -germline_db_J "$REF/igdb/J" -auxiliary_data "$REF/igblast.aux" \
    -organism human -ig_seqtype $SEQTYPE -outfmt 19 -num_threads "$THREADS" -query "$READS" -out "$OUT.igblast.tsv" 2> "$OUT.igblast.log"
fi

if [[ " $TOOLS " == *" mixcr "* ]]; then
  : "${MIXCR_DIR:?set MIXCR_DIR}"
  CP="$MIXCR_DIR/target/mixcr-3.0.13.jar:$MIXCR_DIR/target/deps/*"
  ( cd "$REF" && timed mixcr java -Xmx6g -cp "$CP" com.milaboratory.mixcr.cli.Main \
      align -f -s hs -b genairr -t "$THREADS" --report "$OUT.mixcr.report" "$READS" "$OUT.vdjca" > "$OUT.mixcr.log" 2>&1 \
    && java -cp "$CP" com.milaboratory.mixcr.cli.Main \
      exportAlignments -f -readIds -vHitsWithScore -dHitsWithScore -jHitsWithScore "$OUT.vdjca" "$OUT.mixcr.raw.tsv" >> "$OUT.mixcr.log" 2>&1 ) \
    || { echo "mixcr failed"; tail "$OUT.mixcr.log"; }
  python3 - "$READS" "$OUT.mixcr.raw.tsv" "$OUT" <<'EOF'
import csv, re, sys
reads, raw, out = sys.argv[1:4]
ids = [l[1:].split()[0] for l in open(reads) if l.startswith(">")]
def hits(field):
    return [(m.group(1), float(m.group(2))) for m in re.finditer(r"([^,()]+)\(([-0-9.]+)\)", field or "")]
rows = {}
for r in csv.DictReader(open(raw), delimiter="\t"):
    cols = {k.lower().replace(" ", ""): v for k, v in r.items()}
    key = cols.get("readid") or cols.get("readids")
    rows[int(key.split(",")[0])] = {s: hits(next((v for k, v in cols.items() if k.startswith(f"all{s}hits")), "")) for s in "vdj"}
for variant in ["mixcr", "mixcr-ties"]:
    with open(f"{out}.{variant}.tsv", "w", newline="") as h:
        w = csv.writer(h, delimiter="\t"); w.writerow(["sequence_id", "v_call", "d_call", "j_call"])
        for i, sid in enumerate(ids):
            calls = []
            for s in "vdj":
                hs = rows.get(i, {}).get(s, [])
                if not hs:
                    calls.append(""); continue
                if variant == "mixcr":
                    calls.append(hs[0][0])
                else:
                    top = max(x[1] for x in hs); calls.append(",".join(n for n, sc in hs if sc == top))
            w.writerow([sid] + calls)
EOF
fi
