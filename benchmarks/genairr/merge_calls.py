"""merge_calls.py VJ.tsv D.tsv OUT.tsv: V and J calls from the first file,
D calls from the second."""
import csv
import sys

csv.field_size_limit(1 << 30)
vj = {r["sequence_id"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
with open(sys.argv[3], "w", newline="") as h:
    w = csv.writer(h, delimiter="\t")
    w.writerow(["sequence_id", "v_call", "d_call", "j_call"])
    for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
        b = vj[r["sequence_id"]]
        w.writerow([r["sequence_id"], b["v_call"], r["d_call"], b["j_call"]])
