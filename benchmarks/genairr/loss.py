"""Where accuracy is lost: loss = 1 - |T∩C|/|T∪C|, split by error class.

    python loss.py TRUTH.tsv TOOL=CALLS.tsv ...
"""
import csv
import sys
from collections import Counter, defaultdict

from genairr_bench import call_set, jaccard, read_calls

csv.field_size_limit(1 << 30)


def gene(allele):
    return allele.split("*")[0]


def classify(t, c):
    if t == c:
        return "exact"
    if not c:
        return "missing call"
    if not t:
        return "spurious call"
    if c < t:
        return "subset of tie (under-reported)"
    if c > t:
        return "superset of tie (over-reported)"
    if t & c:
        return "partial overlap"
    if {gene(a) for a in t} & {gene(a) for a in c}:
        return "wrong allele, right gene"
    return "wrong gene"


def main(truth_path, tools):
    truth = list(csv.DictReader(open(truth_path), delimiter="\t"))
    n = len(truth)
    segments = "vdj" if any(r["d_call"] for r in truth) else "vj"
    for name, path in tools:
        calls = read_calls(path)
        print(f"\n### {name}")
        for s in segments:
            loss = defaultdict(float)
            count = Counter()
            for r in truth:
                t = call_set(r[f"{s}_call"])
                c = calls.get(r["sequence_id"], {}).get(s, frozenset())
                k = classify(t, c)
                count[k] += 1
                loss[k] += 1 - jaccard(t, c)
            total = sum(loss.values())
            print(f"{s.upper()}: accuracy {100 * (1 - total / n):.2f} %, loss {100 * total / n:.2f} points")
            for k, v in sorted(loss.items(), key=lambda kv: -kv[1]):
                if k != "exact":
                    print(f"   {k:34} reads {count[k]:6}  loss {100 * v / n:5.2f} pts")


if __name__ == "__main__":
    main(sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]])
