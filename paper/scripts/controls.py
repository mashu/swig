"""Negative controls: calls made on random reads and on IGK reads run
against the IGH database, by swig and IgBLAST.

controls.py IGSIM_DIR
"""
import csv
import sys

csv.field_size_limit(1 << 30)


def counts(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    n = len(rows)
    has = lambda r, k: bool((r.get(k) or "").strip())
    return {
        "reads": n,
        "any": sum(has(r, "v_call") or has(r, "j_call") for r in rows),
        "V": sum(has(r, "v_call") for r in rows),
        "D": sum(has(r, "d_call") for r in rows),
        "J": sum(has(r, "j_call") for r in rows),
        "junction": sum(has(r, "junction") for r in rows),
    }


def main(root):
    for name in ["random", "igk_mid"]:
        for tool in ["swig", "igblast"]:
            c = counts(f"{root}/{name}.{tool}.tsv")
            print(name, tool, " ".join(f"{k}={v}" for k, v in c.items()))


if __name__ == "__main__":
    main(*sys.argv[1:])
