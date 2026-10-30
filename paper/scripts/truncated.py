"""Truncated-read sets from IgSim truth, and their D-call scoring.

make: truncated.py make TRUTH.tsv READS.fasta OUTDIR
score: truncated.py score OUTDIR/SET.truth.tsv CALLS.tsv

Sets (forward-strand reads only):
  dcut      read[:d_end]                 reads with a D: a D at the read's open 3' end
  cdr3start read[v_end:]                 reads with a D: start inside the junction
  vcut30    read[:v_end] + 30 random     D-free by construction (no D or J bases)
  vcut100   read[:v_end] + 100 random
  jcut100   100 random + read[j_start:]  D-free by construction
"""
import csv
import random
import sys


def fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def make(truth_path, reads_path, out):
    reads = fasta(reads_path)
    rows = list(csv.DictReader(open(truth_path), delimiter="\t"))
    rng = random.Random(5)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    sets = {k: [] for k in ["dcut", "cdr3start", "vcut30", "vcut100", "jcut100"]}
    for r in rows:
        if r["rev_comp"] != "F" or not r["v_end"] or not r["j_start"]:
            continue
        s = reads[r["sequence_id"]]
        ve, js = int(r["v_end"]), int(r["j_start"])
        if r["d_call"]:
            sets["dcut"].append((r, s[: int(r["d_end"])]))
            sets["cdr3start"].append((r, s[ve:]))
        sets["vcut30"].append((None, s[:ve] + rand(30)))
        sets["vcut100"].append((None, s[:ve] + rand(100)))
        sets["jcut100"].append((None, rand(100) + s[js:]))
    for name, items in sets.items():
        items = items[:2000]
        with open(f"{out}/{name}.fasta", "w") as f, open(f"{out}/{name}.truth.tsv", "w") as t:
            t.write("sequence_id\td_equivalent\n")
            for i, (r, seq) in enumerate(items):
                f.write(f">{name}{i}\n{seq}\n")
                t.write(f"{name}{i}\t{r['d_equivalent'] if r else ''}\n")
        print(name, len(items))


def score(truth_path, calls_path):
    truth = {r["sequence_id"]: r["d_equivalent"] for r in csv.DictReader(open(truth_path), delimiter="\t")}
    calls = {r["sequence_id"]: r["d_call"] for r in csv.DictReader(open(calls_path), delimiter="\t")}
    n = len(truth)
    with_d = sum(1 for k in truth if calls.get(k))
    right = sum(
        1
        for k, eq in truth.items()
        if eq and calls.get(k) and calls[k].split(",")[0] in eq.split(",")
    )
    has_d = sum(1 for eq in truth.values() if eq)
    print(f"reads {n}  D called {100 * with_d / n:.1f}%  right D {100 * right / has_d if has_d else float('nan'):.1f}%")


if __name__ == "__main__":
    {"make": make, "score": score}[sys.argv[1]](*sys.argv[2:])
