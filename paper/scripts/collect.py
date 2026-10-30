"""Collect swig evaluate outputs into one JSON for the paper's tables.

collect.py IGSIM_DIR OUT.json [SWIG_BINARY]
"""
import json
import re
import subprocess
import sys
from pathlib import Path

SW = "swig"


def evaluate(truth, *calls):
    return subprocess.run([SW, "evaluate", "--truth", str(truth), *map(str, calls)],
                          capture_output=True, text=True).stdout


def parse_single(text):
    out = {}
    m = re.search(r"oriented ([\d.]+)%\s+unassigned (\d+)\s+junction exact ([\d.NaN]+)% \(of (\d+)\), spurious ([\d.NaN]+)% \(of (\d+)\)", text)
    if m:
        out.update(oriented=float(m[1]), unassigned=int(m[2]), junction=float(m[3]), junction_reads=int(m[4]),
                   spurious_junction=float(m[5]), junction_absent=int(m[6]))
    for seg in ["V", "D", "J"]:
        m = re.search(rf"^{seg}\s+(\d+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+(\d+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+([\d.NaN]+)", text, re.M)
        if m:
            keys = ["reads", "called", "equiv", "exact", "listed", "gene", "names", "spurious", "start_err", "end_err", "start_eq", "end_eq"]
            out[seg] = {k: float(v) for k, v in zip(keys, m.groups())}
    def block(header):
        # The rows under a table header, up to the next unindented line.
        m = re.search(rf"^{header}.*\n((?:\s+.*\n)*)", text, re.M)
        return m[1] if m else ""

    rem = {}
    for m in re.finditer(r"^\s+(none|1-3|4-7|8-11|12-15|16\+)\s+(\d+)\s+([\d.NaN]+)\s+([\d.NaN]+)", block("D remnant"), re.M):
        rem[m[1]] = {"reads": int(m[2]), "called": float(m[3]), "correct": float(m[4])}
    out["remnant"] = rem
    c = {}
    for m in re.finditer(r"^\s+(none|\d+-\d+|\d+\+)\s+(\d+)\s+(.*)$", block("C length"), re.M):
        c[m[1]] = {"reads": int(m[2]), "vals": [float(v) for v in m[3].split()]}
    out["c_length"] = c
    m = re.search(r"ECE %\s+(.*)\n\s+Brier\s+(.*)", text)
    if m:
        out["ece"] = [float(x) for x in m[1].split()]
        out["brier"] = [float(x) for x in m[2].split()]
    bins = []
    for m in re.finditer(r"^\s+(\d\.\d{3})-(\d(?:\.\d{3})?)\s+(.*)$", text, re.M):
        vals = m[3].split()
        bins.append({"lo": float(m[1]), "hi": float(m[2]), "vals": [float(v) for v in vals]})
    out["bins"] = bins
    return out


def parse_paired(text):
    out = {}
    for m in re.finditer(r"^(V allele|D allele|D 1-3 nt|D 4-7 nt|D 8\+ nt|no D call|J allele|junction|C allele|no C call)\s+(\d+)\s+([\d.NaN]+)\s+([\d.NaN]+)\s+([-+][\d.]+|NaN)\s+([\d.NaN]+)\s+(\d+)\s+(\d+)\s+(\S+)", text, re.M):
        out[m[1]] = {"reads": int(m[2]), "first": float(m[3]), "second": float(m[4]), "delta": float(m[5]),
                     "ci": float(m[6]), "first_only": int(m[7]), "second_only": int(m[8]), "p": float(m[9])}
    for m in re.finditer(r"^(D allele|D presence)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([-+][\d.]+)\s+([\d.NaN]+)\s+(\S+)$", text, re.M):
        out["brier " + m[1]] = {"reads": int(m[2]), "first": float(m[3]), "second": float(m[4]), "delta": float(m[5]),
                                "ci": m[6], "t": m[7]}
    return out


def main(root, out_path, binary=None):
    global SW
    SW = binary or SW
    root = Path(root)
    res = {"single": {}, "paired": {}}
    for f in sorted(root.glob("*.tsv")):
        name = f.name
        if name.endswith(".truth.tsv") or name.startswith("trial"):
            continue
        set_name = name.split(".")[0]
        truth = root / f"{set_name}.truth.tsv"
        if not truth.exists():
            continue
        res["single"][name] = parse_single(evaluate(truth, f))
        base = root / f"{set_name}.swig.tsv"
        if f != base and base.exists():
            res["paired"][name] = parse_paired(evaluate(truth, base, f))
    Path(out_path).write_text(json.dumps(res, indent=1))
    print(len(res["single"]), len(res["paired"]))


if __name__ == "__main__":
    main(*sys.argv[1:])
