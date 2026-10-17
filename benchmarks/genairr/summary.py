"""Markdown table of V/D/J scores for every run and tool variant present."""
import os
import sys

from genairr_bench import score

RUNS = ["igh_naive", "igh_shm2", "igh_shm5", "igh_shm10", "igh_realistic",
        "igk_shm5", "igk_realistic", "igl_shm5", "igl_realistic", "tcrb_naive", "tcrb_tcr_realistic"]
TOOLS = [("IgBLAST 1.22", "igblast"), ("MiXCR 3.0.13", "mixcr"), ("MiXCR, all top ties", "mixcr-ties"),
         ("swig 0.39.6", "swig"), ("+ positional J", "recall"),
         ("+ D posterior, geometric priors", "final-parametric"), ("+ D posterior, GenAIRR priors", "final-genairr")]


def main(root):
    print("| Run | Tool | V | D | J |")
    print("|---|---|---:|---:|---:|")
    for run in RUNS:
        truth = f"{root}/{run}.truth.tsv"
        if not os.path.exists(truth):
            continue
        tools = [(label, f"{root}/{run}.{suffix}.tsv") for label, suffix in TOOLS if os.path.exists(f"{root}/{run}.{suffix}.tsv")]
        results, n = score(truth, tools)
        best = {s: max(r.get(s, 0) for r in results.values()) for s in "vdj"}
        for label, r in results.items():
            cells = []
            for s in "vdj":
                if s not in r:
                    cells.append("–")
                else:
                    v = f"{100 * r[s]:.2f}"
                    cells.append(f"**{v}**" if r[s] == best[s] else v)
            print(f"| {run} ({n}) | {label} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "runs")
