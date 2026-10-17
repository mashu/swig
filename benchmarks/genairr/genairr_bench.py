"""GenAIRR benchmark for swig, IgBLAST and MiXCR.

Truth is GenAIRR's own call: the tie set of alleles with the most matching
bases at the true allele's germline positions (``v_call``/``d_call``/
``j_call``). A tool's call is scored per segment as

    |truth ∩ call| / |truth ∪ call|      (empty and empty scores 1)

so a single correct allele against a two-allele tie scores 1/2, and so does
a two-allele call against a single truth allele.

    python genairr_bench.py export   CONFIG OUTDIR
    python genairr_bench.py simulate CONFIG SCENARIO N SEED OUTPREFIX
    python genairr_bench.py score    TRUTH.tsv TOOL=CALLS.tsv ...
"""
import csv
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict

csv.field_size_limit(1 << 30)

SCENARIOS = {
    # No SHM, no sequencing artefacts: pure recombination.
    "naive": lambda e: e,
    # S5F somatic hypermutation at two levels.
    "shm2": lambda e: e.mutate(model="s5f", rate=0.02),
    "shm5": lambda e: e.mutate(model="s5f", rate=0.05),
    "shm10": lambda e: e.mutate(model="s5f", rate=0.10),
    # SHM plus library-prep and sequencing artefacts.
    "realistic": lambda e: (
        e.mutate(model="s5f", rate=0.05)
        .primer_trim_5prime(length=(0, 15))
        .primer_trim_3prime(length=(0, 5))
        .polymerase_indels(count=(0, 1), insertion_prob=0.5)
        .pcr_amplify(rate=0.002)
        .sequencing_errors(rate=0.002)
        .ambiguous_base_calls(count=(0, 2))
    ),
    # TCR: no SHM; library-prep and sequencing artefacts only.
    "tcr_realistic": lambda e: (
        e.primer_trim_5prime(length=(0, 15))
        .primer_trim_3prime(length=(0, 5))
        .polymerase_indels(count=(0, 1), insertion_prob=0.5)
        .pcr_amplify(rate=0.002)
        .sequencing_errors(rate=0.002)
        .ambiguous_base_calls(count=(0, 2))
    ),
}


def experiment(config):
    import GenAIRR as ga
    e = ga.Experiment.on(config)
    if "tcr" in config:
        # GenAIRR's TRB set holds ORF/pseudogene anchors it calls curatable.
        e = e.allow_curatable_refdata()
    return e


def alleles(rd, kind):
    size = getattr(rd, f"{kind}_pool_size")()
    get = getattr(rd, f"{kind}_allele")
    out = []
    for i in range(size):
        a = get(i)
        seq = a.seq()
        if isinstance(seq, bytes):
            seq = seq.decode()
        anchor = a.anchor
        # Keep an anchor only on its conserved codon (V Cys; J Trp/Phe):
        # GenAIRR marks ORF/pseudogene alleles with placeholder anchors.
        codon = seq.upper()[anchor:anchor + 3] if anchor is not None else ""
        if kind == "v" and codon not in ("TGT", "TGC"):
            anchor = None
        if kind == "j" and codon not in ("TGG", "TTT", "TTC"):
            anchor = None
        out.append((a.name, seq.upper(), anchor, list(a.subregions or [])))
    return out


def export(config, outdir):
    """Write the GenAIRR reference for swig, IgBLAST and MiXCR."""
    os.makedirs(outdir, exist_ok=True)
    rd = experiment(config).refdata
    has_d = rd.has_d() if callable(rd.has_d) else rd.has_d
    kinds = ["v", "d", "j"] if has_d else ["v", "j"]
    meta = {}
    for kind in kinds:
        records = alleles(rd, kind)
        meta[kind] = records
        with open(os.path.join(outdir, f"{kind.upper()}.plain.fasta"), "w") as plain, \
             open(os.path.join(outdir, f"{kind.upper()}.swig.fasta"), "w") as swig:
            for name, seq, anchor, regions in records:
                plain.write(f">{name}\n{seq}\n")
                keys = []
                if anchor is not None and kind in ("v", "j"):
                    keys.append(f"anchor={anchor}")
                if kind == "v":
                    starts = {r[0]: r[1] for r in regions}
                    for label, key in [("FWR1", "fr1"), ("CDR1", "cdr1"), ("FWR2", "fr2"), ("CDR2", "cdr2"), ("FWR3", "fr3")]:
                        if label in starts:
                            keys.append(f"{key}={starts[label]}")
                swig.write(f">{name} {' '.join(keys)}\n{seq}\n".replace(" \n", "\n"))
    locus = (rd.identity() if callable(rd.identity) else rd.identity)["locus"]
    json.dump({"config": config, "locus": locus, "has_d": has_d,
               "counts": {k: len(v) for k, v in meta.items()}}, open(os.path.join(outdir, "reference.json"), "w"))
    # IgBLAST: auxiliary J data (0-based coding frame start, chain type,
    # 0-based last CDR3 base).
    chain = {"IGH": "JH", "IGK": "JK", "IGL": "JL", "TRB": "JB", "TRA": "JA", "TRD": "JD", "TRG": "JG"}.get(locus, "JH")
    with open(os.path.join(outdir, "igblast.aux"), "w") as aux:
        for name, seq, anchor, _ in meta["j"]:
            if anchor is not None:
                aux.write(f"{name}\t{anchor % 3}\t{chain}\t{anchor - 1}\t0\n")
    # MiXCR (repseqio) library: per-allele anchor points.
    genes = []
    for kind in kinds:
        for name, seq, anchor, regions in meta[kind]:
            points = {}
            if kind == "v":
                starts = {r[0]: r[1] for r in regions}
                for label, point in [("FWR1", "FR1Begin"), ("CDR1", "CDR1Begin"), ("FWR2", "FR2Begin"), ("CDR2", "CDR2Begin"), ("FWR3", "FR3Begin")]:
                    if label in starts:
                        points[point] = starts[label]
                points.setdefault("FR1Begin", 0)
                if anchor is not None:
                    points["CDR3Begin"] = anchor
                points["VEnd"] = len(seq)
            elif kind == "d":
                points = {"DBegin": 0, "DEnd": len(seq)}
            else:
                points = {"JBegin": 0, "FR4End": len(seq)}
                if anchor is not None:
                    points["FR4Begin"] = anchor + 3
            genes.append({
                "baseSequence": f"file://{kind.upper()}.plain.fasta#{name}",
                "name": name, "geneType": kind.upper(), "isFunctional": True,
                "chains": [locus.replace("TR", "TR")], "anchorPoints": points,
            })
    json.dump([{"taxonId": 9606, "speciesNames": ["hs"], "genes": genes}],
              open(os.path.join(outdir, "genairr.raw.json"), "w"), indent=1)


def simulate(config, scenario, n, seed, prefix):
    e = SCENARIOS[scenario](experiment(config).recombine())
    result = e.run_records(n=n, seed=seed, expose_provenance=True)
    fields = ["sequence_id", "v_call", "d_call", "j_call", "truth_v_call", "truth_d_call", "truth_j_call",
              "n_mutations", "mutation_rate", "n_v_mutations", "n_d_mutations", "n_j_mutations",
              "v_trim_5", "v_trim_3", "d_trim_5", "d_trim_3", "j_trim_5", "j_trim_3",
              "d_sequence_start", "d_sequence_end", "np1_length", "np2_length",
              "n_indels", "n_quality_errors", "n_pcr_errors", "end_loss_5_length", "end_loss_3_length",
              "productive", "rev_comp"]
    with open(f"{prefix}.fasta", "w") as fasta, open(f"{prefix}.truth.tsv", "w", newline="") as truth:
        writer = csv.writer(truth, delimiter="\t")
        writer.writerow(fields)
        for rec in result:
            fasta.write(f">{rec['sequence_id']}\n{rec['sequence']}\n")
            writer.writerow(["" if rec.get(f) is None else rec.get(f) for f in fields])


def call_set(value):
    return frozenset(x.strip() for x in (value or "").split(",") if x.strip())


def jaccard(truth, call):
    if not truth and not call:
        return 1.0
    if not truth or not call:
        return 0.0
    return len(truth & call) / len(truth | call)


def read_calls(path):
    """sequence_id -> {v,d,j: frozenset}; AIRR TSV (swig, IgBLAST) or the
    normalised TSV the MiXCR runner writes."""
    out = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["sequence_id"]] = {s: call_set(row.get(f"{s}_call", "")) for s in "vdj"}
    return out


def score(truth_path, tools):
    truth = list(csv.DictReader(open(truth_path), delimiter="\t"))
    has_d = any(r["d_call"] for r in truth)
    segments = "vdj" if has_d else "vj"
    results = {}
    for name, path in tools:
        calls = read_calls(path)
        totals = defaultdict(float)
        for r in truth:
            c = calls.get(r["sequence_id"], {s: frozenset() for s in "vdj"})
            for s in segments:
                totals[s] += jaccard(call_set(r[f"{s}_call"]), c[s])
        results[name] = {s: totals[s] / len(truth) for s in segments}
    return results, len(truth)


def main(argv):
    command = argv[1]
    if command == "export":
        export(argv[2], argv[3])
    elif command == "simulate":
        simulate(argv[2], argv[3], int(argv[4]), int(argv[5]), argv[6])
    elif command == "score":
        tools = [a.split("=", 1) for a in argv[3:]]
        results, n = score(argv[2], tools)
        print(f"{n} reads")
        for name, r in results.items():
            print(f"{name:24} " + "  ".join(f"{s.upper()} {100 * v:6.2f}" for s, v in r.items()))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
