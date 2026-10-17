"""D call by posterior over (allele, offset, 5' trim, 3' trim) (prototype).

The junction is the read between swig's V end and J start. A hypothesis
h = (A, o, t5, t3) puts germline A[t5 : len(A)-t3] at read o+t5; the other
junction bases are N. Likelihood ratio against "all N":

    prod over D bases of  (1-mu)/(1/4) if matched, (mu/3)/(1/4) if not

times trim and NP-length priors. Each h implies GenAIRR's truth set for it:
the alleles with the most matches at the same germline indices over the D
region, widened by GenAIRR's narrowing walk into the N bases. The reported
set maximises the expected Jaccard score under the posterior.

    python d_posterior.py REFDIR SWIG.tsv BASE.tsv OUT.tsv PRIOR [CONFIG]

BASE.tsv supplies the V and J calls (e.g. a recall.py output). PRIOR:
    parametric[:SHAPE]  negative-binomial trims (means 5 and 7) and NP
                        lengths (mean 8), uniform genes; SHAPE 1 (the
                        default) is geometric, the best tested
    generic             flat trims and NP lengths
    genairr[:PARTS]     GenAIRR's own trim, NP-length and gene-usage tables
                        from DataConfig CONFIG; PARTS a subset of
                        trim,np,use
"""
import csv
import math
import sys
from collections import defaultdict

import numpy as np

from recall import fasta

csv.field_size_limit(1 << 30)
MAX_TRIM = 20


def priors(kind, config, d_names):
    """log P(t5 | allele), log P(t3 | allele), log P(np1), log P(np2), log P(allele).
    kind: generic | genairr | genairr:trim,np,use (any subset)."""
    flat = np.zeros(MAX_TRIM + 1)
    parts = set(kind.split(":", 1)[1].split(",")) if ":" in kind else {"trim", "np", "use"}
    if kind.startswith("parametric"):
        # Negative binomial log-pmf with the given mean and shape r; no
        # simulator knowledge, typical values for human IGH.
        def negbin(mean, r, n_max):
            p_ = r / (r + mean)
            ks = np.arange(n_max + 1)
            lg = np.array([math.lgamma(k + r) - math.lgamma(k + 1) - math.lgamma(r) for k in ks])
            out = lg + r * math.log(p_) + ks * math.log(1 - p_)
            return out - np.log(np.exp(out).sum())
        shape = float(kind.split(":", 1)[1]) if ":" in kind else 1.0
        a5 = negbin(5.0, shape, MAX_TRIM); a3 = negbin(7.0, shape, MAX_TRIM)
        npl = negbin(8.0, 2.0, 200)
        t5 = {n: a5 for n in d_names}; t3 = {n: a3 for n in d_names}
        np1 = np2 = lambda n: float(npl[min(n, 200)])
        # Uniform over genes, each gene's mass split among its alleles.
        per_gene = defaultdict(int)
        for n in d_names:
            per_gene[n.split("*")[0]] += 1
        allele = {n: -math.log(per_gene[n.split("*")[0]]) for n in d_names}
        return t5, t3, np1, np2, allele
    if kind == "generic":
        t5 = {n: flat for n in d_names}
        t3 = {n: flat for n in d_names}
        np1 = np2 = lambda n: 0.0
        allele = {n: 0.0 for n in d_names}
        return t5, t3, np1, np2, allele
    import GenAIRR.data as gd
    dc = getattr(gd, config)
    def table(d, gene):
        fam = gene.split("-")[0]
        dist = d.get(fam, {}).get(gene) or next(iter(d.get(fam, {}).values()), {})
        arr = np.full(MAX_TRIM + 1, 1e-6)
        for k, v in (dist or {}).items():
            if 0 <= int(k) <= MAX_TRIM:
                arr[int(k)] += v
        return np.log(arr / arr.sum())
    t5 = {n: table(dc.trim_dicts["D_5"], n.split("*")[0]) for n in d_names}
    t3 = {n: table(dc.trim_dicts["D_3"], n.split("*")[0]) for n in d_names}
    def np_table(key):
        d = dc.NP_lengths[key]
        tot = sum(d.values())
        return lambda n: math.log(max(d.get(n, 0.0), 1e-6) / tot)
    use = dc.gene_use_dict.get("D", {})
    per_gene = defaultdict(int)
    for n in d_names:
        per_gene[n.split("*")[0]] += 1
    allele = {n: math.log(max(use.get(n.split("*")[0], 1e-3), 1e-6) / per_gene[n.split("*")[0]]) for n in d_names}
    np1, np2 = np_table("NP1"), np_table("NP2")
    if "trim" not in parts:
        t5 = {n: flat for n in d_names}; t3 = {n: flat for n in d_names}
    if "np" not in parts:
        np1 = np2 = lambda n: 0.0
    if "use" not in parts:
        allele = {n: 0.0 for n in d_names}
    return t5, t3, np1, np2, allele


def truth_set(M, alleles_list, lengths, s, e, lo_i, hi_i, t5, t3):
    """GenAIRR's tie set. M[b, i] = read base at germline index i of the
    placement matches allele b (False past b's end or outside the read);
    the D region is germline [s, e); N bases cover indices [lo_i, s) and
    [e, hi_i)."""
    scores = M[:, s:e].sum(axis=1).astype(int)
    def walk(indices):
        for i in indices:
            top = scores.max()
            tied = scores == top
            hit = tied & M[:, i]
            n_hit = int(hit.sum())
            if not (0 < n_hit < int(tied.sum())):
                return
            scores[hit] += 1
    walk(range(s - 1, max(lo_i, s - t5) - 1, -1))
    walk(range(e, min(hi_i, e + t3)))
    top = scores.max()
    if top == 0:
        return frozenset(alleles_list)
    return frozenset(alleles_list[k] for k in np.flatnonzero(scores == top))


def call_d(read, v_end, j_start, alleles, mu, pri, top_k=200):
    t5p, t3p, np1p, np2p, allelep = pri
    if j_start - v_end < 1:
        return ""
    names = list(alleles)
    lengths = np.array([len(alleles[n]) for n in names])
    maxlen = int(lengths.max())
    codes = np.full((len(names), maxlen), 0, dtype=np.uint8)
    for k, n in enumerate(names):
        codes[k, :len(alleles[n])] = np.frombuffer(alleles[n].encode(), dtype=np.uint8)
    rd = np.frombuffer(read.encode(), dtype=np.uint8)
    w_match = math.log((1 - mu) / 0.25)
    w_miss = math.log((mu / 3) / 0.25)
    J = j_start - v_end
    np1a = np.array([np1p(n) for n in range(J + 1)])
    np2a = np.array([np2p(n) for n in range(J + 1)])
    T = MAX_TRIM + 1
    t5_idx = np.arange(T)[:, None]
    t3_idx = np.arange(T)[None, :]
    all_sc, all_meta = [], []
    match_at = {}
    for o in range(v_end - maxlen + 1, j_start):
        i0 = max(0, v_end - o); i1 = min(maxlen, j_start - o)
        if i1 <= i0:
            continue
        M = np.zeros((len(names), maxlen), dtype=bool)
        M[:, i0:i1] = (codes[:, i0:i1] == rd[o + i0:o + i1][None, :]) & (np.arange(i0, i1)[None, :] < lengths[:, None])
        match_at[o] = (M, i0, i1)
        W = np.where(M, w_match, w_miss)
        W[:, :i0] = 0.0; W[:, i1:] = 0.0
        CS = np.concatenate([np.zeros((len(names), 1)), np.cumsum(W, axis=1)], axis=1)
        for k, name in enumerate(names):
            L = int(lengths[k])
            if L <= i0:
                continue
            e = L - t3_idx                                   # region end index
            valid = (t5_idx >= i0) & (e <= min(L, i1)) & (e > t5_idx) & (t5_idx <= L) & (t3_idx <= L)
            if not valid.any():
                continue
            ec = np.clip(e, 0, maxlen); tc = np.clip(t5_idx, 0, maxlen)
            sc = CS[k][ec] - CS[k][tc] + t5p[name][:, None] + t3p[name][None, :] + allelep[name]
            n1 = np.clip(o + t5_idx - v_end, 0, J); n2 = np.clip(j_start - (o + e), 0, J)
            sc = sc + np1a[n1] + np2a[n2]
            ii, jj = np.nonzero(valid)
            all_sc.append(sc[ii, jj])
            all_meta.append(np.stack([np.full(ii.size, k), np.full(ii.size, o), ii, jj], axis=1))
    if not all_sc:
        return ""
    scores = np.concatenate(all_sc); meta = np.concatenate(all_meta)
    top = np.argpartition(-scores, min(top_k, scores.size - 1))[:top_k]
    top = top[np.argsort(-scores[top])]
    best = scores[top[0]]
    kept = [(float(scores[t]), names[int(meta[t, 0])], int(meta[t, 1]), int(meta[t, 2]), int(meta[t, 3])) for t in top if scores[t] > best - 12]
    mass = defaultdict(float)
    cache = {}
    for sc, name, o, t5, t3 in kept:
        L = len(alleles[name])
        key = (o, t5, L - t3, t5, t3)
        if key not in cache:
            M, i0, i1 = match_at[o]
            cache[key] = truth_set(M, names, lengths, t5, L - t3, i0, i1, t5, t3)
        mass[cache[key]] += math.exp(sc - best)
    z = sum(mass.values())
    sets = sorted(mass.items(), key=lambda kv: -kv[1])
    candidates = [s for s, _ in sets[:8]]
    candidates += [a | b for k, a in enumerate(candidates) for b in candidates[k + 1:4]]
    def expected(R):
        return sum(m * len(R & T) / len(R | T) for T, m in mass.items()) / z
    R = max(candidates, key=expected)
    return ",".join(sorted(R))


def main(refdir, swig_tsv, base_tsv, out_tsv, kind, config="HUMAN_IGH_OGRDB"):
    D = fasta(f"{refdir}/D.plain.fasta")
    pri = priors(kind, config, list(D))
    base = {r["sequence_id"]: r for r in csv.DictReader(open(base_tsv), delimiter="\t")}
    with open(out_tsv, "w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["sequence_id", "v_call", "d_call", "j_call"])
        for row in csv.DictReader(open(swig_tsv), delimiter="\t"):
            sid = row["sequence_id"]
            read = row["sequence"].upper()
            if row.get("rev_comp") == "T":
                read = read[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
            b = base[sid]
            d = b["d_call"]
            if row["v_call"] and row["j_call"] and row["v_sequence_end"] and row["j_sequence_start"]:
                v_end = int(row["v_sequence_end"]); j_start = int(row["j_sequence_start"]) - 1
                try:
                    mu = min(0.3, max(0.005, 1 - float(row["v_identity"])))
                except ValueError:
                    mu = 0.05
                d = call_d(read, v_end, j_start, D, mu, pri) or d
            w.writerow([sid, b["v_call"], d, b["j_call"]])


if __name__ == "__main__":
    main(*sys.argv[1:])
