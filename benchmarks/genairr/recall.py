"""Positional re-call of swig's V/D/J tie sets (prototype).

GenAIRR's truth is the set of alleles with the most read bases matching at
the true allele's germline index, over the segment's true read region (plus
NP bases while they narrow the tie). swig's tie set comes from local
alignment scores, which stop at the last good column. This re-scores every
allele of the segment at the called placement's germline frame.

    python recall.py REFDIR SWIG.tsv OUT.tsv MODE[,MODE...]

Modes, per segment letter (e.g. v=full,d=core,j=narrow):
    keep    swig's call unchanged
    core    positional tie set over swig's aligned span only
    narrow  core, then GenAIRR's walk: extend one base at a time outwards
            while the base strictly narrows the tie set
    full    core plus every read base up to the neighbouring segment (or the
            germline's end), all counted
    endpost posterior over the segment's true junction-side end (V 3', J 5')
            and over the frames of swig's co-optimal alleles; each end
            implies GenAIRR's tie set (positional over the segment, then the
            narrowing walk); report the set of highest expected Jaccard
"""
import math
import csv
import sys

csv.field_size_limit(1 << 30)


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line.upper())
    return {k: "".join(v) for k, v in out.items()}


def load(refdir):
    refs = {}
    for s in "VDJ":
        try:
            refs[s.lower()] = fasta(f"{refdir}/{s}.plain.fasta")
        except FileNotFoundError:
            refs[s.lower()] = {}
    return refs


def span(row, s):
    try:
        return (int(row[f"{s}_sequence_start"]) - 1, int(row[f"{s}_sequence_end"]),
                int(row[f"{s}_germline_start"]) - 1, int(row[f"{s}_germline_end"]))
    except (KeyError, ValueError):
        return None


class Frame:
    """Read position -> germline index through the segment's CIGAR (AIRR
    extended CIGAR: S, N, =, X, M, I, D). Inserted read bases map to None,
    as GenAIRR skips them; outside the alignment the nearest end's diagonal
    is used."""

    def __init__(self, cigar, q0, q1, g0, g1):
        import re
        self.map = {}
        q, g = None, None
        for n, op in re.findall(r"(\d+)([SNMXDI=])", cigar or ""):
            n = int(n)
            if op == "S" and q is None:
                q = n
            elif op == "N" and g is None:
                g = n
            elif op in "=XM":
                q = q or 0; g = g or 0
                for k in range(n):
                    self.map[q + k] = g + k
                q += n; g += n
            elif op == "I":
                q = q or 0
                for k in range(n):
                    self.map[q + k] = None
                q += n
            elif op == "D":
                g = (g or 0) + n
        if not self.map:
            self.map = {q0 + k: g0 + k for k in range(q1 - q0)}
        self.left = q0 - g0
        self.right = q1 - g1

    def __call__(self, p):
        if p in self.map:
            return self.map[p]
        return p - (self.left if p < min(self.map) else self.right)


def positional(read, alleles, offset, positions):
    frame = offset if callable(offset) else (lambda p, o=offset: p - o)
    scores = {}
    idx = [(p, frame(p)) for p in positions]
    for name, g in alleles.items():
        scores[name] = sum(1 for p, gi in idx if gi is not None and 0 <= gi < len(g) and read[p] == g[gi])
    return scores


def tie(scores):
    top = max(scores.values())
    return {n for n, v in scores.items() if v == top}, top


def narrows(scores, read, alleles, offset, p):
    """GenAIRR's rule: the base at p matches some but not all alleles tied at the max."""
    gi = offset(p) if callable(offset) else p - offset
    top = max(scores.values())
    tied = [n for n, v in scores.items() if v == top]
    hit = [n for n in tied if gi is not None and 0 <= gi < len(alleles[n]) and read[p] == alleles[n][gi]]
    return 0 < len(hit) < len(tied)


def recall_segment(row, s, alleles, mode, bounds):
    calls = [c for c in row[f"{s}_call"].split(",") if c]
    sp = span(row, s)
    if mode == "keep" or not calls or sp is None or calls[0] not in alleles:
        return row[f"{s}_call"]
    q0, q1, g0, g1 = sp
    read = row["sequence_alignment_full"]
    frame = Frame(row.get(f"{s}_cigar"), q0, q1, g0, g1)
    offset = q0 - g0
    primary = alleles[calls[0]]
    lo_bound, hi_bound = bounds
    core = list(range(q0, q1))
    if mode == "full":
        right_end = min(hi_bound, offset + len(primary))
        left_start = max(lo_bound, offset) if s in "dj" else q0
        core = list(range(left_start, q0)) + core + list(range(q1, right_end))
    scores = positional(read, alleles, frame, core)
    if mode == "narrow":
        p = q1
        while p < min(hi_bound, len(read)) and narrows(scores, read, alleles, frame, p):
            gi = frame(p)
            for n, g in alleles.items():
                if 0 <= gi < len(g) and read[p] == g[gi]:
                    scores[n] += 1
            p += 1
        if s in "dj":
            p = q0 - 1
            while p >= max(lo_bound, 0) and narrows(scores, read, alleles, frame, p):
                gi = frame(p)
                for n, g in alleles.items():
                    if 0 <= gi < len(g) and read[p] == g[gi]:
                        scores[n] += 1
                p -= 1
    best, _ = tie(scores)
    # Keep swig's first name first when it is in the set.
    ordered = [c for c in calls if c in best] + sorted(best - set(calls))
    return ",".join(ordered)



def negbin(mean, r, n_max):
    p_ = r / (r + mean)
    out = [math.lgamma(k + r) - math.lgamma(k + 1) - math.lgamma(r) + r * math.log(p_) + k * math.log(1 - p_) for k in range(n_max + 1)]
    z = math.log(sum(math.exp(x) for x in out))
    return [x - z for x in out]


TRIM_PRIOR = {"v": negbin(1.5, 1.0, 60), "j": negbin(4.5, 2.0, 60)}


def frames_of(row, s, alleles, margin=8):
    """(allele, offset) frames: swig's call (its span's offset) and the
    near-tied alternatives within `margin` score points, each at its own
    reported coordinates."""
    sp = span(row, s)
    q0, q1, g0, _ = sp
    calls = [c for c in row[f"{s}_call"].split(",") if c]
    try:
        top = float(row[f"{s}_score"])
    except (KeyError, ValueError):
        top = None
    out = [(c, q0 - g0) for c in calls if c in alleles]
    for alt in (row.get(f"{s}_alternatives") or "").split(";"):
        f = alt.split("|")
        if len(f) >= 7 and f[0] in alleles and top is not None and float(f[1]) >= top - margin:
            out.append((f[0], (int(f[3]) - 1) - (int(f[5]) - 1)))
    seen, uniq = set(), []
    for fr in out:
        if fr not in seen:
            seen.add(fr); uniq.append(fr)
    return uniq


def endpost_segment(row, s, alleles, bounds, mu):
    calls = [c for c in row[f"{s}_call"].split(",") if c]
    sp = span(row, s)
    if not calls or sp is None:
        return row[f"{s}_call"]
    q0, q1, g0, _ = sp
    read = row["sequence_alignment_full"]
    lo_bound, hi_bound = bounds
    hi_bound = min(hi_bound, len(read))
    lm, lx, ln = math.log(1 - mu), math.log(mu / 3), math.log(0.25)
    mass = {}
    for frame, offset in frames_of(row, s, alleles):
        g = alleles[frame]
        hyps = []
        if s == "v":
            start = max(q0, offset)
            g_end_read = offset + len(g)
            base = [lm if read[p] == g[p - offset] else lx for p in range(start, min(g_end_read, hi_bound))]
            prefix = [0.0]
            for x in base:
                prefix.append(prefix[-1] + x)
            for e in range(max(q1 - 3, start + 1), min(hi_bound, g_end_read) + 1):
                ll = prefix[e - start] + ln * (hi_bound - e) + TRIM_PRIOR["v"][min(g_end_read - e, 60)]
                hyps.append((ll, e))
        else:
            end = min(q1, offset + len(g))
            lo = max(lo_bound, offset, 0)
            base = [lm if read[p] == g[p - offset] else lx for p in range(lo, end)]
            suffix = [0.0]
            for x in reversed(base):
                suffix.append(suffix[-1] + x)
            for b in range(lo, min(q0 + 3, end - 1) + 1):
                ll = suffix[end - b] + ln * (b - lo_bound) + TRIM_PRIOR["j"][min(b - offset, 60)]
                hyps.append((ll, b))
        for ll, edge in hyps:
            if s == "v":
                positions = range(max(q0, offset), edge)
                cap = offset + len(g) - edge
                walk = range(edge, min(hi_bound, edge + cap))
            else:
                positions = range(edge, min(q1, offset + len(g)))
                cap = edge - offset
                walk = range(edge - 1, max(lo_bound, edge - cap) - 1, -1)
            scores = positional(read, alleles, offset, positions)
            for pos in walk:
                if not narrows(scores, read, alleles, offset, pos):
                    break
                for n, a in alleles.items():
                    if 0 <= pos - offset < len(a) and read[pos] == a[pos - offset]:
                        scores[n] += 1
            T = frozenset(tie(scores)[0])
            mass[T] = mass.get(T, 0.0) + ll
            mass[T] = mass[T]
            mass.setdefault(("w", T), []).append(ll)
    sets = {k[1]: v for k, v in mass.items() if isinstance(k, tuple) and k and k[0] == "w"}
    if not sets:
        return row[f"{s}_call"]
    best_ll = max(max(v) for v in sets.values())
    weight = {T: sum(math.exp(x - best_ll) for x in v) for T, v in sets.items()}
    z = sum(weight.values())
    cands = sorted(weight, key=lambda t: -weight[t])[:6]
    cands += [a | b for k, a in enumerate(cands) for b in cands[k + 1:3]]
    best = max(cands, key=lambda R: sum(m * len(R & T) / len(R | T) for T, m in weight.items()) / z)
    ordered = [c for c in calls if c in best] + sorted(best - set(calls))
    return ",".join(ordered)


def main(refdir, swig_tsv, out_tsv, modes):
    refs = load(refdir)
    mode = {s: "keep" for s in "vdj"}
    for item in modes.split(","):
        k, v = item.split("=")
        mode[k] = v
    rows = list(csv.DictReader(open(swig_tsv), delimiter="\t"))
    with open(out_tsv, "w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["sequence_id", "v_call", "d_call", "j_call"])
        for row in rows:
            seq = row["sequence"].upper()
            if row.get("rev_comp") == "T":
                seq = seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
            row["sequence_alignment_full"] = seq
            n = len(seq)
            spans = {s: span(row, s) for s in "vdj"}
            v_end = spans["v"][1] if spans["v"] and row["v_call"] else 0
            j_start = spans["j"][0] if spans["j"] and row["j_call"] else n
            d_lo = spans["d"][0] if spans["d"] and row["d_call"] else None
            d_hi = spans["d"][1] if spans["d"] and row["d_call"] else None
            bounds = {
                "v": (0, d_lo if d_lo is not None else j_start),
                "d": (v_end, j_start),
                "j": (d_hi if d_hi is not None else v_end, n),
            }
            out = [row["sequence_id"]]
            try:
                mu = min(0.3, max(0.005, 1 - float(row["v_identity"])))
            except (KeyError, ValueError):
                mu = 0.05
            for s in "vdj":
                if not refs[s]:
                    out.append(row.get(f"{s}_call", ""))
                elif mode[s] == "endpost":
                    out.append(endpost_segment(row, s, refs[s], bounds[s], mu))
                else:
                    out.append(recall_segment(row, s, refs[s], mode[s], bounds[s]))
            w.writerow(out)


if __name__ == "__main__":
    main(*sys.argv[1:5])
