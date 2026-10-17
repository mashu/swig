//! The D posterior prototype (`../d_posterior.py`, `parametric:1.0` priors)
//! in Rust: fast enough to time, and with switches for the ablation in
//! BENCHMARK_GENAIRR_0.39.6.md.
//!
//!     dposterior D.fasta J.fasta SWIG.tsv OUT.tsv
//!
//! Hypotheses (allele, offset, 5′ trim, 3′ trim) explain the read between
//! swig's V end and J start, ungapped, at the read's divergence
//! μ = clamp(1 − v_identity, 0.005, 0.3). Priors: geometric trims (means 5
//! and 7), negative-binomial N lengths (mean 8, shape 2), genes uniform with
//! each gene's mass split among its alleles. The trims are separable given
//! the D's ends, so each (allele, offset) is bounded in O(length) by a
//! Kadane pass over its match mask and evaluated exactly, in bound order,
//! only while its bound can reach the kept hypotheses. Each kept hypothesis
//! implies GenAIRR's tie set for it; OUT.tsv gets one D call per read.
//!
//! Environment:
//!   TRIM=flat         uniform trims
//!   NP=flat | geom    uniform N lengths | geometric (mean 8)
//!   ALLELE=flat       uniform alleles
//!   ONEOFF=1          each allele only at its offset of best local score
//!                     (the engine's DAG holds one D placement per allele)
//!   DECIDE=jaccard    the set of best expected Jaccard score (default)
//!          mapset     GenAIRR's tie set of the best hypothesis
//!          map        the best hypothesis's allele
//!          marg       the allele of largest posterior marginal
//!   NATS, TOPK        keep hypotheses within NATS of the best (default 12),
//!                     at most TOPK of them (default 200)
//!
//! Timings of the D call and of the positional J re-call's core loop go to
//! stderr.
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;

const MAX_TRIM: usize = 20;

fn env_or<T: std::str::FromStr>(name: &str, default: T) -> T {
    std::env::var(name).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn flag(name: &str, value: &str) -> bool {
    std::env::var(name).map(|v| v == value).unwrap_or(false)
}

fn fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in BufReader::new(std::fs::File::open(path).unwrap()).lines() {
        let line = line.unwrap();
        if let Some(name) = line.strip_prefix('>') {
            out.push((name.split_whitespace().next().unwrap().to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    out
}

/// Negative-binomial log-pmf on 0..=n_max, renormalised.
fn negbin(mean: f64, r: f64, n_max: usize) -> Vec<f64> {
    let p = r / (r + mean);
    let raw: Vec<f64> = (0..=n_max)
        .map(|k| {
            let k = k as f64;
            lgamma(k + r) - lgamma(k + 1.0) - lgamma(r) + r * p.ln() + k * (1.0 - p).ln()
        })
        .collect();
    let z = raw.iter().map(|x| x.exp()).sum::<f64>().ln();
    raw.iter().map(|x| x - z).collect()
}

/// ln Γ(x) by Lanczos (x > 0).
fn lgamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9, 676.520_368_121_885_1, -1_259.139_216_722_402_8,
        771.323_428_777_653_1, -176.615_029_162_140_6, 12.507_343_278_686_905,
        -0.138_571_095_265_720_12, 9.984_369_578_019_572e-6, 1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - lgamma(1.0 - x);
    }
    let x = x - 1.0;
    let t = x + 7.5;
    let a = G.iter().enumerate().skip(1).fold(G[0], |a, (i, g)| a + g / (x + i as f64));
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// `x` low bits set.
fn low(x: usize) -> u64 {
    if x >= 64 { u64::MAX } else { (1u64 << x) - 1 }
}

/// A hypothesis: score (nats), allele, offset, 5′ trim, 3′ trim.
type Hyp = (f64, usize, i64, usize, usize);

#[derive(PartialEq)]
struct Ranked(Hyp);
impl Eq for Ranked {}
impl PartialOrd for Ranked {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Ranked {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.0.0.total_cmp(&o.0.0)
    }
}

struct Model {
    /// pm[k][c]: bit i set when allele k has base c at index i.
    pm: Vec<[u64; 4]>,
    names: Vec<String>,
    lens: Vec<usize>,
    maxlen: usize,
    t5: Vec<f64>,
    t3: Vec<f64>,
    np: Vec<f64>,
    allele: Vec<f64>,
}

impl Model {
    fn new(d: Vec<(String, Vec<u8>)>) -> Self {
        assert!(d.len() <= 64, "the tie-set bitmask holds 64 alleles");
        let maxlen = d.iter().map(|(_, s)| s.len()).max().unwrap();
        assert!(maxlen <= 63, "D alleles fit a 64-bit mask");
        let gene = |n: &str| n.split('*').next().unwrap().to_string();
        let mut per_gene: HashMap<String, usize> = HashMap::new();
        for (n, _) in &d {
            *per_gene.entry(gene(n)).or_default() += 1;
        }
        let allele = d
            .iter()
            .map(|(n, _)| if flag("ALLELE", "flat") { 0.0 } else { -(per_gene[&gene(n)] as f64).ln() })
            .collect();
        let pm = d
            .iter()
            .map(|(_, s)| {
                let mut m = [0u64; 4];
                for (i, &b) in s.iter().enumerate() {
                    if let Some(c) = base_index(b) {
                        m[c] |= 1 << i;
                    }
                }
                m
            })
            .collect();
        let trim = |mean| if flag("TRIM", "flat") { vec![0.0; MAX_TRIM + 1] } else { negbin(mean, 1.0, MAX_TRIM) };
        let np = if flag("NP", "flat") {
            vec![0.0; 201]
        } else if flag("NP", "geom") {
            negbin(8.0, 1.0, 200)
        } else {
            negbin(8.0, 2.0, 200)
        };
        Self {
            pm,
            names: d.iter().map(|(n, _)| n.clone()).collect(),
            lens: d.iter().map(|(_, s)| s.len()).collect(),
            maxlen,
            t5: trim(5.0),
            t3: trim(7.0),
            np,
            allele,
        }
    }

    /// The D call (a set of alleles as a bitmask) for the junction
    /// `[v_end, j_start)` of `read`.
    fn call(&self, read: &[u8], v_end: usize, j_start: usize, mu: f64) -> Option<u64> {
        if j_start <= v_end {
            return None;
        }
        let n = self.lens.len();
        let jlen = j_start - v_end;
        let (wm, wx) = (((1.0 - mu) / 0.25).ln(), ((mu / 3.0) / 0.25).ln());
        let np: Vec<f64> = (0..=jlen).map(|k| self.np[k.min(200)]).collect();
        let max = |v: &[f64]| v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let (np_max, t5_max, t3_max) = (max(&np), max(&self.t5), max(&self.t3));
        let lo = v_end as i64 - self.maxlen as i64 + 1;
        // Germline indices of an allele at offset o that fall in the junction.
        let window = |o: i64| ((v_end as i64 - o).max(0) as usize, ((j_start as i64 - o) as usize).min(self.maxlen));
        // Read base masks: bit p of rm[c] is set when read[p] is base c.
        let words = read.len() / 64 + 2;
        let mut rm = [vec![0u64; words], vec![0u64; words], vec![0u64; words], vec![0u64; words]];
        for (p, &b) in read.iter().enumerate() {
            if let Some(c) = base_index(b) {
                rm[c][p / 64] |= 1 << (p % 64);
            }
        }
        let at = |v: &Vec<u64>, o: usize| -> u64 {
            let (w, r) = (o / 64, o % 64);
            if r == 0 { v[w] } else { (v[w] >> r) | (v[w + 1] << (64 - r)) }
        };
        // Match mask of allele k at offset o within its window [i0, min(l, i1)).
        let match_mask = |k: usize, o: i64, i0: usize, i1: usize| -> (u64, u64) {
            let win = low(self.lens[k].min(i1)) & !low(i0);
            let o = o as usize;
            let m = (0..4).fold(0u64, |m, c| m | (self.pm[k][c] & at(&rm[c], o)));
            (m & win, win)
        };
        // cs[x]: the evidence of germline indices < x.
        let fill = |k: usize, o: i64, i0: usize, i1: usize, cs: &mut [f64]| {
            let (m, win) = match_mask(k, o, i0, i1);
            for (x, c) in cs.iter_mut().enumerate().take(self.lens[k] + 1) {
                let hits = (m & low(x)).count_ones() as f64;
                let cols = (win & low(x)).count_ones() as f64;
                *c = wm * hits + wx * (cols - hits);
            }
        };
        // Best contiguous evidence within the window (Kadane over the runs
        // of matches and mismatches of the match mask): a tight upper bound.
        let bound_of = |k: usize, o: i64, i0: usize, i1: usize| -> f64 {
            let (m, win) = match_mask(k, o, i0, i1);
            if m == 0 {
                return wx;
            }
            let (mut pos, end) = (win.trailing_zeros(), 64 - win.leading_zeros());
            let (mut run, mut best) = (0.0f64, f64::NEG_INFINITY);
            while pos < end {
                let rest = m >> pos;
                let (len, w) = if rest & 1 == 1 {
                    (rest.trailing_ones().min(end - pos), wm)
                } else {
                    (rest.trailing_zeros().min(end - pos), wx)
                };
                let v = len as f64 * w;
                run = if w > 0.0 && run <= 0.0 { v } else { run + v };
                best = best.max(run);
                pos += len;
            }
            best
        };
        let mut cs = vec![0.0f64; self.maxlen + 1];
        let mut a_buf = [f64::NEG_INFINITY; MAX_TRIM + 1];
        // The best hypothesis of (k, o) over its trims; with `floor`, also
        // every hypothesis above it, pushed to `out`.
        let mut exact = |k: usize, o: i64, i0: usize, i1: usize, cs: &[f64], floor: f64, out: &mut Vec<Hyp>| -> f64 {
            let l = self.lens[k];
            let e_hi = l.min(i1);
            let e_lo = l.saturating_sub(MAX_TRIM).max(i0 + 1);
            let t5_hi = MAX_TRIM.min(e_hi.saturating_sub(1));
            if i0 > t5_hi {
                return f64::NEG_INFINITY;
            }
            for t5 in i0..=t5_hi {
                a_buf[t5] = -cs[t5] + self.t5[t5] + np[(o + t5 as i64 - v_end as i64).clamp(0, jlen as i64) as usize];
            }
            // 5′ trims by score, best first.
            let mut order = [0usize; MAX_TRIM + 1];
            let n_t5 = t5_hi + 1 - i0;
            for (j, t5) in (i0..=t5_hi).enumerate() {
                order[j] = t5;
            }
            order[..n_t5].sort_unstable_by(|a, b| a_buf[*b].partial_cmp(&a_buf[*a]).unwrap());
            let (mut run, mut next, mut m) = (f64::NEG_INFINITY, i0, f64::NEG_INFINITY);
            for e in e_lo..=e_hi {
                while next < e && next <= t5_hi {
                    run = run.max(a_buf[next]);
                    next += 1;
                }
                if run == f64::NEG_INFINITY {
                    continue;
                }
                let be = cs[e] + self.t3[l - e] + np[(j_start as i64 - (o + e as i64)).clamp(0, jlen as i64) as usize] + self.allele[k];
                m = m.max(run + be);
                for &t5 in &order[..n_t5] {
                    let sc = a_buf[t5] + be;
                    if sc <= floor {
                        break;
                    }
                    if t5 < e {
                        out.push((sc, k, o, t5, l - e));
                    }
                }
            }
            m
        };
        // Bound every (allele, offset).
        let mut bounds: Vec<(f64, usize, i64)> = Vec::new();
        let mut best_o: Vec<Option<(f64, i64)>> = vec![None; n];
        for o in lo..j_start as i64 {
            let (i0, i1) = window(o);
            if i1 <= i0 {
                continue;
            }
            for k in 0..n {
                if self.lens[k] <= i0 {
                    continue;
                }
                let b = bound_of(k, o, i0, i1);
                if best_o[k].is_none_or(|(v, _)| b > v) {
                    best_o[k] = Some((b, o));
                }
                bounds.push((b + t5_max + t3_max + 2.0 * np_max + self.allele[k], k, o));
            }
        }
        if flag("ONEOFF", "1") {
            bounds.retain(|&(_, k, o)| best_o[k].is_some_and(|(_, bo)| bo == o));
        }
        bounds.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        // Evaluate in bound order, keeping the TOPK best hypotheses in a
        // min-heap; stop at the first bound that cannot beat the heap's
        // minimum (once full) or the best minus NATS.
        let (nats, topk) = (env_or("NATS", 12.0), env_or("TOPK", 200usize));
        let mut heap: BinaryHeap<Reverse<Ranked>> = BinaryHeap::with_capacity(topk + 1);
        let mut best = f64::NEG_INFINITY;
        let mut batch = Vec::new();
        for &(bound, k, o) in &bounds {
            let heap_min = if heap.len() >= topk { heap.peek().unwrap().0.0.0 } else { f64::NEG_INFINITY };
            let floor = (best - nats).max(heap_min);
            if bound <= floor {
                break;
            }
            let (i0, i1) = window(o);
            fill(k, o, i0, i1, &mut cs);
            batch.clear();
            best = best.max(exact(k, o, i0, i1, &cs, floor, &mut batch));
            for &h in &batch {
                if heap.len() < topk {
                    heap.push(Reverse(Ranked(h)));
                } else if h.0 > heap.peek().unwrap().0.0.0 {
                    heap.pop();
                    heap.push(Reverse(Ranked(h)));
                }
            }
        }
        if best == f64::NEG_INFINITY {
            return None;
        }
        let mut hyps: Vec<Hyp> = heap.into_iter().map(|Reverse(Ranked(h))| h).filter(|h| h.0 > best - nats).collect();
        hyps.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        if flag("DECIDE", "map") {
            return Some(1u64 << hyps[0].1);
        }
        if flag("DECIDE", "marg") {
            let mut m = vec![0.0f64; n];
            for &(sc, k, ..) in &hyps {
                m[k] += (sc - best).exp();
            }
            let k = (0..n).max_by(|&a, &b| m[a].partial_cmp(&m[b]).unwrap()).unwrap();
            return Some(1u64 << k);
        }
        // GenAIRR's tie set of each hypothesis, from per-offset masks
        // (mk[i], bit b: allele b matches the read at germline index i).
        let mut masks: HashMap<i64, Vec<u64>> = HashMap::new();
        let mut cache: HashMap<(i64, usize, usize, usize), u64> = HashMap::new();
        let mut mass: Vec<(u64, f64)> = Vec::new();
        for &(sc, k, o, t5, t3) in &hyps {
            let (i0, i1) = window(o);
            let mk = masks.entry(o).or_insert_with(|| {
                let per: Vec<u64> = (0..n).map(|b| if self.lens[b] > i0 { match_mask(b, o, i0, i1).0 } else { 0 }).collect();
                (0..self.maxlen).map(|i| (0..n).fold(0u64, |m, b| m | ((per[b] >> i) & 1) << b)).collect()
            });
            let (s0, e0) = (t5, self.lens[k] - t3);
            let t = *cache.entry((o, s0, e0, t3)).or_insert_with(|| {
                let mut scores = [0u32; 64];
                for &m in &mk[s0..e0] {
                    let mut m = m;
                    while m != 0 {
                        scores[m.trailing_zeros() as usize] += 1;
                        m &= m - 1;
                    }
                }
                // The narrowing walk into the N bases: a base counts while it
                // splits the current tie.
                let step = |i: usize, scores: &mut [u32; 64]| -> bool {
                    let top = *scores[..n].iter().max().unwrap();
                    let tied = (0..n).filter(|&b| scores[b] == top).fold(0u64, |m, b| m | 1 << b);
                    let hit = tied & mk[i];
                    if hit == 0 || hit == tied {
                        return false;
                    }
                    let mut h = hit;
                    while h != 0 {
                        scores[h.trailing_zeros() as usize] += 1;
                        h &= h - 1;
                    }
                    true
                };
                for i in (i0.max(s0.saturating_sub(t5))..s0).rev() {
                    if !step(i, &mut scores) {
                        break;
                    }
                }
                for i in e0..i1.min(e0 + t3) {
                    if !step(i, &mut scores) {
                        break;
                    }
                }
                let top = *scores[..n].iter().max().unwrap();
                if top == 0 {
                    low(n)
                } else {
                    (0..n).filter(|&b| scores[b] == top).fold(0u64, |m, b| m | 1 << b)
                }
            });
            let w = (sc - best).exp();
            match mass.iter_mut().find(|(s, _)| *s == t) {
                Some(entry) => entry.1 += w,
                None => mass.push((t, w)),
            }
        }
        if flag("DECIDE", "mapset") {
            return Some(mass[0].0);
        }
        // The set of best expected Jaccard score among the 8 likeliest tie
        // sets and the pairwise unions of the first 4.
        let z: f64 = mass.iter().map(|m| m.1).sum();
        let mut sets = mass.clone();
        sets.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        let mut cands: Vec<u64> = sets.iter().take(8).map(|s| s.0).collect();
        let base = cands.clone();
        for (i, a) in base.iter().enumerate() {
            for b in base.iter().take(4).skip(i + 1) {
                cands.push(a | b);
            }
        }
        let expected = |r: u64| mass.iter().map(|&(t, m)| m * (r & t).count_ones() as f64 / (r | t).count_ones() as f64).sum::<f64>() / z;
        cands.into_iter().map(|r| (expected(r), r)).fold(None, |acc: Option<(f64, u64)>, (e, r)| match acc {
            Some((be, _)) if be >= e => acc,
            _ => Some((e, r)),
        }).map(|(_, r)| r)
    }

    fn names_of(&self, set: u64) -> String {
        let mut v: Vec<&str> = (0..self.names.len()).filter(|&b| set >> b & 1 == 1).map(|b| self.names[b].as_str()).collect();
        v.sort();
        v.join(",")
    }
}

/// The positional J re-call's core: every J allele scored at the called
/// allele's germline index over the called span (indels not handled).
fn positional_j(j: &[(String, Vec<u8>)], read: &[u8], q0: usize, q1: usize, g0: usize) -> u64 {
    let mut best = 0;
    let mut set = 0u64;
    for (b, (_, g)) in j.iter().enumerate() {
        let n = (q0..q1).filter(|&p| p + g0 - q0 < g.len() && read[p] == g[p + g0 - q0]).count();
        if n > best {
            best = n;
            set = 1 << b;
        } else if n == best {
            set |= 1 << b;
        }
    }
    set
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("usage: dposterior D.fasta J.fasta SWIG.tsv OUT.tsv");
        std::process::exit(2);
    }
    let model = Model::new(fasta(&args[1]));
    let jset = fasta(&args[2]);
    let mut lines = BufReader::new(std::fs::File::open(&args[3]).unwrap()).lines();
    let header: Vec<String> = lines.next().unwrap().unwrap().split('\t').map(String::from).collect();
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    let (c_id, c_seq, c_rc, c_ve, c_js, c_vi) =
        (col("sequence_id"), col("sequence"), col("rev_comp"), col("v_sequence_end"), col("j_sequence_start"), col("v_identity"));
    let (c_jge, c_jgs) = (col("j_sequence_end"), col("j_germline_start"));
    let rows: Vec<Vec<String>> = lines.map(|l| l.unwrap().split('\t').map(String::from).collect()).collect();
    let complement = |b: u8| match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    };
    let reads: Vec<Vec<u8>> = rows
        .iter()
        .map(|r| {
            let s: Vec<u8> = r[c_seq].bytes().map(|b| b.to_ascii_uppercase()).collect();
            if r[c_rc] == "T" { s.iter().rev().map(|&b| complement(b)).collect() } else { s }
        })
        .collect();
    let t = Instant::now();
    let calls: Vec<Option<String>> = rows
        .iter()
        .zip(&reads)
        .map(|(r, read)| match (r[c_ve].parse::<usize>(), r[c_js].parse::<usize>()) {
            (Ok(ve), Ok(js)) => {
                let mu = r[c_vi].parse::<f64>().map(|v| (1.0 - v).clamp(0.005, 0.3)).unwrap_or(0.05);
                model.call(read, ve, js - 1, mu).map(|s| model.names_of(s))
            }
            _ => None,
        })
        .collect();
    let d_time = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let mut acc = 0u64;
    for (r, read) in rows.iter().zip(&reads) {
        if let (Ok(q0), Ok(q1), Ok(g0)) = (r[c_js].parse::<usize>(), r[c_jge].parse::<usize>(), r[c_jgs].parse::<usize>()) {
            acc ^= positional_j(&jset, read, q0 - 1, q1, g0 - 1);
        }
    }
    let j_time = t.elapsed().as_secs_f64();
    let per_read = |s: f64| s * 1e6 / rows.len().max(1) as f64;
    eprintln!(
        "{} reads: D posterior {:.3} s ({:.1} us/read), positional J {:.4} s ({:.2} us/read) [{acc:x}]",
        rows.len(), d_time, per_read(d_time), j_time, per_read(j_time)
    );
    let mut w = std::io::BufWriter::new(std::fs::File::create(&args[4]).unwrap());
    writeln!(w, "sequence_id\td_call").unwrap();
    for (r, c) in rows.iter().zip(calls) {
        writeln!(w, "{}\t{}", r[c_id], c.unwrap_or_default()).unwrap();
    }
}
