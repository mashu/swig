//! Swig's double-D screen: VDDJ evidence beside the ordinary call.
//!
//! Screening leaves every AIRR row unchanged. A read whose call places V and
//! J on a locus with D genes (IGH, TRB, TRD) is searched again between the
//! called V's end and J's start for two D segments, and reported in a
//! sparse sidecar table when the pair passes the screen's evidence rules
//! (`public/methods/04_DOUBLE_D_SCREEN.md`):
//!
//! - **span**: in `long-span` mode, the V–J span is at least `min_vj_span`
//!   bases;
//! - **seeds**: each D hit holds an exact `seed_length`-word of a D allele
//!   of the locus, or a `max(6, seed_length − 3)` rescue word whose
//!   ungapped extension (match +2, mismatch −2, best prefix either way)
//!   reaches `seed_length` aligned bases; hits on one allele diagonal count
//!   once;
//! - **gain**: two ordered, non-overlapping hits score at least
//!   `min_score_gain` above the best single hit;
//! - **pseudo-tandem**: the read from D1's start to D2's end, `pseudo_trim`
//!   bases shorter, is more than `pseudo_mismatches` Hamming mismatches from
//!   every substring of a D allele of the locus (the IgScout Δ-distance
//!   test), so no single D explains it.
//!
//! Pairs are tried in evidence order (pair score, gain, matched bases,
//! shorter NP2, names); pairs within two score units of the reported one
//! are listed as alternatives. A hit's call lists every allele with the same
//! score and coordinates.

use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};

use crate::alphabet::{N, decode};
use crate::germline::{AlleleId, GermlineSet};
use crate::locus::Loci;
use crate::report::{AirrRecord, Row, fraction, push_run};
use crate::scoring::Score;

/// When the screen runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DoubleDMode {
    Off,
    /// Every read whose call places V and J.
    All,
    /// Only reads whose V–J span is at least `min_vj_span`.
    LongSpan,
}

impl DoubleDMode {
    pub fn name(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::All => "all",
            Self::LongSpan => "long_span",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "off" => Some(Self::Off),
            "all" => Some(Self::All),
            "long-span" | "long_span" => Some(Self::LongSpan),
            _ => None,
        }
    }
}

/// The screen's settings (see the module documentation).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DoubleDOptions {
    pub mode: DoubleDMode,
    pub min_vj_span: u32,
    /// One of [`SEED_LENGTHS`].
    pub seed_length: u32,
    pub pseudo_trim: u32,
    pub pseudo_mismatches: u32,
    pub min_score_gain: Score,
}

impl Default for DoubleDOptions {
    fn default() -> Self {
        Self {
            mode: DoubleDMode::Off,
            min_vj_span: 40,
            seed_length: 11,
            pseudo_trim: 5,
            pseudo_mismatches: 3,
            min_score_gain: 8,
        }
    }
}

/// The seed lengths the screen indexes.
pub const SEED_LENGTHS: std::ops::RangeInclusive<u32> = 6..=24;

/// Alternatives reported beside the selected pair.
const MAX_ALTERNATIVES: usize = 12;

/// The sidecar's columns. `swig_batch_record_index` is the read's 0-based
/// index in the reads of one [`crate::run::Run::assign`] call (in the
/// WebAssembly ABI, the batch), and in the whole input in
/// [`crate::run::Run::write_airr`] (the command line's `--double-d-out`).
pub const COLUMNS: [&str; 38] = [
    "swig_batch_record_index",
    "sequence_id",
    "locus",
    "v_call",
    "standard_d_call",
    "d_call",
    "d2_call",
    "j_call",
    "d_score",
    "d_identity",
    "d_cigar",
    "d2_score",
    "d2_identity",
    "d2_cigar",
    "d_sequence_start",
    "d_sequence_end",
    "d_germline_start",
    "d_germline_end",
    "d2_sequence_start",
    "d2_sequence_end",
    "d2_germline_start",
    "d2_germline_end",
    "d_sequence_alignment",
    "d_germline_alignment",
    "d2_sequence_alignment",
    "d2_germline_alignment",
    "np2",
    "np3",
    "np2_length",
    "np3_length",
    "swig_double_d_vj_span",
    "swig_double_d_seed_length",
    "swig_double_d_pair_score",
    "swig_double_d_best_single_score",
    "swig_double_d_score_gain",
    "swig_double_d_pseudo_distance",
    "swig_double_d_mode",
    "swig_double_d_alternatives",
];

/// The sidecar's header line, without a line break.
pub fn header() -> String {
    COLUMNS.join("\t")
}

/// A word of A/C/G/T codes as 2-bit digits; `None` if it holds `N`.
fn word(codes: &[u8]) -> Option<u64> {
    codes.iter().try_fold(0u64, |packed, &c| {
        (c < N).then_some((packed << 2) | u64::from(c))
    })
}

type WordIndex = HashMap<u64, Vec<(AlleleId, u32)>>;

/// Every `len`-word of each allele: the allele and the word's offset.
fn index(alleles: &[(AlleleId, &[u8])], len: usize) -> WordIndex {
    let mut words = WordIndex::new();
    for &(allele, codes) in alleles {
        for (offset, window) in codes.windows(len).enumerate() {
            if let Some(w) = word(window) {
                words.entry(w).or_default().push((allele, offset as u32));
            }
        }
    }
    words
}

/// An ungapped D hit; query coordinates are in the whole oriented read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Hit {
    allele: AlleleId,
    q_lo: u32,
    q_hi: u32,
    r_lo: u32,
    r_hi: u32,
    matches: u32,
    score: Score,
}

impl Hit {
    fn len(&self) -> u32 {
        self.q_hi - self.q_lo
    }

    fn same_alignment(&self, other: &Hit) -> bool {
        (self.score, self.q_lo, self.q_hi, self.r_lo, self.r_hi)
            == (other.score, other.q_lo, other.q_hi, other.r_lo, other.r_hi)
    }
}

fn same(a: u8, b: u8) -> bool {
    a < N && a == b
}

/// Steps of the best-scoring prefix of `pairs` (+2 a match, −2 otherwise).
fn best_prefix(pairs: impl Iterator<Item = (u8, u8)>) -> usize {
    let (mut running, mut best, mut steps) = (0, 0, 0);
    for (step, (a, b)) in pairs.enumerate() {
        running += if same(a, b) { 2 } else { -2 };
        if running > best {
            (best, steps) = (running, step + 1);
        }
    }
    steps
}

/// The seed `query[q..q + len]` = `reference[r..r + len]` extended without
/// gaps to its best-scoring prefix on each side, within `query`; `offset`
/// moves query coordinates to the whole read.
fn extend(
    allele: AlleleId,
    query: &[u8],
    reference: &[u8],
    (q, r, len): (usize, usize, usize),
    offset: u32,
) -> Hit {
    let left = best_prefix(
        query[..q]
            .iter()
            .rev()
            .copied()
            .zip(reference[..r].iter().rev().copied()),
    );
    let right = best_prefix(
        query[q + len..]
            .iter()
            .copied()
            .zip(reference[r + len..].iter().copied()),
    );
    let (q_lo, q_hi, r_lo) = (q - left, q + len + right, r - left);
    let matches = query[q_lo..q_hi]
        .iter()
        .zip(&reference[r_lo..])
        .filter(|&(&a, &b)| same(a, b))
        .count() as u32;
    let columns = (q_hi - q_lo) as u32;
    Hit {
        allele,
        q_lo: q_lo as u32 + offset,
        q_hi: q_hi as u32 + offset,
        r_lo: r_lo as u32,
        r_hi: (r + len + right) as u32,
        matches,
        score: 2 * matches as Score - 2 * (columns - matches) as Score,
    }
}

/// The CIGAR of an ungapped hit, as the main rows write theirs: the read's
/// bases before it soft-clipped, the germline's skipped, `=`/`X` columns,
/// then the read's bases after it.
fn cigar(hit: &Hit, query: &[u8], reference: &[u8]) -> String {
    let mut text = String::new();
    push_run(&mut text, 'S', hit.q_lo);
    push_run(&mut text, 'N', hit.r_lo);
    let mut run: Option<(char, u32)> = None;
    let columns = query[hit.q_lo as usize..hit.q_hi as usize]
        .iter()
        .zip(&reference[hit.r_lo as usize..]);
    for (&a, &b) in columns {
        let symbol = if same(a, b) { '=' } else { 'X' };
        run = match run {
            Some((s, n)) if s == symbol => Some((s, n + 1)),
            Some((s, n)) => {
                push_run(&mut text, s, n);
                Some((symbol, 1))
            }
            None => Some((symbol, 1)),
        };
    }
    if let Some((s, n)) = run {
        push_run(&mut text, s, n);
    }
    push_run(&mut text, 'S', query.len() as u32 - hit.q_hi);
    text
}

/// Fewest mismatches between `span` shortened by `trim` bases, at any
/// offset, and any equally long substring of an allele in `alleles`;
/// `None` when the shortened span is shorter than `min_len` or longer than
/// every allele.
fn delta_distance<'a>(
    span: &[u8],
    alleles: impl Iterator<Item = &'a [u8]> + Clone,
    trim: usize,
    min_len: usize,
) -> Option<u32> {
    let len = span
        .len()
        .checked_sub(trim)
        .filter(|&l| l > 0 && l >= min_len)?;
    let mut best: Option<u32> = None;
    for start in 0..=trim {
        let window = &span[start..start + len];
        for allele in alleles.clone() {
            for germline in allele.windows(len) {
                let bound = best.unwrap_or(u32::MAX);
                let mut distance = 0;
                for (a, b) in window.iter().zip(germline) {
                    distance += u32::from(a != b);
                    if distance >= bound {
                        break;
                    }
                }
                best = Some(distance.min(bound));
            }
        }
    }
    best
}

/// Two hits (indices into the read's hits) and their evidence.
#[derive(Clone, Copy, Debug)]
struct Pair {
    first: usize,
    second: usize,
    score: Score,
    gain: Score,
    pseudo: Option<u32>,
}

/// The screen over one D database: its seed indices, built once.
pub struct DoubleDScreen {
    options: DoubleDOptions,
    seeds: WordIndex,
    /// Words `rescue_length` long; empty when that is not shorter.
    rescue: WordIndex,
    rescue_length: usize,
}

impl DoubleDScreen {
    /// Index the alleles of `d` that are at least `seed_length` long.
    pub fn new(d: &GermlineSet, options: DoubleDOptions) -> Self {
        let seed = options.seed_length as usize;
        let rescue_length = seed.saturating_sub(3).max(6);
        let alleles: Vec<(AlleleId, &[u8])> = (0..d.len() as AlleleId)
            .map(|id| (id, d.allele(id).codes.as_slice()))
            .filter(|(_, codes)| codes.len() >= seed)
            .collect();
        Self {
            options,
            seeds: index(&alleles, seed),
            rescue: if rescue_length < seed {
                index(&alleles, rescue_length)
            } else {
                WordIndex::new()
            },
            rescue_length,
        }
    }

    /// The sidecar fields, from `locus` on (without
    /// `swig_batch_record_index` and `sequence_id`), of the read `query`
    /// (codes, oriented as called) whose main row is `record`, with the D
    /// alleles of `d_set`, if the screen finds a VDDJ pair. The line has no
    /// line break.
    pub fn screen(&self, d_set: &GermlineSet, query: &[u8], record: &AirrRecord) -> Option<String> {
        let options = &self.options;
        if options.mode == DoubleDMode::Off {
            return None;
        }
        let [v, d, j] = [&record.v, &record.d, &record.j].map(Option::as_ref);
        let locus = record.locus;
        if locus.is_some_and(|l| !l.has_d()) {
            return None;
        }
        let loci = locus.map_or(Loci::ANY, Loci::of);
        let (v, j) = (v?, j?);
        // The called V's end and J's start bound the window.
        let (start, end) = (v.sequence_end, j.sequence_start.checked_sub(1)?);
        if end <= start || end as usize > query.len() {
            return None;
        }
        let vj_span = end - start;
        if options.mode == DoubleDMode::LongSpan && vj_span < options.min_vj_span {
            return None;
        }
        let seed = options.seed_length as usize;
        let window = &query[start as usize..end as usize];
        if window.len() < 2 * seed {
            return None;
        }

        let mut hits: Vec<Hit> = Vec::new();
        let mut diagonals: HashSet<(AlleleId, i64)> = HashSet::new();
        for (words, len, rescue) in [
            (&self.seeds, seed, false),
            (&self.rescue, self.rescue_length, true),
        ] {
            if words.is_empty() {
                continue;
            }
            for (q, w) in window.windows(len).enumerate() {
                let Some(found) = word(w).and_then(|w| words.get(&w)) else {
                    continue;
                };
                for &(allele, r) in found {
                    let diagonal = (allele, q as i64 - i64::from(r));
                    if !d_set.allele(allele).loci.pairs_with(loci) || diagonals.contains(&diagonal)
                    {
                        continue;
                    }
                    let codes = &d_set.allele(allele).codes;
                    let hit = extend(allele, window, codes, (q, r as usize, len), start);
                    // A rescue word only shortens the seed, not the evidence.
                    if rescue && (hit.len() as usize) < seed {
                        continue;
                    }
                    diagonals.insert(diagonal);
                    hits.push(hit);
                }
            }
        }
        if hits.len() < 2 {
            return None;
        }
        let name = |hit: &Hit| d_set.allele(hit.allele).name.as_str();
        hits.sort_by(|a, b| {
            (a.q_lo, Reverse(a.q_hi), Reverse(a.score), name(a)).cmp(&(
                b.q_lo,
                Reverse(b.q_hi),
                Reverse(b.score),
                name(b),
            ))
        });
        let best_single = hits.iter().map(|h| h.score).max().unwrap_or(0);

        let mut candidates: Vec<Pair> = Vec::new();
        for (first, d1) in hits.iter().enumerate() {
            for (second, d2) in hits.iter().enumerate().skip(first + 1) {
                let score = d1.score + d2.score;
                let gain = score - best_single;
                if d1.q_hi <= d2.q_lo && gain >= options.min_score_gain {
                    candidates.push(Pair {
                        first,
                        second,
                        score,
                        gain,
                        pseudo: None,
                    });
                }
            }
        }
        let order = |p: &Pair| {
            let (d1, d2) = (&hits[p.first], &hits[p.second]);
            (
                Reverse(p.score),
                Reverse(p.gain),
                Reverse(d1.matches + d2.matches),
                d2.q_lo - d1.q_hi,
                (name(d1), name(d2)),
            )
        };
        candidates.sort_by(|a, b| order(a).cmp(&order(b)));

        // Δ-distance is the costly test: try pairs in evidence order only as
        // far as an alternative could reach, once per outer span.
        let alleles = d_set
            .alleles()
            .iter()
            .filter(|a| a.loci.pairs_with(loci))
            .map(|a| a.codes.as_slice());
        let mut pseudo_by_span: HashMap<(u32, u32), Option<u32>> = HashMap::new();
        let mut pairs: Vec<Pair> = Vec::new();
        for mut candidate in candidates {
            if pairs
                .first()
                .is_some_and(|best| candidate.score + 2 < best.score)
            {
                break;
            }
            let (d1, d2) = (&hits[candidate.first], &hits[candidate.second]);
            candidate.pseudo = *pseudo_by_span.entry((d1.q_lo, d2.q_hi)).or_insert_with(|| {
                delta_distance(
                    &query[d1.q_lo as usize..d2.q_hi as usize],
                    alleles.clone(),
                    options.pseudo_trim as usize,
                    seed,
                )
            });
            if candidate
                .pseudo
                .is_none_or(|p| p > options.pseudo_mismatches)
            {
                pairs.push(candidate);
            }
        }
        let selected = *pairs.first()?;

        // Every allele aligned exactly as `hit`, by name.
        let calls = |hit: &Hit| {
            let mut names: Vec<&str> = hits
                .iter()
                .filter(|h| h.same_alignment(hit))
                .flat_map(|h| d_set.allele(h.allele).names())
                .collect();
            names.sort_unstable();
            names.dedup();
            names.join(",")
        };
        let (d1, d2) = (&hits[selected.first], &hits[selected.second]);
        // Read bases as the main row quotes them, IUPAC codes included.
        let text = record.oriented();
        let between = |from: u32, to: u32| &text[from as usize..to.max(from) as usize];
        let np2 = between(d1.q_hi, d2.q_lo);
        let np3 = between(d2.q_hi, end);
        let alternatives: Vec<String> = pairs[1..]
            .iter()
            .take(MAX_ALTERNATIVES)
            .map(|p| {
                format!(
                    "{}+{}|{}|{}|{}",
                    calls(&hits[p.first]),
                    calls(&hits[p.second]),
                    p.score,
                    p.gain,
                    p.pseudo.map_or(-1, i64::from)
                )
            })
            .collect();

        let mut row = Row::with_capacity(1024);
        row.field(locus.map_or("", |l| l.name()));
        row.field(&v.call);
        row.field(d.map_or("", |d| &d.call));
        row.field(&calls(d1));
        row.field(&calls(d2));
        row.field(&j.call);
        for hit in [d1, d2] {
            let reference = &d_set.allele(hit.allele).codes;
            row.number(Some(hit.score));
            row.field(&fraction(f64::from(hit.matches) / f64::from(hit.len())));
            row.field(&cigar(hit, query, reference));
        }
        for hit in [d1, d2] {
            for value in [hit.q_lo + 1, hit.q_hi, hit.r_lo + 1, hit.r_hi] {
                row.number(Some(value));
            }
        }
        for hit in [d1, d2] {
            let reference = &d_set.allele(hit.allele).codes;
            row.field(&text[hit.q_lo as usize..hit.q_hi as usize]);
            row.field(&decode(&reference[hit.r_lo as usize..hit.r_hi as usize]));
        }
        row.field(np2);
        row.field(np3);
        row.number(Some(np2.len()));
        row.number(Some(np3.len()));
        for value in [
            Score::try_from(vj_span).unwrap_or(Score::MAX),
            Score::try_from(options.seed_length).unwrap_or(Score::MAX),
            selected.score,
            best_single,
            selected.gain,
        ] {
            row.number(Some(value));
        }
        row.number(selected.pseudo);
        row.field(options.mode.name());
        row.field(&alternatives.join(";"));
        Some(row.line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::encode;

    #[test]
    fn delta_distance_is_the_fewest_mismatches_of_a_trimmed_span() {
        let d1 = b"ACGTACGTACGTAAAA".as_slice();
        let d2 = b"TTTTGGGGCCCC".as_slice();
        // One allele's substring, one base changed, two foreign bases at
        // the end: trimming 2 aligns the rest at one mismatch.
        let span = b"ACGTACCTACGTGG";
        assert_eq!(delta_distance(span, [d1, d2].into_iter(), 2, 8), Some(1));
        assert_eq!(delta_distance(span, [d1, d2].into_iter(), 0, 8), Some(3));
        // Too short once trimmed, or longer than every allele.
        assert_eq!(delta_distance(b"ACGTAC", [d1].into_iter(), 2, 8), None);
        assert_eq!(delta_distance(&[0; 20], [d2].into_iter(), 2, 8), None);
    }

    #[test]
    fn seeds_extend_to_their_best_prefix_on_each_side() {
        let reference = encode(b"GGGGACGTACGTACGTCCCC");
        // The seed ACGTACGT at query 3 = reference 4, one more matching
        // base on each side, then mismatches only.
        let query = encode(b"TTGACGTACGTATTTT");
        let hit = extend(7, &query, &reference, (3, 4, 8), 100);
        assert_eq!((hit.q_lo, hit.q_hi, hit.r_lo, hit.r_hi), (102, 112, 3, 13));
        assert_eq!((hit.matches, hit.score), (10, 20));
        let local = extend(7, &query, &reference, (3, 4, 8), 0);
        assert_eq!(cigar(&local, &query, &reference), "2S3N10=4S");
    }

    #[test]
    fn modes_read_their_names() {
        for mode in [DoubleDMode::Off, DoubleDMode::All, DoubleDMode::LongSpan] {
            assert_eq!(DoubleDMode::from_name(mode.name()), Some(mode));
        }
        assert_eq!(
            DoubleDMode::from_name("long-span"),
            Some(DoubleDMode::LongSpan)
        );
        assert_eq!(header().split('\t').count(), COLUMNS.len());
    }
}
