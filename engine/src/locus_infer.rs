//! Locus from sequence: a V, J or C allele whose name, `locus=` key and
//! header name no locus takes the loci of its closest relatives of named
//! locus ([`LocusSource::Relative`]).
//!
//! The relatives are first the database's own alleles of the same kind
//! whose locus is named ([`LocusSource::Named`], [`LocusSource::Header`]),
//! then, when those do not decide, an embedded panel of human and mouse
//! family representatives (`data/locus_panel.fa`, 122 V, 118 J and 29 C
//! alleles taken from the shipped IMGT pack by
//! `scripts/build-locus-data.mjs`), so that a database whose V alleles are
//! all custom-named still gets its loci.
//!
//! Each allele is compared with its [`CANDIDATES`] relatives sharing most
//! distinct 5-mers by a local alignment (`AffineScoring(2, −3, −5, −1)`,
//! as annotation transfer uses) in a band of ±[`HALF_BAND`] diagonals
//! around the diagonal most of their shared 5-mers lie on. The allele takes
//! the loci of every relative scoring within [`MARGIN`] of the best, when
//! the best scores at least [`MIN_SCORE`] and at least a fraction of the
//! allele's score against itself: [`MIN_SELF_FRACTION`] for the database's
//! own alleles, [`MIN_PANEL_FRACTION`] for the panel's. Several loci may
//! result: a relative of the shared TRA/TRD V genes gives both. In the
//! design's leave-one-allele-out tests over the IMGT pack, database
//! relatives placed an allele in a wrong locus only for 12 of 1,277 TCR J,
//! and the panel placed none of 1,500 V alleles of other species wrongly
//! (ALGORITHM §3a).
//!
//! D alleles are too short to place by sequence: they, and every allele
//! this step leaves undecided, take the loci the database names (a D those
//! with D genes, [`crate::prepare`]). Nothing runs when every allele names its
//! locus, as in every shipped reference.

use std::sync::OnceLock;

use crate::align::{AlignJob, Aligner, InterSequenceGotoh, InterSequenceWorkspace};
use crate::alphabet;
use crate::band::Band;
use crate::fasta::parse_fasta;
use crate::germline::{GermlineSet, LocusSource, SegmentKind};
use crate::locus::Loci;
use crate::scoring::{AffineScoring, Score, SubstitutionProfile};

const SCORING: AffineScoring = AffineScoring::new(2, -3, -5, -1);
/// Length of the prefilter's k-mers.
const K: usize = 5;
/// Relatives aligned per allele: those sharing most distinct 5-mers.
pub const CANDIDATES: usize = 8;
/// Half-width of the alignment band around the shared 5-mers' diagonal.
pub const HALF_BAND: i32 = 16;
/// Relatives scoring within this of the best give their loci too.
pub const MARGIN: Score = 15;
/// The best relative must score at least this.
pub const MIN_SCORE: Score = 30;
/// A database relative must score at least this fraction of the allele's
/// score against itself.
pub const MIN_SELF_FRACTION: f64 = 0.3;
/// A panel relative, possibly of another species, at least this fraction:
/// above the 30–40 a chance local alignment of two unrelated V sequences
/// reaches, below what relatives of other vertebrates score (5 % of the
/// pack's V alleles of species outside the panel stay below 0.085; at 0.1
/// the panel assigns none of them a wrong locus).
pub const MIN_PANEL_FRACTION: f64 = 0.1;

/// The embedded panel, FASTA with `locus=` and `kind=` keys.
const PANEL_FASTA: &str = include_str!("../data/locus_panel.fa");

/// A relative of named locus, with what the prefilter and the band need.
struct Relative {
    codes: Vec<u8>,
    loci: Loci,
    /// Distinct 5-mers as a 1024-bit set.
    kmers: [u64; 16],
    /// `(5-mer, position)` of every 5-mer, sorted.
    positions: Vec<(u16, u32)>,
}

impl Relative {
    fn new(codes: Vec<u8>, loci: Loci) -> Self {
        let positions = kmer_positions(&codes);
        let mut kmers = [0u64; 16];
        for &(key, _) in &positions {
            kmers[usize::from(key) / 64] |= 1 << (key % 64);
        }
        Self {
            codes,
            loci,
            kmers,
            positions,
        }
    }

    fn shared(&self, other: &Self) -> u32 {
        self.kmers
            .iter()
            .zip(&other.kmers)
            .map(|(a, b)| (a & b).count_ones())
            .sum()
    }
}

/// `(5-mer, position)` of every 5-mer of `codes` without an N, sorted.
fn kmer_positions(codes: &[u8]) -> Vec<(u16, u32)> {
    let mut out = Vec::with_capacity(codes.len());
    let (mut key, mut valid) = (0u16, 0usize);
    for (at, &code) in codes.iter().enumerate() {
        if code >= alphabet::N {
            valid = 0;
            continue;
        }
        key = ((key << 2) | u16::from(code)) & ((1 << (2 * K)) - 1);
        valid += 1;
        if valid >= K {
            out.push((key, (at + 1 - K) as u32));
        }
    }
    out.sort_unstable();
    out
}

/// The diagonal `d = i − j` (target position minus relative position) that
/// most shared 5-mer occurrences lie on.
fn main_diagonal(target: &Relative, relative: &Relative, votes: &mut Vec<u32>) -> i32 {
    let offset = relative.codes.len();
    votes.clear();
    votes.resize(target.codes.len() + relative.codes.len() + 1, 0);
    let (a, b) = (&target.positions, &relative.positions);
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].0.cmp(&b[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                let key = a[i].0;
                let j_end = j + b[j..].iter().take_while(|p| p.0 == key).count();
                while i < a.len() && a[i].0 == key {
                    for p in &b[j..j_end] {
                        votes[a[i].1 as usize + offset - p.1 as usize] += 1;
                    }
                    i += 1;
                }
                j = j_end;
            }
        }
    }
    let best = (0..votes.len())
        .max_by_key(|&d| (votes[d], std::cmp::Reverse(d)))
        .unwrap_or(offset);
    best as i32 - offset as i32
}

/// Relatives of named locus with their alignment profiles.
struct Pool {
    relatives: Vec<Relative>,
    profiles: Vec<SubstitutionProfile>,
}

impl Pool {
    fn new(relatives: Vec<Relative>) -> Self {
        let profiles = relatives
            .iter()
            .map(|r| SubstitutionProfile::build(&SCORING, &r.codes))
            .collect();
        Self {
            relatives,
            profiles,
        }
    }
}

/// The local scores of each target against its [`CANDIDATES`] relatives
/// of `pool` sharing most 5-mers (`skip(target, relative)` excluded), as
/// `(relative, score)`. Every alignment runs in one batch of the SIMD
/// kernel, whose scores equal the scalar kernel's.
fn score_relatives(
    targets: &[&Relative],
    pool: &Pool,
    skip: impl Fn(usize, usize) -> bool,
) -> Vec<Vec<(usize, Score)>> {
    let mut votes = Vec::new();
    // (target, relative, band) of each alignment.
    let mut planned: Vec<(usize, usize, Band)> = Vec::new();
    let mut ranked: Vec<(u32, usize)> = Vec::new();
    for (t, target) in targets.iter().enumerate() {
        ranked.clear();
        ranked.extend(
            pool.relatives
                .iter()
                .enumerate()
                .filter(|&(at, _)| !skip(t, at))
                .map(|(at, r)| (target.shared(r), at))
                .filter(|&(shared, _)| shared > 0),
        );
        ranked.sort_unstable_by_key(|&(shared, at)| (std::cmp::Reverse(shared), at));
        for &(_, at) in ranked.iter().take(CANDIDATES) {
            let diagonal = main_diagonal(target, &pool.relatives[at], &mut votes);
            planned.push((t, at, Band::around(diagonal, HALF_BAND)));
        }
    }
    let jobs: Vec<AlignJob<'_>> = planned
        .iter()
        .map(|&(t, at, band)| AlignJob {
            query: &targets[t].codes,
            profile: &pool.profiles[at],
            gaps: SCORING.gaps,
            band,
        })
        .collect();
    let mut scores = Vec::new();
    InterSequenceGotoh::new().score_batch(
        &jobs,
        &mut InterSequenceWorkspace::default(),
        &mut scores,
    );
    let mut scored: Vec<Vec<(usize, Score)>> = vec![Vec::new(); targets.len()];
    for (&(t, at, _), hit) in planned.iter().zip(&scores) {
        if let Some(hit) = hit {
            scored[t].push((at, hit.score));
        }
    }
    scored
}

/// The loci of each target from `pool`: those of every relative within
/// [`MARGIN`] of the best local score, if the best reaches the target's
/// floor; `None` otherwise.
fn decide(targets: &[(&Relative, Score)], pool: &Pool) -> Vec<Option<Loci>> {
    let relatives: Vec<&Relative> = targets.iter().map(|&(t, _)| t).collect();
    score_relatives(&relatives, pool, |_, _| false)
        .iter()
        .zip(targets)
        .map(|(scored, &(_, floor))| {
            let top = scored.iter().map(|&(_, score)| score).max()?;
            (top >= floor).then(|| {
                scored
                    .iter()
                    .filter(|&&(_, score)| score + MARGIN >= top)
                    .fold(Loci::NONE, |loci, &(at, _)| {
                        loci.union(pool.relatives[at].loci)
                    })
            })
        })
        .collect()
}

/// The panel's relatives of `kind` (V, J or C), parsed once.
fn panel(kind: SegmentKind) -> &'static Pool {
    static PANEL: OnceLock<[Pool; 3]> = OnceLock::new();
    let panel = PANEL.get_or_init(|| {
        let mut kinds: [Vec<Relative>; 3] = Default::default();
        for record in parse_fasta(PANEL_FASTA).expect("the embedded panel is FASTA") {
            let slot = match record.attribute("kind") {
                Some("V") => 0,
                Some("J") => 1,
                Some("C") => 2,
                _ => panic!("panel record {} has no kind", record.id()),
            };
            let loci = record
                .attribute("locus")
                .and_then(Loci::parse_list)
                .expect("panel records name their loci");
            kinds[slot].push(Relative::new(alphabet::encode(&record.sequence), loci));
        }
        kinds.map(Pool::new)
    });
    match kind {
        SegmentKind::V => &panel[0],
        SegmentKind::J => &panel[1],
        SegmentKind::C | SegmentKind::D => &panel[2],
    }
}

/// Give each allele of `set` whose locus is unknown ([`LocusSource::Any`])
/// the loci of its closest relatives of named locus, in the set or the
/// embedded panel ([`LocusSource::Relative`]); D alleles are left alone.
/// Returns how many alleles got loci.
pub fn infer_loci(set: &mut GermlineSet) -> usize {
    let kind = set.kind();
    if kind == SegmentKind::D {
        return 0;
    }
    let ids: Vec<u32> = (0..set.len() as u32)
        .filter(|&id| set.allele(id).locus_source == LocusSource::Any)
        .collect();
    if ids.is_empty() {
        return 0;
    }
    let targets: Vec<Relative> = ids
        .iter()
        .map(|&id| Relative::new(set.allele(id).codes.clone(), Loci::ANY))
        .collect();
    // A floor of `fraction` of each target's score against itself.
    let floors = |fraction: f64| -> Vec<(&Relative, Score)> {
        targets
            .iter()
            .map(|target| {
                let own = 2 * target.codes.iter().filter(|&&c| c < alphabet::N).count();
                let floor = (fraction * own as f64).ceil() as Score;
                (target, MIN_SCORE.max(floor))
            })
            .collect()
    };
    let named = Pool::new(
        set.alleles()
            .iter()
            .filter(|a| matches!(a.locus_source, LocusSource::Named | LocusSource::Header))
            .map(|a| Relative::new(a.codes.clone(), a.loci))
            .collect(),
    );
    let mut loci = decide(&floors(MIN_SELF_FRACTION), &named);
    // The panel decides what the database's own alleles leave open.
    let open: Vec<usize> = (0..ids.len()).filter(|&k| loci[k].is_none()).collect();
    if !open.is_empty() {
        let panel_floors = floors(MIN_PANEL_FRACTION);
        let open_targets: Vec<(&Relative, Score)> = open.iter().map(|&k| panel_floors[k]).collect();
        for (k, found) in open.into_iter().zip(decide(&open_targets, panel(kind))) {
            loci[k] = found;
        }
    }
    let mut inferred = 0;
    for (&id, found) in ids.iter().zip(loci) {
        if let Some(found) = found {
            set.set_loci(id, found, LocusSource::Relative);
            inferred += 1;
        }
    }
    inferred
}

/// Self-normalised local score ([`normalised`]) at which a V allele named
/// without `*` joins a gene: about 2 % divergence.
pub const SAME_GENE: f64 = 0.95;

/// `score` over the smaller of the two sequences' scores against themselves,
/// so a 5′-partial allele still matches its full-length relative.
fn normalised(score: Score, a: &Relative, b: &Relative) -> f64 {
    let own = |r: &Relative| 2 * r.codes.iter().filter(|&&c| c < alphabet::N).count();
    f64::from(score) / own(a).min(own(b)).max(1) as f64
}

/// Genes for the V alleles of `set` whose names hold no `*` (REVIEW §0
/// item 11): such a name is otherwise a gene of its own, and the retrieval
/// budget, gene frames and the gene counts of the priors key on genes. An
/// allele joins the gene of its nearest `*`-named relative when their
/// self-normalised local score reaches [`SAME_GENE`]; the rest are joined
/// to each other at the same level (single linkage), a group taking the
/// name of its first allele as its gene. Names with `*`, as in every
/// shipped reference, keep their genes. `None` when no gene changes.
pub fn group_genes(set: &GermlineSet) -> Option<GermlineSet> {
    if set.kind() != SegmentKind::V {
        return None;
    }
    let alleles = set.alleles();
    let open: Vec<usize> = (0..alleles.len())
        .filter(|&i| !alleles[i].name.contains('*') && alleles[i].gene_key.is_none())
        .collect();
    if open.is_empty() {
        return None;
    }
    let relative = |i: usize| Relative::new(alleles[i].codes.clone(), alleles[i].loci);
    let named: Vec<usize> = (0..alleles.len())
        .filter(|&i| alleles[i].name.contains('*'))
        .collect();
    let named_pool = Pool::new(named.iter().map(|&i| relative(i)).collect());
    let targets: Vec<Relative> = open.iter().map(|&i| relative(i)).collect();
    let target_refs: Vec<&Relative> = targets.iter().collect();
    let mut keys: Vec<Option<String>> = alleles.iter().map(|a| a.gene_key.clone()).collect();
    let mut changed = false;
    let mut left = Vec::new();
    for (k, scored) in score_relatives(&target_refs, &named_pool, |_, _| false)
        .into_iter()
        .enumerate()
    {
        let best = scored
            .iter()
            .map(|&(at, score)| {
                (
                    normalised(score, &targets[k], &named_pool.relatives[at]),
                    at,
                )
            })
            .filter(|&(level, _)| level >= SAME_GENE)
            .max_by(|a, b| a.0.total_cmp(&b.0).then(b.1.cmp(&a.1)));
        match best {
            Some((_, at)) => {
                keys[open[k]] = Some(alleles[named[at]].gene_label().to_string());
                changed = true;
            }
            None => left.push(k),
        }
    }
    // Single linkage among the alleles no named gene took.
    let pool = Pool::new(left.iter().map(|&k| relative(open[k])).collect());
    let refs: Vec<&Relative> = pool.relatives.iter().collect();
    let mut parent: Vec<usize> = (0..left.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for (a, scored) in score_relatives(&refs, &pool, |t, at| t == at)
        .into_iter()
        .enumerate()
    {
        for (b, score) in scored {
            if normalised(score, refs[a], refs[b]) >= SAME_GENE {
                let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
                // The group's root is its first allele.
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
    }
    for k in 0..left.len() {
        let first = root(&mut parent, k);
        if first != k {
            keys[open[left[k]]] = Some(alleles[open[left[first]]].name.clone());
            changed = true;
        }
    }
    changed.then(|| set.with_gene_keys(&keys))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::decode;
    use crate::locus::Locus;
    use crate::testing::{random_codes, rng};

    /// Human IGHV3-23*01, IGKV1-5*01 and IGLV2-14*01.
    const IGHV3_23: &str = "GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAAGA";
    const IGKV1_5: &str = "GACATCCAGATGACCCAGTCTCCTTCCACCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCCAGTCAGAGTATTAGTAGCTGGTTGGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGATGCCTCCAGTTTGGAAAGTGGGGTCCCATCAAGGTTCAGCGGCAGTGGATCTGGGACAGAATTCACTCTCACCATCAGCAGCCTGCAGCCTGATGATTTTGCAACTTATTACTGCCAACAGTATAATAGTTATTCTCC";
    const IGLV2_14: &str = "CAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACTGGAACCAGCAGTGACGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAGCACCCAGGCAAAGCCCCCAAACTCATGATTTATGAGGTCAGTAATCGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCTGGCAACACGGCCTCCCTGACCATCTCTGGGCTCCAGGCTGAGGACGAGGCTGATTATTACTGCAGCTCATATACAAGCAGCAGCACTCTC";

    /// `sequence` with every `step`-th base changed.
    fn mutated(sequence: &str, step: usize) -> String {
        sequence
            .bytes()
            .enumerate()
            .map(|(i, b)| match (i % step == step / 2, b) {
                (true, b'A') => 'C',
                (true, _) => 'A',
                (false, b) => b as char,
            })
            .collect()
    }

    /// Custom names take the locus of their nearest named relative in the
    /// database, then of the panel; an allele related to nothing keeps its
    /// unknown locus; named alleles are untouched.
    #[test]
    fn custom_names_take_their_relatives_loci() {
        let mut next = rng(21);
        let unrelated = decode(&random_codes(&mut next, 290));
        let text = [
            format!(">IGHV3-23*01\n{IGHV3_23}\n"),
            format!(">IGLV2-14*01\n{IGLV2_14}\n"),
            // 5 % divergent from a named allele of the database.
            format!(">novel_h\n{}\n", mutated(IGHV3_23, 20)),
            format!(">novel_l\n{}\n", mutated(IGLV2_14, 20)),
            // A kappa V: no kappa V in the database, the panel's IGKV1-12
            // is its relative.
            format!(">sample_k\n{}\n", mutated(IGKV1_5, 25)),
            format!(">junk\n{unrelated}\n"),
        ]
        .concat();
        let mut set = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
        assert_eq!(infer_loci(&mut set), 3);
        let of = |locus| Loci::of(locus);
        let got: Vec<(Loci, LocusSource)> = set
            .alleles()
            .iter()
            .map(|a| (a.loci, a.locus_source))
            .collect();
        assert_eq!(
            got,
            [
                (of(Locus::Igh), LocusSource::Named),
                (of(Locus::Igl), LocusSource::Named),
                (of(Locus::Igh), LocusSource::Relative),
                (of(Locus::Igl), LocusSource::Relative),
                (of(Locus::Igk), LocusSource::Relative),
                (Loci::ANY, LocusSource::Any),
            ]
        );
        // Nothing to infer: nothing runs.
        let mut named =
            GermlineSet::from_fasta(SegmentKind::V, &format!(">IGHV3-23*01\n{IGHV3_23}\n"))
                .unwrap();
        assert_eq!(infer_loci(&mut named), 0);
    }

    /// The panel places human and mouse J and C families of every locus,
    /// TRAV/DV relatives in both TRA and TRD; D is never placed.
    #[test]
    fn the_panel_places_every_kind_but_d() {
        for (kind, slot) in [
            (SegmentKind::V, 0),
            (SegmentKind::J, 1),
            (SegmentKind::C, 2),
        ] {
            let relatives = &panel(kind).relatives;
            assert!(!relatives.is_empty(), "{kind:?}");
            let loci = relatives.iter().fold(Loci::NONE, |l, r| l.union(r.loci));
            assert_eq!(loci, Loci::ANY, "{kind:?} {slot}");
        }
        let tra_trd = Loci::parse_list("TRA,TRD").unwrap();
        assert!(
            panel(SegmentKind::V)
                .relatives
                .iter()
                .any(|r| r.loci == tra_trd)
        );
        let mut d = GermlineSet::from_fasta(SegmentKind::D, ">d1\nGGTATAACTGGAACGAC\n").unwrap();
        assert_eq!(infer_loci(&mut d), 0);
        assert_eq!(d.allele(0).loci, Loci::ANY);
    }

    /// The band's diagonal follows a shift between the sequences: a
    /// relative missing its first 40 bases still scores its whole overlap.
    #[test]
    fn the_band_follows_the_shared_diagonal() {
        let full = Relative::new(alphabet::encode(IGHV3_23.as_bytes()), Loci::ANY);
        let cut = Relative::new(alphabet::encode(&IGHV3_23.as_bytes()[40..]), Loci::ANY);
        let mut votes = Vec::new();
        assert_eq!(main_diagonal(&full, &cut, &mut votes), 40);
        assert_eq!(main_diagonal(&cut, &full, &mut votes), -40);
        let loci = Loci::of(Locus::Igh);
        let pool = Pool::new(vec![Relative::new(cut.codes.clone(), loci)]);
        let own = 2 * (IGHV3_23.len() - 40) as Score;
        assert_eq!(
            decide(&[(&full, own), (&full, own + 1)], &pool),
            [Some(loci), None]
        );
    }
}
