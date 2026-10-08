//! Germline annotation accuracy on synthetic genes with known boundaries.
//!
//! V genes are assembled from heavy / light framework templates (IMGT layout:
//! heavy FR1 25 aa with Cys23 at index 21, light FR1 26 aa with Cys23 at 22;
//! FR2 17 aa with Trp41 at index 2; FR3 38 / 36 aa ending in Cys104), random
//! CDR1 / CDR2 of realistic lengths drawn from all 20 residues (so decoy C and
//! W occur), random mutation of non-conserved framework residues, random
//! synonymous codons, and optional 5′ truncation that shifts the frame.

use vdj_assign::annotate::{AnchorSource, AnnotateOptions, ChainType, annotate};
use vdj_assign::prelude::*;
use vdj_assign::sim::Rng;
use vdj_assign::translate::{codon, translate};

const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Synonymous codons per residue, as ASCII triplets.
fn codons_for(residue: u8) -> Vec<[u8; 3]> {
    let bases = b"ACGT";
    let mut out = Vec::new();
    for a in 0..4u8 {
        for b in 0..4u8 {
            for c in 0..4u8 {
                if codon(a, b, c) == residue {
                    out.push([bases[a as usize], bases[b as usize], bases[c as usize]]);
                }
            }
        }
    }
    out
}

fn back_translate(protein: &[u8], rng: &mut Rng) -> Vec<u8> {
    protein
        .iter()
        .flat_map(|&aa| {
            let options = codons_for(aa);
            options[rng.below(options.len())]
        })
        .collect()
}

fn random_residues(len: usize, rng: &mut Rng) -> Vec<u8> {
    (0..len)
        .map(|_| RESIDUES[rng.below(RESIDUES.len())])
        .collect()
}

/// Mutate ~12 % of residues except the conserved `keep` indices.
fn mutate(template: &[u8], keep: &[usize], rng: &mut Rng) -> Vec<u8> {
    template
        .iter()
        .enumerate()
        .map(|(i, &aa)| {
            if !keep.contains(&i) && rng.below(100) < 12 {
                RESIDUES[rng.below(20)]
            } else {
                aa
            }
        })
        .collect()
}

struct SyntheticV {
    sequence: Vec<u8>,
    anchor: u32,
    starts: [u32; 5],
}

fn synthetic_v(chain: ChainType, rng: &mut Rng, truncate: usize) -> SyntheticV {
    let (fr1, cys23, cdr1_len, fr2, cdr2_len, fr3, tail) = match chain {
        ChainType::Heavy => (
            &b"QVQLVQSGAEVKKPGASVKVSCKAS"[..],
            21,
            5 + rng.below(8),
            &b"MHWVRQAPGQRLEWMGW"[..],
            3 + rng.below(8),
            &b"KYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYC"[..],
            &b"AR"[..],
        ),
        ChainType::Light => (
            &b"DIQMTQSPSSLSASVGDRVTITCRAS"[..],
            22,
            5 + rng.below(8),
            &b"LAWYQQKPGKAPKLLIY"[..],
            3 + rng.below(5),
            &b"SLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYC"[..],
            &b"QQSYSTP"[..],
        ),
    };
    let fr1 = mutate(fr1, &[cys23, cys23 - 1], rng);
    let fr2 = mutate(fr2, &[2, 5, 7, 8], rng);
    let fr3 = mutate(fr3, &[fr3.len() - 1, fr3.len() - 2, fr3.len() - 3], rng);
    let mut protein = fr1.clone();
    protein.extend(
        random_residues(cdr1_len, rng)
            .iter()
            .map(|&a| if a == b'C' { b'S' } else { a }),
    );
    let fr2_start = protein.len();
    protein.extend(&fr2);
    let cdr2_start = protein.len();
    protein.extend(random_residues(cdr2_len, rng));
    let fr3_start = protein.len();
    protein.extend(&fr3);
    let cys104 = protein.len() - 1;
    protein.extend(tail);
    let mut sequence = back_translate(&protein, rng);
    sequence.extend(b"GA"); // partial codon at the V 3' end
    let shift = |codon_index: usize| (3 * codon_index - truncate) as u32;
    SyntheticV {
        sequence: sequence[truncate..].to_vec(),
        anchor: shift(cys104),
        starts: [
            ((3 - truncate) % 3) as u32,
            shift(fr1.len()),
            shift(fr2_start),
            shift(cdr2_start),
            shift(fr3_start),
        ],
    }
}

fn synthetic_j(chain: ChainType, rng: &mut Rng) -> (Vec<u8>, u32) {
    let (protein, anchor_index) = match chain {
        ChainType::Heavy => (mutate(b"YFDYWGQGTLVTVSS", &[4, 5, 7], rng), 4),
        ChainType::Light => (mutate(b"WTFGQGTKVEIK", &[2, 3, 5], rng), 2),
    };
    let lead = rng.below(3);
    let mut sequence: Vec<u8> = (0..lead).map(|_| b"ACGT"[rng.below(4)]).collect();
    sequence.extend(back_translate(&protein, rng));
    sequence.push(b'G');
    (sequence, (lead + 3 * anchor_index) as u32)
}

fn database(
    chain: ChainType,
    genes: usize,
    rng: &mut Rng,
) -> (GermlineSet, Vec<SyntheticV>, GermlineSet, Vec<u32>) {
    let mut v = GermlineSet::new(SegmentKind::V);
    let mut truth = Vec::new();
    for g in 0..genes {
        let truncate = if g % 7 == 3 { 1 + rng.below(2) } else { 0 };
        let gene = synthetic_v(chain, rng, truncate);
        v.push(&format!("V{g}*01"), &gene.sequence, None).unwrap();
        truth.push(gene);
    }
    let mut j = GermlineSet::new(SegmentKind::J);
    let mut j_truth = Vec::new();
    for g in 0..6 {
        let (sequence, anchor) = synthetic_j(chain, rng);
        j.push(&format!("J{g}*01"), &sequence, None).unwrap();
        j_truth.push(anchor);
    }
    (v, truth, j, j_truth)
}

#[test]
fn recovers_exact_boundaries_for_heavy_and_light_genes() {
    for (chain, seed) in [(ChainType::Heavy, 1), (ChainType::Light, 2)] {
        let mut rng = Rng::new(seed);
        let (mut v, truth, mut j, j_truth) = database(chain, 300, &mut rng);
        let report = annotate(&mut v, Some(&mut j), AnnotateOptions::default());
        assert_eq!(report.chain, Some(chain));
        for (allele, expected) in v.alleles().iter().zip(&truth) {
            assert_eq!(
                allele.anchor,
                Some(expected.anchor),
                "{chain:?} {} anchor",
                allele.name
            );
            assert_eq!(
                allele.regions.values(),
                expected.starts.map(Some),
                "{chain:?} {} regions ({})",
                allele.name,
                String::from_utf8_lossy(&translate(&allele.codes, expected.anchor as usize % 3))
            );
        }
        for (allele, &expected) in j.alleles().iter().zip(&j_truth) {
            assert_eq!(allele.anchor, Some(expected), "{chain:?} {}", allele.name);
        }
        assert_eq!(report.uncertain(), 0, "{}", report.summary());
    }
}

#[test]
fn pseudogene_anchor_is_transferred_from_its_closest_sibling() {
    let mut rng = Rng::new(3);
    let (mut v, truth, mut j, _) = database(ChainType::Heavy, 40, &mut rng);
    // A copy of V0 with Cys104 and Trp41 destroyed and two substitutions elsewhere.
    let mut broken = truth[0].sequence.clone();
    let anchor = truth[0].anchor as usize;
    broken[anchor..anchor + 3].copy_from_slice(b"AGT");
    let trp = truth[0].starts[2] as usize + 6;
    broken[trp..trp + 3].copy_from_slice(b"TTG");
    broken[30] = if broken[30] == b'A' { b'C' } else { b'A' };
    v.push("V0*02", &broken, None).unwrap();
    let report = annotate(&mut v, Some(&mut j), AnnotateOptions::default());
    let rescued = v.alleles().last().unwrap();
    assert_eq!(rescued.anchor, Some(truth[0].anchor));
    assert_eq!(rescued.regions.values(), truth[0].starts.map(Some));
    let annotation = report.v.last().unwrap();
    assert_eq!(
        annotation.source,
        AnchorSource::Transfer {
            from: "V0*01".into()
        }
    );
    assert!(
        annotation
            .notes
            .iter()
            .any(|n| n.contains("conserved residue is S"))
    );
}

/// A 3′-partial allele (IMGT "partial in 3′") ends inside its Cys104 codon.
/// Its closest sibling places the anchor past the end, so it gets no anchor,
/// even when an out-of-frame TGT near the end looks like a weak Cys104.
#[test]
fn three_prime_partial_allele_gets_no_anchor() {
    let mut rng = Rng::new(5);
    let (mut v, truth, mut j, _) = database(ChainType::Heavy, 40, &mut rng);
    let anchor = truth[0].anchor as usize;
    let mut partial = truth[0].sequence[..anchor + 2].to_vec();
    let decoy = anchor - 31;
    partial[decoy..decoy + 3].copy_from_slice(b"TGT");
    v.push("V0*02", &partial, None).unwrap();
    let report = annotate(&mut v, Some(&mut j), AnnotateOptions::default());
    let allele = v.alleles().last().unwrap();
    assert_eq!(allele.anchor, None);
    // Region starts inside the allele still come from the sibling.
    assert_eq!(allele.regions.values(), truth[0].starts.map(Some));
    let annotation = report.v.last().unwrap();
    assert_eq!(annotation.source, AnchorSource::Missing);
    assert!(
        annotation
            .notes
            .iter()
            .any(|n| n.contains("past the 3′ end")),
        "{:?}",
        annotation.notes
    );
}

/// An allele that keeps its Cys104 but lost 21 bases just before it: the
/// sibling alignment stops at the deletion, so the extrapolated anchor lands
/// past the end. That is not evidence of a 3′-partial allele.
#[test]
fn a_deletion_before_cys104_is_not_a_partial_allele() {
    let mut rng = Rng::new(6);
    let (mut v, truth, mut j, _) = database(ChainType::Heavy, 40, &mut rng);
    let anchor = truth[0].anchor as usize;
    let mut deleted = truth[0].sequence[..anchor - 21].to_vec();
    deleted.extend_from_slice(&truth[0].sequence[anchor..]);
    v.push("V0*02", &deleted, None).unwrap();
    let report = annotate(&mut v, Some(&mut j), AnnotateOptions::default());
    let annotation = report.v.last().unwrap();
    assert!(
        !annotation.notes.iter().any(|n| n.contains("3′-partial")),
        "{annotation:?}"
    );
    let allele = v.alleles().last().unwrap();
    assert_eq!(allele.anchor, Some(truth[0].anchor - 21), "{annotation:?}");
}

#[test]
fn header_values_are_never_overridden() {
    let mut rng = Rng::new(4);
    let gene = synthetic_v(ChainType::Heavy, &mut rng, 0);
    let text = format!(
        ">V1*01 anchor={} cdr1=3\n{}\n",
        gene.anchor,
        String::from_utf8(gene.sequence.clone()).unwrap()
    );
    let mut v = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
    let report = annotate(
        &mut v,
        None,
        AnnotateOptions {
            chain: Some(ChainType::Heavy),
        },
    );
    assert_eq!(report.v[0].source, AnchorSource::Header);
    let regions = v.allele(0).regions;
    assert_eq!(regions.cdr1, Some(3));
    assert_eq!(regions.fr2, Some(gene.starts[2]));
    assert_eq!(regions.fr3, Some(gene.starts[4]));
}

#[test]
fn annotation_of_thousands_of_alleles_takes_milliseconds() {
    let mut rng = Rng::new(5);
    let (mut v, _, mut j, _) = database(ChainType::Heavy, 3000, &mut rng);
    let started = std::time::Instant::now();
    annotate(&mut v, Some(&mut j), AnnotateOptions::default());
    let elapsed = started.elapsed();
    eprintln!("annotated 3000 V + 6 J alleles in {elapsed:?}");
    assert!(elapsed.as_millis() < if cfg!(debug_assertions) { 3000 } else { 200 });
}

/// 12,000 adversarial alleles per chain type: anchors, FR2, CDR2 and FR3 must
/// be exact; a CDR1 start may only be wrong where the report flags Cys23 as
/// ambiguous (a mutation-made decoy Cys next to Cys23 in a 5′-truncated allele).
#[test]
fn every_error_is_flagged_across_many_databases() {
    for chain in [ChainType::Heavy, ChainType::Light] {
        let (mut total, mut flagged_errors) = (0, 0);
        for seed in 100..120 {
            let mut rng = Rng::new(seed);
            let (mut v, truth, mut j, j_truth) = database(chain, 300, &mut rng);
            let report = annotate(&mut v, Some(&mut j), AnnotateOptions::default());
            for (id, (allele, expected)) in v.alleles().iter().zip(&truth).enumerate() {
                total += 1;
                assert_eq!(
                    allele.anchor,
                    Some(expected.anchor),
                    "{chain:?} seed {seed} {}",
                    allele.name
                );
                let got = allele.regions.values();
                let want = expected.starts.map(Some);
                for field in [0, 2, 3, 4] {
                    assert_eq!(
                        got[field], want[field],
                        "{chain:?} seed {seed} {} field {field}",
                        allele.name
                    );
                }
                if got[1] != want[1] {
                    assert!(
                        report.v[id]
                            .notes
                            .iter()
                            .any(|n| n.contains("Cys23 ambiguous")),
                        "{chain:?} seed {seed} {}: silent CDR1 error",
                        allele.name
                    );
                    flagged_errors += 1;
                }
            }
            for (allele, &expected) in j.alleles().iter().zip(&j_truth) {
                assert_eq!(allele.anchor, Some(expected));
            }
        }
        eprintln!(
            "{chain:?}: {total} alleles, {flagged_errors} flagged CDR1 errors, no silent errors"
        );
    }
}

#[test]
fn ambiguous_cys23_of_a_truncated_allele_is_resolved_from_its_sibling() {
    let mut rng = Rng::new(6);
    let (mut v, truth, mut j, _) = database(ChainType::Light, 30, &mut rng);
    // A 5′-truncated copy of V0 whose Cys23 gets a decoy Cys right after it.
    let mut sequence = truth[0].sequence.clone();
    let cys23 = truth[0].starts[1] as usize - 12;
    sequence[cys23 + 3..cys23 + 6].copy_from_slice(b"TGC");
    let truncated = sequence[2..].to_vec();
    v.push("V0*02", &truncated, None).unwrap();
    let report = annotate(&mut v, Some(&mut j), AnnotateOptions::default());
    let allele = v.alleles().last().unwrap();
    assert_eq!(allele.regions.cdr1, Some(truth[0].starts[1] - 2));
    let notes = &report.v.last().unwrap().notes;
    assert!(
        notes
            .iter()
            .any(|n| n == "ambiguous Cys23 resolved from V0*01"),
        "{notes:?}"
    );
}

#[test]
fn read_regions_match_the_simulated_rearrangement() {
    use vdj_assign::prepare::prepare;
    use vdj_assign::read_regions::ReadRegions;
    let mut rng = Rng::new(7);
    let mut v = GermlineSet::new(SegmentKind::V);
    let mut v_truth = Vec::new();
    for g in 0..30 {
        let gene = synthetic_v(ChainType::Heavy, &mut rng, 0);
        v.push(&format!("V{g}*01"), &gene.sequence, None).unwrap();
        v_truth.push(gene);
    }
    let mut j = GermlineSet::new(SegmentKind::J);
    let mut j_truth = Vec::new();
    for g in 0..4 {
        let (sequence, anchor) = synthetic_j(ChainType::Heavy, &mut rng);
        j.push(&format!("J{g}*01"), &sequence, None).unwrap();
        j_truth.push((sequence, anchor));
    }
    let mut d = GermlineSet::new(SegmentKind::D);
    let d_seqs: Vec<Vec<u8>> = (0..6)
        .map(|_| (0..17).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    for (g, seq) in d_seqs.iter().enumerate() {
        d.push(&format!("D{g}*01"), seq, None).unwrap();
    }
    let db = prepare(v, Some(d), j, AnnotateOptions::default()).unwrap();
    assert_eq!(db.annotation.uncertain(), 0);
    let assigner =
        DefaultAssigner::from_config(&EngineConfig::default(), db.v, db.d, db.j).unwrap();
    let mut ws = assigner.workspace();
    let (mut in_frame, mut out_of_frame) = (0, 0);
    for read in 0..200 {
        let (vi, ji, di) = (rng.below(30), rng.below(4), rng.below(6));
        let gene = &v_truth[vi];
        let (j_seq, j_anchor) = &j_truth[ji];
        let random = |len: usize, rng: &mut Rng| -> Vec<u8> {
            (0..len).map(|_| b"ACGT"[rng.below(4)]).collect()
        };
        let mut query = gene.sequence.clone();
        query.extend(random(2 + rng.below(6), &mut rng));
        query.extend(&d_seqs[di]);
        query.extend(random(2 + rng.below(6), &mut rng));
        let j_start = query.len() as u32;
        query.extend(j_seq);
        let reverse = read % 2 == 1;
        let input = if reverse {
            query
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                })
                .collect()
        } else {
            query.clone()
        };
        let call = assigner.assign(&encode(&input), &mut ws);
        assert_eq!(call.orientation.is_reverse(), reverse);
        let regions = ReadRegions::from_call(&call, &assigner);
        let [fr1, cdr1, fr2, cdr2, fr3] = gene.starts;
        assert_eq!(regions.fwr1, Some(fr1..cdr1), "read {read}");
        assert_eq!(regions.cdr1, Some(cdr1..fr2), "read {read}");
        assert_eq!(regions.fwr2, Some(fr2..cdr2), "read {read}");
        assert_eq!(regions.cdr2, Some(cdr2..fr3), "read {read}");
        assert_eq!(regions.fwr3, Some(fr3..gene.anchor + 3), "read {read}");
        let w = j_start + j_anchor;
        assert_eq!(regions.junction, Some(gene.anchor..w + 3), "read {read}");
        assert_eq!(regions.cdr3, Some(gene.anchor + 3..w), "read {read}");
        assert_eq!(regions.fwr4, Some(w..query.len() as u32), "read {read}");
        let frame_ok = (w + 3 - gene.anchor).is_multiple_of(3);
        assert_eq!(regions.vj_in_frame, Some(frame_ok));
        let codes = encode(&query);
        let stop = translate(&codes[..query.len()], (gene.anchor % 3) as usize).contains(&b'*');
        assert_eq!(regions.stop_codon, Some(stop), "read {read}");
        assert_eq!(regions.productive, Some(frame_ok && !stop));
        if frame_ok {
            in_frame += 1
        } else {
            out_of_frame += 1
        }
    }
    assert!(
        in_frame > 20 && out_of_frame > 20,
        "both junction frames exercised"
    );
}

/// Gap columns for IMGT codons `positions`: `codons` in order, `...` at the
/// codons in `gaps`.
fn place(
    out: &mut Vec<u8>,
    codons: &[u8],
    positions: std::ops::RangeInclusive<usize>,
    gaps: &[usize],
) {
    let mut codons = codons.chunks(3);
    for position in positions {
        if gaps.contains(&position) {
            out.extend(b"...");
        } else {
            out.extend(codons.next().expect("enough codons"));
        }
    }
    assert!(codons.next().is_none(), "every codon placed");
}

/// The middle `positions − n` codons of a CDR of `n` codons, IMGT style.
fn middle_gaps(positions: std::ops::RangeInclusive<usize>, n: usize) -> Vec<usize> {
    let first = positions.start() + n.div_ceil(2);
    (first..first + (positions.count() - n)).collect()
}

/// IMGT-gapped text of a full-length synthetic gene: FR1 on IMGT codons
/// 1–26 (heavy lacks 10), CDR1 on 27–38, FR2 on 39–55, CDR2 on 56–65 and
/// FR3 on 66–104 (heavy lacks 73; light 73, 81 and 82), CDR gaps in the
/// middle, then the 3′ tail.
fn imgt_gapped(gene: &SyntheticV, chain: ChainType) -> Vec<u8> {
    let [fr1, cdr1, fr2, cdr2, fr3] = gene.starts.map(|s| s as usize);
    assert_eq!(fr1, 0, "full-length gene");
    let s = &gene.sequence;
    let end = gene.anchor as usize + 3;
    let (fr1_gaps, fr3_gaps): (&[usize], &[usize]) = match chain {
        ChainType::Heavy => (&[10], &[73]),
        ChainType::Light => (&[], &[73, 81, 82]),
    };
    let mut out = Vec::new();
    place(&mut out, &s[..cdr1], 1..=26, fr1_gaps);
    let cdr1_gaps = middle_gaps(27..=38, (fr2 - cdr1) / 3);
    place(&mut out, &s[cdr1..fr2], 27..=38, &cdr1_gaps);
    place(&mut out, &s[fr2..cdr2], 39..=55, &[]);
    let cdr2_gaps = middle_gaps(56..=65, (fr3 - cdr2) / 3);
    place(&mut out, &s[cdr2..fr3], 56..=65, &cdr2_gaps);
    place(&mut out, &s[fr3..end], 66..=104, fr3_gaps);
    out.extend(&s[end..]);
    out
}

fn record(name: &str, header: &str, sequence: &[u8]) -> String {
    format!(">{name}{header}\n{}\n", String::from_utf8_lossy(sequence))
}

fn j_set(chain: ChainType, rng: &mut Rng) -> GermlineSet {
    let mut j = GermlineSet::new(SegmentKind::J);
    for g in 0..4 {
        let (sequence, _) = synthetic_j(chain, rng);
        j.push(&format!("J{g}*01"), &sequence, None).unwrap();
    }
    j
}

/// An IMGT-gapped database is numbered from its gap columns: every anchor
/// and region start is exact, including a Cys104 → Ser pseudogene, a
/// 3′-partial record (no anchor), a 5′-partial record that starts inside
/// CDR1 (no FR1 or CDR1, although a motif hit would offer both), a record
/// with a frameshift marker, and a header key that overrides a gap column.
#[test]
fn imgt_gapped_input_numbers_every_region() {
    for (chain, seed) in [(ChainType::Heavy, 31), (ChainType::Light, 32)] {
        let mut rng = Rng::new(seed);
        let genes: Vec<SyntheticV> = (0..40).map(|_| synthetic_v(chain, &mut rng, 0)).collect();
        let gapped: Vec<Vec<u8>> = genes.iter().map(|g| imgt_gapped(g, chain)).collect();
        let mut text = String::new();
        for (g, columns) in gapped.iter().enumerate() {
            text += &record(&format!("V{g}*01"), "", columns);
        }
        // V0*02: Cys104 → Ser.
        let mut pseudogene = gapped[0].clone();
        pseudogene[309..312].copy_from_slice(b"AGC");
        text += &record("V0*02", "", &pseudogene);
        // V1*02: 3′-partial, ending inside FR3.
        text += &record("V1*02", "", &gapped[1][..250]);
        // V2*02: 5′-partial, starting inside CDR1.
        let mut late = gapped[2].clone();
        late[..84].fill(b'.');
        let lost = gapped[2][..84].iter().filter(|&&c| c != b'.').count() as u32;
        text += &record("V2*02", "", &late);
        // V3*02: a base lost in FR3, marked by one gap column (codon 84).
        let mut shifted = gapped[3].clone();
        shifted[250] = b'.';
        text += &record("V3*02", "", &shifted);
        // V4*02: a header FR3 start one codon later than the gap columns'.
        let fr3 = genes[4].starts[4] + 3;
        text += &record("V4*02", &format!(" fr3={fr3}"), &gapped[4]);

        let mut v = GermlineSet::from_fasta(SegmentKind::V, &text).unwrap();
        let mut j = j_set(chain, &mut rng);
        let report = annotate(&mut v, Some(&mut j), AnnotateOptions::default());
        assert!(
            report
                .summary()
                .contains("IMGT numbering read from 45 gapped V alleles"),
            "{}",
            report.summary()
        );
        for (allele, gene) in v.alleles().iter().zip(&genes) {
            assert_eq!(allele.anchor, Some(gene.anchor), "{}", allele.name);
            assert_eq!(
                allele.regions.values(),
                gene.starts.map(Some),
                "{}",
                allele.name
            );
        }
        let by_name = |name: &str| {
            let id = v.alleles().iter().position(|a| a.name == name).unwrap();
            (&v.alleles()[id], &report.v[id])
        };
        for name in ["V0*01", "V0*02", "V2*02", "V3*02", "V4*02"] {
            assert_eq!(by_name(name).1.source, AnchorSource::ImgtGaps, "{name}");
        }
        let (allele, annotation) = by_name("V0*02");
        assert_eq!(allele.anchor, Some(genes[0].anchor));
        assert!(
            annotation.notes.iter().any(|n| n == "IMGT 104 codon is S"),
            "{annotation:?}"
        );

        let (allele, annotation) = by_name("V1*02");
        assert_eq!(allele.anchor, None);
        assert_eq!(annotation.source, AnchorSource::Missing);
        assert_eq!(allele.regions.values(), genes[1].starts.map(Some));
        assert!(
            annotation
                .notes
                .iter()
                .any(|n| n == "3′-partial: the IMGT numbering ends before Cys104"),
            "{annotation:?}"
        );

        let (allele, _) = by_name("V2*02");
        let [_, _, fr2, cdr2, fr3_start] = genes[2].starts;
        assert_eq!(allele.anchor, Some(genes[2].anchor - lost));
        assert_eq!(
            allele.regions.values(),
            [
                None,
                None,
                Some(fr2 - lost),
                Some(cdr2 - lost),
                Some(fr3_start - lost)
            ]
        );

        let (allele, annotation) = by_name("V3*02");
        assert_eq!(allele.anchor, Some(genes[3].anchor - 1));
        assert_eq!(allele.regions.values(), genes[3].starts.map(Some));
        assert!(
            annotation
                .notes
                .iter()
                .any(|n| n == "IMGT frameshift marker in codon 84"),
            "{annotation:?}"
        );

        let (allele, _) = by_name("V4*02");
        assert_eq!(allele.regions.fr3, Some(fr3));
        assert_eq!(allele.anchor, Some(genes[4].anchor));
    }
}

/// A codon alignment written with `.` but not in IMGT columns (one extra
/// gap codon in FR1 moves Cys23, Trp41 and Cys104 off IMGT 23, 41 and
/// 104) is not read as IMGT numbering: annotation equals the ungapped run.
#[test]
fn non_imgt_gapped_alignment_is_ignored() {
    let mut rng = Rng::new(33);
    let genes: Vec<SyntheticV> = (0..40)
        .map(|_| synthetic_v(ChainType::Heavy, &mut rng, 0))
        .collect();
    let (mut shifted, mut plain) = (String::new(), String::new());
    for (g, gene) in genes.iter().enumerate() {
        let mut columns = imgt_gapped(gene, ChainType::Heavy);
        columns.splice(30..30, *b"...");
        shifted += &record(&format!("V{g}*01"), "", &columns);
        plain += &record(&format!("V{g}*01"), "", &gene.sequence);
    }
    let j = j_set(ChainType::Heavy, &mut rng);
    let run = |text: &str| {
        let mut v = GermlineSet::from_fasta(SegmentKind::V, text).unwrap();
        let report = annotate(&mut v, Some(&mut j.clone()), AnnotateOptions::default());
        (v, report)
    };
    let (v, report) = run(&shifted);
    let (expected, expected_report) = run(&plain);
    assert!(
        report.summary().contains(
            "gaps of 40 V alleles ignored, not IMGT numbering (0 of 40 gapped V alleles have C at IMGT 23)"
        ),
        "{}",
        report.summary()
    );
    let pairs = v.alleles().iter().zip(&report.v);
    let expected_pairs = expected.alleles().iter().zip(&expected_report.v);
    for ((allele, annotation), (want, want_annotation)) in pairs.zip(expected_pairs) {
        assert_eq!(allele.anchor, want.anchor, "{}", allele.name);
        assert_eq!(allele.regions, want.regions, "{}", allele.name);
        assert_eq!(annotation.source, want_annotation.source);
        assert_eq!(annotation.notes, want_annotation.notes);
    }
}
