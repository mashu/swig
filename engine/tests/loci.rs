//! One database of several loci, as Swig composes a BCR or TCR study: every
//! rearrangement stays within one locus.

use swig::germline::GermlineSet;
use swig::locus::{Loci, Locus};
use swig::prelude::*;
use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

fn fasta(set: &GermlineSet) -> String {
    let mut text = Vec::new();
    set.write_fasta(&mut text, |_| None).unwrap();
    String::from_utf8(text).unwrap()
}

/// A heavy-chain database (IGH V, D, J) and a light-chain one built the same
/// way from another seed, renamed to IGK, without D.
fn databases() -> (
    SyntheticDatabase,
    SyntheticDatabase,
    GermlineSet,
    GermlineSet,
) {
    let heavy = SyntheticDatabase::generate(&SyntheticSpec::default(), 1);
    let light = SyntheticDatabase::generate(&SyntheticSpec::default(), 2);
    // The written `locus=IGH` would win over the new labels: drop it, so the
    // loci are read from the labels.
    let rename = |set: &GermlineSet, kind, from: &str, to: &str| {
        let text = fasta(set).replace(from, to).replace(" locus=IGH", "");
        GermlineSet::from_fasta(kind, &text).unwrap()
    };
    let light_v = rename(&light.v, SegmentKind::V, "IGHV", "IGKV");
    let light_j = rename(&light.j, SegmentKind::J, "IGHJ", "IGKJ");
    (heavy, light, light_v, light_j)
}

#[test]
fn rearrangements_stay_within_one_locus() {
    let (heavy, _, light_v, light_j) = databases();
    let v = fasta(&heavy.v) + &fasta(&light_v);
    let j = fasta(&heavy.j) + &fasta(&light_j);
    let d = fasta(&heavy.d);
    let texts = GermlineTexts {
        v: &v,
        d: Some(&d),
        j: &j,
        c: None,
    };
    let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
    assert_eq!(
        database.v.loci(),
        Loci::of(Locus::Igh).union(Loci::of(Locus::Igk))
    );
    // Each locus is annotated with its own chain type and models.
    assert_eq!(database.annotation.parts.len(), 2);
    let config = CallingProfile::ROptimized.config();
    let assigner =
        DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap();

    let heavy_reads = Simulator::new(&heavy.v, Some(&heavy.d), &heavy.j, SimParams::mid())
        .unwrap()
        .reads(150, 11)
        .unwrap();
    let light_reads = Simulator::new(&light_v, None, &light_j, SimParams::mid())
        .unwrap()
        .reads(150, 12)
        .unwrap();
    let mut ws = assigner.workspace();
    let mut d_on_light = 0;
    let mut called = 0;
    for (reads, locus) in [(&heavy_reads, Locus::Igh), (&light_reads, Locus::Igk)] {
        for read in reads {
            let call = assigner.assign(&encode(&read.sequence), &mut ws);
            let Some(path) = &call.path else { continue };
            called += 1;
            for segment in &path.segments {
                let placement = &call.placements[segment.placement];
                let set = assigner.segment(placement.kind).unwrap().set();
                for &allele in &placement.alleles {
                    assert_eq!(
                        set.allele(allele).loci,
                        Loci::of(locus),
                        "{} on a {locus} read",
                        set.allele(allele).name
                    );
                }
                d_on_light += usize::from(locus == Locus::Igk && placement.kind == SegmentKind::D);
            }
        }
    }
    assert_eq!(d_on_light, 0);
    assert!(called >= 290, "{called} of 300 reads called");
}

/// A simulated read joins V, D, J and C of one locus: light-chain reads
/// carry no IGHD and no IGHC, heavy-chain reads no IGKV or IGKJ.
#[test]
fn simulated_reads_stay_within_one_locus() {
    let (heavy, _, light_v, light_j) = databases();
    let joined = |kind, a: &GermlineSet, b: &GermlineSet| {
        GermlineSet::from_fasta(kind, &(fasta(a) + &fasta(b))).unwrap()
    };
    let v = joined(SegmentKind::V, &heavy.v, &light_v);
    let j = joined(SegmentKind::J, &heavy.j, &light_j);
    let sim = Simulator::new(&v, Some(&heavy.d), &j, SimParams::mid())
        .unwrap()
        .with_constant(&heavy.c)
        .unwrap();
    let (mut light, mut heavy_d, mut heavy_c) = (0, 0, 0);
    for read in sim.reads(600, 5).unwrap() {
        let locus = j.allele(read.j as u32).loci;
        assert_eq!(v.allele(read.v as u32).loci, locus);
        if locus == Loci::of(Locus::Igk) {
            light += 1;
            assert_eq!((read.d, read.d2, read.c), (None, None, None));
            continue;
        }
        for d in [read.d, read.d2].into_iter().flatten() {
            assert_eq!(heavy.d.allele(d as u32).loci, locus);
            heavy_d += 1;
        }
        heavy_c += usize::from(read.c.is_some());
    }
    assert!(
        light > 200 && heavy_d > 200 && heavy_c > 150,
        "{light} {heavy_d} {heavy_c}"
    );
    // A J that pairs with no V is refused rather than retried forever.
    assert_eq!(
        Simulator::new(&light_v, None, &heavy.j, SimParams::mid()).err(),
        Some(swig::sim::SimError::Unpaired)
    );
}

/// A read that starts inside CDR3 holds no V, but a few bases of another
/// locus' V can still give a chance V hit. Such a hit neither bounds the J
/// window nor chooses the J's locus: the read's own J is found and called.
#[test]
fn a_chance_v_hit_does_not_choose_the_j_locus() {
    let (heavy, _, light_v, light_j) = databases();
    let v = fasta(&heavy.v) + &fasta(&light_v);
    let j = fasta(&heavy.j) + &fasta(&light_j);
    let d = fasta(&heavy.d);
    let texts = GermlineTexts {
        v: &v,
        d: Some(&d),
        j: &j,
        c: None,
    };
    let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
    let config = CallingProfile::ROptimized.config();
    let assigner =
        DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap();
    let mut ws = assigner.workspace();
    let mut chance_v = 0;
    for (k, light) in light_v.alleles().iter().enumerate().take(20) {
        let heavy_j = heavy.j.alleles().get(k % heavy.j.len()).unwrap();
        let tail = &light.codes[light.codes.len() - 15..];
        let read: Vec<u8> = tail.iter().chain(&heavy_j.codes).copied().collect();
        let call = assigner.assign(&read, &mut ws);
        chance_v += usize::from(call.placements.iter().any(|p| p.kind == SegmentKind::V));
        let path = call.path.as_ref().expect("the whole J anchors a call");
        let j_set = assigner.segment(SegmentKind::J).unwrap().set();
        let called = path
            .segments
            .iter()
            .map(|s| &call.placements[s.placement])
            .find(|p| p.kind == SegmentKind::J)
            .expect("a J call");
        assert_eq!(
            j_set.allele(called.primary_allele()).loci,
            Loci::of(Locus::Igh),
            "read {k}"
        );
    }
    assert!(chance_v > 0, "no read got a chance V hit to test against");
}

/// A read holding two loci, as an scFv does (VH–JH, a linker, then VL–JL),
/// anchors a V of each. Its J search covers both, but each locus keeps its
/// own J alignments: four heavy J alleles of one gene outscore a short
/// light J, and must not crowd it out of the J stage's top three.
#[test]
fn each_anchoring_locus_keeps_its_own_j() {
    let (heavy, light, light_v, light_j) = databases();
    let heavy_j = &heavy.j.alleles()[0];
    let anchor = heavy_j.anchor.expect("synthetic J alleles carry anchors");
    let mut siblings = String::new();
    for (k, at) in [1usize, 3, 5].into_iter().enumerate() {
        assert!(at < anchor as usize, "a SNP before the Trp118 codon");
        let mut codes = heavy_j.codes.clone();
        codes[at] = (codes[at] + 1) % 4;
        let name = heavy_j.name.replace("*01", &format!("*{:02}", k + 5));
        siblings += &format!(">{name} anchor={anchor}\n{}\n", decode(&codes));
    }
    let v = fasta(&heavy.v) + &fasta(&light_v);
    let j = fasta(&heavy.j) + &siblings + &fasta(&light_j);
    let d = fasta(&heavy.d);
    let texts = GermlineTexts {
        v: &v,
        d: Some(&d),
        j: &j,
        c: None,
    };
    let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
    let config = CallingProfile::ROptimized.config();
    let assigner =
        DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap();
    let j_set = assigner.segment(SegmentKind::J).unwrap().set();
    let mut ws = assigner.workspace();
    // Whole V genes of one allele each, so the V stage keeps both V's.
    let single = |set: &GermlineSet| {
        let alleles = |gene| set.alleles().iter().filter(|a| a.gene == gene).count();
        set.alleles()
            .iter()
            .filter(|a| alleles(a.gene) == 1)
            .max_by_key(|a| a.codes.len())
            .unwrap()
            .codes
            .clone()
    };
    let (heavy_v, light_v) = (single(&heavy.v), single(&light.v));
    let linker = encode(b"GGTGGAGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG");
    for light_j in light_j.alleles() {
        let read: Vec<u8> = [
            &heavy_v[..],
            &encode(b"GGTACC"),
            &heavy_j.codes,
            &linker,
            &light_v,
            &encode(b"TC"),
            &light_j.codes[..20],
        ]
        .concat();
        let call = assigner.assign(&read, &mut ws);
        let placed = |locus| {
            call.placements.iter().any(|p| {
                p.kind == SegmentKind::J && j_set.allele(p.primary_allele()).loci == Loci::of(locus)
            })
        };
        assert!(placed(Locus::Igh), "no heavy J with {}", light_j.name);
        assert!(placed(Locus::Igk), "{} lost to the heavy J", light_j.name);
    }
}
