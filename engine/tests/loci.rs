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
