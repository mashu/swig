//! Germline input as Swig ships and writes it: `SWIGMETA` coordinates,
//! prepared databases read back, and the loci of alleles, end to end from
//! FASTA text to AIRR rows.

use swig::prelude::*;

/// Human TRAV1-2*01 with Swig's shipped SWIGMETA (Cys104 at 261).
const TRAV1_2: &str = ">TRAV1-2*01 SWIGMETA=0,-1,0,75,75,93,93,144,144,162,162,264,1
GGACAAAACATTGACCAGCCCACTGAGATGACAGCTACGGAAGGTGCCATTGTCCAGATCAACTGCACGTACCAGACATCTGGGTTCAACGGGCTGTTCTGGTACCAGCAACATGCTGGCGAAGCACCCACATTTCTGTCTTACAATGTTCTGGATGGTTTGGAGGAGAAAGGTCGTTTTTCTTCATTCCTTAGTCGGTCTAAAGGGTACAGTTACCTCCTTTTGAAGGAGCTCCAGATGAAAGACTCTGCCTCTTACCTCTGTGCTGTGAGAGA
";

/// Human TRAJ48*01 and TRAJ24*01 with Swig's shipped SWIGMETA, whose
/// anchors (8: `FGNE`, 14: `WGKF`) precede each allele's `FGXG` (29).
const TRAJ: &str = ">TRAJ48*01 SWIGMETA=2,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4
TATCTAACTTTGGAAATGAGAAATTAACCTTTGGGACTGGAACAAGACTCACCATCATACCCA
>TRAJ24*01 SWIGMETA=2,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4
TGACAACTGACAGCTGGGGGAAATTCGAGTTTGGAGCAGGGACCCAGGTTGTGGTCACCCCAG
";

/// The first 59 bases of human TRAC*01.
const TRAC: &str = ">TRAC*01
ATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAG
";

/// Rows of `reads` (FASTA) against the database `texts`, as maps from
/// column name to value.
fn assign(texts: GermlineTexts<'_>, reads: &str) -> Vec<std::collections::HashMap<String, String>> {
    let options = RunOptions::default();
    let database = PreparedDatabase::from_fasta(texts, options.annotate()).unwrap();
    let mut run = Run::new(&options, database).unwrap();
    let mut out = Vec::new();
    run.write_airr(reads.as_bytes(), &mut out, None).unwrap();
    let text = String::from_utf8(out).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

/// Productive TRAV1-2/TRAJ48 and TRAV1-2/TRAJ24 rearrangements spliced to
/// TRAC: the junction ends at each J's `FGXG`, not at the shipped SWIGMETA
/// anchor upstream of it (which cut CDR3 short by five and seven codons).
#[test]
fn junctions_end_at_the_j_motif_a_swigmeta_anchor_misplaces() {
    let v = TRAV1_2.lines().nth(1).unwrap();
    let c = TRAC.lines().nth(1).unwrap();
    let mut reads = String::new();
    for (j, name) in TRAJ.lines().skip(1).step_by(2).zip(["TRAJ48", "TRAJ24"]) {
        reads += &format!(">{name}\n{v}GGAAGTC{}{c}\n", &j[2..]);
    }
    let texts = GermlineTexts {
        v: TRAV1_2,
        d: None,
        j: TRAJ,
        c: Some(TRAC),
    };
    let rows = assign(texts, &reads);
    let expected = [
        ("TRAJ48*01", "CAVREEVSNFGNEKLTF", "FGTGTRLTIIP"),
        ("TRAJ24*01", "CAVREEVTTDSWGKFEF", "FGAGTQVVVTP"),
    ];
    for (row, (j_call, junction, fwr4)) in rows.iter().zip(expected) {
        assert_eq!(row["j_call"], j_call);
        assert_eq!(row["junction_aa"], junction, "{j_call}");
        assert_eq!(row["fwr4_aa"], fwr4, "{j_call}");
        assert_eq!(row["productive"], "T", "{j_call}");
        assert_eq!(row["c_call"], "TRAC*01", "{j_call}");
        // The anchor is the motif's, so the SWIGMETA provenance is not reported.
        assert_eq!(row["j_annotation_source"], "motif", "{j_call}");
        assert_eq!(row["v_annotation_source"], "IMGT-gapped");
    }
}

/// Kappa V alleles labelled without their locus code (`KV…`) in a heavy +
/// kappa database, whose headers read "human kappa light chain": `LIGHT`
/// holds `IGH` but names no locus. The kappa reads keep kappa J calls and
/// get no D, as with the IMGT names.
#[test]
fn header_words_name_no_locus() {
    use swig::germline::GermlineSet;
    use swig::locus::{Loci, Locus};
    use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};
    let heavy = SyntheticDatabase::generate(&SyntheticSpec::default(), 1);
    let light = SyntheticDatabase::generate(&SyntheticSpec::default(), 2);
    let fasta = |set: &GermlineSet, rename: &dyn Fn(&str) -> String| {
        set.alleles()
            .iter()
            .map(|a| format!(">{}\n{}\n", rename(&a.name), decode(&a.codes)))
            .collect::<String>()
    };
    let same = |name: &str| name.to_string();
    let light_v = fasta(&light.v, &|name| {
        format!("{} human kappa light chain", name.replace("IGHV", "KV"))
    });
    let light_j = fasta(&light.j, &|name| name.replace("IGHJ", "IGKJ"));
    let v = fasta(&heavy.v, &same) + &light_v;
    let j = fasta(&heavy.j, &same) + &light_j;
    let d = fasta(&heavy.d, &same);
    let texts = GermlineTexts {
        v: &v,
        d: Some(&d),
        j: &j,
        c: None,
    };
    let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
    let kappa_v = database
        .v
        .alleles()
        .iter()
        .filter(|a| a.name.starts_with("KV"));
    assert!(kappa_v.clone().count() > 0);
    assert!(kappa_v.into_iter().all(|a| a.loci == Loci::ANY));
    let reads_v = GermlineSet::from_fasta(SegmentKind::V, &light_v).unwrap();
    let reads_j = GermlineSet::from_fasta(SegmentKind::J, &light_j).unwrap();
    let reads = Simulator::new(&reads_v, None, &reads_j, SimParams::mid())
        .unwrap()
        .reads(150, 12)
        .unwrap();
    let config = CallingProfile::ROptimized.config();
    let assigner =
        DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap();
    let mut ws = assigner.workspace();
    let (mut kappa_j, mut d_calls) = (0, 0);
    for read in &reads {
        let call = assigner.assign(&encode(&read.sequence), &mut ws);
        let Some(path) = &call.path else { continue };
        for segment in &path.segments {
            let placement = &call.placements[segment.placement];
            let set = assigner.segment(placement.kind).unwrap().set();
            let loci = set.allele(placement.alleles[0]).loci;
            kappa_j +=
                usize::from(placement.kind == SegmentKind::J && loci == Loci::of(Locus::Igk));
            d_calls += usize::from(placement.kind == SegmentKind::D);
        }
    }
    assert_eq!(d_calls, 0);
    assert!(kappa_j >= 145, "{kappa_j} of 150 reads called a kappa J");
}
