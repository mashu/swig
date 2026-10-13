//! Germline input as Swig ships and writes it: `SWIGMETA` coordinates,
//! prepared databases read back, and the loci of alleles, end to end from
//! FASTA text to AIRR rows.

use swig::prelude::*;

/// Human TRAV1-2*01 with Swig's shipped SWIGMETA (Cys104 at 261).
const TRAV1_2: &str = ">TRAV1-2*01 SWIGMETA=0,-1,0,75,75,93,93,144,144,162,162,264,1
GGACAAAACATTGACCAGCCCACTGAGATGACAGCTACGGAAGGTGCCATTGTCCAGATCAACTGCACGTACCAGACATCTGGGTTCAACGGGCTGTTCTGGTACCAGCAACATGCTGGCGAAGCACCCACATTTCTGTCTTACAATGTTCTGGATGGTTTGGAGGAGAAAGGTCGTTTTTCTTCATTCCTTAGTCGGTCTAAAGGGTACAGTTACCTCCTTTTGAAGGAGCTCCAGATGAAAGACTCTGCCTCTTACCTCTGTGCTGTGAGAGA
";

/// IMGT gap runs `(bases before, run)` that put TRAV1-2*01 on IMGT columns
/// (one FR1 gap codon at IMGT 10, CDR1 and CDR2 gaps in the middle, FR3
/// gaps at 73 and 81–84): Cys23, Trp41 and Cys104 at IMGT 23, 41 and 104.
const TRAV1_2_GAPS: [(usize, usize); 5] = [(27, 3), (84, 18), (153, 12), (183, 3), (204, 12)];

/// `sequence` with the gap runs `gaps` inserted.
fn imgt_gapped(sequence: &str, gaps: &[(usize, usize)]) -> String {
    let mut out = String::new();
    let mut at = 0;
    for &(before, run) in gaps {
        out += &sequence[at..before];
        out += &".".repeat(run);
        at = before;
    }
    out + &sequence[at..]
}

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

/// A prepared database read back keeps what its annotation knew: values
/// annotation derived stay as trusted as they were (a transferred anchor is
/// no donor and still counts as uncertain), doubtful starts stay flagged,
/// every allele reports the provenance it had, and what the IMGT numbering
/// of an ungapped copy no longer shows stays known: the region starts
/// before a 5′-partial sequence stay empty, and a 3′-partial one stays
/// without Cys104.
#[test]
fn a_prepared_database_reloads_with_its_annotation() {
    let v = TRAV1_2.lines().nth(1).unwrap();
    let mut broken = v.as_bytes().to_vec();
    broken[99..102].copy_from_slice(b"TTG"); // Trp41
    broken[261..264].copy_from_slice(b"AGT"); // Cys104
    let broken = String::from_utf8(broken).unwrap();
    let gapped = imgt_gapped(v, &TRAV1_2_GAPS);
    // *04 starts inside CDR1 (IMGT column 85), *05 ends inside FR3 (270).
    let late = ".".repeat(84) + &gapped[84..];
    let v_text = format!(
        "{TRAV1_2}>TRAV1-2*02\n{broken}\n>TRAV1-2*03\n{}\n>TRAV1-2*04\n{late}\n>TRAV1-2*05\n{}\n",
        &v[..240],
        &gapped[..270]
    );
    let prepare = |v: &str, j: &str| {
        let texts = GermlineTexts {
            v,
            d: None,
            j,
            c: None,
        };
        PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap()
    };
    let first = prepare(&v_text, TRAJ);
    let labels = |db: &PreparedDatabase| -> Vec<(String, [bool; 5], Option<String>)> {
        db.annotation
            .v
            .iter()
            .chain(&db.annotation.j)
            .zip(db.v.alleles().iter().chain(db.j.alleles()))
            .map(|(a, allele)| (a.source.label(), a.uncertain, allele.source.clone()))
            .collect()
    };
    let partial_notes = |db: &PreparedDatabase| db.annotation.v[4].notes_field();
    let expected = labels(&first);
    assert_eq!(expected[1].0, "transfer:TRAV1-2*01");
    assert_eq!(
        (expected[2].0.as_str(), expected[2].1),
        ("missing", [true; 5])
    );
    assert_eq!(expected[0].2.as_deref(), Some("IMGT-gapped"));
    let lost = 81; // FR1 and the first two CDR1 codons
    let late = first.v.allele(3);
    assert_eq!(expected[3].0, "imgt-gaps");
    assert_eq!(late.anchor, Some(261 - lost));
    let starts = [
        None,
        None,
        Some(93 - lost),
        Some(144 - lost),
        Some(162 - lost),
    ];
    assert_eq!(late.regions.values(), starts);
    assert_eq!(first.v.allele(4).anchor, None);
    let note = "3′-partial: the IMGT numbering ends before Cys104";
    assert_eq!(partial_notes(&first), note);
    let files: std::collections::HashMap<&str, String> = first.files().into_iter().collect();
    let again = prepare(&files["V.fasta"], &files["J.fasta"]);
    assert_eq!(labels(&again), expected);
    assert_eq!(partial_notes(&again), note);
    assert_eq!(again.annotation.uncertain(), first.annotation.uncertain());
    assert_eq!(again.warnings(3), first.warnings(3));
    for (allele, before) in again.v.alleles().iter().zip(first.v.alleles()) {
        assert_eq!(allele.anchor, before.anchor, "{}", allele.name);
        assert_eq!(allele.regions, before.regions, "{}", allele.name);
    }
}
