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
