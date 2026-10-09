//! Degenerate and adversarial reads: the pipeline must neither panic nor
//! produce malformed output, and must call nothing it cannot support.

use vdj_assign::prelude::*;
use vdj_assign::sim::{C_STREAM, Rng, SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

fn database() -> SyntheticDatabase {
    SyntheticDatabase::generate(&SyntheticSpec::default(), 11)
}

fn assigner(db: &SyntheticDatabase, config: &EngineConfig) -> DefaultAssigner {
    DefaultAssigner::from_config(config, db.v.clone(), Some(db.d.clone()), db.j.clone())
        .expect("synthetic database builds")
}

/// [`assigner`] with the C stage.
fn constant_assigner(db: &SyntheticDatabase, config: &EngineConfig) -> DefaultAssigner {
    assigner(db, config)
        .with_constant(db.c.clone(), &config.c)
        .expect("synthetic C builds")
}

fn random(rng: &mut Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4)]).collect()
}

/// Every AIRR record of an assigner with a C stage is well formed: one
/// field per column, coordinates inside the read, and spans consistent with
/// their alignment strings.
fn check_record(record: &AirrRecord, query_len: usize) {
    let columns = AirrRecord::tsv_header_with(true).split('\t').count();
    assert_eq!(record.to_tsv().split('\t').count(), columns);
    for segment in [&record.v, &record.d, &record.d2, &record.j, &record.c]
        .into_iter()
        .flatten()
    {
        assert!(segment.sequence_start >= 1);
        assert!(segment.sequence_end as usize <= query_len);
        assert!(segment.sequence_start <= segment.sequence_end);
        assert_eq!(
            segment.sequence_alignment.len(),
            segment.germline_alignment.len()
        );
        let read_bases = segment
            .sequence_alignment
            .bytes()
            .filter(|&b| b != b'-')
            .count();
        assert_eq!(
            read_bases,
            (segment.sequence_end - segment.sequence_start + 1) as usize
        );
        if let Some(p) = segment.posterior {
            assert!((0.0..=1.0).contains(&p), "posterior {p}");
        }
    }
    if let Some(p) = record.d_presence_posterior {
        assert!((0.0..=1.0 + 1e-12).contains(&p), "D presence {p}");
    }
}

fn degenerate_reads(db: &SyntheticDatabase) -> Vec<(String, Vec<u8>)> {
    let mut rng = Rng::new(17);
    let v = vdj_assign::alphabet::decode(&db.v.allele(0).codes).into_bytes();
    let j = vdj_assign::alphabet::decode(&db.j.allele(0).codes).into_bytes();
    let c = vdj_assign::alphabet::decode(&db.c.allele(0).codes).into_bytes();
    let reverse = |s: &[u8]| vdj_assign::sim::reverse_complement(s);
    let vj = [v.clone(), j.clone()].concat();
    let vjc = [vj.clone(), c.clone()].concat();
    let mut reads = vec![
        ("empty".to_string(), Vec::new()),
        ("one base".into(), b"A".to_vec()),
        ("shorter than any seed".into(), b"ACGTAC".to_vec()),
        ("all N".into(), vec![b'N'; 500]),
        ("IUPAC and lower case".into(), b"acgtRYKMSWbdhvn".repeat(30)),
        ("poly-A".into(), vec![b'A'; 700]),
        ("dinucleotide repeat".into(), b"CA".repeat(400)),
        ("random 2 kb".into(), random(&mut rng, 2000)),
        ("J only".into(), j.clone()),
        ("J only, reverse complement".into(), reverse(&j)),
        ("V only".into(), v.clone()),
        ("V twice (chimera)".into(), [v.clone(), v.clone()].concat()),
        ("J before V".into(), [j.clone(), v.clone()].concat()),
        ("V with an N every 7 bases".into(), {
            let mut s = v.clone();
            s.iter_mut().step_by(7).for_each(|b| *b = b'N');
            s
        }),
        (
            "V, then 20 kb of random sequence".into(),
            [v.clone(), random(&mut rng, 20_000)].concat(),
        ),
        ("C only".into(), c.clone()),
        ("V, J".into(), vj.clone()),
        (
            "V, J, one base".into(),
            [vj.clone(), b"A".to_vec()].concat(),
        ),
        ("V, J, all N".into(), [vj.clone(), vec![b'N'; 60]].concat()),
        ("V, J, C".into(), vjc.clone()),
        ("V, J, C, reverse complement".into(), reverse(&vjc)),
        (
            "V, J, C, then 20 kb of random sequence".into(),
            [vjc.clone(), random(&mut rng, 20_000)].concat(),
        ),
        // Net drift beyond the splice spread: the alignment stops there.
        (
            "V, J, C with 6 bases deleted".into(),
            [vj.clone(), c[..12].to_vec(), c[18..].to_vec()].concat(),
        ),
    ];
    // Every prefix and suffix length of a simulated read through C.
    let params = SimParams {
        include_c: 1.0,
        c_whole_flank: 1.0,
        ..SimParams::mid()
    };
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params)
        .and_then(|sim| sim.with_constant(&db.c))
        .expect("simulator");
    let mut read = sim.read(&mut rng).expect("read");
    sim.append_constant(&mut read, &mut Rng::new(17 ^ C_STREAM));
    assert!(read.c.is_some());
    let read = read.sequence;
    for len in (0..=read.len()).step_by(7) {
        reads.push((format!("prefix {len}"), read[..len].to_vec()));
        reads.push((format!("suffix {len}"), read[read.len() - len..].to_vec()));
    }
    reads
}

#[test]
fn degenerate_reads_are_handled_by_every_configuration() {
    let db = database();
    let reads = degenerate_reads(&db);
    let queries: Vec<Vec<u8>> = reads.iter().map(|(_, s)| encode(s)).collect();
    let mut configs = vec![EngineConfig::default()];
    for profile in CallingProfile::ALL {
        let mut config = profile.config();
        config.grammar = VdjGrammar::TANDEM_D;
        config.posterior = profile == CallingProfile::SensitiveD;
        configs.push(config);
    }
    for config in &configs {
        let assigner = constant_assigner(&db, config);
        let calls = assign_all(&assigner, &queries, 2);
        let mut constant = std::collections::HashMap::new();
        for ((name, sequence), call) in reads.iter().zip(&calls) {
            assert_eq!(call.query.len(), sequence.len(), "{name}");
            let record = AirrRecord::from_call(
                name,
                &String::from_utf8_lossy(sequence),
                call,
                &assigner,
                &[None, None, None],
            );
            check_record(&record, sequence.len());
            constant.insert(name.as_str(), record.c.map(|c| c.call));
        }
        // C needs a J before it and C bases after the J.
        for name in ["C only", "V, J", "V, J, one base", "V, J, all N"] {
            assert_eq!(constant[name], None, "{name}");
        }
        let names = constant["V, J, C"].clone().expect("C is called");
        assert!(names.split(',').any(|n| n == db.c.allele(0).name));
        for name in [
            "V, J, C, reverse complement",
            "V, J, C, then 20 kb of random sequence",
        ] {
            assert_eq!(constant[name].as_ref(), Some(&names), "{name}");
        }
        // Reads without germline content are never called.
        for (name, call) in reads.iter().zip(&calls).filter(|((name, _), _)| {
            [
                "empty",
                "one base",
                "shorter than any seed",
                "all N",
                "poly-A",
                "dinucleotide repeat",
                "random 2 kb",
            ]
            .contains(&name.as_str())
        }) {
            assert!(call.path.is_none(), "{} was called", name.0);
        }
    }
}

#[test]
fn a_lone_j_is_called_on_either_strand() {
    let db = database();
    let assigner = assigner(&db, &EngineConfig::default());
    let j = vdj_assign::alphabet::decode(&db.j.allele(0).codes).into_bytes();
    let mut ws = assigner.workspace();
    for (sequence, reverse) in [
        (j.clone(), false),
        (vdj_assign::sim::reverse_complement(&j), true),
    ] {
        let call = assigner.assign(&encode(&sequence), &mut ws);
        assert_eq!(call.orientation.is_reverse(), reverse);
        let path = call.path.expect("a significant J anchors a call");
        assert_eq!(path.segments.len(), 1);
        assert_eq!(
            call.placements[path.segments[0].placement].kind,
            SegmentKind::J
        );
    }
}

/// Fixed sequences that follow a J in real libraries are not constant
/// regions: after every J allele, no sequencing adapter, homopolymer or
/// splice donor is called C, whatever the C allele and start diagonal.
#[test]
fn fixed_tails_after_a_j_are_not_called_c() {
    let tails: [&[u8]; 7] = [
        b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",  // TruSeq read 1 adapter
        b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",   // TruSeq read 2 adapter
        b"CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",  // Nextera adapter
        &[b'A'; 40],                            // poly-A
        &[b'G'; 40],                            // two-colour no-signal run
        b"GTAAGTGGCCACTCTAGGGCCTTTGTTTTCTGCTA", // J splice donor, intron
        b"GTGAGTCCTCACAACCTCTCTCCTGCTTTAACTCT",
    ];
    let db = database();
    let assigner = constant_assigner(&db, &EngineConfig::default());
    let mut ws = assigner.workspace();
    let v = vdj_assign::alphabet::decode(&db.v.allele(0).codes);
    for j in db.j.alleles() {
        let j_bases = vdj_assign::alphabet::decode(&j.codes);
        for tail in tails {
            let read = format!("{v}GGTAC{j_bases}{}", String::from_utf8_lossy(tail));
            let call = assigner.assign(&encode(read.as_bytes()), &mut ws);
            let path = call.path.as_ref().expect("V and J are called");
            let last = &call.placements[path.segments.last().unwrap().placement];
            assert_eq!(last.kind, SegmentKind::J);
            assert!(call.c.is_none(), "{} + {}", j.name, read.len());
        }
    }
}
