//! Retrieval regressions on reads the full-length presets rarely make, each
//! on the synthetic database with fixed seeds.

use swig::germline::GermlineSet;
use swig::placement::Placement;
use swig::prelude::*;
use swig::sim::{Region, SimParams, SimulatedRead, Simulator, SyntheticDatabase, SyntheticSpec};

fn database() -> SyntheticDatabase {
    SyntheticDatabase::generate(&SyntheticSpec::default(), 11)
}

fn assigner(v: &GermlineSet, db: &SyntheticDatabase, config: &EngineConfig) -> DefaultAssigner {
    DefaultAssigner::from_config(config, v.clone(), Some(db.d.clone()), db.j.clone())
        .expect("synthetic database builds")
}

/// The called placement of `kind`, if the call has one.
fn called(call: &Call, kind: SegmentKind) -> Option<&Placement> {
    let path = call.path.as_ref()?;
    path.segments
        .iter()
        .map(|s| &call.placements[s.placement])
        .find(|p| p.kind == kind)
}

/// Names of the alleles of the called placement of `kind`.
fn names(assigner: &DefaultAssigner, call: &Call, kind: SegmentKind) -> Vec<String> {
    let set = assigner.segment(kind).expect("segment").set();
    called(call, kind).map_or(Vec::new(), |p| {
        p.alleles
            .iter()
            .map(|&a| set.allele(a).name.clone())
            .collect()
    })
}

fn gene(name: &str) -> &str {
    name.split('*').next().unwrap_or(name)
}

/// `read`'s germline-sense bases from `start`, emitted on the read's strand.
fn cut(read: &SimulatedRead, start: usize) -> Vec<u8> {
    let forward = &read.forward[start..];
    if read.reverse_complemented {
        swig::sim::reverse_complement(forward)
    } else {
        forward.to_vec()
    }
}

/// An indel of 9 to 30 bases in V puts the seeds on its two sides on
/// diagonals further apart than one seed cluster spans: the band must
/// span both, and the V must still be called through the indel.
#[test]
fn long_indels_in_v_keep_the_call() {
    let db = database();
    let assigner = assigner(&db.v, &db, &EngineConfig::default());
    let sim = Simulator::new(
        &db.v,
        Some(&db.d),
        &db.j,
        SimParams::easy().with_reverse_complement(0.0),
    )
    .expect("simulator");
    let reads = sim.reads(40, 6).expect("reads");
    let mut filler = swig::sim::Rng::new(9);
    let (mut total, mut right) = (0, 0);
    for read in &reads {
        let span = read.span(Region::V).expect("V");
        let truth = gene(&db.v.allele(read.v as u32).name);
        for len in [9usize, 12, 18, 30] {
            for fraction in [0.3, 0.6] {
                let at = span.start + ((span.end - span.start) as f64 * fraction) as usize;
                let deleted = [&read.forward[..at], &read.forward[at + len..]].concat();
                let inserted: Vec<u8> = read.forward[..at]
                    .iter()
                    .copied()
                    .chain(filler.codes(len).into_iter().map(|c| b"ACGT"[c as usize]))
                    .chain(read.forward[at..].iter().copied())
                    .collect();
                for (query, shift) in [(deleted, -(len as i64)), (inserted, len as i64)] {
                    total += 1;
                    let mut workspace = assigner.workspace();
                    let call = assigner.assign(&encode(&query), &mut workspace);
                    let v_end = called(&call, SegmentKind::V).map(|p| p.alignment.q_end as i64);
                    let gene_right = names(&assigner, &call, SegmentKind::V)
                        .first()
                        .is_some_and(|n| gene(n) == truth);
                    let end_right =
                        v_end.is_some_and(|e| (e - (span.end as i64 + shift)).abs() <= 2);
                    right += usize::from(gene_right && end_right);
                }
            }
        }
    }
    println!("V gene and end right on {right} of {total} reads with a long indel");
    assert!(right * 100 >= total * 98, "{right} of {total}");
}

/// A read that starts inside CDR3 holds no V, but a chance V hit in its
/// 3′ flank or C region can sit downstream of the real J. The J must still
/// be found as often as by a search of the whole read.
#[test]
fn reads_starting_in_cdr3_keep_their_j() {
    let db = database();
    let sim = Simulator::new(
        &db.v,
        Some(&db.d),
        &db.j,
        SimParams::mid().with_reverse_complement(0.0),
    )
    .expect("simulator")
    .with_constant(&db.c)
    .expect("C simulator");
    let reads = sim.reads(1500, 5).expect("reads");
    let queries: Vec<Vec<u8>> = reads
        .iter()
        .map(|read| encode(&cut(read, read.span(Region::V).expect("V").end + 3)))
        .collect();
    let truth: Vec<&str> = reads
        .iter()
        .map(|read| gene(&db.j.allele(read.j as u32).name))
        .collect();
    let right = |config: &EngineConfig| {
        let assigner = assigner(&db.v, &db, config);
        let calls = assign_all(&assigner, &queries, 0);
        calls
            .iter()
            .zip(&truth)
            .filter(|(call, truth)| {
                names(&assigner, call, SegmentKind::J)
                    .first()
                    .is_some_and(|n| gene(n) == **truth)
            })
            .count()
    };
    let default = right(&EngineConfig::default());
    let mut whole = EngineConfig::default();
    whole.j.window.reclaim = u32::MAX;
    let whole = right(&whole);
    println!(
        "J gene right: {default} (whole-read J search {whole}) of {}",
        reads.len()
    );
    assert!(
        default + reads.len() / 200 >= whole,
        "J gene right {default}, whole-read search {whole}"
    );
}
