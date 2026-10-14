//! SwiftIG's plausibility guards (Swig ≤ 0.38), which the engine keeps
//! ([`swig::pipeline::PairGuards`]): implausible and chimeric reads keep
//! their better segment instead of becoming a rearrangement, a D needs a V
//! and a J, a junction needs aligned anchors, and C is called as SwiftIG
//! anchored it, on 12 columns at 90 % identity.

use swig::prelude::*;
use swig::sim::{Rng, SyntheticDatabase, SyntheticSpec};

fn database() -> SyntheticDatabase {
    SyntheticDatabase::generate(&SyntheticSpec::default(), 11)
}

fn random(rng: &mut Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4)]).collect()
}

fn bases(codes: &[u8]) -> Vec<u8> {
    swig::alphabet::decode(codes).into_bytes()
}

/// The AIRR record of `read` under `config`, with the C stage.
fn record(db: &SyntheticDatabase, config: &EngineConfig, read: &[u8]) -> AirrRecord {
    let assigner = db
        .assigner(config)
        .expect("synthetic database builds")
        .with_constant(db.c.clone(), &config.c)
        .expect("synthetic C builds");
    let call = assigner.assign(&encode(read), &mut assigner.workspace());
    AirrRecord::from_call(
        "read",
        &String::from_utf8_lossy(read),
        &call,
        &assigner,
        &[None, None, None],
    )
}

fn unguarded() -> EngineConfig {
    let mut config = EngineConfig::default();
    config.options.guards = None;
    config
}

fn paired(record: &AirrRecord) -> bool {
    record.v.is_some() && record.j.is_some()
}

/// A clean rearrangement keeps its V, D, J, C and junction.
#[test]
fn a_rearrangement_keeps_every_segment() {
    let db = database();
    let (v, d, j, c) = (
        bases(&db.v.allele(0).codes),
        bases(&db.d.allele(0).codes),
        bases(&db.j.allele(0).codes),
        bases(&db.c.allele(0).codes),
    );
    let read = [
        v,
        b"GGTAC".to_vec(),
        d,
        b"CCA".to_vec(),
        j,
        c[..60].to_vec(),
    ]
    .concat();
    let r = record(&db, &EngineConfig::default(), &read);
    assert!(paired(&r) && r.d.is_some() && r.c.is_some());
    assert!(r.regions.junction.is_some() && r.regions.productive.is_some());
}

/// A V fragment that stops 220 bases short of its 3′ end does not pair
/// with a J 20 random bases later: it neither faces the junction nor holds
/// a long split, and its Cys104 is not aligned. Without the guards it does.
#[test]
fn a_v_that_does_not_reach_its_junction_does_not_pair() {
    let db = database();
    let mut rng = Rng::new(3);
    let v = bases(&db.v.allele(0).codes);
    let j = bases(&db.j.allele(0).codes);
    let read = [v[..70].to_vec(), random(&mut rng, 20), j].concat();
    let r = record(&db, &EngineConfig::default(), &read);
    assert!(!paired(&r) && r.d.is_none());
    assert_eq!(r.regions.junction, None);
    assert!(paired(&record(&db, &unguarded(), &read)));
}

/// V and J more than 225 bases apart do not pair; 200 apart, they do.
#[test]
fn v_and_j_beyond_the_junction_span_do_not_pair() {
    let db = database();
    let mut rng = Rng::new(5);
    let v = bases(&db.v.allele(0).codes);
    let j = bases(&db.j.allele(0).codes);
    for (gap, pairs) in [(200, true), (260, false), (400, false)] {
        let read = [v.clone(), random(&mut rng, gap), j.clone()].concat();
        let r = record(&db, &EngineConfig::default(), &read);
        assert_eq!(paired(&r), pairs, "{gap} bases between V and J");
        if !pairs {
            assert!(r.d.is_none() && r.regions.junction.is_none(), "{gap}");
        }
    }
}

/// Eight bases of J after 12 random ones are not a J call (10 columns).
#[test]
fn a_j_of_fewer_than_ten_columns_is_not_called() {
    let db = database();
    let mut rng = Rng::new(7);
    let v = bases(&db.v.allele(0).codes);
    for allele in db.j.alleles() {
        let j = bases(&allele.codes);
        let read = [v.clone(), random(&mut rng, 12), j[..8].to_vec()].concat();
        let r = record(&db, &EngineConfig::default(), &read);
        assert!(r.j.is_none(), "{}", allele.name);
    }
}

/// A D is called only between a V and a J: a read that stops eight bases
/// after its D (too few for a J) or starts eight bases before it (too few
/// for a V) gets no D. Without the guards those reads keep the D.
#[test]
fn a_d_needs_a_v_and_a_j() {
    let db = database();
    let mut rng = Rng::new(9);
    let v = bases(&db.v.allele(0).codes);
    let j = bases(&db.j.allele(0).codes);
    let mut unguarded_d = 0;
    for (k, allele) in db.d.alleles().iter().enumerate().take(8) {
        let d = bases(&allele.codes);
        let v_side = [
            v.clone(),
            random(&mut rng, 6),
            d.clone(),
            random(&mut rng, 8),
        ]
        .concat();
        let j_side = [random(&mut rng, 8), d, random(&mut rng, 6), j.clone()].concat();
        for read in [v_side, j_side] {
            let r = record(&db, &EngineConfig::default(), &read);
            assert!(r.d.is_none(), "D {k}: {:?}", r.d.map(|d| d.call));
            assert!(!paired(&r), "D {k}");
            unguarded_d += usize::from(record(&db, &unguarded(), &read).d.is_some());
        }
    }
    assert!(unguarded_d > 0, "no read tested the rule");
}

/// A V cut six bases before its Cys104 codon pairs with its J, but the
/// junction is not extrapolated through bases the V does not align.
#[test]
fn the_junction_needs_aligned_anchors() {
    let db = database();
    let v = &db.v.allele(0);
    let cys = v.anchor.expect("synthetic V alleles carry anchors") as usize;
    let read = [
        bases(&v.codes[..cys - 6]),
        b"GGTACCAT".to_vec(),
        bases(&db.j.allele(0).codes),
    ]
    .concat();
    let r = record(&db, &EngineConfig::default(), &read);
    assert!(paired(&r));
    assert_eq!(r.regions.junction, None);
    assert_eq!(r.regions.productive, None);
    assert_eq!(r.regions.complete_vdj, None);
    assert!(record(&db, &unguarded(), &read).regions.junction.is_some());
}

/// C is anchored within 7 bases of the J exit, as SwiftIG anchored it:
/// five bases between J and C, a J missing its last four bases or a C
/// missing its first three still give the C.
#[test]
fn c_is_called_up_to_seven_bases_off_the_splice() {
    let db = database();
    let v = bases(&db.v.allele(0).codes);
    let j = bases(&db.j.allele(0).codes);
    let c = bases(&db.c.allele(0).codes);
    let name = &db.c.allele(0).name;
    let n = j.len();
    for (what, read) in [
        (
            "5 bases between J and C",
            [v.clone(), j.clone(), b"GATCA".to_vec(), c[..80].to_vec()].concat(),
        ),
        (
            "J missing 4 bases",
            [v.clone(), j[..n - 4].to_vec(), c[..80].to_vec()].concat(),
        ),
        (
            "C missing 3 bases",
            [v.clone(), j.clone(), c[3..80].to_vec()].concat(),
        ),
    ] {
        let r = record(&db, &EngineConfig::default(), &read);
        let called = r.c.map(|c| c.call).unwrap_or_default();
        assert!(called.split(',').any(|a| a == name), "{what}: {called:?}");
    }
}

/// C is called on at least 12 columns at 90 % identity: ten bases of C, or
/// 80 bases at one substitution in four, are not called.
#[test]
fn c_needs_twelve_columns_at_ninety_percent() {
    let db = database();
    let v = bases(&db.v.allele(0).codes);
    let j = bases(&db.j.allele(0).codes);
    let c = bases(&db.c.allele(0).codes);
    let short = [v.clone(), j.clone(), c[..10].to_vec()].concat();
    assert!(record(&db, &EngineConfig::default(), &short).c.is_none());
    let mut diverged = c[..80].to_vec();
    for base in diverged.iter_mut().skip(1).step_by(4) {
        *base = if *base == b'A' { b'C' } else { b'A' };
    }
    let read = [v, j, diverged].concat();
    assert!(record(&db, &EngineConfig::default(), &read).c.is_none());
    let mut config = EngineConfig::default();
    config.c.accept.min_identity = 0.0;
    assert!(record(&db, &config, &read).c.is_some());
}

/// A V fragment too short for a long split (under 80 columns) does not
/// pair with a J that is not significant on its own, unless both conserved
/// anchors are aligned (SwiftIG's motif rescue): a chimera's 40 V bases do
/// not vouch for a heavily mutated J fragment that starts after its
/// Trp/Phe118 codon. The 8 bases between them differ from the J germline
/// before the fragment, so no chance match aligns the J's anchor.
#[test]
fn a_short_v_does_not_vouch_for_a_chance_j() {
    let db = database();
    let other = |b: &u8| match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    let v = bases(&db.v.allele(0).codes);
    let mut checked = 0;
    for allele in db.j.alleles() {
        let trp = allele.anchor.expect("synthetic J alleles carry anchors") as usize;
        let germline = bases(&allele.codes);
        let spacer: Vec<u8> = germline[trp - 5..trp + 3].iter().map(other).collect();
        let mut j = germline[trp + 3..].to_vec();
        for base in j.iter_mut().skip(2).step_by(4) {
            *base = other(base);
        }
        let read = [v[v.len() - 40..].to_vec(), spacer, j].concat();
        let unguarded = record(&db, &unguarded(), &read);
        if !paired(&unguarded) {
            continue;
        }
        checked += 1;
        let r = record(&db, &EngineConfig::default(), &read);
        assert!(!paired(&r), "{}: {:?}", allele.name, r.j.map(|j| j.call));
        assert_eq!(r.regions.productive, None);
    }
    assert!(checked > 0, "no read paired without the guards");
}
