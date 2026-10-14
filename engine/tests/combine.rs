//! The pair guards with the open-flank cost (Swig 0.39.5): a read that
//! shows both a V and a J keeps every guard of 0.39.4, but a D at the read's
//! open end, after a V with no J behind it or before a J with no V ahead of
//! it, is called and pays for the J or V the read would show beyond it
//! ([`swig::config::OpenFlank`]). Swig 0.39.4 called no such D. These tests
//! use only the API that 0.39.4 has, so they run against it unchanged.

use swig::prelude::*;
use swig::sim::{
    Region, Rng, SimParams, SimulatedRead, Simulator, SyntheticDatabase, SyntheticSpec,
};

fn database() -> SyntheticDatabase {
    SyntheticDatabase::generate(&SyntheticSpec::default(), 11)
}

fn unguarded() -> EngineConfig {
    let mut config = EngineConfig::default();
    config.options.guards = None;
    config
}

/// Germline-sense simulated reads with a D of at least 12 bases.
fn simulated(db: &SyntheticDatabase, count: usize, seed: u64) -> Vec<SimulatedRead> {
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::mid()).unwrap();
    sim.reads(count * 4, seed)
        .unwrap()
        .into_iter()
        .filter(|r| !r.reverse_complemented && r.span(Region::D).is_some_and(|d| d.len() >= 12))
        .take(count)
        .collect()
}

/// The AIRR record of `read` under `config`.
fn record(db: &SyntheticDatabase, config: &EngineConfig, read: &[u8]) -> AirrRecord {
    let assigner = db.assigner(config).expect("synthetic database builds");
    let call = assigner.assign(&encode(read), &mut assigner.workspace());
    AirrRecord::from_call(
        "read",
        &String::from_utf8_lossy(read),
        &call,
        &assigner,
        &[None, None, None],
    )
}

/// Of `reads` (each with its true D allele), those called with a D, and
/// those whose D call names the true allele.
fn d_calls(
    db: &SyntheticDatabase,
    config: &EngineConfig,
    reads: &[(Vec<u8>, usize)],
) -> (usize, usize) {
    let (mut called, mut right) = (0, 0);
    for (read, truth) in reads {
        let r = record(db, config, read);
        if let Some(d) = r.d {
            called += 1;
            let name = &db.d.alleles()[*truth].name;
            right += usize::from(d.call.split(',').any(|a| a == name));
        }
    }
    (called, right)
}

/// A read that ends two bases after its D (its V, N1 and D, and no J)
/// keeps its D: the path V–D pays nothing within the N/P allowance (55 of
/// 60 right). Swig 0.39.4 called no D on any of them (a D only between a V
/// and a J).
#[test]
fn a_read_cut_at_the_end_of_its_d_keeps_it() {
    let db = database();
    let reads: Vec<(Vec<u8>, usize)> = simulated(&db, 60, 21)
        .into_iter()
        .map(|r| {
            let end = (r.span(Region::D).unwrap().end + 2).min(r.forward.len());
            (r.forward[..end].to_vec(), r.d.unwrap())
        })
        .collect();
    assert!(reads.len() >= 50, "{}", reads.len());
    let (called, right) = d_calls(&db, &EngineConfig::default(), &reads);
    assert!(
        right * 10 >= reads.len() * 8,
        "{called} D calls, {right} right of {}",
        reads.len()
    );
    // As many as without the guards (0.39.3's rule with the cost).
    let (_, free) = d_calls(&db, &unguarded(), &reads);
    assert!(
        right + 2 >= free,
        "{right} right D calls, {free} without the guards"
    );
}

/// A read that starts ten bases before its D (no V) keeps its D: the path
/// D–J pays nothing within the allowance (58 of 60 right). Swig 0.39.4
/// called no D on any.
#[test]
fn a_read_starting_after_its_v_keeps_its_d() {
    let db = database();
    let reads: Vec<(Vec<u8>, usize)> = simulated(&db, 60, 22)
        .into_iter()
        .map(|r| {
            let start = r.span(Region::D).unwrap().start.saturating_sub(10);
            (r.forward[start..].to_vec(), r.d.unwrap())
        })
        .collect();
    assert!(reads.len() >= 50, "{}", reads.len());
    let (called, right) = d_calls(&db, &EngineConfig::default(), &reads);
    assert!(
        right * 10 >= reads.len() * 8,
        "{called} D calls, {right} right of {}",
        reads.len()
    );
}

/// A full V followed by 100 random bases (no D, no J) gets few Ds: a
/// chance D far from the V pays for the J the read would show after it.
/// The unguarded engine without the cost (0.39.3) calls one on most such
/// reads, the cost on 3 of 60; Swig 0.39.4 called none, by forbidding
/// every D without a J.
#[test]
fn a_v_followed_by_100_random_bases_gets_few_ds() {
    let db = database();
    let mut rng = Rng::new(23);
    let reads: Vec<Vec<u8>> = simulated(&db, 60, 23)
        .into_iter()
        .map(|r| {
            let v_end = r.span(Region::V).unwrap().end;
            [r.forward[..v_end].to_vec(), rng.bases(100)].concat()
        })
        .collect();
    let called = reads
        .iter()
        .filter(|read| record(&db, &EngineConfig::default(), read).d.is_some())
        .count();
    assert!(
        called * 10 <= reads.len(),
        "{called} D calls of {}",
        reads.len()
    );
}

/// A read whose J fails the pair guards: a full V, five bases and a D,
/// then random bases up to 230 bases after the V, and a J. V and J lie
/// beyond the 225-base junction span, so they do not pair; the unguarded
/// engine calls V, D and J. Under the guards the J gives way and the read
/// keeps its V and D, the D paying the open-flank cost over the bases
/// after it. Swig 0.39.4 kept the V alone.
#[test]
fn a_j_that_does_not_pair_gives_way_to_the_v_and_its_d() {
    let db = database();
    let mut rng = Rng::new(24);
    let v = bases(&db.v.allele(0).codes);
    let j = bases(&db.j.allele(0).codes);
    let mut kept = 0;
    for allele in
        db.d.alleles()
            .iter()
            .filter(|a| a.codes.len() >= 20)
            .take(8)
    {
        let d = bases(&allele.codes);
        let read = [
            v.clone(),
            b"GGTAC".to_vec(),
            d.clone(),
            rng.bases(225 - d.len()),
            j.clone(),
        ]
        .concat();
        let free = record(&db, &unguarded(), &read);
        assert!(free.v.is_some() && free.j.is_some(), "{}", allele.name);
        let r = record(&db, &EngineConfig::default(), &read);
        assert!(r.v.is_some() && r.j.is_none(), "{}", allele.name);
        assert_eq!(r.regions.junction, None, "{}", allele.name);
        if let Some(called) = r.d {
            assert!(
                called.call.split(',').any(|a| a == allele.name),
                "{}",
                allele.name
            );
            kept += 1;
        }
    }
    assert!(kept >= 6, "{kept} of 8 reads kept their D");
}

/// A D at the open end lies where a read without the dropped partner
/// would have searched it: within the D window's reach (60 bases) of the
/// V's end, or of the J's start. Reads of a full V, random bases, a D and
/// a J 240 bases after the V (beyond the 225-base junction span, so the
/// guards drop the pair): with the D 100 to 160 bases after the V and more
/// than 60 bases before the J, the call is the V (or the J) alone, as on
/// the same read without the J and in Swig 0.39.4, not the V with a D that
/// only the dropped J's window found (a chance D in the first 60 random
/// bases, at most 4 of 24, may be); with the D five bases before the J,
/// the call is never the V with that D (D–J, or the V alone). Before the
/// reach was checked, every read of the first kind was called V–D with a
/// D presence of 1.
#[test]
fn an_open_end_d_lies_within_reach_of_its_partner() {
    let db = database();
    let mut rng = Rng::new(31);
    let v = bases(&db.v.allele(0).codes);
    let j = bases(&db.j.allele(0).codes);
    let config = EngineConfig::default();
    let (mut far, mut chance, mut near_j, mut d_j) = (0, 0, 0, 0);
    for allele in
        db.d.alleles()
            .iter()
            .filter(|a| a.codes.len() >= 20)
            .take(8)
    {
        let d = bases(&allele.codes);
        for n1 in [100, 130, 160] {
            let tail = 240 - n1 - d.len();
            let read = [
                v.clone(),
                rng.bases(n1),
                d.clone(),
                rng.bases(tail),
                j.clone(),
            ]
            .concat();
            let r = record(&db, &config, &read);
            assert!(r.v.is_some() || r.j.is_some(), "{} n1 {n1}", allele.name);
            // A chance D in the first 60 random bases may be called, as on
            // the read without the J; the planted D, 100 or more bases on,
            // may not.
            if let Some(called) = &r.d {
                assert!(
                    called.sequence_end as usize <= v.len() + 60 + 2,
                    "{} n1 {n1}: D {} at {}..{} beyond the reach",
                    allele.name,
                    called.call,
                    called.sequence_start,
                    called.sequence_end
                );
                chance += 1;
            } else {
                assert!(
                    r.d_presence_posterior.is_none_or(|p| p < 0.5),
                    "{} n1 {n1}: D presence {:?}",
                    allele.name,
                    r.d_presence_posterior
                );
            }
            far += 1;
        }
        let n1 = 240 - 5 - d.len();
        let read = [v.clone(), rng.bases(n1), d.clone(), rng.bases(5), j.clone()].concat();
        let r = record(&db, &config, &read);
        if r.d.is_some() {
            assert!(
                r.j.is_some() && r.v.is_none(),
                "{}: V with a far D",
                allele.name
            );
            d_j += 1;
        }
        near_j += 1;
    }
    assert!(
        far == 24 && near_j == 8 && chance <= 4,
        "{far} far ({chance} with a chance D near the V), {near_j} near the J"
    );
    eprintln!("{d_j} of {near_j} reads with the D beside a dropped J called D-J");
}

fn bases(codes: &[u8]) -> Vec<u8> {
    swig::alphabet::decode(codes).into_bytes()
}
