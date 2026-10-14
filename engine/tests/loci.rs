//! One database of several loci, as Swig composes a BCR or TCR study: every
//! rearrangement stays within one locus.

use swig::germline::GermlineSet;
use swig::locus::{Loci, Locus};
use swig::prelude::*;
use swig::sim::{Region, SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

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

/// On a read holding two loci, each locus is adjudicated on its own and the
/// call is the best path over both in the model the call maximises. That
/// choice must not depend on whether posteriors are computed, under either
/// call model: every locus is read at the divergence of the whole read.
#[test]
fn the_called_locus_does_not_depend_on_posteriors() {
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
    let assigner = |calls, posterior| {
        let mut config = CallingProfile::ROptimized.config();
        (config.calls, config.posterior) = (calls, posterior);
        DefaultAssigner::from_config(
            &config,
            database.v.clone(),
            database.d.clone(),
            database.j.clone(),
        )
        .unwrap()
    };
    // A mutated heavy chain, a linker and a less mutated light chain cut
    // from its start, so that their paths come close in score.
    let heavy_reads = Simulator::new(&heavy.v, Some(&heavy.d), &heavy.j, SimParams::mid())
        .unwrap()
        .reads(100, 21)
        .unwrap();
    let light_reads = Simulator::new(&light_v, None, &light_j, SimParams::easy())
        .unwrap()
        .reads(100, 22)
        .unwrap();
    let linker = encode(b"GGTGGAGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG");
    let mut reads: Vec<Vec<u8>> = Vec::new();
    for (h, l) in heavy_reads.iter().zip(&light_reads) {
        for cut in (0..200).step_by(25) {
            let light = encode(&l.forward[cut.min(l.forward.len())..]);
            reads.push([encode(&h.forward), linker.clone(), light].concat());
        }
    }
    let locus = |assigner: &DefaultAssigner, read: &[u8]| {
        let call = assigner.assign(read, &mut assigner.workspace());
        let path = call.path?;
        let last = &call.placements[path.segments.last()?.placement];
        let set = assigner.segment(last.kind).unwrap().set();
        Some(set.allele(last.primary_allele()).loci)
    };
    let flips = [CallModel::PerRead, CallModel::FixedScale].map(|calls| {
        let (with, without) = (assigner(calls, true), assigner(calls, false));
        reads
            .iter()
            .filter(|read| locus(&with, read) != locus(&without, read))
            .count()
    });
    assert_eq!(flips, [0, 0], "of {} reads", reads.len());
}

/// On a read holding two chains, heavy first (VH–D–JH, a linker, then
/// VL–JL), the best V may be the light one. The D window must still lie
/// between the heavy chain's own V and J: between the light V and the
/// heavy J it would be empty, and the heavy D would never be searched.
#[test]
fn each_chain_bounds_its_own_d_window() {
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
    let v_set = assigner.segment(SegmentKind::V).unwrap().set();
    let mut ws = assigner.workspace();
    // A mutated heavy chain, so that the lightly mutated light V often
    // scores higher.
    let heavy_reads = Simulator::new(&heavy.v, Some(&heavy.d), &heavy.j, SimParams::mid())
        .unwrap()
        .reads(200, 31)
        .unwrap();
    let light_reads = Simulator::new(&light_v, None, &light_j, SimParams::easy())
        .unwrap()
        .reads(200, 32)
        .unwrap();
    let linker = encode(b"GGTGGAGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG");
    let (mut light_first, mut lost) = (0, Vec::new());
    let mut windows = Vec::new();
    for (k, (h, l)) in heavy_reads.iter().zip(&light_reads).enumerate() {
        let at = |region| h.regions.iter().position(|&r| r == region);
        let Some(d_start) = at(Region::D) else {
            continue;
        };
        let d_end = h.regions.iter().rposition(|&r| r == Region::D).unwrap() + 1;
        if d_end - d_start < 8 {
            continue;
        }
        let read = [encode(&h.forward), linker.clone(), encode(&l.forward)].concat();
        let call = assigner.assign(&read, &mut ws);
        if call.orientation != Orientation::Forward {
            continue;
        }
        // Reads whose best anchoring V is the light one, with a heavy V
        // anchoring too (on others the V stage kept light alleles only).
        let mut anchors = call
            .placements
            .iter()
            .filter(|p| p.kind == SegmentKind::V && p.anchor.is_some())
            .map(|p| v_set.allele(p.primary_allele()).loci);
        if anchors.next() != Some(Loci::of(Locus::Igk))
            || !anchors.any(|loci| loci == Loci::of(Locus::Igh))
        {
            continue;
        }
        light_first += 1;
        assigner.d_windows(call.query.len(), &call.placements, &mut windows);
        if !windows.iter().any(|w| w.start <= d_start && d_end <= w.end) {
            lost.push(k);
        }
    }
    assert!(light_first > 0, "no read ranks the light V first");
    assert!(
        lost.is_empty(),
        "heavy D outside every window: {lost:?} ({light_first} light-first)"
    );
}

/// Every allele label of the shipped references names exactly its locus:
/// the IMGT pack (its labels listed in `tests/data/imgt-labels.tsv` by
/// `scripts/build-locus-data.mjs`), the IgT, IgZ and IgO constant regions of
/// fish and the platypus included, and the kimdb macaque sets.
#[test]
fn every_shipped_label_names_its_locus() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let listed = std::fs::read_to_string(root.join("tests/data/imgt-labels.tsv")).unwrap();
    let mut checked = 0;
    let mut wrong = Vec::new();
    for line in listed.lines() {
        let mut fields = line.split('\t');
        let (locus, _kind, labels) = (
            fields.next().unwrap(),
            fields.next().unwrap(),
            fields.next().unwrap(),
        );
        let locus = Loci::parse_list(locus).unwrap();
        for label in labels.split(' ') {
            let expected = if label.starts_with("TRAV") && label.contains("/DV") {
                Loci::parse_list("TRA,TRD").unwrap()
            } else {
                locus
            };
            if Loci::infer(label, label) != expected {
                wrong.push(label.to_string());
            }
            checked += 1;
        }
    }
    let kimdb = root.join("../public/references/kimdb-1.1");
    for species in ["Macaca_mulatta", "Macaca_fascicularis"] {
        for (kind, file) in [
            (SegmentKind::V, "V.fasta"),
            (SegmentKind::D, "D.fasta"),
            (SegmentKind::J, "J.fasta"),
        ] {
            let text = std::fs::read_to_string(kimdb.join(species).join("IGH").join(file)).unwrap();
            let set = GermlineSet::from_fasta(kind, &text).unwrap();
            for allele in set.alleles() {
                if allele.loci != Loci::of(Locus::Igh) {
                    wrong.push(allele.name.clone());
                }
                checked += 1;
            }
        }
    }
    assert!(checked > 5000, "{checked} labels");
    assert!(wrong.is_empty(), "{} of {checked}: {wrong:?}", wrong.len());
}

/// `text` with every label starting `from` renamed to start `to` and the
/// written loci dropped: alleles whose names give no locus.
fn rename(text: &str, from: &str, to: &str) -> String {
    text.replace(&format!(">{from}"), &format!(">{to}"))
        .lines()
        .map(|line| match line.split_once(" locus=") {
            Some((head, _)) if line.starts_with('>') => head.to_string(),
            _ => line.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

/// The D calls on `reads` assigned against `database`.
fn d_calls(database: PreparedDatabase, reads: &[swig::sim::SimulatedRead]) -> usize {
    let config = CallingProfile::ROptimized.config();
    let assigner =
        DefaultAssigner::from_config(&config, database.v, database.d, database.j).unwrap();
    let mut ws = assigner.workspace();
    reads
        .iter()
        .filter(|read| {
            let call = assigner.assign(&encode(&read.sequence), &mut ws);
            call.path.is_some_and(|path| {
                path.segments
                    .iter()
                    .any(|s| call.placements[s.placement].kind == SegmentKind::D)
            })
        })
        .count()
}

/// D alleles with custom names in a database of a heavy and a light locus
/// (BCR) or of TRB and TRA (TCR) belong to the locus with D genes, not to
/// every locus: light-chain and TRA reads get no D, heavy and TRB reads
/// keep theirs.
#[test]
fn custom_named_d_alleles_stay_out_of_light_and_tra_reads() {
    let (heavy, _, light_v, light_j) = databases();
    for (with_d, without_d) in [("IGH", "IGK"), ("TRB", "TRA")] {
        let v =
            rename(&fasta(&heavy.v), "IGH", with_d) + &rename(&fasta(&light_v), "IGK", without_d);
        let j =
            rename(&fasta(&heavy.j), "IGH", with_d) + &rename(&fasta(&light_j), "IGK", without_d);
        let d = rename(&fasta(&heavy.d), "IGHD", "Dx");
        let texts = GermlineTexts {
            v: &v,
            d: Some(&d),
            j: &j,
            c: None,
        };
        let database = || PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
        let light_reads = Simulator::new(&light_v, None, &light_j, SimParams::mid())
            .unwrap()
            .reads(200, 41)
            .unwrap();
        let heavy_reads = Simulator::new(&heavy.v, Some(&heavy.d), &heavy.j, SimParams::mid())
            .unwrap()
            .reads(100, 42)
            .unwrap();
        assert_eq!(d_calls(database(), &light_reads), 0, "{without_d} reads");
        let heavy_d = d_calls(database(), &heavy_reads);
        assert!(heavy_d > 60, "{heavy_d} of 100 {with_d} reads with a D");
    }
}

/// J alleles with custom names in a heavy + light database: C follows the
/// locus the read's path shares, so a light read with a heavy constant
/// region behind its J gets no heavy C; heavy reads keep their C calls.
#[test]
fn custom_named_j_alleles_keep_c_in_the_read_locus() {
    let (heavy, light, light_v, light_j) = databases();
    let v = fasta(&heavy.v) + &fasta(&light_v);
    let j = rename(&fasta(&heavy.j), "IGHJ", "hj") + &rename(&fasta(&light_j), "IGKJ", "kj");
    let d = fasta(&heavy.d);
    let c = fasta(&heavy.c);
    let texts = GermlineTexts {
        v: &v,
        d: Some(&d),
        j: &j,
        c: Some(&c),
    };
    let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
    let config = CallingProfile::ROptimized.config();
    let assigner = DefaultAssigner::from_config(&config, database.v, database.d, database.j)
        .unwrap()
        .with_constant(database.c.unwrap(), &config.c)
        .unwrap();
    let mut ws = assigner.workspace();
    // Light V–J rearrangements ending at their J's 3′ end, then 150 nt of
    // a heavy C.
    let constant = &heavy.c.alleles()[0].codes[..150];
    let mut light_c = 0;
    for (k, lv) in light.v.alleles().iter().take(30).enumerate() {
        let lj = &light_j.alleles()[k % light_j.len()];
        let read = [&lv.codes[..], &encode(b"GGTC"), &lj.codes, constant].concat();
        let call = assigner.assign(&read, &mut ws);
        assert!(call.path.is_some(), "light read {k}");
        light_c += usize::from(call.c.is_some());
    }
    assert_eq!(light_c, 0, "light reads with a heavy C");
    let sim = Simulator::new(&heavy.v, Some(&heavy.d), &heavy.j, SimParams::mid())
        .unwrap()
        .with_constant(&heavy.c)
        .unwrap();
    let reads = sim.reads(100, 43).unwrap();
    let heavy_c = reads
        .iter()
        .filter(|read| read.c.is_some())
        .filter(|read| {
            assigner
                .assign(&encode(&read.sequence), &mut ws)
                .c
                .is_some()
        })
        .count();
    let reaching = reads.iter().filter(|read| read.c.is_some()).count();
    assert!(heavy_c * 10 >= reaching * 9, "{heavy_c} of {reaching}");
}

/// C alleles of no named locus and no relative in a database of heavy and
/// light C do not join the light locus' gate: `N_C` of the kappa gate
/// counts its own C genes only, as without them.
#[test]
fn a_c_allele_of_no_named_locus_leaves_other_gates_alone() {
    let (heavy, _, light_v, light_j) = databases();
    let v = fasta(&heavy.v) + &fasta(&light_v);
    let j = fasta(&heavy.j) + &fasta(&light_j);
    // Sequences of their own (not merged into the heavy ones as aliases).
    let mutated = |text: String, at: usize| {
        text.lines()
            .map(|line| {
                let mut line = line.to_string();
                if !line.starts_with('>') {
                    let base = if &line[at..=at] == "A" { "C" } else { "A" };
                    line.replace_range(at..=at, base);
                }
                line + "\n"
            })
            .collect::<String>()
    };
    // One kappa C gene: a heavy C renamed.
    let kappa = mutated(rename(&fasta(&heavy.c), "IGHC1-1", "IGKC1-1"), 10)
        .split('>')
        .filter(|r| r.starts_with("IGKC"))
        .map(|r| format!(">{r}"))
        .collect::<String>();
    let gate = |c: &str| {
        let texts = GermlineTexts {
            v: &v,
            d: None,
            j: &j,
            c: Some(c),
        };
        let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
        let stage = swig::constant::ConstantStage::new(
            database.c.unwrap(),
            &CallingProfile::ROptimized.config().c,
        )
        .unwrap();
        stage
            .locus_significance
            .iter()
            .find(|(locus, _)| *locus == Locus::Igk)
            .map(|(_, gate)| gate.ln_hypotheses)
            .unwrap()
    };
    let named = fasta(&heavy.c) + &kappa;
    // Custom-named C alleles related to no C (random sequences): their
    // locus falls back to the database's, IGH and IGK.
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut base = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        ['A', 'C', 'G', 'T'][(state % 4) as usize]
    };
    let custom = (1..=3).fold(named.clone(), |text, k| {
        let sequence: String = (0..330).map(|_| base()).collect();
        text + &format!(">Cx{k}\n{sequence}\n")
    });
    assert_eq!(gate(&custom), gate(&named));
}

/// A database of one locus whose alleles all have custom names: the
/// embedded panel places its V and J in IGH (its D follows them), so the
/// AIRR `locus` column is filled.
#[test]
fn a_one_locus_database_of_custom_names_reports_its_locus() {
    const IGHV3_23: &str = "GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAAGA";
    const IGHD3_10: &str = "GTATTACTATGGTTCGGGGAGTTATTATAAC";
    const IGHJ4: &str = "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG";
    let v = format!(">sample_v1\n{IGHV3_23}\n");
    let d = format!(">sample_d1\n{IGHD3_10}\n");
    let j = format!(">sample_j1\n{IGHJ4}\n");
    let texts = GermlineTexts {
        v: &v,
        d: Some(&d),
        j: &j,
        c: None,
    };
    let options = RunOptions::default();
    let database = PreparedDatabase::from_fasta(texts, options.annotate()).unwrap();
    assert!(
        database.locus_warnings().is_empty(),
        "{:?}",
        database.locus_warnings()
    );
    let mut run = Run::new(&options, database).unwrap();
    let reads: String = (0..5)
        .map(|k| {
            let n1 = &"GGGACTCC"[..2 + k];
            format!(">r{k}\n{IGHV3_23}{n1}{}TC{IGHJ4}\n", &IGHD3_10[2..])
        })
        .collect();
    let mut out = Vec::new();
    run.write_airr(reads.as_bytes(), &mut out, None).unwrap();
    let text = String::from_utf8(out).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let locus = header.iter().position(|&c| c == "locus").unwrap();
    let rows: Vec<&str> = lines.map(|l| l.split('\t').nth(locus).unwrap()).collect();
    assert_eq!(rows, ["IGH"; 5]);
}

/// The records of `tests/data/tcr-other-species.fa` of `kind` (and of
/// `species`, if given) as FASTA: TCR alleles of species outside the
/// locus panel (`scripts/build-locus-data.mjs others`). With `custom`, each
/// label `TR<locus><rest>` becomes `Q<locus>x<rest>_<k>` (`k` its record
/// number, as some labels recur across species), a name of no locus.
fn other_species(kind: &str, species: Option<&str>, custom: bool) -> String {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(root.join("tests/data/tcr-other-species.fa")).unwrap();
    let mut out = String::new();
    for (k, record) in text.split('>').skip(1).enumerate() {
        let (header, sequence) = record.split_once('\n').unwrap();
        let mut fields = header.split(' ');
        let label = fields.next().unwrap();
        let attribute = |key: &str| {
            header
                .split(' ')
                .find_map(|f| f.strip_prefix(key)?.strip_prefix('='))
        };
        if attribute("kind") != Some(kind)
            || species.is_some_and(|s| attribute("species") != Some(s))
        {
            continue;
        }
        let name = if custom {
            let locus = label[2..3].to_ascii_lowercase();
            format!("Q{locus}x{}_{k}", &label[3..])
        } else {
            label.to_string()
        };
        out += &format!(">{name}\n{}", sequence);
    }
    out
}

/// The locus a custom name of [`other_species`] stands for.
fn true_locus(name: &str) -> Locus {
    let code = format!("TR{}", name[1..2].to_ascii_uppercase());
    Locus::from_name(&code).unwrap()
}

/// TCR alleles of species outside the panel, all custom-named: no V, J or
/// C is placed outside its own locus, and none in an IG locus, so the loci
/// the database falls back to are TCR loci too. The cases are those the
/// panel once misplaced: TRDV5 of rabbit, ferret and sheep and ferret TRAV33
/// (TRA and TRD share V genes; a relative of another species does not tell
/// them apart), and short J alleles whose closest panel relatives are of
/// another locus (TRBJ1-2 of cat, camel and sheep near IGLJ).
#[test]
fn custom_tcr_alleles_of_other_species_stay_in_their_own_loci() {
    let [v, d, j, c] = ["V", "D", "J", "C"].map(|kind| other_species(kind, None, true));
    let texts = GermlineTexts {
        v: &v,
        d: Some(&d),
        j: &j,
        c: Some(&c),
    };
    let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
    let tcr = Loci::parse_list("TRA,TRB,TRD,TRG").unwrap();
    let mut wrong = Vec::new();
    for set in [
        Some(&database.v),
        database.d.as_ref(),
        Some(&database.j),
        database.c.as_ref(),
    ]
    .into_iter()
    .flatten()
    {
        for allele in set.alleles() {
            let own = Loci::of(true_locus(&allele.name));
            if !allele.loci.pairs_with(own) || allele.loci.union(tcr) != tcr {
                wrong.push(format!("{} {:?}", allele.name, allele.loci));
            }
        }
    }
    assert!(wrong.is_empty(), "{wrong:?}");
}

/// A TCR database of a species outside the panel (rabbit), every allele
/// custom-named: TRA and TRG reads get no D, as with the IMGT names, and
/// TRD reads keep the D, J and C calls the IMGT names give them, though
/// their V (rabbit TRDV5) is closest to a human TRAV.
#[test]
fn an_all_custom_tcr_database_of_another_species_keeps_every_locus_whole() {
    let rabbit = Some("Oryctolagus_cuniculus_Thorbecke_inbred");
    let [v, d, j, c] = ["V", "D", "J", "C"].map(|kind| other_species(kind, rabbit, false));
    let [cv, cd, cj, cc] = ["V", "D", "J", "C"].map(|kind| other_species(kind, rabbit, true));
    let assigner = |v: &str, d: &str, j: &str, c: &str| {
        let texts = GermlineTexts {
            v,
            d: Some(d),
            j,
            c: Some(c),
        };
        let database = PreparedDatabase::from_fasta(texts, AnnotateOptions::default()).unwrap();
        let config = CallingProfile::ROptimized.config();
        DefaultAssigner::from_config(&config, database.v, database.d, database.j)
            .unwrap()
            .with_constant(database.c.unwrap(), &config.c)
            .unwrap()
    };
    let (named, custom) = (assigner(&v, &d, &j, &c), assigner(&cv, &cd, &cj, &cc));
    let set = |kind, text: &str| GermlineSet::from_fasta(kind, text).unwrap();
    let (sv, sd, sj, sc) = (
        set(SegmentKind::V, &v),
        set(SegmentKind::D, &d),
        set(SegmentKind::J, &j),
        set(SegmentKind::C, &c),
    );
    let reads = Simulator::new(&sv, Some(&sd), &sj, SimParams::mid())
        .unwrap()
        .with_constant(&sc)
        .unwrap()
        .reads(600, 51)
        .unwrap();
    // The read's D presence, J and C, as IMGT labels (custom ones mapped
    // back by their locus letter: `QdxJ1*01_30` → `TRDJ1*01`).
    let calls = |assigner: &DefaultAssigner, read: &[u8]| {
        let call = assigner.assign(read, &mut assigner.workspace());
        let label = |kind, allele| {
            let name = &assigner.segment(kind).unwrap().set().allele(allele).name;
            let name = name.split('_').next().unwrap().replace('x', "");
            match name.strip_prefix('Q') {
                Some(rest) => format!("TR{}{}", rest[..1].to_ascii_uppercase(), &rest[1..]),
                None => name,
            }
        };
        let path = call.path.as_ref();
        let segment = |kind| {
            path.and_then(|p| {
                p.segments
                    .iter()
                    .map(|s| &call.placements[s.placement])
                    .find(|p| p.kind == kind)
            })
        };
        (
            segment(SegmentKind::D).is_some(),
            segment(SegmentKind::J).map(|p| label(SegmentKind::J, p.primary_allele())),
            call.c.as_ref().map(|k| {
                label(
                    SegmentKind::C,
                    call.constant[k.segment.placement].primary_allele(),
                )
            }),
        )
    };
    let (mut no_d_reads, mut d_on_no_d) = (0, 0);
    let (mut trd, mut trd_d_same, mut trd_jc_differ) = (0, 0, Vec::new());
    for read in &reads {
        let locus = sv.allele(read.v as u32).loci;
        let sequence = encode(&read.sequence);
        let got = calls(&custom, &sequence);
        if locus == Loci::of(Locus::Tra) || locus == Loci::of(Locus::Trg) {
            no_d_reads += 1;
            d_on_no_d += usize::from(got.0);
        } else if locus == Loci::of(Locus::Trd) {
            trd += 1;
            let want = calls(&named, &sequence);
            trd_d_same += usize::from(got.0 == want.0);
            if (&got.1, &got.2) != (&want.1, &want.2) {
                trd_jc_differ.push((sv.allele(read.v as u32).name.clone(), got, want));
            }
        }
    }
    assert!(no_d_reads > 250 && trd > 60, "{no_d_reads} {trd}");
    assert_eq!(d_on_no_d, 0, "TRA and TRG reads with a D");
    assert!(trd_jc_differ.is_empty(), "TRD reads: {trd_jc_differ:?}");
    // D presence differs only on a few chance D of a few bases.
    assert!(trd_d_same * 100 >= trd * 90, "{trd_d_same} of {trd}");
}
