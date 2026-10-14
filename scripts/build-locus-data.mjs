#!/usr/bin/env node

// Locus data the Swig engine derives from the shipped IMGT reference pack.
//
//   node scripts/build-locus-data.mjs labels <pack.json.gz> <out.tsv>
//     Every allele label of the pack with its locus and segment kind, one
//     line per locus and kind: `LOCUS<TAB>KIND<TAB>label label …`. The
//     engine's tests check that each label names exactly its pack locus
//     (engine/tests/loci.rs).
//
//   node scripts/build-locus-data.mjs panel <pack.json.gz> <out.fa>
//     The family-representative panel embedded in the engine
//     (engine/data/locus_panel.fa): V, J and C alleles of human and mouse
//     whose locus a custom-named allele is inferred from when the database
//     names none of its relatives (engine/src/locus_infer.rs).
//
//   node scripts/build-locus-data.mjs others <pack.json.gz> <out.fa>
//     TCR alleles of species outside the panel for the engine's tests
//     (engine/tests/data/tcr-other-species.fa): a subset of the rabbit's
//     TCR loci, and the V and J alleles of other species that the panel
//     once placed outside their locus (engine/tests/loci.rs).
//
// The pack is only read.

import fs from "node:fs";
import zlib from "node:zlib";

const [, , mode, packPath, outPath] = process.argv;
if (!["labels", "panel", "others"].includes(mode) || !packPath || !outPath) {
  console.error("Usage: node scripts/build-locus-data.mjs labels|panel|others <pack.json.gz> <out>");
  process.exit(2);
}

const pack = JSON.parse(zlib.gunzipSync(fs.readFileSync(packPath)).toString("utf8"));
const LOCI = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
const KINDS = ["V", "D", "J", "C"];

if (mode === "labels") {
  const groups = new Map();
  for (const species of pack.species) {
    for (const [locus, kinds] of Object.entries(species.loci)) {
      for (const [kind, records] of Object.entries(kinds)) {
        const key = `${locus}\t${kind}`;
        if (!groups.has(key)) groups.set(key, new Set());
        for (const record of records) groups.get(key).add(record[0]);
      }
    }
  }
  const lines = [];
  for (const locus of LOCI) {
    for (const kind of KINDS) {
      const labels = groups.get(`${locus}\t${kind}`);
      if (labels) lines.push(`${locus}\t${kind}\t${[...labels].sort().join(" ")}`);
    }
  }
  fs.writeFileSync(outPath, `${lines.join("\n")}\n`);
} else if (mode === "others") {
  const OTHERS = {
    "Oryctolagus cuniculus_Thorbecke inbred": "TRAV10*01 TRAV12*01 TRAV13-1*01 TRAV14-1*01 TRAV16*01 TRAV17*01 TRDV1*01 TRDV3*01 TRDV5*01 TRBV1*01 TRBV10*01 TRBV12*01 TRBV13*01 TRGV1-1*01 TRGV2*01 TRBD1*01 TRDD1*01 TRDD2*01 TRAJ10*01 TRAJ11*01 TRAJ12*01 TRAJ13*01 TRAJ14*01 TRAJ15*01 TRBJ1-1*01 TRBJ1-2*01 TRBJ1-3*01 TRBJ1-4*01 TRDJ1*01 TRDJ2*01 TRDJ3*01 TRGJ1*01 TRGJ2*01 TRAC*01 TRBC1*01 TRBC2*01 TRDC*01 TRGC*01",
    "Mustela putorius furo": "TRAV33*01 TRDV5*01",
    "Mustela putorius furo_sable": "TRBJ1-5*01",
    "Ovis aries_Rambouillet": "TRDV5*01",
    "Ovis aries_Romanov": "TRBJ1-2*01",
    "Camelus dromedarius": "TRBJ1-2*01",
    "Felis catus_Abyssinian": "TRBJ1-2*01 TRGJ2-2*01 TRGJ3-2*01",
    "Oncorhynchus mykiss": "TRBJ3*01",
  };
  const out = [];
  for (const [speciesName, labels] of Object.entries(OTHERS)) {
    const species = pack.species.find((s) => s.name === speciesName);
    if (!species) throw new Error(`no ${speciesName} in the pack`);
    for (const label of labels.split(" ")) {
      const found = [];
      for (const locus of LOCI) {
        for (const kind of KINDS) {
          for (const record of species.loci[locus]?.[kind] ?? []) {
            if (record[0] === label) found.push([kind, record[1]]);
          }
        }
      }
      if (found.length !== 1) throw new Error(`${speciesName} ${label}: ${found.length} records`);
      const [kind, sequence] = found[0];
      out.push(`>${label} kind=${kind} species=${speciesName.replace(/[ /]/g, "_")}`);
      out.push(sequence.toUpperCase().replace(/[^ACGTN]/g, ""));
    }
  }
  fs.writeFileSync(outPath, `${out.join("\n")}\n`);
} else {
  // The representatives measured in the design (scratchpad design4/locus,
  // "family-representative panel"): per human V family its first allele,
  // every human J and C gene's first allele, and the mouse families and
  // genes that no human representative stands for.
  const PANEL = {
    "Homo sapiens": {
      V: "IGHV1-18*01 IGHV2-10*01 IGHV3-11*01 IGHV4-28*01 IGHV5-10-1*01 IGHV6-1*01 IGHV7-34-1*01 IGHV8-51-1*01 IGKV1-12*01 IGKV1D-12*01 IGKV2-18*01 IGKV2D-18*01 IGKV3-11*01 IGKV3D-11*01 IGKV4-1*01 IGKV5-2*01 IGKV6-21*01 IGKV6D-21*01 IGKV7-3*01 IGLV1-36*01 IGLV10-54*01 IGLV11-55*01 IGLV2-11*01 IGLV3-1*01 IGLV4-10-1*01 IGLV5-37*01 IGLV6-57*01 IGLV7-43*01 IGLV8-61*01 IGLV9-49*01 TRAV1-1*01 TRAV10*01 TRAV11*01 TRAV12-1*01 TRAV13-1*01 TRAV14/DV4*01 TRAV16*01 TRAV17*01 TRAV18*01 TRAV19*01 TRAV2*01 TRAV20*01 TRAV21*01 TRAV22*01 TRAV23/DV6*01 TRAV24*01 TRAV25*01 TRAV26-1*01 TRAV27*01 TRAV29/DV5*01 TRAV3*01 TRAV30*01 TRAV34*01 TRAV35*01 TRAV36/DV7*01 TRAV38-1*01 TRAV39*01 TRAV4*01 TRAV40*01 TRAV41*01 TRAV5*01 TRAV6*01 TRAV7*01 TRAV8-1*01 TRAV9-1*01 TRBV1*01 TRBV10-1*01 TRBV11-1*01 TRBV12-1*01 TRBV13*01 TRBV14*01 TRBV15*01 TRBV16*01 TRBV17*01 TRBV18*01 TRBV19*01 TRBV2*01 TRBV20-1*01 TRBV21-1*01 TRBV23-1*01 TRBV24-1*01 TRBV25-1*01 TRBV26*01 TRBV27*01 TRBV28*01 TRBV29-1*01 TRBV3-1*01 TRBV30*01 TRBV4-1*01 TRBV5-1*01 TRBV6-1*01 TRBV7-1*01 TRBV9*01 TRDV1*01 TRDV2*01 TRDV3*01 TRGV1*01 TRGV10*01 TRGV11*01 TRGV2*01 TRGV3*01 TRGV4*01 TRGV5*01 TRGV5P*01 TRGV8*01 TRGV9*01 TRGVA*01",
      J: "IGHJ1*01 IGHJ2*01 IGHJ3*01 IGHJ4*01 IGHJ5*01 IGHJ6*01 IGKJ1*01 IGKJ2*01 IGKJ3*01 IGKJ4*01 IGKJ5*01 IGLJ-C/OR18*01 IGLJ1*01 IGLJ2*01 IGLJ2A*01 IGLJ2B*01 IGLJ3*01 IGLJ4*01 IGLJ5*01 IGLJ6*01 IGLJ7*01 TRAJ1*01 TRAJ10*01 TRAJ11*01 TRAJ12*01 TRAJ13*01 TRAJ14*01 TRAJ15*01 TRAJ16*01 TRAJ17*01 TRAJ18*01 TRAJ19*01 TRAJ2*01 TRAJ20*01 TRAJ21*01 TRAJ22*01 TRAJ23*01 TRAJ24*01 TRAJ25*01 TRAJ26*01 TRAJ27*01 TRAJ28*01 TRAJ29*01 TRAJ3*01 TRAJ30*01 TRAJ31*01 TRAJ32*01 TRAJ33*01 TRAJ34*01 TRAJ35*01 TRAJ36*01 TRAJ37*01 TRAJ38*01 TRAJ39*01 TRAJ4*01 TRAJ40*01 TRAJ41*01 TRAJ42*01 TRAJ43*01 TRAJ44*01 TRAJ45*01 TRAJ46*01 TRAJ47*01 TRAJ48*01 TRAJ49*01 TRAJ5*01 TRAJ50*01 TRAJ51*01 TRAJ52*01 TRAJ53*01 TRAJ54*01 TRAJ55*01 TRAJ56*01 TRAJ57*01 TRAJ58*01 TRAJ59*01 TRAJ6*01 TRAJ60*01 TRAJ61*01 TRAJ7*01 TRAJ8*01 TRAJ9*01 TRBJ1-1*01 TRBJ1-2*01 TRBJ1-3*01 TRBJ1-4*01 TRBJ1-5*01 TRBJ1-6*01 TRBJ2-1*01 TRBJ2-2*01 TRBJ2-2P*01 TRBJ2-3*01 TRBJ2-4*01 TRBJ2-5*01 TRBJ2-6*01 TRBJ2-7*01 TRDJ1*01 TRDJ2*01 TRDJ3*01 TRDJ4*01 TRGJ1*01 TRGJ2*01 TRGJP*01 TRGJP1*01 TRGJP2*01",
      C: "IGHA1*01 IGHA2*01 IGHD*01 IGHE*01 IGHG1*01 IGHG2*01 IGHG3*01 IGHG4*01 IGHG4A*01 IGHGP*01 IGHM*01 IGKC*01 IGLC1*01 IGLC2*01 IGLC2A*01 IGLC2B*01 IGLC3*01 IGLC6*01 IGLC7*01 IGLJ-C/OR18*01 TRAC*01 TRBC1*01 TRBC2*01 TRDC*01 TRGC1*01 TRGC2*01",
    },
    "Mus musculus": {
      V: "IGHV1-54*02 IGHV2-6-7*02 IGHV5-15*05 IGKV10-96*07 IGKV19-93*02 IGKV4-65*02 TRAV4D-3*04 TRAV6D-6*05 TRGV1*02 TRGV2*02 TRGV3*01 TRGV4*01 TRGV5*01 TRGV6*02 TRGV7*02",
      J: "IGHJ1*02 IGHJ2*01 IGHJ3*01 IGHJ4*01 IGKJ1*01 IGKJ2*01 IGKJ3*01 IGKJ4*01 IGKJ5*01 TRAJ31*02 TRAJ45*02 TRGJ3*01 TRGJ4*01",
      C: "IGHG2B*02 IGHM*03 IGKC*01",
    },
  };
  const out = [];
  for (const [speciesName, kinds] of Object.entries(PANEL)) {
    const species = pack.species.find((s) => s.name === speciesName);
    if (!species) throw new Error(`no ${speciesName} in the pack`);
    for (const kind of ["V", "J", "C"]) {
      for (const label of kinds[kind].split(" ")) {
        const found = [];
        for (const locus of LOCI) {
          for (const record of species.loci[locus]?.[kind] ?? []) {
            if (record[0] === label) found.push([locus, record[1]]);
          }
        }
        if (found.length !== 1) throw new Error(`${speciesName} ${label}: ${found.length} records`);
        const [locus, sequence] = found[0];
        const loci = label.startsWith("TRAV") && label.includes("/DV") ? "TRA,TRD" : locus;
        out.push(`>${label} locus=${loci} kind=${kind} species=${speciesName.replace(/ /g, "_")}`);
        out.push(sequence.toUpperCase().replace(/[^ACGTN]/g, ""));
      }
    }
  }
  fs.writeFileSync(outPath, `${out.join("\n")}\n`);
}
