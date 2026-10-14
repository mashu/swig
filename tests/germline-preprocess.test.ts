import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import zlib from "node:zlib";

import {
  applyIgblastAuxiliaryData,
  applyIgblastDFrameData,
  applyIgblastInternalData,
  prepareIgblastStyleGermlineFasta,
  preprocessGermlineFasta,
  headerLoci,
  namedLoci,
  preprocessGermlineFastaAcrossTiers,
  type MetadataAllele,
} from "../src/germline-preprocess.ts";
import { composeReferenceOverrides, referenceCellKey, segmentAppliesToLocus } from "../src/reference-composition.ts";
import { germlineTemplateTiers } from "../src/reference-pack.ts";
import {
  collectionsFor,
  collectionsForDatabase,
  databaseOptionsFor,
  databasesForCell,
  DEFAULT_DATABASE_ID,
  preferredDatabaseIdFor,
  REFERENCE_DATABASES,
} from "../src/reference-catalog.ts";

const pack = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("../public/references/imgt-202632-7-swig-0.7.json.gz", import.meta.url))).toString());
const human = pack.species.find((entry: { name: string }) => entry.name === "Homo sapiens");
const humanIghV = human.loci.IGH.V as MetadataAllele[];
const humanIghJ = human.loci.IGH.J as MetadataAllele[];

function firstAnnotated(alleles: MetadataAllele[], segment: "V" | "J"): MetadataAllele {
  const allele = alleles.find((entry) => segment === "V"
    ? entry[2]?.slice(2, 12).every((value) => value >= 0)
    : Boolean(entry[2] && entry[2]![0] >= 0 && entry[2]![1] >= 0));
  assert.ok(allele);
  return allele;
}

test("V delineation transfers from the nearest annotated relative without using functionality labels", () => {
  const template = firstAnnotated(humanIghV, "V");
  const changed = `${template[1].slice(0, 145)}${template[1][145] === "A" ? "C" : "A"}${template[1].slice(146)}`;
  const report = preprocessGermlineFasta(`>${template[0]}_S9999 pseudogene\n${changed}\n`, "V", humanIghV, ["IGH"]);
  assert.equal(report.annotated, 1);
  assert.equal(report.transferred, 1);
  assert.match(report.fasta, /SWIGMETA=[^\n]*,3\n/);
});

test("V coordinate transfer retains boundaries across a sequence-level frameshift", () => {
  const template = firstAnnotated(humanIghV, "V");
  const changed = `${template[1].slice(0, 130)}${template[1].slice(131)}`;
  const report = preprocessGermlineFasta(`>${template[0]}_S9998\n${changed}\n`, "V", humanIghV, ["IGH"]);
  assert.equal(report.annotated, 1);
  assert.equal(report.transferred, 1);
});

test("best-guess metadata matching removes identity rejection but preserves structural validation", () => {
  const template = firstAnnotated(humanIghV, "V");
  const anchorStart = template[2]![11] - 3;
  const changed = template[1].split("").map((base, index) => {
    if (index >= anchorStart && index < anchorStart + 3) return base;
    if (index % 3 !== 0) return base;
    return base === "A" ? "C" : "A";
  }).join("");
  const input = `>${template[0]}_S_best_guess\n${changed}\n`;
  const strict = preprocessGermlineFasta(input, "V", [template], ["IGH"], { includeDiagnostics: true });
  assert.equal(strict.annotated, 0);
  assert.ok((strict.diagnostics?.[0].rejectionCounts?.below_identity ?? 0) > 0);

  const guessed = preprocessGermlineFasta(input, "V", [template], ["IGH"], { mode: "best_guess", includeDiagnostics: true });
  assert.equal(guessed.annotated, 1);
  assert.equal(guessed.transferred, 1);
  assert.equal(guessed.diagnostics?.[0].status, "transferred");
  assert.equal(guessed.diagnostics?.[0].template, template[0]);
  assert.ok((guessed.diagnostics?.[0].identity ?? 1) < 0.8);
});

test("a V sequence that does not span the mapped IMGT intervals is not assigned fabricated boundaries", () => {
  const template = firstAnnotated(humanIghV, "V");
  const report = preprocessGermlineFasta(`>${template[0]}_S9997\n${template[1].slice(0, 90)}\n`, "V", humanIghV, ["IGH"]);
  assert.equal(report.annotated, 0);
  assert.equal(report.unannotated, 1);
  assert.doesNotMatch(report.fasta, /SWIGMETA=/);
});

test("J homology resolves multiple motif candidates but still verifies the mapped F/W-G anchor", () => {
  const template = firstAnnotated(humanIghJ, "J");
  const anchor = template[2]![1] + 1;
  assert.ok(anchor > 6);
  const changed = `TTCGGT${template[1].slice(6)}`;
  assert.notEqual(changed.slice(anchor, anchor + 6), "TTCGGT");
  const report = preprocessGermlineFasta(`>${template[0]}_S9999\n${changed}\n`, "J", humanIghJ, ["IGH"]);
  assert.equal(report.annotated, 1);
  assert.equal(report.transferred, 1);
  assert.match(report.fasta, /SWIGMETA=[^\n]*,6\n/);
});

test("IgBLAST-compatible preparation and annotation preserve each file's coordinate convention", () => {
  const v = firstAnnotated(humanIghV, "V");
  const j = firstAnnotated(humanIghJ, "J");
  const d = (human.loci.IGH.D as MetadataAllele[])[0];
  assert.ok(d);
  const dotted = `${v[1].slice(0, 20)}...${v[1].slice(20)}`;
  const normalized = prepareIgblastStyleGermlineFasta(`>IMGT|${v[0]}|Homo sapiens SWIGMETA=0,0,0,0,0,0,0,0,0,0,0,0,5\n${dotted}\n`, "V", ["IGH"]);
  assert.equal(normalized.fasta, `>${v[0]}\n${v[1]}\n`);
  assert.doesNotMatch(normalized.fasta, /SWIGMETA|\./);

  const bounds = v[2]!.slice(2, 12);
  const internal = `${v[0]} ${bounds.map((value, index) => index % 2 === 0 ? value + 1 : value).join(" ")} VH ${v[2]![0]}\n`;
  const vApplied = applyIgblastInternalData(normalized.fasta, internal);
  assert.equal(vApplied.matched, 1);
  assert.match(vApplied.fasta, new RegExp(`SWIGMETA=${v[2]![0]},-1,${bounds.join(",")},5`));

  const auxiliary = `# all coordinates are 0-based\n${j[0]} ${j[2]![0]} JH ${j[2]![1]} 1\n`;
  const jApplied = applyIgblastAuxiliaryData(`>${j[0]}\n${j[1]}\n`, auxiliary);
  assert.equal(jApplied.annotated, 1);
  assert.equal(jApplied.fwr4EndOffsets?.[j[0]], 1);
  assert.match(jApplied.fasta, new RegExp(`SWIGMETA=${j[2]![0]},${j[2]![1]},`));

  const dApplied = applyIgblastDFrameData(`>${d[0]}\n${d[1]}\n`, `${d[0]} 2\n`);
  assert.equal(dApplied.annotated, 1);
  assert.match(dApplied.fasta, /SWIGMETA=2,-1,/);
});

test("germline preprocessing rejects a locus mismatch", () => {
  assert.throws(
    () => preprocessGermlineFasta(">TRBV1*01\nACGTACGTACGT\n", "V", humanIghV, ["IGH"]),
    /outside the selected IGH search space/,
  );
});

test("the KI catalog exposes human IGH/TCR and macaque IGH collections by matching taxon and locus", () => {
  assert.deepEqual(collectionsFor("Homo sapiens", "IGH").map((entry) => entry.id), ["kiarva-human-igh"]);
  for (const locus of ["TRA", "TRB", "TRD", "TRG"] as const) {
    const collections = collectionsFor("Homo sapiens", locus);
    assert.equal(collections.length, 1);
    assert.equal(collections[0].locus, locus);
    assert.ok(collections[0].segments.V);
    assert.ok(collections[0].segments.J);
  }
  assert.equal(collectionsFor("Macaca mulatta_AG07107", "IGH")[0]?.id, "kimdb-macaca_mulatta");
  assert.equal(collectionsFor("Macaca fascicularis", "IGH")[0]?.id, "kimdb-macaca_fascicularis");
  assert.deepEqual(collectionsFor("Mus musculus_C57BL/6", "IGH"), []);
});

test("database choices prefer the compatible KI collection and retain IMGT as fallback", () => {
  const cat = databaseOptionsFor("Felis catus_Abyssinian", "202632-7");
  assert.deepEqual(cat.map((option) => option.id), [DEFAULT_DATABASE_ID]);
  assert.match(cat[0].label, /^IMGT\/GENE-DB 202632-7/);
  assert.deepEqual(databaseOptionsFor("Homo sapiens", "202632-7", "BCR").map((option) => option.id), ["kiarva-human-igh", DEFAULT_DATABASE_ID, "ki-human-tcr"]);
  assert.deepEqual(databaseOptionsFor("Homo sapiens", "202632-7", "TCR").map((option) => option.id), ["ki-human-tcr", DEFAULT_DATABASE_ID, "kiarva-human-igh"]);
  assert.deepEqual(databaseOptionsFor("Macaca mulatta_AG07107", "202632-7", "BCR").map((option) => option.id), ["kimdb-macaca_mulatta", DEFAULT_DATABASE_ID]);
  assert.equal(preferredDatabaseIdFor("Homo sapiens","BCR"),"kiarva-human-igh");
  assert.equal(preferredDatabaseIdFor("Homo sapiens","TCR"),"ki-human-tcr");
  assert.equal(preferredDatabaseIdFor("Macaca fascicularis_Cynomolgus","IGH"),"kimdb-macaca_fascicularis");
  assert.equal(preferredDatabaseIdFor("Felis catus_Abyssinian","BCR"),DEFAULT_DATABASE_ID);
});

test("published databases compose by receptor, locus, and segment", () => {
  const kiarva = REFERENCE_DATABASES.find((database) => database.id === "kiarva-human-igh");
  const tcr = REFERENCE_DATABASES.find((database) => database.id === "ki-human-tcr");
  assert.ok(kiarva);
  assert.ok(tcr);
  assert.deepEqual(collectionsForDatabase(kiarva, "BCR").map((collection) => collection.locus), ["IGH"]);
  assert.deepEqual(collectionsForDatabase(tcr, "TCR").map((collection) => collection.locus), ["TRA", "TRB", "TRD", "TRG"]);
  assert.deepEqual(collectionsForDatabase(tcr, "TRB").map((collection) => collection.locus), ["TRB"]);
  assert.deepEqual(databasesForCell("Homo sapiens", "IGH", "V").map((database) => database.id), ["kiarva-human-igh"]);
  assert.deepEqual(databasesForCell("Homo sapiens", "IGK", "V"), []);
  assert.deepEqual(databasesForCell("Homo sapiens", "TRB", "D").map((database) => database.id), ["ki-human-tcr"]);
});

test("a locus-specific replacement retains IMGT records from the other combined-chain loci", () => {
  const customIghV = ">IGHV_CUSTOM*01\nACGTACGTACGT\n";
  const composed = composeReferenceOverrides(
    ["IGH", "IGK", "IGL"],
    { [referenceCellKey("IGH", "V")]: { text: customIghV } },
    (locus, segment) => (human.loci[locus]?.[segment] ?? []).map((allele: MetadataAllele) => `>${allele[0]}\n${allele[1]}\n`).join(""),
  );
  assert.match(composed.V ?? "", /^>IGHV_CUSTOM\*01/m);
  assert.match(composed.V ?? "", new RegExp(`^>${human.loci.IGK.V[0][0].replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.match(composed.V ?? "", new RegExp(`^>${human.loci.IGL.V[0][0].replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.equal(composed.D, undefined);
  assert.equal(segmentAppliesToLocus("IGK", "D"), false);
  assert.equal(segmentAppliesToLocus("IGH", "D"), true);
});

test("per-cell allele exclusions alter the initial composed assignment FASTA only in that locus and segment", () => {
  const excludedIghV = human.loci.IGH.V[0][0];
  const retainedIghV = human.loci.IGH.V[1][0];
  const representativeIgkV = human.loci.IGK.V[0][0];
  const representativeIghJ = human.loci.IGH.J[0][0];
  const composed = composeReferenceOverrides(
    ["IGH", "IGK", "IGL"],
    {},
    (locus, segment) => (human.loci[locus]?.[segment] ?? []).map((allele: MetadataAllele) => `>${allele[0]}\n${allele[1]}\n`).join(""),
    { [referenceCellKey("IGH", "V")]: [excludedIghV] },
  );
  assert.ok(composed.V);
  assert.doesNotMatch(composed.V!, new RegExp(`^>${excludedIghV.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.match(composed.V!, new RegExp(`^>${retainedIghV.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.match(composed.V!, new RegExp(`^>${representativeIgkV.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.equal(composed.J, undefined);
  const baselineVCount = human.loci.IGH.V.length + human.loci.IGK.V.length + human.loci.IGL.V.length;
  assert.equal((composed.V!.match(/^>/gm) ?? []).length, baselineVCount - 1);
  const baselineJ = (human.loci.IGH.J ?? []).map((allele: MetadataAllele) => `>${allele[0]}\n${allele[1]}\n`).join("");
  assert.match(baselineJ, new RegExp(`^>${representativeIghJ.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
});

test("bundled KIMDB files are intact FASTA exports", () => {
  for (const species of ["Macaca_mulatta", "Macaca_fascicularis"]) {
    for (const segment of ["V", "D", "J"]) {
      const text = fs.readFileSync(new URL(`../public/references/kimdb-1.1/${species}/IGH/${segment}.fasta`, import.meta.url), "utf8");
      assert.ok(text.startsWith(">"));
      assert.ok((text.match(/^>/gm) ?? []).length >= (segment === "V" ? 600 : 14));
    }
  }
});

test("the built IMGT pack includes assembled heavy-chain and TCR constant references", () => {
  const rhesus = pack.species.find((entry: { name: string }) => entry.name === "Macaca mulatta_AG07107");
  const cynomolgus = pack.species.find((entry: { name: string }) => entry.name === "Macaca fascicularis");
  assert.ok(human.loci.IGH.C.length > 50);
  for (const prefix of ["IGHM", "IGHD", "IGHG", "IGHA", "IGHE"]) {
    assert.ok(human.loci.IGH.C.some(([name]: MetadataAllele) => name.startsWith(prefix)), prefix);
  }
  assert.ok(human.loci.TRA.C.length > 0);
  assert.ok(human.loci.TRB.C.length > 0);
  assert.ok(rhesus.loci.IGH.C.length > 0);
  assert.ok(cynomolgus.loci.IGH.C.length > 0);
  assert.ok((rhesus.loci.IGH.C as MetadataAllele[]).every((allele) => allele[1].length >= 300));
});

test("every bundled KIMDB macaque V and J record receives validated in-browser annotation", { timeout: 30_000 }, () => {
  for (const [speciesName, slug] of [["Macaca mulatta_AG07107", "Macaca_mulatta"], ["Macaca fascicularis", "Macaca_fascicularis"]]) {
    const selected = pack.species.find((entry: { name: string }) => entry.name === speciesName);
    assert.ok(selected);
    for (const segment of ["V", "J"] as const) {
      const text = fs.readFileSync(new URL(`../public/references/kimdb-1.1/${slug}/IGH/${segment}.fasta`, import.meta.url), "utf8");
      const report = preprocessGermlineFastaAcrossTiers(text, segment, germlineTemplateTiers(pack, selected, "IGH", segment), ["IGH"]);
      assert.equal(report.annotated, report.count, `${speciesName} ${segment}: ${report.warnings.join("; ")}`);
      assert.equal((report.fasta.match(/SWIGMETA=/g) ?? []).length, report.count);
    }
  }
});

test("the host reads a locus from an allele name exactly as the engine does", () => {
  for (const [label, loci] of [
    ["IGHV1-2*02", ["IGH"]],
    ["IGHT2D*01", ["IGH"]],
    ["IGHZ*01", ["IGH"]],
    ["IGHO*01", ["IGH"]],
    ["IGHG2A*01", ["IGH"]],
    ["IGHEP1*01", ["IGH"]],
    ["IGHD1-1*01", ["IGH"]],
    ["TRAV14/DV4*01", ["TRA", "TRD"]],
    ["Homo_sapiens_IGKV1-5*01", ["IGK"]],
    ["IGHTEST1", []],
    ["IgNAR*01", []],
    ["novel_1", []],
    ["myIGHV1", []],
    ["TRANSFER", []],
  ] as const) {
    assert.deepEqual(namedLoci(label), loci, label);
  }
});

test("a custom-named allele is written with the locus of its template", () => {
  const template = firstAnnotated(humanIghV, "V");
  const changed = `${template[1].slice(0, 145)}${template[1][145] === "A" ? "C" : "A"}${template[1].slice(146)}`;
  const report = preprocessGermlineFasta(
    `>novel_7 IGH sample allele\n${changed}\n>${template[0]}\n${template[1]}\n`,
    "V",
    humanIghV,
    ["IGH"],
  );
  assert.equal(report.transferred, 2);
  assert.match(report.fasta, /^>novel_7 SWIGMETA=[\d,-]+ locus=IGH$/m);
  // A name the engine places carries no key.
  assert.doesNotMatch(report.fasta, new RegExp(`^>${template[0].replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}.*locus=`, "m"));
  // A declared locus is kept, and survives the next template tier.
  const declared = preprocessGermlineFastaAcrossTiers(">sample_d1 locus=TRB,TRD\nGGGACAGGGGGC\n", "D", [], ["TRB", "TRD"]);
  assert.match(declared.fasta, /^>sample_d1 locus=TRB,TRD$/m);
});

test("the locus key of a custom allele follows its template and the header's words, not substrings", () => {
  // Regression: 'LIGHT' holds IGH and 'STRAIN' holds TRA, which the
  // substring test that admits a record (inferLocus) reads as loci.
  const oneChange = (sequence: string, at: number) => `${sequence.slice(0, at)}${sequence[at] === "A" ? "C" : "A"}${sequence.slice(at + 1)}`;
  const bcrV = (["IGH", "IGK", "IGL"] as const).flatMap((locus) => human.loci[locus].V as MetadataAllele[]);
  const igkv = bcrV.find((allele) => allele[0] === "IGKV1-5*01");
  assert.ok(igkv);
  const keyOf = (fasta: string, name: string) => new RegExp(`^>${name}\\b.*?(?: locus=([A-Z,]+))?$`, "m").exec(fasta)?.[1] ?? "";
  for (const [header, key] of [
    ["sample_k1 IGK light chain", "IGK"],
    ["sample_k1 light chain", "IGK"],
    ["sample_k1 IGK", "IGK"],
    // The header and the template disagree: no key, the engine decides.
    ["sample_k1 IGL light chain", ""],
  ] as const) {
    const report = preprocessGermlineFasta(`>${header}\n${oneChange(igkv[1], 145)}\n`, "V", bcrV, ["IGH", "IGK", "IGL"]);
    assert.equal(report.transferred, 1, header);
    assert.equal(keyOf(report.fasta, "sample_k1"), key, header);
  }
  const balbc = pack.species.find((entry: { name: string }) => entry.name === "Mus musculus_BALB/c");
  const tcrV = (["TRA", "TRB", "TRD", "TRG"] as const).flatMap((locus) => (balbc.loci[locus]?.V ?? []) as MetadataAllele[]);
  const trbv = tcrV.find((allele) => allele[0] === "TRBV12-2*02");
  assert.ok(trbv);
  const mouse = preprocessGermlineFasta(`>novel_b1 TRB strain BALB/c\n${oneChange(trbv[1], 145)}\n`, "V", tcrV, ["TRA", "TRB", "TRD", "TRG"]);
  assert.equal(keyOf(mouse.fasta, "novel_b1"), "TRB");
  // No template and no locus word: no key.
  const d = preprocessGermlineFasta(">sample_d2 light chain D\nGGTATAACTGGAACGAC\n", "D", [], ["IGH", "IGK", "IGL"]);
  assert.equal(keyOf(d.fasta, "sample_d2"), "");
  const stated = preprocessGermlineFasta(">sample_d3 TRB\nGGGACAGGGGGC\n", "D", [], ["TRB", "TRD"]);
  assert.equal(keyOf(stated.fasta, "sample_d3"), "TRB");
  for (const [header, loci] of [
    ["x IGK light chain", ["IGK"]],
    ["x light STRAIN transfer", []],
    ["x similar to IGKV1-5*01", ["IGK"]],
    ["x TRA TRD", ["TRA", "TRD"]],
    ["IGK annotation=transfer:IGHV1-2*01", []],
  ] as const) {
    assert.deepEqual(headerLoci(header), loci, header);
  }
});

test("a custom FASTA in a reference cell is given the cell's locus", () => {
  const custom = ">sample_k1\nACGTACGTACGT\n>IGKV1-5*99\nACGTACGTACGA\n>other locus=IGL\nACGTACGTACGC\n";
  const composed = composeReferenceOverrides(
    ["IGH", "IGK", "IGL"],
    { [referenceCellKey("IGK", "V")]: { text: custom } },
    (locus, segment) => (human.loci[locus]?.[segment] ?? []).map((allele: MetadataAllele) => `>${allele[0]}\n${allele[1]}\n`).join(""),
  );
  assert.match(composed.V ?? "", /^>sample_k1 locus=IGK$/m);
  assert.match(composed.V ?? "", /^>IGKV1-5\*99$/m);
  assert.match(composed.V ?? "", /^>other locus=IGL$/m);
  // Pack records are untouched.
  assert.doesNotMatch(composed.V ?? "", /^>IGHV.* locus=/m);
});

test("the IgBLAST-style preparation keeps the locus a custom record's header gives", () => {
  // Regression: the default swig-cli --vdj route wrote the bare name, so a
  // custom TRB D fell back to TRB,TRD and TRD reads were given it.
  const prepared = prepareIgblastStyleGermlineFasta(
    ">myBD1 locus=TRB\nGGGACAGGGGGC\n>myBD2 TRB D gene\nGGGACTAGCGGGGGGG\n>TRBD1*01 locus=TRB\nGGGACAGGGGGC\n>myDD1 TRD\nGAAATAGT\n",
    "D",
  );
  assert.equal(
    prepared.fasta,
    ">myBD1 locus=TRB\nGGGACAGGGGGC\n>myBD2 locus=TRB\nGGGACTAGCGGGGGGG\n>TRBD1*01 locus=TRB\nGGGACAGGGGGC\n>myDD1 locus=TRD\nGAAATAGT\n",
  );
  // The IgBLAST annotation files rewrite the headers and keep the key.
  const framed = applyIgblastDFrameData(prepared.fasta, "myBD1 0\nmyBD2 -1\n");
  assert.match(framed.fasta, /^>myBD1 SWIGMETA=0,-1,[-\d,]+ locus=TRB$/m);
  assert.match(framed.fasta, /^>myBD2 locus=TRB$/m);
  assert.match(framed.fasta, /^>myDD1 locus=TRD$/m);
  const j = firstAnnotated(humanIghJ, "J");
  const jApplied = applyIgblastAuxiliaryData(`>myJ1 IGH\n${j[1]}\n>other_j IGH\n${j[1]}\n`, `myJ1 ${j[2]![0]} JH ${j[2]![1]}\n`);
  assert.match(jApplied.fasta, /^>myJ1 SWIGMETA=[-\d,]+ locus=IGH$/m);
  assert.match(jApplied.fasta, /^>other_j locus=IGH$/m);
  const v = firstAnnotated(humanIghV, "V");
  const bounds = v[2]!.slice(2, 12);
  const vApplied = applyIgblastInternalData(`>myV1 IGH\n${v[1]}\n`, `myV1 ${bounds.map((value, index) => index % 2 === 0 ? value + 1 : value).join(" ")} VH ${v[2]![0]}\n`);
  assert.match(vApplied.fasta, /^>myV1 SWIGMETA=[-\d,]+ locus=IGH$/m);
});

test("later template tiers read a record's original header", () => {
  // Regression: tier 1 re-read tier 0's output, whose header had lost the
  // words that admitted the record, and refused it.
  const bcr = (segment: "V") => germlineTemplateTiers(pack, human, "BCR", segment);
  const igkv18 = (human.loci.IGK.V as MetadataAllele[]).find((allele) => allele[0] === "IGKV1-8*03");
  assert.ok(igkv18);
  const unannotated = `>IGKV99-1*01\n${"ACGT".repeat(75)}\n`;
  // The record's key, "" for none; undefined when the record is missing.
  const keyOf = (fasta: string, name: string) => {
    const found = new RegExp(`^>${name}\\b.*?(?: locus=([A-Z,]+))?$`, "m").exec(fasta);
    return found ? found[1] ?? "" : undefined;
  };
  for (const [header, key] of [
    // The header and the template disagree: no key, the engine decides.
    ["sample_k1 IGL", ""],
    ["sample_k1 IGK", "IGK"],
    ["sample_k1 light chain", "IGK"],
  ] as const) {
    const alone = preprocessGermlineFastaAcrossTiers(`>${header}\n${igkv18[1]}\n`, "V", bcr("V"), ["IGH", "IGK", "IGL"]);
    const report = preprocessGermlineFastaAcrossTiers(`>${header}\n${igkv18[1]}\n${unannotated}`, "V", bcr("V"), ["IGH", "IGK", "IGL"], { includeDiagnostics: true });
    assert.equal(report.count, 2, header);
    assert.equal(report.annotated, 1, header);
    assert.equal(keyOf(report.fasta, "sample_k1"), key, header);
    assert.equal(keyOf(alone.fasta, "sample_k1"), key, header);
    // The first tier's template stays the record's.
    assert.equal(report.diagnostics?.find((entry) => entry.name === "sample_k1")?.template, "IGKV1-8*03", header);
  }
});

test("a custom name with no locus code is placed by its scope or its template", () => {
  const igkv15 = (human.loci.IGK.V as MetadataAllele[]).find((allele) => allele[0] === "IGKV1-5*01");
  assert.ok(igkv15);
  const changed = `${igkv15[1].slice(0, 145)}${igkv15[1][145] === "A" ? "C" : "A"}${igkv15[1].slice(146)}`;
  // The record's key, "" for none; undefined when the record is missing.
  const keyOf = (fasta: string, name: string) => {
    const found = new RegExp(`^>${name}\\b.*?(?: locus=([A-Z,]+))?$`, "m").exec(fasta);
    return found ? found[1] ?? "" : undefined;
  };
  const igkTiers = germlineTemplateTiers(pack, human, "IGK", "V");
  // A reference cell's scope is its one locus: a bare name, or words holding
  // a code inside them ('light' holds IGH), take that locus.
  for (const header of ["novel_1", "novel_1 kappa light chain"]) {
    const cell = preprocessGermlineFastaAcrossTiers(`>${header}\n${changed}\n`, "V", igkTiers, ["IGK"]);
    assert.equal(cell.annotated, 1, header);
    assert.equal(keyOf(cell.fasta, "novel_1"), "IGK", header);
    assert.deepEqual(cell.loci, ["IGK"], header);
  }
  const bareJ = preprocessGermlineFasta(">novel_j\nGTGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC\n", "J", [], ["IGK"]);
  assert.equal(keyOf(bareJ.fasta, "novel_j"), "IGK");
  // A gene word of another locus is still refused.
  assert.throws(
    () => preprocessGermlineFasta(`>novel_1 IGHV-like\n${changed}\n`, "V", igkTiers[0], ["IGK"]),
    /novel_1 belongs to IGH, outside the selected IGK search space/,
  );
  // In a multi-locus scope the template places it.
  const bcrTiers = germlineTemplateTiers(pack, human, "BCR", "V");
  const bcr = preprocessGermlineFastaAcrossTiers(`>novel_1\n${changed}\n`, "V", bcrTiers, ["IGH", "IGK", "IGL"], { includeDiagnostics: true });
  assert.equal(keyOf(bcr.fasta, "novel_1"), "IGK");
  assert.deepEqual(bcr.loci, ["IGK"]);
  assert.equal(bcr.diagnostics?.[0].locus, "IGK");
  // Nothing places it: refused, saying how to place it.
  assert.throws(
    () => preprocessGermlineFastaAcrossTiers(">novel_d\nGGTATAACTGGAACGAC\n", "D", [], ["IGH", "IGK", "IGL"]),
    /novel_d names no IG\/TR locus and no IGH\/IGK\/IGL relative places it; name its locus in the header \(for example locus=IGH\)\./,
  );
  assert.throws(
    () => preprocessGermlineFasta(`>novel_v\n${"ACGT".repeat(75)}\n`, "V", bcrTiers[0], ["IGH", "IGK", "IGL"]),
    /novel_v names no IG\/TR locus/,
  );
});

test("constant-region FASTA may hold IMGT's IgD constant gene and J-C orphons", () => {
  // Regression: IGHD*01 read as a D gene and IGLJ-C/OR18*01 as a J gene.
  const c = ">IGHD*01\nGCACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGTGCAGACACCCAAAGGATAACAGCCCTGTGGTCCTGGCATGCTTGATAACTGGGTACCACCCAACGTCCGTGACTGTCACCTGGTACATGGGGACACAGAGCCAGCCCCAGAGAACCTTCCCTGAGATACAAAGACGGGACAGCTACTACATGACAAGCAGCCAGCTCTCCACCCCCCTCCAGCAGTGGCGCCAAGGCGAGTACAAATGCGTGGTCCAGCACACCGCCAGCAAGAGTAAGAAGGAGATCTTCCGCTGGCCAG\n>IGLJ-C/OR18*01\nGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTAGGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAAGTGACTTCTACCCGGGAGCCGTGACAGTGGCCTGGAAGGCAGATAGCAGCCCCGTCAAGGCGGGAGTGGAGACCACCACACCCTCCAAACAAAGCAACAACAAGTACGCGGCCAGCAGCTACCTGAGCCTGACGCCTGAGCAGTGGAAGTCCCACAGAAGCTACAGCTGCCAGGTCACGCATGAAGGGAGCACCGTGGAGAAGACAGTGGCCCCTACAGAATGTTCA\n";
  assert.equal(preprocessGermlineFasta(c, "C", [], ["IGH", "IGK", "IGL"]).count, 2);
  assert.equal(prepareIgblastStyleGermlineFasta(c, "C").count, 2);
  for (const [header, segment] of [["IGHD1-1*01", "D"], ["IGHJ4*02", "J"], ["IGLJ1*01", "J"]] as const) {
    assert.throws(() => preprocessGermlineFasta(`>${header}\nACGTACGTACGT\n`, "C", [], ["IGH", "IGL"]), new RegExp(`looks like a ${segment} gene, not a C gene`), header);
    assert.throws(() => prepareIgblastStyleGermlineFasta(`>${header}\nACGTACGTACGT\n`, "C"), new RegExp(`looks like a ${segment} gene, not a C gene`), header);
  }
});
