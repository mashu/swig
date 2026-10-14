import assert from "node:assert/strict";
import test from "node:test";

import "fake-indexeddb/auto";
import { AirrResultStore, EMPTY_FILTERS } from "../src/result-store.ts";

const encoder = new TextEncoder();
const header = [
  "sequence_id", "sequence", "locus", "v_call", "d_call", "j_call", "productive",
  "cdr3", "cdr3_aa", "junction_aa", "v_identity", "d_identity", "j_identity",
  "vj_in_frame", "stop_codon", "complete_vdj", "rev_comp",
].join("\t");

function makeBody(start: number, count: number): Uint8Array {
  let body = "";
  for (let offset = 0; offset < count; offset += 1) {
    const ordinal = start + offset;
    const locus = ordinal % 4 < 2 ? "IGH" : ordinal % 4 === 2 ? "IGK" : "IGL";
    body += [
      `read_${ordinal}`,
      "ACGT".repeat(75),
      locus,
      `IGHV${ordinal % 8 + 1}*01`,
      locus === "IGH" ? `IGHD${ordinal % 4 + 1}*01` : "",
      `${locus}J${ordinal % 5 + 1}*01`,
      ordinal % 3 ? "T" : "F",
      "TGTGCCAGAGATCGT",
      "CARD R".replace(" ", ""),
      "CARDRWG",
      "0.98",
      locus === "IGH" ? "1" : "",
      "0.97",
      "T",
      "F",
      "T",
      "F",
    ].join("\t") + "\n";
  }
  return encoder.encode(body);
}

test("50k AIRR records are committed, filtered, retrieved, and exported batchwise", async () => {
  const store = new AirrResultStore();
  const batchSize = 1000;
  for (let start = 0; start < 50_000; start += batchSize) {
    await store.appendBatch(header, makeBody(start, batchSize));
  }
  await store.finalize();

  assert.equal(store.count, 50_000);
  assert.equal(store.summary.assigned, 50_000);
  assert.equal(store.summary.withCdr3, 50_000);
  assert.equal(store.facets().loci.reduce((sum, value) => sum + value.count, 0), 50_000);
  const repertoire = store.repertoire({ locus: "IGH", ambiguity: "fractional" });
  assert.equal(repertoire.records, 25_000);
  assert.equal(repertoire.vCalls.reduce((sum, value) => sum + value.count, 0), 25_000);
  assert.ok(repertoire.vjPairs.length > 0);

  const page = await store.page({ ...EMPTY_FILTERS, locus: "IGK" }, 0, 25);
  assert.equal(page.rows.length, 25);
  assert.ok(page.rows.every((row) => row.locus === "IGK"));
  const detail = await store.detail(page.rows[0]);
  assert.equal(detail.sequence_id, page.rows[0].sequenceId);
  assert.equal(detail.sequence.length, 300);

  let exportedBytes = 0;
  await store.writeAirr(async (part) => {
    exportedBytes += typeof part === "string" ? encoder.encode(part).byteLength : part instanceof Blob ? part.size : part.byteLength;
  });
  assert.equal(exportedBytes, store.outputBytes);
  await store.clear();
});

test("gene and allele filters distinguish unambiguous from multi-hit calls", async () => {
  const store = new AirrResultStore();
  const rows = new TextDecoder().decode(makeBody(0, 3)).trimEnd().split("\n").map((line) => line.split("\t"));
  rows[0][3] = "IGHV1-2*01";
  rows[1][3] = "IGHV1-2*02,IGHV3-23*01";
  rows[2][3] = "IGHV3-23*01";
  await store.appendBatch(header, `${rows.map((row) => row.join("\t")).join("\n")}\n`);
  await store.finalize();

  const geneStrict = await store.page({ ...EMPTY_FILTERS, vCall: "IGHV1-2" }, 0, 10);
  assert.deepEqual(geneStrict.rows.map((row) => row.sequenceId), ["read_0"]);
  const geneAmbiguous = await store.page({ ...EMPTY_FILTERS, vCall: "ighv1-2", vCallIncludeAmbiguous: true }, 0, 10);
  assert.deepEqual(geneAmbiguous.rows.map((row) => row.sequenceId), ["read_0", "read_1"]);
  assert.equal((await store.page({ ...EMPTY_FILTERS, vCall: "IGHV1-2*02" }, 0, 10)).rows.length, 0);
  assert.equal((await store.page({ ...EMPTY_FILTERS, vCall: "IGHV1-2*02", vCallIncludeAmbiguous: true }, 0, 10)).rows.length, 1);
  assert.deepEqual((await store.page({ ...EMPTY_FILTERS, vCall: "IGHV3-23*01, IGHV1-2*02" }, 0, 10)).rows.map((row)=>row.sequenceId),["read_1"]);
  await store.clear();
});

test("constant-region evidence yields isotype facets and repertoire summaries", async () => {
  const store = new AirrResultStore();
  const constantHeader = [
    "sequence_id", "sequence", "locus", "v_call", "d_call", "j_call", "c_call",
    "productive", "cdr3_aa", "v_identity", "d_identity", "j_identity", "c_identity", "c_support",
    "c_sequence_alignment", "vj_in_frame", "stop_codon", "complete_vdj", "rev_comp",
  ].join("\t");
  const constantBody = [
    "ig_m\tACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ1*01\tIGHM*01\tT\tCARDR\t.98\t1\t.97\t.99\t\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_g_short\tACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ1*01\tIGHG1*01\tT\tCARDR\t.98\t1\t.97\t.99\t5e-5\tAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_g_short_anchored\tACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ1*01\tIGHG2*01\tT\tCARDR\t.98\t1\t.97\t.99\t2e-4\tAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_g_short_weak\tACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ1*01\tIGHG4*01\tT\tCARDR\t.98\t1\t.97\t.99\t2e-3\tAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_a\tACGT\tIGH\tIGHV2*01\tIGHD1*01\tIGHJ2*01\tIGHA2*01\tF\tCARDRR\t.95\t1\t.96\t.97\t\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_g_ambiguous\tACGT\tIGH\tIGHV2*01\tIGHD1*01\tIGHJ2*01\tIGHG1*01,IGHG3*01\tT\tCARDRR\t.95\t1\t.96\t.98\t\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
    "ig_e_weak\tACGT\tIGH\tIGHV2*01\tIGHD1*01\tIGHJ2*01\tIGHE*01\tT\tCARDRR\t.95\t1\t.96\t.89\t\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tT\tF\tT\tF",
  ].join("\n") + "\n";
  await store.appendBatch(constantHeader, encoder.encode(constantBody));
  await store.finalize();
  // A short C stretch needs the aligner's anchored gate (c_support ≤ 1e-3).
  assert.deepEqual(store.facets().isotypes.map((item) => item.value).sort(), ["IgA2", "IgG", "IgG1", "IgG2", "IgM"]);
  assert.equal(store.repertoire().isotypes.reduce((sum, item) => sum + item.count, 0), 5);
  const page = await store.page({ ...EMPTY_FILTERS, isotype: "IgM" }, 0, 10);
  assert.equal(page.rows.length, 1);
  assert.equal(page.rows[0].cCall, "IGHM*01");
  await store.clear();
});

test("study metadata is indexed for dataset, sample, donor, cohort, timepoint, and compartment filtering", async () => {
  const store = new AirrResultStore();
  const studyHeader = `${header}\tswig_dataset_id\tsample_id\tsubject_id\tswig_cohort\tswig_timepoint\tswig_compartment`;
  const baseRows = new TextDecoder().decode(makeBody(0, 3)).trimEnd().split("\n");
  const body = [
    `${baseRows[0]}\tdataset_1\tsample_A\tdonor_1\tcase\tday_0\tblood`,
    `${baseRows[1]}\tdataset_2\tsample_B\tdonor_1\tcase\tday_30\tlymph_node`,
    `${baseRows[2]}\tdataset_3\tsample_C\tdonor_2\tcontrol\tday_0\tblood`,
  ].join("\n") + "\n";
  await store.appendBatch(studyHeader, body);
  await store.finalize();
  assert.deepEqual(store.facets().samples.map((item) => item.value).sort(), ["sample_A", "sample_B", "sample_C"]);
  assert.equal(store.facets().subjects.find((item) => item.value === "donor_1")?.count, 2);
  assert.equal((await store.page({ ...EMPTY_FILTERS, subjectId: "donor_1" }, 0, 10)).rows.length, 2);
  assert.equal((await store.page({ ...EMPTY_FILTERS, cohort: "control", timepoint: "day_0" }, 0, 10)).rows[0].sampleId, "sample_C");
  assert.equal((await store.page({ ...EMPTY_FILTERS, compartment: "lymph_node" }, 0, 10)).rows[0].sampleId, "sample_B");
  const updated=[
    {datasetId:"dataset_1",inputName:"one.tsv",sampleId:"shared_sample",subjectId:"corrected_donor",cohort:"case",timepoint:"day_0",compartment:"blood"},
    {datasetId:"dataset_2",inputName:"two.tsv",sampleId:"shared_sample",subjectId:"corrected_donor",cohort:"case",timepoint:"day_30",compartment:"node"},
    {datasetId:"dataset_3",inputName:"three.tsv",sampleId:"sample_C",subjectId:"donor_2",cohort:"control",timepoint:"day_0",compartment:"blood"},
  ];
  await store.updateStudyMetadata(updated);
  assert.equal((await store.page({...EMPTY_FILTERS,subjectId:"corrected_donor"},0,10)).rows.length,2);
  assert.equal(store.facets().samples.find((item)=>item.value==="shared_sample")?.count,2);
  const [updatedDetail]=await store.detailMany([1]);
  assert.equal(updatedDetail.values.sample_id,"shared_sample");
  assert.equal(updatedDetail.values.subject_id,"corrected_donor");
  assert.equal(updatedDetail.values.swig_compartment,"node");
  const scanned:string[]=[];
  await store.scanAirrRows(["sample_id","subject_id"],(rows)=>{for(const row of rows)if(row.ordinal<2)scanned.push(`${row.values.sample_id}/${row.values.subject_id}`);});
  assert.deepEqual(scanned,["shared_sample/corrected_donor","shared_sample/corrected_donor"]);
  let exported="";await store.writeAirr(async(part)=>{exported+=typeof part==="string"?part:part instanceof Blob?await part.text():new TextDecoder().decode(part);});
  const exportedRows=exported.trimEnd().split("\n").map((line)=>line.split("\t"));const sampleColumn=exportedRows[0].indexOf("sample_id");const subjectColumn=exportedRows[0].indexOf("subject_id");
  assert.deepEqual(exportedRows.slice(1,3).map((row)=>[row[sampleColumn],row[subjectColumn]]),[["shared_sample","corrected_donor"],["shared_sample","corrected_donor"]]);
  await store.clear();
});

test("cancelling metadata re-indexing aborts the transaction and retains the previous index", async () => {
  const store = new AirrResultStore();
  const studyHeader = `${header}\tswig_dataset_id\tsample_id\tsubject_id\tswig_cohort\tswig_timepoint\tswig_compartment`;
  const sourceRows = new TextDecoder().decode(makeBody(0, 600)).trimEnd().split("\n");
  const body = `${sourceRows.map((row) => `${row}\tdataset_1\told_sample\told_donor\told_cohort\tday_0\tblood`).join("\n")}\n`;
  await store.appendBatch(studyHeader, body);
  await store.finalize();
  const controller = new AbortController();
  await assert.rejects(store.updateStudyMetadata([
    { datasetId: "dataset_1", inputName: "study.tsv", sampleId: "new_sample", subjectId: "new_donor", cohort: "new_cohort", timepoint: "day_1", compartment: "node" },
  ], (processed) => { if (processed >= 500) controller.abort(); }, controller.signal), (error: unknown) => error instanceof DOMException && error.name === "AbortError");
  assert.equal(store.facets().samples.find((item) => item.value === "old_sample")?.count, 600);
  assert.equal(store.facets().samples.some((item) => item.value === "new_sample"), false);
  const [detail] = await store.detailMany([599]);
  assert.equal(detail.values.sample_id, "old_sample");
  assert.equal(detail.values.subject_id, "old_donor");
  await store.clear();
});

test("saved study metadata is materialized during the first import without changing the AIRR fingerprint", async () => {
  const studyHeader = `${header}\tswig_dataset_id\tsample_id\tsubject_id\tswig_cohort\tswig_timepoint\tswig_compartment`;
  const baseRows = new TextDecoder().decode(makeBody(0, 3)).trimEnd().split("\n");
  const body = [
    `${baseRows[0]}\tdataset_1\told_A\told_donor\told\told\told`,
    `${baseRows[1]}\tdataset_2\told_B\told_donor\told\told\told`,
    `${baseRows[2]}\tdataset_3\told_C\told_donor\told\told\told`,
  ].join("\n") + "\n";
  const datasets = [
    { datasetId: "dataset_1", inputName: "one.tsv", sampleId: "sample_A", subjectId: "donor_1", cohort: "case", timepoint: "day_0", compartment: "blood" },
    { datasetId: "dataset_2", inputName: "two.tsv", sampleId: "sample_B", subjectId: "donor_1", cohort: "case", timepoint: "day_30", compartment: "node" },
    { datasetId: "dataset_3", inputName: "three.tsv", sampleId: "sample_C", subjectId: "donor_2", cohort: "control", timepoint: "day_0", compartment: "blood" },
  ];

  const baseline = new AirrResultStore();
  await baseline.appendBatch(studyHeader, body);
  await baseline.finalize();

  const restored = new AirrResultStore();
  restored.configureStudyMetadataForImport(datasets);
  await restored.appendBatch(studyHeader, body);
  await restored.finalize();

  assert.equal(restored.fingerprint, baseline.fingerprint, "saved metadata must not alter linked-file verification");
  assert.equal((await restored.page({ ...EMPTY_FILTERS, subjectId: "donor_1" }, 0, 10)).rows.length, 2);
  assert.deepEqual(restored.facets().samples.map((item) => item.value).sort(), ["sample_A", "sample_B", "sample_C"]);
  const [detail] = await restored.detailMany([1]);
  assert.equal(detail.values.sample_id, "sample_B");
  assert.equal(detail.values.subject_id, "donor_1");
  assert.equal(detail.values.swig_compartment, "node");
  let exported = "";
  await restored.writeAirr(async (part) => { exported += typeof part === "string" ? part : part instanceof Blob ? await part.text() : new TextDecoder().decode(part); });
  const exportedRows = exported.trimEnd().split("\n").map((line) => line.split("\t"));
  const sampleColumn = exportedRows[0].indexOf("sample_id");
  assert.deepEqual(exportedRows.slice(1).map((row) => row[sampleColumn]), ["sample_A", "sample_B", "sample_C"]);

  await baseline.clear();
  await restored.clear();
});

test("lineage export adds AIRR clone_id values and leaves excluded records blank", async () => {
  const store = new AirrResultStore();
  await store.appendBatch(header, makeBody(0, 3));
  await store.finalize();
  let output = "";
  await store.writeLineageAirr(Int32Array.from([7, 0, 12]), async (part) => {
    output += typeof part === "string" ? part : part instanceof Blob ? await part.text() : new TextDecoder().decode(part);
  });
  const lines = output.trimEnd().split("\n");
  assert.equal(lines[0].split("\t").at(-1), "clone_id");
  assert.equal(lines[1].split("\t").at(-1), "swig_lineage_7");
  assert.equal(lines[2].split("\t").at(-1), "");
  assert.equal(lines[3].split("\t").at(-1), "swig_lineage_12");
  await store.clear();
});

test("lineage export preserves original clone_id and adds an explicit derived merge column", async () => {
  const store = new AirrResultStore();
  await store.appendBatch(header, makeBody(0, 3));
  await store.finalize();
  let output = "";
  const merged = new Map<number, string>([[7, "swig_merged_lineage_1"], [12, "swig_merged_lineage_1"]]);
  await store.writeLineageAirrFormat(Int32Array.from([7, 0, 12]), "tsv", async (part) => {
    output += typeof part === "string" ? part : part instanceof Blob ? await part.text() : new TextDecoder().decode(part);
  }, merged);
  const rows = output.trimEnd().split("\n").map((line) => line.split("\t"));
  const clone = rows[0].indexOf("clone_id");
  const derived = rows[0].indexOf("swig_merged_lineage_id");
  assert.ok(clone >= 0 && derived >= 0);
  assert.deepEqual(rows.slice(1).map((row) => [row[clone], row[derived]]), [
    ["swig_lineage_7", "swig_merged_lineage_1"],
    ["", ""],
    ["swig_lineage_12", "swig_merged_lineage_1"],
  ]);
  await store.clear();
});

test("masked AIRR scans yield only active ordinals and report filtered progress totals", async () => {
  const store = new AirrResultStore();
  await store.appendBatch(header, makeBody(0, 8));
  await store.finalize();
  const mask = Uint8Array.from([0, 1, 0, 0, 1, 0, 0, 1]);
  const observed: number[] = [];
  const progress: Array<[number, number]> = [];
  await store.scanAirrRows(["sequence_id"], (rows) => {
    observed.push(...rows.map((row) => row.ordinal));
  }, {
    includeMask: mask,
    batchSize: 100,
    onProgress: (processed, total) => progress.push([processed, total]),
  });
  assert.deepEqual(observed, [1, 4, 7]);
  assert.deepEqual(progress.at(-1), [3, 3]);
  await store.clear();
});

test("double-D evidence is indexed and exported separately without rewriting the AIRR table", async () => {
  const store = new AirrResultStore();
  const doubleDHeader = [
    "swig_batch_record_index", "sequence_id", "standard_d_call", "d_call", "d2_call",
    "d_identity", "d2_identity", "d_sequence_start", "d_sequence_end",
    "d2_sequence_start", "d2_sequence_end", "np2", "np3", "swig_double_d_mode",
  ].join("\t");
  const doubleDBody = [
    "1", "read_1", "IGHD2*01", "IGHDA*01", "IGHDB*01", "0.99", "0.98",
    "110", "122", "129", "141", "AACCGG", "TT", "all",
  ].join("\t") + "\n";
  await store.appendBatch(header, makeBody(0, 3), { header: doubleDHeader, body: doubleDBody });
  const secondDoubleDBody = [
    "0", "read_3", "IGHD4*01", "IGHDC*01", "IGHDD*01", "0.97", "0.96",
    "112", "124", "130", "143", "AACCC", "T", "long_span",
  ].join("\t") + "\n";
  await store.appendBatch(header, makeBody(3, 2), { header: doubleDHeader, body: secondDoubleDBody });
  await store.finalize();

  assert.equal(store.doubleDCount, 2);
  const doubleDPage = await store.page({ ...EMPTY_FILTERS, hasDoubleD: true }, 0, 10);
  assert.equal(doubleDPage.rows.length, 2);
  assert.equal(doubleDPage.rows[0].ordinal, 1);
  assert.equal(doubleDPage.rows[1].ordinal, 3, "batch-local evidence indexes were not mapped to global AIRR ordinals");
  assert.equal(doubleDPage.rows[0].dCall, "IGHDA*01");
  assert.equal(doubleDPage.rows[0].d2Call, "IGHDB*01");
  assert.equal(doubleDPage.rows[0].dIdentity, 0.99);
  const detail = await store.detail(doubleDPage.rows[0]);
  assert.equal(detail.standard_d_call, "IGHD2*01");
  assert.equal(detail.d_call, "IGHDA*01");
  assert.equal(detail.d2_call, "IGHDB*01");
  assert.equal(detail.swig_batch_record_index, undefined);

  const main = await store.airrBlob();
  const mainText = await main.text();
  assert.ok(!mainText.split("\n", 1)[0].includes("d2_call"));
  const mainHeader = mainText.split("\n", 1)[0].split("\t");
  const mainRows = mainText.trimEnd().split("\n").slice(1).map((line) => line.split("\t"));
  assert.equal(mainRows[1][mainHeader.indexOf("d_call")], "IGHD2*01");

  let evidence = "";
  await store.writeDoubleD(async (part) => {
    evidence += typeof part === "string" ? part : part instanceof Blob ? await part.text() : new TextDecoder().decode(part);
  });
  const evidenceLines = evidence.trimEnd().split("\n");
  assert.ok(evidenceLines[0].startsWith("swig_airr_ordinal\tsequence_id"));
  assert.ok(!evidenceLines[0].includes("swig_batch_record_index"));
  assert.equal(evidenceLines[1].split("\t", 1)[0], "2");
  assert.ok(evidenceLines[1].includes("IGHDB*01"));
  assert.equal(evidenceLines[2].split("\t", 1)[0], "4");
  await store.clear();
});

test("direct output keeps byte offsets usable for on-demand detail", async () => {
  const parts: BlobPart[] = [];
  let file = new File([], "direct.airr.tsv");
  const store = new AirrResultStore({
    handle: { getFile: async () => file },
    writable: {
      write: async (part) => { parts.push(part instanceof Uint8Array ? part.slice().buffer : part); },
      close: async () => { file = new File(parts, "direct.airr.tsv", { type: "text/tab-separated-values" }); },
      abort: async () => { parts.length = 0; },
    },
  });
  await store.appendBatch(header, makeBody(0, 2));
  await store.finalize();
  const page = await store.page({ ...EMPTY_FILTERS }, 0, 2);
  const detail = await store.detail(page.rows[1]);
  assert.equal(detail.sequence_id, "read_1");
  assert.equal(file.size, store.outputBytes);
  assert.equal((await store.airrBlob()).size, file.size);
  await store.clear();
});
