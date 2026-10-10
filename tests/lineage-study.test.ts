import assert from "node:assert/strict";
import test from "node:test";
import { createHash } from "node:crypto";

import { lineageStudyRange, readLineageAirrSlice, validateLineageStudy, type LineageStudyManifest } from "../src/lineage-study.ts";
import type { LineageSummary } from "../src/post-analysis-core.ts";

function summary(id:number):LineageSummary{return {id,representativeOrdinal:0,uniqueMembers:1,abundance:4,locus:"IGH",vCalls:["IGHV1-2"],jCalls:["IGHJ4"],cdr3Length:9,studyScope:"subject",studyGroup:"donor_1",sampleIds:["sample_A"],sampleCounts:[{sampleId:"sample_A",uniqueMembers:1,abundance:4}],subjectIds:["donor_1"],timepoints:[],compartments:[]};}

function manifest(start:number,end:number,row="r1\tACGT\tsample_A\n"):LineageStudyManifest{return {
  schema:1,application:"Swig lineage study",applicationVersion:"0.29.0",createdAt:new Date(0).toISOString(),
  linkedAirr:{name:"study.lineages.airr.tsv",size:end,records:1,headers:["sequence_id","sequence","sample_id"],sha256:"test"},
  analysis:{inputName:"study",species:"Homo sapiens",scope:"IGH",references:{V:">V\nA\n",D:"",J:">J\nA\n",C:"",counts:{V:1,D:0,J:1,C:0},annotation:{V:{annotated:0,total:1},J:{annotated:0,total:1}},loci:["IGH"]},datasets:[{datasetId:"d1",inputName:"reads.fastq",sampleId:"sample_A",subjectId:"donor_1",cohort:"",timepoint:""}],callingProfile:"truth_optimized",minimumIdentity:.6,strand:0,lineage:{scope:"subject",identity:.85,resolution:"gene",ambiguity:"overlap",productiveOnly:true,maxCandidateComparisons:50_000}},
  summaries:[summary(1)],ranges:[{lineageId:1,start,end,records:1,sha256:createHash("sha256").update(row).digest("hex")}],
};}

test("lineage-study manifests validate ordered exact byte ranges",()=>{
  const header="sequence_id\tsequence\tsample_id\n";
  const row="r1\tACGT\tsample_A\n";
  const start=Buffer.byteLength(header),end=start+Buffer.byteLength(row);
  const value=manifest(start,end);
  assert.equal(lineageStudyRange(validateLineageStudy(value),1).start,start);
  assert.throws(()=>validateLineageStudy({...value,ranges:[{...value.ranges[0]},{...value.ranges[0],lineageId:2,start:start-1}]}),/invalid AIRR byte range/);
});

test("lazy lineage loading reads only the selected AIRR byte slice",async()=>{
  const header="sequence_id\tsequence\tsample_id\n";
  const row="r1\tACGT\tsample_A\n";
  const text=header+row;
  const start=Buffer.byteLength(header),end=Buffer.byteLength(text);
  const value=manifest(start,end);
  const file=new Blob([text]) as File;
  const slice=await readLineageAirrSlice(file,value,1);
  assert.equal(slice.header,header.trimEnd());
  assert.equal(slice.body,row);
  assert.equal(slice.records,1);
  await assert.rejects(readLineageAirrSlice(new Blob([header+"r2\tTGCA\tsample_B\n"]) as File,value,1),/does not match the indexed AIRR content/);
});
