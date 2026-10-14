import assert from "node:assert/strict";
import test from "node:test";

import "fake-indexeddb/auto";
import { alignmentText, tableHeader, tableRow, treeNexus } from "../src/export-formats.ts";
import { augmentReferenceFasta, candidateFasta, MissingAlleleAccumulator, requiresCompleteGermlineWarning, type MissingAlleleCandidate } from "../src/germline-evidence.ts";
import { PersonalizedGermlineAccumulator, personalizedGermlineFasta, shmContextMutability } from "../src/personalized-germline.ts";
import { filterReferenceFasta, parseReferenceFasta } from "../src/reference-fasta.ts";
import { DEFAULT_REPERTOIRE_SELECTION, selectRepertoire } from "../src/repertoire-selection.ts";
import { AirrResultStore } from "../src/result-store.ts";
import { decodeSession, encodeSession, packSessionVector, SWIG_SESSION_SCHEMA, unpackSessionVector, type SwigSession } from "../src/session-state.ts";
import { computeShmMetric, ShmAccumulator } from "../src/shm-analysis.ts";

test("tabular, alignment, and tree serializers cover standard interchange formats", () => {
  const fields=["sequence_id","note"];
  assert.equal(tableHeader(fields,"csv"),"sequence_id,note\n");
  assert.equal(tableRow(fields,{sequence_id:"read,1",note:'a "quote"'},"csv"),'"read,1","a ""quote"""\n');
  assert.deepEqual(JSON.parse(tableRow(fields,{sequence_id:"r1",note:"ok"},"jsonl")),{sequence_id:"r1",note:"ok"});
  const records=[{name:"a one",sequence:"AC-G"},{name:"b",sequence:"ACTG"}];
  assert.match(alignmentText(records,"fasta"),/>a one/);
  assert.match(alignmentText(records,"clustal"),/^CLUSTAL W/);
  assert.match(alignmentText(records,"phylip"),/^2 4/);
  assert.match(alignmentText(records,"stockholm"),/^# STOCKHOLM 1.0/);
  assert.match(alignmentText(records,"nexus"),/dimensions ntax=2 nchar=4/);
  assert.match(treeNexus("(a:0.1,b:0.2);"),/\[&R\] \(a:0.1,b:0.2\);/);
});

test("SHM metrics distinguish synonymous and replacement changes and retain multiplicity", () => {
  const row={sequence_id:"r1",sample_id:"day_0",subject_id:"donor_A",swig_timepoint:"day_0",swig_compartment:"blood",v_call:"IGHV1-2*01",j_call:"IGHJ4*02",cdr3:"TGTGCCAAA",cdr3_aa:"CAK",duplicate_count:"3",v_sequence_alignment:"GCCGATGCT",v_germline_alignment:"GCTGCTGCT",v_sequence_start:"1",sequence_frame:"1",cdr1_start:"1",cdr1_end:"3",fwr1_start:"4",fwr1_end:"9"};
  const metric=computeShmMetric(row,0,7);
  assert.ok(metric);assert.equal(metric.vNtMutations,2);assert.equal(metric.synonymous,1);assert.equal(metric.vAaReplacements,1);assert.equal(metric.cdrNtMutations,1);assert.equal(metric.frameworkNtMutations,1);assert.equal(metric.duplicateCount,3);
  const accumulator=new ShmAccumulator({metric:"vNtRate",maxSamplesPerLineage:10});accumulator.add(row,0,7);accumulator.add({...row,sequence_id:"r2",sample_id:"day_30",swig_timepoint:"day_30",swig_compartment:"lymph_node",cdr3:"TGTGCCCCC",cdr3_aa:"CAP",duplicate_count:"2",v_sequence_alignment:"GCCGATGTT"},1,7);const dashboard=accumulator.finish();assert.equal(dashboard.analyzedAbundance,5);assert.equal(dashboard.lineages[0].label,"Lineage 7");assert.equal(dashboard.lineageSamples.length,2);assert.deepEqual(dashboard.lineageSamples.map((group)=>group.sampleId).sort(),["day_0","day_30"]);assert.ok(dashboard.lineages[0].maximum>=dashboard.lineages[0].p95);assert.equal(dashboard.lineageSamples.find((group)=>group.sampleId==="day_30")?.compartment,"lymph_node");assert.deepEqual(dashboard.lowestByLineage?.map((value)=>[value.lineageId,value.ordinal,value.cdr3Nt,value.cdr3Aa]),[[7,0,"TGTGCCAAA","CAK"]]);
});

test("SHM codons follow the read's frame when the V alignment starts inside a codon", () => {
  // Codons start at read positions 1, 4, 7, 10 (sequence_frame 1); the V
  // alignment starts at position 2, so its first two bases end the codon at 1.
  const row={sequence_id:"r1",v_call:"IGHV1-2*01",j_call:"IGHJ4*02",v_sequence_alignment:"AAGCCGATGCT",v_germline_alignment:"AAGCTGCTGCT",v_sequence_start:"2",sequence_frame:"1"};
  const metric=computeShmMetric(row);
  assert.ok(metric);assert.equal(metric.vNtMutations,2);assert.equal(metric.comparedCodons,3);assert.equal(metric.synonymous,1);assert.equal(metric.vAaReplacements,1);
});

test("SHM codons follow the V germline's frame on rows without sequence_frame", () => {
  // IgBLAST rows carry no sequence_frame; the alignment starts at V germline
  // base 2, so its first two bases end germline codon 1.
  const row={sequence_id:"r1",v_call:"IGHV1-2*01",j_call:"IGHJ4*02",v_sequence_alignment:"AAGCCGATGCT",v_germline_alignment:"AAGCTGCTGCT",v_sequence_start:"5",v_germline_start:"2"};
  const metric=computeShmMetric(row);
  assert.ok(metric);assert.equal(metric.comparedCodons,3);assert.equal(metric.synonymous,1);assert.equal(metric.vAaReplacements,1);
});

type GermlineTestRow={row:Record<string,string>;ordinal:number;lineageId:number};
function cdr3For(index:number,length=30+index%4){const encoded=index.toString(4).split("").map(value=>"ACGT"[Number(value)]).join("");return `TGT${encoded}${"A".repeat(length-encoded.length-6)}TGG`;}
function germlineRow(reference:string,index:number,mutations:number[]=[]):Record<string,string>{
  const start=11;const germline=[...reference.slice(start-1)];const query=[...germline];for(const fullPosition of mutations){const local=fullPosition-start;query[local]=query[local]==="A"?"G":"A";}
  const addPureGap=(value:string[])=>[...value.slice(0,20),".",...value.slice(20)].join("");
  return {subject_id:"donor_1",v_call:"IGHV1-2*01",j_call:`IGHJ${index%4+1}*01`,cdr3:cdr3For(index),v_germline_start:String(start),v_sequence_alignment:addPureGap(query),v_germline_alignment:addPureGap(germline)};
}
function runMissingAllele(rows:GermlineTestRow[],reference:string,options:ConstructorParameters<typeof MissingAlleleAccumulator>[0]={}){
  const discovery=new MissingAlleleAccumulator(options);for(const item of rows)discovery.add(item.row,item.ordinal,item.lineageId);const validator=discovery.prepareValidation(`>IGHV1-2*01\n${reference}\n`);for(const item of rows)validator.add(item.row,item.ordinal,item.lineageId);return validator.finish();
}

test("two-pass missing-V screen uses germline coordinates and emits independent-lineage evidence",()=>{
  const reference="ACGT".repeat(55);const rows:GermlineTestRow[]=[];for(let index=0;index<24;index+=1)rows.push({row:germlineRow(reference,index,index<9?[50]:[]),ordinal:index,lineageId:index+1});
  const result=runMissingAllele(rows,reference);assert.equal(result.mode,"lineage");assert.equal(result.validationPasses,2);assert.equal(result.independentUnits,24);assert.equal(result.candidates.length,1);assert.equal(result.candidates[0].independentUnits,9);assert.equal(result.candidates[0].substitutions[0].position,50);assert.equal(result.candidates[0].distinctCdr3s,9);assert.notEqual(result.candidates[0].sequence,reference);
});

test("reference exclusions are exact and preserve aligned FASTA metadata",()=>{
  const source=">IGHV1*01 SWIGMETA=1,2,3,4\nAC-GT\n>IGHV1*02 note\nACCGT\n";
  const filtered=filterReferenceFasta(source,["IGHV1*01","absent*01"]);
  assert.equal(filtered.total,2);assert.equal(filtered.retained,1);assert.equal(filtered.excluded,1);
  assert.deepEqual(filtered.excludedNames,["IGHV1*01"]);assert.deepEqual(filtered.unmatchedExclusions,["absent*01"]);
  assert.equal(parseReferenceFasta(filtered.fasta)[0].header,"IGHV1*02 note");
});

test("selected missing-V candidates export alone and append to an annotated V reference",()=>{
  const reference="ACGT".repeat(55);const rows:GermlineTestRow[]=[];for(let index=0;index<24;index+=1)rows.push({row:germlineRow(reference,index,index<9?[50]:[]),ordinal:index,lineageId:index+1});
  const dashboard=runMissingAllele(rows,reference);const candidate=dashboard.candidates[0];
  assert.ok(candidate);assert.equal(candidateFasta(dashboard,[]),"");assert.match(candidateFasta(dashboard,[candidate.id]),new RegExp(`>${candidate.id}`));
  const original=`>IGHV1-2*01 SWIGMETA=1,2,3,4\n${reference}\n>IGHV9-9*01\n${reference}\n`;
  const augmented=augmentReferenceFasta(original,dashboard,[candidate.id]);
  assert.equal(augmented.originalRecords,2);assert.equal(augmented.addedRecords,1);assert.equal(augmented.inheritedAnnotationRecords,1);
  const parsed=parseReferenceFasta(augmented.fasta);assert.equal(parsed.length,3);assert.equal(parsed[2].sequence,candidate.sequence);
  assert.match(parsed[2].header,/SWIGMETA=1,2,3,4/);assert.match(parsed[2].header,/SWIG_CANDIDATE=missing_allele_hint/);assert.match(parsed[2].header,/PARENT=IGHV1-2\*01/);
});

test("large or numerous possible-missing-V signals trigger the complete-reference warning",()=>{
  const candidate=(id:string,independentUnits:number)=>({id,independentUnits} as MissingAlleleCandidate);
  assert.equal(requiresCompleteGermlineWarning({candidates:[candidate("large",51)]}),true);
  assert.equal(requiresCompleteGermlineWarning({candidates:Array.from({length:6},(_,index)=>candidate(`small-${index}`,8))}),true);
  assert.equal(requiresCompleteGermlineWarning({candidates:[candidate("boundary",50),...Array.from({length:4},(_,index)=>candidate(`small-${index}`,8))]}),false);
});

test("a parent-reference nucleotide anywhere else in a supporting lineage vetoes that lineage",()=>{
  const reference="ACGT".repeat(55);const rows:GermlineTestRow[]=[];let ordinal=0;
  for(let index=0;index<24;index+=1){rows.push({row:germlineRow(reference,index,index<9?[50]:[]),ordinal:ordinal++,lineageId:index+1});if(index<9)rows.push({row:germlineRow(reference,index,[72,96,120]),ordinal:ordinal++,lineageId:index+1});}
  const result=runMissingAllele(rows,reference);assert.equal(result.candidates.length,0);assert.equal(result.referenceVetoedLineagePatterns,9);
});

test("recurrent third nucleotide states suppress a hotspot-like candidate",()=>{
  const reference="ACGT".repeat(55);const rows:GermlineTestRow[]=[];
  for(let index=0;index<24;index+=1){const row=germlineRow(reference,index,index<9?[50]:[]);if(index>=9&&index<11){const query=[...row.v_sequence_alignment];const column=[...row.v_germline_alignment].findIndex((_,columnIndex)=>{let position=10;for(let offset=0;offset<=columnIndex;offset+=1)if(row.v_germline_alignment[offset]!==".")position+=1;return position===50;});query[column]=row.v_sequence_alignment[column]==="T"?"C":"T";row.v_sequence_alignment=query.join("");}rows.push({row,ordinal:index,lineageId:index+1});}
  const result=runMissingAllele(rows,reference);assert.equal(result.candidates.length,0);assert.equal(result.conflictingLineagePatterns,2);
});

test("co-occurring candidate substitutions remain one linked haplotype rather than singleton warnings",()=>{
  const reference="ACGT".repeat(55);const rows:GermlineTestRow[]=[];for(let index=0;index<24;index+=1)rows.push({row:germlineRow(reference,index,index<10?[50,82]:[]),ordinal:index,lineageId:index+1});
  const result=runMissingAllele(rows,reference);assert.equal(result.candidates.length,1);assert.deepEqual(result.candidates[0].substitutions.map(item=>item.position),[50,82]);assert.equal(result.candidates[0].nearGermlineUnits,10);
});

test("clonal expansion cannot inflate missing-V support",()=>{
  const reference="ACGT".repeat(55);const rows:GermlineTestRow[]=[];let ordinal=0;
  for(let copy=0;copy<100;copy+=1)rows.push({row:germlineRow(reference,0,[50]),ordinal:ordinal++,lineageId:1});
  for(let index=1;index<5;index+=1)rows.push({row:germlineRow(reference,index,[50]),ordinal:ordinal++,lineageId:index+1});
  for(let index=5;index<25;index+=1)rows.push({row:germlineRow(reference,index,[]),ordinal:ordinal++,lineageId:index+1});
  const result=runMissingAllele(rows,reference);assert.equal(result.independentUnits,25);assert.equal(result.candidates.length,0);
});

function personalizedRow(parent:string,query:string,id:string,call="IGHV1-1*01"):Record<string,string>{
  return {sequence_id:id,subject_id:"donor_1",locus:"IGH",v_call:call,v_germline_start:"1",v_sequence_alignment:query,v_germline_alignment:parent};
}

test("personalized germline representative is the lowest current-assignment V-SHM row, not a random row",()=>{
  const reference="ACGT".repeat(50);const changed=[...reference];for(const position of [10,20,30,40,50])changed[position-1]=changed[position-1]==="A"?"C":"A";
  const accumulator=new PersonalizedGermlineAccumulator(`>IGHV1-1*01\n${reference}\n`,{minimumAlignedBases:100,minimumNovelSupport:20});
  accumulator.add(personalizedRow(reference,changed.join(""),"high-shm"),0,1);
  accumulator.add(personalizedRow(reference,reference,"low-shm"),1,1);
  assert.deepEqual(accumulator.selectedRepresentativeOrdinals(),[1]);
});

test("personalized germline fit recovers a two-allele expressed set with one vote per lineage",()=>{
  const first="ACGT".repeat(50);const second=[...first];second[49]=second[49]==="A"?"C":"A";const secondSequence=second.join("");
  const fasta=`>IGHV1-1*01 SWIGMETA=1,2,3\n${first}\n>IGHV1-1*02 SWIGMETA=1,2,3\n${secondSequence}\n>IGHV9-9*01\n${first}\n`;
  const accumulator=new PersonalizedGermlineAccumulator(fasta,{minimumAlignedBases:100,minimumNovelSupport:100});
  for(let lineage=1;lineage<=60;lineage+=1){const query=lineage<=30?first:secondSequence;accumulator.add(personalizedRow(first,query,`lineage-${lineage}`),lineage,lineage);for(let copy=0;copy<10;copy+=1)accumulator.add(personalizedRow(first,query,`lineage-${lineage}-copy-${copy}`),1000+lineage*10+copy,lineage);}
  const dashboard=accumulator.finish();const gene=dashboard.pools[0].genes.find((item)=>item.gene==="IGHV1-1");assert.ok(gene);assert.equal(dashboard.representativeLineages,60);assert.deepEqual(gene.activeAlleles.map((item)=>item.names[0]).sort(),["IGHV1-1*01","IGHV1-1*02"]);assert.ok(gene.activeAlleles.every((item)=>Math.abs(item.frequency-.5)<.02));
});

test("personalized germline fit can retain a recurrent same-length novel substitution haplotype",()=>{
  const reference="ACGT".repeat(50);const novel=[...reference];for(const position of [50,80])novel[position-1]=novel[position-1]==="A"?"C":"A";
  const fasta=`>IGHV1-1*01 SWIGMETA=1,2,3\n${reference}\n>IGHV9-9*01\n${reference}\n`;
  const accumulator=new PersonalizedGermlineAccumulator(fasta,{minimumAlignedBases:100,minimumNovelSupport:4});for(let lineage=1;lineage<=24;lineage+=1)accumulator.add(personalizedRow(reference,novel.join(""),`novel-${lineage}`),lineage,lineage);
  const dashboard=accumulator.finish();const allele=dashboard.pools[0].genes.find((item)=>item.gene==="IGHV1-1")?.activeAlleles[0];assert.ok(allele);assert.equal(allele.known,false);assert.deepEqual(allele.substitutions.map((item)=>item.position),[50,80]);
  const personalized=personalizedGermlineFasta(fasta,dashboard,dashboard.pools[0].id);assert.match(personalized,/SWIG_PERSONALIZED=novel/);assert.match(personalized,/>IGHV9-9\*01/);assert.doesNotMatch(personalized,/>IGHV1-1\*01(?:\s|\n)/);
});

test("personalized germline SHM prior distinguishes classic hot, neutral, and cold contexts",()=>{
  assert.ok(shmContextMutability("AACTA",2)>shmContextMutability("GGCAA",2));
  assert.ok(shmContextMutability("GGCAA",2)>shmContextMutability("CTCAA",2));
});

test("double-D positive repertoire selection uses sparse evidence and composes call/CDR3 filters", async () => {
  const store=new AirrResultStore();const header=["sequence_id","sequence","locus","v_call","d_call","j_call","productive","cdr3","cdr3_aa","v_identity","j_identity"].join("\t");
  const body=["r0\tACGTACGT\tIGH\tIGHV1*01\tIGHD1*01\tIGHJ4*01\tT\tTGTAAA\tCAK\t0.98\t0.97","r1\tACGTACGT\tIGH\tIGHV3*01\tIGHD2*01\tIGHJ6*01\tT\tTGTCCC\tCCP\t0.96\t0.95"].join("\n")+"\n";
  const ddHeader=["swig_batch_record_index","sequence_id","d_call","d2_call","swig_double_d_score_gain"].join("\t");const ddBody="1\tr1\tIGHD2*01\tIGHD3*01\t12\n";
  await store.appendBatch(header,body,{header:ddHeader,body:ddBody});await store.finalize();
  const selection=await selectRepertoire(store,{...DEFAULT_REPERTOIRE_SELECTION,doubleD:"positive",vCall:"IGHV3",cdr3Nt:"CCC"});assert.equal(selection.inputRecords,2);assert.equal(selection.retainedRecords,1);assert.equal(selection.mask[1],1);assert.equal(selection.mask[0],0);
  const page=await store.page({sequenceId:"",cdr3:"",locus:"",productive:"",vCall:"",dCall:"",jCall:"",cCall:"",isotype:"",minVIdentity:0,minDIdentity:0,minJIdentity:0,minCIdentity:0,minCdr3AaLength:0,maxCdr3AaLength:0,vjInFrame:"",stopCodon:"",completeVdj:"",revComp:"",hasD:false,hasDoubleD:true,hasCdr3:false},0,10);assert.equal(page.rows.length,1);assert.equal(page.rows[0].d2Call,"IGHD3*01");
  await store.clear();
});

test("repertoire call filters expose an explicit ambiguous-assignment policy", async () => {
  const store=new AirrResultStore();
  const header=["sequence_id","sequence","locus","v_call","j_call","productive","cdr3"].join("\t");
  const body=[
    "single\tACGT\tIGH\tIGHV1-2*01\tIGHJ4*01\tT\tTGTAAA",
    "ambiguous\tACGT\tIGH\tIGHV3-23*01,IGHV1-2*02\tIGHJ4*01\tT\tTGTCCC",
  ].join("\n")+"\n";
  await store.appendBatch(header,body);await store.finalize();
  const strict=await selectRepertoire(store,{...DEFAULT_REPERTOIRE_SELECTION,vCall:"IGHV1-2"});
  assert.equal(strict.retainedRecords,1);assert.equal(strict.mask[0],1);assert.equal(strict.mask[1],0);
  const inclusive=await selectRepertoire(store,{...DEFAULT_REPERTOIRE_SELECTION,vCall:"IGHV1-2",vCallIncludeAmbiguous:true});
  assert.equal(inclusive.retainedRecords,2);assert.equal(inclusive.mask[1],1);
  await store.clear();
});

test("AIRR session fingerprint is invariant to batch boundaries and changes with content", async()=>{
  const header="sequence_id\tsequence\tv_call\tj_call";const rows=["a\tACGT\tV1\tJ1\n","b\tTGCA\tV2\tJ2\n"];
  const first=new AirrResultStore();await first.appendBatch(header,rows.join(""));await first.finalize();
  const second=new AirrResultStore();await second.appendBatch(header,rows[0]);await second.appendBatch(header,rows[1]);await second.finalize();
  const changed=new AirrResultStore();await changed.appendBatch(header,`${rows[0]}b\tTGCC\tV2\tJ2\n`);await changed.finalize();
  assert.equal(first.fingerprint,second.fingerprint);assert.notEqual(first.fingerprint,changed.fingerprint);
  await first.clear();await second.clear();await changed.clear();
});

test("linked sessions gzip round-trip typed analysis vectors without embedding AIRR rows", async () => {
  const mask=new Uint8Array([1,0,1,1]);
  const session:SwigSession={schema:SWIG_SESSION_SCHEMA,application:"Swig",applicationVersion:"0.13.2",savedAt:new Date(0).toISOString(),linkedAirr:{name:"x.airr.tsv",size:100,lastModified:0,records:4,headers:["sequence_id","sequence"],fingerprint:"abcd"},analysis:{inputName:"x.fasta",species:"Homo sapiens",scope:"IGH",workers:2,callingProfile:"igblast_compatible",minimumIdentity:.6,strand:0,references:{V:">v\nACG\n",D:"",J:">j\nACG\n",C:"",counts:{V:1,D:0,J:1,C:0},annotation:{V:{annotated:1,total:1},J:{annotated:1,total:1}},loci:["IGH"]},referenceExclusions:{"IGH:V":["IGHV1-2*01"]}},doubleD:[],postAnalysis:{workingStages:[],activeMask:packSessionVector(mask)}};
  const restored=await decodeSession(await encodeSession(session));assert.equal(restored.linkedAirr.records,4);assert.equal(restored.analysis.callingProfile,"igblast_compatible");assert.deepEqual(restored.analysis.referenceExclusions,{"IGH:V":["IGHV1-2*01"]});assert.equal(JSON.stringify(restored).includes("ACGTACGT"),false);assert.deepEqual([...unpackSessionVector(restored.postAnalysis.activeMask!)],[1,0,1,1]);
  const controller=new AbortController();controller.abort();
  await assert.rejects(encodeSession(session,controller.signal),(error:unknown)=>error instanceof DOMException&&error.name==="AbortError");
});
