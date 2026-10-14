import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { existsSync } from "node:fs";
import { chmod, mkdtemp, readFile, readdir, rm, writeFile } from "node:fs/promises";
import { availableParallelism } from "node:os";
import { join, resolve } from "node:path";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { gzipSync, gunzipSync } from "node:zlib";

function cliInvocation(root,configPath){
  const standalone=process.env.SWIG_CLI_EXECUTABLE;
  return {command:standalone||process.execPath,arguments_:standalone?["run","--config",configPath]:[join(root,"cli/swig-cli.mjs"),"run","--config",configPath]};
}

function runCli(root,configPath,extra=[],env=process.env){
  const {command,arguments_}=cliInvocation(root,configPath);
  return spawnSync(command,[...arguments_,...extra],{encoding:"utf8",timeout:120_000,env});
}

function runRawCli(root,args,env=process.env,input){
  const standalone=process.env.SWIG_CLI_EXECUTABLE;
  return spawnSync(standalone||process.execPath,standalone?args:[join(root,"cli/swig-cli.mjs"),...args],{encoding:"utf8",timeout:240_000,env,input,maxBuffer:256*1024*1024});
}

/**
 * The environment under which the CLI finds a native aligner: a standalone
 * executable carries its own; the repository bundle uses the engine's release
 * build when there is one (CI builds it first with `npm run test:engine`).
 * Null when neither exists, outside CI.
 */
function nativeAlignerEnvironment(root,cache){
  const env={...process.env,SWIG_CLI_CACHE_DIR:cache};
  delete env.SWIG_CLI_ALIGNER;delete env.SWIG_CLI_ALIGNER_PATH;
  if(process.env.SWIG_CLI_EXECUTABLE)return env;
  const binary=join(root,"engine/target/release",process.platform==="win32"?"swig.exe":"swig");
  if(existsSync(binary))return {...env,SWIG_CLI_ALIGNER_PATH:binary};
  if(process.env.CI)throw new Error(`CI tests the native aligner, but ${binary} does not exist: run npm run test:engine first.`);
  return null;
}

test("the CLI pipeline runs annotation through lazy lineage-study export",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-test-"));
  try{
    const output=join(temporary,"out");
    const configPath=join(temporary,"config.json");
    await writeFile(configPath,JSON.stringify({
      inputs:[{path:join(root,"tests/fixtures/cli-smoke.fasta"),sampleId:"sample-A",subjectId:"donor-A"}],
      annotation:{workers:1},
      pipeline:{lineage:{productiveOnly:false}},
      output:{directory:output,prefix:"smoke",airrCompression:"gzip"},
    }));
    const result=runCli(root,configPath);
    assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);
    assert.equal(summary.inputRecords,2);
    assert.equal(summary.annotatedRecords,2);
    assert.equal(summary.retainedRecords,1);
    assert.equal(summary.lineages,1);
    assert.deepEqual(summary.outputs,{annotated:"smoke.annotated.airr.tsv.gz",processed:"smoke.processed.airr.tsv.gz",airrCompression:"gzip"});
    for(const name of [summary.outputs.annotated,summary.outputs.processed]){
      const compressed=await readFile(join(output,name));
      assert.deepEqual([...compressed.subarray(0,2)],[0x1f,0x8b]);
      assert.match(gunzipSync(compressed).toString("utf8"),/^sequence_id\t/);
    }
    const resolved=JSON.parse(await readFile(join(output,"smoke.resolved-config.json"),"utf8"));
    assert.equal(resolved.annotation.assignerStrategy,undefined);
    assert.equal(resolved.annotation.callingProfile,"r_optimized");
    assert.equal(resolved.output.airrCompression,"gzip");
    const airr=await readFile(join(output,"smoke.lineages.airr.tsv"));
    const headers=airr.toString("utf8").split("\n",1)[0].split("\t");
    for(const field of ["v_support","d_support","j_support","c_support"])assert.ok(headers.includes(field));
    const manifest=JSON.parse(gunzipSync(await readFile(join(output,"smoke.swig-lineage-study.json.gz"))).toString("utf8"));
    assert.equal(manifest.linkedAirr.size,airr.byteLength);
    assert.equal(manifest.linkedAirr.sha256,createHash("sha256").update(airr).digest("hex"));
    assert.deepEqual(manifest.summaries[0].sampleCounts,[{sampleId:"sample-A",uniqueMembers:1,abundance:2}]);
    const range=manifest.ranges[0];
    const lineageBytes=airr.subarray(range.start,range.end);
    assert.equal(lineageBytes.toString("utf8").trimEnd().split("\n").length,1);
    assert.equal(range.sha256,createHash("sha256").update(lineageBytes).digest("hex"));
  }finally{
    await rm(temporary,{recursive:true,force:true});
  }
});

test("CLI FASTQ preprocessing matches browser QC, wrapped parsing, and per-dataset reservoir sampling",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-fastq-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const sequence=fasta.split(/\r?\n/).slice(1).join("").split(">read2",1)[0];
    const wrapped=(name,quality)=>`@${name}\n${sequence.slice(0,187)}\n${sequence.slice(187)}\n+\n${quality.slice(0,120)}\n${quality.slice(120)}\n`;
    const fastq=[wrapped("good-1","I".repeat(sequence.length)),wrapped("bad","!".repeat(sequence.length)),wrapped("good-2","I".repeat(sequence.length)),wrapped("good-3","I".repeat(sequence.length))].join("");
    const input=join(temporary,"reads.fastq");const output=join(temporary,"out");const configPath=join(temporary,"config.json");
    await writeFile(input,fastq);
    await writeFile(configPath,JSON.stringify({inputs:[{path:input,sampleId:"sample-QC",subjectId:"donor-QC"}],preprocessing:{fastqFilter:{enabled:true,maximumExpectedErrors:0.1,phredOffset:33,trim3Prime:{enabled:false}},subsample:{enabled:true,size:2,seed:7}},annotation:{workers:1},pipeline:{lineage:{productiveOnly:false}},output:{directory:output,prefix:"qc"}}));
    const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);
    assert.equal(summary.inputRecords,4);
    assert.equal(summary.eligibleRecords,3);
    assert.equal(summary.annotatedRecords,2);
    assert.equal(summary.preprocessing.fastqFilter.stats.recordsEvaluated,4);
    assert.equal(summary.preprocessing.fastqFilter.stats.recordsRejectedExpectedErrors,1);
    assert.deepEqual(summary.preprocessing.subsample,{enabled:true,size:2,seed:7});
    const resolved=JSON.parse(await readFile(join(output,"qc.resolved-config.json"),"utf8"));
    assert.equal(resolved.preprocessing.fastqFilter.maximumExpectedErrors,0.1);
    assert.equal(resolved.preprocessing.subsample.size,2);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI can address separate samples inside one concatenated gzip and preserves member names",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-members-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const records=fasta.trim().split(/(?=>)/).filter(Boolean);
    const first=gzipSync(`${records[0]}\n`),second=gzipSync(`${records[1]}\n`);
    const input=join(temporary,"joined.fasta.gz");const output=join(temporary,"out");const configPath=join(temporary,"config.json");
    await writeFile(input,Buffer.concat([first,second]));
    await writeFile(configPath,JSON.stringify({studyName:"members",inputs:[{path:input,inputName:"day-0.fasta.gz",format:"fasta",gzipRange:{start:0,end:first.length},datasetId:"d0",sampleId:"day-0",subjectId:"donor-A"},{path:input,inputName:"day-30.fasta.gz",format:"fasta",gzipRange:{start:first.length,end:first.length+second.length},datasetId:"d30",sampleId:"day-30",subjectId:"donor-A"}],annotation:{workers:1},pipeline:{lineage:{productiveOnly:false}},output:{directory:output,prefix:"members"}}));
    const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);assert.equal(summary.inputRecords,2);assert.equal(summary.annotatedRecords,2);assert.equal(summary.retainedRecords,2);assert.equal(summary.lineages,1);
    const manifest=JSON.parse(gunzipSync(await readFile(join(output,"members.swig-lineage-study.json.gz"))).toString("utf8"));
    assert.deepEqual(manifest.analysis.datasets.map((dataset)=>dataset.inputName),["day-0.fasta.gz","day-30.fasta.gz"]);
    const airr=await readFile(join(output,"members.annotated.airr.tsv"),"utf8");
    assert.match(airr,/d0::read1/);assert.match(airr,/d30::read2/);assert.match(airr,/swig_source_sequence_id/);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI AIRR mode can either preserve calls or reannotate sequence, including browser-accepted CSV input",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-airr-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const sequence=fasta.split(/\r?\n/)[1];
    const input=join(temporary,"reads.csv");await writeFile(input,`sequence_id,sequence,v_call,j_call,productive\noriginal,${sequence},FAKEV*01,FAKEJ*01,T\n`);
    for(const airrMode of ["preserve","reannotate"]){
      const output=join(temporary,airrMode);const configPath=join(temporary,`${airrMode}.json`);
      await writeFile(configPath,JSON.stringify({inputs:[{path:input,format:"airr",datasetId:"csv",sampleId:"sample",subjectId:"donor"}],annotation:{workers:1,airrMode},pipeline:{collapse:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:output,prefix:airrMode,writeLineageStudy:false}}));
      const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
      const table=await readFile(join(output,`${airrMode}.annotated.airr.tsv`),"utf8");
      assert.match(table,/csv::original/);assert.match(table,/swig_source_sequence_id/);
      if(airrMode==="preserve")assert.match(table,/FAKEV\*01/);else assert.doesNotMatch(table,/FAKEV\*01/);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI CHMMAIRRa uses the configured worker pool without changing output",async()=>{
  const root=resolve(import.meta.dirname,"..");const temporary=await mkdtemp(join(root,"tmp-cli-chmm-"));
  try{
    const input=join(temporary,"reads.airr.tsv");
    const header="sequence_id\tsequence\tlocus\tv_call\tj_call\tproductive\tcdr3\tv_sequence_alignment\tv_germline_alignment\n";
    const body=Array.from({length:501},(_,index)=>`r${index}\tAACA\tIGH\tIGHV1*01\tIGHJ4*01\tT\tTGTA\tAACA\tAAAA\n`).join("");
    await writeFile(input,header+body);
    const run=async(workers)=>{const output=join(temporary,`out-${workers}`),configPath=join(temporary,`config-${workers}.json`);await writeFile(configPath,JSON.stringify({inputs:[{path:input,format:"airr",sampleId:"sample",subjectId:"donor"}],annotation:{workers,airrMode:"preserve"},pipeline:{collapse:{enabled:false},chimera:{enabled:true,segment:"V",model:"BW",msaSource:"uploaded",uploadedMsa:">IGHV1*01\nAAAA\n>IGHV2*01\nCCCC\n",minimumDfr:1,posteriorThreshold:0.99,retainUnevaluated:true},lineage:{enabled:false},shm:{enabled:false}},output:{directory:output,prefix:"chmm",writeLineageStudy:false}}));const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);return{result,table:await readFile(join(output,"chmm.processed.airr.tsv"))};};
    const serial=await run(1),parallel=await run(2);assert.deepEqual(parallel.table,serial.table);assert.match(parallel.result.stderr,/CHMMAIRRa: 2 workers across 3 row batches/);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI FAD denoising runs independent V/J partitions in workers without changing output",async()=>{
  const root=resolve(import.meta.dirname,"..");const temporary=await mkdtemp(join(root,"tmp-cli-fad-"));
  try{
    const encode=(value)=>{let result="";for(let place=0;place<5;place+=1){result="ACGT"[value&3]+result;value>>>=2;}return result;};
    const input=join(temporary,"reads.airr.tsv");
    const header="sequence_id\tsequence\tlocus\tv_call\tj_call\tproductive\tcdr3\tsequence_alignment\n";
    const body=Array.from({length:600},(_,index)=>{const sequence=`ACGTACGT${encode(index)}`,v=index<300?"IGHV1-2*01":"IGHV3-7*01";return `r${index}\t${sequence}\tIGH\t${v}\tIGHJ4*02\tT\t${sequence.slice(-9)}\t${sequence}\n`;}).join("");
    await writeFile(input,header+body);
    const run=async(workers)=>{const output=join(temporary,`out-${workers}`),configPath=join(temporary,`config-${workers}.json`);await writeFile(configPath,JSON.stringify({inputs:[{path:input,format:"airr",sampleId:"sample",subjectId:"donor"}],annotation:{workers,airrMode:"preserve"},pipeline:{collapse:{enabled:true,mode:"fad"},chimera:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:output,prefix:"fad",writeLineageStudy:false}}));const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);return{result,table:await readFile(join(output,"fad.processed.airr.tsv"))};};
    const serial=await run(1),parallel=await run(2);assert.deepEqual(parallel.table,serial.table);assert.match(parallel.result.stderr,/FAD denoising: 2 workers across 2 independent V\/J partitions/);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI executes a browser-exported pasted input embedded in the JSON config",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-inline-"));
  try{
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const output=join(temporary,"out"),configPath=join(temporary,"config.json");
    await writeFile(configPath,JSON.stringify({inputs:[{path:"pasted-sequences.txt",inputName:"pasted-sequences.txt",inline:fasta,format:"fasta",sampleId:"paste",subjectId:"donor"}],annotation:{workers:1},pipeline:{lineage:{productiveOnly:false}},output:{directory:output,prefix:"inline"}}));
    const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
    const summary=JSON.parse(result.stdout);assert.equal(summary.inputRecords,2);assert.equal(summary.annotatedRecords,2);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI prepares local V/J metadata by default and honors the config opt-out",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-reference-prep-"));
  try{
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const human=pack.species.find((entry)=>entry.name==="Homo sapiens");
    const v=human.loci.IGH.V.find((entry)=>entry[2]?.slice(2,12).every((value)=>value>=0));
    const d=human.loci.IGH.D[0];
    const j=human.loci.IGH.J.find((entry)=>entry[2]?.[0]>=0&&entry[2]?.[1]>=0);
    assert.ok(v&&d&&j);
    // Saved as some editors save them, V with a UTF-8 byte-order mark and J with a Latin-1
    // description: both read as Swig Web reads them.
    await writeFile(join(temporary,"custom-v.fasta"),`\uFEFF>${v[0]}\n${v[1]}\n`);
    await writeFile(join(temporary,"custom-j.fasta"),Buffer.concat([Buffer.from(`>${j[0]} caf`),Buffer.from([0xe9]),Buffer.from(`\n${j[1]}\n`)]));
    const input=join(temporary,"read.fasta");
    await writeFile(input,`>custom-read\n${v[1]}${d[1]}${j[1]}\n`);
    for(const prepareMetadata of [true,false]){
      const label=prepareMetadata?"prepared":"raw";
      const output=join(temporary,label),configPath=join(temporary,`${label}.json`);
      const references={species:"Homo sapiens",scope:"IGH",files:{V:"custom-v.fasta",J:"custom-j.fasta"}};
      if(!prepareMetadata)references.prepareMetadata=false;
      await writeFile(configPath,JSON.stringify({inputs:[{path:input,format:"fasta",sampleId:"sample",subjectId:"donor"}],references,annotation:{workers:1},pipeline:{collapse:{enabled:false},lineage:{productiveOnly:false},shm:{enabled:false}},output:{directory:output,prefix:label}}));
      const result=runCli(root,configPath);assert.equal(result.status,0,result.stderr);
      const summary=JSON.parse(result.stdout);
      assert.equal(summary.references.metadataPreparation.enabled,prepareMetadata);
      const manifest=JSON.parse(gunzipSync(await readFile(join(output,`${label}.swig-lineage-study.json.gz`))).toString("utf8"));
      assert.match(manifest.analysis.references.V,/^>/);
      if(prepareMetadata){
        assert.equal(summary.references.metadataPreparation.segments.V.annotated,1);
        assert.equal(summary.references.metadataPreparation.segments.J.annotated,1);
        assert.match(manifest.analysis.references.V,/SWIGMETA=/);
        assert.match(manifest.analysis.references.J,/SWIGMETA=/);
        assert.match(result.stderr,/Preparing custom V reference metadata/);
      }else{
        assert.deepEqual(summary.references.metadataPreparation.segments,{});
        assert.doesNotMatch(manifest.analysis.references.V,/SWIGMETA=/);
        assert.doesNotMatch(manifest.analysis.references.J,/SWIGMETA=/);
        assert.doesNotMatch(result.stderr,/Preparing custom V reference metadata/);
      }
      const resolved=JSON.parse(await readFile(join(output,`${label}.resolved-config.json`),"utf8"));
      assert.equal(resolved.references.prepareMetadata,prepareMetadata);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("pipeline CLI requires an explicit output directory and command-line workers override config",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-output-"));
  try{
    const input=join(root,"tests/fixtures/cli-smoke.fasta");
    const missingPath=join(temporary,"missing-output.json");
    await writeFile(missingPath,JSON.stringify({inputs:[{path:input}],annotation:{workers:1}}));
    const missing=runCli(root,missingPath);
    assert.notEqual(missing.status,0);
    assert.match(missing.stderr,/requires an explicit output directory/);

    const invalidProfilePath=join(temporary,"invalid-profile.json");
    await writeFile(invalidProfilePath,JSON.stringify({inputs:[{path:input}],annotation:{workers:1,callingProfile:"igblast_balanced"},output:{directory:join(temporary,"invalid-profile")}}));
    const invalidProfile=runCli(root,invalidProfilePath);
    assert.notEqual(invalidProfile.status,0);
    assert.match(invalidProfile.stderr,/annotation\.callingProfile must be r_optimized, truth_optimized, igblast_compatible, sensitive_d/);

    const output=join(temporary,"out"),configPath=join(temporary,"config.json");
    await writeFile(configPath,JSON.stringify({inputs:[{path:input}],annotation:{workers:1},pipeline:{collapse:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:output,prefix:"workers"}}));
    const result=runCli(root,configPath,["--workers","2","--airr-compression","gzip"]);
    assert.equal(result.status,0,result.stderr);
    const resolved=JSON.parse(await readFile(join(output,"workers.resolved-config.json"),"utf8"));
    assert.equal(resolved.annotation.workers,2);
    assert.equal(resolved.output.airrCompression,"gzip");
    assert.deepEqual([...(await readFile(join(output,"workers.processed.airr.tsv.gz"))).subarray(0,2)],[0x1f,0x8b]);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("--vdj streams assignment-only, IgBLAST-data, and Swig-annotation modes without downstream state",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-vdj-"));
  try{
    const help=runRawCli(root,["--vdj","--help"]);
    assert.equal(help.status,0,help.stderr);
    assert.match(help.stdout,/r_optimized: 4-nt D exact-run floor; a D pays its prior cost/);
    assert.match(help.stdout,/sensitive_d: the same scores with a D twice as likely a priori/);
    assert.doesNotMatch(help.stdout,/--assigner/);
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const human=pack.species.find((entry)=>entry.name==="Homo sapiens");
    const v=human.loci.IGH.V.find((entry)=>entry[2]?.slice(2,12).every((value)=>value>=0));
    const d=human.loci.IGH.D[0];
    const j=human.loci.IGH.J.find((entry)=>entry[2]?.[0]>=0&&entry[2]?.[1]>=0);
    assert.ok(v&&d&&j);
    const vPath=join(temporary,"V.fasta"),dPath=join(temporary,"D.fasta"),jPath=join(temporary,"J.fasta"),queryPath=join(temporary,"query.fasta");
    await writeFile(vPath,`>IMGT|${v[0]}|Homo sapiens\n${v[1].slice(0,30)}...${v[1].slice(30)}\n`);
    await writeFile(dPath,`>${d[0]}\n${d[1]}\n`);
    await writeFile(jPath,`>${j[0]}\n${j[1]}\n`);
    await writeFile(queryPath,`>read-1\n${v[1]}AACCGG${d[1]}TTG${j[1]}\n`);
    const bounds=v[2].slice(2,12);
    const internalPath=join(temporary,"human.ndm.imgt"),auxPath=join(temporary,"human_gl.aux");
    await writeFile(internalPath,`# FWR/CDR positions are 1-based; frame is 0-based\n${v[0]} ${bounds.map((value,index)=>index%2===0?value+1:value).join(" ")} VH ${v[2][0]}\n`);
    await writeFile(auxPath,`# name, frame, chain, CDR3 stop, trailing non-coding bases; all but intervals are 0-based\n${j[0]} ${j[2][0]} JH ${j[2][1]} 1\n`);
    const common=["--vdj","-query",queryPath,"-germline_db_V",vPath,"-germline_db_D",dPath,"-germline_db_J",jPath,"-outfmt","19","-num_threads","3","--workers","1","--batch-records","1"];
    const readRow=async(path)=>{
      const [header,line]=String(await readFile(path,"utf8")).trimEnd().split("\n");
      const names=header.split("\t"),values=line.split("\t");return Object.fromEntries(names.map((name,index)=>[name,values[index]??""]));
    };

    const plainPath=join(temporary,"plain.airr.tsv");
    const plain=runRawCli(root,[...common,"-out",plainPath]);
    assert.equal(plain.status,0,plain.stderr);assert.match(plain.stderr,/assignments-only; 1 (worker|thread);/);assert.match(plain.stderr,/Swig aligner/);assert.match(plain.stderr,/R-optimized profile/);
    const plainRow=await readRow(plainPath);
    assert.equal(plainRow.v_call,v[0]);assert.equal(plainRow.j_call,j[0]);
    // Without annotation input the aligner annotates the germlines itself.
    assert.ok(plainRow.cdr1);assert.ok(plainRow.cdr3);
    for(const field of ["v_support","d_support","j_support"])assert.ok(plainRow[field]&&Number.isFinite(Number(plainRow[field]))&&Number(plainRow[field])>=0);
    const truthPath=join(temporary,"truth.airr.tsv");
    const truth=runRawCli(root,[...common,"--calling-profile","truth_optimized","-out",truthPath]);
    assert.equal(truth.status,0,truth.stderr);
    const truthRow=await readRow(truthPath);
    for(const field of ["v_support","d_support","j_support"])assert.ok(truthRow[field]&&Number.isFinite(Number(truthRow[field]))&&Number(truthRow[field])>=0);
    assert.equal(plainRow.c_support,"");

    const gzipPath=join(temporary,"plain.airr.tsv.gz");
    const gzipResult=runRawCli(root,[...common,"-out",gzipPath]);
    assert.equal(gzipResult.status,0,gzipResult.stderr);
    const compressed=await readFile(gzipPath);
    assert.deepEqual([...compressed.subarray(0,2)],[0x1f,0x8b]);
    const gzipText=gunzipSync(compressed).toString("utf8");
    const gzipLines=gzipText.split("\n").filter((line)=>line.length>0);
    const gzipHeaders=gzipLines[0].split("\t");
    assert.equal(gzipHeaders[12],"j_call");
    assert.ok(gzipLines.slice(1).every((line)=>line.split("\t").length===gzipHeaders.length));
    assert.equal(gzipLines[1].split("\t")[12],j[0]);

    const c=human.loci.IGH.C[0];
    const cPath=join(temporary,"C.fasta"),cQueryPath=join(temporary,"query-with-c.fasta");
    const replace={A:"C",C:"G",G:"T",T:"A"};
    const cPrefix=c[1].slice(0,72);
    const mismatchingPrefix=`${replace[cPrefix[0]]}${replace[cPrefix[1]]}${cPrefix.slice(2)}`;
    await writeFile(cPath,`>${c[0]}\n${c[1]}\n`);
    await writeFile(cQueryPath,`>read-c\n${v[1]}AACCGG${d[1]}TTG${j[1]}${mismatchingPrefix}\n`);
    const cCommon=[...common];
    cCommon[cCommon.indexOf("-query")+1]=cQueryPath;
    cCommon.splice(cCommon.indexOf("-outfmt"),0,"-c_region_db",cPath);
    const cPermissivePath=join(temporary,"c-default.airr.tsv");
    const cPermissive=runRawCli(root,[...cCommon,"-out",cPermissivePath]);
    assert.equal(cPermissive.status,0,cPermissive.stderr);
    const cPermissiveRow=await readRow(cPermissivePath);
    assert.ok(cPermissiveRow.c_call);assert.ok(Number(cPermissiveRow.c_prefix_identity)<.99);
    const cStrictPath=join(temporary,"c-strict.airr.tsv");
    const cStrict=runRawCli(root,[...cCommon,"--minimum-c-prefix-identity","0.99","-out",cStrictPath]);
    assert.equal(cStrict.status,0,cStrict.stderr);
    const cStrictRow=await readRow(cStrictPath);
    assert.equal(cStrictRow.c_call,"");
    for(const field of ["v_call","d_call","j_call","v_sequence_start","v_sequence_end","j_sequence_start","j_sequence_end"])assert.equal(cStrictRow[field],cPermissiveRow[field]);

    const optimizedPath=join(temporary,"r-optimized.airr.tsv");
    const optimized=runRawCli(root,[...common,"--calling-profile","r_optimized","-out",optimizedPath]);
    assert.equal(optimized.status,0,optimized.stderr);assert.match(optimized.stderr,/R-optimized profile/);
    const optimizedRow=await readRow(optimizedPath);
    assert.equal(optimizedRow.v_call,v[0]);assert.equal(optimizedRow.j_call,j[0]);

    const sensitivePath=join(temporary,"sensitive-d.airr.tsv");
    const sensitive=runRawCli(root,[...common,"--calling-profile","sensitive_d","-out",sensitivePath]);
    assert.equal(sensitive.status,0,sensitive.stderr);assert.match(sensitive.stderr,/Sensitive-D profile/);
    const sensitiveRow=await readRow(sensitivePath);
    assert.equal(sensitiveRow.v_call,v[0]);assert.equal(sensitiveRow.j_call,j[0]);

    // One search replaced the assignment strategies.
    const removedAssigner=runRawCli(root,[...common,"--assigner","riat_mp","-out",join(temporary,"removed-assigner.airr.tsv")]);
    assert.notEqual(removedAssigner.status,0);assert.match(removedAssigner.stderr,/Unsupported --vdj option --assigner/);
    const invalidProfile=runRawCli(root,[...common,"--calling-profile","igblast_balanced","-out",join(temporary,"invalid-profile.airr.tsv")]);
    assert.notEqual(invalidProfile.status,0);assert.match(invalidProfile.stderr,/--calling-profile must be/);

    const igblastPath=join(temporary,"igblast-data.airr.tsv");
    const igblast=runRawCli(root,[...common,"-custom_internal_data",internalPath,"-auxiliary_data",auxPath,"-out",igblastPath]);
    assert.equal(igblast.status,0,igblast.stderr);assert.match(igblast.stderr,/IgBLAST internal data/);assert.match(igblast.stderr,/IgBLAST auxiliary data/);
    const igblastRow=await readRow(igblastPath);
    assert.ok(igblastRow.cdr1);assert.ok(igblastRow.cdr3);assert.equal(igblastRow.region_definition,"IMGT");
    assert.equal(Number(igblastRow.j_sequence_end)-Number(igblastRow.fwr4_end),1);

    const swigPath=join(temporary,"swig.airr.tsv");
    const swig=runRawCli(root,[...common,"--swigannots","-organism","human","-ig_seqtype","Ig","-out",swigPath]);
    assert.equal(swig.status,0,swig.stderr);assert.match(swig.stderr,/Swig metadata preparation/);
    const swigRow=await readRow(swigPath);assert.ok(swigRow.cdr1);assert.ok(swigRow.cdr3);

    const multiQueryPath=join(temporary,"multi-query.fasta"),serialPath=join(temporary,"serial.airr.tsv"),parallelPath=join(temporary,"parallel.airr.tsv");
    const rearrangement=`${v[1]}AACCGG${d[1]}TTG${j[1]}`;
    await writeFile(multiQueryPath,Array.from({length:7},(_,index)=>`>multi-${index+1}\n${rearrangement}\n`).join(""));
    const references=["--vdj","-query",multiQueryPath,"-germline_db_V",vPath,"-germline_db_D",dPath,"-germline_db_J",jPath,"-outfmt","19"];
    const serial=runRawCli(root,[...references,"--workers","1","--batch-records","1","-out",serialPath]);
    const parallel=runRawCli(root,[...references,"--workers","3","--batch-records","0","-out",parallelPath]);
    assert.equal(serial.status,0,serial.stderr);assert.equal(parallel.status,0,parallel.stderr);
    assert.deepEqual(await readFile(parallelPath),await readFile(serialPath));
    assert.match(parallel.stderr,/1,000 records\/batch/);
    assert.deepEqual((await readdir(temporary)).filter((name)=>/summary|resolved-config|processed|lineage/i.test(name)),[]);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("--vdj applies the Sensitive-D prior on top of explicit IgBLAST-style D tuning",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-sensitive-d-"));
  try{
    const referenceRoot=join(root,"public/references/kimdb-1.1");
    const referencePaths={};
    for(const segment of ["V","D","J"]){
      const records=new Map();
      for(const text of await Promise.all(["Macaca_mulatta","Macaca_fascicularis"].map((species)=>
        readFile(join(referenceRoot,species,"IGH",`${segment}.fasta`),"utf8")))){
        for(const record of text.trim().split(/(?=>)/).filter(Boolean)){
          const name=record.slice(1).trim().split(/\s+/,1)[0];
          if(!records.has(name))records.set(name,record);
        }
      }
      referencePaths[segment]=join(temporary,`${segment}.fasta`);
      await writeFile(referencePaths[segment],`${[...records.values()].join("\n")}\n`);
    }
    const common=[
      "--vdj","-query",join(root,"tests/fixtures/sensitive-d-joint-threshold-regression.fasta"),
      "-germline_db_V",referencePaths.V,"-germline_db_D",referencePaths.D,
      "-germline_db_J",referencePaths.J,
      "--workers","1","--batch-records","1","-outfmt","19","-min_D_match","4",
    ];
    const annotate=async(profile)=>{
      const output=join(temporary,`${profile}.airr.tsv`);
      const result=runRawCli(root,[...common,"--calling-profile",profile,"-out",output]);
      assert.equal(result.status,0,result.stderr);
      const [header,line]=String(await readFile(output,"utf8")).trimEnd().split("\n");
      const names=header.split("\t"),values=line.split("\t");
      return Object.fromEntries(names.map((name,index)=>[name,values[index]??""]));
    };
    const optimized=await annotate("r_optimized");
    const sensitive=await annotate("sensitive_d");
    // A D twice as likely a priori shifts the presence posterior, not the call.
    for(const row of [optimized,sensitive]){assert.equal(row.d_call,"IGHD1-7*01");assert.equal(row.d_score,"12");}
    assert.ok(Number(sensitive.d_presence_posterior)>Number(optimized.d_presence_posterior));
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("prepare-reference persists diagnostics and a hash-checked manifest reusable by --vdj",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-prepared-reference-"));
  try{
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const human=pack.species.find((entry)=>entry.name==="Homo sapiens");
    const v=human.loci.IGH.V.find((entry)=>entry[2]?.slice(2,12).every((value)=>value>=0));
    const d=human.loci.IGH.D[0];
    const j=human.loci.IGH.J.find((entry)=>entry[2]?.[0]>=0&&entry[2]?.[1]>=0);
    assert.ok(v&&d&&j);
    const vPath=join(temporary,"V.fasta"),dPath=join(temporary,"D.fasta"),jPath=join(temporary,"J.fasta");
    const auxPath=join(temporary,"human_gl.aux"),queryPath=join(temporary,"query.fasta"),prefix=join(temporary,"prepared","human-igh");
    await writeFile(vPath,`>${v[0]}\n${v[1]}\n`);
    await writeFile(dPath,`>${d[0]}\n${d[1]}\n`);
    await writeFile(jPath,`>${j[0]}\n${j[1]}\n`);
    await writeFile(auxPath,`${j[0]} ${j[2][0]} JH ${j[2][1]} 1\n`);
    await writeFile(queryPath,`>read-1\n${v[1]}AACCGG${d[1]}TTG${j[1]}\n`);
    const prepared=runRawCli(root,["prepare-reference","-germline_db_V",vPath,"-germline_db_D",dPath,"-germline_db_J",jPath,"-auxiliary_data",auxPath,"-organism","human","-ig_seqtype","Ig","--match-mode","permissive","--out-prefix",prefix]);
    assert.equal(prepared.status,0,prepared.stderr);
    assert.match(prepared.stderr,/\[prepare:V\]/);
    const manifestPath=`${prefix}.swig-reference.json`;
    const manifest=JSON.parse(await readFile(manifestPath,"utf8"));
    assert.equal(manifest.application,"Swig prepared germline reference");
    assert.equal(manifest.match.mode,"permissive");
    assert.equal(manifest.complete,true);
    assert.equal(manifest.fwr4EndOffsets[j[0]],1);
    for(const segment of ["V","D","J"])assert.match(await readFile(`${prefix}.${segment}.fasta`,"utf8"),/^>/);
    const diagnostics=await readFile(`${prefix}.annotation-diagnostics.tsv`,"utf8");
    assert.match(diagnostics,/annotation_source/);assert.match(diagnostics,new RegExp(v[0].replace(/[.*+?^${}()|[\]\\]/g,"\\$&")));

    const output=join(temporary,"prepared.airr.tsv");
    const run=runRawCli(root,["--vdj","-query",queryPath,"--prepared-reference",manifestPath,"-out",output,"--workers","1","--batch-records","1"]);
    assert.equal(run.status,0,run.stderr);assert.match(run.stderr,/prepared-reference:V/);
    const [header,line]=String(await readFile(output,"utf8")).trimEnd().split("\n");
    const names=header.split("\t"),values=line.split("\t"),row=Object.fromEntries(names.map((name,index)=>[name,values[index]??""]));
    assert.ok(row.cdr1);assert.ok(row.cdr3);assert.equal(Number(row.j_sequence_end)-Number(row.fwr4_end),1);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("prepare-reference logs allele-level rejection reasons and can require complete metadata",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-prepared-failure-"));
  try{
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const human=pack.species.find((entry)=>entry.name==="Homo sapiens");
    const j=human.loci.IGH.J.find((entry)=>entry[2]?.[0]>=0&&entry[2]?.[1]>=0);assert.ok(j);
    const vPath=join(temporary,"V.fasta"),jPath=join(temporary,"J.fasta"),prefix=join(temporary,"failed","reference");
    await writeFile(vPath,`>IGHV_FAKE*01\n${"A".repeat(320)}\n`);
    await writeFile(jPath,`>${j[0]}\n${j[1]}\n`);
    const result=runRawCli(root,["prepare-reference","-germline_db_V",vPath,"-germline_db_J",jPath,"-organism","human","-ig_seqtype","Ig","--out-prefix",prefix,"--require-complete"]);
    assert.notEqual(result.status,0);assert.match(result.stderr,/\[prepare:V\] unresolved IGHV_FAKE\*01/);assert.match(result.stderr,/below_identity/);assert.match(result.stderr,/inspect .*annotation-diagnostics\.tsv/);
    const diagnostics=await readFile(`${prefix}.annotation-diagnostics.tsv`,"utf8");assert.match(diagnostics,/IGHV_FAKE\*01\tIGH\tunresolved/);assert.match(diagnostics,/below_identity:/);
    const manifest=JSON.parse(await readFile(`${prefix}.swig-reference.json`,"utf8"));assert.equal(manifest.complete,false);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("--vdj refuses to accumulate output when no destination is supplied",()=>{
  const root=resolve(import.meta.dirname,"..");
  const result=runRawCli(root,["--vdj","-germline_db_V","V.fasta","-germline_db_J","J.fasta"]);
  assert.notEqual(result.status,0);assert.match(result.stderr,/requires -out/);
});

test("--vdj refuses a D or J penalty outside the aligner's range",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const refs=join(root,"public/references/kimdb-1.1/Macaca_mulatta/IGH");
  const temporary=await mkdtemp(join(root,"tmp-cli-penalty-"));
  try{
    for(const option of ["-D_penalty","-J_penalty"]){
      const args=["--vdj","-query",join(root,"tests/fixtures/kimdb-perfect-d-regression.fasta"),"-germline_db_V",join(refs,"V.fasta"),"-germline_db_D",join(refs,"D.fasta"),"-germline_db_J",join(refs,"J.fasta"),option,"-1001","-out",join(temporary,"out.tsv")];
      const result=runRawCli(root,args);
      assert.equal(result.status,1,result.stderr);
      assert.match(result.stderr,new RegExp(`${option} must be a negative integer from -1000 to -1`));
      assert.equal(existsSync(join(temporary,"out.tsv")),false);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("--vdj refuses the C-prefix gate without C germlines on either aligner, and names D tuning without D",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const refs=join(root,"public/references/kimdb-1.1/Macaca_mulatta/IGH");
  const temporary=await mkdtemp(join(root,"tmp-cli-databases-"));
  try{
    const env=nativeAlignerEnvironment(root,join(temporary,"cache"))??process.env;
    const query=join(root,"tests/fixtures/kimdb-perfect-d-regression.fasta");
    const vj=["--vdj","-query",query,"-germline_db_V",join(refs,"V.fasta"),"-germline_db_J",join(refs,"J.fasta"),"-outfmt","19","--workers","1"];
    const messages=[];
    for(const aligner of ["native","wasm"]){
      const out=join(temporary,`${aligner}.tsv`);
      const result=runRawCli(root,[...vj,"-germline_db_D",join(refs,"D.fasta"),"--aligner",aligner,"--minimum-c-prefix-identity","0.99","-out",out],env);
      assert.notEqual(result.status,0,`--aligner ${aligner}: ${result.stderr}`);
      assert.match(result.stderr,/--minimum-c-prefix-identity needs a constant-region database: add -c_region_db/);
      assert.doesNotMatch(result.stderr,/did not start/);
      assert.equal(existsSync(out),false);
      messages.push(result.stderr.trim().split("\n").at(-1));
    }
    assert.equal(messages[0],messages[1]);
    // 0 turns the gate off.
    const off=runRawCli(root,[...vj,"-germline_db_D",join(refs,"D.fasta"),"--aligner","wasm","--minimum-c-prefix-identity","0","-out",join(temporary,"off.tsv")],env);
    assert.equal(off.status,0,off.stderr);
    // IgBLAST-style scripts pass D tuning for every chain: run, and say it changes nothing.
    const tuned=runRawCli(root,[...vj,"--aligner","wasm","-min_D_match","5","-out",join(temporary,"tuned.tsv")],env);
    assert.equal(tuned.status,0,tuned.stderr);
    assert.match(tuned.stderr,/Warning: -min_D_match has no effect without -germline_db_D\./);
    const plain=runRawCli(root,[...vj,"--aligner","wasm","-out",join(temporary,"plain.tsv")],env);
    assert.equal(plain.status,0,plain.stderr);
    assert.doesNotMatch(plain.stderr,/no effect/);
    assert.deepEqual(await readFile(join(temporary,"tuned.tsv")),await readFile(join(temporary,"plain.tsv")));
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("a pipeline config is refused when its reference set lacks the germlines an annotation option needs",async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-config-databases-"));
  try{
    let runs=0;
    const run=async(references,annotation)=>{
      const name=`run-${runs++}`,configPath=join(temporary,`${name}.json`);
      await writeFile(configPath,JSON.stringify({inputs:[{path:join(root,"tests/fixtures/cli-smoke.fasta"),sampleId:"a",subjectId:"d"}],references,annotation:{workers:1,...annotation},pipeline:{collapse:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:join(temporary,name),prefix:"run",writeLineageStudy:false}}));
      return runCli(root,configPath);
    };
    // Checked before the aligner is chosen: the same answer on every route.
    for(const aligner of ["auto","native","wasm"]){
      const screened=await run({species:"Homo sapiens",scope:"IGK"},{aligner,doubleD:{mode:"all"}});
      assert.notEqual(screened.status,0,screened.stderr);
      assert.match(screened.stderr,/annotation\.doubleD needs D germlines, but the Homo sapiens IGK reference set has none; set annotation\.doubleD\.mode to "off"\./);
      assert.doesNotMatch(screened.stderr,/runs only on WebAssembly/);
      const gated=await run({species:"Mus musculus",scope:"TRA"},{aligner,minimumConstantPrefixIdentity:0.97});
      assert.notEqual(gated.status,0,gated.stderr);
      assert.match(gated.stderr,/annotation\.minimumConstantPrefixIdentity needs constant-region \(C\) germlines, but the Mus musculus TRA reference set has none; set it to null/);
    }
    // A blank C text is no C database on either route: the native child is not handed an
    // empty C file (which it refuses), as WebAssembly is not.
    const env=nativeAlignerEnvironment(root,join(temporary,"cache"));
    if(env){
      const blank={species:"Homo sapiens",scope:"IGH",prepareMetadata:false,inline:{C:" \n"}};
      const tables=[];
      for(const aligner of ["native","wasm"]){
        const name=`blank-${aligner}`,configPath=join(temporary,`${name}.json`);
        await writeFile(configPath,JSON.stringify({inputs:[{path:join(root,"tests/fixtures/cli-smoke.fasta"),sampleId:"a",subjectId:"d"}],references:blank,annotation:{workers:1,aligner},pipeline:{collapse:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:join(temporary,name),prefix:"run",writeLineageStudy:false}}));
        const result=runCli(root,configPath,[],env);
        assert.equal(result.status,0,`${aligner}: ${result.stderr}`);
        tables.push(await readFile(join(temporary,name,"run.annotated.airr.tsv")));
      }
      assert.deepEqual(tables[0],tables[1]);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("CLI exposes the embedded IMGT data notice",()=>{
  const root=resolve(import.meta.dirname,"..");
  const result=runRawCli(root,["notices"]);
  assert.equal(result.status,0,result.stderr);
  assert.match(result.stdout,/IMGT\/GENE-DB reference data/);
  assert.match(result.stdout,/CC BY 4\.0/);
  assert.match(result.stdout,/Swig modifies the source data/);
  assert.match(result.stdout,/gki010/);
});

test("personalized-germline CLI streams gzip AIRR into the shared engine and exports a recovered allele",async()=>{
 const root=resolve(import.meta.dirname,".."),temporary=await mkdtemp(join(root,"tmp-cli-personalized-"));
 try{
  const reference='ACGT'.repeat(60),bases=[...reference];bases[39]='A';bases[79]='A';const query=bases.join('');
  const fasta=join(temporary,'V.fasta'),airr=join(temporary,'input.airr.tsv.gz'),out=join(temporary,'out');
  await writeFile(fasta,`>IGHV1-1*01\n${reference}\n`);
  const header=['sequence_id','subject_id','clone_id','v_call','v_germline_start','v_sequence_alignment','v_germline_alignment'];
  const rows=Array.from({length:24},(_,i)=>[`read${i}`,'donor',i+1,'IGHV1-1*01',1,query,reference].join('\t'));
  await writeFile(airr,gzipSync(header.join('\t')+'\n'+rows.join('\n')+'\n'));
  const result=runRawCli(root,['personalized-germline','--airr',airr,'--v-reference',fasta,'--out',out,'--max-candidates','96']);assert.equal(result.status,0,result.stderr);
  const dashboard=JSON.parse(await readFile(join(out,'personalized-germline.json'),'utf8'));
  assert.equal(dashboard.options.maximumNovelCandidatesPerGene,96);
  assert.equal(dashboard.representativeLineages,24);
  assert.ok(dashboard.pools.flatMap(p=>p.genes).flatMap(g=>g.activeAlleles).some(a=>!a.known&&a.sequence===query));
  assert.match(await readFile(join(out,'personalized_pool_1.V.fasta'),'utf8'),/SWIG_PERSONALIZED=novel/);
  assert.match(await readFile(join(out,'personalized-germline.tsv'),'utf8'),/expected_lineages/);
  const jointOut=join(temporary,'joint');const before=await readFile(join(out,'personalized-germline.json'),'utf8');
  const joint=runRawCli(root,['joint-germline','--airr',airr,'--v-reference',fasta,'--out',jointOut,'--max-candidates','16','--iterations','100']);assert.equal(joint.status,0,joint.stderr);
  const independent=JSON.parse(await readFile(join(jointOut,'joint-germline.json'),'utf8'));
  assert.equal(independent.mode,'joint-inherited-somatic-split');assert.equal(independent.options.maximumCandidates,16);
  assert.equal(independent.results.length,1);assert.equal(independent.results[0].trainingLineages+independent.results[0].testLineages,24);
  assert.equal(await readFile(join(out,'personalized-germline.json'),'utf8'),before);

 }finally{await rm(temporary,{recursive:true,force:true});}
});

/** Seeded synthetic rearrangements from a few human IGH alleles, `count` reads, 2% substitutions. */
function syntheticRearrangements(pack,count){
  const human=pack.species.find((entry)=>entry.name==="Homo sapiens");
  const V=human.loci.IGH.V.filter((entry)=>entry[2]?.slice(2,12).every((value)=>value>=0)).slice(0,12);
  const D=human.loci.IGH.D.slice(0,6),J=human.loci.IGH.J.slice(0,6);
  let state=12345;
  const random=()=>((state=Math.imul(state^(state>>>15),2246822519)+0x6d2b79f5|0)>>>0)/2**32;
  const pick=(items)=>items[Math.floor(random()*items.length)];
  const bases=(n)=>Array.from({length:n},()=>"ACGT"[Math.floor(random()*4)]).join("");
  const reads=[];
  for(let index=0;index<count;index+=1){
    const d=pick(D)[1];
    let read=`${pick(V)[1]}${bases(Math.floor(random()*6))}${d.slice(Math.floor(random()*4),d.length-Math.floor(random()*4))}${bases(Math.floor(random()*6))}${pick(J)[1]}`;
    read=[...read].map((base)=>random()<0.02?"ACGT"[Math.floor(random()*4)]:base).join("");
    reads.push(`>read_${index}\n${read}\n`);
  }
  return {reads:reads.join(""),references:{V:V.map(([name,sequence])=>`>${name}\n${sequence}\n`).join(""),D:D.map(([name,sequence])=>`>${name}\n${sequence}\n`).join(""),J:J.map(([name,sequence])=>`>${name}\n${sequence}\n`).join("")}};
}

test("--vdj runs the native aligner by default and writes the same AIRR as WebAssembly",async(t)=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-native-"));
  try{
    const env=nativeAlignerEnvironment(root,join(temporary,"cache"));
    if(!env){t.skip("no native aligner: build engine/ with cargo build --release");return;}
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    // More reads than the engine assigns per batch (16,384), so rows stream back while records are still sent.
    const {reads,references}=syntheticRearrangements(pack,17_500);
    const paths={query:join(temporary,"reads.fasta"),V:join(temporary,"V.fasta"),D:join(temporary,"D.fasta"),J:join(temporary,"J.fasta")};
    await writeFile(paths.query,reads);for(const segment of ["V","D","J"])await writeFile(paths[segment],references[segment]);
    const germlines=["--vdj","-germline_db_V",paths.V,"-germline_db_D",paths.D,"-germline_db_J",paths.J,"-outfmt","19"];
    const common=[...germlines,"--batch-records","700"];
    const nativePath=join(temporary,"native.tsv"),wasmPath=join(temporary,"wasm.tsv"),streamedPath=join(temporary,"streamed.tsv.gz");
    const native=runRawCli(root,[...common,"-query",paths.query,"-out",nativePath],env);
    assert.equal(native.status,0,native.stderr);
    assert.match(native.stderr,/native aligner; assignments-only; (\d+) threads?;/);
    assert.equal(Number(native.stderr.match(/(\d+) threads?;/)[1]),availableParallelism(),"the default uses every core");
    const wasm=runRawCli(root,[...common,"--aligner","wasm","--workers","2","-query",paths.query,"-out",wasmPath],env);
    assert.equal(wasm.status,0,wasm.stderr);assert.match(wasm.stderr,/WebAssembly aligner; assignments-only; 2 workers;/);
    const nativeTable=await readFile(nativePath);
    assert.equal(nativeTable.toString("utf8").trimEnd().split("\n").length,17_501);
    assert.deepEqual(nativeTable,await readFile(wasmPath));
    // Standard input and gzip output take the same route.
    const streamed=runRawCli(root,[...common,"--aligner","native","--workers","2","-query","-","-out",streamedPath],env,reads);
    assert.equal(streamed.status,0,streamed.stderr);
    assert.deepEqual(gunzipSync(await readFile(streamedPath)),nativeTable);
    // 17,500 = 17,000 + 500: the engine writes the first batch's last 616 rows only at the end of
    // its input, so the short last batch must reach it before the first batch is awaited.
    const tailPath=join(temporary,"tail.tsv");
    const tail=runRawCli(root,[...germlines,"--batch-records","17000","--workers","2","-query",paths.query,"-out",tailPath],env);
    assert.equal(tail.status,0,tail.stderr);assert.match(tail.stderr,/native aligner; assignments-only; 2 threads; 17,000 records\/batch;/);
    assert.deepEqual(await readFile(tailPath),nativeTable);
    // A standalone executable repairs an extracted copy that lost its execute bit.
    if(process.env.SWIG_CLI_EXECUTABLE&&process.platform!=="win32"){
      const cache=join(temporary,"cache","swig-cli");
      for(const version of await readdir(cache))await chmod(join(cache,version,"swig"),0o644);
      const repaired=runRawCli(root,[...common,"--workers","1","-query",paths.query,"-out",join(temporary,"repaired.tsv")],env);
      assert.equal(repaired.status,0,repaired.stderr);assert.match(repaired.stderr,/native aligner;/);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("native and WebAssembly routes name unnamed reads by input position",async(t)=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-native-unnamed-"));
  try{
    const env=nativeAlignerEnvironment(root,join(temporary,"cache"));
    if(!env){t.skip("no native aligner: build engine/ with cargo build --release");return;}
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const {reads,references}=syntheticRearrangements(pack,1_200);
    // Every third read has no identifier; each WebAssembly worker would otherwise number its own.
    let index=0;
    // ...and read_7's sequence carries a line " >x7": a sequence line to the aligner, which
    // must not turn into a header (an extra row) once the CLI drops its whitespace.
    const unnamed=reads.replace(/^>read_\d+$/gm,(header)=>(index++%3===0?">":header)).replace(/^(>read_7\n)/m,"$1 >x7\n");
    const paths={query:join(temporary,"reads.fasta"),V:join(temporary,"V.fasta"),D:join(temporary,"D.fasta"),J:join(temporary,"J.fasta")};
    await writeFile(paths.query,unnamed);for(const segment of ["V","D","J"])await writeFile(paths[segment],references[segment]);
    const run=(aligner)=>{
      const out=join(temporary,`${aligner}.tsv`);
      const result=runRawCli(root,["--vdj","--aligner",aligner,"--workers","2","--batch-records","100","-query",paths.query,"-germline_db_V",paths.V,"-germline_db_D",paths.D,"-germline_db_J",paths.J,"-outfmt","19","-out",out],env);
      assert.equal(result.status,0,result.stderr);
      return readFile(out,"utf8");
    };
    const native=await run("native"),wasm=await run("wasm");
    assert.equal(native,wasm);
    const ids=native.trimEnd().split("\n").slice(1).map((row)=>row.split("\t")[0]);
    assert.equal(ids.length,1_200);assert.equal(new Set(ids).size,1_200);
    assert.deepEqual([ids[0],ids[1],ids[3],ids[1_197]],["sequence_1","read_1","sequence_4","sequence_1198"]);
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("the pipeline annotates FASTA, FASTQ and AIRR datasets identically on both aligners",async(t)=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-native-run-"));
  try{
    const env=nativeAlignerEnvironment(root,join(temporary,"cache"));
    if(!env){t.skip("no native aligner: build engine/ with cargo build --release");return;}
    const fasta=await readFile(join(root,"tests/fixtures/cli-smoke.fasta"),"utf8");
    const records=fasta.trim().split(/(?=>)/).map((record)=>{const [header,...lines]=record.split(/\r?\n/);return {id:header.slice(1).trim().split(/\s+/)[0],sequence:lines.join("").trim()};});
    const fastaPath=join(temporary,"reads.fasta"),fastqPath=join(temporary,"reads.fastq"),airrPath=join(temporary,"reads.airr.tsv");
    // An unnamed record in FASTA and AIRR, and an AIRR row whose sequence is blank (no record).
    await writeFile(fastaPath,`${fasta.trimEnd()}\n>\n${records[0].sequence}\n`);
    await writeFile(fastqPath,records.map((r)=>`@${r.id}\n${r.sequence}\n+\n${"I".repeat(r.sequence.length)}\n`).join(""));
    await writeFile(airrPath,`sequence_id\tsequence\n${records.map((r)=>`${r.id}\t${r.sequence}`).join("\n")}\nblank\t  \n\t${records[1].sequence}\n`);
    // Names padded as the aligner trims them, no sequence_id column, rows longer than the header.
    const loosePath=join(temporary,"loose.airr.tsv");
    await writeFile(loosePath,`v_call\t sequence \n${records.map((r)=>`IGHV1-2*02\t${r.sequence}\textra`).join("\n")}\n`);
    // CSV whose values hold a tab, which must not add a column once the row is tab-separated.
    const csvPath=join(temporary,"tabbed.airr.csv");
    await writeFile(csvPath,`sequence_id,note,sequence\n${records.map((r)=>`${r.id},a\tb,${r.sequence}`).join("\n")}\n`);
    let runs=0;
    const annotate=async(aligner,doubleD="off",{references,environment=env}={})=>{
      const name=`run-${runs++}`,output=join(temporary,name),configPath=join(temporary,`${name}.json`);
      await writeFile(configPath,JSON.stringify({inputs:[{path:fastaPath,sampleId:"a",subjectId:"d"},{path:fastqPath,sampleId:"b",subjectId:"d"},{path:airrPath,format:"airr",sampleId:"c",subjectId:"d"},{path:csvPath,format:"airr",sampleId:"f",subjectId:"d"},{path:loosePath,format:"airr",sampleId:"e",subjectId:"d"}],...(references?{references}:{}),annotation:{workers:2,batchRecords:1,airrMode:"reannotate",aligner,doubleD:{mode:doubleD}},pipeline:{collapse:{enabled:false},lineage:{enabled:false},shm:{enabled:false}},output:{directory:output,prefix:"run",writeLineageStudy:false}}));
      const result=runCli(root,configPath,[],environment);
      return {result,table:result.status===0?await readFile(join(output,"run.annotated.airr.tsv")):null};
    };
    const native=await annotate("native"),wasm=await annotate("wasm");
    // One record per batch: AIRR repeats its header in every batch, which the native child reads once.
    assert.equal(native.result.status,0,native.result.stderr);assert.match(native.result.stderr,/Using the native aligner \(2 threads; 1 records\/batch\)/);
    assert.equal(wasm.result.status,0,wasm.result.stderr);assert.match(wasm.result.stderr,/Starting the WebAssembly aligner pool \(2 workers; 1 records\/batch\)/);
    assert.deepEqual(native.table,wasm.table);
    const ids=native.table.toString("utf8").trimEnd().split("\n").slice(1).map((row)=>row.split("\t")[0]);
    assert.equal(ids.length,records.length*5+2);
    assert.equal(ids.filter((id)=>id.endsWith(`::sequence_${records.length+1}`)).length,2);
    assert.deepEqual(ids.slice(-records.length),records.map((_,index)=>`loose::sequence_${index+1}`));
    // The double-D sidecar runs only on WebAssembly: auto switches, native refuses.
    const screened=await annotate("auto","all");
    assert.equal(screened.result.status,0,screened.result.stderr);assert.match(screened.result.stderr,/Using the WebAssembly aligner: the double-D screen runs only on it/);
    const refused=await annotate("native","all");
    assert.notEqual(refused.result.status,0);assert.match(refused.result.stderr,/annotation\.aligner "native": the double-D screen runs only on WebAssembly/);
    const required=await annotate(undefined,"all",{environment:{...env,SWIG_CLI_ALIGNER:"native"}});
    assert.notEqual(required.result.status,0);assert.match(required.result.stderr,/SWIG_CLI_ALIGNER=native: the double-D screen runs only on WebAssembly/);
    // A germline error stops either aligner (no fallback from native), and a WebAssembly pool
    // whose workers fail to start still lets the process exit.
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const badV=join(temporary,"bad-V.fasta");
    await writeFile(badV,`>IGHV-empty\n${syntheticRearrangements(pack,1).references.V}`);
    for(const aligner of ["auto","wasm"]){
      const bad=await annotate(aligner,"off",{references:{prepareMetadata:false,files:{V:badV}}});
      assert.notEqual(bad.result.status,null,`--aligner ${aligner} with a bad germline must exit`);assert.notEqual(bad.result.status,0);
      assert.match(bad.result.stderr,/V database: record 'IGHV-empty' has no sequence/);
      assert.doesNotMatch(bad.result.stderr,/did not start/);
    }
  }finally{await rm(temporary,{recursive:true,force:true});}
});

test("an unusable native aligner falls back to WebAssembly unless --aligner native insists",{skip:process.platform==="win32"},async()=>{
  const root=resolve(import.meta.dirname,"..");
  const temporary=await mkdtemp(join(root,"tmp-cli-native-fallback-"));
  try{
    const fake=join(temporary,"swig");await writeFile(fake,"#!/bin/sh\necho 'swig 0.0.0'\n");await chmod(fake,0o755);
    const env={...process.env,SWIG_CLI_CACHE_DIR:join(temporary,"cache"),SWIG_CLI_ALIGNER_PATH:fake};delete env.SWIG_CLI_ALIGNER;
    const fasta=join(root,"tests/fixtures/cli-smoke.fasta");
    const pack=JSON.parse(gunzipSync(await readFile(join(root,"public/references/imgt-202632-7-swig-0.7.json.gz"))).toString("utf8"));
    const {references}=syntheticRearrangements(pack,1);
    for(const segment of ["V","D","J"])await writeFile(join(temporary,`${segment}.fasta`),references[segment]);
    const args=(aligner)=>["--vdj","--aligner",aligner,"-query",fasta,"-germline_db_V",join(temporary,"V.fasta"),"-germline_db_D",join(temporary,"D.fasta"),"-germline_db_J",join(temporary,"J.fasta"),"--workers","1","-outfmt","19","-out",join(temporary,`${aligner}.tsv`)];
    const auto=runRawCli(root,args("auto"),env);
    assert.equal(auto.status,0,auto.stderr);
    const insisted=runRawCli(root,args("native"),env);
    if(process.env.SWIG_CLI_EXECUTABLE){
      // The executable's own binary still serves after the bad override.
      assert.match(auto.stderr,/native aligner;/);assert.equal(insisted.status,0,insisted.stderr);
    }else{
      assert.match(auto.stderr,/Using the WebAssembly aligner \(.*reports 'swig 0\.0\.0'/);assert.match(auto.stderr,/WebAssembly aligner;/);
      assert.notEqual(insisted.status,0);assert.match(insisted.stderr,/--aligner native: .*reports 'swig 0\.0\.0'/);
    }
    const invalid=runRawCli(root,args("gpu"),env);
    assert.notEqual(invalid.status,0);assert.match(invalid.stderr,/--aligner must be auto, native, wasm/);
    // A binary that passes the version check but cannot start (an illegal instruction, say).
    const version=JSON.parse(await readFile(join(root,"package.json"),"utf8")).version;
    const crashing=join(temporary,"crashing-swig");
    await writeFile(crashing,`#!/bin/sh\nif [ "$1" = "--version" ]; then echo 'swig ${version}'; exit 0; fi\necho 'Illegal instruction' >&2\nexit 132\n`);await chmod(crashing,0o755);
    const crashEnv={...env,SWIG_CLI_ALIGNER_PATH:crashing};
    const fellBack=runRawCli(root,args("auto"),crashEnv);
    assert.equal(fellBack.status,0,fellBack.stderr);
    assert.match(fellBack.stderr,/The native aligner did not start \(The native Swig aligner exited with status 132: Illegal instruction\); using WebAssembly\./);
    assert.match(fellBack.stderr,/WebAssembly aligner;/);
    const crashed=runRawCli(root,args("native"),crashEnv);
    assert.notEqual(crashed.status,0);assert.match(crashed.stderr,/swig-cli: The native Swig aligner exited with status 132: Illegal instruction/);
    // A binary that answers every record and then writes one row too many: the extra row comes
    // after the last batch was answered, and the run must still fail rather than drop a record.
    const surplus=join(temporary,"surplus-swig");
    await writeFile(surplus,`#!/bin/sh\nif [ "$1" = "--version" ]; then echo 'swig ${version}'; exit 0; fi\necho 'swig: database prepared in 1 ms' >&2\nprintf 'sequence_id\\tsequence\\n'\nwhile IFS= read -r line; do case "$line" in '>'*) printf '%s\\tACGT\\n' "\${line#>}";; esac; done\nprintf 'extra\\tACGT\\n'\n`);await chmod(surplus,0o755);
    const overflow=runRawCli(root,args("native"),{...env,SWIG_CLI_ALIGNER_PATH:surplus});
    assert.notEqual(overflow.status,0,overflow.stderr);assert.match(overflow.stderr,/swig-cli: The native aligner returned more rows than records\./);
    // A binary that exits with status 0 before its input ended has not assigned the records.
    const quitter=join(temporary,"quitting-swig");
    await writeFile(quitter,`#!/bin/sh\nif [ "$1" = "--version" ]; then echo 'swig ${version}'; exit 0; fi\necho 'swig: database prepared in 1 ms' >&2\nexit 0\n`);await chmod(quitter,0o755);
    const quit=runRawCli(root,args("native"),{...env,SWIG_CLI_ALIGNER_PATH:quitter});
    assert.notEqual(quit.status,0,quit.stderr);assert.match(quit.stderr,/swig-cli: The native aligner exited before its input ended\./);
    // A write error on the output (a full disk) ends the run with a message, not an unhandled event.
    if(existsSync("/dev/full")){
      const full=runRawCli(root,[...args("wasm").slice(0,-1),"/dev/full"],env);
      assert.notEqual(full.status,0);assert.match(full.stderr,/swig-cli: ENOSPC/);assert.doesNotMatch(full.stderr,/Unhandled 'error' event/);
    }
    // SWIG_CLI_ALIGNER: empty means unset; a bad value is named as such.
    const unset=runRawCli(root,args("auto").filter((value,index,all)=>value!=="--aligner"&&all[index-1]!=="--aligner"),{...env,SWIG_CLI_ALIGNER:""});
    assert.equal(unset.status,0,unset.stderr);
    const badEnvironment=runRawCli(root,args("auto").filter((value,index,all)=>value!=="--aligner"&&all[index-1]!=="--aligner"),{...env,SWIG_CLI_ALIGNER:"gpu"});
    assert.notEqual(badEnvironment.status,0);assert.match(badEnvironment.stderr,/SWIG_CLI_ALIGNER must be auto, native, wasm/);
    const configPath=join(temporary,"gpu.json");
    await writeFile(configPath,JSON.stringify({inputs:[{path:fasta,sampleId:"a",subjectId:"d"}],annotation:{workers:1,aligner:"gpu"},output:{directory:join(temporary,"gpu")}}));
    const badConfig=runCli(root,configPath,[],env);
    assert.notEqual(badConfig.status,0);assert.match(badConfig.stderr,/annotation\.aligner must be auto, native, wasm/);
  }finally{await rm(temporary,{recursive:true,force:true});}
});
