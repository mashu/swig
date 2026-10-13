import { runUnifiedGermline } from "./unified-germline.mjs";
import { runPersonalizedGermline } from "./personalized-germline.mjs";
import { createReadStream, createWriteStream, readFileSync, writeFileSync } from "node:fs";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import { basename, dirname, extname, isAbsolute, join, resolve } from "node:path";
import { createInterface } from "node:readline";
import { fileURLToPath } from "node:url";
import { Worker as NodeWorker } from "node:worker_threads";
import { createGunzip, createGzip, gzipSync, gunzipSync } from "node:zlib";
import { once } from "node:events";
import { availableParallelism } from "node:os";
import { Readable } from "node:stream";
import { pipeline } from "node:stream/promises";

import { SparseEvidenceAccumulator } from "../src/allele-refinement/evidence.ts";
import { refinedCall } from "../src/allele-refinement/export.ts";
import { toRefinementInputRow } from "../src/allele-refinement/input.ts";
import { fitSparseAlleleModel } from "../src/allele-refinement/model.ts";
import { buildReferenceAlleleGraph } from "../src/allele-refinement/reference-graph.ts";
import { CALLING_PROFILES, isCallingProfile } from "../src/aligner-engine.ts";
import { DEFAULT_CLI_CONFIG, normalizeCliConfig } from "../src/pipeline-config.ts";
import { MissingAlleleAccumulator } from "../src/germline-evidence.ts";
import {
  applyIgblastAuxiliaryData,
  applyIgblastDFrameData,
  applyIgblastInternalData,
  prepareIgblastStyleGermlineFasta,
  preprocessGermlineFastaAcrossTiers,
  resolveGermlineMatchOptions,
} from "../src/germline-preprocess.ts";
import {
  DenoiseAccumulator,
  assignLineages,
  createExactDedupPlan,
  deduplicate,
  finishExactDedupPlan,
  prepareReferenceMsa,
  runDenoisePartitionJob,
  runExactDedupJob,
} from "../src/post-analysis-core.ts";
import { airrRowToPostAnalysisRecord, denoiseVdjSequence } from "../src/post-analysis-record.ts";
import { compileReferences, germlineTemplateTiers, lociForScope } from "../src/reference-pack.ts";
import { repertoireRowMatches, validateRepertoireSelection } from "../src/repertoire-selection.ts";
import { ShmAccumulator } from "../src/shm-analysis.ts";
import {
  addFastqQualityFilterStats,
  airrRecords,
  airrRowNamer,
  canonicalFastq,
  emptyFastqQualityFilterStats,
  fastaRecords,
  fastqRecords,
  filterFastqRecord,
  namedFastaRecord,
  namedFastqRecord,
  seededRandom,
  validateFastqQualityFilter,
} from "../src/sequence-stream.ts";
import { annotateAirrBatch, annotateDoubleDBatch, stableDatasetSeed } from "../src/study-design.ts";
import { ALIGNER_MODES, NativePool, resolveNativeAligner } from "./native-aligner.mjs";

const VERSION="0.39.3";
const CLI_STREAM_HIGH_WATER_MARK=8*1024*1024;
const CLI_GZIP_CHUNK_SIZE=1024*1024;
const CLI_DIRECTORY=dirname(fileURLToPath(import.meta.url));

export function defaultCliAssets(){
  const directory=join(CLI_DIRECTORY,"assets");
  return {wasmPath:join(directory,"swig.wasm"),referencePackPath:join(directory,"imgt-reference-pack.json.gz")};
}

function usage(){
  return `swig-cli ${VERSION}\n\n`+
    `Run a complete non-phylogenetic Swig pipeline:\n`+
    `  swig-cli run reads.fastq.gz --out swig-output\n`+
    `  swig-cli run --config swig.config.json [--out DIRECTORY] [--workers N] [--airr-compression none|gzip]\n`+
    `      [--aligner auto|native|wasm]\n\n`+
    `Experimental joint inherited/SHM model: swig-cli joint-germline --help\n\nInfer personalized V alleles from lineage-assigned AIRR:\n  swig-cli personalized-germline --airr processed.airr.tsv.gz --v-reference V.fasta --out germline\n\nRun only streaming V(D)J assignment (AIRR outfmt 19):\n`+
    `  swig-cli --vdj -query reads.fasta -germline_db_V V.fasta -germline_db_D D.fasta \\\n`+
    `    -germline_db_J J.fasta -out calls.airr.tsv\n\n`+
    `Prepare custom germlines once and reuse their inferred annotations:\n`+
    `  swig-cli prepare-reference -germline_db_V V.fasta -germline_db_D D.fasta \\\n`+
    `    -germline_db_J J.fasta -organism human -ig_seqtype Ig --out-prefix refs/custom\n\n`+
    `Display bundled-data attribution and license:\n`+
    `  swig-cli notices\n\n`+
    `Create an editable config:\n`+
    `  swig-cli init swig.config.json\n\n`+
    `Single-input metadata options:\n`+
    `  --sample SAMPLE_ID  --donor SUBJECT_ID  --dataset DATASET_ID\n\n`+
    `Samples with the same subjectId/--donor are treated as the same donor.\n`+
    `The aligner runs natively where this executable carries the platform's binary (--aligner auto,\n`+
    `the default), else as WebAssembly; --aligner wasm forces WebAssembly, --aligner native fails\n`+
    `instead of falling back. Double-D screening always uses WebAssembly.\n`+
    `Lineage phylogenetics is intentionally not run by swig-cli.`;
}

function prepareReferenceUsage(){
  return `swig-cli ${VERSION} prepare-reference\n\n`+
    `Infer, validate, and persist reusable SWIGMETA germline annotations.\n\n`+
    `Required:\n`+
    `  -germline_db_V FASTA  -germline_db_J FASTA  --out-prefix PREFIX\n\n`+
    `Optional references and exact metadata:\n`+
    `  -germline_db_D FASTA  -c_region_db FASTA\n`+
    `  -custom_internal_data FILE  -auxiliary_data FILE  -d_frame_data FILE\n`+
    `  -organism NAME (default human)  -ig_seqtype Ig|TCR (default Ig)\n\n`+
    `Matching controls:\n`+
    `  --match-mode strict|permissive|best-guess  (default strict)\n`+
    `  --best-guess             Alias for --match-mode best-guess; disables identity floors\n`+
    `  --nearest-candidates N   Non-gene candidates aligned after named candidates fail\n`+
    `  --v-same-gene-min-identity X  --v-nearest-min-identity X\n`+
    `  --j-same-gene-min-identity X  --j-nearest-min-identity X\n`+
    `  --require-complete       Exit nonzero if any V/J record remains unresolved\n\n`+
    `Outputs are PREFIX.V/D/J/C.fasta, PREFIX.swig-reference.json, and\n`+
    `PREFIX.annotation-diagnostics.tsv. The manifest can be passed directly to\n`+
    `swig-cli --vdj with --prepared-reference.`;
}

function vdjUsage(){
  return `swig-cli ${VERSION} --vdj\n\n`+
    `Low-overhead, streaming V(D)J assignment by the Swig aligner with IgBLAST-style option names.\n\n`+
    `Required:\n`+
    `  -out AIRR_TSV, plus either --prepared-reference MANIFEST or\n`+
    `  -germline_db_V FASTA and -germline_db_J FASTA\n\n`+
    `Input and optional references:\n`+
    `  -query FASTA            Query FASTA or '-' for stdin (default '-')\n`+
    `  -germline_db_D FASTA    D germline FASTA\n`+
    `  -c_region_db FASTA      Constant-region FASTA\n\n`+
    `Annotation modes (default: the aligner's own annotation from conserved motifs and annotated\n`+
    `relatives of the ungapped FASTA; v_annotation_source and j_annotation_source name each source):\n`+
    `  -custom_internal_data FILE  IgBLAST V .ndm.imgt data (1-based inclusive intervals)\n`+
    `  -auxiliary_data FILE        IgBLAST J .aux data (0-based frame/CDR3 stop)\n`+
    `  -d_frame_data FILE          IgBLAST D frame-one starts\n`+
    `  --swigannots                Infer/validate metadata as in Swig Web\n\n`+
    `  --prepared-reference FILE   Reuse a prepare-reference manifest and its FASTAs\n`+
    `  --match-mode MODE           Metadata transfer only: strict, permissive, best-guess\n`+
    `  --best-guess                Disable metadata-transfer identity floors (not read mapping)\n`+
    `Execution:\n`+
    `  -num_threads N          Exact thread count; --workers N overrides it\n`+
    `  --workers N             Exact threads (or WebAssembly workers); 0 or omitted uses every core\n`+
    `  --aligner MODE          auto (default): the native aligner this executable carries, else\n`+
    `                          WebAssembly; native: fail instead of falling back; wasm: WebAssembly\n`+
    `  --batch-records N       Records per bounded batch; 0/omitted selects 2000, 1000,\n`+
    `                          or 500 according to worker count\n`+
    `  --calling-profile NAME  r_optimized (default), sensitive_d, truth_optimized or igblast_compatible\n`+
    `    r_optimized: 4-nt D exact-run floor; a D pays its prior cost, derived from the D genes of\n`+
    `                 the read's locus (10.55 score units for human IGH)\n`+
    `    sensitive_d: the same scores with a D twice as likely a priori (9.45 for human IGH)\n`+
    `  C is anchored at the called J's end and called at an anchored E-value of at most 1e-3\n`+
    `  --minimum-c-prefix-identity X  Optional 0..1 identity gate over the complete covered C prefix\n`+
    `  --minimum-identity X    Drop V, D and J hits whose local alignment identity is below X (default 0.6)\n`+
    `  -min_D_match N  -D_penalty S  -J_penalty S  -min_J_length N   override the profile's D/J settings\n`+
    `  -num_alignments_D/J N   accepted for IgBLAST compatibility; every admissible D and the best J\n`+
    `                          alignments enter adjudication, and near-tied alleles are in *_alternatives\n`+
    `  -strand both|plus|minus -outfmt 19 -organism NAME -ig_seqtype Ig|TCR\n\n`+
    `The germline options take FASTA files, not makeblastdb binary prefixes. Output is Swig AIRR,\n`+
    `not IgBLAST pairwise/tabular formatting. The output path is mandatory and is written incrementally;\n`+
    `a path ending in .gz is gzip-compressed.`;
}

function thirdPartyNotices(){
  return `Bundled IMGT/GENE-DB reference data\n\n`+
    `Source: IMGT/GENE-DB release 202632-7, retrieved 2026-08-08.\n`+
    `Copyright © 1995-2026 IMGT®, the international ImMunoGeneTics information system®.\n`+
    `Attribution: IMGT®, the international ImMunoGeneTics information system®, https://www.imgt.org/, Institute of Human Genetics, Université de Montpellier and CNRS.\n`+
    `License: CC BY 4.0, https://creativecommons.org/licenses/by/4.0/\n`+
    `Terms: https://www.imgt.org/about/termsofuse.php\n`+
    `Citation: Giudicelli V, Chaume D, Lefranc M-P. Nucleic Acids Research. 2005;33:D593-D597. https://doi.org/10.1093/nar/gki010\n\n`+
    `Swig modifies the source data by selecting and reorganizing IG/TR V/D/J/C records, normalizing and ungapping nucleotide sequences, deriving compact coordinate metadata, selecting one source sequence per allele identifier, and joining selected coding IGH/TR constant exons. Membrane-only and untranslated constant exons are omitted. IMGT, Université de Montpellier, and CNRS do not endorse Swig or warrant the modified pack or its use.`;
}

function argumentValue(args,name){const index=args.indexOf(name);if(index>=0)return args[index+1];const inline=args.find((value)=>value.startsWith(`${name}=`));return inline?.slice(name.length+1);}
function hasFlag(args,name){return args.includes(name);}
function positional(args){const value=[];for(let i=0;i<args.length;i+=1){if(args[i].startsWith("--")){if(!args[i].includes("=")&&!["--help","--version"].includes(args[i]))i+=1;continue;}value.push(args[i]);}return value;}
function cleanCell(value){return String(value??"").replace(/[\t\r\n]+/g," ");}

function parseIntegerOption(value,label,{minimum=0,allowZero=true}={}){
  if(value===undefined||!/^\d+$/.test(value))throw new Error(`${label} requires an integer value.`);
  const parsed=Number(value);
  if(!Number.isSafeInteger(parsed)||parsed<minimum||(!allowZero&&parsed===0))throw new Error(`${label} has an invalid value: ${value}.`);
  return parsed;
}

function parseFiniteOption(value,label){
  const parsed=Number(value);if(value===undefined||!Number.isFinite(parsed))throw new Error(`${label} requires a numeric value.`);return parsed;
}

function parseVdjArguments(rawArgs,context="vdj"){
  const aliases=new Map([
    ["--query","-query"],["--output","-out"],["--out","-out"],
    ["--germline-db-v","-germline_db_V"],["--germline-db-d","-germline_db_D"],
    ["--germline-db-j","-germline_db_J"],["--c-region-db","-c_region_db"],
    ["--min-c-prefix-identity","--minimum-c-prefix-identity"],
    ["--out-prefix","-out"],
  ]);
  const valued=new Set([
    "-query","-out","-germline_db_V","-germline_db_D","-germline_db_J","-c_region_db",
    "-custom_internal_data","-auxiliary_data","-d_frame_data","-organism","-domain_system",
    "-ig_seqtype","-strand","-outfmt","-num_threads","--workers","--batch-records",
    "--minimum-identity","--minimum-c-prefix-identity","--calling-profile","-min_D_match","-min_J_length",
    "-num_alignments_D","-num_alignments_J","-D_penalty","-J_penalty",
    "--prepared-reference","--match-mode","--nearest-candidates",
    "--v-same-gene-min-identity","--v-nearest-min-identity",
    "--j-same-gene-min-identity","--j-nearest-min-identity","--aligner",
  ]);
  const flags=new Set(["--swigannots","-show_translation","--best-guess","--require-complete"]);
  const options={};
  for(let index=0;index<rawArgs.length;index+=1){
    let token=rawArgs[index];
    if(["--vdj","--precompute_aux","--precompute-aux"].includes(token))continue;
    if(["-h","-help","--help"].includes(token)){options.help=true;continue;}
    if(["-version","--version"].includes(token)){options.version=true;continue;}
    const equals=token.indexOf("=");let inline;
    if(equals>0){inline=token.slice(equals+1);token=token.slice(0,equals);}
    token=aliases.get(token)??token;
    if(flags.has(token)){
      if(inline!==undefined)throw new Error(`${token} does not take a value.`);
      options[token]=true;continue;
    }
    if(!valued.has(token))throw new Error(`Unsupported ${context==="prepare"?"prepare-reference":"--vdj"} option ${token}.\n\n${context==="prepare"?prepareReferenceUsage():vdjUsage()}`);
    const value=inline!==undefined?inline:rawArgs[++index];
    if(value===undefined)throw new Error(`${token} requires a value.`);
    options[token]=value;
  }
  return options;
}

function germlineMatchOptions(options,{diagnostics=false}={}){
  let mode=String(options["--match-mode"]??"strict").replaceAll("-","_");
  if(options["--best-guess"]){
    if(options["--match-mode"]&&mode!=="best_guess")throw new Error("--best-guess conflicts with a different --match-mode value.");
    mode="best_guess";
  }
  const optionalIdentity=(name)=>options[name]===undefined?undefined:parseFiniteOption(options[name],name);
  return resolveGermlineMatchOptions({
    mode,
    nearestCandidates:options["--nearest-candidates"]===undefined?undefined:parseIntegerOption(options["--nearest-candidates"],"--nearest-candidates",{minimum:1,allowZero:false}),
    vSameGeneMinIdentity:optionalIdentity("--v-same-gene-min-identity"),
    vNearestMinIdentity:optionalIdentity("--v-nearest-min-identity"),
    jSameGeneMinIdentity:optionalIdentity("--j-same-gene-min-identity"),
    jNearestMinIdentity:optionalIdentity("--j-nearest-min-identity"),
    includeDiagnostics:diagnostics,
  });
}

function resolveFrom(base,value){return isAbsolute(value)?value:resolve(base,value);}

function detectFormat(path,requested="auto"){
  if(requested!=="auto")return requested;
  const plain=path.replace(/\.gz$/i,"");
  const extension=extname(plain).toLowerCase();
  if([".fa",".fasta",".fna",".fas"].includes(extension))return "fasta";
  if([".fq",".fastq"].includes(extension))return "fastq";
  if([".tsv",".airr"].includes(extension)||plain.toLowerCase().endsWith(".airr.tsv"))return "airr";
  throw new Error(`Cannot infer the input format for ${path}; set format to fasta, fastq, or airr in the config.`);
}

function inputLines(input){
  const range=input.gzipRange;
  if(typeof input.inline==="string"){
    if(range)throw new Error(`${input.path} cannot combine inline input with gzipRange.`);
    return createInterface({input:Readable.from([input.inline]),crlfDelay:Infinity});
  }
  if(input.path==="-"&&range)throw new Error("Standard input cannot be combined with gzipRange.");
  if(range&&(!Number.isSafeInteger(range.start)||!Number.isSafeInteger(range.end)||range.start<0||range.end<=range.start))throw new Error(`${input.path} has an invalid gzipRange.`);
  const readOptions=range
    ? {start:range.start,end:range.end-1,highWaterMark:CLI_GZIP_CHUNK_SIZE}
    : {highWaterMark:CLI_GZIP_CHUNK_SIZE};
  const raw=input.path==="-"?process.stdin:createReadStream(input.path,readOptions);
  const stream=/\.gz$/i.test(input.path)||range?raw.pipe(createGunzip({chunkSize:CLI_GZIP_CHUNK_SIZE})):raw;
  return createInterface({input:stream,crlfDelay:Infinity});
}

async function* sequenceBatches(input,batchRecords,preprocessing,airrMode,datasetIndex,state){
  const format=detectFormat(input.path,input.format);
  const lines=inputLines(input);
  const filter=preprocessing.fastqFilter;
  validateFastqQualityFilter(filter);
  state.fastqFilter=emptyFastqQualityFilterStats(filter.enabled,filter.enabled&&format==="fastq");
  // A record without an identifier is named sequence_N by its position N in the input before it
  // is batched, as Swig Web names it (unnamedSequenceId).
  const parsed=async function*(){
    if(format==="airr"){
      let namer=null;
      for await(const record of airrRecords(lines)){
        state.inputRecords+=1;state.eligibleRecords+=1;
        if(filter.enabled){state.fastqFilter.recordsRetained+=1;state.fastqFilter.recordsPassedThrough+=1;}
        const delimiter=record.header.includes("\t")?"\t":",";
        let header=delimiter==="\t"?record.header:record.header.split(delimiter).join("\t");
        // A tab inside a CSV value would add a column once the row is tab-separated.
        let row=delimiter==="\t"?record.row:record.row.split(delimiter).map((value)=>value.replace(/\t/g," ")).join("\t");
        if(airrMode!=="preserve"){
          namer??=airrRowNamer(header,"\t");
          row=namer.name(row,state.inputRecords);header=namer.header;
        }
        yield {ordinal:state.inputRecords-1,text:`${row}\n`,header};
      }
      return;
    }
    if(format==="fasta"){
      for await(const record of fastaRecords(lines)){
        state.inputRecords+=1;state.eligibleRecords+=1;
        if(filter.enabled){state.fastqFilter.recordsRetained+=1;state.fastqFilter.recordsPassedThrough+=1;}
        yield {ordinal:state.inputRecords-1,text:namedFastaRecord(record,state.inputRecords)};
      }
      return;
    }
    for await(const read of fastqRecords(lines)){
      state.inputRecords+=1;
      const record=namedFastqRecord(read,state.inputRecords);
      const text=filter.enabled?filterFastqRecord(record,filter,state.fastqFilter):canonicalFastq(record);
      if(text!==null){state.eligibleRecords+=1;yield {ordinal:state.inputRecords-1,text};}
    }
  };
  let selected=parsed();
  if(preprocessing.subsample.enabled){
    const size=preprocessing.subsample.size;
    const random=seededRandom(stableDatasetSeed(preprocessing.subsample.seed,datasetIndex));
    const reservoir=[];
    for await(const record of selected){
      if(reservoir.length<size)reservoir.push(record);
      else{const replacement=Math.floor(random()*state.eligibleRecords);if(replacement<size)reservoir[replacement]=record;}
    }
    reservoir.sort((left,right)=>left.ordinal-right.ordinal);
    selected=(async function*(){yield* reservoir;})();
  }
  let batch=[];
  const emit=()=>{
    const header=format==="airr"?batch[0]?.header:undefined;
    const body=batch.map((record)=>record.text).join("");
    const result=format==="airr"&&airrMode==="preserve"
      ? {format,count:batch.length,header,body}
      : {format,count:batch.length,text:header?`${header}\n${body}`:body};
    batch=[];
    return result;
  };
  for await(const record of selected){batch.push(record);state.selectedRecords+=1;if(batch.length>=batchRecords)yield emit();}
  if(batch.length)yield emit();
  if(!state.inputRecords)throw new Error(`No sequence records were found in ${input.path}.`);
  if(!state.eligibleRecords&&state.fastqFilter.applicable)throw new Error(`The FASTQ quality filter rejected all ${state.inputRecords.toLocaleString()} reads in ${input.path}.`);
}

/**
 * A Bun worker for the CLI entry `name` (cli-src/<name>.mjs). A standalone executable embeds it as
 * <name>.js and finds it from a relative specifier. Given a file URL there instead, Bun 1.3 on
 * Windows finds the embedded entry as B:\~BUN\root\<name>.js but then looks that spelling up
 * verbatim and reads it from disk (ENOENT).
 */
function bunWorker(name){
  const standalone=/\$bunfs|~BUN|%7EBUN/i.test(import.meta.url);
  return new globalThis.Worker(standalone?`./${name}.mjs`:new URL(`./${name}.js`,import.meta.url));
}

class WasmPool {
  constructor(size,init){this.size=size;this.init=init;this.workers=[];this.available=[];this.queued=[];this.pending=new Map();this.nextId=1;this.failed=null;}
  async start(){
    for(let index=0;index<this.size;index+=1){
      const webWorker=Boolean(process.versions.bun&&globalThis.Worker);
      const worker=webWorker?bunWorker("swig-worker"):new NodeWorker(new URL("./swig-worker.mjs",import.meta.url));
      const receive=(message)=>{
        const pending=this.pending.get(message.id);if(!pending)return;
        this.pending.delete(message.id);
        const {id,error,...result}=message;
        if(error)pending.reject(new Error(error));else pending.resolve(result);
        if(pending.release)this.release(worker);
      };
      const fail=(error)=>{
        const failure=error instanceof Error?error:new Error(error?.message??String(error));this.failed=failure;
        for(const [id,pending] of this.pending){this.pending.delete(id);pending.reject(failure);}
        while(this.queued.length)this.queued.shift().reject(failure);
      };
      if(webWorker){worker.addEventListener("message",(event)=>receive(event.data));worker.addEventListener("error",fail);}
      else{worker.on("message",receive);worker.on("error",fail);}
      this.workers.push(worker);
    }
    // One compilation shared by every worker: per-worker compiles compete with the workers for the
    // CPU (under Bun 1.3, four workers ran no faster than two).
    const module=await WebAssembly.compile(readFileSync(this.init.wasmPath));
    await Promise.all(this.workers.map((worker)=>this.request(worker,{type:"init",...this.init,module})));
    this.available.push(...this.workers);
  }
  request(worker,message,release=false){
    const id=this.nextId++;
    return new Promise((resolvePromise,reject)=>{this.pending.set(id,{resolve:resolvePromise,reject,worker,release});worker.postMessage({id,...message});});
  }
  dispatch(worker,job){
    const id=this.nextId++;
    this.pending.set(id,{resolve:job.resolve,reject:job.reject,worker,release:true});
    worker.postMessage({id,type:"annotate",...job.message});
  }
  release(worker){const job=this.queued.shift();if(job)this.dispatch(worker,job);else this.available.push(worker);}
  run(message){
    const promise=this.failed?Promise.reject(this.failed):new Promise((resolvePromise,reject)=>{
      const job={message,resolve:resolvePromise,reject};const worker=this.available.shift();
      if(worker)this.dispatch(worker,job);else this.queued.push(job);
    });
    // The caller awaits batches later, in order, and stops at the first failure: the batches
    // queued behind it (rejected when the pool closes) must not be reported as unhandled.
    promise.catch(()=>{});
    return promise;
  }
  async close(){
    const error=new Error("The aligner's worker pool closed before queued work completed.");
    while(this.queued.length)this.queued.shift().reject(error);
    for(const worker of this.workers)await worker.terminate();
    this.workers=[];this.available=[];
  }
}

async function runPostAnalysisTasks(tasks,requestedWorkers,onCompleted){
  if(!tasks.length)return [];
  const workerCount=Math.max(1,Math.min(Math.floor(requestedWorkers)||1,tasks.length));
  if(workerCount===1)return tasks.map((task)=>{const result=task.kind==="denoise"?runDenoisePartitionJob(task.job):runExactDedupJob(task.job);onCompleted?.(result);return result;});
  return new Promise((resolvePromise,reject)=>{
    const workers=[];const results=new Array(tasks.length);let next=0,completed=0,settled=false;
    const close=()=>{for(const worker of workers)void worker.terminate();};
    const fail=(error)=>{if(settled)return;settled=true;close();reject(error instanceof Error?error:new Error(error?.message??String(error)));};
    const dispatch=(worker)=>{if(next>=tasks.length)return;const id=next++;worker.postMessage({id,...tasks[id]});};
    for(let index=0;index<workerCount;index+=1){
      const webWorker=Boolean(process.versions.bun&&globalThis.Worker);
      const worker=webWorker?bunWorker("post-analysis-worker"):new NodeWorker(new URL("./post-analysis-worker.mjs",import.meta.url));
      const receive=(message)=>{if(settled)return;if(message.error||!message.result){fail(new Error(message.error||"A post-analysis worker returned no result."));return;}results[message.id]=message.result;onCompleted?.(message.result);completed+=1;if(completed===tasks.length){settled=true;close();resolvePromise(results);}else dispatch(worker);};
      if(webWorker){worker.addEventListener("message",(event)=>receive(event.data));worker.addEventListener("error",fail);}
      else{worker.on("message",receive);worker.on("error",fail);}
      workers.push(worker);dispatch(worker);
    }
  });
}

class ChmmPool {
  constructor(size,init){this.size=size;this.init=init;this.workers=[];this.available=[];this.queued=[];this.pending=new Map();this.nextId=1;this.failed=null;}
  async start(){
    for(let index=0;index<this.size;index+=1){
      const webWorker=Boolean(process.versions.bun&&globalThis.Worker);
      const worker=webWorker?bunWorker("chmmairra-worker"):new NodeWorker(new URL("./chmmairra-worker.mjs",import.meta.url));
      const receive=(message)=>{const pending=this.pending.get(message.id);if(!pending)return;this.pending.delete(message.id);if(message.error)pending.reject(new Error(message.error));else pending.resolve(message.result);if(pending.release)this.release(worker);};
      const fail=(error)=>{const failure=error instanceof Error?error:new Error(error?.message??String(error));this.failed=failure;for(const [id,pending] of this.pending){this.pending.delete(id);pending.reject(failure);}while(this.queued.length)this.queued.shift().reject(failure);};
      if(webWorker){worker.addEventListener("message",(event)=>receive(event.data));worker.addEventListener("error",fail);}else{worker.on("message",receive);worker.on("error",fail);}
      this.workers.push(worker);
    }
    await Promise.all(this.workers.map((worker)=>this.request(worker,{type:"init",...this.init})));
    this.available.push(...this.workers);
  }
  request(worker,message,release=false){const id=this.nextId++;return new Promise((resolvePromise,reject)=>{this.pending.set(id,{resolve:resolvePromise,reject,release});worker.postMessage({id,...message});});}
  dispatch(worker,job){const id=this.nextId++;this.pending.set(id,{resolve:job.resolve,reject:job.reject,release:true});worker.postMessage({id,type:"batch",rows:job.rows});}
  release(worker){const job=this.queued.shift();if(job)this.dispatch(worker,job);else this.available.push(worker);}
  run(rows){if(this.failed)return Promise.reject(this.failed);return new Promise((resolvePromise,reject)=>{const job={rows,resolve:resolvePromise,reject};const worker=this.available.shift();if(worker)this.dispatch(worker,job);else this.queued.push(job);});}
  async close(){const error=new Error("CHMMAIRRa worker pool closed before queued work completed.");while(this.queued.length)this.queued.shift().reject(error);await Promise.all(this.workers.map((worker)=>Promise.resolve(worker.terminate())));this.workers=[];this.available=[];}
}

function workerInitialization(wasmPath,references,settings){
  return {wasmPath,references:{V:references.V??"",D:references.D??"",J:references.J??"",C:references.C??""},settings};
}

function workerAnnotation(text,count){
  return {text,count};
}

function automaticBatchRecords(workers){return workers<=2?2_000:workers<=4?1_000:500;}

function parseTable(headerText,bodyText){
  const headers=headerText.replace(/\r$/,"").split("\t");
  const rows=[];
  for(const line of bodyText.split(/\r?\n/)){if(!line)continue;const fields=line.split("\t");const values={};headers.forEach((header,index)=>{values[header]=fields[index]??"";});rows.push(values);}
  return {headers,rows};
}

function trimAuxiliaryFwr4(headerText,bodyText,endOffsets,jLengths){
  if(!Object.keys(endOffsets).length)return bodyText;
  const headers=headerText.split("\t");const at=new Map(headers.map((header,index)=>[header,index]));
  const required=["j_call","j_sequence_start","j_germline_start","j_sequence_alignment","j_germline_alignment","fwr4","fwr4_aa","fwr4_end"];
  if(required.some((header)=>!at.has(header)))return bodyText;
  const lines=[];
  for(const line of bodyText.split(/\r?\n/)){
    if(!line)continue;
    const values=line.split("\t");const call=(values[at.get("j_call")]??"").split(",",1)[0];
    const offset=endOffsets[call],referenceLength=jLengths.get(call);
    if(offset===undefined||!referenceLength||offset===0){lines.push(line);continue;}
    const queryStart=Number(values[at.get("j_sequence_start")]);
    const referenceStart=Number(values[at.get("j_germline_start")]);
    const queryAlignment=values[at.get("j_sequence_alignment")]??"";
    const germlineAlignment=values[at.get("j_germline_alignment")]??"";
    const currentEnd=Number(values[at.get("fwr4_end")]);
    const targetReferenceEnd=referenceLength-offset;
    if(!Number.isFinite(queryStart)||!Number.isFinite(referenceStart)||!Number.isFinite(currentEnd)||targetReferenceEnd<referenceStart-1){lines.push(line);continue;}
    let queryPosition=queryStart-1,referencePosition=referenceStart-1;
    const columns=Math.min(queryAlignment.length,germlineAlignment.length);
    for(let column=0;column<columns;column+=1){
      if(referencePosition>=targetReferenceEnd)break;
      if(queryAlignment[column]!=="-")queryPosition+=1;
      if(germlineAlignment[column]!=="-")referencePosition+=1;
    }
    if(referencePosition<targetReferenceEnd||queryPosition>=currentEnd){lines.push(line);continue;}
    const trim=currentEnd-queryPosition;
    const nucleotide=values[at.get("fwr4")]??"";
    const kept=Math.max(0,nucleotide.length-trim);
    values[at.get("fwr4")]=nucleotide.slice(0,kept);
    values[at.get("fwr4_aa")]=(values[at.get("fwr4_aa")]??"").slice(0,Math.floor(kept/3));
    values[at.get("fwr4_end")]=String(queryPosition);
    lines.push(values.join("\t"));
  }
  return lines.length?`${lines.join("\n")}\n`:"";
}

function serializeRowBody(headers,rows){return `${rows.map((row)=>headers.map((header)=>cleanCell(row[header])).join("\t")).join("\n")}${rows.length?"\n":""}`;}
function addHeaders(headers,names){for(const name of names)if(!headers.includes(name))headers.push(name);}

// A user's file reads as Swig Web reads it (TextDecoder): a UTF-8 byte-order mark is dropped,
// a byte that is no UTF-8 reads as U+FFFD. Buffer's toString keeps the mark before the first header.
function decodeText(bytes){return new TextDecoder().decode(bytes);}

async function loadReferences(config,base,referencePackPath){
  const packed=JSON.parse(gunzipSync(readFileSync(referencePackPath)).toString("utf8"));
  const species=packed.species.find((entry)=>entry.name===config.references.species);
  if(!species)throw new Error(`The bundled reference pack has no exact species named ${config.references.species}.`);
  const overrides={...config.references.inline};
  const preparation={};
  const allowedLoci=lociForScope(species,config.references.scope);
  if(!allowedLoci.length)throw new Error(`The bundled reference pack has no ${config.references.scope} loci for ${config.references.species}.`);
  for(const segment of ["V","D","J","C"]){
    const path=config.references.files?.[segment];if(!path)continue;
    const raw=decodeText(await readFile(resolveFrom(base,path)));
    if(config.references.prepareMetadata){
      process.stderr.write(`Preparing custom ${segment} reference metadata from ${basename(path)}…\n`);
      const report=preprocessGermlineFastaAcrossTiers(raw,segment,germlineTemplateTiers(packed,species,config.references.scope,segment),allowedLoci);
      overrides[segment]=report.fasta;
      const {fasta,...summary}=report;preparation[segment]={file:path,...summary};
      const detail=segment==="V"?`${report.annotated.toLocaleString()} with validated FWR/CDR metadata`:segment==="J"?`${report.annotated.toLocaleString()} with validated frame/CDR3-anchor metadata`:"validated and normalized";
      process.stderr.write(`Prepared ${report.count.toLocaleString()} ${segment} record${report.count===1?"":"s"}; ${detail}.\n`);
      for(const warning of report.warnings)process.stderr.write(`Reference warning: ${warning}\n`);
    }else overrides[segment]=raw;
  }
  const references=compileReferences(species,config.references.scope,overrides);
  if(!references.V.trim()||!references.J.trim())throw new Error("The selected reference composition must contain V and J records.");
  return {references,preparation};
}

async function readMaybeCompressedText(path,label){
  if(path==="-")throw new Error(`${label} must be a FASTA/data file; standard input is reserved for -query.`);
  let bytes;
  try{bytes=await readFile(resolve(path));}
  catch(error){
    if(error?.code==="ENOENT")throw new Error(`${label} was not found at ${path}. The IgBLAST-style germline options require source FASTA files, not makeblastdb prefixes.`);
    throw error;
  }
  if(bytes[0]===0x1f&&bytes[1]===0x8b)bytes=gunzipSync(bytes);
  return decodeText(bytes);
}

async function readVdjFasta(path,label){
  const text=await readMaybeCompressedText(path,label);
  if(!text.trimStart().startsWith(">"))throw new Error(`${label} must point to a nucleotide FASTA file. Binary BLAST database files/prefixes are not decoded by swig-cli.`);
  return text;
}

function vdjScope(options){
  const value=String(options["-ig_seqtype"]??"Ig").toUpperCase();
  if(value==="IG")return "BCR";
  if(value==="TCR")return "TCR";
  throw new Error("-ig_seqtype must be Ig or TCR.");
}

function vdjSpecies(pack,requested){
  const aliases={
    human:"Homo sapiens",
    mouse:"Mus musculus_C57BL/6",
    rat:"Rattus norvegicus_BN; Sprague-Dawley",
    rabbit:"Oryctolagus cuniculus",
    rhesus_monkey:"Macaca mulatta_AG07107",
  };
  const value=String(requested??"human");
  const target=aliases[value.toLowerCase()]??value;
  const species=pack.species.find((entry)=>entry.name.toLowerCase()===target.toLowerCase());
  if(!species)throw new Error(`The embedded Swig reference pack has no species matching -organism ${value}. Use one of human, mouse, rat, rabbit, rhesus_monkey, or an exact pack species name.`);
  return species;
}

function fastaLengths(fasta){
  const result=new Map();let name="",sequence="";
  const finish=()=>{if(name)result.set(name,sequence.replace(/\s/g,"").length);};
  for(const line of fasta.split(/\r?\n/)){
    if(line.startsWith(">")){finish();name=line.slice(1).trim().split(/\s+/,1)[0];sequence="";}
    else if(name)sequence+=line.trim();
  }
  finish();return result;
}

function unmatchedMessage(kind,application){
  const preview=application.unmatched.slice(0,8).join(", ");
  return `${kind} did not match ${application.unmatched.length.toLocaleString()} selected germline identifier${application.unmatched.length===1?"":"s"}${preview?`: ${preview}${application.unmatched.length>8?", …":""}`:""}.`;
}

function sha256(value){return createHash("sha256").update(value).digest("hex");}

function diagnosticTsv(diagnostics){
  const header=["segment","name","locus","status","annotation_source","taxonomic_tier","match_kind","template","identity","attempted_candidates","best_candidate","best_identity","rejections"];
  const rows=diagnostics.map((item)=>[
    item.segment,item.name,item.locus,item.status,item.source,item.taxonomicTier??"",item.matchKind??"",item.template??"",
    item.identity===undefined?"":item.identity.toFixed(6),item.attemptedCandidates,item.bestCandidate??"",
    item.bestIdentity===undefined?"":item.bestIdentity.toFixed(6),
    Object.entries(item.rejectionCounts??{}).sort(([left],[right])=>left.localeCompare(right)).map(([reason,count])=>`${reason}:${count}`).join(";"),
  ].map(cleanCell).join("\t"));
  return `${header.join("\t")}\n${rows.join("\n")}${rows.length?"\n":""}`;
}

function referencePreparationSummary(report){
  const statusCounts={};
  for(const diagnostic of report.diagnostics??[])statusCounts[diagnostic.status]=(statusCounts[diagnostic.status]??0)+1;
  const {fasta,diagnostics,...summary}=report;
  return {...summary,statusCounts};
}

function logReferenceFailures(segment,report){
  for(const diagnostic of report.diagnostics??[]){
    if(diagnostic.status!=="unresolved")continue;
    const best=diagnostic.bestCandidate
      ? `; best candidate ${diagnostic.bestCandidate}${diagnostic.bestIdentity===undefined?"":` at ${(diagnostic.bestIdentity*100).toFixed(1)}% identity`}`
      : "; no annotated candidate";
    const failures=Object.entries(diagnostic.rejectionCounts??{}).sort((left,right)=>right[1]-left[1]||left[0].localeCompare(right[0])).map(([reason,count])=>`${reason}=${count}`).join(", ");
    process.stderr.write(`[prepare:${segment}] unresolved ${diagnostic.name}${best}${failures?`; rejected: ${failures}`:""}.\n`);
  }
}

async function runPrepareReference(rawArgs,assets){
  const options=parseVdjArguments(rawArgs,"prepare");
  if(options.help){process.stdout.write(`${prepareReferenceUsage()}\n`);return;}
  if(options.version){process.stdout.write(`${VERSION}\n`);return;}
  const outputValue=options["-out"];
  if(!outputValue||outputValue==="-")throw new Error("prepare-reference requires --out-prefix with a filesystem path.");
  const paths={V:options["-germline_db_V"],D:options["-germline_db_D"],J:options["-germline_db_J"],C:options["-c_region_db"]};
  if(!paths.V||!paths.J)throw new Error("prepare-reference requires -germline_db_V and -germline_db_J source FASTA files.");
  const packBytes=readFileSync(assets.referencePackPath);
  const pack=JSON.parse(gunzipSync(packBytes).toString("utf8"));
  const species=vdjSpecies(pack,options["-organism"]);
  const scope=vdjScope(options);
  const allowedLoci=lociForScope(species,scope);
  if(!allowedLoci.length)throw new Error(`The embedded reference pack has no ${scope} loci for ${species.name}.`);
  const match=germlineMatchOptions(options,{diagnostics:true});
  const outputPrefix=resolve(String(outputValue));
  await mkdir(dirname(outputPrefix),{recursive:true});
  process.stderr.write(`[prepare] ${species.name} ${scope}; loci ${allowedLoci.join(",")}; ${match.mode.replaceAll("_","-")} matching; ${match.nearestCandidates} nearest candidates.\n`);
  process.stderr.write(`[prepare] V identity floors same-gene=${match.vSameGeneMinIdentity.toFixed(3)}, nearest=${match.vNearestMinIdentity.toFixed(3)}; J same-gene=${match.jSameGeneMinIdentity.toFixed(3)}, nearest=${match.jNearestMinIdentity.toFixed(3)}.\n`);
  if(match.mode==="best_guess")process.stderr.write("[prepare] Best-guess mode disables identity rejection; coordinate, frame, and conserved-anchor validation remain mandatory.\n");

  const exactData={
    V:options["-custom_internal_data"]?await readMaybeCompressedText(options["-custom_internal_data"],"-custom_internal_data"):null,
    J:options["-auxiliary_data"]?await readMaybeCompressedText(options["-auxiliary_data"],"-auxiliary_data"):null,
    D:options["-d_frame_data"]?await readMaybeCompressedText(options["-d_frame_data"],"-d_frame_data"):null,
  };
  if(exactData.D&&!paths.D)throw new Error("-d_frame_data requires -germline_db_D.");
  const references={V:"",D:"",J:"",C:""};
  const files={};
  const segments={};
  const allDiagnostics=[];
  let fwr4EndOffsets={};
  for(const segment of ["V","D","J","C"]){
    const path=paths[segment];
    if(!path)continue;
    const raw=await readVdjFasta(path,segment==="C"?"-c_region_db":`-germline_db_${segment}`);
    let input=raw;
    process.stderr.write(`[prepare:${segment}] Reading and validating ${basename(path)}.\n`);
    if(segment==="V"&&exactData.V){
      const application=applyIgblastInternalData(input,exactData.V);
      input=application.fasta;
      process.stderr.write(`[prepare:V] IgBLAST internal data matched ${application.matched.toLocaleString()}/${application.total.toLocaleString()} records; unmatched records continue to homology transfer.\n`);
    }
    if(segment==="J"&&exactData.J){
      const application=applyIgblastAuxiliaryData(input,exactData.J);
      input=application.fasta;fwr4EndOffsets=application.fwr4EndOffsets??{};
      process.stderr.write(`[prepare:J] IgBLAST auxiliary data matched ${application.matched.toLocaleString()}/${application.total.toLocaleString()} records and annotated ${application.annotated.toLocaleString()} CDR3 stops; unmatched records continue to homology/motif inference.\n`);
    }
    if(segment==="D"&&exactData.D){
      const application=applyIgblastDFrameData(input,exactData.D);
      input=application.fasta;
      process.stderr.write(`[prepare:D] IgBLAST D-frame data matched ${application.matched.toLocaleString()}/${application.total.toLocaleString()} records.\n`);
    }
    const started=performance.now();
    const report=preprocessGermlineFastaAcrossTiers(input,segment,germlineTemplateTiers(pack,species,scope,segment),allowedLoci,match);
    const seconds=(performance.now()-started)/1000;
    references[segment]=report.fasta;
    allDiagnostics.push(...(report.diagnostics??[]));
    logReferenceFailures(segment,report);
    const outputPath=`${outputPrefix}.${segment}.fasta`;
    files[segment]={path:basename(outputPath),sha256:sha256(report.fasta),records:report.count,annotated:report.annotated};
    segments[segment]={source:basename(path),sourceSha256:sha256(raw),...referencePreparationSummary(report),seconds:Number(seconds.toFixed(3))};
    process.stderr.write(`[prepare:${segment}] ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} annotated; ${report.unannotated.toLocaleString()} unresolved; ${seconds.toFixed(2)} s.\n`);
  }

  const diagnosticsText=diagnosticTsv(allDiagnostics);
  const diagnosticsPath=`${outputPrefix}.annotation-diagnostics.tsv`;
  const manifestPath=`${outputPrefix}.swig-reference.json`;
  const incomplete=(segments.V?.unannotated??0)+(segments.J?.unannotated??0);
  const manifest={
    schema:1,
    application:"Swig prepared germline reference",
    applicationVersion:VERSION,
    createdAt:new Date().toISOString(),
    complete:incomplete===0,
    species:species.name,
    scope,
    loci:allowedLoci,
    referencePack:{source:pack.source,release:pack.release,retrieved:pack.retrieved,sha256:sha256(packBytes)},
    match:{...match,includeDiagnostics:undefined},
    files,
    diagnostics:{path:basename(diagnosticsPath),sha256:sha256(diagnosticsText),records:allDiagnostics.length},
    segments,
    fwr4EndOffsets,
  };
  for(const segment of Object.keys(files))await writeFile(join(dirname(outputPrefix),files[segment].path),references[segment]);
  await writeFile(diagnosticsPath,diagnosticsText);
  await writeFile(manifestPath,`${JSON.stringify(manifest,null,2)}\n`);
  process.stderr.write(`[prepare] Wrote ${Object.keys(files).length} prepared FASTA${Object.keys(files).length===1?"":"s"}, ${basename(diagnosticsPath)}, and ${basename(manifestPath)}.\n`);
  process.stdout.write(`${JSON.stringify({manifest:manifestPath,diagnostics:diagnosticsPath,complete:manifest.complete,segments:Object.fromEntries(Object.entries(segments).map(([segment,value])=>[segment,{count:value.count,annotated:value.annotated,unannotated:value.unannotated,seconds:value.seconds}]))},null,2)}\n`);
  if(options["--require-complete"]&&incomplete)throw new Error(`${incomplete.toLocaleString()} V/J germline record${incomplete===1?" remains":"s remain"} unresolved; inspect ${diagnosticsPath}.`);
}

async function loadPreparedReference(path){
  const manifestPath=resolve(path);
  let manifest;
  try{manifest=JSON.parse(await readFile(manifestPath,"utf8"));}
  catch(error){throw new Error(`Could not read prepared-reference manifest ${path}: ${error instanceof Error?error.message:String(error)}`);}
  if(manifest?.schema!==1||manifest?.application!=="Swig prepared germline reference")throw new Error(`${path} is not a supported Swig prepared-reference manifest.`);
  const base=dirname(manifestPath);const references={V:"",D:"",J:"",C:""};
  for(const segment of ["V","D","J","C"]){
    const entry=manifest.files?.[segment];if(!entry)continue;
    const relative=typeof entry==="string"?entry:entry.path;
    if(typeof relative!=="string"||!relative)throw new Error(`Prepared-reference manifest has no valid ${segment} FASTA path.`);
    const fasta=await readVdjFasta(resolve(base,relative),`prepared ${segment} reference`);
    const expected=typeof entry==="object"?entry.sha256:undefined;
    if(expected&&sha256(fasta)!==expected)throw new Error(`Prepared ${segment} FASTA failed its SHA-256 check: ${resolve(base,relative)}.`);
    references[segment]=fasta;
    process.stderr.write(`[prepared-reference:${segment}] Loaded and verified ${basename(relative)}.\n`);
  }
  if(!references.V||!references.J)throw new Error("Prepared-reference manifest must provide V and J FASTAs.");
  if(manifest.complete===false)process.stderr.write("Reference warning: this prepared-reference manifest contains unresolved V/J metadata; assignments remain available but some region annotations may be blank.\n");
  return {references,mode:"prepared-reference",fwr4EndOffsets:manifest.fwr4EndOffsets??{},jLengths:fastaLengths(references.J)};
}

async function prepareVdjReferences(options,assets){
  const preparedPath=options["--prepared-reference"];
  if(preparedPath){
    const conflicts=["-germline_db_V","-germline_db_D","-germline_db_J","-c_region_db","-custom_internal_data","-auxiliary_data","-d_frame_data","--swigannots","--match-mode","--best-guess","--nearest-candidates","--v-same-gene-min-identity","--v-nearest-min-identity","--j-same-gene-min-identity","--j-nearest-min-identity"].filter((name)=>options[name]);
    if(conflicts.length)throw new Error(`--prepared-reference cannot be combined with ${conflicts.join(", ")}.`);
    return loadPreparedReference(preparedPath);
  }
  const paths={V:options["-germline_db_V"],D:options["-germline_db_D"],J:options["-germline_db_J"],C:options["-c_region_db"]};
  if(!paths.V||!paths.J)throw new Error("--vdj requires -germline_db_V and -germline_db_J source FASTA files.");
  const raw={
    V:await readVdjFasta(paths.V,"-germline_db_V"),
    D:paths.D?await readVdjFasta(paths.D,"-germline_db_D"):"",
    J:await readVdjFasta(paths.J,"-germline_db_J"),
    C:paths.C?await readVdjFasta(paths.C,"-c_region_db"):"",
  };
  const useSwig=Boolean(options["--swigannots"]);
  const internalPath=options["-custom_internal_data"];
  const auxiliaryPath=options["-auxiliary_data"];
  if(useSwig&&(internalPath||auxiliaryPath))throw new Error("--swigannots cannot be combined with -custom_internal_data or -auxiliary_data; choose one annotation source.");
  const scope=vdjScope(options);
  const match=germlineMatchOptions(options);
  let pack;let species;let allowedLoci;
  const needPack=useSwig||Boolean(auxiliaryPath&&!internalPath);
  if(needPack){
    pack=JSON.parse(gunzipSync(readFileSync(assets.referencePackPath)).toString("utf8"));
    species=vdjSpecies(pack,options["-organism"]);
    allowedLoci=lociForScope(species,scope);
    if(!allowedLoci.length)throw new Error(`The embedded reference pack has no ${scope} loci for ${species.name}.`);
  }

  const references={V:"",D:"",J:"",C:""};
  let mode="assignments-only";
  if(useSwig){
    mode="swig-metadata";
    for(const segment of ["V","D","J","C"]){
      if(!raw[segment])continue;
      const report=preprocessGermlineFastaAcrossTiers(raw[segment],segment,germlineTemplateTiers(pack,species,scope,segment),allowedLoci,match);
      references[segment]=report.fasta;
      process.stderr.write(`Swig metadata preparation: ${segment} ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} annotated.\n`);
      for(const warning of report.warnings)process.stderr.write(`Reference warning: ${warning}\n`);
    }
  }else{
    for(const segment of ["V","D","J","C"]){
      if(raw[segment])references[segment]=prepareIgblastStyleGermlineFasta(raw[segment],segment).fasta;
    }
  }

  if(internalPath){
    mode="igblast-data";
    const application=applyIgblastInternalData(references.V,await readMaybeCompressedText(internalPath,"-custom_internal_data"));
    if(application.matched!==application.total)throw new Error(`${unmatchedMessage("-custom_internal_data",application)} IgBLAST requires custom internal data for every selected V sequence.`);
    references.V=application.fasta;
    process.stderr.write(`IgBLAST internal data: ${application.matched.toLocaleString()}/${application.total.toLocaleString()} V records annotated.\n`);
  }else if(auxiliaryPath){
    // IgBLAST normally obtains V domains from its organism-specific internal
    // directory when only -auxiliary_data is supplied. The standalone Swig
    // equivalent uses the fixed embedded pack and the web metadata transfer.
    mode="igblast-data";
    const report=preprocessGermlineFastaAcrossTiers(raw.V,"V",germlineTemplateTiers(pack,species,scope,"V"),allowedLoci,match);
    references.V=report.fasta;
    process.stderr.write(`Embedded ${species.name} V metadata: ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} records annotated.\n`);
    for(const warning of report.warnings)process.stderr.write(`Reference warning: ${warning}\n`);
  }

  let fwr4EndOffsets={};
  if(auxiliaryPath){
    mode="igblast-data";
    const application=applyIgblastAuxiliaryData(references.J,await readMaybeCompressedText(auxiliaryPath,"-auxiliary_data"));
    if(!application.matched)throw new Error("-auxiliary_data has no exact identifiers matching the selected J FASTA.");
    references.J=application.fasta;fwr4EndOffsets=application.fwr4EndOffsets??{};
    process.stderr.write(`IgBLAST auxiliary data: ${application.annotated.toLocaleString()}/${application.total.toLocaleString()} J records have CDR3-stop annotations.\n`);
    if(application.unmatched.length)process.stderr.write(`Reference warning: ${unmatchedMessage("-auxiliary_data",application)}\n`);
  }

  const dFramePath=options["-d_frame_data"];
  if(dFramePath){
    if(!references.D)throw new Error("-d_frame_data requires -germline_db_D.");
    const application=applyIgblastDFrameData(references.D,await readMaybeCompressedText(dFramePath,"-d_frame_data"));
    if(!application.matched)throw new Error("-d_frame_data has no exact identifiers matching the selected D FASTA.");
    references.D=application.fasta;
    process.stderr.write(`IgBLAST D-frame data: ${application.annotated.toLocaleString()}/${application.total.toLocaleString()} D records annotated.\n`);
    if(application.unmatched.length)process.stderr.write(`Reference warning: ${unmatchedMessage("-d_frame_data",application)}\n`);
  }
  return {references,mode,fwr4EndOffsets,jLengths:fastaLengths(references.J)};
}

function compactAlleleResult(result){
  return {version:result.version,options:result.options,totalRecords:result.totalRecords,activeRecords:result.activeRecords,runAt:result.runAt,warnings:result.warnings,segments:Object.fromEntries(Object.entries(result.segments).map(([segment,value])=>[segment,{segment:value.segment,models:value.models,modeledRows:value.modeledRows,changedMapRows:value.changedMapRows,skippedRows:value.skippedRows,matrixNonZeros:value.matrixNonZeros,truncatedRows:value.truncatedRows,exactDuplicateLabels:value.exactDuplicateLabels}]))};
}

function fitAlleles(rows,references,options){
  const segments={};
  for(const segment of options.segments){
    const fasta=references[segment];if(!fasta?.trim())continue;
    const graph=buildReferenceAlleleGraph(fasta,segment,options.neighbourRadius);
    const accumulator=new SparseEvidenceAccumulator(rows.length,graph,options);
    rows.forEach((values,ordinal)=>accumulator.add(toRefinementInputRow({ordinal,values},segment)));
    segments[segment]=fitSparseAlleleModel(accumulator.finish(),graph,options,rows.length);
  }
  const warnings=[];
  if(options.model==="active-set")warnings.push("The hurdle model estimates repertoire-active usage, not literal genomic presence.");
  if(options.weighting==="abundance")warnings.push("Abundance weighting allows clonal expansion to influence the fitted mixture.");
  return {version:1,options,totalRecords:rows.length,activeRecords:rows.length,segments,runAt:new Date().toISOString(),warnings};
}

function applyAlleles(rows,result,policy,threshold,headers){
  for(const segment of ["V","D","J"]){
    const prefix=segment.toLowerCase();
    if(!result.segments[segment])continue;
    addHeaders(headers,[`swig_original_${prefix}_call`,`swig_repertoire_${prefix}_call`]);
    rows.forEach((row,ordinal)=>{const call=refinedCall(result,segment,ordinal,policy,threshold);if(!call)return;row[`swig_original_${prefix}_call`]=row[`${prefix}_call`]??"";row[`swig_repertoire_${prefix}_call`]=call;row[`${prefix}_call`]=call;});
  }
}

async function runChimera(rows,activeMask,config,scope,headers,references,workers){
  const selectedSegment=config.segment.toUpperCase();
  const msaText=config.uploadedMsa?.trim()||(config.msaSource==="selected"?references[selectedSegment]?.trim():"");
  if(!msaText)throw new Error("CLI chimera filtering requires either pipeline.chimera.uploadedMsa or an aligned selected-segment reference.");
  try{prepareReferenceMsa(msaText);}
  catch(error){
    if(config.msaSource==="selected"&&!config.uploadedMsa?.trim())throw new Error(`The selected ${selectedSegment} references are not already aligned. Export this pipeline from Swig Web to embed its Kalign MSA, or set pipeline.chimera.uploadedMsa explicitly. ${error instanceof Error?error.message:String(error)}`);
    throw error;
  }
  const segment=config.segment.toLowerCase();
  const method=config.model==="auto"?(scope==="TCR"||String(scope).startsWith("TR")?"DB":"BW"):config.model;
  const options={method,priorProbability:config.priorProbability,baseMutationProbability:config.baseMutationProbability,mutationRates:config.mutationRates,mutationSwitchProbability:config.mutationSwitchProbability,detailed:config.detailed};
  let evaluated=0,flagged=0,unevaluated=0;
  addHeaders(headers,["swig_chimera_probability","swig_chimera_status"]);
  if(config.detailed)addHeaders(headers,["swig_chimera_starting_reference","swig_chimera_recombinations"]);
  const batches=[];let batch=[];
  for(let ordinal=0;ordinal<rows.length;ordinal+=1){
    if(!activeMask[ordinal])continue;
    const row=rows[ordinal],call=row[`${segment}_call`]??"",sequence=row[`${segment}_sequence_alignment`]??"",germline=row[`${segment}_germline_alignment`]??"";
    batch.push({ordinal,call,sequenceAlignment:sequence,germlineAlignment:germline});if(batch.length>=250){batches.push(batch);batch=[];}
  }
  if(batch.length)batches.push(batch);
  if(!batches.length)return {evaluated,flagged,unevaluated,threshold:config.posteriorThreshold};
  const workerCount=Math.max(1,Math.min(Math.floor(workers)||1,batches.length));
  process.stderr.write(`CHMMAIRRa: ${workerCount} worker${workerCount===1?"":"s"} across ${batches.length.toLocaleString()} row batches.\n`);
  const pool=new ChmmPool(workerCount,{msa:msaText,options,minDfr:config.minimumDfr});
  try{
    await pool.start();const outputs=await Promise.all(batches.map((rows)=>pool.run(rows)));
    for(const output of outputs)for(const result of output.results){
      const row=rows[result.ordinal];row.swig_chimera_status=result.status;
      if(result.status==="evaluated"){
        row.swig_chimera_probability=String(result.probability);evaluated+=1;
        if(config.detailed){row.swig_chimera_starting_reference=result.startingReference;row.swig_chimera_recombinations=result.recombinations.map((event)=>`${event.left}->${event.right}@${event.position}`).join(";");}
        if(result.probability>=config.posteriorThreshold){activeMask[result.ordinal]=0;flagged+=1;}
      }else{unevaluated+=1;if(!config.retainUnevaluated)activeMask[result.ordinal]=0;}
    }
  }finally{await pool.close();}
  return {evaluated,flagged,unevaluated,threshold:config.posteriorThreshold};
}

function writeChunk(stream,chunk){if(stream.errored)return Promise.reject(stream.errored);if(stream.write(chunk))return Promise.resolve();return once(stream,"drain");}

function createCliWriteStream(path){return createWriteStream(path,{highWaterMark:CLI_STREAM_HIGH_WATER_MARK});}

function createVdjOutput(path){
  const file=createCliWriteStream(path);
  // A failed run abandons its output: late write errors must not be reported over the real one.
  const discard=(stream)=>{stream.on("error",()=>{});stream.destroy();};
  if(!/\.gz$/i.test(path)){
    // A write error (a full disk, say) is reported by the next write or by finish(), not as an
    // unhandled 'error' event.
    file.on("error",()=>{});
    return {stream:file,opened:once(file,"open"),finish:()=>finishWritable(file),destroy:()=>discard(file)};
  }
  const gzip=createGzip();
  const completed=pipeline(gzip,file);
  return {
    stream:gzip,
    opened:once(file,"open"),
    finish:async()=>{gzip.end();await completed;},
    destroy:()=>{discard(gzip);discard(file);void completed.catch(()=>{});},
  };
}

function airrOutputName(prefix,kind,compression){
  return `${prefix}.${kind}.airr.tsv${compression==="gzip"?".gz":""}`;
}

async function finishWritable(stream){if(stream.errored)throw stream.errored;const completed=once(stream,"finish");stream.end();await completed;}

async function writeRowsFile(path,headers,rows,include=()=>true){
  const output=createVdjOutput(path);const stream=output.stream;
  try{
    await output.opened;
    await writeChunk(stream,`${headers.join("\t")}\n`);
    let batch=[];
    for(let ordinal=0;ordinal<rows.length;ordinal+=1){
      if(!include(rows[ordinal],ordinal))continue;
      batch.push(rows[ordinal]);
      if(batch.length>=2_000){await writeChunk(stream,serializeRowBody(headers,batch));batch=[];}
    }
    if(batch.length)await writeChunk(stream,serializeRowBody(headers,batch));
    await output.finish();
  }catch(error){output.destroy();throw error;}
}

async function writeLineageStudy(path,manifestPath,headers,rows,activeMask,lineages,config,references,shm){
  const stream=createCliWriteStream(path);const hash=createHash("sha256");let offset=0;
  const append=async(text)=>{const chunk=Buffer.from(text,"utf8");hash.update(chunk);offset+=chunk.byteLength;await writeChunk(stream,chunk);};
  await append(`${headers.join("\t")}\n`);
  const retainedIds=new Set(lineages.summaries.map((summary)=>summary.id));
  const indexed=rows.map((row,ordinal)=>({row,ordinal,lineageId:lineages.assignments[ordinal]??0})).filter((item)=>activeMask[item.ordinal]&&retainedIds.has(item.lineageId)).sort((left,right)=>left.lineageId-right.lineageId||left.ordinal-right.ordinal);
  const ranges=[];let current=0,start=offset,count=0,rangeHash=createHash("sha256");
  const finishRange=()=>{if(current>0)ranges.push({lineageId:current,start,end:offset,records:count,sha256:rangeHash.digest("hex")});};
  for(const item of indexed){
    if(item.lineageId!==current){finishRange();current=item.lineageId;start=offset;count=0;rangeHash=createHash("sha256");}
    const text=`${headers.map((header)=>cleanCell(item.row[header])).join("\t")}\n`;
    rangeHash.update(text,"utf8");await append(text);count+=1;
  }
  finishRange();stream.end();await once(stream,"finish");
  const lowestShmByLineage=new Map((shm?.lowestByLineage??[]).map((value)=>[value.lineageId,value]));
  const shmSummaries=shm?.lineages.flatMap((group)=>{const match=/^Lineage\s+(\d+)$/.exec(group.label);if(!match)return[];const lineageId=Number(match[1]),lowest=lowestShmByLineage.get(lineageId);return[{lineageId,mean:group.mean,p95:group.p95??0,ordinal:lowest?.ordinal,cdr3Nt:lowest?.cdr3Nt,cdr3Aa:lowest?.cdr3Aa}];})??[];
  const manifest={schema:1,application:"Swig lineage study",applicationVersion:VERSION,createdAt:new Date().toISOString(),linkedAirr:{name:basename(path),size:offset,records:indexed.length,headers,sha256:hash.digest("hex")},analysis:{inputName:config.studyName,species:config.references.species,scope:config.references.scope,references,datasets:config.inputs.map((input)=>({datasetId:input.datasetId,inputName:input.inputName||basename(input.path),sampleId:input.sampleId,subjectId:input.subjectId,cohort:input.cohort,timepoint:input.timepoint,compartment:input.compartment})),callingProfile:config.annotation.callingProfile,minimumIdentity:config.annotation.minimumIdentity,strand:config.annotation.strand,lineage:{scope:config.pipeline.lineage.scope,identity:config.pipeline.lineage.identity,resolution:config.pipeline.lineage.resolution,ambiguity:config.pipeline.lineage.ambiguity,productiveOnly:config.pipeline.lineage.productiveOnly,maxCandidateComparisons:config.pipeline.lineage.maxCandidateComparisons}},summaries:lineages.summaries,ranges,shm:shm?{metric:shm.metric,summaries:shmSummaries}:undefined};
  await writeFile(manifestPath,gzipSync(JSON.stringify(manifest)));
  return manifest;
}

async function runPipeline(config,base,assets){
  if(!config.inputs.length)throw new Error("No input datasets were specified.");
  const outputDirectory=resolveFrom(base,config.output.directory);await mkdir(outputDirectory,{recursive:true});
  const prefix=config.output.prefix;
  const loadedReferences=await loadReferences(config,base,assets.referencePackPath);
  const references=loadedReferences.references;
  if(config.annotation.airrMode==="preserve"&&config.annotation.doubleD.mode!=="off"&&config.inputs.some((input)=>detectFormat(input.path,input.format)==="airr"))throw new Error("Double-D screening of AIRR input requires annotation.airrMode = \"reannotate\".");
  const needsAnnotation=config.inputs.some((input)=>detectFormat(input.path,input.format)!=="airr"||config.annotation.airrMode==="reannotate");
  const annotationBatchRecords=config.annotation.batchRecords||automaticBatchRecords(config.annotation.workers);
  const settings={callingProfile:config.annotation.callingProfile,minimumIdentity:config.annotation.minimumIdentity,minimumConstantPrefixIdentity:config.annotation.minimumConstantPrefixIdentity,strand:config.annotation.strand,doubleD:config.annotation.doubleD};
  const init=workerInitialization(assets.wasmPath,references,settings);
  const mode=config.annotation.aligner;
  let nativeBinary=needsAnnotation?chooseNativeAligner(mode,assets,config.annotation.doubleD.mode!=="off"?"the double-D screen":"",config.annotation.alignerLabel):null;
  let pool=null;
  const wasmPool=async()=>{
    if(!pool){
      process.stderr.write(`Starting the WebAssembly aligner pool (${config.annotation.workers} worker${config.annotation.workers===1?"":"s"}; ${annotationBatchRecords.toLocaleString()} records/batch).\n`);
      pool=new WasmPool(config.annotation.workers,init);
      try{await pool.start();}
      catch(error){const failed=pool;pool=null;await failed.close().catch(()=>{});throw error;}
    }
    return pool;
  };
  if(nativeBinary)process.stderr.write(`Using the native aligner (${config.annotation.workers} thread${config.annotation.workers===1?"":"s"}; ${annotationBatchRecords.toLocaleString()} records/batch).\n`);
  else if(needsAnnotation)await wasmPool();
  const annotatedName=airrOutputName(prefix,"annotated",config.output.airrCompression);
  const annotatedOutput=config.output.writeAnnotatedAirr?createVdjOutput(join(outputDirectory,annotatedName)):null;
  const annotatedStream=annotatedOutput?.stream??null;
  let annotatedOutputHeaders=null;
  const rows=[];const headers=[];let annotatedRecords=0;let inputRecords=0;let eligibleRecords=0;
  let fastqFilterStats=emptyFastqQualityFilterStats(config.preprocessing.fastqFilter.enabled,false);
  try{
    if(annotatedOutput)await annotatedOutput.opened;
    for(let datasetIndex=0;datasetIndex<config.inputs.length;datasetIndex+=1){
      const input=config.inputs[datasetIndex];
      const format=detectFormat(input.path,input.format);
      const preserveAirr=format==="airr"&&config.annotation.airrMode==="preserve";
      process.stderr.write(`${preserveAirr?"Loading existing AIRR calls from":"Annotating"} ${input.inputName||basename(input.path)} (${input.sampleId}; donor ${input.subjectId})…\n`);
      const preprocessingState={inputRecords:0,eligibleRecords:0,selectedRecords:0,fastqFilter:emptyFastqQualityFilterStats(false,false)};
      const pending=[];
      let datasetPool=null;
      if(needsAnnotation&&!preserveAirr){
        if(nativeBinary){
          // One child per dataset: the engine detects FASTA, FASTQ or AIRR once per stream.
          const native=new NativePool(nativeBinary,config.annotation.workers,init,annotationBatchRecords);
          try{await native.start();datasetPool=native;}
          catch(error){
            await native.close().catch(()=>{});
            if(mode==="native"||/database/i.test(error.message))throw error;
            process.stderr.write(`The native aligner did not start (${error.message}); using WebAssembly.\n`);
            nativeBinary=null;
          }
        }
        datasetPool??=await wasmPool();
      }
      const window=datasetPool instanceof NativePool?datasetPool.window:Math.max(2,config.annotation.workers*2);
      const consume=async(item)=>{
        const result=await item.promise;
        const manifest={datasetId:input.datasetId,inputName:input.inputName||basename(input.path),sampleId:input.sampleId,subjectId:input.subjectId,cohort:input.cohort,timepoint:input.timepoint,compartment:input.compartment};
        const annotated=annotateAirrBatch(result.header,result.body,manifest);
        const table=parseTable(annotated.header,annotated.body);
        const doubleDAnnotated=result.doubleDHeader?annotateDoubleDBatch(result.doubleDHeader,result.doubleDBody,manifest):null;
        const doubleD=doubleDAnnotated?parseTable(doubleDAnnotated.header,doubleDAnnotated.body):null;
        const ddById=new Map((doubleD?.rows??[]).map((row)=>[row.sequence_id,row]));
        const batchHeaders=[...table.headers];
        if(doubleD)addHeaders(batchHeaders,doubleD.headers);
        for(const row of table.rows){
          const dd=ddById.get(row.sequence_id);if(dd)Object.assign(row,dd);
          rows.push(row);
        }
        addHeaders(headers,batchHeaders);
        if(annotatedStream){
          if(!annotatedOutputHeaders){annotatedOutputHeaders=batchHeaders;await writeChunk(annotatedStream,`${annotatedOutputHeaders.join("\t")}\n`);}
          else{
            const newHeaders=batchHeaders.filter((header)=>!annotatedOutputHeaders.includes(header));
            if(newHeaders.length)throw new Error(`Streaming annotated AIRR output cannot add columns after its header was written (${newHeaders.join(", ")}). Use matching AIRR schemas for all preserved inputs or reannotate them.`);
          }
          await writeChunk(annotatedStream,serializeRowBody(annotatedOutputHeaders,table.rows));
        }
        annotatedRecords+=table.rows.length;
      };
      try{
        for await(const batch of sequenceBatches(input,annotationBatchRecords,config.preprocessing,config.annotation.airrMode,datasetIndex,preprocessingState)){
          // Await the oldest batch only once another is ready to queue: then every batch held is
          // full, which the native window counts on (a short last batch waits for finishInput).
          if(pending.length>=window)await consume(pending.shift());
          const promise=batch.format==="airr"&&config.annotation.airrMode==="preserve"
            ? Promise.resolve({direct:true,header:batch.header,body:batch.body,doubleDHeader:"",doubleDBody:""})
            : datasetPool.run(workerAnnotation(batch.text,batch.count));
          pending.push({promise,count:batch.count});
          if(datasetPool instanceof NativePool)await datasetPool.ready();
        }
        if(datasetPool instanceof NativePool)datasetPool.finishInput();
        while(pending.length)await consume(pending.shift());
        if(datasetPool instanceof NativePool)await datasetPool.finish();
      }finally{
        if(datasetPool instanceof NativePool)await datasetPool.close();
      }
      inputRecords+=preprocessingState.inputRecords;eligibleRecords+=preprocessingState.eligibleRecords;
      fastqFilterStats=addFastqQualityFilterStats(fastqFilterStats,preprocessingState.fastqFilter);
    }
  }catch(error){annotatedOutput?.destroy();throw error;}
  finally{await pool?.close();}
  if(annotatedOutput)await annotatedOutput.finish();
  rows.forEach((row,ordinal)=>{if(!row.sequence_id)row.sequence_id=`swig_${ordinal+1}`;});

  let alleleResult=null;
  if(config.pipeline.alleleRefinement.enabled){
    process.stderr.write("Fitting repertoire-level allele model…\n");
    alleleResult=fitAlleles(rows,references,config.pipeline.alleleRefinement);
    applyAlleles(rows,alleleResult,config.pipeline.alleleRefinement.reassignmentPolicy,config.pipeline.alleleRefinement.applyMinimumPosterior,headers);
    await writeFile(join(outputDirectory,`${prefix}.allele-models.json`),JSON.stringify(compactAlleleResult(alleleResult),null,2));
  }

  const records=rows.map((values,ordinal)=>airrRowToPostAnalysisRecord({ordinal,values}));
  let collapseResult=null;let activeMask=new Uint8Array(rows.length);activeMask.fill(1);
  if(config.pipeline.collapse.enabled){
    process.stderr.write(`${config.pipeline.collapse.mode==="exact"?"Collapsing exact duplicates":"Denoising reads"}…\n`);
    if(config.pipeline.collapse.mode==="exact"){
      const shardCount=records.length>=10_000?config.annotation.workers:1;
      if(shardCount===1)collapseResult=deduplicate(records,config.pipeline.collapse.key,config.pipeline.collapse.unresolvedPolicy,config.pipeline.collapse.scope,config.pipeline.collapse.respectConstantCall);
      else{
        const plan=createExactDedupPlan(records,config.pipeline.collapse.key,config.pipeline.collapse.unresolvedPolicy,config.pipeline.collapse.scope,config.pipeline.collapse.respectConstantCall,shardCount);
        const workerCount=Math.max(1,Math.min(config.annotation.workers,plan.jobs.length));
        process.stderr.write(`Exact collapse: ${workerCount} worker${workerCount===1?"":"s"} across ${plan.jobs.length.toLocaleString()} deterministic key shards.\n`);
        const results=await runPostAnalysisTasks(plan.jobs.map((job)=>({kind:"exact",job})),workerCount);
        collapseResult=finishExactDedupPlan(plan,results);
      }
    }
    else{
      const options={...config.pipeline.collapse.denoise,mode:config.pipeline.collapse.mode,scope:config.pipeline.collapse.scope,respectConstantCall:config.pipeline.collapse.respectConstantCall,unresolvedPolicy:config.pipeline.collapse.unresolvedPolicy};
      const accumulator=new DenoiseAccumulator(records,options);rows.forEach((values,ordinal)=>accumulator.add(ordinal,denoiseVdjSequence({ordinal,values})));
      const jobs=accumulator.preparePartitionJobs();const variants=jobs.reduce((total,job)=>total+job.variants.length,0);
      const requested=variants>=500?config.annotation.workers:1;const workerCount=Math.max(1,Math.min(requested,jobs.length||1));
      process.stderr.write(`${config.pipeline.collapse.mode.toUpperCase()} denoising: ${workerCount} worker${workerCount===1?"":"s"} across ${jobs.length.toLocaleString()} independent V/J partition${jobs.length===1?"":"s"}.\n`);
      const results=await runPostAnalysisTasks(jobs.map((job)=>({kind:"denoise",job})),workerCount);
      collapseResult=accumulator.finishWithPartitionResults(results);
    }
    activeMask=Uint8Array.from(collapseResult.counts,(count)=>count>0?1:0);addHeaders(headers,["duplicate_count"]);
    rows.forEach((row,ordinal)=>{if(collapseResult.counts[ordinal])row.duplicate_count=String(collapseResult.counts[ordinal]);});
  }

  let chimeraSummary=null;
  if(config.pipeline.chimera.enabled){process.stderr.write("Filtering candidate chimeras…\n");chimeraSummary=await runChimera(rows,activeMask,config.pipeline.chimera,config.references.scope,headers,references,config.annotation.workers);}

  if(config.pipeline.selection.enabled){
    const errors=validateRepertoireSelection(config.pipeline.selection);if(errors.length)throw new Error(errors.join(" "));
    process.stderr.write("Applying explicit repertoire selection…\n");
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal]&&!repertoireRowMatches(row,config.pipeline.selection))activeMask[ordinal]=0;});
  }

  let lineages=null;
  if(config.pipeline.lineage.enabled){
    process.stderr.write("Assigning lineages…\n");
    const doubleDMask=Uint8Array.from(rows,(row)=>row.d2_call?1:0);
    lineages=assignLineages(records,{identity:config.pipeline.lineage.identity,callResolution:config.pipeline.lineage.resolution,ambiguity:config.pipeline.lineage.ambiguity,productiveOnly:config.pipeline.lineage.productiveOnly,requireSameLocus:true,maxCandidateComparisons:config.pipeline.lineage.maxCandidateComparisons,scope:config.pipeline.lineage.scope},collapseResult??undefined,activeMask,doubleDMask);
    addHeaders(headers,["clone_id"]);rows.forEach((row,ordinal)=>{row.clone_id=lineages.assignments[ordinal]>0?String(lineages.assignments[ordinal]):"";});
  }

  let shm=null;
  if(config.pipeline.shm.enabled){
    process.stderr.write("Calculating SHM summaries…\n");
    const accumulator=new ShmAccumulator({metric:config.pipeline.shm.metric,maxSamplesPerLineage:2000});
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal])accumulator.add(row,ordinal,lineages?.assignments[ordinal]??0,"All selected");});shm=accumulator.finish();
    await writeFile(join(outputDirectory,`${prefix}.shm.json`),JSON.stringify(shm,null,2));
  }

  let missingAlleles=null;
  if(config.pipeline.missingAlleles.enabled){
    if(!lineages)throw new Error("Missing-allele screening requires lineage assignment.");
    process.stderr.write("Screening for possible missing V alleles…\n");
    const accumulator=new MissingAlleleAccumulator(config.pipeline.missingAlleles);
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal])accumulator.add(row,ordinal,lineages.assignments[ordinal]??0);});
    const validator=accumulator.prepareValidation(references.V);
    rows.forEach((row,ordinal)=>{if(activeMask[ordinal])validator.add(row,ordinal,lineages.assignments[ordinal]??0);});missingAlleles=validator.finish();
    await writeFile(join(outputDirectory,`${prefix}.missing-v-alleles.json`),JSON.stringify(missingAlleles,null,2));
  }

  addHeaders(headers,["swig_retained"]);rows.forEach((row,ordinal)=>{row.swig_retained=activeMask[ordinal]?"T":"F";});
  const processedName=airrOutputName(prefix,"processed",config.output.airrCompression);
  await writeRowsFile(join(outputDirectory,processedName),headers,rows,(_,ordinal)=>Boolean(activeMask[ordinal]));

  let lineageStudy=null;
  if(config.output.writeLineageStudy&&lineages){
    process.stderr.write("Writing lazy lineage-study bundle…\n");
    lineageStudy=await writeLineageStudy(join(outputDirectory,`${prefix}.lineages.airr.tsv`),join(outputDirectory,`${prefix}.swig-lineage-study.json.gz`),headers,rows,activeMask,lineages,config,references,shm);
  }

  const retained=activeMask.reduce((sum,value)=>sum+(value?1:0),0);
  const summary={application:"swig-cli",version:VERSION,completedAt:new Date().toISOString(),inputRecords,eligibleRecords,annotatedRecords,retainedRecords:retained,lineages:lineages?.lineageCount??0,outputs:{annotated:config.output.writeAnnotatedAirr?annotatedName:null,processed:processedName,airrCompression:config.output.airrCompression},references:{metadataPreparation:{enabled:config.references.prepareMetadata,segments:loadedReferences.preparation}},preprocessing:{subsample:config.preprocessing.subsample,fastqFilter:{options:config.preprocessing.fastqFilter,stats:fastqFilterStats}},collapse:collapseResult?{mode:collapseResult.mode,inputRecords:collapseResult.inputRecords,inputAbundance:collapseResult.inputAbundance,uniqueRecords:collapseResult.uniqueRecords,collapsedRecords:collapseResult.collapsedRecords,warnings:collapseResult.warnings}:null,chimera:chimeraSummary,alleleRefinement:alleleResult?compactAlleleResult(alleleResult):null,shm:shm?{analyzedRecords:shm.analyzedRecords,analyzedAbundance:shm.analyzedAbundance,skippedRecords:shm.skippedRecords,metric:shm.metric}:null,missingAlleles:missingAlleles?{candidates:missingAlleles.candidates.length,warnings:missingAlleles.warnings}:null,lineageStudy:lineageStudy?{manifest:`${prefix}.swig-lineage-study.json.gz`,airr:`${prefix}.lineages.airr.tsv`,records:lineageStudy.linkedAirr.records}:null};
  await writeFile(join(outputDirectory,`${prefix}.summary.json`),JSON.stringify(summary,null,2));
  await writeFile(join(outputDirectory,`${prefix}.resolved-config.json`),JSON.stringify(config,null,2));
  process.stdout.write(`${JSON.stringify(summary,null,2)}\n`);
}

function vdjTuning(options){
  const tuning={};
  if(options["-min_D_match"]!==undefined)tuning.minimumDMatch=parseIntegerOption(options["-min_D_match"],"-min_D_match",{minimum:4});
  if(options["-min_J_length"]!==undefined)tuning.minimumJLength=parseIntegerOption(options["-min_J_length"],"-min_J_length",{minimum:0});
  for(const [option,key] of [["-D_penalty","dMismatch"],["-J_penalty","jMismatch"]]){
    if(options[option]===undefined)continue;
    const value=parseFiniteOption(options[option],option);
    // The aligner bounds mismatch scores at -1000 (engine/src/run.rs MIN_MISMATCH).
    if(!Number.isInteger(value)||value>=0||value<-1000)throw new Error(`${option} must be a negative integer from -1000 to -1.`);
    tuning[key]=value;
  }
  for(const option of ["-num_alignments_D","-num_alignments_J"]){
    if(options[option]!==undefined)parseIntegerOption(options[option],option,{minimum:1,allowZero:false});
  }
  return tuning;
}

/**
 * The aligner mode, "auto", "native" or "wasm" (a non-empty SWIG_CLI_ALIGNER sets the default),
 * and the setting that chose it, for messages. `source` names `value`.
 */
function alignerMode(value,source="--aligner"){
  const origin=value==null&&process.env.SWIG_CLI_ALIGNER?"SWIG_CLI_ALIGNER":source;
  const mode=String(value??(process.env.SWIG_CLI_ALIGNER||"auto")).toLowerCase();
  if(!ALIGNER_MODES.includes(mode))throw new Error(`${origin} must be ${ALIGNER_MODES.join(", ")}.`);
  const label=origin==="SWIG_CLI_ALIGNER"?`SWIG_CLI_ALIGNER=${mode}`:origin==="--aligner"?`--aligner ${mode}`:`${origin} "${mode}"`;
  return {mode,label};
}

/**
 * The native aligner binary for `mode`, or null for WebAssembly (saying why on stderr in auto
 * mode). `unsupported` names a setting only the WebAssembly route carries.
 */
function chooseNativeAligner(mode,assets,unsupported="",label=`--aligner ${mode}`){
  if(mode==="wasm")return null;
  if(unsupported){
    if(mode==="native")throw new Error(`${label}: ${unsupported} runs only on WebAssembly.`);
    process.stderr.write(`Using the WebAssembly aligner: ${unsupported} runs only on it.\n`);
    return null;
  }
  const notes=[];const binary=resolveNativeAligner(assets,VERSION,notes);
  if(binary)return binary;
  const reason=notes.length?notes.join("; "):"this build carries no native aligner";
  if(mode==="native")throw new Error(`${label}: ${reason}.`);
  process.stderr.write(`Using the WebAssembly aligner (${reason}).\n`);
  return null;
}

async function runVdj(rawArgs,assets){
  const options=parseVdjArguments(rawArgs);
  if(options.help){process.stdout.write(`${vdjUsage()}\n`);return;}
  if(options.version){process.stdout.write(`${VERSION}\n`);return;}
  const outputValue=options["-out"];
  if(!outputValue)throw new Error("--vdj requires -out (or --out) so AIRR rows can be streamed to an explicit destination.");
  const outfmt=String(options["-outfmt"]??"19").trim();
  if(outfmt!=="19")throw new Error("swig-cli --vdj currently emits only AIRR rearrangement format; use -outfmt 19.");
  const domain=String(options["-domain_system"]??"imgt").toLowerCase();
  if(domain!=="imgt")throw new Error("swig-cli --vdj supports only -domain_system imgt; Kabat coordinates must not be labeled as IMGT/AIRR annotations.");
  const strandValue=String(options["-strand"]??"both").toLowerCase();
  const strand={both:0,plus:1,minus:2}[strandValue];
  if(strand===undefined)throw new Error("-strand must be both, plus, or minus.");
  const minimumIdentity=options["--minimum-identity"]===undefined?0.6:parseFiniteOption(options["--minimum-identity"],"--minimum-identity");
  if(minimumIdentity<0||minimumIdentity>1)throw new Error("--minimum-identity must be between 0 and 1.");
  const threadValue=options["--workers"]??options["-num_threads"];
  let workers=threadValue===undefined?0:parseIntegerOption(threadValue,options["--workers"]!==undefined?"--workers":"-num_threads",{minimum:0});
  if(options["--workers"]===undefined&&threadValue!==undefined&&workers===0)throw new Error("-num_threads must be at least 1; use --workers 0 for automatic selection.");
  if(workers===0)workers=Math.max(1,availableParallelism());
  const {mode,label:alignerLabel}=alignerMode(options["--aligner"]);
  const requestedBatchRecords=options["--batch-records"]===undefined?0:parseIntegerOption(options["--batch-records"],"--batch-records",{minimum:0});
  const batchRecords=requestedBatchRecords||automaticBatchRecords(workers);
  const callingProfile=String(options["--calling-profile"]??"r_optimized");
  if(!isCallingProfile(callingProfile))throw new Error(`--calling-profile must be ${CALLING_PROFILES.join(", ")}.`);
  const queryValue=String(options["-query"]??"-");
  const queryPath=queryValue==="-"?"-":resolve(queryValue);
  const outputPath=String(outputValue)==="-"?"-":resolve(String(outputValue));
  if(outputPath!=="-")await mkdir(dirname(outputPath),{recursive:true});
  const prepared=await prepareVdjReferences(options,assets);
  const tuning=vdjTuning(options);
  const minimumConstantPrefixIdentity=options["--minimum-c-prefix-identity"]===undefined?null:parseFiniteOption(options["--minimum-c-prefix-identity"],"--minimum-c-prefix-identity");
  if(minimumConstantPrefixIdentity!==null&&(minimumConstantPrefixIdentity<0||minimumConstantPrefixIdentity>1))throw new Error("--minimum-c-prefix-identity must be between 0 and 1.");
  const init=workerInitialization(assets.wasmPath,prepared.references,{callingProfile,minimumIdentity,minimumConstantPrefixIdentity,strand,tuning});
  const nativeBinary=chooseNativeAligner(mode,assets,"",alignerLabel);
  let pool=nativeBinary?new NativePool(nativeBinary,workers,init,batchRecords):new WasmPool(workers,init);
  try{await pool.start();}
  catch(error){
    await pool.close().catch(()=>{});
    // Germline errors would stop WebAssembly too; anything else falls back to it.
    if(!nativeBinary||mode==="native"||/database/i.test(error.message))throw error;
    process.stderr.write(`The native aligner did not start (${error.message}); using WebAssembly.\n`);
    pool=new WasmPool(workers,init);
    try{await pool.start();}
    catch(fallbackError){await pool.close().catch(()=>{});throw fallbackError;}
  }
  const native=pool instanceof NativePool;
  const window=native?pool.window:Math.max(2,workers*2);
  // Rows go straight to the output as bytes unless -auxiliary_data rewrites FWR4.
  const bytes=native&&!Object.keys(prepared.fwr4EndOffsets).length;
  const outputHandle=outputPath==="-"?null:createVdjOutput(outputPath);
  const output=outputHandle?.stream??process.stdout;
  let outputHeader=null;let records=0;let completed=false;
  try{
    if(outputHandle)await outputHandle.opened;
    const profileLabel=callingProfile==="r_optimized"?"R-optimized":callingProfile==="sensitive_d"?"Sensitive-D":callingProfile==="igblast_compatible"?"IgBLAST-agreement":"truth-optimized";
    process.stderr.write(`Streaming Swig aligner V(D)J assignments (${native?"native aligner":"WebAssembly aligner"}; ${prepared.mode}; ${workers} ${native?"thread":"worker"}${workers===1?"":"s"}; ${batchRecords.toLocaleString()} records/batch; ${profileLabel} profile) to ${outputPath}.\n`);
    const state={inputRecords:0,eligibleRecords:0,selectedRecords:0,fastqFilter:emptyFastqQualityFilterStats(false,false)};
    const pending=[];
    const consume=async(item)=>{
      const result=await item.promise;
      if(outputHeader===null){outputHeader=result.header;await writeChunk(output,`${outputHeader}\n`);}
      else if(result.header!==outputHeader)throw new Error("The aligner changed the AIRR schema between V(D)J batches.");
      const body=trimAuxiliaryFwr4(result.header,result.body,prepared.fwr4EndOffsets,prepared.jLengths);
      await writeChunk(output,body);records+=result.count;
    };
    const preprocessing={fastqFilter:{...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter,trim3Prime:{...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter.trim3Prime}},subsample:{...DEFAULT_CLI_CONFIG.preprocessing.subsample,enabled:false}};
    const input={path:queryPath,format:"fasta"};
    for await(const batch of sequenceBatches(input,batchRecords,preprocessing,"reannotate",0,state)){
      // Await the oldest batch only once another is ready to queue: then every batch held is
      // full, which the native window counts on (a short last batch waits for finishInput).
      if(pending.length>=window)await consume(pending.shift());
      pending.push({promise:pool.run({...workerAnnotation(batch.text,batch.count),bytes})});
      if(native)await pool.ready();
    }
    if(native)pool.finishInput();
    while(pending.length)await consume(pending.shift());
    if(native)await pool.finish();
    if(outputHandle)await outputHandle.finish();
    completed=true;
  }finally{
    if(!completed&&outputHandle)outputHandle.destroy();
    await pool.close();
  }
  process.stderr.write(`Completed ${records.toLocaleString()} streaming V(D)J assignment${records===1?"":"s"}.\n`);
}

export async function runCli(assets=defaultCliAssets()){
  const args=process.argv.slice(2);
  if(args.includes("--precompute_aux")||args.includes("--precompute-aux")){await runPrepareReference(args,assets);return;}
  if(args.includes("--vdj")){await runVdj(args,assets);return;}
  const command=args[0]&&!args[0].startsWith("-")?args[0]:"run";const rest=command===args[0]?args.slice(1):args;
  if(command==="joint-germline"){await runUnifiedGermline(rest);return;}
  if(command==="personalized-germline"){await runPersonalizedGermline(rest);return;}
  if(command==="prepare-reference"){await runPrepareReference(rest,assets);return;}
  if(hasFlag(args,"--help")||command==="help"){process.stdout.write(`${usage()}\n`);return;}
  if(hasFlag(args,"--version")||command==="version"){process.stdout.write(`${VERSION}\n`);return;}
  if(hasFlag(args,"--notices")||command==="notices"){process.stdout.write(`${thirdPartyNotices()}\n`);return;}
  if(command==="init"){
    const target=positional(rest)[0]??"swig.config.json";
    const example=normalizeCliConfig({...DEFAULT_CLI_CONFIG,inputs:[{path:"reads.fastq.gz",sampleId:"sample-1",subjectId:"donor-1"}]});
    writeFileSync(target,`${JSON.stringify(example,null,2)}\n`,{flag:"wx"});process.stdout.write(`Created ${target}\n`);return;
  }
  if(command!=="run")throw new Error(`Unknown command ${command}.\n\n${usage()}`);
  const configPath=argumentValue(rest,"--config");let base=process.cwd();let raw={};
  if(configPath){const absolute=resolve(configPath);base=dirname(absolute);raw=JSON.parse(await readFile(absolute,"utf8"));}
  const inputs=positional(rest);
  if(inputs.length){raw.inputs=inputs.map((path,index)=>({path,datasetId:index?undefined:argumentValue(rest,"--dataset"),sampleId:index?undefined:argumentValue(rest,"--sample"),subjectId:index?undefined:argumentValue(rest,"--donor")}));}
  const commandOutput=argumentValue(rest,"--out");
  const configuredOutput=typeof raw.output?.directory==="string"&&raw.output.directory.trim();
  if(!commandOutput&&!configuredOutput)throw new Error("The pipeline CLI requires an explicit output directory in output.directory or --out so results can be streamed to disk.");
  if(commandOutput)raw.output={...(raw.output??{}),directory:commandOutput};
  const commandWorkers=argumentValue(rest,"--workers");
  if(commandWorkers!==undefined)raw.annotation={...(raw.annotation??{}),workers:parseIntegerOption(commandWorkers,"--workers",{minimum:0})};
  const commandAligner=argumentValue(rest,"--aligner");
  if(commandAligner!==undefined)raw.annotation={...(raw.annotation??{}),aligner:commandAligner};
  const commandCompression=argumentValue(rest,"--airr-compression");
  if(commandCompression!==undefined){if(!["none","gzip"].includes(commandCompression))throw new Error("--airr-compression must be none or gzip.");raw.output={...(raw.output??{}),airrCompression:commandCompression};}
  const config=normalizeCliConfig(raw);
  // An explicit --aligner or annotation.aligner wins over SWIG_CLI_ALIGNER.
  const aligner=alignerMode(raw.annotation?.aligner,commandAligner!==undefined?"--aligner":"annotation.aligner");
  config.annotation.aligner=aligner.mode;config.annotation.alignerLabel=aligner.label;
  if(!isCallingProfile(config.annotation.callingProfile)){
    throw new Error(`annotation.callingProfile must be ${CALLING_PROFILES.join(", ")}.`);
  }
  if(config.annotation.workers===0)config.annotation.workers=Math.max(1,Math.min(8,availableParallelism()));
  config.inputs=config.inputs.map((input)=>({...input,path:typeof input.inline==="string"?input.path:resolveFrom(base,input.path)}));
  await runPipeline(config,base,assets);
}
