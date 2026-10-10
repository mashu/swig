import assert from "node:assert/strict";
import test from "node:test";

import { cliConfigFromBrowser, cliStateFromPostAnalysis, normalizeCliConfig } from "../src/pipeline-config.ts";
import { DEFAULT_MISSING_ALLELE_OPTIONS } from "../src/germline-evidence.ts";
import { DEFAULT_REPERTOIRE_SELECTION } from "../src/repertoire-selection.ts";
import type { PostAnalysisSessionSnapshot } from "../src/session-state.ts";
import type { DatasetManifestEntry } from "../src/study-design.ts";
import type { CompiledReferences } from "../src/reference-pack.ts";
import { DEFAULT_PIPELINE_PLAN } from "../src/study-design.ts";

const references:CompiledReferences={
  V:">IGHV1-2*01\nACGT\n",
  D:">IGHD1-1*01\nAC\n",
  J:">IGHJ4*02\nTGCA\n",
  C:"",
  counts:{V:1,D:1,J:1,C:0},
  annotation:{V:{annotated:0,total:1},J:{annotated:0,total:1}},
  loci:["IGH"],
};

test("small CLI configs receive analysis defaults and explicit donor grouping",()=>{
  const config=normalizeCliConfig({inputs:[
    {path:"day 0.fastq.gz",sampleId:"day_0",subjectId:"donor A"},
    {path:"day 30.fastq.gz",sampleId:"day_30",subjectId:"donor A"},
  ]});
  assert.equal(config.inputs[0].subjectId,"donor-A");
  assert.equal(config.inputs[1].subjectId,"donor-A");
  assert.equal(config.pipeline.collapse.enabled,true);
  assert.equal(config.pipeline.selection.enabled,false);
  assert.equal(config.pipeline.lineage.enabled,true);
  assert.equal(config.pipeline.lineage.scope,"subject");
  assert.equal(config.pipeline.shm.enabled,true);
  assert.equal(config.annotation.workers,0,"zero means choose the host default at CLI startup");
  assert.equal(config.annotation.callingProfile,"r_optimized");
  assert.equal(config.annotation.minimumConstantPrefixIdentity,null);
  assert.equal(config.output.airrCompression,"none");
  assert.equal(normalizeCliConfig({annotation:{minimumConstantPrefixIdentity:.98},output:{airrCompression:"gzip"}}).annotation.minimumConstantPrefixIdentity,.98);
  assert.equal(normalizeCliConfig({output:{airrCompression:"gzip"}}).output.airrCompression,"gzip");
  assert.equal(normalizeCliConfig({annotation:{callingProfile:"r_optimized"}}).annotation.callingProfile,"r_optimized");
  assert.equal(normalizeCliConfig({annotation:{callingProfile:"sensitive_d"}}).annotation.callingProfile,"sensitive_d");
  assert.equal(config.references.prepareMetadata,true);
  assert.equal(normalizeCliConfig({references:{prepareMetadata:false}}).references.prepareMetadata,false);
});

test("browser CLI export preserves exact references, methods, and sample-to-donor map",()=>{
  const config=cliConfigFromBrowser({
    studyName:"vaccine study",
    studyDesign:"longitudinal",
    datasets:[
      {datasetId:"d0",inputName:"day0.member.fastq.gz",inputPath:"study.fastq.gz",inputFormat:"fastq",gzipRange:{start:0,end:1234},sampleId:"day_0",subjectId:"person_1",cohort:"vaccinated",timepoint:"0",compartment:"blood"},
      {datasetId:"d30",inputName:"day30.fastq.gz",sampleId:"day_30",subjectId:"person_1",cohort:"vaccinated",timepoint:"30",compartment:"lymph"},
    ],
    references,
    species:"Homo sapiens",
    scope:"IGH",
    workers:3,
    callingProfile:"truth_optimized",
    minimumIdentity:0.7,
    minimumConstantPrefixIdentity:0.98,
    strand:0,
    fastqFilter:{enabled:true,maximumExpectedErrors:0.2,phredOffset:33,trim3Prime:{enabled:true,windowSize:6,minimumMeanPhred:25,minimumLength:100}},
    subsample:{enabled:true,size:12_345,seed:19},
    doubleD:{mode:"off",minimumVjSpan:40,seedLength:11,pseudoTrim:5,maximumPseudoMismatches:3,minimumScoreGain:8},
    pipeline:{...DEFAULT_PIPELINE_PLAN,enabled:true,selection:{...DEFAULT_PIPELINE_PLAN.selection,enabled:true,sampleId:"day_30"}},
  });
  assert.equal(config.output.prefix,"vaccine-study");
  assert.equal(config.inputs[0].subjectId,config.inputs[1].subjectId);
  assert.equal(config.inputs[0].path,"study.fastq.gz");
  assert.equal(config.inputs[0].inputName,"day0.member.fastq.gz");
  assert.deepEqual(config.inputs[0].gzipRange,{start:0,end:1234});
  assert.equal(config.references.inline?.V,references.V);
  assert.equal(config.references.inline?.J,references.J);
  assert.equal(config.pipeline.selection.enabled,true);
  assert.equal(config.pipeline.selection.sampleId,"day_30");
  assert.equal(config.annotation.workers,3);
  assert.equal(config.annotation.airrMode,"reannotate");
  assert.equal(config.annotation.minimumConstantPrefixIdentity,0.98);
  assert.equal(config.preprocessing.fastqFilter.trim3Prime.minimumMeanPhred,25);
  assert.deepEqual(config.preprocessing.subsample,{enabled:true,size:12_345,seed:19});
});

test("completed interactive settings map back to executable CLI fields without losing advanced options",()=>{
  const post={
    workingStages:[],
    collapse:{mode:"indel",options:{dedupKey:"cdr3",collapseScope:"cohort",respectConstantCall:false,denoiseErrorRate:0.012,denoiseAlpha:0.02,denoiseResolution:"gene",denoiseAmbiguity:"top",minimumParentCount:7,denoiseAmbiguousPolicy:"retain",denoiseUnresolvedPolicy:"retain",fadNeighborThreshold:1.5,fadMethod:1,expectedZeroErrorFraction:0.9,maximumDenoiseDistance:2,maximumEditDistance:4,minimumIndelParentRatio:3,denoiseCandidateCap:1234}},
    chimera:{options:{segment:"J",method:"DB",priorProbability:0.08,baseMutationProbability:0.04,mutationRates:"0.01,0.03",mutationSwitchProbability:0.02,minDfr:4,detailed:true,retainUnevaluated:false,uploadedMsaName:"j.msa.fasta"},msa:">IGHJ1*01\nAC-G\n>IGHJ2*01\nACTG\n",filterThreshold:0.81},
    selection:{options:{...DEFAULT_REPERTOIRE_SELECTION,motif:"CAR",motifSyntax:"regex"}},
    lineage:{options:{identity:0.91,resolution:"allele",ambiguity:"strict",productiveOnly:false,candidateCap:4321,lineageScope:"global"}},
    shm:{metric:"synonymous",dashboard:{}},
    missingAlleles:{options:{...DEFAULT_MISSING_ALLELE_OPTIONS,minimumIndependentUnits:11},dashboard:{}},
  } as unknown as PostAnalysisSessionSnapshot;
  // The browser may have started in interactive mode; only stages with a
  // committed result should become enabled in its Results-page export.
  const state=cliStateFromPostAnalysis({...DEFAULT_PIPELINE_PLAN,enabled:false},post);
  const config=cliConfigFromBrowser({studyName:"interactive",studyDesign:"custom",datasets:[{datasetId:"d",inputName:"reads.airr.tsv",sampleId:"s",subjectId:"p",cohort:"",timepoint:""}],references,species:"Homo sapiens",scope:"IGH",workers:2,callingProfile:"truth_optimized",minimumIdentity:0.6,strand:0,doubleD:{mode:"off",minimumVjSpan:40,seedLength:11,pseudoTrim:5,maximumPseudoMismatches:3,minimumScoreGain:8},...state});
  assert.equal(config.pipeline.collapse.mode,"indel");
  assert.equal(config.pipeline.collapse.scope,"cohort");
  assert.equal(config.pipeline.collapse.denoise.maximumEditDistance,4);
  assert.equal(config.pipeline.chimera.segment,"J");
  assert.equal(config.pipeline.chimera.posteriorThreshold,0.81);
  assert.equal(config.pipeline.chimera.minimumDfr,4);
  assert.equal(config.pipeline.chimera.detailed,true);
  assert.equal(config.pipeline.chimera.uploadedMsaName,"j.msa.fasta");
  assert.equal(config.pipeline.selection.motif,"CAR");
  assert.equal(config.pipeline.lineage.maxCandidateComparisons,4321);
  assert.equal(config.pipeline.lineage.scope,"global");
  assert.equal(config.pipeline.shm.metric,"synonymous");
  assert.equal(config.pipeline.missingAlleles.minimumIndependentUnits,11);
});

test("pre-run export embeds a pasted dataset instead of emitting an unusable text-file path",()=>{
  const pasted={datasetId:"paste",inputName:"pasted-sequences.txt",inputPath:"pasted-sequences.txt",inputFormat:"fasta",sampleId:"sample",subjectId:"donor",cohort:"",timepoint:"",source:">read\nACGT\n"} as DatasetManifestEntry & {source:string};
  const config=cliConfigFromBrowser({studyName:"paste",studyDesign:"independent",datasets:[pasted],references,species:"Homo sapiens",scope:"IGH",workers:1,callingProfile:"truth_optimized",minimumIdentity:0.6,strand:0,doubleD:{mode:"off",minimumVjSpan:40,seedLength:11,pseudoTrim:5,maximumPseudoMismatches:3,minimumScoreGain:8},pipeline:{...DEFAULT_PIPELINE_PLAN,enabled:true}});
  assert.equal(config.inputs[0].format,"fasta");
  assert.equal(config.inputs[0].inline,">read\nACGT\n");
});

for(const callingProfile of ["r_optimized","sensitive_d","truth_optimized","igblast_compatible"] as const){
  test(`browser CLI export preserves the ${callingProfile} calling profile`,()=>{
    const config=cliConfigFromBrowser({studyName:"profile",studyDesign:"independent",datasets:[],references,species:"Macaca mulatta",scope:"IGH",workers:2,callingProfile,minimumIdentity:0.6,strand:1,doubleD:{mode:"off",minimumVjSpan:40,seedLength:11,pseudoTrim:5,maximumPseudoMismatches:3,minimumScoreGain:8},pipeline:{...DEFAULT_PIPELINE_PLAN,enabled:false}});
    assert.equal(config.annotation.callingProfile,callingProfile);
    assert.equal("assignerStrategy" in config.annotation,false);
  });
}

test("assignment defaults to R-optimized; a saved assignment strategy is dropped, its profile kept",()=>{
  const defaults=normalizeCliConfig({});
  assert.equal(defaults.annotation.callingProfile,"r_optimized");
  assert.equal("assignerStrategy" in defaults.annotation,false);
  const saved=normalizeCliConfig({annotation:{assignerStrategy:"riat_mp",callingProfile:"truth_optimized"}});
  assert.equal("assignerStrategy" in saved.annotation,false);
  assert.equal(saved.annotation.callingProfile,"truth_optimized");
  assert.equal(normalizeCliConfig({annotation:{assignerStrategy:"riat_mp"}}).annotation.callingProfile,"r_optimized");
});
