import { DEFAULT_ALLELE_REFINEMENT_OPTIONS, type AlleleRefinementOptions } from "./allele-refinement/types.ts";
import { DEFAULT_MISSING_ALLELE_OPTIONS, type MissingAlleleOptions } from "./germline-evidence.ts";
import type { DenoiseOptions } from "./post-analysis-core.ts";
import type { CompiledReferences, ScopeKey } from "./reference-pack.ts";
import { DEFAULT_REPERTOIRE_SELECTION, type RepertoireSelectionOptions } from "./repertoire-selection.ts";
import { DEFAULT_PIPELINE_PLAN, type DatasetManifestEntry, type PipelinePlan, type StudyDesign } from "./study-design.ts";
import type { CallingProfile, DoubleDScreenOptions } from "./aligner-engine.ts";
import { DEFAULT_FASTQ_QUALITY_FILTER, type FastqQualityFilterOptions } from "./sequence-stream.ts";
import type { PostAnalysisSessionSnapshot } from "./session-state.ts";

export const SWIG_CLI_CONFIG_SCHEMA = 1 as const;

export type CliInputFormat = "auto" | "fasta" | "fastq" | "airr";

export interface CliInputDataset {
  path: string;
  /** Human-readable source/member name retained in lineage-study metadata. */
  inputName?: string;
  /** Small pasted inputs can be embedded directly; when present no file is opened. */
  inline?: string;
  format?: CliInputFormat;
  /** Optional compressed-byte range selecting one complete gzip member. */
  gzipRange?: { start: number; end: number };
  datasetId?: string;
  sampleId?: string;
  /** Samples sharing this value are treated as belonging to the same donor. */
  subjectId?: string;
  cohort?: string;
  timepoint?: string;
  compartment?: string;
}

export interface CliReferenceConfig {
  species: string;
  scope: ScopeKey;
  /** Prepare and validate local `files` references with the browser's metadata-transfer pipeline. */
  prepareMetadata: boolean;
  /** Exact references exported by Swig Web. These take precedence over built-ins. */
  inline?: Partial<Record<"V" | "D" | "J" | "C", string>>;
  /** Optional local FASTA replacements, resolved relative to the config file. */
  files?: Partial<Record<"V" | "D" | "J" | "C", string>>;
}

export interface CliAnnotationConfig {
  workers: number;
  batchRecords: number;
  callingProfile: CallingProfile;
  minimumIdentity: number;
  /** Optional extra gate over the complete query-covered C-reference prefix. */
  minimumConstantPrefixIdentity: number | null;
  strand: 0 | 1 | 2;
  /** Re-run the aligner on AIRR `sequence`, or trust the existing AIRR calls. */
  airrMode: "reannotate" | "preserve";
  /** swig-cli only: the native aligner where the executable carries it ("auto"), or one route. */
  aligner?: "auto" | "native" | "wasm";
  doubleD: DoubleDScreenOptions;
}

export interface CliPreprocessingConfig {
  fastqFilter: FastqQualityFilterOptions;
  /** Exact seeded reservoir sampling is applied independently to each input dataset. */
  subsample: { enabled: boolean; size: number; seed: number };
}

export interface CliPipelineConfig {
  collapse: PipelinePlan["collapse"] & { denoise: DenoiseOptions };
  chimera: PipelinePlan["chimera"] & {
    priorProbability: number;
    baseMutationProbability: number;
    mutationRates: number[];
    mutationSwitchProbability: number;
    minimumDfr: number;
    /** Precompute Viterbi parent/breakpoint labels during the repertoire scan. */
    detailed: boolean;
  };
  selection: RepertoireSelectionOptions & { enabled: boolean };
  alleleRefinement: AlleleRefinementOptions & {
    enabled: boolean;
    reassignmentPolicy: "best" | "confidence";
    applyMinimumPosterior: number;
  };
  lineage: PipelinePlan["lineage"] & { maxCandidateComparisons: number };
  shm: PipelinePlan["shm"];
  missingAlleles: MissingAlleleOptions & { enabled: boolean };
}

export interface SwigCliConfig {
  schema: typeof SWIG_CLI_CONFIG_SCHEMA;
  application: "swig-cli";
  studyName: string;
  studyDesign: StudyDesign;
  inputs: CliInputDataset[];
  references: CliReferenceConfig;
  preprocessing: CliPreprocessingConfig;
  annotation: CliAnnotationConfig;
  pipeline: CliPipelineConfig;
  output: {
    directory: string;
    prefix: string;
    airrCompression: "none" | "gzip";
    writeAnnotatedAirr: boolean;
    writeLineageStudy: boolean;
  };
}

export type PartialSwigCliConfig = Partial<Omit<SwigCliConfig, "inputs" | "references" | "preprocessing" | "annotation" | "pipeline" | "output">> & {
  inputs?: CliInputDataset[];
  references?: Partial<CliReferenceConfig>;
  preprocessing?: Partial<Omit<CliPreprocessingConfig, "fastqFilter" | "subsample">> & {
    fastqFilter?: Partial<Omit<FastqQualityFilterOptions, "trim3Prime">> & { trim3Prime?: Partial<FastqQualityFilterOptions["trim3Prime"]> };
    subsample?: Partial<CliPreprocessingConfig["subsample"]>;
  };
  annotation?: Partial<CliAnnotationConfig> & { doubleD?: Partial<DoubleDScreenOptions> };
  pipeline?: Partial<{
    collapse: Partial<Omit<CliPipelineConfig["collapse"],"denoise">> & { denoise?: Partial<DenoiseOptions> };
    chimera: Partial<CliPipelineConfig["chimera"]>;
    selection: Partial<CliPipelineConfig["selection"]>;
    alleleRefinement: Partial<CliPipelineConfig["alleleRefinement"]>;
    lineage: Partial<CliPipelineConfig["lineage"]>;
    shm: Partial<CliPipelineConfig["shm"]>;
    missingAlleles: Partial<CliPipelineConfig["missingAlleles"]>;
  }>;
  output?: Partial<SwigCliConfig["output"]>;
};

const DEFAULT_DENOISE: DenoiseOptions = {
  mode: "conservative",
  errorRate: 0.00473,
  alpha: 0.01,
  callResolution: "allele",
  ambiguity: "strict",
  minimumParentCount: 2,
  ambiguousPolicy: "exclude",
  unresolvedPolicy: "discard",
  fadNeighborThreshold: 1,
  fadMethod: 2,
  expectedZeroErrorFraction: 1,
  maximumHammingDistance: 1,
  maximumEditDistance: 2,
  minimumIndelParentRatio: 2,
  maxCandidatesPerVariant: 50_000,
  scope: "sample",
  respectConstantCall: true,
};

export const DEFAULT_CLI_CONFIG: SwigCliConfig = {
  schema: SWIG_CLI_CONFIG_SCHEMA,
  application: "swig-cli",
  studyName: "swig-study",
  studyDesign: "independent",
  inputs: [],
  references: { species: "Homo sapiens", scope: "BCR", prepareMetadata: true },
  preprocessing: {
    fastqFilter: { ...DEFAULT_FASTQ_QUALITY_FILTER, trim3Prime: { ...DEFAULT_FASTQ_QUALITY_FILTER.trim3Prime } },
    subsample: { enabled: false, size: 10_000, seed: 1 },
  },
  annotation: {
    workers: 0,
    batchRecords: 0,
    callingProfile: "r_optimized",
    minimumIdentity: 0.6,
    minimumConstantPrefixIdentity: null,
    strand: 0,
    airrMode: "preserve",
    doubleD: { mode: "off", minimumVjSpan: 40, seedLength: 11, pseudoTrim: 5, maximumPseudoMismatches: 3, minimumScoreGain: 8 },
  },
  pipeline: {
    collapse: { ...DEFAULT_PIPELINE_PLAN.collapse, enabled: true, denoise: { ...DEFAULT_DENOISE } },
    chimera: {
      ...DEFAULT_PIPELINE_PLAN.chimera,
      enabled: false,
      priorProbability: 0.05,
      baseMutationProbability: 0.05,
      mutationRates: [0,0.0179,0.0357,0.0536,0.0714,0.0893,0.1071,0.125,0.1429,0.1607,0.1786,0.1964,0.2143,0.2321,0.25],
      mutationSwitchProbability: 0,
      minimumDfr: 1,
      detailed: false,
    },
    selection: { ...DEFAULT_REPERTOIRE_SELECTION, enabled: false },
    alleleRefinement: { ...DEFAULT_ALLELE_REFINEMENT_OPTIONS, segments: [...DEFAULT_ALLELE_REFINEMENT_OPTIONS.segments], enabled: false, reassignmentPolicy: "confidence", applyMinimumPosterior: 0.8 },
    lineage: { ...DEFAULT_PIPELINE_PLAN.lineage, enabled: true, maxCandidateComparisons: 50_000 },
    shm: { ...DEFAULT_PIPELINE_PLAN.shm, enabled: true },
    missingAlleles: { ...DEFAULT_MISSING_ALLELE_OPTIONS, enabled: false },
  },
  output: { directory: "swig-output", prefix: "swig", airrCompression: "none", writeAnnotatedAirr: true, writeLineageStudy: true },
};

function finite(value: unknown, fallback: number): number {
  const parsed=Number(value);
  return Number.isFinite(parsed)?parsed:fallback;
}

function safeIdentifier(value: string, fallback: string): string {
  return value.trim().replace(/[^A-Za-z0-9_.-]+/g,"-").replace(/^-+|-+$/g,"")||fallback;
}

export function normalizeCliConfig(value: PartialSwigCliConfig): SwigCliConfig {
  const collapse=value.pipeline?.collapse??{};
  const chimera=value.pipeline?.chimera??{};
  const selection=value.pipeline?.selection??{};
  const allele=value.pipeline?.alleleRefinement??{};
  const lineage=value.pipeline?.lineage??{};
  const shm=value.pipeline?.shm??{};
  const missing=value.pipeline?.missingAlleles??{};
  const fastqFilter=value.preprocessing?.fastqFilter??{};
  const trim3Prime=fastqFilter.trim3Prime??{};
  const subsample=value.preprocessing?.subsample??{};
  const inputs=(value.inputs??[]).map((input,index)=>{
    const stem=input.path.split(/[\\/]/).pop()?.replace(/\.(?:fasta?|fna|fas|fastq|fq|airr\.tsv|tsv)(?:\.gz)?$/i,"")||`dataset-${index+1}`;
    const datasetId=safeIdentifier(input.datasetId??stem,`dataset-${index+1}`);
    const sampleId=safeIdentifier(input.sampleId??datasetId,`sample-${index+1}`);
    return {...input,format:input.format??"auto",datasetId,sampleId,subjectId:safeIdentifier(input.subjectId??sampleId,`subject-${index+1}`),cohort:input.cohort??"",timepoint:input.timepoint??"",compartment:input.compartment??""};
  });
  // Swig 0.38 configs also name a SwiftIG assigner strategy; the Rust aligner has one search.
  const {assignerStrategy:_strategy,...given}=(value.annotation??{}) as Partial<CliAnnotationConfig>&{assignerStrategy?:unknown};
  const annotation={...DEFAULT_CLI_CONFIG.annotation,...given,doubleD:{...DEFAULT_CLI_CONFIG.annotation.doubleD,...given.doubleD}};
  annotation.workers=Math.max(0,Math.floor(finite(annotation.workers,0)));
  annotation.batchRecords=Math.max(0,Math.floor(finite(annotation.batchRecords,0)));
  annotation.minimumIdentity=Math.max(0,Math.min(1,finite(annotation.minimumIdentity,0.6)));
  annotation.minimumConstantPrefixIdentity=value.annotation?.minimumConstantPrefixIdentity===null||value.annotation?.minimumConstantPrefixIdentity===undefined
    ? null
    : Math.max(0,Math.min(1,finite(value.annotation.minimumConstantPrefixIdentity,0)));
  const normalizedSubsample={...DEFAULT_CLI_CONFIG.preprocessing.subsample,...subsample};
  normalizedSubsample.enabled=Boolean(normalizedSubsample.enabled);
  normalizedSubsample.size=Math.max(1,Math.floor(finite(normalizedSubsample.size,10_000)));
  normalizedSubsample.seed=Math.trunc(finite(normalizedSubsample.seed,1));
  return {
    ...DEFAULT_CLI_CONFIG,
    ...value,
    schema:SWIG_CLI_CONFIG_SCHEMA,
    application:"swig-cli",
    inputs,
    references:{...DEFAULT_CLI_CONFIG.references,...value.references,prepareMetadata:value.references?.prepareMetadata!==false,inline:value.references?.inline?{...value.references.inline}:undefined,files:value.references?.files?{...value.references.files}:undefined},
    preprocessing:{
      fastqFilter:{...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter,...fastqFilter,trim3Prime:{...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter.trim3Prime,...trim3Prime}},
      subsample:normalizedSubsample,
    },
    annotation,
    pipeline:{
      collapse:{...DEFAULT_CLI_CONFIG.pipeline.collapse,...collapse,denoise:{...DEFAULT_DENOISE,...collapse.denoise,mode:(collapse.mode??collapse.denoise?.mode??DEFAULT_DENOISE.mode) as DenoiseOptions["mode"],scope:collapse.scope??collapse.denoise?.scope??DEFAULT_DENOISE.scope,respectConstantCall:collapse.respectConstantCall??collapse.denoise?.respectConstantCall??true}},
      chimera:{...DEFAULT_CLI_CONFIG.pipeline.chimera,...chimera,mutationRates:[...(chimera.mutationRates??DEFAULT_CLI_CONFIG.pipeline.chimera.mutationRates)]},
      selection:{...DEFAULT_CLI_CONFIG.pipeline.selection,...selection},
      alleleRefinement:{...DEFAULT_CLI_CONFIG.pipeline.alleleRefinement,...allele,segments:[...(allele.segments??DEFAULT_CLI_CONFIG.pipeline.alleleRefinement.segments)]},
      lineage:{...DEFAULT_CLI_CONFIG.pipeline.lineage,...lineage},
      shm:{...DEFAULT_CLI_CONFIG.pipeline.shm,...shm},
      missingAlleles:{...DEFAULT_CLI_CONFIG.pipeline.missingAlleles,...missing},
    },
    output:{...DEFAULT_CLI_CONFIG.output,...value.output,airrCompression:value.output?.airrCompression==="gzip"?"gzip":"none"},
  };
}

export interface BrowserCliExport {
  studyName: string;
  studyDesign: StudyDesign;
  datasets: DatasetManifestEntry[];
  references: CompiledReferences;
  species: string;
  scope: ScopeKey;
  workers: number;
  callingProfile: CallingProfile;
  minimumIdentity: number;
  minimumConstantPrefixIdentity?: number | null;
  strand: 0|1|2;
  fastqFilter?: FastqQualityFilterOptions;
  subsample?: CliPreprocessingConfig["subsample"];
  doubleD: DoubleDScreenOptions;
  pipeline: PipelinePlan;
  collapseOptions?: Partial<DenoiseOptions>;
  chimeraOptions?: Partial<CliPipelineConfig["chimera"]>;
  selectionOptions?: Partial<RepertoireSelectionOptions>;
  alleleOptions?: Partial<AlleleRefinementOptions> & { reassignmentPolicy?: "best"|"confidence"; applyMinimumPosterior?: number };
  lineageOptions?: Partial<CliPipelineConfig["lineage"]>;
  missingAlleleOptions?: Partial<MissingAlleleOptions>;
}

export interface BrowserPostCliState {
  pipeline: PipelinePlan;
  collapseOptions?: Partial<DenoiseOptions>;
  chimeraOptions?: Partial<CliPipelineConfig["chimera"]>;
  selectionOptions?: Partial<RepertoireSelectionOptions>;
  alleleOptions?: Partial<AlleleRefinementOptions> & { reassignmentPolicy?: "best"|"confidence"; applyMinimumPosterior?: number };
  lineageOptions?: Partial<CliPipelineConfig["lineage"]>;
  missingAlleleOptions?: Partial<MissingAlleleOptions>;
}

/** Translate completed interactive state back to the same explicit CLI fields used by an up-front pipeline export. */
export function cliStateFromPostAnalysis(base: PipelinePlan, post: PostAnalysisSessionSnapshot): BrowserPostCliState {
  const pipeline:PipelinePlan={
    ...DEFAULT_PIPELINE_PLAN,...base,
    collapse:{...DEFAULT_PIPELINE_PLAN.collapse,...base.collapse},
    chimera:{...DEFAULT_PIPELINE_PLAN.chimera,...base.chimera},
    selection:{...DEFAULT_PIPELINE_PLAN.selection,...base.selection},
    alleleRefinement:{...DEFAULT_PIPELINE_PLAN.alleleRefinement,...base.alleleRefinement,segments:[...(base.alleleRefinement?.segments??DEFAULT_PIPELINE_PLAN.alleleRefinement.segments)]},
    lineage:{...DEFAULT_PIPELINE_PLAN.lineage,...base.lineage},
    shm:{...DEFAULT_PIPELINE_PLAN.shm,...base.shm},
    missingAlleles:{...DEFAULT_PIPELINE_PLAN.missingAlleles,...base.missingAlleles},
  };
  if(!base.enabled){
    pipeline.collapse.enabled=false;pipeline.chimera.enabled=false;pipeline.selection.enabled=false;pipeline.alleleRefinement.enabled=false;pipeline.lineage.enabled=false;pipeline.shm.enabled=false;pipeline.missingAlleles.enabled=false;
  }
  let collapseOptions:Partial<DenoiseOptions>|undefined;
  if(post.collapse){
    const o=post.collapse.options;
    pipeline.collapse.enabled=true;
    if(["exact","fad","conservative","indel"].includes(post.collapse.mode))pipeline.collapse.mode=post.collapse.mode;
    if(["sequence","trimmed","cdr3","rearrangement"].includes(String(o.dedupKey)))pipeline.collapse.key=o.dedupKey as PipelinePlan["collapse"]["key"];
    if(["dataset","sample","subject","cohort","global"].includes(String(o.collapseScope)))pipeline.collapse.scope=o.collapseScope as PipelinePlan["collapse"]["scope"];
    if(o.denoiseUnresolvedPolicy==="discard"||o.denoiseUnresolvedPolicy==="retain")pipeline.collapse.unresolvedPolicy=o.denoiseUnresolvedPolicy;
    if(typeof o.respectConstantCall==="boolean")pipeline.collapse.respectConstantCall=o.respectConstantCall;
    if(post.collapse.mode!=="exact")collapseOptions={
      mode:post.collapse.mode,errorRate:Number(o.denoiseErrorRate),alpha:Number(o.denoiseAlpha),callResolution:o.denoiseResolution as DenoiseOptions["callResolution"],ambiguity:o.denoiseAmbiguity as DenoiseOptions["ambiguity"],minimumParentCount:Number(o.minimumParentCount),ambiguousPolicy:o.denoiseAmbiguousPolicy as DenoiseOptions["ambiguousPolicy"],unresolvedPolicy:pipeline.collapse.unresolvedPolicy,fadNeighborThreshold:Number(o.fadNeighborThreshold),fadMethod:Number(o.fadMethod) as 1|2,expectedZeroErrorFraction:Number(o.expectedZeroErrorFraction),maximumHammingDistance:Number(o.maximumDenoiseDistance),maximumEditDistance:Number(o.maximumEditDistance),minimumIndelParentRatio:Number(o.minimumIndelParentRatio),maxCandidatesPerVariant:Number(o.denoiseCandidateCap),scope:pipeline.collapse.scope,respectConstantCall:pipeline.collapse.respectConstantCall,
    };
  }
  let chimeraOptions:Partial<CliPipelineConfig["chimera"]>|undefined;
  if(post.chimera){
    const o=post.chimera.options;
    pipeline.chimera.enabled=true;
    if(o.segment==="V"||o.segment==="J")pipeline.chimera.segment=o.segment;
    if(o.method==="BW"||o.method==="DB")pipeline.chimera.model=o.method;
    pipeline.chimera.posteriorThreshold=post.chimera.filterThreshold;
    if(typeof o.retainUnevaluated==="boolean")pipeline.chimera.retainUnevaluated=o.retainUnevaluated;
    pipeline.chimera.msaSource="upload";pipeline.chimera.uploadedMsa=post.chimera.msa??"";pipeline.chimera.uploadedMsaName=String(o.uploadedMsaName??`${pipeline.chimera.segment.toLowerCase()}-reference-msa.fasta`);
    const rates=Array.isArray(o.mutationRates)?o.mutationRates.map(Number):String(o.mutationRates??"").split(/[\s,;]+/).map(Number).filter(Number.isFinite);
    chimeraOptions={priorProbability:Number(o.priorProbability),baseMutationProbability:Number(o.baseMutationProbability),mutationRates:rates,mutationSwitchProbability:Number(o.mutationSwitchProbability),minimumDfr:Number(o.minDfr),detailed:Boolean(o.detailed),posteriorThreshold:post.chimera.filterThreshold,retainUnevaluated:pipeline.chimera.retainUnevaluated,msaSource:"upload",uploadedMsa:pipeline.chimera.uploadedMsa,uploadedMsaName:pipeline.chimera.uploadedMsaName,segment:pipeline.chimera.segment,model:pipeline.chimera.model,enabled:true};
  }
  const selectionOptions=post.selection?.options;
  if(selectionOptions)pipeline.selection.enabled=true;
  const alleleOptions=post.alleleRefinement?{...post.alleleRefinement.options,reassignmentPolicy:post.alleleRefinement.reassignmentPolicy,applyMinimumPosterior:post.alleleRefinement.applyMinimumPosterior}:undefined;
  if(alleleOptions)pipeline.alleleRefinement.enabled=true;
  let lineageOptions:Partial<CliPipelineConfig["lineage"]>|undefined;
  if(post.lineage){const o=post.lineage.options;pipeline.lineage.enabled=true;if(typeof o.identity==="number")pipeline.lineage.identity=o.identity;if(o.resolution==="gene"||o.resolution==="allele")pipeline.lineage.resolution=o.resolution;if(o.ambiguity==="overlap"||o.ambiguity==="top"||o.ambiguity==="strict")pipeline.lineage.ambiguity=o.ambiguity;if(typeof o.productiveOnly==="boolean")pipeline.lineage.productiveOnly=o.productiveOnly;if(["dataset","sample","subject","cohort","global"].includes(String(o.lineageScope)))pipeline.lineage.scope=o.lineageScope as PipelinePlan["lineage"]["scope"];lineageOptions={...pipeline.lineage,maxCandidateComparisons:Number(o.candidateCap)};}
  if(post.shm){pipeline.shm.enabled=true;pipeline.shm.metric=post.shm.metric;}
  const missingAlleleOptions=post.missingAlleles?.options;
  if(missingAlleleOptions){pipeline.missingAlleles.enabled=true;pipeline.lineage.enabled=true;}
  return {pipeline,collapseOptions,chimeraOptions,selectionOptions,alleleOptions,lineageOptions,missingAlleleOptions};
}

/** Build the portable, human-editable configuration emitted by Swig Web. */
export function cliConfigFromBrowser(run: BrowserCliExport): SwigCliConfig {
  const config=normalizeCliConfig({
    studyName:run.studyName,
    studyDesign:run.studyDesign,
    inputs:run.datasets.map((dataset)=>{const source=(dataset as DatasetManifestEntry & {source?:unknown}).source;return {path:dataset.inputPath??dataset.inputName,inputName:dataset.inputName,inline:typeof source==="string"?source:undefined,format:dataset.inputFormat??"auto",gzipRange:dataset.gzipRange?{...dataset.gzipRange}:undefined,datasetId:dataset.datasetId,sampleId:dataset.sampleId,subjectId:dataset.subjectId,cohort:dataset.cohort,timepoint:dataset.timepoint,compartment:dataset.compartment??""};}),
    references:{species:run.species,scope:run.scope,inline:{V:run.references.V,D:run.references.D,J:run.references.J,C:run.references.C}},
    preprocessing:{fastqFilter:run.fastqFilter,subsample:run.subsample},
    annotation:{workers:run.workers,callingProfile:run.callingProfile,minimumIdentity:run.minimumIdentity,minimumConstantPrefixIdentity:run.minimumConstantPrefixIdentity??null,strand:run.strand,airrMode:"reannotate",doubleD:run.doubleD},
    pipeline:{
      collapse:{...run.pipeline.collapse,denoise:run.collapseOptions},
      chimera:{...run.pipeline.chimera,...run.chimeraOptions},
      selection:{...run.pipeline.selection,...run.selectionOptions},
      alleleRefinement:{...run.pipeline.alleleRefinement,...run.alleleOptions},
      lineage:{...run.pipeline.lineage,...run.lineageOptions},
      shm:{...run.pipeline.shm},
      missingAlleles:{...run.pipeline.missingAlleles,...run.missingAlleleOptions},
    },
    output:{prefix:safeIdentifier(run.studyName,"swig")},
  });
  return config;
}
