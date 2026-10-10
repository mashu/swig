import { overallAssignmentProgress } from "./assignment-progress";
import {
  ChangeEvent,
  DragEvent,
  useEffect,
  useLayoutEffect,
  useMemo,
  useRef,
  useState,
} from "react";

import { AlignmentViewer } from "./alignment-view";
import { CommitNumberInput } from "./commit-number-input";
import { CommitTextInput } from "./commit-text-input";
import { FacetPicker, uniqueFacetItems } from "./facet-picker";
import { callFacetItems } from "./call-facets";
import { annotationCoverage, type GermlinePreprocessReport } from "./germline-preprocess";
import { preprocessGermlinesInWorker } from "./germline-preprocess-client";
import { RepertoireDashboard } from "./repertoire-charts";
import { PostAnalysisWorkbench, type PostAnalysisSessionHandle } from "./post-analysis";
import { tableExtension, type TableExportFormat } from "./export-formats";
import { compressedExtension, compressedMime, compressedName, outputWriter, type OutputCompression } from "./gzip-output";
import {
  allelesToFasta,
  availableScopes,
  compileReferences,
  germlineTemplateTiers,
  loadReferencePack,
  lociForScope,
  makeDemoFasta,
  type ReferencePack,
  type ReferenceSpecies,
  type CompiledReferences,
  type LocusKey,
  type ScopeKey,
  type SegmentKey,
} from "./reference-pack";
import {
  composeReferenceOverrides,
  referenceCellKey,
  segmentAppliesToLocus,
} from "./reference-composition";
import {
  collectionsForDatabase,
  databasesForCell,
  databaseOptionsFor,
  DEFAULT_DATABASE_ID,
  loadCollectionSegment,
  preferredDatabaseIdFor,
  type ReferenceDatabase,
} from "./reference-catalog";
import {
  AirrResultStore,
  EMPTY_FILTERS,
  inferIsotype,
  type AirrOutputHandle,
  type AirrOutputWritable,
  type AirrIndexRecord,
  type DoubleDEvidenceRecord,
  type DirectAirrOutput,
  type ResultFacets,
  type ResultFilters,
  type ResultPage,
} from "./result-store";
import { callingProfileOrDefault, type CallingProfile } from "./aligner-engine";
import {
  runAligner,
  withAnalysisWebLock,
  type DoubleDScreenMode,
  type DoubleDScreenOptions,
} from "./aligner-runtime";
import {
  initializingAssignmentTelemetry,
  type AssignmentTelemetry,
} from "./assignment-telemetry";
import {
  addFastqQualityFilterStats,
  DEFAULT_FASTQ_QUALITY_FILTER,
  emptyFastqQualityFilterStats,
  streamSequenceBatches,
  estimateSequenceCharacters,
  type FastqQualityFilterOptions,
  type FastqQualityFilterStats,
  type SequenceSource,
} from "./sequence-stream";
import {
  gzipMemberFile,
  gzipMemberSource,
  inspectGzipMembers,
  suggestedGzipMemberName,
  type GzipMemberInspection,
} from "./gzip-members";
import { prepareReferenceMsa } from "./post-analysis-core";
import { runKalignTask } from "./biowasm-task-runtime";
import { cliConfigFromBrowser, cliStateFromPostAnalysis } from "./pipeline-config";
import { decodeLineageStudy, type LineageStudyManifest } from "./lineage-study";
import { LineageStudyPage } from "./lineage-study-page";
import { ReferenceAlleleExclusionEditor } from "./reference-allele-exclusion";
import { filterReferenceFasta } from "./reference-fasta";
import { createSampleColorMap, sampleColor, sampleIds, type SampleColorMap } from "./sample-colors";
import { decodeSession, encodeSession, linkedAirrMatches, sessionBaseName, SWIG_SESSION_SCHEMA, type PostAnalysisSessionSnapshot, type SwigSession } from "./session-state";
import {
  candidatesFromDirectoryPicker,
  collectDroppedInput,
  donorForFlatRoot,
  inferDirectoryDonors,
  type InputFileCandidate,
} from "./directory-input";
import {
  activeProjectRun,
  appendProjectLog,
  attachProjectDirectory,
  loadActiveProjectFiles,
  prepareProjectRun,
  projectDirectoriesSupported,
  saveProjectCheckpoint,
  selectProjectDirectory,
  writeProjectDatasetManifest,
  type PreparedProjectRun,
  type ProjectWorkspace,
} from "./project-directory";
import {
  annotateAirrBatch,
  annotateDoubleDBatch,
  DATASET_SCOPE_LABELS,
  ALLELE_POOL_SCOPE_LABELS,
  DEFAULT_PIPELINE_PLAN,
  stableDatasetSeed,
  studyScopeDefaults,
  type DatasetManifestEntry,
  type DatasetScope,
  type PipelinePlan,
  type StudyDesign,
} from "./study-design";

type AppPage = "home" | "analyze" | "results" | "lineage-study";
type WebAnalysisMode = "vdj" | "repertoire" | "lineage";
type InputFormat = "FASTA" | "FASTQ" | "AIRR TSV";
type InputSource = "upload" | "paste";
type OutputStorageMode = "auto" | "browser" | "disk";
type AirrRow = Record<string, string>;

interface InputData {
  name: string;
  source: SequenceSource;
  count: number | null;
  format: InputFormat;
  formatCode: 1 | 2 | 3;
  size: number;
  inputPath?: string;
  inputFormat?: "auto" | "fasta" | "fastq" | "airr";
  gzipRange?: { start: number; end: number };
  gzipMemberCount?: number;
  gzipMembersMerged?: boolean;
}

interface DatasetInput extends InputData, DatasetManifestEntry {}

interface ReferenceOverride {
  name: string;
  text: string;
  count: number;
  size: number;
  report: GermlinePreprocessReport;
  sourceKind: "database" | "upload";
  sourceDatabaseId?: string;
}

type ReferenceCellMap = Record<string, ReferenceOverride>;
type ReferenceAlleleExclusionMap = Record<string, string[]>;

interface ResultSession {
  id: number;
  webMode: WebAnalysisMode;
  store: AirrResultStore;
  total: number;
  seconds: number;
  inputName: string;
  datasets: DatasetManifestEntry[];
  studyDesign: StudyDesign;
  pipeline: PipelinePlan;
  species: string;
  scope: ScopeKey;
  facets: ResultFacets;
  summary: { assigned: number; productive: number; withCdr3: number };
  workers: number;
  outputBytes: number;
  streamedDirectly: boolean;
  inputTotal: number;
  subsampleSize: number | null;
  subsampleSeed: number | null;
  fastqFilter: FastqQualityFilterOptions;
  fastqFilterStats: FastqQualityFilterStats;
  references: CompiledReferences;
  referenceExclusions: ReferenceAlleleExclusionMap;
  callingProfile: CallingProfile;
  minimumIdentity: number;
  minimumConstantPrefixIdentity: number | null;
  strand: 0 | 1 | 2;
  doubleD: DoubleDScreenOptions;
  doubleDCount: number;
  sampleColors: SampleColorMap;
  postAnalysis?: PostAnalysisSessionSnapshot;
  restored?: boolean;
  project?: ProjectWorkspace;
  projectStatus?: string;
}

const APP_VERSION = "0.39.0";
const SEGMENTS: SegmentKey[] = ["V", "D", "J", "C"];
const PAGE_SIZE = 50;
const MAX_INLINE_COUNT_BYTES = 2 * 1024 * 1024;
const DIRECT_OUTPUT_COUNT = 100_000;
const DIRECT_OUTPUT_INPUT_BYTES = 32 * 1024 * 1024;

function isAbortError(error: unknown): boolean {
  return typeof error === "object" && error !== null && "name" in error && (error as { name?: unknown }).name === "AbortError";
}

interface SaveFileHandle extends AirrOutputHandle {
  createWritable: () => Promise<AirrOutputWritable>;
}

type SavePicker = (options: {
  suggestedName: string;
  types: Array<{ description: string; accept: Record<string, string[]> }>;
}) => Promise<SaveFileHandle>;

let downloadWorkerRegistration: Promise<ServiceWorkerRegistration> | null = null;

function browserWorkerLimit(): number {
  if (typeof navigator === "undefined") return 4;
  const cores = Math.max(1, navigator.hardwareConcurrency || 4);
  return Math.max(1, Math.min(16, cores - 1 || 1));
}

function recommendedWorkerCount(): number {
  const hardwareLimit = browserWorkerLimit();
  const memory = (navigator as Navigator & { deviceMemory?: number }).deviceMemory;
  const memoryLimit = memory === undefined ? 8 : memory <= 2 ? 2 : memory <= 4 ? 4 : memory <= 8 ? 8 : 12;
  return Math.max(1, Math.min(hardwareLimit, memoryLimit));
}

const FIELD_HELP: Record<string,string> = {
  "species / strain":"Selects the organism whose germline references and locus definitions are used for assignment.",
  "database":"Applies a published reference preset. Individual locus and segment sources can still be replaced below.",
  "receptor":"Chooses immunoglobulin (BCR/IG) or T-cell receptor (TCR/TR) locus definitions.",
  "chain / locus":"Limits assignment to one chain or searches all compatible chains for the selected receptor.",
  "calling profile":"Selects the segment scoring tuples, the D exact-run floor and the D prior of the aligner's path model.",
  "search strand":"Controls whether forward, reverse-complement, or both query orientations are tested.",
  "parallel compute workers":"Sets the browser worker pool used for V(D)J assignment and independently parallel post-analysis partitions.",
  "airr results destination":"Chooses whether result batches remain in browser storage or are written incrementally to a file.",
  "minimum alignment identity":"Rejects segment alignments below this aligned-base identity fraction.",
  "maximum expected errors":"Filters a FASTQ read when the sum of base-error probabilities exceeds this value.",
  "quality encoding":"Interprets FASTQ quality characters using the selected Phred offset.",
  "collapse boundary":"Defines the largest study unit across which two records may be compared for collapse or denoising.",
  "study boundary":"Defines the study unit within which records may be assigned to the same lineage.",
  "cdr3 identity":"Sets the minimum nucleotide identity used to connect equal-length CDR3s by single linkage.",
  "double-d / vddj screening":"Enables an optional second-D evidence screen without changing standard calls when disabled.",
  "posterior threshold":"Sets the CHMMAIRRa posterior probability used to flag candidate chimeras.",
  "sample membership":"Filters lineages by whether they contain one, several, or specifically selected samples.",
};

function fieldHelpText(label:string,detail:string):string {
  const normalized=label.trim().replace(/\s+/g," ").toLowerCase();
  if(FIELD_HELP[normalized])return FIELD_HELP[normalized];
  if(normalized.includes("identity"))return `Sets the identity criterion for ${label.toLowerCase()}.`;
  if(normalized.includes("minimum"))return `Sets the minimum accepted value for ${label.replace(/^Minimum\s+/i,"").toLowerCase()}.`;
  if(normalized.includes("maximum"))return `Sets the maximum accepted value for ${label.replace(/^Maximum\s+/i,"").toLowerCase()}.`;
  if(normalized.includes("include")||normalized.includes("retain")||normalized.includes("require")||normalized.includes("keep"))return `Controls whether ${label.toLowerCase()} is included in this step.`;
  if(normalized.includes("method")||normalized.includes("mode")||normalized.includes("strategy"))return `Selects the ${label.toLowerCase()} used for this action.`;
  if(normalized.includes("call"))return `Filters or groups records using the selected ${label.toLowerCase()}.`;
  return detail ? `${label}. ${detail}` : `Controls ${label.toLowerCase()} for this action.`;
}

function addFieldHelp(root:ParentNode):void {
  root.querySelectorAll<HTMLLabelElement>("label").forEach((label)=>{
    if(label.dataset.fieldHelp||label.title)return;
    const control=label.querySelector<HTMLElement>("input,select,textarea");
    if(!control)return;
    const heading=label.querySelector<HTMLElement>(":scope > span");
    const emphasized=heading?.querySelector<HTMLElement>("strong,b");
    const raw=(emphasized?.textContent||heading?.textContent||control.getAttribute("aria-label")||"").trim();
    if(!raw)return;
    const detail=label.querySelector<HTMLElement>(":scope > small, :scope > span > small")?.textContent?.trim()||"";
    const help=label.dataset.help||fieldHelpText(raw,detail);
    label.dataset.fieldHelp=help;
    label.title=help;
    if(!control.title)control.title=help;
  });
  root.querySelectorAll<HTMLElement>("[role='radio'],.mode-toggle button,.receptor-selector button").forEach((control)=>{
    if(control.title)return;
    const name=(control.querySelector("strong,b")?.textContent||control.textContent||"").trim().replace(/\s+/g," ");
    const detail=control.querySelector("small")?.textContent?.trim()||"";
    if(name)control.title=detail?`${name}. ${detail}`:`Select ${name}.`;
  });
}

function outputName(inputName: string): string {
  return `${inputName.replace(/(\.gz)?\.[^.]+$/, "") || "swig"}.airr.tsv`;
}

function formattedOutputName(inputName:string,format:TableExportFormat):string{return `${inputName.replace(/(\.gz)?\.[^.]+$/, "")||"swig"}.airr${tableExtension(format)}`;}

function doubleDOutputName(inputName: string): string {
  return `${inputName.replace(/(\.gz)?\.[^.]+$/, "") || "swig"}.double-d-evidence.tsv`;
}

function formattedDoubleDOutputName(inputName:string,format:TableExportFormat):string{return `${inputName.replace(/(\.gz)?\.[^.]+$/, "")||"swig"}.double-d-evidence${tableExtension(format)}`;}

function likelyLargeInput(input: InputData, selectedCount?: number | null): boolean {
  if (selectedCount) return selectedCount >= DIRECT_OUTPUT_COUNT;
  if (input.count !== null) return input.count >= DIRECT_OUTPUT_COUNT;
  const gzipThreshold = 8 * 1024 * 1024;
  return input.size >= (input.name.toLowerCase().endsWith(".gz") ? gzipThreshold : DIRECT_OUTPUT_INPUT_BYTES);
}

function savePicker(): SavePicker | undefined {
  return (window as Window & { showSaveFilePicker?: SavePicker }).showSaveFilePicker;
}

function registerDownloadWorker(): Promise<ServiceWorkerRegistration> {
  if (!("serviceWorker" in navigator) || !window.isSecureContext) {
    return Promise.reject(new Error("Streaming downloads require a secure browser context."));
  }
  downloadWorkerRegistration ??= navigator.serviceWorker.register(
    `${import.meta.env.BASE_URL}download-worker.js`,
    { scope: import.meta.env.BASE_URL },
  );
  return downloadWorkerRegistration;
}

const LOCUS_LABELS: Record<string, string> = {
  BCR: "All BCR chains",
  TCR: "All TCR chains",
  IGH: "Heavy · IGH",
  IGK: "Kappa light · IGK",
  IGL: "Lambda light · IGL",
  TRA: "Alpha · TRA",
  TRB: "Beta · TRB",
  TRD: "Delta · TRD",
  TRG: "Gamma · TRG",
};

function bytes(value: number): string {
  if (value < 1024) return `${value} B`;
  if (value < 1024 * 1024) return `${(value / 1024).toFixed(1)} KB`;
  if (value < 1024 * 1024 * 1024) return `${(value / (1024 * 1024)).toFixed(1)} MB`;
  return `${(value / (1024 * 1024 * 1024)).toFixed(2)} GB`;
}

function percentage(value: number, total: number): string {
  return total ? `${Math.round(value / total * 100)}%` : "0%";
}

function friendlySpecies(value: string): string {
  return value.replaceAll("_", " · ");
}

function callingProfileLabel(value: CallingProfile): string {
  if (value === "igblast_compatible") return "IgBLAST-agreement";
  if (value === "r_optimized") return "R-optimized";
  if (value === "sensitive_d") return "Sensitive-D";
  return "Truth-optimized";
}

function favoriteSpecies(species: ReferenceSpecies[]): ReferenceSpecies[] {
  const preferred = [
    "Homo sapiens",
    "Mus musculus_C57BL/6",
    "Macaca mulatta_AG07107",
    "Canis lupus familiaris_boxer",
    "Felis catus_Abyssinian",
    "Bos taurus_Hereford",
    "Danio rerio_Tuebingen",
  ];
  return [...species].sort((a, b) => {
    const aIndex = preferred.indexOf(a.name);
    const bIndex = preferred.indexOf(b.name);
    if (aIndex >= 0 || bIndex >= 0) {
      if (aIndex < 0) return 1;
      if (bIndex < 0) return -1;
      return aIndex - bIndex;
    }
    return a.name.localeCompare(b.name);
  });
}

function formatFromName(name: string): { format: InputFormat; formatCode: 1 | 2 | 3 } | null {
  const plain = name.toLowerCase().replace(/\.gz$/, "");
  if (/\.(fa|fasta|fna|fas)$/.test(plain)) return { format: "FASTA", formatCode: 1 };
  if (/\.(fq|fastq)$/.test(plain)) return { format: "FASTQ", formatCode: 2 };
  if (/\.(tsv|csv)$/.test(plain)) return { format: "AIRR TSV", formatCode: 3 };
  return null;
}

function cliInputFormat(formatCode: 1|2|3): "fasta"|"fastq"|"airr" {
  return formatCode===1?"fasta":formatCode===2?"fastq":"airr";
}

function detectFormat(text: string): { format: InputFormat; formatCode: 1 | 2 | 3 } {
  const first = text.trimStart()[0];
  if (first === ">") return { format: "FASTA", formatCode: 1 };
  if (first === "@") return { format: "FASTQ", formatCode: 2 };
  const header = text.split(/\r?\n/, 1)[0]?.split(/\t|,/g) ?? [];
  if (header.includes("sequence")) return { format: "AIRR TSV", formatCode: 3 };
  throw new Error("Expected FASTA, FASTQ, or an AIRR table with a sequence column.");
}

function countTextRecords(text: string, formatCode: number): number {
  if (formatCode === 1) return (text.match(/^>/gm) ?? []).length;
  if (formatCode === 2) {
    const lines = text.split(/\r?\n/);
    let count = 0;
    let position = 0;
    while (position < lines.length) {
      while (position < lines.length && !lines[position]) position += 1;
      if (position >= lines.length) break;
      if (!lines[position].startsWith("@")) break;
      count += 1;
      position += 1;
      let sequenceLength = 0;
      while (position < lines.length && !lines[position].startsWith("+")) {
        sequenceLength += lines[position].replace(/\s/g, "").length;
        position += 1;
      }
      position += 1;
      let qualityLength = 0;
      while (position < lines.length && qualityLength < sequenceLength) {
        qualityLength += lines[position].length;
        position += 1;
      }
    }
    return count;
  }
  return Math.max(0, text.trimEnd().split(/\r?\n/).filter(Boolean).length - 1);
}

function inspectText(name: string, text: string, size = text.length): InputData {
  const detected = detectFormat(text);
  const count = countTextRecords(text, detected.formatCode);
  if (!count) throw new Error(`No ${detected.format} sequence records were found.`);
  return { name, source: text, count, size, ...detected, inputPath:name, inputFormat:cliInputFormat(detected.formatCode) };
}

async function inspectFile(file: File): Promise<InputData> {
  const namedFormat = formatFromName(file.name);
  if (file.name.toLowerCase().endsWith(".gz")) {
    if (!namedFormat) throw new Error("Name gzip inputs with .fasta.gz, .fastq.gz, or .tsv.gz.");
    return { name: file.name, source: file, count: null, size: file.size, ...namedFormat };
  }
  const sample = await file.slice(0, Math.min(file.size, MAX_INLINE_COUNT_BYTES)).text();
  const detected = namedFormat ?? detectFormat(sample);
  const count = file.size <= MAX_INLINE_COUNT_BYTES ? countTextRecords(sample, detected.formatCode) : null;
  if (count === 0) throw new Error(`No ${detected.format} sequence records were found.`);
  return { name: file.name, source: file, count, size: file.size, ...detected };
}

function inputStem(name: string): string {
  return name.replace(/\.gz$/i, "").replace(/\.(fa|fasta|fna|fas|fq|fastq|tsv|csv|txt)$/i, "").replace(/[^A-Za-z0-9_.-]+/g, "_") || "sample";
}

function datasetInput(input: InputData, ordinal: number): DatasetInput {
  const stem = inputStem(input.name);
  return {
    ...input,
    datasetId: `dataset_${ordinal}`,
    inputName: input.name,
    // Keep independent libraries separated even when files from different
    // directories share the same basename. Technical replicates are merged
    // only by an explicit shared sample ID or the technical-replicate preset.
    sampleId: `${stem}_${ordinal}`,
    subjectId: `subject_${ordinal}`,
    cohort: "cohort_1",
    timepoint: "",
    compartment: "",
    records: input.count,
  };
}

interface DirectoryDatasetInput extends DatasetInput {
  directoryRoot?: string;
  nestedDirectoryDonor?: string;
}

type GzipImportChoice = "separate" | "merge" | "cancel";

interface PendingGzipImport {
  inputName: string;
  format: InputFormat;
  members: GzipMemberInspection[];
  suggestedNames: string[];
}

function copyPipeline(value?: PipelinePlan): PipelinePlan {
  const source = value ?? DEFAULT_PIPELINE_PLAN;
  const legacyChimera = source.chimera as PipelinePlan["chimera"] & { excludedAlleles?: unknown };
  const { excludedAlleles: _discardedLegacyChimeraExclusions, ...chimera } = legacyChimera;
  return {
    ...DEFAULT_PIPELINE_PLAN,
    ...source,
    collapse: { ...DEFAULT_PIPELINE_PLAN.collapse, ...source.collapse },
    chimera: { ...DEFAULT_PIPELINE_PLAN.chimera, ...chimera },
    selection: { ...DEFAULT_PIPELINE_PLAN.selection, ...source.selection },
    alleleRefinement: { ...DEFAULT_PIPELINE_PLAN.alleleRefinement, ...source.alleleRefinement, segments: [...(source.alleleRefinement?.segments ?? DEFAULT_PIPELINE_PLAN.alleleRefinement.segments)] },
    lineage: { ...DEFAULT_PIPELINE_PLAN.lineage, ...source.lineage },
    shm: { ...DEFAULT_PIPELINE_PLAN.shm, ...source.shm },
    missingAlleles: { ...DEFAULT_PIPELINE_PLAN.missingAlleles, ...source.missingAlleles },
  };
}

function copyFastqQualityFilter(value?: FastqQualityFilterOptions): FastqQualityFilterOptions {
  const source = value ?? DEFAULT_FASTQ_QUALITY_FILTER;
  return {
    ...DEFAULT_FASTQ_QUALITY_FILTER,
    ...source,
    trim3Prime: { ...DEFAULT_FASTQ_QUALITY_FILTER.trim3Prime, ...source.trim3Prime },
  };
}

function fastqQualityResultText(session: ResultSession): string {
  if (!session.fastqFilter.enabled) return "";
  const stats = session.fastqFilterStats;
  const fastqRetained = Math.max(0, stats.recordsRetained - stats.recordsPassedThrough);
  const rejected = stats.recordsRejectedExpectedErrors + stats.recordsRejectedMinimumLength;
  const applied = stats.recordsEvaluated
    ? `FASTQ QC ${stats.recordsEvaluated.toLocaleString()} → ${fastqRetained.toLocaleString()} (${rejected.toLocaleString()} rejected${stats.recordsTrimmed ? `; ${stats.recordsTrimmed.toLocaleString()} end-trimmed` : ""})`
    : "FASTQ QC enabled; no FASTQ records present";
  return stats.recordsPassedThrough
    ? ` · ${applied} · ${stats.recordsPassedThrough.toLocaleString()} non-FASTQ records passed through`
    : ` · ${applied}`;
}

async function readUploadedText(file: File, signal?: AbortSignal): Promise<string> {
  if (signal?.aborted) throw new DOMException("File reading was cancelled.", "AbortError");
  if (!file.name.toLowerCase().endsWith(".gz")) {
    const text = await file.text();
    if (signal?.aborted) throw new DOMException("File reading was cancelled.", "AbortError");
    return text;
  }
  if (!("DecompressionStream" in globalThis)) {
    throw new Error("This browser cannot decompress gzip files. Decompress the file first.");
  }
  return new Response(file.stream().pipeThrough(new DecompressionStream("gzip"), { signal })).text();
}

function downloadBlob(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

async function sessionArtifact(
  session: ResultSession,
  datasets: DatasetManifestEntry[],
  sampleColors: SampleColorMap,
  postAnalysis: PostAnalysisSessionSnapshot,
  signal?: AbortSignal,
): Promise<SwigSession> {
  if(signal?.aborted)throw new DOMException("Session saving was cancelled.","AbortError");
  const projectRun = session.project ? activeProjectRun(session.project) : null;
  return {
    schema: SWIG_SESSION_SCHEMA,
    application: "Swig",
    applicationVersion: APP_VERSION,
    savedAt: new Date().toISOString(),
    linkedAirr: {
      name: projectRun?.airrPath.split("/").pop() || outputName(session.inputName),
      size: session.outputBytes,
      lastModified: 0,
      records: session.store.count,
      headers: [...session.store.airrHeaders],
      fingerprint: session.store.fingerprint,
    },
    analysis: {
      webMode: session.webMode,
      inputName: session.inputName,
      species: session.species,
      scope: session.scope,
      workers: session.workers,
      callingProfile: session.callingProfile,
      minimumIdentity: session.minimumIdentity,
      minimumConstantPrefixIdentity: session.minimumConstantPrefixIdentity,
      strand: session.strand,
      subsample: { enabled:session.subsampleSize!==null,size:session.subsampleSize??10_000,seed:session.subsampleSeed??1 },
      fastqFilter: copyFastqQualityFilter(session.fastqFilter),
      fastqFilterStats: { ...session.fastqFilterStats },
      references: session.references,
      referenceExclusions: Object.fromEntries(Object.entries(session.referenceExclusions).map(([key, names]) => [key, [...names]])),
      doubleD: { ...session.doubleD },
      datasets,
      studyDesign: session.studyDesign,
      pipeline: session.pipeline,
      sampleColors,
    },
    doubleD: await session.store.doubleDRecords(signal),
    postAnalysis,
  };
}

function Brand() {
  return (
    <span className="brand">
      <span className="brand-glyph" aria-hidden="true"><i /><i /><i /></span>
      <span><b>SWIG</b><small>SwiftIG · WebAssembly</small></span>
    </span>
  );
}

function AppHeader({ page, hasResults, hasLineageStudy, projectName, onNavigate, onLoadSession, onLoadLineageStudy, onOpenProject }: {
  page: AppPage;
  hasResults: boolean;
  hasLineageStudy:boolean;
  projectName?: string;
  onNavigate: (page: AppPage) => void;
  onLoadSession: () => void;
  onLoadLineageStudy:()=>void;
  onOpenProject: () => void;
}) {
  return (
    <header className="app-header">
      <button className="brand-button" type="button" onClick={() => onNavigate("home")}><Brand /></button>
      <nav aria-label="Primary navigation">
        <button className={page === "home" ? "active" : ""} type="button" onClick={() => onNavigate("home")}>Overview</button>
        <button className={page === "analyze" ? "active" : ""} type="button" onClick={() => onNavigate("analyze")}>Analyze</button>
        {hasResults && <button className={page === "results" ? "active" : ""} type="button" onClick={() => onNavigate("results")}>Results</button>}
        {hasLineageStudy&&<button className={page==="lineage-study"?"active":""} type="button" onClick={()=>onNavigate("lineage-study")}>Lineage study</button>}
        <button type="button" onClick={onLoadSession}>Load full session</button>
        <button type="button" onClick={onLoadLineageStudy}>Load lineage study</button>
        <button type="button" onClick={onOpenProject}>{projectName ? `Project · ${projectName}` : "Load project directory"}</button>
      </nav>
      <span className="local-badge"><i /> Query records remain in this browser</span>
    </header>
  );
}

function WorkflowStepper({ active }: { active: 1 | 2 | 3 }) {
  return (
    <ol className="workflow-stepper" aria-label="Analysis workflow">
      {["Choose data", "Run annotation", "Explore results"].map((label, index) => {
        const value = index + 1;
        return <li key={label} className={value === active ? "active" : value < active ? "complete" : ""}><span>{value < active ? "✓" : value}</span><b>{label}</b></li>;
      })}
    </ol>
  );
}

function LandingPage({ onChoose }: { onChoose: (mode: WebAnalysisMode) => void }) {
  return (
    <main className="workflow-landing-page">
      <section className="workflow-choice-grid" aria-label="Choose a Swig workflow">
        <button type="button" onClick={() => onChoose("vdj")}>
          <span>01</span>
          <strong>VDJ assignment, annotation, and visualization</strong>
          <i>→</i>
        </button>
        <button type="button" onClick={() => onChoose("repertoire")}>
          <span>02</span>
          <strong>End-to-end repertoire analysis</strong>
          <i>→</i>
        </button>
        <button type="button" onClick={() => onChoose("lineage")}>
          <span>03</span>
          <strong>Single lineage analysis</strong>
          <i>→</i>
        </button>
      </section>
      <p className="workflow-local-statement"><i /> Swig runs on your machine. Your sequence data are not uploaded.</p>
    </main>
  );
}

function ReferenceCellControl({ speciesName, locus, segment, builtInFasta, reference, excluded, busy, pendingSourceId, onSelect, onFile, onEditAlleles }: {
  speciesName: string;
  locus: LocusKey;
  segment: SegmentKey;
  builtInFasta: string;
  reference?: ReferenceOverride;
  excluded: string[];
  busy: boolean;
  pendingSourceId?: string;
  onSelect: (locus: LocusKey, segment: SegmentKey, sourceId: string) => void;
  onFile: (locus: LocusKey, segment: SegmentKey, file: File) => void;
  onEditAlleles: (locus: LocusKey, segment: SegmentKey) => void;
}) {
  const input = useRef<HTMLInputElement>(null);
  const databases = databasesForCell(speciesName, locus, segment);
  const value = pendingSourceId ?? (reference?.sourceKind === "upload" ? "upload" : reference?.sourceDatabaseId ?? DEFAULT_DATABASE_ID);
  const sourceFasta = reference?.text ?? builtInFasta;
  const filtered = useMemo(() => filterReferenceFasta(sourceFasta, excluded), [excluded, sourceFasta]);
  const coverage = useMemo(() => segment === "V" || segment === "J" ? annotationCoverage(filtered.fasta, segment) : undefined, [filtered.fasta, segment]);
  return (
    <div className={`composition-cell ${reference?.sourceKind ?? "imgt"} ${excluded.length ? "has-exclusions" : ""} ${busy ? "busy" : ""}`}>
      <select aria-label={`${locus} ${segment} reference source`} aria-busy={busy} value={value} onChange={(event) => onSelect(locus, segment, event.target.value)}>
        <option value={DEFAULT_DATABASE_ID}>IMGT/GENE-DB{filtered.total || reference ? "" : " · no records"}</option>
        {databases.map((database) => <option value={database.id} key={database.id}>{database.name}</option>)}
        {reference?.sourceKind === "upload" && <option value="upload">Loaded file · {reference.name}</option>}
      </select>
      <div className="composition-cell-meta">
        <b>{busy ? `Preparing ${value === DEFAULT_DATABASE_ID ? "IMGT" : databases.find((database) => database.id === value)?.name ?? "reference"}…` : `${filtered.retained.toLocaleString()} active allele${filtered.retained === 1 ? "" : "s"}`}</b>
        {(segment === "V" || segment === "J") && coverage && <em className={coverage.annotated === coverage.total ? "complete" : coverage.annotated ? "partial" : "missing"} title={reference?.report.warnings.join("\n")}>{coverage.annotated.toLocaleString()}/{coverage.total.toLocaleString()} {segment === "V" ? "regions" : "anchors"}</em>}
      </div>
      <input ref={input} className="visually-hidden" type="file" accept=".fa,.fasta,.fna,.fas,.txt,.gz" onChange={(event) => {
        const file = event.target.files?.[0];
        if (file) onFile(locus, segment, file);
        event.target.value = "";
      }} />
      <div className="composition-cell-actions"><button className="cell-alleles" type="button" title={`Include or exclude individual ${locus} ${segment} alleles before assignment`} disabled={busy || !filtered.total} onClick={() => onEditAlleles(locus, segment)}>{excluded.length ? `${excluded.length} excluded` : "Exclude alleles…"}</button><button className="cell-download" type="button" title={`Download the active ${locus} ${segment} FASTA, including metadata and current allele exclusions`} aria-label={`Download ${locus} ${segment} database`} disabled={busy || !filtered.retained} onClick={() => downloadBlob(new Blob([filtered.fasta], { type: "text/plain;charset=utf-8" }), `${inputStem(speciesName)}.${locus}.${segment}.fasta`)}>↓</button><button className="cell-upload" type="button" disabled={busy} onClick={() => input.current?.click()}>{reference?.sourceKind === "upload" ? "Replace FASTA" : "Load FASTA"}</button></div>
    </div>
  );
}

function ReferenceCompositionMatrix({
  species,
  scope,
  references,
  exclusions,
  busyCells,
  pendingSources,
  onSelect,
  onFile,
  onEditAlleles,
}: {
  species: ReferenceSpecies;
  scope: ScopeKey;
  references: ReferenceCellMap;
  exclusions: ReferenceAlleleExclusionMap;
  busyCells: Set<string>;
  pendingSources: Record<string,string>;
  onSelect: (locus: LocusKey, segment: SegmentKey, sourceId: string) => void;
  onFile: (locus: LocusKey, segment: SegmentKey, file: File) => void;
  onEditAlleles: (locus: LocusKey, segment: SegmentKey) => void;
}) {
  const loci = lociForScope(species, scope);
  return (
    <div className="composition-matrix" role="table" aria-label="Germline reference composition">
      <div className="composition-row composition-head" role="row"><span role="columnheader">Locus</span>{SEGMENTS.map((segment) => <span className={`segment-${segment.toLowerCase()}`} role="columnheader" key={segment}>{segment}<small>{segment === "V" ? "variable" : segment === "D" ? "diversity" : segment === "J" ? "joining" : "constant"}</small></span>)}</div>
      {loci.map((locus) => <div className="composition-row" role="row" key={locus}>
        <div className="composition-locus" role="rowheader"><strong>{locus}</strong><small>{LOCUS_LABELS[locus]}</small></div>
        {SEGMENTS.map((segment) => {
          if (!segmentAppliesToLocus(locus, segment)) return <div className="composition-na" role="cell" key={segment}><span>Not used</span></div>;
          const alleles = species.loci[locus]?.[segment] ?? [];
          const key = referenceCellKey(locus, segment);
          return <div role="cell" key={segment}><ReferenceCellControl speciesName={species.name} locus={locus} segment={segment} builtInFasta={allelesToFasta(alleles)} reference={references[key]} excluded={exclusions[key] ?? []} busy={busyCells.has(key)} pendingSourceId={pendingSources[key]} onSelect={onSelect} onFile={onFile} onEditAlleles={onEditAlleles} /></div>;
        })}
      </div>)}
    </div>
  );
}

function CompositionSummary({
  databases,
  hasUploads,
  excludedAlleles,
  busy,
  release,
  onCancel,
}: {
  databases: ReferenceDatabase[];
  hasUploads: boolean;
  excludedAlleles: number;
  busy: boolean;
  release: string;
  onCancel: () => void;
}) {
  const sources = [`IMGT/GENE-DB ${release || "reference pack"}`, ...databases.map((database) => database.name), ...(hasUploads ? ["local FASTA"] : [])];
  return (
    <section className={`database-summary ${databases.length || hasUploads || excludedAlleles ? "alternative" : ""}`} aria-label="Reference composition summary">
      <div>
        <span>Reference composition</span>
        <strong>{sources.join(" + ")}</strong>
        <p>{busy ? "Downloading and validating published germline FASTA in this browser…" : `Apply a database above, then refine any locus/segment independently in the matrix. Unchanged cells retain IMGT.${excludedAlleles ? ` ${excludedAlleles.toLocaleString()} exact allele${excludedAlleles === 1 ? " is" : "s are"} currently removed from the FASTA used for initial assignment.` : ""}`}</p>
      </div>
      {busy && <button className="post-cancel" type="button" onClick={onCancel}>Cancel reference preparation</button>}
      <div className="database-links"><a href="https://doi.org/10.1093/nar/gki010" target="_blank" rel="noreferrer">IMGT/GENE-DB citation ↗</a><a href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="noreferrer">IMGT data · CC BY 4.0 ↗</a>{databases.flatMap((database) => [<a href={database.sourceUrl} target="_blank" rel="noreferrer" key={`${database.id}:source`}>{database.name} source ↗</a>, <a href={database.citationUrl} target="_blank" rel="noreferrer" key={`${database.id}:citation`}>Citation ↗</a>, ...(database.terms ? [<a href={database.terms.url} target="_blank" rel="noreferrer" key={`${database.id}:terms`}>{database.terms.label} ↗</a>] : [])])}</div>
    </section>
  );
}

function formatQueriesPerSecond(value: number | null | undefined): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return "—";
  if (value >= 100) return Math.round(value).toLocaleString();
  return value.toLocaleString(undefined, { minimumFractionDigits: 1, maximumFractionDigits: 1 });
}

function assignmentPhaseLabel(telemetry: AssignmentTelemetry | null): string {
  if (!telemetry || telemetry.phase === "initializing") return "Loading workers";
  if (telemetry.phase === "streaming") return "Streaming input";
  if (telemetry.phase === "draining") return "Draining final batches";
  return "Finalizing index";
}

export function AnalysisProgress({ stage, value, telemetry, onCancel, lockState, hidden }: {
  stage: string;
  value: number;
  telemetry: AssignmentTelemetry | null;
  onCancel: () => void;
  lockState: "unsupported" | "waiting" | "held";
  hidden: boolean;
}) {
  const phases = [
    { label: "Load WASM", threshold: 0.08 },
    { label: "Replicate germlines", threshold: 0.14 },
    { label: "Stream + annotate", threshold: 0.96 },
    { label: "Finalize index", threshold: 1 },
  ];
  const percent = Math.min(100, Math.max(0, Math.round(value * 100)));
  const activeWorkers = telemetry ? `${telemetry.activeWorkers} / ${telemetry.totalWorkers}` : "—";
  const phase = assignmentPhaseLabel(telemetry);
  const phaseDetail = !telemetry || telemetry.phase === "initializing"
    ? "Replicating per-worker germline indexes"
    : telemetry.phase === "draining"
    ? `${telemetry.recordsOutstanding.toLocaleString()} records still outstanding`
    : telemetry.phase === "finalizing"
      ? "All assignment batches committed"
      : `${telemetry.recordsCommitted.toLocaleString()} records committed`;
  const workerDetail = telemetry?.phase === "initializing"
    ? "loading WASM and germline indexes"
    : "currently assigning batches";
  return (
    <section className="progress-stage" aria-live="polite" data-assignment-phase={telemetry?.phase ?? "initializing"}>
      <div className="progress-orbit"><span><b>{percent}</b><small>%</small></span><i style={{ "--progress": `${percent * 3.6}deg` } as React.CSSProperties} /></div>
      <div className="progress-copy"><p className="eyebrow"><span>SwiftIG is running locally</span></p><div className="progress-current-operation"><span>Current operation</span><h2 title={stage}>{stage}</h2></div><p>The page will move to Results automatically when the local AIRR index is ready.</p><dl className="assignment-telemetry" aria-label="Live assignment telemetry"><div><dt>Queries / s</dt><dd data-query-rate={telemetry?.queriesPerSecond ?? "pending"}>{formatQueriesPerSecond(telemetry?.queriesPerSecond)}</dd><small>2 s EMA · AIRR committed</small></div><div><dt>Active workers</dt><dd data-active-workers={telemetry?.activeWorkers ?? 0}>{activeWorkers}</dd><small>{workerDetail}</small></div><div><dt>Assignment state</dt><dd>{phase}</dd><small>{phaseDetail}</small></div></dl><p className={`background-run-state ${lockState}`}>{lockState === "held" ? `${hidden ? "Background tab · " : ""}active-run Web Lock held; Chrome Energy Saver freeze is inhibited while this analysis is running.` : lockState === "waiting" ? "Waiting for another Swig analysis tab to release the active-run lock." : "This browser does not expose Web Locks; keep the tab active for long runs."}</p><div className="main-progress"><i style={{ width: `${percent}%` }} /></div><ol>{phases.map((phase, index) => {
        const previous = index ? phases[index - 1].threshold : 0;
        const state = value >= phase.threshold ? "complete" : value >= previous ? "active" : "";
        return <li className={state} key={phase.label}><span>{state === "complete" ? "✓" : index + 1}</span>{phase.label}</li>;
      })}</ol><button className="cancel-button" type="button" onClick={onCancel}>Cancel analysis</button></div>
    </section>
  );
}

interface AlternativeHit {
  call: string;
  score: number;
  identity: number;
  queryStart: number;
  queryEnd: number;
  germlineStart: number;
  germlineEnd: number;
}

function splitCalls(value: string): string[] {
  return value.split(",").map((call) => call.trim()).filter(Boolean);
}

function parseAlternatives(value: string): AlternativeHit[] {
  if (!value) return [];
  return value.split(";").flatMap((entry) => {
    const [call, score, identity, queryStart, queryEnd, germlineStart, germlineEnd] = entry.split("|");
    const parsed = { call, score: Number(score), identity: Number(identity), queryStart: Number(queryStart), queryEnd: Number(queryEnd), germlineStart: Number(germlineStart), germlineEnd: Number(germlineEnd) };
    return call && Number.isFinite(parsed.score) ? [parsed] : [];
  });
}

function formatAlignmentSupport(value: string): string {
  if (!value) return "";
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) return value;
  if (parsed === 0) return "E < 5e−324";
  return `E ${parsed.toExponential(2).replace("e-", "e−")}`;
}

function ResultDetail({ row, onClose }: { row: AirrRow; onClose: () => void }) {
  const [mode, setMode] = useState<"nt" | "aa">("nt");
  const sequenceLength = row.sequence?.length ?? 0;
  const displayedSegments = row.d2_call ? ["v", "d", "d2", "j", "c"] : ["v", "d", "j", "c"];
  const mapSegments = displayedSegments.flatMap((segment) => {
    const start = Number(row[`${segment}_sequence_start`]);
    const end = Number(row[`${segment}_sequence_end`]);
    if (!row[`${segment}_call`] || !start || !end || !sequenceLength) return [];
    return [{ segment, call: row[`${segment}_call`], left: (start - 1) / sequenceLength * 100, width: Math.max(1.5, (end - start + 1) / sequenceLength * 100) }];
  });
  const isotype = row.isotype || inferIsotype(
    row.c_call,
    row.c_sequence_alignment,
    row.c_identity ? Number(row.c_identity) : null,
    row.c_support ? Number(row.c_support) : null,
  );
  const uncertainSegments = ["v", "d", "j", "c"].flatMap((segment) => {
    const selected = splitCalls(row[`${segment}_call`] || "");
    const alternatives = parseAlternatives(row[`${segment}_alternatives`] || "");
    return selected.length > 1 || alternatives.length ? [{ segment, selected, alternatives }] : [];
  });

  return (
    <article className="detail-panel">
      <header className="detail-header">
        <div><span className="section-kicker">Selected rearrangement</span><h2>{row.sequence_id}</h2><div className="detail-tags"><span>{row.locus || "unassigned"}</span>{row.sample_id&&<span>sample · {row.sample_id}</span>}{row.subject_id&&<span>donor · {row.subject_id}</span>}{row.swig_timepoint&&<span>{row.swig_timepoint}</span>}{row.swig_compartment&&<span>{row.swig_compartment}</span>}<span className={row.productive === "T" ? "good" : "warn"}>{row.productive === "T" ? "Productive" : "Non-productive"}</span>{row.d2_call && <span>VDDJ screen-supported</span>}{isotype && <span className="isotype-tag">{isotype}</span>}{row.rev_comp === "T" && <span>Reverse complement</span>}</div></div>
        <button className="close-detail" type="button" onClick={onClose}>Close <span>×</span></button>
      </header>

      <section className="rearrangement-overview">
        <div className="overview-heading"><h3>Rearrangement map</h3><span>{sequenceLength.toLocaleString()} nt</span></div>
        <div className="rearrangement-track">
          <i className="baseline" />
          {mapSegments.map((segment) => <span key={segment.segment} className={`map-segment ${segment.segment}`} style={{ left: `${segment.left}%`, width: `${segment.width}%` }} title={segment.call}>{segment.segment.toUpperCase()}</span>)}
        </div>
        <div className="call-grid">
          {displayedSegments.map((segment) => {
            const selected = splitCalls(row[`${segment}_call`] || "");
            const alternatives = parseAlternatives(row[`${segment}_alternatives`] || "");
            const support = formatAlignmentSupport(row[`${segment}_support`] || "");
            return <div key={segment}><span>{segment.toUpperCase()} call</span><strong>{selected[0] || "—"}</strong><small>{row[`${segment}_identity`] ? `${(Number(row[`${segment}_identity`]) * 100).toFixed(1)}% identity` : "not assigned"}{support ? ` · ${support}` : ""}{selected.length > 1 ? ` · ${selected.length} co-optimal` : alternatives.length ? ` · ${alternatives.length} near-tied` : ""}</small></div>;
          })}
        </div>
      </section>

      {row.d2_call && <section className="uncertainty-panel double-d-evidence">
        <div className="uncertainty-heading"><div><span className="section-kicker">Opt-in double-D evidence</span><h3>{row.d_call} → {row.d2_call}</h3></div><p>Two ordered, non-overlapping D seeds passed the configured evidence rule and could not be explained by the allowed single-D pseudo-tandem distance. The standard SwiftIG call was {row.standard_d_call || "unassigned"}.</p></div>
        <div className="post-stat-grid compact"><article><span>V–J search span</span><strong>{row.swig_double_d_vj_span} nt</strong></article><article><span>Exact seed</span><strong>{row.swig_double_d_seed_length} nt / D</strong></article><article><span>Two-D score gain</span><strong>{row.swig_double_d_score_gain}</strong></article><article><span>Single-D Δ-distance</span><strong>{row.swig_double_d_pseudo_distance || "not explainable"}</strong></article></div>
        <div className="junction-panel"><div><span className="section-kicker">D1 → D2 insertion</span><code>{row.np2 || "zero length"}</code></div><dl><div><dt>NP2</dt><dd>{row.np2_length || "0"} nt</dd></div><div><dt>D2 → J / NP3</dt><dd>{row.np3_length || "0"} nt</dd></div><div><dt>Screen mode</dt><dd>{row.swig_double_d_mode?.replace("_", " ")}</dd></div></dl></div>
        {row.swig_double_d_alternatives && <p className="scientific-note"><span>i</span>Near-tied ordered pairs: <code>{row.swig_double_d_alternatives}</code></p>}
      </section>}

      {uncertainSegments.length > 0 && <section className="uncertainty-panel">
        <div className="uncertainty-heading"><div><span className="section-kicker">Call uncertainty</span><h3>Co-optimal and near-tied hits</h3></div><p>Comma-separated calls have the same optimal score. Additional candidates are retained within SwiftIG’s configured uncertainty window. Full alignment strings are stored only for the selected call to bound per-record output size.</p></div>
        <div className="uncertainty-stack">{uncertainSegments.map(({ segment, selected, alternatives }) => <article key={segment}>
          <header><span className={`segment-symbol segment-${segment}`}>{segment.toUpperCase()}</span><div><strong>{selected.length > 1 ? `${selected.length} co-optimal calls` : "Selected + alternate evidence"}</strong><small>Ranked alignment evidence for this segment</small></div></header>
          <div className="cooptimal-calls">{selected.map((call, index) => <span key={call} className={index === 0 ? "primary" : ""}>{call}{index === 0 ? " · reported first" : " · tied"}</span>)}</div>
          {alternatives.length > 0 && <div className="alternate-table"><div className="alternate-head"><span>Alternative</span><span>Score</span><span>Identity</span><span>Query</span><span>Germline</span></div>{alternatives.map((hit) => <div key={`${hit.call}-${hit.queryStart}`}><strong>{hit.call}</strong><span>{hit.score}</span><span>{(hit.identity * 100).toFixed(1)}%</span><span>{hit.queryStart}–{hit.queryEnd}</span><span>{hit.germlineStart}–{hit.germlineEnd}</span></div>)}</div>}
        </article>)}</div>
      </section>}

      <section className="junction-panel">
        <div><span className="section-kicker">Junction evidence</span><h3>{row.junction_aa || "No translated junction"}</h3><code>{row.junction || "—"}</code></div>
        <dl><div><dt>CDR3 AA</dt><dd>{row.cdr3_aa || "—"}</dd></div><div><dt>Length</dt><dd>{row.junction_aa_length ? `${row.junction_aa_length} aa` : "—"}</dd></div><div><dt>In frame</dt><dd>{row.vj_in_frame || "—"}</dd></div><div><dt>Stop codon</dt><dd>{row.stop_codon || "—"}</dd></div></dl>
      </section>

      <AlignmentViewer row={row} mode={mode} onMode={setMode} />

      <details className="raw-fields">
        <summary>All AIRR fields <span>{Object.values(row).filter(Boolean).length} populated</span></summary>
        <div>{Object.entries(row).filter(([, value]) => value !== "").map(([key, value]) => <dl key={key}><dt>{key}</dt><dd>{value}</dd></dl>)}</div>
      </details>
    </article>
  );
}

function DoubleDExplorer({session,onInspect}:{session:ResultSession;onInspect:(ordinal:number)=>void}){
  const [records,setRecords]=useState<DoubleDEvidenceRecord[]>([]);const [loading,setLoading]=useState(true);const [d1,setD1]=useState("");const [d2,setD2]=useState("");const [sequenceId,setSequenceId]=useState("");const [cdr3,setCdr3]=useState("");const [minimumGain,setMinimumGain]=useState(0);const [minimumSpan,setMinimumSpan]=useState(0);const [selected,setSelected]=useState<AirrRow|null>(null);const [selectedOrdinal,setSelectedOrdinal]=useState<number|null>(null);const [mode,setMode]=useState<"nt"|"aa">("nt");const detailRef=useRef<HTMLElement>(null);
  useEffect(()=>{let cancelled=false;setLoading(true);void(async()=>{const evidence=await session.store.doubleDRecords();const indexed=await session.store.indexRecords(evidence.map((record)=>record.ordinal));const byOrdinal=new Map(indexed.map((record)=>[record.ordinal,record]));const hydrated=evidence.map((record)=>{const index=byOrdinal.get(record.ordinal);return {...record,values:{...record.values,cdr3:record.values.cdr3||index?.cdr3||"",cdr3_aa:record.values.cdr3_aa||index?.cdr3Aa||""}};});if(!cancelled){setRecords(hydrated);setLoading(false);}})().catch(()=>{if(!cancelled){setRecords([]);setLoading(false);}});return()=>{cancelled=true;};},[session]);
  const filtered=useMemo(()=>records.filter(record=>{const values=record.values;const normalizedCdr3=cdr3.toUpperCase();return(!d1||String(values.d_call||"").toUpperCase().includes(d1.toUpperCase()))&&(!d2||String(values.d2_call||"").toUpperCase().includes(d2.toUpperCase()))&&(!sequenceId||String(values.sequence_id||"").toLowerCase().includes(sequenceId.toLowerCase()))&&(!normalizedCdr3||(values.cdr3||"").toUpperCase().includes(normalizedCdr3)||(values.cdr3_aa||"").toUpperCase().includes(normalizedCdr3))&&Number(values.swig_double_d_score_gain||0)>=minimumGain&&Number(values.swig_double_d_vj_span||0)>=minimumSpan;}),[records,d1,d2,sequenceId,cdr3,minimumGain,minimumSpan]);
  async function open(record:DoubleDEvidenceRecord){setSelectedOrdinal(record.ordinal);const [detail]=await session.store.detailMany([record.ordinal]);if(!detail)return;setSelected(detail.values);window.requestAnimationFrame(()=>{detailRef.current?.scrollIntoView({behavior:window.matchMedia("(prefers-reduced-motion: reduce)").matches?"auto":"smooth",block:"start"});detailRef.current?.focus({preventScroll:true});});}
  const pairCounts=useMemo(()=>{const map=new Map<string,number>();for(const record of filtered){const key=`${record.values.d_call||"D1 —"} → ${record.values.d2_call||"D2 —"}`;map.set(key,(map.get(key)??0)+1);}return [...map.entries()].sort((a,b)=>b[1]-a[1]).slice(0,12);},[filtered]);
  return <section className="double-d-explorer embedded-double-d-explorer">
        <header className="repertoire-heading"><div><span className="section-kicker">Opt-in VDDJ evidence</span><h2>Double-D call explorer</h2><p>Inspect ordered D1/D2 evidence. Calls remain separate from the standard AIRR annotation path.</p></div><span className="aggregate-badge">{filtered.length.toLocaleString()} / {records.length.toLocaleString()} supported calls</span></header>
        <div className="double-d-filter-grid"><label><span>Sequence ID</span><CommitTextInput value={sequenceId} onCommit={setSequenceId} placeholder="contains…"/></label><label><span>CDR3 nt or AA</span><CommitTextInput value={cdr3} onCommit={setCdr3} placeholder="CARDR / TGTGCC…"/></label><label><span>D1 call</span><CommitTextInput value={d1} onCommit={setD1} placeholder="IGHD…"/></label><label><span>D2 call</span><CommitTextInput value={d2} onCommit={setD2} placeholder="IGHD…"/></label><label><span>Minimum score gain</span><CommitNumberInput min="0" value={minimumGain} onCommit={setMinimumGain}/></label><label><span>Minimum V–J span</span><CommitNumberInput min="0" value={minimumSpan} onCommit={setMinimumSpan}/></label></div>
        {pairCounts.length?<div className="double-d-pair-summary">{pairCounts.map(([pair,count])=><article key={pair}><code>{pair}</code><strong>{count.toLocaleString()}</strong></article>)}</div>:null}
        <div className="lineage-table-wrap"><table><thead><tr><th>Sequence</th><th>Locus</th><th>CDR3</th><th>D1</th><th>D2</th><th>Score gain</th><th>V–J span</th><th>NP2</th><th/></tr></thead><tbody>{filtered.slice(0,1000).map(record=><tr key={record.ordinal} className={selectedOrdinal===record.ordinal?"selected":""} onClick={()=>void open(record)}><td><strong>{record.values.sequence_id||`#${record.ordinal+1}`}</strong></td><td>{record.values.locus||"—"}</td><td className="cdr3-table-cell"><code title={record.values.cdr3_aa||record.values.cdr3}>{record.values.cdr3_aa||record.values.cdr3||"—"}</code>{record.values.cdr3_aa&&record.values.cdr3?<small title={record.values.cdr3}>{record.values.cdr3}</small>:null}</td><td><code>{record.values.d_call||"—"}</code></td><td><code>{record.values.d2_call||"—"}</code></td><td>{record.values.swig_double_d_score_gain||"—"}</td><td>{record.values.swig_double_d_vj_span||"—"}</td><td><code>{record.values.np2||"—"}</code></td><td><button type="button">Open alignment →</button></td></tr>)}</tbody></table>{loading?<div className="detail-loading">Loading sparse double-D evidence…</div>:!filtered.length?<div className="empty-results"><span>∅</span><h3>No supported call matches these filters.</h3></div>:null}</div>
      {selected?<section ref={detailRef} className="double-d-alignment-detail" tabIndex={-1}><header><div><span className="section-kicker">Selected VDDJ record</span><h3>{selected.sequence_id}</h3><p>{selected.v_call} → {selected.d_call} → {selected.d2_call} → {selected.j_call}</p><code className="selected-cdr3">CDR3 {selected.cdr3_aa||selected.cdr3||"—"}{selected.cdr3_aa&&selected.cdr3?` · ${selected.cdr3}`:""}</code></div><div className="result-actions"><button type="button" onClick={()=>onInspect(selectedOrdinal??0)}>Open complete AIRR record</button><button type="button" onClick={()=>setSelected(null)}>Close alignment</button></div></header><div className="post-stat-grid compact"><article><span>Two-D score gain</span><strong>{selected.swig_double_d_score_gain||"—"}</strong></article><article><span>V–J search span</span><strong>{selected.swig_double_d_vj_span||"—"} nt</strong></article><article><span>D1→D2 insertion</span><strong>{selected.np2_length||0} nt</strong></article><article><span>D2→J insertion</span><strong>{selected.np3_length||0} nt</strong></article></div><AlignmentViewer row={selected} mode={mode} onMode={setMode}/></section>:null}
  </section>;
}

type ResultsView = "repertoire" | "sequences" | "post";
type ResultsUtilityPanel = "study" | "palette" | null;

function ResultsPage({ session, onNewAnalysis }: { session: ResultSession; onNewAnalysis: () => void }) {
  const initialView: ResultsView = session.webMode === "lineage" ? "post" : session.webMode === "vdj" ? "sequences" : session.pipeline.enabled ? "post" : session.total <= 3 ? "sequences" : "repertoire";
  const [view, setView] = useState<ResultsView>(initialView);
  const [openedViews, setOpenedViews] = useState<Set<ResultsView>>(() => new Set([
    initialView,
    ...(session.postAnalysis || session.pipeline.enabled ? ["post" as const] : []),
  ]));
  const [filters, setFilters] = useState<ResultFilters>({ ...EMPTY_FILTERS });
  const [page, setPage] = useState(0);
  const [results, setResults] = useState<ResultPage>({ rows: [], hasMore: false, totalMatches: session.total, scanned: 0 });
  const [searching, setSearching] = useState(true);
  const [scanCount, setScanCount] = useState(0);
  const [selected, setSelected] = useState<AirrIndexRecord | null>(null);
  const [detail, setDetail] = useState<AirrRow | null>(null);
  const [sequenceWorkspace,setSequenceWorkspace]=useState<"records"|"detail"|"double-d">(session.total<=3?"detail":"records");
  const [downloading, setDownloading] = useState(false);
  const [doubleDDownloading, setDoubleDDownloading] = useState(false);
  const [downloadError, setDownloadError] = useState("");
  const [downloadFormat,setDownloadFormat]=useState<TableExportFormat>("tsv");
  const [downloadCompression,setDownloadCompression]=useState<OutputCompression>("none");
  const [savingSession,setSavingSession]=useState(false);
  const [cliConfigExporting,setCliConfigExporting]=useState(false);
  const [projectSaveStatus,setProjectSaveStatus]=useState(session.project?session.projectStatus||"Project state current":"");
  const [datasets,setDatasets]=useState<DatasetManifestEntry[]>(()=>session.datasets.map((dataset)=>({...dataset})));
  const [metadataDraft,setMetadataDraft]=useState<DatasetManifestEntry[]>(()=>session.datasets.map((dataset)=>({...dataset})));
  const [facets,setFacets]=useState<ResultFacets>(session.facets);
  const [metadataRevision,setMetadataRevision]=useState(0);
  const [metadataBusy,setMetadataBusy]=useState(false);
  const [metadataProgress,setMetadataProgress]=useState({processed:0,total:session.total});
  const [metadataStatus,setMetadataStatus]=useState("");
  const [sampleColors,setSampleColors]=useState<SampleColorMap>(()=>createSampleColorMap(session.datasets,session.sampleColors));
  const [utilityPanel,setUtilityPanel]=useState<ResultsUtilityPanel>(null);
  const postSessionRef=useRef<PostAnalysisSessionHandle|null>(null);
  const metadataAbortRef=useRef<AbortController|null>(null);
  const exportAbortRef=useRef<AbortController|null>(null);
  const sessionSaveAbortRef=useRef<AbortController|null>(null);
  const cliConfigExportAbortRef=useRef<AbortController|null>(null);
  const searchAbortRef=useRef<AbortController|null>(null);
  const autoOpened = useRef(false);
  const detailRef = useRef<HTMLElement>(null);
  const scrollToDetail = useRef(false);
  const resultsTabsRef=useRef<HTMLElement>(null);
  const viewScrollPositionsRef=useRef<Partial<Record<ResultsView,number>>>({});
  const pendingViewScrollRef=useRef<number|null>(null);
  const projectSaveTimerRef=useRef<number|null>(null);
  const projectSaveChainRef=useRef<Promise<void>>(Promise.resolve());
  const projectMetadataReadyRef=useRef(false);

  function tabDocumentTop():number{
    let top=0;
    let element:HTMLElement|null=resultsTabsRef.current;
    while(element){top+=element.offsetTop;element=element.offsetParent as HTMLElement|null;}
    return top;
  }

  function activateView(next:ResultsView,restoreScroll=true){
    if(next===view)return;
    viewScrollPositionsRef.current[view]=window.scrollY;
    setOpenedViews((current)=>{if(current.has(next))return current;const updated=new Set(current);updated.add(next);return updated;});
    pendingViewScrollRef.current=restoreScroll?(viewScrollPositionsRef.current[next]??Math.max(0,tabDocumentTop()-8)):null;
    setView(next);
  }

  useLayoutEffect(()=>{
    const target=pendingViewScrollRef.current;
    if(target===null)return;
    pendingViewScrollRef.current=null;
    let second=0;
    const first=window.requestAnimationFrame(()=>{second=window.requestAnimationFrame(()=>window.scrollTo({top:target,left:0,behavior:"auto"}));});
    return()=>{window.cancelAnimationFrame(first);if(second)window.cancelAnimationFrame(second);};
  },[view]);

  useEffect(()=>{
    if(!utilityPanel)return;
    const close=(event:KeyboardEvent)=>{if(event.key==="Escape")setUtilityPanel(null);};
    document.addEventListener("keydown",close);
    return()=>document.removeEventListener("keydown",close);
  },[utilityPanel]);

  async function persistProjectState(reason:string){
    if(!session.project)return;
    setProjectSaveStatus("Writing project checkpoint…");setDownloadError("");
    const postAnalysis=postSessionRef.current?await postSessionRef.current.snapshot():session.postAnalysis??{workingStages:[]};
    const artifact=await sessionArtifact(session,datasets,sampleColors,postAnalysis);
    await saveProjectCheckpoint(session.project,artifact,await encodeSession(artifact),reason);
    session.postAnalysis=postAnalysis;
    setProjectSaveStatus(`Project checkpoint ${activeProjectRun(session.project)?.checkpointCount??""} written`);
  }

  function saveProjectNow(reason:string){
    if(!session.project)return;
    if(projectSaveTimerRef.current!==null){window.clearTimeout(projectSaveTimerRef.current);projectSaveTimerRef.current=null;}
    projectSaveChainRef.current=projectSaveChainRef.current.then(()=>persistProjectState(reason)).catch((error)=>{
      setProjectSaveStatus("Project checkpoint failed");
      setDownloadError(error instanceof Error?error.message:String(error));
    });
  }

  function scheduleProjectSave(reason:string,delay=900){
    if(!session.project)return;
    if(projectSaveTimerRef.current!==null)window.clearTimeout(projectSaveTimerRef.current);
    setProjectSaveStatus("Project state changed · checkpoint pending");
    projectSaveTimerRef.current=window.setTimeout(()=>{
      projectSaveTimerRef.current=null;
      saveProjectNow(reason);
    },delay);
  }

  useEffect(()=>()=>{if(projectSaveTimerRef.current!==null)window.clearTimeout(projectSaveTimerRef.current);metadataAbortRef.current?.abort();exportAbortRef.current?.abort();sessionSaveAbortRef.current?.abort();cliConfigExportAbortRef.current?.abort();searchAbortRef.current?.abort();},[]);

  useEffect(()=>{
    if(!session.project)return;
    if(!projectMetadataReadyRef.current){projectMetadataReadyRef.current=true;return;}
    scheduleProjectSave("study_metadata_or_palette_changed");
  },[datasets,sampleColors]);

  useEffect(() => {
    if (session.webMode === "lineage") {
      setSearching(false);
      return;
    }
    let cancelled = false;
    const controller = new AbortController();
    searchAbortRef.current = controller;
    setSearching(true);
    setScanCount(0);
    const delay = filters.sequenceId.trim() || filters.cdr3.trim() ? 180 : 0;
    const timer = window.setTimeout(() => {
      session.store.page(filters, page * PAGE_SIZE, PAGE_SIZE, (scanned) => {
        if (!cancelled) setScanCount(scanned);
      }, controller.signal).then((next) => {
        if (cancelled) return;
        setResults(next);
        setSearching(false);
      }).catch(() => { if (!cancelled) setSearching(false); });
    }, delay);
    return () => {
      cancelled = true;
      controller.abort();
      if(searchAbortRef.current===controller)searchAbortRef.current=null;
      window.clearTimeout(timer);
    };
  }, [filters, page, session]);

  useEffect(() => {
    if (session.total <= 3 && !autoOpened.current && !selected && results.rows[0]) {
      autoOpened.current = true;
      setSelected(results.rows[0]);
    }
  }, [results.rows, selected, session.total]);

  useEffect(() => {
    let cancelled = false;
    if (!selected) {
      setDetail(null);
      return;
    }
    session.store.detail(selected).then((row) => { if (!cancelled) setDetail(row); });
    return () => { cancelled = true; };
  }, [selected, session]);

  useEffect(() => {
    if (!selected || !scrollToDetail.current) return;
    scrollToDetail.current = false;
    const frame = window.requestAnimationFrame(() => {
      detailRef.current?.scrollIntoView({
        behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth",
        block: "start",
      });
      detailRef.current?.focus({ preventScroll: true });
    });
    return () => window.cancelAnimationFrame(frame);
  }, [selected]);

  function openRecord(row: AirrIndexRecord) {
    setSequenceWorkspace("detail");
    if (selected?.ordinal === row.ordinal) {
      window.requestAnimationFrame(() => detailRef.current?.scrollIntoView({
        behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth",
        block: "start",
      }));
      return;
    }
    setDetail(null);
    scrollToDetail.current = true;
    setSelected(row);
  }

  async function inspectOrdinal(ordinal: number) {
    const [record] = await session.store.indexRecords([ordinal]);
    if (!record) return;
    activateView("sequences",false);
    setSequenceWorkspace("detail");
    if(selected?.ordinal===record.ordinal){
      window.requestAnimationFrame(()=>detailRef.current?.scrollIntoView({behavior:window.matchMedia("(prefers-reduced-motion: reduce)").matches?"auto":"smooth",block:"start"}));
      return;
    }
    setDetail(null);
    scrollToDetail.current = true;
    setSelected(record);
  }

  function updateFilter<K extends keyof ResultFilters>(key: K, value: ResultFilters[K]) {
    setFilters((current) => ({ ...current, [key]: value }));
    setPage(0);
    setSelected(null);
  }

  function updateMetadataDraft(datasetId:string,field:"sampleId"|"subjectId"|"cohort"|"timepoint"|"compartment",value:string){
    setMetadataDraft((current)=>current.map((dataset)=>dataset.datasetId===datasetId?{...dataset,[field]:value}:dataset));
    setMetadataStatus("");
  }

  async function applyStudyMetadata(){
    const controller=new AbortController();metadataAbortRef.current=controller;
    setMetadataBusy(true);setMetadataStatus("");setDownloadError("");setMetadataProgress({processed:0,total:session.total});
    try{
      const updated=metadataDraft.map((dataset)=>({...dataset,sampleId:dataset.sampleId.trim(),subjectId:dataset.subjectId.trim(),cohort:dataset.cohort.trim(),timepoint:dataset.timepoint.trim(),compartment:(dataset.compartment??"").trim()}));
      const nextFacets=await session.store.updateStudyMetadata(updated,(processed,total)=>setMetadataProgress({processed,total}),controller.signal);
      setDatasets(updated);setMetadataDraft(updated.map((dataset)=>({...dataset})));setFacets(nextFacets);
      setSampleColors((current)=>createSampleColorMap(updated,current));
      setFilters({...EMPTY_FILTERS});setPage(0);setSelected(null);setDetail(null);
      session.datasets=updated;session.facets=nextFacets;session.postAnalysis=undefined;
      postSessionRef.current=null;setMetadataRevision((value)=>value+1);
      setMetadataStatus("Study metadata applied. V(D)J calls are unchanged; downstream filters, collapse, chimera, lineage, SHM, queries, alignments, and trees were cleared and must be rerun.");
    }catch(error){if(isAbortError(error))setMetadataStatus("Metadata re-indexing cancelled. The previous index and all downstream state were retained.");else setDownloadError(error instanceof Error?error.message:String(error));}
    finally{if(metadataAbortRef.current===controller)metadataAbortRef.current=null;setMetadataBusy(false);}
  }

  async function downloadAll() {
    const controller=new AbortController();exportAbortRef.current=controller;
    setDownloading(true);
    setDownloadError("");
    try {
      const name = compressedName(formattedOutputName(session.inputName,downloadFormat),downloadCompression);
      const picker = savePicker();
      if (picker) {
        const handle = await picker.call(window, {
          suggestedName: name,
          types: [{ description: downloadCompression==="gzip"?"Gzip-compressed AIRR rearrangement table":"AIRR rearrangement table", accept: { [compressedMime("text/plain",downloadCompression)]: [compressedExtension(tableExtension(downloadFormat),downloadCompression)] } }],
        });
        const writable = outputWriter(await handle.createWritable(),downloadCompression);
        try {
          await session.store.writeAirrFormat(downloadFormat,(part) => writable.write(part),undefined,controller.signal);
          await writable.close();
        } catch (error) {
          await writable.abort(error);
          throw error;
        }
      } else if(downloadFormat==="tsv"&&!session.streamedDirectly&&!session.store.hasStudyMetadataOverrides) {
        try {
          await registerDownloadWorker();
          const anchor = document.createElement("a");
          anchor.href = session.store.streamingDownloadUrl(import.meta.env.BASE_URL, name, downloadCompression);
          anchor.download = name;
          anchor.click();
        } catch {
          if(downloadCompression==="none")downloadBlob(await session.store.airrBlob(controller.signal), name);
          else throw new Error("This browser could not create the streaming gzip download. Use Chrome/Edge on HTTPS so Swig can stream through Save As.");
        }
      } else if(downloadFormat==="tsv"&&downloadCompression==="none") {
        if (session.streamedDirectly || session.store.hasStudyMetadataOverrides) {
          downloadBlob(await session.store.airrBlob(controller.signal), name);
          return;
        }
      } else {
        if(session.outputBytes>256*1024*1024)throw new Error("This converted or gzip export is too large for a memory-backed browser download. Use Chrome/Edge on HTTPS so Swig can stream it through Save As.");
        const parts:BlobPart[]=[];
        const writable=outputWriter({write:async(part)=>{parts.push(part instanceof Uint8Array?part.slice().buffer:part);}},downloadCompression);
        try{await session.store.writeAirrFormat(downloadFormat,writable.write,undefined,controller.signal);await writable.close();}
        catch(error){await writable.abort(error);throw error;}
        downloadBlob(new Blob(parts,{type:compressedMime("text/plain;charset=utf-8",downloadCompression)}),name);
      }
    } catch (error) {
      if (!(error instanceof DOMException && error.name === "AbortError")) {
        setDownloadError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      if(exportAbortRef.current===controller)exportAbortRef.current=null;
      setDownloading(false);
    }
  }

  async function downloadDoubleD() {
    const controller=new AbortController();exportAbortRef.current=controller;
    setDoubleDDownloading(true);
    setDownloadError("");
    try {
      const name = compressedName(formattedDoubleDOutputName(session.inputName,downloadFormat),downloadCompression);
      const picker = savePicker();
      if (picker) {
        const handle = await picker.call(window, {
          suggestedName: name,
          types: [{ description: downloadCompression==="gzip"?"Gzip-compressed Swig double-D evidence table":"Swig double-D evidence table", accept: { [compressedMime("text/plain",downloadCompression)]: [compressedExtension(tableExtension(downloadFormat),downloadCompression)] } }],
        });
        const writable = outputWriter(await handle.createWritable(),downloadCompression);
        try {
          await session.store.writeDoubleDFormat(downloadFormat,(part) => writable.write(part),controller.signal);
          await writable.close();
        } catch (error) {
          await writable.abort(error);
          throw error;
        }
      } else {
        const parts: BlobPart[] = [];
        const writable=outputWriter({write:async(part)=>{parts.push(part instanceof Uint8Array?part.slice().buffer:part);}},downloadCompression);
        try{await session.store.writeDoubleDFormat(downloadFormat,writable.write,controller.signal);await writable.close();}
        catch(error){await writable.abort(error);throw error;}
        downloadBlob(new Blob(parts, { type: compressedMime("text/plain;charset=utf-8",downloadCompression) }), name);
      }
    } catch (error) {
      if (!(error instanceof DOMException && error.name === "AbortError")) {
        setDownloadError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      if(exportAbortRef.current===controller)exportAbortRef.current=null;
      setDoubleDDownloading(false);
    }
  }

  async function saveAnalysisSession(){
    const controller=new AbortController();sessionSaveAbortRef.current=controller;
    setSavingSession(true);setDownloadError("");
    try{
      const postAnalysis=postSessionRef.current?await postSessionRef.current.snapshot(controller.signal):session.postAnalysis??{workingStages:[]};
      const artifact=await sessionArtifact(session,datasets,sampleColors,postAnalysis,controller.signal);
      downloadBlob(await encodeSession(artifact,controller.signal),sessionBaseName(session.inputName));
    }catch(error){if(!isAbortError(error))setDownloadError(error instanceof Error?error.message:String(error));}finally{if(sessionSaveAbortRef.current===controller)sessionSaveAbortRef.current=null;setSavingSession(false);}
  }

  async function exportCliConfig(){
    if(cliConfigExporting){cliConfigExportAbortRef.current?.abort();return;}
    const controller=new AbortController();cliConfigExportAbortRef.current=controller;setCliConfigExporting(true);
    setDownloadError("");
    try{
      const post=postSessionRef.current?await postSessionRef.current.snapshot(controller.signal):session.postAnalysis??{workingStages:[]};
      const state=cliStateFromPostAnalysis(session.pipeline,post);
      const config=cliConfigFromBrowser({
        studyName:session.inputName,studyDesign:session.studyDesign,datasets,references:session.references,species:session.species,scope:session.scope,workers:session.workers,callingProfile:session.callingProfile,minimumIdentity:session.minimumIdentity,minimumConstantPrefixIdentity:session.minimumConstantPrefixIdentity,strand:session.strand,fastqFilter:session.fastqFilter,subsample:{enabled:session.subsampleSize!==null,size:session.subsampleSize??10_000,seed:session.subsampleSeed??1},doubleD:session.doubleD,pipeline:state.pipeline,
        collapseOptions:state.collapseOptions,chimeraOptions:state.chimeraOptions,selectionOptions:state.selectionOptions,alleleOptions:state.alleleOptions,lineageOptions:state.lineageOptions,missingAlleleOptions:state.missingAlleleOptions,
      });
      config.output.airrCompression=downloadCompression;
      if(config.pipeline.chimera.enabled&&config.pipeline.chimera.msaSource==="selected"&&!config.pipeline.chimera.uploadedMsa.trim()){
        const msa=await runKalignTask(session.references[config.pipeline.chimera.segment],controller.signal);
        prepareReferenceMsa(msa);config.pipeline.chimera.msaSource="upload";config.pipeline.chimera.uploadedMsa=msa;config.pipeline.chimera.uploadedMsaName=`${config.pipeline.chimera.segment.toLowerCase()}-reference-msa.fasta`;
      }
      if(controller.signal.aborted)throw new DOMException("Config export cancelled.","AbortError");
      downloadBlob(new Blob([`${JSON.stringify(config,null,2)}\n`],{type:"application/json"}),`${session.inputName.replace(/\.[^.]+$/,"")||"swig"}.swig-cli.json`);
    }catch(error){if(!isAbortError(error))setDownloadError(error instanceof Error?error.message:String(error));}
    finally{if(cliConfigExportAbortRef.current===controller)cliConfigExportAbortRef.current=null;setCliConfigExporting(false);}
  }

  const callSuggestions=useMemo(()=>({
    vCall:callFacetItems(facets.vCalls),dCall:callFacetItems(facets.dCalls),
    jCall:callFacetItems(facets.jCalls),cCall:callFacetItems(facets.cCalls),
  }),[facets.cCalls,facets.dCalls,facets.jCalls,facets.vCalls]);
  const filtered = Object.entries(filters).some(([key, value]) => key.endsWith("IncludeAmbiguous") ? false : key.startsWith("min") ? Number(value) > 0 : Boolean(value));
  const pageStart = page * PAGE_SIZE + (results.rows.length ? 1 : 0);
  const pageEnd = page * PAGE_SIZE + results.rows.length;

  function resultsSidebarTools(){
    return <div className="results-rail-tools">
      <div className="results-rail-summary">
        <span>Current analysis</span>
        <strong>{session.total.toLocaleString()} records</strong>
        <small>{friendlySpecies(session.species)} · {LOCUS_LABELS[session.scope]}</small>
      </div>
      {session.webMode === "repertoire" && <div className="results-rail-section">
        <span className="results-rail-section-label">Study</span>
        <button type="button" className={utilityPanel==="study"?"active":""} onClick={()=>setUtilityPanel("study")}><strong>Study design</strong><small>Samples, donors and timepoints</small></button>
        <button type="button" className={utilityPanel==="palette"?"active":""} onClick={()=>setUtilityPanel("palette")}><strong>Sample colors</strong><small>Shared figure palette</small></button>
      </div>}
      <div className="results-rail-section results-rail-export">
        <span className="results-rail-section-label">Files</span>
        <label><span>Table format</span><select value={downloadFormat} onChange={(event)=>setDownloadFormat(event.target.value as TableExportFormat)}><option value="tsv">AIRR TSV</option><option value="csv">CSV</option><option value="jsonl">JSON Lines</option></select></label>
        <label><span>Compression</span><select value={downloadCompression} onChange={(event)=>setDownloadCompression(event.target.value as OutputCompression)}><option value="none">None</option><option value="gzip">Gzip (.gz)</option></select></label>
        <button className={downloading?"post-cancel":"primary"} type="button" onClick={()=>downloading?exportAbortRef.current?.abort():void downloadAll()} disabled={doubleDDownloading}>{downloading?"Cancel results export":`Download ${downloadFormat.toUpperCase()}`}</button>
        {session.doubleDCount>0&&<button className={doubleDDownloading?"post-cancel":undefined} type="button" onClick={()=>doubleDDownloading?exportAbortRef.current?.abort():void downloadDoubleD()} disabled={downloading}>{doubleDDownloading?"Cancel evidence export":`Double-D evidence · ${session.doubleDCount.toLocaleString()}`}</button>}
        <button className={savingSession?"post-cancel":undefined} type="button" onClick={()=>savingSession?sessionSaveAbortRef.current?.abort():void saveAnalysisSession()}>{savingSession?"Cancel session save":"Save portable session"}</button>
        {session.webMode === "repertoire" && <button className={cliConfigExporting?"post-cancel":undefined} type="button" onClick={()=>void exportCliConfig()}><strong>{cliConfigExporting?"Cancel CLI config export":"Export CLI config"}</strong><small>{cliConfigExporting?"Stop MSA preparation and discard the export":"Inputs, QC, references, methods, and donor map"}</small></button>}
        {session.project&&<button type="button" onClick={()=>saveProjectNow("manual_checkpoint")}>Save project checkpoint</button>}
      </div>
      {session.project&&<small className="results-rail-project">{projectSaveStatus}</small>}
      {downloadError&&<small className="results-rail-error" role="alert">{downloadError}</small>}
      <button className="results-new-analysis" type="button" onClick={onNewAnalysis}>New analysis</button>
    </div>;
  }

  const repertoireOverview=<section className="repertoire-run-summary" aria-label="Run summary">
    <header><div><span className="section-kicker">Run summary</span><h2>{session.total.toLocaleString()} rearrangements</h2><p>{session.restored?"Saved state restored":`Completed in ${session.seconds.toFixed(2)} s · ${Math.round(session.total/Math.max(session.seconds,0.001)).toLocaleString()} reads/s`}</p></div><strong className="run-state">Complete</strong></header>
    <div className="result-summary app-result-summary">
      <article><span>V + J assigned</span><strong>{session.summary.assigned.toLocaleString()}</strong><small>{percentage(session.summary.assigned,session.total)} of input</small></article>
      <article><span>Productive</span><strong>{session.summary.productive.toLocaleString()}</strong><small>{percentage(session.summary.productive,session.total)} of input</small></article>
      <article><span>CDR3 called</span><strong>{session.summary.withCdr3.toLocaleString()}</strong><small>{percentage(session.summary.withCdr3,session.total)} of input</small></article>
      <article><span>Loci observed</span><strong>{facets.loci.length}</strong><small>{facets.loci.map((item)=>item.value).join(" · ")||"none"}</small></article>
      {session.doubleD.mode!=="off"&&<article><span>Supported Double-D</span><strong>{session.doubleDCount.toLocaleString()}</strong><small>opt-in evidence screen</small></article>}
    </div>
    <details className="run-technical-details"><summary>Run settings</summary><div><span>{datasets.length.toLocaleString()} dataset{datasets.length===1?"":"s"}</span><span>{callingProfileLabel(session.callingProfile)} calling</span><span>{session.workers} WASM worker{session.workers===1?"":"s"}</span><span>{bytes(session.outputBytes)} AIRR</span>{session.minimumConstantPrefixIdentity!==null&&<span>covered C-prefix identity ≥ {Math.round(session.minimumConstantPrefixIdentity*100)}%</span>}{session.subsampleSize&&<span>subsampled from {session.inputTotal.toLocaleString()} input records · seed {session.subsampleSeed}</span>}{fastqQualityResultText(session)&&<span>{fastqQualityResultText(session).replace(/^ · /,"")}</span>}{session.doubleD.mode!=="off"&&<span>Double-D: {session.doubleD.mode==="all"?"all eligible junctions":`V–J span ≥ ${session.doubleD.minimumVjSpan} nt`}</span>}</div></details>
    {session.project&&<div className="repertoire-project-status"><strong>{session.project.root.name}</strong><span>{projectSaveStatus}</span></div>}
  </section>;

  const focusedOverview = session.webMode === "repertoire" ? null : <section className="focused-results-overview">
    <div><span className="section-kicker">{session.webMode === "lineage" ? "Single lineage ready" : "VDJ annotation complete"}</span><h1>{session.total.toLocaleString()} sequence{session.total === 1 ? "" : "s"}</h1><p>{session.summary.assigned.toLocaleString()} V+J assigned · {session.summary.withCdr3.toLocaleString()} CDR3 called · {session.seconds.toFixed(2)} s · {callingProfileLabel(session.callingProfile)}</p></div>
    <div className="focused-results-actions"><label><span>Download format</span><select value={downloadFormat} onChange={(event)=>setDownloadFormat(event.target.value as TableExportFormat)}><option value="tsv">AIRR TSV</option><option value="csv">CSV</option><option value="jsonl">JSON Lines</option></select></label><label><span>Compression</span><select value={downloadCompression} onChange={(event)=>setDownloadCompression(event.target.value as OutputCompression)}><option value="none">None</option><option value="gzip">Gzip (.gz)</option></select></label><button className={downloading?"post-cancel":"post-primary"} type="button" onClick={()=>downloading?exportAbortRef.current?.abort():void downloadAll()}>{downloading?"Cancel export":`Download ${downloadFormat.toUpperCase()}${downloadCompression==="gzip"?".GZ":""}`}</button><button className={savingSession?"post-cancel":undefined} type="button" onClick={()=>savingSession?sessionSaveAbortRef.current?.abort():void saveAnalysisSession()}>{savingSession?"Cancel session save":"Save session"}</button><button type="button" onClick={onNewAnalysis}>New analysis</button></div>
    {downloadError && <p className="run-error" role="alert">{downloadError}</p>}
  </section>;

  return (
    <main className={`results-page results-application-page results-mode-${session.webMode}`}>
      {focusedOverview}
      <section className="results-hero" hidden aria-hidden="true">
        <WorkflowStepper active={3} />
        <div className="results-title"><div><p className="eyebrow"><span>{session.restored?"Saved state restored":`Analysis complete · ${session.seconds.toFixed(2)} s · ${Math.round(session.total / Math.max(session.seconds, 0.001)).toLocaleString()} reads/s`}</span></p><h1>{session.total.toLocaleString()} analyzed<br /><em>rearrangements.</em></h1><p>{datasets.length.toLocaleString()} dataset{datasets.length===1?"":"s"} · {friendlySpecies(session.species)} · {LOCUS_LABELS[session.scope]} · {callingProfileLabel(session.callingProfile)} calling · {session.workers} WASM worker{session.workers===1?"":"s"} · {bytes(session.outputBytes)} AIRR{session.subsampleSize?` · exact random sample per dataset from ${session.inputTotal.toLocaleString()} input records (base seed ${session.subsampleSeed})`:""}{fastqQualityResultText(session)}{session.doubleD.mode!=="off"?` · double-D screen ${session.doubleD.mode==="all"?"all eligible junctions":`V–J spans ≥ ${session.doubleD.minimumVjSpan} nt`}`:""}</p>{session.project&&<span className="project-state-line"><b>Project directory · {session.project.root.name}</b><small>{projectSaveStatus}</small></span>}</div><div className="results-actions"><label className="compact-export-format"><span>Table format</span><select value={downloadFormat} onChange={(event)=>setDownloadFormat(event.target.value as TableExportFormat)}><option value="tsv">AIRR TSV</option><option value="csv">CSV</option><option value="jsonl">JSON Lines</option></select></label><button className={downloading?"post-cancel":"download-primary"} type="button" onClick={()=>downloading?exportAbortRef.current?.abort():void downloadAll()} disabled={doubleDDownloading}>{downloading?"Cancel results export":<>Download {downloadFormat.toUpperCase()}<span>↓</span></>}</button>{session.doubleDCount>0&&<button className={doubleDDownloading?"post-cancel":undefined} type="button" onClick={()=>doubleDDownloading?exportAbortRef.current?.abort():void downloadDoubleD()} disabled={downloading}>{doubleDDownloading?"Cancel evidence export":`Double-D evidence (${session.doubleDCount.toLocaleString()})`}</button>}{session.project&&<button type="button" onClick={()=>saveProjectNow("manual_checkpoint")}>Save project checkpoint</button>}<button className={savingSession?"post-cancel":undefined} type="button" onClick={()=>savingSession?sessionSaveAbortRef.current?.abort():void saveAnalysisSession()}>{savingSession?"Cancel session save":"Save portable session"}</button><button type="button" onClick={onNewAnalysis}>New analysis</button>{downloadError&&<small role="alert">{downloadError}</small>}</div></div>
        <div className="result-summary">
          <article><span>V + J assigned</span><strong>{session.summary.assigned.toLocaleString()}</strong><small>{percentage(session.summary.assigned, session.total)} of input</small></article>
          <article><span>Productive</span><strong>{session.summary.productive.toLocaleString()}</strong><small>{percentage(session.summary.productive, session.total)} of input</small></article>
          <article><span>CDR3 called</span><strong>{session.summary.withCdr3.toLocaleString()}</strong><small>{percentage(session.summary.withCdr3, session.total)} of input</small></article>
          <article><span>Loci observed</span><strong>{facets.loci.length}</strong><small>{facets.loci.map((item) => item.value).join(" · ") || "none"}</small></article>
          {session.doubleD.mode !== "off" && <article><span>Supported double-D</span><strong>{session.doubleDCount.toLocaleString()}</strong><small>opt-in screen · separate evidence table</small></article>}
        </div>
        <details className="study-metadata-editor"><summary><span>Study design</span><strong>Edit sample, donor, timepoint, and compartment metadata</strong><small>Corrections re-index downstream grouping without rerunning V(D)J assignment.</small></summary><div className="study-metadata-table"><div className="study-metadata-head"><span>Dataset</span><span>Sample</span><span>Donor / subject</span><span>Cohort</span><span>Timepoint</span><span>Compartment / tissue</span></div>{metadataDraft.map((dataset)=><div className="study-metadata-row" key={dataset.datasetId}><span><strong>{dataset.inputName}</strong><small>{dataset.datasetId}</small></span><input aria-label={`${dataset.inputName} sample ID`} value={dataset.sampleId} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"sampleId",event.target.value)}/><input aria-label={`${dataset.inputName} subject ID`} value={dataset.subjectId} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"subjectId",event.target.value)}/><input aria-label={`${dataset.inputName} cohort`} value={dataset.cohort} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"cohort",event.target.value)}/><input aria-label={`${dataset.inputName} timepoint`} value={dataset.timepoint} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"timepoint",event.target.value)}/><input aria-label={`${dataset.inputName} compartment`} value={dataset.compartment??""} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"compartment",event.target.value)}/></div>)}</div><div className="study-metadata-actions">{metadataBusy?<button className="post-cancel" type="button" onClick={()=>metadataAbortRef.current?.abort()}>Cancel re-indexing</button>:<button className="post-primary" type="button" onClick={()=>void applyStudyMetadata()}>Apply metadata + clear downstream state</button>}<button type="button" disabled={metadataBusy} onClick={()=>{setMetadataDraft(datasets.map((dataset)=>({...dataset})));setMetadataStatus("");}}>Discard edits</button>{metadataBusy&&<span>{metadataProgress.processed.toLocaleString()} / {metadataProgress.total.toLocaleString()}</span>}</div>{metadataStatus&&<p className="scientific-note"><span>i</span>{metadataStatus}</p>}</details>
        <details className="sample-palette-editor"><summary><span>Study palette</span><strong>Sample colors</strong><small>One association is reused in sample-level figures and phylogenies.</small></summary><div>{sampleIds(datasets).map((sample)=><label key={sample}><input type="color" value={sampleColor(sample,sampleColors)} onChange={(event)=>setSampleColors((current)=>({...current,[sample]:event.target.value}))}/><span><i style={{background:sampleColor(sample,sampleColors)}}/>{sample}</span></label>)}<button type="button" onClick={()=>setSampleColors(createSampleColorMap(datasets))}>Reset palette</button></div></details>
      </section>

      {session.webMode === "repertoire" && <nav ref={resultsTabsRef} className="results-view-tabs" aria-label="Results view" role="tablist">
        <button id="results-tab-repertoire" role="tab" aria-selected={view === "repertoire"} aria-controls="results-panel-repertoire" className={view === "repertoire" ? "active" : ""} type="button" onClick={() => activateView("repertoire")}><b>01</b><span>Repertoire<small>Composition and figures</small></span></button>
        <button id="results-tab-sequences" role="tab" aria-selected={view === "sequences"} aria-controls="results-panel-sequences" className={view === "sequences" ? "active" : ""} type="button" onClick={() => activateView("sequences")}><b>02</b><span>Sequences<small>Filter and inspect calls</small></span></button>
        <button id="results-tab-post" role="tab" aria-selected={view === "post"} aria-controls="results-panel-post" className={view === "post" ? "active" : ""} type="button" onClick={() => activateView("post")}><b>03</b><span>Post-analysis<small>Lineages, SHM and trees</small></span></button>
      </nav>}

      {session.webMode === "repertoire" && openedViews.has("repertoire") && <div id="results-panel-repertoire" role="tabpanel" aria-labelledby="results-tab-repertoire" className="results-view-panel" hidden={view !== "repertoire"}><RepertoireDashboard key={metadataRevision} store={session.store} loci={facets.loci} inputName={session.inputName} samples={facets.samples} sampleColors={sampleColors} sidebarTools={resultsSidebarTools()} overview={repertoireOverview} /></div>}
      {session.webMode !== "lineage" && openedViews.has("sequences") && <div id="results-panel-sequences" role="tabpanel" aria-labelledby={session.webMode === "repertoire" ? "results-tab-sequences" : undefined} className="results-view-panel" hidden={view !== "sequences"}>
      <div className="sequence-context-workspace contextual-workspace">
        <nav className="context-rail" aria-label="Sequence result panels">
          <div className="context-rail-heading"><span>Sequences</span><small>{session.total.toLocaleString()} AIRR records</small></div>
          <button type="button" className={sequenceWorkspace==="records"?"active":""} onClick={()=>setSequenceWorkspace("records")}><b>01</b><span>{session.webMode === "vdj" ? "Assigned records" : "Records + filters"}<small>{results.totalMatches===null?`${pageEnd.toLocaleString()}+ matches`:`${results.totalMatches.toLocaleString()} matches`}{filtered?" · filtered":""}</small></span></button>
          <button type="button" disabled={!selected} className={sequenceWorkspace==="detail"?"active":""} onClick={()=>setSequenceWorkspace("detail")}><b>02</b><span>Record detail<small>{selected?.sequenceId||"Select a record"}</small></span></button>
          {session.doubleDCount>0&&<button type="button" className={`context-rail-secondary ${sequenceWorkspace==="double-d"?"active":""}`} onClick={()=>setSequenceWorkspace("double-d")}><span>Double-D explorer<small>{session.doubleDCount.toLocaleString()} supported calls</small></span></button>}
          {resultsSidebarTools()}
        </nav>
        <div className="context-main sequence-context-main">

      <section className={`explorer-shell sequence-workspace-${sequenceWorkspace}`}>
        <div className="sequence-records-stack" hidden={sequenceWorkspace!=="records"}>
        {session.webMode === "repertoire" && <aside className="filter-panel compact-filter-panel">
          <div className="filter-heading"><div><span className="section-kicker">Local query</span><h2>Filter records</h2><p>Filters update the table below.</p></div>{filtered && <button type="button" onClick={() => { setFilters({ ...EMPTY_FILTERS }); setPage(0); }}>Clear all</button>}</div>
          <div className="sequence-filter-primary-grid">
            <label className="filter-field"><span>Sequence ID contains</span><CommitTextInput type="search" value={filters.sequenceId} placeholder="e.g. clonotype_104" onCommit={(value) => updateFilter("sequenceId", value)} /></label>
            <label className="filter-field"><span>CDR3 substring <small>nt or AA</small></span><CommitTextInput className="monospace" type="search" value={filters.cdr3} placeholder="CARDR / TGTGCC…" onCommit={(value) => updateFilter("cdr3", value)} /></label>
            {datasets.length>1&&<FacetPicker label="Sample" value={filters.sampleId} items={facets.samples} placeholder="Any sample" onChange={(value)=>updateFilter("sampleId",value)}/>}
            <label className="filter-field"><span>Locus</span><select value={filters.locus} onChange={(event) => updateFilter("locus", event.target.value)}><option value="">Any locus</option>{facets.loci.map((item) => <option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label>
            {(["vCall","jCall"] as const).map((key)=>{const segment=key[0].toUpperCase();const ambiguityKey=`${key}IncludeAmbiguous` as const;return <div className="call-filter" key={key}><FacetPicker label={`${segment} gene or allele`} value={filters[key]} items={callSuggestions[key]} allowCustom placeholder={`Any ${segment} call`} help={`Choose an observed ${segment} gene or allele, or type a value. Search narrows this menu without scanning the dataset.`} onChange={(value)=>updateFilter(key,value)}/><label className="check-filter compact"><input type="checkbox" checked={filters[ambiguityKey]} onChange={(event)=>updateFilter(ambiguityKey,event.target.checked)}/><span>Include multi-call assignments containing this {segment}</span></label></div>})}
            <label className="filter-field"><span>Productivity</span><select value={filters.productive} onChange={(event) => updateFilter("productive", event.target.value)}><option value="">Either</option><option value="T">Productive</option><option value="F">Non-productive</option></select></label>
            <label className="filter-field"><span>Isotype / constant class</span><select value={filters.isotype} onChange={(event)=>updateFilter("isotype",event.target.value)}><option value="">Any isotype</option>{facets.isotypes.map((item)=><option value={item.value} key={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label>
          </div>
          <details className="filter-advanced inline-filter-more"><summary>More filters</summary><div>
            {datasets.length>1&&<section className="filter-section"><h3>Study metadata</h3><div className="study-filter-grid"><FacetPicker label="Dataset" value={filters.datasetId} items={facets.datasets} placeholder="Any dataset" onChange={(value)=>updateFilter("datasetId",value)}/><FacetPicker label="Donor / subject" value={filters.subjectId} items={facets.subjects} placeholder="Any donor" onChange={(value)=>updateFilter("subjectId",value)}/><FacetPicker label="Cohort" value={filters.cohort} items={facets.cohorts} placeholder="Any cohort" onChange={(value)=>updateFilter("cohort",value)}/><FacetPicker label="Timepoint" value={filters.timepoint} items={facets.timepoints} placeholder="Any timepoint" onChange={(value)=>updateFilter("timepoint",value)}/><FacetPicker label="Compartment / tissue" value={filters.compartment} items={facets.compartments} placeholder="Any compartment" onChange={(value)=>updateFilter("compartment",value)}/></div></section>}
            <section className="filter-section"><h3>Additional germline calls</h3><div className="sequence-filter-secondary-grid">{(["dCall","cCall"] as const).map((key)=>{const segment=key[0].toUpperCase();const ambiguityKey=`${key}IncludeAmbiguous` as const;return <div className="call-filter" key={key}><FacetPicker label={`${segment} gene or allele`} value={filters[key]} items={callSuggestions[key]} allowCustom placeholder={`Any ${segment} call`} help={`Choose an observed ${segment} gene or allele, or type a value. Search narrows this menu without scanning the dataset.`} onChange={(value)=>updateFilter(key,value)}/><label className="check-filter compact"><input type="checkbox" checked={filters[ambiguityKey]} onChange={(event)=>updateFilter(ambiguityKey,event.target.checked)}/><span>Include multi-call assignments containing this {segment}</span></label></div>})}</div></section>
            <section className="filter-section identity-qc-section"><h3>Identity, junction and QC</h3>
            <label className="identity-filter"><span>Minimum V identity <b>{Math.round(filters.minVIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minVIdentity} onChange={(event) => updateFilter("minVIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum D identity <b>{Math.round(filters.minDIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minDIdentity} onChange={(event) => updateFilter("minDIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum J identity <b>{Math.round(filters.minJIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minJIdentity} onChange={(event) => updateFilter("minJIdentity", Number(event.target.value))} /></label>
            <label className="identity-filter"><span>Minimum C identity <b>{Math.round(filters.minCIdentity * 100)}%</b></span><input type="range" min="0" max="1" step="0.01" value={filters.minCIdentity} onChange={(event) => updateFilter("minCIdentity", Number(event.target.value))} /></label>
            <div className="length-filters">
              <label><span>CDR3 AA min</span><CommitNumberInput min="0" max="250" value={filters.minCdr3AaLength} blankWhenZero placeholder="Any" onCommit={(value) => updateFilter("minCdr3AaLength", value)} /></label>
              <label><span>CDR3 AA max</span><CommitNumberInput min="0" max="250" value={filters.maxCdr3AaLength} blankWhenZero placeholder="Any" onCommit={(value) => updateFilter("maxCdr3AaLength", value)} /></label>
            </div>
            <div className="qc-filter-grid">
              <label><span>VJ frame</span><select value={filters.vjInFrame} onChange={(event) => updateFilter("vjInFrame", event.target.value)}><option value="">Either</option><option value="T">In frame</option><option value="F">Out of frame</option></select></label>
              <label><span>Stop codon</span><select value={filters.stopCodon} onChange={(event) => updateFilter("stopCodon", event.target.value)}><option value="">Either</option><option value="F">Absent</option><option value="T">Present</option></select></label>
              <label><span>Completeness</span><select value={filters.completeVdj} onChange={(event) => updateFilter("completeVdj", event.target.value)}><option value="">Either</option><option value="T">Complete V(D)J</option><option value="F">Partial</option></select></label>
              <label><span>Orientation</span><select value={filters.revComp} onChange={(event) => updateFilter("revComp", event.target.value)}><option value="">Either</option><option value="F">Forward</option><option value="T">Reverse-comp.</option></select></label>
            </div>
            <label className="check-filter"><input type="checkbox" checked={filters.hasD} onChange={(event) => updateFilter("hasD", event.target.checked)} /><span>Require a D assignment</span></label>
            <label className="check-filter"><input type="checkbox" checked={filters.hasCdr3} onChange={(event) => updateFilter("hasCdr3", event.target.checked)} /><span>Require a CDR3 call</span></label>
            </section>
            {session.doubleDCount>0&&<section className="filter-section niche-filter-section"><h3>Rare VDDJ evidence</h3><label className="check-filter"><input type="checkbox" checked={filters.hasDoubleD} onChange={(event) => updateFilter("hasDoubleD", event.target.checked)} /><span>Only supported Double-D calls</span></label><button type="button" onClick={()=>setSequenceWorkspace("double-d")}>Open Double-D evidence explorer ({session.doubleDCount.toLocaleString()})</button></section>}
            <p className="index-note"><span>i</span> Gene and allele filters use browser-local token indexes. A complete comma-separated call set matches that exact ambiguous assignment; sequence substrings scan indexed candidates on demand.</p>
          </div></details>
        </aside>}

        <div className="result-browser">
          <header className="browser-heading"><div><span className="section-kicker">AIRR records</span><h2>{searching ? "Searching local index…" : results.totalMatches !== null ? `${results.totalMatches.toLocaleString()} matching records` : `${(pageEnd + (results.hasMore ? 1 : 0)).toLocaleString()}+ matching records`}</h2><p>{searching && scanCount ? `${scanCount.toLocaleString()} candidates scanned` : results.rows.length ? `Showing ${pageStart.toLocaleString()}–${pageEnd.toLocaleString()}` : "Adjust filters to broaden the query."}</p></div>{searching&&<button className="post-cancel" type="button" onClick={()=>searchAbortRef.current?.abort()}>Cancel search</button>}<span className="scale-mode">{session.total <= 3 ? "detail mode" : session.total >= 100000 ? "large-run index" : "paged index"}</span></header>
          <div className={`results-table-wrap ${searching ? "loading" : ""}`}>
            <table className="results-table">
              <colgroup><col className="column-sequence"/>{datasets.length>1&&<col className="column-sample"/>}<col className="column-locus"/><col className="column-v"/><col className="column-d"/><col className="column-j"/><col className="column-isotype"/><col className="column-cdr3"/><col className="column-productive"/><col className="column-open"/></colgroup>
              <thead><tr><th>Sequence</th>{datasets.length>1&&<th>Sample</th>}<th>Locus</th><th>V call</th><th>D call</th><th>J call</th><th>Isotype</th><th>CDR3</th><th>Productive</th><th /></tr></thead>
              <tbody>{results.rows.map((row) => <tr className={selected?.ordinal === row.ordinal ? "selected" : ""} key={row.ordinal} tabIndex={0} onClick={() => openRecord(row)} onKeyDown={(event) => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); openRecord(row); } }}>
                <td><strong title={row.sequenceId}>{row.sequenceId}</strong><small>#{(row.ordinal + 1).toLocaleString()}</small></td>
                {datasets.length>1&&<td><strong className="sample-colored-label"><i style={{background:sampleColor(row.sampleId,sampleColors)}}/>{row.sampleId||"—"}</strong><small>{row.timepoint||row.subjectId||""}</small></td>}
                <td><span className="locus-pill">{row.locus || "—"}</span></td>
                <td title={row.vCall}>{row.vCall || <i>—</i>}</td><td title={[row.dCall, row.d2Call].filter(Boolean).join(" → ")}>{row.dCall || <i>—</i>}{row.d2Call && <small className="d2-table-call">→ {row.d2Call}</small>}</td><td title={row.jCall}>{row.jCall || <i>—</i>}</td><td>{row.isotype || <i>—</i>}</td>
                <td className="cdr3-table-cell"><code title={row.cdr3Aa || row.cdr3}>{row.cdr3Aa || row.cdr3 || "—"}</code>{row.cdr3Aa && row.cdr3 ? <small title={row.cdr3}>{row.cdr3}</small> : null}</td>
                <td><span className={`productive-dot ${row.productive === "T" ? "yes" : "no"}`} />{row.productive === "T" ? "Yes" : row.productive === "F" ? "No" : "—"}</td>
                <td><button type="button" aria-label={`Open ${row.sequenceId}`}>→</button></td>
              </tr>)}</tbody>
            </table>
            {!searching && !results.rows.length && <div className="empty-results"><span>∅</span><h3>No records match these filters.</h3><p>Try clearing an allele, identity threshold, or substring.</p></div>}
          </div>
          <footer className="table-pagination"><button type="button" disabled={!page || searching} onClick={() => setPage((value) => Math.max(0, value - 1))}>← Previous</button><span>Page {(page + 1).toLocaleString()}</span><button type="button" disabled={!results.hasMore || searching} onClick={() => setPage((value) => value + 1)}>Next →</button></footer>
        </div>
        </div>
      </section>

      {sequenceWorkspace==="detail"&&selected?<section ref={detailRef} className="detail-shell" tabIndex={-1} aria-label={`Details for ${selected.sequenceId}`}>{detail ? <ResultDetail row={detail} onClose={() => {setSelected(null);setSequenceWorkspace("records");}} /> : <div className="detail-loading">Loading selected AIRR record…</div>}</section>:sequenceWorkspace==="detail"?<div className="method-placeholder"><span>↳</span><h3>Select an AIRR record</h3><p>Open a row from Records to inspect its calls and V(D)J alignment.</p><button type="button" onClick={()=>setSequenceWorkspace("records")}>Browse records</button></div>:null}
      {sequenceWorkspace==="double-d"&&<DoubleDExplorer key={metadataRevision} session={session} onInspect={(ordinal)=>void inspectOrdinal(ordinal)}/>}
        </div>
      </div>
      </div>}
      {(session.webMode === "lineage" || openedViews.has("post")) && <div id="results-panel-post" role="tabpanel" aria-labelledby={session.webMode === "repertoire" ? "results-tab-post" : undefined} className="results-view-panel" hidden={view !== "post"}><PostAnalysisWorkbench key={metadataRevision} store={session.store} references={session.references} scope={session.scope} loci={facets.loci} resultFacets={facets} inputName={session.inputName} workers={session.workers} callingProfile={session.callingProfile} minimumIdentity={session.minimumIdentity} strand={session.strand} datasets={datasets} sampleColors={sampleColors} defaultCollapseScope={session.pipeline.collapse.scope} defaultLineageScope={session.pipeline.lineage.scope} doubleDCount={session.doubleDCount} autoPipeline={metadataRevision===0&&session.webMode === "repertoire"&&session.pipeline.enabled&&!session.postAnalysis?session.pipeline:null} sidebarTools={resultsSidebarTools()} onInspect={(ordinal) => { if (session.webMode !== "lineage") void inspectOrdinal(ordinal); }} onSessionChange={(reason)=>scheduleProjectSave(reason)} sessionHandleRef={postSessionRef} initialSession={session.postAnalysis??null} directLineage={session.webMode === "lineage"} /></div>}
      {utilityPanel&&<div className="results-tool-drawer-backdrop" role="presentation" onMouseDown={(event)=>{if(event.target===event.currentTarget)setUtilityPanel(null);}}><section className="results-tool-drawer" role="dialog" aria-modal="true" aria-labelledby="results-tool-drawer-title">
        <header><div><span className="section-kicker">Analysis settings</span><h2 id="results-tool-drawer-title">{utilityPanel==="study"?"Study design":"Sample colors"}</h2><p>{utilityPanel==="study"?"Edit grouping metadata. Applying changes re-indexes downstream grouping without rerunning V(D)J assignment.":"The same sample-to-color association is used throughout figures and lineage trees."}</p></div><button type="button" onClick={()=>setUtilityPanel(null)} aria-label="Close analysis settings">×</button></header>
        {utilityPanel==="study"?<div className="results-drawer-study"><div className="study-metadata-table"><div className="study-metadata-head"><span>Dataset</span><span>Sample</span><span>Donor / subject</span><span>Cohort</span><span>Timepoint</span><span>Compartment / tissue</span></div>{metadataDraft.map((dataset)=><div className="study-metadata-row" key={dataset.datasetId}><span><strong>{dataset.inputName}</strong><small>{dataset.datasetId}</small></span><input aria-label={`${dataset.inputName} sample ID`} value={dataset.sampleId} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"sampleId",event.target.value)}/><input aria-label={`${dataset.inputName} subject ID`} value={dataset.subjectId} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"subjectId",event.target.value)}/><input aria-label={`${dataset.inputName} cohort`} value={dataset.cohort} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"cohort",event.target.value)}/><input aria-label={`${dataset.inputName} timepoint`} value={dataset.timepoint} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"timepoint",event.target.value)}/><input aria-label={`${dataset.inputName} compartment`} value={dataset.compartment??""} onChange={(event)=>updateMetadataDraft(dataset.datasetId,"compartment",event.target.value)}/></div>)}</div><div className="study-metadata-actions">{metadataBusy?<button className="post-cancel" type="button" onClick={()=>metadataAbortRef.current?.abort()}>Cancel re-indexing</button>:<button className="post-primary" type="button" onClick={()=>void applyStudyMetadata()}>Apply metadata + clear downstream state</button>}<button type="button" disabled={metadataBusy} onClick={()=>{setMetadataDraft(datasets.map((dataset)=>({...dataset})));setMetadataStatus("");}}>Discard edits</button>{metadataBusy&&<span>{metadataProgress.processed.toLocaleString()} / {metadataProgress.total.toLocaleString()}</span>}</div>{metadataStatus&&<p className="scientific-note"><span>i</span>{metadataStatus}</p>}</div>:<div className="results-drawer-palette"><div>{sampleIds(datasets).map((sample)=><label key={sample}><input type="color" value={sampleColor(sample,sampleColors)} onChange={(event)=>setSampleColors((current)=>({...current,[sample]:event.target.value}))}/><span><i style={{background:sampleColor(sample,sampleColors)}}/>{sample}</span></label>)}</div><button type="button" onClick={()=>setSampleColors(createSampleColorMap(datasets))}>Reset palette</button></div>}
      </section></div>}
    </main>
  );
}

export default function SwigApp() {
  const [page, setPage] = useState<AppPage>("home");
  const [webMode, setWebMode] = useState<WebAnalysisMode>("repertoire");
  const [pack, setPack] = useState<ReferencePack | null>(null);
  const [packError, setPackError] = useState("");
  const [speciesName, setSpeciesName] = useState("Homo sapiens");
  const [scope, setScope] = useState<ScopeKey>("BCR");
  const [inputSource, setInputSource] = useState<InputSource>("upload");
  const [fileInputs, setFileInputs] = useState<DirectoryDatasetInput[]>([]);
  const [pendingDirectoryInputs,setPendingDirectoryInputs]=useState<{inputs:DirectoryDatasetInput[];flatRoots:string[]}|null>(null);
  const [pendingGzipImport,setPendingGzipImport]=useState<PendingGzipImport|null>(null);
  const [inputInspecting,setInputInspecting]=useState(false);
  const [pasteText, setPasteText] = useState("");
  const [pasteMetadata, setPasteMetadata] = useState<DatasetManifestEntry>({ datasetId: "dataset_paste", inputName: "pasted-sequences.txt", sampleId: "sample_1", subjectId: "subject_1", cohort: "cohort_1", timepoint: "", compartment: "", records: null });
  const [studyName, setStudyName] = useState("swig-study");
  const [studyDesign, setStudyDesign] = useState<StudyDesign>("longitudinal");
  const [pipeline, setPipeline] = useState<PipelinePlan>(() => copyPipeline());
  const [inputError, setInputError] = useState("");
  const [cellReferences, setCellReferences] = useState<ReferenceCellMap>({});
  const [alleleExclusions, setAlleleExclusions] = useState<ReferenceAlleleExclusionMap>({});
  const [editingReferenceCell, setEditingReferenceCell] = useState<{ locus: LocusKey; segment: SegmentKey } | null>(null);
  const [busyCells, setBusyCells] = useState<Set<string>>(new Set());
  const [pendingCellSources,setPendingCellSources]=useState<Record<string,string>>({});
  const [databaseBusy, setDatabaseBusy] = useState(false);
  const [pendingDatabaseId,setPendingDatabaseId]=useState<string|null>(null);
  const [minimumIdentity, setMinimumIdentity] = useState(0.6);
  const [constantPrefixIdentityEnabled, setConstantPrefixIdentityEnabled] = useState(false);
  const [minimumConstantPrefixIdentity, setMinimumConstantPrefixIdentity] = useState(0.97);
  const [callingProfile, setCallingProfile] = useState<CallingProfile>("r_optimized");
  const [strand, setStrand] = useState<0 | 1 | 2>(0);
  const [workerCount, setWorkerCount] = useState(recommendedWorkerCount);
  const [outputStorage, setOutputStorage] = useState<OutputStorageMode>("auto");
  const [subsampleEnabled, setSubsampleEnabled] = useState(false);
  const [subsampleSize, setSubsampleSize] = useState(10_000);
  const [subsampleSeed, setSubsampleSeed] = useState(1);
  const [fastqFilter, setFastqFilter] = useState<FastqQualityFilterOptions>(() => copyFastqQualityFilter());
  const [doubleDMode, setDoubleDMode] = useState<DoubleDScreenMode>("off");
  const [doubleDMinimumSpan, setDoubleDMinimumSpan] = useState(40);
  const [doubleDSeedLength, setDoubleDSeedLength] = useState(11);
  const [doubleDPseudoTrim, setDoubleDPseudoTrim] = useState(5);
  const [doubleDMaximumPseudoMismatches, setDoubleDMaximumPseudoMismatches] = useState(3);
  const [doubleDMinimumScoreGain, setDoubleDMinimumScoreGain] = useState(8);
  const [outputPrompt, setOutputPrompt] = useState(false);
  const [configExporting,setConfigExporting]=useState(false);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState({ stage: "Preparing analysis", value: 0 });
  const [assignmentTelemetry, setAssignmentTelemetry] = useState<AssignmentTelemetry | null>(null);
  const [analysisLockState, setAnalysisLockState] = useState<"unsupported" | "waiting" | "held">("unsupported");
  const [pageHidden, setPageHidden] = useState(typeof document !== "undefined" && document.hidden);
  const [runError, setRunError] = useState("");
  const [session, setSession] = useState<ResultSession | null>(null);
  const [projectWorkspace,setProjectWorkspace]=useState<ProjectWorkspace|null>(null);
  const [projectStatus,setProjectStatus]=useState("");
  const [projectBusy,setProjectBusy]=useState(false);
  const inputRef = useRef<HTMLInputElement>(null);
  const directoryInputRef=useRef<HTMLInputElement>(null);
  const sessionInputRef=useRef<HTMLInputElement>(null);
  const linkedAirrInputRef=useRef<HTMLInputElement>(null);
  const lineageStudyInputRef=useRef<HTMLInputElement>(null);
  const lineageStudyAirrInputRef=useRef<HTMLInputElement>(null);
  const [pendingLoadedSession,setPendingLoadedSession]=useState<SwigSession|null>(null);
  const [pendingLineageStudy,setPendingLineageStudy]=useState<LineageStudyManifest|null>(null);
  const [lineageStudyAirr,setLineageStudyAirr]=useState<File|null>(null);
  const [sessionLoadError,setSessionLoadError]=useState("");
  const [sessionLoadProgress,setSessionLoadProgress]=useState({records:0,total:0,stage:""});
  const [loadingSession,setLoadingSession]=useState(false);
  const sessionLoadAbortRef=useRef<AbortController|null>(null);
  const [leavePrompt,setLeavePrompt]=useState(false);
  const abortRef = useRef<AbortController | null>(null);
  const configExportAbortRef=useRef<AbortController|null>(null);
  const databaseRequestRef = useRef(0);
  const cellRequestRef = useRef<Record<string, number>>({});
  const referenceCacheRef = useRef<Map<string, ReferenceOverride>>(new Map());
  const referenceAbortRef=useRef<Map<string,AbortController>>(new Map());
  const databaseAbortRef=useRef<AbortController|null>(null);
  const preferredDatabaseContextRef = useRef("");
  const nextDatasetIdRef = useRef(1);
  const gzipChoiceResolverRef=useRef<((choice:GzipImportChoice)=>void)|null>(null);
  const inputImportBusyRef=useRef(false);
  const inputInspectAbortRef=useRef<AbortController|null>(null);
  const liveWorkRef=useRef(false);
  const historyGuardInstalledRef=useRef(false);
  const confirmedLeaveRef=useRef(false);

  const hasLiveWork = running
    || Boolean(session)
    || fileInputs.length > 0
    || pasteText.trim().length > 0
    || Boolean(projectWorkspace)
    || Boolean(pendingLoadedSession)
    || Boolean(pendingLineageStudy)
    || Boolean(lineageStudyAirr)
    || loadingSession;

  useEffect(() => {
    loadReferencePack().then(setPack).catch((error) => setPackError(error instanceof Error ? error.message : String(error)));
    if ("serviceWorker" in navigator && window.isSecureContext) void registerDownloadWorker().catch(() => undefined);
  }, []);

  useEffect(()=>{
    const root=document.querySelector<HTMLElement>(".site-shell");
    if(!root)return;
    addFieldHelp(root);
    const observer=new MutationObserver(()=>addFieldHelp(root));
    observer.observe(root,{childList:true,subtree:true});
    return()=>observer.disconnect();
  },[]);

  useEffect(() => {
    const updateVisibility = () => setPageHidden(document.hidden);
    document.addEventListener("visibilitychange", updateVisibility);
    return () => document.removeEventListener("visibilitychange", updateVisibility);
  }, []);

  useEffect(() => { liveWorkRef.current = hasLiveWork; }, [hasLiveWork]);

  useEffect(() => {
    if (!hasLiveWork) return;
    const warnBeforeLeaving = (event: BeforeUnloadEvent) => {
      if (confirmedLeaveRef.current) return;
      event.preventDefault();
      event.returnValue = "";
    };
    window.addEventListener("beforeunload", warnBeforeLeaving);
    return () => window.removeEventListener("beforeunload", warnBeforeLeaving);
  }, [hasLiveWork]);

  useEffect(() => {
    const interceptHistoryDeparture = () => {
      if (confirmedLeaveRef.current || !liveWorkRef.current) return;
      const previous = history.state && typeof history.state === "object" ? history.state : {};
      history.pushState({ ...previous, swigLeaveGuard: true }, "", location.href);
      historyGuardInstalledRef.current = true;
      setLeavePrompt(true);
    };
    window.addEventListener("popstate", interceptHistoryDeparture);
    return () => window.removeEventListener("popstate", interceptHistoryDeparture);
  }, []);

  useEffect(() => {
    if (hasLiveWork && !historyGuardInstalledRef.current) {
      const previous = history.state && typeof history.state === "object" ? history.state : {};
      history.pushState({ ...previous, swigLeaveGuard: true }, "", location.href);
      historyGuardInstalledRef.current = true;
    } else if (!hasLiveWork && historyGuardInstalledRef.current) {
      historyGuardInstalledRef.current = false;
      history.back();
    }
  }, [hasLiveWork]);

  function confirmPageDeparture() {
    confirmedLeaveRef.current = true;
    setLeavePrompt(false);
    history.go(-2);
    window.setTimeout(() => { confirmedLeaveRef.current = false; }, 1200);
  }

  useEffect(() => {
    return () => {
      abortRef.current?.abort();
      configExportAbortRef.current?.abort();
      if (session) void session.store.clear();
    };
  }, [session]);

  const speciesList = useMemo(() => favoriteSpecies(pack?.species ?? []), [pack]);
  const species = useMemo(() => speciesList.find((candidate) => candidate.name === speciesName) ?? speciesList[0], [speciesList, speciesName]);
  const scopes = useMemo(() => species ? availableScopes(species) : [], [species]);
  const activeScope = scopes.includes(scope) ? scope : scopes[0] ?? "BCR";
  const referenceOverrides = useMemo(() => species ? composeReferenceOverrides(
    lociForScope(species, activeScope),
    cellReferences,
    (locus, segment) => allelesToFasta(species.loci[locus]?.[segment] ?? []),
    alleleExclusions,
  ) : {}, [activeScope, alleleExclusions, cellReferences, species]);
  const compiled = useMemo(() => species ? compileReferences(species, activeScope, referenceOverrides) : null, [activeScope, referenceOverrides, species]);
  const databaseOptions = useMemo(() => databaseOptionsFor(species?.name ?? "", pack?.release ?? "", activeScope), [activeScope, pack?.release, species?.name]);
  const activeReferenceEntries = useMemo(() => {
    if (!species) return [] as Array<[string, ReferenceOverride]>;
    const activeLoci = new Set(lociForScope(species, activeScope));
    return Object.entries(cellReferences).filter(([key]) => activeLoci.has(key.split(":", 1)[0] as LocusKey));
  }, [activeScope, cellReferences, species]);
  const activeAlleleExclusions = useMemo(() => {
    if (!species) return {} as ReferenceAlleleExclusionMap;
    const activeLoci = new Set(lociForScope(species, activeScope));
    return Object.fromEntries(Object.entries(alleleExclusions).filter(([key, names]) => names.length && activeLoci.has(key.split(":", 1)[0] as LocusKey)));
  }, [activeScope, alleleExclusions, species]);
  const activeAlleleExclusionCount = useMemo(() => Object.values(activeAlleleExclusions).reduce((total, names) => total + names.length, 0), [activeAlleleExclusions]);
  const usedDatabaseIds = useMemo(() => new Set(activeReferenceEntries.flatMap(([, reference]) => reference.sourceDatabaseId ? [reference.sourceDatabaseId] : [])), [activeReferenceEntries]);
  const usedDatabases = useMemo(() => databaseOptions.flatMap((option) => option.database && usedDatabaseIds.has(option.database.id) ? [option.database] : []), [databaseOptions, usedDatabaseIds]);
  const hasUploadedReferences = activeReferenceEntries.some(([, reference]) => reference.sourceKind === "upload");
  const compositionMode = useMemo(() => {
    if (!activeReferenceEntries.length) return DEFAULT_DATABASE_ID;
    const sourceIds = new Set(activeReferenceEntries.map(([, reference]) => reference.sourceKind === "database" ? reference.sourceDatabaseId : "upload"));
    if (sourceIds.size !== 1 || sourceIds.has("upload")) return "mixed";
    const sourceId=[...sourceIds][0];
    const database=databaseOptions.find((option)=>option.id===sourceId)?.database;
    if(!database)return "mixed";
    const targetKeys=new Set(collectionsForDatabase(database,activeScope).flatMap((collection)=>SEGMENTS.filter((segment)=>collection.segments[segment]).map((segment)=>referenceCellKey(collection.locus,segment))));
    const activeKeys=new Set(activeReferenceEntries.map(([key])=>key));
    return targetKeys.size===activeKeys.size&&[...targetKeys].every((key)=>activeKeys.has(key))?sourceId||"mixed":"mixed";
  }, [activeReferenceEntries,activeScope,databaseOptions]);
  const pendingDatabaseForScope = pendingDatabaseId && databaseOptions.some((option)=>option.id===pendingDatabaseId) ? pendingDatabaseId : null;
  const displayedDatabaseId = pendingDatabaseForScope ?? compositionMode;
  const pendingDatabase = pendingDatabaseForScope && pendingDatabaseForScope !== DEFAULT_DATABASE_ID
    ? databaseOptions.find((option) => option.id === pendingDatabaseForScope)?.database
    : undefined;
  const summaryDatabases = pendingDatabase ? [pendingDatabase] : usedDatabases;
  const databaseLabel = `${activeReferenceEntries.length ? ["IMGT", ...usedDatabases.map((database) => database.name), ...(hasUploadedReferences ? ["local FASTA"] : [])].join(" + ") : `IMGT/GENE-DB ${pack?.release ?? "reference pack"}`}${activeAlleleExclusionCount ? ` · ${activeAlleleExclusionCount.toLocaleString()} allele${activeAlleleExclusionCount === 1 ? "" : "s"} excluded` : ""}`;
  const editingReferenceKey = editingReferenceCell ? referenceCellKey(editingReferenceCell.locus, editingReferenceCell.segment) : "";
  const editingReferenceOverride = editingReferenceKey ? cellReferences[editingReferenceKey] : undefined;
  const editingReferenceFasta = editingReferenceCell && species
    ? editingReferenceOverride?.text ?? allelesToFasta(species.loci[editingReferenceCell.locus]?.[editingReferenceCell.segment] ?? [])
    : "";
  const editingReferenceSource = editingReferenceCell
    ? editingReferenceOverride?.name ?? `IMGT/GENE-DB ${pack?.release ?? "reference pack"}`
    : "";
  const receptor = activeScope.startsWith("IG") || activeScope === "BCR" ? "BCR" : "TCR";
  const receptorScopes = scopes.filter((value) => receptor === "BCR" ? value === "BCR" || value.startsWith("IG") : value === "TCR" || value.startsWith("TR"));
  const pasteInput = useMemo(() => {
    if (!pasteText.trim()) return null;
    try { return inspectText("pasted-sequences.txt", pasteText); } catch { return null; }
  }, [pasteText]);
  const activeDatasets = useMemo<DatasetInput[]>(() => {
    if (inputSource === "upload") return fileInputs;
    return pasteInput ? [{ ...pasteInput, ...pasteMetadata, inputName: pasteInput.name, records: pasteInput.count }] : [];
  }, [fileInputs, inputSource, pasteInput, pasteMetadata]);
  const activeInput = activeDatasets[0] ?? null;
  const pipelineSelectionFacets=useMemo(()=>({
    datasets:activeDatasets.map((dataset)=>({value:dataset.datasetId,count:dataset.records??undefined})),
    samples:uniqueFacetItems(activeDatasets.map((dataset)=>dataset.sampleId)),
    subjects:uniqueFacetItems(activeDatasets.map((dataset)=>dataset.subjectId)),
    cohorts:uniqueFacetItems(activeDatasets.map((dataset)=>dataset.cohort)),
    timepoints:uniqueFacetItems(activeDatasets.map((dataset)=>dataset.timepoint)),
    compartments:uniqueFacetItems(activeDatasets.map((dataset)=>dataset.compartment??"")),
    loci:(compiled?.loci??[]).map((value)=>({value})),
  }),[activeDatasets,compiled?.loci]);
  const activeInputName = activeDatasets.length > 1 ? (studyName.trim() || "swig-study") : activeInput?.name ?? "swig";
  const knownInputCount = activeDatasets.every((input) => input.count !== null)
    ? activeDatasets.reduce((sum, input) => sum + (input.count ?? 0), 0)
    : null;
  const fastqDatasetCount = activeDatasets.filter((input) => input.formatCode === 2).length;
  const nonFastqDatasetCount = activeDatasets.length - fastqDatasetCount;

  function selectSpecies(nextName:string){
    const next=speciesList.find((item)=>item.name===nextName);
    if(!next||next.name===species?.name)return;
    const nextScopes=availableScopes(next);
    const nextScope=nextScopes.includes(activeScope)?activeScope:nextScopes[0]??"BCR";
    const preferred=preferredDatabaseIdFor(next.name,nextScope);
    databaseRequestRef.current+=1;
    cellRequestRef.current={};
    preferredDatabaseContextRef.current="";
    setRunError("");
    setSpeciesName(next.name);
    setScope(nextScope);
    setCellReferences({});
    setAlleleExclusions({});
    setEditingReferenceCell(null);
    setBusyCells(new Set());
    setPendingCellSources({});
    setDatabaseBusy(false);
    setPendingDatabaseId(preferred===DEFAULT_DATABASE_ID?null:preferred);
    referenceCacheRef.current.clear();
  }

  function selectReceptor(next:"BCR"|"TCR"){
    if(next===receptor)return;
    const preferred=preferredDatabaseIdFor(species?.name??"",next);
    resetReferenceContext();
    preferredDatabaseContextRef.current="";
    setPendingDatabaseId(preferred===DEFAULT_DATABASE_ID?null:preferred);
    setScope(next);
  }

  useEffect(() => {
    if (!pack || !species) return;
    const receptorScope: ScopeKey = activeScope === "BCR" || activeScope.startsWith("IG") ? "BCR" : "TCR";
    const context = `${species.name}|${receptorScope}`;
    if (preferredDatabaseContextRef.current === context) return;
    preferredDatabaseContextRef.current = context;
    const preferred = preferredDatabaseIdFor(species.name, receptorScope);
    if (preferred !== DEFAULT_DATABASE_ID) {
      setPendingDatabaseId(preferred);
      void applyReferenceDatabase(preferred);
    } else setPendingDatabaseId(null);
  }, [activeScope, pack, species]);

  function navigate(next: AppPage) {
    if (next === "results" && !session) next = "analyze";
    if(next==="lineage-study"&&(!pendingLineageStudy||!lineageStudyAirr))next="analyze";
    setPage(next);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function chooseWebAnalysisMode(next: WebAnalysisMode) {
    if (next !== webMode) {
      cancelInputInspection();
      setFileInputs([]);
      setPasteText("");
      setInputSource("upload");
      setInputError("");
      setRunError("");
      setPendingDirectoryInputs(null);
      setProjectWorkspace(null);
      setProjectStatus("");
    }
    setWebMode(next);
    if (next !== "repertoire") {
      setPipeline((current) => ({ ...current, enabled: false }));
      setSubsampleEnabled(false);
      setFastqFilter(copyFastqQualityFilter());
    }
    navigate("analyze");
  }

  async function acceptSessionFile(file:File){
    setSessionLoadError("");
    try{const loaded=await decodeSession(file);setPendingLoadedSession(loaded);setSessionLoadProgress({records:0,total:loaded.linkedAirr.records,stage:"Session manifest read; select the linked AIRR table."});}
    catch(error){setPendingLoadedSession(null);setSessionLoadError(error instanceof Error?error.message:String(error));}
  }

  async function acceptLineageStudyFile(file:File){
    setSessionLoadError("");
    try{
      const manifest=await decodeLineageStudy(file);
      setPendingLineageStudy(manifest);setLineageStudyAirr(null);
    }catch(error){setPendingLineageStudy(null);setSessionLoadError(error instanceof Error?error.message:String(error));}
  }

  function restoreLineageStudyAirr(file:File){
    const manifest=pendingLineageStudy;if(!manifest)return;
    if(file.size!==manifest.linkedAirr.size){setSessionLoadError(`This AIRR file has ${file.size.toLocaleString()} bytes; the lineage-study manifest expects ${manifest.linkedAirr.size.toLocaleString()}. Choose the uncompressed CLI lineage AIRR file.`);return;}
    setSessionLoadError("");setLineageStudyAirr(file);setPage("lineage-study");window.scrollTo({top:0});
  }

  async function restoreSavedSession(saved:SwigSession,file:File,restoredProject?:ProjectWorkspace){
    let store:AirrResultStore|undefined;
    const controller=new AbortController();sessionLoadAbortRef.current=controller;
    setLoadingSession(true);setSessionLoadError("");setSessionLoadProgress({records:0,total:saved.linkedAirr.records,stage:"Streaming linked AIRR table into the local index"});
    try{
      const restoredDatasets=(saved.analysis.datasets?.length?saved.analysis.datasets:[{datasetId:"legacy",inputName:saved.analysis.inputName,sampleId:"sample_1",subjectId:"subject_1",cohort:"",timepoint:"",compartment:"",records:saved.linkedAirr.records}]).map((dataset)=>({...dataset,compartment:dataset.compartment??""}));
      store=new AirrResultStore();store.configureStudyMetadataForImport(restoredDatasets);let records=0;
      for await(const batch of streamSequenceBatches({source:file,format:3,batchSize:2000,signal:controller.signal,onProgress:(value)=>setSessionLoadProgress({records:value.recordsRead,total:saved.linkedAirr.records,stage:"Reading and indexing linked AIRR records"})})){
        const newline=batch.text.indexOf("\n");const header=batch.text.slice(0,newline).replace(/\r$/,"");const body=batch.text.slice(newline+1);
        if(!header.includes("\t"))throw new Error("Session loading requires the original AIRR TSV or TSV.gz file, not a comma-separated conversion.");
        await store.appendBatch(header,body);records+=batch.count;setSessionLoadProgress({records,total:saved.linkedAirr.records,stage:"Building browser-local AIRR indexes"});
      }
      if(controller.signal.aborted)throw new DOMException("Session restoration was cancelled.","AbortError");
      await store.finalize();
      if(controller.signal.aborted)throw new DOMException("Session restoration was cancelled.","AbortError");
      const mismatches=linkedAirrMatches(saved,{name:file.name,size:file.size,lastModified:file.lastModified,records:store.count,headers:[...store.airrHeaders],fingerprint:store.fingerprint});
      if(mismatches.length)throw new Error(`This is not the AIRR table linked by the session: ${mismatches.join("; ")}.`);
      if(saved.doubleD.length)await store.importDoubleDRecords(saved.doubleD);
      if(controller.signal.aborted)throw new DOMException("Session restoration was cancelled.","AbortError");
      const dd={mode:"off",minimumVjSpan:40,seedLength:11,pseudoTrim:5,maximumPseudoMismatches:3,minimumScoreGain:8,...saved.analysis.doubleD} as DoubleDScreenOptions;
      setSessionLoadProgress({records:store.count,total:store.count,stage:"Local AIRR indexes ready"});
      const restoredWebMode=saved.analysis.webMode??"repertoire";
      setWebMode(restoredWebMode);
      setSession({id:Date.now(),webMode:restoredWebMode,store,total:store.count,seconds:0,inputName:saved.analysis.inputName,datasets:restoredDatasets,studyDesign:saved.analysis.studyDesign??"longitudinal",pipeline:copyPipeline(saved.analysis.pipeline),species:saved.analysis.species,scope:saved.analysis.scope,facets:store.facets(),summary:store.summary,workers:saved.analysis.workers,outputBytes:store.outputBytes,streamedDirectly:false,inputTotal:store.count,subsampleSize:saved.analysis.subsample?.enabled?saved.analysis.subsample.size:null,subsampleSeed:saved.analysis.subsample?.enabled?saved.analysis.subsample.seed:null,fastqFilter:copyFastqQualityFilter(saved.analysis.fastqFilter),fastqFilterStats:saved.analysis.fastqFilterStats??emptyFastqQualityFilterStats(Boolean(saved.analysis.fastqFilter?.enabled),false),references:saved.analysis.references,referenceExclusions:Object.fromEntries(Object.entries(saved.analysis.referenceExclusions??{}).map(([key,names])=>[key,[...names]])),callingProfile:callingProfileOrDefault(saved.analysis.callingProfile),minimumIdentity:saved.analysis.minimumIdentity,minimumConstantPrefixIdentity:saved.analysis.minimumConstantPrefixIdentity??null,strand:saved.analysis.strand,doubleD:dd,doubleDCount:store.doubleDCount,sampleColors:createSampleColorMap(restoredDatasets,saved.analysis.sampleColors),postAnalysis:saved.postAnalysis,restored:true,project:restoredProject});
      if(restoredProject){setProjectWorkspace(restoredProject);const run=activeProjectRun(restoredProject);if(run)await appendProjectLog(restoredProject,run,"project_opened",{records:store.count});setProjectStatus(`Restored ${run?.id??"active run"}`);}
      setPendingLoadedSession(null);setSessionLoadProgress({records:store.count,total:store.count,stage:"Session restored"});setPage("results");window.scrollTo({top:0});
    }catch(error){if(store)await store.clear();if(isAbortError(error))setSessionLoadProgress({records:0,total:saved.linkedAirr.records,stage:"Session restoration cancelled; no state was changed"});else setSessionLoadError(error instanceof Error?error.message:String(error));}finally{if(sessionLoadAbortRef.current===controller)sessionLoadAbortRef.current=null;setLoadingSession(false);}
  }

  function cancelSessionRestore(){
    sessionLoadAbortRef.current?.abort();
    setSessionLoadProgress((current)=>({...current,stage:"Cancelling session restoration…"}));
  }

  async function restoreLinkedAirr(file:File){
    const saved=pendingLoadedSession;if(!saved)return;
    await restoreSavedSession(saved,file);
  }

  async function chooseProjectDirectory(){
    setProjectBusy(true);setProjectStatus("");setSessionLoadError("");
    try{
      const root=await selectProjectDirectory();
      const attached=await attachProjectDirectory(root,APP_VERSION);
      setProjectWorkspace(attached.workspace);
      if(attached.existing&&activeProjectRun(attached.workspace)){
        setProjectStatus("Reading project manifest and active run…");
        const files=await loadActiveProjectFiles(attached.workspace);
        await restoreSavedSession(await decodeSession(files.sessionFile),files.airrFile,attached.workspace);
      }else{
        setProjectStatus(`Project directory ready · ${root.name}`);
        setPage("analyze");window.scrollTo({top:0});
      }
    }catch(error){
      if(!(error instanceof DOMException&&error.name==="AbortError")){
        const message=error instanceof Error?error.message:String(error);
        setProjectStatus(message);setSessionLoadError(message);
      }
    }finally{setProjectBusy(false);}
  }

  function requestGzipImportChoice(pending:PendingGzipImport):Promise<GzipImportChoice>{
    return new Promise((resolve)=>{
      gzipChoiceResolverRef.current=resolve;
      setPendingGzipImport(pending);
    });
  }

  function resolveGzipImportChoice(choice:GzipImportChoice){
    const resolve=gzipChoiceResolverRef.current;
    gzipChoiceResolverRef.current=null;
    setPendingGzipImport(null);
    resolve?.(choice);
  }

  async function acceptInputCandidates(candidates: InputFileCandidate[]) {
    if (!candidates.length) return;
    if (webMode !== "repertoire" && candidates.length !== 1) {
      setInputError("This focused workflow accepts one input file. Choose a single FASTA, FASTQ, or AIRR file.");
      return;
    }
    if(inputImportBusyRef.current){setInputError("Finish the current gzip-member choice before adding more files.");return;}
    inputImportBusyRef.current=true;
    const controller=new AbortController();inputInspectAbortRef.current=controller;
    const startingDatasetId=nextDatasetIdRef.current;
    setInputInspecting(true);
    setInputError("");
    setRunError("");
    const donorPlan=inferDirectoryDonors(candidates);
    const donorByPath=new Map(donorPlan.suggestions.map((entry)=>[entry.relativePath,entry.donor]));
    const accepted: DirectoryDatasetInput[] = [];
    const failures: string[] = [];
    let ignored=0;
    let promoteFocusedImportToRepertoire=false;
    try{
      for (const candidate of candidates) {
        if(controller.signal.aborted)throw new DOMException("Input inspection was cancelled.","AbortError");
        if(candidate.fromDirectory&&(!/\.(fa|fasta|fna|fas|fq|fastq|tsv|csv|txt)(\.gz)?$/i.test(candidate.file.name)||candidate.file.name.startsWith("."))){ignored+=1;continue;}
        try {
          let inspectedInputs:InputData[]=[];
          const namedFormat=formatFromName(candidate.file.name);
          if(candidate.file.name.toLowerCase().endsWith(".gz")){
            if(!namedFormat)throw new Error("Name gzip inputs with .fasta.gz, .fastq.gz, or .tsv.gz.");
            const members=await inspectGzipMembers(candidate.file,namedFormat.formatCode,controller.signal);
            if(members.length>1){
              const source=gzipMemberSource(candidate.file,members);
              const merged:InputData={name:candidate.file.name,source,count:null,size:candidate.file.size,...namedFormat,inputPath:candidate.relativePath,inputFormat:cliInputFormat(namedFormat.formatCode),gzipMemberCount:members.length,gzipMembersMerged:true};
              if(members.every((member)=>member.startsAtRecordBoundary)){
                const seen=new Set<string>();
                const suggestedNames=members.map((member,index)=>{
                  let name=suggestedGzipMemberName(candidate.file.name,member,index);
                  if(seen.has(name.toLowerCase())){
                    const gzip=name.match(/\.gz$/i)?.[0]??"";
                    const base=gzip?name.slice(0,-gzip.length):name;
                    name=`${base}.member-${index+1}${gzip}`;
                  }
                  seen.add(name.toLowerCase());
                  return name;
                });
                const choice=await requestGzipImportChoice({inputName:candidate.relativePath,format:namedFormat.format,members,suggestedNames});
                if(choice==="cancel")continue;
                if(choice==="separate"&&webMode!=="repertoire")promoteFocusedImportToRepertoire=true;
                inspectedInputs=choice==="separate"
                  ? members.map((member,index)=>{
                    const memberFile=gzipMemberFile(candidate.file,member,suggestedNames[index]);
                    return {name:memberFile.name,source:memberFile,count:null,size:memberFile.size,...namedFormat,inputPath:candidate.relativePath,inputFormat:cliInputFormat(namedFormat.formatCode),gzipRange:{start:member.start,end:member.end},gzipMemberCount:1,gzipMembersMerged:false};
                  })
                  : [merged];
              }else inspectedInputs=[merged];
            }else inspectedInputs=[await inspectFile(candidate.file)];
          }else inspectedInputs=[await inspectFile(candidate.file)];

          const slash=candidate.relativePath.lastIndexOf("/");
          const directoryPrefix=candidate.fromDirectory&&slash>=0?candidate.relativePath.slice(0,slash+1):"";
          for(const inspected of inspectedInputs){
            if(candidate.fromDirectory)inspected.name=`${directoryPrefix}${inspected.name}`;
            inspected.inputPath??=candidate.relativePath;
            inspected.inputFormat??=cliInputFormat(inspected.formatCode);
            const ordinal = nextDatasetIdRef.current++;
            const input:DirectoryDatasetInput={...datasetInput(inspected,ordinal),directoryRoot:candidate.fromDirectory?candidate.rootName:undefined,nestedDirectoryDonor:donorByPath.get(candidate.relativePath)??undefined};
            if(input.nestedDirectoryDonor)input.subjectId=input.nestedDirectoryDonor;
            accepted.push(input);
          }
        } catch (error) {
          if(isAbortError(error))throw error;
          failures.push(`${candidate.relativePath}: ${error instanceof Error ? error.message : String(error)}`);
        }
      }
      const acceptingAsRepertoire=webMode === "repertoire"||promoteFocusedImportToRepertoire;
      if(promoteFocusedImportToRepertoire)setWebMode("repertoire");
      const acceptedFlatRoots=acceptingAsRepertoire ? donorPlan.flatRoots.filter((root)=>accepted.some((input)=>input.directoryRoot===root&&!input.nestedDirectoryDonor)) : [];
      if(accepted.length){
        if(acceptedFlatRoots.length)setPendingDirectoryInputs({inputs:accepted,flatRoots:acceptedFlatRoots});
        else setFileInputs((current)=>acceptingAsRepertoire ? [...current,...accepted] : accepted.slice(0,1));
      }
      if(failures.length||(!accepted.length&&ignored))setInputError(`${failures.slice(0,8).join(" ")}${failures.length>8?` ${failures.length-8} additional file(s) could not be read.`:""}${ignored?` ${ignored.toLocaleString()} non-dataset file(s) were ignored.`:""}`.trim());
    }catch(error){
      if(isAbortError(error)){nextDatasetIdRef.current=startingDatasetId;setInputError("");}
      else setInputError(error instanceof Error?error.message:String(error));
    }finally{
      if(inputInspectAbortRef.current===controller)inputInspectAbortRef.current=null;
      inputImportBusyRef.current=false;
      setInputInspecting(false);
    }
  }

  function cancelInputInspection(){
    inputInspectAbortRef.current?.abort();
    if(gzipChoiceResolverRef.current)resolveGzipImportChoice("cancel");
  }

  async function acceptInputFiles(files: File[]) {
    await acceptInputCandidates(files.map((file)=>({file,relativePath:file.name,rootName:file.name,fromDirectory:false})));
  }

  function commitDirectoryInputs(sameDonor:boolean){
    const pending=pendingDirectoryInputs;if(!pending)return;
    const flatRoots=new Set(pending.flatRoots);
    const resolved=pending.inputs.map((input)=>{
      if(input.nestedDirectoryDonor||!input.directoryRoot||!flatRoots.has(input.directoryRoot)||!sameDonor)return input;
      return {...input,subjectId:donorForFlatRoot(input.directoryRoot)};
    });
    setFileInputs((current)=>[...current,...resolved]);
    setPendingDirectoryInputs(null);
  }

  function updateDataset(datasetId: string, patch: Partial<DatasetManifestEntry>) {
    setFileInputs((current) => current.map((input) => input.datasetId === datasetId ? { ...input, ...patch } : input));
  }

  function applyStudyDesign(next: StudyDesign) {
    setStudyDesign(next);
    const defaults = studyScopeDefaults(next);
    setPipeline((current) => ({
      ...current,
      collapse: { ...current.collapse, scope: defaults.collapse },
      alleleRefinement: { ...current.alleleRefinement, scope: defaults.lineage },
      lineage: { ...current.lineage, scope: defaults.lineage },
    }));
    if (next === "technical") {
      setFileInputs((current) => current.map((input) => ({ ...input, sampleId: "sample_1", subjectId: "subject_1" })));
    } else if (next === "independent" || next === "cohort") {
      setFileInputs((current) => current.map((input, index) => ({ ...input, sampleId: `${inputStem(input.inputName)}_${index + 1}`, subjectId: `subject_${index + 1}` })));
    } else if (next === "longitudinal") {
      // A previous technical-replicate preset may have merged sample IDs. Restore
      // file-level specimens while preserving explicitly grouped donor IDs.
      setFileInputs((current) => current.map((input, index) => ({ ...input, sampleId: `${inputStem(input.inputName)}_${index + 1}` })));
    }
  }

  function updatePipelineAlleleSegment(segment: "V" | "D" | "J", selected: boolean) {
    setPipeline((current) => {
      const segments = selected
        ? [...new Set([...current.alleleRefinement.segments, segment])]
        : current.alleleRefinement.segments.filter((value) => value !== segment);
      return { ...current, alleleRefinement: { ...current.alleleRefinement, segments } };
    });
  }

  async function prepareReferenceOverride(
    segment: SegmentKey,
    text: string,
    name: string,
    sourceKind: "database" | "upload",
    referenceScope: ScopeKey,
    sourceDatabaseId?: string,
    signal?: AbortSignal,
  ): Promise<ReferenceOverride> {
    if (!pack || !species) throw new Error("The baseline reference pack is not ready.");
    const allowedLoci = lociForScope(species, referenceScope);
    if (!allowedLoci.length) throw new Error("Choose a supported locus before adding germlines.");
    const report = await preprocessGermlinesInWorker(
      text,
      segment,
      germlineTemplateTiers(pack, species, referenceScope, segment),
      allowedLoci,
      signal,
    );
    return {
      name,
      text: report.fasta,
      count: report.count,
      size: new Blob([report.fasta]).size,
      report,
      sourceKind,
      sourceDatabaseId,
    };
  }

  function markCellsBusy(keys: string[], busy: boolean) {
    setBusyCells((current) => {
      const next = new Set(current);
      for (const key of keys) busy ? next.add(key) : next.delete(key);
      return next;
    });
  }

  async function databaseCellReference(
    database: ReferenceDatabase,
    locus: LocusKey,
    segment: SegmentKey,
    signal?: AbortSignal,
  ): Promise<ReferenceOverride> {
    if (!species) throw new Error("Choose a species before loading a database.");
    const cacheKey = `${species.name}:${database.id}:${locus}:${segment}`;
    const cached = referenceCacheRef.current.get(cacheKey);
    if (cached) return cached;
    const collection = collectionsForDatabase(database, locus).find((candidate) => candidate.segments[segment]);
    if (!collection) throw new Error(`${database.name} does not provide ${locus} ${segment} records.`);
    const text = await loadCollectionSegment(collection, segment, signal);
    const prepared = await prepareReferenceOverride(
      segment,
      text,
      `${database.name} · ${locus} ${segment}`,
      "database",
      locus,
      database.id,
      signal,
    );
    referenceCacheRef.current.set(cacheKey, prepared);
    return prepared;
  }

  async function acceptReferenceFile(locus: LocusKey, segment: SegmentKey, file: File) {
    const contextRequest = databaseRequestRef.current;
    const key = referenceCellKey(locus, segment);
    const cellRequest = (cellRequestRef.current[key] ?? 0) + 1;
    cellRequestRef.current[key] = cellRequest;
    referenceAbortRef.current.get(key)?.abort();
    const controller=new AbortController();referenceAbortRef.current.set(key,controller);
    markCellsBusy([key], true);
    try {
      setRunError("");
      const text = await readUploadedText(file,controller.signal);
      if (text.trimStart()[0] !== ">") throw new Error(`${locus} ${segment} references must be FASTA.`);
      const override = await prepareReferenceOverride(segment, text, file.name, "upload", locus, undefined, controller.signal);
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      setCellReferences((current) => ({ ...current, [key]: override }));
      setAlleleExclusions((current) => {
        if (!current[key]) return current;
        const next = { ...current };
        delete next[key];
        return next;
      });
    } catch (error) {
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      if(!isAbortError(error))setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      if (databaseRequestRef.current === contextRequest && cellRequestRef.current[key] === cellRequest) {
        if(referenceAbortRef.current.get(key)===controller)referenceAbortRef.current.delete(key);
        markCellsBusy([key], false);
        setPendingCellSources((current)=>{const next={...current};delete next[key];return next;});
      }
    }
  }

  async function selectCellSource(locus: LocusKey, segment: SegmentKey, sourceId: string) {
    const key = referenceCellKey(locus, segment);
    const cellRequest = (cellRequestRef.current[key] ?? 0) + 1;
    cellRequestRef.current[key] = cellRequest;
    if (sourceId === "upload") return;
    setPendingCellSources((current)=>({...current,[key]:sourceId}));
    setRunError("");
    if (sourceId === DEFAULT_DATABASE_ID) {
      setCellReferences((current) => {
        const next = { ...current };
        delete next[key];
        return next;
      });
      setAlleleExclusions((current) => {
        if (!current[key]) return current;
        const next = { ...current };
        delete next[key];
        return next;
      });
      markCellsBusy([key], false);
      setPendingCellSources((current)=>{const next={...current};delete next[key];return next;});
      return;
    }
    const database = databaseOptions.find((option) => option.id === sourceId)?.database;
    if (!database) {
      setPendingCellSources((current)=>{const next={...current};delete next[key];return next;});
      setRunError("That database is not available for the selected species.");
      return;
    }
    const contextRequest = databaseRequestRef.current;
    referenceAbortRef.current.get(key)?.abort();
    const controller=new AbortController();referenceAbortRef.current.set(key,controller);
    markCellsBusy([key], true);
    try {
      const reference = await databaseCellReference(database, locus, segment, controller.signal);
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      setCellReferences((current) => ({ ...current, [key]: reference }));
      setAlleleExclusions((current) => {
        if (!current[key]) return current;
        const next = { ...current };
        delete next[key];
        return next;
      });
    } catch (error) {
      if (databaseRequestRef.current !== contextRequest || cellRequestRef.current[key] !== cellRequest) return;
      if(!isAbortError(error))setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      if (databaseRequestRef.current === contextRequest && cellRequestRef.current[key] === cellRequest) {
        if(referenceAbortRef.current.get(key)===controller)referenceAbortRef.current.delete(key);
        markCellsBusy([key], false);
        setPendingCellSources((current)=>{const next={...current};delete next[key];return next;});
      }
    }
  }

  async function applyReferenceDatabase(nextId: string) {
    if (nextId === "mixed") return;
    const requestId = ++databaseRequestRef.current;
    databaseAbortRef.current?.abort();
    const controller=new AbortController();databaseAbortRef.current=controller;
    cellRequestRef.current = {};
    setRunError("");
    setPendingDatabaseId(nextId);
    if (nextId === DEFAULT_DATABASE_ID) {
      if(databaseAbortRef.current===controller)databaseAbortRef.current=null;
      setCellReferences({});
      setAlleleExclusions({});
      setEditingReferenceCell(null);
      setBusyCells(new Set());
      setPendingCellSources({});
      setDatabaseBusy(false);
      setPendingDatabaseId(null);
      return;
    }
    const database = databaseOptions.find((option) => option.id === nextId)?.database;
    if (!database) {
      if(databaseAbortRef.current===controller)databaseAbortRef.current=null;
      setPendingDatabaseId(null);
      setRunError("That database is not available for the selected species.");
      return;
    }
    let targetScope = activeScope;
    let collections = collectionsForDatabase(database, targetScope);
    if (!collections.length) {
      targetScope = database.defaultScope;
      collections = collectionsForDatabase(database, targetScope);
    }
    if (!collections.length || !scopes.includes(targetScope)) {
      if(databaseAbortRef.current===controller)databaseAbortRef.current=null;
      setPendingDatabaseId(null);
      setRunError(`${database.name} has no records compatible with this species and search space.`);
      return;
    }
    const currentReceptor = activeScope === "BCR" || activeScope.startsWith("IG") ? "BCR" : "TCR";
    const targetReceptor = targetScope === "BCR" || targetScope.startsWith("IG") ? "BCR" : "TCR";
    const replaceComposition = currentReceptor !== targetReceptor;
    const targets = collections.flatMap((collection) => SEGMENTS.flatMap((segment) => collection.segments[segment] ? [{ locus: collection.locus, segment }] : []));
    const keys = targets.map(({ locus, segment }) => referenceCellKey(locus, segment));
    setDatabaseBusy(true);
    setPendingCellSources(Object.fromEntries(keys.map((key)=>[key,nextId])));
    markCellsBusy(keys, true);
    try {
      const preparedEntries = await Promise.all(targets.map(async ({ locus, segment }) => ({
        key: referenceCellKey(locus, segment),
        reference: await databaseCellReference(database, locus, segment, controller.signal),
      })));
      if (databaseRequestRef.current !== requestId) return;
      if (targetScope !== activeScope) setScope(targetScope);
      setCellReferences((current) => ({ ...(replaceComposition ? {} : current), ...Object.fromEntries(preparedEntries.map((entry) => [entry.key, entry.reference])) }));
      setAlleleExclusions((current) => {
        if (replaceComposition) return {};
        const next = { ...current };
        for (const key of keys) delete next[key];
        return next;
      });
      setEditingReferenceCell(null);
    } catch (error) {
      if (databaseRequestRef.current !== requestId) return;
      if(!isAbortError(error))setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      if (databaseRequestRef.current === requestId) {
        if(databaseAbortRef.current===controller)databaseAbortRef.current=null;
        setDatabaseBusy(false);
        markCellsBusy(keys, false);
        setPendingCellSources({});
        setPendingDatabaseId(null);
      }
    }
  }

  function cancelReferencePreparation(){
    databaseAbortRef.current?.abort();
    referenceAbortRef.current.forEach((controller)=>controller.abort());
  }

  async function exportConfiguredCliConfig(){
    if(!pipeline.enabled||!compiled||!species||!activeDatasets.length)return;
    if(configExporting){configExportAbortRef.current?.abort();return;}
    const invalidDataset=activeDatasets.find((input)=>!input.sampleId.trim()||!input.subjectId.trim());
    if(invalidDataset){setRunError(`Dataset ${invalidDataset.inputName} needs both a sample ID and donor/subject ID.`);return;}
    const controller=new AbortController();configExportAbortRef.current=controller;setConfigExporting(true);setRunError("");
    try{
      const exportPipeline=copyPipeline(pipeline);
      if(exportPipeline.chimera.enabled&&exportPipeline.chimera.msaSource==="selected"){
        const msa=await runKalignTask(compiled[exportPipeline.chimera.segment],controller.signal);
        prepareReferenceMsa(msa);exportPipeline.chimera.msaSource="upload";exportPipeline.chimera.uploadedMsa=msa;exportPipeline.chimera.uploadedMsaName=`${exportPipeline.chimera.segment.toLowerCase()}-reference-msa.fasta`;
      }
      if(controller.signal.aborted)throw new DOMException("Config export cancelled.","AbortError");
      const config=cliConfigFromBrowser({
        studyName:activeInputName,studyDesign,datasets:activeDatasets,references:compiled,species:species.name,scope:activeScope,workers:workerCount,callingProfile,minimumIdentity,minimumConstantPrefixIdentity:constantPrefixIdentityEnabled?minimumConstantPrefixIdentity:null,strand,fastqFilter:copyFastqQualityFilter(fastqFilter),subsample:{enabled:subsampleEnabled,size:Math.max(1,Math.floor(subsampleSize)),seed:Math.trunc(subsampleSeed)},doubleD:{mode:doubleDMode,minimumVjSpan:Math.round(doubleDMinimumSpan),seedLength:Math.round(doubleDSeedLength),pseudoTrim:Math.round(doubleDPseudoTrim),maximumPseudoMismatches:Math.round(doubleDMaximumPseudoMismatches),minimumScoreGain:Math.round(doubleDMinimumScoreGain)},pipeline:exportPipeline,
      });
      downloadBlob(new Blob([`${JSON.stringify(config,null,2)}\n`],{type:"application/json"}),`${activeInputName.replace(/\.[^.]+$/,"")||"swig"}.swig-cli.json`);
    }catch(error){if(!isAbortError(error))setRunError(error instanceof Error?error.message:String(error));}
    finally{if(configExportAbortRef.current===controller)configExportAbortRef.current=null;setConfigExporting(false);}
  }

  function resetReferenceContext(restorePreferred = false) {
    databaseAbortRef.current?.abort();databaseAbortRef.current=null;
    referenceAbortRef.current.forEach((controller)=>controller.abort());referenceAbortRef.current.clear();
    databaseRequestRef.current += 1;
    cellRequestRef.current = {};
    setCellReferences({});
    setAlleleExclusions({});
    setEditingReferenceCell(null);
    setBusyCells(new Set());
    setPendingCellSources({});
    setDatabaseBusy(false);
    setPendingDatabaseId(null);
    referenceCacheRef.current.clear();
    if (restorePreferred && species) {
      const preferred = preferredDatabaseIdFor(species.name, activeScope);
      if (preferred !== DEFAULT_DATABASE_ID) {
        setPendingDatabaseId(preferred);
        void applyReferenceDatabase(preferred);
      }
    }
  }

  function chooseDemo() {
    if (!species) {
      navigate("analyze");
      return;
    }
    const text = makeDemoFasta(species, activeScope);
    if (!text.trim()) {
      setInputError("No complete demo locus is available for this reference selection.");
      navigate("analyze");
      return;
    }
    setPasteText(text);
    setInputSource("paste");
    setInputError("");
    navigate("analyze");
  }

  function requestRun() {
    if (!activeInput || !activeDatasets.length || !compiled || !species) return;
    const pipelineActive = webMode === "repertoire" && pipeline.enabled;
    const subsampleActive = webMode === "repertoire" && subsampleEnabled;
    const fastqFilterActive = webMode === "repertoire" && fastqFilter.enabled;
    const invalidDataset = activeDatasets.find((input) => !input.sampleId.trim() || !input.subjectId.trim());
    if (invalidDataset) {
      setRunError(`Dataset ${invalidDataset.inputName} needs both a sample ID and donor/subject ID.`);
      return;
    }
    if (pipelineActive && pipeline.lineage.enabled && (pipeline.lineage.identity < 0 || pipeline.lineage.identity > 1)) {
      setRunError("Pipeline lineage identity must be between 0 and 1.");
      return;
    }
    if (pipelineActive && pipeline.alleleRefinement.enabled && !pipeline.alleleRefinement.segments.length) {
      setRunError("Select at least one V, D, or J segment for pipeline allele pooling.");
      return;
    }
    if (pipelineActive && pipeline.alleleRefinement.enabled && pipeline.alleleRefinement.reassignmentPolicy === "confidence" && (pipeline.alleleRefinement.applyMinimumPosterior < 0 || pipeline.alleleRefinement.applyMinimumPosterior > 1)) {
      setRunError("Pipeline allele-pooling posterior threshold must be between 0 and 1.");
      return;
    }
    if (pipelineActive && pipeline.chimera.enabled && (pipeline.chimera.posteriorThreshold < 0 || pipeline.chimera.posteriorThreshold > 1)) {
      setRunError("Pipeline chimera posterior threshold must be between 0 and 1.");
      return;
    }
    if (pipelineActive && pipeline.chimera.enabled && pipeline.chimera.msaSource === "upload" && !pipeline.chimera.uploadedMsa) {
      setRunError("Choose an aligned FASTA reference MSA for the pipeline CHMMAIRRa stage, or build it from the selected references.");
      return;
    }
    if (pipelineActive && pipeline.chimera.enabled && pipeline.chimera.msaSource === "upload") {
      try { prepareReferenceMsa(pipeline.chimera.uploadedMsa); }
      catch (error) { setRunError(error instanceof Error ? error.message : String(error)); return; }
    }
    if (databaseBusy || busyCells.size) {
      setRunError("Wait for reference validation to finish.");
      return;
    }
    if (!compiled.counts.V || !compiled.counts.J) {
      setRunError("This reference selection requires both V and J germline records.");
      return;
    }
    if (!compiled.annotation.V.annotated) {
      setRunError("None of the selected V records has a validated FWR/CDR delineation. Provide IMGT-gapped records, SWIGMETA/AIRR-C metadata, or a sufficiently close annotated V reference.");
      return;
    }
    if (!compiled.annotation.J.annotated) {
      setRunError("None of the selected J records has a validated coding frame and conserved F/W–G anchor. Junction and CDR3 calls cannot be made from this set.");
      return;
    }
    if (subsampleActive && (!Number.isFinite(subsampleSize) || subsampleSize < 1)) {
      setRunError("Random subsample size must be at least one record.");
      return;
    }
    if (fastqFilterActive) {
      if (!Number.isFinite(fastqFilter.maximumExpectedErrors) || fastqFilter.maximumExpectedErrors < 0) {
        setRunError("Maximum FASTQ expected errors must be a non-negative number.");
        return;
      }
      const trim = fastqFilter.trim3Prime;
      if (trim.enabled && (
        !Number.isFinite(trim.windowSize) || trim.windowSize < 1
        || !Number.isFinite(trim.minimumMeanPhred) || trim.minimumMeanPhred < 0
        || !Number.isFinite(trim.minimumLength) || trim.minimumLength < 1
      )) {
        setRunError("FASTQ 3′ trimming requires a positive window and retained length, and a non-negative mean Phred threshold.");
        return;
      }
    }
    if (doubleDMode !== "off") {
      if (!compiled.counts.D) {
        setRunError("Double-D screening requires at least one D germline in the composed reference set.");
        return;
      }
      const integerParameters = [doubleDMinimumSpan, doubleDSeedLength, doubleDPseudoTrim, doubleDMaximumPseudoMismatches, doubleDMinimumScoreGain];
      if (integerParameters.some((value) => !Number.isFinite(value) || value < 0) || doubleDSeedLength < 4) {
        setRunError("Double-D screen parameters must be non-negative integers; the exact seed must be at least 4 nt.");
        return;
      }
    }
    if (session && !window.confirm("Starting this analysis will replace the current in-browser result. Cancel to save the current session first.")) return;
    const selectedCount = subsampleActive ? Math.floor(subsampleSize) : null;
    if(projectWorkspace){void run("browser");return;}
    const wantsDisk = outputStorage === "disk" || (outputStorage === "auto" && activeDatasets.some((input) => likelyLargeInput(input, selectedCount)));
    if (wantsDisk && !savePicker() && outputStorage === "disk") {
      setRunError("Direct-to-disk streaming is unavailable in this browser. Use Auto/browser storage or a Chromium-based browser.");
      return;
    }
    if (wantsDisk && savePicker()) {
      setOutputPrompt(true);
      return;
    }
    void run("browser");
  }

  async function run(outputDestination: "browser" | "disk") {
    if (!activeInput || !activeDatasets.length || !compiled || !species) return;
    setOutputPrompt(false);
    const datasetSnapshot = activeDatasets.map((input) => ({ ...input }));
    const runName = activeInputName;
    const runSubsampleEnabled = webMode === "repertoire" && subsampleEnabled;
    const fastqFilterSnapshot = webMode === "repertoire" ? copyFastqQualityFilter(fastqFilter) : copyFastqQualityFilter();

    let directOutput: DirectAirrOutput | undefined;
    let preparedProject:PreparedProjectRun|undefined;
    if(projectWorkspace){
      try{
        setProjectStatus("Creating numbered project run…");
        preparedProject=await prepareProjectRun(projectWorkspace,runName,APP_VERSION);
        await writeProjectDatasetManifest(preparedProject.workspace,preparedProject.run,datasetSnapshot);
        directOutput={handle:preparedProject.airrHandle,writable:{
          write:(data)=>preparedProject!.writable.write(data instanceof Uint8Array?new Uint8Array(data).buffer:data),
          close:()=>preparedProject!.writable.close(),
          abort:()=>preparedProject!.writable.abort(),
        }};
        setProjectWorkspace(preparedProject.workspace);
        setProjectStatus(`Writing ${preparedProject.run.airrPath}`);
      }catch(error){await preparedProject?.writable.abort().catch(()=>undefined);setRunError(error instanceof Error?error.message:String(error));return;}
    }else if (outputDestination === "disk") {
      const picker = savePicker();
      if (!picker) {
        setRunError("Direct-to-disk streaming is unavailable in this browser. Use Auto/browser storage or a Chromium-based browser.");
        return;
      }
      try {
        const handle = await picker.call(window, {
          suggestedName: outputName(runName),
          types: [{ description: "Swig AIRR analysis output", accept: { "text/tab-separated-values": [".tsv"] } }],
        });
        directOutput = { handle, writable: await handle.createWritable() };
      } catch (error) {
        if (!(error instanceof DOMException && error.name === "AbortError")) {
          setRunError(error instanceof Error ? error.message : String(error));
        }
        return;
      }
    }

    if (session) await session.store.clear();
    setSession(null);
    const store = new AirrResultStore(directOutput);
    const controller = new AbortController();
    abortRef.current = controller;
    setRunning(true);
    setRunError("");
    setProgress({ stage: "Preparing analysis", value: 0.01 });
    setAssignmentTelemetry(initializingAssignmentTelemetry(workerCount));
    const started = performance.now();
    void navigator.storage?.persist?.().catch(() => false);
    const doubleD: DoubleDScreenOptions = {
      mode: doubleDMode,
      minimumVjSpan: Math.round(doubleDMinimumSpan),
      seedLength: Math.round(doubleDSeedLength),
      pseudoTrim: Math.round(doubleDPseudoTrim),
      maximumPseudoMismatches: Math.round(doubleDMaximumPseudoMismatches),
      minimumScoreGain: Math.round(doubleDMinimumScoreGain),
    };
    try {
      const result = await withAnalysisWebLock(controller.signal, setAnalysisLockState, async () => {
        let count = 0;
        let inputRecords = 0;
        let workers = 1;
        let fastqFilterStats = emptyFastqQualityFilterStats(fastqFilterSnapshot.enabled, false);
        // Never mix exact read counts with a one-read placeholder for unknown inputs.
        const allCounted = datasetSnapshot.every((input) => input.count !== null);
        const estimates = allCounted ? datasetSnapshot.map((input) => input.count!)
          : await Promise.all(datasetSnapshot.map((input) => estimateSequenceCharacters(input.source)));
        const weights = estimates.map((estimate, index) => runSubsampleEnabled
          ? Math.min(datasetSnapshot[index].count ?? (allCounted ? estimate : estimate / 500), Math.floor(subsampleSize))
          : estimate);
        const totalWeight = weights.reduce((sum, value) => sum + value, 0) || datasetSnapshot.length;
        let completedWeight = 0;
        for (let datasetIndex = 0; datasetIndex < datasetSnapshot.length; datasetIndex += 1) {
          const input = datasetSnapshot[datasetIndex];
          const manifest: DatasetManifestEntry = {
            datasetId: input.datasetId,
            inputName: input.inputName,
            sampleId: input.sampleId,
            subjectId: input.subjectId,
            cohort: input.cohort,
            timepoint: input.timepoint,
            compartment: input.compartment,
            records: input.count,
          };
          const weight = weights[datasetIndex];
          setAssignmentTelemetry(initializingAssignmentTelemetry(workerCount));
          const completed = await runAligner({
            query: input.source,
            format: input.formatCode,
            references: compiled,
            settings: {
              callingProfile,
              minimumIdentity,
              minimumConstantPrefixIdentity: constantPrefixIdentityEnabled ? minimumConstantPrefixIdentity : null,
              strand,
              doubleD,
            },
            workers: workerCount,
            countHint: input.count,
            subsample: runSubsampleEnabled ? { size: Math.floor(subsampleSize), seed: stableDatasetSeed(subsampleSeed, datasetIndex) } : undefined,
            fastqFilter: fastqFilterSnapshot,
            signal: controller.signal,
            onProgress: (stage, value) => setProgress({
              stage: datasetSnapshot.length > 1 ? `${input.sampleId} · ${stage}` : stage,
              value: overallAssignmentProgress(completedWeight, weight, totalWeight, value),
            }),
            onTelemetry: setAssignmentTelemetry,
            onBatch: async (batch) => {
              const annotated = annotateAirrBatch(batch.header, batch.body, manifest);
              const doubleDBatch = batch.doubleDHeader !== undefined && batch.doubleDBody !== undefined
                ? annotateDoubleDBatch(batch.doubleDHeader, batch.doubleDBody, manifest)
                : undefined;
              await store.appendBatch(annotated.header, annotated.body, doubleDBatch);
              setProgress((current) => ({
                ...current,
                stage: `${input.sampleId} · committed ${batch.processed.toLocaleString()}${batch.total === null ? "" : ` / ${batch.total.toLocaleString()}`} AIRR records${store.doubleDCount ? ` · ${store.doubleDCount.toLocaleString()} double-D supported` : ""}`,
              }));
            },
          });
          count += completed.count;
          inputRecords += completed.inputRecords;
          workers = Math.max(workers, completed.workers);
          fastqFilterStats = addFastqQualityFilterStats(fastqFilterStats, completed.fastqFilter);
          completedWeight += weight;
        }
        setAssignmentTelemetry((current) => current ? {
          ...current,
          phase: "finalizing",
          activeWorkers: 0,
          recordsOutstanding: 0,
        } : current);
        await store.finalize();
        return { count, inputRecords, workers, fastqFilterStats };
      });
      setProgress({ stage: "Results ready", value: 1 });
      const resultSession:ResultSession={
        id: Date.now(),
        webMode,
        store,
        total: result.count,
        seconds: (performance.now() - started) / 1000,
        inputName: runName,
        datasets: datasetSnapshot.map((input) => ({ datasetId: input.datasetId, inputName: input.inputName, inputPath:input.inputPath, inputFormat:input.inputFormat, gzipRange:input.gzipRange?{...input.gzipRange}:undefined, sampleId: input.sampleId, subjectId: input.subjectId, cohort: input.cohort, timepoint: input.timepoint, compartment: input.compartment, records: input.count })),
        studyDesign,
        pipeline: webMode === "repertoire" ? copyPipeline(pipeline) : { ...copyPipeline(pipeline), enabled: false },
        species: species.name,
        scope: activeScope,
        facets: store.facets(),
        summary: store.summary,
        workers: result.workers,
        outputBytes: store.outputBytes,
        streamedDirectly: store.streamedDirectly,
        inputTotal: result.inputRecords,
        subsampleSize: runSubsampleEnabled ? Math.floor(subsampleSize) : null,
        subsampleSeed: runSubsampleEnabled ? Math.trunc(subsampleSeed) : null,
        fastqFilter: fastqFilterSnapshot,
        fastqFilterStats: result.fastqFilterStats,
        references: compiled,
        referenceExclusions: Object.fromEntries(Object.entries(activeAlleleExclusions).map(([key, names]) => [key, [...names]])),
        callingProfile,
        minimumIdentity,
        minimumConstantPrefixIdentity: constantPrefixIdentityEnabled ? minimumConstantPrefixIdentity : null,
        strand,
        doubleD,
        doubleDCount: store.doubleDCount,
        sampleColors: createSampleColorMap(datasetSnapshot),
        project:preparedProject?.workspace,
      };
      if(preparedProject){
        try{
          await appendProjectLog(preparedProject.workspace,preparedProject.run,"annotation_complete",{records:result.count,inputRecords:result.inputRecords,workers:result.workers,seconds:resultSession.seconds,outputBytes:store.outputBytes,fastqFilter:fastqFilterSnapshot,fastqFilterStats:result.fastqFilterStats});
          const artifact=await sessionArtifact(resultSession,resultSession.datasets,resultSession.sampleColors,{workingStages:[]});
          await saveProjectCheckpoint(preparedProject.workspace,artifact,await encodeSession(artifact),"annotation_complete");
          resultSession.projectStatus=`Project checkpoint ${preparedProject.run.checkpointCount} written`;
          setProjectStatus(`${preparedProject.run.id} · annotation checkpoint written`);
        }catch(error){
          resultSession.projectStatus="AIRR result written; initial state checkpoint failed";
          setProjectStatus(error instanceof Error?error.message:String(error));
        }
      }
      setSession(resultSession);
      setPage("results");
      window.scrollTo({ top: 0 });
    } catch (error) {
      await store.abort();
      await store.clear();
      if(preparedProject)await appendProjectLog(preparedProject.workspace,preparedProject.run,error instanceof DOMException&&error.name==="AbortError"?"analysis_cancelled":"analysis_failed",{message:error instanceof Error?error.message:String(error)}).catch(()=>undefined);
      if (!(error instanceof DOMException && error.name === "AbortError")) {
        setRunError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      abortRef.current = null;
      setRunning(false);
      setAnalysisLockState("unsupported");
    }
  }

  const hasBcr = scopes.includes("BCR");
  const hasTcr = scopes.includes("TCR");
  const focusedWebMode = webMode !== "repertoire";
  const analysisTitle = webMode === "vdj" ? "VDJ assignment and visualization" : webMode === "lineage" ? "Single lineage analysis" : "End-to-end repertoire analysis";

  return (
    <div className="site-shell">
      <AppHeader page={page} hasResults={Boolean(session)} hasLineageStudy={Boolean(pendingLineageStudy&&lineageStudyAirr)} projectName={projectWorkspace?.root.name} onNavigate={navigate} onLoadSession={()=>sessionInputRef.current?.click()} onLoadLineageStudy={()=>lineageStudyInputRef.current?.click()} onOpenProject={()=>void chooseProjectDirectory()} />
      <input ref={sessionInputRef} className="visually-hidden" type="file" accept=".swig-session,.json,.gz" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(file)void acceptSessionFile(file);}}/>
      <input ref={linkedAirrInputRef} className="visually-hidden" type="file" accept=".tsv,.tsv.gz,.gz" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(file)void restoreLinkedAirr(file);}}/>
      <input ref={lineageStudyInputRef} className="visually-hidden" type="file" accept=".json,.gz,.swig-lineage-study" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(file)void acceptLineageStudyFile(file);}}/>
      <input ref={lineageStudyAirrInputRef} className="visually-hidden" type="file" accept=".tsv,.airr.tsv" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(file)restoreLineageStudyAirr(file);}}/>
      {page === "home" && <LandingPage onChoose={chooseWebAnalysisMode} />}

      {page === "analyze" && (
        <main className="analysis-page">
          <section className="analysis-intro"><WorkflowStepper active={running ? 2 : 1} /><div><p className="eyebrow"><span>{analysisTitle}</span></p><h1>{running ? "Calling rearrangements…" : webMode === "vdj" ? "Annotate and inspect one input." : webMode === "lineage" ? "Annotate one pre-defined lineage." : "Configure an annotation run."}</h1><p>{running ? "SwiftIG is processing bounded batches and writing AIRR records into a browser-local index." : webMode === "lineage" ? "Every input sequence will enter the alignment workbench as one lineage; clustering and lineage selection are bypassed." : focusedWebMode ? "Upload or paste sequences, choose germlines, and run V(D)J/CDR3 annotation." : "Provide sequences, select the biological search space, and specify any germline replacements."}</p></div></section>

          {running ? <AnalysisProgress stage={progress.stage} value={progress.value} telemetry={assignmentTelemetry} onCancel={() => abortRef.current?.abort()} lockState={analysisLockState} hidden={pageHidden} /> : (
            <div className="analysis-layout single-action-layout">
              <div className="analysis-forms">
                <section className="analysis-card input-card">
                  <header><span className="card-number">01</span><div><h2>{focusedWebMode ? "Input sequences" : "Datasets and study structure"}</h2><p>{focusedWebMode ? "Load one file or paste sequences directly. No sample naming or study metadata are required." : "Load one or more datasets, a directory tree, or paste one dataset directly."}</p></div>{activeInput && <span className="ready-tag">{focusedWebMode ? "input ready" : `${activeDatasets.length} dataset${activeDatasets.length===1?"":"s"}`}</span>}</header>
                  {!focusedWebMode && <div className="project-directory-panel"><div><span>Optional project directory</span><strong>{projectWorkspace?projectWorkspace.root.name:"No project directory selected"}</strong><small>{projectWorkspace?projectStatus||"AIRR output, state checkpoints, and an event log will be written here.":projectDirectoriesSupported()?"Select a directory to write numbered runs, intermediate state, and logs automatically. Selecting an existing Swig project restores its active run.":"Unavailable in this browser; portable session files remain available."}</small></div><button type="button" disabled={projectBusy||!projectDirectoriesSupported()} onClick={()=>void chooseProjectDirectory()}>{projectBusy?"Reading directory…":projectWorkspace?"Change / load project":"Select / load project"}</button>{projectWorkspace&&<button type="button" onClick={()=>{setProjectWorkspace(null);setProjectStatus("");}}>Use browser storage only</button>}</div>}
                  <div className="source-tabs" role="tablist"><button className={inputSource === "upload" ? "active" : ""} type="button" onClick={() => setInputSource("upload")}>{focusedWebMode ? "Upload file" : "Load dataset(s)"}</button><button className={inputSource === "paste" ? "active" : ""} type="button" onClick={() => setInputSource("paste")}>{focusedWebMode ? "Paste sequences" : "Paste one dataset"}</button></div>
                  <input ref={inputRef} className="visually-hidden" type="file" multiple={!focusedWebMode} accept=".fa,.fasta,.fna,.fas,.fq,.fastq,.tsv,.csv,.txt,.gz" onChange={(event: ChangeEvent<HTMLInputElement>) => {
                    const files = [...(event.target.files ?? [])];
                    if (files.length) void acceptInputFiles(files);
                    event.target.value = "";
                  }} />
                  <input ref={directoryInputRef} className="visually-hidden" type="file" multiple {...{webkitdirectory:"",directory:""}} onChange={(event:ChangeEvent<HTMLInputElement>)=>{const files=[...(event.target.files??[])];event.target.value="";if(files.length)void acceptInputCandidates(candidatesFromDirectoryPicker(files));}}/>
                  {inputSource === "upload" ? (
                    fileInputs.length ? (
                      focusedWebMode ? (
                    <div className="focused-input-summary"><div><span>Selected input</span><strong>{fileInputs[0].name}</strong><small>{fileInputs[0].count === null ? "Records counted while streaming" : `${fileInputs[0].count.toLocaleString()} records`} · {fileInputs[0].format} · {bytes(fileInputs[0].size)}{fileInputs[0].gzipMembersMerged&&fileInputs[0].gzipMemberCount&&fileInputs[0].gzipMemberCount>1?` · ${fileInputs[0].gzipMemberCount} gzip members merged`:""}</small></div><button type="button" onClick={() => inputRef.current?.click()}>Replace file</button><button type="button" onClick={() => setFileInputs([])}>Remove</button></div>
                  ) : (
                    <div className="dataset-import-stack" onDragOver={(event:DragEvent)=>event.preventDefault()} onDrop={(event:DragEvent)=>{event.preventDefault();void collectDroppedInput(event.dataTransfer).then((candidates)=>acceptInputCandidates(candidates)).catch((error)=>setInputError(error instanceof Error?error.message:String(error)));}}>
                      <div className="dataset-import-heading"><div><strong>{fileInputs.length.toLocaleString()} dataset{fileInputs.length === 1 ? "" : "s"}</strong><span>{inputInspecting?"Inspecting gzip member boundaries…":knownInputCount === null ? "Record totals will be counted during analysis" : `${knownInputCount.toLocaleString()} total records`}</span></div><div>{inputInspecting?<button className="post-cancel" type="button" onClick={cancelInputInspection}>Cancel inspection</button>:<><button type="button" onClick={() => inputRef.current?.click()}>＋ Add files</button><button type="button" onClick={()=>directoryInputRef.current?.click()}>＋ Add directory</button></>}</div></div>
                      <div className="dataset-manifest-table" role="table" aria-label="Dataset metadata">
                        <div className="dataset-manifest-head" role="row"><span>Input</span><span>Sample ID</span><span>Donor / subject</span><span>Cohort</span><span>Timepoint</span><span>Compartment / tissue</span><span /></div>
                        {fileInputs.map((input) => <div className="dataset-manifest-row" role="row" key={input.datasetId}>
                          <div><strong title={input.name}>{input.name}</strong><small>{input.datasetId} · {input.count === null ? "stream counted" : `${input.count.toLocaleString()} records`} · {input.format} · {bytes(input.size)}{input.gzipMembersMerged&&input.gzipMemberCount&&input.gzipMemberCount>1?` · ${input.gzipMemberCount} gzip members merged`:""}</small></div>
                          <label><span className="visually-hidden">Sample ID for {input.name}</span><input value={input.sampleId} onChange={(event) => updateDataset(input.datasetId, { sampleId: event.target.value })} /></label>
                          <label><span className="visually-hidden">Donor or subject for {input.name}</span><input value={input.subjectId} onChange={(event) => updateDataset(input.datasetId, { subjectId: event.target.value })} /></label>
                          <label><span className="visually-hidden">Cohort for {input.name}</span><input value={input.cohort} onChange={(event) => updateDataset(input.datasetId, { cohort: event.target.value })} /></label>
                          <label><span className="visually-hidden">Timepoint for {input.name}</span><input value={input.timepoint} placeholder="e.g. day_0" onChange={(event) => updateDataset(input.datasetId, { timepoint: event.target.value })} /></label>
                          <label><span className="visually-hidden">Compartment or tissue for {input.name}</span><input value={input.compartment} placeholder="e.g. blood" onChange={(event) => updateDataset(input.datasetId, { compartment: event.target.value })} /></label>
                          <button type="button" aria-label={`Remove ${input.name}`} onClick={() => setFileInputs((current) => current.filter((candidate) => candidate.datasetId !== input.datasetId))}>×</button>
                        </div>)}
                      </div>
                      <div className="dataset-add-dropzone">Drop additional files or a directory anywhere in this dataset panel.</div>
                    </div>
                      )
                    ) : (
                    <div className="input-dropzone" role="button" aria-disabled={inputInspecting} tabIndex={0} onKeyDown={(event)=>{if(!inputInspecting&&(event.key==="Enter"||event.key===" ")){event.preventDefault();inputRef.current?.click();}}} onClick={() => {if(!inputInspecting)inputRef.current?.click();}} onDragOver={(event: DragEvent) => event.preventDefault()} onDrop={(event: DragEvent) => {
                      event.preventDefault();
                      if(inputInspecting)return;
                      void collectDroppedInput(event.dataTransfer).then((candidates)=>acceptInputCandidates(candidates)).catch((error)=>setInputError(error instanceof Error?error.message:String(error)));
                    }}><span>＋</span><strong>{inputInspecting?"Inspecting gzip member boundaries…":focusedWebMode?"Drop one input file here":"Drop files or an entire directory here"}</strong><small>.fasta(.gz) · .fastq(.gz) · AIRR .tsv(.gz)</small><i>{inputInspecting?"No file is committed until inspection finishes":focusedWebMode?"Choose file":"Choose files"}</i>{inputInspecting?<button className="post-cancel" type="button" onClick={(event)=>{event.stopPropagation();cancelInputInspection();}}>Cancel inspection</button>:!focusedWebMode&&<button type="button" onClick={(event)=>{event.stopPropagation();directoryInputRef.current?.click();}}>Choose directory</button>}</div>
                    )
                  ) : (
                    <div className="paste-input"><textarea spellCheck={false} value={pasteText} onChange={(event) => setPasteText(event.target.value)} placeholder={">sequence_1\nCAGGTGCAGCTGGTG...\n\n—or—\n\nsequence_id\tsequence\nread_1\tCAGGTGCAGCTGGTG..."} /><footer><span>{pasteInput ? `${pasteInput.count?.toLocaleString()} ${pasteInput.format} records detected` : pasteText.trim() ? "Waiting for valid FASTA, FASTQ, or AIRR…" : "Nothing pasted yet"}</span><button type="button" onClick={chooseDemo}>Insert demo</button></footer>{pasteInput && !focusedWebMode && <div className="paste-metadata"><label><span>Sample ID</span><input value={pasteMetadata.sampleId} onChange={(event)=>setPasteMetadata((current)=>({...current,sampleId:event.target.value}))}/></label><label><span>Donor / subject</span><input value={pasteMetadata.subjectId} onChange={(event)=>setPasteMetadata((current)=>({...current,subjectId:event.target.value}))}/></label><label><span>Cohort</span><input value={pasteMetadata.cohort} onChange={(event)=>setPasteMetadata((current)=>({...current,cohort:event.target.value}))}/></label><label><span>Timepoint</span><input value={pasteMetadata.timepoint} onChange={(event)=>setPasteMetadata((current)=>({...current,timepoint:event.target.value}))}/></label><label><span>Compartment / tissue</span><input value={pasteMetadata.compartment} onChange={(event)=>setPasteMetadata((current)=>({...current,compartment:event.target.value}))}/></label></div>}</div>
                  )}
                  {!focusedWebMode && activeDatasets.length > 0 && <div className="study-design-panel">
                    <header><div><span>Study behavior</span><strong>Define which biological boundaries downstream methods may cross.</strong></div>{activeDatasets.length > 1 && <label><span>Study / output name</span><input value={studyName} onChange={(event)=>setStudyName(event.target.value)} /></label>}</header>
                    <div className="study-design-options">
                      {([
                        ["independent", "Independent samples", "Collapse and lineages stay within each sample."],
                        ["cohort", "Cohort study", "Samples stay independent; cohort metadata supports comparison and plotting."],
                        ["longitudinal", "Longitudinal / compartmental", "Collapse stays within sample; lineages may span timepoints or compartments for the same donor."],
                        ["technical", "Technical replicates", "Files sharing a sample ID may collapse together; lineages stay within that sample."],
                        ["custom", "Custom boundaries", "Choose collapse and lineage scopes independently in pipeline or post-analysis."],
                      ] as Array<[StudyDesign,string,string]>).map(([value,label,description])=><label className={studyDesign===value?"selected":""} key={value}><input type="radio" checked={studyDesign===value} onChange={()=>applyStudyDesign(value)}/><span><strong>{label}</strong><small>{description}</small></span></label>)}
                    </div>
                    {studyDesign === "longitudinal" && <p className="scientific-note"><span>i</span>Give all timepoints and compartments from one donor the same donor/subject ID. Keep distinct biological specimens as separate sample IDs. Lineage assignment is donor-bounded by default.</p>}
                  </div>}
                  {inputError && <p className="inline-error" role="alert">{inputError}</p>}
                </section>

                <section className="analysis-card reference-card">
                  <header><span className="card-number">02</span><div><h2>Biological search space</h2><p>Choose species, receptor, locus, and reference database. Customize individual segments only when needed.</p></div>{(Object.keys(cellReferences).length > 0 || Object.keys(alleleExclusions).length > 0) && <button className="reset-button" type="button" onClick={() => resetReferenceContext(true)}>Restore recommended sources</button>}</header>
                  <div className="reference-selectors">
                    <label><span>Species / strain</span><select value={species?.name??""} disabled={!pack} onChange={(event)=>selectSpecies(event.target.value)}>{speciesList.map((item)=><option value={item.name} key={item.name}>{friendlySpecies(item.name)}</option>)}</select></label>
                    <label className={databaseBusy?"database-selector preparing":"database-selector"}><span>Database</span><select value={displayedDatabaseId} disabled={!pack} onChange={(event) => void applyReferenceDatabase(event.target.value)}>{databaseOptions.map((option) => <option value={option.id} key={option.id}>{option.label}</option>)}{displayedDatabaseId === "mixed" && <option value="mixed" disabled>Mixed sources · configured below</option>}</select><small>{databaseBusy?"Selection applied; preparing reference metadata…":displayedDatabaseId==="mixed"?"Mixed by locus or segment":"Ready"}</small></label>
                    <div className="receptor-selector"><span>Receptor</span><div><button className={receptor === "BCR" ? "active" : ""} type="button" disabled={!hasBcr} onClick={() => selectReceptor("BCR")}>BCR <small>IG</small></button><button className={receptor === "TCR" ? "active" : ""} type="button" disabled={!hasTcr} onClick={() => selectReceptor("TCR")}>TCR <small>TR</small></button></div></div>
                    <label><span>Chain / locus</span><select value={activeScope} onChange={(event) => setScope(event.target.value as ScopeKey)}>{receptorScopes.map((value) => <option value={value} key={value}>{LOCUS_LABELS[value]}</option>)}</select></label>
                  </div>
                  {packError && <p className="inline-error" role="alert">{packError}</p>}
                  <CompositionSummary databases={summaryDatabases} hasUploads={hasUploadedReferences} excludedAlleles={activeAlleleExclusionCount} busy={databaseBusy || Boolean(busyCells.size)} release={pack?.release ?? ""} onCancel={cancelReferencePreparation} />
                  <div className="reference-downloads" aria-label="Download assignment databases"><span>Download active FASTA</span>{(["all", ...SEGMENTS] as const).map((segment) => {
                    const fasta = compiled ? (segment === "all" ? SEGMENTS.map((key) => compiled[key]).filter((text) => text.trim()).join("\n") : compiled[segment]) : "";
                    return <button type="button" key={segment} disabled={databaseBusy || Boolean(busyCells.size) || !fasta.trim()} title="Download the exact selected assignment references with metadata and allele exclusions" onClick={() => downloadBlob(new Blob([fasta], { type: "text/plain;charset=utf-8" }), `${inputStem(species?.name || "references")}.${activeScope}.${segment}.fasta`)}>{segment === "all" ? "All" : segment} ↓</button>;
                  })}</div>
                  <details className="reference-composition-details">
                    <summary><span><b>Customize individual loci, V/D/J/C sources, or allele inclusion</b><small>Replace any segment with another published collection or a local FASTA; exclude individual alleles before assignment.</small></span></summary>
                    {species && <ReferenceCompositionMatrix species={species} scope={activeScope} references={cellReferences} exclusions={alleleExclusions} busyCells={busyCells} pendingSources={pendingCellSources} onSelect={(locus, segment, sourceId) => void selectCellSource(locus, segment, sourceId)} onFile={(locus, segment, file) => void acceptReferenceFile(locus, segment, file)} onEditAlleles={(locus,segment)=>setEditingReferenceCell({locus,segment})} />}
                    <p className="reference-footnote"><span>i</span><b>{compiled?.loci.join(" + ") || "No locus"}</b> · {databaseLabel}. The database selector applies a preset only where compatible; every matrix cell remains independently selectable, replaceable from a local file, and filterable by exact allele identifier. Exclusions are applied to the composed FASTA before SwiftIG builds its indexes; the source database/file is not modified. Replacing a cell's source clears that cell's exclusions. In-browser preprocessing assigns V FWR/CDR boundaries and J frame/junction-anchor metadata by transfer from validated, locus-matched IMGT relatives; metadata are retained only when mapped intervals and conserved anchors validate. D and C records are validated and indexed but do not have FWR/CDR boundaries.</p>
                  </details>
                </section>

                {editingReferenceCell && <div className="output-modal-backdrop reference-allele-modal-backdrop" role="presentation"><section className="output-modal reference-allele-modal" role="dialog" aria-modal="true" aria-labelledby="reference-allele-modal-title">
                  <button className="output-modal-close" type="button" onClick={()=>setEditingReferenceCell(null)} aria-label="Close allele editor">×</button>
                  <span className="output-direction">ASSIGNMENT DATABASE · {editingReferenceCell.locus} {editingReferenceCell.segment}</span>
                  <h2 id="reference-allele-modal-title">Include or exclude individual alleles</h2>
                  <p><b>{editingReferenceSource}</b>. Checked records are excluded from the exact {editingReferenceCell.locus} {editingReferenceCell.segment} FASTA passed to V(D)J assignment. Other loci and segments are unchanged.</p>
                  <ReferenceAlleleExclusionEditor
                    key={editingReferenceKey}
                    fasta={editingReferenceFasta}
                    excluded={alleleExclusions[editingReferenceKey] ?? []}
                    onChange={(excluded) => setAlleleExclusions((current) => {
                      const next = { ...current };
                      if (excluded.length) next[editingReferenceKey] = excluded;
                      else delete next[editingReferenceKey];
                      return next;
                    })}
                    label={`${editingReferenceCell.locus} ${editingReferenceCell.segment}`}
                    description="Search by exact FASTA identifier or a substring, then exclude individual matches or the complete filtered result. Exclusions affect assignment, alternative-hit reporting, alignments, CDR/junction calls, and downstream analyses that use this run's references. The original source is retained unchanged."
                    variant="panel"
                  />
                  <div className="reference-allele-modal-actions"><small>{(alleleExclusions[editingReferenceKey] ?? []).length.toLocaleString()} excluded · {filterReferenceFasta(editingReferenceFasta, alleleExclusions[editingReferenceKey] ?? []).retained.toLocaleString()} retained</small><button className="output-save-primary" type="button" onClick={()=>setEditingReferenceCell(null)}>Apply to assignment setup</button></div>
                </section></div>}

                <section className="analysis-card settings-card">
                  <header><span className="card-number">03</span><div><h2>Assignment options</h2><p>{focusedWebMode ? "Advanced mapping, calling, and Double-D controls remain available without adding downstream repertoire stages." : "Configure optional read preprocessing and rare-rearrangement screening. Calibrated compute settings are available below."}</p></div></header>
                  {!focusedWebMode && <details className={`fastq-quality-control ${fastqFilter.enabled ? "active" : ""}`}>
                    <summary><span><b>FASTQ quality filter</b><small>Expected errors and optional 3′ quality trimming before V(D)J assignment.</small></span><i>{fastqFilter.enabled ? `Enabled · EE ≤ ${fastqFilter.maximumExpectedErrors}` : "Off"}</i></summary>
                    <div className="fastq-quality-body">
                      <label className="fastq-quality-switch"><input type="checkbox" checked={fastqFilter.enabled} onChange={(event)=>setFastqFilter((current)=>({...current,enabled:event.target.checked}))}/><span><b>Enable FASTQ filtering</b><small>Each retained read must have Σ 10<sup>−Q/10</sup> ≤ the maximum expected errors. Processing is streaming and linear in the number of quality scores.</small></span></label>
                      {fastqFilter.enabled && <>
                        <div className="fastq-quality-fields">
                          <label><span>Maximum expected errors</span><CommitNumberInput min="0" step="0.001" value={fastqFilter.maximumExpectedErrors} onCommit={(maximumExpectedErrors)=>setFastqFilter((current)=>({...current,maximumExpectedErrors}))}/><small>Keep when EE ≤ threshold</small></label>
                          <label><span>Quality encoding</span><select value={fastqFilter.phredOffset} onChange={(event)=>setFastqFilter((current)=>({...current,phredOffset:Number(event.target.value) as 33|64}))}><option value={33}>Phred+33 · standard</option><option value={64}>Phred+64 · legacy</option></select><small>Invalid characters stop the run</small></label>
                        </div>
                        <label className="fastq-trim-switch"><input type="checkbox" checked={fastqFilter.trim3Prime.enabled} onChange={(event)=>setFastqFilter((current)=>({...current,trim3Prime:{...current.trim3Prime,enabled:event.target.checked}}))}/><span><b>Trim low-quality 3′ ends first</b><small>Remove trailing bases until the terminal window reaches the selected mean Phred score. EE is then tested only on the retained bases.</small></span></label>
                        {fastqFilter.trim3Prime.enabled && <div className="fastq-quality-fields trim">
                          <label><span>Window bases</span><CommitNumberInput min="1" max="100" step="1" value={fastqFilter.trim3Prime.windowSize} onCommit={(windowSize)=>setFastqFilter((current)=>({...current,trim3Prime:{...current.trim3Prime,windowSize}}))}/><small>Terminal running window</small></label>
                          <label><span>Minimum mean Phred</span><CommitNumberInput min="0" max="93" step="1" value={fastqFilter.trim3Prime.minimumMeanPhred} onCommit={(minimumMeanPhred)=>setFastqFilter((current)=>({...current,trim3Prime:{...current.trim3Prime,minimumMeanPhred}}))}/><small>Stop trimming when met</small></label>
                          <label><span>Minimum retained length</span><CommitNumberInput min="1" step="1" value={fastqFilter.trim3Prime.minimumLength} onCommit={(minimumLength)=>setFastqFilter((current)=>({...current,trim3Prime:{...current.trim3Prime,minimumLength}}))}/><small>Shorter trimmed reads are discarded</small></label>
                        </div>}
                        {fastqDatasetCount === 0 ? <p className="fastq-format-warning" role="status"><span>!</span>No selected input is FASTQ. FASTA and AIRR records will pass through this step unchanged because they contain no Phred scores.</p> : nonFastqDatasetCount > 0 ? <p className="fastq-format-warning" role="status"><span>i</span>{fastqDatasetCount} FASTQ dataset{fastqDatasetCount===1?"":"s"} will be filtered; {nonFastqDatasetCount} FASTA/AIRR dataset{nonFastqDatasetCount===1?"":"s"} will pass through unchanged.</p> : <p className="fastq-order-note"><span>Order</span>3′ trim → minimum retained length → expected-error filter → optional random subsample → V(D)J assignment.</p>}
                      </>}
                    </div>
                  </details>}
                  {!focusedWebMode && <div className={`subsample-control ${subsampleEnabled ? "active" : ""}`}><label className="subsample-switch"><input type="checkbox" checked={subsampleEnabled} onChange={(event) => setSubsampleEnabled(event.target.checked)} /><span><b>Analyze a random subsample</b><small>Exact reservoir sampling scans the full stream but retains only the requested records {activeDatasets.length > 1 ? "from each dataset" : ""} in memory and output.{fastqFilter.enabled ? " Sampling is from reads retained by FASTQ QC; non-FASTQ records pass through that step." : ""}</small></span></label>{subsampleEnabled && <div className="subsample-fields"><label><span>{activeDatasets.length > 1 ? "Records per dataset" : "Records to analyze"}</span><CommitNumberInput min="1" step="1000" value={subsampleSize} onCommit={setSubsampleSize} /></label><label><span>Base random seed</span><CommitNumberInput step="1" value={subsampleSeed} onCommit={setSubsampleSeed} /></label></div>}</div>}
                  <details className="advanced-settings progressive-settings">
                    <summary title="Calling profile, strand, workers, output storage, and identity floor"><span><b>Advanced options</b><small>Calling profile, strand, workers, output storage, and identity floor.</small></span></summary>
                    <div className="advanced-settings-grid">
                    <label><span>Calling profile</span><select value={callingProfile} onChange={(event) => setCallingProfile(event.target.value as CallingProfile)}><option value="r_optimized">R-optimized · default</option><option value="sensitive_d">Sensitive-D</option><option value="truth_optimized">Truth-optimized</option><option value="igblast_compatible">IgBLAST-agreement</option></select></label>
                    <label><span>Search strand</span><select value={strand} onChange={(event) => setStrand(Number(event.target.value) as 0 | 1 | 2)}><option value={0}>Both orientations</option><option value={1}>Plus only</option><option value={2}>Minus only</option></select></label>
                    <label><span>Parallel compute workers</span><CommitNumberInput min="1" max={browserWorkerLimit()} step="1" value={workerCount} onCommit={(value)=>setWorkerCount(Math.max(1,Math.min(browserWorkerLimit(),Math.round(value))))}/><small>{recommendedWorkerCount()} recommended on this device · {browserWorkerLimit()} maximum</small></label>
                    {!focusedWebMode && <label><span>AIRR results destination</span><select disabled={Boolean(projectWorkspace)} value={projectWorkspace?"project":outputStorage} onChange={(event) => setOutputStorage(event.target.value as OutputStorageMode)}>{projectWorkspace&&<option value="project">Project directory · save while analyzing</option>}<option value="auto">Auto · ask to save large results</option><option value="browser">Browser · compressed local index</option><option value="disk">File · save while analyzing</option></select></label>}
                    <label className="minimum-slider"><span>Minimum alignment identity <b>{Math.round(minimumIdentity * 100)}%</b></span><input type="range" min="0.45" max="0.9" step="0.01" value={minimumIdentity} onChange={(event) => setMinimumIdentity(Number(event.target.value))} /><small>V/D/J floor. C is called when its 5′-anchored alignment is significant (E ≤ 10⁻³ over the locus’ C genes), including short exact tracts.</small></label>
                    <label><span>Covered C-prefix identity gate</span><select value={constantPrefixIdentityEnabled?"on":"off"} onChange={(event)=>setConstantPrefixIdentityEnabled(event.target.value==="on")}><option value="off">Off · evidence-calibrated C calling</option><option value="on">On · require a user-set C-prefix identity</option></select><small>Optional. A failed C gate clears only C fields; V/D/J calls and boundaries are retained.</small></label>
                    {constantPrefixIdentityEnabled&&<label className="minimum-slider"><span>Minimum covered C-prefix identity <b>{Math.round(minimumConstantPrefixIdentity*100)}%</b></span><input type="range" min="0.9" max="1" step="0.01" value={minimumConstantPrefixIdentity} onChange={(event)=>setMinimumConstantPrefixIdentity(Number(event.target.value))}/><small>Counts the complete covered prefix from C-reference base 1 through the end of the matching tract, including mismatches and indels that local alignment could otherwise skip.</small></label>}
                    <p className="scientific-note calling-profile-note"><span>i</span>The Swig aligner retrieves V by seed votes over gene frames, aligns J and D exhaustively in their windows, and chooses the topology (VJ or VDJ), the alleles and both cut points of every junction together, exactly, in one placement DAG; each read is called under its own divergence, and the same model gives the posterior probabilities. J, D and C stay within the V's locus.</p>
                    <p className="scientific-note calling-profile-note"><span>i</span>{callingProfile === "igblast_compatible" ? "IgBLAST-agreement scoring: D +2/−4/−11/−1 with a 5-nt exact-run floor, J +2/−4/−13/−1, and no D-presence cost beyond beating the no-D path. It is not an IgBLAST implementation." : callingProfile === "sensitive_d" ? "R-optimized scoring with a D twice as likely a priori: the D-presence cost falls by ln 2 nats (9.45 instead of 10.55 score units for human IGH), so more short or mutated D are reported." : callingProfile === "r_optimized" ? "Default. V +2/−3/−9/−1, D +2/−3/−13/−1 with a 4-nt exact-run floor, J +2/−3/−17/−2. A D is included only when it beats the no-D path by its prior cost, derived from the D genes of the read's locus (one of N genes at one of 25 junction positions, even odds: 10.55 score units for human IGH)." : "V +2/−3/−5/−1, a 6-nt D exact-run floor and no D-presence cost beyond beating the no-D path."}</p>
                    </div>
                  </details>
                  <div className="double-d-mode-control"><label><span>Double-D / VDDJ screening</span><select value={doubleDMode} onChange={(event) => setDoubleDMode(event.target.value as DoubleDScreenMode)}><option value="off">Off · standard V(D)J only</option><option value="long_span">Long inter-V/J spans only</option><option value="all">All eligible D-bearing junctions</option></select></label></div>
                  {doubleDMode !== "off" && <details className="double-d-parameters">
                      <summary title="Conservative evidence thresholds for the opt-in VDDJ screen"><span><b>Evidence options</b><small>Conservative evidence thresholds for the opt-in VDDJ screen.</small></span></summary>
                      <div><span className="section-kicker">Rare rearrangement screen</span><h3>Two ordered D segments</h3><p>This evidence screen runs after the ordinary V(D)J call. It does not rewrite the main AIRR TSV; supported D1/D2 calls are retained in a separate evidence table and merged into the interactive detail viewer.</p></div>
                      <div className="double-d-fields">
                        {doubleDMode === "long_span" && <label><span>Minimum inter-V/J span</span><CommitNumberInput min="0" max="10000" step="1" value={doubleDMinimumSpan} onCommit={setDoubleDMinimumSpan} /><small>nt between the baseline V end and J start</small></label>}
                        <label><span>Exact seed per D</span><CommitNumberInput min="6" max="24" step="1" value={doubleDSeedLength} onCommit={setDoubleDSeedLength} /><small>nt; 11 follows the IgScout tandem-D screen</small></label>
                        <label><span>Single-D Δ trim</span><CommitNumberInput min="0" max="24" step="1" value={doubleDPseudoTrim} onCommit={setDoubleDPseudoTrim} /><small>nt removed for the pseudo-tandem test</small></label>
                        <label><span>Maximum pseudo mismatches</span><CommitNumberInput min="0" max="24" step="1" value={doubleDMaximumPseudoMismatches} onCommit={setDoubleDMaximumPseudoMismatches} /><small>reject a pair at or below this distance</small></label>
                        <label><span>Minimum two-D score gain</span><CommitNumberInput min="0" max="1000" step="1" value={doubleDMinimumScoreGain} onCommit={setDoubleDMinimumScoreGain} /><small>over the best supported single D</small></label>
                      </div>
                    </details>}
                </section>

                {!focusedWebMode && <section className={`analysis-card pipeline-card ${pipeline.enabled ? "active" : ""}`}>
                  <header><span className="card-number">04</span><div><h2>Execution mode</h2><p>Stop after annotation, or run the selected post-analysis stages sequentially.</p></div><div className="mode-toggle"><button className={!pipeline.enabled ? "active" : ""} type="button" disabled={configExporting} onClick={()=>setPipeline((current)=>({...current,enabled:false}))}>Interactive</button><button className={pipeline.enabled ? "active" : ""} type="button" disabled={configExporting} onClick={()=>setPipeline((current)=>({...current,enabled:true}))}>Pipeline</button></div></header>
                  {!pipeline.enabled ? <div className="pipeline-interactive-note"><span>Manual post-analysis</span><strong>Annotation finishes at Results.</strong><p>Collapse, chimera filtering, selection, lineage assignment, SHM, missing-allele evidence, queries, alignments, and trees remain available as explicit steps.</p></div> : <>
                    <div className="pipeline-stage-grid">
                      <article className={pipeline.collapse.enabled?"enabled":""}><label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.collapse.enabled} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,enabled:event.target.checked}}))}/><span><b>2 · Collapse / denoise</b><small>Uses policy-selected V/D/J calls when pooling is enabled.</small></span></label>{pipeline.collapse.enabled&&<div className="pipeline-fields"><label><span>Method</span><select value={pipeline.collapse.mode} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,mode:event.target.value as PipelinePlan["collapse"]["mode"]}}))}><option value="exact">Exact deduplication</option><option value="fad">FAD-compatible denoising</option><option value="conservative">Conservative indexed model</option><option value="indel">Indel-aware method D</option></select></label>{pipeline.collapse.mode==="exact"&&<label><span>Exact key</span><select value={pipeline.collapse.key} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,key:event.target.value as PipelinePlan["collapse"]["key"]}}))}><option value="sequence">Full input sequence</option><option value="trimmed">V–J-trimmed sequence</option><option value="cdr3">CDR3 nucleotide</option><option value="rearrangement">V/J calls + CDR3</option></select></label>}<label><span>Boundary</span><select value={pipeline.collapse.scope} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,scope:event.target.value as DatasetScope}}))}>{Object.entries(DATASET_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label><label><span>Unusable V/J or trim</span><select value={pipeline.collapse.unresolvedPolicy} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,unresolvedPolicy:event.target.value as "discard"|"retain"}}))}><option value="discard">Discard · default</option><option value="retain">Retain unchanged</option></select></label><label className="check-line"><input type="checkbox" checked={pipeline.collapse.respectConstantCall??true} onChange={(event)=>setPipeline((current)=>({...current,collapse:{...current.collapse,respectConstantCall:event.target.checked}}))}/><span>Separate C-gene/isotype calls</span></label></div>}</article>
                      <article className={pipeline.chimera.enabled?"enabled":""}>
                        <label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.chimera.enabled} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,enabled:event.target.checked}}))}/><span><b>3 · CHMMAIRRa filter</b><small>Consumes the collapsed working set.</small></span></label>
                        <div className="pipeline-fields"><label><span>Segment</span><select value={pipeline.chimera.segment} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,segment:event.target.value as "V"|"J"}}))}><option value="V">V</option><option value="J">J</option></select></label><label><span>Model</span><select value={pipeline.chimera.model} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,model:event.target.value as "auto"|"BW"|"DB"}}))}><option value="auto">Auto · IG BW / TR DB</option><option value="BW">Baum–Welch</option><option value="DB">Discretized Bayesian</option></select></label><label><span>Exclude posterior ≥</span><CommitNumberInput min="0" max="1" step="0.01" value={pipeline.chimera.posteriorThreshold} onCommit={(posteriorThreshold)=>setPipeline((current)=>({...current,chimera:{...current.chimera,posteriorThreshold}}))}/></label><label className="check-line"><input type="checkbox" checked={pipeline.chimera.retainUnevaluated} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,retainUnevaluated:event.target.checked}}))}/><span>Retain unevaluated</span></label><label><span>Reference MSA</span><select value={pipeline.chimera.msaSource} onChange={(event)=>setPipeline((current)=>({...current,chimera:{...current.chimera,msaSource:event.target.value as "selected"|"upload"}}))}><option value="selected">Build from assignment references</option><option value="upload">Use aligned FASTA from file</option></select></label>{pipeline.chimera.msaSource==="upload"&&<label className="pipeline-file-field"><span>Aligned FASTA</span><input type="file" accept=".fa,.fasta,.fas,.aln,.txt" onChange={(event)=>{const file=event.target.files?.[0];event.target.value="";if(!file)return;void file.text().then((text)=>{prepareReferenceMsa(text);setPipeline((current)=>({...current,chimera:{...current.chimera,uploadedMsa:text,uploadedMsaName:file.name}}));setRunError("");}).catch((error)=>setRunError(error instanceof Error?error.message:String(error)));}}/><small>{pipeline.chimera.uploadedMsaName||"No MSA selected"}</small></label>}</div>
                      </article>
                      <article className={pipeline.selection.enabled?"enabled":""}>
                        <label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.selection.enabled} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,enabled:event.target.checked}}))}/><span><b>4 · Repertoire selection</b><small>Commit study, call, CDR3, and QC filters before lineage assignment.</small></span></label>
                        <div className="pipeline-fields pipeline-selection-fields">
                          <FacetPicker label="Dataset / library" value={pipeline.selection.datasetId} items={pipelineSelectionFacets.datasets} multiple placeholder="Any dataset" onChange={(datasetId)=>setPipeline((current)=>({...current,selection:{...current.selection,datasetId}}))}/>
                          <FacetPicker label="Sample" value={pipeline.selection.sampleId} items={pipelineSelectionFacets.samples} multiple placeholder="Any sample" onChange={(sampleId)=>setPipeline((current)=>({...current,selection:{...current.selection,sampleId}}))}/>
                          <FacetPicker label="Donor / subject" value={pipeline.selection.subjectId} items={pipelineSelectionFacets.subjects} multiple placeholder="Any donor" onChange={(subjectId)=>setPipeline((current)=>({...current,selection:{...current.selection,subjectId}}))}/>
                          <FacetPicker label="Cohort" value={pipeline.selection.cohort} items={pipelineSelectionFacets.cohorts} multiple placeholder="Any cohort" onChange={(cohort)=>setPipeline((current)=>({...current,selection:{...current.selection,cohort}}))}/>
                          <FacetPicker label="Timepoint" value={pipeline.selection.timepoint} items={pipelineSelectionFacets.timepoints} multiple placeholder="Any timepoint" onChange={(timepoint)=>setPipeline((current)=>({...current,selection:{...current.selection,timepoint}}))}/>
                          <FacetPicker label="Compartment / tissue" value={pipeline.selection.compartment} items={pipelineSelectionFacets.compartments} multiple placeholder="Any compartment" onChange={(compartment)=>setPipeline((current)=>({...current,selection:{...current.selection,compartment}}))}/>
                          <FacetPicker label="Locus" value={pipeline.selection.locus} items={pipelineSelectionFacets.loci} multiple placeholder="Any locus" onChange={(locus)=>setPipeline((current)=>({...current,selection:{...current.selection,locus}}))}/>
                          <label><span>V call</span><input value={pipeline.selection.vCall} placeholder="gene or allele" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,vCall:event.target.value}}))}/></label>
                          <label className="check-line"><input type="checkbox" checked={pipeline.selection.vCallIncludeAmbiguous} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,vCallIncludeAmbiguous:event.target.checked}}))}/><span>Include ambiguous V calls</span></label>
                          <label><span>J call</span><input value={pipeline.selection.jCall} placeholder="gene or allele" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,jCall:event.target.value}}))}/></label>
                          <label className="check-line"><input type="checkbox" checked={pipeline.selection.jCallIncludeAmbiguous} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,jCallIncludeAmbiguous:event.target.checked}}))}/><span>Include ambiguous J calls</span></label>
                          <label><span>CDR3 nt</span><input value={pipeline.selection.cdr3Nt} placeholder="substring" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,cdr3Nt:event.target.value}}))}/></label>
                          <label><span>CDR3 aa</span><input value={pipeline.selection.cdr3Aa} placeholder="substring" onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,cdr3Aa:event.target.value}}))}/></label>
                          <label><span>Productive</span><select value={pipeline.selection.productive} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,productive:event.target.value as "any"|"yes"|"no"}}))}><option value="any">Either</option><option value="yes">Productive</option><option value="no">Non-productive</option></select></label>
                          <label><span>CDR3 assigned</span><select value={pipeline.selection.hasCdr3} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,hasCdr3:event.target.value as "any"|"yes"|"no"}}))}><option value="any">Either</option><option value="yes">Required</option><option value="no">Absent</option></select></label>
                          <label><span>Double-D</span><select value={pipeline.selection.doubleD} onChange={(event)=>setPipeline((current)=>({...current,selection:{...current.selection,doubleD:event.target.value as "any"|"positive"|"negative"}}))}><option value="any">Either</option><option value="positive">Positive only</option><option value="negative">Exclude positive</option></select></label>
                        </div>
                      </article>
                      <article className={`pipeline-allele-stage ${pipeline.alleleRefinement.enabled?"enabled":""}`}><label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.alleleRefinement.enabled} onChange={(event)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,enabled:event.target.checked}}))}/><span><b>1 · Repertoire allele pooling</b><small>Runs first so reassigned V/D/J calls feed collapse and every downstream stage.</small></span></label>{pipeline.alleleRefinement.enabled&&<div className="pipeline-fields"><label><span>Repertoire model</span><select value={pipeline.alleleRefinement.model} onChange={(event)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,model:event.target.value as PipelinePlan["alleleRefinement"]["model"]}}))}><option value="dirichlet">Continuous Dirichlet · legacy default</option><option value="active-set">Fast hurdle active set</option></select></label><label><span>Pooling boundary</span><select value={pipeline.alleleRefinement.scope} onChange={(event)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,scope:event.target.value as DatasetScope}}))}>{Object.entries(ALLELE_POOL_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label><label><span>Record weighting</span><select value={pipeline.alleleRefinement.weighting} onChange={(event)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,weighting:event.target.value as "unique"|"abundance"}}))}><option value="unique">Unique records · default</option><option value="abundance">duplicate_count abundance</option></select></label><label><span>Reassignment policy</span><select value={pipeline.alleleRefinement.reassignmentPolicy} onChange={(event)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,reassignmentPolicy:event.target.value as PipelinePlan["alleleRefinement"]["reassignmentPolicy"]}}))}><option value="confidence">Best posterior if confidence passes · default</option><option value="best">Best posterior for every modeled record</option></select></label>{pipeline.alleleRefinement.reassignmentPolicy==="confidence"&&<label><span>Minimum MAP confidence</span><CommitNumberInput min="0" max="1" step="0.01" value={pipeline.alleleRefinement.applyMinimumPosterior} onCommit={(applyMinimumPosterior)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,applyMinimumPosterior}}))}/></label>}<label><span>Zero-SHM neighbour odds</span><CommitNumberInput min="0" max="0.5" step="0.001" value={pipeline.alleleRefinement.baselineNeighbourOdds} onCommit={(baselineNeighbourOdds)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,baselineNeighbourOdds}}))}/></label><label><span>SHM sensitivity</span><CommitNumberInput min="0" max="10" step="0.1" value={pipeline.alleleRefinement.shmLeakageSensitivity} onCommit={(shmLeakageSensitivity)=>setPipeline((current)=>({...current,alleleRefinement:{...current.alleleRefinement,shmLeakageSensitivity}}))}/></label><fieldset className="pipeline-segment-checks"><legend>Segments</legend>{(["V","J","D"] as const).map((segment)=><label key={segment}><input type="checkbox" checked={pipeline.alleleRefinement.segments.includes(segment)} onChange={(event)=>updatePipelineAlleleSegment(segment,event.target.checked)}/><span>{segment}{segment==="D"?" · experimental":""}</span></label>)}</fieldset></div>}</article>
                      <article className={pipeline.lineage.enabled?"enabled":""}><label className="pipeline-stage-switch"><input type="checkbox" checked={pipeline.lineage.enabled} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,enabled:event.target.checked}}))}/><span><b>5 · Lineage assignment</b><small>May span samples only within the selected boundary.</small></span></label><div className="pipeline-fields"><label><span>Boundary</span><select value={pipeline.lineage.scope} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,scope:event.target.value as DatasetScope}}))}>{Object.entries(DATASET_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label><label><span>CDR3 nt identity</span><CommitNumberInput min="0.7" max="1" step="0.01" value={pipeline.lineage.identity} onCommit={(identity)=>setPipeline((current)=>({...current,lineage:{...current.lineage,identity}}))}/></label><label><span>V/J level</span><select value={pipeline.lineage.resolution} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,resolution:event.target.value as "gene"|"allele"}}))}><option value="gene">Gene</option><option value="allele">Allele</option></select></label><label><span>Multiple-call policy</span><select value={pipeline.lineage.ambiguity} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,ambiguity:event.target.value as PipelinePlan["lineage"]["ambiguity"]}}))}><option value="overlap">Any V/J overlap</option><option value="top">Top call only</option><option value="strict">Identical call sets</option></select></label><label className="check-line"><input type="checkbox" checked={pipeline.lineage.productiveOnly} onChange={(event)=>setPipeline((current)=>({...current,lineage:{...current.lineage,productiveOnly:event.target.checked}}))}/><span>Productive only</span></label></div></article>
                      <article className={pipeline.shm.enabled||pipeline.missingAlleles.enabled?"enabled":""}><div className="pipeline-diagnostic-switches"><label><input type="checkbox" checked={pipeline.shm.enabled} onChange={(event)=>setPipeline((current)=>({...current,shm:{...current.shm,enabled:event.target.checked}}))}/><span><b>6 · SHM summary</b><small>On the final working set.</small></span></label><label><input type="checkbox" checked={pipeline.missingAlleles.enabled} onChange={(event)=>setPipeline((current)=>({...current,missingAlleles:{enabled:event.target.checked},lineage:event.target.checked?{...current.lineage,enabled:true}:current.lineage}))}/><span><b>7 · Missing-allele hints</b><small>Enabling this also enables lineage assignment.</small></span></label></div>{pipeline.shm.enabled&&<div className="pipeline-fields"><label><span>SHM metric</span><select value={pipeline.shm.metric} onChange={(event)=>setPipeline((current)=>({...current,shm:{...current.shm,metric:event.target.value as PipelinePlan["shm"]["metric"]}}))}><option value="vNtRate">V nucleotide mutation rate</option><option value="vNtMutations">V nucleotide mutation count</option><option value="vAaRate">V amino-acid replacement rate</option><option value="vAaReplacements">V amino-acid replacement count</option><option value="synonymous">V synonymous mutation count</option><option value="cdrNtRate">CDR nucleotide mutation rate</option><option value="frameworkNtRate">Framework nucleotide mutation rate</option></select></label></div>}</article>
                    </div>
                    <p className="scientific-note"><span>i</span>The pipeline runs annotation → optional allele pooling and policy reassignment → collapse → chimera exclusion → repertoire selection → lineage assignment → post-lineage analyses in that order. Every stage receives the retained set from the previous stage. Targeted queries, lineage alignments, and phylogenies remain on-demand because they require a selected target.</p>
                    <div className="pipeline-config-export">
                      <div><strong>Run this setup with swig-cli</strong><small>Exports the same complete pipeline config available after analysis, including input QC, subsampling, references, annotation, and every enabled stage. If CHMMAIRRa uses the selected references, its MSA is built and embedded now.</small></div>
                      <button className={configExporting?"post-cancel":""} type="button" disabled={!configExporting&&(!activeDatasets.length||!compiled||databaseBusy||Boolean(busyCells.size)||inputInspecting)} onClick={()=>void exportConfiguredCliConfig()}>{configExporting?"Cancel config export":"Export CLI config"}</button>
                    </div>
                  </>}
                </section>}
              </div>

              <section className="run-summary launch-panel">
                <span className="section-kicker">Start analysis</span><h2>{activeInput ? focusedWebMode ? "Input ready" : `${activeDatasets.length.toLocaleString()} dataset${activeDatasets.length === 1 ? "" : "s"} ready` : "Data required"}</h2><p>{activeInput ? `${activeInputName} · ${friendlySpecies(species?.name ?? "")} · ${LOCUS_LABELS[activeScope]}${compiled ? ` · ${compiled.counts.V}/${compiled.counts.D}/${compiled.counts.J}/${compiled.counts.C} V/D/J/C references` : ""}` : "Load or paste sequence data above before starting."}</p>
                {runError && <p className="run-error" role="alert">{runError}</p>}
                <button className="analyze-button" type="button" disabled={!activeInput || !compiled || Boolean(packError) || databaseBusy || Boolean(busyCells.size) || inputInspecting} onClick={requestRun}><span>{inputInspecting?"Inspecting input…":databaseBusy || busyCells.size ? "Validating references…" : webMode === "vdj" ? "Run VDJ annotation" : webMode === "lineage" ? "Annotate + open lineage" : pipeline.enabled ? "Run annotation + pipeline" : "Analyze with SwiftIG"}</span><b>→</b></button>
                <p className="privacy-copy"><span>i</span> Query sequences, locally loaded germlines, and AIRR results are processed in this browser; Swig does not transmit them. A remotely hosted alternative database is requested from its named provider only when selected.</p>
              </section>
            </div>
          )}
        </main>
      )}

      {page === "results" && session && <ResultsPage session={session} onNewAnalysis={() => navigate("analyze")} />}
      {page==="lineage-study"&&pendingLineageStudy&&lineageStudyAirr&&<LineageStudyPage manifest={pendingLineageStudy} airrFile={lineageStudyAirr} onClose={()=>{setPendingLineageStudy(null);setLineageStudyAirr(null);setPage("analyze");}}/>}
      {leavePrompt&&<div className="output-modal-backdrop leave-guard-backdrop" role="presentation"><section className="output-modal leave-guard-modal" role="alertdialog" aria-modal="true" aria-labelledby="leave-guard-title" aria-describedby="leave-guard-description" onKeyDown={(event)=>{if(event.key==="Escape")setLeavePrompt(false);}}>
        <span className="output-direction">UNSAVED BROWSER STATE</span>
        <h2 id="leave-guard-title">Leave Swig and lose the current work?</h2>
        <p id="leave-guard-description">This page still holds loaded data, analysis results, or an active run. Leaving can discard information that has not been exported or saved as a portable session.</p>
        <div className="output-modal-actions"><button className="output-save-primary" type="button" autoFocus onClick={()=>setLeavePrompt(false)}><span>Stay on this page</span><small>Keep every in-browser result and setting</small><b>Cancel</b></button><button className="leave-guard-confirm" type="button" onClick={confirmPageDeparture}><span>Leave anyway</span><small>Discard unsaved browser state</small></button></div>
      </section></div>}
      {pendingGzipImport&&<div className="output-modal-backdrop" role="presentation"><section className="output-modal gzip-member-modal" role="dialog" aria-modal="true" aria-labelledby="gzip-member-title">
        <button className="output-modal-close" type="button" onClick={cancelInputInspection} aria-label="Cancel this input import">×</button>
        <span className="output-direction">CONCATENATED GZIP · SAMPLE STRUCTURE</span>
        <h2 id="gzip-member-title"><em>{pendingGzipImport.inputName}</em> contains {pendingGzipImport.members.length.toLocaleString()} independently readable {pendingGzipImport.format} files.</h2>
        <p>Each gzip member begins at a complete {pendingGzipImport.format} record boundary. Choose whether those members represent one biological sample or separate samples.</p>{webMode!=="repertoire"&&<p className="scientific-note"><span>i</span>Separate samples need editable sample names and metadata. Choosing that option will move this input into the end-to-end repertoire workflow; merging keeps the current focused workflow.</p>}
        <div className="gzip-member-summary">{pendingGzipImport.members.map((member,index)=><span key={`${member.start}-${member.end}`}><b>{index+1} · {pendingGzipImport.suggestedNames[index]}</b><small>{bytes(member.compressedBytes)} compressed{member.uncompressedBytes===null?"":` · ${bytes(member.uncompressedBytes)} decoded`}</small><code title={member.firstLine}>{member.firstRecordId||member.firstLine||"empty member"}</code></span>)}</div>
        <div className="output-modal-actions"><button className="output-save-primary" type="button" onClick={()=>resolveGzipImportChoice("separate")}><span>Import as {pendingGzipImport.members.length} separate samples</span><small>Independent sample, donor, cohort, timepoint, and tissue fields</small><b>Separate →</b></button><button type="button" onClick={()=>resolveGzipImportChoice("merge")}><span>Merge into one sample</span><small>Decode every member in order as one continuous dataset</small></button><button className="post-cancel" type="button" onClick={cancelInputInspection}>Cancel import</button></div>
        <p className="output-safety"><span>i</span>All metadata remain editable in the dataset table. Neither choice drops records or changes sequence content.</p>
      </section></div>}
      {pendingDirectoryInputs&&<div className="output-modal-backdrop" role="presentation"><section className="output-modal directory-donor-modal" role="dialog" aria-modal="true" aria-labelledby="directory-donor-title"><button className="output-modal-close" type="button" onClick={()=>setPendingDirectoryInputs(null)} aria-label="Cancel directory loading">×</button><span className="output-direction">DIRECTORY · DONOR METADATA</span><h2 id="directory-donor-title">Do files directly inside {pendingDirectoryInputs.flatRoots.length===1?<em>{pendingDirectoryInputs.flatRoots[0]}</em>:"these directories"} come from the same donor?</h2><p>{pendingDirectoryInputs.inputs.length.toLocaleString()} compatible dataset file{pendingDirectoryInputs.inputs.length===1?" was":"s were"} found. This choice only initializes <b>Donor / subject</b>; every value remains editable in the dataset table.</p><div className="directory-donor-summary">{pendingDirectoryInputs.flatRoots.map((root)=><span key={root}><b>{root}</b><small>{pendingDirectoryInputs.inputs.filter((input)=>input.directoryRoot===root&&!input.nestedDirectoryDonor).length} direct file(s)</small></span>)}</div><div className="output-modal-actions"><button className="output-save-primary" type="button" onClick={()=>commitDirectoryInputs(true)}><span>Same donor within each directory</span><b>Group →</b></button><button type="button" onClick={()=>commitDirectoryInputs(false)}><span>Different donor for each file</span><small>Keep separate initial donor IDs</small></button></div><p className="output-safety"><span>i</span>Nested directories are assigned automatically from the first folder beneath the selected root.</p></section></div>}
      {pendingLineageStudy&&!lineageStudyAirr&&<div className="output-modal-backdrop" role="presentation"><section className="output-modal session-load-modal" role="dialog" aria-modal="true" aria-labelledby="lineage-study-load-title"><button className="output-modal-close" type="button" onClick={()=>{setPendingLineageStudy(null);setSessionLoadError("");}}>×</button><span className="output-direction">LINEAGE STUDY · LAZY LINKED AIRR</span><h2 id="lineage-study-load-title">Select the lineage-sorted AIRR file.</h2><p>Swig will load the {pendingLineageStudy.summaries.length.toLocaleString()} lineage summaries now, but it will not scan or index the linked AIRR table. Only the byte range for a lineage you open is read.</p><div className="output-flow"><div><span>Lineage manifest</span><strong>{pendingLineageStudy.analysis.inputName}</strong><small>{pendingLineageStudy.summaries.length.toLocaleString()} selectable lineages</small></div><b>+</b><div className="destination"><span>Required uncompressed file</span><strong>{pendingLineageStudy.linkedAirr.name}</strong><small>{pendingLineageStudy.linkedAirr.size.toLocaleString()} bytes · {pendingLineageStudy.linkedAirr.records.toLocaleString()} lineage rows</small></div></div><button className="output-save-primary" type="button" onClick={()=>lineageStudyAirrInputRef.current?.click()}><span>Choose lineage AIRR TSV</span><b>Open lazily →</b></button>{sessionLoadError?<p className="run-error" role="alert">{sessionLoadError}</p>:null}</section></div>}
      {(pendingLoadedSession||(sessionLoadError&&!pendingLineageStudy)||loadingSession)&&<div className="output-modal-backdrop" role="presentation"><section className="output-modal session-load-modal" role="dialog" aria-modal="true" aria-labelledby="session-load-title"><button className="output-modal-close" type="button" onClick={()=>{if(loadingSession)cancelSessionRestore();else{setPendingLoadedSession(null);setSessionLoadError("");}}}>×</button><span className="output-direction">{pendingLoadedSession?"SESSION · LINKED AIRR DATA":"PROJECT DIRECTORY · RESTORE"}</span><h2 id="session-load-title">{pendingLoadedSession?"Select the AIRR table linked to this session.":loadingSession?"Restoring the active project run.":"Saved state could not be loaded."}</h2>{pendingLoadedSession?<><p>The session stores references, options, masks, counts, lineage assignments, plots, and sparse double-D evidence. It deliberately does not duplicate the main AIRR table.</p><div className="output-flow"><div><span>Saved analysis</span><strong>{pendingLoadedSession.analysis.inputName}</strong><small>{pendingLoadedSession.linkedAirr.records.toLocaleString()} records · fingerprint {pendingLoadedSession.linkedAirr.fingerprint.slice(0,12)}…</small></div><b>+</b><div className="destination"><span>Required linked file</span><strong>{pendingLoadedSession.linkedAirr.name}</strong><small>AIRR TSV or TSV.gz; content is verified before restoration</small></div></div>{!loadingSession&&<button className="output-save-primary" type="button" onClick={()=>linkedAirrInputRef.current?.click()}><span>Choose linked AIRR TSV</span><b>Open →</b></button>}</>:loadingSession?<p>The active run's AIRR table is being verified and indexed from the selected directory.</p>:null}{loadingSession&&<div className="post-progress"><div><span>{sessionLoadProgress.stage}</span><span className="post-progress-actions"><strong>{sessionLoadProgress.total?`${Math.min(100,sessionLoadProgress.records/sessionLoadProgress.total*100).toFixed(1)}%`:"working"}</strong><button type="button" onClick={cancelSessionRestore}>Cancel</button></span></div><progress max={Math.max(1,sessionLoadProgress.total)} value={sessionLoadProgress.records}/><small>{sessionLoadProgress.records.toLocaleString()} / {sessionLoadProgress.total.toLocaleString()} records</small></div>}{sessionLoadError?<p className="run-error" role="alert">{sessionLoadError}</p>:null}</section></div>}
      {outputPrompt && activeInput && <div className="output-modal-backdrop" role="presentation"><section className="output-modal" role="dialog" aria-modal="true" aria-labelledby="output-dialog-title">
        <button className="output-modal-close" type="button" onClick={() => setOutputPrompt(false)} aria-label="Cancel output selection">×</button>
        <span className="output-direction">OUTPUT · SAVE RESULTS</span>
        <h2 id="output-dialog-title">Choose where Swig will <em>write the AIRR output.</em></h2>
        <p>The next system window is a <b>Save As</b> dialog. You are naming a new results file—this is not another sequence import.</p>
        <div className="output-flow"><div><span>Input already selected</span><strong>{activeDatasets.length === 1 ? activeInput.name : `${activeDatasets.length} datasets · ${activeInputName}`}</strong></div><b>→</b><div className="destination"><span>New output file</span><strong>{outputName(activeInputName)}</strong><small>Written incrementally during analysis</small></div></div>
        <div className="output-modal-actions"><button className="output-save-primary" type="button" onClick={() => void run("disk")}><span>Choose output file &amp; start</span><b>Save AIRR →</b></button><button type="button" onClick={() => void run("browser")}><span>Keep output in browser instead</span><small>Compressed local index; download after the run</small></button></div>
        <p className="output-safety"><span>i</span> Query sequences remain in this browser and are not transmitted by Swig.</p>
      </section></div>}
      {page !== "home" && <footer className="site-footer"><Brand /><p>Swig {APP_VERSION} · SwiftIG WebAssembly interface · research software · validate study-critical calls independently.</p><div><a href="./METHODS_INDEX.md" target="_blank" rel="noreferrer">Methods ↗</a><a href="https://github.com/MurrellGroup/swiftig" target="_blank" rel="noreferrer">Source ↗</a><a href="./references/README.md" target="_blank" rel="noreferrer">IMGT data notice ↗</a><a href="https://docs.airr-community.org/" target="_blank" rel="noreferrer">AIRR ↗</a></div></footer>}
    </div>
  );
}
