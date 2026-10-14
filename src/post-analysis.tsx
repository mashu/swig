import {
  useEffect,
  useMemo,
  useRef,
  useState,
  type ChangeEvent,
  type MutableRefObject,
  type ReactNode,
} from "react";

import type { FastTreeRun } from "./biowasm-runtime";
import { projectCodonAlignment, translateCodonSequence } from "./alignment-model";
import { runCodonAwareKalignTask, runFastTreeTask, runKalignTask } from "./biowasm-task-runtime";
import { createAlivibeMsaJob, runAlivibeMsaTask } from "./alivibe-msa-runtime";
import {
  ALIVIBE_BRIDGE_VERSION,
  ALIVIBE_SOURCE_REVISION,
  assertAlivibeInitialLoad,
  assertAlivibeRoundTripTarget,
  getAlivibeBridge,
  loadAlivibeNucleotideFasta,
  readAlivibeNucleotideFasta,
  type AlivibeEditorWindow,
} from "./alivibe-roundtrip";
import { inspectAlignment, validateCorrectedAlignment } from "./alignment-provenance";
import { chimeraVisiblePositions, classifyChimeraQuerySite } from "./chimera-view-model";
import { CommitNumberInput } from "./commit-number-input";
import { CommitTextInput, CommitTextarea } from "./commit-text-input";
import { FacetPicker, uniqueFacetItems } from "./facet-picker";
import { callFacetItems } from "./call-facets";
import {
  runChmmairra,
  runChmmairraDetail,
  writeChmmairra,
  writeChmmairraTsv,
  type ChmmDashboard,
  type ChmmDetail,
  type ChmmRunOptions,
  type ChmmSegment,
} from "./chmmairra-runtime";
import { alignmentExtension, alignmentText, tableExtension, tableHeader, tableRow, treeNexus, type AlignmentExportFormat, type TableExportFormat } from "./export-formats";
import { compressedExtension, compressedMime, compressedName, outputWriter, type OutputCompression } from "./gzip-output";
import { MissingAlleleAccumulator, DEFAULT_MISSING_ALLELE_OPTIONS, type MissingAlleleDashboard, type MissingAlleleOptions } from "./germline-evidence";
import { DEFAULT_PERSONALIZED_GERMLINE_OPTIONS, type PersonalizedGermlineDashboard, type PersonalizedGermlineOptions } from "./personalized-germline";
import { UnifiedGermlineRuntime } from "./unified-germline-runtime";
import { DEFAULT_UNIFIED_OPTIONS, type UnifiedDashboard, type UnifiedOptions } from "./unified-germline";
import { UnifiedGermlineResultsPanel } from "./unified-germline-panel";
import { PersonalizedGermlineRuntime } from "./personalized-germline-runtime";
import { GERMLINE_OUTGROUP, alignedSequenceFrameOffset, inferLineageGermline, isProductiveLineageRow, lineageInputFasta, lineageTemplateChoices, projectLineageTemplate, restoreLineageOuterFlanks, retainLineageTemplate, type LineageGermlineMethod } from "./lineage-alignment";
import {
  buildLineageGermlineSketchIndex,
  scoreGermlineCandidate,
  screenLineageGermlineCandidates,
  type GermlineNeighbourScore,
  type LineageGermlineSketchIndex,
} from "./lineage-neighbours";
import { LineageTreeViewer } from "./lineage-tree-viewer";
import { inferAlignedReadingFrame, translateAlignedNucleotides, type AlignmentFrameOffset } from "./lineage-phylogeny";
import type { CallingProfile } from "./aligner-engine";
import { withAnalysisWebLock } from "./aligner-runtime";
import {
  canonicalizeTree,
  collapseShortInternalBranches,
  parseNewick,
  rootOnOutgroup,
  serializeNewick,
} from "./phylogeny";
import {
  parseFasta,
  prepareReferenceMsa,
  lineageDoubleDMatches,
  type AmbiguityPolicy,
  type CallResolution,
  type CollapseMode,
  type DedupKey,
  type LineageSummary,
  type LineageNeighbourHit,
  type LineageNeighbourResult,
  type LineageDoubleDFilter,
  type QueryHit,
  type QueryConstraint,
  type QueryMetric,
  type QueryTarget,
} from "./post-analysis-core";
import {
  PostAnalysisRuntime,
  type DedupDashboard,
  type LineageDashboard,
} from "./post-analysis-runtime";
import { inferQueryAssignments, type InferredQueryAssignment } from "./query-inference-runtime";
import type { CompiledReferences, ScopeKey } from "./reference-pack";
import type { AirrDetailRow, AirrIndexRecord, AirrResultStore, FacetValue, ResultFacets } from "./result-store";
import { ColoredSequence, sequenceColor } from "./sequence-colors";
import { MissingAlleleResultsPanel, PersonalizedGermlineResultsPanel, ShmResultsPanel } from "./post-analysis-extensions";
import { DEFAULT_REPERTOIRE_SELECTION, selectRepertoire, validateRepertoireSelection, type RepertoireSelectionOptions, type RepertoireSelectionResult } from "./repertoire-selection";
import { ShmAccumulator, type ShmDashboard, type ShmMetricKey } from "./shm-analysis";
import { packSessionVector, unpackSessionVector, type PostAnalysisSessionSnapshot } from "./session-state";
import { DATASET_SCOPE_LABELS, datasetScopeValue, type DatasetManifestEntry, type DatasetScope, type PipelinePlan } from "./study-design";
import { sampleColor, type SampleColorMap } from "./sample-colors";
import { PhyloUcaPanel, type PhyloUcaPanelState } from "./phylo-uca/panel";
import { prepareUploadedLineageTree } from "./uploaded-lineage-tree";
import { AlleleRefinementPanel } from "./allele-refinement/panel";
import { AlleleRefinementRuntime } from "./allele-refinement/runtime";
import { refinedCall, refineDetailRows, modelSummaryTable, writeRefinedAirr, writeRefinementSidecar } from "./allele-refinement/export";
import { restoreAlleleRefinement, saveAlleleRefinement } from "./allele-refinement/serialization";
import { DEFAULT_ALLELE_REFINEMENT_OPTIONS, type AlleleReassignmentPolicy, type AlleleRefinementOptions, type AlleleRefinementResult } from "./allele-refinement/types";

interface Props {
  store: AirrResultStore;
  references: CompiledReferences;
  scope: ScopeKey;
  loci: FacetValue[];
  resultFacets: ResultFacets;
  inputName: string;
  workers: number;
  callingProfile: CallingProfile;
  minimumIdentity: number;
  strand: 0 | 1 | 2;
  datasets?: DatasetManifestEntry[];
  sampleColors?: SampleColorMap;
  defaultCollapseScope?: DatasetScope;
  defaultLineageScope?: DatasetScope;
  doubleDCount?: number;
  autoPipeline?: PipelinePlan | null;
  sidebarTools?: ReactNode;
  onInspect: (ordinal: number) => void;
  /** Signals a committed or user-visible state change for optional project checkpoints. */
  onSessionChange?: (reason: string) => void;
  sessionHandleRef?: MutableRefObject<PostAnalysisSessionHandle | null>;
  initialSession?: PostAnalysisSessionSnapshot | null;
  /** Treat the complete assigned input as one preselected lineage. */
  directLineage?: boolean;
}

export interface PostAnalysisSessionHandle {
  snapshot: (signal?: AbortSignal) => Promise<PostAnalysisSessionSnapshot>;
}

interface SaveFileHandle {
  createWritable: () => Promise<{ write: (value: string | Blob | Uint8Array) => Promise<void>; close: () => Promise<void>; abort?: (reason?: unknown) => Promise<void> }>;
}

interface ChartDatum {
  label: string;
  value: number;
}

interface WorkingSetStage {
  id: "dedup" | "chimera" | "selection";
  label: string;
  input: number;
  retained: number;
  discarded: number;
  detail: string;
}

interface TreeSnapshot extends FastTreeRun {
  rootedNewick: string;
  stableNewick: string;
  collapsedEdges: number;
  collapseThreshold: number;
  source: string;
  frameOffset: AlignmentFrameOffset;
  origin?: "fasttree" | "uploaded";
  /** Uploaded trees containing the synthetic germline guide are display-only. */
  observedOnly?: boolean;
}

type TreeViewMode = "stable" | "rooted" | "raw";
type PostModuleId = "dedup" | "chimera" | "selection" | "alleles" | "lineage" | "diagnostics" | "workbench" | "query";
type PostWorkspaceId = "overview" | PostModuleId;

const POST_MODULE_ORDER: PostModuleId[] = ["alleles", "dedup", "chimera", "selection", "lineage", "diagnostics", "workbench", "query"];

interface EditedAlignmentState {
  key: string;
  lineageIds: number[];
  productiveOnly?: boolean;
  fasta: string;
  source: string;
  frameOffset?: AlignmentFrameOffset;
  savedAt: string;
}

interface AlivibeRoundTripSession {
  token: string;
  popup: AlivibeEditorWindow;
  baseline: string;
  baselineFingerprint: string;
  lineageIds: number[];
  groupKey: string;
  frameOffset: AlignmentFrameOffset;
}

interface LineageMergeState {
  id: string;
  label: string;
  originalLineageIds: number[];
  createdAt: string;
}

interface CombinedNeighbourHit {
  lineageId: number;
  cdr3?: LineageNeighbourHit;
  germline?: GermlineNeighbourScore;
}

interface QueryLineageHit {
  lineageId: number;
  bestScore: number;
  bestDistance: number;
  bestOrdinal: number;
  matchedSequences: number;
  matchedQueries: number;
}

function validAlignmentFrameOffset(value: unknown): AlignmentFrameOffset | undefined {
  return value === 0 || value === 1 || value === 2 ? value : undefined;
}

type LineageSortKey = "id" | "studyGroup" | "abundance" | "uniqueMembers" | "doubleDPositiveMembers" | "sampleCount" | "locus" | "vCalls" | "jCalls" | "cdr3Length" | "shmMean" | "shmMaximum" | "shmP95";
type LineageSampleFilterMode = "any" | "multiple" | "single" | "selected-any" | "selected-all" | "selected-only";
type LineageShmFilterStatistic = "mean" | "maximum" | "p95";

function lineageGroupKey(lineageIds: Iterable<number>): string {
  return [...new Set(lineageIds)].filter((value)=>value>0).sort((left,right)=>left-right).join("+");
}

function lineageAlignmentKey(lineageIds: Iterable<number>, productiveOnly = false): string {
  const group = lineageGroupKey(lineageIds);
  return group && productiveOnly ? `${group}|productive` : group;
}

function stratifiedLineageRows(rows: AirrDetailRow[], lineageByOrdinal: ReadonlyMap<number, number>, limit: number): AirrDetailRow[] {
  if (rows.length <= limit) return rows;
  const groups = new Map<number, AirrDetailRow[]>();
  for (const row of rows) {
    const lineageId = lineageByOrdinal.get(row.record.ordinal) ?? 0;
    const values = groups.get(lineageId);
    if (values) values.push(row); else groups.set(lineageId, [row]);
  }
  const ordered = [...groups.values()];
  const offsets = new Uint32Array(ordered.length);
  const selected: AirrDetailRow[] = [];
  while (selected.length < limit) {
    let added = false;
    for (let group = 0; group < ordered.length && selected.length < limit; group += 1) {
      const row = ordered[group][offsets[group]++];
      if (!row) continue;
      selected.push(row);
      added = true;
    }
    if (!added) break;
  }
  return selected;
}

function baseName(name: string): string {
  return name.replace(/(\.airr)?\.(tsv|csv|txt|fa|fasta|fastq)(\.gz)?$/i, "") || "swig";
}

function downloadBlob(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1_000);
}

async function saveStream(
  name: string,
  description: string,
  extension: string,
  produce: (writer: { write: (value: string | Blob | Uint8Array) => Promise<void> }) => Promise<void>,
  compression: OutputCompression = "none",
) {
  const outputName=compressedName(name,compression);
  const outputExtension=compressedExtension(extension,compression);
  const picker = (window as Window & { showSaveFilePicker?: (options: unknown) => Promise<SaveFileHandle> }).showSaveFilePicker;
  if (picker) {
    const handle = await picker.call(window, {
      suggestedName: outputName,
      types: [{ description:compression==="gzip"?`Gzip-compressed ${description}`:description, accept: { [compressedMime("text/plain",compression)]: [outputExtension] } }],
    });
    const writable = outputWriter(await handle.createWritable(),compression);
    try {
      await produce(writable);
      await writable.close();
    } catch (error) {
      await writable.abort(error);
      throw error;
    }
    return;
  }
  const parts: BlobPart[] = [];
  const writable=outputWriter({write:async(value)=>{parts.push(value instanceof Uint8Array?value.slice().buffer:value);}},compression);
  try{await produce(writable);await writable.close();}
  catch(error){await writable.abort(error);throw error;}
  downloadBlob(new Blob(parts, { type: compressedMime("text/plain;charset=utf-8",compression) }), outputName);
}

function downloadText(value: string, name: string, type = "text/plain;charset=utf-8") {
  downloadBlob(new Blob([value], { type }), name);
}

function saveSvg(svg: SVGSVGElement | null, name: string) {
  if (!svg) return;
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  downloadText(new XMLSerializer().serializeToString(clone), name, "image/svg+xml;charset=utf-8");
}

function formatNumber(value: number): string {
  return value >= 1_000_000 ? `${(value / 1_000_000).toFixed(1)}m` : value >= 1_000 ? `${(value / 1_000).toFixed(1)}k` : value.toLocaleString();
}

function BarChart({ title, subtitle, data, color, name, controls }: {
  title: string;
  subtitle: string;
  data: ChartDatum[];
  color: string;
  name: string;
  controls?: ReactNode;
}) {
  const svgRef = useRef<SVGSVGElement>(null);
  const visible = data.slice(0, 24);
  const maximum = Math.max(1, ...visible.map((item) => item.value));
  const height = Math.max(210, visible.length * 25 + 52);
  return <article className="post-chart-card">
    <header><div><span className="section-kicker">Repertoire summary</span><h3>{title}</h3><p>{subtitle}</p></div><div className="post-chart-actions">{controls}<button type="button" onClick={() => {let value=tableHeader(["label","value"],"csv");for(const item of data)value+=tableRow(["label","value"],{label:item.label,value:item.value},"csv");downloadText(value,name.replace(/\.svg$/i,".csv"),"text/csv;charset=utf-8");}}>Data CSV ↓</button><button type="button" onClick={() => saveSvg(svgRef.current, name)}>SVG ↓</button></div></header>
    <div className="post-chart-scroll"><svg ref={svgRef} role="img" aria-label={title} viewBox={`0 0 760 ${height}`}>
      <rect width="760" height={height} fill="#fffdf8" />
      {visible.map((item, index) => {
        const y = 34 + index * 25;
        const width = item.value / maximum * 450;
        return <g key={`${item.label}-${index}`}>
          <text x="8" y={y + 14} fontFamily="Inter,Arial,sans-serif" fontSize="11" fill="#31413c">{item.label.length > 25 ? `${item.label.slice(0, 23)}…` : item.label}</text>
          <rect x="190" y={y} width="450" height="16" rx="3" fill="#e7ebe5" />
          <rect x="190" y={y} width={Math.max(1, width)} height="16" rx="3" fill={color} />
          <text x={Math.min(704, 200 + width)} y={y + 13} fontFamily="Inter,Arial,sans-serif" fontSize="10" fontWeight="700" fill="#16231f">{formatNumber(item.value)}</text>
        </g>;
      })}
      {!visible.length && <text x="380" y="110" textAnchor="middle" fontFamily="Inter,Arial,sans-serif" fontSize="13" fill="#6d7975">No values</text>}
    </svg></div>
  </article>;
}

function AlignmentPreview({ fasta, mode, frameOffset }: { fasta: string; mode: "nt" | "aa"; frameOffset: AlignmentFrameOffset }) {
  const [showNames, setShowNames] = useState(true);
  const records = parseFasta(fasta, true);
  return <div className="lineage-alignment-preview-shell">
    <div className="lineage-plot-display-controls"><label className="check-line"><input type="checkbox" checked={showNames} onChange={(event) => setShowNames(event.target.checked)} /><span>Show sequence names</span></label><small>Hidden-name space is reassigned to the alignment viewport.</small></div>
    <div className={`lineage-alignment-preview${showNames ? "" : " names-hidden"}`}>
      <div className="alignment-ruler">{showNames && <span>Name</span>}<span>{mode === "nt" ? "Nucleotide alignment" : "Codon translation"} · showing {Math.min(80, records.length)} of {records.length}</span></div>
      {records.slice(0, 80).map((record) => <div className={record.name === GERMLINE_OUTGROUP ? "germline-row" : ""} key={record.name}>{showNames && <strong title={record.name}>{record.name}</strong>}<ColoredSequence sequence={mode === "nt" ? record.sequence : translateAlignedNucleotides(record.sequence, frameOffset)} alphabet={mode} /></div>)}
    </div>
  </div>;
}

function ChimeraHighlighter({ detail, name, onInspect }: { detail: ChmmDetail; name: string; onInspect: (ordinal: number) => void }) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [colorMode, setColorMode] = useState<"nucleotide" | "highlighter">("highlighter");
  const [eventIndex, setEventIndex] = useState(0);
  useEffect(() => setEventIndex(0), [detail.ordinal]);
  const cellWidth = 12;
  const labelWidth = 210;
  const rowHeight = 27;
  const top = 45;
  const selectedEvent = detail.recombinations[Math.min(eventIndex, Math.max(0, detail.recombinations.length - 1))];
  const parentIndex = (parentName: string | undefined, fallback: number) => {
    const found = detail.parents.findIndex((parent) => parent.name === parentName);
    return found >= 0 ? found : Math.min(fallback, Math.max(0, detail.parents.length - 1));
  };
  const parentAIndex = parentIndex(selectedEvent?.left, 0);
  let parentBIndex = parentIndex(selectedEvent?.right, parentAIndex === 0 ? 1 : 0);
  if (parentBIndex === parentAIndex && detail.parents.length > 1) parentBIndex = parentAIndex === 0 ? 1 : 0;
  const parentA = detail.parents[parentAIndex] ?? { name: "Parent A unavailable", sequence: "" };
  const parentB = detail.parents[parentBIndex] ?? { name: "Parent B unavailable", sequence: "" };
  const query = detail.threadedObservation;
  const visiblePositions = chimeraVisiblePositions(query, parentA.sequence, parentB.sequence);
  const hiddenTripleGaps = query.length - visiblePositions.length;
  const rows = [
    { kind: "parent" as const, label: `PARENT A · ${parentA.name}`, sequence: parentA.sequence },
    { kind: "query" as const, label: `QUERY · ${detail.sequenceId}`, sequence: query },
    { kind: "parent" as const, label: `PARENT B · ${parentB.name}`, sequence: parentB.sequence },
  ];
  const alignmentBottom = top + (rows.length - 1) * rowHeight + rowHeight - 16;
  const breakpointLabelStart = alignmentBottom + 22;
  const displayBoundary = (position: number) => {
    const original = Math.max(0, position - 1);
    const index = visiblePositions.findIndex((value) => value >= original);
    return index < 0 ? visiblePositions.length : index;
  };
  const highlighterColor = (position: number) => {
    const category = classifyChimeraQuerySite(query[position] ?? "-", parentA.sequence[position] ?? "-", parentB.sequence[position] ?? "-");
    return category === "parent_a" ? "#d84a4a" : category === "parent_b" ? "#4277c7" : "#cfd3d1";
  };
  const width = labelWidth + visiblePositions.length * cellWidth + 18;
  const height = breakpointLabelStart + Math.max(1, detail.recombinations.length) * 10 + 24;
  const informative = [...detail.threadedObservation].filter((base) => base !== "-" && base !== "N").length;
  return <section className="chimera-detail" tabIndex={-1}>
    <header><div><span className="section-kicker">On-demand Viterbi reconstruction · record #{(detail.ordinal + 1).toLocaleString()}</span><h4>{detail.sequenceId}</h4><p>{detail.segment} call {detail.call || "unassigned"} · posterior {(detail.probability * 100).toFixed(3)}% · DFR {detail.dfr} · {informative} informative MSA columns · {hiddenTripleGaps} triple-gap sites hidden</p></div><div className="post-chart-actions"><button type="button" onClick={() => onInspect(detail.ordinal)}>Open AIRR record</button><button type="button" onClick={() => saveSvg(svgRef.current, name)}>SVG ↓</button></div></header>
    <div className="chimera-viz-controls"><div className="mode-toggle"><button className={colorMode === "nucleotide" ? "active" : ""} type="button" onClick={() => setColorMode("nucleotide")}>Nucleotide colors</button><button className={colorMode === "highlighter" ? "active" : ""} type="button" onClick={() => setColorMode("highlighter")}>Parent-match highlighter</button></div>{detail.recombinations.length > 1 && <label><span>Parent pair / breakpoint</span><select value={eventIndex} onChange={(event) => setEventIndex(Number(event.target.value))}>{detail.recombinations.map((event, index) => <option value={index} key={`${event.position}-${index}`}>{index + 1}. {event.left} → {event.right} at MSA {event.position}</option>)}</select></label>}<p>Tile colors are computed from nucleotide identity only; the Viterbi path is never used for coloring.</p></div>
    {colorMode === "highlighter" && <div className="chimera-breakpoint-strip"><span className="parent-a-key"><i />red · query matches parent A only</span><span className="parent-b-key"><i />blue · query matches parent B only</span><span className="neutral-key"><i />gray · matches both or neither</span></div>}
    <div className="chimera-alignment-scroll"><svg ref={svgRef} width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label="Viterbi parent highlighter nucleotide alignment">
      <rect width={width} height={height} fill="#fbfaf5" />
      <text x="12" y="22" fontFamily="Inter,Arial,sans-serif" fontSize="10" fontWeight="700" fill="#34433e">CHMMAIRRa parent-match alignment · {colorMode === "nucleotide" ? "standard nucleotide palette" : "literal query/parent matches"}</text>
      {visiblePositions.filter((_, index) => index % 10 === 0).map((position, tick) => <g key={position}><line x1={labelWidth + tick * 10 * cellWidth} x2={labelWidth + tick * 10 * cellWidth} y1={top - 11} y2={alignmentBottom} stroke={tick % 5 === 0 ? "#9aa7a1" : "#d9ddd8"} strokeWidth={tick % 5 === 0 ? 1 : 0.5} /><text x={labelWidth + tick * 10 * cellWidth + 2} y={top - 16} fontFamily="ui-monospace,monospace" fontSize="8" fill="#6b7974">{position + 1}</text></g>)}
      {rows.map((row, rowIndex) => {
        const y = top + rowIndex * rowHeight;
        return <g key={`${row.kind}-${row.label}`}>
          <rect x="0" y={y - 14} width={labelWidth - 5} height={rowHeight - 2} fill={row.kind === "query" ? "#172622" : "#e9ece7"} />
          <text x="10" y={y + 3} fontFamily="ui-monospace,monospace" fontSize="8" fontWeight="700" fill={row.kind === "query" ? "#f4f6f3" : "#384842"}>{row.label.length > 31 ? `${row.label.slice(0, 29)}…` : row.label}<title>{row.label}</title></text>
          {visiblePositions.map((position, displayIndex) => {
            const base = row.sequence[position] ?? "-";
            const fill = colorMode === "nucleotide" ? sequenceColor(base, "nt") : row.kind === "query" ? highlighterColor(position) : "#f7f6f1";
            const dark = fill === "#d84a4a" || fill === "#4277c7";
            return <g key={position}>
              <rect x={labelWidth + displayIndex * cellWidth} y={y - 14} width={cellWidth} height={rowHeight - 2} fill={fill} />
              <text x={labelWidth + displayIndex * cellWidth + cellWidth / 2} y={y + 4} textAnchor="middle" fontFamily="ui-monospace,monospace" fontSize="10" fontWeight="700" fill={dark ? "#ffffff" : "#17231f"}>{base}</text>
            </g>;
          })}
        </g>;
      })}
      {detail.recombinations.map((event, index) => {
        const x = labelWidth + displayBoundary(event.position) * cellWidth;
        const selected = index === eventIndex || detail.recombinations.length === 1;
        return <g key={`${event.position}-${index}`}><line x1={x} x2={x} y1={top - 12} y2={alignmentBottom + 7} stroke={selected ? "#7c2435" : "#8c9692"} strokeWidth={selected ? 2.4 : 1.2} strokeDasharray="4 3" /><text x={x + 4} y={breakpointLabelStart + index * 10} fontFamily="ui-monospace,monospace" fontSize="8" fontWeight={selected ? "700" : "500"} fill={selected ? "#84253a" : "#66736e"}>Viterbi break {event.position}: {event.left.replace(/\s.*/, "")} → {event.right.replace(/\s.*/, "")}</text></g>;
      })}
      <text x="12" y={height - 11} fontFamily="Inter,Arial,sans-serif" fontSize="8" fill="#65736e">Pure A/query/B gap columns are hidden. Dashed lines are Viterbi breakpoints; tile colors are independent of the inferred path.</text>
    </svg></div>
  </section>;
}

function parseQueries(text: string): string[] {
  if (text.trimStart().startsWith(">")) return parseFasta(text).map((record) => record.sequence).filter(Boolean);
  return text.split(/\r?\n/).map((line) => line.trim()).filter((line) => line && !line.startsWith("#")).map((line) => line.replace(/\s/g, ""));
}

function isAbortError(error: unknown): boolean {
  return typeof error === "object" && error !== null && "name" in error && (error as { name?: unknown }).name === "AbortError";
}

export function PostAnalysisWorkbench({ store, references, scope, loci, resultFacets, inputName, workers, callingProfile, minimumIdentity, strand, datasets = [], sampleColors = {}, defaultCollapseScope = "sample", defaultLineageScope = "sample", doubleDCount = 0, autoPipeline, sidebarTools, onInspect, onSessionChange, sessionHandleRef, initialSession, directLineage = false }: Props) {
  const runtime = useMemo(() => new PostAnalysisRuntime(store, workers), [store, workers]);
  const alleleRuntime = useMemo(() => new AlleleRefinementRuntime(), [store]);
  const unifiedGermlineRuntime = useMemo(() => new UnifiedGermlineRuntime(), [store]);
  const personalizedGermlineRuntime = useMemo(() => new PersonalizedGermlineRuntime(), [store]);
  const postLockAbortRef = useRef<AbortController | null>(null);
  const cancellationRecoveryRef = useRef<Promise<void> | null>(null);
  const pipelineActiveRef = useRef(false);
  const postRestoreActiveRef = useRef(false);
  useEffect(() => () => {
    postLockAbortRef.current?.abort();
    runtime.terminate();
    alleleRuntime.terminate();
    personalizedGermlineRuntime.terminate();
    unifiedGermlineRuntime.terminate();
  }, [alleleRuntime, personalizedGermlineRuntime, unifiedGermlineRuntime, runtime]);
  const [busy, setBusy] = useState("");
  const [cancelling, setCancelling] = useState(false);
  const [postLockState, setPostLockState] = useState<"unsupported" | "waiting" | "held">("unsupported");
  const [progress, setProgress] = useState<{ processed: number; total: number; unit?: string }>({ processed: 0, total: store.count });
  const [error, setError] = useState("");
  const [activeWorkspace,setActiveWorkspace]=useState<PostWorkspaceId>(directLineage ? "workbench" : "alleles");
  const openModules=useMemo(()=>new Set<PostModuleId>(activeWorkspace === "overview" ? [] : [activeWorkspace]),[activeWorkspace]);
  const [skippedModules,setSkippedModules]=useState<Set<PostModuleId>>(
    ()=>new Set(initialSession?.skippedModules ?? (autoPipeline?.selection.enabled?[]:["selection"])),
  );

  const [dedupKey, setDedupKey] = useState<DedupKey>("trimmed");
  const [collapseMode, setCollapseMode] = useState<CollapseMode>("exact");
  const [collapseScope, setCollapseScope] = useState<DatasetScope>(defaultCollapseScope);
  const [respectConstantCall, setRespectConstantCall] = useState(true);
  const [dedup, setDedup] = useState<DedupDashboard | null>(null);
  const [denoiseErrorRate, setDenoiseErrorRate] = useState(0.00473);
  const [denoiseAlpha, setDenoiseAlpha] = useState(0.01);
  const [denoiseResolution, setDenoiseResolution] = useState<CallResolution>("allele");
  const [denoiseAmbiguity, setDenoiseAmbiguity] = useState<"top" | "strict">("strict");
  const [minimumParentCount, setMinimumParentCount] = useState(2);
  const [denoiseAmbiguousPolicy, setDenoiseAmbiguousPolicy] = useState<"exclude" | "retain">("exclude");
  const [denoiseUnresolvedPolicy, setDenoiseUnresolvedPolicy] = useState<"discard" | "retain">("discard");
  const [fadNeighborThreshold, setFadNeighborThreshold] = useState(1);
  const [fadMethod, setFadMethod] = useState<1 | 2>(2);
  const [expectedZeroErrorFraction, setExpectedZeroErrorFraction] = useState(1);
  const [maximumDenoiseDistance, setMaximumDenoiseDistance] = useState(1);
  const [maximumEditDistance, setMaximumEditDistance] = useState(2);
  const [minimumIndelParentRatio, setMinimumIndelParentRatio] = useState(2);
  const [denoiseCandidateCap, setDenoiseCandidateCap] = useState(50_000);
  const [workingMask, setWorkingMask] = useState<Uint8Array | null>(null);
  const [workingStages, setWorkingStages] = useState<WorkingSetStage[]>([]);
  const [selectionDraft, setSelectionDraft] = useState<RepertoireSelectionOptions>({ ...DEFAULT_REPERTOIRE_SELECTION });
  const [selectionPreview, setSelectionPreview] = useState<RepertoireSelectionResult | null>(null);
  const [selectionBaseMask, setSelectionBaseMask] = useState<Uint8Array | null>(null);
  const [selectionApplied, setSelectionApplied] = useState(false);
  const [alleleOptions, setAlleleOptions] = useState<AlleleRefinementOptions>({ ...DEFAULT_ALLELE_REFINEMENT_OPTIONS, segments: [...DEFAULT_ALLELE_REFINEMENT_OPTIONS.segments] });
  const [alleleRefinement, setAlleleRefinement] = useState<AlleleRefinementResult | null>(null);
  const [alleleApplied, setAlleleApplied] = useState(false);
  const [alleleReassignmentPolicy, setAlleleReassignmentPolicy] = useState<AlleleReassignmentPolicy>("confidence");
  const [alleleApplyMinimumPosterior, setAlleleApplyMinimumPosterior] = useState(0.8);
  const [alleleProgress, setAlleleProgress] = useState<{ processed: number; total: number; phase: string } | null>(null);
  const selectionFacets=useMemo(()=>({
    datasets:uniqueFacetItems(datasets.map((dataset)=>dataset.datasetId)),
    samples:uniqueFacetItems(datasets.map((dataset)=>dataset.sampleId)),
    subjects:uniqueFacetItems(datasets.map((dataset)=>dataset.subjectId)),
    cohorts:uniqueFacetItems(datasets.map((dataset)=>dataset.cohort)),
    timepoints:uniqueFacetItems(datasets.map((dataset)=>dataset.timepoint)),
    compartments:uniqueFacetItems(datasets.map((dataset)=>dataset.compartment??"")),
    loci:loci.map((item)=>({value:item.value,count:item.count})),
    vCalls:callFacetItems(resultFacets.vCalls),
    dCalls:callFacetItems(resultFacets.dCalls),
    jCalls:callFacetItems(resultFacets.jCalls),
    cCalls:callFacetItems(resultFacets.cCalls),
    isotypes:resultFacets.isotypes,
  }),[datasets,loci,resultFacets.cCalls,resultFacets.dCalls,resultFacets.isotypes,resultFacets.jCalls,resultFacets.vCalls]);
  const [exportFormat, setExportFormat] = useState<TableExportFormat>("tsv");
  const [exportCompression, setExportCompression] = useState<OutputCompression>("none");
  const [alignmentExportFormat, setAlignmentExportFormat] = useState<AlignmentExportFormat>("fasta");

  const [shmMetric, setShmMetric] = useState<ShmMetricKey>("vNtRate");
  const [shmStratum, setShmStratum] = useState<"all" | "locus" | "v_call" | "isotype" | "sample_id" | "subject_id" | "swig_cohort" | "swig_timepoint" | "swig_compartment">("all");
  const [shmSampleCap, setShmSampleCap] = useState(2000);
  const [shmDashboard, setShmDashboard] = useState<ShmDashboard | null>(null);
  const [shmSampleOrder, setShmSampleOrder] = useState<string[]>(()=>[...new Set(datasets.map((dataset)=>dataset.sampleId).filter(Boolean))]);
  const [missingAlleleOptions, setMissingAlleleOptions] = useState<MissingAlleleOptions>({ ...DEFAULT_MISSING_ALLELE_OPTIONS });
  const [missingAlleles, setMissingAlleles] = useState<MissingAlleleDashboard | null>(null);
  const [selectedMissingAlleleIds, setSelectedMissingAlleleIds] = useState<Set<string>>(new Set());
  const [unifiedOptions,setUnifiedOptions]=useState<UnifiedOptions>({...DEFAULT_UNIFIED_OPTIONS});
  const [unifiedGermline,setUnifiedGermline]=useState<UnifiedDashboard|null>(null);
  const [personalizedGermlineOptions, setPersonalizedGermlineOptions] = useState<PersonalizedGermlineOptions>({ ...DEFAULT_PERSONALIZED_GERMLINE_OPTIONS });
  const [personalizedGermline, setPersonalizedGermline] = useState<PersonalizedGermlineDashboard | null>(null);

  const [identity, setIdentity] = useState(0.85);
  const [lineageScope, setLineageScope] = useState<DatasetScope>(defaultLineageScope);
  const [resolution, setResolution] = useState<CallResolution>("gene");
  const [ambiguity, setAmbiguity] = useState<AmbiguityPolicy>("overlap");
  const [productiveOnly, setProductiveOnly] = useState(true);
  const [candidateCap, setCandidateCap] = useState(50_000);
  const [lineages, setLineages] = useState<LineageDashboard | null>(null);
  const [lineageSort,setLineageSort]=useState<{key:LineageSortKey;direction:"asc"|"desc"}>({key:"abundance",direction:"desc"});
  const [lineageSearch,setLineageSearch]=useState("");
  const [lineageLocusFilter,setLineageLocusFilter]=useState("");
  const [lineageSpecificSample,setLineageSpecificSample]=useState("");
  const [lineageSampleFilterMode,setLineageSampleFilterMode]=useState<LineageSampleFilterMode>("any");
  const [lineageDoubleDFilter,setLineageDoubleDFilter]=useState<LineageDoubleDFilter>("any");
  const [lineageSelectedSamples,setLineageSelectedSamples]=useState<Set<string>>(new Set());
  const [lineageSampleSearch,setLineageSampleSearch]=useState("");
  const [lineageMinAbundance,setLineageMinAbundance]=useState(0);
  const [lineageMinUnique,setLineageMinUnique]=useState(0);
  const [lineageMinSamples,setLineageMinSamples]=useState(0);
  const [lineageMinCdr3Length,setLineageMinCdr3Length]=useState(0);
  const [lineageMaxCdr3Length,setLineageMaxCdr3Length]=useState(0);
  const [lineageShmFilterStatistic,setLineageShmFilterStatistic]=useState<LineageShmFilterStatistic>("mean");
  const [lineageMinimumShm,setLineageMinimumShm]=useState(0);
  const [lineageTableSelectedIds,setLineageTableSelectedIds]=useState<Set<number>>(new Set());
  const [geneMetric, setGeneMetric] = useState<"abundance" | "lineages">("abundance");
  const [chartColor, setChartColor] = useState("#08796f");
  const [topGenes, setTopGenes] = useState(15);

  const [selectedLineage, setSelectedLineage] = useState<LineageSummary | null>(null);
  const [selectedLineageIds, setSelectedLineageIds] = useState<number[]>([]);
  const [lineageRows, setLineageRows] = useState<AirrDetailRow[]>([]);
  const [alignmentProductiveOnly, setAlignmentProductiveOnly] = useState(Boolean(initialSession?.alignmentProductiveOnly));
  const [lineageMultiplicity, setLineageMultiplicity] = useState<Map<number, number>>(new Map());
  const [originalLineageByOrdinal, setOriginalLineageByOrdinal] = useState<Map<number,number>>(new Map());
  const [lineageTotal, setLineageTotal] = useState(0);
  const workbenchRef = useRef<HTMLElement>(null);
  const [alignment, setAlignment] = useState("");
  const [alignmentMode, setAlignmentMode] = useState<"nt" | "aa">("nt");
  const [alignmentFrameOffset, setAlignmentFrameOffset] = useState<AlignmentFrameOffset>(0);
  const [alignmentMethod, setAlignmentMethod] = useState<"quick" | "kalign" | "codon">("quick");
  const [lineageGermlineMethod, setLineageGermlineMethod] = useState<LineageGermlineMethod>(() => initialSession?.lineageGermlineMethod ?? "closest");
  const [lineageTemplateSelections, setLineageTemplateSelections] = useState<Record<string, number>>(() => initialSession?.lineageTemplateSelections ?? {});
  const [alignmentLimit, setAlignmentLimit] = useState(200);
  const [alignmentSource, setAlignmentSource] = useState("");
  const [alignmentEdited, setAlignmentEdited] = useState(false);
  const [editedAlignments, setEditedAlignments] = useState<Map<string,EditedAlignmentState>>(new Map());
  const [lineageMerges,setLineageMerges]=useState<LineageMergeState[]>([]);
  const alignmentRevisionRef = useRef(0);
  const [treeRun, setTreeRun] = useState<TreeSnapshot | null>(null);
  const [treeViewMode, setTreeViewMode] = useState<TreeViewMode>("rooted");
  const [treeError, setTreeError] = useState("");
  const [treeModel, setTreeModel] = useState<"gtr" | "jc">("gtr");
  const [treeFast, setTreeFast] = useState(false);
  const [phyloUcaState, setPhyloUcaState] = useState<PhyloUcaPanelState | null>(() => initialSession?.phyloUca ?? null);
  const treeResultRef = useRef<HTMLDivElement>(null);
  const alivibeSessionRef = useRef<AlivibeRoundTripSession | null>(null);
  const alignmentRef = useRef(alignment);
  const [alignmentEditorStatus, setAlignmentEditorStatus] = useState("");
  const [alignmentEditorError, setAlignmentEditorError] = useState("");
  const alignmentInfo = useMemo(() => {
    if (!alignment) return null;
    try { return inspectAlignment(alignment); } catch { return null; }
  }, [alignment]);
  const productiveLineageRows = useMemo(()=>lineageRows.filter(isProductiveLineageRow),[lineageRows]);
  const workbenchLineageRows = alignmentProductiveOnly ? productiveLineageRows : lineageRows;
  const excludedNonProductiveRows = lineageRows.length - productiveLineageRows.length;
  const selectedGroupKey = useMemo(()=>lineageAlignmentKey(selectedLineageIds,alignmentProductiveOnly),[alignmentProductiveOnly,selectedLineageIds]);
  const requestedTemplateOrdinal = lineageTemplateSelections[selectedGroupKey];
  const templateOrdinal = workbenchLineageRows.some((row) => row.record.ordinal === requestedTemplateOrdinal) ? requestedTemplateOrdinal : undefined;
  const templateChoices = useMemo(() => lineageTemplateChoices(workbenchLineageRows, { references }), [workbenchLineageRows, references]);
  const selectedGroupKeyRef = useRef(selectedGroupKey);
  useEffect(() => { alignmentRef.current = alignment; }, [alignment]);
  useEffect(() => { selectedGroupKeyRef.current = selectedGroupKey; }, [selectedGroupKey]);
  const savedEditedAlignment = editedAlignments.get(selectedGroupKey);
  const mergedIdByOriginal = useMemo(()=>{
    const map=new Map<number,string>();
    lineageMerges.forEach((merge)=>merge.originalLineageIds.forEach((lineageId)=>map.set(lineageId,merge.id)));
    return map;
  },[lineageMerges]);

  function clearAlignmentArtifacts() {
    alignmentRevisionRef.current += 1;
    setAlignment("");
    setAlignmentSource("");
    setAlignmentEdited(false);
    setAlignmentFrameOffset(0);
    setTreeRun(null);
    setTreeError("");
    setPhyloUcaState(null);
  }

  function installAlignment(next: string, source: string, edited = false, lineageIds: number[] = selectedLineageIds, frameOffset?: AlignmentFrameOffset, productiveOnly = alignmentProductiveOnly) {
    const inspected = inspectAlignment(next);
    const resolvedFrameOffset = frameOffset ?? inferAlignedReadingFrame(inspected.records.map((record) => record.sequence)).offset;
    alignmentRevisionRef.current += 1;
    setAlignment(inspected.fasta);
    setAlignmentSource(source);
    setAlignmentEdited(edited);
    setAlignmentFrameOffset(resolvedFrameOffset);
    setTreeRun(null);
    setTreeError("");
    setPhyloUcaState(null);
    const key=lineageAlignmentKey(lineageIds,productiveOnly);
    if(edited&&key){
      const entry:EditedAlignmentState={key,lineageIds:[...new Set(lineageIds)].sort((a,b)=>a-b),productiveOnly,fasta:inspected.fasta,source,frameOffset:resolvedFrameOffset,savedAt:new Date().toISOString()};
      setEditedAlignments((current)=>{const updated=new Map(current);updated.set(key,entry);return updated;});
    }
    return inspected;
  }

  function changeAlignmentFrameOffset(frameOffset: AlignmentFrameOffset) {
    setAlignmentFrameOffset(frameOffset);
    setTreeRun((current) => current ? { ...current, frameOffset } : current);
    if (!alignmentEdited || !selectedGroupKey) return;
    setEditedAlignments((current) => {
      const entry = current.get(selectedGroupKey);
      if (!entry) return current;
      const updated = new Map(current);
      updated.set(selectedGroupKey, { ...entry, frameOffset });
      return updated;
    });
  }

  function downloadCurrentAlignment() {
    if(!alignment||!selectedLineage)return;
    const records=parseFasta(alignment,true).map((record)=>({name:record.name,sequence:record.sequence}));
    const text=alignmentText(records,alignmentExportFormat);
    downloadText(text,`${baseName(inputName)}.lineage-${selectedLineage.id}.alignment${alignmentExtension(alignmentExportFormat)}`);
  }

  const isTcr = scope === "TCR" || String(scope).startsWith("TR");
  const [chmmSegment, setChmmSegment] = useState<ChmmSegment>("V");
  const [chmmMethod, setChmmMethod] = useState<"BW" | "DB">(isTcr ? "DB" : "BW");
  const [chmmSource, setChmmSource] = useState<"selected" | "upload">("selected");
  const [uploadedMsa, setUploadedMsa] = useState("");
  const [uploadedMsaName, setUploadedMsaName] = useState("");
  const [preparedMsa, setPreparedMsa] = useState("");
  const [chmmThreshold, setChmmThreshold] = useState(0.95);
  const [chmmPrior, setChmmPrior] = useState(0.05);
  const [chmmMinDfr, setChmmMinDfr] = useState(1);
  const [chmmDetailed, setChmmDetailed] = useState(false);
  const [mutationRates, setMutationRates] = useState(isTcr ? "0.005" : "0,0.0179,0.0357,0.0536,0.0714,0.0893,0.1071,0.125,0.1429,0.1607,0.1786,0.1964,0.2143,0.2321,0.25");
  const [chmm, setChmm] = useState<ChmmDashboard | null>(null);
  const [chmmTopIndex, setChmmTopIndex] = useState<Map<number, AirrIndexRecord>>(new Map());
  const [chmmRun, setChmmRun] = useState<{ msa: string; options: ChmmRunOptions; inputMask: Uint8Array | null } | null>(null);
  const [chmmFilterThreshold, setChmmFilterThreshold] = useState(0.95);
  useEffect(() => {
    let cancelled = false;
    if (!chmm?.top.length) {
      setChmmTopIndex(new Map());
      return () => { cancelled = true; };
    }
    void store.indexRecords(chmm.top.slice(0, 12).map((record) => record.ordinal)).then((records) => {
      if (!cancelled) setChmmTopIndex(new Map(records.map((record) => [record.ordinal, record])));
    });
    return () => { cancelled = true; };
  }, [chmm, store]);
  const [retainUnevaluated, setRetainUnevaluated] = useState(true);
  const [chimeraDetail, setChimeraDetail] = useState<ChmmDetail | null>(null);
  const chimeraDetailRef = useRef<HTMLDivElement>(null);

  const [queryText, setQueryText] = useState("");
  const [queryTarget, setQueryTarget] = useState<QueryTarget>("cdr3_nt");
  const [queryMetric, setQueryMetric] = useState<QueryMetric>("hamming");
  const [queryIdentity, setQueryIdentity] = useState(0.85);
  const [queryLimit, setQueryLimit] = useState(250);
  const [queryLocus, setQueryLocus] = useState("");
  const [queryV, setQueryV] = useState("");
  const [queryJ, setQueryJ] = useState("");
  const [queryConstraintMode, setQueryConstraintMode] = useState<"manual" | "infer">("manual");
  const [queryResultMode,setQueryResultMode]=useState<"sequences"|"lineages">("sequences");
  const [queryInference, setQueryInference] = useState<InferredQueryAssignment[]>([]);
  const [queryHits, setQueryHits] = useState<QueryHit[]>([]);
  const [queryRecords, setQueryRecords] = useState<Map<number, AirrIndexRecord>>(new Map());
  const [queryLineageHits,setQueryLineageHits]=useState<QueryLineageHit[]>([]);
  const [expanded, setExpanded] = useState<{ ordinals: number[]; comparisons: number; capped: boolean } | null>(null);
  const restoredSessionRef = useRef(false);
  const directLineageInitializedRef = useRef(false);
  const pipelineRunRef = useRef(false);
  const sessionChangeReadyRef = useRef(false);
  const sessionChangeCallbackRef = useRef(onSessionChange);
  const [pipelineReport, setPipelineReport] = useState<string[]>([]);
  const [neighbourCdr3Identity,setNeighbourCdr3Identity]=useState(()=>Math.max(0.7,identity-0.05));
  const [neighbourGermlineIdentity,setNeighbourGermlineIdentity]=useState(0.97);
  const [neighbourMethod,setNeighbourMethod]=useState<"cdr3"|"germline"|"either">("either");
  const [neighbourCdr3Metric,setNeighbourCdr3Metric]=useState<"hamming"|"edit">("hamming");
  const [neighbourCallPolicy,setNeighbourCallPolicy]=useState<"both"|"either"|"ignore">("both");
  const [neighbourScope,setNeighbourScope]=useState<DatasetScope>(defaultLineageScope);
  const [neighbourLimit,setNeighbourLimit]=useState(50);
  const [neighbourResult,setNeighbourResult]=useState<LineageNeighbourResult|null>(null);
  const [germlineNeighbourScores,setGermlineNeighbourScores]=useState<GermlineNeighbourScore[]>([]);
  const [germlineSketchIndex,setGermlineSketchIndex]=useState<LineageGermlineSketchIndex|null>(null);
  const [selectedNeighbourIds,setSelectedNeighbourIds]=useState<Set<number>>(new Set());
  const combinedNeighbourHits = useMemo(() => {
    const hits = new Map<number, CombinedNeighbourHit>();
    for (const cdr3 of neighbourResult?.hits ?? []) hits.set(cdr3.lineageId, { lineageId: cdr3.lineageId, cdr3 });
    for (const germline of germlineNeighbourScores) {
      const current = hits.get(germline.lineageId) ?? { lineageId: germline.lineageId };
      current.germline = germline;
      hits.set(germline.lineageId, current);
    }
    return [...hits.values()].sort((left, right) =>
      Math.max(right.cdr3?.cdr3Identity ?? 0, right.germline?.germlineIdentity ?? 0) - Math.max(left.cdr3?.cdr3Identity ?? 0, left.germline?.germlineIdentity ?? 0) ||
      left.lineageId - right.lineageId,
    ).slice(0, neighbourLimit);
  }, [neighbourResult, germlineNeighbourScores, neighbourLimit]);

  useEffect(()=>{
    const currentSamples=[...new Set(datasets.map((dataset)=>dataset.sampleId).filter(Boolean))];
    setShmSampleOrder((current)=>{
      const retained=current.filter((sample)=>currentSamples.includes(sample));
      const added=currentSamples.filter((sample)=>!retained.includes(sample));
      const next=[...retained,...added];
      return next.length===current.length&&next.every((sample,index)=>sample===current[index])?current:next;
    });
    setLineageSelectedSamples((current)=>new Set([...current].filter((sample)=>currentSamples.includes(sample))));
  },[datasets]);

  const nextUnskippedModule=(from:PostModuleId, skipped=skippedModules):PostWorkspaceId=>{
    const start=POST_MODULE_ORDER.indexOf(from);
    for(let index=start+1;index<POST_MODULE_ORDER.length;index+=1){
      const candidate=POST_MODULE_ORDER[index];
      if(!skipped.has(candidate)&&(candidate!=="workbench"||selectedLineage))return candidate;
    }
    return "overview";
  };
  const toggleModule=(module:PostModuleId)=>setActiveWorkspace(module);
  const advanceModule=(_from:PostModuleId,to:PostModuleId)=>setActiveWorkspace(skippedModules.has(to)?nextUnskippedModule(to):to);
  const openModule=(module:PostModuleId)=>setActiveWorkspace(module);
  const skipModule=(module:PostModuleId)=>{
    const next=new Set(skippedModules);next.add(module);setSkippedModules(next);setActiveWorkspace(nextUnskippedModule(module,next));
  };
  const includeModule=(module:PostModuleId)=>setSkippedModules((current)=>{const next=new Set(current);next.delete(module);return next;});
  const skipStepButton=(module:PostModuleId)=>skippedModules.has(module)
    ? <button className="module-skip-step skipped" type="button" onClick={()=>includeModule(module)}>Include step</button>
    : <button className="module-skip-step" type="button" onClick={()=>skipModule(module)}>Skip step</button>;
  const moduleClass=(module:PostModuleId,base:string)=>`${base}${openModules.has(module)?" is-open":" is-collapsed"}${skippedModules.has(module)?" is-skipped":""}`;
  const lineageGermline = useMemo(() => {
    if (!workbenchLineageRows.length) return null;
    try { return inferLineageGermline(workbenchLineageRows, lineageGermlineMethod, { references, templateOrdinal }); } catch { return null; }
  }, [workbenchLineageRows, lineageGermlineMethod, references, templateOrdinal]);

  useEffect(() => {
    if (!alignment || !workbenchLineageRows.length || lineageGermlineMethod !== "closest") return;
    try {
      const repaired = restoreLineageOuterFlanks(alignment, workbenchLineageRows, { references, templateOrdinal });
      if (!repaired.changed) return;
      installAlignment(repaired.fasta, `${alignmentSource} · restored outer germline flanks`, alignmentEdited, selectedLineageIds, ((alignmentFrameOffset + repaired.addedLeftColumns) % 3) as AlignmentFrameOffset);
      setAlignmentEditorStatus("Restored missing V 5′ / J 3′ germline flanks in the guide. Curated interior columns and read nucleotides were retained; rerun tree/UCA inference.");
    } catch {
      // Curation may have removed the anchor's query bases. Do not guess a map.
      setAlignmentEditorError("The current MSA cannot safely project this template's AIRR coordinates. Rebuild the MSA or choose an intact member to restore missing germline flanks.");
    }
  }, [alignment, workbenchLineageRows, references, templateOrdinal, lineageGermlineMethod]);

  useEffect(() => {
    if (directLineage || !autoPipeline?.enabled || initialSession || pipelineRunRef.current) return;
    pipelineRunRef.current = true;
    void (async () => {
      const report: string[] = [];
      pipelineActiveRef.current = true;
      runtime.beginTransaction();
      setError("");
      setPipelineReport([]);
      try {
        await runInActiveLock(async () => {
          setBusy("Pipeline · indexing AIRR records");
          await runtime.ensureIndexed((processed, total) => setProgress({ processed, total }), postLockAbortRef.current?.signal);
          let activeMask: Uint8Array | null = null;
          let collapseResult: DedupDashboard | null = null;
          const stages: WorkingSetStage[] = [];

          // Assignment pooling must precede every stage whose partitions can
          // depend on V/D/J calls (notably rearrangement-key collapse).
          let pipelineAlleleResult: AlleleRefinementResult | null = null;
          const pipelineAllelePolicy = autoPipeline.alleleRefinement.reassignmentPolicy ?? "confidence";
          const pipelineAlleleThreshold = autoPipeline.alleleRefinement.applyMinimumPosterior;
          const overlayPipelineAlleleCalls = (row: { ordinal: number; values: Record<string,string> }) => {
            if (!pipelineAlleleResult) return;
            const v = refinedCall(pipelineAlleleResult, "V", row.ordinal, pipelineAllelePolicy, pipelineAlleleThreshold);
            const d = refinedCall(pipelineAlleleResult, "D", row.ordinal, pipelineAllelePolicy, pipelineAlleleThreshold);
            const j = refinedCall(pipelineAlleleResult, "J", row.ordinal, pipelineAllelePolicy, pipelineAlleleThreshold);
            if (v) row.values.v_call = v;
            if (d) row.values.d_call = d;
            if (j) row.values.j_call = j;
          };
          if (autoPipeline.alleleRefinement.enabled) {
            const options: AlleleRefinementOptions = {
              ...DEFAULT_ALLELE_REFINEMENT_OPTIONS,
              model: autoPipeline.alleleRefinement.model ?? "dirichlet",
              scope: autoPipeline.alleleRefinement.scope,
              segments: [...autoPipeline.alleleRefinement.segments],
              weighting: autoPipeline.alleleRefinement.weighting,
              baselineNeighbourOdds: autoPipeline.alleleRefinement.baselineNeighbourOdds,
              shmLeakageSensitivity: autoPipeline.alleleRefinement.shmLeakageSensitivity,
            };
            setAlleleOptions(options);
            setAlleleReassignmentPolicy(pipelineAllelePolicy);
            setAlleleApplyMinimumPosterior(pipelineAlleleThreshold);
            setBusy("Pipeline · pooling germline evidence before downstream partitioning");
            pipelineAlleleResult = await alleleRuntime.run(store, references, options, null, (next) => {
              setAlleleProgress(next);
              setBusy(`Pipeline · ${next.phase}`);
              setProgress({ processed: next.processed, total: next.total, unit: "records or independent donor pools" });
            }, postLockAbortRef.current?.signal);
            setAlleleProgress(null);
            setAlleleRefinement(pipelineAlleleResult);
            setAlleleApplied(true);
            await runtime.setRepertoireCallOverrides(pipelineAlleleResult, pipelineAllelePolicy, pipelineAlleleThreshold);
            report.push(`Allele pooling first fitted ${Object.values(pipelineAlleleResult.segments).reduce((sum, segment) => sum + (segment?.models.length ?? 0), 0).toLocaleString()} independent ${options.scope === "subject" ? "donor" : options.scope} / locus / segment ${options.model === "active-set" ? "hurdle active-set" : "Dirichlet"} models and applied ${pipelineAllelePolicy === "best" ? "the posterior MAP call to every modeled record" : `posterior MAP calls at confidence ≥ ${pipelineAlleleThreshold}`}.`);
            if (!autoPipeline.lineage.enabled) openModule("alleles");
          }

          if (autoPipeline.collapse.enabled) {
            setCollapseMode(autoPipeline.collapse.mode);
            setDedupKey(autoPipeline.collapse.key);
            setCollapseScope(autoPipeline.collapse.scope);
            setDenoiseUnresolvedPolicy(autoPipeline.collapse.unresolvedPolicy);
            setRespectConstantCall(autoPipeline.collapse.respectConstantCall ?? true);
            setBusy(`Pipeline · ${autoPipeline.collapse.mode === "exact" ? "exact collapse" : `${autoPipeline.collapse.mode} denoising`}`);
            collapseResult = autoPipeline.collapse.mode === "exact"
              ? await runtime.deduplicate(autoPipeline.collapse.key, autoPipeline.collapse.unresolvedPolicy, autoPipeline.collapse.scope, autoPipeline.collapse.respectConstantCall ?? true)
              : await runtime.denoise({
                mode: autoPipeline.collapse.mode,
                errorRate: 0.00473,
                alpha: 0.01,
                callResolution: "allele",
                ambiguity: "strict",
                minimumParentCount: 2,
                ambiguousPolicy: autoPipeline.collapse.mode === "fad" ? "exclude" : "retain",
                unresolvedPolicy: autoPipeline.collapse.unresolvedPolicy,
                fadNeighborThreshold: 1,
                fadMethod: 2,
                expectedZeroErrorFraction: 1,
                maximumHammingDistance: 1,
                maximumEditDistance: 2,
                minimumIndelParentRatio: 2,
                maxCandidatesPerVariant: 50_000,
                scope: autoPipeline.collapse.scope,
                respectConstantCall: autoPipeline.collapse.respectConstantCall ?? true,
              }, (processed, total, phase) => {
                const stage = phase === "variants" ? "indexed variant denoising" : phase === "finalize" ? "representative materialization" : "streaming VDJ sequences";
                setBusy(`Pipeline · ${autoPipeline.collapse.mode} denoising · ${stage}`);
                setProgress({ processed, total, unit: phase === "variants" ? "unique sequence variants processed" : phase === "finalize" ? "representative-state operations completed" : "AIRR records streamed into the denoising index" });
              }, postLockAbortRef.current?.signal);
            setDedup(collapseResult);
            const applied = await runtime.applyDedupFilter();
            activeMask = applied.mask;
            stages.push({
              id: "dedup",
              label: `${collapseResult.mode === "exact" ? "Exact collapse" : `${collapseResult.mode} denoising`} · ${DATASET_SCOPE_LABELS[autoPipeline.collapse.scope]}`,
              input: collapseResult.inputRecords,
              retained: applied.retained,
              discarded: collapseResult.collapsedRecords,
              detail: `${collapseResult.algorithm}. Multiplicity retained in duplicate_count.`,
            });
            report.push(`Collapse retained ${applied.retained.toLocaleString()} representatives.`);
          }

          if (autoPipeline.chimera.enabled) {
            const segment = autoPipeline.chimera.segment;
            setChmmSegment(segment);
            setBusy(autoPipeline.chimera.msaSource === "upload" ? `Pipeline · validating loaded ${segment} reference MSA` : `Pipeline · building ${segment} reference MSA`);
            const referenceSource = autoPipeline.chimera.msaSource === "upload" ? autoPipeline.chimera.uploadedMsa : references[segment];
            const msa = autoPipeline.chimera.msaSource === "upload" ? referenceSource : await runKalignTask(referenceSource, postLockAbortRef.current?.signal);
            prepareReferenceMsa(msa);
            setPreparedMsa(msa);
            const tcr = scope === "TCR" || String(scope).startsWith("TR");
            const method = autoPipeline.chimera.model === "auto" ? (tcr ? "DB" : "BW") : autoPipeline.chimera.model;
            const options: ChmmRunOptions = {
              segment,
              method,
              priorProbability: 0.05,
              baseMutationProbability: 0.05,
              mutationRates: method === "DB" ? (tcr ? [0.005] : [0,0.0179,0.0357,0.0536,0.0714,0.0893,0.1071,0.125,0.1429,0.1607,0.1786,0.1964,0.2143,0.2321,0.25]) : [],
              mutationSwitchProbability: 0,
              detailed: false,
              minDfr: 1,
              threshold: autoPipeline.chimera.posteriorThreshold,
              workers,
            };
            setChmmMethod(method);
            setChmmThreshold(autoPipeline.chimera.posteriorThreshold);
            setChmmFilterThreshold(autoPipeline.chimera.posteriorThreshold);
            setRetainUnevaluated(autoPipeline.chimera.retainUnevaluated);
            setBusy(`Pipeline · CHMMAIRRa ${segment}`);
            const dashboard = await runChmmairra(store, msa, options, activeMask ?? undefined, (processed, total) => setProgress({ processed, total }), postLockAbortRef.current?.signal);
            setChmm(dashboard);
            setChmmRun({ msa, options, inputMask: activeMask });
            const next = new Uint8Array(store.count);
            let retained = 0;
            let unevaluatedRetained = 0;
            for (let ordinal = 0; ordinal < store.count; ordinal += 1) {
              if (activeMask && !activeMask[ordinal]) continue;
              const probability = dashboard.probabilities[ordinal];
              const keep = Number.isFinite(probability) ? probability < autoPipeline.chimera.posteriorThreshold : autoPipeline.chimera.retainUnevaluated;
              if (!keep) continue;
              next[ordinal] = 1;
              retained += 1;
              if (!Number.isFinite(probability)) unevaluatedRetained += 1;
            }
            await runtime.setActiveMask(next);
            activeMask = next;
            stages.push({ id: "chimera", label: `${segment} chimera posterior < ${autoPipeline.chimera.posteriorThreshold}`, input: dashboard.inputRecords, retained, discarded: dashboard.inputRecords - retained, detail: autoPipeline.chimera.retainUnevaluated ? `${unevaluatedRetained.toLocaleString()} unevaluated records retained.` : "Unevaluated records excluded." });
            report.push(`CHMMAIRRa retained ${retained.toLocaleString()} records.`);
          }

          if (autoPipeline.selection.enabled) {
            setBusy("Pipeline · applying repertoire selection");
            const selectionOptions: RepertoireSelectionOptions = {
              ...DEFAULT_REPERTOIRE_SELECTION,
              datasetId: autoPipeline.selection.datasetId,
              sampleId: autoPipeline.selection.sampleId,
              subjectId: autoPipeline.selection.subjectId,
              cohort: autoPipeline.selection.cohort,
              timepoint: autoPipeline.selection.timepoint,
              compartment: autoPipeline.selection.compartment,
              locus: autoPipeline.selection.locus,
              vCall: autoPipeline.selection.vCall,
              vCallIncludeAmbiguous: autoPipeline.selection.vCallIncludeAmbiguous,
              jCall: autoPipeline.selection.jCall,
              jCallIncludeAmbiguous: autoPipeline.selection.jCallIncludeAmbiguous,
              cdr3Nt: autoPipeline.selection.cdr3Nt,
              cdr3Aa: autoPipeline.selection.cdr3Aa,
              productive: autoPipeline.selection.productive,
              hasCdr3: autoPipeline.selection.hasCdr3,
              doubleD: autoPipeline.selection.doubleD,
            };
            const baseMask = activeMask?.slice() ?? null;
            const selected = await selectRepertoire(store, selectionOptions, activeMask ?? undefined, (processed) => setProgress({ processed, total: activeMask ? activeMask.reduce((sum, value) => sum + value, 0) : store.count }), postLockAbortRef.current?.signal);
            await runtime.setActiveMask(selected.mask);
            activeMask = selected.mask;
            setSelectionDraft(selectionOptions);
            setSelectionBaseMask(baseMask);
            setSelectionPreview(selected);
            setSelectionApplied(true);
            stages.push({ id: "selection", label: `Repertoire selection · ${selected.summary}`, input: selected.inputRecords, retained: selected.retainedRecords, discarded: selected.discardedRecords, detail: "Configured before annotation and committed automatically by pipeline mode." });
            report.push(`Repertoire selection retained ${selected.retainedRecords.toLocaleString()} records.`);
          }

          let lineageResult: LineageDashboard | null = null;
          if (autoPipeline.lineage.enabled) {
            setIdentity(autoPipeline.lineage.identity);
            setResolution(autoPipeline.lineage.resolution);
            setAmbiguity(autoPipeline.lineage.ambiguity);
            setProductiveOnly(autoPipeline.lineage.productiveOnly);
            setLineageScope(autoPipeline.lineage.scope);
            setBusy("Pipeline · assigning lineages");
            lineageResult = await runtime.assignLineages({
              identity: autoPipeline.lineage.identity,
              callResolution: autoPipeline.lineage.resolution,
              ambiguity: autoPipeline.lineage.ambiguity,
              productiveOnly: autoPipeline.lineage.productiveOnly,
              requireSameLocus: true,
              maxCandidateComparisons: 50_000,
              scope: autoPipeline.lineage.scope,
            }, Boolean(collapseResult));
            setLineages(lineageResult);
            openModule("lineage");
            report.push(`Assigned ${lineageResult.lineageCount.toLocaleString()} lineages within ${DATASET_SCOPE_LABELS[autoPipeline.lineage.scope].toLowerCase()}.`);
          }

          if (autoPipeline.shm.enabled) {
            setShmMetric(autoPipeline.shm.metric);
            setBusy("Pipeline · calculating SHM distributions");
            const assignments = lineageResult ? await runtime.lineageAssignments() : null;
            const counts = collapseResult ? await runtime.dedupCounts() : null;
            const accumulator = new ShmAccumulator({ metric: autoPipeline.shm.metric, maxSamplesPerLineage: 2000 });
            const fields = ["sequence_id","v_call","j_call","locus","isotype","sample_id","subject_id","swig_cohort","swig_timepoint","swig_compartment","duplicate_count","v_sequence_alignment","v_germline_alignment","sequence_alignment","germline_alignment","v_sequence_start","v_germline_start","sequence_frame","v_frame","cdr1_start","cdr1_end","cdr2_start","cdr2_end","fwr1_start","fwr1_end","fwr2_start","fwr2_end","fwr3_start","fwr3_end"];
            await store.scanAirrRows(fields, async (rows) => { for (const row of rows) { overlayPipelineAlleleCalls(row); if (counts?.[row.ordinal]) row.values.duplicate_count = String(counts[row.ordinal]); accumulator.add(row.values, row.ordinal, assignments?.[row.ordinal] ?? 0, "All selected"); } }, { batchSize: 1500, includeMask: activeMask ?? undefined, onProgress: (processed, total) => setProgress({ processed, total }), signal: postLockAbortRef.current?.signal });
            const dashboard = accumulator.finish();
            setShmDashboard(dashboard);
            report.push(`SHM summarized ${dashboard.analyzedRecords.toLocaleString()} records.`);
          }

          if (autoPipeline.missingAlleles.enabled) {
            if (!lineageResult) throw new Error("Pipeline missing-allele screening requires lineage assignment to be enabled.");
            setBusy("Pipeline · discovering recurrent V haplotypes (pass 1 of 2)");
            const assignments = await runtime.lineageAssignments();
            const accumulator = new MissingAlleleAccumulator(DEFAULT_MISSING_ALLELE_OPTIONS);
            const fields = ["subject_id","v_call","j_call","cdr3","junction","v_germline_start","v_sequence_alignment","v_germline_alignment"];
            await store.scanAirrRows(fields, async (rows) => { for (const row of rows) { overlayPipelineAlleleCalls(row); accumulator.add(row.values, row.ordinal, assignments[row.ordinal] ?? 0); } }, { batchSize: 1500, includeMask: activeMask ?? undefined, onProgress: (processed, total) => setProgress({ processed, total: total * 2 }), signal: postLockAbortRef.current?.signal });
            setBusy("Pipeline · screening every lineage member (pass 2 of 2)");
            const validator = accumulator.prepareValidation(references.V);
            await store.scanAirrRows(fields, async (rows) => { for (const row of rows) { overlayPipelineAlleleCalls(row); validator.add(row.values, row.ordinal, assignments[row.ordinal] ?? 0); } }, { batchSize: 1500, includeMask: activeMask ?? undefined, onProgress: (processed, total) => setProgress({ processed: total + processed, total: total * 2 }), signal: postLockAbortRef.current?.signal });
            const dashboard = validator.finish();
            setMissingAlleles(dashboard);
            setSelectedMissingAlleleIds(new Set());
            report.push(`Missing-allele screen produced ${dashboard.candidates.length.toLocaleString()} candidate${dashboard.candidates.length === 1 ? "" : "s"}.`);
          }

          setWorkingMask(activeMask);
          setWorkingStages(stages);
        });
        runtime.commitTransaction();
        pipelineActiveRef.current = false;
        setPipelineReport(report.length ? report : ["Annotation completed; no automatic post-analysis stage was selected."]);
      } catch (pipelineError) {
        if (!cancellationRecoveryRef.current) {
          if (runtime.requiresRecoveryForCancellation()) await runtime.rollbackTransaction((processed, total) => setProgress({ processed, total, unit: "AIRR records restored" }));
          else runtime.commitTransaction();
          setError(pipelineError instanceof Error ? pipelineError.message : String(pipelineError));
          setPipelineReport([...report, "Pipeline stopped before all selected stages completed."]);
        }
      } finally {
        if (!cancellationRecoveryRef.current) setBusy("");
        pipelineActiveRef.current = false;
      }
    })();
  }, [alleleRuntime, autoPipeline, directLineage, initialSession, references, runtime, scope, store, workers]);

  useEffect(() => {
    if(!sessionHandleRef)return;
    const handle:PostAnalysisSessionHandle={snapshot:async(signal)=>{
      if(signal?.aborted)throw new DOMException("Session saving was cancelled.","AbortError");
      const activeMask=await runtime.activeMask();
      if(signal?.aborted)throw new DOMException("Session saving was cancelled.","AbortError");
      const collapse=dedup?await (async()=>{const state=await runtime.dedupState();return {mode:dedup.mode,options:{dedupKey,collapseMode,collapseScope,respectConstantCall,denoiseErrorRate,denoiseAlpha,denoiseResolution,denoiseAmbiguity,minimumParentCount,denoiseAmbiguousPolicy,denoiseUnresolvedPolicy,fadNeighborThreshold,fadMethod,expectedZeroErrorFraction,maximumDenoiseDistance,maximumEditDistance,minimumIndelParentRatio,denoiseCandidateCap},counts:packSessionVector(state.counts),representatives:packSessionVector(state.representatives),dashboard:{...dedup}};})():undefined;
      if(signal?.aborted)throw new DOMException("Session saving was cancelled.","AbortError");
      const lineage=lineages?{options:{identity,resolution,ambiguity,productiveOnly,candidateCap,lineageScope},assignments:packSessionVector(await runtime.lineageAssignments()),dashboard:{...lineages}}:undefined;
      if(signal?.aborted)throw new DOMException("Session saving was cancelled.","AbortError");
      const chimera=chmm&&chmmRun?{options:{...chmmRun.options,chmmSource,uploadedMsaName,mutationRates,retainUnevaluated},msa:chmmRun.msa,dashboard:Object.fromEntries(Object.entries(chmm).filter(([key])=>key!=="probabilities"&&key!=="dfr")),filterThreshold:chmmFilterThreshold,probabilities:packSessionVector(chmm.probabilities),dfr:packSessionVector(chmm.dfr),retainedMask:chmmRun.inputMask?packSessionVector(chmmRun.inputMask):undefined}:undefined;
      return {skippedModules:[...skippedModules],workingStages:[...workingStages],activeMask:activeMask?packSessionVector(activeMask):undefined,collapse,chimera,selection:selectionApplied||selectionPreview?{options:{...selectionDraft},mask:selectionPreview?packSessionVector(selectionPreview.mask):undefined,baseMask:selectionBaseMask?packSessionVector(selectionBaseMask):undefined}:undefined,alleleRefinement:alleleRefinement?saveAlleleRefinement(alleleRefinement,alleleApplied,alleleReassignmentPolicy,alleleApplyMinimumPosterior):undefined,lineage,selectedLineageIds:[...selectedLineageIds],lineageGermlineMethod,lineageTemplateSelections:{...lineageTemplateSelections},
        alignmentFrameOffset,alignmentProductiveOnly,
        query:{queryText,queryTarget,queryMetric,queryIdentity,queryLimit,queryLocus,queryV,queryJ,queryConstraintMode,queryResultMode,queryInference,queryHits,queryLineageHits,expanded},
        editedAlignments:[...editedAlignments.values()].map((entry)=>({...entry,lineageIds:[...entry.lineageIds]})),
        lineageMerges:lineageMerges.map((merge)=>({...merge,originalLineageIds:[...merge.originalLineageIds]})),
        tree:treeRun?{rawNewick:treeRun.newick,rootedNewick:treeRun.rootedNewick,stableNewick:treeRun.stableNewick,source:treeRun.source,lineageIds:[...selectedLineageIds],run:{...treeRun}}:undefined,phyloUca:phyloUcaState??undefined,
        shm:shmDashboard?{metric:shmMetric,dashboard:shmDashboard,sampleOrder:[...shmSampleOrder]}:undefined,missingAlleles:missingAlleles?{options:missingAlleleOptions,dashboard:missingAlleles,selectedCandidateIds:[...selectedMissingAlleleIds]}:undefined,unifiedGermline:unifiedGermline?{options:unifiedOptions,dashboard:unifiedGermline}:undefined,personalizedGermline:personalizedGermline?{options:personalizedGermlineOptions,dashboard:personalizedGermline}:undefined} satisfies PostAnalysisSessionSnapshot;
    }};
    sessionHandleRef.current=handle;
    return()=>{if(sessionHandleRef.current===handle)sessionHandleRef.current=null;};
  });

  useEffect(()=>{sessionChangeCallbackRef.current=onSessionChange;},[onSessionChange]);

  useEffect(()=>{
    if(!sessionChangeReadyRef.current){sessionChangeReadyRef.current=true;return;}
    const reason=treeRun?"phylogeny_changed":editedAlignments.size?"edited_alignment_changed":alignment?"lineage_alignment_changed":personalizedGermline?"personalized_germline_changed":missingAlleles?"missing_allele_screen_changed":shmDashboard?"shm_changed":lineages?"lineages_changed":chmm?"chimera_state_changed":dedup?"collapse_state_changed":selectionApplied||selectionPreview?"repertoire_selection_changed":"post_analysis_state_changed";
    sessionChangeCallbackRef.current?.(reason);
  },[
    alignment,alignmentFrameOffset,alignmentProductiveOnly,alleleApplied,alleleReassignmentPolicy,alleleApplyMinimumPosterior,alleleOptions,alleleRefinement,chmm,dedup,editedAlignments,expanded,lineageGermlineMethod,lineageTemplateSelections,lineageMerges,lineages,respectConstantCall,
    missingAlleleOptions,missingAlleles,personalizedGermlineOptions,personalizedGermline,unifiedOptions,unifiedGermline,queryConstraintMode,queryHits,queryIdentity,queryJ,queryLimit,
    queryLocus,queryMetric,queryResultMode,queryTarget,queryText,queryV,selectedLineageIds,selectionApplied,skippedModules,
    selectedMissingAlleleIds,selectionPreview,shmDashboard,shmMetric,shmSampleOrder,treeRun,phyloUcaState,workingStages,
  ]);

  useEffect(() => {
    if (!directLineage || directLineageInitializedRef.current) return;
    directLineageInitializedRef.current = true;
    let cancelled = false;
    void (async () => {
      setBusy("Loading the uploaded lineage");
      setError("");
      try {
        const rowLimit = Math.min(store.count, 500);
        const ordinals = Array.from({ length: rowLimit }, (_, index) => rowLimit === store.count
          ? index
          : Math.floor(index * (store.count - 1) / Math.max(1, rowLimit - 1)));
        const rows = await store.detailMany(ordinals);
        if (cancelled) return;
        if (!rows.length) throw new Error("The assigned input contains no records to open as a lineage.");

        const valueCounts = (values: string[]) => [...values.reduce((counts, value) => {
          if (value) counts.set(value, (counts.get(value) ?? 0) + 1);
          return counts;
        }, new Map<string, number>()).entries()].sort((left, right) => right[1] - left[1] || left[0].localeCompare(right[0]));
        const locus = valueCounts(rows.map((row) => row.values.locus || row.record.locus))[0]?.[0] ?? "";
        const sampleCounts = valueCounts(rows.map((row) => row.values.sample_id || row.record.sampleId)).map(([sampleId, count]) => ({ sampleId, uniqueMembers: count, abundance: count }));
        const samples = sampleCounts.map((item) => item.sampleId);
        const summary: LineageSummary = {
          id: 1,
          representativeOrdinal: rows[0].record.ordinal,
          uniqueMembers: store.count,
          abundance: store.count,
          locus,
          vCalls: [...new Set(rows.flatMap((row) => (row.values.v_call || row.record.vCall).split(",").map((value) => value.trim()).filter(Boolean)))].slice(0, 12),
          jCalls: [...new Set(rows.flatMap((row) => (row.values.j_call || row.record.jCall).split(",").map((value) => value.trim()).filter(Boolean)))].slice(0, 12),
          cdr3Length: Math.max(0, Math.round(rows.reduce((sum, row) => sum + (row.values.cdr3 || row.record.cdr3).length, 0) / rows.length)),
          studyScope: "global",
          studyGroup: "uploaded lineage",
          sampleIds: samples,
          sampleCounts,
          subjectIds: [...new Set(rows.map((row) => row.values.subject_id || row.record.subjectId).filter(Boolean))],
          timepoints: [...new Set(rows.map((row) => row.values.swig_timepoint || row.record.timepoint).filter(Boolean))],
          compartments: [...new Set(rows.map((row) => row.values.swig_compartment || row.record.compartment).filter(Boolean))],
          doubleDPositiveMembers: rows.filter((row) => Boolean(row.values.d2_call || row.record.d2Call)).length,
          doubleDPositiveAbundance: rows.filter((row) => Boolean(row.values.d2_call || row.record.d2Call)).length,
        };
        const lineageByOrdinal = new Map(rows.map((row) => [row.record.ordinal, 1] as const));
        const multiplicity = new Map(rows.map((row) => {
          const imported = Number(row.values.duplicate_count);
          return [row.record.ordinal, Number.isFinite(imported) && imported > 0 ? Math.floor(imported) : 1] as const;
        }));
        setSelectedLineage(summary);
        setSelectedLineageIds([1]);
        setLineageRows(rows);
        setLineageTotal(store.count);
        setOriginalLineageByOrdinal(lineageByOrdinal);
        setLineageMultiplicity(multiplicity);
        setActiveWorkspace("workbench");

        const restoredEdit = initialSession?.editedAlignments?.find((entry) => entry.lineageIds.includes(1) && Boolean(entry.productiveOnly) === alignmentProductiveOnly);
        const savedFrameOffset = validAlignmentFrameOffset(initialSession?.alignmentFrameOffset);
        if (restoredEdit) {
          installAlignment(restoredEdit.fasta, restoredEdit.source, true, [1], validAlignmentFrameOffset(restoredEdit.frameOffset) ?? savedFrameOffset, alignmentProductiveOnly);
          const key=lineageAlignmentKey([1],alignmentProductiveOnly);
          setEditedAlignments(new Map([[key, { ...restoredEdit, key, lineageIds: [1], productiveOnly:alignmentProductiveOnly, frameOffset: validAlignmentFrameOffset(restoredEdit.frameOffset) }]]));
        } else if (initialSession?.tree?.run && typeof initialSession.tree.run.alignmentFasta === "string") {
          installAlignment(initialSession.tree.run.alignmentFasta, initialSession.tree.source || "Saved tree input", false, [1], validAlignmentFrameOffset(initialSession.tree.run.frameOffset) ?? savedFrameOffset);
        } else if ((alignmentProductiveOnly ? rows.filter(isProductiveLineageRow) : rows).length >= 2) {
          const workspaceRows = alignmentProductiveOnly ? rows.filter(isProductiveLineageRow) : rows;
          const requested = initialSession?.lineageTemplateSelections?.[lineageAlignmentKey([1], alignmentProductiveOnly)];
          const selected = inferLineageGermline(workspaceRows, lineageGermlineMethod, { references, templateOrdinal: workspaceRows.some((row) => row.record.ordinal === requested) ? requested : undefined });
          const selectedRows = retainLineageTemplate(workspaceRows, stratifiedLineageRows(workspaceRows, lineageByOrdinal, Math.min(200, workspaceRows.length)), selected.selectedOrdinal);
          const input = lineageInputFasta(selectedRows, lineageGermlineMethod, { references, templateOrdinal: selected.selectedOrdinal });
          const records = parseFasta(input.fasta, true);
          const maximum = Math.max(...records.map((record) => record.sequence.length));
          installAlignment(records.map((record) => `>${record.name}\n${record.sequence.padEnd(maximum, "-")}`).join("\n") + "\n", "AIRR-anchored reference quick view", false, [1], input.alignmentFrameOffset);
        } else {
          setError("This input contains only one sequence. V(D)J annotation succeeded, but an MSA needs at least two biological sequences.");
        }
        if (initialSession?.tree?.run) setTreeRun(initialSession.tree.run as unknown as TreeSnapshot);
        if (initialSession?.phyloUca) setPhyloUcaState(initialSession.phyloUca);
      } catch (loadError) {
        if (!cancelled) setError(loadError instanceof Error ? loadError.message : String(loadError));
      } finally {
        if (!cancelled) setBusy("");
      }
    })();
    return () => { cancelled = true; };
  }, [directLineage, initialSession, lineageGermlineMethod, store]);

  useEffect(()=>{
    if(directLineage||!initialSession||restoredSessionRef.current)return;restoredSessionRef.current=true;
    void (async()=>{
      setBusy("Restoring saved post-analysis state");setError("");
      postRestoreActiveRef.current=true;runtime.beginTransaction();
      try{
        await runInActiveLock(async(signal)=>{
        const active=initialSession.activeMask?unpackSessionVector(initialSession.activeMask) as Uint8Array:null;
        if(initialSession.lineageGermlineMethod==="closest"||initialSession.lineageGermlineMethod==="consensus")setLineageGermlineMethod(initialSession.lineageGermlineMethod);
        setAlignmentProductiveOnly(Boolean(initialSession.alignmentProductiveOnly));
        const collapse=initialSession.collapse;const lineage=initialSession.lineage;
        const dedupState=collapse?.counts&&collapse.representatives&&collapse.dashboard?{dashboard:collapse.dashboard as unknown as DedupDashboard,counts:unpackSessionVector(collapse.counts) as Uint32Array,representatives:unpackSessionVector(collapse.representatives) as Int32Array}:undefined;
        const lineageState=lineage?.assignments&&lineage.dashboard?{dashboard:lineage.dashboard as unknown as LineageDashboard,assignments:unpackSessionVector(lineage.assignments) as Int32Array}:undefined;
        let restoredAllele:AlleleRefinementResult|null=null;let restoredAlleleApplied=false;let restoredAllelePolicy:AlleleReassignmentPolicy="confidence";let restoredAlleleThreshold=0.8;
        if(initialSession.alleleRefinement){restoredAllele=restoreAlleleRefinement(initialSession.alleleRefinement);restoredAllele.options={...DEFAULT_ALLELE_REFINEMENT_OPTIONS,...restoredAllele.options,segments:[...restoredAllele.options.segments]};restoredAlleleApplied=Boolean(initialSession.alleleRefinement.applied);restoredAllelePolicy=initialSession.alleleRefinement.reassignmentPolicy??"confidence";restoredAlleleThreshold=initialSession.alleleRefinement.applyMinimumPosterior??0.8;setAlleleRefinement(restoredAllele);setAlleleOptions(restoredAllele.options);setAlleleApplied(restoredAlleleApplied);setAlleleReassignmentPolicy(restoredAllelePolicy);setAlleleApplyMinimumPosterior(restoredAlleleThreshold);if(restoredAlleleApplied)await runtime.setRepertoireCallOverrides(restoredAllele,restoredAllelePolicy,restoredAlleleThreshold);}
        const restoredRuntimeState=await runtime.restoreState({activeMask:active,dedup:dedupState,lineages:lineageState});
        setWorkingMask(active);setWorkingStages(initialSession.workingStages as WorkingSetStage[]);
        if(collapse?.dashboard){setDedup(collapse.dashboard as unknown as DedupDashboard);setCollapseMode(collapse.mode);const o=collapse.options; if(typeof o.dedupKey==="string")setDedupKey(o.dedupKey as DedupKey);if(typeof o.collapseScope==="string")setCollapseScope(o.collapseScope as DatasetScope);if(typeof o.respectConstantCall==="boolean")setRespectConstantCall(o.respectConstantCall);if(typeof o.denoiseUnresolvedPolicy==="string")setDenoiseUnresolvedPolicy(o.denoiseUnresolvedPolicy as "discard"|"retain");}
        if(initialSession.selection){setSelectionDraft({...DEFAULT_REPERTOIRE_SELECTION,...initialSession.selection.options});setSelectionApplied(initialSession.workingStages.some(stage=>stage.id==="selection"));if(initialSession.selection.baseMask)setSelectionBaseMask(unpackSessionVector(initialSession.selection.baseMask) as Uint8Array);if(initialSession.selection.mask){const mask=unpackSessionVector(initialSession.selection.mask) as Uint8Array;let retained=0;for(const value of mask)retained+=value?1:0;setSelectionPreview({mask,inputRecords:initialSession.workingStages.find(stage=>stage.id==="selection")?.input??store.count,retainedRecords:retained,discardedRecords:(initialSession.workingStages.find(stage=>stage.id==="selection")?.input??store.count)-retained,summary:"restored saved selection"});}}
        const rawRestoredLineages=restoredRuntimeState.lineages??(lineage?.dashboard?lineage.dashboard as unknown as LineageDashboard:null);
        const restoredLineages=rawRestoredLineages?{...rawRestoredLineages,summaries:rawRestoredLineages.summaries.map((summary)=>({...summary,studyScope:summary.studyScope??"global",studyGroup:summary.studyGroup||"complete study",sampleIds:summary.sampleIds??[],sampleCounts:summary.sampleCounts??(summary.sampleIds??[]).map((sampleId)=>({sampleId,uniqueMembers:(summary.sampleIds??[]).length===1?summary.uniqueMembers:0,abundance:(summary.sampleIds??[]).length===1?summary.abundance:0})),subjectIds:summary.subjectIds??[],timepoints:summary.timepoints??[],compartments:summary.compartments??[],doubleDPositiveMembers:summary.doubleDPositiveMembers??0,doubleDPositiveAbundance:summary.doubleDPositiveAbundance??0}))}:null;
        if(restoredLineages){setLineages(restoredLineages);const o=lineage?.options??{};if(typeof o.identity==="number")setIdentity(o.identity);if(typeof o.resolution==="string")setResolution(o.resolution as CallResolution);if(typeof o.ambiguity==="string")setAmbiguity(o.ambiguity as AmbiguityPolicy);if(typeof o.productiveOnly==="boolean")setProductiveOnly(o.productiveOnly);if(typeof o.lineageScope==="string")setLineageScope(o.lineageScope as DatasetScope);}
        const restoredMerges=(initialSession.lineageMerges??[]).map((merge)=>({...merge,originalLineageIds:[...new Set(merge.originalLineageIds)].filter((value)=>value>0).sort((a,b)=>a-b)}));
        setLineageMerges(restoredMerges);
        const restoredEdited=new Map<string,EditedAlignmentState>();
        for(const entry of initialSession.editedAlignments??[]){const productiveOnly=Boolean(entry.productiveOnly);const key=lineageAlignmentKey(entry.lineageIds,productiveOnly);if(key)restoredEdited.set(key,{...entry,key,productiveOnly,lineageIds:[...new Set(entry.lineageIds)].sort((a,b)=>a-b),frameOffset:validAlignmentFrameOffset(entry.frameOffset)});}
        if(initialSession.alignment&&/(corrected|manual|alivibe|edited)/i.test(initialSession.alignment.source)){const lineageIds=initialSession.alignment.selectedLineageId?[initialSession.alignment.selectedLineageId]:[];const productiveOnly=Boolean(initialSession.alignmentProductiveOnly);const key=lineageAlignmentKey(lineageIds,productiveOnly);if(key&&!restoredEdited.has(key))restoredEdited.set(key,{key,lineageIds,productiveOnly,fasta:initialSession.alignment.fasta,source:initialSession.alignment.source,frameOffset:validAlignmentFrameOffset(initialSession.alignment.frameOffset),savedAt:new Date().toISOString()});}
        setEditedAlignments(restoredEdited);
        const chimera=initialSession.chimera;if(chimera?.dashboard&&chimera.probabilities&&chimera.dfr&&chimera.msa){const dashboard={...chimera.dashboard,probabilities:unpackSessionVector(chimera.probabilities) as Float32Array,dfr:unpackSessionVector(chimera.dfr) as Uint16Array} as unknown as ChmmDashboard;const rawOptions=chimera.options;const options=rawOptions as unknown as ChmmRunOptions;const inputMask=chimera.retainedMask?unpackSessionVector(chimera.retainedMask) as Uint8Array:null;setChmm(dashboard);setChmmRun({msa:chimera.msa,options,inputMask});setPreparedMsa(chimera.msa);setChmmFilterThreshold(chimera.filterThreshold);setChmmSegment(options.segment);if(rawOptions.chmmSource==="selected"||rawOptions.chmmSource==="upload")setChmmSource(rawOptions.chmmSource);if(typeof rawOptions.uploadedMsaName==="string")setUploadedMsaName(rawOptions.uploadedMsaName);if(rawOptions.chmmSource==="upload")setUploadedMsa(chimera.msa);}
        if(initialSession.unifiedGermline){setUnifiedOptions({...DEFAULT_UNIFIED_OPTIONS,...initialSession.unifiedGermline.options});setUnifiedGermline(initialSession.unifiedGermline.dashboard?.mode==="joint-inherited-somatic-split"?initialSession.unifiedGermline.dashboard:null);}
        if(initialSession.shm){setShmMetric(initialSession.shm.metric);setShmDashboard(initialSession.shm.dashboard);if(initialSession.shm.sampleOrder?.length)setShmSampleOrder([...initialSession.shm.sampleOrder]);}if(initialSession.missingAlleles){setMissingAlleleOptions({...DEFAULT_MISSING_ALLELE_OPTIONS,...initialSession.missingAlleles.options,unit:"lineage"});setMissingAlleles(initialSession.missingAlleles.dashboard?.validationPasses===2?initialSession.missingAlleles.dashboard:null);setSelectedMissingAlleleIds(new Set(initialSession.missingAlleles.selectedCandidateIds??[]));}if(initialSession.personalizedGermline){setPersonalizedGermlineOptions({...DEFAULT_PERSONALIZED_GERMLINE_OPTIONS,...initialSession.personalizedGermline.options});setPersonalizedGermline(initialSession.personalizedGermline.dashboard?.contextModel==="hs5f-hurdle-profiled-junction"?initialSession.personalizedGermline.dashboard:null);}
        const q=initialSession.query??{};if(typeof q.queryText==="string")setQueryText(q.queryText);if(q.queryResultMode==="lineages"||q.queryResultMode==="sequences")setQueryResultMode(q.queryResultMode);
        if(Array.isArray(q.queryLineageHits))setQueryLineageHits(q.queryLineageHits as QueryLineageHit[]);
        const restoredHits=Array.isArray(q.queryHits)?q.queryHits as QueryHit[]:[];setQueryHits(restoredHits);const restoredExpansion=q.expanded as NonNullable<typeof expanded>|undefined;if(restoredExpansion)setExpanded(restoredExpansion);
        const restoredQueryOrdinals=[...new Set([...(restoredExpansion?.ordinals??[]),...restoredHits.map(hit=>hit.ordinal)])].slice(0,500);
        if(restoredQueryOrdinals.length){const rows=await store.indexRecords(restoredQueryOrdinals);setQueryRecords(new Map(rows.map(row=>[row.ordinal,row])));}
        const restoredSelectedIds=[...new Set(initialSession.selectedLineageIds??initialSession.tree?.lineageIds??(initialSession.alignment?.selectedLineageId?[initialSession.alignment.selectedLineageId]:[]))].filter((value)=>value>0).sort((a,b)=>a-b);
        if(restoredSelectedIds.length){
          openModule("workbench");
          const selectedId=restoredSelectedIds[0];
          const summary=restoredLineages?.summaries.find(item=>item.id===selectedId);
          if(summary){
            const limitPerLineage=restoredSelectedIds.length===1?500:Math.max(20,Math.floor(1000/restoredSelectedIds.length));const memberGroups=await runtime.lineageMembersMany(restoredSelectedIds,limitPerLineage);const rows=refineDetailRows(await store.detailMany(memberGroups.flatMap((group)=>group.ordinals)),restoredAllele,restoredAllelePolicy,restoredAlleleThreshold,restoredAlleleApplied);
            const counts=initialSession.workingStages.some(stage=>stage.id==="dedup")?await runtime.dedupCounts():null;
            const originalByOrdinal=new Map<number,number>();memberGroups.forEach((group)=>group.ordinals.forEach((ordinal)=>originalByOrdinal.set(ordinal,group.lineageId)));
            setSelectedLineage(summary);setSelectedLineageIds(restoredSelectedIds);setLineageRows(rows);setLineageTotal(memberGroups.reduce((sum,group)=>sum+group.total,0));setOriginalLineageByOrdinal(originalByOrdinal);setLineageMultiplicity(new Map(rows.map(row=>{const imported=Number(row.values.duplicate_count);const count=counts?.[row.record.ordinal]||(Number.isFinite(imported)&&imported>0?imported:1);return [row.record.ordinal,Math.max(1,Math.floor(count))] as const;})));
          }
          const restoredProductiveOnly=Boolean(initialSession.alignmentProductiveOnly);const key=lineageAlignmentKey(restoredSelectedIds,restoredProductiveOnly);const restored=restoredEdited.get(key);
          const savedFrameOffset=validAlignmentFrameOffset(initialSession.alignmentFrameOffset);
          if(restored)installAlignment(restored.fasta,restored.source,true,restored.lineageIds,restored.frameOffset??savedFrameOffset,restoredProductiveOnly);
          else if(initialSession.tree?.run&&typeof initialSession.tree.run.alignmentFasta==="string")installAlignment(initialSession.tree.run.alignmentFasta,initialSession.tree.source||"Saved tree input",false,restoredSelectedIds,validAlignmentFrameOffset(initialSession.tree.run.frameOffset)??savedFrameOffset,restoredProductiveOnly);
        }
        if(initialSession.tree?.run)setTreeRun(initialSession.tree.run as unknown as TreeSnapshot);if(initialSession.phyloUca)setPhyloUcaState(initialSession.phyloUca);
        if(signal.aborted)throw new DOMException("Session restoration was cancelled.","AbortError");
        });runtime.commitTransaction();postRestoreActiveRef.current=false;
      }catch(restoreError){if(!cancellationRecoveryRef.current){if(runtime.requiresRecoveryForCancellation())await runtime.rollbackTransaction((processed,total)=>setProgress({processed,total,unit:"AIRR records restored"}));else runtime.commitTransaction();setError(restoreError instanceof Error?restoreError.message:String(restoreError));}}finally{if(!cancellationRecoveryRef.current)setBusy("");postRestoreActiveRef.current=false;}
    })();
  },[directLineage,initialSession,runtime,store.count]);

  async function cancelActiveOperation() {
    const controller = postLockAbortRef.current;
    if (!controller || cancellationRecoveryRef.current) return;
    const recoverRuntime = runtime.requiresRecoveryForCancellation();
    const cancelledPipeline = pipelineActiveRef.current;
    const cancelledRestore = postRestoreActiveRef.current;
    setCancelling(true);
    setBusy(recoverRuntime ? "Cancelling · restoring the last committed analysis state" : "Cancelling active computation");
    setError("");
    const recovery = (async () => {
      controller.abort();
      alleleRuntime.cancel();
      personalizedGermlineRuntime.cancel();
      unifiedGermlineRuntime.cancel();
      if (recoverRuntime) {
        await runtime.cancelAndRestore((processed, total) => setProgress({ processed, total, unit: "AIRR records restored" }), true);
      } else {
        runtime.commitTransaction();
      }
      if (cancelledPipeline || cancelledRestore) {
        setAlleleRefinement(null);
        setAlleleApplied(false);
        clearDownstreamStageState();
        invalidateAssignmentDependentAnalyses();
        if (cancelledPipeline) setPipelineReport(["Pipeline cancelled. No pipeline result or working-set change was committed."]);
      }
    })();
    cancellationRecoveryRef.current = recovery;
    try {
      await recovery;
    } catch (recoveryError) {
      setError(`Cancellation succeeded, but restoring the previous worker state failed: ${recoveryError instanceof Error ? recoveryError.message : String(recoveryError)}`);
    } finally {
      if (cancellationRecoveryRef.current === recovery) cancellationRecoveryRef.current = null;
      setCancelling(false);
      setBusy("");
    }
  }

  async function runInActiveLock<T>(action: (signal: AbortSignal) => Promise<T>): Promise<T> {
    const controller = new AbortController();
    postLockAbortRef.current = controller;
    try {
      return await withAnalysisWebLock(controller.signal, setPostLockState, () => action(controller.signal));
    } finally {
      if (postLockAbortRef.current === controller) postLockAbortRef.current = null;
      setPostLockState("unsupported");
    }
  }

  async function operation<T>(label: string, action: (signal: AbortSignal) => Promise<T>): Promise<T | undefined> {
    runtime.beginTransaction();
    setBusy(label);
    setError("");
    try {
      const result = await runInActiveLock(async (signal) => {
        await runtime.ensureIndexed((processed, total) => setProgress({ processed, total }), signal);
        return action(signal);
      });
      runtime.commitTransaction();
      return result;
    } catch (operationError) {
      if (!cancellationRecoveryRef.current) {
        if (runtime.requiresRecoveryForCancellation()) {
          setBusy("Restoring the last committed analysis state");
          await runtime.rollbackTransaction((processed, total) => setProgress({ processed, total, unit: "AIRR records restored" }));
        } else runtime.commitTransaction();
      }
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
      return undefined;
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function runDedup() {
    setProgress({ processed: 0, total: store.count });
    const label = collapseMode === "exact" ? "Deduplicating AIRR records" : collapseMode === "fad" ? "Running FAD-compatible denoising" : collapseMode === "indel" ? "Running indel-aware denoising" : "Running conservative error-model denoising";
    const result = await operation(label, async (signal) => {
      // Repertoire allele pooling is the upstream assignment stage. Collapse
      // therefore uses whichever policy-selected calls are currently applied.
      return collapseMode === "exact" ? runtime.deduplicate(dedupKey, denoiseUnresolvedPolicy, collapseScope, respectConstantCall) : runtime.denoise({
      mode: collapseMode,
      errorRate: denoiseErrorRate,
      alpha: denoiseAlpha,
      callResolution: denoiseResolution,
      ambiguity: denoiseAmbiguity,
      minimumParentCount,
      ambiguousPolicy: denoiseAmbiguousPolicy,
      unresolvedPolicy: denoiseUnresolvedPolicy,
      fadNeighborThreshold,
      fadMethod,
      expectedZeroErrorFraction,
      maximumHammingDistance: maximumDenoiseDistance,
      maximumEditDistance,
      minimumIndelParentRatio,
      maxCandidatesPerVariant: denoiseCandidateCap,
      scope: collapseScope,
      respectConstantCall,
    }, (processed, total, phase) => {
      const stage = phase === "variants" ? "Denoising indexed sequence variants" : phase === "finalize" ? "Materializing representatives and multiplicities" : "Streaming VDJ sequences into the denoising index";
      setBusy(stage);
      setProgress({ processed, total, unit: phase === "variants" ? "unique sequence variants processed" : phase === "finalize" ? "representative-state operations completed" : "AIRR records streamed into the denoising index" });
      }, signal);
    });
    if (result) {
      setDedup(result);
      setWorkingMask(null);
      setWorkingStages([]);
      setSelectionPreview(null);
      setSelectionBaseMask(null);
      setSelectionApplied(false);
      setChmm(null);
      setChmmRun(null);
      setChimeraDetail(null);
      invalidateAssignmentDependentAnalyses();
    }
  }

  async function applyDedupFilter() {
    if (!dedup) return;
    const result = await operation("Applying deduplicated representatives to the downstream working set", () => runtime.applyDedupFilter());
    if (!result) return;
    setWorkingMask(result.mask);
    setWorkingStages([{
      id: "dedup",
      label: dedup.mode === "exact" ? `Exact collapse · ${dedup.key}` : dedup.mode === "fad" ? "FAD denoising" : dedup.mode === "indel" ? "Indel-aware denoising" : "Conservative denoising",
      input: dedup.inputRecords,
      retained: result.retained,
      discarded: dedup.collapsedRecords,
      detail: `${dedup.algorithm}. Collapsed abundance remains in duplicate_count and lineage weights.${dedup.unresolvedRecords ? ` ${dedup.unresolvedRecords.toLocaleString()} ineligible records ${denoiseUnresolvedPolicy === "discard" ? "were excluded" : "remain unchanged"}.` : ""}`,
    }]);
    setSelectionPreview(null);
    setSelectionBaseMask(null);
    setSelectionApplied(false);
    setChmm(null);
    setChmmRun(null);
    setChimeraDetail(null);
    invalidateAssignmentDependentAnalyses();
    advanceModule("dedup","chimera");
  }

  async function resetWorkingSet() {
    const result = await operation("Restoring the complete downstream working set", () => runtime.setActiveMask(null));
    if (!result) return;
    setWorkingMask(null);
    setWorkingStages([]);
    setSelectionPreview(null);
    setSelectionBaseMask(null);
    setSelectionApplied(false);
    setChmm(null);
    setChmmRun(null);
    setChimeraDetail(null);
    invalidateAssignmentDependentAnalyses();
    setActiveWorkspace("dedup");
  }

  function invalidateAssignmentDependentAnalyses() {
    setLineages(null);
    setLineageTableSelectedIds(new Set());
    setSelectedLineage(null);
    setSelectedLineageIds([]);
    setLineageRows([]);
    setAlignmentProductiveOnly(false);
    setLineageMultiplicity(new Map());
    setOriginalLineageByOrdinal(new Map());
    setEditedAlignments(new Map());
    setLineageMerges([]);
    clearNeighbourResults(true);
    clearAlignmentArtifacts();
    setQueryHits([]);
    setQueryLineageHits([]);
    setQueryRecords(new Map());
    setExpanded(null);
    setShmDashboard(null);
    setMissingAlleles(null);
    setSelectedMissingAlleleIds(new Set());
    setPersonalizedGermline(null);
    setUnifiedGermline(null);
  }

  function clearDownstreamStageState() {
    setDedup(null);
    setWorkingMask(null);
    setWorkingStages([]);
    setSelectionPreview(null);
    setSelectionBaseMask(null);
    setSelectionApplied(false);
    setChmm(null);
    setChmmRun(null);
    setChimeraDetail(null);
  }

  function discardAppliedAllelePolicy() {
    if (!alleleApplied) return;
    setAlleleApplied(false);
    clearDownstreamStageState();
    invalidateAssignmentDependentAnalyses();
    void (async () => {
      try {
        await runtime.setActiveMask(null);
        await runtime.setRepertoireCallOverrides(null);
      } catch (runtimeError) {
        setError(runtimeError instanceof Error ? runtimeError.message : String(runtimeError));
      }
    })();
  }

  function overlayRefinedCalls(row: { ordinal: number; values: Record<string,string> }) {
    if (!alleleRefinement || !alleleApplied) return;
    const v = refinedCall(alleleRefinement, "V", row.ordinal, alleleReassignmentPolicy, alleleApplyMinimumPosterior);
    const d = refinedCall(alleleRefinement, "D", row.ordinal, alleleReassignmentPolicy, alleleApplyMinimumPosterior);
    const j = refinedCall(alleleRefinement, "J", row.ordinal, alleleReassignmentPolicy, alleleApplyMinimumPosterior);
    if (v) row.values.v_call = v;
    if (d) row.values.d_call = d;
    if (j) row.values.j_call = j;
  }

  async function previewSelection() {
    const validation = validateRepertoireSelection(selectionDraft);
    if (validation.length) { setError(validation.join(" ")); return; }
    const result = await operation("Previewing the repertoire selection", async (signal) => {
      const current = await runtime.activeMask();
      const base = selectionApplied ? selectionBaseMask : current;
      const snapshot = base ? base.slice() : null;
      const preview = await selectRepertoire(store, selectionDraft, snapshot ?? undefined, (processed) => setProgress({ processed, total: snapshot ? snapshot.reduce((sum, value) => sum + value, 0) : store.count }), signal);
      return { preview, base: snapshot };
    });
    if (!result) return;
    setSelectionPreview(result.preview);
    setSelectionBaseMask(result.base);
  }

  async function applySelection() {
    if (!selectionPreview) return;
    const result = await operation("Committing the repertoire selection", () => runtime.setActiveMask(selectionPreview.mask));
    if (!result) return;
    setWorkingMask(selectionPreview.mask);
    setSelectionApplied(true);
    setWorkingStages((stages) => [...stages.filter((stage) => stage.id !== "selection"), {
      id: "selection", label: `Repertoire selection · ${selectionPreview.summary}`, input: selectionPreview.inputRecords,
      retained: selectionPreview.retainedRecords, discarded: selectionPreview.discardedRecords,
      detail: "Explicitly committed selection; downstream lineage, SHM, reference evidence, and sequence queries inherit this mask.",
    }]);
    invalidateAssignmentDependentAnalyses();
    advanceModule("selection","lineage");
  }

  async function removeSelection() {
    const result = await operation("Removing the committed repertoire selection", () => runtime.setActiveMask(selectionBaseMask));
    if (!result) return;
    setWorkingMask(selectionBaseMask);
    setWorkingStages((stages) => {
      const index=stages.findIndex((stage)=>stage.id==="selection");
      return index<0?stages:stages.slice(0,index);
    });
    setChmm(null);setChmmRun(null);setChimeraDetail(null);
    setSelectionApplied(false);
    setSelectionPreview(null);
    setSelectionBaseMask(null);
    invalidateAssignmentDependentAnalyses();
  }

  async function runAlleleRefinement() {
    setAlleleProgress({ processed: 0, total: store.count, phase: "Preparing repertoire allele evidence" });
    const result = await operation("Pooling ambiguous germline assignments across the repertoire", async () => {
      // This first-stage fit always sees the complete assigned input. A later
      // collapse or selection mask must not alter donor-level allele evidence.
      const fitted = await alleleRuntime.run(store, references, alleleOptions, null, (next) => {
        setAlleleProgress(next);
        setBusy(next.phase);
        setProgress({ processed: next.processed, total: next.total, unit: "records or independent repertoire pools" });
      }, postLockAbortRef.current?.signal);
      if (alleleApplied) {
        await runtime.setActiveMask(null);
        await runtime.setRepertoireCallOverrides(null);
      }
      return fitted;
    });
    setAlleleProgress(null);
    if (!result) return;
    if (alleleApplied) clearDownstreamStageState();
    setAlleleRefinement(result);
    setAlleleApplied(false);
    invalidateAssignmentDependentAnalyses();
    openModule("alleles");
  }

  async function applyAlleleRefinement() {
    if (!alleleRefinement) return;
    const result = await operation("Applying repertoire posterior calls to downstream analyses", async () => {
      await runtime.setActiveMask(null);
      await runtime.setRepertoireCallOverrides(alleleRefinement, alleleReassignmentPolicy, alleleApplyMinimumPosterior);
      return true;
    });
    if (!result) return;
    clearDownstreamStageState();
    setAlleleApplied(true);
    invalidateAssignmentDependentAnalyses();
    openModule("alleles");
  }

  async function resetAlleleRefinement() {
    const result = await operation("Restoring original AIRR germline calls downstream", async () => {
      await runtime.setActiveMask(null);
      await runtime.setRepertoireCallOverrides(null);
      return true;
    });
    if (!result) return;
    clearDownstreamStageState();
    setAlleleApplied(false);
    invalidateAssignmentDependentAnalyses();
    openModule("alleles");
  }

  function downloadAlleleModel() {
    if (!alleleRefinement) return;
    const extension = tableExtension(exportFormat);
    downloadText(modelSummaryTable(alleleRefinement, exportFormat), `${baseName(inputName)}.repertoire-allele-model${extension}`);
  }

  async function downloadAlleleSidecar() {
    if (!alleleRefinement) return;
    setBusy("Writing per-record repertoire allele posterior sidecar");
    setError("");
    try {
      await runInActiveLock(async (signal) => {
        const mask = await runtime.activeMask();
        const extension = tableExtension(exportFormat);
        await saveStream(`${baseName(inputName)}.repertoire-allele-posteriors${extension}`, "Sparse per-record repertoire allele posterior", extension, (writer) => writeRefinementSidecar(store, alleleRefinement, alleleReassignmentPolicy, alleleApplyMinimumPosterior, exportFormat, writer.write, mask ?? undefined, signal),exportCompression);
      });
    } catch (operationError) { if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError)); }
    finally { if (!cancellationRecoveryRef.current) setBusy(""); }
  }

  async function downloadRefinedAirr() {
    if (!alleleRefinement) return;
    setBusy("Writing AIRR table with repertoire-refined germline calls");
    setError("");
    try {
      await runInActiveLock(async (signal) => {
        const mask = await runtime.activeMask();
        const extension = tableExtension(exportFormat);
        await saveStream(`${baseName(inputName)}.repertoire-refined.airr${extension}`, "AIRR table with policy-selected repertoire allele calls", extension, (writer) => writeRefinedAirr(store, alleleRefinement, alleleReassignmentPolicy, alleleApplyMinimumPosterior, exportFormat, writer.write, mask ?? undefined, signal),exportCompression);
      });
    } catch (operationError) { if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError)); }
    finally { if (!cancellationRecoveryRef.current) setBusy(""); }
  }

  async function downloadActivePopulation() {
    setBusy("Writing the current selected population");setError("");
    try { await runInActiveLock(async(signal)=>{const mask=await runtime.activeMask();const extension=tableExtension(exportFormat);await saveStream(`${baseName(inputName)}.selected.airr${extension}`,"Selected AIRR population",extension,(writer)=>store.writeAirrFormat(exportFormat,writer.write,mask??undefined,signal),exportCompression);}); }
    catch(operationError){if(!isAbortError(operationError))setError(operationError instanceof Error?operationError.message:String(operationError));}finally{if(!cancellationRecoveryRef.current)setBusy("");}
  }

  async function runShmAnalysis() {
    const result=await operation("Calculating somatic hypermutation on the current working set",async(signal)=>{
      const mask=await runtime.activeMask();const assignments=lineages?await analysisLineageAssignments():null;const counts=workingStages.some((stage)=>stage.id==="dedup")?await runtime.dedupCounts():null;
      const accumulator=new ShmAccumulator({metric:shmMetric,maxSamplesPerLineage:shmSampleCap});
      const fields=["sequence_id","v_call","j_call","locus","isotype","sample_id","subject_id","swig_cohort","swig_timepoint","swig_compartment","duplicate_count","v_sequence_alignment","v_germline_alignment","sequence_alignment","germline_alignment","v_sequence_start","v_germline_start","sequence_frame","v_frame","cdr1_start","cdr1_end","cdr2_start","cdr2_end","fwr1_start","fwr1_end","fwr2_start","fwr2_end","fwr3_start","fwr3_end"];
      await store.scanAirrRows(fields,async(rows)=>{for(const row of rows){overlayRefinedCalls(row);if(counts?.[row.ordinal])row.values.duplicate_count=String(counts[row.ordinal]);const stratum=shmStratum==="all"?"All selected":row.values[shmStratum]||"Unassigned";accumulator.add(row.values,row.ordinal,assignments?.[row.ordinal]??0,stratum);}}, {batchSize:1500,includeMask:mask??undefined,onProgress:(processed,total)=>setProgress({processed,total}),signal});
      return accumulator.finish();
    });
    if(result)setShmDashboard(result);
  }

  async function runMissingAlleleAnalysis() {
    if(!lineages){setError("Missing-V screening requires lineage assignments on the current selected population. Assign lineages first so clonal descendants cannot count as independent evidence.");return;}
    const result=await operation("Two-pass screening for linked V germline discrepancies",async(signal)=>{
      const mask=await runtime.activeMask();const assignments=await analysisLineageAssignments();const accumulator=new MissingAlleleAccumulator(missingAlleleOptions);
      const fields=["subject_id","v_call","j_call","cdr3","junction","v_germline_start","v_sequence_alignment","v_germline_alignment"];
      setBusy("Missing-V screen · discovery representatives (pass 1 of 2)");
      await store.scanAirrRows(fields,async(rows)=>{for(const row of rows){overlayRefinedCalls(row);accumulator.add(row.values,row.ordinal,assignments[row.ordinal]??0);}},{batchSize:1500,includeMask:mask??undefined,onProgress:(processed,total)=>setProgress({processed,total:total*2}),signal});
      const validator=accumulator.prepareValidation(references.V);
      setBusy("Missing-V screen · all-member reference veto (pass 2 of 2)");
      await store.scanAirrRows(fields,async(rows)=>{for(const row of rows){overlayRefinedCalls(row);validator.add(row.values,row.ordinal,assignments[row.ordinal]??0);}},{batchSize:1500,includeMask:mask??undefined,onProgress:(processed,total)=>setProgress({processed:total+processed,total:total*2}),signal});
      return validator.finish();
    });
    if(result){setMissingAlleles(result);setSelectedMissingAlleleIds(new Set());}
  }

  async function runPersonalizedGermlineAnalysis() {
    if(!lineages){setError("Personalized germline inference requires lineage assignments on the current selected population. Assign lineages first so every clone contributes exactly one observation.");return;}
    const result=await operation("Inferring a lineage-weighted personalized V germline set",async(signal)=>{
      const mask=await runtime.activeMask();
      const assignments=await analysisLineageAssignments();
      await personalizedGermlineRuntime.begin(references.V,personalizedGermlineOptions);
      const fields=["sequence_id","sequence","rev_comp","subject_id","locus","v_call","v_sequence_start","v_sequence_end","d_sequence_start","j_sequence_start","v_germline_start","v_sequence_alignment","v_germline_alignment"];
      setBusy("Personalized V germline · selecting the lowest-current-SHM member of every lineage");
      await store.scanAirrRows(fields,async(rows)=>{
        for(const row of rows)overlayRefinedCalls(row);
        await personalizedGermlineRuntime.ingest(rows,assignments);
      },{batchSize:1500,includeMask:mask??undefined,onProgress:(processed,total)=>setProgress({processed,total,unit:"AIRR records screened"}),signal});
      return personalizedGermlineRuntime.finish((next)=>{setBusy(next.phase);setProgress({processed:next.processed,total:next.total,unit:"V-gene fits"});});
    });
    if(result)setPersonalizedGermline(result);
  }

  async function runUnifiedGermlineAnalysis() {
    if(!lineages){setError("Personalized germline inference requires lineage assignments on the current selected population. Assign lineages first so every clone contributes exactly one observation.");return;}
    const result=await operation("Testing a joint inherited/SHM model on held-out lineages",async(signal)=>{
      const mask=await runtime.activeMask();
      const assignments=await analysisLineageAssignments();
      await unifiedGermlineRuntime.begin(references.V,unifiedOptions);
      const fields=["sequence_id","sequence","rev_comp","subject_id","locus","v_call","v_sequence_start","v_sequence_end","d_sequence_start","j_sequence_start","v_germline_start","v_sequence_alignment","v_germline_alignment"];
      setBusy("Joint inherited/SHM model · selecting the lowest-current-SHM member of every lineage");
      await store.scanAirrRows(fields,async(rows)=>{
        for(const row of rows)overlayRefinedCalls(row);
        await unifiedGermlineRuntime.ingest(rows,assignments);
      },{batchSize:1500,includeMask:mask??undefined,onProgress:(processed,total)=>setProgress({processed,total,unit:"AIRR records screened"}),signal});
      return unifiedGermlineRuntime.finish((next)=>{setBusy(next.phase);setProgress({processed:next.processed,total:next.total,unit:"V-gene fits"});});
    });
    if(result)setUnifiedGermline(result);
  }

  async function downloadDeduplicated() {
    setBusy("Writing deduplicated AIRR table");
    setError("");
    try {
      await runInActiveLock(async (signal) => {
        const counts = await runtime.dedupCounts();
        const suffix = dedup?.mode === "exact" ? "deduplicated" : "denoised";
        const extension = tableExtension(exportFormat);
        await saveStream(`${baseName(inputName)}.${suffix}.airr${extension}`, "Collapsed AIRR rearrangement table with multiplicity", extension, async (writer) => store.writeDeduplicatedAirrFormat(counts, exportFormat, writer.write, signal),exportCompression);
      });
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function runLineages() {
    const result = await operation("Assigning lineages", () => runtime.assignLineages({
      identity,
      callResolution: resolution,
      ambiguity,
      productiveOnly,
      requireSameLocus: true,
      maxCandidateComparisons: candidateCap,
      scope: lineageScope,
    }, workingStages.some((stage) => stage.id === "dedup")));
    if (result) {
      setLineages(result);
      setLineageTableSelectedIds(new Set());
      setNeighbourScope(lineageScope);
      setNeighbourCdr3Identity(Math.max(0.5, identity - 0.05));
      setSelectedLineage(null);
      setLineageRows([]);
      setAlignmentProductiveOnly(false);
      setLineageMultiplicity(new Map());
      setSelectedLineageIds([]);
      setOriginalLineageByOrdinal(new Map());
      setEditedAlignments(new Map());
      setLineageMerges([]);
      setQueryLineageHits([]);
      clearNeighbourResults(true);
      clearAlignmentArtifacts();
      // Keep the assignment controls and lineage table visible: the next
      // normal action is to inspect a lineage, not to jump to diagnostics.
      openModule("lineage");
    }
  }

  async function downloadLineages() {
    setBusy("Writing AIRR table with clone identifiers");
    setError("");
    try {
      await runInActiveLock(async (signal) => {
        const assignments = await runtime.lineageAssignments();
        const extension = tableExtension(exportFormat);
        await saveStream(`${baseName(inputName)}.lineages.airr${extension}`, "AIRR rearrangement table with original and merged lineage identifiers", extension, async (writer) => store.writeLineageAirrFormat(assignments, exportFormat, writer.write, mergedIdByOriginal, signal),exportCompression);
      });
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  function expandedMergedLineageIds(values: Iterable<number>): number[] {
    const result = new Set<number>();
    for (const lineageId of values) {
      const merge = lineageMerges.find((entry) => entry.originalLineageIds.includes(lineageId));
      if (merge) merge.originalLineageIds.forEach((value) => result.add(value));
      else if (lineageId > 0) result.add(lineageId);
    }
    return [...result].sort((left, right) => left - right);
  }

  async function analysisLineageAssignments(): Promise<Int32Array> {
    const original = await runtime.lineageAssignments();
    if (!lineageMerges.length) return original;
    const result = original.slice();
    let maximum = 0;
    for (const value of original) maximum = Math.max(maximum, value);
    const synthetic = new Map<string, number>();
    lineageMerges.forEach((merge, index) => synthetic.set(merge.id, maximum + index + 1));
    for (let ordinal = 0; ordinal < result.length; ordinal += 1) {
      const merged = mergedIdByOriginal.get(original[ordinal]);
      if (merged) result[ordinal] = synthetic.get(merged) ?? original[ordinal];
    }
    return result;
  }

  function clearNeighbourResults(clearIndex = false) {
    setNeighbourResult(null);
    setGermlineNeighbourScores([]);
    setSelectedNeighbourIds(new Set());
    if (clearIndex) setGermlineSketchIndex(null);
  }

  async function loadLineageGroup(summary: LineageSummary, requestedIds: number[]) {
    const lineageIds = expandedMergedLineageIds(requestedIds);
    if (!lineageIds.length) return;
    setBusy(lineageIds.length === 1 ? `Loading lineage ${lineageIds[0]}` : `Loading ${lineageIds.length} lineages together`);
    setError("");
    try {
      const loaded = await runInActiveLock(async (signal) => {
        const limitPerLineage = lineageIds.length === 1 ? 500 : Math.max(20, Math.floor(1_000 / lineageIds.length));
        const memberGroups = await runtime.lineageMembersMany(lineageIds, limitPerLineage);
        const ordinals = memberGroups.flatMap((group) => group.ordinals);
        const rows = refineDetailRows(await store.detailMany(ordinals), alleleRefinement, alleleReassignmentPolicy, alleleApplyMinimumPosterior, alleleApplied);
        const lineageByOrdinal = new Map<number, number>();
        memberGroups.forEach((group) => group.ordinals.forEach((ordinal) => lineageByOrdinal.set(ordinal, group.lineageId)));
        const deduplicationApplied = workingStages.some((stage) => stage.id === "dedup");
        const counts = deduplicationApplied ? await runtime.dedupCounts() : null;
        if (signal.aborted) throw new DOMException("Lineage loading was cancelled.", "AbortError");
        const multiplicity = new Map(rows.map((row) => {
          const imported = Number(row.values.duplicate_count);
          const value = counts?.[row.record.ordinal] || (Number.isFinite(imported) && imported > 0 ? imported : 1);
          return [row.record.ordinal, Math.max(1, Math.floor(value))] as const;
        }));
        return { memberGroups, rows, lineageByOrdinal, multiplicity };
      });
      setSelectedLineage(summary);
      setSelectedLineageIds(lineageIds);
      setLineageRows(loaded.rows);
      setAlignmentProductiveOnly(false);
      setOriginalLineageByOrdinal(loaded.lineageByOrdinal);
      setLineageMultiplicity(loaded.multiplicity);
      setLineageTotal(loaded.memberGroups.reduce((sum, group) => sum + group.total, 0));
      clearAlignmentArtifacts();
      const restored = editedAlignments.get(lineageAlignmentKey(lineageIds, false));
      if (restored) installAlignment(restored.fasta, restored.source, true, restored.lineageIds, restored.frameOffset, false);
      setAlignmentEditorStatus("");
      setAlignmentEditorError("");
      clearNeighbourResults();
      openModule("workbench");
      window.requestAnimationFrame(() => workbenchRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function openLineage(summary: LineageSummary) {
    await loadLineageGroup(summary, [summary.id]);
  }

  function setLineageTableSelected(lineageId: number, selected: boolean) {
    setLineageTableSelectedIds((current) => {
      const next = new Set(current);
      if (selected) next.add(lineageId); else next.delete(lineageId);
      return next;
    });
  }

  async function openSelectedLineagesTogether() {
    const lineageIds = [...lineageTableSelectedIds].filter((value) => value > 0).sort((left, right) => left - right);
    if (!lineageIds.length || !lineages) return;
    const anchor = lineages.summaries.find((summary) => summary.id === lineageIds[0]);
    if (!anchor) return;
    await loadLineageGroup(anchor, lineageIds);
  }

  function setNeighbourSelected(lineageId: number, selected: boolean) {
    setSelectedNeighbourIds((current) => {
      const next = new Set(current);
      if (selected) next.add(lineageId); else next.delete(lineageId);
      return next;
    });
  }

  async function searchLineageNeighbours() {
    if (!lineages || !selectedLineage || !selectedLineageIds.length) return;
    setBusy("Searching for neighbouring lineages");
    setError("");
    try {
      await runInActiveLock(async (signal) => {
      let nextNeighbourResult: LineageNeighbourResult | null = null;
      let nextGermlineScores: GermlineNeighbourScore[] = [];
      let nextSketchIndex = germlineSketchIndex;
      if (neighbourMethod === "cdr3" || neighbourMethod === "either") {
        const cdr3 = await runtime.lineageNeighbours({
          identity,
          callResolution: resolution,
          ambiguity,
          productiveOnly,
          requireSameLocus: true,
          maxCandidateComparisons: candidateCap,
          scope: neighbourScope,
          sourceLineageIds: selectedLineageIds,
          minimumIdentity: neighbourCdr3Identity,
          maximumResults: Math.max(neighbourLimit * 5, 100),
          metric: neighbourCdr3Metric,
          callPolicy: neighbourCallPolicy,
        }, workingStages.some((stage) => stage.id === "dedup"));
        nextNeighbourResult = cdr3;
      }

      if (neighbourMethod === "germline" || neighbourMethod === "either") {
        const assignments = await runtime.lineageAssignments();
        const activeMask = await runtime.activeMask();
        let index = nextSketchIndex;
        if (!index) {
          setBusy("Indexing one compact germline sketch per lineage");
          index = await buildLineageGermlineSketchIndex(store, assignments, lineages.lineageCount, neighbourScope, activeMask, (processed, total) => setProgress({ processed, total }), signal);
          nextSketchIndex = index;
        }
        const sourceGroups = new Map<string, { lineageId: number; studyGroup: string; rows: AirrDetailRow[] }>();
        for (const row of lineageRows) {
          const lineageId = originalLineageByOrdinal.get(row.record.ordinal) ?? 0;
          if (!lineageId) continue;
          const studyGroup = datasetScopeValue({
            datasetId: row.values.swig_dataset_id || row.record.datasetId,
            sampleId: row.values.sample_id || row.record.sampleId,
            subjectId: row.values.subject_id || row.record.subjectId,
            cohort: row.values.swig_cohort || row.record.cohort,
          }, neighbourScope);
          const key = `${lineageId}\u0000${studyGroup}`;
          const group = sourceGroups.get(key);
          if (group) group.rows.push(row); else sourceGroups.set(key, { lineageId, studyGroup, rows: [row] });
        }
        const screens = new Map<number, { screen: ReturnType<typeof screenLineageGermlineCandidates>[number]; sourceLineageId: number; sourceRows: AirrDetailRow[] }>();
        for (const { lineageId: sourceLineageId, studyGroup, rows } of sourceGroups.values()) {
          const source = inferLineageGermline(rows, lineageGermlineMethod);
          const first = rows[0];
          const candidates = screenLineageGermlineCandidates(index, source, selectedLineageIds, first.values.locus || first.record.locus, studyGroup, Math.min(2_000, Math.max(500, neighbourLimit * 20)));
          for (const screen of candidates) {
            const previous = screens.get(screen.lineageId);
            if (!previous || screen.sketchSimilarity > previous.screen.sketchSimilarity) screens.set(screen.lineageId, { screen, sourceLineageId, sourceRows: rows });
          }
        }
        // Cheap representative-level exact scoring removes weak MinHash hits;
        // only the best shortlist needs multi-member germline reconstruction.
        const representativeRows = refineDetailRows(await store.detailMany([...screens.values()].map((entry) => entry.screen.representativeOrdinal)), alleleRefinement, alleleReassignmentPolicy, alleleApplyMinimumPosterior, alleleApplied);
        const representativeByOrdinal = new Map(representativeRows.map((row) => [row.record.ordinal, row]));
        const ranked = [...screens.values()].map((entry) => {
          const representative = representativeByOrdinal.get(entry.screen.representativeOrdinal);
          const score = representative ? scoreGermlineCandidate(entry.sourceRows, [representative], entry.screen, 1, Math.max(0.7, neighbourGermlineIdentity - 0.12), entry.sourceLineageId, lineageGermlineMethod) : null;
          return score ? { entry, score: score.germlineIdentity } : null;
        }).filter((value): value is NonNullable<typeof value> => Boolean(value))
          .sort((left, right) => right.score - left.score || right.entry.screen.sketchSimilarity - left.entry.screen.sketchSimilarity)
          .slice(0, Math.min(400, Math.max(100, neighbourLimit * 6)));
        const memberGroups = await runtime.lineageMembersMany(ranked.map((value) => value.entry.screen.lineageId), 250);
        const memberByLineage = new Map(memberGroups.map((group) => [group.lineageId, group]));
        const candidateRows = refineDetailRows(await store.detailMany(memberGroups.flatMap((group) => group.ordinals)), alleleRefinement, alleleReassignmentPolicy, alleleApplyMinimumPosterior, alleleApplied);
        const rowsByLineage = new Map<number, AirrDetailRow[]>();
        for (const row of candidateRows) {
          const lineageId = assignments[row.record.ordinal];
          const values = rowsByLineage.get(lineageId);
          if (values) values.push(row); else rowsByLineage.set(lineageId, [row]);
        }
        const exact: GermlineNeighbourScore[] = [];
        for (const { entry } of ranked) {
          const group = memberByLineage.get(entry.screen.lineageId);
          const rows = rowsByLineage.get(entry.screen.lineageId) ?? [];
          if (!group || !rows.length) continue;
          const score = scoreGermlineCandidate(entry.sourceRows, rows, entry.screen, group.total, neighbourGermlineIdentity, entry.sourceLineageId, lineageGermlineMethod);
          if (score) exact.push(score);
        }
        exact.sort((left, right) => right.germlineIdentity - left.germlineIdentity || right.sketchSimilarity - left.sketchSimilarity || left.lineageId - right.lineageId);
        nextGermlineScores = exact.slice(0, Math.max(neighbourLimit * 3, neighbourLimit));
      }
      if (signal.aborted) throw new DOMException("Neighbour searching was cancelled.", "AbortError");
      setNeighbourResult(nextNeighbourResult);
      setGermlineNeighbourScores(nextGermlineScores);
      setSelectedNeighbourIds(new Set());
      if (nextSketchIndex !== germlineSketchIndex) setGermlineSketchIndex(nextSketchIndex);
      });
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  function createMergedLineage(lineageIdsValue: Iterable<number>): LineageMergeState | null {
    const requested = expandedMergedLineageIds(lineageIdsValue);
    if (requested.length < 2) return null;
    const touching = lineageMerges.filter((merge) => merge.originalLineageIds.some((lineageId) => requested.includes(lineageId)));
    const all = new Set(requested);
    touching.forEach((merge) => merge.originalLineageIds.forEach((lineageId) => all.add(lineageId)));
    const retainedId = touching[0]?.id;
    let serial = 1;
    while (!retainedId && lineageMerges.some((merge) => merge.id === `swig_merged_lineage_${serial}`)) serial += 1;
    const id = retainedId ?? `swig_merged_lineage_${serial}`;
    const merged: LineageMergeState = {
      id,
      label: touching[0]?.label ?? `Merged lineage ${serial}`,
      originalLineageIds: [...all].sort((left, right) => left - right),
      createdAt: touching[0]?.createdAt ?? new Date().toISOString(),
    };
    setLineageMerges((current) => [...current.filter((entry) => !touching.some((touch) => touch.id === entry.id)), merged]);
    setShmDashboard(null);
    setMissingAlleles(null);
    setPersonalizedGermline(null);
    setUnifiedGermline(null);
    return merged;
  }

  function removeLineageMerge(id: string) {
    setLineageMerges((current) => current.filter((entry) => entry.id !== id));
    setShmDashboard(null);
    setMissingAlleles(null);
    setPersonalizedGermline(null);
    setUnifiedGermline(null);
  }

  async function viewSelectedNeighbourGroup(merge = false) {
    if (!selectedLineage) return;
    const lineageIds = expandedMergedLineageIds([...selectedLineageIds, ...selectedNeighbourIds]);
    if (lineageIds.length < 2) return;
    if (merge) createMergedLineage(lineageIds);
    await loadLineageGroup(selectedLineage, lineageIds);
  }

  function toggleAlignmentProductiveFilter() {
    const next = !alignmentProductiveOnly;
    setAlignmentProductiveOnly(next);
    clearAlignmentArtifacts();
    setAlignmentEditorError("");
    setAlignmentEditorStatus(next
      ? `Dropped ${excludedNonProductiveRows.toLocaleString()} non-productive or unresolved row${excludedNonProductiveRows===1?"":"s"} from the lineage alignment workspace. Re-run the desired MSA; repertoire data and lineage assignments were not changed.`
      : "Restored the complete loaded lineage population. Re-run the desired MSA.");
  }

  function chooseLineageTemplate(value: string) {
    const ordinal = value === "auto" ? undefined : Number(value);
    setLineageTemplateSelections((current) => {
      const next = { ...current };
      if (ordinal === undefined) delete next[selectedGroupKey]; else next[selectedGroupKey] = ordinal;
      return next;
    });
    setAlignmentEditorError("");
    if (alignment) {
      try {
        const projected = projectLineageTemplate(alignment, workbenchLineageRows, { references, templateOrdinal: ordinal });
        installAlignment(projected.fasta, `${alignmentSource} · selected AIRR germline template`, alignmentEdited, selectedLineageIds, ((alignmentFrameOffset + projected.addedLeftColumns) % 3) as AlignmentFrameOffset);
        setAlignmentEditorStatus("Applied the selected member's trimming and full outer reference flanks to this MSA. Tree and UCA results were invalidated; run them again.");
        return;
      } catch (projectionError) {
        clearAlignmentArtifacts();
        setAlignmentEditorStatus(projectionError instanceof Error ? projectionError.message : String(projectionError));
        return;
      }
    }
    clearAlignmentArtifacts();
    setAlignmentEditorStatus("Template choice saved. Build the MSA to apply its trimming and full outer reference flanks.");
  }

  async function runAlignment() {
    if (workbenchLineageRows.length < 2) {
      setError("A lineage MSA needs at least two biological rows after the current productivity filter.");
      return;
    }
    setBusy(alignmentMethod === "quick" ? "Preparing AIRR-anchored alignment" : alignmentMethod === "codon" ? "Running codon-aware Kalign WASM" : "Running Kalign WASM");
    setError("");
    try {
      const next = await runInActiveLock(async (signal) => {
        const rows = retainLineageTemplate(workbenchLineageRows, stratifiedLineageRows(workbenchLineageRows, originalLineageByOrdinal, Math.max(2, alignmentLimit)), lineageGermline?.selectedOrdinal);
        const input = lineageInputFasta(rows, lineageGermlineMethod, { references, templateOrdinal: lineageGermline?.selectedOrdinal });
        if (alignmentMethod === "quick") {
          const records = parseFasta(input.fasta, true);
          const maximum = Math.max(...records.map((record) => record.sequence.length));
          const fasta = records.map((record) => `>${record.name}\n${record.sequence.padEnd(maximum, "-")}`).join("\n") + "\n";
          return { fasta, frameOffset: input.alignmentFrameOffset };
        }
        if (alignmentMethod === "codon") {
          return { fasta: await runCodonAwareKalignTask(input.fasta, input.frames, signal), frameOffset: 0 as AlignmentFrameOffset };
        }
        const fasta = await runKalignTask(input.fasta, signal);
        const anchor = parseFasta(fasta, true).find((record) => record.name === input.frameAnchorName);
        const derivedFrame = anchor
          ? alignedSequenceFrameOffset(anchor.sequence, input.frameAnchorUngappedOffset)
          : input.alignmentFrameOffset;
        return { fasta, frameOffset: derivedFrame };
      });
      installAlignment(next.fasta, alignmentMethod === "quick" ? "AIRR-anchored reference quick view" : alignmentMethod === "codon" ? "Codon-aware Kalign 3.3.1" : "Nucleotide Kalign 3.3.1", false, selectedLineageIds, next.frameOffset);
      setAlignmentEditorStatus("");
      setAlignmentEditorError("");
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function runDirectAlivibeMsa() {
    if (workbenchLineageRows.length < 2) {
      setError("A lineage MSA needs at least two biological rows after the current productivity filter.");
      return;
    }
    setBusy("Running Alivibe-compatible POA MSA WASM");
    setError("");
    try {
      const next = await runInActiveLock(async (signal) => {
        const rows = retainLineageTemplate(workbenchLineageRows, stratifiedLineageRows(workbenchLineageRows, originalLineageByOrdinal, Math.max(2, alignmentLimit)), lineageGermline?.selectedOrdinal);
        const input = lineageInputFasta(rows, lineageGermlineMethod, { references, templateOrdinal: lineageGermline?.selectedOrdinal });
        const records = parseFasta(input.fasta, true);
        const aligned = await runAlivibeMsaTask(
          records.map((record) => record.sequence.replaceAll("-", "")),
          signal,
          3,
          "nucleotide",
        );
        const fasta = records.map((record, index) => `>${record.name}\n${aligned[index]}`).join("\n") + "\n";
        const anchor = records.findIndex((record) => record.name === input.frameAnchorName);
        const derivedFrame = anchor >= 0
          ? alignedSequenceFrameOffset(aligned[anchor], input.frameAnchorUngappedOffset)
          : input.alignmentFrameOffset;
        return { fasta, frameOffset: derivedFrame };
      });
      installAlignment(
        next.fasta,
        "Alivibe-compatible POA MSA · nucleotide N wildcard · WASM · three refinement passes",
        false,
        selectedLineageIds,
        next.frameOffset,
      );
      setAlignmentEditorStatus("");
      setAlignmentEditorError("");
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function runDirectAlivibeAaMsa() {
    if (workbenchLineageRows.length < 2) {
      setError("A lineage MSA needs at least two biological rows after the current productivity filter.");
      return;
    }
    setBusy("Running frame-preserving amino-acid Alivibe-compatible MSA WASM");
    setError("");
    try {
      const next = await runInActiveLock(async (signal) => {
        const rows = retainLineageTemplate(workbenchLineageRows, stratifiedLineageRows(workbenchLineageRows, originalLineageByOrdinal, Math.max(2, alignmentLimit)), lineageGermline?.selectedOrdinal);
        const input = lineageInputFasta(rows, lineageGermlineMethod, { references, templateOrdinal: lineageGermline?.selectedOrdinal });
        const records = parseFasta(input.fasta, true);
        const aminoAcids = records.map((record, index) => translateCodonSequence(record.sequence, input.frames[index] ?? 0));
        const alignedAminoAcids = await runAlivibeMsaTask(aminoAcids, signal, 3, "amino-acid");
        const projected = records.map((record, index) => projectCodonAlignment(
          record.sequence,
          alignedAminoAcids[index],
          input.frames[index] ?? 0,
        ));
        const maximum = Math.max(...projected.map((sequence) => sequence.length));
        const fasta = records.map((record, index) => `>${record.name}\n${projected[index].padEnd(maximum, "-")}`).join("\n") + "\n";
        return { fasta, frameOffset: 0 as AlignmentFrameOffset };
      });
      installAlignment(
        next.fasta,
        "Alivibe-compatible AA POA MSA · X wildcard · codon back-translation · WASM · three refinement passes",
        false,
        selectedLineageIds,
        next.frameOffset,
      );
      setAlignmentEditorStatus("");
      setAlignmentEditorError("");
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function inferTree() {
    if (!alignment) return;
    setBusy("Running FastTree WASM and rooting on the N-masked germline");
    setTreeError("");
    const revision = alignmentRevisionRef.current;
    const alignmentSnapshot = alignment;
    const sourceSnapshot = alignmentSource;
    try {
      const execution = await runInActiveLock((signal) => runFastTreeTask(alignmentSnapshot, treeModel, treeFast, signal));
      if (revision !== alignmentRevisionRef.current) throw new Error("The alignment changed while FastTree was running. Run the tree again on the current alignment.");
      const rooted = rootOnOutgroup(parseNewick(execution.newick), GERMLINE_OUTGROUP);
      const canonicalRooted = canonicalizeTree(rooted);
      const stable = collapseShortInternalBranches(canonicalRooted);
      setTreeRun({
        ...execution,
        rootedNewick: `${serializeNewick(canonicalRooted)};`,
        stableNewick: `${serializeNewick(stable.root)};`,
        collapsedEdges: stable.collapsedEdges,
        collapseThreshold: stable.threshold,
        source: sourceSnapshot,
        frameOffset: alignmentFrameOffset,
        origin: "fasttree",
        observedOnly: false,
      });
      setTreeViewMode("rooted");
      window.requestAnimationFrame(() => treeResultRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
    } catch (operationError) {
      if (!isAbortError(operationError)) {
        setTreeError(operationError instanceof Error ? operationError.message : String(operationError));
        window.requestAnimationFrame(() => treeResultRef.current?.scrollIntoView({ behavior: "smooth", block: "center" }));
      }
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function acceptUploadedTree(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file || !alignment) return;
    setTreeError("");
    const revision = alignmentRevisionRef.current;
    try {
      const text = await file.text();
      if (revision !== alignmentRevisionRef.current) throw new Error("The alignment changed while the tree file was being read. Choose the tree again for the current alignment.");
      const inspected = inspectAlignment(alignment);
      const observedNames = inspected.records.filter((record) => record.name !== GERMLINE_OUTGROUP).map((record) => record.name);
      const uploaded = prepareUploadedLineageTree(text, observedNames, true);
      setTreeRun({
        newick: uploaded.rawNewick,
        rootedNewick: uploaded.canonicalNewick,
        stableNewick: uploaded.stableNewick,
        inputFasta: "",
        alignmentFasta: inspected.fasta,
        command: "Uploaded Newick · FastTree was not run",
        rows: observedNames.length,
        columns: inspected.columns,
        fingerprint: inspected.fingerprint,
        collapsedEdges: uploaded.collapsedEdges,
        collapseThreshold: uploaded.collapseThreshold,
        source: `Uploaded tree · ${file.name}`,
        frameOffset: alignmentFrameOffset,
        origin: "uploaded",
        observedOnly: uploaded.observedOnly,
      });
      setPhyloUcaState(null);
      setTreeViewMode("rooted");
      window.requestAnimationFrame(() => treeResultRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
    } catch (uploadError) {
      setTreeError(uploadError instanceof Error ? uploadError.message : String(uploadError));
      window.requestAnimationFrame(() => treeResultRef.current?.scrollIntoView({ behavior: "smooth", block: "center" }));
    }
  }

  function importEditedAlignment(
    text: string,
    source = "Alivibe-corrected alignment",
    transferredFrameOffset?: AlignmentFrameOffset,
    options: { baseline?: string; lineageIds?: number[]; requireExactSerialization?: boolean } = {},
  ) {
    const baseline = options.baseline ?? alignment;
    const lineageIds = options.lineageIds ?? selectedLineageIds;
    if (!baseline) throw new Error("Create a lineage alignment before importing a correction.");
    const inspected = validateCorrectedAlignment(baseline, text);
    if (options.requireExactSerialization && inspected.fasta !== text) {
      throw new Error("Swig would have to normalize the nucleotide FASTA returned by Alivibe. The return was refused instead of silently changing it.");
    }
    // A corrected alignment is derived from the current MSA. Preserve its
    // explicit biological phase unless the versioned Alivibe bridge returns a
    // user-selected replacement. Gap placement is evidence about alignment
    // quality, not permission to redefine the coding frame.
    const frameOffset = transferredFrameOffset ?? alignmentFrameOffset;
    installAlignment(inspected.fasta, source, true, lineageIds, frameOffset);
    setAlignmentEditorError("");
    setAlignmentEditorStatus(`Accepted corrected nucleotide alignment: ${inspected.rows.toLocaleString()} rows × ${inspected.columns.toLocaleString()} columns${inspected.removedRows.length?` · ${inspected.removedRows.length.toLocaleString()} biological row${inspected.removedRows.length===1?"":"s"} deleted`:""}${inspected.removedNucleotides?` · ${inspected.removedNucleotides.toLocaleString()} nucleotide character${inspected.removedNucleotides===1?"":"s"} deleted`:""} · AA reading frame starts at nucleotide column ${frameOffset + 1} · fingerprint ${inspected.fingerprint}. FastTree will use these exact retained rows and columns.`);
  }

  function returnFromAlivibe(session: AlivibeRoundTripSession) {
    if (alivibeSessionRef.current?.token !== session.token) {
      throw new Error("This Alivibe editor belongs to an expired Swig round trip. Close it and reopen the current alignment.");
    }
    if (session.popup.closed) throw new Error("The Alivibe editor has been closed.");
    const currentAlignment = alignmentRef.current;
    if (!currentAlignment) throw new Error("The originating Swig alignment is no longer loaded.");
    assertAlivibeRoundTripTarget(
      { groupKey: session.groupKey, alignmentFingerprint: session.baselineFingerprint },
      { groupKey: selectedGroupKeyRef.current, alignmentFingerprint: inspectAlignment(currentAlignment).fingerprint },
    );
    const returned = readAlivibeNucleotideFasta(session.popup);
    importEditedAlignment(returned.fasta, `Alivibe-corrected alignment · ${returned.sourceRevision.slice(0, 12)}`, returned.frameOffset, {
      baseline: session.baseline,
      lineageIds: session.lineageIds,
      requireExactSerialization: true,
    });
    alivibeSessionRef.current = null;
    window.focus();
    session.popup.close();
  }

  function openAlivibeEditor() {
    if (!alignment) return;
    const existing = alivibeSessionRef.current;
    if (existing && !existing.popup.closed) {
      existing.popup.focus();
      setAlignmentEditorError("An Alivibe round trip is already open. Return or close that editor before opening another alignment.");
      return;
    }
    setAlignmentEditorError("");
    setAlignmentEditorStatus("Opening the bundled Alivibe editor…");
    const base = import.meta.env.BASE_URL.endsWith("/") ? import.meta.env.BASE_URL : `${import.meta.env.BASE_URL}/`;
    const editorUrl = new URL(`${base}tools/alivibe.html`, window.location.origin);
    editorUrl.searchParams.set("swigBridge", String(ALIVIBE_BRIDGE_VERSION));
    editorUrl.searchParams.set("source", ALIVIBE_SOURCE_REVISION.slice(0, 12));
    const token = globalThis.crypto?.randomUUID?.() ?? `${Date.now()}-${Math.random().toString(36).slice(2)}`;
    const popup = window.open(editorUrl.href, `swig-alivibe-${token}`, "popup,width=1500,height=920") as AlivibeEditorWindow | null;
    if (!popup) {
      setAlignmentEditorError("The browser blocked the Alivibe window. Allow pop-ups for this page and try again.");
      return;
    }
    const baseline = alignment;
    const baselineFingerprint = inspectAlignment(baseline).fingerprint;
    const lineageIds = [...new Set(selectedLineageIds)].sort((left, right) => left - right);
    const groupKey = selectedGroupKey;
    const session: AlivibeRoundTripSession = { token, popup, baseline, baselineFingerprint, lineageIds, groupKey, frameOffset: alignmentFrameOffset };
    alivibeSessionRef.current = session;
    let attempts = 0;
    const timer = window.setInterval(() => {
      attempts += 1;
      if (popup.closed) {
        window.clearInterval(timer);
        if (alivibeSessionRef.current?.token === token) alivibeSessionRef.current = null;
        setAlignmentEditorStatus("Alivibe closed without a direct return. Its downloaded nucleotide FASTA can still be loaded with Import corrected FASTA.");
        return;
      }
      try {
        const bridge = getAlivibeBridge(popup);
        if (!bridge) {
          if (attempts < 240) return;
          throw new Error("The bundled Alivibe bridge did not become ready. Close the editor and try again.");
        }
        bridge.installMsaRunner((sequences, scoringMode) => createAlivibeMsaJob(sequences, 3, scoringMode));
        const loaded = loadAlivibeNucleotideFasta(popup, session.baseline, session.frameOffset);
        assertAlivibeInitialLoad(session.baseline, loaded);
        const controls = popup.document.getElementById("controls");
        if (!controls) return;
        if (!popup.document.getElementById("swig-return-control")) {
          const group = popup.document.createElement("div");
          group.id = "swig-return-control";
          group.className = "control-group";
          const label = popup.document.createElement("label");
          label.textContent = "Swig round trip";
          const button = popup.document.createElement("button");
          button.type = "button";
          button.textContent = "Return alignment to Swig";
          button.className = "active";
          button.onclick = () => {
            try {
              returnFromAlivibe(session);
            } catch (importError) {
              setAlignmentEditorError(importError instanceof Error ? importError.message : String(importError));
              window.focus();
            }
          };
          group.append(label, button);
          controls.prepend(group);
        }
        setAlignmentEditorStatus(`Exact nucleotide alignment loaded in bundled Alivibe ${ALIVIBE_SOURCE_REVISION.slice(0, 12)} with reading frame ${session.frameOffset + 1}. Edit it, then press Return alignment to Swig in Alivibe’s toolbar.`);
        window.clearInterval(timer);
      } catch (openError) {
        setAlignmentEditorError(openError instanceof Error ? openError.message : String(openError));
        setAlignmentEditorStatus("Alivibe round trip stopped before editing. No Swig alignment was changed.");
        if (alivibeSessionRef.current?.token === token) alivibeSessionRef.current = null;
        window.clearInterval(timer);
      }
    }, 250);
  }

  function importFromAlivibe() {
    setAlignmentEditorError("");
    try {
      const session = alivibeSessionRef.current;
      if (!session || session.popup.closed) throw new Error("No live Alivibe round trip is available. Open the current alignment in Alivibe, or load an exported nucleotide FASTA file.");
      returnFromAlivibe(session);
    } catch (importError) {
      setAlignmentEditorError(importError instanceof Error ? importError.message : String(importError));
    }
  }

  async function acceptEditedAlignment(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file) return;
    try {
      importEditedAlignment(await file.text(), `Corrected alignment imported from ${file.name}`);
    } catch (importError) {
      setAlignmentEditorError(importError instanceof Error ? importError.message : String(importError));
    }
  }

  async function acceptMsa(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file) return;
    try {
      const text = await file.text();
      prepareReferenceMsa(text);
      setUploadedMsa(text);
      setUploadedMsaName(file.name);
      setChmmSource("upload");
      setPreparedMsa("");
      setChmm(null);
      setChmmRun(null);
      setChimeraDetail(null);
      setError("");
    } catch (operationError) {
      setError(operationError instanceof Error ? operationError.message : String(operationError));
    }
  }

  async function runChmmAnalysis() {
    setBusy(chmmSource === "selected" ? `Building ${chmmSegment} reference MSA with Kalign WASM` : "Validating loaded reference MSA");
    setError("");
    try {
      await runInActiveLock(async (signal) => {
        const referenceSource = chmmSource === "upload" ? uploadedMsa : references[chmmSegment];
        let msa = chmmSource === "upload" ? referenceSource : preparedMsa;
        if (!msa) msa = await runKalignTask(referenceSource, signal);
        prepareReferenceMsa(msa);
        setBusy(`Running CHMMAIRRa ${chmmSegment} model with ${Math.min(workers, 16)} workers`);
        const rates = mutationRates.split(/[\s,;]+/).map(Number).filter((value) => Number.isFinite(value) && value >= 0 && value < 1);
        if (chmmMethod === "DB" && !rates.length) throw new Error("Provide at least one mutation-rate state for the discretized Bayesian model.");
        const options: ChmmRunOptions = {
          segment: chmmSegment,
          method: chmmMethod,
          priorProbability: chmmPrior,
          baseMutationProbability: 0.05,
          mutationRates: rates,
          mutationSwitchProbability: 0,
          detailed: chmmDetailed,
          minDfr: chmmMinDfr,
          threshold: chmmThreshold,
          workers,
        };
        // The worker owns the canonical cumulative filter. Reading it here avoids
        // launching an HMM scan from a stale React render after a filter commit.
        const inputMask = await runtime.activeMask();
        const result = await runChmmairra(store, msa, options, inputMask ?? undefined, (processed, total) => setProgress({ processed, total }), signal);
        setPreparedMsa(msa);
        setChmm(result);
        setChmmRun({ msa, options, inputMask });
        setChmmFilterThreshold(chmmThreshold);
        setChimeraDetail(null);
      });
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function applyChimeraFilter() {
    if (!chmm || !chmmRun) return;
    const next = new Uint8Array(store.count);
    let retained = 0;
    let unevaluatedRetained = 0;
    for (let ordinal = 0; ordinal < store.count; ordinal += 1) {
      if (chmmRun.inputMask && !chmmRun.inputMask[ordinal]) continue;
      const probability = chmm.probabilities[ordinal];
      const keep = Number.isFinite(probability) ? probability < chmmFilterThreshold : retainUnevaluated;
      if (!keep) continue;
      next[ordinal] = 1;
      retained += 1;
      if (!Number.isFinite(probability)) unevaluatedRetained += 1;
    }
    const result = await operation("Applying CHMMAIRRa exclusion to the downstream working set", () => runtime.setActiveMask(next));
    if (!result) return;
    setWorkingMask(next);
    setWorkingStages((stages) => [...stages.filter((stage) => stage.id !== "chimera"), {
      id: "chimera",
      label: `${chmm.segment} chimera posterior < ${chmmFilterThreshold}`,
      input: chmm.inputRecords,
      retained,
      discarded: chmm.inputRecords - retained,
      detail: retainUnevaluated ? `${unevaluatedRetained.toLocaleString()} unevaluated records retained.` : "Unevaluated records excluded.",
    }]);
    invalidateAssignmentDependentAnalyses();
    setTreeError("");
    advanceModule("chimera","selection");
  }

  async function openChimera(ordinal: number) {
    if (!chmmRun) {
      setError("Run CHMMAIRRa before opening a Viterbi parent reconstruction.");
      return;
    }
    setBusy("Running the selected record's detailed CHMMAIRRa Viterbi path");
    setError("");
    try {
      const detail = await runInActiveLock((signal) => runChmmairraDetail(store, chmmRun.msa, chmmRun.options, ordinal, signal));
      setChimeraDetail(detail);
      window.requestAnimationFrame(() => chimeraDetailRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }));
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function downloadChmm() {
    if (!chmm) return;
    setBusy("Writing CHMMAIRRa result table");
    try {
      await runInActiveLock(async(signal)=>{const extension=tableExtension(exportFormat);await saveStream(`${baseName(inputName)}.chmmairra-${chmm.segment.toLowerCase()}${extension}`, "CHMMAIRRa result table", extension, (writer) => writeChmmairra(store, chmm, exportFormat, writer, signal),exportCompression);});
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function runQuery() {
    const queries = parseQueries(queryText);
    if (!queries.length) {
      setError("Paste at least one nucleotide/amino-acid sequence or FASTA record.");
      return;
    }
    if(queryResultMode==="lineages"&&!lineages){setError("Lineage result mode requires lineage assignments on the current working set. Assign lineages first, then run this search.");return;}
    const metric = queryTarget === "trimmed" ? "sketch" : queryMetric;
    const response = await operation(queryConstraintMode === "infer" ? "Assigning query V/J calls with the Swig aligner, then searching" : "Searching the assigned repertoire", async (signal) => {
      let searchQueries = queries;
      let queryConstraints: QueryConstraint[] | undefined;
      let inferredAssignments: InferredQueryAssignment[] = [];
      if (queryConstraintMode === "infer") {
        const inferred = await inferQueryAssignments(queries, queryTarget, references, callingProfile, minimumIdentity, strand, workers, signal);
        inferredAssignments = inferred;
        const usable = inferred.filter((assignment) => assignment.searchSequence && (queryV || assignment.vCall) && (queryJ || assignment.jCall));
        if (!usable.length) {
          throw new Error("The aligner could not infer both V and J for any query. Use complete rearranged nucleotide sequences, lower the main identity floor in a new run, or provide manual V/J overrides.");
        }
        searchQueries = usable.map((assignment) => assignment.searchSequence);
        queryConstraints = usable.map((assignment) => ({
          locus: queryLocus || assignment.locus || undefined,
          vCall: queryV || assignment.vCall || undefined,
          jCall: queryJ || assignment.jCall || undefined,
        }));
      }
      const hits = await runtime.query(searchQueries, {
        target: queryTarget,
        metric,
        identity: queryIdentity,
        maxResults: queryLimit,
        locus: queryConstraintMode === "manual" ? queryLocus || undefined : undefined,
        vCall: queryConstraintMode === "manual" ? queryV || undefined : undefined,
        jCall: queryConstraintMode === "manual" ? queryJ || undefined : undefined,
        queryConstraints,
        callResolution: resolution,
        ambiguity,
        productiveOnly,
        resultMode: queryResultMode,
      });
      return { hits, inferredAssignments };
    });
    if (!response) return;
    const { hits, inferredAssignments } = response;
    setQueryInference(inferredAssignments);
    setQueryHits(hits);
    setExpanded(null);
    if(queryResultMode==="lineages"){
      setQueryLineageHits(hits.flatMap((hit)=>hit.lineageId?[{lineageId:hit.lineageId,bestScore:hit.score,bestDistance:hit.distance,bestOrdinal:hit.ordinal,matchedSequences:hit.matchedSequences??1,matchedQueries:hit.matchedQueries??1}]:[]));
    }else setQueryLineageHits([]);
    const indexRows = await store.indexRecords([...new Set(hits.map((hit) => hit.ordinal))]);
    setQueryRecords(new Map(indexRows.map((record) => [record.ordinal, record])));
  }

  async function openQueryLineage(hit:QueryLineageHit){
    if(!lineages)return;
    let summary=lineages.summaries.find((item)=>item.id===hit.lineageId);
    if(!summary){
      const members=await runtime.lineageMembers(hit.lineageId,0,1);
      const [record]=await store.indexRecords([hit.bestOrdinal,...members.ordinals]);
      if(!record){setError(`Could not load lineage ${hit.lineageId}.`);return;}
      summary={id:hit.lineageId,representativeOrdinal:record.ordinal,uniqueMembers:members.total,abundance:members.total,locus:record.locus,vCalls:record.vCall.split(",").map((value)=>value.trim()).filter(Boolean),jCalls:record.jCall.split(",").map((value)=>value.trim()).filter(Boolean),cdr3Length:record.cdr3.length,studyScope:lineageScope,studyGroup:datasetScopeValue(record,lineageScope),sampleIds:record.sampleId?[record.sampleId]:[],sampleCounts:record.sampleId?[{sampleId:record.sampleId,uniqueMembers:members.total,abundance:members.total}]:[],subjectIds:record.subjectId?[record.subjectId]:[],timepoints:record.timepoint?[record.timepoint]:[],compartments:record.compartment?[record.compartment]:[]};
    }
    await openLineage(summary);
  }

  async function previewQueryInference() {
    const queries = parseQueries(queryText);
    if (!queries.length) {
      setError("Paste at least one complete rearranged nucleotide sequence or FASTA record.");
      return;
    }
    setBusy("Assigning seed V/J calls with the main aligner configuration");
    setError("");
    try {
      setQueryInference(await runInActiveLock((signal) => inferQueryAssignments(queries, queryTarget, references, callingProfile, minimumIdentity, strand, workers, signal)));
    } catch (operationError) {
      if (!isAbortError(operationError)) setError(operationError instanceof Error ? operationError.message : String(operationError));
    } finally {
      if (!cancellationRecoveryRef.current) setBusy("");
    }
  }

  async function expandMatches() {
    const seeds = expanded?.ordinals.length ? expanded.ordinals : [...new Set(queryHits.map((hit) => hit.ordinal))];
    if (!seeds.length) {
      setError("Run a sequence query before expanding its matched set.");
      return;
    }
    const result = await operation("Expanding a single-linkage CDR3 neighborhood", () => runtime.expand(seeds, {
      identity: queryIdentity,
      callResolution: resolution,
      ambiguity,
      productiveOnly,
      requireSameLocus: true,
      maxCandidateComparisons: candidateCap,
      maxResults: Math.max(queryLimit, 10_000),
      scope: lineageScope,
    }));
    if (!result) return;
    setExpanded(result);
    const indexRows = await store.indexRecords(result.ordinals.slice(0, 500));
    setQueryRecords(new Map(indexRows.map((record) => [record.ordinal, record])));
  }

  const workingCount = workingStages.length ? workingStages[workingStages.length - 1].retained : store.count;
  const guidedAction = !skippedModules.has("dedup")&&!dedup ? "run-dedup"
    : !skippedModules.has("dedup")&&!workingStages.some((stage)=>stage.id==="dedup") ? "apply-dedup"
    : !skippedModules.has("chimera")&&!chmm ? "run-chimera"
    : !skippedModules.has("chimera")&&!workingStages.some((stage)=>stage.id==="chimera") ? "apply-chimera"
    : !skippedModules.has("selection")&&!selectionPreview ? "preview-selection"
    : !skippedModules.has("selection")&&selectionPreview&&!selectionApplied ? "apply-selection"
    : !skippedModules.has("lineage")&&!lineages ? "run-lineages"
    : !skippedModules.has("diagnostics")&&lineages&&!shmDashboard ? "run-shm"
    : "complete";
  const guidedClass=(action:string,base="post-primary")=>`${base}${guidedAction===action?" guided-next":""}`;
  useEffect(()=>{
    const module:PostModuleId|undefined=guidedAction==="run-dedup"||guidedAction==="apply-dedup"?"dedup":guidedAction==="run-chimera"||guidedAction==="apply-chimera"?"chimera":guidedAction==="preview-selection"||guidedAction==="apply-selection"?"selection":guidedAction==="run-lineages"?"lineage":undefined;
    if(module)openModule(module);
  },[guidedAction]);
  const chmmFilterThresholdValid = Number.isFinite(chmmFilterThreshold) && chmmFilterThreshold >= 0 && chmmFilterThreshold <= 1;
  const chmmFilterPreview = useMemo(() => {
    if (!chmm || !chmmRun || !Number.isFinite(chmmFilterThreshold) || chmmFilterThreshold < 0 || chmmFilterThreshold > 1) return null;
    let retained = 0;
    let excluded = 0;
    let unevaluated = 0;
    for (let ordinal = 0; ordinal < store.count; ordinal += 1) {
      if (chmmRun.inputMask && !chmmRun.inputMask[ordinal]) continue;
      const probability = chmm.probabilities[ordinal];
      if (!Number.isFinite(probability)) {
        unevaluated += 1;
        if (retainUnevaluated) retained += 1;
        else excluded += 1;
      } else if (probability >= chmmFilterThreshold) excluded += 1;
      else retained += 1;
    }
    return { retained, excluded, unevaluated };
  }, [chmm, chmmFilterThreshold, chmmRun, retainUnevaluated, store.count]);
  const displayedQueryOrdinals = expanded?.ordinals ?? queryHits.map((hit) => hit.ordinal);
  const displayedQueryRows = [...new Set(displayedQueryOrdinals)].slice(0, 500).flatMap((ordinal) => queryRecords.get(ordinal) ?? []);
  const lineageShmById=useMemo(()=>{
    const map=new Map<number,{mean:number;maximum:number;p95:number}>();
    for(const group of shmDashboard?.lineages??[]){const match=/^Lineage\s+(\d+)$/.exec(group.label);if(match)map.set(Number(match[1]),{mean:group.mean,maximum:group.maximum??0,p95:group.p95??0});}
    return map;
  },[shmDashboard]);
  const lineageLowestCdr3ById=useMemo(()=>new Map((shmDashboard?.lowestByLineage??[]).map((value)=>[value.lineageId,value])),[shmDashboard]);
  const lineageRepresentativeCdr3=(summary:LineageSummary)=>{
    const lowest=lineageLowestCdr3ById.get(summary.id);
    return {nt:lowest?.cdr3Nt||summary.representativeCdr3Nt||"",aa:lowest?.cdr3Aa||summary.representativeCdr3Aa||"",lowestShm:Boolean(lowest)};
  };
  const sortedLineageSummaries=useMemo(()=>{
    const selected=lineageSelectedSamples;
    const search=lineageSearch.trim().toLowerCase();
    const summaries=(lineages?.summaries??[]).filter((summary)=>{
      const samples=summary.sampleIds??[];
      const lowestCdr3=lineageLowestCdr3ById.get(summary.id);
      if(search&&!`${summary.id} ${summary.studyGroup} ${summary.locus} ${summary.vCalls.join(" ")} ${summary.jCalls.join(" ")} ${summary.representativeCdr3Nt??""} ${summary.representativeCdr3Aa??""} ${lowestCdr3?.cdr3Nt??""} ${lowestCdr3?.cdr3Aa??""}`.toLowerCase().includes(search))return false;
      if(lineageLocusFilter&&summary.locus!==lineageLocusFilter)return false;
      if(lineageMinAbundance&&summary.abundance<lineageMinAbundance)return false;
      if(lineageMinUnique&&summary.uniqueMembers<lineageMinUnique)return false;
      if(lineageMinSamples&&samples.length<lineageMinSamples)return false;
      if(lineageMinCdr3Length&&summary.cdr3Length<lineageMinCdr3Length)return false;
      if(lineageMaxCdr3Length&&summary.cdr3Length>lineageMaxCdr3Length)return false;
      if(!lineageDoubleDMatches(summary,lineageDoubleDFilter))return false;
      if(lineageSpecificSample&&!samples.includes(lineageSpecificSample))return false;
      if(lineageSampleFilterMode==="multiple"&&samples.length<2)return false;
      if(lineageSampleFilterMode==="single"&&samples.length!==1)return false;
      if(lineageSampleFilterMode==="selected-any"&&(!selected.size||![...selected].some((sample)=>samples.includes(sample))))return false;
      if(lineageSampleFilterMode==="selected-all"&&(!selected.size||![...selected].every((sample)=>samples.includes(sample))))return false;
      if(lineageSampleFilterMode==="selected-only"&&(!selected.size||!samples.length||samples.some((sample)=>!selected.has(sample))))return false;
      if(lineageMinimumShm){const statistic=lineageShmById.get(summary.id)?.[lineageShmFilterStatistic];if(statistic===undefined||statistic<lineageMinimumShm)return false;}
      return true;
    });
    const text=(summary:LineageSummary,key:LineageSortKey)=>key==="vCalls"?summary.vCalls.join(", "):key==="jCalls"?summary.jCalls.join(", "):key==="studyGroup"?summary.studyGroup:key==="locus"?summary.locus:"";
    const numeric=(summary:LineageSummary,key:LineageSortKey)=>key==="sampleCount"?(summary.sampleIds??[]).length:key==="doubleDPositiveMembers"?(summary.doubleDPositiveMembers??0):key==="shmMean"?(lineageShmById.get(summary.id)?.mean??Number.NEGATIVE_INFINITY):key==="shmMaximum"?(lineageShmById.get(summary.id)?.maximum??Number.NEGATIVE_INFINITY):key==="shmP95"?(lineageShmById.get(summary.id)?.p95??Number.NEGATIVE_INFINITY):Number(summary[key as "id"|"abundance"|"uniqueMembers"|"cdr3Length"]);
    const numericKeys=new Set<LineageSortKey>(["id","abundance","uniqueMembers","doubleDPositiveMembers","sampleCount","cdr3Length","shmMean","shmMaximum","shmP95"]);
    summaries.sort((left,right)=>{const comparison=numericKeys.has(lineageSort.key)?numeric(left,lineageSort.key)-numeric(right,lineageSort.key):text(left,lineageSort.key).localeCompare(text(right,lineageSort.key),undefined,{numeric:true,sensitivity:"base"});return (lineageSort.direction==="asc"?comparison:-comparison)||right.abundance-left.abundance||left.id-right.id;});
    return summaries;
  },[lineages,lineageSort,lineageSearch,lineageLocusFilter,lineageMinAbundance,lineageMinUnique,lineageMinSamples,lineageMinCdr3Length,lineageMaxCdr3Length,lineageDoubleDFilter,lineageSpecificSample,lineageSampleFilterMode,lineageSelectedSamples,lineageMinimumShm,lineageShmFilterStatistic,lineageShmById,lineageLowestCdr3ById]);
  const sortLineages=(key:LineageSortKey)=>setLineageSort((current)=>current.key===key?{key,direction:current.direction==="desc"?"asc":"desc"}:{key,direction:"desc"});
  const sortIndicator=(key:LineageSortKey)=>lineageSort.key===key?(lineageSort.direction==="desc"?" ↓":" ↑"):"";
  const lineageSampleOptions=useMemo(()=>[...new Set([...datasets.map((dataset)=>dataset.sampleId),...(lineages?.summaries??[]).flatMap((summary)=>summary.sampleIds??[])].filter(Boolean))].sort((left,right)=>left.localeCompare(right,undefined,{numeric:true,sensitivity:"base"})),[datasets,lineages]);
  const visibleLineageSampleOptions=useMemo(()=>{
    const query=lineageSampleSearch.trim().toLowerCase();
    if(query)return lineageSampleOptions.filter((sample)=>sample.toLowerCase().includes(query)).slice(0,200);
    const selected=lineageSampleOptions.filter((sample)=>lineageSelectedSamples.has(sample));
    return [...selected,...lineageSampleOptions.filter((sample)=>!lineageSelectedSamples.has(sample)).slice(0,Math.max(0,40-selected.length))];
  },[lineageSampleOptions,lineageSampleSearch,lineageSelectedSamples]);
  const lineageLocusOptions=useMemo(()=>[...new Set((lineages?.summaries??[]).map((summary)=>summary.locus).filter(Boolean))].sort(),[lineages]);
  const lineageShmLabel=(lineageId:number,key:LineageShmFilterStatistic)=>{const value=lineageShmById.get(lineageId)?.[key];if(value===undefined)return "—";return shmMetric.toLowerCase().includes("rate")?`${(value*100).toFixed(1)}%`:value.toFixed(value<10?2:0);};
  const vChart = lineages?.vUsage.slice(0, topGenes).map((item) => ({ label: item.call, value: item[geneMetric] })) ?? [];
  const jChart = lineages?.jUsage.slice(0, topGenes).map((item) => ({ label: item.call, value: item[geneMetric] })) ?? [];
  const railClass=(module:PostModuleId)=>`${activeWorkspace===module?"active":""}${skippedModules.has(module)?" skipped":""}`.trim();

  return <section className={`post-analysis-shell${directLineage ? " direct-lineage-shell" : ""}`}>
    <header className="post-analysis-heading"><div><span className="section-kicker">{directLineage ? "Single lineage analysis" : "Post-assignment analyses"}</span><h2>{directLineage ? "Alignment, tree and UCA inference" : "Repertoire structure and lineage analysis"}</h2><p>{directLineage ? "Every sequence in the uploaded input is treated as a member of one lineage. No clustering, repertoire filtering, or lineage selection is performed." : "Exact collapse or denoising and chimera exclusion modify an explicit cumulative working set. CHMMAIRRa, lineage assignment, repertoire querying, and expansion consume that set; alignment and tree inference consume the selected lineage."}</p><a href="./METHODS_INDEX.md" target="_blank" rel="noreferrer">Complete methods index ↗</a></div><div className="local-method-note"><span>Data handling</span><strong>Browser-local</strong><small>Input, germlines, and results are not submitted to an analysis server.</small></div></header>

    <div className="post-context-workspace contextual-workspace">
      {!directLineage && <nav className="context-rail post-context-rail" aria-label="Post-analysis sections">
        <div className="context-rail-heading"><span>Post-analysis</span><small>{workingCount.toLocaleString()} active records</small></div>
        <button type="button" className={activeWorkspace==="overview"?"active":""} onClick={()=>setActiveWorkspace("overview")}><b>00</b><span>Overview<small>Working set + exports</small></span></button>
        <button type="button" className={railClass("alleles")} onClick={()=>setActiveWorkspace("alleles")}><b>01</b><span>Allele pooling<small>{skippedModules.has("alleles")?"Skipped":alleleRefinement?`${alleleRefinement.activeRecords.toLocaleString()} modeled`:"Optional · runs first"}</small></span></button>
        <button type="button" className={railClass("dedup")} onClick={()=>setActiveWorkspace("dedup")}><b>02</b><span>Collapse<small>{skippedModules.has("dedup")?"Skipped":dedup?`${dedup.uniqueRecords.toLocaleString()} representatives`:"Not run"}</small></span></button>
        <button type="button" className={railClass("chimera")} onClick={()=>setActiveWorkspace("chimera")}><b>03</b><span>Chimera<small>{skippedModules.has("chimera")?"Skipped":chmm?`${chmm.evaluated.toLocaleString()} evaluated`:"Optional"}</small></span></button>
        <button type="button" className={railClass("selection")} onClick={()=>setActiveWorkspace("selection")}><b>04</b><span>Selection<small>{skippedModules.has("selection")?"Skipped by default":selectionApplied?"Committed":selectionPreview?"Preview ready":"Configure filters"}</small></span></button>
        <button type="button" className={railClass("lineage")} onClick={()=>setActiveWorkspace("lineage")}><b>05</b><span>Lineages<small>{skippedModules.has("lineage")?"Skipped":lineages?`${lineages.summaries.length.toLocaleString()} assigned`:"Not run"}</small></span></button>
        <button type="button" className={railClass("diagnostics")} onClick={()=>setActiveWorkspace("diagnostics")}><b>06</b><span>Post-lineage<small>{skippedModules.has("diagnostics")?"Skipped":"SHM + allele evidence"}</small></span></button>
        <button type="button" disabled={!selectedLineage} className={railClass("workbench")} onClick={()=>setActiveWorkspace("workbench")}><b>07</b><span>Workbench<small>{skippedModules.has("workbench")?"Skipped":selectedLineage?`Lineage ${selectedLineage.id}`:"Select a lineage"}</small></span></button>
        <button type="button" className={railClass("query")} onClick={()=>setActiveWorkspace("query")}><b>08</b><span>Query<small>{skippedModules.has("query")?"Skipped":"Search + expand"}</small></span></button>
        {sidebarTools}
      </nav>}

      <div className="context-main post-context-main">
    {activeWorkspace==="overview"&&<section className="post-overview-panel">
    <div className="post-overview-intro"><strong>Post-analysis workspace</strong><span>Each stage consumes the retained set from the preceding stage. Computation and intermediate state remain in this browser.</span></div>

    {autoPipeline?.enabled&&<section className={`pipeline-run-banner ${busy?"running":"complete"}`}><div><span className="section-kicker">Pipeline execution</span><h3>{busy?busy:"Selected repertoire-scale stages complete"}</h3><p>{pipelineReport.length?pipelineReport.join(" "):"The configured stages are running in order; each stage receives the retained set from the previous stage."}</p></div><strong>{busy?"RUNNING":"COMPLETE"}</strong></section>}

    <div className="post-method-map"><article><b>01</b><span>Alleles</span><strong>Optional donor-level evidence pooling</strong></article><article><b>02</b><span>Collapse</span><strong>Exact deduplication or denoising</strong></article><article><b>03</b><span>QC</span><strong>Optional CHMMAIRRa</strong></article><article><b>04</b><span>Select</span><strong>Skipped by default; opt in to filtering</strong></article><article><b>05</b><span>Repertoire</span><strong>Assign lineages</strong></article><article><b>06</b><span>Post-lineage</span><strong>SHM + allele evidence</strong></article><article><b>07</b><span>On demand</span><strong>Align + infer tree</strong></article><article><b>08</b><span>Targeted</span><strong>Query + expand</strong></article></div>

    {datasets.length>0&&<section className="study-policy-panel"><header><div><span className="section-kicker">Multi-dataset policy</span><h3>{datasets.length.toLocaleString()} dataset{datasets.length===1?"":"s"} · {new Set(datasets.map((item)=>item.sampleId)).size.toLocaleString()} sample{new Set(datasets.map((item)=>item.sampleId)).size===1?"":"s"} · {new Set(datasets.map((item)=>item.subjectId)).size.toLocaleString()} donor{new Set(datasets.map((item)=>item.subjectId)).size===1?"":"s"}</h3></div><p>Metadata are stored on every AIRR row. Scope controls below are candidate-generation boundaries, not display-only labels.</p></header><div>{datasets.slice(0,12).map((item)=><article key={item.datasetId}><strong>{item.sampleId}</strong><span>{item.subjectId}</span><small>{[item.cohort,item.timepoint,item.compartment,item.inputName].filter(Boolean).join(" · ")}</small></article>)}</div>{datasets.length>12&&<small>+ {(datasets.length-12).toLocaleString()} additional datasets</small>}</section>}

    <section className="working-set-panel" aria-label="Downstream working-set pipeline">
      <header><div><span className="section-kicker">Cumulative downstream filter</span><h3>{workingCount.toLocaleString()} of {store.count.toLocaleString()} records active</h3><p>Nothing is deleted from the AIRR result. Applying a stage excludes rows from later computation; resetting restores the complete input.</p></div><button type="button" disabled={Boolean(busy) || !workingStages.length} onClick={() => void resetWorkingSet()}>Reset to all records</button></header>
      <div className="working-set-flow"><article className="source"><span>Assigned input</span><strong>{store.count.toLocaleString()}</strong><small>records</small></article>{workingStages.map((stage,index) => <article key={`${stage.id}-${index}`} className={stage.id}><span>{stage.label}</span><strong>{stage.retained.toLocaleString()}</strong><small>retained · {stage.discarded.toLocaleString()} excluded at this step</small><p>{stage.detail}</p></article>)}{!workingStages.length && <article className="pass-through"><span>No applied filter</span><strong>All records</strong><small>pass downstream</small></article>}<article className="consumers"><span>Current consumers</span><strong>lineages · SHM · allele hints · query</strong><small>alignment/tree follow the selected lineage</small></article></div>
    </section>

    <section className="post-export-center"><div><span className="section-kicker">Export center</span><h3>Machine-readable analysis outputs</h3><p>The selected table is streamed from browser storage. Gzip uses a bounded native compression stream when the browser supports it.</p></div><label><span>Tabular format</span><select value={exportFormat} onChange={(event)=>setExportFormat(event.target.value as TableExportFormat)}><option value="tsv">AIRR TSV</option><option value="csv">CSV</option><option value="jsonl">JSON Lines</option></select></label><label><span>Compression</span><select value={exportCompression} onChange={(event)=>setExportCompression(event.target.value as OutputCompression)}><option value="none">None</option><option value="gzip">Gzip (.gz)</option></select></label><button type="button" disabled={Boolean(busy)} onClick={()=>void downloadActivePopulation()}>Download current population</button></section>
    </section>}

    {busy && <div className="post-progress" role="status"><div><span>{busy}</span><span className="post-progress-actions"><strong>{progress.total ? `${Math.min(100, progress.processed / progress.total * 100).toFixed(1)}%` : "working"}</strong><button type="button" disabled={cancelling || !postLockAbortRef.current} onClick={()=>void cancelActiveOperation()}>{cancelling ? "Restoring…" : "Cancel"}</button></span></div><progress max={Math.max(1, progress.total)} value={progress.processed} /><small>{progress.processed.toLocaleString()} / {progress.total.toLocaleString()} {progress.unit ?? "AIRR records indexed or scanned"} · {postLockState === "held" ? "background-run lock held" : postLockState === "waiting" ? "waiting for another Swig tab" : "Web Locks unavailable"}</small></div>}
    {error && <div className="post-error" role="alert"><strong>Post-analysis stopped</strong><span>{error}</span><button type="button" onClick={() => setError("")}>Dismiss</button></div>}

    <section className={moduleClass("dedup","post-module dedup-module")}>
      <header><div className="module-number">02</div><div><span className="section-kicker">Abundance preservation after optional call reassignment</span><h3>Collapse exact duplicates or denoise read errors</h3><p>When allele pooling is applied, rearrangement-key partitions use the policy-selected V/D/J calls. Every retained representative carries the sum of its source multiplicities in <code>duplicate_count</code>.</p></div><a href="./methods/06_COLLAPSE_AND_DENOISING.md" target="_blank" rel="noreferrer">Method details ↗</a>{skipStepButton("dedup")}<button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("dedup")} onClick={()=>toggleModule("dedup")}>{openModules.has("dedup")?"Collapse ↑":"Expand ↓"}</button></header>
      <div className="collapse-mode-grid" role="radiogroup" aria-label="Collapse method">
        <button type="button" role="radio" aria-checked={collapseMode === "exact"} className={collapseMode === "exact" ? "selected" : ""} onClick={() => { setCollapseMode("exact"); setDedup(null); }}><b>A</b><span><strong>Exact deduplication</strong><small>Collapse identical keys only. No error model.</small></span></button>
        <button type="button" role="radio" aria-checked={collapseMode === "fad"} className={collapseMode === "fad" ? "selected" : ""} onClick={() => { setCollapseMode("fad"); setDenoiseAmbiguousPolicy("exclude"); setDedup(null); }}><b>B</b><span><strong>FAD-compatible denoising</strong><small>Published 6-mer distance and abundance/Poisson rule.</small></span></button>
        <button type="button" role="radio" aria-checked={collapseMode === "conservative"} className={collapseMode === "conservative" ? "selected" : ""} onClick={() => { setCollapseMode("conservative"); setDenoiseAmbiguousPolicy("retain"); setDedup(null); }}><b>C</b><span><strong>Exact-neighbor error model</strong><small>Experimental; conservative, indexed Hamming candidates.</small></span></button>
        <button type="button" role="radio" aria-checked={collapseMode === "indel"} className={collapseMode === "indel" ? "selected" : ""} onClick={() => { setCollapseMode("indel"); setDenoiseAmbiguousPolicy("retain"); setDedup(null); }}><b>D</b><span><strong>Indel-aware error model</strong><small>Complete 1–2 edit index; abundance-directed indel collapse.</small></span></button>
      </div>
      <div className="scope-policy-row"><label><span>Collapse boundary</span><select value={collapseScope} onChange={(event)=>{setCollapseScope(event.target.value as DatasetScope);setDedup(null);}}>{Object.entries(DATASET_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label><p><b>{DATASET_SCOPE_LABELS[collapseScope]}:</b> exact matches and denoising parents are never compared across this boundary. Files sharing a sample ID are treated as technical libraries of the same specimen.</p></div>
      <label className="check-line constant-collapse-policy"><input type="checkbox" checked={respectConstantCall} onChange={(event)=>{setRespectConstantCall(event.target.checked);setDedup(null);}} /><span><strong>Keep different constant-gene/isotype calls separate</strong><small>Default. The top C call is normalized to gene level and used only as a partition; constant-region tail sequence and tail length are not part of the VDJ comparison key. Records without a C call form a separate unassigned partition.</small></span></label>
      {collapseMode === "exact" ? <>
        <div className="module-controls">
          <label><span>Identity key</span><select value={dedupKey} onChange={(event) => setDedupKey(event.target.value as DedupKey)}><option value="sequence">Full input sequence</option><option value="trimmed">VDJ-aligned sequence</option><option value="cdr3">Locus + CDR3 nucleotide</option><option value="rearrangement">Locus + V/J calls + CDR3</option></select></label>
          <label><span>Missing selected key fields</span><select value={denoiseUnresolvedPolicy} onChange={(event) => { setDenoiseUnresolvedPolicy(event.target.value as "discard" | "retain"); setDedup(null); }}><option value="discard">Discard from downstream · default</option><option value="retain">Retain unchanged</option></select></label>
          <button className={guidedClass("run-dedup")} type="button" disabled={Boolean(busy)} onClick={() => void runDedup()}>Run exact deduplication</button>
        </div>
        {(dedupKey === "sequence" || dedupKey === "trimmed") && <p className="scientific-note"><span>i</span>Sequence-key modes compare normalized length plus a 128-bit fingerprint so complete sequence payloads do not remain in memory. Existing <code>duplicate_count</code> values are summed rather than reset. Large runs use deterministic key shards across up to {workers} workers; all copies of one key remain in one shard and the first AIRR row remains the representative.</p>}
      </> : <div className="denoise-config">
        <div className="denoise-essential"><label><span>Missing trimmed sequence or V/J</span><select value={denoiseUnresolvedPolicy} onChange={(event) => { setDenoiseUnresolvedPolicy(event.target.value as "discard" | "retain"); setDedup(null); }}><option value="discard">Discard from downstream · default</option><option value="retain">Retain unchanged</option></select></label></div>
        <details className="post-advanced denoise-advanced"><summary>Advanced {collapseMode === "fad" ? "FAD" : collapseMode === "indel" ? "indel-aware" : "error-model"} settings</summary><div>
        <div className="control-grid three">
          <label><span>V/J call level</span><select value={denoiseResolution} onChange={(event) => setDenoiseResolution(event.target.value as CallResolution)}><option value="allele">Allele</option><option value="gene">Gene</option></select></label>
          <label><span>Uncertain assignments</span><select value={denoiseAmbiguity} onChange={(event) => setDenoiseAmbiguity(event.target.value as "top" | "strict")}><option value="strict">Exact sorted call set</option><option value="top">Top V and J calls</option></select></label>
          <label><span>Ambiguous N bases</span><select value={denoiseAmbiguousPolicy} onChange={(event) => setDenoiseAmbiguousPolicy(event.target.value as "exclude" | "retain")}><option value="exclude">Exclude (published FAD behavior)</option><option value="retain">Keep as exact-only representatives</option></select></label>
          <label><span>Per-base error rate</span><CommitNumberInput min="0.000001" max="0.2" step="0.00001" value={denoiseErrorRate} onCommit={setDenoiseErrorRate} /></label>
          <label><span>Significance α</span><CommitNumberInput min="0.000001" max="0.5" step="0.001" value={denoiseAlpha} onCommit={setDenoiseAlpha} /></label>
          <label><span>Minimum parent abundance</span><CommitNumberInput min="1" max="1000000" step="1" value={minimumParentCount} onCommit={setMinimumParentCount} /></label>
        </div>
        {collapseMode === "fad" ? <div className="control-grid three denoise-specific">
          <label><span>FAD decision</span><select value={fadMethod} onChange={(event) => setFadMethod(Number(event.target.value) as 1 | 2)}><option value={2}>Method 2 · Poisson</option><option value={1}>Method 1 · abundance only</option></select></label>
          <label><span>Corrected 6-mer radius</span><CommitNumberInput min="0" max="5" step="0.25" value={fadNeighborThreshold} onCommit={setFadNeighborThreshold} /></label>
          <label><span>Expected zero-error fraction</span><CommitNumberInput min="0.0001" max="1" step="0.01" value={expectedZeroErrorFraction} onCommit={setExpectedZeroErrorFraction} /></label>
        </div> : collapseMode === "conservative" ? <div className="control-grid three denoise-specific"><label><span>Maximum Hamming errors</span><select value={maximumDenoiseDistance} onChange={(event) => setMaximumDenoiseDistance(Number(event.target.value))}><option value={1}>1 substitution</option><option value={2}>2 substitutions</option><option value={3}>3 substitutions</option></select></label></div> : <div className="control-grid three denoise-specific">
          <label><span>Maximum edit distance</span><select value={maximumEditDistance} onChange={(event) => setMaximumEditDistance(Number(event.target.value))}><option value={1}>1 edit</option><option value={2}>2 edits</option></select></label>
          <label><span>Minimum indel parent:child ratio</span><CommitNumberInput min="1.01" max="1000000" step="0.1" value={minimumIndelParentRatio} onCommit={setMinimumIndelParentRatio} /></label>
        </div>}
        <details className="post-advanced"><summary>Performance guard</summary><div className="control-grid three"><label><span>Candidate cap / variant</span><CommitNumberInput min="100" max="1000000" step="1000" value={denoiseCandidateCap} onCommit={setDenoiseCandidateCap} /></label></div></details>
        <div className="algorithm-note"><strong>{collapseMode === "fad" ? "FAD-compatible indexed implementation" : collapseMode === "indel" ? "Complete bounded-edit implementation" : "Conservative exact-neighbor implementation"}</strong><span>{collapseMode === "fad" ? "Partition locus/V/J → exact dereplication → complete corrected-6-mer radius index → published template test → exact nearest accepted centroid within each partition." : collapseMode === "indel" ? "Partition locus/V/J → exact dereplication → complete length-aware d+1 segment join → exact banded edit profiling. Low-abundance indel paths use the parent:child ratio; substitution-only paths retain method C’s Poisson test. No distant-centroid assignment is performed." : "Partition locus/V/J → exact dereplication → d+1 block index → exact Hamming verification → sequence-specific Poisson error test. Isolated singletons are retained; no read is assigned to a distant centroid."} Independent study/C/locus/V/J partitions are dynamically scheduled across up to {workers} workers; no partition is split internally, so serial decisions and tie-breaks are unchanged.</span></div>
        <p className="scientific-note"><span>i</span>Trimmed VDJ sequences are streamed from the AIRR store into a two-bit packed arena. Temporary neighbor profiles are released one V/J partition at a time. The 0.00473 default reproduces the linked workflow’s MiSeq median error-rate setting; it should be changed for a different assay.</p>
        </div></details>
        {collapseMode === "indel" && <p className="scientific-note warning"><span>!</span>Method D deliberately treats a one- or two-edit path containing an indel as likely technical error once the parent:child abundance ratio is met. True biological length variants can therefore be merged; set a stricter ratio or use method C when that assumption is inappropriate.</p>}
        <button className={guidedClass("run-dedup","post-primary denoise-run")} type="button" disabled={Boolean(busy)} onClick={() => void runDedup()}>{collapseMode === "fad" ? "Run FAD-compatible denoising" : collapseMode === "indel" ? "Run indel-aware denoising" : "Run exact-neighbor denoising"}</button>
      </div>}
      {dedup && <div className="module-result"><div className="post-stat-grid"><article><span>Input rows / abundance</span><strong>{dedup.inputRecords.toLocaleString()} / {dedup.inputAbundance.toLocaleString()}</strong></article><article><span>Retained representatives</span><strong>{dedup.uniqueRecords.toLocaleString()}</strong></article><article><span>Collapsed or excluded rows</span><strong>{dedup.collapsedRecords.toLocaleString()}</strong></article><article><span>Largest multiplicity</span><strong>{(dedup.largestGroups[0]?.count ?? 1).toLocaleString()}</strong></article>{dedup.mode !== "exact" && <><article><span>V/J partitions</span><strong>{dedup.partitions.toLocaleString()}</strong></article><article><span>Verified candidates</span><strong>{dedup.candidateComparisons.toLocaleString()}</strong></article></>}{dedup.unresolvedRecords > 0 && <article><span>{denoiseUnresolvedPolicy === "discard" ? "Ineligible discarded" : "Ineligible retained"}</span><strong>{dedup.unresolvedRecords.toLocaleString()}</strong></article>}{dedup.mode === "indel" && <><article><span>Indel variants merged</span><strong>{dedup.indelMergedVariants.toLocaleString()}</strong></article><article><span>Substitution variants merged</span><strong>{dedup.substitutionMergedVariants.toLocaleString()}</strong></article></>}</div><div className="denoise-provenance"><strong>{dedup.algorithm}</strong><span>{dedup.mode === "exact" ? `Key: ${dedup.key} · ineligible ${denoiseUnresolvedPolicy}` : `${denoiseResolution}-level V/J partitioning · ${denoiseAmbiguity} assignment policy · ineligible ${denoiseUnresolvedPolicy}`} · up to {workers} compute workers</span></div>{dedup.warnings.map((warning) => <div key={warning} className="scientific-note warning"><span>!</span><p>{warning}</p></div>)}<div className="filter-commit"><div><span>Downstream action</span><strong>Retain representatives; preserve collapsed abundance as counts</strong><p>This action changes the working set from {dedup.inputRecords.toLocaleString()} to {dedup.uniqueRecords.toLocaleString()} rows and invalidates downstream results.</p></div><button className={guidedClass("apply-dedup")} type="button" disabled={Boolean(busy) || workingStages.some((stage) => stage.id === "dedup")} onClick={() => void applyDedupFilter()}>{workingStages.some((stage) => stage.id === "dedup") ? "Applied downstream" : `Apply ${dedup.uniqueRecords.toLocaleString()} representatives`}</button></div><div className="result-actions"><button type="button" onClick={() => void downloadDeduplicated()}>Download collapsed AIRR + multiplicity</button></div></div>}
    </section>

    <section className={moduleClass("chimera","post-module chmm-module")}>
      <header><div className="module-number amber">03</div><div><span className="section-kicker">Optional PCR-chimera model</span><h3>CHMMAIRRa after V(D)J assignment</h3><p>The browser port threads each AIRR local V or J alignment onto a reference MSA, then evaluates the CHMMera posterior. V is the manuscript default; D is not modeled.</p></div><a href="./methods/07_CHIMERA_INFERENCE.md" target="_blank" rel="noreferrer">Method details ↗</a>{skipStepButton("chimera")}<button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("chimera")} onClick={()=>toggleModule("chimera")}>{openModules.has("chimera")?"Collapse ↑":"Expand ↓"}</button></header>
      <div className="chmm-grid">
        <div className="chmm-config">
          <div className="control-grid three"><label><span>Segment</span><select value={chmmSegment} onChange={(event) => { setChmmSegment(event.target.value as ChmmSegment); setPreparedMsa(""); setChmm(null); setChmmRun(null); setChimeraDetail(null); }}><option value="V">V (recommended)</option><option value="J">J (optional)</option></select></label><label><span>Model</span><select value={chmmMethod} onChange={(event) => setChmmMethod(event.target.value as "BW" | "DB")}><option value="BW">Baum–Welch · IG default</option><option value="DB">Discretized Bayesian · TCR default</option></select></label><label><span>Posterior threshold</span><CommitNumberInput min="0" max="1" step="0.01" value={chmmThreshold} onCommit={setChmmThreshold} /></label></div>
          <details className="post-advanced chimera-reference-settings"><summary>Reference MSA and advanced model parameters</summary><div><fieldset className="msa-source"><legend>Reference multiple-sequence alignment</legend><label className={chmmSource === "selected" ? "selected" : ""}><input type="radio" checked={chmmSource === "selected"} onChange={() => { setChmmSource("selected"); setPreparedMsa(""); setChmm(null); setChmmRun(null); setChimeraDetail(null); }} /><span><strong>Build from this run’s assignment {chmmSegment} references</strong><small>Kalign 3.3.1 WASM; uses the exact IMGT/KI/local-file composition and allele exclusions already applied before V(D)J assignment.</small></span></label><label className={chmmSource === "upload" ? "selected" : ""}><input type="radio" checked={chmmSource === "upload"} onChange={() => { setChmmSource("upload"); setPreparedMsa(""); setChmm(null); setChmmRun(null); setChimeraDetail(null); }} /><span><strong>Use an aligned FASTA MSA from file</strong><small>{uploadedMsaName || "Every record must have equal aligned length and names matching AIRR calls."}</small></span><input className="file-inline" type="file" accept=".fa,.fasta,.fas,.aln,.txt" onChange={(event) => void acceptMsa(event)} /></label></fieldset><div className="control-grid three"><label><span>Chimera prior</span><CommitNumberInput min="0.00001" max="0.5" step="0.01" value={chmmPrior} onCommit={setChmmPrior} /></label><label><span>Minimum DFR</span><CommitNumberInput min="0" max="100" step="1" value={chmmMinDfr} onCommit={setChmmMinDfr} /></label><label><span>DB mutation rates</span><CommitTextInput value={mutationRates} onCommit={setMutationRates} /></label><label className="check-line"><input type="checkbox" checked={chmmDetailed} onChange={(event) => setChmmDetailed(event.target.checked)} /><span>Precompute breakpoint labels during the repertoire scan (the full path remains on-demand)</span></label></div></div></details>
          <div className="scientific-note warning"><span>!</span><p>Reference completeness matters: an absent true V/J allele can produce a false switch signal. Loaded MSAs are never silently supplemented. Low-DFR records are reported as unevaluated rather than forced through the model.</p></div>
          <div className="result-actions"><button className={guidedClass("run-chimera","post-primary amber")} type="button" disabled={Boolean(busy) || (chmmSource === "upload" && !uploadedMsa)} onClick={() => void runChmmAnalysis()}>Run CHMMAIRRa on {workingCount.toLocaleString()}</button>{preparedMsa && <button type="button" onClick={() => downloadText(preparedMsa, `${baseName(inputName)}.${chmmSegment.toLowerCase()}-reference-msa.fasta`)}>Download reference MSA</button>}</div>
        </div>
        <div className="chmm-result-panel">
          {chmm ? <>
            <div className="post-stat-grid compact"><article><span>Working-set input</span><strong>{chmm.inputRecords.toLocaleString()}</strong></article><article><span>Evaluated</span><strong>{chmm.evaluated.toLocaleString()}</strong></article><article><span>Posterior ≥ {chmm.threshold}</span><strong>{chmm.flagged.toLocaleString()}</strong></article><article><span>Below DFR</span><strong>{chmm.lowDfr.toLocaleString()}</strong></article><article><span>Missing reference</span><strong>{chmm.missingReference.toLocaleString()}</strong></article>{chmm.upstreamExcluded > 0 && <article><span>Excluded upstream</span><strong>{chmm.upstreamExcluded.toLocaleString()}</strong></article>}</div>
            <div className="chmm-filter-commit">
              <header><div><span className="section-kicker">Explicit downstream exclusion</span><h4>Filter chimeras from later steps</h4><p>Change the threshold without rerunning the HMM, inspect the resulting counts, then apply the filter.</p></div></header>
              <div className="chmm-filter-controls"><label><span>Exclude posterior ≥</span><CommitNumberInput min="0" max="1" step="0.01" value={chmmFilterThreshold} onCommit={setChmmFilterThreshold} /></label><label className="check-line"><input type="checkbox" checked={retainUnevaluated} onChange={(event) => setRetainUnevaluated(event.target.checked)} /><span>Retain unevaluated records</span></label></div>
              {chmmFilterPreview && <div className="filter-preview"><div><span>Retained downstream</span><strong>{chmmFilterPreview.retained.toLocaleString()}</strong></div><div><span>Excluded as chimera</span><strong>{chmmFilterPreview.excluded.toLocaleString()}</strong></div><div><span>Unevaluated in input</span><strong>{chmmFilterPreview.unevaluated.toLocaleString()}</strong></div></div>}
              {!chmmFilterThresholdValid && <div className="inline-field-error">Enter a posterior threshold from 0 to 1.</div>}
              <button className={guidedClass("apply-chimera","post-primary amber")} type="button" disabled={Boolean(busy) || !chmmFilterThresholdValid} onClick={() => void applyChimeraFilter()}>Apply chimera filter downstream</button>
            </div>
            <BarChart title={`${chmm.segment} chimera posterior`} subtitle="Evaluated rearrangements by posterior interval" data={chmm.histogram.map((item) => ({ label: item.label, value: item.count }))} color="#d49a19" name={`${baseName(inputName)}.chmmairra-posterior.svg`} />
            {chmm.top.length > 0 && <div className="chmm-top-list"><header><strong>Highest posteriors</strong><span>Click to run the detailed Viterbi parent view</span></header>{chmm.top.slice(0, 12).map((record) => {const indexed=chmmTopIndex.get(record.ordinal);return <button type="button" key={record.ordinal} onClick={() => void openChimera(record.ordinal)}><b>#{(record.ordinal + 1).toLocaleString()}</b><span>{(record.probability * 100).toFixed(2)}%</span><strong title={indexed?.sequenceId}>{indexed?.sequenceId||"sequence loading…"}</strong><code title={indexed?.cdr3Aa||indexed?.cdr3}>{indexed?.cdr3Aa||indexed?.cdr3||"CDR3 —"}</code><small>DFR {record.dfr}{record.recombinations.length ? ` · ${record.recombinations.map((event) => `${event.left}→${event.right}@${event.position}`).join("; ")}` : " · Viterbi path on click"}</small></button>;})}</div>}
            <button type="button" onClick={() => void downloadChmm()}>Download CHMMAIRRa TSV</button>
          </> : <div className="method-placeholder"><span>HMM</span><h4>No CHMMAIRRa run</h4><p>{isTcr ? "TCR mode defaults to a fixed 0.005 mutation-rate state (DB)." : "IG mode defaults to per-reference Baum–Welch mutation estimates."}</p></div>}
        </div>
      </div>
      {chimeraDetail && <div ref={chimeraDetailRef}><ChimeraHighlighter detail={chimeraDetail} name={`${baseName(inputName)}.record-${chimeraDetail.ordinal + 1}.chmmairra-viterbi.svg`} onInspect={onInspect} /></div>}
    </section>

    <section className={moduleClass("selection","post-module selection-module")}>
      <header><div className="module-number dark">04</div><div><span className="section-kicker">Composable repertoire population</span><h3>Select the records used downstream</h3><p>Combine assignment, CDR3, motif, quality, SHM, and double-D evidence filters. This step is skipped by default; include it when you need an explicit filtered population.</p></div><a href="./methods/05_STORAGE_DASHBOARD_SELECTION.md" target="_blank" rel="noreferrer">Method details ↗</a>{skipStepButton("selection")}<button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("selection")} onClick={()=>toggleModule("selection")}>{openModules.has("selection")?"Collapse ↑":"Expand ↓"}</button></header>
      <div className="control-grid four selection-call-grid selection-common-grid">
        <label><span>Sequence ID contains</span><CommitTextInput value={selectionDraft.sequenceId} onCommit={(sequenceId)=>{setSelectionDraft(value=>({...value,sequenceId}));setSelectionPreview(null);}} placeholder="one or more IDs" /></label>
        {selectionFacets.samples.length>0&&<FacetPicker label="Sample" value={selectionDraft.sampleId} items={selectionFacets.samples} multiple placeholder="Any sample" onChange={(sampleId)=>{setSelectionDraft(value=>({...value,sampleId}));setSelectionPreview(null);}}/>}
        <FacetPicker label="Locus" value={selectionDraft.locus} items={selectionFacets.loci} multiple placeholder="Any locus" onChange={(locus)=>{setSelectionDraft(value=>({...value,locus}));setSelectionPreview(null);}}/>
        <div className="call-filter selection-call-picker"><FacetPicker label="V gene or allele" value={selectionDraft.vCall} items={selectionFacets.vCalls} multiple allowCustom placeholder="Any V call" help="Choose one or more observed V genes or alleles, or type a value. Selected values are combined with OR." onChange={(vCall)=>{setSelectionDraft(value=>({...value,vCall}));setSelectionPreview(null);}}/><label className="check-filter compact"><input type="checkbox" checked={selectionDraft.vCallIncludeAmbiguous} onChange={(event)=>{setSelectionDraft(value=>({...value,vCallIncludeAmbiguous:event.target.checked}));setSelectionPreview(null);}}/><span>Include multi-call assignments containing a selected V</span></label></div>
        <div className="call-filter selection-call-picker"><FacetPicker label="J gene or allele" value={selectionDraft.jCall} items={selectionFacets.jCalls} multiple allowCustom placeholder="Any J call" help="Choose one or more observed J genes or alleles, or type a value. Selected values are combined with OR." onChange={(jCall)=>{setSelectionDraft(value=>({...value,jCall}));setSelectionPreview(null);}}/><label className="check-filter compact"><input type="checkbox" checked={selectionDraft.jCallIncludeAmbiguous} onChange={(event)=>{setSelectionDraft(value=>({...value,jCallIncludeAmbiguous:event.target.checked}));setSelectionPreview(null);}}/><span>Include multi-call assignments containing a selected J</span></label></div>
        <label><span>CDR3 nucleotide contains</span><CommitTextInput value={selectionDraft.cdr3Nt} onCommit={(cdr3Nt)=>{setSelectionDraft(value=>({...value,cdr3Nt}));setSelectionPreview(null);}} placeholder="literal substring" /></label>
        <label><span>CDR3 amino acid contains</span><CommitTextInput value={selectionDraft.cdr3Aa} onCommit={(cdr3Aa)=>{setSelectionDraft(value=>({...value,cdr3Aa}));setSelectionPreview(null);}} placeholder="literal substring" /></label>
      </div>
      <details className="post-advanced"><summary>Study metadata filters</summary><div className="control-grid four">
        <FacetPicker label="Dataset" value={selectionDraft.datasetId} items={selectionFacets.datasets} multiple placeholder="Any dataset" onChange={(datasetId)=>{setSelectionDraft(value=>({...value,datasetId}));setSelectionPreview(null);}}/>
        <FacetPicker label="Donor / subject" value={selectionDraft.subjectId} items={selectionFacets.subjects} multiple placeholder="Any donor" onChange={(subjectId)=>{setSelectionDraft(value=>({...value,subjectId}));setSelectionPreview(null);}}/>
        <FacetPicker label="Cohort" value={selectionDraft.cohort} items={selectionFacets.cohorts} multiple placeholder="Any cohort" onChange={(cohort)=>{setSelectionDraft(value=>({...value,cohort}));setSelectionPreview(null);}}/>
        <FacetPicker label="Timepoint" value={selectionDraft.timepoint} items={selectionFacets.timepoints} multiple placeholder="Any timepoint" onChange={(timepoint)=>{setSelectionDraft(value=>({...value,timepoint}));setSelectionPreview(null);}}/>
        <FacetPicker label="Compartment / tissue" value={selectionDraft.compartment} items={selectionFacets.compartments} multiple placeholder="Any compartment" onChange={(compartment)=>{setSelectionDraft(value=>({...value,compartment}));setSelectionPreview(null);}}/>
      </div></details>
      <details className="post-advanced"><summary>Additional D, constant-gene, isotype, and rare-event filters</summary><div className="control-grid four">
        <div className="call-filter selection-call-picker"><FacetPicker label="D1 gene or allele" value={selectionDraft.d1Call} items={selectionFacets.dCalls} multiple allowCustom placeholder="Any D1 call" onChange={(d1Call)=>{setSelectionDraft(value=>({...value,d1Call}));setSelectionPreview(null);}}/><label className="check-filter compact"><input type="checkbox" checked={selectionDraft.d1CallIncludeAmbiguous} onChange={(event)=>{setSelectionDraft(value=>({...value,d1CallIncludeAmbiguous:event.target.checked}));setSelectionPreview(null);}}/><span>Include multi-call assignments containing a selected D1</span></label></div>
        <label><span>D2 call contains</span><CommitTextInput value={selectionDraft.d2Call} onCommit={(d2Call)=>{setSelectionDraft(value=>({...value,d2Call}));setSelectionPreview(null);}} placeholder="double-D only" /></label>
        <div className="call-filter selection-call-picker"><FacetPicker label="Constant gene or allele" value={selectionDraft.cCall} items={selectionFacets.cCalls} multiple allowCustom placeholder="Any constant call" onChange={(cCall)=>{setSelectionDraft(value=>({...value,cCall}));setSelectionPreview(null);}}/><label className="check-filter compact"><input type="checkbox" checked={selectionDraft.cCallIncludeAmbiguous} onChange={(event)=>{setSelectionDraft(value=>({...value,cCallIncludeAmbiguous:event.target.checked}));setSelectionPreview(null);}}/><span>Include multi-call assignments containing a selected C</span></label></div>
        <FacetPicker label="Isotype / constant class" value={selectionDraft.isotype} items={selectionFacets.isotypes} multiple allowCustom placeholder="Any isotype" onChange={(isotype)=>{setSelectionDraft(value=>({...value,isotype}));setSelectionPreview(null);}}/>
        {doubleDCount>0&&<label><span>Double-D evidence</span><select value={selectionDraft.doubleD} onChange={(event)=>{setSelectionDraft(value=>({...value,doubleD:event.target.value as RepertoireSelectionOptions["doubleD"]}));setSelectionPreview(null);}}><option value="any">Any record</option><option value="positive">Supported Double-D only</option><option value="negative">Exclude supported Double-D</option></select></label>}
      </div></details>
      <details className="post-advanced motif-filter"><summary>Sequence motif filter</summary><div><label><span>Motif(s)</span><CommitTextarea value={selectionDraft.motif} onCommit={(motif)=>{setSelectionDraft(value=>({...value,motif}));setSelectionPreview(null);}} placeholder="one motif per line; leave blank to disable" /></label><div className="control-grid three"><label><span>Target</span><select value={selectionDraft.motifTarget} onChange={(event)=>{setSelectionDraft(value=>({...value,motifTarget:event.target.value as RepertoireSelectionOptions["motifTarget"]}));setSelectionPreview(null);}}><option value="sequence">Full sequence / alignment</option><option value="cdr3_nt">CDR3 nucleotide</option><option value="cdr3_aa">CDR3 amino acid</option><option value="junction_aa">Junction amino acid</option></select></label><label><span>Syntax</span><select value={selectionDraft.motifSyntax} onChange={(event)=>{setSelectionDraft(value=>({...value,motifSyntax:event.target.value as RepertoireSelectionOptions["motifSyntax"]}));setSelectionPreview(null);}}><option value="substring">Literal substring</option><option value="iupac">IUPAC nucleotide</option><option value="regex">Regular expression</option></select></label><label><span>Multiple motifs</span><select value={selectionDraft.motifMode} onChange={(event)=>{setSelectionDraft(value=>({...value,motifMode:event.target.value as "any"|"all"}));setSelectionPreview(null);}}><option value="any">Match any</option><option value="all">Match all</option></select></label></div></div></details>
      <details className="post-advanced"><summary>Quality, length, identity, and mutation filters</summary><div className="control-grid four">
        {([['productive','Productive'],['completeVdj','Complete V(D)J'],['vjInFrame','VJ in frame'],['stopCodon','Stop codon'],['hasD','D assigned'],['hasCdr3','CDR3 assigned']] as Array<[keyof RepertoireSelectionOptions,string]>).map(([field,label])=><label key={field}><span>{label}</span><select value={String(selectionDraft[field])} onChange={(event)=>{setSelectionDraft(value=>({...value,[field]:event.target.value}));setSelectionPreview(null);}}><option value="any">Any</option><option value="yes">Yes</option><option value="no">No</option></select></label>)}
        <label><span>Minimum CDR3 nt length</span><CommitNumberInput min="0" value={selectionDraft.minCdr3NtLength} onCommit={(minCdr3NtLength)=>{setSelectionDraft(value=>({...value,minCdr3NtLength}));setSelectionPreview(null);}} /></label><label><span>Maximum CDR3 nt length</span><CommitNumberInput min="0" value={selectionDraft.maxCdr3NtLength} onCommit={(maxCdr3NtLength)=>{setSelectionDraft(value=>({...value,maxCdr3NtLength}));setSelectionPreview(null);}} /></label>
        <label><span>Minimum CDR3 aa length</span><CommitNumberInput min="0" value={selectionDraft.minCdr3AaLength} onCommit={(minCdr3AaLength)=>{setSelectionDraft(value=>({...value,minCdr3AaLength}));setSelectionPreview(null);}} /></label><label><span>Maximum CDR3 aa length</span><CommitNumberInput min="0" value={selectionDraft.maxCdr3AaLength} onCommit={(maxCdr3AaLength)=>{setSelectionDraft(value=>({...value,maxCdr3AaLength}));setSelectionPreview(null);}} /></label>
        <label><span>Minimum V identity</span><CommitNumberInput min="0" max="1" step="0.01" value={selectionDraft.minVIdentity} onCommit={(minVIdentity)=>{setSelectionDraft(value=>({...value,minVIdentity}));setSelectionPreview(null);}} /></label><label><span>Minimum J identity</span><CommitNumberInput min="0" max="1" step="0.01" value={selectionDraft.minJIdentity} onCommit={(minJIdentity)=>{setSelectionDraft(value=>({...value,minJIdentity}));setSelectionPreview(null);}} /></label>
        <label><span>Minimum V mutation fraction</span><CommitNumberInput min="0" max="1" step="0.01" value={selectionDraft.minVMutation} onCommit={(minVMutation)=>{setSelectionDraft(value=>({...value,minVMutation}));setSelectionPreview(null);}} /></label><label><span>Maximum V mutation fraction</span><CommitNumberInput min="0" max="1" step="0.01" value={selectionDraft.maxVMutation} onCommit={(maxVMutation)=>{setSelectionDraft(value=>({...value,maxVMutation}));setSelectionPreview(null);}} /></label>
      </div></details>
      <div className="selection-commit"><div><span>Inherited input</span><strong>{(selectionPreview?.inputRecords??workingCount).toLocaleString()} records</strong><small>{workingStages.filter(stage=>stage.id!=="selection").map(stage=>stage.label).join(" → ")||"Complete assigned input"}</small></div><span className="selection-arrow">→</span><div><span>Preview retained</span><strong>{selectionPreview?selectionPreview.retainedRecords.toLocaleString():"Run preview"}</strong><small>{selectionPreview?`${selectionPreview.discardedRecords.toLocaleString()} would be excluded`:"Editing fields does not rescan automatically"}</small></div><div className="result-actions"><button className={guidedAction==="preview-selection"?"guided-next":""} type="button" disabled={Boolean(busy)} onClick={()=>void previewSelection()}>Preview count</button><button className={guidedClass("apply-selection","post-primary dark")} type="button" disabled={Boolean(busy)||!selectionPreview} onClick={()=>void applySelection()}>{selectionApplied?"Re-apply selection":"Apply selection downstream"}</button>{selectionApplied?<button type="button" onClick={()=>void removeSelection()}>Remove selection stage</button>:null}</div></div>
      <p className="scientific-note"><span>i</span>Fields are committed only on Enter or blur. Previewing performs the scan; typing does not repeatedly resolve a million-row dataset. All non-empty conditions compose with logical AND, while comma/newline-separated call filters and motifs use their stated any/all rule.</p>
    </section>

    <section className={moduleClass("alleles","post-module allele-refinement-module")}>
      <header><div className="module-number amber">01</div><div><span className="section-kicker">Optional first-stage repertoire assignment model</span><h3>Resolve ambiguous germline calls by pooling repertoire evidence</h3><p>Choose the continuous Dirichlet mixture or the fast exact-zero hurdle model before collapse or filtering, so policy-selected V/D/J calls define every downstream partition. The default is one independent fit per donor, combining that donor's samples without crossing participant IDs.</p></div><a href="./REPERTOIRE_ALLELE_REFINEMENT.md" target="_blank" rel="noreferrer">Method details ↗</a>{skipStepButton("alleles")}<button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("alleles")} onClick={()=>toggleModule("alleles")}>{openModules.has("alleles")?"Collapse ↑":"Expand ↓"}</button></header>
      <AlleleRefinementPanel references={references} options={alleleOptions} onOptionsChange={(next)=>{setAlleleOptions(next);setAlleleRefinement(null);discardAppliedAllelePolicy();}} result={alleleRefinement} applied={alleleApplied} reassignmentPolicy={alleleReassignmentPolicy} onReassignmentPolicyChange={(policy)=>{setAlleleReassignmentPolicy(policy);discardAppliedAllelePolicy();}} applyMinimumPosterior={alleleApplyMinimumPosterior} onApplyMinimumPosteriorChange={(value)=>{setAlleleApplyMinimumPosterior(Math.max(0,Math.min(1,value)));discardAppliedAllelePolicy();}} busy={Boolean(busy)} progress={alleleProgress} onRun={()=>void runAlleleRefinement()} onCancel={()=>void cancelActiveOperation()} onApply={()=>void applyAlleleRefinement()} onReset={()=>void resetAlleleRefinement()} onDownloadModel={downloadAlleleModel} onDownloadSidecar={()=>void downloadAlleleSidecar()} onDownloadAirr={()=>void downloadRefinedAirr()} />
    </section>

    <section className={moduleClass("lineage","post-module lineage-module")}>
      <header><div className="module-number">05</div><div><span className="section-kicker">Repertoire-scale clonal grouping</span><h3>Assign lineages from CDR3 nucleotide distance</h3><p>Default: same locus, overlapping V/J gene assignments, exact CDR3 nucleotide length, and single-linkage at ≥85% identity. The threshold is a starting point and remains dataset-adjustable.</p></div><a href="./methods/08_LINEAGE_ASSIGNMENT.md" target="_blank" rel="noreferrer">Method details ↗</a>{skipStepButton("lineage")}<button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("lineage")} onClick={()=>toggleModule("lineage")}>{openModules.has("lineage")?"Collapse ↑":"Expand ↓"}</button></header>
      <div className="lineage-config">
        <div className="control-grid three"><label><span>Study boundary</span><select value={lineageScope} onChange={(event)=>{setLineageScope(event.target.value as DatasetScope);setLineages(null);}}>{Object.entries(DATASET_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label><label><span>CDR3 identity</span><div className="range-number"><input type="range" min="0.7" max="1" step="0.01" value={identity} onChange={(event) => setIdentity(Number(event.target.value))} /><b>{Math.round(identity * 100)}%</b></div></label><label className="check-line"><input type="checkbox" checked={productiveOnly} onChange={(event) => setProductiveOnly(event.target.checked)} /><span>Productive only</span></label></div>
        <details className="post-advanced"><summary>Advanced call matching and performance settings</summary><div className="control-grid three"><label><span>Call level</span><select value={resolution} onChange={(event) => setResolution(event.target.value as CallResolution)}><option value="gene">Gene</option><option value="allele">Allele</option></select></label><label><span>Ambiguous calls</span><select value={ambiguity} onChange={(event) => setAmbiguity(event.target.value as AmbiguityPolicy)}><option value="overlap">Any assignment overlaps</option><option value="top">Top call only</option><option value="strict">Exact call sets</option></select></label><label><span>Candidate cap / record</span><CommitNumberInput min="100" max="1000000" step="1000" value={candidateCap} onCommit={setCandidateCap} /></label></div><div className="algorithm-note"><strong>Exact accelerated single-linkage · {DATASET_SCOPE_LABELS[lineageScope]}</strong><span>Partition by study boundary → locus → V/J calls → CDR3 length → d+1 exact blocks → verify normalized Hamming distance → union-find components.</span></div></details>
        <div className="current-step-input"><span>Input inherited from applied filters</span><strong>{workingCount.toLocaleString()} active records</strong><small>{workingStages.length ? workingStages.map((stage) => stage.label).join(" → ") : "No upstream exclusion applied"}</small></div><button className={guidedClass("run-lineages")} type="button" disabled={Boolean(busy)} onClick={() => void runLineages()}>Assign lineages on current set</button>
      </div>
      {lineages && <div className="lineage-results"><div className="post-stat-grid"><article><span>Original lineages</span><strong>{lineages.lineageCount.toLocaleString()}</strong></article><article><span>Assigned records</span><strong>{lineages.assignedRecords.toLocaleString()}</strong></article><article><span>Largest lineage</span><strong>{(lineages.summaries[0]?.abundance ?? 0).toLocaleString()}</strong></article><article><span>Exact comparisons</span><strong>{lineages.candidateComparisons.toLocaleString()}</strong></article></div>{lineages.truncatedCandidates > 0 && <div className="scientific-note warning"><span>!</span><p>{lineages.truncatedCandidates.toLocaleString()} records reached the candidate cap. Increase it and rerun before treating components as complete.</p></div>}<div className="result-actions"><button type="button" onClick={() => void downloadLineages()}>Download AIRR + original/merged lineage IDs</button><button className="post-primary dark" type="button" disabled={!lineageTableSelectedIds.size||Boolean(busy)} onClick={()=>void openSelectedLineagesTogether()}>Open selected together ({lineageTableSelectedIds.size})</button>{lineageTableSelectedIds.size>0&&<button type="button" onClick={()=>setLineageTableSelectedIds(new Set())}>Clear selection</button>}</div>{lineageMerges.length > 0 && <div className="lineage-merge-register"><header><div><span className="section-kicker">Derived merge register</span><strong>{lineageMerges.length.toLocaleString()} merged lineage{lineageMerges.length === 1 ? "" : "s"}</strong></div><small>Original clone_id values never change. The derived swig_merged_lineage_id column is added to exports.</small></header>{lineageMerges.map((merge) => <div key={merge.id}><code>{merge.id}</code><span>{merge.originalLineageIds.map((id) => `L${id}`).join(" + ")}</span><button type="button" onClick={() => removeLineageMerge(merge.id)}>Remove merge</button></div>)}</div>}<div className="chart-customizer"><label><span>Gene chart metric</span><select value={geneMetric} onChange={(event) => setGeneMetric(event.target.value as "abundance" | "lineages")}><option value="abundance">Sequence abundance</option><option value="lineages">Lineage count</option></select></label><label><span>Top genes</span><CommitNumberInput min="5" max="24" value={topGenes} onCommit={setTopGenes} /></label><label><span>Figure color</span><input type="color" value={chartColor} onChange={(event) => setChartColor(event.target.value)} /></label></div><div className="post-chart-grid"><BarChart title="Lineage abundance distribution" subtitle="Lineage count in each abundance interval" data={lineages.sizeHistogram.map((item) => ({ label: item.label, value: item.count }))} color={chartColor} name={`${baseName(inputName)}.lineage-size-distribution.svg`} /><BarChart title="Largest lineages" subtitle="Abundance retained after deduplication" data={lineages.summaries.slice(0, 20).map((item) => ({ label: `Lineage ${item.id} · ${item.studyGroup}`, value: item.abundance }))} color={chartColor} name={`${baseName(inputName)}.largest-lineages.svg`} /><BarChart title="V germline use by lineage" subtitle={geneMetric === "abundance" ? "Sequence abundance across lineage representatives" : "Number of lineages represented"} data={vChart} color={chartColor} name={`${baseName(inputName)}.lineage-v-use.svg`} /><BarChart title="J germline use by lineage" subtitle={geneMetric === "abundance" ? "Sequence abundance across lineage representatives" : "Number of lineages represented"} data={jChart} color={chartColor} name={`${baseName(inputName)}.lineage-j-use.svg`} /></div><div className="lineage-table-wrap"><table><thead><tr><th>Select</th><th><button type="button" onClick={()=>sortLineages("id")}>Original lineage{sortIndicator("id")}</button></th><th>Merged lineage</th><th><button type="button" onClick={()=>sortLineages("studyGroup")}>{DATASET_SCOPE_LABELS[lineageScope]}{sortIndicator("studyGroup")}</button></th><th><button type="button" onClick={()=>sortLineages("abundance")}>Abundance{sortIndicator("abundance")}</button></th><th><button type="button" onClick={()=>sortLineages("uniqueMembers")}>Unique{sortIndicator("uniqueMembers")}</button></th><th><button type="button" onClick={()=>sortLineages("locus")}>Locus{sortIndicator("locus")}</button></th><th><button type="button" onClick={()=>sortLineages("vCalls")}>V calls{sortIndicator("vCalls")}</button></th><th><button type="button" onClick={()=>sortLineages("jCalls")}>J calls{sortIndicator("jCalls")}</button></th><th>Representative CDR3</th><th><button type="button" onClick={()=>sortLineages("cdr3Length")}>CDR3 nt{sortIndicator("cdr3Length")}</button></th><th /></tr></thead><tbody>{sortedLineageSummaries.slice(0, 250).map((summary) => {const representative=lineageRepresentativeCdr3(summary);return <tr key={summary.id} className={`${selectedLineageIds.includes(summary.id)?"selected":""}${lineageTableSelectedIds.has(summary.id)?" multi-selected":""}`.trim()} onClick={() => void openLineage(summary)}><td onClick={(event)=>event.stopPropagation()}><input aria-label={`Select lineage ${summary.id}`} type="checkbox" checked={lineageTableSelectedIds.has(summary.id)} onChange={(event)=>setLineageTableSelected(summary.id,event.target.checked)}/></td><td><strong>{summary.id}</strong></td><td>{mergedIdByOriginal.get(summary.id) || "—"}</td><td>{summary.studyGroup}</td><td>{summary.abundance.toLocaleString()}</td><td>{summary.uniqueMembers.toLocaleString()}</td><td>{summary.locus}</td><td>{summary.vCalls.join(", ")}</td><td>{summary.jCalls.join(", ")}</td><td className="lineage-cdr3-cell"><strong>{representative.aa||representative.nt||"—"}</strong>{representative.aa&&representative.nt&&<small>{representative.nt}</small>}{representative.lowestShm&&<small>lowest V-SHM member</small>}</td><td>{summary.cdr3Length} nt</td><td><button type="button" onClick={(event)=>{event.stopPropagation();void openLineage(summary);}}>Open →</button></td></tr>;})}</tbody></table></div></div>}
      {lineages && <section className="lineage-explorer-controls">
        <header>
          <div>
            <span className="section-kicker">Lineage explorer</span>
            <strong>{sortedLineageSummaries.length.toLocaleString()} of {lineages.summaries.length.toLocaleString()} indexed summaries shown</strong>
            <small>Filters compose with logical AND. Up to the 10,000 largest lineage summaries are retained for interactive interrogation.</small>
          </div>
          <button type="button" onClick={() => {
            setLineageSearch("");
            setLineageLocusFilter("");
            setLineageSpecificSample("");
            setLineageSampleFilterMode("any");
            setLineageDoubleDFilter("any");
            setLineageSelectedSamples(new Set());
            setLineageSampleSearch("");
            setLineageMinAbundance(0);
            setLineageMinUnique(0);
            setLineageMinSamples(0);
            setLineageMinCdr3Length(0);
            setLineageMaxCdr3Length(0);
            setLineageMinimumShm(0);
          }}>Clear lineage filters</button>
        </header>
        <div className="lineage-filter-grid">
          <label><span>Lineage, V, J, locus, or study group</span><CommitTextInput type="search" value={lineageSearch} onCommit={setLineageSearch} placeholder="contains…"/></label>
          <label><span>Locus</span><select value={lineageLocusFilter} onChange={(event) => setLineageLocusFilter(event.target.value)}><option value="">Any locus</option>{lineageLocusOptions.map((locus) => <option value={locus} key={locus}>{locus}</option>)}</select></label>
          <label><span>Contains reads from sample</span><select value={lineageSpecificSample} onChange={(event)=>setLineageSpecificSample(event.target.value)}><option value="">Any sample</option>{lineageSampleOptions.map((sample)=><option value={sample} key={sample}>{sample}</option>)}</select><small>Other samples may also occur in the lineage.</small></label>
          <label><span>Sample membership</span><select value={lineageSampleFilterMode} onChange={(event) => setLineageSampleFilterMode(event.target.value as LineageSampleFilterMode)}><option value="any">Any sample pattern</option><option value="multiple">Multiple samples</option><option value="single">One sample only</option><option value="selected-any">Any selected sample</option><option value="selected-all">Every selected sample</option><option value="selected-only">No samples outside selection</option></select></label>
          {doubleDCount > 0 && <label><span>Double-D membership</span><select value={lineageDoubleDFilter} onChange={(event) => setLineageDoubleDFilter(event.target.value as LineageDoubleDFilter)}><option value="any">Any Double-D status</option><option value="present">At least one positive member</option><option value="all">Every active member positive</option><option value="absent">No positive members</option></select><small>Supported VDDJ calls among active lineage representatives.</small></label>}
          <label><span>Minimum abundance</span><CommitNumberInput min="0" value={lineageMinAbundance} blankWhenZero placeholder="Any" onCommit={setLineageMinAbundance}/></label>
          <label><span>Minimum unique sequences</span><CommitNumberInput min="0" value={lineageMinUnique} blankWhenZero placeholder="Any" onCommit={setLineageMinUnique}/></label>
          <label><span>Minimum sample breadth</span><CommitNumberInput min="0" value={lineageMinSamples} blankWhenZero placeholder="Any" onCommit={setLineageMinSamples}/></label>
          <label><span>CDR3 nt minimum length</span><CommitNumberInput min="0" value={lineageMinCdr3Length} blankWhenZero placeholder="Any" onCommit={setLineageMinCdr3Length}/></label>
          <label><span>CDR3 nt maximum length</span><CommitNumberInput min="0" value={lineageMaxCdr3Length} blankWhenZero placeholder="Any" onCommit={setLineageMaxCdr3Length}/></label>
          {shmDashboard && <>
            <label><span>SHM filter statistic</span><select value={lineageShmFilterStatistic} onChange={(event) => setLineageShmFilterStatistic(event.target.value as LineageShmFilterStatistic)}><option value="mean">Weighted mean</option><option value="maximum">Maximum</option><option value="p95">Upper 95% quantile</option></select></label>
            <label><span>Minimum SHM value</span><CommitNumberInput min="0" step="0.001" value={lineageMinimumShm} blankWhenZero placeholder="Any" onCommit={setLineageMinimumShm}/><small>{shmMetric.toLowerCase().includes("rate") ? "Enter a fraction; 0.05 = 5%." : "Mutation count"}</small></label>
          </>}
        </div>
        {lineageSampleOptions.length > 1 && <fieldset className="lineage-sample-picker"><legend>Specific samples</legend><div className="lineage-sample-search"><CommitTextInput type="search" value={lineageSampleSearch} onCommit={setLineageSampleSearch} placeholder={`Search ${lineageSampleOptions.length.toLocaleString()} samples`}/><small>{lineageSelectedSamples.size.toLocaleString()} selected · showing {visibleLineageSampleOptions.length.toLocaleString()}</small></div>{visibleLineageSampleOptions.map((sample) => <label key={sample}><input type="checkbox" checked={lineageSelectedSamples.has(sample)} onChange={(event) => setLineageSelectedSamples((current) => { const next = new Set(current); if (event.target.checked) next.add(sample); else next.delete(sample); return next; })}/><i style={{background: sampleColor(sample, sampleColors)}}/><span>{sample}</span></label>)}</fieldset>}
        <div className="lineage-multi-open-bar"><span><strong>{lineageTableSelectedIds.size.toLocaleString()}</strong> lineage{lineageTableSelectedIds.size===1?"":"s"} selected</span><button className="post-primary dark" type="button" disabled={!lineageTableSelectedIds.size||Boolean(busy)} onClick={()=>void openSelectedLineagesTogether()}>Open together</button>{lineageTableSelectedIds.size>0&&<button type="button" onClick={()=>setLineageTableSelectedIds(new Set())}>Clear</button>}</div>
        <div className="lineage-table-wrap">
          <table>
            <thead><tr>
              <th>Select</th>
              <th><button type="button" onClick={() => sortLineages("id")}>Lineage{sortIndicator("id")}</button></th>
              <th>Merged</th>
              <th><button type="button" onClick={() => sortLineages("studyGroup")}>{DATASET_SCOPE_LABELS[lineageScope]}{sortIndicator("studyGroup")}</button></th>
              <th><button type="button" onClick={() => sortLineages("abundance")}>Abundance{sortIndicator("abundance")}</button></th>
              <th><button type="button" onClick={() => sortLineages("uniqueMembers")}>Unique{sortIndicator("uniqueMembers")}</button></th>
              {doubleDCount > 0 && <th><button type="button" onClick={() => sortLineages("doubleDPositiveMembers")}>Double-D{sortIndicator("doubleDPositiveMembers")}</button></th>}
              <th><button type="button" onClick={() => sortLineages("sampleCount")}>Samples{sortIndicator("sampleCount")}</button></th>
              <th><button type="button" onClick={() => sortLineages("locus")}>Locus{sortIndicator("locus")}</button></th>
              <th><button type="button" onClick={() => sortLineages("vCalls")}>V calls{sortIndicator("vCalls")}</button></th>
              <th><button type="button" onClick={() => sortLineages("jCalls")}>J calls{sortIndicator("jCalls")}</button></th>
              <th>Representative CDR3</th>
              <th><button type="button" onClick={() => sortLineages("cdr3Length")}>CDR3 nt{sortIndicator("cdr3Length")}</button></th>
              {shmDashboard && <><th><button type="button" onClick={() => sortLineages("shmMean")}>Mean SHM{sortIndicator("shmMean")}</button></th><th><button type="button" onClick={() => sortLineages("shmP95")}>SHM upper q95{sortIndicator("shmP95")}</button></th></>}
              <th/>
            </tr></thead>
            <tbody>{sortedLineageSummaries.slice(0, 500).map((summary) => {const representative=lineageRepresentativeCdr3(summary);return <tr key={summary.id} className={`${selectedLineageIds.includes(summary.id)?"selected":""}${lineageTableSelectedIds.has(summary.id)?" multi-selected":""}`.trim()} onClick={() => void openLineage(summary)}>
              <td onClick={(event)=>event.stopPropagation()}><input aria-label={`Select lineage ${summary.id}`} type="checkbox" checked={lineageTableSelectedIds.has(summary.id)} onChange={(event)=>setLineageTableSelected(summary.id,event.target.checked)}/></td>
              <td><strong>{summary.id}</strong></td>
              <td>{mergedIdByOriginal.get(summary.id) || "—"}</td>
              <td>{summary.studyGroup}</td>
              <td>{summary.abundance.toLocaleString()}</td>
              <td>{summary.uniqueMembers.toLocaleString()}</td>
              {doubleDCount > 0 && <td><strong>{(summary.doubleDPositiveMembers ?? 0).toLocaleString()} / {summary.uniqueMembers.toLocaleString()}</strong><small>{(summary.doubleDPositiveAbundance ?? 0).toLocaleString()} weighted abundance</small></td>}
              <td className="lineage-sample-counts"><strong>{(summary.sampleIds ?? []).length} sample{(summary.sampleIds ?? []).length===1?"":"s"}</strong>{(summary.sampleCounts??[]).map((sample)=><small key={sample.sampleId} title={`${sample.uniqueMembers.toLocaleString()} unique sequence${sample.uniqueMembers===1?"":"s"}`}>{sample.sampleId}: {sample.abundance.toLocaleString()} read{sample.abundance===1?"":"s"}</small>)}</td>
              <td>{summary.locus}</td>
              <td>{summary.vCalls.join(", ")}</td>
              <td>{summary.jCalls.join(", ")}</td>
              <td className="lineage-cdr3-cell"><strong>{representative.aa||representative.nt||"—"}</strong>{representative.aa&&representative.nt&&<small>{representative.nt}</small>}{representative.lowestShm&&<small>lowest V-SHM member</small>}</td>
              <td>{summary.cdr3Length} nt</td>
              {shmDashboard && <><td>{lineageShmLabel(summary.id, "mean")}</td><td>{lineageShmLabel(summary.id, "p95")}</td></>}
              <td><button type="button" onClick={(event)=>{event.stopPropagation();void openLineage(summary);}}>Open →</button></td>
            </tr>;})}</tbody>
          </table>
          {!sortedLineageSummaries.length && <div className="empty-results"><span>∅</span><h3>No indexed lineage matches these filters.</h3></div>}
        </div>
      </section>}
    </section>

    <section className={moduleClass("diagnostics","post-module downstream-viz-module")}>
      <header><div className="module-number amber">06</div><div><span className="section-kicker">Post-lineage analysis</span><h3>Somatic hypermutation and reference-set evidence</h3><p>These are separate downstream analyses of the current population. SHM summarizes mutation distributions; the missing-V screen tests linked haplotypes; personalized inference fits an expressed V set with one observation per lineage.</p></div><a href="./methods/09_SHM_AND_GERMLINE_DIAGNOSTICS.md" target="_blank" rel="noreferrer">Method details ↗</a>{skipStepButton("diagnostics")}<button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("diagnostics")} onClick={()=>toggleModule("diagnostics")}>{openModules.has("diagnostics")?"Collapse ↑":"Expand ↓"}</button></header>
      <div className="diagnostic-grid">
        <article><span className="section-kicker">SHM</span><h4>Mutation distributions</h4><div className="control-grid two"><label><span>Measure</span><select value={shmMetric} onChange={(event)=>{setShmMetric(event.target.value as ShmMetricKey);setShmDashboard(null);}}><option value="vNtRate">V nucleotide rate</option><option value="vNtMutations">V nucleotide count</option><option value="vAaRate">V amino-acid replacement rate</option><option value="vAaReplacements">V amino-acid replacement count</option><option value="synonymous">Synonymous codon count</option><option value="cdrNtRate">CDR1/2 nucleotide rate</option><option value="frameworkNtRate">Framework nucleotide rate</option></select></label><label><span>Stratify by</span><select value={shmStratum} onChange={(event)=>{setShmStratum(event.target.value as typeof shmStratum);setShmDashboard(null);}}><option value="all">No additional stratum</option><option value="sample_id">Sample</option><option value="subject_id">Donor / subject</option><option value="swig_cohort">Cohort</option><option value="swig_timepoint">Timepoint</option><option value="swig_compartment">Compartment / tissue</option><option value="locus">Locus</option><option value="v_call">V call</option><option value="isotype">Isotype</option></select></label></div><details className="post-advanced"><summary>Advanced plot sampling</summary><div className="control-grid"><label><span>Plot sample / lineage</span><CommitNumberInput min="50" max="10000" step="50" value={shmSampleCap} onCommit={(value)=>{setShmSampleCap(value);setShmDashboard(null);}} /></label></div></details><button className={guidedClass("run-shm")} type="button" disabled={Boolean(busy)} onClick={()=>void runShmAnalysis()}>Calculate SHM on {workingCount.toLocaleString()} records</button><p>Plot memory is bounded per lineage; scalar counts still cover every selected row.</p></article>
        <article><span className="section-kicker">Reference warning</span><h4>Possible missing V alleles</h4>
          <details className="post-advanced"><summary>Advanced missing-allele evidence thresholds</summary><div className="control-grid three">
            <div className="fixed-method-value"><span>Independent unit</span><strong>Assigned lineage</strong><small>Method invariant; record-level evidence is not permitted.</small></div>
            <label><span>Minimum independent support</span><CommitNumberInput min="3" max="1000" value={missingAlleleOptions.minimumIndependentUnits} onCommit={(minimumIndependentUnits)=>{setMissingAlleleOptions(value=>({...value,minimumIndependentUnits}));setMissingAlleles(null);}} /></label>
            <label><span>Minimum jointly covered lineages</span><CommitNumberInput min="5" max="100000" value={missingAlleleOptions.minimumCoveredUnits} onCommit={(minimumCoveredUnits)=>{setMissingAlleleOptions(value=>({...value,minimumCoveredUnits}));setMissingAlleles(null);}} /></label>
            <label><span>Minimum candidate fraction</span><CommitNumberInput min="0.01" max="1" step="0.01" value={missingAlleleOptions.minimumAlleleFraction} onCommit={(minimumAlleleFraction)=>{setMissingAlleleOptions(value=>({...value,minimumAlleleFraction}));setMissingAlleles(null);}} /></label>
            <label><span>Maximum discovery SHM</span><CommitNumberInput min="0" max="0.5" step="0.01" value={missingAlleleOptions.maximumShmRate} onCommit={(maximumShmRate)=>{setMissingAlleleOptions(value=>({...value,maximumShmRate}));setMissingAlleles(null);}} /></label>
            <label><span>Maximum linked SNPs</span><CommitNumberInput min="1" max="20" value={missingAlleleOptions.maximumCandidateSnps} onCommit={(maximumCandidateSnps)=>{setMissingAlleleOptions(value=>({...value,maximumCandidateSnps}));setMissingAlleles(null);}} /></label>
          </div><div className="control-grid three">
            <label><span>Minimum aligned V bases</span><CommitNumberInput min="30" max="1000" value={missingAlleleOptions.minimumAlignedBases} onCommit={(minimumAlignedBases)=>{setMissingAlleleOptions(value=>({...value,minimumAlignedBases}));setMissingAlleles(null);}} /></label>
            <label><span>Minimum linked fraction</span><CommitNumberInput min="0.5" max="1" step="0.01" value={missingAlleleOptions.minimumLinkedFraction} onCommit={(minimumLinkedFraction)=>{setMissingAlleleOptions(value=>({...value,minimumLinkedFraction}));setMissingAlleles(null);}} /></label>
            <label><span>Maximum screening-tail probability</span><CommitNumberInput min="0.000000000001" max="0.1" step="0.000001" value={missingAlleleOptions.maximumPValue} onCommit={(maximumPValue)=>{setMissingAlleleOptions(value=>({...value,maximumPValue}));setMissingAlleles(null);}} /></label>
            <label><span>Minimum distinct J calls</span><CommitNumberInput min="1" max="100" value={missingAlleleOptions.minimumDistinctJCalls} onCommit={(minimumDistinctJCalls)=>{setMissingAlleleOptions(value=>({...value,minimumDistinctJCalls}));setMissingAlleles(null);}} /></label>
            <label><span>Minimum distinct CDR3 lengths</span><CommitNumberInput min="1" max="100" value={missingAlleleOptions.minimumDistinctJunctionLengths} onCommit={(minimumDistinctJunctionLengths)=>{setMissingAlleleOptions(value=>({...value,minimumDistinctJunctionLengths}));setMissingAlleles(null);}} /></label>
            <label><span>Minimum distinct CDR3 sequences</span><CommitNumberInput min="1" max="10000" value={missingAlleleOptions.minimumDistinctCdr3s} onCommit={(minimumDistinctCdr3s)=>{setMissingAlleleOptions(value=>({...value,minimumDistinctCdr3s}));setMissingAlleles(null);}} /></label>
            <label><span>Near-germline supporting lineages</span><CommitNumberInput min="1" max="1000" value={missingAlleleOptions.minimumNearGermlineUnits} onCommit={(minimumNearGermlineUnits)=>{setMissingAlleleOptions(value=>({...value,minimumNearGermlineUnits}));setMissingAlleles(null);}} /><small>At most two other V mismatches in the discovery representative.</small></label>
            <label><span>Maximum other-alternate fraction</span><CommitNumberInput min="0" max="1" step="0.01" value={missingAlleleOptions.maximumOtherAlternateFraction} onCommit={(maximumOtherAlternateFraction)=>{setMissingAlleleOptions(value=>({...value,maximumOtherAlternateFraction}));setMissingAlleles(null);}} /><small>Third nucleotide states among covered lineages.</small></label>
          </div></details>
          <div className="algorithm-note"><strong>Two-pass linked-haplotype rule</strong><span>Pass 1 uses one lowest-SHM representative per lineage to propose nucleotide sets that co-occur on the same molecule. Pass 2 rescans every retained member: any covered parent-reference base at any proposed site vetoes that lineage. Candidates must also span distinct CDR3 sequences, CDR3 lengths, and J calls.</span></div>
          <button className="post-primary amber" type="button" disabled={Boolean(busy)} onClick={()=>void runMissingAlleleAnalysis()}>Run two-pass screen on {workingCount.toLocaleString()} records</button><p>This produces referral candidates only; it never edits the selected database.</p>
        </article>
        <article><span className="section-kicker">Personalized reference</span><h4>Infer expressed V allele set</h4>
          <details className="post-advanced"><summary>Advanced personalized-germline settings</summary><div className="control-grid three">
            <div className="fixed-method-value"><span>Lineage representative</span><strong>Lowest current V-SHM</strong><small>Directly from the current V alignment; no circularity correction.</small></div>
            <label><span>Minimum aligned V bases</span><CommitNumberInput min="30" max="1000" value={personalizedGermlineOptions.minimumAlignedBases} onCommit={(minimumAlignedBases)=>{setPersonalizedGermlineOptions(value=>({...value,minimumAlignedBases}));setPersonalizedGermline(null);}} /></label>
            <label><span>Known-allele SNP radius</span><CommitNumberInput min="0" max="30" value={personalizedGermlineOptions.maximumKnownAlleleSnps} onCommit={(maximumKnownAlleleSnps)=>{setPersonalizedGermlineOptions(value=>({...value,maximumKnownAlleleSnps:Math.floor(maximumKnownAlleleSnps)}));setPersonalizedGermline(null);}} /></label>
            <label><span>Maximum novel SNPs</span><CommitNumberInput min="1" max="20" value={personalizedGermlineOptions.maximumNovelSnps} onCommit={(maximumNovelSnps)=>{setPersonalizedGermlineOptions(value=>({...value,maximumNovelSnps:Math.floor(maximumNovelSnps)}));setPersonalizedGermline(null);}} /></label>
            <label><span>Novel candidate budget / gene</span><CommitNumberInput min="1" max="4096" value={personalizedGermlineOptions.maximumNovelCandidatesPerGene} onCommit={(maximumNovelCandidatesPerGene)=>{setPersonalizedGermlineOptions(value=>({...value,maximumNovelCandidatesPerGene:Math.floor(maximumNovelCandidatesPerGene)}));setPersonalizedGermline(null);}} /><small>Increase after a candidate-budget warning; larger searches take longer.</small></label>
            <label><span>Minimum novel support</span><CommitNumberInput min="2" max="10000" value={personalizedGermlineOptions.minimumNovelSupport} onCommit={(minimumNovelSupport)=>{setPersonalizedGermlineOptions(value=>({...value,minimumNovelSupport:Math.floor(minimumNovelSupport)}));setPersonalizedGermline(null);}} /></label>
            <label><span>Minimum novel fraction</span><CommitNumberInput min="0" max="1" step="0.001" value={personalizedGermlineOptions.minimumNovelFraction} onCommit={(minimumNovelFraction)=>{setPersonalizedGermlineOptions(value=>({...value,minimumNovelFraction}));setPersonalizedGermline(null);}} /><small>Default 0: no fraction cutoff; support and model evidence still apply.</small></label>
            <label><span>Extra log-evidence gain</span><CommitNumberInput min="0" max="100" step="0.5" value={personalizedGermlineOptions.minimumLogEvidenceGain} onCommit={(minimumLogEvidenceGain)=>{setPersonalizedGermlineOptions(value=>({...value,minimumLogEvidenceGain}));setPersonalizedGermline(null);}} /><small>Required after frequency and sequence-search penalties.</small></label>
            <label><span>Sequencing error floor</span><CommitNumberInput min="0.000001" max="0.1" step="0.0001" value={personalizedGermlineOptions.sequencingErrorRate} onCommit={(sequencingErrorRate)=>{setPersonalizedGermlineOptions(value=>({...value,sequencingErrorRate}));setPersonalizedGermline(null);}} /></label>
            <label><span>Active-allele computational guard / gene</span><CommitNumberInput min="1" max="30" value={personalizedGermlineOptions.maximumActiveAllelesPerGene} onCommit={(maximumActiveAllelesPerGene)=>{setPersonalizedGermlineOptions(value=>({...value,maximumActiveAllelesPerGene:Math.floor(maximumActiveAllelesPerGene)}));setPersonalizedGermline(null);}} /></label>
          </div></details>
          <div className="algorithm-note"><strong>Lowest-SHM lineage row → calibrated SHM evidence → sparse allele set</strong><span>Each lineage contributes its lowest V-mismatch-rate member under the current assignment. Same-length, alignment-compatible alternatives compete across nearby genes. Novel haplotypes must beat a linked SHM null with an unmutated component, flexible hotspot rates and sequence-search penalties. Boundary extensions must also beat gene-specific trimming and junction-start bias. Final-two-base changes remain unresolved. Subjects never share evidence.</span></div>
          <button className="post-primary amber" type="button" disabled={Boolean(busy)} onClick={()=>void runPersonalizedGermlineAnalysis()}>Infer personalized V set from {lineages?.lineageCount.toLocaleString()??"assigned"} lineages</button><p>This creates a downloadable candidate reference; it does not change current calls.</p>
        </article>
        <article><span className="section-kicker">Independent experimental method</span><h4>Joint inherited / SHM model</h4><a href="./methods/10_JOINT_INHERITED_SHM_MODEL.md" target="_blank" rel="noreferrer">Model, evidence rule and limitations ↗</a>
          <p>Propose sequences on one half of the lineages, then test inheritance against a flexible somatic explanation on the other half. Both explanations use the same mutation and V-end likelihood. Results are separate from the existing personalized-reference method.</p>
          <details className="post-advanced"><summary>Joint-model settings</summary><div className="control-grid three">
            <label><span>Training candidate budget / component</span><CommitNumberInput min="1" max="512" value={unifiedOptions.maximumCandidates} onCommit={maximumCandidates=>{setUnifiedOptions(v=>({...v,maximumCandidates:Math.floor(maximumCandidates)}));setUnifiedGermline(null);}} /></label>
            <label><span>Maximum fitting iterations</span><CommitNumberInput min="50" max="2000" value={unifiedOptions.maximumIterations} onCommit={maximumIterations=>{setUnifiedOptions(v=>({...v,maximumIterations:Math.floor(maximumIterations)}));setUnifiedGermline(null);}} /><small>Unfinished optimization lowers reported evidence conservatively.</small></label>
          </div></details>
          <div className="algorithm-note"><strong>Experimental evidence comparison</strong><span>This model can abstain when inherited and somatic explanations are indistinguishable. Final two V bases remain unresolved. Scores are not calibrated real-data false-positive probabilities. Large competitive components can take several minutes.</span></div>
          <button className="post-primary amber" type="button" disabled={Boolean(busy)} onClick={()=>void runUnifiedGermlineAnalysis()}>Run joint inherited / SHM model</button>
        </article>
      </div>
      {shmDashboard?<ShmResultsPanel dashboard={shmDashboard} name={baseName(inputName)} color={chartColor} sampleColors={sampleColors} stratum={shmStratum} datasets={datasets} sampleOrder={shmSampleOrder} onSampleOrderChange={setShmSampleOrder}/>:null}
      {missingAlleles?<MissingAlleleResultsPanel dashboard={missingAlleles} name={baseName(inputName)} referenceFasta={references.V} selectedIds={selectedMissingAlleleIds} onSelectedIdsChange={setSelectedMissingAlleleIds}/>:null}
      {unifiedGermline?<UnifiedGermlineResultsPanel dashboard={unifiedGermline} name={baseName(inputName)}/>:null}
      {personalizedGermline?<PersonalizedGermlineResultsPanel dashboard={personalizedGermline} name={baseName(inputName)} referenceFasta={references.V}/>:null}
    </section>

    <section className={moduleClass("query","post-module query-module")}>
      <header><div className="module-number dark">08</div><div><span className="section-kicker">Targeted repertoire search</span><h3>Query sequences, then expand the matched set</h3><p>Paste one or more sequences or FASTA records. V/J constraints can be supplied directly or inferred per seed with the same references, assignment strategy, and calling profile as the main analysis.</p></div><a href="./methods/10_RETRIEVAL_AND_NEIGHBOURS.md" target="_blank" rel="noreferrer">Method details ↗</a>{skipStepButton("query")}<button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("query")} onClick={()=>toggleModule("query")}>{openModules.has("query")?"Collapse ↑":"Expand ↓"}</button></header>
      <div className="query-layout">
        <div className="query-input">
          <label><span>Query sequence(s) or FASTA</span><CommitTextarea value={queryText} onCommit={setQueryText} placeholder=">seed_1\nTGTGCGAGAGAT…\n>seed_2\nTGTGCGAGGGAT…" /></label>
          <div className="control-grid three">
            <label><span>Target</span><select value={queryTarget} onChange={(event) => setQueryTarget(event.target.value as QueryTarget)}><option value="cdr3_nt">CDR3 nucleotide</option><option value="cdr3_aa">CDR3 amino acid</option><option value="trimmed">VDJ-aligned k-mer sketch</option></select></label>
            <label><span>Return mode</span><select value={queryResultMode} onChange={(event)=>{setQueryResultMode(event.target.value as "sequences"|"lineages");setQueryLineageHits([]);setExpanded(null);}}><option value="sequences">Matching sequences</option><option value="lineages" disabled={!lineages}>Best-matching lineages</option></select><small>{lineages?"Lineages rank by their best member hit.":"Assign lineages to enable lineage mode."}</small></label>
            <label><span>Metric</span><select disabled={queryTarget === "trimmed"} value={queryTarget === "trimmed" ? "sketch" : queryMetric} onChange={(event) => setQueryMetric(event.target.value as QueryMetric)}><option value="exact">Exact</option><option value="substring">Substring</option><option value="hamming">Hamming</option><option value="edit">Banded edit distance</option><option value="sketch">MinHash k-mer estimate</option></select></label>
            <label><span>Identity / similarity</span><CommitNumberInput min="0" max="1" step="0.01" value={queryIdentity} onCommit={setQueryIdentity} /></label>
            <label><span>V/J constraint source</span><select value={queryConstraintMode} onChange={(event) => { setQueryConstraintMode(event.target.value as "manual" | "infer"); setQueryInference([]); }}><option value="manual">Manual / unconstrained</option><option value="infer">Infer per query with the aligner</option></select></label>
            <label><span>{queryConstraintMode === "infer" ? "Locus override" : "Locus constraint"}</span><select value={queryLocus} onChange={(event) => setQueryLocus(event.target.value)}><option value="">{queryConstraintMode === "infer" ? "Use each inferred locus" : "Any locus"}</option>{loci.map((item) => <option value={item.value} key={item.value}>{item.value}</option>)}</select></label>
            <label><span>{queryConstraintMode === "infer" ? "V override" : "V constraint"}</span><CommitTextInput value={queryV} onCommit={setQueryV} placeholder={queryConstraintMode === "infer" ? "otherwise inferred per seed" : "optional IGHV…"} /></label>
            <label><span>{queryConstraintMode === "infer" ? "J override" : "J constraint"}</span><CommitTextInput value={queryJ} onCommit={setQueryJ} placeholder={queryConstraintMode === "infer" ? "otherwise inferred per seed" : "optional IGHJ…"} /></label>
            <label><span>Maximum initial hits</span><CommitNumberInput min="1" max="10000" value={queryLimit} onCommit={setQueryLimit} /></label>
          </div>
          {queryInference.length > 0 && <div className="query-inference"><header><strong>Aligner seed assignments</strong><span>{queryInference.filter((item) => item.assigned).length} / {queryInference.length} with both V and J</span></header>{queryInference.slice(0, 20).map((assignment) => <div className={assignment.assigned ? "assigned" : "unassigned"} key={assignment.queryIndex}><b>Seed {assignment.queryIndex + 1}</b><span>{assignment.locus || "locus —"}</span><code>{assignment.vCall || "V —"}</code><code>{assignment.jCall || "J —"}</code><small>{assignment.searchSequence ? `${assignment.searchSequence.length} ${queryTarget === "cdr3_aa" ? "aa" : "nt"} search target` : "no target derived"}</small></div>)}</div>}
          <div className="current-step-input"><span>Search universe inherited from applied filters</span><strong>{workingCount.toLocaleString()} active records</strong><small>Initial search and every single-linkage expansion both use this same working set.</small></div>
          <div className="result-actions">{queryConstraintMode === "infer" && <button type="button" disabled={Boolean(busy)} onClick={() => void previewQueryInference()}>Preview inferred V/J</button>}<button className="post-primary dark" type="button" disabled={Boolean(busy)} onClick={() => void runQuery()}>{queryConstraintMode === "infer" ? `Assign seeds + rank ${queryResultMode}` : `Search + rank ${queryResultMode}`}</button><button type="button" disabled={Boolean(busy) || !queryHits.length} onClick={() => void expandMatches()}>Single-linkage expand sequence set</button></div>
          <p className="scientific-note"><span>i</span>{queryConstraintMode === "infer" ? `Input is interpreted as rearranged nucleotide sequence. The Swig aligner uses this run’s composed references, ${callingProfile === "igblast_compatible" ? "IgBLAST-agreement" : callingProfile === "r_optimized" ? "R-optimized" : callingProfile === "sensitive_d" ? "Sensitive-D" : "truth-optimized"} calling profile, ${Math.round(minimumIdentity * 100)}% identity floor, and ${strand === 0 ? "both strands" : strand === 1 ? "plus strand" : "minus strand"}; ambiguous V/J calls remain ambiguity-aware and are applied per seed. Non-empty override fields replace the inferred value.` : queryTarget === "trimmed" ? "VDJ searches lazily build a packed index of eight independent 7-mer MinHash values per record. Scores are approximate and intended for candidate retrieval." : "Hamming search requires equal length. Edit distance is banded by the selected identity. Expansion always uses equal-length CDR3 nucleotide Hamming edges."}</p>
        </div>
        <div className="query-results"><header><div><span className="section-kicker">Matched set</span><h4>{queryResultMode==="lineages"?`${queryLineageHits.length.toLocaleString()} matching lineages`:expanded ? `${expanded.ordinals.length.toLocaleString()} expanded records` : `${queryHits.length.toLocaleString()} initial hits`}</h4><p>{queryResultMode==="lineages"?"Each lineage score is its best match to any member sequence.":expanded ? `${expanded.comparisons.toLocaleString()} exact edge checks${expanded.capped ? " · result cap reached" : " · fixed point reached"}` : "Ranked by similarity, then distance."}</p></div></header><div className="query-result-list">{queryResultMode==="lineages"?queryLineageHits.slice(0,100).map((hit)=><button type="button" key={hit.lineageId} onClick={()=>void openQueryLineage(hit)}><span><strong>Lineage {hit.lineageId}</strong><small>{hit.matchedSequences.toLocaleString()} matched member{hit.matchedSequences===1?"":"s"} · {hit.matchedQueries.toLocaleString()} query seed{hit.matchedQueries===1?"":"s"}</small></span><code>best AIRR row #{(hit.bestOrdinal+1).toLocaleString()}</code><b>{(hit.bestScore*100).toFixed(1)}% · open →</b></button>):displayedQueryRows.slice(0, 100).map((record) => { const hit = queryHits.find((value) => value.ordinal === record.ordinal); return <button type="button" key={record.ordinal} onClick={() => onInspect(record.ordinal)}><span><strong>{record.sequenceId}</strong><small>{record.locus} · {record.vCall || "V—"} · {record.jCall || "J—"}</small></span><code>{record.cdr3 || record.cdr3Aa || "—"}</code><b>{hit ? `${(hit.score * 100).toFixed(1)}%` : "expanded"}</b></button>; })}{(queryResultMode==="lineages"?!queryLineageHits.length:!displayedQueryRows.length) && <div className="method-placeholder small"><span>⌕</span><h4>No query results</h4><p>Provide a sequence and search the assigned repertoire.</p></div>}</div></div>
      </div>
    </section>

    {selectedLineage && <section ref={workbenchRef} className={moduleClass("workbench","post-module lineage-workbench")} tabIndex={-1}>
      <header><div className="module-number dark">{directLineage ? "01" : "07"}</div><div><span className="section-kicker">{directLineage ? "Complete uploaded input · one lineage" : selectedLineageIds.length > 1 ? `Combined view · ${selectedLineageIds.length} original lineages` : `Selected lineage ${selectedLineage.id}`} · {selectedLineage.studyGroup}</span><h3>{directLineage ? "Lineage alignment, phylogeny and UCA" : "Lineage neighbours, alignment and rooted phylogeny"}</h3><p>{lineageTotal.toLocaleString()} input rows · {selectedLineage.locus}{directLineage ? " · lineage clustering was bypassed" : ` · original assignments ${selectedLineageIds.map((id) => `L${id}`).join(", ")}`}</p></div><a href="./methods/11_ALIGNMENT_AND_PHYLOGENY.md" target="_blank" rel="noreferrer">Method details ↗</a>{!directLineage && <>{skipStepButton("workbench")}<button type="button" onClick={() => { setSelectedLineage(null); setSelectedLineageIds([]); setLineageRows([]); setAlignmentProductiveOnly(false); setLineageMultiplicity(new Map()); setOriginalLineageByOrdinal(new Map()); clearNeighbourResults(); setActiveWorkspace("lineage"); }}>Close lineage</button><button className="module-collapse-toggle" type="button" aria-expanded={openModules.has("workbench")} onClick={()=>toggleModule("workbench")}>{openModules.has("workbench")?"Collapse ↑":"Expand ↓"}</button></>}</header>
      <div className="member-strip"><div><strong>{workbenchLineageRows.length.toLocaleString()}</strong><span>active workspace rows</span><small>{alignmentProductiveOnly?`${excludedNonProductiveRows.toLocaleString()} non-productive or unresolved loaded rows excluded`:lineageTotal > lineageRows.length ? `${lineageRows.length} stratified members loaded from ${lineageTotal.toLocaleString()}; analysis remains bounded` : "complete selected lineage working set"}</small></div><div className="member-pills">{workbenchLineageRows.slice(0, 18).map((row) => {const sample=row.values.sample_id||row.record.sampleId||"";const original=originalLineageByOrdinal.get(row.record.ordinal);const cdr3=row.values.cdr3_aa||row.values.cdr3||"CDR3 —";return <button type="button" key={row.record.ordinal} onClick={() => onInspect(row.record.ordinal)} title={`${sample||"sample unassigned"}${original?` · original lineage ${original}`:""} · ${cdr3}`}><i style={{background:sampleColor(sample,sampleColors)}}/><span><b>{row.record.sequenceId}</b><small>{cdr3}</small></span></button>;})}</div></div>
      {!directLineage && <div className="lineage-neighbour-explorer">
        <header><div><span className="section-kicker">Exploratory boundary review</span><h4>Lineage neighbours</h4><p>Find separate assigned lineages with CDR3 links below the clustering cutoff, similar inferred lineage germlines, or either criterion. Search is read-only until you explicitly merge.</p></div><span className="neighbour-source-chip">Source {selectedLineageIds.map((id)=>`L${id}`).join(" + ")}</span></header>
        <div className="control-grid five">
          <label><span>Evidence route</span><select value={neighbourMethod} onChange={(event)=>{setNeighbourMethod(event.target.value as typeof neighbourMethod);clearNeighbourResults();}}><option value="either">CDR3 or inferred germline</option><option value="cdr3">CDR3 only</option><option value="germline">Inferred germline only</option></select></label>
          <label><span>Search boundary</span><select value={neighbourScope} onChange={(event)=>{setNeighbourScope(event.target.value as DatasetScope);clearNeighbourResults(true);}}>{Object.entries(DATASET_SCOPE_LABELS).map(([value,label])=><option value={value} key={value}>{label}</option>)}</select></label>
          <label><span>CDR3 distance</span><select value={neighbourCdr3Metric} onChange={(event)=>{setNeighbourCdr3Metric(event.target.value as typeof neighbourCdr3Metric);clearNeighbourResults();}}><option value="hamming">Equal-length Hamming · assignment-compatible</option><option value="edit">Indel-aware edit distance · length variation</option></select></label>
          <label><span>V/J evidence gate</span><select value={neighbourCallPolicy} onChange={(event)=>{setNeighbourCallPolicy(event.target.value as typeof neighbourCallPolicy);clearNeighbourResults();}}><option value="both">Require compatible V and J · assignment-compatible</option><option value="either">Require compatible V or J · noisy-call tolerant</option><option value="ignore">Ignore V/J · same locus and boundary only</option></select></label>
          <label><span>Minimum CDR3 identity</span><CommitNumberInput min="0.5" max="1" step="0.01" value={neighbourCdr3Identity} onCommit={setNeighbourCdr3Identity}/><small>Assignment cutoff: {(identity*100).toFixed(0)}%</small></label>
          <label><span>Minimum inferred germline identity</span><CommitNumberInput min="0.5" max="1" step="0.01" value={neighbourGermlineIdentity} onCommit={setNeighbourGermlineIdentity}/></label>
          <label><span>Maximum displayed</span><CommitNumberInput min="1" max="250" value={neighbourLimit} onCommit={setNeighbourLimit}/></label>
        </div>
        <div className="algorithm-note"><strong>Indexed shortlist → exact verification</strong><span>{neighbourCdr3Metric==="edit"?"Indel-aware CDR3 mode uses a study/locus-partitioned BK tree and exact banded Levenshtein verification, so insertions and deletions—and therefore different CDR3 lengths—can match.":"Hamming mode uses the same exact d+1 block index and equal-length verification as lineage assignment."} {neighbourCallPolicy==="both"?"Both V and J must remain compatible.":neighbourCallPolicy==="either"?"Either V or J compatibility is sufficient, making the search tolerant of one noisy germline call.":"V/J calls do not gate candidates; locus and the selected study boundary remain mandatory."} Germline mode independently reconstructs shortlisted lineage germlines with the selected {lineageGermlineMethod === "closest" ? "closest-member" : "equal-weight consensus"} method. No route modifies lineage IDs.</span></div>
        <div className="result-actions"><button className="post-primary dark" type="button" disabled={Boolean(busy)} onClick={()=>void searchLineageNeighbours()}>Search neighbouring lineages</button>{combinedNeighbourHits.length>0&&<><button type="button" disabled={!selectedNeighbourIds.size} onClick={()=>void viewSelectedNeighbourGroup(false)}>View selected together</button><button type="button" disabled={!selectedNeighbourIds.size} onClick={()=>void viewSelectedNeighbourGroup(true)}>Merge + view together</button></>}</div>
        {(neighbourResult||germlineNeighbourScores.length>0)&&<><div className="neighbour-search-stats"><span>{combinedNeighbourHits.length.toLocaleString()} displayed candidates</span>{neighbourResult&&<span>{neighbourResult.candidateComparisons.toLocaleString()} exact CDR3 comparisons</span>}{germlineSketchIndex&&<span>{germlineSketchIndex.representedLineages.toLocaleString()} lineage germline sketches</span>}<span>{selectedNeighbourIds.size.toLocaleString()} selected</span></div><div className="lineage-table-wrap neighbour-table"><table><thead><tr><th>Select</th><th>Original lineage</th><th>Best source</th><th>CDR3 identity</th><th>Edit / length Δ</th><th>V/J agreement</th><th>Inferred germline identity</th><th>Members / abundance</th><th>Boundary</th></tr></thead><tbody>{combinedNeighbourHits.map((hit)=><tr key={hit.lineageId} className={selectedNeighbourIds.has(hit.lineageId)?"selected":""}><td><input aria-label={`Select lineage ${hit.lineageId}`} type="checkbox" checked={selectedNeighbourIds.has(hit.lineageId)} onChange={(event)=>setNeighbourSelected(hit.lineageId,event.target.checked)}/></td><td><strong>L{hit.lineageId}</strong>{mergedIdByOriginal.get(hit.lineageId)&&<small>{mergedIdByOriginal.get(hit.lineageId)}</small>}</td><td>L{hit.cdr3?.sourceLineageId??hit.germline?.sourceLineageId??"—"}</td><td>{hit.cdr3?`${(hit.cdr3.cdr3Identity*100).toFixed(2)}%`:"—"}</td><td>{hit.cdr3?`${hit.cdr3.cdr3Distance??0} ${hit.cdr3.cdr3Metric==="edit"?"edit":"mismatch"}${hit.cdr3.cdr3LengthDelta?` · ${hit.cdr3.cdr3LengthDelta>0?"+":""}${hit.cdr3.cdr3LengthDelta} nt`:" · same length"}`:"—"}</td><td>{hit.cdr3?.callAgreement==="both"?"V + J":hit.cdr3?.callAgreement==="v"?"V only":hit.cdr3?.callAgreement==="j"?"J only":hit.cdr3?"none":"—"}</td><td>{hit.germline?`${(hit.germline.germlineIdentity*100).toFixed(2)}%`:"—"}</td><td>{hit.cdr3?`${hit.cdr3.uniqueMembers.toLocaleString()} / ${hit.cdr3.abundance.toLocaleString()}`:hit.germline?`${hit.germline.candidateTotalRows.toLocaleString()} / —`:"—"}</td><td>{hit.cdr3?.studyGroup||"same selected boundary"}</td></tr>)}</tbody></table></div></>}
        {(neighbourResult&&neighbourResult.hits.length===0&&!germlineNeighbourScores.length)&&<div className="method-placeholder small"><span>∅</span><h4>No candidate met the selected neighbour criteria</h4><p>Lower the exploratory threshold or broaden the search boundary; original assignments remain unchanged.</p></div>}
      </div>}
      <div className={`alignment-population-filter${alignmentProductiveOnly?" active":""}`}>
        <div><span className="section-kicker">Alignment population</span><strong>{workbenchLineageRows.length.toLocaleString()} of {lineageRows.length.toLocaleString()} loaded rows included</strong><small>{alignmentProductiveOnly?"Only AIRR rows explicitly marked productive are used for germline construction, MSA, tree metadata, and UCA inference.":"All loaded rows currently contribute to the lineage workspace."} The repertoire and lineage assignment are never modified.</small></div>
        <button type="button" disabled={Boolean(busy)||(!alignmentProductiveOnly&&excludedNonProductiveRows===0)} onClick={toggleAlignmentProductiveFilter}>{alignmentProductiveOnly?"Restore all loaded rows":excludedNonProductiveRows?`Drop non-productive (${excludedNonProductiveRows.toLocaleString()})`:"All loaded rows productive"}</button>
      </div>
      {lineageGermline&&<div className="lineage-germline-method">
        <div>
          <span className="section-kicker">Lineage root construction</span>
          <label><span>Germline / UCA method</span><select value={lineageGermlineMethod} onChange={(event)=>{const method=event.target.value as LineageGermlineMethod;setLineageGermlineMethod(method);clearAlignmentArtifacts();clearNeighbourResults(true);}}><option value="closest">Closest member by matched V + J identity · default</option><option value="consensus">Equal-weight member consensus · alternative</option></select></label>
          {lineageGermlineMethod === "closest" && <label><span>Trimming template</span><select aria-label="Lineage germline trimming template" disabled={Boolean(busy)} value={templateOrdinal ?? "auto"} onChange={(event) => chooseLineageTemplate(event.target.value)}><option value="auto">Automatic · complete alignments preferred</option>{templateChoices.map((choice) => <option key={choice.ordinal} value={choice.ordinal}>{choice.sequenceId} · row {choice.ordinal + 1} · {choice.complete ? "complete" : "partial / unverified"}{choice.identity >= 0 ? ` · ${(choice.identity * 100).toFixed(2)}% V/J` : ""}{choice.doubleDResolved ? " · VDDJ" : ""}</option>)}</select></label>}
          <strong>{lineageGermlineMethod==="closest"?`Template member: ${lineageGermline.selectedSequenceId||`AIRR row ${(lineageGermline.selectedOrdinal??0)+1}`}`:"Equal-weight AIRR-anchored member voting"}</strong>
          <p>{lineageGermlineMethod==="closest"?"Automatic selection favors alignments reaching the V 5′ and J 3′ reference ends, then ranks by equal-weight V/J identity, combined identity, coverage, and AIRR order. The chosen member supplies junction-facing trimming; missing outer coverage is restored in the guide from its exact called references and validated AIRR coordinates. Read gaps remain missing data. A manual choice overrides the ranking and Double-D preference. With supported Double-D evidence, automatic selection prefers safely projected VDDJ members. Junction N sites remain unresolved in the tree guide.":"Each loaded representative contributes one equal vote at every AIRR-anchored V/D/J germline column, including validated outer reference flanks. For a lineage with supported Double-D evidence, only safely projected V–D1–D2–J members vote. A germline base requires ≥80% agreement; unresolved junction bases stay N in the tree root, while the comparison UCA may fill them with ≥60% unweighted query consensus."}</p>
          {lineageGermlineMethod === "closest" && <small>Outer reference bases restored: V 5′ {lineageGermline.restoredVPrefix ?? 0} · J 3′ {lineageGermline.restoredJSuffix ?? 0}. Completion requires exact called references and consistent AIRR coordinates. {alignment && <button type="button" disabled={Boolean(busy)} onClick={() => chooseLineageTemplate(templateOrdinal === undefined ? "auto" : String(templateOrdinal))}>Apply template to current MSA</button>}</small>}
          {lineageGermline.doubleDPositiveRows>0&&<div className={`double-d-root-status${lineageGermline.doubleDTemplate?" resolved":" unresolved"}`} role="status">
            <b>{lineageGermline.doubleDTemplate?"V–D1–D2–J root":"Double-D root not projected"}</b>
            <span>{lineageGermline.doubleDTemplate
              ? `${lineageGermline.doubleDRowsUsed.toLocaleString()} VDDJ-aware ${lineageGermlineMethod==="closest"?"member supplied":"members supplied"} the root${lineageGermline.selectedDCall&&lineageGermline.selectedD2Call?` · ${lineageGermline.selectedDCall} → ${lineageGermline.selectedD2Call}`:""}. D1 and D2 reference bases are retained; NP1, NP2, and NP3 remain N.`
              : `${lineageGermline.doubleDIncompleteRows.toLocaleString()} positive row(s) lacked a complete coordinate/alignment projection. The displayed root remains the baseline V(D)J template; inspect or rerun the Double-D evidence before inferring a tree.`}</span>
          </div>}
        </div>
        <div><b>{lineageGermline.rowsUsed.toLocaleString()}</b><span>{lineageGermlineMethod==="closest"?"member used":"members voted"}</span></div>
        <div><b>{lineageGermline.knownColumns.toLocaleString()}</b><span>reference-resolved columns</span></div>
        <div><b>{lineageGermline.inferredColumns.toLocaleString()}</b><span>N sites filled in comparison UCA</span></div>
        <div><b>{lineageGermlineMethod==="closest"&&lineageGermline.selectedVjIdentity!==undefined?`${(lineageGermline.selectedVjIdentity*100).toFixed(2)}%`:lineageGermline.conflictingColumns.toLocaleString()}</b><span>{lineageGermlineMethod==="closest"?"equal-weight V/J identity":"conflicting columns retained as N"}</span></div>
      </div>}
      <div className="alignment-controls"><label><span>Alignment method</span><select value={alignmentMethod} onChange={(event) => setAlignmentMethod(event.target.value as "quick" | "kalign" | "codon")}><option value="quick">Ref-anchored quick view · default</option><option value="kalign">Kalign 3.3.1 WASM · nucleotide</option><option value="codon">Kalign 3.3.1 WASM · codon-aware</option></select></label><label><span>Maximum sequences</span><CommitNumberInput min="2" max="500" value={alignmentLimit} onCommit={setAlignmentLimit} /></label><button className="post-primary" type="button" disabled={Boolean(busy)||workbenchLineageRows.length<2} onClick={() => void runAlignment()}>{alignmentMethod==="quick"?"Prepare quick view":"Align selected lineage"}</button><button type="button" title="Align nucleotides directly. N receives the full nucleotide match score; ordinary gap penalties remain in force." disabled={Boolean(busy)||workbenchLineageRows.length<2} onClick={() => void runDirectAlivibeMsa()}>MSA lineage · Alivibe WASM</button><button type="button" title="Translate each row in its AIRR-derived frame, align amino acids with X as an unknown-residue wildcard, then project every amino-acid gap back as three nucleotide gaps." disabled={Boolean(busy)||workbenchLineageRows.length<2} onClick={() => void runDirectAlivibeAaMsa()}>AA MSA · Alivibe WASM</button>{savedEditedAlignment&&!alignmentEdited&&<button type="button" onClick={()=>installAlignment(savedEditedAlignment.fasta,savedEditedAlignment.source,true,savedEditedAlignment.lineageIds,savedEditedAlignment.frameOffset,savedEditedAlignment.productiveOnly)}>Restore saved manual edit</button>}{alignment && <><label><span>AA reading frame</span><select value={alignmentFrameOffset} onChange={(event)=>changeAlignmentFrameOffset(Number(event.target.value) as AlignmentFrameOffset)}><option value="0">Start at nucleotide column 1</option><option value="1">Start at nucleotide column 2</option><option value="2">Start at nucleotide column 3</option></select></label><label><span>Alignment export</span><select value={alignmentExportFormat} onChange={(event)=>setAlignmentExportFormat(event.target.value as AlignmentExportFormat)}><option value="fasta">Aligned FASTA</option><option value="clustal">Clustal</option><option value="phylip">Relaxed PHYLIP</option><option value="stockholm">Stockholm</option><option value="nexus">NEXUS</option></select></label><button type="button" onClick={downloadCurrentAlignment}>Download alignment ↓</button></>}</div>
      {alignment && <>
        <div className={`alignment-editor-transfer${alignmentEdited?" edited-saved":""}`}><div><span className="section-kicker">Optional manual correction</span><h4>Round trip through Alivibe</h4><p>The nucleotide and AA Alivibe WASM buttons above perform their MSA entirely in Swig. Open this editor only when you want to inspect or curate that alignment manually. Gap edits, deleted alignment columns, deleted nucleotide characters, and removal of bad biological rows are accepted. Added or renamed rows and nucleotide substitutions are rejected. Keep <code>{GERMLINE_OUTGROUP}</code> so rooting remains reproducible. The bundled editor returns the complete ordered records from the same nucleotide state used by Alivibe’s NT viewer and NT export; its AA frame is transferred separately. Deleted biological rows are also removed from the AIRR metadata used by downstream UCA preparation.</p>{alignmentEdited?<strong className="session-preserved-badge">Manual alignment + AA frame · included in Save session</strong>:<small>Generated alignments are reproducible and omitted from sessions unless a tree needs their exact input. Any manually returned/imported correction is preserved.</small>}</div><div className="result-actions"><button type="button" onClick={openAlivibeEditor}>Open + load in Alivibe ↗</button><button type="button" onClick={importFromAlivibe}>Read live Alivibe NT view</button><label className="file-button">Import corrected FASTA<input type="file" accept=".fa,.fasta,.fas,.aln,.txt" onChange={(event) => void acceptEditedAlignment(event)} /></label>{savedEditedAlignment&&<button type="button" onClick={()=>{setEditedAlignments((current)=>{const next=new Map(current);next.delete(selectedGroupKey);return next;});if(alignmentEdited)setAlignmentEdited(false);}}>Discard saved manual edit</button>}</div>{alignmentEditorStatus && <p className="editor-status">{alignmentEditorStatus}</p>}{alignmentEditorError && <div className="inline-method-error" role="alert">{alignmentEditorError}</div>}<small>The live return never reads the system clipboard and is bound to the lineage/alignment from which the editor was opened. For a downloaded corrected FASTA, verify the adjacent AA reading-frame control because FASTA itself does not encode codon phase. Sequence data remain in the browser.</small></div>
        <div className="alignment-view-controls"><div className="mode-toggle"><button className={alignmentMode === "nt" ? "active" : ""} type="button" onClick={() => setAlignmentMode("nt")}>Nucleotide</button><button className={alignmentMode === "aa" ? "active" : ""} type="button" onClick={() => setAlignmentMode("aa")}>Amino acid</button></div><span>{alignmentInfo?.rows.toLocaleString() ?? parseFasta(alignment, true).length.toLocaleString()} aligned rows · {alignmentInfo?.columns.toLocaleString() ?? "—"} columns · AA frame starts at nucleotide column {alignmentFrameOffset + 1} · {alignmentSource || "alignment"} · fingerprint {alignmentInfo?.fingerprint ?? "—"}</span></div>
        <AlignmentPreview fasta={alignment} mode={alignmentMode} frameOffset={alignmentFrameOffset} />
        <div ref={treeResultRef} className="tree-operation-region">
          <div className="tree-controls"><div><span className="section-kicker">On-demand tree</span><h4>Infer or upload a lineage tree</h4><p>FastTree uses the exact current nucleotide alignment and its synthetic germline guide. An uploaded Newick may contain that guide for ordinary display, or may be observed-only for optional reuse by UCA inference.</p></div><label><span>FastTree model</span><select value={treeModel} onChange={(event) => setTreeModel(event.target.value as "gtr" | "jc")}><option value="gtr">GTR</option><option value="jc">Jukes–Cantor</option></select></label><label className="check-line"><input type="checkbox" checked={treeFast} onChange={(event) => setTreeFast(event.target.checked)} /><span>Fastest heuristic</span></label><button className="post-primary dark" type="button" disabled={Boolean(busy)} onClick={() => void inferTree()}>{busy.startsWith("Running FastTree") ? "Inferring tree…" : "Infer with FastTree"}</button><label className="file-button tree-upload-button">Upload Newick<input type="file" accept=".nwk,.newick,.tree,.tre,.txt" onChange={(event) => void acceptUploadedTree(event)} /></label></div>
          {treeError && <div className="inline-method-error tree-error" role="alert"><strong>Tree could not be loaded or inferred</strong><span>{treeError}</span><button type="button" onClick={() => setTreeError("")}>Dismiss</button></div>}
          {treeRun && <>
            <div className="tree-provenance"><div><span className="section-kicker">{treeRun.origin === "uploaded" ? "Validated upload" : "Executed input"}</span><strong>{treeRun.rows.toLocaleString()} biological rows × {treeRun.columns.toLocaleString()} columns</strong><small>{treeRun.source} · alignment fingerprint {treeRun.fingerprint}</small><code>{treeRun.command}</code>{treeRun.origin === "uploaded" && <em className={treeRun.observedOnly ? "tree-uca-eligible" : "tree-display-only"}>{treeRun.observedOnly ? "Observed-only · eligible as an optional UCA tree" : "Contains germline guide · display only; UCA will infer a fresh observed-only tree"}</em>}</div><div className="result-actions"><button type="button" onClick={() => downloadText(treeRun.alignmentFasta, `${baseName(inputName)}.lineage-${selectedLineage.id}.tree-alignment.fasta`)}>Named alignment FASTA ↓</button>{treeRun.inputFasta && <button type="button" onClick={() => downloadText(treeRun.inputFasta, `${baseName(inputName)}.lineage-${selectedLineage.id}.fasttree-exact-input.fasta`)}>Exact numeric input ↓</button>}<button type="button" onClick={()=>downloadText(treeViewMode==="stable"?treeRun.stableNewick:treeViewMode==="rooted"?treeRun.rootedNewick:treeRun.newick,`${baseName(inputName)}.lineage-${selectedLineage.id}.${treeViewMode}.nwk`)}>Newick ↓</button><button type="button" onClick={()=>downloadText(treeNexus(treeViewMode==="stable"?treeRun.stableNewick:treeViewMode==="rooted"?treeRun.rootedNewick:treeRun.newick,`lineage_${selectedLineage.id}`),`${baseName(inputName)}.lineage-${selectedLineage.id}.${treeViewMode}.nex`)}>NEXUS ↓</button></div></div>
            <div className="tree-output-switch"><div className="mode-toggle"><button className={treeViewMode === "rooted" ? "active" : ""} type="button" onClick={() => setTreeViewMode("rooted")}>{treeRun.origin === "uploaded" ? "Uploaded · ordered" : "Rooted · resolved"}</button><button className={treeViewMode === "stable" ? "active" : ""} type="button" onClick={() => setTreeViewMode("stable")}>{treeRun.origin === "uploaded" ? "Uploaded · floor-collapsed" : "Rooted · floor-collapsed"}</button><button className={treeViewMode === "raw" ? "active" : ""} type="button" onClick={() => setTreeViewMode("raw")}>{treeRun.origin === "uploaded" ? "Uploaded · original order" : "Raw FastTree"}</button></div><span>{treeRun.collapsedEdges.toLocaleString()} numerical-floor internal edges can optionally be displayed as polytomies; the complete resolved tree is the default.</span></div>
            <LineageTreeViewer
              newick={treeViewMode === "stable" ? treeRun.stableNewick : treeViewMode === "rooted" ? treeRun.rootedNewick : treeRun.newick}
              alignmentFasta={treeRun.alignmentFasta}
              rows={workbenchLineageRows}
              multiplicityByOrdinal={lineageMultiplicity}
              sampleColors={sampleColors}
              lineageByOrdinal={originalLineageByOrdinal}
              variant={treeViewMode}
              collapsedEdges={treeRun.collapsedEdges}
              collapseThreshold={treeRun.collapseThreshold}
              mode={alignmentMode}
              onModeChange={setAlignmentMode}
              frameOffset={alignmentFrameOffset}
              isTcr={isTcr}
              name={`${baseName(inputName)}.lineage-${selectedLineage.id}.${treeViewMode === "stable" ? "rooted-floor-collapsed" : treeViewMode === "rooted" ? "rooted-resolved" : "raw-fasttree"}-tree`}
            />
          </>}
        </div>
        <PhyloUcaPanel
          alignment={alignment}
          templateOrdinal={lineageGermlineMethod === "closest" ? lineageGermline?.selectedOrdinal : undefined}
          lineageRows={workbenchLineageRows}
          lineageIds={selectedLineageIds}
          lineageLabel={selectedLineageIds.length > 1 ? `lineages-${selectedLineageIds.join("-")}` : `lineage-${selectedLineage.id}`}
          locus={selectedLineage.locus || workbenchLineageRows[0]?.values.locus || ""}
          references={references}
          inputName={inputName}
          frameOffset={alignmentFrameOffset}
          isTcr={isTcr}
          sampleColors={sampleColors}
          multiplicityByOrdinal={lineageMultiplicity}
          lineageByOrdinal={originalLineageByOrdinal}
          observedTreeNewick={treeRun?.origin === "uploaded" && treeRun.observedOnly ? treeRun.newick : undefined}
          observedTreeSource={treeRun?.origin === "uploaded" && treeRun.observedOnly ? treeRun.source : undefined}
          initialState={phyloUcaState}
          onStateChange={setPhyloUcaState}
        />
      </>}
    </section>}
      </div>
    </div>
  </section>;
}
