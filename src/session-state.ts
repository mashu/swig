import type { UnifiedDashboard, UnifiedOptions } from "./unified-germline";
import type { CompiledReferences, ScopeKey } from "./reference-pack";
import type { RepertoireSelectionOptions } from "./repertoire-selection";
import type { MissingAlleleDashboard, MissingAlleleOptions } from "./germline-evidence";
import type { PersonalizedGermlineDashboard, PersonalizedGermlineOptions } from "./personalized-germline";
import type { ShmDashboard, ShmMetricKey } from "./shm-analysis";
import type { CallingProfile } from "./aligner-engine";
import type { DatasetManifestEntry, PipelinePlan, StudyDesign } from "./study-design";
import type { SampleColorMap } from "./sample-colors";
import type { LineageGermlineMethod } from "./lineage-alignment";
import type { AlignmentFrameOffset } from "./lineage-phylogeny";
import type { FastqQualityFilterOptions, FastqQualityFilterStats } from "./sequence-stream";
import type { PhyloUcaSavedState } from "./phylo-uca/types";
import type { SavedAlleleRefinement } from "./allele-refinement/types";

export const SWIG_SESSION_SCHEMA = 1 as const;

export interface LinkedAirrBinding {
  name: string;
  size: number;
  lastModified: number;
  records: number;
  headers: string[];
  fingerprint: string;
}

export interface SessionVector {
  type: "u8" | "u16" | "u32" | "i32" | "f32";
  length: number;
  base64: string;
}

export interface CollapseReplayState {
  mode: "exact" | "fad" | "conservative" | "indel";
  options: Record<string, string | number | boolean>;
  counts?: SessionVector;
  representatives?: SessionVector;
  representativeMask?: SessionVector;
  dashboard?: Record<string, unknown>;
}

export interface ChimeraReplayState {
  options: Record<string, string | number | boolean>;
  msa?: string;
  dashboard?: Record<string, unknown>;
  filterThreshold: number;
  probabilities?: SessionVector;
  dfr?: SessionVector;
  retainedMask?: SessionVector;
}

export interface LineageReplayState {
  options: Record<string, string | number | boolean>;
  assignments?: SessionVector;
  dashboard?: Record<string, unknown>;
}

export interface PostAnalysisSessionSnapshot {
  /** Explicitly bypassed post-analysis cards. Selection is skipped by default for new workspaces. */
  skippedModules?: Array<"alleles" | "dedup" | "chimera" | "selection" | "lineage" | "diagnostics" | "workbench" | "query">;
  workingStages: Array<{ id: "dedup" | "chimera" | "selection"; label: string; input: number; retained: number; discarded: number; detail: string }>;
  activeMask?: SessionVector;
  collapse?: CollapseReplayState;
  chimera?: ChimeraReplayState;
  selection?: { options: RepertoireSelectionOptions; mask?: SessionVector; baseMask?: SessionVector };
  lineage?: LineageReplayState;
  /** Sparse repertoire-level allele posterior and whether its MAP calls were applied downstream. */
  alleleRefinement?: SavedAlleleRefinement;
  /** Original lineage IDs currently opened together in the workbench. */
  selectedLineageIds?: number[];
  lineageGermlineMethod?: LineageGermlineMethod;
  /** AIRR template ordinal per lineage group/productivity view; absent means automatic. */
  lineageTemplateSelections?: Record<string, number>;
  query?: Record<string, unknown>;
  alignment?: { fasta: string; source: string; selectedLineageId?: number; frameOffset?: AlignmentFrameOffset };
  /** Shared nucleotide-column offset used for codon translation in the current lineage MSA. */
  alignmentFrameOffset?: AlignmentFrameOffset;
  /** True when the lineage workbench excludes AIRR rows not explicitly marked productive. */
  alignmentProductiveOnly?: boolean;
  /** Only manual/corrected alignments are retained; generated alignments are reproducible from AIRR rows. */
  editedAlignments?: Array<{ key: string; lineageIds: number[]; productiveOnly?: boolean; fasta: string; source: string; frameOffset?: AlignmentFrameOffset; savedAt: string }>;
  lineageMerges?: Array<{ id: string; label: string; originalLineageIds: number[]; createdAt: string }>;
  tree?: { rawNewick: string; rootedNewick: string; stableNewick: string; source: string; lineageIds?: number[]; run?: Record<string, unknown> };
  /** Full result is retained because it depends on a searched placement posterior, not only the ordinary lineage tree. */
  phyloUca?: PhyloUcaSavedState;
  shm?: { metric: ShmMetricKey; dashboard: ShmDashboard; sampleOrder?: string[] };
  missingAlleles?: { options: MissingAlleleOptions; dashboard: MissingAlleleDashboard; selectedCandidateIds?: string[] };
  unifiedGermline?: {options:UnifiedOptions;dashboard:UnifiedDashboard};
  personalizedGermline?: { options: PersonalizedGermlineOptions; dashboard: PersonalizedGermlineDashboard };
}

export interface SwigSession {
  schema: typeof SWIG_SESSION_SCHEMA;
  application: "Swig";
  applicationVersion: string;
  savedAt: string;
  linkedAirr: LinkedAirrBinding;
  analysis: {
    /** Browser entry route; omitted by sessions created before the focused workflows existed. */
    webMode?: "vdj" | "repertoire" | "lineage";
    inputName: string;
    species: string;
    scope: ScopeKey;
    workers: number;
    callingProfile?: CallingProfile;
    minimumIdentity: number;
    minimumConstantPrefixIdentity?: number | null;
    strand: 0 | 1 | 2;
    subsample?: { enabled: boolean; size: number; seed: number };
    fastqFilter?: FastqQualityFilterOptions;
    fastqFilterStats?: FastqQualityFilterStats;
    references: CompiledReferences;
    /** Exact per-locus/segment FASTA identifiers removed before assignment. */
    referenceExclusions?: Record<string, string[]>;
    doubleD?: Record<string, string | number | boolean>;
    datasets?: DatasetManifestEntry[];
    studyDesign?: StudyDesign;
    pipeline?: PipelinePlan;
    sampleColors?: SampleColorMap;
  };
  doubleD: Array<{ ordinal: number; values: Record<string, string> }>;
  postAnalysis: PostAnalysisSessionSnapshot;
}

function bytesToBase64(bytes: Uint8Array): string {
  let output = "";
  for (let offset = 0; offset < bytes.length; offset += 0x8000) output += String.fromCharCode(...bytes.subarray(offset, offset + 0x8000));
  return btoa(output);
}

function base64ToBytes(value: string): Uint8Array {
  const binary = atob(value); const result = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) result[index] = binary.charCodeAt(index);
  return result;
}

export function packSessionVector(value: Uint8Array | Uint16Array | Uint32Array | Int32Array | Float32Array): SessionVector {
  const type = value instanceof Uint8Array ? "u8" : value instanceof Uint16Array ? "u16" : value instanceof Uint32Array ? "u32" : value instanceof Int32Array ? "i32" : "f32";
  return { type, length: value.length, base64: bytesToBase64(new Uint8Array(value.buffer, value.byteOffset, value.byteLength)) };
}

export function unpackSessionVector(value: SessionVector): Uint8Array | Uint16Array | Uint32Array | Int32Array | Float32Array {
  const bytes = base64ToBytes(value.base64); const copy = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
  const result = value.type === "u8" ? new Uint8Array(copy) : value.type === "u16" ? new Uint16Array(copy) : value.type === "u32" ? new Uint32Array(copy) : value.type === "i32" ? new Int32Array(copy) : new Float32Array(copy);
  if (result.length !== value.length) throw new Error("The saved session contains a malformed analysis vector.");
  return result;
}

export function validateSession(value: unknown): SwigSession {
  if (!value || typeof value !== "object") throw new Error("This file is not a Swig session.");
  const session = value as Partial<SwigSession>;
  if (session.schema !== SWIG_SESSION_SCHEMA || session.application !== "Swig" || !session.linkedAirr || !session.analysis || !session.postAnalysis) {
    throw new Error("This session uses an unsupported or incomplete Swig schema.");
  }
  return session as SwigSession;
}

function sessionAbortError(): DOMException {
  return new DOMException("Session saving was cancelled.", "AbortError");
}

export async function encodeSession(session: SwigSession, signal?: AbortSignal): Promise<Blob> {
  if (signal?.aborted) throw sessionAbortError();
  const json = new TextEncoder().encode(JSON.stringify(session));
  if (signal?.aborted) throw sessionAbortError();
  if ("CompressionStream" in globalThis) {
    const reader = new Blob([json]).stream().pipeThrough(new CompressionStream("gzip")).getReader();
    const chunks: BlobPart[] = [];
    const abort = () => { void reader.cancel(sessionAbortError()); };
    signal?.addEventListener("abort", abort, { once: true });
    try {
      while (true) {
        if (signal?.aborted) throw sessionAbortError();
        const { done, value } = await reader.read();
        if (done) break;
        chunks.push(value.buffer.slice(value.byteOffset, value.byteOffset + value.byteLength) as ArrayBuffer);
      }
    } finally {
      signal?.removeEventListener("abort", abort);
      reader.releaseLock();
    }
    if (signal?.aborted) throw sessionAbortError();
    return new Blob(chunks, { type: "application/gzip" });
  }
  if (signal?.aborted) throw sessionAbortError();
  return new Blob([json], { type: "application/json" });
}

export async function decodeSession(file: Blob): Promise<SwigSession> {
  const bytes = new Uint8Array(await file.arrayBuffer());
  let text: string;
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
    if (!("DecompressionStream" in globalThis)) throw new Error("This browser cannot read gzip-compressed Swig sessions.");
    text = await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text();
  } else text = new TextDecoder().decode(bytes);
  try { return validateSession(JSON.parse(text)); } catch (error) {
    if (error instanceof SyntaxError) throw new Error("The selected file is not valid JSON or gzip-compressed JSON.");
    throw error;
  }
}

export function linkedAirrMatches(session: SwigSession, binding: LinkedAirrBinding): string[] {
  const expected=session.linkedAirr;const errors:string[]=[];
  if(binding.records!==expected.records)errors.push(`record count ${binding.records.toLocaleString()} does not match ${expected.records.toLocaleString()}`);
  if(binding.fingerprint!==expected.fingerprint)errors.push("content fingerprint differs");
  if(binding.headers.join("\t")!==expected.headers.join("\t"))errors.push("AIRR columns differ");
  return errors;
}

export function sessionBaseName(name:string):string{return (name.replace(/\.swig-session(?:\.json)?(?:\.gz)?$/i,"").replace(/(\.airr)?\.(tsv|csv|txt)(\.gz)?$/i,"")||"swig")+".swig-session.json.gz";}
