import type { LineageSummary } from "./post-analysis-core.ts";
import type { CompiledReferences, ScopeKey } from "./reference-pack.ts";
import type { ShmMetricKey } from "./shm-analysis.ts";
import type { DatasetManifestEntry, DatasetScope } from "./study-design.ts";
import type { CallingProfile } from "./aligner-engine.ts";

export const SWIG_LINEAGE_STUDY_SCHEMA = 1 as const;

export interface LineageStudyRange {
  lineageId: number;
  /** UTF-8 byte offsets in the linked, uncompressed AIRR TSV. */
  start: number;
  end: number;
  records: number;
  /** SHA-256 over this exact row slice, allowing lazy integrity checks. */
  sha256: string;
}

export interface LineageStudyShmSummary {
  lineageId: number;
  mean: number;
  p95: number;
  /** Exact lowest-V-SHM member retained during the original streaming scan. */
  ordinal?: number;
  cdr3Nt?: string;
  cdr3Aa?: string;
}

export interface LineageStudyManifest {
  schema: typeof SWIG_LINEAGE_STUDY_SCHEMA;
  application: "Swig lineage study";
  applicationVersion: string;
  createdAt: string;
  linkedAirr: {
    name: string;
    size: number;
    records: number;
    headers: string[];
    /** SHA-256 over the exact linked AIRR bytes. */
    sha256: string;
  };
  analysis: {
    inputName: string;
    species: string;
    scope: ScopeKey;
    references: CompiledReferences;
    datasets: DatasetManifestEntry[];
    callingProfile: CallingProfile;
    minimumIdentity: number;
    strand: 0|1|2;
    lineage: {
      scope: DatasetScope;
      identity: number;
      resolution: "gene"|"allele";
      ambiguity: "overlap"|"top"|"strict";
      productiveOnly: boolean;
      maxCandidateComparisons: number;
    };
  };
  summaries: LineageSummary[];
  ranges: LineageStudyRange[];
  shm?: { metric: ShmMetricKey; summaries: LineageStudyShmSummary[] };
}

export function validateLineageStudy(value: unknown): LineageStudyManifest {
  if(!value||typeof value!=="object")throw new Error("This file is not a Swig lineage-study manifest.");
  const manifest=value as Partial<LineageStudyManifest>;
  if(manifest.schema!==SWIG_LINEAGE_STUDY_SCHEMA||manifest.application!=="Swig lineage study"||!manifest.linkedAirr||!manifest.analysis||!Array.isArray(manifest.summaries)||!Array.isArray(manifest.ranges))throw new Error("This lineage-study file uses an unsupported or incomplete schema.");
  const ranges=[...manifest.ranges].sort((left,right)=>left.start-right.start);
  let previous=0;const ids=new Set<number>();
  for(const range of ranges){
    if(!Number.isInteger(range.lineageId)||range.lineageId<=0||ids.has(range.lineageId)||!Number.isInteger(range.start)||!Number.isInteger(range.end)||range.start<previous||range.end<range.start||range.end>manifest.linkedAirr.size||!/^[0-9a-f]{64}$/i.test(range.sha256))throw new Error("The lineage-study manifest contains an invalid AIRR byte range.");
    ids.add(range.lineageId);
    previous=range.end;
  }
  return manifest as LineageStudyManifest;
}

export async function decodeLineageStudy(file: Blob): Promise<LineageStudyManifest> {
  const bytes=new Uint8Array(await file.arrayBuffer());
  let text:string;
  if(bytes[0]===0x1f&&bytes[1]===0x8b){
    if(!("DecompressionStream" in globalThis))throw new Error("This browser cannot open compressed lineage-study manifests.");
    text=await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text();
  }else text=new TextDecoder().decode(bytes);
  return validateLineageStudy(JSON.parse(text));
}

export async function sha256Hex(blob: Blob): Promise<string> {
  const digest=await crypto.subtle.digest("SHA-256",await blob.arrayBuffer());
  return [...new Uint8Array(digest)].map((value)=>value.toString(16).padStart(2,"0")).join("");
}

export function lineageStudyRange(manifest: LineageStudyManifest, lineageId: number): LineageStudyRange {
  const range=manifest.ranges.find((candidate)=>candidate.lineageId===lineageId);
  if(!range)throw new Error(`Lineage ${lineageId} has no linked AIRR byte range.`);
  return range;
}

/** Read only the selected lineage. The linked AIRR file is deliberately never indexed in full. */
export async function readLineageAirrSlice(file: File, manifest: LineageStudyManifest, lineageId: number): Promise<{header:string;body:string;records:number}> {
  if(file.size!==manifest.linkedAirr.size)throw new Error(`The selected AIRR file has ${file.size.toLocaleString()} bytes; the manifest expects ${manifest.linkedAirr.size.toLocaleString()}.`);
  const range=lineageStudyRange(manifest,lineageId);
  const source=file.slice(range.start,range.end);
  const [body,digest]=await Promise.all([source.text(),sha256Hex(source)]);
  if(digest.toLowerCase()!==range.sha256.toLowerCase())throw new Error(`Lineage ${lineageId} does not match the indexed AIRR content. Select the exact linked AIRR file emitted with this study.`);
  const records=body.split(/\r?\n/).filter(Boolean).length;
  if(records!==range.records)throw new Error(`Lineage ${lineageId} contains ${records.toLocaleString()} rows; the manifest expects ${range.records.toLocaleString()}.`);
  return {header:manifest.linkedAirr.headers.join("\t"),body,records};
}
