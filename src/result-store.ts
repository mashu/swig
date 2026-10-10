import { tableHeader, tableRow, type TableExportFormat } from "./export-formats.ts";
import type { DatasetManifestEntry } from "./study-design.ts";

export interface AirrIndexRecord {
  ordinal: number;
  chunk: number;
  line: number;
  sequenceId: string;
  datasetId: string;
  sampleId: string;
  subjectId: string;
  cohort: string;
  timepoint: string;
  compartment: string;
  locus: string;
  vCall: string;
  dCall: string;
  d2Call: string;
  jCall: string;
  cCall: string;
  isotype: string;
  productive: string;
  cdr3: string;
  cdr3Aa: string;
  junctionAa: string;
  vIdentity: number | null;
  dIdentity: number | null;
  jIdentity: number | null;
  cIdentity: number | null;
  cdr3AaLength: number | null;
  vjInFrame: string;
  stopCodon: string;
  completeVdj: string;
  revComp: string;
}

export interface ResultFilters {
  sequenceId: string;
  cdr3: string;
  datasetId: string;
  sampleId: string;
  subjectId: string;
  cohort: string;
  timepoint: string;
  compartment: string;
  locus: string;
  productive: string;
  vCall: string;
  vCallIncludeAmbiguous: boolean;
  dCall: string;
  dCallIncludeAmbiguous: boolean;
  jCall: string;
  jCallIncludeAmbiguous: boolean;
  cCall: string;
  cCallIncludeAmbiguous: boolean;
  isotype: string;
  minVIdentity: number;
  minDIdentity: number;
  minJIdentity: number;
  minCIdentity: number;
  minCdr3AaLength: number;
  maxCdr3AaLength: number;
  vjInFrame: string;
  stopCodon: string;
  completeVdj: string;
  revComp: string;
  hasD: boolean;
  hasDoubleD: boolean;
  hasCdr3: boolean;
}

export interface FacetValue {
  value: string;
  count: number;
}

export interface ResultFacets {
  datasets: FacetValue[];
  samples: FacetValue[];
  subjects: FacetValue[];
  cohorts: FacetValue[];
  timepoints: FacetValue[];
  compartments: FacetValue[];
  loci: FacetValue[];
  productive: FacetValue[];
  vCalls: FacetValue[];
  dCalls: FacetValue[];
  jCalls: FacetValue[];
  cCalls: FacetValue[];
  isotypes: FacetValue[];
}

export interface RepertoireOptions {
  locus?: string;
  productiveOnly?: boolean;
  ambiguity?: "top" | "fractional";
}

export interface RepertoirePair {
  v: string;
  j: string;
  count: number;
}

export interface RepertoireSnapshot {
  records: number;
  vCalls: FacetValue[];
  dCalls: FacetValue[];
  jCalls: FacetValue[];
  cCalls: FacetValue[];
  isotypes: FacetValue[];
  cdr3Lengths: FacetValue[];
  vIdentityBins: FacetValue[];
  vjPairs: RepertoirePair[];
}

export interface ResultPage {
  rows: AirrIndexRecord[];
  hasMore: boolean;
  totalMatches: number | null;
  scanned: number;
}

type AirrRow = Record<string, string>;

export interface AirrScanRow {
  ordinal: number;
  values: AirrRow;
}

export interface AirrDetailRow {
  record: AirrIndexRecord;
  values: AirrRow;
}

interface PackedIndexRecord {
  o: number;
  c: number;
  n: number;
  i: string;
  e: string;
  m: string;
  b: string;
  g: string;
  t: string;
  w: string;
  l: string;
  v: string;
  va: string[];
  d: string;
  da: string[];
  h: string;
  j: string;
  ja: string[];
  k: string;
  ka: string[];
  y: string;
  p: string;
  r: string;
  a: string;
  u: string;
  vi: number | null;
  di: number | null;
  ji: number | null;
  ci: number | null;
  z: number | null;
  f: string;
  s: string;
  q: string;
  x: string;
}

interface ChunkRecord {
  index: number;
  storage: "indexed" | "external";
  data?: Blob;
  compressed?: boolean;
  start?: number;
  length?: number;
}

interface ManifestRecord {
  key: "manifest";
  headerLine: string;
  chunks: number;
  records: number;
  outputBytes: number;
}

export interface DoubleDEvidenceRecord {
  ordinal: number;
  values: Record<string, string>;
}

export interface DoubleDBatch {
  header: string;
  body: string | Uint8Array;
}

export interface AirrOutputWritable {
  write: (data: string | Blob | Uint8Array) => Promise<void>;
  close: () => Promise<void>;
  abort?: () => Promise<void>;
}

export interface AirrOutputHandle {
  getFile: () => Promise<File>;
}

export interface DirectAirrOutput {
  handle: AirrOutputHandle;
  writable: AirrOutputWritable;
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();

function updateHash(states: Uint32Array, value: string) {
  const bytes = encoder.encode(value);
  for (const byte of bytes) {
    states[0] = Math.imul(states[0] ^ byte, 0x01000193) >>> 0;
    states[1] = Math.imul((states[1] + byte + 0x9e3779b9) ^ (states[1] >>> 13), 0x85ebca6b) >>> 0;
    states[2] = Math.imul(states[2] ^ (byte + (states[2] << 6) + (states[2] >>> 2)), 0xc2b2ae35) >>> 0;
    states[3] = Math.imul((states[3] ^ byte) + 0x27d4eb2d, 0x165667b1) >>> 0;
  }
}

function cooperativeYield(): Promise<void> {
  const scheduler = (globalThis as typeof globalThis & { scheduler?: { yield?: () => Promise<void> } }).scheduler;
  if (scheduler?.yield) return scheduler.yield();
  return new Promise((resolve) => {
    const channel = new MessageChannel();
    channel.port1.onmessage = () => {
      channel.port1.close();
      channel.port2.close();
      resolve();
    };
    channel.port2.postMessage(null);
  });
}
const MAX_FALLBACK_BLOB_BYTES = 256 * 1024 * 1024;

function requestResult<T>(request: IDBRequest<T>): Promise<T> {
  return new Promise((resolve, reject) => {
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error ?? new Error("The local result index failed."));
  });
}

function transactionDone(transaction: IDBTransaction): Promise<void> {
  return new Promise((resolve, reject) => {
    transaction.oncomplete = () => resolve();
    transaction.onerror = () => reject(transaction.error ?? new Error("The local result index failed."));
    transaction.onabort = () => reject(transaction.error ?? new Error("The local result index was interrupted."));
  });
}

function bump(map: Map<string, number>, value: string) {
  if (!value) return;
  map.set(value, (map.get(value) ?? 0) + 1);
}

function bumpBy(map: Map<string, number>, value: string, count: number) {
  if (!value || !count) return;
  map.set(value, (map.get(value) ?? 0) + count);
}

function calls(value: string): string[] {
  return value.split(",").map((call) => call.trim()).filter(Boolean);
}

function callTokens(value: string): string[] {
  const tokens = new Set<string>();
  for (const call of calls(value)) {
    const normalized = call.toUpperCase();
    tokens.add(normalized);
    tokens.add(normalized.replace(/\*.*$/, ""));
  }
  return [...tokens].filter(Boolean);
}

function matchesCall(value: string, query: string, includeAmbiguous: boolean): boolean {
  const assignments = calls(value);
  if (!query.trim()) return true;
  const requestedAssignments=calls(query).map((call)=>call.toUpperCase()).sort();
  if(requestedAssignments.length>1){
    const observed=assignments.map((call)=>call.toUpperCase()).sort();
    return observed.length===requestedAssignments.length&&observed.every((call,index)=>call===requestedAssignments[index]);
  }
  if (!includeAmbiguous && assignments.length > 1) return false;
  const normalized = query.trim().toUpperCase();
  return assignments.some((assignment) => {
    const allele = assignment.toUpperCase();
    return allele === normalized || allele.replace(/\*.*$/, "") === normalized;
  });
}

function topCall(value: string): string {
  return calls(value)[0] ?? "";
}

function addCallCounts(top: Map<string, number>, fractional: Map<string, number>, value: string) {
  const values = calls(value);
  if (!values.length) return;
  bump(top, values[0]);
  const weight = 1 / values.length;
  values.forEach((call) => bumpBy(fractional, call, weight));
}

function isotypeLabel(call: string): string {
  const gene = call.toUpperCase().replace(/\*.*$/, "");
  const heavy = gene.match(/IGH([MDE])(?:\d+)?$/) ?? gene.match(/IGH(G|A)(\d+)?$/);
  if (heavy) return `Ig${heavy[1]}${heavy[2] ?? ""}`;
  if (gene.startsWith("IGKC")) return "Igκ constant";
  if (gene.startsWith("IGLC")) return "Igλ constant";
  if (gene.startsWith("TR") && gene.includes("C")) return gene;
  return gene;
}

export function inferIsotype(
  call: string,
  alignedSequence: string,
  identity: number | null,
  support: number | null = null,
): string {
  const alignedBases = alignedSequence.replace(/[-.\s]/g, "").length;
  const longLocalEvidence = alignedBases >= 30;
  const shortAnchoredEvidence = alignedBases >= 12 && support !== null && support <= 1e-4;
  if (!call || (identity ?? 0) < 0.90 || (!longLocalEvidence && !shortAnchoredEvidence)) return "";
  const labels = [...new Set(calls(call).map(isotypeLabel).filter(Boolean))];
  if (labels.length <= 1) return labels[0] ?? "";
  if (labels.every((label) => /^IgG\d*$/.test(label))) return "IgG";
  if (labels.every((label) => /^IgA\d*$/.test(label))) return "IgA";
  return labels.join("/");
}

interface RepertoireBucket {
  records: number;
  vTop: Map<string, number>;
  vFractional: Map<string, number>;
  dTop: Map<string, number>;
  dFractional: Map<string, number>;
  jTop: Map<string, number>;
  jFractional: Map<string, number>;
  cTop: Map<string, number>;
  cFractional: Map<string, number>;
  isotypes: Map<string, number>;
  cdr3Lengths: Map<string, number>;
  vIdentityBins: Map<string, number>;
  vjPairsTop: Map<string, number>;
  vjPairsFractional: Map<string, number>;
}

function makeRepertoireBucket(): RepertoireBucket {
  return {
    records: 0,
    vTop: new Map(), vFractional: new Map(),
    dTop: new Map(), dFractional: new Map(),
    jTop: new Map(), jFractional: new Map(),
    cTop: new Map(), cFractional: new Map(),
    isotypes: new Map(), cdr3Lengths: new Map(), vIdentityBins: new Map(),
    vjPairsTop: new Map(), vjPairsFractional: new Map(),
  };
}

function facet(map: Map<string, number>): FacetValue[] {
  return [...map.entries()]
    .map(([value, count]) => ({ value, count }))
    .sort((a, b) => b.count - a.count || a.value.localeCompare(b.value));
}

function numeric(value: string): number | null {
  const parsed = Number(value);
  return value !== "" && Number.isFinite(parsed) ? parsed : null;
}

function normalizedStudyMetadata(datasets: readonly DatasetManifestEntry[]): Map<string, DatasetManifestEntry> {
  const byDataset = new Map<string, DatasetManifestEntry>();
  for (const dataset of datasets) {
    const datasetId = dataset.datasetId.trim();
    if (!datasetId) throw new Error("Every dataset must retain a non-empty dataset ID.");
    if (byDataset.has(datasetId)) throw new Error(`Dataset ID ${datasetId} occurs more than once.`);
    byDataset.set(datasetId, {
      ...dataset,
      datasetId,
      sampleId: dataset.sampleId.trim(),
      subjectId: dataset.subjectId.trim(),
      cohort: dataset.cohort.trim(),
      timepoint: dataset.timepoint.trim(),
      compartment: (dataset.compartment ?? "").trim(),
    });
  }
  return byDataset;
}

function packRecord(record: AirrIndexRecord): PackedIndexRecord {
  return {
    o: record.ordinal, c: record.chunk, n: record.line, i: record.sequenceId,
    e: record.datasetId, m: record.sampleId, b: record.subjectId, g: record.cohort, t: record.timepoint, w: record.compartment,
    l: record.locus, v: record.vCall, va: callTokens(record.vCall), d: record.dCall, da: callTokens(record.dCall), h: record.d2Call, j: record.jCall, ja: callTokens(record.jCall),
    k: record.cCall, ka: callTokens(record.cCall), y: record.isotype,
    p: record.productive, r: record.cdr3, a: record.cdr3Aa, u: record.junctionAa,
    vi: record.vIdentity, di: record.dIdentity, ji: record.jIdentity, ci: record.cIdentity,
    z: record.cdr3AaLength, f: record.vjInFrame, s: record.stopCodon,
    q: record.completeVdj, x: record.revComp,
  };
}

function unpackRecord(record: PackedIndexRecord): AirrIndexRecord {
  return {
    ordinal: record.o, chunk: record.c, line: record.n, sequenceId: record.i,
    datasetId: record.e ?? "", sampleId: record.m ?? "", subjectId: record.b ?? "", cohort: record.g ?? "", timepoint: record.t ?? "", compartment: record.w ?? "",
    locus: record.l, vCall: record.v, dCall: record.d, d2Call: record.h ?? "", jCall: record.j,
    cCall: record.k, isotype: record.y,
    productive: record.p, cdr3: record.r, cdr3Aa: record.a, junctionAa: record.u,
    vIdentity: record.vi, dIdentity: record.di, jIdentity: record.ji, cIdentity: record.ci,
    cdr3AaLength: record.z, vjInFrame: record.f, stopCodon: record.s,
    completeVdj: record.q, revComp: record.x,
  };
}

export const EMPTY_FILTERS: ResultFilters = {
  sequenceId: "",
  cdr3: "",
  datasetId: "",
  sampleId: "",
  subjectId: "",
  cohort: "",
  timepoint: "",
  compartment: "",
  locus: "",
  productive: "",
  vCall: "",
  vCallIncludeAmbiguous: false,
  dCall: "",
  dCallIncludeAmbiguous: false,
  jCall: "",
  jCallIncludeAmbiguous: false,
  cCall: "",
  cCallIncludeAmbiguous: false,
  isotype: "",
  minVIdentity: 0,
  minDIdentity: 0,
  minJIdentity: 0,
  minCIdentity: 0,
  minCdr3AaLength: 0,
  maxCdr3AaLength: 0,
  vjInFrame: "",
  stopCodon: "",
  completeVdj: "",
  revComp: "",
  hasD: false,
  hasDoubleD: false,
  hasCdr3: false,
};

export class AirrResultStore {
  readonly databaseName = `swig-results-${Date.now()}-${crypto.randomUUID()}`;
  private readonly database: Promise<IDBDatabase>;
  private readonly facetMaps = {
    datasets: new Map<string, number>(),
    samples: new Map<string, number>(),
    subjects: new Map<string, number>(),
    cohorts: new Map<string, number>(),
    timepoints: new Map<string, number>(),
    compartments: new Map<string, number>(),
    loci: new Map<string, number>(),
    productive: new Map<string, number>(),
    vCalls: new Map<string, number>(),
    dCalls: new Map<string, number>(),
    jCalls: new Map<string, number>(),
    cCalls: new Map<string, number>(),
    isotypes: new Map<string, number>(),
  };
  private readonly repertoireBuckets = new Map<string, RepertoireBucket>();
  private headerLine = "";
  private headers: string[] = [];
  private nextChunk = 0;
  private nextOrdinal = 0;
  private assigned = 0;
  private productive = 0;
  private withCdr3 = 0;
  private readonly directOutput?: DirectAirrOutput;
  private outputByteCount = 0;
  private doubleDHeaderLine = "";
  private doubleDHeaders: string[] = [];
  private doubleDRecordCount = 0;
  private readonly fingerprintState = new Uint32Array([0x811c9dc5, 0x9e3779b9, 0x85ebca6b, 0xc2b2ae35]);
  private finalized = false;
  private studyMetadataOverrides = new Map<string, DatasetManifestEntry>();

  constructor(directOutput?: DirectAirrOutput) {
    this.directOutput = directOutput;
    this.database = new Promise((resolve, reject) => {
      const request = indexedDB.open(this.databaseName, 1);
      request.onupgradeneeded = () => {
        const database = request.result;
        database.createObjectStore("chunks", { keyPath: "index" });
        database.createObjectStore("meta", { keyPath: "key" });
        database.createObjectStore("doubleD", { keyPath: "ordinal" });
        const records = database.createObjectStore("records", { keyPath: "o" });
        records.createIndex("datasetId", "e");
        records.createIndex("sampleId", "m");
        records.createIndex("subjectId", "b");
        records.createIndex("cohort", "g");
        records.createIndex("timepoint", "t");
        records.createIndex("compartment", "w");
        records.createIndex("locus", "l");
        records.createIndex("productive", "p");
        records.createIndex("vCall", "v");
        records.createIndex("vCallToken", "va", { multiEntry: true });
        records.createIndex("dCall", "d");
        records.createIndex("dCallToken", "da", { multiEntry: true });
        records.createIndex("d2Call", "h");
        records.createIndex("jCall", "j");
        records.createIndex("jCallToken", "ja", { multiEntry: true });
        records.createIndex("cCall", "k");
        records.createIndex("cCallToken", "ka", { multiEntry: true });
        records.createIndex("isotype", "y");
      };
      request.onsuccess = () => resolve(request.result);
      request.onerror = () => reject(request.error ?? new Error("Could not create the local result index."));
    });
  }

  get count(): number {
    return this.nextOrdinal;
  }

  get airrHeaders(): string[] {
    return this.headers;
  }

  get outputBytes(): number {
    return this.outputByteCount;
  }

  get streamedDirectly(): boolean {
    return Boolean(this.directOutput);
  }

  get doubleDCount(): number {
    return this.doubleDRecordCount;
  }

  get hasStudyMetadataOverrides(): boolean {
    return this.studyMetadataOverrides.size > 0;
  }

  /**
   * Install saved dataset metadata before the first AIRR batch is indexed.
   * Import can then materialize the final searchable values in its existing
   * write, avoiding a second read-modify-write pass over every record.
   */
  configureStudyMetadataForImport(datasets: readonly DatasetManifestEntry[]): void {
    if (this.nextOrdinal || this.nextChunk || this.finalized) {
      throw new Error("Saved study metadata must be configured before AIRR records are imported.");
    }
    this.studyMetadataOverrides = normalizedStudyMetadata(datasets);
  }

  get fingerprint(): string {
    return [...this.fingerprintState].map((value) => value.toString(16).padStart(8, "0")).join("");
  }

  get summary() {
    return { assigned: this.assigned, productive: this.productive, withCdr3: this.withCdr3 };
  }

  facets(): ResultFacets {
    return {
      datasets: facet(this.facetMaps.datasets),
      samples: facet(this.facetMaps.samples),
      subjects: facet(this.facetMaps.subjects),
      cohorts: facet(this.facetMaps.cohorts),
      timepoints: facet(this.facetMaps.timepoints),
      compartments: facet(this.facetMaps.compartments),
      loci: facet(this.facetMaps.loci),
      productive: facet(this.facetMaps.productive),
      vCalls: facet(this.facetMaps.vCalls),
      dCalls: facet(this.facetMaps.dCalls),
      jCalls: facet(this.facetMaps.jCalls),
      cCalls: facet(this.facetMaps.cCalls),
      isotypes: facet(this.facetMaps.isotypes),
    };
  }

  private overlayStudyMetadata(values: AirrRow, datasetId = values.swig_dataset_id): AirrRow {
    if (!this.studyMetadataOverrides.size) return values;
    const dataset = this.studyMetadataOverrides.get(datasetId);
    if (!dataset) return values;
    return {
      ...values,
      sample_id: dataset.sampleId,
      subject_id: dataset.subjectId,
      swig_cohort: dataset.cohort,
      swig_timepoint: dataset.timepoint,
      swig_compartment: dataset.compartment ?? "",
    };
  }

  /**
   * Re-index editable study metadata without touching sequence assignments.
   * Chunk payloads stay immutable; AIRR scans/details/exports receive a small
   * dataset-keyed overlay so a million-row table is not duplicated in storage.
   */
  async updateStudyMetadata(
    datasets: readonly DatasetManifestEntry[],
    onProgress?: (processed: number, total: number) => void,
    signal?: AbortSignal,
  ): Promise<ResultFacets> {
    if (signal?.aborted) throw new DOMException("Study-metadata re-indexing was cancelled.", "AbortError");
    const byDataset = normalizedStudyMetadata(datasets);
    const nextSamples = new Map<string, number>();
    const nextSubjects = new Map<string, number>();
    const nextCohorts = new Map<string, number>();
    const nextTimepoints = new Map<string, number>();
    const nextCompartments = new Map<string, number>();
    const database = await this.database;
    const transaction = database.transaction("records", "readwrite");
    const abort = () => transaction.abort();
    signal?.addEventListener("abort", abort, { once: true });
    const records = transaction.objectStore("records");
    let processed = 0;
    try {
      const scan = new Promise<void>((resolve, reject) => {
      const request = records.openCursor();
      request.onerror = () => reject(request.error ?? new Error("Could not update study metadata in the local index."));
      request.onsuccess = () => {
        if (signal?.aborted) {
          transaction.abort();
          return;
        }
        const cursor = request.result;
        if (!cursor) {
          resolve();
          return;
        }
        const current = unpackRecord(cursor.value as PackedIndexRecord);
        const dataset = byDataset.get(current.datasetId);
        const updated = dataset ? {
          ...current,
          sampleId: dataset.sampleId,
          subjectId: dataset.subjectId,
          cohort: dataset.cohort,
          timepoint: dataset.timepoint,
          compartment: dataset.compartment ?? "",
        } : current;
        cursor.update(packRecord(updated));
        bump(nextSamples, updated.sampleId);
        bump(nextSubjects, updated.subjectId);
        bump(nextCohorts, updated.cohort);
        bump(nextTimepoints, updated.timepoint);
        bump(nextCompartments, updated.compartment);
        processed += 1;
        if (processed % 500 === 0) onProgress?.(processed, this.count);
        if (signal?.aborted) return;
        cursor.continue();
      };
      });
      // Wait for both cursor exhaustion and the transaction commit. On an
      // explicit abort some IndexedDB implementations do not dispatch a
      // cursor error, so transactionDone must be observed concurrently rather
      // than only after the cursor promise resolves.
      await Promise.all([scan, transactionDone(transaction)]);
    } catch (error) {
      if (signal?.aborted) throw new DOMException("Study-metadata re-indexing was cancelled.", "AbortError");
      throw error;
    } finally {
      signal?.removeEventListener("abort", abort);
    }
    this.studyMetadataOverrides = byDataset;
    this.facetMaps.samples.clear();
    this.facetMaps.subjects.clear();
    this.facetMaps.cohorts.clear();
    this.facetMaps.timepoints.clear();
    this.facetMaps.compartments.clear();
    nextSamples.forEach((count, value) => this.facetMaps.samples.set(value, count));
    nextSubjects.forEach((count, value) => this.facetMaps.subjects.set(value, count));
    nextCohorts.forEach((count, value) => this.facetMaps.cohorts.set(value, count));
    nextTimepoints.forEach((count, value) => this.facetMaps.timepoints.set(value, count));
    nextCompartments.forEach((count, value) => this.facetMaps.compartments.set(value, count));
    onProgress?.(processed, this.count);
    return this.facets();
  }

  repertoire(options: RepertoireOptions = {}): RepertoireSnapshot {
    const key = `${options.productiveOnly ? "productive" : "all"}|${options.locus || "*"}`;
    const bucket = this.repertoireBuckets.get(key) ?? makeRepertoireBucket();
    const fractional = options.ambiguity === "fractional";
    return {
      records: bucket.records,
      vCalls: facet(fractional ? bucket.vFractional : bucket.vTop),
      dCalls: facet(fractional ? bucket.dFractional : bucket.dTop),
      jCalls: facet(fractional ? bucket.jFractional : bucket.jTop),
      cCalls: facet(fractional ? bucket.cFractional : bucket.cTop),
      isotypes: facet(bucket.isotypes),
      cdr3Lengths: facet(bucket.cdr3Lengths).sort((a, b) => Number(a.value) - Number(b.value)),
      vIdentityBins: facet(bucket.vIdentityBins).sort((a, b) => Number(a.value) - Number(b.value)),
      vjPairs: [...(fractional ? bucket.vjPairsFractional : bucket.vjPairsTop).entries()]
        .map(([value, count]) => {
          const [v, j] = value.split("\u0000");
          return { v, j, count };
        })
        .sort((a, b) => b.count - a.count || a.v.localeCompare(b.v) || a.j.localeCompare(b.j)),
    };
  }

  private accumulateRepertoire(record: AirrIndexRecord) {
    const keys = [`all|*`];
    if (record.locus) keys.push(`all|${record.locus}`);
    if (record.productive === "T") {
      keys.push("productive|*");
      if (record.locus) keys.push(`productive|${record.locus}`);
    }
    for (const key of keys) {
      let bucket = this.repertoireBuckets.get(key);
      if (!bucket) {
        bucket = makeRepertoireBucket();
        this.repertoireBuckets.set(key, bucket);
      }
      bucket.records += 1;
      addCallCounts(bucket.vTop, bucket.vFractional, record.vCall);
      addCallCounts(bucket.dTop, bucket.dFractional, record.dCall);
      addCallCounts(bucket.jTop, bucket.jFractional, record.jCall);
      addCallCounts(bucket.cTop, bucket.cFractional, record.cCall);
      bump(bucket.isotypes, record.isotype);
      if (record.cdr3AaLength !== null) bump(bucket.cdr3Lengths, String(record.cdr3AaLength));
      if (record.vIdentity !== null) bump(bucket.vIdentityBins, String(Math.round(record.vIdentity * 100)));
      const vValues = calls(record.vCall);
      const jValues = calls(record.jCall);
      if (vValues.length && jValues.length) {
        bump(bucket.vjPairsTop, `${vValues[0]}\u0000${jValues[0]}`);
        const weight = 1 / (vValues.length * jValues.length);
        for (const v of vValues) for (const j of jValues) {
          bumpBy(bucket.vjPairsFractional, `${v}\u0000${j}`, weight);
        }
      }
    }
  }

  async appendBatch(headerLine: string, body: string | Uint8Array, doubleD?: DoubleDBatch): Promise<void> {
    if (this.finalized) throw new Error("Cannot append to a finalized AIRR result store.");
    let normalizedHeader = headerLine.replace(/\r$/, "");
    let bodyText = typeof body === "string" ? body : decoder.decode(body);
    const incomingHeaders = normalizedHeader.split("\t");
    if (!incomingHeaders.includes("isotype")) {
      const incomingPositions = Object.fromEntries(incomingHeaders.map((name, index) => [name, index]));
      bodyText = bodyText.split("\n").map((line) => {
        const clean = line.replace(/\r$/, "");
        if (!clean) return "";
        const values = clean.split("\t");
        const at = (name: string) => values[incomingPositions[name]] ?? "";
        return `${clean}\t${inferIsotype(at("c_call"), at("c_sequence_alignment"), numeric(at("c_identity")), numeric(at("c_support")))}`;
      }).join("\n");
      normalizedHeader += "\tisotype";
    }
    if (!this.headerLine) {
      this.headerLine = normalizedHeader;
      this.headers = this.headerLine.split("\t");
      updateHash(this.fingerprintState, `${this.headerLine}\n`);
      const headerBytes = encoder.encode(`${this.headerLine}\n`);
      this.outputByteCount += headerBytes.byteLength;
      if (this.directOutput) await this.directOutput.writable.write(headerBytes);
    } else if (this.headerLine !== normalizedHeader) {
      throw new Error("The aligner returned inconsistent AIRR columns between batches.");
    }

    const bodyBytes = encoder.encode(bodyText);
    const lines = bodyText.split("\n").map((line) => line.replace(/\r$/, "")).filter(Boolean);
    if (!lines.length) return;
    for (const line of lines) updateHash(this.fingerprintState, `${line}\n`);
    const doubleDByLine = new Map<number, Record<string, string>>();
    if (doubleD) {
      const incomingHeader = doubleD.header.replace(/\r$/, "");
      if (!this.doubleDHeaderLine) {
        this.doubleDHeaderLine = incomingHeader;
        this.doubleDHeaders = incomingHeader.split("\t");
      } else if (incomingHeader !== this.doubleDHeaderLine) {
        throw new Error("The aligner returned inconsistent double-D evidence columns between batches.");
      }
      const doubleDText = typeof doubleD.body === "string" ? doubleD.body : decoder.decode(doubleD.body);
      for (const rawLine of doubleDText.split("\n")) {
        const line = rawLine.replace(/\r$/, "");
        if (!line) continue;
        const values = line.split("\t");
        const row = Object.fromEntries(this.doubleDHeaders.map((name, index) => [name, values[index] ?? ""]));
        const relative = Number(row.swig_batch_record_index);
        if (!Number.isInteger(relative) || relative < 0 || relative >= lines.length || doubleDByLine.has(relative)) {
          throw new Error("The aligner returned an invalid double-D batch record index.");
        }
        delete row.swig_batch_record_index;
        doubleDByLine.set(relative, row);
      }
    }
    const chunkIndex = this.nextChunk;
    let chunk: ChunkRecord;
    if (this.directOutput) {
      const start = this.outputByteCount;
      await this.directOutput.writable.write(bodyBytes);
      chunk = { index: chunkIndex, storage: "external", start, length: bodyBytes.byteLength };
    } else {
      let data: Blob;
      let compressed = false;
      if ("CompressionStream" in globalThis) {
        const raw = new Blob([bodyBytes.slice().buffer]);
        data = await new Response(raw.stream().pipeThrough(new CompressionStream("gzip"))).blob();
        compressed = true;
      } else {
        data = new Blob([bodyBytes.slice().buffer]);
      }
      chunk = { index: chunkIndex, storage: "indexed", data, compressed };
    }
    this.outputByteCount += bodyBytes.byteLength;

    const database = await this.database;
    const transaction = database.transaction(["chunks", "records", "doubleD"], "readwrite");
    const chunks = transaction.objectStore("chunks");
    const records = transaction.objectStore("records");
    const doubleDRecords = transaction.objectStore("doubleD");
    chunks.put(chunk);

    const positions = Object.fromEntries(this.headers.map((name, index) => [name, index]));
    const at = (values: string[], name: string) => values[positions[name]] ?? "";
    for (let line = 0; line < lines.length; line += 1) {
      const values = lines[line].split("\t");
      const datasetId = at(values, "swig_dataset_id");
      const savedDataset = this.studyMetadataOverrides.get(datasetId);
      const record: AirrIndexRecord = {
        ordinal: this.nextOrdinal++,
        chunk: chunkIndex,
        line,
        sequenceId: at(values, "sequence_id"),
        datasetId,
        sampleId: savedDataset?.sampleId ?? at(values, "sample_id"),
        subjectId: savedDataset?.subjectId ?? at(values, "subject_id"),
        cohort: savedDataset?.cohort ?? at(values, "swig_cohort"),
        timepoint: savedDataset?.timepoint ?? at(values, "swig_timepoint"),
        compartment: savedDataset?.compartment ?? at(values, "swig_compartment"),
        locus: at(values, "locus"),
        vCall: at(values, "v_call"),
        dCall: doubleDByLine.get(line)?.d_call || at(values, "d_call"),
        d2Call: doubleDByLine.get(line)?.d2_call ?? "",
        jCall: at(values, "j_call"),
        cCall: at(values, "c_call"),
        isotype: at(values, "isotype") || inferIsotype(
          at(values, "c_call"),
          at(values, "c_sequence_alignment"),
          numeric(at(values, "c_identity")),
          numeric(at(values, "c_support")),
        ),
        productive: at(values, "productive"),
        cdr3: at(values, "cdr3"),
        cdr3Aa: at(values, "cdr3_aa"),
        junctionAa: at(values, "junction_aa"),
        vIdentity: numeric(at(values, "v_identity")),
        dIdentity: numeric(doubleDByLine.get(line)?.d_identity ?? at(values, "d_identity")),
        jIdentity: numeric(at(values, "j_identity")),
        cIdentity: numeric(at(values, "c_identity")),
        cdr3AaLength: at(values, "cdr3_aa") ? at(values, "cdr3_aa").length : null,
        vjInFrame: at(values, "vj_in_frame"),
        stopCodon: at(values, "stop_codon"),
        completeVdj: at(values, "complete_vdj"),
        revComp: at(values, "rev_comp"),
      };
      records.put(packRecord(record));
      const doubleDValues = doubleDByLine.get(line);
      if (doubleDValues) {
        doubleDRecords.put({ ordinal: record.ordinal, values: doubleDValues } satisfies DoubleDEvidenceRecord);
        this.doubleDRecordCount += 1;
      }
      if (record.vCall && record.jCall) this.assigned += 1;
      if (record.productive === "T") this.productive += 1;
      if (record.cdr3 || record.cdr3Aa) this.withCdr3 += 1;
      bump(this.facetMaps.datasets, record.datasetId);
      bump(this.facetMaps.samples, record.sampleId);
      bump(this.facetMaps.subjects, record.subjectId);
      bump(this.facetMaps.cohorts, record.cohort);
      bump(this.facetMaps.timepoints, record.timepoint);
      bump(this.facetMaps.compartments, record.compartment);
      bump(this.facetMaps.loci, record.locus);
      bump(this.facetMaps.productive, record.productive);
      bump(this.facetMaps.vCalls, record.vCall);
      bump(this.facetMaps.dCalls, record.dCall);
      bump(this.facetMaps.jCalls, record.jCall);
      bump(this.facetMaps.cCalls, record.cCall);
      bump(this.facetMaps.isotypes, record.isotype);
      this.accumulateRepertoire(record);
    }
    try {
      await transactionDone(transaction);
      this.nextChunk += 1;
    } catch (error) {
      if (error instanceof DOMException && error.name === "QuotaExceededError") {
        throw new Error("Browser storage is full. Re-run with ‘stream AIRR to disk’ so output is never retained in browser storage.");
      }
      throw error;
    }
  }

  async page(
    filters: ResultFilters,
    offset: number,
    limit: number,
    onScan?: (scanned: number) => void,
    signal?: AbortSignal,
  ): Promise<ResultPage> {
    const database = await this.database;
    const transaction = database.transaction("records", "readonly");
    const store = transaction.objectStore("records");
    const normalizedId = filters.sequenceId.trim().toLowerCase();
    const normalizedCdr3 = filters.cdr3.trim().toUpperCase();
    const filtered = Boolean(
      normalizedId || normalizedCdr3 || filters.locus || filters.productive ||
      filters.datasetId || filters.sampleId || filters.subjectId || filters.cohort || filters.timepoint || filters.compartment ||
      filters.vCall || filters.dCall || filters.jCall || filters.cCall || filters.isotype ||
      filters.minVIdentity || filters.minDIdentity || filters.minJIdentity || filters.minCIdentity || filters.minCdr3AaLength ||
      filters.maxCdr3AaLength || filters.vjInFrame || filters.stopCodon ||
      filters.completeVdj || filters.revComp || filters.hasD || filters.hasDoubleD || filters.hasCdr3,
    );

    if (!filtered) {
      if (offset >= this.count) return { rows: [], hasMore: false, totalMatches: this.count, scanned: 0 };
      const upper = Math.min(this.count - 1, offset + limit - 1);
      const packed = await requestResult(store.getAll(IDBKeyRange.bound(offset, upper), limit)) as PackedIndexRecord[];
      const rows = packed.map(unpackRecord);
      return { rows, hasMore: offset + rows.length < this.count, totalMatches: this.count, scanned: rows.length };
    }

    let source: IDBObjectStore | IDBIndex = store;
    let range: IDBKeyRange | undefined;
    const exactCandidates: Array<[keyof ResultFilters, string]> = [
      ["datasetId", "datasetId"], ["sampleId", "sampleId"], ["subjectId", "subjectId"],
      ["cohort", "cohort"], ["timepoint", "timepoint"], ["compartment", "compartment"],
      ["isotype", "isotype"],
      ["locus", "locus"], ["productive", "productive"],
    ];
    if (filters.hasDoubleD) {
      source = store.index("d2Call");
      range = IDBKeyRange.lowerBound("", true);
    } else {
      const callCandidate = ([
        [filters.vCall, "vCallToken"], [filters.dCall, "dCallToken"],
        [filters.jCall, "jCallToken"], [filters.cCall, "cCallToken"],
      ] as Array<[string, string]>).find(([value]) => value.trim()&&!value.includes(","));
      if (callCandidate) {
        source = store.index(callCandidate[1]);
        range = IDBKeyRange.only(callCandidate[0].trim().toUpperCase());
      } else for (const [filterName, indexName] of exactCandidates) {
        const value = filters[filterName];
        if (typeof value === "string" && value) {
          source = store.index(indexName);
          range = IDBKeyRange.only(value);
          break;
        }
      }
    }

    const matches = (record: AirrIndexRecord) => {
      if (normalizedId && !record.sequenceId.toLowerCase().includes(normalizedId)) return false;
      if (normalizedCdr3 && !record.cdr3.toUpperCase().includes(normalizedCdr3) &&
        !record.cdr3Aa.toUpperCase().includes(normalizedCdr3)) return false;
      if (filters.datasetId && record.datasetId !== filters.datasetId) return false;
      if (filters.sampleId && record.sampleId !== filters.sampleId) return false;
      if (filters.subjectId && record.subjectId !== filters.subjectId) return false;
      if (filters.cohort && record.cohort !== filters.cohort) return false;
      if (filters.timepoint && record.timepoint !== filters.timepoint) return false;
      if (filters.compartment && record.compartment !== filters.compartment) return false;
      if (filters.locus && record.locus !== filters.locus) return false;
      if (filters.productive && record.productive !== filters.productive) return false;
      if (!matchesCall(record.vCall, filters.vCall, filters.vCallIncludeAmbiguous)) return false;
      if (!matchesCall(record.dCall, filters.dCall, filters.dCallIncludeAmbiguous)) return false;
      if (!matchesCall(record.jCall, filters.jCall, filters.jCallIncludeAmbiguous)) return false;
      if (!matchesCall(record.cCall, filters.cCall, filters.cCallIncludeAmbiguous)) return false;
      if (filters.isotype && record.isotype !== filters.isotype) return false;
      if (filters.minVIdentity && (record.vIdentity ?? 0) < filters.minVIdentity) return false;
      if (filters.minDIdentity && (record.dIdentity ?? 0) < filters.minDIdentity) return false;
      if (filters.minJIdentity && (record.jIdentity ?? 0) < filters.minJIdentity) return false;
      if (filters.minCIdentity && (record.cIdentity ?? 0) < filters.minCIdentity) return false;
      if (filters.minCdr3AaLength && (record.cdr3AaLength ?? 0) < filters.minCdr3AaLength) return false;
      if (filters.maxCdr3AaLength && (record.cdr3AaLength ?? Number.POSITIVE_INFINITY) > filters.maxCdr3AaLength) return false;
      if (filters.vjInFrame && record.vjInFrame !== filters.vjInFrame) return false;
      if (filters.stopCodon && record.stopCodon !== filters.stopCodon) return false;
      if (filters.completeVdj && record.completeVdj !== filters.completeVdj) return false;
      if (filters.revComp && record.revComp !== filters.revComp) return false;
      if (filters.hasD && !record.dCall) return false;
      if (filters.hasDoubleD && !record.d2Call) return false;
      if (filters.hasCdr3 && !record.cdr3 && !record.cdr3Aa) return false;
      return true;
    };

    return new Promise((resolve, reject) => {
      const rows: AirrIndexRecord[] = [];
      let scanned = 0;
      let matched = 0;
      let resolved = false;
      const request = source.openCursor(range);
      request.onerror = () => reject(request.error ?? new Error("Could not search the local result index."));
      request.onsuccess = () => {
        if (signal?.aborted) {
          if (!resolved) {
            resolved = true;
            resolve({ rows: [], hasMore: false, totalMatches: null, scanned });
          }
          return;
        }
        const cursor = request.result;
        if (!cursor) {
          if (!resolved) resolve({ rows, hasMore: false, totalMatches: matched, scanned });
          return;
        }
        scanned += 1;
        if (scanned % 2500 === 0) onScan?.(scanned);
        const record = unpackRecord(cursor.value as PackedIndexRecord);
        if (matches(record)) {
          if (matched >= offset && rows.length < limit) rows.push(record);
          matched += 1;
          if (rows.length === limit && matched > offset + limit) {
            resolved = true;
            resolve({ rows, hasMore: true, totalMatches: null, scanned });
            return;
          }
        }
        cursor.continue();
      };
    });
  }

  async detail(record: AirrIndexRecord): Promise<AirrRow> {
    const database = await this.database;
    const transaction = database.transaction("chunks", "readonly");
    const chunk = await requestResult(transaction.objectStore("chunks").get(record.chunk)) as ChunkRecord | undefined;
    if (!chunk) throw new Error("That AIRR result is no longer present in the local index.");
    const text = await this.chunkText(chunk);
    const line = text.split("\n")[record.line]?.replace(/\r$/, "") ?? "";
    const values = line.split("\t");
    const row = Object.fromEntries(this.headers.map((header, index) => [header, values[index] ?? ""]));
    const doubleDTransaction = database.transaction("doubleD", "readonly");
    const doubleD = await requestResult(doubleDTransaction.objectStore("doubleD").get(record.ordinal)) as DoubleDEvidenceRecord | undefined;
    return this.overlayStudyMetadata(doubleD ? { ...row, ...doubleD.values } : row);
  }

  async indexRecords(ordinals: readonly number[]): Promise<AirrIndexRecord[]> {
    if (!ordinals.length) return [];
    const database = await this.database;
    const transaction = database.transaction("records", "readonly");
    const store = transaction.objectStore("records");
    const packed = await Promise.all(ordinals.map((ordinal) => requestResult(store.get(ordinal)) as Promise<PackedIndexRecord | undefined>));
    return packed.flatMap((record) => record ? [unpackRecord(record)] : []);
  }

  async detailMany(ordinals: readonly number[]): Promise<AirrDetailRow[]> {
    const records = await this.indexRecords(ordinals);
    const byChunk = new Map<number, AirrIndexRecord[]>();
    for (const record of records) {
      const values = byChunk.get(record.chunk);
      if (values) values.push(record);
      else byChunk.set(record.chunk, [record]);
    }
    const result = new Map<number, AirrDetailRow>();
    const database = await this.database;
    for (const [chunkIndex, chunkRecords] of byChunk) {
      const transaction = database.transaction("chunks", "readonly");
      const chunk = await requestResult(transaction.objectStore("chunks").get(chunkIndex)) as ChunkRecord | undefined;
      if (!chunk) continue;
      const lines = (await this.chunkText(chunk)).split("\n");
      for (const record of chunkRecords) {
        const line = lines[record.line]?.replace(/\r$/, "") ?? "";
        const values = line.split("\t");
        result.set(record.ordinal, {
          record,
          values: this.overlayStudyMetadata(Object.fromEntries(this.headers.map((header, index) => [header, values[index] ?? ""]))),
        });
      }
    }
    if (result.size) {
      const transaction = database.transaction("doubleD", "readonly");
      const doubleDStore = transaction.objectStore("doubleD");
      await Promise.all([...result.entries()].map(async ([ordinal, detail]) => {
        const doubleD = await requestResult(doubleDStore.get(ordinal)) as DoubleDEvidenceRecord | undefined;
        if (doubleD) detail.values = { ...detail.values, ...doubleD.values };
      }));
    }
    return ordinals.flatMap((ordinal) => result.get(ordinal) ?? []);
  }

  async scanAirrRows(
    fields: readonly string[],
    onBatch: (rows: AirrScanRow[]) => void | Promise<void>,
    options: {
      batchSize?: number;
      onProgress?: (processed: number, total: number) => void;
      signal?: AbortSignal;
      includeMask?: Uint8Array;
    } = {},
  ): Promise<void> {
    const selected = [...new Set(fields)];
    const positions = selected.map((field) => this.headers.indexOf(field));
    const metadataDatasetPosition = this.headers.indexOf("swig_dataset_id");
    const batchSize = Math.max(100, options.batchSize ?? 2_000);
    if (options.includeMask && options.includeMask.length !== this.count) {
      throw new Error("The AIRR scan mask does not match the result-store record count.");
    }
    let includedTotal = this.count;
    if (options.includeMask) {
      includedTotal = 0;
      for (const value of options.includeMask) includedTotal += value ? 1 : 0;
    }
    const database = await this.database;
    let ordinal = 0;
    let included = 0;
    let batch: AirrScanRow[] = [];
    for (let index = 0; index < this.nextChunk; index += 1) {
      if (options.signal?.aborted) throw new DOMException("Post-analysis was cancelled.", "AbortError");
      const transaction = database.transaction("chunks", "readonly");
      const chunk = await requestResult(transaction.objectStore("chunks").get(index)) as ChunkRecord | undefined;
      if (!chunk) continue;
      const lines = (await this.chunkText(chunk)).split("\n");
      for (const rawLine of lines) {
        const line = rawLine.replace(/\r$/, "");
        if (!line) continue;
        const rowOrdinal = ordinal++;
        if (options.includeMask && !options.includeMask[rowOrdinal]) continue;
        const values = line.split("\t");
        const row: AirrRow = {};
        selected.forEach((field, fieldIndex) => {
          const position = positions[fieldIndex];
          row[field] = position >= 0 ? values[position] ?? "" : "";
        });
        batch.push({ ordinal: rowOrdinal, values: this.overlayStudyMetadata(row, metadataDatasetPosition >= 0 ? values[metadataDatasetPosition] ?? "" : "") });
        included += 1;
        if (batch.length >= batchSize) {
          await onBatch(batch);
          batch = [];
          options.onProgress?.(included, includedTotal);
          await cooperativeYield();
          if (options.signal?.aborted) throw new DOMException("Post-analysis was cancelled.", "AbortError");
        }
      }
    }
    if (batch.length) await onBatch(batch);
    options.onProgress?.(included, includedTotal);
  }

  private async chunkBlob(chunk: ChunkRecord): Promise<Blob> {
    if (chunk.storage === "external") {
      if (!this.directOutput || chunk.start === undefined || chunk.length === undefined) {
        throw new Error("The streamed AIRR output file is no longer available.");
      }
      const file = await this.directOutput.handle.getFile();
      return file.slice(chunk.start, chunk.start + chunk.length);
    }
    if (!chunk.data) throw new Error("An indexed AIRR output batch is missing.");
    if (!chunk.compressed) return chunk.data;
    if (!("DecompressionStream" in globalThis)) {
      throw new Error("This browser cannot decompress the local AIRR result batch.");
    }
    return new Response(chunk.data.stream().pipeThrough(new DecompressionStream("gzip"))).blob();
  }

  private async chunkText(chunk: ChunkRecord): Promise<string> {
    return (await this.chunkBlob(chunk)).text();
  }

  async finalize(): Promise<void> {
    if (this.finalized) return;
    const database = await this.database;
    const transaction = database.transaction("meta", "readwrite");
    transaction.objectStore("meta").put({
      key: "manifest",
      headerLine: this.headerLine,
      chunks: this.nextChunk,
      records: this.nextOrdinal,
      outputBytes: this.outputByteCount,
    } satisfies ManifestRecord);
    await transactionDone(transaction);
    await this.directOutput?.writable.close();
    this.finalized = true;
  }

  async abort(): Promise<void> {
    if (!this.finalized) {
      try {
        await this.directOutput?.writable.abort?.();
      } catch {
        // The browser may already have discarded an interrupted output stream.
      }
    }
  }

  streamingDownloadUrl(baseUrl: string, name: string, compression: "none"|"gzip" = "none"): string {
    const url = new URL(`${baseUrl}__swig_download__`, globalThis.location?.href);
    url.searchParams.set("database", this.databaseName);
    url.searchParams.set("name", name);
    if(compression==="gzip")url.searchParams.set("compression","gzip");
    return url.href;
  }

  async airrBlob(signal?: AbortSignal): Promise<Blob> {
    if (signal?.aborted) throw new DOMException("AIRR export was cancelled.", "AbortError");
    if (this.directOutput && !this.studyMetadataOverrides.size) {
      const file = await this.directOutput.handle.getFile();
      if (signal?.aborted) throw new DOMException("AIRR export was cancelled.", "AbortError");
      return file;
    }
    if (this.outputByteCount > MAX_FALLBACK_BLOB_BYTES) {
      throw new Error("This AIRR table is too large for a memory-backed download. Enable the streaming download worker or re-run with direct-to-disk output.");
    }
    const parts: BlobPart[] = [];
    await this.writeAirr(async (part) => { parts.push(part instanceof Uint8Array ? part.slice().buffer : part); }, signal);
    return new Blob(parts, { type: "text/tab-separated-values;charset=utf-8" });
  }

  async writeAirr(
    write: (part: string | Blob | Uint8Array) => Promise<void>,
    signal?: AbortSignal,
  ): Promise<void> {
    if (signal?.aborted) throw new DOMException("AIRR export was cancelled.", "AbortError");
    if (this.directOutput && !this.studyMetadataOverrides.size) {
      await write(await this.directOutput.handle.getFile());
      if (signal?.aborted) throw new DOMException("AIRR export was cancelled.", "AbortError");
      return;
    }
    if (this.studyMetadataOverrides.size) {
      await write(`${this.headerLine}\n`);
      await this.scanAirrRows(this.headers, async (rows) => {
        let body = "";
        for (const row of rows) body += tableRow(this.headers, row.values, "tsv");
        if (body) await write(body);
      }, { batchSize: 2_000, signal });
      return;
    }
    await write(`${this.headerLine}\n`);
    const database = await this.database;
    for (let index = 0; index < this.nextChunk; index += 1) {
      if (signal?.aborted) throw new DOMException("AIRR export was cancelled.", "AbortError");
      const transaction = database.transaction("chunks", "readonly");
      const chunk = await requestResult(transaction.objectStore("chunks").get(index)) as ChunkRecord | undefined;
      if (chunk) await write(await this.chunkBlob(chunk));
    }
  }

  async writeAirrFormat(
    format: TableExportFormat,
    write: (part: string | Blob | Uint8Array) => Promise<void>,
    includeMask?: Uint8Array,
    signal?: AbortSignal,
  ): Promise<void> {
    if (format === "tsv" && !includeMask) return this.writeAirr(write, signal);
    if (includeMask && includeMask.length !== this.count) throw new Error("The export mask does not match the result-store record count.");
    const header = tableHeader(this.headers, format);
    if (header) await write(header);
    await this.scanAirrRows(this.headers, async (rows) => {
      let body = "";
      for (const row of rows) body += tableRow(this.headers, row.values, format);
      if (body) await write(body);
    }, { batchSize: 2000, includeMask, signal });
  }

  async doubleDRecords(signal?: AbortSignal): Promise<DoubleDEvidenceRecord[]> {
    if (!this.doubleDRecordCount) return [];
    if (signal?.aborted) throw new DOMException("Double-D export was cancelled.", "AbortError");
    const database = await this.database;
    const transaction = database.transaction("doubleD", "readonly");
    const records = await requestResult(transaction.objectStore("doubleD").getAll()) as DoubleDEvidenceRecord[];
    if (signal?.aborted) throw new DOMException("Double-D export was cancelled.", "AbortError");
    if (!this.studyMetadataOverrides.size) return records;
    return records.map((record) => ({ ...record, values: this.overlayStudyMetadata(record.values) }));
  }

  async doubleDMask(): Promise<Uint8Array> {
    const mask = new Uint8Array(this.count);
    for (const record of await this.doubleDRecords()) if (record.values.d2_call) mask[record.ordinal] = 1;
    return mask;
  }

  async importDoubleDRecords(records: DoubleDEvidenceRecord[]): Promise<void> {
    if (!records.length) return;
    for (const evidence of records) if (!Number.isInteger(evidence.ordinal) || evidence.ordinal < 0 || evidence.ordinal >= this.count) {
      throw new Error("The session contains a double-D ordinal outside the linked AIRR table.");
    }
    const indexed = new Map((await this.indexRecords(records.map((record) => record.ordinal))).map((record) => [record.ordinal, record]));
    const database = await this.database;
    const transaction = database.transaction(["doubleD", "records"], "readwrite");
    const evidenceStore = transaction.objectStore("doubleD");
    const recordStore = transaction.objectStore("records");
    const headers = new Set<string>();
    for (const evidence of records) {
      Object.keys(evidence.values).forEach((header) => headers.add(header));
      evidenceStore.put(evidence);
      const current = indexed.get(evidence.ordinal);
      if (current) recordStore.put(packRecord({ ...current, dCall: evidence.values.d_call || current.dCall, d2Call: evidence.values.d2_call ?? "", dIdentity: numeric(evidence.values.d_identity) ?? current.dIdentity }));
    }
    await transactionDone(transaction);
    this.doubleDHeaders = [...headers];
    this.doubleDHeaderLine = this.doubleDHeaders.join("\t");
    this.doubleDRecordCount = records.length;
  }

  async writeDoubleD(write: (part: string | Blob | Uint8Array) => Promise<void>): Promise<void> {
    return this.writeDoubleDFormat("tsv", write);
  }

  async writeDoubleDFormat(
    format: TableExportFormat,
    write: (part: string | Blob | Uint8Array) => Promise<void>,
    signal?: AbortSignal,
  ): Promise<void> {
    if (!this.doubleDHeaderLine || !this.doubleDRecordCount) {
      throw new Error("This analysis has no supported double-D calls to export.");
    }
    const exportedHeaders = this.doubleDHeaders.filter((header) => header !== "swig_batch_record_index");
    const fields = ["swig_airr_ordinal", ...exportedHeaders];
    const header = tableHeader(fields, format);
    if (header) await write(header);
    const records = await this.doubleDRecords(signal);
    let body = "";
    for (const record of records) {
      if (signal?.aborted) throw new DOMException("Double-D export was cancelled.", "AbortError");
      body += tableRow(fields, { swig_airr_ordinal: record.ordinal + 1, ...record.values }, format);
      if (body.length >= 1_000_000) {
        await write(body);
        body = "";
        await cooperativeYield();
      }
    }
    if (body) await write(body);
  }

  async doubleDTsv(): Promise<string> {
    const parts: string[] = [];
    await this.writeDoubleDFormat("tsv", async (part) => { parts.push(typeof part === "string" ? part : await new Response(part instanceof Blob ? part : new Blob([part.slice().buffer])).text()); });
    return parts.join("");
  }

  async writeDeduplicatedAirr(
    counts: Uint32Array,
    write: (part: string | Blob | Uint8Array) => Promise<void>,
  ): Promise<void> {
    return this.writeDeduplicatedAirrFormat(counts, "tsv", write);
  }

  async writeDeduplicatedAirrFormat(
    counts: Uint32Array,
    format: TableExportFormat,
    write: (part: string | Blob | Uint8Array) => Promise<void>,
    signal?: AbortSignal,
  ): Promise<void> {
    if (counts.length < this.count) throw new Error("The duplicate-count vector does not cover every AIRR record.");
    const fields = [...this.headers, "duplicate_count"];
    const header = tableHeader(fields, format);
    if (header) await write(header);
    await this.scanAirrRows(this.headers, async (rows) => {
      let body = "";
      for (const row of rows) if (counts[row.ordinal]) body += tableRow(fields, { ...row.values, duplicate_count: counts[row.ordinal] }, format);
      if (body) await write(body);
    }, { batchSize: 2000, signal });
  }

  async writeLineageAirr(
    assignments: Int32Array,
    write: (part: string | Blob | Uint8Array) => Promise<void>,
  ): Promise<void> {
    return this.writeLineageAirrFormat(assignments, "tsv", write);
  }

  async writeLineageAirrFormat(
    assignments: Int32Array,
    format: TableExportFormat,
    write: (part: string | Blob | Uint8Array) => Promise<void>,
    mergedLineageByOriginal?: ReadonlyMap<number, string>,
    signal?: AbortSignal,
  ): Promise<void> {
    if (assignments.length < this.count) throw new Error("The lineage-assignment vector does not cover every AIRR record.");
    const clonePosition = this.headers.indexOf("clone_id");
    const fields = clonePosition >= 0 ? [...this.headers] : [...this.headers, "clone_id"];
    if (mergedLineageByOriginal?.size && !fields.includes("swig_merged_lineage_id")) fields.push("swig_merged_lineage_id");
    const header = tableHeader(fields, format);
    if (header) await write(header);
    await this.scanAirrRows(this.headers, async (rows) => {
      let body = "";
      for (const row of rows) {
        const lineage = assignments[row.ordinal];
        const cloneId = lineage > 0 ? `swig_lineage_${lineage}` : "";
        const mergedLineageId = lineage > 0 ? mergedLineageByOriginal?.get(lineage) ?? "" : "";
        body += tableRow(fields, { ...row.values, clone_id: cloneId, swig_merged_lineage_id: mergedLineageId }, format);
      }
      if (body) await write(body);
    }, { batchSize: 2000, signal });
  }

  async clear(): Promise<void> {
    await this.abort();
    const database = await this.database;
    database.close();
    await new Promise<void>((resolve) => {
      const request = indexedDB.deleteDatabase(this.databaseName);
      request.onsuccess = () => resolve();
      request.onerror = () => resolve();
      request.onblocked = () => resolve();
    });
  }
}
