import type { AirrResultStore, AirrScanRow } from "./result-store";
import type {
  CollapseMode,
  DedupKey,
  DenoiseOptions,
  ExpansionOptions,
  LineageOptions,
  LineageNeighbourOptions,
  LineageNeighbourResult,
  LineageSummary,
  QueryHit,
  QueryOptions,
} from "./post-analysis-core";
import type { DatasetScope } from "./study-design";
import type { AlleleReassignmentPolicy, AlleleRefinementResult, RefinementSegment } from "./allele-refinement/types";
import { denoiseVdjSequence } from "./post-analysis-record";

export interface DedupDashboard {
  mode: CollapseMode;
  key: DedupKey;
  algorithm: string;
  inputRecords: number;
  inputAbundance: number;
  uniqueRecords: number;
  collapsedRecords: number;
  largestGroups: Array<{ ordinal: number; count: number }>;
  partitions: number;
  candidateComparisons: number;
  indelMergedVariants: number;
  substitutionMergedVariants: number;
  excludedAmbiguous: number;
  unresolvedRecords: number;
  warnings: string[];
}

export interface LineageDashboard {
  summaries: LineageSummary[];
  lineageCount: number;
  sizeHistogram: Array<{ label: string; count: number }>;
  vUsage: Array<{ call: string; lineages: number; abundance: number }>;
  jUsage: Array<{ call: string; lineages: number; abundance: number }>;
  assignedRecords: number;
  unassignedRecords: number;
  candidateComparisons: number;
  truncatedCandidates: number;
}

type WorkerResult = Record<string, unknown>;

interface CallOverridePayload {
  reassignmentPolicy: AlleleReassignmentPolicy;
  minimumPosterior: number;
  v?: { labels: string[]; mapNode: Int32Array; probability: Float32Array };
  j?: { labels: string[]; mapNode: Int32Array; probability: Float32Array };
}

interface RuntimeCheckpoint {
  activeMask: Uint8Array | null;
  dedup?: { dashboard: DedupDashboard; counts: Uint32Array; representatives: Int32Array };
  lineages?: { dashboard: LineageDashboard; assignments: Int32Array };
  callOverrides: CallOverridePayload | null;
}

function abortError(message = "Post-analysis was cancelled."): DOMException {
  return new DOMException(message, "AbortError");
}

function cloneCheckpoint(checkpoint: RuntimeCheckpoint): RuntimeCheckpoint {
  return {
    activeMask: checkpoint.activeMask?.slice() ?? null,
    dedup: checkpoint.dedup ? {
      dashboard: checkpoint.dedup.dashboard,
      counts: checkpoint.dedup.counts.slice(),
      representatives: checkpoint.dedup.representatives.slice(),
    } : undefined,
    lineages: checkpoint.lineages ? {
      dashboard: checkpoint.lineages.dashboard,
      assignments: checkpoint.lineages.assignments.slice(),
    } : undefined,
    callOverrides: checkpoint.callOverrides ? {
      reassignmentPolicy: checkpoint.callOverrides.reassignmentPolicy,
      minimumPosterior: checkpoint.callOverrides.minimumPosterior,
      v: checkpoint.callOverrides.v ? {
        labels: [...checkpoint.callOverrides.v.labels],
        mapNode: checkpoint.callOverrides.v.mapNode.slice(),
        probability: checkpoint.callOverrides.v.probability.slice(),
      } : undefined,
      j: checkpoint.callOverrides.j ? {
        labels: [...checkpoint.callOverrides.j.labels],
        mapNode: checkpoint.callOverrides.j.mapNode.slice(),
        probability: checkpoint.callOverrides.j.probability.slice(),
      } : undefined,
    } : null,
  };
}

export class PostAnalysisRuntime {
  private worker!: Worker;
  private readonly pending = new Map<number, { resolve: (value: unknown) => void; reject: (error: Error) => void; onProgress?: (progress: { processed: number; total: number; phase: string }) => void }>();
  private nextId = 1;
  private indexed = false;
  private indexing: Promise<void> | null = null;
  private sketched = false;
  private sketching: Promise<void> | null = null;
  private generation = 0;
  private closed = false;
  private checkpoint: RuntimeCheckpoint = { activeMask: null, callOverrides: null };
  private transactionCheckpoint: RuntimeCheckpoint | null = null;
  private transactionDirty = false;
  private recovery: Promise<void> | null = null;

  constructor(private readonly store: AirrResultStore, private readonly parallelWorkers = 1) {
    this.installWorker();
  }

  private installWorker() {
    const generation = this.generation;
    this.worker = new Worker(new URL("./post-analysis-worker.ts", import.meta.url), { type: "module" });
    this.worker.onmessage = (event: MessageEvent<{ id: number; result?: unknown; error?: string; progress?: { processed: number; total: number; phase: string } }>) => {
      if (generation !== this.generation) return;
      const pending = this.pending.get(event.data.id);
      if (!pending) return;
      if (event.data.progress) {
        pending.onProgress?.(event.data.progress);
        return;
      }
      this.pending.delete(event.data.id);
      if (event.data.error) pending.reject(new Error(event.data.error));
      else pending.resolve(event.data.result);
    };
    this.worker.onerror = (event) => {
      if (generation !== this.generation) return;
      const error = new Error(event.message || "The post-analysis worker stopped unexpectedly.");
      this.pending.forEach((pending) => pending.reject(error));
      this.pending.clear();
    };
  }

  private replaceWorker(error: Error) {
    this.generation += 1;
    this.worker.terminate();
    this.pending.forEach((pending) => pending.reject(error));
    this.pending.clear();
    this.indexed = false;
    this.indexing = null;
    this.sketched = false;
    this.sketching = null;
    if (!this.closed) this.installWorker();
  }

  private request<T>(message: Record<string, unknown>, onProgress?: (progress: { processed: number; total: number; phase: string }) => void): Promise<T> {
    if (this.closed) return Promise.reject(new Error("The post-analysis runtime is closed."));
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(id, { resolve: resolve as (value: unknown) => void, reject, onProgress });
      this.worker.postMessage({ id, ...message });
    });
  }

  ensureIndexed(onProgress?: (processed: number, total: number) => void, signal?: AbortSignal): Promise<void> {
    if (this.indexed) return Promise.resolve();
    if (this.indexing) return this.indexing;
    const generation = this.generation;
    const indexing = (async () => {
      if (signal?.aborted || generation !== this.generation) throw abortError();
      const doubleDMask = await this.store.doubleDMask();
      if (signal?.aborted || generation !== this.generation) throw abortError();
      await this.request({ type: "init", total: this.store.count, doubleDMask });
      const fields = [
        "sequence_id", "sequence", "sequence_alignment", "locus", "v_call", "j_call", "c_call",
        "cdr3", "cdr3_aa", "productive", "duplicate_count", "swig_dataset_id", "sample_id",
        "subject_id", "swig_cohort", "swig_timepoint", "swig_compartment",
      ];
      await this.store.scanAirrRows(fields, async (rows) => {
        if (signal?.aborted || generation !== this.generation) throw abortError();
        await this.request({ type: "ingest", rows: rows.map(toWorkerRow) });
      }, { batchSize: 2_000, onProgress, signal });
      if (signal?.aborted || generation !== this.generation) throw abortError();
      this.indexed = true;
    })();
    this.indexing = indexing;
    void indexing.finally(() => { if (this.indexing === indexing) this.indexing = null; }).catch(() => undefined);
    return this.indexing;
  }

  async deduplicate(key: DedupKey, unresolvedPolicy: "discard" | "retain" = "discard", scope: DatasetScope = "global", respectConstantCall = true): Promise<DedupDashboard> {
    await this.ensureIndexed();
    const dashboard = await this.request<DedupDashboard>({ type: "dedup", key, unresolvedPolicy, scope, respectConstantCall, workers: this.parallelWorkers });
    const state = await this.request<{ counts: Uint32Array; representatives: Int32Array }>({ type: "dedupState" });
    this.checkpoint = { ...this.checkpoint, activeMask: null, dedup: { dashboard, counts: state.counts, representatives: state.representatives }, lineages: undefined };
    this.transactionDirty ||= Boolean(this.transactionCheckpoint);
    return dashboard;
  }

  async denoise(options: DenoiseOptions, onProgress?: (processed: number, total: number, phase?: "ingest" | "variants" | "finalize") => void, signal?: AbortSignal): Promise<DedupDashboard> {
    await this.ensureIndexed(onProgress, signal);
    await this.request({ type: "denoiseInit", options });
    await this.store.scanAirrRows(["sequence_alignment", "sequence", "rev_comp", "v_sequence_start", "j_sequence_end"], async (rows) => {
      if (signal?.aborted) throw new DOMException("Denoising was cancelled.", "AbortError");
      await this.request({
        type: "denoiseIngest",
        rows: rows.map((row) => ({ ordinal: row.ordinal, sequence: denoiseVdjSequence(row) })),
      });
    }, { batchSize: 2_000, onProgress: onProgress ? (processed, total) => onProgress(processed, total, "ingest") : undefined, signal });
    const dashboard = await this.request<DedupDashboard>({ type: "denoiseFinish", workers: this.parallelWorkers }, (progress) => {
      const phase = progress.phase === "finalize" ? "finalize" : "variants";
      onProgress?.(progress.processed, progress.total, phase);
    });
    const state = await this.request<{ counts: Uint32Array; representatives: Int32Array }>({ type: "dedupState" });
    this.checkpoint = { ...this.checkpoint, activeMask: null, dedup: { dashboard, counts: state.counts, representatives: state.representatives }, lineages: undefined };
    this.transactionDirty ||= Boolean(this.transactionCheckpoint);
    return dashboard;
  }

  async applyDedupFilter(): Promise<{ mask: Uint8Array; retained: number }> {
    await this.ensureIndexed();
    const result = await this.request<{ mask: Uint8Array; retained: number }>({ type: "applyDedupFilter" });
    this.checkpoint = { ...this.checkpoint, activeMask: result.mask.slice(), lineages: undefined };
    this.transactionDirty ||= Boolean(this.transactionCheckpoint);
    return result;
  }

  async setActiveMask(mask: Uint8Array | null): Promise<{ retained: number }> {
    await this.ensureIndexed();
    const result = await this.request<{ retained: number }>({ type: "setActiveMask", mask });
    this.checkpoint = { ...this.checkpoint, activeMask: mask?.slice() ?? null, lineages: undefined };
    this.transactionDirty ||= Boolean(this.transactionCheckpoint);
    return result;
  }

  async activeMask(): Promise<Uint8Array | null> {
    await this.ensureIndexed();
    const result = await this.request<{ mask: Uint8Array | null }>({ type: "activeMask" });
    return result.mask;
  }

  async setRepertoireCallOverrides(
    refinement: AlleleRefinementResult | null,
    reassignmentPolicy: AlleleReassignmentPolicy = "confidence",
    minimumPosterior = 0.8,
  ): Promise<{ changedV: number; changedJ: number; policy: AlleleReassignmentPolicy; threshold: number }> {
    await this.ensureIndexed();
    const payload = (segment: RefinementSegment) => {
      const result = refinement?.segments[segment];
      return result ? {
        labels: result.nodes.map((node) => node.names.join(",")),
        mapNode: result.mapNode,
        probability: result.mapProbability,
      } : undefined;
    };
    const callOverrides: CallOverridePayload = {
      reassignmentPolicy,
      minimumPosterior,
      v: payload("V"),
      j: payload("J"),
    };
    const result = await this.request<{ changedV: number; changedJ: number; policy: AlleleReassignmentPolicy; threshold: number }>({
      type: "setCallOverrides",
      ...callOverrides,
    });
    this.checkpoint = { ...this.checkpoint, callOverrides: refinement ? cloneCheckpoint({ activeMask: null, callOverrides }).callOverrides : null, lineages: undefined };
    this.transactionDirty ||= Boolean(this.transactionCheckpoint);
    return result;
  }

  async assignLineages(options: LineageOptions, useDedup: boolean): Promise<LineageDashboard> {
    await this.ensureIndexed();
    const dashboard = await this.request<LineageDashboard>({ type: "lineages", options, useDedup });
    const state = await this.request<{ assignments: Int32Array }>({ type: "lineageAssignments" });
    this.checkpoint = { ...this.checkpoint, lineages: { dashboard, assignments: state.assignments } };
    this.transactionDirty ||= Boolean(this.transactionCheckpoint);
    return dashboard;
  }

  async query(queries: string[], options: QueryOptions): Promise<QueryHit[]> {
    await this.ensureIndexed();
    if (options.target === "trimmed") await this.ensureSketches();
    const result = await this.request<{ hits: QueryHit[] }>({ type: "query", queries, options });
    return result.hits;
  }

  private ensureSketches(): Promise<void> {
    if (this.sketched) return Promise.resolve();
    if (this.sketching) return this.sketching;
    const generation = this.generation;
    const sketching = (async () => {
      if (generation !== this.generation) throw abortError();
      await this.request({ type: "initSketches" });
      await this.store.scanAirrRows(["sequence_alignment", "sequence"], async (rows) => {
        if (generation !== this.generation) throw abortError();
        await this.request({
          type: "ingestSketches",
          rows: rows.map((row) => ({ ordinal: row.ordinal, sequence: row.values.sequence_alignment || row.values.sequence })),
        });
      }, { batchSize: 2_000 });
      if (generation !== this.generation) throw abortError();
      this.sketched = true;
    })();
    this.sketching = sketching;
    void sketching.finally(() => { if (this.sketching === sketching) this.sketching = null; }).catch(() => undefined);
    return this.sketching;
  }

  async expand(seedOrdinals: number[], options: ExpansionOptions): Promise<{ ordinals: number[]; comparisons: number; capped: boolean }> {
    await this.ensureIndexed();
    return this.request({ type: "expand", seedOrdinals, options });
  }

  async lineageMembers(lineageId: number, offset = 0, limit = 500): Promise<{ ordinals: number[]; total: number }> {
    return this.request({ type: "lineageMembers", lineageId, offset, limit });
  }

  async lineageMembersMany(lineageIds: number[], limitPerLineage = 500): Promise<Array<{ lineageId: number; ordinals: number[]; total: number }>> {
    const result = await this.request<{ members: Array<{ lineageId: number; ordinals: number[]; total: number }> }>({
      type: "lineageMembersMany",
      lineageIds,
      limitPerLineage,
    });
    return result.members;
  }

  async lineageNeighbours(options: LineageNeighbourOptions, useDedup: boolean): Promise<LineageNeighbourResult> {
    await this.ensureIndexed();
    return this.request<LineageNeighbourResult>({ type: "lineageNeighbours", options, useDedup });
  }

  async lineageAssignments(): Promise<Int32Array> {
    const result = await this.request<{ assignments: Int32Array }>({ type: "lineageAssignments" });
    return result.assignments;
  }

  async dedupMembers(representative: number, offset = 0, limit = 500): Promise<{ ordinals: number[]; total: number }> {
    return this.request({ type: "dedupMembers", representative, offset, limit });
  }

  async dedupCounts(): Promise<Uint32Array> {
    const result = await this.request<{ counts: Uint32Array }>({ type: "dedupCounts" });
    return result.counts;
  }

  async dedupState(): Promise<{ counts: Uint32Array; representatives: Int32Array }> {
    return this.request({ type: "dedupState" });
  }

  async restoreState(state: {
    activeMask?: Uint8Array | null;
    dedup?: { dashboard: DedupDashboard; counts: Uint32Array; representatives: Int32Array };
    lineages?: { dashboard: LineageDashboard; assignments: Int32Array };
  }): Promise<{ lineages: LineageDashboard | null }> {
    await this.ensureIndexed();
    const result = await this.request<{ lineages: LineageDashboard | null }>({ type: "restoreState", ...state });
    this.checkpoint = {
      ...this.checkpoint,
      activeMask: state.activeMask?.slice() ?? null,
      dedup: state.dedup ? { dashboard: state.dedup.dashboard, counts: state.dedup.counts.slice(), representatives: state.dedup.representatives.slice() } : undefined,
      lineages: state.lineages ? { dashboard: state.lineages.dashboard, assignments: state.lineages.assignments.slice() } : undefined,
    };
    this.transactionDirty ||= Boolean(this.transactionCheckpoint);
    return result;
  }

  beginTransaction() {
    if (!this.transactionCheckpoint) {
      this.transactionCheckpoint = cloneCheckpoint(this.checkpoint);
      this.transactionDirty = false;
    }
  }

  commitTransaction() {
    this.transactionCheckpoint = null;
    this.transactionDirty = false;
  }

  requiresRecoveryForCancellation(): boolean {
    return this.pending.size > 0 || Boolean(this.indexing) || Boolean(this.sketching) || this.transactionDirty;
  }

  async rollbackTransaction(onProgress?: (processed: number, total: number) => void): Promise<void> {
    await this.cancelAndRestore(onProgress, true);
  }

  cancelAndRestore(onProgress?: (processed: number, total: number) => void, rollbackTransaction = true): Promise<void> {
    if (this.closed) return Promise.resolve();
    if (this.recovery) return this.recovery;
    if (rollbackTransaction && this.transactionCheckpoint) {
      this.checkpoint = this.transactionCheckpoint;
      this.transactionCheckpoint = null;
      this.transactionDirty = false;
    }
    const target = cloneCheckpoint(this.checkpoint);
    this.replaceWorker(abortError());
    this.recovery = (async () => {
      await this.ensureIndexed(onProgress);
      if (target.callOverrides) await this.request({ type: "setCallOverrides", ...target.callOverrides });
      await this.request({
        type: "restoreState",
        activeMask: target.activeMask,
        dedup: target.dedup,
        lineages: target.lineages,
      });
      this.checkpoint = target;
    })().finally(() => { this.recovery = null; });
    return this.recovery;
  }

  terminate() {
    this.closed = true;
    this.generation += 1;
    this.worker.terminate();
    const error = new Error("Post-analysis was closed.");
    this.pending.forEach((pending) => pending.reject(error));
    this.pending.clear();
  }
}

function toWorkerRow(row: AirrScanRow): WorkerResult {
  return {
    ordinal: row.ordinal,
    sequence_id: row.values.sequence_id,
    swig_dataset_id: row.values.swig_dataset_id,
    sample_id: row.values.sample_id,
    subject_id: row.values.subject_id,
    swig_cohort: row.values.swig_cohort,
    swig_timepoint: row.values.swig_timepoint,
    swig_compartment: row.values.swig_compartment,
    sequence: row.values.sequence,
    sequence_alignment: row.values.sequence_alignment,
    locus: row.values.locus,
    v_call: row.values.v_call,
    j_call: row.values.j_call,
    c_call: row.values.c_call,
    cdr3: row.values.cdr3,
    cdr3_aa: row.values.cdr3_aa,
    productive: row.values.productive,
    duplicate_count: row.values.duplicate_count,
  };
}
