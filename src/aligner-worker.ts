/// <reference lib="webworker" />

import { assignmentCompletion } from "./assignment-progress";
import {
  emptyFastqQualityFilterStats,
  streamSequenceBatches,
  type FastqQualityFilterOptions,
  type FastqQualityFilterStats,
  type SequenceBatch,
  type SequenceFormat,
  type SequenceSource,
  sequenceSourceSize,
} from "./sequence-stream";
import type { AlignerSettings, GermlineTexts } from "./aligner-engine";
import {
  ASSIGNMENT_TELEMETRY_INTERVAL_MS,
  initializingAssignmentTelemetry,
  updateQueriesPerSecondEma,
  type AssignmentPhase,
  type AssignmentTelemetry,
} from "./assignment-telemetry";

interface StartRequest {
  type: "start";
  id: number;
  query: SequenceSource;
  format: SequenceFormat;
  references: GermlineTexts;
  settings: AlignerSettings;
  workers: number;
  countHint: number | null;
  subsample?: { size: number; seed: number };
  fastqFilter?: FastqQualityFilterOptions;
}

interface AckRequest {
  type: "ack";
  id: number;
  batch: number;
}

interface ComputeBatch {
  batch: number;
  header: string;
  body: ArrayBuffer;
  count: number;
  milliseconds: number;
  doubleDHeader?: string;
  doubleDBody?: ArrayBuffer;
  doubleDCount?: number;
}

interface ComputeSlot {
  index: number;
  worker: Worker;
  busy: boolean;
  /** The aligner's database and run warnings, from loading the germlines. */
  warnings: string[];
  resolve?: (batch: ComputeBatch) => void;
  reject?: (error: Error) => void;
}

const acknowledgements = new Map<string, () => void>();
const encoder = new TextEncoder();

function postProgress(id: number, stage: string, value: number) {
  self.postMessage({ id, type: "progress", stage, value });
}

function postInitializationTelemetry(id: number, totalWorkers: number, activeWorkers: number) {
  self.postMessage({
    id,
    type: "telemetry",
    telemetry: {
      ...initializingAssignmentTelemetry(totalWorkers),
      activeWorkers: Math.max(0, Math.min(totalWorkers, activeWorkers)),
    },
  });
}

function waitForAcknowledgement(id: number, batch: number): Promise<void> {
  return new Promise((resolve) => acknowledgements.set(`${id}:${batch}`, resolve));
}

async function compileCore(): Promise<WebAssembly.Module> {
  const url = `${import.meta.env.BASE_URL}swig.wasm`;
  try {
    return await WebAssembly.compileStreaming(fetch(url));
  } catch {
    const response = await fetch(url);
    if (!response.ok) throw new Error("The Swig aligner's WebAssembly module could not be loaded.");
    return WebAssembly.compile(await response.arrayBuffer());
  }
}

function effectiveWorkerCount(request: StartRequest): number {
  const requested = Math.max(1, Math.min(16, Math.floor(request.workers || 1)));
  const effectiveCount = request.subsample
    ? Math.min(request.countHint ?? request.subsample.size, request.subsample.size)
    : request.countHint;
  if (effectiveCount !== null && effectiveCount <= 3) return 1;
  if (effectiveCount !== null && effectiveCount <= 500) return Math.min(2, requested);
  return requested;
}

function batchSize(countHint: number | null): number {
  if (countHint !== null && countHint <= 10) return 1;
  if (countHint !== null && countHint <= 1000) return 100;
  return 1000;
}

async function initializeSlot(
  index: number,
  module: WebAssembly.Module,
  references: GermlineTexts,
  settings: AlignerSettings,
): Promise<ComputeSlot> {
  const worker = new Worker(new URL("./aligner-compute-worker.ts", import.meta.url), { type: "module" });
  const slot: ComputeSlot = { index, worker, busy: false, warnings: [] };
  return new Promise((resolve, reject) => {
    const fail = (error: Error) => {
      worker.terminate();
      reject(error);
    };
    worker.onerror = (event) => fail(new Error(event.message || `Compute worker ${index + 1} stopped unexpectedly.`));
    worker.onmessage = (event) => {
      const message = event.data as { type: "ready" | "error"; records?: number; warnings?: string[]; message?: string };
      if (message.type === "error") {
        fail(new Error(message.message || `Compute worker ${index + 1} could not initialize.`));
        return;
      }
      slot.warnings = message.warnings ?? [];
      resolve(slot);
    };
    worker.postMessage({ type: "initialize", worker: index, module, references, settings });
  });
}

function installComputeHandler(slot: ComputeSlot, onFatal: (error: Error) => void) {
  slot.worker.onerror = (event) => {
    const error = new Error(event.message || `Compute worker ${slot.index + 1} stopped unexpectedly.`);
    slot.reject?.(error);
    onFatal(error);
  };
  slot.worker.onmessage = (event) => {
    const message = event.data as ComputeBatch & { type: "batch" | "error"; message?: string };
    if (message.type === "error") {
      const error = new Error(message.message || `Compute worker ${slot.index + 1} failed.`);
      slot.reject?.(error);
      slot.resolve = undefined;
      slot.reject = undefined;
      onFatal(error);
      return;
    }
    slot.resolve?.(message);
    slot.resolve = undefined;
    slot.reject = undefined;
  };
}

function annotate(slot: ComputeSlot, batch: SequenceBatch, request: StartRequest): Promise<ComputeBatch> {
  const query = encoder.encode(batch.text);
  return new Promise((resolve, reject) => {
    slot.resolve = resolve;
    slot.reject = reject;
    slot.worker.postMessage({
      type: "annotate",
      batch: batch.index,
      query: query.buffer,
      count: batch.count,
      doubleD: Boolean(request.settings.doubleD && request.settings.doubleD.mode !== "off"),
    }, [query.buffer]);
  });
}

async function handleRequest(request: StartRequest) {
  const slots: ComputeSlot[] = [];
  let telemetryTimer: number | undefined;
  try {
    const workerCount = effectiveWorkerCount(request);
    postProgress(request.id, `Loading WebAssembly for ${workerCount} parallel worker${workerCount === 1 ? "" : "s"}`, 0.03);
    postInitializationTelemetry(request.id, workerCount, 0);
    const module = await compileCore();
    postProgress(request.id, `Indexing germlines in ${workerCount} worker${workerCount === 1 ? "" : "s"}`, 0.08);
    // Each compute worker receives its own aligner instance and germline index.
    // The initializations are launched together, which is the expected brief
    // all-core burst at the start of every browser dataset.
    let activeInitializers = workerCount;
    postInitializationTelemetry(request.id, workerCount, activeInitializers);
    const initializationTasks = Array.from(
      { length: workerCount },
      async (_, index) => {
        try {
          return await initializeSlot(index, module, request.references, request.settings);
        } finally {
          activeInitializers -= 1;
          postInitializationTelemetry(request.id, workerCount, activeInitializers);
        }
      },
    );
    slots.push(...await Promise.all(initializationTasks));

    let fatalError: Error | null = null;
    let wakeAvailability: (() => void) | null = null;
    let resolveFinished!: () => void;
    const finished = new Promise<void>((resolve) => {
      resolveFinished = resolve;
    });
    const pending = new Map<number, { result: ComputeBatch; slot: ComputeSlot }>();
    const maximumOutstanding = workerCount * 2;
    let nextCommit = 0;
    let parsed = 0;
    let committed = 0;
    let acknowledged = 0;
    let inputDone = false;
    let finalizing = false;
    let flushing = false;
    let lastProgress = 0.14;
    let bytesRead = 0;
    let totalBytes = sequenceSourceSize(request.query);
    let inputFraction = 0;
    let inputRecords = 0;
    let eligibleRecords = 0;
    let fastqFilterStats: FastqQualityFilterStats = emptyFastqQualityFilterStats(
      Boolean(request.fastqFilter?.enabled),
      Boolean(request.fastqFilter?.enabled && request.format === 2),
    );
    let telemetryRate: number | null = null;
    let telemetryPreviousTime = performance.now();
    let telemetryPreviousCommitted = 0;

    const assignmentPhase = (): AssignmentPhase => finalizing
      ? "finalizing"
      : inputDone ? "draining" : "streaming";

    const telemetrySnapshot = (): AssignmentTelemetry => ({
      phase: assignmentPhase(),
      activeWorkers: slots.reduce((count, slot) => count + (slot.busy ? 1 : 0), 0),
      totalWorkers: workerCount,
      queriesPerSecond: telemetryRate,
      recordsParsed: parsed,
      recordsCommitted: acknowledged,
      recordsOutstanding: Math.max(0, parsed - acknowledged),
    });

    const publishTelemetry = (updateRate: boolean) => {
      const now = performance.now();
      const elapsed = now - telemetryPreviousTime;
      if (updateRate && elapsed > 0) {
        telemetryRate = updateQueriesPerSecondEma(
          telemetryRate,
          acknowledged - telemetryPreviousCommitted,
          elapsed,
        );
        telemetryPreviousTime = now;
        telemetryPreviousCommitted = acknowledged;
      }
      self.postMessage({ id: request.id, type: "telemetry", telemetry: telemetrySnapshot() });
    };

    publishTelemetry(false);
    telemetryTimer = self.setInterval(
      () => publishTelemetry(true),
      ASSIGNMENT_TELEMETRY_INTERVAL_MS,
    );

    const fail = (error: Error) => {
      if (fatalError) return;
      fatalError = error;
      wakeAvailability?.();
      wakeAvailability = null;
      resolveFinished();
    };
    slots.forEach((slot) => installComputeHandler(slot, fail));

    const report = (stage: string) => {
      const completion = assignmentCompletion({ acknowledged, parsed, inputDone, inputRecords, eligibleRecords, countHint: request.countHint, subsampleSize: request.subsample?.size, inputFraction });
      lastProgress = Math.max(lastProgress, 0.14 + completion * 0.82);
      postProgress(request.id, stage, Math.min(0.96, lastProgress));
    };

    const release = (slot: ComputeSlot) => {
      slot.busy = false;
      wakeAvailability?.();
      wakeAvailability = null;
    };

    const maybeFinish = () => {
      if (inputDone && committed === parsed && slots.every((slot) => !slot.busy) && !pending.size) {
        resolveFinished();
      }
    };

    const flush = async () => {
      if (flushing || fatalError) return;
      flushing = true;
      try {
        while (pending.has(nextCommit)) {
          const { result, slot } = pending.get(nextCommit)!;
          pending.delete(nextCommit);
          committed += result.count;
          const acknowledgement = waitForAcknowledgement(request.id, nextCommit);
          self.postMessage({
            id: request.id,
            type: "batch",
            batch: nextCommit,
            header: result.header,
            body: result.body,
            count: result.count,
            processed: committed,
            total: inputDone ? parsed : null,
            milliseconds: result.milliseconds,
            doubleDHeader: result.doubleDHeader,
            doubleDBody: result.doubleDBody,
            doubleDCount: result.doubleDCount ?? 0,
          }, [result.body, ...(result.doubleDBody ? [result.doubleDBody] : [])]);
          await acknowledgement;
          acknowledged += result.count;
          nextCommit += 1;
          wakeAvailability?.();
          wakeAvailability = null;
          report(inputDone
            ? `Annotated ${committed.toLocaleString()} of ${parsed.toLocaleString()} sequences`
            : `Streamed + annotated ${committed.toLocaleString()} sequences`);
        }
      } catch (error) {
        fail(error instanceof Error ? error : new Error(String(error)));
      } finally {
        flushing = false;
        maybeFinish();
      }
    };

    const acquire = async (batchIndex: number): Promise<ComputeSlot> => {
      while (true) {
        if (fatalError) throw fatalError;
        const available = slots.find((slot) => !slot.busy);
        if (available && batchIndex - nextCommit < maximumOutstanding) {
          available.busy = true;
          return available;
        }
        await new Promise<void>((resolve) => { wakeAvailability = resolve; });
      }
    };

    postProgress(request.id, `Ready · ${workerCount} parallel worker${workerCount === 1 ? "" : "s"}`, 0.14);
    for await (const batch of streamSequenceBatches({
      source: request.query,
      format: request.format,
      batchSize: batchSize(request.subsample
        ? Math.min(request.countHint ?? request.subsample.size, request.subsample.size)
        : request.countHint),
      subsample: request.subsample,
      fastqFilter: request.fastqFilter,
      onProgress: (state) => {
        inputFraction = state.totalCharacters ? Math.min(1, state.charactersRead / state.totalCharacters) : 0;
        bytesRead = state.bytesRead;
        totalBytes = state.totalBytes;
        inputRecords = state.recordsRead;
        eligibleRecords = state.recordsEligible;
        fastqFilterStats = state.fastqFilter;
        if (state.recordsRead && state.recordsRead % 1000 === 0) {
          if (state.fastqFilter.applicable) {
            const rejected = state.fastqFilter.recordsRejectedExpectedErrors + state.fastqFilter.recordsRejectedMinimumLength;
            report(`FASTQ quality filter · ${state.recordsRead.toLocaleString()} scanned · ${state.recordsEligible.toLocaleString()} retained · ${rejected.toLocaleString()} rejected${request.subsample ? ` · sampling ${state.recordsSelected.toLocaleString()}` : ""}`);
          } else {
            report(request.subsample
              ? `Random subsample · ${state.recordsRead.toLocaleString()} scanned · retaining ${state.recordsSelected.toLocaleString()}`
              : `Streaming input · ${state.recordsRead.toLocaleString()} sequences parsed`);
          }
        }
      },
    })) {
      if (fatalError) throw fatalError;
      parsed += batch.count;
      const slot = await acquire(batch.index);
      void annotate(slot, batch, request).then((result) => {
        pending.set(batch.index, { result, slot });
        release(slot);
        void flush();
      }).catch(fail);
    }
    inputDone = true;
    publishTelemetry(false);
    bytesRead = totalBytes;
    if (fastqFilterStats.applicable) {
      const rejected = fastqFilterStats.recordsRejectedExpectedErrors + fastqFilterStats.recordsRejectedMinimumLength;
      report(`FASTQ filter complete · ${inputRecords.toLocaleString()} scanned · ${fastqFilterStats.recordsRetained.toLocaleString()} retained · ${rejected.toLocaleString()} rejected · finishing ${parsed.toLocaleString()} sequences`);
    } else if (request.fastqFilter?.enabled) {
      report(`FASTQ filter bypassed for non-FASTQ input · finishing ${parsed.toLocaleString()} sequences`);
    } else {
      report(`Input complete · finishing ${parsed.toLocaleString()} sequences`);
    }
    maybeFinish();
    await finished;
    if (fatalError) throw fatalError;
    finalizing = true;
    publishTelemetry(false);
    postProgress(request.id, "Finalizing the local AIRR index", 0.98);
    self.postMessage({
      id: request.id,
      type: "result",
      count: committed,
      total: parsed,
      inputRecords,
      workers: workerCount,
      fastqFilter: fastqFilterStats,
      // Every worker loads the same germlines and options: the first one's warnings are all of them.
      warnings: slots[0]?.warnings ?? [],
    });
  } catch (error) {
    self.postMessage({
      id: request.id,
      type: "error",
      message: error instanceof Error ? error.message : String(error),
    });
  } finally {
    if (telemetryTimer !== undefined) self.clearInterval(telemetryTimer);
    slots.forEach((slot) => slot.worker.terminate());
  }
}

self.addEventListener("message", (event: MessageEvent<StartRequest | AckRequest>) => {
  const message = event.data;
  if (message.type === "ack") {
    const key = `${message.id}:${message.batch}`;
    acknowledgements.get(key)?.();
    acknowledgements.delete(key);
    return;
  }
  void handleRequest(message);
});

export {};
