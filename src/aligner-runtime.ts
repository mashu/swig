import type { AlignerSettings } from "./aligner-engine";
import type { CompiledReferences } from "./reference-pack";
import type { FastqQualityFilterOptions, FastqQualityFilterStats, SequenceSource } from "./sequence-stream";
import type { AssignmentTelemetry } from "./assignment-telemetry";

export interface ResultBatch {
  header: string;
  body: Uint8Array;
  count: number;
  processed: number;
  total: number | null;
  doubleDHeader?: string;
  doubleDBody?: Uint8Array;
  doubleDCount?: number;
}

export type { CallingProfile, DoubleDScreenMode, DoubleDScreenOptions } from "./aligner-engine";

export interface RunOptions {
  query: SequenceSource;
  format: 1 | 2 | 3;
  references: CompiledReferences;
  settings: AlignerSettings;
  workers: number;
  countHint?: number | null;
  subsample?: { size: number; seed: number };
  fastqFilter?: FastqQualityFilterOptions;
  onProgress?: (stage: string, value: number) => void;
  onTelemetry?: (telemetry: AssignmentTelemetry) => void;
  onBatch?: (batch: ResultBatch) => void | Promise<void>;
  signal?: AbortSignal;
}

export interface RunResult {
  count: number;
  total: number;
  inputRecords: number;
  workers: number;
  fastqFilter: FastqQualityFilterStats;
  /** The aligner's database and run warnings. */
  warnings: string[];
}

let requestId = 0;

interface AnalysisLockManager {
  request: <T>(
    name: string,
    options: { mode: "exclusive"; signal: AbortSignal },
    callback: (lock: unknown) => Promise<T>,
  ) => Promise<T>;
}

/**
 * Chrome's Energy Saver freeze eligibility excludes a browsing context group
 * that owns a Web Lock. The lock also prevents two tabs on the same origin
 * from concurrently saturating the machine and competing for result storage.
 */
export async function withAnalysisWebLock<T>(
  signal: AbortSignal,
  onState: (state: "unsupported" | "waiting" | "held") => void,
  action: () => Promise<T>,
): Promise<T> {
  const locks = (navigator as Navigator & { locks?: AnalysisLockManager }).locks;
  if (!locks) {
    onState("unsupported");
    return action();
  }
  onState("waiting");
  return locks.request("swig-active-analysis", { mode: "exclusive", signal }, async () => {
    onState("held");
    return action();
  });
}

export function runAligner(options: RunOptions): Promise<RunResult> {
  const id = ++requestId;
  const worker = new Worker(new URL("./aligner-worker.ts", import.meta.url), { type: "module" });
  return new Promise((resolve, reject) => {
    let finished = false;
    const fail = (error: Error) => {
      if (finished) return;
      finished = true;
      options.signal?.removeEventListener("abort", abort);
      worker.terminate();
      reject(error);
    };
    const abort = () => fail(new DOMException("Analysis cancelled.", "AbortError"));
    if (options.signal?.aborted) {
      abort();
      return;
    }
    options.signal?.addEventListener("abort", abort, { once: true });

    worker.onmessage = (event) => {
      const message = event.data as {
        id: number;
        type: "progress" | "telemetry" | "batch" | "result" | "error";
        stage?: string;
        value?: number;
        batch?: number;
        header?: string;
        body?: ArrayBuffer;
        count?: number;
        processed?: number;
        total?: number | null;
        workers?: number;
        inputRecords?: number;
        fastqFilter?: FastqQualityFilterStats;
        message?: string;
        doubleDHeader?: string;
        doubleDBody?: ArrayBuffer;
        doubleDCount?: number;
        telemetry?: AssignmentTelemetry;
        warnings?: string[];
      };
      if (message.id !== id || finished) return;
      if (message.type === "progress") {
        options.onProgress?.(message.stage ?? "Working", message.value ?? 0);
        return;
      }
      if (message.type === "telemetry") {
        if (message.telemetry) options.onTelemetry?.(message.telemetry);
        return;
      }
      if (message.type === "batch") {
        Promise.resolve(options.onBatch?.({
          header: message.header ?? "",
          body: new Uint8Array(message.body ?? new ArrayBuffer(0)),
          count: message.count ?? 0,
          processed: message.processed ?? 0,
          total: message.total ?? null,
          doubleDHeader: message.doubleDHeader,
          doubleDBody: message.doubleDBody ? new Uint8Array(message.doubleDBody) : undefined,
          doubleDCount: message.doubleDCount ?? 0,
        })).then(() => {
          worker.postMessage({ type: "ack", id, batch: message.batch ?? 0 });
        }).catch((error) => fail(error instanceof Error ? error : new Error(String(error))));
        return;
      }
      if (message.type === "error") {
        fail(new Error(message.message ?? "The Swig aligner failed."));
        return;
      }
      finished = true;
      options.signal?.removeEventListener("abort", abort);
      worker.terminate();
      resolve({
        count: message.count ?? 0,
        total: message.total ?? 0,
        inputRecords: message.inputRecords ?? message.total ?? 0,
        workers: message.workers ?? 1,
        fastqFilter: message.fastqFilter ?? {
          enabled: false,
          applicable: false,
          recordsEvaluated: 0,
          recordsRetained: 0,
          recordsPassedThrough: 0,
          recordsRejectedExpectedErrors: 0,
          recordsRejectedMinimumLength: 0,
          recordsTrimmed: 0,
          basesTrimmed: 0,
        },
        warnings: message.warnings ?? [],
      });
    };
    worker.onerror = (event) => fail(new Error(event.message || "The aligner worker stopped unexpectedly."));
    worker.postMessage({
      type: "start",
      id,
      query: options.query,
      format: options.format,
      references: {
        V: options.references.V,
        D: options.references.D,
        J: options.references.J,
        C: options.references.C,
      },
      settings: options.settings,
      workers: options.workers,
      countHint: options.countHint ?? null,
      subsample: options.subsample,
      fastqFilter: options.fastqFilter,
    });
  });
}
