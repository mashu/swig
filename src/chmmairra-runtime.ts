import { chmmJunctionFlank, type ChmmOptions } from "./post-analysis-core";
import type { AirrOutputWritable, AirrResultStore } from "./result-store";
import { tableHeader, tableRow, type TableExportFormat } from "./export-formats";

export type ChmmSegment = "V" | "J";

interface WorkerResult {
  ordinal: number;
  probability: number;
  dfr: number;
  startingReference: string;
  recombinations: Array<{ position: number; left: string; right: string }>;
  status: "evaluated" | "low_dfr" | "missing_alignment" | "missing_reference" | "error";
  error?: string;
}

interface ChmmClient {
  worker: Worker;
  request: <T>(message: Record<string, unknown>) => Promise<T>;
  terminate: (reason?: Error) => void;
}

export interface ChmmRunOptions extends ChmmOptions {
  segment: ChmmSegment;
  minDfr: number;
  threshold: number;
  workers: number;
}

export interface ChmmFlaggedRecord {
  ordinal: number;
  probability: number;
  dfr: number;
  startingReference: string;
  recombinations: Array<{ position: number; left: string; right: string }>;
}

export interface ChmmDashboard {
  segment: ChmmSegment;
  threshold: number;
  inputRecords: number;
  upstreamExcluded: number;
  evaluated: number;
  flagged: number;
  lowDfr: number;
  missingAlignment: number;
  missingReference: number;
  errors: number;
  histogram: Array<{ label: string; count: number }>;
  top: ChmmFlaggedRecord[];
  probabilities: Float32Array;
  dfr: Uint16Array;
}

export interface ChmmDetail {
  ordinal: number;
  sequenceId: string;
  segment: ChmmSegment;
  call: string;
  probability: number;
  dfr: number;
  startingReference: string;
  recombinations: Array<{ position: number; left: string; right: string }>;
  threadedObservation: string;
  parents: Array<{ name: string; sequence: string }>;
}

function modelOptions(options: ChmmRunOptions): ChmmOptions {
  return {
    method: options.method,
    priorProbability: options.priorProbability,
    baseMutationProbability: options.baseMutationProbability,
    mutationRates: options.mutationRates,
    mutationSwitchProbability: options.mutationSwitchProbability,
    detailed: options.detailed,
    tracePath: false,
  };
}

function makeClient(): ChmmClient {
  const worker = new Worker(new URL("./chmmairra-worker.ts", import.meta.url), { type: "module" });
  const pending = new Map<number, { resolve: (value: unknown) => void; reject: (error: Error) => void }>();
  let nextId = 1;
  worker.onmessage = (event: MessageEvent<{ id: number; result?: unknown; error?: string }>) => {
    const request = pending.get(event.data.id);
    if (!request) return;
    pending.delete(event.data.id);
    if (event.data.error) request.reject(new Error(event.data.error));
    else request.resolve(event.data.result);
  };
  worker.onerror = (event) => {
    const error = new Error(event.message || "A CHMMAIRRa worker stopped unexpectedly.");
    pending.forEach((request) => request.reject(error));
    pending.clear();
  };
  return {
    worker,
    request: <T>(message: Record<string, unknown>) => {
      const id = nextId++;
      return new Promise<T>((resolve, reject) => {
        pending.set(id, { resolve: resolve as (value: unknown) => void, reject });
        worker.postMessage({ id, ...message });
      });
    },
    terminate: (reason = new Error("The CHMMAIRRa worker was closed.")) => {
      worker.terminate();
      pending.forEach((request) => request.reject(reason));
      pending.clear();
    },
  };
}

export async function runChmmairra(
  store: AirrResultStore,
  msa: string,
  options: ChmmRunOptions,
  workingMask?: Uint8Array,
  onProgress?: (processed: number, total: number) => void,
  signal?: AbortSignal,
): Promise<ChmmDashboard> {
  let inputRecords = store.count;
  if (workingMask) {
    if (workingMask.length !== store.count) throw new Error("The CHMMAIRRa working set does not match the AIRR record count.");
    inputRecords = 0;
    for (const value of workingMask) inputRecords += value ? 1 : 0;
  }
  const workerCount = Math.max(1, Math.min(options.workers, 16, inputRecords || 1));
  const clients = Array.from({ length: workerCount }, makeClient);
  const workerOptions = modelOptions(options);
  try {
    await Promise.all(clients.map((client) => client.request({ type: "init", msa, options: workerOptions, minDfr: options.minDfr })));
    const probabilities = new Float32Array(store.count);
    probabilities.fill(Number.NaN);
    const dfr = new Uint16Array(store.count);
    const histogram = new Uint32Array(10);
    const counters = { evaluated: 0, flagged: 0, lowDfr: 0, missingAlignment: 0, missingReference: 0, errors: 0 };
    const top: ChmmFlaggedRecord[] = [];
    const smaller = (left: ChmmFlaggedRecord, right: ChmmFlaggedRecord) => left.probability < right.probability || (left.probability === right.probability && left.dfr < right.dfr);
    const retainTop = (value: ChmmFlaggedRecord) => {
      if (top.length < 500) {
        top.push(value);
        let index = top.length - 1;
        while (index > 0) {
          const parent = Math.floor((index - 1) / 2);
          if (!smaller(top[index], top[parent])) break;
          [top[index], top[parent]] = [top[parent], top[index]];
          index = parent;
        }
        return;
      }
      if (!smaller(top[0], value)) return;
      top[0] = value;
      let index = 0;
      while (true) {
        const left = index * 2 + 1;
        const right = left + 1;
        let smallest = index;
        if (left < top.length && smaller(top[left], top[smallest])) smallest = left;
        if (right < top.length && smaller(top[right], top[smallest])) smallest = right;
        if (smallest === index) break;
        [top[index], top[smallest]] = [top[smallest], top[index]];
        index = smallest;
      }
    };
    const slots: Array<Promise<void>> = clients.map(() => Promise.resolve());
    let slot = 0;
    const consume = (results: WorkerResult[]) => {
      for (const result of results) {
        dfr[result.ordinal] = Math.min(65535, result.dfr);
        if (result.status === "evaluated") {
          probabilities[result.ordinal] = result.probability;
          counters.evaluated += 1;
          histogram[Math.min(9, Math.floor(result.probability * 10))] += 1;
          if (result.probability >= options.threshold) {
            counters.flagged += 1;
            retainTop({ ordinal: result.ordinal, probability: result.probability, dfr: result.dfr, startingReference: result.startingReference, recombinations: result.recombinations });
          }
        } else if (result.status === "low_dfr") counters.lowDfr += 1;
        else if (result.status === "missing_alignment") counters.missingAlignment += 1;
        else if (result.status === "missing_reference") counters.missingReference += 1;
        else counters.errors += 1;
      }
    };
    const prefix = options.segment.toLowerCase();
    await store.scanAirrRows(
      [`${prefix}_call`, `${prefix}_sequence_alignment`, `${prefix}_germline_alignment`, "sequence", "rev_comp", `${prefix}_sequence_start`, `${prefix}_sequence_end`],
      async (rows) => {
        if (!rows.length) return;
        const current = slot;
        slot = (slot + 1) % clients.length;
        await slots[current];
        const workerRows = rows.map((row) => ({
          ordinal: row.ordinal,
          call: row.values[`${prefix}_call`],
          sequenceAlignment: row.values[`${prefix}_sequence_alignment`],
          germlineAlignment: row.values[`${prefix}_germline_alignment`],
          segment: options.segment,
          flank: chmmJunctionFlank(options.segment, row.values.sequence, row.values.rev_comp, row.values[`${prefix}_sequence_start`], row.values[`${prefix}_sequence_end`]),
        }));
        slots[current] = clients[current].request<{ results: WorkerResult[] }>({ type: "batch", rows: workerRows }).then(({ results }) => consume(results));
      },
      { batchSize: 250, onProgress, signal, includeMask: workingMask },
    );
    await Promise.all(slots);
    return {
      segment: options.segment,
      threshold: options.threshold,
      inputRecords,
      upstreamExcluded: store.count - inputRecords,
      ...counters,
      histogram: Array.from(histogram, (count, index) => ({ label: `${(index / 10).toFixed(1)}–${((index + 1) / 10).toFixed(1)}`, count })),
      top: top.sort((a, b) => b.probability - a.probability || b.dfr - a.dfr || a.ordinal - b.ordinal),
      probabilities,
      dfr,
    };
  } finally {
    const reason = signal?.aborted ? new DOMException("CHMMAIRRa was cancelled.", "AbortError") : undefined;
    clients.forEach((client) => client.terminate(reason));
  }
}

export async function runChmmairraDetail(
  store: AirrResultStore,
  msa: string,
  options: ChmmRunOptions,
  ordinal: number,
  signal?: AbortSignal,
): Promise<ChmmDetail> {
  if (signal?.aborted) throw new DOMException("CHMMAIRRa detail reconstruction was cancelled.", "AbortError");
  const [detail] = await store.detailMany([ordinal]);
  if (!detail) throw new Error("That CHMMAIRRa record is no longer present in the local result index.");
  const prefix = options.segment.toLowerCase();
  const call = detail.values[`${prefix}_call`] ?? "";
  const client = makeClient();
  const abort = () => client.terminate(new DOMException("CHMMAIRRa detail reconstruction was cancelled.", "AbortError"));
  signal?.addEventListener("abort", abort, { once: true });
  try {
    await client.request({ type: "init", msa, options: modelOptions(options), minDfr: options.minDfr });
    const result = await client.request<Omit<ChmmDetail, "sequenceId" | "segment" | "call">>({
      type: "detail",
      row: {
        ordinal,
        call,
        sequenceAlignment: detail.values[`${prefix}_sequence_alignment`] ?? "",
        germlineAlignment: detail.values[`${prefix}_germline_alignment`] ?? "",
        segment: options.segment,
        flank: chmmJunctionFlank(options.segment, detail.values.sequence, detail.values.rev_comp, detail.values[`${prefix}_sequence_start`], detail.values[`${prefix}_sequence_end`]),
      },
    });
    if (signal?.aborted) throw new DOMException("CHMMAIRRa detail reconstruction was cancelled.", "AbortError");
    return {
      ...result,
      sequenceId: detail.values.sequence_id || detail.record.sequenceId,
      segment: options.segment,
      call,
    };
  } finally {
    signal?.removeEventListener("abort", abort);
    client.terminate();
  }
}

export async function writeChmmairraTsv(
  store: AirrResultStore,
  dashboard: ChmmDashboard,
  writable: Pick<AirrOutputWritable, "write">,
  signal?: AbortSignal,
): Promise<void> {
  return writeChmmairra(store, dashboard, "tsv", writable, signal);
}

export async function writeChmmairra(
  store: AirrResultStore,
  dashboard: ChmmDashboard,
  format: TableExportFormat,
  writable: Pick<AirrOutputWritable, "write">,
  signal?: AbortSignal,
): Promise<void> {
  const segment = dashboard.segment.toLowerCase();
  const fields = ["sequence_id", `${segment}_chimera_probability`, `${segment}_chimeric`, `${segment}_distance_from_reference`];
  const header = tableHeader(fields, format);
  if (header) await writable.write(header);
  await store.scanAirrRows(["sequence_id"], async (rows) => {
    let body = "";
    for (const row of rows) {
      const probability = dashboard.probabilities[row.ordinal];
      body += tableRow(fields, {
        sequence_id: row.values.sequence_id,
        [`${segment}_chimera_probability`]: Number.isNaN(probability) ? "" : probability.toFixed(8),
        [`${segment}_chimeric`]: Number.isNaN(probability) ? "" : probability >= dashboard.threshold ? "T" : "F",
        [`${segment}_distance_from_reference`]: dashboard.dfr[row.ordinal],
      }, format);
    }
    await writable.write(body);
  }, { signal });
}
