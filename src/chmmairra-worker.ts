/// <reference lib="webworker" />

import {
  chmmairraDistanceFromReference,
  extendChmmAlignment,
  prepareReferenceMsa,
  runChmm,
  threadSequenceToMsa,
  type ChmmOptions,
  type ReferenceMsa,
} from "./post-analysis-core";

interface ChmmRow {
  ordinal: number;
  call: string;
  sequenceAlignment: string;
  germlineAlignment: string;
  /** "V" or "J", with `flank` the read bases on the junction side of the local alignment. */
  segment?: "V" | "J";
  flank?: string;
}

function extended(row: ChmmRow): ChmmRow {
  if (!row.segment || !row.flank || !msa) return row;
  return { ...row, ...extendChmmAlignment(row.segment, row.sequenceAlignment, row.germlineAlignment, row.call, row.flank, msa) };
}

interface ChmmWorkerResult {
  ordinal: number;
  probability: number;
  dfr: number;
  startingReference: string;
  recombinations: Array<{ position: number; left: string; right: string }>;
  status: "evaluated" | "low_dfr" | "missing_alignment" | "missing_reference" | "error";
  error?: string;
}

type Request =
  | { id: number; type: "init"; msa: string; options: ChmmOptions; minDfr: number }
  | { id: number; type: "batch"; rows: ChmmRow[] }
  | { id: number; type: "detail"; row: ChmmRow };

const worker = self as unknown as DedicatedWorkerGlobalScope;
let msa: ReferenceMsa | null = null;
let options: ChmmOptions | null = null;
let minDfr = 1;
const cache = new Map<string, Omit<ChmmWorkerResult, "ordinal">>();

worker.onmessage = (event: MessageEvent<Request>) => {
  const request = event.data;
  try {
    if (request.type === "init") {
      msa = prepareReferenceMsa(request.msa);
      options = request.options;
      minDfr = request.minDfr;
      cache.clear();
      worker.postMessage({ id: request.id, result: { references: msa.names.length, length: msa.length } });
      return;
    }
    if (!msa || !options) throw new Error("CHMMAIRRa worker is not initialized.");
    if (request.type === "detail") {
      const row = extended(request.row);
      if (!row.call || !row.sequenceAlignment || !row.germlineAlignment) {
        throw new Error("This record does not contain the selected segment alignment.");
      }
      const threadedObservation = threadSequenceToMsa(row.sequenceAlignment, row.germlineAlignment, row.call, msa);
      const result = runChmm(msa, threadedObservation, row.sequenceAlignment, row.germlineAlignment, {
        ...options,
        detailed: true,
        tracePath: true,
      });
      if (!result.referencePath) throw new Error("The CHMMAIRRa Viterbi trace was not produced.");
      const used = [...new Set(result.referencePath)].sort((left, right) => {
        const leftFirst = result.referencePath!.indexOf(left);
        const rightFirst = result.referencePath!.indexOf(right);
        return leftFirst - rightFirst;
      });
      worker.postMessage({ id: request.id, result: {
        ordinal: row.ordinal,
        probability: result.probability,
        dfr: result.dfr,
        startingReference: result.startingReference,
        recombinations: result.recombinations,
        threadedObservation,
        parents: used.map((reference) => ({ name: msa!.names[reference], sequence: msa!.sequences[reference] })),
      } });
      return;
    }
    const results: ChmmWorkerResult[] = [];
    for (const input of request.rows) {
      if (!input.call || !input.sequenceAlignment || !input.germlineAlignment) {
        results.push({ ordinal: input.ordinal, probability: Number.NaN, dfr: 0, startingReference: "", recombinations: [], status: "missing_alignment" });
        continue;
      }
      const row = extended(input);
      const dfr = chmmairraDistanceFromReference(row.sequenceAlignment, row.germlineAlignment);
      if (dfr < minDfr) {
        results.push({ ordinal: row.ordinal, probability: Number.NaN, dfr, startingReference: "", recombinations: [], status: "low_dfr" });
        continue;
      }
      const key = `${row.call}\u0000${row.sequenceAlignment}\u0000${row.germlineAlignment}`;
      const cached = cache.get(key);
      if (cached) {
        results.push({ ordinal: row.ordinal, ...cached });
        continue;
      }
      try {
        const threaded = threadSequenceToMsa(row.sequenceAlignment, row.germlineAlignment, row.call, msa);
        const result = runChmm(msa, threaded, row.sequenceAlignment, row.germlineAlignment, options);
        const value: Omit<ChmmWorkerResult, "ordinal"> = { ...result, status: "evaluated" };
        cache.set(key, value);
        results.push({ ordinal: row.ordinal, ...value });
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        const status = message.includes("missing") ? "missing_reference" : "error";
        const value: Omit<ChmmWorkerResult, "ordinal"> = { probability: Number.NaN, dfr, startingReference: "", recombinations: [], status, error: message };
        cache.set(key, value);
        results.push({ ordinal: row.ordinal, ...value });
      }
    }
    worker.postMessage({ id: request.id, result: { results } });
  } catch (error) {
    worker.postMessage({ id: request.id, error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
