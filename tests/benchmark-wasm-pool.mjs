import assert from "node:assert/strict";
import crypto from "node:crypto";
import fs from "node:fs";
import os from "node:os";
import { performance } from "node:perf_hooks";
import { Worker } from "node:worker_threads";
import zlib from "node:zlib";

const TOTAL = Number(process.env.SWIG_BENCH_TOTAL || 50_000);
const BATCH_SIZE = Number(process.env.SWIG_BENCH_BATCH || 1000);
const requestedWorkers = Number(process.env.SWIG_BENCH_WORKERS || 0);
const WORKERS = requestedWorkers || Math.max(1, Math.min(8, os.availableParallelism() - 1));
const PROFILE = process.env.SWIG_BENCH_PROFILE || "r_optimized";
const encoder = new TextEncoder();

const pack = JSON.parse(zlib.gunzipSync(
  fs.readFileSync(new URL("../public/references/imgt-202632-7-swig-0.7.json.gz", import.meta.url)),
));
const human = pack.species.find((entry) => entry.name === "Homo sapiens");
assert.ok(human?.loci?.IGH, "Human IGH references are unavailable");
const locus = human.loci.IGH;

function fasta(records) {
  return records.map(([name, sequence, metadata]) => `>${name}${metadata ? ` SWIGMETA=${metadata.join(",")}` : ""}\n${sequence}\n`).join("");
}

const references = {
  V: fasta(locus.V),
  D: fasta(locus.D || []),
  J: fasta(locus.J),
  C: "",
};
const templateV = locus.V.find((allele) => allele[2]?.slice(2, 12).every((value) => value >= 0)) ?? locus.V[0];
const templateJ = locus.J.find((allele) => allele[2]?.[0] >= 0 && allele[2]?.[1] >= 0) ?? locus.J[0];
const template = `${templateV[1]}AACCGG${locus.D?.[0]?.[1] || ""}TTG${templateJ[1]}`;

class Client {
  constructor(index) {
    this.index = index;
    this.worker = new Worker(new URL("./benchmark-wasm-worker.mjs", import.meta.url), {
      workerData: {
        profile: PROFILE,
        wasmPath: process.env.SWIG_BENCH_WASM || "",
      },
    });
    this.pending = new Map();
    this.worker.on("message", (message) => {
      if (message.type === "ready") {
        this.readyResolve(message.genes);
        return;
      }
      const pending = this.pending.get(message.id);
      if (!pending) return;
      this.pending.delete(message.id);
      if (message.type === "error") pending.reject(new Error(message.message));
      else pending.resolve(message);
    });
    this.worker.on("error", (error) => {
      this.readyReject?.(error);
      for (const pending of this.pending.values()) pending.reject(error);
    });
    this.ready = new Promise((resolve, reject) => {
      this.readyResolve = resolve;
      this.readyReject = reject;
    });
    this.worker.postMessage({ type: "initialize", references });
  }

  annotate(id, query, count) {
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.worker.postMessage({ type: "annotate", id, query: query.buffer, count }, [query.buffer]);
    });
  }

  terminate() {
    return this.worker.terminate();
  }
}

const clients = Array.from({ length: WORKERS }, (_, index) => new Client(index));
const genes = await Promise.all(clients.map((client) => client.ready));
assert.ok(genes.every((count) => count === genes[0]));

let nextRecord = 0;
let completed = 0;
let outputBytes = 0;
const outputHash = crypto.createHash("sha256");
const baselineRss = process.resourceUsage().maxRSS * 1024;

// Every read differs (3% substitutions, seeded by its index), so neither
// worker's repeated-sequence cache can answer it.
function variant(index) {
  let state = (index * 2654435761 + 1) >>> 0;
  const random = () => (state = (Math.imul(state, 1664525) + 1013904223) >>> 0) / 2 ** 32;
  const bases = [...template];
  for (let position = 0; position < bases.length; position += 1) {
    if (random() < 0.03) bases[position] = "ACGT"[Math.floor(random() * 4)];
  }
  return bases.join("");
}

function makeBatch(start, count) {
  let text = "";
  for (let offset = 0; offset < count; offset += 1) {
    text += `>stress_${start + offset}\n${variant(start + offset)}\n`;
  }
  return encoder.encode(text);
}

const started = performance.now();
await Promise.all(clients.map(async (client) => {
  while (true) {
    const start = nextRecord;
    if (start >= TOTAL) return;
    const count = Math.min(BATCH_SIZE, TOTAL - start);
    nextRecord += count;
    const result = await client.annotate(start / BATCH_SIZE, makeBatch(start, count), count);
    assert.equal(result.count, count);
    completed += result.count;
    outputBytes += result.output.byteLength;
    outputHash.update(new Uint8Array(result.output));
  }
}));
const seconds = (performance.now() - started) / 1000;
const peakRss = process.resourceUsage().maxRSS * 1024;
await Promise.all(clients.map((client) => client.terminate()));

assert.equal(completed, TOTAL);
const summary = {
  records: completed,
  workers: WORKERS,
  profile: PROFILE,
  germlineAllelesPerWorker: genes[0],
  batchSize: BATCH_SIZE,
  seconds: Number(seconds.toFixed(3)),
  readsPerSecond: Math.round(completed / seconds),
  airrOutputMiB: Number((outputBytes / 1024 / 1024).toFixed(1)),
  outputSha256: outputHash.digest("hex"),
  baselineRssMiB: Number((baselineRss / 1024 / 1024).toFixed(1)),
  peakRssMiB: Number((peakRss / 1024 / 1024).toFixed(1)),
  boundedInFlightBatches: WORKERS,
};
process.stdout.write(`${JSON.stringify(summary, null, 2)}\n`);
