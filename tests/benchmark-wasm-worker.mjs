import fs from "node:fs";
import { parentPort, workerData } from "node:worker_threads";

import { AlignerEngine } from "../src/aligner-engine.ts";

let engine;

async function initialize(references) {
  const wasmBytes = fs.readFileSync(
    workerData?.wasmPath || new URL("../public/swig.wasm", import.meta.url),
  );
  const module = await WebAssembly.compile(wasmBytes);
  const settings = {
    callingProfile: workerData?.profile || "r_optimized",
    minimumIdentity: 0.6,
    strand: 0,
  };
  const created = await AlignerEngine.create(module, settings, references);
  engine = created.engine;
  return created.records;
}

function annotate(message) {
  const batch = engine.annotate(new Uint8Array(message.query), message.count, false);
  const output = new Uint8Array(batch.header.length + 1 + batch.body.byteLength);
  output.set(new TextEncoder().encode(`${batch.header}\n`));
  output.set(batch.body, batch.header.length + 1);
  parentPort.postMessage({ type: "batch", id: message.id, count: batch.count, output: output.buffer }, [output.buffer]);
}

parentPort.on("message", (message) => {
  Promise.resolve(message.type === "initialize" ? initialize(message.references) : annotate(message)).then((genes) => {
    if (message.type === "initialize") parentPort.postMessage({ type: "ready", genes });
  }).catch((error) => parentPort.postMessage({ type: "error", id: message.id, message: error.message }));
});
