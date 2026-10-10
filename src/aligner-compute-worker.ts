/// <reference lib="webworker" />

import { AlignerEngine, type AlignerSettings, type GermlineTexts } from "./aligner-engine";

interface InitializeRequest {
  type: "initialize";
  worker: number;
  module: WebAssembly.Module;
  references: GermlineTexts;
  settings: AlignerSettings;
}

interface AnnotateRequest {
  type: "annotate";
  batch: number;
  query: ArrayBuffer;
  count: number;
  doubleD: boolean;
}

let engine: AlignerEngine | null = null;

async function initialize(request: InitializeRequest) {
  const created = await AlignerEngine.create(request.module, request.settings, request.references);
  engine = created.engine;
  self.postMessage({ type: "ready", worker: request.worker, records: created.records });
}

function annotate(request: AnnotateRequest) {
  if (!engine) throw new Error("The aligner's compute worker is not initialized.");
  const started = performance.now();
  const result = engine.annotate(new Uint8Array(request.query), request.count, request.doubleD);
  const message: Record<string, unknown> = {
    type: "batch",
    batch: request.batch,
    header: result.header,
    body: result.body.buffer,
    count: result.count,
    milliseconds: performance.now() - started,
  };
  const transfers: ArrayBuffer[] = [result.body.buffer as ArrayBuffer];
  if (result.doubleDBody) {
    message.doubleDHeader = result.doubleDHeader;
    message.doubleDBody = result.doubleDBody.buffer;
    message.doubleDCount = result.doubleDCount;
    transfers.push(result.doubleDBody.buffer as ArrayBuffer);
  }
  self.postMessage(message, transfers);
}

self.addEventListener("message", (event: MessageEvent<InitializeRequest | AnnotateRequest>) => {
  const request = event.data;
  Promise.resolve().then(() => request.type === "initialize" ? initialize(request) : annotate(request)).catch((error) => {
    self.postMessage({
      type: "error",
      batch: request.type === "annotate" ? request.batch : -1,
      message: error instanceof Error ? error.message : String(error),
    });
  });
});

export {};
