#!/usr/bin/env bun

// Bun embeds these immutable assets in the standalone executable's virtual
// filesystem. End users receive one native binary and need neither Bun nor
// Node nor a separately installed WASM runtime.
import referencePackPath from "../public/references/imgt-202632-7-swig-0.7.json.gz" with { type: "file" };
import wasmPath from "../public/swig.wasm" with { type: "file" };
import { runCli } from "./swig-cli.mjs";

runCli({wasmPath,referencePackPath}).catch((error)=>{
  process.stderr.write(`swig-cli: ${error instanceof Error?error.message:String(error)}\n`);
  process.exitCode=1;
});
