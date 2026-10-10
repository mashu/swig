#!/usr/bin/env bun

// Bun embeds these immutable assets in the standalone executable's virtual
// filesystem. End users receive one native binary and need neither Bun nor
// Node nor a separately installed WASM runtime. native/swig is the platform's
// native aligner (staged by the release workflow, or locally by
// scripts/stage-native-aligner.mjs); swig-cli extracts it to the user cache and
// runs it, and keeps the WebAssembly module as the fallback.
import referencePackPath from "../public/references/imgt-202632-7-swig-0.7.json.gz" with { type: "file" };
import wasmPath from "../public/swig.wasm" with { type: "file" };
import nativeAlignerPath from "../native/swig" with { type: "file" };
import { runCli } from "./swig-cli.mjs";

runCli({wasmPath,referencePackPath,nativeAlignerPath}).catch((error)=>{
  process.stderr.write(`swig-cli: ${error instanceof Error?error.message:String(error)}\n`);
  process.exitCode=1;
});
