import { chmod, copyFile, mkdir, rm } from "node:fs/promises";
import { build } from "rolldown";

await rm(new URL("../cli",import.meta.url),{recursive:true,force:true});
await build({
  input:{"swig-cli":new URL("../cli-src/swig-cli-node.mjs",import.meta.url).pathname,"swig-worker":new URL("../cli-src/swig-worker.mjs",import.meta.url).pathname,"post-analysis-worker":new URL("../cli-src/post-analysis-worker.mjs",import.meta.url).pathname,"chmmairra-worker":new URL("../cli-src/chmmairra-worker.mjs",import.meta.url).pathname},
  external:/^node:/,
  output:{
    dir:new URL("../cli",import.meta.url).pathname,
    format:"es",
    entryFileNames:"[name].mjs",
    chunkFileNames:"chunks/[name]-[hash].mjs",
  },
});
const assets=new URL("../cli/assets/",import.meta.url);
await mkdir(assets,{recursive:true});
await copyFile(new URL("../public/swig.wasm",import.meta.url),new URL("swig.wasm",assets));
await copyFile(new URL("../public/references/imgt-202632-7-swig-0.7.json.gz",import.meta.url),new URL("imgt-reference-pack.json.gz",assets));
await chmod(new URL("../cli/swig-cli.mjs",import.meta.url),0o755);
