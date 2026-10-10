// standalone.js: the script of index.html, the page that runs opened from
// disk or from any static host. wasm/build.sh replaces each `//@ name` line
// with that file, so the page fetches nothing: a file: page may not fetch a
// module or start a worker from a file, but it may compile bytes it holds
// and start a worker from a blob.
'use strict';

/** Both WebAssembly modules, base64. */
const VDJ_WASM = {
//@ vdj_assign_simd.wasm
//@ vdj_assign.wasm
};

/** vdj-assign.js, as a function so a worker can run it too. */
function vdjAssignApi() {
//@ vdj-assign.js
  return { simdSupported, compile, load, VdjAssign, assignInWorkers };
}

/** worker.js, given the API. */
function vdjAssignWorker({ load }) {
//@ worker.js
}

/** A blob that runs the API and worker.js: every worker starts from it. */
const VDJ_WORKER = URL.createObjectURL(
  new Blob([`(${vdjAssignWorker})((${vdjAssignApi})());\n`], { type: 'text/javascript' }),
);

(async () => {
  const api = vdjAssignApi();
  const { load } = api;
  /** The module this engine runs, compiled from the embedded bytes. */
  const compile = ({ simd = api.simdSupported() } = {}) => {
    const text = atob(simd ? VDJ_WASM.simd : VDJ_WASM.baseline);
    const bytes = new Uint8Array(text.length);
    for (let i = 0; i < text.length; i++) bytes[i] = text.charCodeAt(i);
    return WebAssembly.compile(bytes);
  };
  const assignInWorkers = (options) =>
    api.assignInWorkers({ ...options, createWorker: () => new Worker(VDJ_WORKER) });
//@ app.js
})();
