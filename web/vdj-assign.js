// vdj-assign.js: the JavaScript API of the vdj_assign WebAssembly module.
//
// No dependencies and no build step. Runs in browsers, Web Workers and
// Node.js 18 or later. The module computes what `vdj-assign` computes on the
// command line, byte for byte: the same options by the same names, the same
// database preparation, the same AIRR rows.
//
//   import { load } from './vdj-assign.js';
//   const vdj = await load();            // the SIMD build where the engine has SIMD
//   const run = vdj.createRun({ v, d, j, options: { profile: 'r', 'tandem-d': true } });
//   const tsv = run.assign(reads);       // FASTA or FASTQ text in, AIRR TSV out
//   run.free();
//
// In a browser, assignInWorkers() spreads one run over Web Workers.
// In Node.js, pass the module's bytes: load({ bytes: fs.readFileSync(path) }).

const SIMD_FILE = 'vdj_assign_simd.wasm';
const BASELINE_FILE = 'vdj_assign.wasm';

// The smallest module using a SIMD instruction: it validates exactly where
// the engine runs WebAssembly SIMD.
const SIMD_PROBE = Uint8Array.of(
  0, 97, 115, 109, 1, 0, 0, 0, 1, 5, 1, 96, 0, 1, 123, 3, 2, 1, 0, 10, 10, 1, 8, 0, 65, 0, 253, 15,
  253, 98, 11,
);

/** Whether this engine runs WebAssembly SIMD (every current browser does). */
export function simdSupported() {
  return WebAssembly.validate(SIMD_PROBE);
}

/**
 * Compile the module: the SIMD build when the engine has SIMD, else the
 * baseline build, fetched next to this file unless `url` says otherwise. A
 * compiled module can be posted to workers, which then only instantiate it.
 */
export async function compile({ simd = simdSupported(), url } = {}) {
  url ??= new URL(simd ? SIMD_FILE : BASELINE_FILE, import.meta.url);
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${url}: ${response.status} ${response.statusText}`);
  if (response.headers.get('content-type') === 'application/wasm') {
    return WebAssembly.compileStreaming(response);
  }
  return WebAssembly.compile(await response.arrayBuffer());
}

/**
 * Instantiate the module from a compiled `module`, its `bytes`, or by
 * compiling it (see compile()).
 */
export async function load({ module, bytes, ...options } = {}) {
  module ??= bytes ? await WebAssembly.compile(bytes) : await compile(options);
  const instance = await WebAssembly.instantiate(module, {});
  return new VdjAssign(instance.exports);
}

/** One instance of the module. */
export class VdjAssign {
  #x;
  #encoder = new TextEncoder();
  #decoder = new TextDecoder();
  #broken = null;

  constructor(exports) {
    this.#x = exports;
    exports.vdj_init();
    const info = this.#call(() => {
      exports.vdj_version();
      return JSON.parse(this.#output());
    });
    /** The crate version. */
    this.version = info.version;
    /** The alignment kernel a run uses by default: simd128 or scalar. */
    this.backend = info.backend;
  }

  /**
   * Demonstration data: a synthetic IGH-like database and `reads` reads
   * simulated under IgSim preset `preset` (easy, mid, train or hard), half
   * of them reverse-complemented. Returns { 'V.fasta', 'D.fasta', 'J.fasta',
   * 'C.fasta', 'reads.fasta', 'truth.tsv' }. The synthetic database is a
   * test fixture, not a real locus: most of its V genes carry stop codons.
   */
  demo({ seed = 1, reads = 1000, preset = 'mid' } = {}) {
    return this.#call(() => {
      this.#check(this.#withText(preset, (ptr, len) => this.#x.vdj_demo(seed, reads, ptr, len)));
      return JSON.parse(this.#output());
    });
  }

  /**
   * The prepared (deduplicated, annotated) database of `germlines`, as
   * `vdj-assign prepare` writes it: { 'V.fasta', 'D.fasta', 'J.fasta',
   * 'C.fasta', 'annotation.tsv' } (D and C when given).
   */
  prepare(germlines) {
    return this.#call(() => {
      const builder = this.#builder(germlines);
      try {
        this.#check(this.#x.vdj_builder_prepare(builder));
        return JSON.parse(this.#output());
      } finally {
        this.#x.vdj_builder_free(builder);
      }
    });
  }

  /**
   * A run over germline FASTA texts { v, d, j, c } (d and c optional) with
   * `options` named as on the command line: { profile: 'igblast',
   * 'tandem-d': true, 'c-evalue': 0.01 }. Flags take `true`; `false`, null
   * and undefined leave an option out. Free the run with free().
   */
  createRun(germlines) {
    return this.#call(() => {
      const x = this.#x;
      const handle = x.vdj_run_new(this.#builder(germlines));
      if (!handle) throw new Error(this.#output());
      const { summary, warnings, backend, header } = JSON.parse(this.#output());
      const vdj = this;
      let live = handle;
      const alive = () => {
        if (!live) throw new Error('the run was freed');
        return live;
      };
      return {
        /** The database summary line. */
        summary,
        /** What the command line would warn about the database. */
        warnings,
        /** The alignment kernel. */
        backend,
        /** The AIRR TSV header, without a line break. */
        header,
        /**
         * A session over `reads` (FASTA or FASTQ, as text or bytes) that
         * assigns block `worker`, `worker + workers`, … of `block` reads:
         * one session per worker shares a run. next() returns { block,
         * rows } (rows: AIRR lines, no header) or null at the end.
         */
        session(reads, { worker = 0, workers = 1, block = 256 } = {}) {
          return vdj.#session(alive(), reads, worker, workers, block);
        },
        /** AIRR TSV of every read, header included. */
        assign(reads) {
          const session = this.session(reads);
          try {
            const parts = [header, '\n'];
            for (let next; (next = session.next()); ) parts.push(next.rows);
            return parts.join('');
          } finally {
            session.free();
          }
        },
        free() {
          if (live) x.vdj_run_free(live);
          live = 0;
        },
      };
    });
  }

  #session(run, reads, worker, workers, block) {
    return this.#call(() => {
      const x = this.#x;
      // The session takes the buffer, also when it fails.
      const [ptr, len] = this.#buffer(reads);
      const handle = x.vdj_session_new(run, ptr, len, worker, workers, block);
      if (!handle) throw new Error(this.#output());
      let live = handle;
      const alive = () => {
        if (!live) throw new Error('the session was freed');
        return live;
      };
      return {
        next: () =>
          this.#call(() => {
            const index = x.vdj_session_next(alive());
            if (index === -1) return null;
            if (index < 0) throw new Error(this.#output());
            return { block: index, rows: this.#output() };
          }),
        /** The fraction of the reads consumed so far. */
        get progress() {
          return x.vdj_session_progress(alive());
        },
        free() {
          if (live) x.vdj_session_free(live);
          live = 0;
        },
      };
    });
  }

  /** A builder holding `germlines` and their options. */
  #builder({ v, d, j, c, options = {} }) {
    const x = this.#x;
    const builder = x.vdj_builder_new();
    const set = (key, value) => {
      const status = this.#withText(key, (keyPtr, keyLen) =>
        value === true
          ? x.vdj_builder_set(builder, keyPtr, keyLen, 0, 0, 0)
          : this.#withText(String(value), (ptr, len) =>
              x.vdj_builder_set(builder, keyPtr, keyLen, ptr, len, 1),
            ),
      );
      this.#check(status);
    };
    try {
      for (const [key, text] of Object.entries({ v, d, j, c })) {
        if (text != null) set(key, text);
      }
      for (const [key, value] of Object.entries(options)) {
        if (value != null && value !== false) set(key, value);
      }
      return builder;
    } catch (error) {
      x.vdj_builder_free(builder);
      throw error;
    }
  }

  /** Run `f`; a trap leaves the instance unusable, with the panic message. */
  #call(f) {
    if (this.#broken) throw this.#broken;
    try {
      return f();
    } catch (error) {
      if (!(error instanceof WebAssembly.RuntimeError)) throw error;
      let message = error.message;
      try {
        this.#x.vdj_panic_message();
        message = this.#output() || message;
      } catch {
        // The instance cannot even report its panic.
      }
      this.#broken = new Error(`vdj_assign stopped: ${message}`);
      throw this.#broken;
    }
  }

  /** The module's output buffer as text. */
  #output() {
    const x = this.#x;
    const bytes = new Uint8Array(x.memory.buffer, x.vdj_output_ptr(), x.vdj_output_len());
    return this.#decoder.decode(bytes);
  }

  #check(status) {
    if (status !== 0) throw new Error(this.#output());
  }

  /** `text` (a string or bytes) copied into a module buffer: [ptr, len]. */
  #buffer(text) {
    const bytes = typeof text === 'string' ? this.#encoder.encode(text) : new Uint8Array(text);
    const ptr = this.#x.vdj_alloc(bytes.length);
    // After the allocation: it may have grown, and so replaced, the memory.
    new Uint8Array(this.#x.memory.buffer, ptr, bytes.length).set(bytes);
    return [ptr, bytes.length];
  }

  /** f(ptr, len) over a temporary copy of `text`. */
  #withText(text, f) {
    const [ptr, len] = this.#buffer(text);
    try {
      return f(ptr, len);
    } finally {
      this.#x.vdj_free(ptr, len);
    }
  }
}

/**
 * Assign `reads` with `germlines` (as for createRun) on `workers` Web Workers
 * (browsers only), each instantiating the module once and taking every
 * `workers`-th block of reads. Calls onInfo(run) once a worker has built
 * the run and onProgress(fraction) as blocks finish. Resolves to { tsv,
 * reads, called, info }: the AIRR TSV in input order, header included.
 */
export async function assignInWorkers({
  germlines,
  reads,
  workers = globalThis.navigator?.hardwareConcurrency || 1,
  block = 256,
  module,
  onInfo,
  onProgress,
}) {
  module ??= await compile();
  const bytes = typeof reads === 'string' ? new TextEncoder().encode(reads) : new Uint8Array(reads);
  const blocks = [];
  const progress = new Array(workers).fill(0);
  const threads = [];
  let info;
  const report = () => onProgress?.(progress.reduce((sum, p) => sum + p, 0) / workers);
  try {
    await Promise.all(
      Array.from(
        { length: workers },
        (_, worker) =>
          new Promise((resolve, reject) => {
            const thread = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
            threads.push(thread);
            thread.onmessage = ({ data }) => {
              if (data.type === 'ready') {
                if (!info) onInfo?.((info = data.info));
              } else if (data.type === 'rows') {
                blocks[data.block] = data.rows;
                progress[worker] = data.progress;
                report();
              } else if (data.type === 'done') {
                progress[worker] = 1;
                report();
                resolve();
              } else if (data.type === 'error') {
                reject(new Error(data.message));
              }
            };
            thread.onerror = (event) => reject(new Error(event.message || 'worker failed'));
            thread.postMessage({ module, germlines, reads: bytes, worker, workers, block });
          }),
      ),
    );
  } finally {
    for (const thread of threads) thread.terminate();
  }
  const rows = blocks.join('');
  const columns = info.header.split('\t');
  const [v, j] = [columns.indexOf('v_call'), columns.indexOf('j_call')];
  let count = 0;
  let called = 0;
  for (const line of rows.split('\n')) {
    if (!line) continue;
    count += 1;
    // A call holds an anchor, a V or a J, so its v_call or j_call is set.
    const fields = line.split('\t');
    if (fields[v] || fields[j]) called += 1;
  }
  return { tsv: `${info.header}\n${rows}`, reads: count, called, info };
}
