// The native Swig aligner (engine/, the `swig` binary) as a child process,
// behind the same interface as the WebAssembly worker pool: start(),
// run({text,count}) → {header,body,count}, close().
//
// One child per input stream. The CLI writes its prepared FASTA/FASTQ/AIRR
// records to the child's stdin and reads the AIRR TSV from its stdout. The
// engine writes exactly one row per record, in input order, so rows are handed
// back to the queued batches by counting line breaks. It assigns ENGINE_BATCH
// records at a time and writes a batch's rows only once it has read that many
// (or the input ended), so the caller keeps more than ENGINE_BATCH records in
// flight past the batch it awaits (`window`), and calls finishInput() after its
// last batch, before it awaits the batches still queued.

import { spawn, spawnSync } from "node:child_process";
import { createHash, randomUUID } from "node:crypto";
import { chmodSync, existsSync, lstatSync, mkdirSync, mkdtempSync, readFileSync, renameSync, rmSync, writeFileSync } from "node:fs";
import { homedir, tmpdir } from "node:os";
import { join } from "node:path";
import { alignerOptions } from "../src/aligner-engine.ts";

/** `swig::run::BATCH` (engine/src/run.rs): records the engine reads before it writes rows. */
export const ENGINE_BATCH = 16_384;

/** The aligner a run uses: "auto" (native when available, else WebAssembly), "native" or "wasm". */
export const ALIGNER_MODES = ["auto", "native", "wasm"];

function executableName() {
  return process.platform === "win32" ? "swig.exe" : "swig";
}

/**
 * Where an embedded binary is extracted to, most preferred first. The last
 * resort is the shared temporary directory, so on POSIX it gets a directory of
 * its own per user (`private`: checked before anything in it is run).
 */
function cacheRoots() {
  const roots = [];
  if (process.env.SWIG_CLI_CACHE_DIR) roots.push({ path: process.env.SWIG_CLI_CACHE_DIR, private: false });
  if (process.platform === "win32") {
    if (process.env.LOCALAPPDATA) roots.push({ path: process.env.LOCALAPPDATA, private: false });
  } else if (process.platform === "darwin") {
    roots.push({ path: join(homedir(), "Library", "Caches"), private: false });
  } else {
    roots.push({ path: process.env.XDG_CACHE_HOME || join(homedir(), ".cache"), private: false });
  }
  if (typeof process.getuid === "function") roots.push({ path: join(tmpdir(), `swig-cli-${process.getuid()}`), private: true });
  else roots.push({ path: tmpdir(), private: false });
  return roots;
}

/** Create `path` for this user only, and refuse it if another user could change what is in it. */
function privateDirectory(path) {
  mkdirSync(path, { recursive: true, mode: 0o700 });
  const status = lstatSync(path);
  if (!status.isDirectory() || status.uid !== process.getuid() || (status.mode & 0o077) !== 0) {
    throw new Error(`${path} is not a directory private to this user`);
  }
}

function sha256(bytes) {
  return createHash("sha256").update(bytes).digest("hex");
}

/**
 * A binary embedded in a Bun executable cannot be executed in place: copy it
 * once per version and content hash to a cache directory (atomically, so
 * concurrent first runs agree) and reuse the copy after checking its hash.
 */
function extract(sourcePath, version) {
  const bytes = readFileSync(sourcePath);
  const digest = sha256(bytes);
  const failures = [];
  for (const root of cacheRoots()) {
    const directory = join(root.path, "swig-cli", `${version}-${digest.slice(0, 16)}`);
    const target = join(directory, executableName());
    try {
      if (root.private) privateDirectory(root.path);
      if (existsSync(target) && sha256(readFileSync(target)) === digest) return target;
      mkdirSync(directory, { recursive: true });
      const temporary = `${target}.${randomUUID()}.tmp`;
      writeFileSync(temporary, bytes, { mode: 0o755 });
      chmodSync(temporary, 0o755);
      try {
        renameSync(temporary, target);
      } catch (error) {
        // Another process may have won the race (Windows refuses to replace a running file).
        rmSync(temporary, { force: true });
        if (!(existsSync(target) && sha256(readFileSync(target)) === digest)) throw error;
      }
      return target;
    } catch (error) {
      failures.push(`${directory}: ${error instanceof Error ? error.message : String(error)}`);
    }
  }
  throw new Error(`the embedded native aligner could not be extracted (${failures.join("; ")})`);
}

/** Why `binary` cannot serve this CLI version, or null. */
function probe(binary, version) {
  const result = spawnSync(binary, ["--version"], { encoding: "utf8", timeout: 15_000, windowsHide: true });
  if (result.error) return `${binary}: ${result.error.message}`;
  if (result.status !== 0) {
    // The loader's reason, e.g. a missing shared library.
    const detail = (result.stderr ?? "").trim().split(/\r?\n/).at(-1);
    return `${binary} --version exited with ${result.status === null ? `signal ${result.signal}` : `status ${result.status}`}${detail ? `: ${detail}` : ""}`;
  }
  const reported = result.stdout.trim();
  return reported === `swig ${version}` ? null : `${binary} reports '${reported}', not 'swig ${version}'`;
}

/**
 * The native aligner to run, or null with the reasons appended to `notes`:
 * SWIG_CLI_ALIGNER_PATH, then the binary embedded in a standalone executable
 * (`assets.nativeAlignerPath`). It must report this CLI's version.
 */
export function resolveNativeAligner(assets, version, notes = []) {
  const candidates = [];
  if (process.env.SWIG_CLI_ALIGNER_PATH) candidates.push(() => process.env.SWIG_CLI_ALIGNER_PATH);
  if (assets.nativeAlignerPath) {
    // Bun serves embedded files from a virtual file system ($bunfs, or ~BUN on Windows).
    const embedded = /\$bunfs|~BUN/.test(assets.nativeAlignerPath);
    candidates.push(() => (embedded ? extract(assets.nativeAlignerPath, version) : assets.nativeAlignerPath));
  }
  for (const candidate of candidates) {
    let path;
    try {
      path = candidate();
    } catch (error) {
      notes.push(error instanceof Error ? error.message : String(error));
      continue;
    }
    if (!existsSync(path)) {
      notes.push(`${path}: not found`);
      continue;
    }
    const problem = probe(path, version);
    if (!problem) return path;
    notes.push(problem);
  }
  return null;
}

/** A refused batch: the caller may never await it (ready() throws first), so it is marked handled. */
function refuse(error) {
  const refused = Promise.reject(error);
  refused.catch(() => {});
  return refused;
}

export class NativePool {
  /**
   * @param binary the native `swig` executable
   * @param threads its `--threads`
   * @param init {references: {V, D, J, C}, settings}, as for the WebAssembly pool
   * @param batchRecords records per run() call, which sets the in-flight window
   */
  constructor(binary, threads, init, batchRecords) {
    this.binary = binary;
    this.threads = threads;
    this.init = init;
    this.batchRecords = batchRecords;
    this.jobs = [];
    this.header = null;
    this.failed = null;
    this.stderr = "";
    this.airrHeader = null;
    this.carry = null;
    this.ended = false;
    this.drain = null;
    this.records = 0;
    this.rows = 0;
  }

  /**
   * Batches the caller may hold before it awaits the oldest one. The caller
   * awaits the oldest only when a further batch is ready to queue, so every
   * batch held is full and (window - 1) × batchRecords > ENGINE_BATCH records
   * follow the oldest: enough to close the engine batch that holds its rows.
   */
  get window() {
    return Math.floor(ENGINE_BATCH / Math.max(1, this.batchRecords)) + 2;
  }

  async start() {
    this.directory = mkdtempSync(join(tmpdir(), "swig-cli-native-"));
    const args = [];
    for (const segment of ["V", "D", "J", "C"]) {
      const text = this.init.references[segment];
      if (!text) continue;
      const path = join(this.directory, `${segment}.fasta`);
      writeFileSync(path, text);
      args.push(`--${segment.toLowerCase()}`, path);
    }
    args.push("--threads", String(this.threads));
    for (const [name, value] of alignerOptions(this.init.settings)) {
      if (name.startsWith("double-d")) throw new Error("the native route does not carry the double-D screen's sidecar");
      args.push(`--${name}`, value);
    }
    args.push("-");
    this.child = spawn(this.binary, args, { stdio: ["pipe", "pipe", "pipe"], windowsHide: true });
    this.exited = new Promise((resolve) => this.child.on("close", (code, signal) => resolve({ code, signal })));
    this.child.on("error", (error) => this.fail(error));
    this.child.stdin.on("error", (error) => {
      if (error.code !== "EPIPE") this.fail(error);
    });
    this.child.stdout.on("data", (chunk) => this.receive(chunk));
    this.child.stderr.setEncoding("utf8");
    let readyResolve;
    let readyReject;
    const ready = new Promise((resolve, reject) => {
      readyResolve = resolve;
      readyReject = reject;
    });
    this.child.stderr.on("data", (text) => {
      this.stderr = (this.stderr + text).slice(-65_536);
      if (/database prepared/.test(this.stderr)) readyResolve();
    });
    this.exited.then(({ code, signal }) => {
      const message = this.exitMessage(code, signal);
      if (message) {
        const error = new Error(message);
        readyReject(error);
        this.fail(error);
      } else if (this.jobs.length) {
        this.fail(
          new Error(
            `The native aligner exited with ${this.jobs.length} batch(es) unanswered (${this.rows.toLocaleString()} rows for ${this.records.toLocaleString()} records).`,
          ),
        );
      }
      readyResolve();
    });
    await ready;
  }

  exitMessage(code, signal) {
    if (code === 0) return null;
    const lines = this.stderr.split(/\r?\n/).filter((line) => line && !/^swig: (database prepared|warning:)/.test(line));
    // An argument error is followed by the usage text, and a panic by a note:
    // report the engine's own "swig: " line, or the panic message.
    const error = lines.findLast((line) => line.startsWith("swig: "));
    const panic = lines.findIndex((line) => / panicked at /.test(line));
    const detail = error ? error.slice("swig: ".length) : panic >= 0 ? lines.slice(panic, panic + 2).join(" ") : (lines.at(-1) ?? "");
    return `The native Swig aligner ${signal ? `was stopped by ${signal}` : `exited with status ${code}`}${detail ? `: ${detail}` : "."}`;
  }

  fail(error) {
    if (this.failed) return;
    this.failed = error;
    for (const job of this.jobs.splice(0)) job.reject(error);
  }

  /** Rows from stdout, handed to the oldest unanswered batch by line count. */
  receive(chunk) {
    const buffer = this.carry ? Buffer.concat([this.carry, chunk]) : chunk;
    this.carry = null;
    let start = 0;
    if (this.header === null) {
      const newline = buffer.indexOf(10);
      if (newline < 0) {
        this.carry = buffer;
        return;
      }
      this.header = buffer.toString("utf8", 0, newline).replace(/\r$/, "");
      start = newline + 1;
    }
    while (start < buffer.length) {
      const job = this.jobs[0];
      if (!job) {
        this.fail(new Error("The native aligner returned more rows than records."));
        return;
      }
      let end = -1;
      let position = start;
      while (job.remaining > 0) {
        const newline = buffer.indexOf(10, position);
        if (newline < 0) break;
        job.remaining -= 1;
        this.rows += 1;
        position = newline + 1;
        end = position;
      }
      if (end < 0) {
        // A partial row: keep it for the next chunk.
        job.parts.push(buffer.subarray(start));
        return;
      }
      job.parts.push(buffer.subarray(start, end));
      start = end;
      if (job.remaining === 0) {
        this.jobs.shift();
        const body = job.parts.length === 1 ? job.parts[0] : Buffer.concat(job.parts);
        job.resolve({ header: this.header, body: job.bytes ? body : body.toString("utf8"), count: job.count });
      }
    }
  }

  /**
   * Queue `count` records (`text`). AIRR input repeats its header line in
   * every batch; the child reads it once per stream. With `bytes`, the rows
   * come back as a Buffer.
   */
  run({ text, count, bytes = false }) {
    if (this.failed) return refuse(this.failed);
    if (this.ended) return refuse(new Error("The native aligner's input is closed."));
    if (!text.startsWith(">") && !text.startsWith("@")) {
      const newline = text.indexOf("\n");
      const header = text.slice(0, newline);
      if (this.airrHeader === null) this.airrHeader = header;
      else if (header !== this.airrHeader) return refuse(new Error("The AIRR input changed its columns between batches."));
      else text = text.slice(newline + 1);
    }
    this.records += count;
    const promise = new Promise((resolve, reject) => {
      this.jobs.push({ count, remaining: count, parts: [], resolve, reject, bytes });
      if (!this.child.stdin.write(text)) {
        // Resolves on drain, or when the pipe closes because the child died.
        this.drain ??= new Promise((done) => {
          const finish = () => {
            this.child.stdin.off("drain", finish);
            this.child.stdin.off("close", finish);
            this.drain = null;
            done();
          };
          this.child.stdin.on("drain", finish);
          this.child.stdin.on("close", finish);
        });
      }
    });
    // The caller awaits batches later, in order; a failure must not go unhandled meanwhile.
    promise.catch(() => {});
    return promise;
  }

  /** Back-pressure for the producer: waits for stdin to drain, and throws the child's failure. */
  async ready() {
    if (this.drain) await this.drain;
    if (this.child.stdin.destroyed && !this.ended) await this.exited;
    if (this.failed) throw this.failed;
  }

  /** No more records: the child assigns what it holds and exits. */
  finishInput() {
    if (!this.ended) {
      this.ended = true;
      this.child?.stdin.end();
    }
  }

  /**
   * After the caller has taken every batch: waits for the child to exit and
   * throws unless it wrote exactly one row per record. A surplus row arrives
   * after the last batch was answered, and nothing else would report it.
   */
  async finish() {
    this.finishInput();
    await this.exited;
    if (this.failed) throw this.failed;
    if (this.rows !== this.records) {
      throw new Error(`The native aligner wrote ${this.rows.toLocaleString()} rows for ${this.records.toLocaleString()} records.`);
    }
  }

  async close() {
    this.finishInput();
    if (this.child && this.child.exitCode === null && this.child.signalCode === null && this.jobs.length) this.child.kill();
    if (this.exited) await this.exited;
    if (this.directory) rmSync(this.directory, { recursive: true, force: true });
  }
}
