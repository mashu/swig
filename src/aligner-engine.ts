// The Swig aligner's WebAssembly module (engine/, Rust; public/swig.wasm)
// as one JavaScript object per worker. The browser's compute workers and
// the command line's workers share this wrapper and its option mapping.

export type CallingProfile = "truth_optimized" | "igblast_compatible" | "r_optimized" | "sensitive_d";

export const CALLING_PROFILES: readonly CallingProfile[] = [
  "r_optimized",
  "truth_optimized",
  "igblast_compatible",
  "sensitive_d",
];

export type DoubleDScreenMode = "off" | "all" | "long_span";

export interface DoubleDScreenOptions {
  mode: DoubleDScreenMode;
  minimumVjSpan: number;
  seedLength: number;
  pseudoTrim: number;
  maximumPseudoMismatches: number;
  minimumScoreGain: number;
}

/** IgBLAST-style overrides of the profile's D and J settings (`swig-cli --vdj`). */
export interface AlignerTuning {
  minimumDMatch?: number;
  dMismatch?: number;
  jMismatch?: number;
  minimumJLength?: number;
}

export interface AlignerSettings {
  callingProfile: CallingProfile;
  /** V, D and J hits below this identity are dropped (0–1). */
  minimumIdentity: number;
  /** Report C only where its germline from the first base reaches this identity. */
  minimumConstantPrefixIdentity?: number | null;
  /** 0 both strands, 1 the read as given, 2 its reverse complement. */
  strand: 0 | 1 | 2;
  doubleD?: DoubleDScreenOptions;
  tuning?: AlignerTuning;
}

export interface GermlineTexts {
  V: string;
  D: string;
  J: string;
  C: string;
}

export interface AnnotatedBatch {
  header: string;
  body: Uint8Array;
  count: number;
  doubleDHeader?: string;
  doubleDBody?: Uint8Array;
  doubleDCount?: number;
}

export function isCallingProfile(value: unknown): value is CallingProfile {
  return typeof value === "string" && (CALLING_PROFILES as readonly string[]).includes(value);
}

/** A saved calling profile, or the default where it names none the aligner has. */
export function callingProfileOrDefault(value: unknown): CallingProfile {
  return isCallingProfile(value) ? value : "r_optimized";
}

/** The settings as the engine's options, by the names of the native `swig` command line. */
export function alignerOptions(settings: AlignerSettings): Array<[string, string]> {
  const options: Array<[string, string]> = [
    ["profile", settings.callingProfile],
    ["strand", settings.strand === 1 ? "forward" : settings.strand === 2 ? "reverse" : "both"],
    ["min-identity", String(settings.minimumIdentity)],
  ];
  if (settings.minimumConstantPrefixIdentity != null) {
    options.push(["c-prefix-identity", String(settings.minimumConstantPrefixIdentity)]);
  }
  const doubleD = settings.doubleD;
  if (doubleD && doubleD.mode !== "off") {
    options.push(
      ["double-d", doubleD.mode],
      ["double-d-min-span", String(Math.round(doubleD.minimumVjSpan))],
      ["double-d-seed-length", String(Math.round(doubleD.seedLength))],
      ["double-d-pseudo-trim", String(Math.round(doubleD.pseudoTrim))],
      ["double-d-pseudo-mismatches", String(Math.round(doubleD.maximumPseudoMismatches))],
      ["double-d-min-gain", String(Math.round(doubleD.minimumScoreGain))],
    );
  }
  const tuning = settings.tuning ?? {};
  if (tuning.minimumDMatch != null) options.push(["min-d-match", String(Math.round(tuning.minimumDMatch))]);
  if (tuning.dMismatch != null) options.push(["d-mismatch", String(Math.round(tuning.dMismatch))]);
  if (tuning.jMismatch != null) options.push(["j-mismatch", String(Math.round(tuning.jMismatch))]);
  if (tuning.minimumJLength != null) options.push(["min-j-length", String(Math.round(tuning.minimumJLength))]);
  return options;
}

interface AlignerExports extends WebAssembly.Exports {
  memory: WebAssembly.Memory;
  swig_alloc: (length: number) => number;
  swig_free: (pointer: number, length: number) => void;
  swig_set_option: (name: number, nameLength: number, value: number, valueLength: number) => number;
  swig_init_database: (
    v: number, vLength: number, d: number, dLength: number,
    j: number, jLength: number, c: number, cLength: number,
  ) => number;
  swig_annotate: (reads: number, length: number) => number;
  swig_result_ptr: () => number;
  swig_result_len: () => number;
  swig_double_d_result_ptr: () => number;
  swig_double_d_result_len: () => number;
  swig_double_d_count: () => number;
  swig_error_ptr: () => number;
  swig_error_len: () => number;
  swig_panic_message: () => void;
  swig_version: () => void;
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();

/** Split `table` (header line, then rows) into its header and rows. */
function splitTable(table: Uint8Array, what: string): { header: string; body: Uint8Array } {
  const newline = table.indexOf(10);
  if (newline < 0) throw new Error(`The Swig aligner returned an invalid ${what}.`);
  return {
    header: decoder.decode(table.subarray(0, newline)).replace(/\r$/, ""),
    body: table.slice(newline + 1),
  };
}

export class AlignerEngine {
  /** The database summary and warnings from loading the germlines. */
  readonly summary: { summary: string; warnings: string[] };
  private readonly exports: AlignerExports;

  private constructor(exports: AlignerExports, summary: { summary: string; warnings: string[] }) {
    this.exports = exports;
    this.summary = summary;
  }

  /** One engine over `germlines` with `settings`, from the compiled module. */
  static async create(
    module: WebAssembly.Module,
    settings: AlignerSettings,
    germlines: GermlineTexts,
  ): Promise<{ engine: AlignerEngine; records: number }> {
    const instance = await WebAssembly.instantiate(module, {});
    const exports = instance.exports as AlignerExports;
    const call = <T>(action: () => T): T => {
      try {
        return action();
      } catch (error) {
        // A panic traps the instance; its message is still readable.
        exports.swig_panic_message();
        const message = AlignerEngine.read(exports, exports.swig_error_ptr(), exports.swig_error_len());
        throw new Error(message || (error instanceof Error ? error.message : String(error)));
      }
    };
    for (const [name, value] of alignerOptions(settings)) {
      const nameBytes = AlignerEngine.put(exports, encoder.encode(name));
      const valueBytes = AlignerEngine.put(exports, encoder.encode(value));
      try {
        if (call(() => exports.swig_set_option(...nameBytes, ...valueBytes)) !== 0) {
          throw new Error(AlignerEngine.error(exports));
        }
      } finally {
        exports.swig_free(...nameBytes);
        exports.swig_free(...valueBytes);
      }
    }
    const texts = [germlines.V, germlines.D, germlines.J, germlines.C].map((text) =>
      AlignerEngine.put(exports, encoder.encode(text || "")));
    try {
      const records = call(() => exports.swig_init_database(...texts.flat() as [
        number, number, number, number, number, number, number, number,
      ]));
      if (records < 0) throw new Error(AlignerEngine.error(exports));
      const summary = JSON.parse(AlignerEngine.read(exports, exports.swig_result_ptr(), exports.swig_result_len())) as {
        summary: string;
        warnings: string[];
      };
      return { engine: new AlignerEngine(exports, summary), records };
    } finally {
      texts.forEach((buffer) => exports.swig_free(...buffer));
    }
  }

  /** The engine's version and alignment kernel. */
  version(): { version: string; backend: string } {
    this.exports.swig_version();
    return JSON.parse(AlignerEngine.read(this.exports, this.exports.swig_result_ptr(), this.exports.swig_result_len()));
  }

  /**
   * Assign `count` reads (FASTA, FASTQ or an AIRR table with its header line):
   * the AIRR TSV header and one row per read in input order, and the double-D
   * sidecar when the screen is on.
   */
  annotate(reads: Uint8Array, count: number, doubleD: boolean): AnnotatedBatch {
    const exports = this.exports;
    const input = AlignerEngine.put(exports, reads);
    let returned: number;
    try {
      returned = exports.swig_annotate(...input);
    } catch (error) {
      exports.swig_panic_message();
      throw new Error(AlignerEngine.error(exports) || (error instanceof Error ? error.message : String(error)));
    } finally {
      exports.swig_free(...input);
    }
    if (returned < 0) throw new Error(AlignerEngine.error(exports));
    if (returned !== count) {
      throw new Error(`The Swig aligner returned ${returned} records for a ${count}-record input batch.`);
    }
    const table = new Uint8Array(exports.memory.buffer, exports.swig_result_ptr(), exports.swig_result_len());
    const batch: AnnotatedBatch = { ...splitTable(table, "AIRR table"), count };
    if (doubleD) {
      const sidecar = new Uint8Array(
        exports.memory.buffer,
        exports.swig_double_d_result_ptr(),
        exports.swig_double_d_result_len(),
      );
      const { header, body } = splitTable(sidecar, "double-D evidence table");
      batch.doubleDHeader = header;
      batch.doubleDBody = body;
      batch.doubleDCount = exports.swig_double_d_count();
    }
    return batch;
  }

  private static put(exports: AlignerExports, bytes: Uint8Array): [number, number] {
    const pointer = exports.swig_alloc(bytes.byteLength);
    if (!pointer) throw new Error("The Swig aligner ran out of WebAssembly memory.");
    new Uint8Array(exports.memory.buffer, pointer, bytes.byteLength).set(bytes);
    return [pointer, bytes.byteLength];
  }

  private static read(exports: AlignerExports, pointer: number, length: number): string {
    return decoder.decode(new Uint8Array(exports.memory.buffer, pointer, length));
  }

  private static error(exports: AlignerExports): string {
    return AlignerEngine.read(exports, exports.swig_error_ptr(), exports.swig_error_len()) ||
      "The Swig aligner could not complete the request.";
  }
}
