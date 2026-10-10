export type SequenceFormat = 1 | 2 | 3;

export interface GzipMemberRange {
  start: number;
  end: number;
}

/**
 * A gzip file whose independently compressed members must be decoded in
 * sequence. Browsers currently reject bytes after the first gzip member, so
 * the upload inspector records the validated member boundaries explicitly.
 */
export interface GzipMemberSource {
  kind: "gzip-members";
  file: File;
  members: GzipMemberRange[];
}

export type SequenceSource = string | File | GzipMemberSource;

export function isGzipMemberSource(source: SequenceSource): source is GzipMemberSource {
  return typeof source !== "string" && !(source instanceof File) && source.kind === "gzip-members";
}

export function sequenceSourceSize(source: SequenceSource): number {
  if (typeof source === "string") return source.length;
  if (source instanceof File) return source.size;
  return source.members.reduce((total, member) => total + Math.max(0, member.end - member.start), 0);
}

/** Cheap workload proxy. Gzip ISIZE is modulo 2^32, so wrapped sizes use a rough expansion estimate. */
export async function estimateSequenceCharacters(source: SequenceSource): Promise<number> {
  if (typeof source === "string") return source.length;
  const grouped = isGzipMemberSource(source);
  const file = grouped ? source.file : source;
  if (!grouped && !file.name.toLowerCase().endsWith(".gz")) return file.size;
  const ranges = grouped ? source.members : [{ start: 0, end: file.size }];
  let total = 0;
  for (const range of ranges) {
    if (range.end - range.start < 18) { total += Math.max(1, range.end - range.start) * 4; continue; }
    const trailer = await file.slice(range.end - 4, range.end).arrayBuffer();
    const size = new DataView(trailer).getUint32(0, true);
    total += size > 0 && range.end - range.start < 512 * 1024 * 1024 ? size : (range.end - range.start) * 4;
  }
  return total;
}

export interface SequenceBatch {
  index: number;
  text: string;
  count: number;
  format: SequenceFormat;
}

export interface SequenceStreamProgress {
  charactersRead: number;
  totalCharacters: number;
  bytesRead: number;
  totalBytes: number;
  recordsRead: number;
  recordsEligible: number;
  recordsSelected: number;
  maxBatchCharacters: number;
  maxCarryCharacters: number;
  fastqFilter: FastqQualityFilterStats;
}

export interface SequenceSubsample {
  size: number;
  seed: number;
}

export interface FastqEndTrimOptions {
  enabled: boolean;
  windowSize: number;
  minimumMeanPhred: number;
  minimumLength: number;
}

export interface FastqQualityFilterOptions {
  enabled: boolean;
  maximumExpectedErrors: number;
  phredOffset: 33 | 64;
  trim3Prime: FastqEndTrimOptions;
}

export interface FastqQualityFilterStats {
  enabled: boolean;
  applicable: boolean;
  recordsEvaluated: number;
  recordsRetained: number;
  recordsPassedThrough: number;
  recordsRejectedExpectedErrors: number;
  recordsRejectedMinimumLength: number;
  recordsTrimmed: number;
  basesTrimmed: number;
}

export const DEFAULT_FASTQ_QUALITY_FILTER: FastqQualityFilterOptions = {
  enabled: false,
  maximumExpectedErrors: 0.01,
  phredOffset: 33,
  trim3Prime: {
    enabled: false,
    windowSize: 4,
    minimumMeanPhred: 20,
    minimumLength: 50,
  },
};

export function emptyFastqQualityFilterStats(
  enabled = false,
  applicable = false,
): FastqQualityFilterStats {
  return {
    enabled,
    applicable,
    recordsEvaluated: 0,
    recordsRetained: 0,
    recordsPassedThrough: 0,
    recordsRejectedExpectedErrors: 0,
    recordsRejectedMinimumLength: 0,
    recordsTrimmed: 0,
    basesTrimmed: 0,
  };
}

export function addFastqQualityFilterStats(
  target: FastqQualityFilterStats,
  source: FastqQualityFilterStats,
): FastqQualityFilterStats {
  return {
    enabled: target.enabled || source.enabled,
    applicable: target.applicable || source.applicable,
    recordsEvaluated: target.recordsEvaluated + source.recordsEvaluated,
    recordsRetained: target.recordsRetained + source.recordsRetained,
    recordsPassedThrough: target.recordsPassedThrough + source.recordsPassedThrough,
    recordsRejectedExpectedErrors: target.recordsRejectedExpectedErrors + source.recordsRejectedExpectedErrors,
    recordsRejectedMinimumLength: target.recordsRejectedMinimumLength + source.recordsRejectedMinimumLength,
    recordsTrimmed: target.recordsTrimmed + source.recordsTrimmed,
    basesTrimmed: target.basesTrimmed + source.basesTrimmed,
  };
}

export interface SequenceStreamOptions {
  source: SequenceSource;
  format: SequenceFormat;
  batchSize?: number;
  subsample?: SequenceSubsample;
  fastqFilter?: FastqQualityFilterOptions;
  signal?: AbortSignal;
  onProgress?: (progress: SequenceStreamProgress) => void;
}

const STRING_CHUNK_SIZE = 256 * 1024;

function abortIfNeeded(signal?: AbortSignal) {
  if (signal?.aborted) throw new DOMException("Analysis cancelled.", "AbortError");
}

async function* decodedChunks(
  source: SequenceSource,
  signal: AbortSignal | undefined,
  onBytes: (bytesRead: number, totalBytes: number) => void,
): AsyncGenerator<string> {
  if (typeof source === "string") {
    const total = source.length;
    for (let offset = 0; offset < total; offset += STRING_CHUNK_SIZE) {
      abortIfNeeded(signal);
      const next = Math.min(offset + STRING_CHUNK_SIZE, total);
      yield source.slice(offset, next);
      onBytes(next, total);
      await Promise.resolve();
    }
    return;
  }

  const concatenatedGzip = isGzipMemberSource(source);
  const gzipMembers = concatenatedGzip ? source.members : null;
  const file = concatenatedGzip ? source.file : source;
  const ranges = gzipMembers ?? [{ start: 0, end: file.size }];
  if (!ranges.length) throw new Error("The concatenated gzip source has no members.");
  let previousEnd = -1;
  for (const range of ranges) {
    if (
      !Number.isSafeInteger(range.start)
      || !Number.isSafeInteger(range.end)
      || range.start < 0
      || range.end <= range.start
      || range.end > file.size
      || range.start < previousEnd
    ) throw new Error("The concatenated gzip source has invalid member boundaries.");
    previousEnd = range.end;
  }

  const total = ranges.reduce((sum, range) => sum + range.end - range.start, 0);
  let bytesRead = 0;
  const decoder = new TextDecoder();
  for (const range of ranges) {
    const measured = file.slice(range.start, range.end).stream().pipeThrough(new TransformStream<Uint8Array, Uint8Array>({
      transform(chunk, controller) {
        abortIfNeeded(signal);
        bytesRead += chunk.byteLength;
        onBytes(bytesRead, total);
        controller.enqueue(chunk);
      },
    }));
    const bytes = gzipMembers || file.name.toLowerCase().endsWith(".gz")
      ? measured.pipeThrough(makeGzipDecoder())
      : measured;
    const reader = bytes.getReader();
    try {
      while (true) {
        abortIfNeeded(signal);
        const { done, value } = await reader.read();
        if (done) break;
        const text = decoder.decode(value, { stream: true });
        if (text) yield text;
      }
    } finally {
      reader.releaseLock();
    }
  }
  const tail = decoder.decode();
  if (tail) yield tail;
}

function makeGzipDecoder(): TransformStream<Uint8Array, Uint8Array> {
  if (!("DecompressionStream" in globalThis)) {
    throw new Error("This browser cannot stream gzip input. Decompress the file first.");
  }
  return new DecompressionStream("gzip") as TransformStream<Uint8Array, Uint8Array>;
}

async function* lines(
  source: SequenceSource,
  signal: AbortSignal | undefined,
  onBytes: (bytesRead: number, totalBytes: number) => void,
  onCarry: (characters: number) => void,
  onConsumed?: (characters: number) => void,
): AsyncGenerator<string> {
  let carry = "";
  for await (const chunk of decodedChunks(source, signal, onBytes)) {
    carry += chunk;
    let start = 0;
    while (true) {
      const end = carry.indexOf("\n", start);
      if (end < 0) break;
      onConsumed?.(end + 1 - start);
      yield carry.slice(start, end).replace(/\r$/, "");
      start = end + 1;
    }
    carry = carry.slice(start);
    onCarry(carry.length);
  }
  if (carry) { onConsumed?.(carry.length); yield carry.replace(/\r$/, ""); }
}

export async function* fastaRecords(source: AsyncIterable<string>): AsyncGenerator<string> {
  let header = "";
  let sequence: string[] = [];
  for await (const line of source) {
    if (line.startsWith(">")) {
      if (header) {
        const joined = sequence.join("");
        if (!joined) throw new Error(`Empty FASTA record: ${header.slice(1).trim() || "unnamed sequence"}.`);
        yield `${header}\n${joined}\n`;
      }
      header = line;
      sequence = [];
      continue;
    }
    if (!line.trim()) continue;
    if (!header) throw new Error("Expected a FASTA header beginning with '>'.");
    sequence.push(line.replace(/\s/g, ""));
  }
  if (header) {
    const joined = sequence.join("");
    if (!joined) throw new Error(`Empty FASTA record: ${header.slice(1).trim() || "unnamed sequence"}.`);
    yield `${header}\n${joined}\n`;
  }
}

export interface FastqRecord {
  header: string;
  sequence: string;
  plus: string;
  quality: string;
}

export async function* fastqRecords(source: AsyncIterable<string>): AsyncGenerator<FastqRecord> {
  let header = "";
  let plus = "";
  let sequence: string[] = [];
  let quality: string[] = [];
  let sequenceLength = 0;
  let qualityLength = 0;
  let state: "header" | "sequence" | "quality" = "header";

  for await (const line of source) {
    if (state === "header") {
      if (!line.trim()) continue;
      if (!line.startsWith("@")) throw new Error("Expected a FASTQ header beginning with '@'.");
      header = line;
      plus = "";
      sequence = [];
      quality = [];
      sequenceLength = 0;
      qualityLength = 0;
      state = "sequence";
      continue;
    }
    if (state === "sequence") {
      if (line.startsWith("+")) {
        if (!sequenceLength) throw new Error(`Empty FASTQ sequence: ${header.slice(1).trim() || "unnamed sequence"}.`);
        plus = line;
        state = "quality";
      } else {
        const normalized = line.replace(/\s/g, "");
        sequence.push(normalized);
        sequenceLength += normalized.length;
      }
      continue;
    }

    quality.push(line);
    qualityLength += line.length;
    if (qualityLength > sequenceLength) {
      throw new Error(`FASTQ sequence/quality length mismatch: ${header.slice(1).trim() || "unnamed sequence"}.`);
    }
    if (qualityLength === sequenceLength) {
      yield { header, sequence: sequence.join(""), plus, quality: quality.join("") };
      state = "header";
    }
  }

  if (state === "sequence") throw new Error(`The FASTQ input ended before a '+' line for ${header.slice(1).trim() || "a sequence"}.`);
  if (state === "quality") throw new Error(`The FASTQ quality string is truncated for ${header.slice(1).trim() || "a sequence"}.`);
}

export async function* airrRecords(source: AsyncIterable<string>): AsyncGenerator<{ header: string; row: string }> {
  let header = "";
  let delimiter = "\t";
  let sequenceColumn = -1;
  for await (const line of source) {
    if (!header) {
      if (!line.trim()) continue;
      header = line;
      delimiter = header.includes("\t") ? "\t" : ",";
      sequenceColumn = header.split(delimiter).indexOf("sequence");
      if (sequenceColumn < 0) throw new Error("AIRR input requires a 'sequence' column.");
      continue;
    }
    if (!line.trim()) continue;
    const values = line.split(delimiter);
    // A blank sequence is no record, as for the aligner, which counts on one row per record.
    if (values[sequenceColumn]?.trim()) yield { header, row: line };
  }
}

interface StreamRecord {
  ordinal: number;
  text: string;
  header?: string;
}

function expectedErrorTable(offset: 33 | 64): Float64Array {
  const table = new Float64Array(127);
  for (let code = offset; code < table.length; code += 1) {
    table[code] = 10 ** (-(code - offset) / 10);
  }
  return table;
}

const PHRED33_EXPECTED_ERRORS = expectedErrorTable(33);
const PHRED64_EXPECTED_ERRORS = expectedErrorTable(64);

export function canonicalFastq(record: FastqRecord, end = record.sequence.length): string {
  if (end === record.sequence.length) {
    return `${record.header}\n${record.sequence}\n${record.plus}\n${record.quality}\n`;
  }
  return `${record.header}\n${record.sequence.slice(0, end)}\n${record.plus}\n${record.quality.slice(0, end)}\n`;
}

/**
 * Apply the FASTQ filter without constructing an output record for rejected
 * reads. The full-read expected-error sum is accumulated once; when 3' trim
 * removes a base, its contribution is subtracted before the threshold test.
 */
export function filterFastqRecord(
  record: FastqRecord,
  options: FastqQualityFilterOptions,
  stats: FastqQualityFilterStats,
): string | null {
  stats.recordsEvaluated += 1;
  const offset = options.phredOffset;
  const errorTable = offset === 64 ? PHRED64_EXPECTED_ERRORS : PHRED33_EXPECTED_ERRORS;
  const quality = record.quality;
  const trim = options.trim3Prime;
  const windowSize = Math.max(1, Math.floor(trim.windowSize));
  let end = quality.length;
  const initialWindowStart = Math.max(0, end - windowSize);
  let windowStart = initialWindowStart;
  let windowPhredSum = 0;
  let expectedErrors = 0;

  for (let index = 0; index < quality.length; index += 1) {
    const code = quality.charCodeAt(index);
    if (code < offset || code >= errorTable.length) {
      const name = record.header.slice(1).trim() || "unnamed sequence";
      throw new Error(`FASTQ quality character outside Phred+${offset} range in ${name} at base ${index + 1}.`);
    }
    expectedErrors += errorTable[code];
    if (trim.enabled && index >= initialWindowStart) windowPhredSum += code - offset;
  }

  if (trim.enabled) {
    while (
      end > 0
      && windowPhredSum < trim.minimumMeanPhred * (end - windowStart)
    ) {
      const removedCode = quality.charCodeAt(end - 1);
      expectedErrors -= errorTable[removedCode];
      windowPhredSum -= removedCode - offset;
      end -= 1;
      const nextWindowStart = Math.max(0, end - windowSize);
      if (nextWindowStart < windowStart) {
        const addedCode = quality.charCodeAt(nextWindowStart);
        windowPhredSum += addedCode - offset;
      }
      windowStart = nextWindowStart;
    }
    if (end < Math.max(1, Math.floor(trim.minimumLength))) {
      stats.recordsRejectedMinimumLength += 1;
      return null;
    }
  }

  if (Math.max(0, expectedErrors) > options.maximumExpectedErrors) {
    stats.recordsRejectedExpectedErrors += 1;
    return null;
  }
  const trimmedBases = quality.length - end;
  if (trimmedBases) {
    stats.recordsTrimmed += 1;
    stats.basesTrimmed += trimmedBases;
  }
  stats.recordsRetained += 1;
  return canonicalFastq(record, end);
}

export function seededRandom(seed: number): () => number {
  let value = Math.trunc(seed) >>> 0;
  return () => {
    value = (value + 0x6d2b79f5) >>> 0;
    let mixed = value;
    mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
    mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
    return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
  };
}

export function validateFastqQualityFilter(options: FastqQualityFilterOptions): void {
  if (!options.enabled) return;
  if (!Number.isFinite(options.maximumExpectedErrors) || options.maximumExpectedErrors < 0) {
    throw new Error("Maximum FASTQ expected errors must be a non-negative number.");
  }
  if (options.phredOffset !== 33 && options.phredOffset !== 64) {
    throw new Error("FASTQ quality encoding must be Phred+33 or Phred+64.");
  }
  if (options.trim3Prime.enabled && (
    !Number.isFinite(options.trim3Prime.windowSize) || options.trim3Prime.windowSize < 1
    || !Number.isFinite(options.trim3Prime.minimumMeanPhred) || options.trim3Prime.minimumMeanPhred < 0
    || !Number.isFinite(options.trim3Prime.minimumLength) || options.trim3Prime.minimumLength < 1
  )) {
    throw new Error("FASTQ 3' trimming requires a positive window and retained length, and a non-negative mean Phred threshold.");
  }
}

export async function* streamSequenceBatches(options: SequenceStreamOptions): AsyncGenerator<SequenceBatch> {
  const batchSize = Math.max(1, Math.floor(options.batchSize ?? 1000));
  const subsampleSize = options.subsample
    ? Math.max(1, Math.floor(options.subsample.size))
    : 0;
  let bytesRead = 0;
  let totalBytes = sequenceSourceSize(options.source);
  const totalCharacters = await estimateSequenceCharacters(options.source);
  let charactersRead = 0;
  let recordsRead = 0;
  let recordsEligible = 0;
  let recordsSelected = 0;
  let maxBatchCharacters = 0;
  let maxCarryCharacters = 0;
  const filterOptions = options.fastqFilter ?? DEFAULT_FASTQ_QUALITY_FILTER;
  validateFastqQualityFilter(filterOptions);
  const fastqFilter = emptyFastqQualityFilterStats(
    filterOptions.enabled,
    filterOptions.enabled && options.format === 2,
  );
  const report = () => options.onProgress?.({
    charactersRead, totalCharacters,
    bytesRead,
    totalBytes,
    recordsRead,
    recordsEligible,
    recordsSelected,
    maxBatchCharacters,
    maxCarryCharacters,
    fastqFilter: { ...fastqFilter },
  });
  const inputLines = lines(
    options.source,
    options.signal,
    (read, total) => {
      bytesRead = read;
      totalBytes = total;
      report();
    },
    (characters) => {
      maxCarryCharacters = Math.max(maxCarryCharacters, characters);
    },
    (characters) => { charactersRead += characters; },
  );

  const parsedRecords = async function* (): AsyncGenerator<StreamRecord> {
    if (options.format === 3) {
      for await (const record of airrRecords(inputLines)) {
        abortIfNeeded(options.signal);
        recordsRead += 1;
        recordsEligible += 1;
        if (filterOptions.enabled) {
          fastqFilter.recordsRetained += 1;
          fastqFilter.recordsPassedThrough += 1;
        }
        recordsSelected = subsampleSize ? Math.min(recordsEligible, subsampleSize) : recordsEligible;
        if (recordsRead % 1000 === 0) report();
        yield { ordinal: recordsRead - 1, text: `${record.row}\n`, header: record.header };
      }
      return;
    }
    if (options.format === 1) {
      for await (const record of fastaRecords(inputLines)) {
        abortIfNeeded(options.signal);
        recordsRead += 1;
        recordsEligible += 1;
        if (filterOptions.enabled) {
          fastqFilter.recordsRetained += 1;
          fastqFilter.recordsPassedThrough += 1;
        }
        recordsSelected = subsampleSize ? Math.min(recordsEligible, subsampleSize) : recordsEligible;
        if (recordsRead % 1000 === 0) report();
        yield { ordinal: recordsRead - 1, text: record };
      }
      return;
    }
    for await (const record of fastqRecords(inputLines)) {
      abortIfNeeded(options.signal);
      recordsRead += 1;
      const text = filterOptions.enabled
        ? filterFastqRecord(record, filterOptions, fastqFilter)
        : canonicalFastq(record);
      if (text !== null) recordsEligible += 1;
      recordsSelected = subsampleSize ? Math.min(recordsEligible, subsampleSize) : recordsEligible;
      if (recordsRead % 1000 === 0) report();
      if (text !== null) yield { ordinal: recordsRead - 1, text };
    }
  };

  let selected: AsyncIterable<StreamRecord> = parsedRecords();
  if (subsampleSize) {
    const random = seededRandom(options.subsample?.seed ?? 1);
    const reservoir: StreamRecord[] = [];
    for await (const record of selected) {
      if (reservoir.length < subsampleSize) {
        reservoir.push(record);
      } else {
        const replacement = Math.floor(random() * recordsEligible);
        if (replacement < subsampleSize) reservoir[replacement] = record;
      }
      recordsSelected = Math.min(recordsEligible, subsampleSize);
    }
    reservoir.sort((a, b) => a.ordinal - b.ordinal);
    selected = (async function* () { yield* reservoir; })();
    report();
  }

  let batchIndex = 0;
  let records: StreamRecord[] = [];
  const emit = (): SequenceBatch => {
    const header = options.format === 3 ? records[0]?.header : undefined;
    const body = records.map((record) => record.text).join("");
    const text = header ? `${header}\n${body}` : body;
    maxBatchCharacters = Math.max(maxBatchCharacters, text.length);
    const batch = { index: batchIndex++, text, count: records.length, format: options.format };
    records = [];
    report();
    return batch;
  };

  for await (const record of selected) {
    abortIfNeeded(options.signal);
    records.push(record);
    recordsSelected = subsampleSize ? Math.min(recordsEligible, subsampleSize) : recordsEligible;
    if (records.length === batchSize) yield emit();
  }
  if (records.length) yield emit();

  bytesRead = totalBytes;
  recordsSelected = subsampleSize ? Math.min(recordsEligible, subsampleSize) : recordsEligible;
  report();
  if (!recordsRead) throw new Error("No sequence records were found in the input.");
  if (!recordsEligible && fastqFilter.applicable) {
    throw new Error(`The FASTQ quality filter rejected all ${recordsRead.toLocaleString()} input reads.`);
  }
}
