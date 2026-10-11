import assert from "node:assert/strict";
import test from "node:test";
import { gzipSync } from "node:zlib";

import {
  DEFAULT_FASTQ_QUALITY_FILTER,
  isBlank,
  streamSequenceBatches,
  type FastqQualityFilterOptions,
  type SequenceSource,
  type SequenceStreamProgress,
} from "../src/sequence-stream.ts";
import {
  gzipMemberFile,
  gzipMemberSource,
  inspectGzipMembers,
  suggestedGzipMemberName,
} from "../src/gzip-members.ts";

async function collect(source: SequenceSource, format: 1 | 2 | 3, batchSize = 2) {
  const batches = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source,
    format,
    batchSize,
    onProgress: (next) => { progress = next; },
  })) batches.push(batch);
  return { batches, progress: progress as SequenceStreamProgress | null };
}

async function sampledIds(source: string, size: number, seed: number) {
  const ids: string[] = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source,
    format: 1,
    batchSize: 3,
    subsample: { size, seed },
    onProgress: (next) => { progress = next; },
  })) {
    ids.push(...[...batch.text.matchAll(/^>([^\n]+)/gm)].map((match) => match[1]));
  }
  return { ids, progress: progress as SequenceStreamProgress | null };
}

function qualityFilter(overrides: Partial<FastqQualityFilterOptions> = {}): FastqQualityFilterOptions {
  return {
    ...DEFAULT_FASTQ_QUALITY_FILTER,
    enabled: true,
    ...overrides,
    trim3Prime: {
      ...DEFAULT_FASTQ_QUALITY_FILTER.trim3Prime,
      ...overrides.trim3Prime,
    },
  };
}

async function filteredFastq(
  source: string | File,
  filter: FastqQualityFilterOptions,
  subsample?: { size: number; seed: number },
) {
  const batches = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source,
    format: 2,
    batchSize: 1000,
    fastqFilter: filter,
    subsample,
    onProgress: (next) => { progress = next; },
  })) batches.push(batch);
  return { batches, progress: progress as SequenceStreamProgress | null };
}

test("incremental FASTA parser normalizes wrapped records into bounded batches", async () => {
  const result = await collect(">one notes\nAC GT\nTG\n>two\nAAAA\n>three\nCCCC\n", 1);
  assert.deepEqual(result.batches.map((batch) => batch.count), [2, 1]);
  assert.equal(result.batches[0].text, ">one notes\nACGTTG\n>two\nAAAA\n");
  assert.equal(result.progress?.recordsRead, 3);
});

test("incremental FASTQ parser supports wrapped sequence and quality lines", async () => {
  const result = await collect(
    "@one notes\nACGT\nTG\n+one\nIIII\nII\n@two\nAAAA\n+\n@@@@\n",
    2,
    10,
  );
  assert.equal(result.batches[0].count, 2);
  assert.equal(result.batches[0].text, "@one notes\nACGTTG\n+one\nIIIIII\n@two\nAAAA\n+\n@@@@\n");
});

test("FASTQ expected-error threshold is inclusive and rejects the next lower quality", async () => {
  const sequence = "A".repeat(100);
  const result = await filteredFastq(
    `@q40\n${sequence}\n+\n${"I".repeat(100)}\n@q39\n${sequence}\n+\n${"H".repeat(100)}\n`,
    qualityFilter({ maximumExpectedErrors: 0.01 }),
  );
  assert.equal(result.batches.length, 1);
  assert.match(result.batches[0].text, /^@q40$/m);
  assert.doesNotMatch(result.batches[0].text, /^@q39$/m);
  assert.equal(result.progress?.recordsRead, 2);
  assert.equal(result.progress?.recordsEligible, 1);
  assert.equal(result.progress?.fastqFilter.recordsEvaluated, 2);
  assert.equal(result.progress?.fastqFilter.recordsRejectedExpectedErrors, 1);
});

test("3-prime trimming happens before expected-error filtering", async () => {
  const source = "@tail\nACGTACGTACGT\n+\nIIIIIIII!!!!\n";
  const withoutTrim = await assert.rejects(
    async () => filteredFastq(source, qualityFilter({ maximumExpectedErrors: 0.001 })),
    /rejected all 1 input reads/,
  );
  assert.equal(withoutTrim, undefined);

  const trimmed = await filteredFastq(source, qualityFilter({
    maximumExpectedErrors: 0.001,
    trim3Prime: { enabled: true, windowSize: 4, minimumMeanPhred: 40, minimumLength: 8 },
  }));
  assert.equal(trimmed.batches[0].text, "@tail\nACGTACGT\n+\nIIIIIIII\n");
  assert.equal(trimmed.progress?.fastqFilter.recordsTrimmed, 1);
  assert.equal(trimmed.progress?.fastqFilter.basesTrimmed, 4);
  assert.equal(trimmed.progress?.fastqFilter.recordsRejectedExpectedErrors, 0);
});

test("3-prime trimming rejects reads shorter than the configured retained length", async () => {
  await assert.rejects(
    async () => filteredFastq(
      "@short_after_trim\nACGTACGTACGT\n+\nIIIIIIII!!!!\n",
      qualityFilter({
        maximumExpectedErrors: 10,
        trim3Prime: { enabled: true, windowSize: 4, minimumMeanPhred: 20, minimumLength: 11 },
      }),
    ),
    /rejected all 1 input reads/,
  );
});

test("enabled FASTQ filtering passes FASTA through unchanged and reports the bypass", async () => {
  const batches = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source: ">one\nACGT\n>two\nTGCA\n",
    format: 1,
    batchSize: 10,
    fastqFilter: qualityFilter(),
    onProgress: (next) => { progress = next; },
  })) batches.push(batch);
  assert.equal(batches[0].text, ">one\nACGT\n>two\nTGCA\n");
  assert.equal(progress?.fastqFilter.applicable, false);
  assert.equal(progress?.fastqFilter.recordsEvaluated, 0);
  assert.equal(progress?.fastqFilter.recordsPassedThrough, 2);
  assert.equal(progress?.fastqFilter.recordsRetained, 2);
});

test("reservoir subsampling draws only from FASTQ reads that pass quality filtering", async () => {
  const source = Array.from({ length: 20 }, (_, index) => {
    const quality = index % 2 === 0 ? "IIII" : "!!!!";
    return `@read_${index}\nACGT\n+\n${quality}\n`;
  }).join("");
  const result = await filteredFastq(source, qualityFilter({ maximumExpectedErrors: 0.1 }), { size: 5, seed: 19 });
  const ids = result.batches.flatMap((batch) => [...batch.text.matchAll(/^@read_(\d+)$/gm)].map((match) => Number(match[1])));
  assert.equal(ids.length, 5);
  assert.ok(ids.every((id) => id % 2 === 0));
  assert.equal(result.progress?.recordsRead, 20);
  assert.equal(result.progress?.recordsEligible, 10);
  assert.equal(result.progress?.recordsSelected, 5);
});

test("gzip FASTQ is decompressed and parsed without whole-file text conversion", async () => {
  const source = "@gzip_one\nACGTACGT\n+\nIIIIIIII\n@gzip_two\nTGCATGCA\n+\nHHHHHHHH\n";
  const file = new File([gzipSync(source)], "reads.fastq.gz", { type: "application/gzip" });
  const result = await collect(file, 2, 1);
  assert.deepEqual(result.batches.map((batch) => batch.count), [1, 1]);
  assert.equal(result.progress?.bytesRead, file.size);
  assert.ok(result.batches[1].text.startsWith("@gzip_two\n"));
});

test("concatenated gzip datasets are detected, named, split, or decoded as one bounded source", async () => {
  const first = ">SRR2126754.1 1 length=4\nACGT\n>SRR2126754.2 2 length=4\nTGCA\n";
  const second = ">SRR2126755.1 1 length=4\nAAAA\n>SRR2126755.2 2 length=4\nCCCC\n";
  const firstGzip = gzipSync(first);
  const secondGzip = gzipSync(second);
  const file = new File([firstGzip, secondGzip], "SRR2126754.fasta.gz", { type: "application/gzip" });

  const members = await inspectGzipMembers(file, 1);
  assert.equal(members.length, 2);
  assert.deepEqual(members.map((member) => member.start), [0, firstGzip.byteLength]);
  assert.deepEqual(members.map((member) => member.end), [firstGzip.byteLength, file.size]);
  assert.deepEqual(members.map((member) => member.firstRecordId), ["SRR2126754.1", "SRR2126755.1"]);
  assert.ok(members.every((member) => member.startsAtRecordBoundary));
  assert.deepEqual(members.map((member, index) => suggestedGzipMemberName(file.name, member, index)), [
    "SRR2126754.fasta.gz",
    "SRR2126755.fasta.gz",
  ]);

  const merged = await collect(gzipMemberSource(file, members), 1, 3);
  assert.deepEqual(merged.batches.map((batch) => batch.count), [3, 1]);
  assert.equal(merged.progress?.bytesRead, file.size);
  assert.match(merged.batches.map((batch) => batch.text).join(""), /^>SRR2126754\.1/m);
  assert.match(merged.batches.map((batch) => batch.text).join(""), /^>SRR2126755\.2/m);

  const secondFile = gzipMemberFile(file, members[1], "SRR2126755.fasta.gz");
  const separate = await collect(secondFile, 1, 10);
  assert.equal(separate.batches[0].count, 2);
  assert.ok(separate.batches[0].text.startsWith(">SRR2126755.1 1 length=4\n"));
});

test("block gzip members that continue a record are merged without being presented as separate datasets", async () => {
  const first = gzipSync(">one\nAC");
  const second = gzipSync("GT\n>two\nTTAA\n");
  const file = new File([first, second], "blocked.fasta.gz", { type: "application/gzip" });
  const members = await inspectGzipMembers(file, 1);
  assert.equal(members.length, 2);
  assert.deepEqual(members.map((member) => member.startsAtRecordBoundary), [true, false]);
  const merged = await collect(gzipMemberSource(file, members), 1, 10);
  assert.equal(merged.batches[0].text, ">one\nACGT\n>two\nTTAA\n");
});

test("AIRR input retains its header in every compute batch", async () => {
  const source = "sequence_id\tsequence\tlocus\na\tACGT\tIGH\nb\tTGCA\tTRB\nc\t\tIGH\nd\tAAAA\tIGK\n";
  const result = await collect(source, 3, 2);
  assert.deepEqual(result.batches.map((batch) => batch.count), [2, 1]);
  assert.ok(result.batches.every((batch) => batch.text.startsWith("sequence_id\tsequence\tlocus\n")));
  assert.ok(!result.batches.some((batch) => batch.text.includes("\nc\t")));
});

test("blank means the aligner's Unicode whitespace, and FASTA records keep the aligner's boundaries", async () => {
  assert.ok(isBlank("") && isBlank(" \t\u00a0") && isBlank("\u0085") && isBlank(undefined));
  assert.ok(!isBlank("\ufeff") && !isBlank(" a "));
  // " >b" is a sequence line to the aligner; without its whitespace it must not become a header.
  const result = await collect(">a\nACGT\n >b\nTTGG\n>c\nGGCC\n", 1, 10);
  assert.equal(result.batches[0].count, 2);
  assert.equal(result.batches[0].text, ">a\nACGTbTTGG\n>c\nGGCC\n");
  // Column names are matched trimmed, and a sequence of only U+0085 is blank.
  const airr = await collect("sequence_id\tsequence \na\tACGT\nnel\t\u0085\n", 3, 10);
  assert.deepEqual(airr.batches.map((batch) => batch.count), [1]);
});

test("AIRR rows whose sequence is only whitespace are no records, as for the aligner", async () => {
  // The aligner writes one row per record; a row it skips must not be counted in a batch.
  const source = "sequence_id\tsequence\na\tACGT\nblank\t   \nb\tTGCA\n";
  const result = await collect(source, 3, 2);
  assert.deepEqual(result.batches.map((batch) => batch.count), [2]);
  assert.ok(!result.batches.some((batch) => batch.text.includes("blank")));
});

test("seeded reservoir subsampling is exact, reproducible, and input-order stable", async () => {
  const source = Array.from({ length: 100 }, (_, index) => `>read_${index}\nACGT${index % 10}\n`).join("");
  const first = await sampledIds(source, 10, 73);
  const again = await sampledIds(source, 10, 73);
  const other = await sampledIds(source, 10, 74);
  assert.deepEqual(first.ids, again.ids);
  assert.notDeepEqual(first.ids, other.ids);
  assert.equal(first.ids.length, 10);
  assert.deepEqual(first.ids, [...first.ids].sort((a, b) => Number(a.slice(5)) - Number(b.slice(5))));
  assert.equal(first.progress?.recordsRead, 100);
  assert.equal(first.progress?.recordsSelected, 10);
});

test("AIRR reservoir sampling keeps the header and scans all source rows", async () => {
  const source = `sequence_id\tsequence\n${Array.from({ length: 25 }, (_, index) => `r${index}\tACGT`).join("\n")}\n`;
  const batches = [];
  let progress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source,
    format: 3,
    batchSize: 4,
    subsample: { size: 7, seed: 9 },
    onProgress: (next) => { progress = next; },
  })) batches.push(batch);
  assert.deepEqual(batches.map((batch) => batch.count), [4, 3]);
  assert.ok(batches.every((batch) => batch.text.startsWith("sequence_id\tsequence\n")));
  assert.equal(progress?.recordsRead, 25);
  assert.equal(progress?.recordsSelected, 7);
});

test("50k gzip FASTA stress path and 10k reservoir sample remain batch-bounded", async () => {
  const sequence = "ACGT".repeat(75);
  const records = Array.from({ length: 50_000 }, (_, index) => `>read_${index}\n${sequence}\n`).join("");
  const compressed = gzipSync(records, { level: 1 });
  const file = new File([compressed], "stress.fasta.gz", { type: "application/gzip" });
  let count = 0;
  let batches = 0;
  let finalProgress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source: file,
    format: 1,
    batchSize: 1000,
    onProgress: (next) => { finalProgress = next; },
  })) {
    count += batch.count;
    batches += 1;
    assert.ok(batch.count <= 1000);
    assert.ok(batch.text.length < 400_000);
  }
  assert.equal(count, 50_000);
  assert.equal(batches, 50);
  assert.equal(finalProgress?.bytesRead, file.size);
  assert.equal(finalProgress?.recordsRead, 50_000);
  assert.ok((finalProgress?.maxBatchCharacters ?? Infinity) < 400_000);
  assert.ok((finalProgress?.maxCarryCharacters ?? Infinity) < 10_000);

  let sampledCount = 0;
  let sampledBatches = 0;
  let sampledProgress: SequenceStreamProgress | null = null;
  for await (const batch of streamSequenceBatches({
    source: file,
    format: 1,
    batchSize: 1000,
    subsample: { size: 10_000, seed: 2024 },
    onProgress: (next) => { sampledProgress = next; },
  })) {
    sampledCount += batch.count;
    sampledBatches += 1;
    assert.ok(batch.count <= 1000);
  }
  assert.equal(sampledCount, 10_000);
  assert.equal(sampledBatches, 10);
  assert.equal(sampledProgress?.recordsRead, 50_000);
  assert.equal(sampledProgress?.recordsSelected, 10_000);
  assert.ok((sampledProgress?.maxBatchCharacters ?? Infinity) < 400_000);
});

test("50k-read FASTQ expected-error filtering remains linear and batch-bounded", async () => {
  const sequence = "ACGT".repeat(30);
  const highQuality = "I".repeat(sequence.length);
  const lowQuality = "!".repeat(sequence.length);
  const source = Array.from({ length: 50_000 }, (_, index) => `@q_${index}\n${sequence}\n+\n${index % 10 ? highQuality : lowQuality}\n`).join("");
  const started = performance.now();
  const result = await filteredFastq(source, qualityFilter({ maximumExpectedErrors: 0.02 }));
  const elapsed = performance.now() - started;
  assert.equal(result.progress?.recordsRead, 50_000);
  assert.equal(result.progress?.recordsEligible, 45_000);
  assert.equal(result.progress?.fastqFilter.recordsRejectedExpectedErrors, 5_000);
  assert.equal(result.batches.length, 45);
  assert.ok(result.batches.every((batch) => batch.count <= 1000));
  assert.ok((result.progress?.maxBatchCharacters ?? Infinity) < 300_000);
  assert.ok(elapsed < 8_000, `FASTQ filter took ${elapsed.toFixed(0)} ms`);
});
