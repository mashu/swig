import assert from "node:assert/strict";
import test from "node:test";

import {
  assignLineages,
  boundedEditProfile,
  chmmairraDistanceFromReference,
  chmmJunctionFlank,
  extendChmmAlignment,
  createExactDedupPlan,
  DenoiseAccumulator,
  deduplicate,
  expandSingleLinkage,
  lineageDoubleDMatches,
  minHashSketch,
  prepareReferenceMsa,
  poissonStrictUpperTail,
  queryRecords,
  finishExactDedupPlan,
  runDenoisePartitionJob,
  runExactDedupJob,
  runChmm,
  sequenceFingerprint,
  threadSequenceToMsa,
  type PostAnalysisRecord,
} from "../src/post-analysis-core.ts";
import { denoiseVdjSequence } from "../src/post-analysis-record.ts";
import { isReverseComplement, reverseComplement } from "../src/nucleotide.ts";

function record(
  ordinal: number,
  cdr3Nt: string,
  vCall = "IGHV1-2*01",
  jCall = "IGHJ4*02",
  sequence = `AAA${cdr3Nt}TTT`,
): PostAnalysisRecord {
  return {
    ordinal,
    sequenceId: `read_${ordinal}`,
    locus: "IGH",
    vCall,
    jCall,
    cdr3Nt,
    cdr3Aa: "CARDR",
    productive: true,
    sequenceFingerprint: sequenceFingerprint(sequence),
    trimmedFingerprint: sequenceFingerprint(sequence),
    trimmedSketch: minHashSketch(sequence),
  };
}

test("deduplication retains abundance and representative membership", () => {
  const records = [record(0, "TGTGCCAAA"), record(1, "TGTGCCAAA"), record(2, "TGTGCCTAA")];
  const result = deduplicate(records, "sequence");
  assert.equal(result.uniqueRecords, 2);
  assert.equal(result.collapsedRecords, 1);
  assert.deepEqual([...result.counts], [2, 0, 1]);
  assert.deepEqual([...result.representatives], [0, 0, 2]);
});

test("exact deduplication sums an existing AIRR duplicate_count", () => {
  const records = [record(0, "TGTGCCAAA"), record(1, "TGTGCCAAA")];
  records[0].inputCount = 7;
  records[1].inputCount = 3;
  const result = deduplicate(records, "sequence");
  assert.equal(result.inputAbundance, 10);
  assert.deepEqual([...result.counts], [10, 0]);
});

test("collapse partitions on C gene by default without using constant-tail sequence length", () => {
  const first = record(0, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", "AAAACCCCTAIL");
  const second = record(1, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", "AAAACCCCTAIL-LONGER");
  first.cCall = "IGHM*01";
  second.cCall = "IGHG1*03";
  first.trimmedFingerprint = sequenceFingerprint("AAAACCCC");
  second.trimmedFingerprint = sequenceFingerprint("AAAACCCC");
  const separated = deduplicate([first, second], "trimmed");
  assert.equal(separated.uniqueRecords, 2);
  second.cCall = "IGHM*02";
  const sameConstantGene = deduplicate([first, second], "trimmed");
  assert.equal(sameConstantGene.uniqueRecords, 1);
  second.cCall = "IGHG1*03";
  const ignored = deduplicate([first, second], "trimmed", "discard", "global", false);
  assert.equal(ignored.uniqueRecords, 1);
});

test("exact collapse discards unusable keys by default and can retain them unchanged", () => {
  const records = [record(0, "TGTGCCAAA"), record(1, "TGTGCCAAA"), record(2, "")];
  records[2].trimmedFingerprint = sequenceFingerprint("");
  records[2].inputCount = 9;
  const discarded = deduplicate(records, "trimmed");
  assert.equal(discarded.unresolvedRecords, 1);
  assert.equal(discarded.uniqueRecords, 1);
  assert.deepEqual([...discarded.counts], [2, 0, 0]);
  assert.equal(discarded.representatives[2], -1);
  const retained = deduplicate(records, "trimmed", "retain");
  assert.equal(retained.unresolvedRecords, 1);
  assert.equal(retained.uniqueRecords, 2);
  assert.equal(retained.counts[2], 9);
  assert.equal(retained.representatives[2], 2);
});

test("partition-worker exact collapse and denoising are byte-for-byte deterministic", () => {
  const records: PostAnalysisRecord[] = [];
  const sequences: string[] = [];
  for (let partition = 0; partition < 4; partition += 1) {
    const parent = `${"ACGT".repeat(14)}${"ACGT"[partition]}`;
    const child = `${parent.slice(0, 21)}${parent[21] === "A" ? "C" : "A"}${parent.slice(22)}`;
    for (let copy = 0; copy < 6; copy += 1) {
      const sequence = copy < 4 ? parent : child;
      const value = record(records.length, `TGTGCCAA${partition}`, `IGHV${partition + 1}-1*01`, "IGHJ4*02", sequence);
      value.inputCount = copy === 0 ? 30 : 1;
      records.push(value);
      sequences.push(sequence);
    }
  }
  const serialExact = deduplicate(records, "trimmed", "discard", "global", true);
  const exactPlan = createExactDedupPlan(records, "trimmed", "discard", "global", true, 4);
  const parallelExact = finishExactDedupPlan(exactPlan, exactPlan.jobs.map(runExactDedupJob).reverse());
  assert.deepEqual([...parallelExact.representatives], [...serialExact.representatives]);
  assert.deepEqual([...parallelExact.counts], [...serialExact.counts]);
  assert.deepEqual(parallelExact.largestGroups, serialExact.largestGroups);

  const options = {
    mode: "conservative" as const, errorRate: 0.00473, alpha: 0.01, callResolution: "allele" as const,
    ambiguity: "strict" as const, minimumParentCount: 2, ambiguousPolicy: "retain" as const,
    fadNeighborThreshold: 1, fadMethod: 2 as const, expectedZeroErrorFraction: 1,
    maximumHammingDistance: 1, maximumEditDistance: 2, minimumIndelParentRatio: 2,
    maxCandidatesPerVariant: 10_000,
  };
  const serialAccumulator = new DenoiseAccumulator(records, options);
  const partitionedAccumulator = new DenoiseAccumulator(records, options);
  sequences.forEach((sequence, ordinal) => { serialAccumulator.add(ordinal, sequence); partitionedAccumulator.add(ordinal, sequence); });
  const serialDenoise = serialAccumulator.finish();
  const jobs = partitionedAccumulator.preparePartitionJobs();
  assert.equal(jobs.length, 4);
  const partitionedDenoise = partitionedAccumulator.finishWithPartitionResults(jobs.map(runDenoisePartitionJob).reverse());
  assert.deepEqual([...partitionedDenoise.representatives], [...serialDenoise.representatives]);
  assert.deepEqual([...partitionedDenoise.counts], [...serialDenoise.counts]);
  assert.equal(partitionedDenoise.candidateComparisons, serialDenoise.candidateComparisons);
  assert.deepEqual(partitionedDenoise.largestGroups, serialDenoise.largestGroups);
});

test("Poisson strict upper tail matches known probabilities used by FAD method 2", () => {
  assert.ok(Math.abs(poissonStrictUpperTail(0, 0.1) - 0.095162581964) < 1e-11);
  assert.ok(Math.abs(poissonStrictUpperTail(1, 0.1) - 0.00467884016044) < 1e-11);
});

test("FAD-compatible denoising preserves multiplicity and partitions by V/J calls", () => {
  const parent = "ACGT".repeat(15);
  const child = `${parent.slice(0, 29)}A${parent.slice(30)}`;
  const records = [record(0, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", parent), record(1, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", child)];
  records[0].inputCount = 100;
  const accumulator = new DenoiseAccumulator(records, {
    mode: "fad",
    errorRate: 0.00473,
    alpha: 0.01,
    callResolution: "allele",
    ambiguity: "strict",
    minimumParentCount: 2,
    ambiguousPolicy: "exclude",
    fadNeighborThreshold: 1,
    fadMethod: 2,
    expectedZeroErrorFraction: 1,
    maximumHammingDistance: 1,
    maximumEditDistance: 2,
    minimumIndelParentRatio: 2,
    maxCandidatesPerVariant: 10_000,
  });
  records.forEach((value) => accumulator.add(value.ordinal, value.ordinal ? child : parent));
  const result = accumulator.finish();
  assert.equal(result.mode, "fad");
  assert.equal(result.partitions, 1);
  assert.equal(result.uniqueRecords, 1);
  assert.deepEqual([...result.counts], [101, 0]);
  assert.deepEqual([...result.representatives], [0, 0]);
});

test("denoising reports variant and finalization progress and respects C-gene partitions", () => {
  const parent = "ACGT".repeat(15);
  const child = `${parent.slice(0, 29)}A${parent.slice(30)}`;
  const records = [record(0, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", parent), record(1, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", child)];
  records[0].inputCount = 100;
  records[0].cCall = "IGHM*01";
  records[1].cCall = "IGHG1*01";
  const options = {
    mode: "indel" as const, errorRate: 0.00473, alpha: 0.01, callResolution: "allele" as const,
    ambiguity: "strict" as const, minimumParentCount: 2, ambiguousPolicy: "retain" as const,
    fadNeighborThreshold: 1, fadMethod: 2 as const, expectedZeroErrorFraction: 1,
    maximumHammingDistance: 1, maximumEditDistance: 2, minimumIndelParentRatio: 2,
    maxCandidatesPerVariant: 10_000,
  };
  const separated = new DenoiseAccumulator(records, options);
  records.forEach((value) => separated.add(value.ordinal, value.ordinal ? child : parent));
  const phases = new Set<string>();
  const separatedResult = separated.finish((_processed, _total, phase) => phases.add(phase));
  assert.equal(separatedResult.uniqueRecords, 2);
  assert.deepEqual([...phases].sort(), ["finalize", "variants"]);
  const combined = new DenoiseAccumulator(records, { ...options, respectConstantCall: false });
  records.forEach((value) => combined.add(value.ordinal, value.ordinal ? child : parent));
  assert.equal(combined.finish().uniqueRecords, 1);
});

test("conservative denoising merges plausible one-base errors but retains isolated singletons", () => {
  const parent = "ACGT".repeat(15);
  const child = `${parent.slice(0, 29)}A${parent.slice(30)}`;
  const isolated = "T".repeat(parent.length);
  const records = [
    record(0, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", parent),
    record(1, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", child),
    record(2, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", isolated),
  ];
  records[0].inputCount = 100;
  const accumulator = new DenoiseAccumulator(records, {
    mode: "conservative",
    errorRate: 0.00473,
    alpha: 0.01,
    callResolution: "allele",
    ambiguity: "strict",
    minimumParentCount: 2,
    ambiguousPolicy: "retain",
    fadNeighborThreshold: 1,
    fadMethod: 2,
    expectedZeroErrorFraction: 1,
    maximumHammingDistance: 1,
    maximumEditDistance: 2,
    minimumIndelParentRatio: 2,
    maxCandidatesPerVariant: 10_000,
  });
  [parent, child, isolated].forEach((sequence, ordinal) => accumulator.add(ordinal, sequence));
  const result = accumulator.finish();
  assert.equal(result.mode, "conservative");
  assert.equal(result.uniqueRecords, 2);
  assert.deepEqual([...result.counts], [101, 0, 1]);
  assert.deepEqual([...result.representatives], [0, 0, 2]);
});

test("denoising discards unresolved records by default and can retain them explicitly", () => {
  const usable = record(0, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", "ACGT".repeat(12));
  const unresolved = record(1, "TGTGCCAAA", "", "IGHJ4*02", "TGCA".repeat(12));
  const options = {
    mode: "conservative" as const,
    errorRate: 0.00473,
    alpha: 0.01,
    callResolution: "allele" as const,
    ambiguity: "strict" as const,
    minimumParentCount: 2,
    ambiguousPolicy: "retain" as const,
    fadNeighborThreshold: 1,
    fadMethod: 2 as const,
    expectedZeroErrorFraction: 1,
    maximumHammingDistance: 1,
    maximumEditDistance: 2,
    minimumIndelParentRatio: 2,
    maxCandidatesPerVariant: 10_000,
  };
  const discarded = new DenoiseAccumulator([usable, unresolved], options);
  discarded.add(0, "ACGT".repeat(12));
  discarded.add(1, "TGCA".repeat(12));
  const discardedResult = discarded.finish();
  assert.equal(discardedResult.unresolvedRecords, 1);
  assert.equal(discardedResult.uniqueRecords, 1);
  assert.deepEqual([...discardedResult.counts], [1, 0]);
  assert.equal(discardedResult.representatives[1], -1);

  const retained = new DenoiseAccumulator([usable, unresolved], { ...options, unresolvedPolicy: "retain" });
  retained.add(0, "ACGT".repeat(12));
  retained.add(1, "TGCA".repeat(12));
  const retainedResult = retained.finish();
  assert.equal(retainedResult.uniqueRecords, 2);
  assert.deepEqual([...retainedResult.counts], [1, 1]);
  assert.equal(retainedResult.representatives[1], 1);
});

test("bounded edit profiling distinguishes substitutions, insertions, deletions and mixed paths", () => {
  assert.deepEqual(boundedEditProfile("ACGT", "AGGT", 1), { distance: 1, substitutions: 1, insertions: 0, deletions: 0 });
  assert.deepEqual(boundedEditProfile("ACGT", "ACGGT", 1), { distance: 1, substitutions: 0, insertions: 1, deletions: 0 });
  assert.deepEqual(boundedEditProfile("ACGGT", "ACGT", 1), { distance: 1, substitutions: 0, insertions: 0, deletions: 1 });
  assert.deepEqual(boundedEditProfile("ACGT", "ATGGT", 2), { distance: 2, substitutions: 1, insertions: 1, deletions: 0 });
  assert.equal(boundedEditProfile("ACGT", "ATGGGT", 1), null);
});

test("method D completely finds one-base insertions and deletions at every sequence boundary", () => {
  const parent = "ACGTCAGTGCATGACCTGATCGTACGATGCCTAGTCGATGCTACGTAGCTGACGTAGCATGCTAGCATCGATGC";
  const sequences = [parent];
  for (let position = 0; position <= parent.length; position += 1) {
    const inserted = "ACGT"[(position + 1) % 4];
    sequences.push(`${parent.slice(0, position)}${inserted}${parent.slice(position)}`);
  }
  for (let position = 0; position < parent.length; position += 1) sequences.push(`${parent.slice(0, position)}${parent.slice(position + 1)}`);
  const records = sequences.map((sequence, ordinal) => record(ordinal, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", sequence));
  records[0].inputCount = 1_000;
  const accumulator = new DenoiseAccumulator(records, {
    mode: "indel",
    errorRate: 0.00473,
    alpha: 0.01,
    callResolution: "allele",
    ambiguity: "strict",
    minimumParentCount: 2,
    ambiguousPolicy: "retain",
    fadNeighborThreshold: 1,
    fadMethod: 2,
    expectedZeroErrorFraction: 1,
    maximumHammingDistance: 1,
    maximumEditDistance: 1,
    minimumIndelParentRatio: 2,
    maxCandidatesPerVariant: 10_000,
  });
  sequences.forEach((sequence, ordinal) => accumulator.add(ordinal, sequence));
  const result = accumulator.finish();
  assert.equal(result.mode, "indel");
  assert.equal(result.uniqueRecords, 1);
  assert.equal(result.counts[0], 1_000 + sequences.length - 1);
  assert.ok(result.indelMergedVariants > 0);
  assert.equal(result.substitutionMergedVariants, 0);
});

test("method D verifies two-edit paths and preserves similarly abundant indel variants", () => {
  const parent = "ACGT".repeat(20);
  const mixed = `${parent.slice(0, 23)}A${parent.slice(23, 48)}${parent[48] === "A" ? "C" : "A"}${parent.slice(49)}`;
  const abundantIndel = `${parent.slice(0, 12)}G${parent.slice(12)}`;
  const otherPartition = `${parent.slice(0, 30)}T${parent.slice(30)}`;
  const records = [
    record(0, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", parent),
    record(1, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", mixed),
    record(2, "TGTGCCAAA", "IGHV1-2*01", "IGHJ4*02", abundantIndel),
    record(3, "TGTGCCAAA", "IGHV1-2*01", "IGHJ6*02", otherPartition),
  ];
  records[0].inputCount = 100;
  records[2].inputCount = 60;
  const accumulator = new DenoiseAccumulator(records, {
    mode: "indel",
    errorRate: 0.00473,
    alpha: 0.01,
    callResolution: "allele",
    ambiguity: "strict",
    minimumParentCount: 2,
    ambiguousPolicy: "retain",
    fadNeighborThreshold: 1,
    fadMethod: 2,
    expectedZeroErrorFraction: 1,
    maximumHammingDistance: 1,
    maximumEditDistance: 2,
    minimumIndelParentRatio: 2,
    maxCandidatesPerVariant: 10_000,
  });
  [parent, mixed, abundantIndel, otherPartition].forEach((sequence, ordinal) => accumulator.add(ordinal, sequence));
  const result = accumulator.finish();
  assert.equal(result.uniqueRecords, 3);
  assert.equal(result.representatives[1], 0);
  assert.equal(result.representatives[2], 2);
  assert.equal(result.representatives[3], 3);
  assert.equal(result.indelMergedVariants, 1);
});

test("deduplicated active representatives retain collapsed abundance in lineage summaries", () => {
  const records = [record(0, "TGTGCCAAA"), record(1, "TGTGCCAAA"), record(2, "TGTGCCTAA")];
  const collapsed = deduplicate(records, "sequence");
  const active = Uint8Array.from([1, 0, 0]);
  const result = assignLineages(records, {
    identity: 0.9,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
  }, collapsed, active);
  assert.equal(result.lineageCount, 1);
  assert.equal(result.assignedRecords, 2);
  assert.equal(result.summaries[0].uniqueMembers, 1);
  assert.equal(result.summaries[0].abundance, 2);
  assert.equal(result.summaries[0].representativeCdr3Nt, "TGTGCCAAA");
  assert.equal(result.summaries[0].representativeCdr3Aa, "CARDR");
});

test("lineages use V/J ambiguity overlap, exact CDR3 length, and bounded Hamming distance", () => {
  const records = [
    record(0, "TGTGCCAAAA", "IGHV1-2*01,IGHV1-3*01"),
    record(1, "TGTGCCAAAT", "IGHV1-3*02"),
    record(2, "TGTGCCAAATG", "IGHV1-3*02"),
    record(3, "CCCCCCCCCC", "IGHV1-3*02"),
  ];
  const result = assignLineages(records, {
    identity: 0.9,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
  });
  assert.equal(result.lineageCount, 3);
  assert.equal(result.assignments[0], result.assignments[1]);
  assert.notEqual(result.assignments[1], result.assignments[2]);
  assert.notEqual(result.assignments[1], result.assignments[3]);
  assert.equal(result.summaries[0].abundance, 2);
});

test("lineage summaries retain exact active Double-D membership for explorer filtering", () => {
  const records = [record(0, "TGTGCCAAAA"), record(1, "TGTGCCAAAT"), record(2, "TGTGCCAATA")];
  records[0].inputCount = 3;
  records[1].inputCount = 2;
  records[2].inputCount = 4;
  const result = assignLineages(records, {
    identity: 0.8,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
  }, undefined, undefined, Uint8Array.from([1, 0, 1]));
  const summary = result.summaries[0];
  assert.equal(summary.uniqueMembers, 3);
  assert.equal(summary.doubleDPositiveMembers, 2);
  assert.equal(summary.doubleDPositiveAbundance, 7);
  assert.equal(lineageDoubleDMatches(summary, "any"), true);
  assert.equal(lineageDoubleDMatches(summary, "present"), true);
  assert.equal(lineageDoubleDMatches(summary, "all"), false);
  assert.equal(lineageDoubleDMatches(summary, "absent"), false);
});

test("sequence query and single-linkage expansion recover a transitive CDR3 neighborhood", () => {
  const records = [
    record(0, "AAAAAAAAAA"),
    record(1, "AAAAAAAATA"),
    record(2, "AAAAAATATA"),
    record(3, "CCCCCCCCCC"),
  ];
  const hits = queryRecords(records, ["AAAAAAAAAA"], {
    target: "cdr3_nt",
    metric: "hamming",
    identity: 0.9,
    maxResults: 100,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
  });
  assert.deepEqual(hits.map((hit) => hit.ordinal), [0, 1]);
  const lineageHits=queryRecords(records,["AAAAAAAAAA"],{
    target:"cdr3_nt",metric:"hamming",identity:0.8,maxResults:100,callResolution:"gene",ambiguity:"overlap",productiveOnly:true,resultMode:"lineages",
  },undefined,undefined,Int32Array.from([1,1,2,3]),3);
  assert.deepEqual(lineageHits.map((hit)=>[hit.lineageId,hit.ordinal,hit.matchedSequences,hit.matchedQueries]),[[1,0,2,1],[2,2,1,1]]);
  const expanded = expandSingleLinkage(records, [0], {
    identity: 0.9,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
    maxResults: 100,
  });
  assert.deepEqual(expanded.ordinals.sort((a, b) => a - b), [0, 1, 2]);
});

test("the cumulative active mask is honored by lineage assignment, query, and expansion", () => {
  const records = [
    record(0, "AAAAAAAAAA"),
    record(1, "AAAAAAAATA"),
    record(2, "AAAAAATATA"),
    record(3, "CCCCCCCCCC"),
  ];
  const active = Uint8Array.from([1, 0, 1, 0]);
  const lineages = assignLineages(records, {
    identity: 0.9,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
  }, undefined, active);
  assert.equal(lineages.assignedRecords, 2);
  assert.equal(lineages.lineageCount, 2, "an excluded intermediate must not bridge active records");
  assert.equal(lineages.assignments[1], 0);
  assert.equal(lineages.assignments[3], 0);

  const hits = queryRecords(records, ["AAAAAAAAAA"], {
    target: "cdr3_nt",
    metric: "hamming",
    identity: 0.8,
    maxResults: 100,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
  }, undefined, active);
  assert.deepEqual(hits.map((hit) => hit.ordinal), [0, 2]);

  const expanded = expandSingleLinkage(records, [0, 1], {
    identity: 0.9,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 10_000,
    maxResults: 100,
  }, active);
  assert.deepEqual(expanded.ordinals, [0], "inactive seeds and bridge records must stay excluded");
});

test("per-query inferred V/J constraints preserve seed-specific V-J combinations", () => {
  const records = [
    record(0, "AAAAAAAAAA", "IGHV1*01", "IGHJ4*01"),
    record(1, "AAAAAAAAAA", "IGHV2*01", "IGHJ6*01"),
    record(2, "AAAAAAAAAA", "IGHV1*01", "IGHJ6*01"),
  ];
  const hits = queryRecords(records, ["AAAAAAAAAA", "AAAAAAAAAA"], {
    target: "cdr3_nt",
    metric: "exact",
    identity: 1,
    maxResults: 100,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    queryConstraints: [
      { locus: "IGH", vCall: "IGHV1*01", jCall: "IGHJ4*01" },
      { locus: "IGH", vCall: "IGHV2*01", jCall: "IGHJ6*01" },
    ],
  });
  assert.deepEqual(hits.map((hit) => [hit.ordinal, hit.queryIndex]), [[0, 0], [1, 1]]);
});

test("CHMMera BW port separates the published chimera and non-chimera example", () => {
  const msa = prepareReferenceMsa(">ref1\nACGTACGTACGT\n>ref2\nACCACCACCAAT\n");
  const options = {
    method: "BW" as const,
    priorProbability: 0.02,
    baseMutationProbability: 0.05,
    mutationRates: [0.005],
    mutationSwitchProbability: 0,
    detailed: true,
    tracePath: true,
  };
  const chimera = runChmm(msa, "ACGTACACCAAT", "ACGTACACCAAT", "ACGTACGTACGT", options);
  const nonChimera = runChmm(msa, "ACCACCACCAGT", "ACCACCACCAGT", "ACCACCACCAAT", options);
  assert.ok(chimera.probability > 0.95, `${chimera.probability}`);
  assert.ok(nonChimera.probability < 0.1, `${nonChimera.probability}`);
  assert.equal(chimera.recombinations.length, 1);
  assert.deepEqual([...chimera.referencePath!].slice(0, 4), [0, 0, 0, 0]);
  assert.deepEqual([...chimera.referencePath!].slice(-4), [1, 1, 1, 1]);
});

test("CHMMera discretized Bayesian mode reports detailed parent switches", () => {
  const msa = prepareReferenceMsa(">left\nAAAAAAAAAAAA\n>right\nCCCCCCCCCCCC\n");
  const result = runChmm(msa, "CCCCCCAAAAAA", "CCCCCCAAAAAA", "AAAAAAAAAAAA", {
    method: "DB",
    priorProbability: 0.05,
    baseMutationProbability: 0.005,
    mutationRates: [0.005],
    mutationSwitchProbability: 0,
    detailed: true,
    tracePath: true,
  });
  assert.ok(result.probability > 0.95, `${result.probability}`);
  assert.ok(result.recombinations.length >= 1);
  assert.notEqual(result.startingReference, result.recombinations[0].right);
  assert.equal(result.referencePath?.length, 12);
  assert.equal(new Set(result.referencePath).size, 2);
});

test("AIRR local alignments are threaded onto an uploaded reference MSA", () => {
  const msa = prepareReferenceMsa(">IGHV1*01\nAA--CCGG\n>IGHV2*01\nAATTCCGG\n");
  const threaded = threadSequenceToMsa("TCGG", "CCGG", "IGHV1*01", msa);
  assert.equal(threaded, "NN--TCGG");
});

test("CHMMAIRRa threading maps ambiguous observations to N and DFR counts gaps and N mismatches", () => {
  const msa = prepareReferenceMsa(">IGHV1*01\nAA--CCGG\n>IGHV2*01\nAATTCCGG\n");
  assert.equal(threadSequenceToMsa("RCGG", "CCGG", "IGHV1*01", msa), "NN--NCGG");
  assert.equal(threadSequenceToMsa("CTCGG", "C-CGG", "IGHV1*01", msa), "NN--CCGG", "query insertions opposite local germline gaps are dropped");
  assert.equal(threadSequenceToMsa("C-GG", "CCGG", "IGHV1*01", msa), "NN--C-GG", "query deletions are retained as non-informative gaps");
  assert.equal(chmmairraDistanceFromReference("AC-N", "ACAN"), 1);
  assert.equal(chmmairraDistanceFromReference("ACNN", "ACAA"), 2);
});

test("50k post-analysis records retain bounded compact state and avoid all-pairs lineage comparisons", () => {
  let state = 0x12345678;
  const randomNt = () => {
    let value = "";
    for (let index = 0; index < 36; index += 1) {
      state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
      value += "ACGT"[state >>> 30];
    }
    return value;
  };
  const records: PostAnalysisRecord[] = [];
  for (let ordinal = 0; ordinal < 50_000; ordinal += 1) {
    const duplicate = ordinal % 5 === 4;
    const cdr3 = duplicate ? records[ordinal - 1].cdr3Nt : randomNt();
    records.push(record(ordinal, cdr3, `IGHV${ordinal % 20 + 1}*01`, `IGHJ${ordinal % 6 + 1}*01`, `AAA${cdr3}TTT`));
  }
  const collapsed = deduplicate(records, "sequence");
  assert.equal(collapsed.inputRecords, 50_000);
  assert.equal(collapsed.uniqueRecords, 40_000);
  const lineages = assignLineages(records, {
    identity: 0.85,
    callResolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
    requireSameLocus: true,
    maxCandidateComparisons: 50_000,
  }, collapsed);
  assert.equal(lineages.assignedRecords, 50_000);
  assert.equal(lineages.truncatedCandidates, 0);
  assert.ok(lineages.candidateComparisons < 50_000 * 100, `${lineages.candidateComparisons}`);
  assert.ok(lineages.lineageCount <= collapsed.uniqueRecords);
});

test("the denoising V-to-J slice reads AIRR coordinates on the aligned strand", () => {
  // The aligned read is TTACGTACGTCC; V starts at 3 and J ends at 10.
  const forward = { ordinal: 0, values: { sequence: "TTACGTACGTCC", rev_comp: "F", v_sequence_start: "3", j_sequence_end: "10" } };
  const reverse = { ordinal: 1, values: { sequence: "GGACGTACGTAA", rev_comp: "T", v_sequence_start: "3", j_sequence_end: "10" } };
  assert.equal(denoiseVdjSequence(forward), "ACGTACGT");
  assert.equal(denoiseVdjSequence(reverse), "ACGTACGT");
  // IUPAC codes are complemented (R↔Y), and both strands come out upper-case.
  const iupac = { ordinal: 2, values: { sequence: "GGACGTRCGTAA", rev_comp: "1", v_sequence_start: "3", j_sequence_end: "10" } };
  assert.equal(denoiseVdjSequence(iupac), "ACGYACGT");
  assert.equal(denoiseVdjSequence({ ordinal: 3, values: { ...forward.values, sequence: "ttacgtacgtcc" } }), "ACGTACGT");
});

test("reverse complements keep IUPAC codes and accept every AIRR spelling of rev_comp", () => {
  assert.equal(reverseComplement("acgtRYKMBVDHSWNU-"), "NANWSDHBVKMRYACGT");
  for (const value of ["T", "true", "TRUE", "1", " t "]) assert.ok(isReverseComplement(value), value);
  for (const value of ["F", "false", "0", "", undefined]) assert.ok(!isReverseComplement(value), String(value));
});

test("CHMMAIRRa extends a local alignment toward the junction by its best-scoring flank", () => {
  const msa = prepareReferenceMsa(">IGHV1*01\nAA--CCGGTACA\n>IGHV2*01\nAATTCCGGTTCA\n>IGHJ1*01\nAA--CCGGTACA\n");
  // V: the aligner stopped before a trailing mismatch the next bases outweigh (X then ==).
  assert.deepEqual(extendChmmAlignment("V", "AACCGG", "AACCGG", "IGHV1*01", "AACAGGG", msa), { sequenceAlignment: "AACCGGAACA", germlineAlignment: "AACCGGTACA" });
  // A flank that never scores positive (junction bases) is not taken.
  assert.deepEqual(extendChmmAlignment("V", "AACCGG", "AACCGG", "IGHV1*01", "CGTGCC", msa), { sequenceAlignment: "AACCGG", germlineAlignment: "AACCGG" });
  // One matching base extends a V.
  assert.deepEqual(extendChmmAlignment("V", "AACCGG", "AACCGG", "IGHV1*01", "TGGG", msa), { sequenceAlignment: "AACCGGT", germlineAlignment: "AACCGGT" });
  // J: extended 5'-ward with the read bases before the alignment.
  assert.deepEqual(extendChmmAlignment("J", "GGTACA", "GGTACA", "IGHJ1*01", "TACC", msa), { sequenceAlignment: "ACCGGTACA", germlineAlignment: "ACCGGTACA" });
  // A J extension needs a score of 2: one chance-matching base is not taken.
  assert.deepEqual(extendChmmAlignment("J", "GGTACA", "GGTACA", "IGHJ1*01", "TTGC", msa), { sequenceAlignment: "GGTACA", germlineAlignment: "GGTACA" });
  assert.equal(chmmJunctionFlank("V", "TTACGTAA", "F", 1, 4), "GTAA");
  assert.equal(chmmJunctionFlank("J", "TTACGTAA", "F", 5, 8), "TTAC");
  assert.equal(chmmJunctionFlank("V", "AAACCCGG", "T", 1, 4), "GTTT", "rev_comp rows use the reverse complement of sequence");
  assert.equal(chmmJunctionFlank("V", "AAACCCGG", "1", 1, 4), "GTTT");
  assert.equal(threadSequenceToMsa("AACCGGAACA", "AACCGGTACA", "IGHV1*01", msa), "AA--CCGGAACA");
});
