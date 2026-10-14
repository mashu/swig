import { datasetScopeKey, datasetScopeValue, type DatasetScope } from "./study-design.ts";

export type CallResolution = "gene" | "allele";
export type AmbiguityPolicy = "overlap" | "top" | "strict";
export type DedupKey = "sequence" | "trimmed" | "cdr3" | "rearrangement";
export type CollapseMode = "exact" | "fad" | "conservative" | "indel";

export interface PostAnalysisRecord {
  ordinal: number;
  sequenceId: string;
  datasetId?: string;
  sampleId?: string;
  subjectId?: string;
  cohort?: string;
  timepoint?: string;
  compartment?: string;
  locus: string;
  vCall: string;
  jCall: string;
  /** Immutable AIRR calls retained when repertoire-level posterior calls are reset. */
  originalVCall?: string;
  originalJCall?: string;
  /** AIRR constant-region call. Used as a gene-level collapse partition when requested. */
  cCall?: string;
  cdr3Nt: string;
  cdr3Aa: string;
  productive: boolean;
  sequenceFingerprint: string;
  trimmedFingerprint: string;
  trimmedSketch?: Uint32Array;
  inputCount?: number;
}

export interface DedupResult {
  mode: CollapseMode;
  key: DedupKey;
  algorithm: string;
  inputRecords: number;
  inputAbundance: number;
  uniqueRecords: number;
  collapsedRecords: number;
  representatives: Int32Array;
  counts: Uint32Array;
  largestGroups: Array<{ ordinal: number; count: number }>;
  partitions: number;
  candidateComparisons: number;
  indelMergedVariants: number;
  substitutionMergedVariants: number;
  excludedAmbiguous: number;
  unresolvedRecords: number;
  warnings: string[];
}

export interface DenoiseOptions {
  mode: "fad" | "conservative" | "indel";
  errorRate: number;
  alpha: number;
  callResolution: CallResolution;
  ambiguity: "top" | "strict";
  minimumParentCount: number;
  ambiguousPolicy: "exclude" | "retain";
  /** Records that cannot enter a V/J-partitioned denoising model. */
  unresolvedPolicy?: "discard" | "retain";
  /** FAD corrected 6-mer distance radius. */
  fadNeighborThreshold: number;
  /** FAD method 1 is abundance-only; method 2 uses its Poisson decision. */
  fadMethod: 1 | 2;
  expectedZeroErrorFraction: number;
  /** Exact Hamming radius for the conservative error model. */
  maximumHammingDistance: number;
  /** Complete bounded Levenshtein radius for indel-aware method D. */
  maximumEditDistance: number;
  /** Required abundance ratio for an indel-containing child to collapse. */
  minimumIndelParentRatio: number;
  maxCandidatesPerVariant: number;
  /** Hard boundary for candidate generation; defaults to global for legacy API calls. */
  scope?: DatasetScope;
  /** Keep records with different C-gene/isotype assignments in separate partitions. */
  respectConstantCall?: boolean;
}

export interface LineageOptions {
  identity: number;
  callResolution: CallResolution;
  ambiguity: AmbiguityPolicy;
  productiveOnly: boolean;
  requireSameLocus: boolean;
  maxCandidateComparisons: number;
  /** Lineages cannot cross this study boundary. */
  scope?: DatasetScope;
}

export interface LineageSummary {
  id: number;
  representativeOrdinal: number;
  uniqueMembers: number;
  abundance: number;
  locus: string;
  vCalls: string[];
  jCalls: string[];
  cdr3Length: number;
  /** Representative CDR3 available without loading lineage members. */
  representativeCdr3Nt?: string;
  representativeCdr3Aa?: string;
  studyScope: DatasetScope;
  studyGroup: string;
  /** Study metadata represented among active members of this lineage summary. */
  sampleIds: string[];
  /** Per-sample unique representatives and duplicate_count-weighted read abundance. */
  sampleCounts: Array<{ sampleId: string; uniqueMembers: number; abundance: number }>;
  subjectIds: string[];
  timepoints: string[];
  compartments: string[];
  /** Active unique lineage representatives with supported VDDJ evidence. */
  doubleDPositiveMembers?: number;
  /** Multiplicity-weighted abundance represented by those positive members. */
  doubleDPositiveAbundance?: number;
}

export type LineageDoubleDFilter = "any" | "present" | "all" | "absent";

export function lineageDoubleDMatches(summary: LineageSummary, mode: LineageDoubleDFilter): boolean {
  const positives = Math.max(0, summary.doubleDPositiveMembers ?? 0);
  if (mode === "present") return positives > 0;
  if (mode === "all") return positives > 0 && positives === summary.uniqueMembers;
  if (mode === "absent") return positives === 0;
  return true;
}

export interface LineageResult {
  assignments: Int32Array;
  summaries: LineageSummary[];
  lineageCount: number;
  sizeHistogram: Array<{ label: string; count: number }>;
  vUsage: Array<{ call: string; lineages: number; abundance: number }>;
  jUsage: Array<{ call: string; lineages: number; abundance: number }>;
  assignedRecords: number;
  unassignedRecords: number;
  candidateComparisons: number;
  truncatedCandidates: number;
}

export interface LineageNeighbourOptions extends LineageOptions {
  sourceLineageIds: number[];
  /** Lower, exploratory CDR3 identity boundary; normally below lineage assignment. */
  minimumIdentity: number;
  maximumResults: number;
  /** Equal-length Hamming preserves the assignment metric; edit admits indels. */
  metric?: "hamming" | "edit";
  /** How much potentially noisy V/J evidence is allowed to gate a CDR3 hit. */
  callPolicy?: "both" | "either" | "ignore";
}

export interface LineageNeighbourHit {
  lineageId: number;
  sourceLineageId: number;
  sourceOrdinal: number;
  candidateOrdinal: number;
  cdr3Identity: number;
  cdr3Distance?: number;
  cdr3Metric?: "hamming" | "edit";
  cdr3LengthDelta?: number;
  callAgreement?: "both" | "v" | "j" | "none";
  uniqueMembers: number;
  abundance: number;
  locus: string;
  vCalls: string[];
  jCalls: string[];
  cdr3Length: number;
  studyGroup: string;
}

export interface LineageNeighbourResult {
  hits: LineageNeighbourHit[];
  indexedRecords: number;
  sourceRecords: number;
  candidateComparisons: number;
  truncatedSourceRecords: number;
}

export type QueryTarget = "cdr3_nt" | "cdr3_aa" | "trimmed";
export type QueryMetric = "exact" | "substring" | "hamming" | "edit" | "sketch";

export interface QueryConstraint {
  locus?: string;
  vCall?: string;
  jCall?: string;
}

export interface QueryOptions {
  target: QueryTarget;
  metric: QueryMetric;
  identity: number;
  maxResults: number;
  locus?: string;
  vCall?: string;
  jCall?: string;
  callResolution: CallResolution;
  ambiguity: AmbiguityPolicy;
  productiveOnly: boolean;
  queryConstraints?: QueryConstraint[];
  /** Rank individual sequence hits or one exact best-member hit per lineage. */
  resultMode?: "sequences" | "lineages";
}

export interface QueryHit {
  ordinal: number;
  queryIndex: number;
  score: number;
  distance: number;
  matched: string;
  lineageId?: number;
  matchedSequences?: number;
  matchedQueries?: number;
}

const HASH_SEEDS = [
  0x811c9dc5, 0x9e3779b9, 0x85ebca6b, 0xc2b2ae35,
  0x27d4eb2f, 0x165667b1, 0xd3a2646c, 0xfd7046c5,
] as const;

function popcount32(value: number): number {
  value -= (value >>> 1) & 0x55555555;
  value = (value & 0x33333333) + ((value >>> 2) & 0x33333333);
  return (((value + (value >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
}

function mix32(value: number): number {
  value ^= value >>> 16;
  value = Math.imul(value, 0x7feb352d);
  value ^= value >>> 15;
  value = Math.imul(value, 0x846ca68b);
  value ^= value >>> 16;
  return value >>> 0;
}

export function hashSequence(value: string, seed: number = HASH_SEEDS[0]): number {
  let hash = seed >>> 0;
  for (let index = 0; index < value.length; index += 1) {
    hash ^= value.charCodeAt(index);
    hash = Math.imul(hash, 0x01000193);
  }
  return mix32(hash);
}

export function sequenceFingerprint(value: string): string {
  const normalized = normalizeNt(value);
  return `${normalized.length}:${[0, 2, 4, 6].map((index) => hashSequence(normalized, HASH_SEEDS[index]).toString(36)).join(":")}`;
}

export function minHashSketch(value: string, k = 7): Uint32Array {
  const normalized = normalizeNt(value).replaceAll("N", "");
  const sketch = new Uint32Array(HASH_SEEDS.length);
  sketch.fill(0xffffffff);
  const width = Math.max(1, Math.min(k, normalized.length));
  if (!normalized.length) return sketch;
  for (let index = 0; index <= normalized.length - width; index += 1) {
    const kmer = normalized.slice(index, index + width);
    for (let seed = 0; seed < HASH_SEEDS.length; seed += 1) {
      const hash = hashSequence(kmer, HASH_SEEDS[seed]);
      if (hash < sketch[seed]) sketch[seed] = hash;
    }
  }
  return sketch;
}

export function sketchSimilarity(left: Uint32Array, right: Uint32Array): number {
  const length = Math.min(left.length, right.length);
  if (!length) return 0;
  let equal = 0;
  let informative = 0;
  for (let index = 0; index < length; index += 1) {
    if (left[index] === 0xffffffff && right[index] === 0xffffffff) continue;
    informative += 1;
    if (left[index] === right[index]) equal += 1;
  }
  return informative ? equal / informative : 0;
}

export function normalizeNt(value: string): string {
  return value.toUpperCase().replaceAll("U", "T").replace(/[^ACGTN]/g, "");
}

export function normalizeAa(value: string): string {
  return value.toUpperCase().replace(/[^A-Z*]/g, "");
}

export function normalizeCall(call: string, resolution: CallResolution): string {
  const clean = call.trim().toUpperCase();
  return resolution === "gene" ? clean.replace(/\*.*$/, "") : clean;
}

export function callSet(value: string, resolution: CallResolution, ambiguity: AmbiguityPolicy): string[] {
  const values = value.split(",").map((call) => normalizeCall(call, resolution)).filter(Boolean);
  const selected = ambiguity === "top" ? values.slice(0, 1) : [...new Set(values)].sort();
  return selected;
}

function callsCompatible(
  left: string,
  right: string,
  resolution: CallResolution,
  ambiguity: AmbiguityPolicy,
): boolean {
  const a = callSet(left, resolution, ambiguity);
  const b = callSet(right, resolution, ambiguity);
  if (!a.length || !b.length) return false;
  if (ambiguity === "strict") return a.length === b.length && a.every((value, index) => value === b[index]);
  const values = new Set(a);
  return b.some((value) => values.has(value));
}

export function hammingDistanceWithin(left: string, right: string, maximum: number, ambiguousN = true): number {
  if (left.length !== right.length) return maximum + 1;
  let distance = 0;
  for (let index = 0; index < left.length; index += 1) {
    if (left[index] !== right[index] && (!ambiguousN || (left[index] !== "N" && right[index] !== "N"))) {
      distance += 1;
      if (distance > maximum) return distance;
    }
  }
  return distance;
}

export function bandedEditDistance(left: string, right: string, maximum: number): number {
  if (Math.abs(left.length - right.length) > maximum) return maximum + 1;
  let previous = new Int32Array(right.length + 1);
  let current = new Int32Array(right.length + 1);
  for (let column = 0; column <= right.length; column += 1) previous[column] = column;
  for (let row = 1; row <= left.length; row += 1) {
    current.fill(maximum + 1);
    current[0] = row;
    const start = Math.max(1, row - maximum);
    const end = Math.min(right.length, row + maximum);
    let rowMinimum = maximum + 1;
    for (let column = start; column <= end; column += 1) {
      const substitution = left[row - 1] === right[column - 1] ? 0 : 1;
      current[column] = Math.min(
        previous[column] + 1,
        current[column - 1] + 1,
        previous[column - 1] + substitution,
      );
      rowMinimum = Math.min(rowMinimum, current[column]);
    }
    if (rowMinimum > maximum) return maximum + 1;
    [previous, current] = [current, previous];
  }
  return previous[right.length];
}

function constantCallPartition(record: PostAnalysisRecord, enabled: boolean): string {
  if (!enabled) return "";
  const topGene = callSet(record.cCall ?? "", "gene", "top")[0] ?? "__C_UNASSIGNED__";
  return `\u0000C:${topGene}`;
}

function dedupKey(record: PostAnalysisRecord, key: DedupKey, scope: DatasetScope, respectConstantCall: boolean): string {
  const prefix = `${datasetScopeKey(record, scope)}${constantCallPartition(record, respectConstantCall)}\u0000`;
  if (key === "sequence") return `${prefix}${record.sequenceFingerprint}`;
  if (key === "trimmed") return `${prefix}${record.trimmedFingerprint}`;
  if (key === "cdr3") return `${prefix}${record.locus}\u0000${record.cdr3Nt}`;
  return `${prefix}${record.locus}\u0000${record.vCall}\u0000${record.jCall}\u0000${record.cdr3Nt}`;
}

function largestCountGroups(counts: Uint32Array, limit = 100): Array<{ ordinal: number; count: number }> {
  const heap: number[] = [];
  const worse = (left: number, right: number) => counts[left] < counts[right] || (counts[left] === counts[right] && left > right);
  const rise = (start: number) => {
    let index = start;
    while (index > 0) {
      const parent = Math.floor((index - 1) / 2);
      if (!worse(heap[index], heap[parent])) break;
      [heap[index], heap[parent]] = [heap[parent], heap[index]];
      index = parent;
    }
  };
  const sink = () => {
    let index = 0;
    while (true) {
      const left = index * 2 + 1;
      const right = left + 1;
      let next = index;
      if (left < heap.length && worse(heap[left], heap[next])) next = left;
      if (right < heap.length && worse(heap[right], heap[next])) next = right;
      if (next === index) break;
      [heap[index], heap[next]] = [heap[next], heap[index]];
      index = next;
    }
  };
  for (let ordinal = 0; ordinal < counts.length; ordinal += 1) {
    if (counts[ordinal] <= 1) continue;
    if (heap.length < limit) {
      heap.push(ordinal);
      rise(heap.length - 1);
    } else if (worse(heap[0], ordinal)) {
      heap[0] = ordinal;
      sink();
    }
  }
  return heap.sort((left, right) => counts[right] - counts[left] || left - right).map((ordinal) => ({ ordinal, count: counts[ordinal] }));
}

function hasUsableDedupKey(record: PostAnalysisRecord, key: DedupKey): boolean {
  if (key === "sequence") return Number(record.sequenceFingerprint.split(":", 1)[0]) > 0;
  if (key === "trimmed") return Number(record.trimmedFingerprint.split(":", 1)[0]) > 0;
  if (key === "cdr3") return Boolean(record.locus && record.cdr3Nt);
  return Boolean(record.locus && record.vCall && record.jCall && record.cdr3Nt);
}

export function deduplicate(
  records: PostAnalysisRecord[],
  key: DedupKey,
  unresolvedPolicy: "discard" | "retain" = "discard",
  scope: DatasetScope = "global",
  respectConstantCall = true,
): DedupResult {
  const representatives = new Int32Array(records.length);
  representatives.fill(-1);
  const counts = new Uint32Array(records.length);
  const seen = new Map<string, number>();
  let inputAbundance = 0;
  let unresolvedRecords = 0;
  let retainedUnresolved = 0;
  for (let index = 0; index < records.length; index += 1) {
    const weight = Math.max(1, Math.floor(records[index].inputCount ?? 1));
    inputAbundance += weight;
    if (!hasUsableDedupKey(records[index], key)) {
      unresolvedRecords += 1;
      if (unresolvedPolicy === "retain") {
        representatives[index] = index;
        counts[index] = weight;
        retainedUnresolved += 1;
      }
      continue;
    }
    const value = dedupKey(records[index], key, scope, respectConstantCall);
    const previous = seen.get(value);
    if (previous === undefined) {
      seen.set(value, index);
      representatives[index] = index;
      counts[index] = weight;
    } else {
      representatives[index] = previous;
      counts[previous] += weight;
    }
  }
  const largestGroups = largestCountGroups(counts);
  const uniqueRecords = seen.size + retainedUnresolved;
  return {
    mode: "exact",
    key,
    algorithm: "Exact key collapse",
    inputRecords: records.length,
    inputAbundance,
    uniqueRecords,
    collapsedRecords: records.length - uniqueRecords,
    representatives,
    counts,
    largestGroups,
    partitions: 1,
    candidateComparisons: 0,
    indelMergedVariants: 0,
    substitutionMergedVariants: 0,
    excludedAmbiguous: 0,
    unresolvedRecords,
    warnings: unresolvedRecords ? [
      `${unresolvedRecords.toLocaleString()} records without a usable ${key} key were ${unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`,
    ] : [],
  };
}

/** A deterministic shard of an exact-collapse job. Identical keys are always
 * routed to the same shard, while records inside each shard remain in AIRR
 * order so the serial implementation's first-record representative is kept. */
export interface ExactDedupJob {
  id: number;
  ordinals: Int32Array;
  weights: Uint32Array;
  keys: string[];
}

export interface ExactDedupJobResult {
  id: number;
  representatives: Int32Array;
  representativeOrdinals: Int32Array;
  representativeCounts: Uint32Array;
  uniqueRecords: number;
}

export interface ExactDedupPlan {
  key: DedupKey;
  recordCount: number;
  inputAbundance: number;
  unresolvedPolicy: "discard" | "retain";
  unresolvedOrdinals: Int32Array;
  unresolvedWeights: Uint32Array;
  jobs: ExactDedupJob[];
}

export function createExactDedupPlan(
  records: PostAnalysisRecord[],
  key: DedupKey,
  unresolvedPolicy: "discard" | "retain" = "discard",
  scope: DatasetScope = "global",
  respectConstantCall = true,
  requestedShards = 1,
): ExactDedupPlan {
  const shardCount = Math.max(1, Math.min(Math.floor(requestedShards) || 1, records.length || 1));
  const shardOrdinals = Array.from({ length: shardCount }, () => [] as number[]);
  const shardWeights = Array.from({ length: shardCount }, () => [] as number[]);
  const shardKeys = Array.from({ length: shardCount }, () => [] as string[]);
  const unresolvedOrdinals: number[] = [];
  const unresolvedWeights: number[] = [];
  let inputAbundance = 0;
  for (let ordinal = 0; ordinal < records.length; ordinal += 1) {
    const weight = Math.max(1, Math.floor(records[ordinal].inputCount ?? 1));
    inputAbundance += weight;
    if (!hasUsableDedupKey(records[ordinal], key)) {
      unresolvedOrdinals.push(ordinal);
      unresolvedWeights.push(weight);
      continue;
    }
    const value = dedupKey(records[ordinal], key, scope, respectConstantCall);
    const shard = hashSequence(value) % shardCount;
    shardOrdinals[shard].push(ordinal);
    shardWeights[shard].push(weight);
    shardKeys[shard].push(value);
  }
  const jobs: ExactDedupJob[] = [];
  for (let shard = 0; shard < shardCount; shard += 1) {
    if (!shardOrdinals[shard].length) continue;
    jobs.push({
      id: shard,
      ordinals: Int32Array.from(shardOrdinals[shard]),
      weights: Uint32Array.from(shardWeights[shard]),
      keys: shardKeys[shard],
    });
  }
  return {
    key,
    recordCount: records.length,
    inputAbundance,
    unresolvedPolicy,
    unresolvedOrdinals: Int32Array.from(unresolvedOrdinals),
    unresolvedWeights: Uint32Array.from(unresolvedWeights),
    jobs,
  };
}

export function runExactDedupJob(job: ExactDedupJob): ExactDedupJobResult {
  if (job.ordinals.length !== job.weights.length || job.ordinals.length !== job.keys.length) {
    throw new Error("The exact-collapse worker received inconsistent shard vectors.");
  }
  const representatives = new Int32Array(job.ordinals.length);
  const seen = new Map<string, number>();
  const countByRepresentative = new Map<number, number>();
  for (let index = 0; index < job.ordinals.length; index += 1) {
    const previous = seen.get(job.keys[index]);
    const representative = previous ?? job.ordinals[index];
    if (previous === undefined) seen.set(job.keys[index], representative);
    representatives[index] = representative;
    countByRepresentative.set(representative, (countByRepresentative.get(representative) ?? 0) + job.weights[index]);
  }
  const representativeOrdinals = Int32Array.from(countByRepresentative.keys());
  const representativeCounts = Uint32Array.from(representativeOrdinals, (ordinal) => countByRepresentative.get(ordinal) ?? 0);
  return { id: job.id, representatives, representativeOrdinals, representativeCounts, uniqueRecords: seen.size };
}

export function finishExactDedupPlan(plan: ExactDedupPlan, results: ExactDedupJobResult[]): DedupResult {
  const representatives = new Int32Array(plan.recordCount);
  representatives.fill(-1);
  const counts = new Uint32Array(plan.recordCount);
  const resultById = new Map(results.map((result) => [result.id, result] as const));
  let uniqueRecords = 0;
  for (const job of plan.jobs) {
    const result = resultById.get(job.id);
    if (!result || result.representatives.length !== job.ordinals.length) throw new Error(`Exact-collapse shard ${job.id} did not return a complete result.`);
    uniqueRecords += result.uniqueRecords;
    for (let index = 0; index < job.ordinals.length; index += 1) representatives[job.ordinals[index]] = result.representatives[index];
    for (let index = 0; index < result.representativeOrdinals.length; index += 1) counts[result.representativeOrdinals[index]] = result.representativeCounts[index];
  }
  if (plan.unresolvedPolicy === "retain") {
    uniqueRecords += plan.unresolvedOrdinals.length;
    for (let index = 0; index < plan.unresolvedOrdinals.length; index += 1) {
      const ordinal = plan.unresolvedOrdinals[index];
      representatives[ordinal] = ordinal;
      counts[ordinal] = plan.unresolvedWeights[index];
    }
  }
  const unresolvedRecords = plan.unresolvedOrdinals.length;
  return {
    mode: "exact",
    key: plan.key,
    algorithm: "Exact key collapse",
    inputRecords: plan.recordCount,
    inputAbundance: plan.inputAbundance,
    uniqueRecords,
    collapsedRecords: plan.recordCount - uniqueRecords,
    representatives,
    counts,
    largestGroups: largestCountGroups(counts),
    partitions: 1,
    candidateComparisons: 0,
    indelMergedVariants: 0,
    substitutionMergedVariants: 0,
    excludedAmbiguous: 0,
    unresolvedRecords,
    warnings: unresolvedRecords ? [
      `${unresolvedRecords.toLocaleString()} records without a usable ${plan.key} key were ${plan.unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`,
    ] : [],
  };
}

interface PackedDnaLocation {
  chunk: number;
  offset: number;
  length: number;
}

/**
 * Append-only two-bit DNA storage. Denoising may need the complete trimmed VDJ
 * sequence, but retaining one JavaScript string per unique read is prohibitive
 * for large repertoires. Fixed chunks also avoid a transient 2x allocation
 * when an expanding typed array is copied.
 */
class PackedDnaArena {
  private readonly chunkWords = 262_144;
  private readonly chunks: Uint32Array[] = [];
  private readonly used: number[] = [];

  append(sequence: string): PackedDnaLocation {
    const words = Math.ceil(sequence.length / 16);
    let chunk = this.chunks.length - 1;
    if (chunk < 0 || this.used[chunk] + words > this.chunkWords) {
      chunk = this.chunks.length;
      this.chunks.push(new Uint32Array(Math.max(this.chunkWords, words)));
      this.used.push(0);
    }
    const offset = this.used[chunk];
    this.used[chunk] += words;
    const target = this.chunks[chunk];
    for (let index = 0; index < sequence.length; index += 1) {
      const code = sequence[index] === "C" ? 1 : sequence[index] === "G" ? 2 : sequence[index] === "T" ? 3 : 0;
      target[offset + (index >>> 4)] |= code << ((index & 15) * 2);
    }
    return { chunk, offset, length: sequence.length };
  }

  base(location: PackedDnaLocation, index: number): string {
    const code = (this.chunks[location.chunk][location.offset + (index >>> 4)] >>> ((index & 15) * 2)) & 3;
    return code === 0 ? "A" : code === 1 ? "C" : code === 2 ? "G" : "T";
  }

  equals(location: PackedDnaLocation, sequence: string): boolean {
    if (location.length !== sequence.length) return false;
    for (let index = 0; index < sequence.length; index += 1) if (this.base(location, index) !== sequence[index]) return false;
    return true;
  }

  decode(location: PackedDnaLocation): string {
    let sequence = "";
    // Joining moderate blocks avoids quadratic string growth for long reads.
    for (let start = 0; start < location.length; start += 1_024) {
      const values: string[] = [];
      const end = Math.min(location.length, start + 1_024);
      for (let index = start; index < end; index += 1) values.push(this.base(location, index));
      sequence += values.join("");
    }
    return sequence;
  }
}

interface DenoiseVariant {
  location: PackedDnaLocation;
  partition: string;
  representative: number;
  count: number;
  target: number;
}

interface DenoisePartitionStats {
  comparisons: number;
  truncated: number;
  indelMergedVariants: number;
  substitutionMergedVariants: number;
}

interface KmerProfile {
  codes: Uint16Array;
  counts: Uint16Array;
  hashes: Uint32Array;
}

function denoisePartition(record: PostAnalysisRecord, options: DenoiseOptions): string | null {
  const v = callSet(record.vCall, options.callResolution, options.ambiguity);
  const j = callSet(record.jCall, options.callResolution, options.ambiguity);
  if (!v.length || !j.length) return null;
  return `${datasetScopeKey(record, options.scope ?? "global")}${constantCallPartition(record, options.respectConstantCall ?? true)}\u0000${record.locus || "?"}\u0000${v.join("+")}\u0000${j.join("+")}`;
}

function kmerProfile(sequence: string, blockCount: number, k = 6): KmerProfile {
  const space = 1 << (2 * k);
  const counts = new Uint16Array(space);
  const touched: number[] = [];
  const mask = space - 1;
  let code = 0;
  for (let index = 0; index < sequence.length; index += 1) {
    const base = sequence[index] === "C" ? 1 : sequence[index] === "G" ? 2 : sequence[index] === "T" ? 3 : 0;
    code = ((code << 2) | base) & mask;
    if (index < k - 1) continue;
    if (!counts[code]) touched.push(code);
    if (counts[code] < 0xffff) counts[code] += 1;
  }
  touched.sort((left, right) => left - right);
  const codes = Uint16Array.from(touched);
  const sparseCounts = Uint16Array.from(touched, (value) => counts[value]);
  const hashes = new Uint32Array(blockCount);
  for (let block = 0; block < blockCount; block += 1) hashes[block] = mix32(0x9e3779b9 ^ block);
  for (let index = 0; index < codes.length; index += 1) {
    const block = codes[index] % blockCount;
    hashes[block] = mix32(hashes[block] ^ mix32(Math.imul(codes[index] + 1, 0x85ebca6b) ^ sparseCounts[index]));
  }
  return { codes, counts: sparseCounts, hashes };
}

function kmerSquaredDistance(left: KmerProfile, right: KmerProfile, ceiling = Number.POSITIVE_INFINITY): number {
  let a = 0;
  let b = 0;
  let distance = 0;
  while (a < left.codes.length || b < right.codes.length) {
    if (b >= right.codes.length || (a < left.codes.length && left.codes[a] < right.codes[b])) {
      distance += left.counts[a] ** 2;
      a += 1;
    } else if (a >= left.codes.length || right.codes[b] < left.codes[a]) {
      distance += right.counts[b] ** 2;
      b += 1;
    } else {
      distance += (left.counts[a] - right.counts[b]) ** 2;
      a += 1;
      b += 1;
    }
    if (distance > ceiling) return distance;
  }
  return distance;
}

interface KmerVpNode {
  point: number;
  radius: number;
  inside: KmerVpNode | null;
  outside: KmerVpNode | null;
}

/** Exact metric index for FAD's final nearest-centroid assignment. */
function buildKmerVpTree(points: number[], profiles: Map<number, KmerProfile>, compared: () => void): KmerVpNode | null {
  if (!points.length) return null;
  const point = points[points.length - 1];
  if (points.length === 1) return { point, radius: 0, inside: null, outside: null };
  const distances = points.slice(0, -1).map((candidate) => {
    compared();
    return { candidate, distance: Math.sqrt(kmerSquaredDistance(profiles.get(point)!, profiles.get(candidate)!)) };
  }).sort((left, right) => left.distance - right.distance || left.candidate - right.candidate);
  const middle = Math.floor(distances.length / 2);
  const radius = distances[middle]?.distance ?? 0;
  return {
    point,
    radius,
    inside: buildKmerVpTree(distances.slice(0, middle).map((value) => value.candidate), profiles, compared),
    outside: buildKmerVpTree(distances.slice(middle).map((value) => value.candidate), profiles, compared),
  };
}

function nearestKmerPoint(
  tree: KmerVpNode,
  query: number,
  profiles: Map<number, KmerProfile>,
  abundance: (point: number) => number,
  compared: () => void,
): number {
  let best = tree.point;
  let bestDistance = Number.POSITIVE_INFINITY;
  const visit = (node: KmerVpNode | null) => {
    if (!node) return;
    compared();
    // Exact distance is required for triangle-inequality pruning. An
    // early-aborted partial distance could select the wrong VP-tree branch.
    const distance = Math.sqrt(kmerSquaredDistance(profiles.get(query)!, profiles.get(node.point)!));
    if (distance < bestDistance || (distance === bestDistance && abundance(node.point) > abundance(best))) {
      best = node.point;
      bestDistance = distance;
    }
    if (!node.inside && !node.outside) return;
    if (distance < node.radius) {
      visit(node.inside);
      if (distance + bestDistance >= node.radius) visit(node.outside);
    } else {
      visit(node.outside);
      if (distance - bestDistance <= node.radius) visit(node.inside);
    }
  };
  visit(tree);
  return best;
}

function logGamma(value: number): number {
  const coefficients = [
    0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.50734327868691, -0.1385710952657201, 9.984369578019572e-6,
    1.505632735149312e-7,
  ];
  if (value < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value)) - logGamma(1 - value);
  let x = coefficients[0];
  const z = value - 1;
  for (let index = 1; index < coefficients.length; index += 1) x += coefficients[index] / (z + index);
  const t = z + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

function regularizedGammaP(shape: number, value: number): number {
  if (!(shape > 0) || value < 0 || Number.isNaN(value)) return Number.NaN;
  if (value === 0) return 0;
  const logScale = -value + shape * Math.log(value) - logGamma(shape);
  if (value < shape + 1) {
    let term = 1 / shape;
    let sum = term;
    let denominator = shape;
    for (let iteration = 0; iteration < 10_000; iteration += 1) {
      denominator += 1;
      term *= value / denominator;
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 1e-14) break;
    }
    return Math.min(1, Math.max(0, sum * Math.exp(logScale)));
  }
  let b = value + 1 - shape;
  const floor = 1e-300;
  let c = 1 / floor;
  let d = 1 / Math.max(floor, b);
  let fraction = d;
  for (let iteration = 1; iteration < 10_000; iteration += 1) {
    const coefficient = -iteration * (iteration - shape);
    b += 2;
    d = coefficient * d + b;
    if (Math.abs(d) < floor) d = floor;
    c = b + coefficient / c;
    if (Math.abs(c) < floor) c = floor;
    d = 1 / d;
    const delta = d * c;
    fraction *= delta;
    if (Math.abs(delta - 1) < 1e-14) break;
  }
  const q = Math.exp(logScale) * fraction;
  return Math.min(1, Math.max(0, 1 - q));
}

/** P(Poisson(lambda) > observed), matching Distributions.jl ccdf. */
export function poissonStrictUpperTail(observed: number, lambda: number): number {
  if (!(lambda > 0)) return 0;
  if (observed < 0) return 1;
  return regularizedGammaP(Math.floor(observed) + 1, lambda);
}

function alternativeCount(length: number, distance: number): number {
  if (distance <= 0) return 1;
  let combinations = 1;
  for (let index = 1; index <= distance; index += 1) combinations *= (length - distance + index) / index;
  return Math.max(1, combinations * 3 ** distance);
}

function sequenceBlocks(sequence: string, count: number): string[] {
  const result: string[] = [];
  for (let index = 0; index < count; index += 1) {
    const start = Math.floor(index * sequence.length / count);
    const end = Math.floor((index + 1) * sequence.length / count);
    result.push(sequence.slice(start, end));
  }
  return result;
}

interface IndexedEditSegment {
  index: number;
  start: number;
  length: number;
  value: string;
}

function indexedEditSegments(sequence: string, count: number): IndexedEditSegment[] {
  const result: IndexedEditSegment[] = [];
  for (let index = 0; index < count; index += 1) {
    const start = Math.floor(index * sequence.length / count);
    const end = Math.floor((index + 1) * sequence.length / count);
    result.push({ index, start, length: end - start, value: sequence.slice(start, end) });
  }
  return result;
}

export interface BoundedEditProfile {
  distance: number;
  substitutions: number;
  insertions: number;
  deletions: number;
}

/**
 * Allocation-bounded Ukkonen-style edit profiler. The narrow diagonal band is
 * O(length × maximum), and tie-breaking chooses fewer indels when two optimal
 * paths exist so the aggressive indel rule is not triggered by an arbitrary
 * alignment of equal-length sequences.
 */
function createBoundedEditProfiler(maximum: number): (parent: string, child: string) => BoundedEditProfile | null {
  const width = maximum * 2 + 3;
  const infinity = maximum + 1;
  let previousCost = new Int16Array(width);
  let previousSubstitutions = new Int16Array(width);
  let previousInsertions = new Int16Array(width);
  let previousDeletions = new Int16Array(width);
  let currentCost = new Int16Array(width);
  let currentSubstitutions = new Int16Array(width);
  let currentInsertions = new Int16Array(width);
  let currentDeletions = new Int16Array(width);

  return (parent: string, child: string): BoundedEditProfile | null => {
    if (Math.abs(parent.length - child.length) > maximum) return null;
    previousCost.fill(infinity);
    previousSubstitutions.fill(0);
    previousInsertions.fill(0);
    previousDeletions.fill(0);
    for (let column = 0; column <= Math.min(child.length, maximum); column += 1) {
      const offset = column + maximum + 1;
      previousCost[offset] = column;
      previousInsertions[offset] = column;
    }

    for (let row = 1; row <= parent.length; row += 1) {
      currentCost.fill(infinity);
      currentSubstitutions.fill(0);
      currentInsertions.fill(0);
      currentDeletions.fill(0);
      const firstColumn = Math.max(0, row - maximum);
      const lastColumn = Math.min(child.length, row + maximum);
      for (let column = firstColumn; column <= lastColumn; column += 1) {
        const offset = column - row + maximum + 1;
        let bestCost = infinity;
        let bestSubstitutions = 0;
        let bestInsertions = 0;
        let bestDeletions = 0;
        const consider = (cost: number, substitutions: number, insertions: number, deletions: number) => {
          const indels = insertions + deletions;
          const bestIndels = bestInsertions + bestDeletions;
          if (cost < bestCost ||
            (cost === bestCost && indels < bestIndels) ||
            (cost === bestCost && indels === bestIndels && substitutions < bestSubstitutions) ||
            (cost === bestCost && indels === bestIndels && substitutions === bestSubstitutions && deletions < bestDeletions)) {
            bestCost = cost;
            bestSubstitutions = substitutions;
            bestInsertions = insertions;
            bestDeletions = deletions;
          }
        };

        if (column > 0 && previousCost[offset] <= maximum) {
          const mismatch = parent[row - 1] === child[column - 1] ? 0 : 1;
          consider(
            previousCost[offset] + mismatch,
            previousSubstitutions[offset] + mismatch,
            previousInsertions[offset],
            previousDeletions[offset],
          );
        }
        if (previousCost[offset + 1] <= maximum) {
          consider(
            previousCost[offset + 1] + 1,
            previousSubstitutions[offset + 1],
            previousInsertions[offset + 1],
            previousDeletions[offset + 1] + 1,
          );
        }
        if (column > 0 && currentCost[offset - 1] <= maximum) {
          consider(
            currentCost[offset - 1] + 1,
            currentSubstitutions[offset - 1],
            currentInsertions[offset - 1] + 1,
            currentDeletions[offset - 1],
          );
        }
        currentCost[offset] = bestCost;
        currentSubstitutions[offset] = bestSubstitutions;
        currentInsertions[offset] = bestInsertions;
        currentDeletions[offset] = bestDeletions;
      }
      [previousCost, currentCost] = [currentCost, previousCost];
      [previousSubstitutions, currentSubstitutions] = [currentSubstitutions, previousSubstitutions];
      [previousInsertions, currentInsertions] = [currentInsertions, previousInsertions];
      [previousDeletions, currentDeletions] = [currentDeletions, previousDeletions];
    }
    const offset = child.length - parent.length + maximum + 1;
    const distance = previousCost[offset];
    if (distance > maximum) return null;
    return {
      distance,
      substitutions: previousSubstitutions[offset],
      insertions: previousInsertions[offset],
      deletions: previousDeletions[offset],
    };
  };
}

export function boundedEditProfile(parent: string, child: string, maximum: number): BoundedEditProfile | null {
  if (!Number.isInteger(maximum) || maximum < 0 || maximum > 8) throw new Error("The bounded edit maximum must be an integer from 0 to 8.");
  return createBoundedEditProfiler(maximum)(parent, child);
}

export interface DenoisePartitionJobVariant {
  sequence: string;
  representative: number;
  count: number;
}

export interface DenoisePartitionJob {
  id: number;
  /** Global variant ordinals retained by the accumulator for deterministic merge. */
  variantIndices: Int32Array;
  variants: DenoisePartitionJobVariant[];
  options: DenoiseOptions;
}

export interface DenoisePartitionJobResult extends DenoisePartitionStats {
  id: number;
  /** Local target index for every local variant. */
  targets: Int32Array;
  variantWork: number;
}

/**
 * Pure partition kernel shared by browser workers, CLI worker threads, and the
 * serial fallback. A V/J partition is an exact independence boundary for all
 * three denoisers; splitting inside it would change parent/centroid choices.
 */
export function runDenoisePartitionJob(job: DenoisePartitionJob): DenoisePartitionJobResult {
  const variants = job.variants;
  if (job.variantIndices.length !== variants.length) throw new Error(`Denoising partition ${job.id} has inconsistent variant vectors.`);
  const targets = Int32Array.from(variants, (_, index) => index);
  const ordered = Array.from({ length: variants.length }, (_, index) => index)
    .sort((left, right) => variants[right].count - variants[left].count || variants[left].representative - variants[right].representative);
  const options = job.options;

  if (options.mode === "fad") {
    const maximumSquared = Math.max(0, Math.floor(12 * options.fadNeighborThreshold + 1e-9));
    const blockCount = Math.max(1, maximumSquared + 1);
    const profiles = new Map<number, KmerProfile>();
    for (let index = 0; index < variants.length; index += 1) profiles.set(index, kmerProfile(variants[index].sequence, blockCount));
    const accepted: number[] = [];
    const blockIndex = new Map<string, number[]>();
    const candidates = new Set<number>();
    let comparisons = 0;
    let truncated = 0;
    const addAccepted = (variantIndex: number) => {
      accepted.push(variantIndex);
      const profile = profiles.get(variantIndex)!;
      profile.hashes.forEach((hash, block) => {
        const key = `${block}:${hash}`;
        const values = blockIndex.get(key);
        if (values) values.push(variantIndex);
        else blockIndex.set(key, [variantIndex]);
      });
    };
    for (const variantIndex of ordered.filter((value) => variants[value].count >= options.minimumParentCount)) {
      const profile = profiles.get(variantIndex)!;
      candidates.clear();
      profile.hashes.forEach((hash, block) => {
        if (candidates.size >= options.maxCandidatesPerVariant) return;
        for (const candidate of blockIndex.get(`${block}:${hash}`) ?? []) {
          if (candidates.size >= options.maxCandidatesPerVariant) break;
          candidates.add(candidate);
        }
      });
      if (candidates.size >= options.maxCandidatesPerVariant) truncated += 1;
      const neighbors: Array<{ index: number; distance: number }> = [];
      for (const candidate of candidates) {
        comparisons += 1;
        const distance = kmerSquaredDistance(profile, profiles.get(candidate)!, maximumSquared);
        if (distance <= maximumSquared) neighbors.push({ index: candidate, distance });
      }
      if (!neighbors.length) {
        addAccepted(variantIndex);
        continue;
      }
      neighbors.sort((left, right) => variants[right.index].count - variants[left.index].count || left.distance - right.distance || variants[left.index].representative - variants[right.index].representative);
      const parent = neighbors[0].index;
      const child = variants[variantIndex];
      const lambda = variants[parent].count / Math.max(Number.MIN_VALUE, options.expectedZeroErrorFraction) * options.errorRate;
      const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * (child.sequence.length || 1));
      if (options.fadMethod === 2 && adjusted < options.alpha) addAccepted(variantIndex);
      else targets[variantIndex] = parent;
    }
    if (!accepted.length && ordered.length) addAccepted(ordered[0]);
    const acceptedSet = new Set(accepted);
    const vpTree = buildKmerVpTree(accepted, profiles, () => { comparisons += 1; });
    for (const variantIndex of ordered) {
      targets[variantIndex] = acceptedSet.has(variantIndex)
        ? variantIndex
        : vpTree
          ? nearestKmerPoint(vpTree, variantIndex, profiles, (point) => variants[point].count, () => { comparisons += 1; })
          : variantIndex;
    }
    return {
      id: job.id,
      targets,
      comparisons,
      truncated,
      indelMergedVariants: 0,
      substitutionMergedVariants: 0,
      variantWork: variants.length + variants.reduce((count, variant) => count + (variant.count >= options.minimumParentCount ? 1 : 0), 0),
    };
  }

  if (options.mode === "indel") {
    const distanceLimit = options.maximumEditDistance;
    const blockCount = distanceLimit + 1;
    const parentIndex = new Map<string, number[]>();
    const shortParentsByLength = new Map<number, number[]>();
    const profileEdit = createBoundedEditProfiler(distanceLimit);
    let comparisons = 0;
    let truncated = 0;
    let indelMergedVariants = 0;
    let substitutionMergedVariants = 0;
    const candidates = new Set<number>();
    const addParent = (variantIndex: number) => {
      if (variants[variantIndex].count < options.minimumParentCount) return;
      const sequence = variants[variantIndex].sequence;
      const segments = indexedEditSegments(sequence, blockCount);
      if (segments.some((segment) => segment.length === 0)) {
        const values = shortParentsByLength.get(sequence.length);
        if (values) values.push(variantIndex);
        else shortParentsByLength.set(sequence.length, [variantIndex]);
        return;
      }
      for (const segment of segments) {
        const key = `${sequence.length}:${segment.index}:${segment.value}`;
        const values = parentIndex.get(key);
        if (values) values.push(variantIndex);
        else parentIndex.set(key, [variantIndex]);
      }
    };
    for (const variantIndex of ordered) {
      const child = variants[variantIndex];
      const sequence = child.sequence;
      candidates.clear();
      let capped = false;
      const addCandidates = (values: number[]) => {
        for (const candidate of values) {
          if (candidates.has(candidate)) continue;
          if (candidates.size >= options.maxCandidatesPerVariant) {
            capped = true;
            return;
          }
          candidates.add(candidate);
        }
      };
      const minimumParentLength = Math.max(1, sequence.length - distanceLimit);
      const maximumParentLength = sequence.length + distanceLimit;
      for (let parentLength = minimumParentLength; parentLength <= maximumParentLength && !capped; parentLength += 1) {
        addCandidates(shortParentsByLength.get(parentLength) ?? []);
        for (let segmentIndex = 0; segmentIndex < blockCount && !capped; segmentIndex += 1) {
          const parentStart = Math.floor(segmentIndex * parentLength / blockCount);
          const parentEnd = Math.floor((segmentIndex + 1) * parentLength / blockCount);
          const segmentLength = parentEnd - parentStart;
          if (!segmentLength) continue;
          const firstStart = Math.max(0, parentStart - distanceLimit);
          const lastStart = Math.min(sequence.length - segmentLength, parentStart + distanceLimit);
          for (let queryStart = firstStart; queryStart <= lastStart && !capped; queryStart += 1) {
            addCandidates(parentIndex.get(`${parentLength}:${segmentIndex}:${sequence.slice(queryStart, queryStart + segmentLength)}`) ?? []);
          }
        }
      }
      if (capped) truncated += 1;
      let best = -1;
      let bestProfile: BoundedEditProfile | null = null;
      for (const candidate of candidates) {
        const parent = variants[candidate];
        if (parent.count <= child.count) continue;
        comparisons += 1;
        const profile = profileEdit(parent.sequence, sequence);
        if (!profile || profile.distance < 1) continue;
        const indels = profile.insertions + profile.deletions;
        let plausible = false;
        if (indels > 0) plausible = parent.count / child.count >= options.minimumIndelParentRatio;
        else {
          const exactErrorProbability = (options.errorRate / 3) ** profile.substitutions * (1 - options.errorRate) ** Math.max(0, sequence.length - profile.substitutions);
          const lambda = parent.count * exactErrorProbability;
          const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, profile.substitutions));
          plausible = adjusted >= options.alpha;
        }
        if (!plausible) continue;
        const bestParent = best >= 0 ? variants[best] : null;
        const isBetter = !bestProfile ||
          profile.distance < bestProfile.distance ||
          (profile.distance === bestProfile.distance && profile.substitutions < bestProfile.substitutions) ||
          (profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count > (bestParent?.count ?? -1)) ||
          (profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count === bestParent?.count && parent.representative < (bestParent?.representative ?? Number.POSITIVE_INFINITY));
        if (isBetter) {
          best = candidate;
          bestProfile = profile;
        }
      }
      if (best >= 0 && bestProfile) {
        targets[variantIndex] = best;
        if (bestProfile.insertions + bestProfile.deletions > 0) indelMergedVariants += 1;
        else substitutionMergedVariants += 1;
      } else {
        targets[variantIndex] = variantIndex;
        addParent(variantIndex);
      }
    }
    return { id: job.id, targets, comparisons, truncated, indelMergedVariants, substitutionMergedVariants, variantWork: variants.length };
  }

  const distanceLimit = options.maximumHammingDistance;
  const blockCount = distanceLimit + 1;
  const parentIndex = new Map<string, number[]>();
  let comparisons = 0;
  let truncated = 0;
  let substitutionMergedVariants = 0;
  const candidates = new Set<number>();
  const addParent = (variantIndex: number) => {
    if (variants[variantIndex].count < options.minimumParentCount) return;
    const sequence = variants[variantIndex].sequence;
    sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
      const key = `${sequence.length}:${blockIndex}:${block}`;
      const values = parentIndex.get(key);
      if (values) values.push(variantIndex);
      else parentIndex.set(key, [variantIndex]);
    });
  };
  for (const variantIndex of ordered) {
    const child = variants[variantIndex];
    const sequence = child.sequence;
    candidates.clear();
    sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
      if (candidates.size >= options.maxCandidatesPerVariant) return;
      for (const candidate of parentIndex.get(`${sequence.length}:${blockIndex}:${block}`) ?? []) {
        if (candidates.size >= options.maxCandidatesPerVariant) break;
        candidates.add(candidate);
      }
    });
    if (candidates.size >= options.maxCandidatesPerVariant) truncated += 1;
    let best = -1;
    let bestLambda = -1;
    let bestDistance = Number.POSITIVE_INFINITY;
    for (const candidate of candidates) {
      comparisons += 1;
      const distance = hammingDistanceWithin(sequence, variants[candidate].sequence, distanceLimit, false);
      if (distance < 1 || distance > distanceLimit) continue;
      const parentCount = variants[candidate].count;
      const exactErrorProbability = (options.errorRate / 3) ** distance * (1 - options.errorRate) ** Math.max(0, sequence.length - distance);
      const lambda = parentCount * exactErrorProbability;
      const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, distance));
      if (adjusted >= options.alpha && (lambda > bestLambda || (lambda === bestLambda && distance < bestDistance))) {
        best = candidate;
        bestLambda = lambda;
        bestDistance = distance;
      }
    }
    if (best >= 0) {
      targets[variantIndex] = best;
      substitutionMergedVariants += 1;
    } else {
      targets[variantIndex] = variantIndex;
      addParent(variantIndex);
    }
  }
  return { id: job.id, targets, comparisons, truncated, indelMergedVariants: 0, substitutionMergedVariants, variantWork: variants.length };
}

/**
 * Streaming builder used by the post-analysis worker. Exact variants are
 * dereplicated while batches are scanned, sequences live in a compact 2-bit
 * arena, and only one V/J partition's temporary neighbor index exists at a
 * time during finalization.
 */
export class DenoiseAccumulator {
  private readonly arena = new PackedDnaArena();
  private readonly variants: DenoiseVariant[] = [];
  private readonly variantsByFingerprint = new Map<string, number[]>();
  private readonly variantByOrdinal: Int32Array;
  private readonly processed: Uint8Array;
  private readonly standalone: number[] = [];
  private readonly retainedAmbiguous = new Map<string, { representative: number; members: number[]; count: number }>();
  private excludedAmbiguous = 0;
  private unresolvedRecords = 0;
  private readonly records: PostAnalysisRecord[];
  private readonly options: DenoiseOptions;
  private partitionJobsCache: DenoisePartitionJob[] | null = null;

  constructor(records: PostAnalysisRecord[], options: DenoiseOptions) {
    this.records = records;
    this.options = options;
    this.variantByOrdinal = new Int32Array(records.length);
    this.variantByOrdinal.fill(-1);
    this.processed = new Uint8Array(records.length);
    if (!(options.errorRate > 0 && options.errorRate < 1)) throw new Error("The denoising error rate must be between 0 and 1.");
    if (!(options.alpha > 0 && options.alpha < 1)) throw new Error("The denoising alpha must be between 0 and 1.");
    if (options.fadNeighborThreshold < 0) throw new Error("The FAD neighbor threshold cannot be negative.");
    if (options.maximumHammingDistance < 1 || options.maximumHammingDistance > 4) throw new Error("The conservative Hamming radius must be from 1 to 4.");
    if (!Number.isInteger(options.maximumEditDistance) || options.maximumEditDistance < 1 || options.maximumEditDistance > 2) throw new Error("The indel-aware edit radius must be 1 or 2.");
    if (!(options.minimumIndelParentRatio > 1)) throw new Error("The indel parent:child abundance ratio must be greater than 1.");
  }

  add(ordinal: number, rawSequence: string) {
    if (this.partitionJobsCache) throw new Error("Cannot add reads after denoising partition finalization has started.");
    if (ordinal < 0 || ordinal >= this.records.length || this.processed[ordinal]) return;
    this.processed[ordinal] = 1;
    const record = this.records[ordinal];
    const sequence = normalizeNt(rawSequence);
    const partition = denoisePartition(record, this.options);
    if (!sequence || !partition) {
      if (this.options.unresolvedPolicy === "retain") this.standalone.push(ordinal);
      this.unresolvedRecords += 1;
      return;
    }
    if (sequence.includes("N")) {
      if (this.options.ambiguousPolicy === "exclude") {
        this.excludedAmbiguous += 1;
        return;
      }
      const key = `${partition}\u0000${sequence}`;
      const existing = this.retainedAmbiguous.get(key);
      const weight = Math.max(1, Math.floor(record.inputCount ?? 1));
      if (existing) {
        existing.members.push(ordinal);
        existing.count += weight;
      } else this.retainedAmbiguous.set(key, { representative: ordinal, members: [ordinal], count: weight });
      return;
    }
    const fingerprint = `${partition}\u0000${sequenceFingerprint(sequence)}`;
    const candidates = this.variantsByFingerprint.get(fingerprint) ?? [];
    let variantIndex = candidates.find((candidate) => this.arena.equals(this.variants[candidate].location, sequence));
    const weight = Math.max(1, Math.floor(record.inputCount ?? 1));
    if (variantIndex === undefined) {
      variantIndex = this.variants.length;
      this.variants.push({ location: this.arena.append(sequence), partition, representative: ordinal, count: weight, target: variantIndex });
      candidates.push(variantIndex);
      this.variantsByFingerprint.set(fingerprint, candidates);
    } else this.variants[variantIndex].count += weight;
    this.variantByOrdinal[ordinal] = variantIndex;
  }

  preparePartitionJobs(): DenoisePartitionJob[] {
    if (this.partitionJobsCache) return this.partitionJobsCache;
    for (let ordinal = 0; ordinal < this.records.length; ordinal += 1) {
      if (!this.processed[ordinal]) {
        if (this.options.unresolvedPolicy === "retain") this.standalone.push(ordinal);
        this.unresolvedRecords += 1;
      }
    }
    const partitions = new Map<string, number[]>();
    this.variants.forEach((variant, index) => {
      const values = partitions.get(variant.partition);
      if (values) values.push(index);
      else partitions.set(variant.partition, [index]);
    });
    this.partitionJobsCache = [...partitions.values()].map((group, id) => ({
      id,
      variantIndices: Int32Array.from(group),
      variants: group.map((index) => ({
        sequence: this.arena.decode(this.variants[index].location),
        representative: this.variants[index].representative,
        count: this.variants[index].count,
      })),
      options: this.options,
    }));
    return this.partitionJobsCache;
  }

  finish(onProgress?: (processed: number, total: number, phase: "variants" | "finalize") => void): DedupResult {
    const results = this.preparePartitionJobs().map(runDenoisePartitionJob);
    return this.finishWithPartitionResults(results, onProgress);
  }

  finishWithPartitionResults(
    results: DenoisePartitionJobResult[],
    onProgress?: (processed: number, total: number, phase: "variants" | "finalize") => void,
    reportVariantProgress = true,
  ): DedupResult {
    const jobs = this.preparePartitionJobs();
    const resultById = new Map(results.map((result) => [result.id, result] as const));
    let candidateComparisons = 0;
    let truncated = 0;
    let indelMergedVariants = 0;
    let substitutionMergedVariants = 0;
    const variantWork = Math.max(1, jobs.reduce((total, job) => total + (this.options.mode === "fad"
      ? job.variants.length + job.variants.reduce((count, variant) => count + (variant.count >= this.options.minimumParentCount ? 1 : 0), 0)
      : job.variants.length), 0));
    let processedVariantWork = 0;
    if (reportVariantProgress) onProgress?.(0, variantWork, "variants");
    for (const job of jobs) {
      const result = resultById.get(job.id);
      if (!result || result.targets.length !== job.variantIndices.length) throw new Error(`Denoising partition ${job.id} did not return a complete result.`);
      for (let localIndex = 0; localIndex < job.variantIndices.length; localIndex += 1) {
        const localTarget = result.targets[localIndex];
        if (localTarget < 0 || localTarget >= job.variantIndices.length) throw new Error(`Denoising partition ${job.id} returned an invalid target.`);
        this.variants[job.variantIndices[localIndex]].target = job.variantIndices[localTarget];
      }
      candidateComparisons += result.comparisons;
      truncated += result.truncated;
      indelMergedVariants += result.indelMergedVariants;
      substitutionMergedVariants += result.substitutionMergedVariants;
      processedVariantWork += result.variantWork;
      if (reportVariantProgress) onProgress?.(Math.min(processedVariantWork, variantWork), variantWork, "variants");
    }
    if (reportVariantProgress && !jobs.length) onProgress?.(variantWork, variantWork, "variants");

    const representatives = new Int32Array(this.records.length);
    representatives.fill(-1);
    const counts = new Uint32Array(this.records.length);
    const finalizeTotal = Math.max(1, this.records.length * 3 + this.variants.length + this.standalone.length);
    let finalized = 0;
    const finalizeStride = Math.max(1, Math.floor(finalizeTotal / 500));
    const reportFinalize = (increment = 1) => {
      finalized += increment;
      if (finalized >= finalizeTotal || finalized % finalizeStride === 0) onProgress?.(Math.min(finalized, finalizeTotal), finalizeTotal, "finalize");
    };
    onProgress?.(0, finalizeTotal, "finalize");
    for (let ordinal = 0; ordinal < this.records.length; ordinal += 1) {
      const variantIndex = this.variantByOrdinal[ordinal];
      if (variantIndex >= 0) {
        const target = this.variants[this.variants[variantIndex].target];
        representatives[ordinal] = target.representative;
      }
      reportFinalize();
    }
    for (const variant of this.variants) {
      const target = this.variants[variant.target];
      counts[target.representative] += variant.count;
      reportFinalize();
    }
    for (const group of this.retainedAmbiguous.values()) {
      counts[group.representative] = group.count;
      group.members.forEach((ordinal) => { representatives[ordinal] = group.representative; });
    }
    for (const ordinal of this.standalone) {
      const weight = Math.max(1, Math.floor(this.records[ordinal].inputCount ?? 1));
      counts[ordinal] = weight;
      representatives[ordinal] = ordinal;
      reportFinalize();
    }
    let uniqueRecords = 0;
    for (const count of counts) {
      if (count > 0) uniqueRecords += 1;
      reportFinalize();
    }
    const largestGroups = largestCountGroups(counts);
    const warnings: string[] = [];
    if (this.excludedAmbiguous) warnings.push(`${this.excludedAmbiguous.toLocaleString()} records containing ambiguous nucleotide symbols were excluded to match the selected policy.`);
    if (this.unresolvedRecords) warnings.push(
      `${this.unresolvedRecords.toLocaleString()} records without a usable trimmed sequence or both V/J calls were ${this.options.unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`,
    );
    if (truncated) warnings.push(`${truncated.toLocaleString()} variants reached the candidate cap; increase it before treating this denoising result as complete.`);
    let inputAbundance = 0;
    for (const record of this.records) {
      inputAbundance += Math.max(1, Math.floor(record.inputCount ?? 1));
      reportFinalize();
    }
    return {
      mode: this.options.mode,
      key: "trimmed",
      algorithm: this.options.mode === "fad"
        ? `FAD-compatible corrected 6-mer / method ${this.options.fadMethod}`
        : this.options.mode === "indel"
          ? "Indel-aware bounded edit model"
          : "Conservative exact-neighbor error model",
      inputRecords: this.records.length,
      inputAbundance,
      uniqueRecords,
      collapsedRecords: this.records.length - uniqueRecords,
      representatives,
      counts,
      largestGroups,
      partitions: jobs.length,
      candidateComparisons,
      indelMergedVariants,
      substitutionMergedVariants,
      excludedAmbiguous: this.excludedAmbiguous,
      unresolvedRecords: this.unresolvedRecords,
      warnings,
    };
  }

  private processFadPartition(group: number[], onVariant?: () => void): DenoisePartitionStats {
    const maximumSquared = Math.max(0, Math.floor(12 * this.options.fadNeighborThreshold + 1e-9));
    const blockCount = Math.max(1, maximumSquared + 1);
    const profiles = new Map<number, KmerProfile>();
    const sequences = new Map<number, string>();
    for (const index of group) {
      const sequence = this.arena.decode(this.variants[index].location);
      sequences.set(index, sequence);
      profiles.set(index, kmerProfile(sequence, blockCount));
    }
    const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
    const accepted: number[] = [];
    const index = new Map<string, number[]>();
    const candidates = new Set<number>();
    let comparisons = 0;
    let truncated = 0;
    const addAccepted = (variantIndex: number) => {
      accepted.push(variantIndex);
      const profile = profiles.get(variantIndex)!;
      profile.hashes.forEach((hash, block) => {
        const key = `${block}:${hash}`;
        const values = index.get(key);
        if (values) values.push(variantIndex);
        else index.set(key, [variantIndex]);
      });
    };
    for (const variantIndex of ordered.filter((value) => this.variants[value].count >= this.options.minimumParentCount)) {
      const profile = profiles.get(variantIndex)!;
      candidates.clear();
      profile.hashes.forEach((hash, block) => {
        if (candidates.size >= this.options.maxCandidatesPerVariant) return;
        for (const candidate of index.get(`${block}:${hash}`) ?? []) {
          if (candidates.size >= this.options.maxCandidatesPerVariant) break;
          candidates.add(candidate);
        }
      });
      if (candidates.size >= this.options.maxCandidatesPerVariant) truncated += 1;
      const neighbors: Array<{ index: number; distance: number }> = [];
      for (const candidate of candidates) {
        comparisons += 1;
        const distance = kmerSquaredDistance(profile, profiles.get(candidate)!, maximumSquared);
        if (distance <= maximumSquared) neighbors.push({ index: candidate, distance });
      }
      if (!neighbors.length) {
        addAccepted(variantIndex);
        onVariant?.();
        continue;
      }
      neighbors.sort((left, right) => this.variants[right.index].count - this.variants[left.index].count || left.distance - right.distance || this.variants[left.index].representative - this.variants[right.index].representative);
      const parent = neighbors[0].index;
      const child = this.variants[variantIndex];
      const parentCount = this.variants[parent].count;
      const lambda = parentCount / Math.max(Number.MIN_VALUE, this.options.expectedZeroErrorFraction) * this.options.errorRate;
      const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * (sequences.get(variantIndex)?.length ?? 1));
      if (this.options.fadMethod === 2 && adjusted < this.options.alpha) addAccepted(variantIndex);
      else child.target = parent;
      onVariant?.();
    }
    if (!accepted.length && ordered.length) addAccepted(ordered[0]);
    // The published FAD implementation assigns every non-template variant to
    // its globally nearest accepted corrected-k-mer centroid, even when it is
    // outside the template-selection radius. V/J partitioning bounds this
    // exact scan and preserves compatibility with that behavior.
    const acceptedSet = new Set(accepted);
    const vpTree = buildKmerVpTree(accepted, profiles, () => { comparisons += 1; });
    for (const variantIndex of ordered) {
      if (acceptedSet.has(variantIndex)) {
        this.variants[variantIndex].target = variantIndex;
        onVariant?.();
        continue;
      }
      this.variants[variantIndex].target = vpTree ? nearestKmerPoint(vpTree, variantIndex, profiles, (point) => this.variants[point].count, () => { comparisons += 1; }) : variantIndex;
      onVariant?.();
    }
    return { comparisons, truncated, indelMergedVariants: 0, substitutionMergedVariants: 0 };
  }

  private processIndelPartition(group: number[], onVariant?: () => void): DenoisePartitionStats {
    const distanceLimit = this.options.maximumEditDistance;
    const blockCount = distanceLimit + 1;
    const sequences = new Map(group.map((index) => [index, this.arena.decode(this.variants[index].location)]));
    const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
    const parentIndex = new Map<string, number[]>();
    const shortParentsByLength = new Map<number, number[]>();
    const profileEdit = createBoundedEditProfiler(distanceLimit);
    let comparisons = 0;
    let truncated = 0;
    let indelMergedVariants = 0;
    let substitutionMergedVariants = 0;
    const candidates = new Set<number>();

    const addParent = (variantIndex: number) => {
      if (this.variants[variantIndex].count < this.options.minimumParentCount) return;
      const sequence = sequences.get(variantIndex)!;
      const segments = indexedEditSegments(sequence, blockCount);
      if (segments.some((segment) => segment.length === 0)) {
        const values = shortParentsByLength.get(sequence.length);
        if (values) values.push(variantIndex);
        else shortParentsByLength.set(sequence.length, [variantIndex]);
        return;
      }
      for (const segment of segments) {
        const key = `${sequence.length}:${segment.index}:${segment.value}`;
        const values = parentIndex.get(key);
        if (values) values.push(variantIndex);
        else parentIndex.set(key, [variantIndex]);
      }
    };

    for (const variantIndex of ordered) {
      const child = this.variants[variantIndex];
      const sequence = sequences.get(variantIndex)!;
      candidates.clear();
      let capped = false;
      const addCandidates = (values: number[]) => {
        for (const candidate of values) {
          if (candidates.has(candidate)) continue;
          if (candidates.size >= this.options.maxCandidatesPerVariant) {
            capped = true;
            return;
          }
          candidates.add(candidate);
        }
      };

      // Complete bounded-edit join: a parent split into d+1 disjoint segments
      // must retain at least one exact segment under at most d edits. Probe the
      // segment at every legal length and shifted start; verification below is
      // exact and allocation-bounded.
      const minimumParentLength = Math.max(1, sequence.length - distanceLimit);
      const maximumParentLength = sequence.length + distanceLimit;
      for (let parentLength = minimumParentLength; parentLength <= maximumParentLength && !capped; parentLength += 1) {
        addCandidates(shortParentsByLength.get(parentLength) ?? []);
        for (let segmentIndex = 0; segmentIndex < blockCount && !capped; segmentIndex += 1) {
          const parentStart = Math.floor(segmentIndex * parentLength / blockCount);
          const parentEnd = Math.floor((segmentIndex + 1) * parentLength / blockCount);
          const segmentLength = parentEnd - parentStart;
          if (!segmentLength) continue;
          const firstStart = Math.max(0, parentStart - distanceLimit);
          const lastStart = Math.min(sequence.length - segmentLength, parentStart + distanceLimit);
          for (let queryStart = firstStart; queryStart <= lastStart && !capped; queryStart += 1) {
            addCandidates(parentIndex.get(`${parentLength}:${segmentIndex}:${sequence.slice(queryStart, queryStart + segmentLength)}`) ?? []);
          }
        }
      }
      if (capped) truncated += 1;

      let best = -1;
      let bestProfile: BoundedEditProfile | null = null;
      for (const candidate of candidates) {
        const parent = this.variants[candidate];
        if (parent.count <= child.count) continue;
        comparisons += 1;
        const profile = profileEdit(sequences.get(candidate)!, sequence);
        if (!profile || profile.distance < 1) continue;
        const indels = profile.insertions + profile.deletions;
        let plausible = false;
        if (indels > 0) {
          plausible = parent.count / child.count >= this.options.minimumIndelParentRatio;
        } else {
          const exactErrorProbability = (this.options.errorRate / 3) ** profile.substitutions * (1 - this.options.errorRate) ** Math.max(0, sequence.length - profile.substitutions);
          const lambda = parent.count * exactErrorProbability;
          const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, profile.substitutions));
          plausible = adjusted >= this.options.alpha;
        }
        if (!plausible) continue;
        const bestParent = best >= 0 ? this.variants[best] : null;
        const isBetter = !bestProfile ||
          profile.distance < bestProfile.distance ||
          (profile.distance === bestProfile.distance && profile.substitutions < bestProfile.substitutions) ||
          (profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count > (bestParent?.count ?? -1)) ||
          (profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count === bestParent?.count && parent.representative < (bestParent?.representative ?? Number.POSITIVE_INFINITY));
        if (isBetter) {
          best = candidate;
          bestProfile = profile;
        }
      }
      if (best >= 0 && bestProfile) {
        child.target = best;
        if (bestProfile.insertions + bestProfile.deletions > 0) indelMergedVariants += 1;
        else substitutionMergedVariants += 1;
      } else {
        child.target = variantIndex;
        addParent(variantIndex);
      }
      onVariant?.();
    }
    return { comparisons, truncated, indelMergedVariants, substitutionMergedVariants };
  }

  private processConservativePartition(group: number[], onVariant?: () => void): DenoisePartitionStats {
    const distanceLimit = this.options.maximumHammingDistance;
    const blockCount = distanceLimit + 1;
    const sequences = new Map(group.map((index) => [index, this.arena.decode(this.variants[index].location)]));
    const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
    const parentIndex = new Map<string, number[]>();
    let comparisons = 0;
    let truncated = 0;
    let substitutionMergedVariants = 0;
    const candidates = new Set<number>();
    const addParent = (variantIndex: number) => {
      if (this.variants[variantIndex].count < this.options.minimumParentCount) return;
      const sequence = sequences.get(variantIndex)!;
      sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
        const key = `${sequence.length}:${blockIndex}:${block}`;
        const values = parentIndex.get(key);
        if (values) values.push(variantIndex);
        else parentIndex.set(key, [variantIndex]);
      });
    };
    for (const variantIndex of ordered) {
      const child = this.variants[variantIndex];
      const sequence = sequences.get(variantIndex)!;
      candidates.clear();
      sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
        if (candidates.size >= this.options.maxCandidatesPerVariant) return;
        for (const candidate of parentIndex.get(`${sequence.length}:${blockIndex}:${block}`) ?? []) {
          if (candidates.size >= this.options.maxCandidatesPerVariant) break;
          candidates.add(candidate);
        }
      });
      if (candidates.size >= this.options.maxCandidatesPerVariant) truncated += 1;
      let best = -1;
      let bestLambda = -1;
      let bestDistance = Number.POSITIVE_INFINITY;
      for (const candidate of candidates) {
        comparisons += 1;
        const distance = hammingDistanceWithin(sequence, sequences.get(candidate)!, distanceLimit, false);
        if (distance < 1 || distance > distanceLimit) continue;
        const parentCount = this.variants[candidate].count;
        const exactErrorProbability = (this.options.errorRate / 3) ** distance * (1 - this.options.errorRate) ** Math.max(0, sequence.length - distance);
        const lambda = parentCount * exactErrorProbability;
        const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, distance));
        // A child is collapsed only when its abundance is statistically
        // compatible with sequencing error from a more abundant accepted read.
        if (adjusted >= this.options.alpha && (lambda > bestLambda || (lambda === bestLambda && distance < bestDistance))) {
          best = candidate;
          bestLambda = lambda;
          bestDistance = distance;
        }
      }
      if (best >= 0) {
        child.target = best;
        substitutionMergedVariants += 1;
      }
      else {
        child.target = variantIndex;
        addParent(variantIndex);
      }
      onVariant?.();
    }
    return { comparisons, truncated, indelMergedVariants: 0, substitutionMergedVariants };
  }
}

class UnionFind {
  readonly parent: Int32Array;
  private readonly rank: Uint8Array;

  constructor(size: number) {
    this.parent = new Int32Array(size);
    this.rank = new Uint8Array(size);
    for (let index = 0; index < size; index += 1) this.parent[index] = index;
  }

  find(value: number): number {
    let root = value;
    while (this.parent[root] !== root) root = this.parent[root];
    while (this.parent[value] !== value) {
      const next = this.parent[value];
      this.parent[value] = root;
      value = next;
    }
    return root;
  }

  union(left: number, right: number) {
    let a = this.find(left);
    let b = this.find(right);
    if (a === b) return;
    if (this.rank[a] < this.rank[b]) [a, b] = [b, a];
    this.parent[b] = a;
    if (this.rank[a] === this.rank[b]) this.rank[a] += 1;
  }
}

function blocks(value: string, count: number): Array<{ index: number; value: string }> {
  const result: Array<{ index: number; value: string }> = [];
  for (let index = 0; index < count; index += 1) {
    const start = Math.floor(index * value.length / count);
    const end = Math.floor((index + 1) * value.length / count);
    result.push({ index, value: value.slice(start, end) });
  }
  return result;
}

function compatibleRecords(left: PostAnalysisRecord, right: PostAnalysisRecord, options: LineageOptions): boolean {
  if (datasetScopeKey(left, options.scope ?? "global") !== datasetScopeKey(right, options.scope ?? "global")) return false;
  if (options.requireSameLocus && left.locus !== right.locus) return false;
  return callsCompatible(left.vCall, right.vCall, options.callResolution, options.ambiguity) &&
    callsCompatible(left.jCall, right.jCall, options.callResolution, options.ambiguity);
}

function recordIndexTokens(record: PostAnalysisRecord, options: LineageOptions): string[] {
  const v = callSet(record.vCall, options.callResolution, options.ambiguity);
  const j = callSet(record.jCall, options.callResolution, options.ambiguity);
  if (!v.length || !j.length) return [];
  const scope = `${datasetScopeKey(record, options.scope ?? "global")}\u0001`;
  if (options.ambiguity === "strict") return [`${scope}${v.join("+")}\u0001${j.join("+")}`];
  return v.flatMap((vCall) => j.map((jCall) => `${scope}${vCall}\u0001${jCall}`));
}

export function assignLineages(
  records: PostAnalysisRecord[],
  options: LineageOptions,
  dedup?: DedupResult,
  activeMask?: Uint8Array,
  doubleDMask?: Uint8Array,
): LineageResult {
  const union = new UnionFind(records.length);
  const bucket = new Map<string, number[]>();
  const exact = new Map<string, number>();
  let candidateComparisons = 0;
  let truncatedCandidates = 0;
  let assignedRecords = 0;
  let activeAbundance = 0;

  for (let index = 0; index < records.length; index += 1) {
    const record = records[index];
    const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
    if (!weight) continue;
    activeAbundance += weight;
    const cdr3 = normalizeNt(record.cdr3Nt);
    const tokens = recordIndexTokens(record, options);
    if (!cdr3 || !tokens.length || (options.productiveOnly && !record.productive)) continue;
    assignedRecords += weight;
    const distanceLimit = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, distanceLimit + 1));
    const prefix = `${options.requireSameLocus ? record.locus : "*"}\u0000${cdr3.length}\u0000`;
    const exactSignature = `${prefix}${tokens.join("|")}\u0000${cdr3}`;
    const exactRepresentative = exact.get(exactSignature);
    if (exactRepresentative !== undefined) {
      union.union(index, exactRepresentative);
      continue;
    }
    exact.set(exactSignature, index);

    const candidates = new Set<number>();
    for (const token of tokens) {
      for (const block of blocks(cdr3, blockCount)) {
        const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
        for (const candidate of bucket.get(key) ?? []) {
          candidates.add(candidate);
          if (candidates.size >= options.maxCandidateComparisons) break;
        }
        if (candidates.size >= options.maxCandidateComparisons) break;
      }
      if (candidates.size >= options.maxCandidateComparisons) break;
    }
    if (candidates.size >= options.maxCandidateComparisons) truncatedCandidates += 1;
    for (const candidate of candidates) {
      candidateComparisons += 1;
      const other = records[candidate];
      if (!compatibleRecords(record, other, options)) continue;
      if (hammingDistanceWithin(cdr3, normalizeNt(other.cdr3Nt), distanceLimit, false) <= distanceLimit) {
        union.union(index, candidate);
      }
    }
    for (const token of tokens) {
      for (const block of blocks(cdr3, blockCount)) {
        const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
        const values = bucket.get(key);
        if (values) values.push(index);
        else bucket.set(key, [index]);
      }
    }
  }

  const eligible = new Uint8Array(records.length);
  const uniqueByRoot = new Uint32Array(records.length);
  const abundanceByRoot = new Float64Array(records.length);
  const representativeByRoot = new Int32Array(records.length);
  representativeByRoot.fill(-1);
  for (let index = 0; index < records.length; index += 1) {
    const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
    if (!weight) continue;
    const record = records[index];
    if (!record.cdr3Nt || !recordIndexTokens(record, options).length || (options.productiveOnly && !record.productive)) continue;
    const root = union.find(index);
    eligible[index] = 1;
    uniqueByRoot[root] += 1;
    abundanceByRoot[root] += weight;
    if (representativeByRoot[root] < 0) representativeByRoot[root] = index;
  }

  const assignments = new Int32Array(records.length);
  const rootIds = new Int32Array(records.length);
  const top: LineageSummary[] = [];
  const summaryLimit = 10_000;
  const sizeBins = new Uint32Array(7);
  const vUsageMap = new Map<string, { lineages: number; abundance: number }>();
  const jUsageMap = new Map<string, { lineages: number; abundance: number }>();
  let lineageCount = 0;
  const smaller = (left: LineageSummary, right: LineageSummary) =>
    left.abundance < right.abundance ||
    (left.abundance === right.abundance && left.uniqueMembers < right.uniqueMembers) ||
    (left.abundance === right.abundance && left.uniqueMembers === right.uniqueMembers && left.id > right.id);
  const retainTop = (summary: LineageSummary) => {
    if (top.length < summaryLimit) {
      top.push(summary);
      let index = top.length - 1;
      while (index > 0) {
        const parent = Math.floor((index - 1) / 2);
        if (!smaller(top[index], top[parent])) break;
        [top[index], top[parent]] = [top[parent], top[index]];
        index = parent;
      }
      return;
    }
    if (!smaller(top[0], summary)) return;
    top[0] = summary;
    let index = 0;
    while (true) {
      const left = index * 2 + 1;
      const right = left + 1;
      let smallest = index;
      if (left < top.length && smaller(top[left], top[smallest])) smallest = left;
      if (right < top.length && smaller(top[right], top[smallest])) smallest = right;
      if (smallest === index) break;
      [top[index], top[smallest]] = [top[smallest], top[index]];
      index = smallest;
    }
  };
  const addUsage = (map: Map<string, { lineages: number; abundance: number }>, call: string, abundance: number) => {
    const value = map.get(call);
    if (value) {
      value.lineages += 1;
      value.abundance += abundance;
    } else map.set(call, { lineages: 1, abundance });
  };
  for (let root = 0; root < records.length; root += 1) {
    if (!uniqueByRoot[root]) continue;
    const id = ++lineageCount;
    rootIds[root] = id;
    const representativeOrdinal = representativeByRoot[root];
    const representative = records[representativeOrdinal];
    const abundance = abundanceByRoot[root];
    const bin = abundance === 1 ? 0 : abundance <= 3 ? 1 : abundance <= 9 ? 2 : abundance <= 24 ? 3 : abundance <= 99 ? 4 : abundance <= 499 ? 5 : 6;
    sizeBins[bin] += 1;
    const vCalls = callSet(representative.vCall, options.callResolution, options.ambiguity);
    const jCalls = callSet(representative.jCall, options.callResolution, options.ambiguity);
    vCalls.forEach((call) => addUsage(vUsageMap, call, abundance));
    jCalls.forEach((call) => addUsage(jUsageMap, call, abundance));
    retainTop({
      id,
      representativeOrdinal,
      uniqueMembers: uniqueByRoot[root],
      abundance,
      locus: representative.locus,
      vCalls,
      jCalls,
      cdr3Length: representative.cdr3Nt.length,
      representativeCdr3Nt: representative.cdr3Nt,
      representativeCdr3Aa: representative.cdr3Aa,
      studyScope: options.scope ?? "global",
      studyGroup: datasetScopeValue(representative, options.scope ?? "global"),
      sampleIds: [],
      sampleCounts: [],
      subjectIds: [],
      timepoints: [],
      compartments: [],
      doubleDPositiveMembers: 0,
      doubleDPositiveAbundance: 0,
    });
  }
  for (let index = 0; index < records.length; index += 1) {
    if (eligible[index]) assignments[index] = rootIds[union.find(index)];
  }
  if (dedup) {
    for (let index = 0; index < dedup.representatives.length; index += 1) {
      const representative = dedup.representatives[index];
      assignments[index] = representative >= 0 ? assignments[representative] : 0;
    }
  }
  // Enrich only the bounded set retained for interactive browsing. This keeps
  // memory proportional to the table cap while preserving exact sample
  // membership for every summary the user can filter or open.
  const summaryById = new Map(top.map((summary) => [summary.id, summary]));
  const samplesById = new Map<number, Set<string>>();
  const sampleCountsById = new Map<number, Map<string, { representatives: Set<number>; abundance: number }>>();
  const subjectsById = new Map<number, Set<string>>();
  const timepointsById = new Map<number, Set<string>>();
  const compartmentsById = new Map<number, Set<string>>();
  const addValue = (map: Map<number, Set<string>>, id: number, value?: string) => {
    if (!value) return;
    const values = map.get(id);
    if (values) values.add(value); else map.set(id, new Set([value]));
  };
  for (let index = 0; index < records.length; index += 1) {
    const representative=dedup?dedup.representatives[index]:index;
    if(representative<0||(activeMask&&!activeMask[representative]))continue;
    const sourceWeight=Math.max(1,Math.floor(records[index].inputCount??1));
    const id = assignments[representative];
    if (!id || !summaryById.has(id)) continue;
    const record = records[index];
    addValue(samplesById, id, record.sampleId);
    if (record.sampleId) {
      let bySample = sampleCountsById.get(id);
      if (!bySample) {
        bySample = new Map();
        sampleCountsById.set(id, bySample);
      }
      const counts = bySample.get(record.sampleId) ?? { representatives: new Set<number>(), abundance: 0 };
      counts.representatives.add(representative);
      counts.abundance += sourceWeight;
      bySample.set(record.sampleId, counts);
    }
    addValue(subjectsById, id, record.subjectId);
    addValue(timepointsById, id, record.timepoint);
    addValue(compartmentsById, id, record.compartment);
  }
  for(let index=0;index<records.length;index+=1){
    const activeWeight=dedup?dedup.counts[index]:Math.max(1,Math.floor(records[index].inputCount??1));
    if(!activeWeight||(activeMask&&!activeMask[index])||!doubleDMask?.[index])continue;
    const id=assignments[index];
    if(summaryById.has(id)){
      const summary = summaryById.get(id)!;
      summary.doubleDPositiveMembers = (summary.doubleDPositiveMembers ?? 0) + 1;
      summary.doubleDPositiveAbundance = (summary.doubleDPositiveAbundance ?? 0) + activeWeight;
    }
  }
  for (const summary of top) {
    summary.sampleIds = [...(samplesById.get(summary.id) ?? [])].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
    summary.sampleCounts = [...(sampleCountsById.get(summary.id) ?? [])]
      .map(([sampleId, counts]) => ({ sampleId, uniqueMembers: counts.representatives.size, abundance: counts.abundance }))
      .sort((left, right) => right.abundance - left.abundance || left.sampleId.localeCompare(right.sampleId, undefined, { numeric: true }));
    summary.subjectIds = [...(subjectsById.get(summary.id) ?? [])].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
    summary.timepoints = [...(timepointsById.get(summary.id) ?? [])].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
    summary.compartments = [...(compartmentsById.get(summary.id) ?? [])].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
  }
  const usage = (map: Map<string, { lineages: number; abundance: number }>) => [...map.entries()]
    .map(([call, value]) => ({ call, ...value }))
    .sort((a, b) => b.abundance - a.abundance || b.lineages - a.lineages || a.call.localeCompare(b.call));

  return {
    assignments,
    summaries: top.sort((a, b) => b.abundance - a.abundance || b.uniqueMembers - a.uniqueMembers || a.id - b.id),
    lineageCount,
    sizeHistogram: ["1", "2–3", "4–9", "10–24", "25–99", "100–499", "500+"].map((label, index) => ({ label, count: sizeBins[index] })),
    vUsage: usage(vUsageMap),
    jUsage: usage(jUsageMap),
    assignedRecords,
    unassignedRecords: activeAbundance - assignedRecords,
    candidateComparisons,
    truncatedCandidates,
  };
}

interface Cdr3BkNode {
  sequence: string;
  ordinals: number[];
  children: Map<number, Cdr3BkNode>;
}

interface Cdr3BkPartition {
  root: Cdr3BkNode | null;
  minimumLength: number;
  maximumLength: number;
}

function exactEditDistance(left: string, right: string): number {
  return bandedEditDistance(left, right, Math.max(left.length, right.length));
}

function insertCdr3Bk(partition: Cdr3BkPartition, sequence: string, ordinal: number) {
  partition.minimumLength = Math.min(partition.minimumLength, sequence.length);
  partition.maximumLength = Math.max(partition.maximumLength, sequence.length);
  if (!partition.root) {
    partition.root = { sequence, ordinals: [ordinal], children: new Map() };
    return;
  }
  let node = partition.root;
  while (true) {
    const distance = exactEditDistance(sequence, node.sequence);
    if (!distance) {
      node.ordinals.push(ordinal);
      return;
    }
    const child = node.children.get(distance);
    if (child) node = child;
    else {
      node.children.set(distance, { sequence, ordinals: [ordinal], children: new Map() });
      return;
    }
  }
}

function searchCdr3Bk(
  node: Cdr3BkNode | null,
  query: string,
  radius: number,
  visit: (node: Cdr3BkNode, distance: number) => boolean,
): boolean {
  if (!node) return true;
  const distance = exactEditDistance(query, node.sequence);
  if (distance <= radius && !visit(node, distance)) return false;
  const minimum = Math.max(0, distance - radius);
  const maximum = distance + radius;
  for (const [edge, child] of node.children) {
    if (edge < minimum || edge > maximum) continue;
    if (!searchCdr3Bk(child, query, radius, visit)) return false;
  }
  return true;
}

function neighbourHardPartition(record: PostAnalysisRecord, options: LineageOptions): string {
  return `${datasetScopeKey(record, options.scope ?? "global")}\u0001${options.requireSameLocus ? record.locus : "*"}`;
}

function neighbourCallAgreement(
  left: PostAnalysisRecord,
  right: PostAnalysisRecord,
  options: LineageNeighbourOptions,
): { accepted: boolean; label: "both" | "v" | "j" | "none" } {
  const v = callsCompatible(left.vCall, right.vCall, options.callResolution, options.ambiguity);
  const j = callsCompatible(left.jCall, right.jCall, options.callResolution, options.ambiguity);
  const label = v && j ? "both" : v ? "v" : j ? "j" : "none";
  const policy = options.callPolicy ?? "both";
  return { accepted: policy === "ignore" || (policy === "either" ? v || j : v && j), label };
}

function maximumEditSearchRadius(sourceLength: number, partition: Cdr3BkPartition, identity: number): number {
  let maximum = 0;
  if (!Number.isFinite(partition.minimumLength)) return maximum;
  for (let length = partition.minimumLength; length <= partition.maximumLength; length += 1) {
    const allowed = Math.floor((1 - identity) * Math.max(sourceLength, length) + 1e-9);
    if (Math.abs(sourceLength - length) <= allowed) maximum = Math.max(maximum, allowed);
  }
  return maximum;
}

/**
 * Read-only search across already-assigned lineage boundaries. The default
 * retains assignment-compatible equal-length Hamming search. Indel-aware mode
 * uses an exact BK-tree shortlist followed by pair-specific banded Levenshtein
 * verification; V/J evidence can remain strict, require either call, or be
 * ignored while study boundary and locus remain hard constraints.
 */
export function findLineageNeighbours(
  records: PostAnalysisRecord[],
  assignments: Int32Array,
  options: LineageNeighbourOptions,
  dedup?: DedupResult,
  activeMask?: Uint8Array,
): LineageNeighbourResult {
  if (assignments.length < records.length) throw new Error("Lineage assignments do not cover the indexed repertoire.");
  const sourceIds = new Set(options.sourceLineageIds.filter((value) => value > 0));
  if (!sourceIds.size) throw new Error("Choose at least one assigned lineage before searching for neighbours.");
  if (!(options.minimumIdentity >= 0 && options.minimumIdentity <= 1)) throw new Error("Neighbour CDR3 identity must be between zero and one.");
  const metric = options.metric ?? "hamming";
  const callPolicy = options.callPolicy ?? "both";
  const lineageCount = assignments.reduce((maximum, value) => Math.max(maximum, value), 0);
  const uniqueMembers = new Uint32Array(lineageCount + 1);
  const abundance = new Float64Array(lineageCount + 1);
  const bucket = new Map<string, number[]>();
  const bkPartitions = new Map<string, Cdr3BkPartition>();
  const sourceRecords: number[] = [];
  let indexedRecords = 0;

  for (let index = 0; index < records.length; index += 1) {
    if (activeMask && !activeMask[index]) continue;
    const lineageId = assignments[index];
    if (!(lineageId > 0)) continue;
    const record = records[index];
    const cdr3 = normalizeNt(record.cdr3Nt);
    const tokens = callPolicy === "both" ? recordIndexTokens(record, options) : ["*"];
    if (!cdr3 || !tokens.length || (options.productiveOnly && !record.productive)) continue;
    const weight = dedup ? dedup.counts[index] : Math.max(1, Math.floor(record.inputCount ?? 1));
    if (dedup && !weight) continue;
    uniqueMembers[lineageId] += 1;
    abundance[lineageId] += weight;
    indexedRecords += 1;
    if (sourceIds.has(lineageId)) sourceRecords.push(index);
    const hard = neighbourHardPartition(record, options);
    if (metric === "edit") {
      let partition = bkPartitions.get(hard);
      if (!partition) {
        partition = { root: null, minimumLength: Number.POSITIVE_INFINITY, maximumLength: 0 };
        bkPartitions.set(hard, partition);
      }
      insertCdr3Bk(partition, cdr3, index);
      continue;
    }
    const distanceLimit = Math.floor((1 - options.minimumIdentity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, distanceLimit + 1));
    for (const token of tokens) {
      for (const block of blocks(cdr3, blockCount)) {
        const key = `${hard}\u0000${cdr3.length}\u0000${token}\u0000${block.index}\u0000${block.value}`;
        const values = bucket.get(key);
        if (values) values.push(index); else bucket.set(key, [index]);
      }
    }
  }

  const best = new Map<number, LineageNeighbourHit>();
  let candidateComparisons = 0;
  let truncatedSourceRecords = 0;
  const agreementRank = (value: LineageNeighbourHit["callAgreement"]) => value === "both" ? 2 : value === "v" || value === "j" ? 1 : 0;
  const consider = (sourceIndex: number, candidateIndex: number, distance: number) => {
    const candidateLineage = assignments[candidateIndex];
    if (candidateLineage <= 0 || sourceIds.has(candidateLineage)) return;
    candidateComparisons += 1;
    const source = records[sourceIndex];
    const candidate = records[candidateIndex];
    if (neighbourHardPartition(source, options) !== neighbourHardPartition(candidate, options)) return;
    const agreement = neighbourCallAgreement(source, candidate, options);
    if (!agreement.accepted) return;
    const sourceCdr3 = normalizeNt(source.cdr3Nt);
    const candidateCdr3 = normalizeNt(candidate.cdr3Nt);
    const denominator = metric === "edit" ? Math.max(sourceCdr3.length, candidateCdr3.length) : sourceCdr3.length;
    const distanceLimit = Math.floor((1 - options.minimumIdentity) * denominator + 1e-9);
    if (distance > distanceLimit) return;
    const cdr3Identity = denominator ? 1 - distance / denominator : 0;
    const previous = best.get(candidateLineage);
    if (previous && (previous.cdr3Identity > cdr3Identity ||
      (previous.cdr3Identity === cdr3Identity && agreementRank(previous.callAgreement) >= agreementRank(agreement.label)))) return;
    best.set(candidateLineage, {
      lineageId: candidateLineage,
      sourceLineageId: assignments[sourceIndex],
      sourceOrdinal: source.ordinal,
      candidateOrdinal: candidate.ordinal,
      cdr3Identity,
      cdr3Distance: distance,
      cdr3Metric: metric,
      cdr3LengthDelta: candidateCdr3.length - sourceCdr3.length,
      callAgreement: agreement.label,
      uniqueMembers: uniqueMembers[candidateLineage],
      abundance: abundance[candidateLineage],
      locus: candidate.locus,
      vCalls: callSet(candidate.vCall, options.callResolution, options.ambiguity),
      jCalls: callSet(candidate.jCall, options.callResolution, options.ambiguity),
      cdr3Length: candidateCdr3.length,
      studyGroup: datasetScopeValue(candidate, options.scope ?? "global"),
    });
  };

  for (const sourceIndex of sourceRecords) {
    const source = records[sourceIndex];
    const cdr3 = normalizeNt(source.cdr3Nt);
    const hard = neighbourHardPartition(source, options);
    if (metric === "edit") {
      const partition = bkPartitions.get(hard);
      const radius = partition ? maximumEditSearchRadius(cdr3.length, partition, options.minimumIdentity) : 0;
      let candidates = 0;
      const completed = searchCdr3Bk(partition?.root ?? null, cdr3, radius, (node, distance) => {
        for (const candidateIndex of node.ordinals) {
          const candidateLineage = assignments[candidateIndex];
          if (candidateLineage <= 0 || sourceIds.has(candidateLineage)) continue;
          consider(sourceIndex, candidateIndex, distance);
          candidates += 1;
          if (candidates >= options.maxCandidateComparisons) return false;
        }
        return true;
      });
      if (!completed) truncatedSourceRecords += 1;
      continue;
    }
    const distanceLimit = Math.floor((1 - options.minimumIdentity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, distanceLimit + 1));
    const candidates = new Set<number>();
    const tokens = callPolicy === "both" ? recordIndexTokens(source, options) : ["*"];
    for (const token of tokens) {
      for (const block of blocks(cdr3, blockCount)) {
        const key = `${hard}\u0000${cdr3.length}\u0000${token}\u0000${block.index}\u0000${block.value}`;
        for (const candidateIndex of bucket.get(key) ?? []) {
          const candidateLineage = assignments[candidateIndex];
          if (candidateLineage <= 0 || sourceIds.has(candidateLineage)) continue;
          candidates.add(candidateIndex);
          if (candidates.size >= options.maxCandidateComparisons) break;
        }
        if (candidates.size >= options.maxCandidateComparisons) break;
      }
      if (candidates.size >= options.maxCandidateComparisons) break;
    }
    if (candidates.size >= options.maxCandidateComparisons) truncatedSourceRecords += 1;
    for (const candidateIndex of candidates) {
      const distance = hammingDistanceWithin(cdr3, normalizeNt(records[candidateIndex].cdr3Nt), distanceLimit, false);
      consider(sourceIndex, candidateIndex, distance);
    }
  }
  const hits = [...best.values()]
    .sort((left, right) => right.cdr3Identity - left.cdr3Identity || agreementRank(right.callAgreement) - agreementRank(left.callAgreement) || right.abundance - left.abundance || left.lineageId - right.lineageId)
    .slice(0, Math.max(1, options.maximumResults));
  return { hits, indexedRecords, sourceRecords: sourceRecords.length, candidateComparisons, truncatedSourceRecords };
}

function queryValue(record: PostAnalysisRecord, target: QueryTarget): string {
  if (target === "cdr3_nt") return normalizeNt(record.cdr3Nt);
  if (target === "cdr3_aa") return normalizeAa(record.cdr3Aa);
  return "";
}

function matchesOneQueryConstraint(record: PostAnalysisRecord, constraint: QueryConstraint | undefined, options: QueryOptions): boolean {
  if (constraint?.locus && record.locus !== constraint.locus) return false;
  if (constraint?.vCall && !callsCompatible(record.vCall, constraint.vCall, options.callResolution, options.ambiguity)) return false;
  if (constraint?.jCall && !callsCompatible(record.jCall, constraint.jCall, options.callResolution, options.ambiguity)) return false;
  return true;
}

function matchesQueryConstraints(record: PostAnalysisRecord, options: QueryOptions): boolean {
  if (options.productiveOnly && !record.productive) return false;
  return matchesOneQueryConstraint(record, options, options);
}

export function queryRecords(
  records: PostAnalysisRecord[],
  queries: string[],
  options: QueryOptions,
  packedSketches?: Uint32Array,
  activeMask?: Uint8Array,
  lineageAssignments?: Int32Array,
  lineageCount = 0,
): QueryHit[] {
  const normalizedQueries = queries.map((query) => options.target === "cdr3_aa" ? normalizeAa(query) : normalizeNt(query));
  const querySketches = options.target === "trimmed" ? normalizedQueries.map((query) => minHashSketch(query)) : [];
  const hits: QueryHit[] = [];
  const lineageMode = options.resultMode === "lineages";
  if (lineageMode && (!lineageAssignments || !(lineageCount > 0))) throw new Error("Lineage query mode requires lineage assignments on the current working set.");
  const bestScores = lineageMode ? new Float64Array(lineageCount + 1).fill(-1) : null;
  const bestDistances = lineageMode ? new Float64Array(lineageCount + 1).fill(Number.POSITIVE_INFINITY) : null;
  const bestOrdinals = lineageMode ? new Int32Array(lineageCount + 1).fill(-1) : null;
  const bestQueries = lineageMode ? new Int32Array(lineageCount + 1).fill(-1) : null;
  const matchedSequenceCounts = lineageMode ? new Uint32Array(lineageCount + 1) : null;
  const lastMatchedOrdinals = lineageMode ? new Int32Array(lineageCount + 1).fill(-1) : null;
  const matchedQueryMasks = lineageMode && normalizedQueries.length <= 32 ? new Uint32Array(lineageCount + 1) : null;
  const matchedQuerySets = lineageMode && normalizedQueries.length > 32 ? new Map<number, Set<number>>() : null;
  const bestMatched = lineageMode ? new Map<number, string>() : null;
  const accept = (hit: QueryHit) => {
    if (!lineageMode) { hits.push(hit); return; }
    const lineageId = lineageAssignments![hit.ordinal] ?? 0;
    if (!(lineageId > 0) || lineageId > lineageCount) return;
    if (lastMatchedOrdinals![lineageId] !== hit.ordinal) {
      lastMatchedOrdinals![lineageId] = hit.ordinal;
      matchedSequenceCounts![lineageId] += 1;
    }
    if (matchedQueryMasks) matchedQueryMasks[lineageId] |= (1 << hit.queryIndex) >>> 0;
    else {
      let values = matchedQuerySets!.get(lineageId);
      if (!values) { values = new Set(); matchedQuerySets!.set(lineageId, values); }
      values.add(hit.queryIndex);
    }
    if (hit.score > bestScores![lineageId] || (hit.score === bestScores![lineageId] && (hit.distance < bestDistances![lineageId] || (hit.distance === bestDistances![lineageId] && hit.ordinal < bestOrdinals![lineageId])))) {
      bestScores![lineageId] = hit.score;
      bestDistances![lineageId] = hit.distance;
      bestOrdinals![lineageId] = hit.ordinal;
      bestQueries![lineageId] = hit.queryIndex;
      bestMatched!.set(lineageId, hit.matched);
    }
  };
  for (let recordIndex = 0; recordIndex < records.length; recordIndex += 1) {
    if (activeMask && !activeMask[recordIndex]) continue;
    const record = records[recordIndex];
    if (!matchesQueryConstraints(record, options)) continue;
    for (let queryIndex = 0; queryIndex < normalizedQueries.length; queryIndex += 1) {
      if (!matchesOneQueryConstraint(record, options.queryConstraints?.[queryIndex], options)) continue;
      const query = normalizedQueries[queryIndex];
      if (!query) continue;
      if (options.target === "trimmed") {
        const sketch = record.trimmedSketch ?? packedSketches?.subarray(record.ordinal * HASH_SEEDS.length, (record.ordinal + 1) * HASH_SEEDS.length);
        if (!sketch) continue;
        const score = sketchSimilarity(sketch, querySketches[queryIndex]);
        if (score >= options.identity) accept({ ordinal: record.ordinal, queryIndex, score, distance: Math.round((1 - score) * 1000), matched: "VDJ k-mer sketch" });
        continue;
      }
      const target = queryValue(record, options.target);
      if (!target) continue;
      let distance = Number.POSITIVE_INFINITY;
      let score = 0;
      if (options.metric === "exact") {
        if (target !== query) continue;
        distance = 0;
        score = 1;
      } else if (options.metric === "substring") {
        if (!target.includes(query) && !query.includes(target)) continue;
        distance = Math.abs(target.length - query.length);
        score = Math.min(target.length, query.length) / Math.max(target.length, query.length);
      } else if (options.metric === "hamming") {
        if (target.length !== query.length) continue;
        const maximum = Math.floor((1 - options.identity) * target.length + 1e-9);
        distance = hammingDistanceWithin(target, query, maximum, options.target !== "cdr3_aa");
        if (distance > maximum) continue;
        score = 1 - distance / target.length;
      } else {
        const maximum = Math.floor((1 - options.identity) * Math.max(target.length, query.length) + 1e-9);
        distance = bandedEditDistance(target, query, maximum);
        if (distance > maximum) continue;
        score = 1 - distance / Math.max(target.length, query.length);
      }
      accept({ ordinal: record.ordinal, queryIndex, score, distance, matched: target });
    }
  }
  if (lineageMode) {
    for (let lineageId = 1; lineageId <= lineageCount; lineageId += 1) {
      if (bestOrdinals![lineageId] < 0) continue;
      const mask = matchedQueryMasks?.[lineageId] ?? 0;
      const matchedQueries = matchedQueryMasks ? popcount32(mask) : matchedQuerySets!.get(lineageId)?.size ?? 0;
      hits.push({
        ordinal: bestOrdinals![lineageId],
        queryIndex: bestQueries![lineageId],
        score: bestScores![lineageId],
        distance: bestDistances![lineageId],
        matched: bestMatched!.get(lineageId) ?? "",
        lineageId,
        matchedSequences: matchedSequenceCounts![lineageId],
        matchedQueries,
      });
    }
  }
  return hits.sort((a, b) => b.score - a.score || a.distance - b.distance || a.ordinal - b.ordinal).slice(0, options.maxResults);
}

export interface ExpansionOptions extends LineageOptions {
  maxResults: number;
}

export function expandSingleLinkage(
  records: PostAnalysisRecord[],
  seedOrdinals: number[],
  options: ExpansionOptions,
  activeMask?: Uint8Array,
): { ordinals: number[]; comparisons: number; capped: boolean } {
  const byLength = new Map<number, Map<string, number[]>>();
  for (let index = 0; index < records.length; index += 1) {
    if (activeMask && !activeMask[index]) continue;
    const record = records[index];
    const cdr3 = normalizeNt(record.cdr3Nt);
    if (!cdr3 || (options.productiveOnly && !record.productive)) continue;
    const maximum = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, maximum + 1));
    let indexMap = byLength.get(cdr3.length);
    if (!indexMap) {
      indexMap = new Map();
      byLength.set(cdr3.length, indexMap);
    }
    for (const block of blocks(cdr3, blockCount)) {
      const key = `${block.index}\u0000${block.value}`;
      const values = indexMap.get(key);
      if (values) values.push(index);
      else indexMap.set(key, [index]);
    }
  }

  const visited = new Set(seedOrdinals.filter((ordinal) => ordinal >= 0 && ordinal < records.length && (!activeMask || Boolean(activeMask[ordinal]))));
  const queue = [...visited];
  let comparisons = 0;
  for (let head = 0; head < queue.length && visited.size < options.maxResults; head += 1) {
    const current = queue[head];
    const source = records[current];
    const cdr3 = normalizeNt(source.cdr3Nt);
    const maximum = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
    const blockCount = Math.max(1, Math.min(cdr3.length, maximum + 1));
    const candidates = new Set<number>();
    const indexMap = byLength.get(cdr3.length);
    for (const block of blocks(cdr3, blockCount)) {
      for (const candidate of indexMap?.get(`${block.index}\u0000${block.value}`) ?? []) candidates.add(candidate);
    }
    for (const candidate of candidates) {
      if (visited.has(candidate)) continue;
      comparisons += 1;
      if (!compatibleRecords(source, records[candidate], options)) continue;
      if (hammingDistanceWithin(cdr3, normalizeNt(records[candidate].cdr3Nt), maximum, false) <= maximum) {
        visited.add(candidate);
        queue.push(candidate);
        if (visited.size >= options.maxResults) break;
      }
    }
  }
  return { ordinals: [...visited], comparisons, capped: visited.size >= options.maxResults };
}

export interface FastaRecord {
  name: string;
  sequence: string;
}

export function parseFasta(text: string, aligned = false): FastaRecord[] {
  const records: FastaRecord[] = [];
  let name = "";
  let sequence: string[] = [];
  const commit = () => {
    if (!name) return;
    const joined = sequence.join("").toUpperCase().replaceAll(".", "-").replace(/\s/g, "");
    const normalized = aligned ? joined.replace(/[^ACGTNRYKMSWBDHV-]/g, "N") : normalizeNt(joined);
    records.push({ name, sequence: normalized });
  };
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      commit();
      name = line.slice(1).trim().split(/\s+/, 1)[0] ?? "";
      sequence = [];
    } else if (name) {
      sequence.push(line);
    }
  }
  commit();
  return records;
}

export interface ReferenceMsa {
  names: string[];
  sequences: string[];
  lookup: Map<string, number>;
  length: number;
}

function nameAliases(name: string): string[] {
  const fields = name.split("|").filter(Boolean);
  const gene = fields.find((field) => /(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)[VJ]/i.test(field));
  return [...new Set([name, gene ?? "", name.split(/\s+/, 1)[0] ?? ""].filter(Boolean))];
}

export function prepareReferenceMsa(text: string): ReferenceMsa {
  const records = parseFasta(text, true);
  if (records.length < 2) throw new Error("CHMMAIRRa requires at least two aligned reference sequences.");
  const length = records[0].sequence.length;
  if (!length || records.some((record) => record.sequence.length !== length)) {
    throw new Error("Reference MSA records must all have the same aligned length.");
  }
  const lookup = new Map<string, number>();
  records.forEach((record, index) => nameAliases(record.name).forEach((alias) => lookup.set(alias, index)));
  return { names: records.map((record) => record.name), sequences: records.map((record) => record.sequence), lookup, length };
}

function topCallName(value: string): string {
  return value.split(",", 1)[0]?.trim() ?? "";
}

function chmmReferenceIndex(callName: string, msa: ReferenceMsa): number | undefined {
  let referenceIndex = msa.lookup.get(callName);
  if (referenceIndex === undefined) {
    const gene = callName.replace(/\*.*$/, "");
    const match = [...msa.lookup.entries()].find(([name]) => name.replace(/\*.*$/, "") === gene);
    referenceIndex = match?.[1];
  }
  return referenceIndex;
}

const CHMM_COMPLEMENT: Record<string, string> = { A: "T", C: "G", G: "C", T: "A" };

/**
 * The read bases beside a segment's local alignment on its junction side (after the V
 * alignment, before the J alignment), in the alignment's orientation. `start`/`end` are the
 * AIRR 1-based `{v,j}_sequence_start`/`_end`, which refer to the reverse complement of
 * `sequence` when `rev_comp` is true.
 */
export function chmmJunctionFlank(
  segment: "V" | "J",
  sequence: string | undefined,
  revComp: string | undefined,
  start: string | number | undefined,
  end: string | number | undefined,
  limit = 64,
): string {
  const raw = (sequence ?? "").toUpperCase();
  if (!raw) return "";
  const oriented = /^(T|TRUE)$/i.test(String(revComp ?? "").trim())
    ? [...raw].reverse().map((base) => CHMM_COMPLEMENT[base] ?? "N").join("")
    : raw;
  if (segment === "V") {
    const last = Number(end);
    return Number.isInteger(last) && last > 0 ? oriented.slice(last, last + limit) : "";
  }
  const first = Number(start);
  return Number.isInteger(first) && first > 1 ? oriented.slice(Math.max(0, first - 1 - limit), first - 1) : "";
}

function bestExtension(read: string, germline: string): number {
  let score = 0;
  let best = 0;
  let length = 0;
  for (let index = 0; index < Math.min(read.length, germline.length); index += 1) {
    score += read[index] === germline[index] && /[ACGT]/.test(read[index]) ? 1 : -1;
    if (score > best) {
      best = score;
      length = index + 1;
    }
  }
  return length;
}

/**
 * Extend a V alignment toward its germline 3' end (a J alignment toward its 5' start) by the
 * best-scoring gapless prefix of the junction-side read bases, scored +1/-1. The aligner may
 * stop a V (start a J) before trailing (leading) mismatches it attributes to the junction;
 * a PCR chimera's last informative columns can lie there. Unchanged when the reference or
 * the local germline cannot be located (threading then reports why).
 */
export function extendChmmAlignment(
  segment: "V" | "J",
  sequenceAlignment: string,
  germlineAlignment: string,
  call: string,
  flank: string,
  msa: ReferenceMsa,
): { sequenceAlignment: string; germlineAlignment: string } {
  const unchanged = { sequenceAlignment, germlineAlignment };
  if (!flank) return unchanged;
  const referenceIndex = chmmReferenceIndex(topCallName(call), msa);
  if (referenceIndex === undefined) return unchanged;
  const reference = msa.sequences[referenceIndex].replaceAll("-", "");
  const local = germlineAlignment.toUpperCase().replaceAll(".", "-").replaceAll("-", "");
  const at = reference.indexOf(local);
  if (at < 0) return unchanged;
  const read = flank.toUpperCase();
  if (segment === "V") {
    const germline = reference.slice(at + local.length);
    const length = bestExtension(read, germline);
    return { sequenceAlignment: sequenceAlignment + read.slice(0, length), germlineAlignment: germlineAlignment + germline.slice(0, length) };
  }
  const reverse = (value: string) => [...value].reverse().join("");
  const germline = reverse(reference.slice(0, at));
  const readReversed = reverse(read);
  const length = bestExtension(readReversed, germline);
  return {
    sequenceAlignment: reverse(readReversed.slice(0, length)) + sequenceAlignment,
    germlineAlignment: reverse(germline.slice(0, length)) + germlineAlignment,
  };
}

export function threadSequenceToMsa(
  sequenceAlignment: string,
  germlineAlignment: string,
  call: string,
  msa: ReferenceMsa,
): string {
  const callName = topCallName(call);
  const referenceIndex = chmmReferenceIndex(callName, msa);
  if (referenceIndex === undefined) throw new Error(`Reference MSA is missing ${callName || "the assigned allele"}.`);
  const gappedFullReference = msa.sequences[referenceIndex];
  const degappedReference = gappedFullReference.replaceAll("-", "");
  const localGermline = germlineAlignment.toUpperCase().replaceAll(".", "-");
  // CHMMera maps every non-ACGT observation other than a gap to its
  // non-informative N state. Apply that conversion before threading so IUPAC
  // query symbols cannot become accidental informative matches.
  const localQuery = sequenceAlignment.toUpperCase().replaceAll(".", "-").replace(/[^ACGT-]/g, "N");
  const degappedLocal = localGermline.replaceAll("-", "");
  const matchStart = degappedReference.indexOf(degappedLocal);
  if (matchStart < 0) throw new Error(`The ${callName} local germline alignment cannot be located in its MSA reference.`);
  const withoutInsertions: string[] = [];
  for (let index = 0; index < localGermline.length; index += 1) {
    if (localGermline[index] !== "-") withoutInsertions.push(localQuery[index] ?? "N");
  }
  const query = [..."N".repeat(matchStart), ...withoutInsertions];
  const threaded = new Array<string>(msa.length).fill("-");
  let position = 0;
  for (let index = 0; index < msa.length; index += 1) {
    if (position >= query.length) break;
    if (gappedFullReference[index] !== "-") threaded[index] = query[position++] ?? "N";
  }
  return threaded.join("");
}

export interface ChmmOptions {
  method: "BW" | "DB";
  priorProbability: number;
  baseMutationProbability: number;
  mutationRates: number[];
  mutationSwitchProbability: number;
  detailed: boolean;
  tracePath?: boolean;
}

export interface ChmmResult {
  probability: number;
  dfr: number;
  startingReference: string;
  recombinations: Array<{ position: number; left: string; right: string }>;
  referencePath?: Int32Array;
}

export function chmmairraDistanceFromReference(query: string, germline: string): number {
  const normalizedQuery = query.toUpperCase().replaceAll(".", "-");
  const normalizedGermline = germline.toUpperCase().replaceAll(".", "-");
  let mismatches = 0;
  // This intentionally follows CHMMAIRRa.jl's add_DFR_column!: zip the local
  // alignments and count every unequal character, including gap and N sites.
  for (let index = 0; index < Math.min(normalizedQuery.length, normalizedGermline.length); index += 1) {
    if (normalizedQuery[index] !== normalizedGermline[index]) mismatches += 1;
  }
  return mismatches;
}

function emission(reference: string, observation: string, mutation: number): number {
  if (observation === "-" || observation === "N") return 1;
  return reference === observation ? 1 - mutation : mutation / 3;
}

function approximateEmissions(msa: ReferenceMsa, observation: string, mutations: Float64Array): Float64Array {
  const values = new Float64Array(msa.sequences.length * msa.length);
  for (let state = 0; state < msa.sequences.length; state += 1) {
    for (let site = 0; site < msa.length; site += 1) {
      values[state * msa.length + site] = emission(msa.sequences[state][site], observation[site], mutations[state]);
    }
  }
  return values;
}

export function estimateMutationProbabilities(msa: ReferenceMsa, observation: string, prior: number, baseMutation: number): Float64Array {
  const states = msa.sequences.length;
  const length = msa.length;
  const mutations = new Float64Array(states);
  mutations.fill(baseMutation);
  if (states < 2 || length < 1) return mutations;
  const emissions = approximateEmissions(msa, observation, mutations);
  const alpha = new Float64Array(states * length);
  const beta = new Float64Array(states * length);
  const scale = new Float64Array(length);
  scale[0] = 1;
  const switchProbability = prior / length;
  const same = 1 - switchProbability;
  const different = switchProbability / (states - 1);
  for (let state = 0; state < states; state += 1) alpha[state * length] = emissions[state * length] / states;
  for (let site = 0; site < length - 1; site += 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += alpha[state * length + site];
    let nextTotal = 0;
    for (let state = 0; state < states; state += 1) {
      const value = ((total - alpha[state * length + site]) * different + alpha[state * length + site] * same) * emissions[state * length + site + 1];
      alpha[state * length + site + 1] = value;
      nextTotal += value;
    }
    scale[site + 1] = nextTotal ? 1 / nextTotal : 1;
    for (let state = 0; state < states; state += 1) alpha[state * length + site + 1] *= scale[site + 1];
  }
  for (let state = 0; state < states; state += 1) beta[state * length + length - 1] = scale[length - 1];
  for (let site = length - 2; site >= 0; site -= 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += beta[state * length + site + 1] * emissions[state * length + site + 1];
    for (let state = 0; state < states; state += 1) {
      const own = beta[state * length + site + 1] * emissions[state * length + site + 1];
      beta[state * length + site] = (different * (total - own) + same * own) * scale[site];
    }
  }
  const mutated = new Float64Array(states);
  const sameCount = new Float64Array(states);
  mutated.fill(2);
  sameCount.fill(10);
  for (let site = 0; site < length; site += 1) {
    if (observation[site] === "-" || observation[site] === "N") continue;
    let normalization = 0;
    for (let state = 0; state < states; state += 1) normalization += alpha[state * length + site] * beta[state * length + site];
    if (!normalization) continue;
    for (let state = 0; state < states; state += 1) {
      const posterior = alpha[state * length + site] * beta[state * length + site] / normalization;
      if (observation[site] === msa.sequences[state][site]) sameCount[state] += posterior;
      else mutated[state] += posterior;
    }
  }
  for (let state = 0; state < states; state += 1) mutations[state] = mutated[state] / (mutated[state] + sameCount[state]);
  return mutations;
}

function approximateProbability(msa: ReferenceMsa, observation: string, prior: number, mutations: Float64Array): number {
  const states = msa.sequences.length;
  if (states < 2) return 0;
  const emissions = approximateEmissions(msa, observation, mutations);
  const nonChimeric = new Float64Array(states);
  const chimeric = new Float64Array(states);
  const switchProbability = prior / msa.length;
  const same = 1 - switchProbability;
  const different = switchProbability / (states - 1);
  for (let state = 0; state < states; state += 1) nonChimeric[state] = emissions[state * msa.length] / states;
  for (let site = 0; site < msa.length - 1; site += 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
    let nextTotal = 0;
    for (let state = 0; state < states; state += 1) {
      const emissionValue = emissions[state * msa.length + site + 1];
      const own = nonChimeric[state] + chimeric[state];
      chimeric[state] = ((total - own) * different + chimeric[state] * same) * emissionValue;
      nonChimeric[state] = nonChimeric[state] * same * emissionValue;
      nextTotal += chimeric[state] + nonChimeric[state];
    }
    const scaling = nextTotal ? 1 / nextTotal : 1;
    for (let state = 0; state < states; state += 1) {
      chimeric[state] *= scaling;
      nonChimeric[state] *= scaling;
    }
  }
  let chimericTotal = 0;
  let total = 0;
  for (let state = 0; state < states; state += 1) {
    chimericTotal += chimeric[state];
    total += chimeric[state] + nonChimeric[state];
  }
  return total ? chimericTotal / total : 0;
}

function fullProbability(msa: ReferenceMsa, observation: string, options: ChmmOptions): number {
  const references = msa.sequences.length;
  const rates = options.mutationRates.length ? options.mutationRates : [0.005];
  const rateCount = rates.length;
  const states = references * rateCount;
  if (references < 2) return 0;
  const nonChimeric = new Float64Array(states);
  const chimeric = new Float64Array(states);
  const nextNonChimeric = new Float64Array(states);
  const nextChimeric = new Float64Array(states);
  const switchProbability = options.priorProbability / msa.length;
  const self = 1 - switchProbability - options.mutationSwitchProbability;
  const differentReference = switchProbability / ((references - 1) * rateCount);
  const differentMutation = rateCount === 1 ? 0 : options.mutationSwitchProbability / (rateCount - 1);
  const emissionAt = (state: number, site: number) => {
    const reference = Math.floor(state / rateCount);
    const rate = state % rateCount;
    return emission(msa.sequences[reference][site], observation[site], rates[rate]);
  };
  for (let state = 0; state < states; state += 1) nonChimeric[state] = emissionAt(state, 0) / states;
  for (let site = 0; site < msa.length - 1; site += 1) {
    let total = 0;
    for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
    let nextTotal = 0;
    for (let reference = 0; reference < references; reference += 1) {
      let nonReference = 0;
      let chimReference = 0;
      for (let rate = 0; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        nonReference += nonChimeric[state];
        chimReference += chimeric[state];
      }
      for (let rate = 0; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        const emit = emissionAt(state, site + 1);
        nextNonChimeric[state] = (nonChimeric[state] * self + (nonReference - nonChimeric[state]) * differentMutation) * emit;
        nextChimeric[state] = (chimeric[state] * self + (chimReference - chimeric[state]) * differentMutation + (total - chimReference - nonReference) * differentReference) * emit;
        nextTotal += nextNonChimeric[state] + nextChimeric[state];
      }
    }
    const scaling = nextTotal ? 1 / nextTotal : 1;
    for (let state = 0; state < states; state += 1) {
      nonChimeric[state] = nextNonChimeric[state] * scaling;
      chimeric[state] = nextChimeric[state] * scaling;
    }
  }
  let chimericTotal = 0;
  let total = 0;
  for (let state = 0; state < states; state += 1) {
    chimericTotal += chimeric[state];
    total += chimeric[state] + nonChimeric[state];
  }
  return total ? chimericTotal / total : 0;
}

function fullViterbi(msa: ReferenceMsa, observation: string, options: ChmmOptions): Pick<ChmmResult, "startingReference" | "recombinations" | "referencePath"> {
  const references = msa.sequences.length;
  const rates = options.mutationRates.length ? options.mutationRates : [0.005];
  const rateCount = rates.length;
  const states = references * rateCount;
  const length = msa.length;
  if (references < 2) return {
    startingReference: msa.names[0] ?? "",
    recombinations: [],
    referencePath: options.tracePath ? new Int32Array(length) : undefined,
  };
  const previous = new Float64Array(states);
  const next = new Float64Array(states);
  const from = new Int32Array(states * length);
  const switchProbability = options.priorProbability / length;
  const logSelf = Math.log(1 - switchProbability - options.mutationSwitchProbability);
  const logDifferentReference = Math.log(switchProbability / ((references - 1) * rateCount));
  const logDifferentMutation = rateCount === 1 ? Number.NEGATIVE_INFINITY : Math.log(options.mutationSwitchProbability / (rateCount - 1));
  const emissionAt = (state: number, site: number) => emission(msa.sequences[Math.floor(state / rateCount)][site], observation[site], rates[state % rateCount]);
  for (let state = 0; state < states; state += 1) {
    previous[state] = Math.log(emissionAt(state, 0) / states);
    from[state * length] = state;
  }
  for (let site = 0; site < length - 1; site += 1) {
    const bestInReference = new Int32Array(references);
    for (let reference = 0; reference < references; reference += 1) {
      let best = reference * rateCount;
      for (let rate = 1; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        if (previous[state] > previous[best]) best = state;
      }
      bestInReference[reference] = best;
    }
    let bestReference = 0;
    let secondReference = references > 1 ? 1 : 0;
    if (previous[bestInReference[secondReference]] > previous[bestInReference[bestReference]]) [bestReference, secondReference] = [secondReference, bestReference];
    for (let reference = 2; reference < references; reference += 1) {
      if (previous[bestInReference[reference]] > previous[bestInReference[bestReference]]) {
        secondReference = bestReference;
        bestReference = reference;
      } else if (previous[bestInReference[reference]] > previous[bestInReference[secondReference]]) secondReference = reference;
    }
    for (let reference = 0; reference < references; reference += 1) {
      const outside = bestInReference[reference === bestReference ? secondReference : bestReference];
      for (let rate = 0; rate < rateCount; rate += 1) {
        const state = reference * rateCount + rate;
        let selected = state;
        let score = previous[state] + logSelf;
        const referenceScore = previous[outside] + logDifferentReference;
        if (referenceScore > score) {
          selected = outside;
          score = referenceScore;
        }
        if (rateCount > 1) {
          let otherRate = reference * rateCount + (rate === 0 ? 1 : 0);
          for (let candidateRate = 0; candidateRate < rateCount; candidateRate += 1) {
            const candidate = reference * rateCount + candidateRate;
            if (candidate !== state && previous[candidate] > previous[otherRate]) otherRate = candidate;
          }
          const mutationScore = previous[otherRate] + logDifferentMutation;
          if (mutationScore > score) {
            selected = otherRate;
            score = mutationScore;
          }
        }
        from[state * length + site + 1] = selected;
        next[state] = score + Math.log(emissionAt(state, site + 1));
      }
    }
    previous.set(next);
  }
  let current = 0;
  for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
  const recombinations: Array<{ position: number; left: string; right: string }> = [];
  const referencePath = options.tracePath ? new Int32Array(length) : undefined;
  if (referencePath) referencePath[length - 1] = Math.floor(current / rateCount);
  for (let site = length - 1; site >= 1; site -= 1) {
    const parent = from[current * length + site];
    const leftReference = Math.floor(parent / rateCount);
    const rightReference = Math.floor(current / rateCount);
    if (leftReference !== rightReference) recombinations.push({ position: site + 1, left: msa.names[leftReference], right: msa.names[rightReference] });
    current = parent;
    if (referencePath) referencePath[site - 1] = Math.floor(current / rateCount);
  }
  return { startingReference: msa.names[Math.floor(current / rateCount)], recombinations: recombinations.reverse(), referencePath };
}

function approximateViterbi(msa: ReferenceMsa, observation: string, prior: number, mutations: Float64Array, tracePath = false): Pick<ChmmResult, "startingReference" | "recombinations" | "referencePath"> {
  const states = msa.sequences.length;
  const length = msa.length;
  const emissions = approximateEmissions(msa, observation, mutations);
  const previous = new Float64Array(states);
  const next = new Float64Array(states);
  const from = new Int32Array(states * length);
  const same = Math.log(1 - prior / length);
  const different = Math.log((prior / length) / (states - 1));
  for (let state = 0; state < states; state += 1) {
    previous[state] = Math.log(emissions[state * length] / states);
    from[state * length] = state;
  }
  for (let site = 0; site < length - 1; site += 1) {
    let best = 0;
    for (let state = 1; state < states; state += 1) if (previous[state] > previous[best]) best = state;
    for (let state = 0; state < states; state += 1) {
      const stay = previous[state] + same;
      const change = previous[best] + different;
      const remain = state === best || stay > change;
      from[state * length + site + 1] = remain ? state : best;
      next[state] = (remain ? stay : change) + Math.log(emissions[state * length + site + 1]);
    }
    previous.set(next);
  }
  let current = 0;
  for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
  const reverse: Array<{ position: number; left: string; right: string }> = [];
  const referencePath = tracePath ? new Int32Array(length) : undefined;
  if (referencePath) referencePath[length - 1] = current;
  for (let site = length - 1; site >= 1; site -= 1) {
    const parent = from[current * length + site];
    if (parent !== current) reverse.push({ position: site + 1, left: msa.names[parent], right: msa.names[current] });
    current = parent;
    if (referencePath) referencePath[site - 1] = current;
  }
  return { startingReference: msa.names[current], recombinations: reverse.reverse(), referencePath };
}

export function runChmm(
  msa: ReferenceMsa,
  threadedObservation: string,
  localSequenceAlignment: string,
  localGermlineAlignment: string,
  options: ChmmOptions,
): ChmmResult {
  if (threadedObservation.length !== msa.length) throw new Error("Threaded query and reference MSA lengths differ.");
  const dfr = chmmairraDistanceFromReference(localSequenceAlignment, localGermlineAlignment);
  if (options.method === "DB") {
    const path = options.detailed ? fullViterbi(msa, threadedObservation, options) : { startingReference: "", recombinations: [] };
    return {
      probability: fullProbability(msa, threadedObservation, options),
      dfr,
      ...path,
    };
  }
  const mutations = estimateMutationProbabilities(msa, threadedObservation, options.priorProbability, options.baseMutationProbability);
  const probability = approximateProbability(msa, threadedObservation, options.priorProbability, mutations);
  const path = options.detailed ? approximateViterbi(msa, threadedObservation, options.priorProbability, mutations, options.tracePath) : { startingReference: "", recombinations: [] };
  return { probability, dfr, ...path };
}
