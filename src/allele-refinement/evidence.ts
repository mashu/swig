import { datasetScopeValue } from "../study-design.ts";
import type {
  AlleleRefinementOptions,
  ReferenceAlleleGraph,
  RefinementInputRow,
  SparseEvidenceMatrix,
} from "./types.ts";

interface AlternativeEvidence {
  call: string;
  score: number | null;
}

function callTokens(value: string): string[] {
  return [...new Set(value.split(",").map((call) => call.trim()).filter(Boolean))];
}

export function parseAlternativeEvidence(value: string): AlternativeEvidence[] {
  return value.split(";").map((entry) => {
    const fields = entry.split("|");
    const score = Number(fields[1]);
    return { call: fields[0]?.trim() ?? "", score: Number.isFinite(score) ? score : null };
  }).filter((entry) => entry.call);
}

function modelScopeValue(row: RefinementInputRow, options: AlleleRefinementOptions): string {
  return datasetScopeValue({
    datasetId: row.datasetId,
    sampleId: row.sampleId,
    subjectId: row.subjectId,
  }, options.scope);
}

export interface SparseEvidenceRow {
  ordinal: number;
  groupKey: string;
  entries: Array<{ node: number; weight: number }>;
  weight: number;
  localTop: number;
  localTopProbability: number;
  truncated: boolean;
}

/**
 * Per-SNP evidence odds assigned to an unreported reference neighbour.
 *
 * If the true and reported alleles differ at one nucleotide, a read-level
 * substitution rate mu gives an equal-substitution likelihood ratio of
 * mu / [3(1-mu)] for changing the true diagnostic base into the reported
 * base rather than retaining it. The baseline term deliberately represents
 * unmodelled assignment uncertainty, not a sequencing-error probability.
 */
export function adaptiveNeighbourOdds(
  readShm: number | null,
  options: AlleleRefinementOptions,
): number {
  const baseline = Math.max(0, Math.min(0.999999, options.baselineNeighbourOdds));
  const maximumShm = Math.max(0, Math.min(0.95, options.maximumShm));
  const mu = Math.max(0, Math.min(maximumShm, readShm ?? 0));
  const sensitivity = Math.max(0, options.shmLeakageSensitivity);
  const somaticOdds = mu > 0 ? sensitivity * mu / (3 * (1 - mu)) : 0;
  const cap = Math.max(0, Math.min(0.999999, options.maximumNeighbourOdds));
  return Math.min(cap, baseline + somaticOdds);
}

/**
 * Converts a literal aligner call into a sparse, database-aware evidence row.
 * Co-optimal calls start at exactly equal weight. Unreported nucleotide
 * neighbours receive geometric leakage; this is an explicit evidence kernel,
 * not a claim that the affine alignment score is a calibrated read likelihood.
 */
export function buildSparseEvidenceRow(
  row: RefinementInputRow,
  graph: ReferenceAlleleGraph,
  options: AlleleRefinementOptions,
): SparseEvidenceRow | null {
  const direct = new Map<number, number>();
  const selected = callTokens(row.call);
  for (const call of selected) {
    const node = graph.callToNode.get(call);
    if (node !== undefined) direct.set(node, 1);
  }
  const temperature = Math.max(1e-6, options.alternativeScoreTemperature);
  for (const alternative of parseAlternativeEvidence(row.alternatives)) {
    const node = graph.callToNode.get(alternative.call);
    if (node === undefined || direct.has(node)) continue;
    const relative = alternative.score !== null && row.score !== null
      ? Math.min(1, Math.exp((alternative.score - row.score) / temperature))
      : options.unscoredAlternativeWeight;
    direct.set(node, Math.max(direct.get(node) ?? 0, relative));
  }
  if (!direct.size) return null;

  const expanded = new Map(direct);
  const substitutionOdds = adaptiveNeighbourOdds(row.shm, options);
  const baselineOdds = Math.max(0, Math.min(options.maximumNeighbourOdds, options.baselineNeighbourOdds));
  if ((substitutionOdds > 0 || baselineOdds > 0) && options.neighbourRadius > 0) {
    for (const [source, sourceWeight] of direct) {
      for (const neighbour of graph.neighbours[source]) {
        if (neighbour.distance > options.neighbourRadius) break;
        const perEditOdds = neighbour.substitutionOnly ? substitutionOdds : baselineOdds;
        const leaked = sourceWeight * perEditOdds ** neighbour.distance;
        if (leaked > (expanded.get(neighbour.index) ?? 0)) expanded.set(neighbour.index, leaked);
      }
    }
  }
  let entries = [...expanded.entries()].map(([node, weight]) => ({ node, weight }))
    .filter((entry) => entry.weight > 0)
    .sort((left, right) => right.weight - left.weight || left.node - right.node);
  const retainedCap = Math.max(1, Math.floor(options.maxCandidatesPerRow), direct.size);
  const truncated = entries.length > retainedCap;
  if (truncated) entries = entries.slice(0, retainedCap);
  const sum = entries.reduce((total, entry) => total + entry.weight, 0);
  if (!(sum > 0)) return null;
  entries.forEach((entry) => { entry.weight /= sum; });
  const top = entries[0];
  const scopeValue = modelScopeValue(row, options);
  const groupKey = `${scopeValue}\u0000${row.locus}\u0000${graph.segment}`;
  return {
    ordinal: row.ordinal,
    groupKey,
    entries,
    weight: options.weighting === "abundance" ? Math.max(1, row.abundance) : 1,
    localTop: top.node,
    localTopProbability: top.weight,
    truncated,
  };
}

export function sparseEvidenceMatrix(rows: readonly SparseEvidenceRow[]): SparseEvidenceMatrix {
  const groups = [...new Set(rows.map((row) => row.groupKey))].sort();
  const groupIndex = new Map(groups.map((key, index) => [key, index] as const));
  const rowOffsets = new Uint32Array(rows.length + 1);
  let nonZeros = 0;
  rows.forEach((row, index) => {
    rowOffsets[index] = nonZeros;
    nonZeros += row.entries.length;
  });
  rowOffsets[rows.length] = nonZeros;
  const columns = new Uint32Array(nonZeros);
  const logEvidence = new Float32Array(nonZeros);
  const ordinals = new Uint32Array(rows.length);
  const weights = new Float32Array(rows.length);
  const rowGroups = new Uint32Array(rows.length);
  const localTop = new Int32Array(rows.length);
  const localTopProbability = new Float32Array(rows.length);
  let offset = 0;
  rows.forEach((row, index) => {
    ordinals[index] = row.ordinal;
    weights[index] = row.weight;
    rowGroups[index] = groupIndex.get(row.groupKey) ?? 0;
    localTop[index] = row.localTop;
    localTopProbability[index] = row.localTopProbability;
    for (const entry of row.entries) {
      columns[offset] = entry.node;
      logEvidence[offset] = Math.log(Math.max(Number.MIN_VALUE, entry.weight));
      offset += 1;
    }
  });
  return {
    rowOffsets,
    columns,
    logEvidence,
    ordinals,
    weights,
    groupKeys: groups,
    rowGroups,
    localTop,
    localTopProbability,
    skippedRows: 0,
    truncatedRows: rows.reduce((total, row) => total + Number(row.truncated), 0),
  };
}

/** Fixed-row/dynamically-grown-NNZ builder used by the browser worker. */
export class SparseEvidenceAccumulator {
  private readonly graph: ReferenceAlleleGraph;
  private readonly options: AlleleRefinementOptions;
  private readonly rowOffsets: Uint32Array;
  private readonly ordinals: Uint32Array;
  private readonly weights: Float32Array;
  private readonly rowGroups: Uint32Array;
  private readonly localTop: Int32Array;
  private readonly localTopProbability: Float32Array;
  private columns = new Uint32Array(4096);
  private logEvidence = new Float32Array(4096);
  private rowCount = 0;
  private nonZeros = 0;
  private readonly groupIndex = new Map<string, number>();
  private readonly groupKeys: string[] = [];
  skippedRows = 0;
  truncatedRows = 0;

  constructor(
    maximumRows: number,
    graph: ReferenceAlleleGraph,
    options: AlleleRefinementOptions,
  ) {
    this.graph = graph;
    this.options = options;
    this.rowOffsets = new Uint32Array(maximumRows + 1);
    this.ordinals = new Uint32Array(maximumRows);
    this.weights = new Float32Array(maximumRows);
    this.rowGroups = new Uint32Array(maximumRows);
    this.localTop = new Int32Array(maximumRows);
    this.localTop.fill(-1);
    this.localTopProbability = new Float32Array(maximumRows);
  }

  private reserve(required: number) {
    if (required <= this.columns.length) return;
    let capacity = this.columns.length;
    while (capacity < required) capacity = Math.max(capacity + 4096, Math.ceil(capacity * 1.6));
    const columns = new Uint32Array(capacity);
    const evidence = new Float32Array(capacity);
    columns.set(this.columns);
    evidence.set(this.logEvidence);
    this.columns = columns;
    this.logEvidence = evidence;
  }

  add(input: RefinementInputRow): boolean {
    const row = buildSparseEvidenceRow(input, this.graph, this.options);
    if (!row) {
      this.skippedRows += 1;
      return false;
    }
    if (this.rowCount >= this.ordinals.length) throw new Error("Allele-refinement input exceeded the declared AIRR record count.");
    this.reserve(this.nonZeros + row.entries.length);
    let group = this.groupIndex.get(row.groupKey);
    if (group === undefined) {
      group = this.groupKeys.length;
      this.groupKeys.push(row.groupKey);
      this.groupIndex.set(row.groupKey, group);
    }
    this.rowOffsets[this.rowCount] = this.nonZeros;
    this.ordinals[this.rowCount] = row.ordinal;
    this.weights[this.rowCount] = row.weight;
    this.rowGroups[this.rowCount] = group;
    this.localTop[this.rowCount] = row.localTop;
    this.localTopProbability[this.rowCount] = row.localTopProbability;
    for (const entry of row.entries) {
      this.columns[this.nonZeros] = entry.node;
      this.logEvidence[this.nonZeros] = Math.log(Math.max(Number.MIN_VALUE, entry.weight));
      this.nonZeros += 1;
    }
    this.rowCount += 1;
    this.rowOffsets[this.rowCount] = this.nonZeros;
    if (row.truncated) this.truncatedRows += 1;
    return true;
  }

  finish(): SparseEvidenceMatrix {
    return {
      // Views avoid temporarily duplicating the complete sparse matrix before
      // fitting. The accumulator is released by the worker immediately after
      // this handoff.
      rowOffsets: this.rowOffsets.subarray(0, this.rowCount + 1),
      columns: this.columns.subarray(0, this.nonZeros),
      logEvidence: this.logEvidence.subarray(0, this.nonZeros),
      ordinals: this.ordinals.subarray(0, this.rowCount),
      weights: this.weights.subarray(0, this.rowCount),
      groupKeys: [...this.groupKeys],
      rowGroups: this.rowGroups.subarray(0, this.rowCount),
      localTop: this.localTop.subarray(0, this.rowCount),
      localTopProbability: this.localTopProbability.subarray(0, this.rowCount),
      skippedRows: this.skippedRows,
      truncatedRows: this.truncatedRows,
    };
  }
}
