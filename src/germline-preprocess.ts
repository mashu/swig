export type GermlineSegment = "V" | "D" | "J" | "C";
export type GermlineLocus = "IGH" | "IGK" | "IGL" | "TRA" | "TRB" | "TRD" | "TRG";

// frame, J CDR3 stop, FWR1/CDR1/FWR2/CDR2/FWR3 start/end pairs, source.
export type CompactMetadata = [
  frame: number,
  cdr3Stop: number,
  fwr1Start: number,
  fwr1End: number,
  cdr1Start: number,
  cdr1End: number,
  fwr2Start: number,
  fwr2End: number,
  cdr2Start: number,
  cdr2End: number,
  fwr3Start: number,
  fwr3End: number,
  source: number,
];

export type MetadataAllele = [name: string, sequence: string, metadata?: CompactMetadata];

export type GermlineMatchMode = "strict" | "permissive" | "best_guess";

export interface GermlineMatchOptions {
  /** Named presets; explicit thresholds below override the selected preset. */
  mode?: GermlineMatchMode;
  vSameGeneMinIdentity?: number;
  vNearestMinIdentity?: number;
  jSameGeneMinIdentity?: number;
  jNearestMinIdentity?: number;
  /** Number of k-mer-ranked non-gene candidates to align after named candidates fail. */
  nearestCandidates?: number;
  /** Include one structured diagnostic record per input allele. */
  includeDiagnostics?: boolean;
  /** Internal provenance index used by the progressively broadened tier search. */
  taxonomicTier?: number;
}

export interface ResolvedGermlineMatchOptions {
  mode: GermlineMatchMode;
  vSameGeneMinIdentity: number;
  vNearestMinIdentity: number;
  jSameGeneMinIdentity: number;
  jNearestMinIdentity: number;
  nearestCandidates: number;
  includeDiagnostics: boolean;
}

export type GermlineDiagnosticStatus =
  | "retained"
  | "imgt"
  | "transferred"
  | "motif"
  | "normalized"
  | "unresolved";

export interface GermlineRecordDiagnostic {
  segment: GermlineSegment;
  name: string;
  locus: GermlineLocus;
  status: GermlineDiagnosticStatus;
  source: string;
  template?: string;
  identity?: number;
  matchKind?: "same_allele" | "same_gene" | "nearest";
  taxonomicTier?: number;
  attemptedCandidates: number;
  bestCandidate?: string;
  bestIdentity?: number;
  rejectionCounts?: Record<string, number>;
}

export interface GermlinePreprocessReport {
  fasta: string;
  count: number;
  annotated: number;
  unannotated: number;
  exactImgt: number;
  transferred: number;
  motifValidated: number;
  ambiguousBases: number;
  loci: GermlineLocus[];
  warnings: string[];
  diagnostics?: GermlineRecordDiagnostic[];
}

export interface GermlineNormalizationReport {
  fasta: string;
  count: number;
  ambiguousBases: number;
  loci: GermlineLocus[];
}

export interface IgblastAnnotationApplication {
  fasta: string;
  matched: number;
  annotated: number;
  total: number;
  unmatched: string[];
  /** IgBLAST .aux column five, keyed by exact J identifier. */
  fwr4EndOffsets?: Record<string, number>;
}

interface ParsedRecord {
  header: string;
  name: string;
  rawSequence: string;
  sequence: string;
  locus: GermlineLocus;
  metadata?: CompactMetadata;
}

interface Alignment {
  query: string;
  reference: string;
  identity: number;
}

interface TransferResult {
  metadata?: CompactMetadata;
  template?: string;
  identity?: number;
  matchKind?: "same_allele" | "same_gene" | "nearest";
  attemptedCandidates: number;
  bestCandidate?: string;
  bestIdentity?: number;
  rejectionCounts: Record<string, number>;
}

const LOCI: GermlineLocus[] = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
const IMGT_V_GAPPED_ENDS = [78, 114, 165, 195, 312] as const;
const EMPTY_BOUNDS = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1] as const;
const SOURCE_IMGT_GAPPED = 1;
const SOURCE_AIRR_C = 2;
const SOURCE_TRANSFERRED_IMGT = 3;
const SOURCE_VALIDATED_J_MOTIF = 4;
const SOURCE_PROVIDED = 5;
const SOURCE_TRANSFERRED_J = 6;

const MATCH_PRESETS: Record<GermlineMatchMode, Omit<ResolvedGermlineMatchOptions, "mode" | "includeDiagnostics">> = {
  strict: {
    vSameGeneMinIdentity: 0.80,
    vNearestMinIdentity: 0.72,
    jSameGeneMinIdentity: 0.75,
    jNearestMinIdentity: 0.68,
    nearestCandidates: 12,
  },
  permissive: {
    vSameGeneMinIdentity: 0.65,
    vNearestMinIdentity: 0.55,
    jSameGeneMinIdentity: 0.60,
    jNearestMinIdentity: 0.50,
    nearestCandidates: 32,
  },
  best_guess: {
    // Identity is deliberately not a rejection criterion in this mode. Motif,
    // coordinate, and frame checks remain hard invariants: emitting internally
    // contradictory SWIGMETA would be worse than reporting an unresolved allele.
    vSameGeneMinIdentity: 0,
    vNearestMinIdentity: 0,
    jSameGeneMinIdentity: 0,
    jNearestMinIdentity: 0,
    nearestCandidates: 64,
  },
};

function identityOption(value: number | undefined, fallback: number, label: string): number {
  const resolved = value ?? fallback;
  if (!Number.isFinite(resolved) || resolved < 0 || resolved > 1) {
    throw new Error(`${label} must be between 0 and 1.`);
  }
  return resolved;
}

export function resolveGermlineMatchOptions(options: GermlineMatchOptions = {}): ResolvedGermlineMatchOptions {
  const mode = options.mode ?? "strict";
  const preset = MATCH_PRESETS[mode];
  if (!preset) throw new Error(`Unsupported germline metadata match mode: ${String(mode)}.`);
  const nearestCandidates = options.nearestCandidates ?? preset.nearestCandidates;
  if (!Number.isSafeInteger(nearestCandidates) || nearestCandidates < 1 || nearestCandidates > 10_000) {
    throw new Error("nearestCandidates must be an integer between 1 and 10000.");
  }
  return {
    mode,
    vSameGeneMinIdentity: identityOption(options.vSameGeneMinIdentity, preset.vSameGeneMinIdentity, "vSameGeneMinIdentity"),
    vNearestMinIdentity: identityOption(options.vNearestMinIdentity, preset.vNearestMinIdentity, "vNearestMinIdentity"),
    jSameGeneMinIdentity: identityOption(options.jSameGeneMinIdentity, preset.jSameGeneMinIdentity, "jSameGeneMinIdentity"),
    jNearestMinIdentity: identityOption(options.jNearestMinIdentity, preset.jNearestMinIdentity, "jNearestMinIdentity"),
    nearestCandidates,
    includeDiagnostics: Boolean(options.includeDiagnostics),
  };
}

const V_CHAIN_LOCI: Record<string, GermlineLocus> = {
  VH: "IGH", VK: "IGK", VL: "IGL", VA: "TRA", VB: "TRB", VD: "TRD", VG: "TRG",
};
const J_CHAIN_LOCI: Record<string, GermlineLocus> = {
  JH: "IGH", JK: "IGK", JL: "IGL", JA: "TRA", JB: "TRB", JD: "TRD", JG: "TRG",
};

export const METADATA_SOURCE_LABELS: Record<number, string> = {
  [SOURCE_IMGT_GAPPED]: "IMGT-gapped delineation",
  [SOURCE_AIRR_C]: "AIRR-C annotation",
  [SOURCE_TRANSFERRED_IMGT]: "validated IMGT-boundary transfer",
  [SOURCE_VALIDATED_J_MOTIF]: "frame-validated J motif",
  [SOURCE_PROVIDED]: "provided annotation",
  [SOURCE_TRANSFERRED_J]: "validated J-anchor transfer",
};

function positiveModulo(value: number, modulus: number): number {
  return ((value % modulus) + modulus) % modulus;
}

function compactMetadata(
  frame = -1,
  cdr3Stop = -1,
  bounds: readonly number[] = EMPTY_BOUNDS,
  source = 0,
): CompactMetadata {
  return [frame, cdr3Stop, ...bounds, source] as CompactMetadata;
}

function inferLocus(value: string): GermlineLocus | null {
  const upper = value.toUpperCase();
  return LOCI.find((locus) => upper.includes(locus)) ?? null;
}

function inferSegment(name: string): GermlineSegment | null {
  const match = /^(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)([VDJC])/.exec(name.toUpperCase());
  return match ? match[1] as GermlineSegment : null;
}

function germlineName(header: string): string {
  const identifier = header.trim().split(/\s+/, 1)[0] ?? "";
  const fields = identifier.split("|");
  return fields.find((field) => inferLocus(field) && field.includes("*"))
    ?? fields.find((field) => inferLocus(field))
    ?? identifier;
}

function normalizeSequence(raw: string): { sequence: string; ambiguous: number } {
  const compact = raw.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
  if (/[^ACGTNRYKMSWBDHV.\-]/.test(compact)) {
    throw new Error("Germline FASTA contains characters outside the IUPAC nucleotide alphabet.");
  }
  const withoutGaps = compact.replace(/[.\-]/g, "");
  const ambiguous = (withoutGaps.match(/[^ACGT]/g) ?? []).length;
  return { sequence: withoutGaps.replace(/[^ACGT]/g, "N"), ambiguous };
}

function normalizeIndexSequence(raw: string): { sequence: string; ambiguous: number } {
  const compact = raw.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
  if (/[^ACGTNRYKMSWBDHV.\-]/.test(compact)) {
    throw new Error("Germline FASTA contains characters outside the IUPAC nucleotide alphabet.");
  }
  const sequence = compact.replace(/[.\-]/g, "");
  return { sequence, ambiguous: (sequence.match(/[^ACGT]/g) ?? []).length };
}

function parseFasta(text: string): Array<{ header: string; rawSequence: string }> {
  const records: Array<{ header: string; rawSequence: string }> = [];
  let header = "";
  let sequence: string[] = [];
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (header) records.push({ header, rawSequence: sequence.join("") });
      header = line.slice(1).trim();
      sequence = [];
    } else if (header && line.trim()) {
      sequence.push(line.trim());
    } else if (!header && line.trim()) {
      throw new Error("Expected a FASTA header beginning with '>'.");
    }
  }
  if (header) records.push({ header, rawSequence: sequence.join("") });
  if (!records.length) throw new Error("No germline FASTA records were found.");
  return records;
}

function parseImgtFields(header: string): string[] | null {
  const identifier = header.split(/\s+/, 1)[0] ?? "";
  const fields = identifier.split("|");
  return fields.length >= 8 && fields[1] && fields[4] ? fields : null;
}

function imgtFrame(fields: string[] | null): number {
  const value = Number.parseInt(fields?.[7] ?? "", 10);
  return value >= 1 && value <= 3 ? value - 1 : -1;
}

function metadataFromHeader(header: string): CompactMetadata | undefined {
  const match = /(?:^|\s)SWIGMETA=([\d,\-]+)/i.exec(header);
  if (!match) return undefined;
  const values = match[1].split(",").map(Number);
  if (values.length !== 13 || values.some((value) => !Number.isInteger(value))) return undefined;
  return values as CompactMetadata;
}

function metadataHeader(metadata?: CompactMetadata): string {
  return metadata ? ` SWIGMETA=${metadata.join(",")}` : "";
}

function hasRegionMetadata(metadata?: CompactMetadata): boolean {
  return Boolean(metadata && metadata.slice(2, 12).every((value) => value >= 0));
}

function hasJMetadata(metadata?: CompactMetadata): boolean {
  return Boolean(metadata && metadata[0] >= 0 && metadata[1] >= 0);
}

function validateMetadata(metadata: CompactMetadata, sequenceLength: number, segment: GermlineSegment): boolean {
  if (metadata[0] < -1 || metadata[0] > 2 || metadata[1] < -1 || metadata[1] >= sequenceLength) return false;
  if (segment === "V" && hasRegionMetadata(metadata)) {
    const bounds = metadata.slice(2, 12);
    if (bounds.some((value, index) => value < 0 || value > sequenceLength || (index && value < bounds[index - 1]))) return false;
    for (let index = 0; index < bounds.length; index += 2) {
      if (bounds[index + 1] <= bounds[index]) return false;
    }
  }
  return true;
}

function nearestFrameCysEnd(sequence: string, predictedEnd: number, frame: number, maximumDistance = 60): number | undefined {
  const candidates: number[] = [];
  const start = Math.max(0, predictedEnd - maximumDistance);
  const stop = Math.min(sequence.length - 2, predictedEnd + maximumDistance);
  for (let index = start; index <= stop; index += 1) {
    const codon = sequence.slice(index, index + 3);
    if (codon !== "TGT" && codon !== "TGC") continue;
    if (frame >= 0 && positiveModulo(index - frame, 3) !== 0) continue;
    candidates.push(index + 3);
  }
  return candidates.sort((left, right) => Math.abs(left - predictedEnd) - Math.abs(right - predictedEnd) || right - left)[0];
}

function imgtVMetadata(rawSequence: string, fields: string[] | null): CompactMetadata | undefined {
  if (fields?.[4]?.toUpperCase() !== "V-REGION") return undefined;
  const gapped = rawSequence.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
  if (gapped.length < IMGT_V_GAPPED_ENDS.at(-1)!) return undefined;
  const sequence = gapped.replace(/[.\-]/g, "");
  const ends = IMGT_V_GAPPED_ENDS.map((end) => gapped.slice(0, end).replace(/[.\-]/g, "").length);
  const anchorEnd = nearestFrameCysEnd(sequence, ends[4], imgtFrame(fields));
  if (anchorEnd && anchorEnd > ends[3]) ends[4] = anchorEnd;
  const bounds = [0, ends[0], ends[0], ends[1], ends[1], ends[2], ends[2], ends[3], ends[3], ends[4]];
  const result = compactMetadata(imgtFrame(fields), -1, bounds, SOURCE_IMGT_GAPPED);
  return validateMetadata(result, sequence.length, "V") ? result : undefined;
}

function isJAnchor(sequence: string, index: number): boolean {
  const aromatic = sequence.slice(index, index + 3);
  return (aromatic === "TGG" || aromatic === "TTT" || aromatic === "TTC")
    && /^GG[ACGT]$/.test(sequence.slice(index + 3, index + 6));
}

function jMetadata(sequence: string, fields: string[] | null): CompactMetadata | undefined {
  const providedFrame = imgtFrame(fields);
  const candidates: number[] = [];
  const limit = Math.min(sequence.length - 6, 54);
  for (let index = Math.max(0, providedFrame); index <= limit; index += providedFrame >= 0 ? 3 : 1) {
    if (isJAnchor(sequence, index)) candidates.push(index);
  }
  if (!candidates.length || (providedFrame < 0 && candidates.length !== 1)) return undefined;
  const anchor = candidates[0];
  const frame = providedFrame >= 0 ? providedFrame : anchor % 3;
  if (anchor % 3 !== frame) return undefined;
  return compactMetadata(frame, anchor - 1, EMPTY_BOUNDS, SOURCE_VALIDATED_J_MOTIF);
}

function globalAlignment(query: string, reference: string): Alignment {
  const columns = reference.length + 1;
  const rows = query.length + 1;
  const scores = new Int16Array(columns);
  const next = new Int16Array(columns);
  const trace = new Uint8Array(rows * columns);
  for (let column = 1; column < columns; column += 1) {
    scores[column] = -3 * column;
    trace[column] = 2;
  }
  for (let row = 1; row < rows; row += 1) {
    next[0] = -3 * row;
    trace[row * columns] = 1;
    for (let column = 1; column < columns; column += 1) {
      const diagonal = scores[column - 1] + (query[row - 1] === reference[column - 1] ? 2 : -2);
      const up = scores[column] - 3;
      const left = next[column - 1] - 3;
      const best = Math.max(diagonal, up, left);
      next[column] = best;
      trace[row * columns + column] = best === diagonal ? 0 : best === up ? 1 : 2;
    }
    scores.set(next);
  }
  let row = query.length;
  let column = reference.length;
  const alignedQuery: string[] = [];
  const alignedReference: string[] = [];
  let matches = 0;
  while (row || column) {
    const direction = trace[row * columns + column];
    if (row && column && direction === 0) {
      const queryBase = query[--row];
      const referenceBase = reference[--column];
      alignedQuery.push(queryBase);
      alignedReference.push(referenceBase);
      if (queryBase === referenceBase) matches += 1;
    } else if (row && (direction === 1 || !column)) {
      alignedQuery.push(query[--row]);
      alignedReference.push("-");
    } else {
      alignedQuery.push("-");
      alignedReference.push(reference[--column]);
    }
  }
  alignedQuery.reverse();
  alignedReference.reverse();
  return {
    query: alignedQuery.join(""),
    reference: alignedReference.join(""),
    identity: matches / Math.max(1, alignedQuery.length),
  };
}

function canonicalAllele(name: string): string {
  return name.replace(/_S\d+$/i, "").toUpperCase();
}

function canonicalGene(name: string): string {
  return canonicalAllele(name).split("*", 1)[0];
}

function kmerSet(sequence: string, size = 9): Set<string> {
  const output = new Set<string>();
  for (let index = 0; index + size <= sequence.length; index += 1) {
    const kmer = sequence.slice(index, index + size);
    if (!kmer.includes("N")) output.add(kmer);
  }
  return output;
}

// Reference alleles are immutable tuples from the fixed pack. Reusing their
// k-mer sets avoids rebuilding the same hundreds of sets for every custom
// allele, which was the dominant cost of large CLI metadata-preparation jobs.
const TEMPLATE_KMER_CACHE = new WeakMap<MetadataAllele, Set<string>>();

function templateKmerSet(template: MetadataAllele): Set<string> {
  const cached = TEMPLATE_KMER_CACHE.get(template);
  if (cached) return cached;
  const kmers = kmerSet(template[1]);
  TEMPLATE_KMER_CACHE.set(template, kmers);
  return kmers;
}

function nearestTemplates(
  query: string,
  templates: MetadataAllele[],
  limit = 12,
): MetadataAllele[] {
  const queryKmers = kmerSet(query);
  return templates
    .map((template) => {
      const templateKmers = templateKmerSet(template);
      let shared = 0;
      for (const kmer of templateKmers) if (queryKmers.has(kmer)) shared += 1;
      const union = queryKmers.size + templateKmers.size - shared;
      const similarity = union ? shared / union : 0;
      return { template, similarity, lengthDelta: Math.abs(query.length - template[1].length) };
    })
    .sort((a, b) => b.similarity - a.similarity || a.lengthDelta - b.lengthDelta)
    .slice(0, limit)
    .map(({ template }) => template);
}

function templateCandidateGroups(
  queryName: string,
  templates: MetadataAllele[],
  eligible: (metadata?: CompactMetadata) => boolean = hasRegionMetadata,
): {
  preferred: MetadataAllele[];
  preferredKind: "same_allele" | "same_gene";
  remaining: MetadataAllele[];
} {
  // Functional/pseudogene labels are deliberately irrelevant here. Any locus-matched
  // template with a complete IMGT delineation may transfer coordinates after the
  // sequence-level validation below.
  const delineated = templates.filter((template) => eligible(template[2]));
  const alleleName = canonicalAllele(queryName);
  const geneName = canonicalGene(queryName);
  const exactAllele = delineated.filter(([name]) => canonicalAllele(name) === alleleName);
  const exactGene = delineated.filter(([name]) => canonicalGene(name) === geneName);
  const preferred = exactAllele.length ? exactAllele : exactGene;
  const preferredKeys = new Set(preferred.map(([name, sequence]) => `${name}\u0000${sequence}`));
  return {
    preferred,
    preferredKind: exactAllele.length ? "same_allele" : "same_gene",
    remaining: delineated.filter(([name, sequence]) => !preferredKeys.has(`${name}\u0000${sequence}`)),
  };
}

function mapReferenceCoordinates(alignment: Alignment, referenceLength: number): number[] {
  const mapped = new Array<number>(referenceLength + 1).fill(-1);
  let queryPosition = 0;
  let referencePosition = 0;
  mapped[0] = 0;
  for (let column = 0; column < alignment.query.length; column += 1) {
    const queryBase = alignment.query[column];
    const referenceBase = alignment.reference[column];
    if (referenceBase === "-") {
      if (queryBase !== "-") queryPosition += 1;
      mapped[referencePosition] = queryPosition;
      continue;
    }
    if (mapped[referencePosition] < 0) mapped[referencePosition] = queryPosition;
    if (queryBase !== "-") queryPosition += 1;
    referencePosition += 1;
    mapped[referencePosition] = queryPosition;
  }
  return mapped;
}

function incrementRejection(result: TransferResult, reason: string): void {
  result.rejectionCounts[reason] = (result.rejectionCounts[reason] ?? 0) + 1;
}

function updateBest(result: TransferResult, template: MetadataAllele, identity: number): void {
  if (result.bestIdentity === undefined || identity > result.bestIdentity) {
    result.bestIdentity = identity;
    result.bestCandidate = template[0];
  }
}

function matchKind(name: string, templateName: string, fallback: "same_allele" | "same_gene" | "nearest"):
"same_allele" | "same_gene" | "nearest" {
  if (canonicalAllele(templateName) === canonicalAllele(name)) return "same_allele";
  if (canonicalGene(templateName) === canonicalGene(name)) return "same_gene";
  return fallback;
}

function transferMetadata(
  name: string,
  sequence: string,
  templates: MetadataAllele[],
  options: ResolvedGermlineMatchOptions,
): TransferResult {
  const result: TransferResult = { attemptedCandidates: 0, rejectionCounts: {} };
  const groups = templateCandidateGroups(name, templates);
  const candidates = [
    { values: groups.preferred, fallback: groups.preferredKind },
    // Compute the broad k-mer ranking only if every named candidate failed.
    { values: null, fallback: "nearest" as const },
  ];
  for (const group of candidates) {
    const values = group.values ?? nearestTemplates(sequence, groups.remaining, options.nearestCandidates);
    if (!values.length) continue;
    let selected: { metadata: CompactMetadata; identity: number; template: string; matchKind: "same_allele" | "same_gene" | "nearest" } | undefined;
    for (const template of values) {
      const alignment = globalAlignment(sequence, template[1]);
      const named = canonicalGene(template[0]) === canonicalGene(name);
      const relation = matchKind(name, template[0], group.fallback);
      result.attemptedCandidates += 1;
      updateBest(result, template, alignment.identity);
      if (alignment.identity < (named ? options.vSameGeneMinIdentity : options.vNearestMinIdentity)) {
        incrementRejection(result, "below_identity");
        continue;
      }
      const templateMetadata = template[2]!;
      const templateBounds = templateMetadata.slice(2, 12);
      const mapped = mapReferenceCoordinates(alignment, template[1].length);
      const bounds = templateBounds.map((boundary) => mapped[boundary]);
      if (bounds.some((value, index) => value < 0 || value > sequence.length || (index && value < bounds[index - 1]))) {
        incrementRejection(result, "unmapped_or_nonmonotonic_boundary");
        continue;
      }
      let valid = true;
      for (let index = 0; index < bounds.length; index += 2) {
        if (bounds[index + 1] <= bounds[index]) valid = false;
      }
      if (!valid) {
        incrementRejection(result, "empty_region");
        continue;
      }
      const templateFrame = templateMetadata[0] >= 0 ? templateMetadata[0] : templateBounds[0] % 3;
      const frame = positiveModulo(bounds[0] + templateFrame - templateBounds[0], 3);
      const anchorEnd = nearestFrameCysEnd(sequence, bounds[9], frame, 24);
      if (!anchorEnd || anchorEnd <= bounds[8]) {
        incrementRejection(result, "missing_frame_consistent_v_anchor");
        continue;
      }
      bounds[9] = anchorEnd;
      const metadata = compactMetadata(frame, -1, bounds, SOURCE_TRANSFERRED_IMGT);
      if (!validateMetadata(metadata, sequence.length, "V")) {
        incrementRejection(result, "invalid_projected_metadata");
        continue;
      }
      if (!selected || alignment.identity > selected.identity) {
        selected = { metadata, identity: alignment.identity, template: template[0], matchKind: relation };
      }
    }
    if (selected) return { ...result, ...selected };
  }
  if (!result.attemptedCandidates) incrementRejection(result, "no_annotated_template_candidates");
  return result;
}

function transferJMetadata(
  name: string,
  sequence: string,
  templates: MetadataAllele[],
  options: ResolvedGermlineMatchOptions,
): TransferResult {
  const result: TransferResult = { attemptedCandidates: 0, rejectionCounts: {} };
  const groups = templateCandidateGroups(name, templates, hasJMetadata);
  const candidates = [
    { values: groups.preferred, fallback: groups.preferredKind },
    { values: null, fallback: "nearest" as const },
  ];
  for (const group of candidates) {
    const values = group.values ?? nearestTemplates(sequence, groups.remaining, options.nearestCandidates);
    if (!values.length) continue;
    let selected: { metadata: CompactMetadata; identity: number; template: string; matchKind: "same_allele" | "same_gene" | "nearest" } | undefined;
    for (const template of values) {
      const alignment = globalAlignment(sequence, template[1]);
      const named = canonicalGene(template[0]) === canonicalGene(name);
      const relation = matchKind(name, template[0], group.fallback);
      result.attemptedCandidates += 1;
      updateBest(result, template, alignment.identity);
      if (alignment.identity < (named ? options.jSameGeneMinIdentity : options.jNearestMinIdentity)) {
        incrementRejection(result, "below_identity");
        continue;
      }
      const metadata = template[2]!;
      const referenceAnchor = metadata[1] + 1;
      if (referenceAnchor < 0 || referenceAnchor + 6 > template[1].length) {
        incrementRejection(result, "invalid_template_j_anchor");
        continue;
      }
      const mapped = mapReferenceCoordinates(alignment, template[1].length);
      const anchor = mapped[referenceAnchor];
      const anchorEnd = mapped[referenceAnchor + 6];
      if (anchor < 0 || anchorEnd - anchor !== 6) {
        incrementRejection(result, "incomplete_j_anchor_projection");
        continue;
      }
      if (!isJAnchor(sequence, anchor)) {
        incrementRejection(result, "target_j_motif_mismatch");
        continue;
      }
      const transferred = compactMetadata(anchor % 3, anchor - 1, EMPTY_BOUNDS, SOURCE_TRANSFERRED_J);
      if (!validateMetadata(transferred, sequence.length, "J")) {
        incrementRejection(result, "invalid_projected_metadata");
        continue;
      }
      if (!selected || alignment.identity > selected.identity) {
        selected = { metadata: transferred, identity: alignment.identity, template: template[0], matchKind: relation };
      }
    }
    if (selected) return { ...result, ...selected };
  }
  if (!result.attemptedCandidates) incrementRejection(result, "no_annotated_template_candidates");
  return result;
}

function annotationPresent(segment: GermlineSegment, metadata?: CompactMetadata): boolean {
  if (segment === "V") return hasRegionMetadata(metadata);
  if (segment === "J") return hasJMetadata(metadata);
  return true;
}

function integerField(value: string | undefined, context: string): number {
  if (value === undefined || !/^-?\d+$/.test(value)) throw new Error(`${context} must be an integer.`);
  return Number.parseInt(value, 10);
}

function annotationDataLines(text: string): Array<{ fields: string[]; line: number }> {
  const result: Array<{ fields: string[]; line: number }> = [];
  text.split(/\r?\n/).forEach((raw, index) => {
    const line = raw.trim();
    if (!line || line.startsWith("#")) return;
    result.push({ fields: line.split(/\s+/), line: index + 1 });
  });
  return result;
}

function assertChainLocus(
  name: string,
  locus: GermlineLocus,
  chain: string,
  chainLoci: Record<string, GermlineLocus>,
  source: string,
): void {
  const declared = chainLoci[chain.toUpperCase()];
  if (declared && declared !== locus) {
    throw new Error(`${source} declares ${name} as ${chain} (${declared}), but its germline identifier is ${locus}.`);
  }
}

/**
 * Perform the FASTA preparation done by IgBLAST's `edit_imgt_file.pl` where it
 * is relevant to the aligner: canonicalize the germline identifier and remove IMGT
 * alignment gaps before indexing. Existing SWIGMETA is deliberately removed.
 */
export function prepareIgblastStyleGermlineFasta(
  text: string,
  segment: GermlineSegment,
  allowedLoci: readonly GermlineLocus[] = LOCI,
): GermlineNormalizationReport {
  const inputRecords = parseFasta(text);
  const seen = new Set<string>();
  const loci = new Set<GermlineLocus>();
  let ambiguousBases = 0;
  const fasta = inputRecords.map((input) => {
    const name = germlineName(input.header);
    if (!name) throw new Error(`A ${segment} germline record has an empty identifier.`);
    if (seen.has(name)) throw new Error(`Duplicate ${segment} germline identifier: ${name}.`);
    seen.add(name);
    const locus = inferLocus(`${name} ${input.header}`);
    if (!locus) throw new Error(`${name} does not identify one of the supported IG/TR loci.`);
    if (!allowedLoci.includes(locus)) {
      throw new Error(`${name} belongs to ${locus}, outside the selected ${allowedLoci.join("/")} search space.`);
    }
    const detectedSegment = inferSegment(name);
    if (detectedSegment && detectedSegment !== segment) {
      throw new Error(`${name} looks like a ${detectedSegment} gene, not a ${segment} gene.`);
    }
    let normalized: { sequence: string; ambiguous: number };
    try { normalized = normalizeIndexSequence(input.rawSequence); }
    catch (error) { throw new Error(`${name}: ${error instanceof Error ? error.message : String(error)}`); }
    if (!normalized.sequence) throw new Error(`${name} has an empty nucleotide sequence.`);
    ambiguousBases += normalized.ambiguous;
    loci.add(locus);
    return `>${name}\n${normalized.sequence}\n`;
  }).join("");
  return { fasta, count: inputRecords.length, ambiguousBases, loci: [...loci].sort() };
}

/** Apply the exact, 1-based inclusive V-domain coordinates in IgBLAST `.ndm.imgt` data. */
export function applyIgblastInternalData(
  fasta: string,
  data: string,
): IgblastAnnotationApplication {
  type Entry = { bounds: number[]; chain: string; frame: number };
  const entries = new Map<string, Entry>();
  for (const record of annotationDataLines(data)) {
    if (record.fields.length !== 13) {
      throw new Error(`Invalid IgBLAST internal-data record at line ${record.line}; expected 13 fields.`);
    }
    const bounds: number[] = [];
    for (let index = 0; index < 10; index += 2) {
      const start = integerField(record.fields[index + 1], `IgBLAST internal-data line ${record.line} start`);
      const stop = integerField(record.fields[index + 2], `IgBLAST internal-data line ${record.line} stop`);
      if (start < 1 || stop < start) {
        throw new Error(`Invalid 1-based FWR/CDR interval at IgBLAST internal-data line ${record.line}.`);
      }
      bounds.push(start - 1, stop);
    }
    const frame = integerField(record.fields[12], `IgBLAST internal-data line ${record.line} coding frame`);
    if (frame < -1 || frame > 2) throw new Error(`Invalid coding frame at IgBLAST internal-data line ${record.line}.`);
    entries.set(record.fields[0], { bounds, chain: record.fields[11], frame });
  }
  if (!entries.size) throw new Error("The IgBLAST internal-data file contains no annotation records.");

  const prepared = prepareIgblastStyleGermlineFasta(fasta, "V");
  let matched = 0;
  const unmatched: string[] = [];
  const output = parseFasta(prepared.fasta).map((record) => {
    const name = germlineName(record.header);
    const normalized = normalizeIndexSequence(record.rawSequence);
    const locus = inferLocus(name)!;
    const entry = entries.get(name);
    if (!entry) { unmatched.push(name); return `>${name}\n${normalized.sequence}\n`; }
    assertChainLocus(name, locus, entry.chain, V_CHAIN_LOCI, "IgBLAST internal data");
    const metadata = compactMetadata(entry.frame, -1, entry.bounds, SOURCE_PROVIDED);
    if (!validateMetadata(metadata, normalized.sequence.length, "V")) {
      throw new Error(`IgBLAST internal-data coordinates for ${name} exceed or contradict its ungapped V sequence.`);
    }
    matched += 1;
    return `>${name}${metadataHeader(metadata)}\n${normalized.sequence}\n`;
  }).join("");
  return { fasta: output, matched, annotated: matched, total: prepared.count, unmatched };
}

/**
 * Apply IgBLAST J auxiliary data. Frame and CDR3-stop coordinates are 0-based;
 * the optional fifth field is retained separately for FWR4 end trimming.
 */
export function applyIgblastAuxiliaryData(
  fasta: string,
  data: string,
): IgblastAnnotationApplication {
  type Entry = { frame: number; chain: string; cdr3Stop: number; fwr4EndOffset?: number };
  const entries = new Map<string, Entry>();
  for (const record of annotationDataLines(data)) {
    if (record.fields.length < 3 || record.fields.length > 5) {
      throw new Error(`Invalid IgBLAST auxiliary record at line ${record.line}; expected 3 to 5 fields.`);
    }
    const frame = integerField(record.fields[1], `IgBLAST auxiliary line ${record.line} coding frame`);
    if (frame < -1 || frame > 2) throw new Error(`Invalid coding frame at IgBLAST auxiliary line ${record.line}.`);
    const cdr3Stop = record.fields.length >= 4
      ? integerField(record.fields[3], `IgBLAST auxiliary line ${record.line} CDR3 stop`) : -1;
    if (cdr3Stop < -1) throw new Error(`Invalid CDR3 stop at IgBLAST auxiliary line ${record.line}.`);
    const fwr4EndOffset = record.fields.length === 5
      ? integerField(record.fields[4], `IgBLAST auxiliary line ${record.line} FWR4 end offset`) : undefined;
    if (fwr4EndOffset !== undefined && fwr4EndOffset < 0) {
      throw new Error(`Invalid FWR4 end offset at IgBLAST auxiliary line ${record.line}.`);
    }
    // IgBLAST itself stores these in maps, so a repeated identifier follows
    // the last record. Its distributed human file contains such a duplicate.
    entries.set(record.fields[0], { frame, chain: record.fields[2], cdr3Stop, fwr4EndOffset });
  }
  if (!entries.size) throw new Error("The IgBLAST auxiliary file contains no annotation records.");

  const prepared = prepareIgblastStyleGermlineFasta(fasta, "J");
  let matched = 0;
  let annotated = 0;
  const unmatched: string[] = [];
  const fwr4EndOffsets: Record<string, number> = {};
  const output = parseFasta(prepared.fasta).map((record) => {
    const name = germlineName(record.header);
    const normalized = normalizeIndexSequence(record.rawSequence);
    const locus = inferLocus(name)!;
    const entry = entries.get(name);
    if (!entry) { unmatched.push(name); return `>${name}\n${normalized.sequence}\n`; }
    assertChainLocus(name, locus, entry.chain, J_CHAIN_LOCI, "IgBLAST auxiliary data");
    const metadata = compactMetadata(entry.frame, entry.cdr3Stop, EMPTY_BOUNDS, SOURCE_PROVIDED);
    if (!validateMetadata(metadata, normalized.sequence.length, "J")) {
      throw new Error(`IgBLAST auxiliary coordinates for ${name} exceed its ungapped J sequence.`);
    }
    matched += 1;
    if (entry.cdr3Stop >= 0) annotated += 1;
    if (entry.fwr4EndOffset !== undefined) fwr4EndOffsets[name] = entry.fwr4EndOffset;
    return `>${name}${metadataHeader(metadata)}\n${normalized.sequence}\n`;
  }).join("");
  return { fasta: output, matched, annotated, total: prepared.count, unmatched, fwr4EndOffsets };
}

/** Apply the exact 0-based coding-frame-one starts in IgBLAST `-d_frame_data`. */
export function applyIgblastDFrameData(
  fasta: string,
  data: string,
): IgblastAnnotationApplication {
  const entries = new Map<string, number>();
  for (const record of annotationDataLines(data)) {
    if (record.fields.length !== 2) {
      throw new Error(`Invalid IgBLAST D-frame record at line ${record.line}; expected 2 fields.`);
    }
    const frame = integerField(record.fields[1], `IgBLAST D-frame line ${record.line} start`);
    if (frame < -1 || frame > 2) throw new Error(`Invalid D-frame start at line ${record.line}.`);
    entries.set(record.fields[0], frame);
  }
  if (!entries.size) throw new Error("The IgBLAST D-frame file contains no annotation records.");

  const prepared = prepareIgblastStyleGermlineFasta(fasta, "D");
  let matched = 0;
  let annotated = 0;
  const unmatched: string[] = [];
  const output = parseFasta(prepared.fasta).map((record) => {
    const name = germlineName(record.header);
    const normalized = normalizeIndexSequence(record.rawSequence);
    const frame = entries.get(name);
    if (frame === undefined) { unmatched.push(name); return `>${name}\n${normalized.sequence}\n`; }
    matched += 1;
    if (frame < 0) return `>${name}\n${normalized.sequence}\n`;
    annotated += 1;
    const metadata = compactMetadata(frame, -1, EMPTY_BOUNDS, SOURCE_PROVIDED);
    return `>${name}${metadataHeader(metadata)}\n${normalized.sequence}\n`;
  }).join("");
  return { fasta: output, matched, annotated, total: prepared.count, unmatched };
}

export function preprocessGermlineFasta(
  text: string,
  segment: GermlineSegment,
  templates: MetadataAllele[] = [],
  allowedLoci: readonly GermlineLocus[] = LOCI,
  match: GermlineMatchOptions = {},
): GermlinePreprocessReport {
  const matchOptions = resolveGermlineMatchOptions(match);
  const inputRecords = parseFasta(text);
  const records: ParsedRecord[] = [];
  const diagnostics: GermlineRecordDiagnostic[] = [];
  const seen = new Set<string>();
  const loci = new Set<GermlineLocus>();
  const warnings: string[] = [];
  let exactImgt = 0;
  let transferred = 0;
  let motifValidated = 0;
  let ambiguousBases = 0;

  for (const input of inputRecords) {
    const name = germlineName(input.header);
    if (!name) throw new Error(`A ${segment} germline record has an empty identifier.`);
    if (seen.has(name)) throw new Error(`Duplicate ${segment} germline identifier: ${name}.`);
    seen.add(name);
    const locus = inferLocus(`${name} ${input.header}`);
    if (!locus) throw new Error(`${name} does not identify one of the supported IG/TR loci.`);
    if (!allowedLoci.includes(locus)) {
      throw new Error(`${name} belongs to ${locus}, outside the selected ${allowedLoci.join("/")} search space.`);
    }
    const detectedSegment = inferSegment(name);
    if (detectedSegment && detectedSegment !== segment) {
      throw new Error(`${name} looks like a ${detectedSegment} gene, not a ${segment} gene.`);
    }
    let normalized: { sequence: string; ambiguous: number };
    try { normalized = normalizeSequence(input.rawSequence); }
    catch (error) { throw new Error(`${name}: ${error instanceof Error ? error.message : String(error)}`); }
    if (!normalized.sequence) throw new Error(`${name} has an empty nucleotide sequence.`);
    ambiguousBases += normalized.ambiguous;
    const fields = parseImgtFields(input.header);
    let metadata = metadataFromHeader(input.header);
    let status: GermlineDiagnosticStatus = metadata ? "retained" : "unresolved";
    let transfer: TransferResult | undefined;
    if (metadata && !validateMetadata(metadata, normalized.sequence.length, segment)) {
      throw new Error(`${name} contains invalid SWIGMETA coordinates.`);
    }
    if (metadata) metadata = [...metadata.slice(0, 12), metadata[12] || SOURCE_PROVIDED] as CompactMetadata;
    if (!metadata && segment === "V") {
      metadata = imgtVMetadata(input.rawSequence, fields);
      if (metadata) status = "imgt";
    }
    if (!metadata && segment === "V" && templates.length) {
      transfer = transferMetadata(name, normalized.sequence, templates, matchOptions);
      metadata = transfer.metadata;
      if (metadata) status = "transferred";
    }
    if (!metadata && segment === "J" && templates.length) {
      transfer = transferJMetadata(name, normalized.sequence, templates, matchOptions);
      metadata = transfer.metadata;
      if (metadata) status = "transferred";
    }
    if (!metadata && segment === "J") {
      metadata = jMetadata(normalized.sequence, fields);
      if (metadata) status = "motif";
    }
    if (!metadata && segment === "D") {
      const frame = imgtFrame(fields);
      if (frame >= 0) {
        metadata = compactMetadata(frame, -1, EMPTY_BOUNDS, SOURCE_IMGT_GAPPED);
        status = "imgt";
      }
    }
    if (!metadata && (segment === "D" || segment === "C")) status = "normalized";
    if (metadata?.[12] === SOURCE_IMGT_GAPPED) exactImgt += 1;
    if (metadata?.[12] === SOURCE_TRANSFERRED_IMGT || metadata?.[12] === SOURCE_TRANSFERRED_J) transferred += 1;
    if (metadata?.[12] === SOURCE_VALIDATED_J_MOTIF) motifValidated += 1;
    if (!annotationPresent(segment, metadata) && warnings.length < 8) {
      warnings.push(`${name}: ${segment === "V" ? "no validated IMGT FWR/CDR delineation" : segment === "J" ? "no unique frame-consistent F/W-G J motif" : "no segment metadata"}.`);
    }
    loci.add(locus);
    records.push({ header: input.header, name, rawSequence: input.rawSequence, sequence: normalized.sequence, locus, metadata });
    if (matchOptions.includeDiagnostics) {
      diagnostics.push({
        segment,
        name,
        locus,
        status,
        source: metadataSource(metadata),
        template: transfer?.template,
        identity: transfer?.identity,
        matchKind: transfer?.matchKind,
        taxonomicTier: transfer?.metadata ? match.taxonomicTier : undefined,
        attemptedCandidates: transfer?.attemptedCandidates ?? 0,
        bestCandidate: transfer?.bestCandidate,
        bestIdentity: transfer?.bestIdentity,
        rejectionCounts: transfer && Object.keys(transfer.rejectionCounts).length ? transfer.rejectionCounts : undefined,
      });
    }
  }

  const annotated = records.filter((record) => annotationPresent(segment, record.metadata)).length;
  const fasta = records.map((record) => `>${record.name}${metadataHeader(record.metadata)}\n${record.sequence}\n`).join("");
  return {
    fasta,
    count: records.length,
    annotated,
    unannotated: records.length - annotated,
    exactImgt,
    transferred,
    motifValidated,
    ambiguousBases,
    loci: [...loci].sort(),
    warnings,
    diagnostics: matchOptions.includeDiagnostics ? diagnostics : undefined,
  };
}

function mergeDiagnosticHistory(
  previous: GermlineRecordDiagnostic | undefined,
  current: GermlineRecordDiagnostic,
): GermlineRecordDiagnostic {
  if (!previous) return current;
  // A metadata record written by an earlier tier is parsed as "retained" on
  // the next pass. Preserve the original transfer source and template.
  if (current.status === "retained" && previous.status !== "unresolved") return previous;
  const rejectionCounts = { ...(previous.rejectionCounts ?? {}) };
  for (const [reason, count] of Object.entries(current.rejectionCounts ?? {})) {
    rejectionCounts[reason] = (rejectionCounts[reason] ?? 0) + count;
  }
  const previousBest = previous.bestIdentity ?? -1;
  const currentBest = current.bestIdentity ?? -1;
  const latest = current.status !== "unresolved" ? current : previous;
  return {
    ...latest,
    attemptedCandidates: previous.attemptedCandidates + current.attemptedCandidates,
    bestCandidate: currentBest > previousBest ? current.bestCandidate : previous.bestCandidate,
    bestIdentity: Math.max(previousBest, currentBest) >= 0 ? Math.max(previousBest, currentBest) : undefined,
    rejectionCounts: Object.keys(rejectionCounts).length ? rejectionCounts : undefined,
  };
}

/**
 * Apply the same progressively broadened template search used by the browser
 * worker. Existing valid SWIGMETA is retained, while records still lacking V
 * or J metadata are offered to each successive template tier.
 */
export function preprocessGermlineFastaAcrossTiers(
  text: string,
  segment: GermlineSegment,
  templateTiers: readonly (readonly MetadataAllele[])[],
  allowedLoci: readonly GermlineLocus[] = LOCI,
  match: GermlineMatchOptions = {},
): GermlinePreprocessReport {
  let report: GermlinePreprocessReport | undefined;
  const diagnosticHistory = new Map<string, GermlineRecordDiagnostic>();
  const tiers = templateTiers.length ? templateTiers : [[]];
  for (let tier = 0; tier < tiers.length; tier += 1) {
    report = preprocessGermlineFasta(report?.fasta ?? text, segment, [...tiers[tier]], allowedLoci, {
      ...match,
      includeDiagnostics: Boolean(match.includeDiagnostics),
      taxonomicTier: tier,
    });
    for (const diagnostic of report.diagnostics ?? []) {
      diagnosticHistory.set(diagnostic.name, mergeDiagnosticHistory(diagnosticHistory.get(diagnostic.name), diagnostic));
    }
    if (segment !== "V" && segment !== "J") break;
    if (report.annotated === report.count) break;
  }
  if (match.includeDiagnostics) report!.diagnostics = [...diagnosticHistory.values()];
  return report!;
}

export function annotationCoverage(fasta: string, segment: GermlineSegment): { annotated: number; total: number } {
  const records = parseFasta(fasta);
  let annotated = 0;
  for (const record of records) {
    const normalized = normalizeSequence(record.rawSequence);
    const metadata = metadataFromHeader(record.header);
    if (metadata && validateMetadata(metadata, normalized.sequence.length, segment) && annotationPresent(segment, metadata)) annotated += 1;
  }
  return { annotated, total: records.length };
}

export function alleleMetadataHeader(allele: MetadataAllele): string {
  return `>${allele[0]}${metadataHeader(allele[2])}\n${allele[1]}\n`;
}

export function metadataSource(metadata?: CompactMetadata): string {
  return metadata ? METADATA_SOURCE_LABELS[metadata[12]] ?? "unclassified annotation" : "none";
}
