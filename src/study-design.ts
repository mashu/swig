export type StudyDesign = "independent" | "cohort" | "longitudinal" | "technical" | "custom";

/**
 * Dataset is one imported file/library. Sample is one biological specimen;
 * technical replicate datasets may deliberately share a sample identifier.
 */
export type DatasetScope = "dataset" | "sample" | "subject" | "cohort" | "global";

export interface DatasetManifestEntry {
  datasetId: string;
  inputName: string;
  /** CLI-resolvable source path; may differ from inputName for one member of a concatenated gzip. */
  inputPath?: string;
  inputFormat?: "auto" | "fasta" | "fastq" | "airr";
  /** Compressed-byte range for one complete member of a concatenated gzip source. */
  gzipRange?: { start: number; end: number };
  sampleId: string;
  subjectId: string;
  cohort: string;
  timepoint: string;
  /** Optional anatomical compartment, tissue, or sampling site. */
  compartment?: string;
  records?: number | null;
}

export interface PipelinePlan {
  enabled: boolean;
  collapse: {
    enabled: boolean;
    mode: "exact" | "fad" | "conservative" | "indel";
    key: "sequence" | "trimmed" | "cdr3" | "rearrangement";
    scope: DatasetScope;
    unresolvedPolicy: "discard" | "retain";
    /** Partition on the normalized top C-gene/isotype call, independent of tail length. */
    respectConstantCall: boolean;
  };
  chimera: {
    enabled: boolean;
    segment: "V" | "J";
    model: "auto" | "BW" | "DB";
    posteriorThreshold: number;
    retainUnevaluated: boolean;
    msaSource: "selected" | "upload";
    uploadedMsa: string;
    uploadedMsaName: string;
  };
  selection: {
    enabled: boolean;
    datasetId: string;
    sampleId: string;
    subjectId: string;
    cohort: string;
    timepoint: string;
    compartment: string;
    locus: string;
    vCall: string;
    vCallIncludeAmbiguous: boolean;
    jCall: string;
    jCallIncludeAmbiguous: boolean;
    cdr3Nt: string;
    cdr3Aa: string;
    productive: "any" | "yes" | "no";
    hasCdr3: "any" | "yes" | "no";
    doubleD: "any" | "positive" | "negative";
  };
  alleleRefinement: {
    enabled: boolean;
    model: "dirichlet" | "active-set";
    scope: DatasetScope;
    segments: Array<"V" | "D" | "J">;
    weighting: "unique" | "abundance";
    baselineNeighbourOdds: number;
    shmLeakageSensitivity: number;
    /** Apply every posterior MAP call, or retain original calls below a confidence gate. */
    reassignmentPolicy: "best" | "confidence";
    applyMinimumPosterior: number;
  };
  lineage: {
    enabled: boolean;
    scope: DatasetScope;
    identity: number;
    resolution: "gene" | "allele";
    ambiguity: "overlap" | "top" | "strict";
    productiveOnly: boolean;
  };
  shm: {
    enabled: boolean;
    metric: "vNtRate" | "vNtMutations" | "vAaRate" | "vAaReplacements" | "synonymous" | "cdrNtRate" | "frameworkNtRate";
  };
  missingAlleles: {
    enabled: boolean;
  };
}

export const DEFAULT_PIPELINE_PLAN: PipelinePlan = {
  enabled: false,
  collapse: {
    enabled: true,
    mode: "exact",
    key: "trimmed",
    scope: "sample",
    unresolvedPolicy: "discard",
    respectConstantCall: true,
  },
  chimera: {
    enabled: false,
    segment: "V",
    model: "auto",
    posteriorThreshold: 0.95,
    retainUnevaluated: true,
    msaSource: "selected",
    uploadedMsa: "",
    uploadedMsaName: "",
  },
  selection: {
    enabled: false,
    datasetId: "",
    sampleId: "",
    subjectId: "",
    cohort: "",
    timepoint: "",
    compartment: "",
    locus: "",
    vCall: "",
    vCallIncludeAmbiguous: false,
    jCall: "",
    jCallIncludeAmbiguous: false,
    cdr3Nt: "",
    cdr3Aa: "",
    productive: "any",
    hasCdr3: "any",
    doubleD: "any",
  },
  alleleRefinement: {
    enabled: false,
    model: "dirichlet",
    scope: "subject",
    segments: ["V", "J"],
    weighting: "unique",
    baselineNeighbourOdds: 0.01,
    shmLeakageSensitivity: 1,
    reassignmentPolicy: "confidence",
    applyMinimumPosterior: 0.8,
  },
  lineage: {
    enabled: true,
    scope: "subject",
    identity: 0.85,
    resolution: "gene",
    ambiguity: "overlap",
    productiveOnly: true,
  },
  shm: {
    enabled: true,
    metric: "vNtRate",
  },
  missingAlleles: {
    enabled: false,
  },
};

export const DATASET_SCOPE_LABELS: Record<DatasetScope, string> = {
  dataset: "Loaded dataset / library",
  sample: "Biological sample",
  subject: "Donor / subject",
  cohort: "Cohort",
  global: "Entire study",
};

/** Explicit pooling language: cross-participant scopes are deliberate overrides. */
export const ALLELE_POOL_SCOPE_LABELS: Record<DatasetScope, string> = {
  dataset: "Loaded dataset / library only",
  sample: "Biological sample only",
  subject: "Donor / subject · default; combine samples",
  cohort: "Cohort · cross-donor override",
  global: "Entire study · cross-donor override",
};

export function studyScopeDefaults(design: StudyDesign): { collapse: DatasetScope; lineage: DatasetScope } {
  if (design === "longitudinal") return { collapse: "sample", lineage: "subject" };
  if (design === "technical") return { collapse: "sample", lineage: "sample" };
  return { collapse: "sample", lineage: "sample" };
}

export function datasetScopeKey(record: {
  datasetId?: string;
  sampleId?: string;
  subjectId?: string;
  cohort?: string;
}, scope: DatasetScope = "global"): string {
  return `${scope}:${datasetScopeValue(record, scope)}`;
}

export function datasetScopeValue(record: {
  datasetId?: string;
  sampleId?: string;
  subjectId?: string;
  cohort?: string;
}, scope: DatasetScope = "global"): string {
  if (scope === "global") return "complete study";
  if (scope === "dataset") return record.datasetId || "legacy";
  if (scope === "sample") return record.sampleId || record.datasetId || "legacy";
  if (scope === "subject") return record.subjectId || record.sampleId || record.datasetId || "legacy";
  return record.cohort || record.subjectId || record.sampleId || record.datasetId || "legacy";
}

const METADATA_FIELDS = [
  "swig_dataset_id",
  "sample_id",
  "subject_id",
  "swig_cohort",
  "swig_timepoint",
  "swig_compartment",
  "swig_source_sequence_id",
] as const;

function normalizedLines(body: string | Uint8Array): string[] {
  const text = typeof body === "string" ? body : new TextDecoder().decode(body);
  return text.split("\n").map((line) => line.replace(/\r$/, "")).filter(Boolean);
}

function replaceOrAppend(values: string[], positions: Map<string, number>, field: string, value: string) {
  const position = positions.get(field);
  if (position === undefined) values.push(value);
  else values[position] = value;
}

/** Adds study metadata to an aligner AIRR batch without materializing the input file. */
export function annotateAirrBatch(
  headerLine: string,
  body: string | Uint8Array,
  dataset: DatasetManifestEntry,
): { header: string; body: string } {
  const headers = headerLine.replace(/\r$/, "").split("\t");
  const positions = new Map(headers.map((field, index) => [field, index]));
  for (const field of METADATA_FIELDS) if (!positions.has(field)) headers.push(field);
  const sourceSequencePosition = positions.get("sequence_id");
  const bodyText = normalizedLines(body).map((line) => {
    const values = line.split("\t");
    const sourceSequenceId = sourceSequencePosition === undefined ? "" : values[sourceSequencePosition] ?? "";
    if (sourceSequencePosition !== undefined) values[sourceSequencePosition] = `${dataset.datasetId}::${sourceSequenceId || "record"}`;
    replaceOrAppend(values, positions, "swig_dataset_id", dataset.datasetId);
    replaceOrAppend(values, positions, "sample_id", dataset.sampleId);
    replaceOrAppend(values, positions, "subject_id", dataset.subjectId);
    replaceOrAppend(values, positions, "swig_cohort", dataset.cohort);
    replaceOrAppend(values, positions, "swig_timepoint", dataset.timepoint ?? "");
    replaceOrAppend(values, positions, "swig_compartment", dataset.compartment ?? "");
    replaceOrAppend(values, positions, "swig_source_sequence_id", sourceSequenceId);
    return values.join("\t");
  }).join("\n");
  return { header: headers.join("\t"), body: bodyText ? `${bodyText}\n` : "" };
}

/** Mirrors metadata and collision-safe IDs into the sparse double-D table. */
export function annotateDoubleDBatch(
  headerLine: string,
  body: string | Uint8Array,
  dataset: DatasetManifestEntry,
): { header: string; body: string } {
  const headers = headerLine.replace(/\r$/, "").split("\t");
  const positions = new Map(headers.map((field, index) => [field, index]));
  for (const field of METADATA_FIELDS.slice(0, 6)) if (!positions.has(field)) headers.push(field);
  const sequencePosition = positions.get("sequence_id");
  const bodyText = normalizedLines(body).map((line) => {
    const values = line.split("\t");
    if (sequencePosition !== undefined) values[sequencePosition] = `${dataset.datasetId}::${values[sequencePosition] || "record"}`;
    replaceOrAppend(values, positions, "swig_dataset_id", dataset.datasetId);
    replaceOrAppend(values, positions, "sample_id", dataset.sampleId);
    replaceOrAppend(values, positions, "subject_id", dataset.subjectId);
    replaceOrAppend(values, positions, "swig_cohort", dataset.cohort);
    replaceOrAppend(values, positions, "swig_timepoint", dataset.timepoint ?? "");
    replaceOrAppend(values, positions, "swig_compartment", dataset.compartment ?? "");
    return values.join("\t");
  }).join("\n");
  return { header: headers.join("\t"), body: bodyText ? `${bodyText}\n` : "" };
}

export function stableDatasetSeed(baseSeed: number, datasetIndex: number): number {
  return (Math.trunc(baseSeed) + Math.imul(datasetIndex + 1, 0x6d2b79f5)) | 0;
}
