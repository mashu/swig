import { normalizeNt, sequenceFingerprint, type PostAnalysisRecord } from "./post-analysis-core.ts";
import { orientedSequence } from "./nucleotide.ts";

export interface OrdinalAirrRow {
  ordinal: number;
  values: Record<string,string>;
}

/** Shared browser-worker/CLI conversion into the compact scientific record. */
export function airrRowToPostAnalysisRecord(row: OrdinalAirrRow, intern: (value:string)=>string=(value)=>value): PostAnalysisRecord {
  const value=(name:string)=>row.values[name]??"";
  const sequence=normalizeNt(value("sequence"));
  return {
    ordinal:row.ordinal,
    sequenceId:value("sequence_id"),
    datasetId:intern(value("swig_dataset_id")),
    sampleId:intern(value("sample_id")),
    subjectId:intern(value("subject_id")),
    cohort:intern(value("swig_cohort")),
    timepoint:intern(value("swig_timepoint")),
    compartment:intern(value("swig_compartment")),
    locus:intern(value("locus")),
    vCall:intern(value("v_call")),
    jCall:intern(value("j_call")),
    originalVCall:intern(value("v_call")),
    originalJCall:intern(value("j_call")),
    cCall:intern(value("c_call")),
    cdr3Nt:normalizeNt(value("cdr3")),
    cdr3Aa:value("cdr3_aa").toUpperCase().replace(/[^A-Z*]/g,""),
    productive:/^(?:T|TRUE)$/i.test(value("productive")),
    sequenceFingerprint:sequenceFingerprint(sequence),
    trimmedFingerprint:sequenceFingerprint(normalizeNt(value("sequence_alignment")||value("sequence"))),
    inputCount:Math.max(1,Math.floor(Number(value("duplicate_count"))||1)),
  };
}

/**
 * Exact V-to-J slice used by every denoising front end. AIRR coordinates are
 * one-based and refer to the read as aligned: its reverse complement when
 * `rev_comp` is true.
 */
export function denoiseVdjSequence(row: OrdinalAirrRow): string {
  const raw=orientedSequence(row.values.sequence,row.values.rev_comp);
  const start=Math.floor(Number(row.values.v_sequence_start));
  const end=Math.floor(Number(row.values.j_sequence_end));
  if(raw&&Number.isFinite(start)&&Number.isFinite(end)&&start>=1&&end>=start&&end<=raw.length)return raw.slice(start-1,end);
  return row.values.sequence_alignment||raw;
}
