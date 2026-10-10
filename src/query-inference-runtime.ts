import type { QueryTarget } from "./post-analysis-core";
import type { CompiledReferences } from "./reference-pack";
import { AirrResultStore } from "./result-store";
import type { CallingProfile } from "./aligner-engine";
import { runAligner } from "./aligner-runtime";

export interface InferredQueryAssignment {
  queryIndex: number;
  sequenceId: string;
  locus: string;
  vCall: string;
  jCall: string;
  cdr3Nt: string;
  cdr3Aa: string;
  trimmed: string;
  searchSequence: string;
  assigned: boolean;
}

function fastaForQueries(queries: string[]): string {
  return queries.map((sequence, index) => `>swig_seed_${index + 1}\n${sequence.toUpperCase().replace(/[^ACGTN]/g, "N")}\n`).join("");
}

function targetSequence(target: QueryTarget, values: Record<string, string>, fallback: string): string {
  if (target === "cdr3_nt") return values.cdr3 || fallback;
  if (target === "cdr3_aa") return values.cdr3_aa || "";
  return values.sequence_alignment || values.sequence || fallback;
}

export async function inferQueryAssignments(
  queries: string[],
  target: QueryTarget,
  references: CompiledReferences,
  callingProfile: CallingProfile,
  minimumIdentity: number,
  strand: 0 | 1 | 2,
  workers: number,
  signal?: AbortSignal,
): Promise<InferredQueryAssignment[]> {
  const store = new AirrResultStore();
  try {
    await runAligner({
      query: fastaForQueries(queries),
      format: 1,
      references,
      settings: { callingProfile, minimumIdentity, strand },
      workers: Math.max(1, Math.min(workers, queries.length, 8)),
      countHint: queries.length,
      signal,
      onBatch: (batch) => store.appendBatch(batch.header, batch.body),
    });
    await store.finalize();
    const rows = new Map<number, Record<string, string>>();
    await store.scanAirrRows(
      ["sequence_id", "sequence", "sequence_alignment", "locus", "v_call", "j_call", "cdr3", "cdr3_aa"],
      (batch) => batch.forEach((row) => rows.set(row.ordinal, row.values)),
      { batchSize: 100 },
    );
    return queries.map((query, queryIndex) => {
      const values = rows.get(queryIndex) ?? {};
      const vCall = values.v_call ?? "";
      const jCall = values.j_call ?? "";
      return {
        queryIndex,
        sequenceId: values.sequence_id || `swig_seed_${queryIndex + 1}`,
        locus: values.locus ?? "",
        vCall,
        jCall,
        cdr3Nt: values.cdr3 ?? "",
        cdr3Aa: values.cdr3_aa ?? "",
        trimmed: values.sequence_alignment ?? "",
        searchSequence: targetSequence(target, values, query),
        assigned: Boolean(vCall && jCall),
      };
    });
  } finally {
    await store.clear().catch(() => undefined);
  }
}
