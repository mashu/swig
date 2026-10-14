import type { LocusKey, SegmentKey } from "./reference-pack";
import { filterReferenceFasta } from "./reference-fasta.ts";
import { namedLoci } from "./germline-preprocess.ts";

const SEGMENTS: SegmentKey[] = ["V", "D", "J", "C"];

export function referenceCellKey(locus: LocusKey, segment: SegmentKey): string {
  return `${locus}:${segment}`;
}

export function segmentAppliesToLocus(locus: LocusKey, segment: SegmentKey): boolean {
  return segment !== "D" || locus === "IGH" || locus === "TRB" || locus === "TRD";
}

/** The allele label the engine reads from a FASTA header (engine `GermlineSet::allele_label`). */
function engineAlleleLabel(header: string): string {
  const id = header.trim().split(/\s+/, 1)[0] ?? "";
  const fields = id.split("|");
  return fields.find((field) => field.includes("*")) ?? fields[0] ?? id;
}

/**
 * Custom FASTA of one reference cell, each record whose name places no locus
 * given the cell's (`locus=`), so the engine pairs, prices and reports it in
 * that locus instead of inferring one. Records with a `locus=` key or a name
 * the engine reads a locus from are left as they are.
 */
export function withCellLocus(text: string, locus: LocusKey): string {
  return text.split(/(\r?\n)/).map((line) => {
    if (!line.startsWith(">")) return line;
    if (/(?:^|\s)locus=/.test(line.slice(1))) return line;
    if (namedLoci(engineAlleleLabel(line.slice(1))).length) return line;
    return `${line.replace(/\s+$/, "")} locus=${locus}`;
  }).join("");
}

export function composeReferenceOverrides(
  loci: LocusKey[],
  cells: Record<string, { text: string }>,
  baseline: (locus: LocusKey, segment: SegmentKey) => string,
  exclusions: Record<string, readonly string[]> = {},
): Partial<Record<SegmentKey, string>> {
  const result: Partial<Record<SegmentKey, string>> = {};
  for (const segment of SEGMENTS) {
    if (!loci.some((locus) => {
      const key = referenceCellKey(locus, segment);
      return Boolean(cells[key] || exclusions[key]?.length);
    })) continue;
    result[segment] = loci.map((locus) => {
      const key = referenceCellKey(locus, segment);
      const custom = cells[key]?.text;
      const source = custom === undefined ? baseline(locus, segment) : withCellLocus(custom, locus);
      return exclusions[key]?.length ? filterReferenceFasta(source, exclusions[key]).fasta : source;
    }).filter(Boolean).join("\n");
  }
  return result;
}
