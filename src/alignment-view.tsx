import { useMemo } from "react";
import {
  ALIGNMENT_SEGMENTS,
  REGIONS,
  biologicalSegmentAlignment,
  buildTrackFeatures,
  type AirrRow,
  type AlignmentMode,
} from "./alignment-model";
import { ColoredSequence } from "./sequence-colors";

function matchLine(query: string, reference: string): string {
  return [...query].map((value, index) => {
    const other = reference[index];
    if (value === "-" || other === "-") return " ";
    return value === other ? "│" : "·";
  }).join("");
}

function alignmentBlocks(query: string, reference: string, width = 72) {
  const output: Array<{ query: string; reference: string; match: string; offset: number }> = [];
  for (let index = 0; index < Math.max(query.length, reference.length); index += width) {
    const queryBlock = query.slice(index, index + width);
    const referenceBlock = reference.slice(index, index + width);
    output.push({ query: queryBlock, reference: referenceBlock, match: matchLine(queryBlock, referenceBlock), offset: index });
  }
  return output;
}

function percent(value: string): string {
  const parsed = Number(value);
  return Number.isFinite(parsed) && value !== "" ? `${(parsed * 100).toFixed(1)}%` : "—";
}

function coordinate(value: number, mode: AlignmentMode, frameOneBased: number): string {
  if (mode === "nt" || !frameOneBased) return value.toLocaleString();
  const frame = frameOneBased - 1;
  if (value - 1 < frame) return "pre-frame";
  return `${Math.floor((value - 1 - frame) / 3) + 1} aa`;
}

function FeatureTracks({ row, mode }: { row: AirrRow; mode: AlignmentMode }) {
  const tracks = useMemo(() => buildTrackFeatures(row), [row]);
  const length = row.sequence?.length ?? 0;
  const frame = Number(row.sequence_frame);
  const middle = Math.max(1, Math.ceil(length / 2));
  return (
    <section className="layered-tracks" aria-label="Query feature tracks">
      <div className="track-ruler"><span>{coordinate(1, mode, frame)}</span><i /><span>{coordinate(middle, mode, frame)}</span><i /><span>{coordinate(length, mode, frame)}</span></div>
      <div className="feature-track query-track"><b>Query</b><div><i /></div></div>
      <div className="feature-track region-track"><b>IMGT regions</b><div>{tracks.regions.map((item) => <span className={item.key.startsWith("cdr") ? "cdr" : "fwr"} key={item.key} style={{ left: `${item.left}%`, width: `${item.width}%` }} title={`${item.label}: ${item.start}–${item.end} nt`}>{item.label}</span>)}</div></div>
      <div className="feature-track germline-track"><b>Germline hits</b><div>{tracks.segments.map((item) => <span className={item.key} key={item.key} style={{ left: `${item.left}%`, width: `${item.width}%` }} title={`${row[`${item.key}_call`]}: ${item.start}–${item.end} nt`}>{item.label}</span>)}</div></div>
      {!tracks.regions.length && <p className="track-warning">No validated FWR/CDR coordinates were reported for the selected V call.</p>}
    </section>
  );
}

function RegionSequences({ row, mode }: { row: AirrRow; mode: AlignmentMode }) {
  return (
    <div className="layer-region-grid">
      {REGIONS.map((region) => {
        const value = row[mode === "aa" ? `${region}_aa` : region] || "";
        const start = Number(row[`${region}_start`]);
        const end = Number(row[`${region}_end`]);
        return <article className={region.startsWith("cdr") ? "cdr" : "fwr"} key={region}><span>{region.toUpperCase()}</span>{value ? <ColoredSequence sequence={value} alphabet={mode} /> : <code>—</code>}<small>{start && end ? `${coordinate(start, mode, Number(row.sequence_frame))}–${coordinate(end, mode, Number(row.sequence_frame))}` : "not assigned"}</small></article>;
      })}
    </div>
  );
}

function CompositeAlignment({ row, mode }: { row: AirrRow; mode: AlignmentMode }) {
  const query = mode === "aa" ? row.sequence_alignment_aa : row.sequence_alignment;
  const reference = mode === "aa" ? row.germline_alignment_aa : row.germline_alignment;
  const rendered = alignmentBlocks(query || "", reference || "");
  return (
    <article className="composite-alignment">
      <header><div><span className="section-kicker">{row.d2_call ? "Baseline single-D composite" : "Combined V(D)J layer"}</span><h4>Query ↔ composite germline</h4></div><small>{mode === "aa" ? `Shared biological frame ${row.sequence_frame ? `+${row.sequence_frame}` : "not available"}` : "AIRR stitched nucleotide alignment"}</small></header>
      {rendered.length ? rendered.map((block) => <div className="alignment-block" key={block.offset}>
        <div><span>query</span><ColoredSequence sequence={block.query} alphabet={mode} /></div>
        <div className="match-row"><span /><code>{block.match}</code></div>
        <div><span>V(D)J ref</span><ColoredSequence sequence={block.reference} alphabet={mode} /></div>
      </div>) : <p className="empty-alignment">{mode === "aa" ? "A shared coding frame could not be established for this record." : "No stitched V(D)J alignment was reported."}</p>}
    </article>
  );
}

function SegmentAlignment({ row, segment, mode }: {
  row: AirrRow;
  segment: typeof ALIGNMENT_SEGMENTS[number];
  mode: AlignmentMode;
}) {
  const call = row[`${segment.key}_call`];
  const nucleotideQuery = row[`${segment.key}_sequence_alignment`] ?? "";
  const nucleotideReference = row[`${segment.key}_germline_alignment`] ?? "";
  const protein = useMemo(() => mode === "aa" ? biologicalSegmentAlignment(
    nucleotideQuery,
    nucleotideReference,
    Number(row[`${segment.key}_sequence_start`]),
    Number(row.sequence_frame),
  ) : null, [mode, nucleotideQuery, nucleotideReference, row, segment.key]);
  if (!call) return null;
  const query = mode === "aa" ? protein?.query ?? "" : nucleotideQuery;
  const reference = mode === "aa" ? protein?.reference ?? "" : nucleotideReference;
  const rendered = alignmentBlocks(query, reference);
  return (
    <article className="alignment-card" style={{ "--segment-color": segment.color } as React.CSSProperties}>
      <header>
        <span className="segment-chip">{segment.label}</span>
        <div><strong>{call}</strong><small>{percent(row[`${segment.key}_identity`])} identity · {row[`${segment.key}_cigar`] || "no CIGAR"}</small></div>
        <dl>
          <div><dt>Query</dt><dd>{row[`${segment.key}_sequence_start`] || "—"}–{row[`${segment.key}_sequence_end`] || "—"}</dd></div>
          <div><dt>Germline</dt><dd>{row[`${segment.key}_germline_start`] || "—"}–{row[`${segment.key}_germline_end`] || "—"}</dd></div>
          <div><dt>Score</dt><dd>{row[`${segment.key}_score`] || "—"}</dd></div>
        </dl>
      </header>
      {rendered.length ? rendered.map((block) => <div className="alignment-block" key={block.offset}>
        <div><span>query</span><ColoredSequence sequence={block.query} alphabet={mode} /></div>
        <div className="match-row"><span /><code>{block.match}</code></div>
        <div><span>{segment.label.toLowerCase()} ref</span><ColoredSequence sequence={block.reference} alphabet={mode} /></div>
      </div>) : <p className="empty-alignment">{mode === "aa" ? "No frame-consistent translated alignment is available." : "No aligned sequence was reported for this segment."}</p>}
    </article>
  );
}

export function AlignmentViewer({ row, mode, onMode }: {
  row: AirrRow;
  mode: AlignmentMode;
  onMode: (mode: AlignmentMode) => void;
}) {
  const hasRegions = Boolean(row.region_definition && row.fwr1_start && row.fwr3_end);
  return (
    <section className="alignment-viewer">
      <div className="alignment-toolbar">
        <div><span className="section-kicker">Layered sequence annotation</span><h3>Query, IMGT regions, and {row.d2_call ? "V(D)(D)J" : "V(D)J"} germlines</h3></div>
        <div className="alignment-actions"><button className="copy-fasta" type="button" onClick={() => void navigator.clipboard.writeText(`>${row.sequence_id}\n${row.sequence}\n`)}>Copy query FASTA</button><div className="mode-toggle" role="group" aria-label="Alignment alphabet"><button className={mode === "nt" ? "active" : ""} type="button" onClick={() => onMode("nt")}>Nucleotide</button><button className={mode === "aa" ? "active" : ""} type="button" onClick={() => onMode("aa")}>Amino acid</button></div></div>
      </div>
      <div className={`annotation-status ${hasRegions ? "valid" : "incomplete"}`}><span>{hasRegions ? "Validated region map" : "Region map unavailable"}</span><p>{hasRegions ? `${row.region_definition} boundaries · V: ${row.v_annotation_source || "unspecified"} · J: ${row.j_annotation_source || "unspecified"}${row.sequence_frame ? ` · query frame +${row.sequence_frame}` : ""}` : "No V-region delineation passed validation for the selected allele. A closest-relative transfer is accepted only when the sequence spans the mapped IMGT intervals."}</p></div>
      {mode === "aa" && <p className="alignment-explainer">All translated layers use the single rearrangement frame anchored by the mapped V cysteine and J F/W-G motif. Frames are not optimized independently per segment.</p>}
      {row.d2_call && <p className="alignment-explainer">The composite layer is the unchanged standard single-D call. The D1 and D2 cards and feature track are the opt-in double-D screening model; the main AIRR export is not rewritten.</p>}
      <FeatureTracks row={row} mode={mode} />
      <RegionSequences row={row} mode={mode} />
      <CompositeAlignment row={row} mode={mode} />
      <div className="segment-evidence-heading"><span className="section-kicker">Per-segment alignments</span><h4>Selected germline calls</h4></div>
      <div className="alignment-stack">{ALIGNMENT_SEGMENTS.map((segment) => <SegmentAlignment key={segment.key} row={row} segment={segment} mode={mode} />)}</div>
    </section>
  );
}
