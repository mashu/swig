import { useState } from "react";

import { CommitNumberInput } from "../commit-number-input.tsx";
import type { CompiledReferences } from "../reference-pack.ts";
import { ALLELE_POOL_SCOPE_LABELS, type DatasetScope } from "../study-design.ts";
import { AlleleAssignmentShiftChart, ReferenceKernelInspector } from "./diagnostic-views.tsx";
import type { AlleleReassignmentPolicy, AlleleRefinementOptions, AlleleRefinementResult, RefinementSegment, SegmentRefinementResult } from "./types.ts";
import { adaptiveNeighbourOdds } from "./evidence.ts";

interface Props {
  references: CompiledReferences;
  options: AlleleRefinementOptions;
  onOptionsChange: (options: AlleleRefinementOptions) => void;
  result: AlleleRefinementResult | null;
  applied: boolean;
  reassignmentPolicy: AlleleReassignmentPolicy;
  onReassignmentPolicyChange: (policy: AlleleReassignmentPolicy) => void;
  applyMinimumPosterior: number;
  onApplyMinimumPosteriorChange: (value: number) => void;
  busy: boolean;
  progress: { processed: number; total: number; phase: string } | null;
  onRun: () => void;
  onCancel: () => void;
  onApply: () => void;
  onReset: () => void;
  onDownloadModel: () => void;
  onDownloadSidecar: () => void;
  onDownloadAirr: () => void;
}

const SEGMENT_LABELS: Record<RefinementSegment, string> = {
  V: "V alleles",
  D: "D alleles · experimental",
  J: "J alleles",
};

export function AlleleRefinementPanel({
  references, options, onOptionsChange, result, applied, reassignmentPolicy, onReassignmentPolicyChange,
  applyMinimumPosterior, onApplyMinimumPosteriorChange, busy, progress, onRun, onCancel, onApply, onReset,
  onDownloadModel, onDownloadSidecar, onDownloadAirr,
}: Props) {
  const [advancedOpen, setAdvancedOpen] = useState(false);
  const update = <K extends keyof AlleleRefinementOptions>(key: K, value: AlleleRefinementOptions[K]) => onOptionsChange({ ...options, [key]: value });
  const toggleSegment = (segment: RefinementSegment) => {
    const selected = options.segments.includes(segment)
      ? options.segments.filter((value) => value !== segment)
      : [...options.segments, segment];
    update("segments", selected);
  };
  const segmentResults = Object.values(result?.segments ?? {}).filter((segment): segment is SegmentRefinementResult => Boolean(segment));
  const modeledRows = segmentResults.reduce((sum, segment) => sum + segment.modeledRows, 0);
  const changedRows = segmentResults.reduce((sum, segment) => sum + segment.changedMapRows, 0);
  const nonZeros = segmentResults.reduce((sum, segment) => sum + segment.matrixNonZeros, 0);
  const models = segmentResults.flatMap((segment) => segment.models);
  const progressFraction = progress ? Math.max(0, Math.min(1, progress.processed / Math.max(1, progress.total))) : 0;
  const leakageExamples = [0, 0.05, 0.1].map((shm) => {
    const odds = adaptiveNeighbourOdds(shm, options);
    return `${Math.round(shm * 100)}% SHM → ${(100 * odds / (1 + odds)).toFixed(2)}%`;
  }).join(" · ");

  return <>
    <div className="allele-refinement-config">
      <div className="control-grid four">
        <label title="The legacy continuous Dirichlet model remains available unchanged. The hurdle model gives excluded alleles exact zero usage and uses fast active-set pruning."><span>Repertoire model</span><select value={options.model} onChange={(event) => update("model", event.target.value as AlleleRefinementOptions["model"])}><option value="dirichlet">Continuous Dirichlet mixture · legacy default</option><option value="active-set">Fast hurdle active set · exact zeroes</option></select></label>
        <label title="Default: one independent fit per donor, combining that donor's samples, timepoints, and compartments. Cohort and entire-study modes are explicit cross-participant overrides."><span>Pooling boundary</span><select value={options.scope} onChange={(event) => update("scope", event.target.value as DatasetScope)}>{Object.entries(ALLELE_POOL_SCOPE_LABELS).map(([value, label]) => <option value={value} key={value}>{label}</option>)}</select></label>
        <fieldset className="allele-segment-picker"><legend title="V and J are enabled by default. D inference is less identifiable because D segments are short and heavily trimmed.">Reference segments</legend>{(["V", "J", "D"] as RefinementSegment[]).map((segment) => <label key={segment}><input type="checkbox" checked={options.segments.includes(segment)} onChange={() => toggleSegment(segment)} /><span>{SEGMENT_LABELS[segment]}</span></label>)}</fieldset>
        <label title="Unique weighting prevents expanded clones or PCR abundance from acting as independent genotype evidence. Abundance weighting is available for deliberate usage estimation."><span>Record weighting</span><select value={options.weighting} onChange={(event) => update("weighting", event.target.value as AlleleRefinementOptions["weighting"])}><option value="unique">One vote per active record · default</option><option value="abundance">Weight by duplicate_count</option></select></label>
      </div>
      <p className={`allele-pooling-boundary-note ${options.scope === "cohort" || options.scope === "global" ? "override" : ""}`}><span>{options.scope === "cohort" || options.scope === "global" ? "Cross-donor override" : "Participant-safe default"}</span>{options.scope === "subject" ? "Each donor is fitted independently; all samples, timepoints, and compartments carrying that donor ID contribute to the same fit. Evidence never crosses participants." : options.scope === "cohort" || options.scope === "global" ? "This setting deliberately pools evidence across participant IDs. Use it only when that is scientifically intended." : "Fits are narrower than donor level and therefore cannot cross participant IDs."}</p>
      <div className="allele-refinement-action"><div><strong>{result ? `${result.activeRecords.toLocaleString()} active records modeled` : "No repertoire posterior fitted"}</strong><small>{options.model === "active-set" ? "The hurdle fit assigns exact zero usage to excluded alleles; retained alleles use the same sparse evidence kernel and per-record posterior projection." : "The original continuous Dirichlet mixture is selected. Explicit co-optimal calls enter with equal local weight."}</small></div>{progress && busy ? <button className="post-cancel" type="button" onClick={onCancel}>Cancel fit</button> : <button className="post-primary" type="button" disabled={busy || !options.segments.length} onClick={onRun}>{result ? "Refit repertoire allele model" : "Fit repertoire allele model"}</button>}</div>
      <details className="post-advanced" onToggle={(event) => setAdvancedOpen(event.currentTarget.open)}><summary>Advanced evidence-kernel and model settings</summary><div><div className="control-grid four allele-refinement-advanced">
        <label title="Irreducible relative evidence odds for a one-nucleotide reference neighbour when estimated SHM is zero. This is an assignment-model leakage term, not a sequencing-error estimate."><span>Zero-SHM neighbour odds</span><CommitNumberInput min="0" max="0.5" step="0.001" value={options.baselineNeighbourOdds} onCommit={(value) => update("baselineNeighbourOdds", value)} /></label>
        <label title="Multiplier on the mechanistic mu/[3(1-mu)] contribution for a specific nucleotide substitution at a diagnostic allele position. Set to zero for SHM-independent leakage."><span>SHM sensitivity</span><CommitNumberInput min="0" max="10" step="0.1" value={options.shmLeakageSensitivity} onCommit={(value) => update("shmLeakageSensitivity", value)} /></label>
        <label title="Upper bound on one-SNP neighbour evidence odds for highly mutated or poorly matched reads."><span>Maximum neighbour odds</span><CommitNumberInput min="0.001" max="0.99" step="0.01" value={options.maximumNeighbourOdds} onCommit={(value) => update("maximumNeighbourOdds", value)} /></label>
        <label title="The best-reference V mutation fraction is clamped at this value before the SHM leakage term is calculated."><span>SHM estimate cap</span><CommitNumberInput min="0" max="0.9" step="0.01" value={options.maximumShm} onCommit={(value) => update("maximumShm", value)} /></label>
        <label title="Maximum bounded nucleotide edit distance from any reported or retained alternative reference. Candidate rows remain sparse."><span>Neighbour edit radius</span><CommitNumberInput min="0" max="5" step="1" value={options.neighbourRadius} onCommit={(value) => update("neighbourRadius", Math.floor(value))} /></label>
        {options.model === "dirichlet" && <label title="Symmetric Dirichlet pseudo-count for every locus-matched reference node. Nodes absent from all sparse candidate rows retain only this prior mass and are not materialized per record."><span>Dirichlet alpha / allele</span><CommitNumberInput min="0.000001" max="100" step="0.01" value={options.alphaPerAllele} onCommit={(value) => update("alphaPerAllele", value)} /></label>}
        {options.model === "active-set" && <>
          <label title="Prior probability that each database allele contributes measurable rearrangements in this fitted donor/locus/segment pool. This is repertoire activity, not genomic presence."><span>Active prior / allele</span><CommitNumberInput min="0.000001" max="0.999999" step="0.01" value={options.activeSetPriorActiveFraction} onCommit={(value) => update("activeSetPriorActiveFraction", Math.max(0.000001, Math.min(0.999999, value)))} /></label>
          <label title="An allele is retained in the exact non-zero active set when its approximate posterior inclusion probability reaches this value. Remapping confidence remains a separate downstream decision."><span>Inclusion posterior threshold</span><CommitNumberInput min="0" max="1" step="0.01" value={options.activeSetInclusionThreshold} onCommit={(value) => update("activeSetInclusionThreshold", Math.max(0, Math.min(1, value)))} /></label>
          <label title="Shape of the positive-frequency gamma slab. Values below one retain substantial density in the low-frequency tail; smaller values make that tail heavier."><span>Positive-use tail shape</span><CommitNumberInput min="0.05" max="10" step="0.05" value={options.activeSetTailShape} onCommit={(value) => update("activeSetTailShape", Math.max(0.05, value))} /></label>
          <label title="Lowest positive frequency represented by numerical quadrature. This is not a minimum retained count or biological frequency cutoff."><span>Tail integration floor</span><CommitNumberInput min="0.000000000001" max="0.05" step="0.000001" value={options.activeSetFrequencyFloor} onCommit={(value) => update("activeSetFrequencyFloor", Math.max(1e-12, value))} /></label>
          <label title="Number of log-spaced positive-frequency points in each one-dimensional inclusion test. Twelve is normally sufficient because each test is smooth."><span>Frequency quadrature points</span><CommitNumberInput min="4" max="64" step="1" value={options.activeSetQuadraturePoints} onCommit={(value) => update("activeSetQuadraturePoints", Math.max(4, Math.min(64, Math.floor(value))))} /></label>
          <label title="Maximum parallel backward-pruning passes. Each pass refits frequencies only over the surviving sparse active set."><span>Active-set pruning sweeps</span><CommitNumberInput min="1" max="50" step="1" value={options.activeSetMaxSweeps} onCommit={(value) => update("activeSetMaxSweeps", Math.max(1, Math.min(50, Math.floor(value))))} /></label>
        </>}
        <label title="Scale converting aligner score differences for retained alternative hits into relative local evidence. Smaller values penalize lower scores more strongly."><span>Alternative score temperature</span><CommitNumberInput min="0.001" max="100" step="0.1" value={options.alternativeScoreTemperature} onCommit={(value) => update("alternativeScoreTemperature", value)} /></label>
        <label title="Fallback local weight for a retained alternative that lacks a numeric alignment score."><span>Unscored alternative weight</span><CommitNumberInput min="0" max="1" step="0.01" value={options.unscoredAlternativeWeight} onCommit={(value) => update("unscoredAlternativeWeight", value)} /></label>
        <label title="Hard memory guard for unusually dense reference neighbourhoods. Candidates are retained in descending local evidence order."><span>Candidate cap / record</span><CommitNumberInput min="1" max="10000" step="1" value={options.maxCandidatesPerRow} onCommit={(value) => update("maxCandidatesPerRow", Math.floor(value))} /></label>
        <label title={options.model === "active-set" ? "Maximum sparse EM frequency updates within each active-set refit." : "Maximum coordinate-ascent updates for each independent pooling group."}><span>Maximum iterations</span><CommitNumberInput min="1" max="10000" step="10" value={options.maxIterations} onCommit={(value) => update("maxIterations", Math.floor(value))} /></label>
        <label title={options.model === "active-set" ? "Stop a sparse active-frequency refit when the maximum absolute mixture-weight change falls below this value." : "Stop when the maximum relative change in any Dirichlet parameter falls below this value."}><span>Convergence tolerance</span><CommitNumberInput min="0.000000001" max="0.1" step="0.000001" value={options.convergenceTolerance} onCommit={(value) => update("convergenceTolerance", value)} /></label>
      </div><div className="algorithm-note"><strong>{options.model === "active-set" ? "Fast hurdle active set" : "Implied one-SNP neighbour probability"}</strong><span>{options.model === "active-set" ? `Sparse EM → one-dimensional inclusion tests → protected parallel pruning. Excluded alleles have exact zero usage. Local evidence examples: ${leakageExamples}.` : `${leakageExamples}. Displayed probabilities are odds/(1+odds); distance-d evidence uses the per-edit odds to the dth power.`}</span></div>{advancedOpen && <ReferenceKernelInspector references={references} options={options} />}</div></details>
    </div>
    {progress && <div className="allele-refinement-progress" role="status"><div><strong>{progress.phase}</strong><span>{progress.processed.toLocaleString()} / {progress.total.toLocaleString()}</span></div><progress max="1" value={progressFraction} /></div>}
    {result && <div className="allele-refinement-results">
      <div className="post-stat-grid"><article><span>Independent models</span><strong>{models.length.toLocaleString()}</strong></article><article><span>Modeled segment-rows</span><strong>{modeledRows.toLocaleString()}</strong></article><article><span>Repertoire MAP differs</span><strong>{changedRows.toLocaleString()}</strong></article><article><span>{options.model === "active-set" ? "Retained active alleles" : "Sparse nonzeros"}</span><strong>{options.model === "active-set" ? models.reduce((sum, model) => sum + (model.activeAlleles ?? 0), 0).toLocaleString() : nonZeros.toLocaleString()}</strong></article></div>
      <div className="allele-apply-row"><div><strong>{applied ? "The selected reassignment policy is active downstream" : "Original AIRR calls remain active downstream"}</strong><small>{reassignmentPolicy === "best" ? "Best posterior: every modeled record receives its posterior MAP allele; unmodeled records retain their immutable original call." : `Confidence gated: a modeled record receives its posterior MAP allele only at ≥ ${(applyMinimumPosterior * 100).toFixed(0)}% posterior confidence; otherwise its immutable original AIRR call is retained.`} Applying, restoring, or changing an active policy clears later collapse/filter results so they can be recomputed from the correct calls.</small></div><label title="Choose whether every modeled record receives its posterior MAP call or whether low-confidence records retain their original AIRR call."><span>Reassignment policy</span><select value={reassignmentPolicy} onChange={(event) => onReassignmentPolicyChange(event.target.value as AlleleReassignmentPolicy)}><option value="confidence">Best posterior if confidence passes · default</option><option value="best">Best posterior for every modeled record</option></select></label>{reassignmentPolicy === "confidence" && <label title="The maximum per-record posterior probability required before its MAP allele replaces the original AIRR call."><span>Minimum MAP confidence</span><CommitNumberInput min="0" max="1" step="0.01" value={applyMinimumPosterior} onCommit={onApplyMinimumPosteriorChange} /></label>}<button type="button" className="post-primary" disabled={busy} onClick={onApply}>{reassignmentPolicy === "best" ? "Apply best posterior calls" : "Apply confidence-gated calls"}</button>{applied && <button type="button" disabled={busy} onClick={onReset}>Restore original calls</button>}</div>
      {result.warnings.map((warning) => <p className="scientific-note warning" key={warning}><span>!</span>{warning}</p>)}
      <div className="result-actions"><button type="button" onClick={onDownloadModel}>Model summary</button><button type="button" onClick={onDownloadSidecar}>Full sparse posterior sidecar</button><button type="button" onClick={onDownloadAirr}>AIRR with policy-selected calls</button></div>
      {models.length > 0 && <AlleleAssignmentShiftChart results={segmentResults} reassignmentPolicy={reassignmentPolicy} minimumPosterior={applyMinimumPosterior} weighting={options.weighting} />}
    </div>}
  </>;
}
