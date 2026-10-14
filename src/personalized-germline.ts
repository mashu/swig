import { boundaryLogLikelihood, downstreamJunction, TERMINAL_WINDOW, INITIAL_JUNCTION_BASES } from "./shm-model/boundary.ts";
import { hs5fRate } from "./shm-model/hs5f.ts";
import { substitutionCompatible } from "./shm-model/alignment.ts";
import { AlignedRateCalibration, baseAt, discoveryCost, testHaplotype } from "./shm-model/discovery.ts";
import { parseReferenceFasta, serializeReferenceFasta, type ReferenceFastaRecord } from "./reference-fasta.ts";
import { orientedSequence } from "./nucleotide.ts";

export interface PersonalizedGermlineOptions {
  minimumAlignedBases: number;
  /** Known alleles are compared only inside the called gene and this Hamming radius. */
  maximumKnownAlleleSnps: number;
  maximumNovelSnps: number;
  minimumNovelSupport: number;
  minimumNovelFraction: number;
  maximumNovelCandidatesPerGene: number;
  /** Additional gain beyond the BIC penalty required for each greedy addition. */
  minimumLogEvidenceGain: number;
  sequencingErrorRate: number;
  maximumActiveAllelesPerGene: number;
  maximumIterations: number;
  convergenceTolerance: number;
}

export const DEFAULT_PERSONALIZED_GERMLINE_OPTIONS: PersonalizedGermlineOptions = {
  minimumAlignedBases: 120,
  maximumKnownAlleleSnps: 8,
  maximumNovelSnps: 6,
  minimumNovelSupport: 4,
  minimumNovelFraction: 0,
  maximumNovelCandidatesPerGene: 64,
  minimumLogEvidenceGain: 0,
  sequencingErrorRate: 0.001,
  maximumActiveAllelesPerGene: 8,
  maximumIterations: 100,
  convergenceTolerance: 1e-7,
};

export interface PersonalizedSubstitution {
  position: number;
  reference: string;
  alternate: string;
}

export interface PersonalizedGermlineAllele {
  id: string;
  names: string[];
  sequence: string;
  known: boolean;
  parentAllele: string;
  substitutions: PersonalizedSubstitution[];
  frequency: number;
  expectedLineages: number;
  localBestLineages: number;
  directSupport: number;
  directCoverage: number;
  selectionGain: number | null;
  discoveryGain?: number | null;
}

export interface PersonalizedGermlineGeneResult {
  gene: string;
  representativeLineages: number;
  testedCandidates: number;
  proposedNovelCandidates: number;
  activeAlleles: PersonalizedGermlineAllele[];
  relativeLogLikelihood: number;
  iterations: number;
  converged: boolean;
  proposalEvidence?: Array<{id:string;gain:number;support:number;coverage:number;baseline:string;selected:boolean;reason?:string}>;
}

export interface PersonalizedGermlinePool {
  id: string;
  subjectId: string;
  locus: string;
  representativeLineages: number;
  genes: PersonalizedGermlineGeneResult[];
}

export interface PersonalizedGermlineDashboard {
  version: 1;
  mode: "lowest-current-v-shm-per-lineage";
  contextModel: "hs5f-hurdle-profiled-junction";
  options: PersonalizedGermlineOptions;
  inputRecords: number;
  eligibleRecords: number;
  assignedLineages: number;
  representativeLineages: number;
  skippedLineages: number;
  pools: PersonalizedGermlinePool[];
  testedCandidates: number;
  activeKnownAlleles: number;
  activeNovelAlleles: number;
  proposalTruncations: number;
  warnings: string[];
}

export interface ReferenceNode {
  index: number;
  locus: string;
  gene: string;
  names: string[];
  headers: string[];
  sequence: string;
}

export interface Observation {
  ordinal: number;
  lineageId: number;
  subjectId: string;
  locus: string;
  gene: string;
  parent: ReferenceNode;
  shmRate: number;
  alignedBases: number;
  positions: Uint16Array;
  query: string;
  terminal?: {start:number;query:string};
  junction?: string;
}

interface Candidate {
  id: string;
  names: string[];
  sequence: string;
  known: boolean;
  parent: ReferenceNode;
  substitutions: PersonalizedSubstitution[];
  directSupport: number;
  directCoverage: number;
  searchCost?: number;
  proposalPriority?: number;
}

interface CandidateProposal {
  candidate: Candidate;
  signature: string;
}

interface MixtureFit {
  theta: Float64Array;
  mixture: Float64Array;
  logLikelihood: number;
  counts: Float64Array;
  iterations: number;
  converged: boolean;
}

const BASES = /^[ACGT]$/;
const W = new Set(["A", "T"]);
const R = new Set(["A", "G"]);
const Y = new Set(["C", "T"]);
const S = new Set(["C", "G"]);

function cleanAlignment(value: string): string {
  return value.toUpperCase().replaceAll("U", "T").replaceAll(".", "-").replace(/[^ACGTN-]/g, "N");
}

function cleanReference(value: string): string {
  return value.toUpperCase().replaceAll("U", "T").replace(/[.\-]/g, "").replace(/[^ACGTN]/g, "N");
}

export function personalizedAlleleGene(value: string): string {
  const token = value.trim().split(/\s+/, 1)[0] ?? "";
  const star = token.indexOf("*");
  return star >= 0 ? token.slice(0, star) : token;
}

function referenceLocus(value: string): string {
  return value.toUpperCase().match(/^(IGH|IGK|IGL|TRA|TRB|TRD|TRG)/)?.[1] ?? "";
}

function callTokens(value: string): string[] {
  return [...new Set(value.split(",").map((item) => item.trim()).filter(Boolean))];
}

function hammingDistance(left: string, right: string, maximum = Number.POSITIVE_INFINITY): number {
  if (left.length !== right.length) return Number.POSITIVE_INFINITY;
  let distance = 0;
  for (let index = 0; index < left.length; index += 1) {
    if (left[index] === right[index]) continue;
    distance += 1;
    if (distance > maximum) return distance;
  }
  return distance;
}

/**
 * Legacy motif annotation helper; personalized inference now uses HS5F. It uses the classic AID WRCY/RGYW hot spots,
 * SYC/GRS cold spots, and WA/TW polymerase-eta hot spots inside a 5-mer
 * window. It is deliberately labelled as a motif prior, not as S5F.
 */
export function shmContextMutability(sequence: string, zeroBasedPosition: number): number {
  const base = sequence[zeroBasedPosition] ?? "N";
  const minus2 = sequence[zeroBasedPosition - 2] ?? "N";
  const minus1 = sequence[zeroBasedPosition - 1] ?? "N";
  const plus1 = sequence[zeroBasedPosition + 1] ?? "N";
  const plus2 = sequence[zeroBasedPosition + 2] ?? "N";
  if (base === "C" && W.has(minus2) && R.has(minus1) && Y.has(plus1)) return 5;
  if (base === "G" && R.has(minus1) && Y.has(plus1) && W.has(plus2)) return 5;
  if (base === "C" && S.has(minus2) && Y.has(minus1)) return 0.25;
  if (base === "G" && R.has(plus1) && S.has(plus2)) return 0.25;
  if (base === "A" && W.has(minus1)) return 2;
  if (base === "T" && W.has(plus1)) return 2;
  return 1;
}

function referenceNodes(fasta: string): { nodes: ReferenceNode[]; byName: Map<string, ReferenceNode> } {
  const byKey = new Map<string, ReferenceNode>();
  const byName = new Map<string, ReferenceNode>();
  for (const record of parseReferenceFasta(fasta)) {
    const sequence = cleanReference(record.sequence);
    const locus = referenceLocus(record.name);
    const gene = personalizedAlleleGene(record.name);
    if (!sequence || !locus || !gene) continue;
    const key = `${locus}\u0000${gene}\u0000${sequence}`;
    let node = byKey.get(key);
    if (!node) {
      node = { index: byKey.size, locus, gene, names: [], headers: [], sequence };
      byKey.set(key, node);
    }
    if (!node.names.includes(record.name)) node.names.push(record.name);
    node.headers.push(record.header);
    byName.set(record.name, node);
  }
  const nodes = [...byKey.values()];
  nodes.forEach((node, index) => {
    node.index = index;
    node.names.sort((left, right) => left.localeCompare(right, undefined, { numeric: true }));
  });
  return { nodes, byName };
}

function parseObservation(
  row: Record<string, string>,
  ordinal: number,
  lineageId: number,
  byName: ReadonlyMap<string, ReferenceNode>,
  options: PersonalizedGermlineOptions,
): Observation | null {
  if (!(lineageId > 0)) return null;
  let parent = callTokens(row.v_call ?? "").map((call) => byName.get(call)).find((value): value is ReferenceNode => Boolean(value));
  if (!parent) return null;
  const queryAlignment = cleanAlignment(row.v_sequence_alignment ?? "");
  const germlineAlignment = cleanAlignment(row.v_germline_alignment ?? "");
  if (!queryAlignment || !germlineAlignment) return null;
  const reportedStart = Math.floor(Number(row.v_germline_start));
  const alignmentCost=(node:ReferenceNode)=>{
    let p=Math.max(0,reportedStart-1),cost=0;
    for(const b of germlineAlignment){if(b==="-")continue;if(BASES.test(b)&&node.sequence[p]!==b)cost++;p++;}
    return cost;
  };
  parent=callTokens(row.v_call??"").map(name=>byName.get(name)).filter((node):node is ReferenceNode=>Boolean(node)).sort((a,b)=>alignmentCost(a)-alignmentCost(b))[0]??parent;
  let position = Number.isFinite(reportedStart) && reportedStart > 0 ? reportedStart - 1 : 0;
  const positions: number[] = [];
  const rawPositions = new Map<number,number>();
  let rawPosition=Math.max(0,Math.floor(Number(row.v_sequence_start)||1)-1);
  const query: string[] = [];
  let currentAligned = 0;
  let currentMismatches = 0;
  const length = Math.min(queryAlignment.length, germlineAlignment.length);
  for (let column = 0; column < length; column += 1) {
    const germlineBase = germlineAlignment[column];
    if (germlineBase !== "-") position += 1;
    const queryBase = queryAlignment[column];
    if (queryBase !== "-") rawPosition += 1;
    if (germlineBase !== "-" && queryBase !== "-") rawPositions.set(position,rawPosition);
    if (BASES.test(queryBase) && BASES.test(germlineBase)) {
      currentAligned += 1;
      if (queryBase !== germlineBase) currentMismatches += 1;
    }
    if (germlineBase === "-" || !BASES.test(queryBase) || position < 1 || position > parent.sequence.length || !BASES.test(parent.sequence[position - 1])) continue;
    positions.push(position);
    query.push(queryBase);
  }
  if (positions.length < options.minimumAlignedBases) return null;
  const cutoff=parent.sequence.length-TERMINAL_WINDOW;
  const raw=orientedSequence(row.sequence,row.rev_comp);
  const anchor=rawPositions.get(cutoff);
  // Downstream D/J templated bases cannot be evidence for a novel V end.
  const downstreamStart=Number(row.d_sequence_start)>0?Number(row.d_sequence_start)-1:Number(row.j_sequence_start)>0?Number(row.j_sequence_start)-1:raw.length;
  const terminal=raw&&anchor!==undefined?{start:cutoff+1,query:raw.slice(anchor,Math.max(anchor,Math.min(anchor+TERMINAL_WINDOW,downstreamStart)))}:undefined;
  // Reconstruct both matching and nonmatching downstream bases from the same anchor.
  // Without raw sequence, terminal bases are missing evidence, not committed V matches.
  while(positions.length&&positions[positions.length-1]>cutoff){positions.pop();query.pop();}
  if(terminal)for(let i=0;i<terminal.query.length;i++)if(BASES.test(terminal.query[i])){positions.push(terminal.start+i);query.push(terminal.query[i]);}

  const subjectId = (row.subject_id || "unassigned-subject").trim() || "unassigned-subject";
  return {
    ordinal,
    lineageId,
    subjectId,
    locus: parent.locus,
    gene: parent.gene,
    parent,
    shmRate: currentAligned ? currentMismatches / currentAligned : 1,
    alignedBases: positions.length,
    positions: Uint16Array.from(positions),
    query: query.join(""),
    terminal,
    junction:raw&&anchor!==undefined?downstreamJunction(raw,anchor+TERMINAL_WINDOW,Number(row.d_sequence_start),Number(row.j_sequence_start)):undefined,
  };
}

function queryBaseAt(observation: Observation, position: number): string | null {
  let low = 0;
  let high = observation.positions.length - 1;
  while (low <= high) {
    const middle = (low + high) >>> 1;
    const found = observation.positions[middle];
    if (found === position) return observation.query[middle] ?? null;
    if (found < position) low = middle + 1;
    else high = middle - 1;
  }
  return null;
}

function candidateSubstitutions(parent: ReferenceNode, sequence: string): PersonalizedSubstitution[] {
  const result: PersonalizedSubstitution[] = [];
  for (let index = 0; index < Math.min(parent.sequence.length, sequence.length); index += 1) {
    if (parent.sequence[index] !== sequence[index]) result.push({ position: index + 1, reference: parent.sequence[index], alternate: sequence[index] });
  }
  return result;
}

function candidateId(parent: ReferenceNode, substitutions: readonly PersonalizedSubstitution[]): string {
  return `${parent.names[0]}__SWIGP_${substitutions.map((item) => `${item.reference}${item.position}${item.alternate}`).join("_")}`.replace(/[^A-Za-z0-9_.|*+\-]/g, "_");
}

function proposalEvidence(observations: readonly Observation[], substitutions: readonly PersonalizedSubstitution[], counts:ReadonlyMap<Observation,{k:number;w:number;n:number}>, calibration:AlignedRateCalibration): { support: number; coverage: number; priority:number } {
  let support = 0;
  let coverage = 0,priority=0;
  const interior=substitutions.filter(c=>c.position<=observations[0].parent.sequence.length-TERMINAL_WINDOW);
  const omittedWeight=interior.reduce((sum,c)=>sum+hs5fRate(observations[0].parent.sequence,c.position-1),0);
  const exposureCache=new Map<string,number>();
  for (const observation of observations) {
    let complete = true;
    let alternate = true;
    for (const substitution of substitutions) {
      const base = queryBaseAt(observation, substitution.position);
      if (!base) {
        complete = false;
        alternate = false;
        break;
      }
      if (base !== substitution.alternate) alternate = false;
    }
    if (complete) coverage += 1;
    if (complete && alternate){
      support += 1;const c=counts.get(observation)!;
      const k=c.k-interior.length,w=Math.max(0,c.w-omittedWeight),n=c.n-interior.length,key=k+'|'+w+'|'+n;
      let naive=exposureCache.get(key);if(naive===undefined){naive=calibration.mutationExposureFromCounts(observation,k,w,n).naive;exposureCache.set(key,naive);}
      priority+=naive;
    }
  }
  return { support, coverage, priority };
}

export function proposeNovelCandidates(
  observations: readonly Observation[],
  options: PersonalizedGermlineOptions,
  calibration:AlignedRateCalibration,
): { proposals: CandidateProposal[]; truncated: boolean } {
  const proposals = new Map<string, CandidateProposal>();
  const byParent = new Map<number, Observation[]>();
  for (const observation of observations) {
    const values = byParent.get(observation.parent.index);
    if (values) values.push(observation);
    else byParent.set(observation.parent.index, [observation]);
  }
  for (const items of byParent.values()) {
    const parent = items[0].parent;
    const counts=new Map(items.map(o=>{
      let k=0,w=0,n=0;for(let i=0;i<o.positions.length;i++){const p=o.positions[i];if(p>parent.sequence.length-TERMINAL_WINDOW)continue;k+=Number(o.query[i]!==parent.sequence[p-1]);w+=hs5fRate(parent.sequence,p-1);n++;}
      return [o,{k,w,n}] as const;
    }));
    const coverage = new Uint32Array(parent.sequence.length + 1);
    const support = new Map<string, number>();
    for (const observation of items) {
      for (let offset = 0; offset < observation.positions.length; offset += 1) {
        const position = observation.positions[offset];
        coverage[position] += 1;
        const alternate = observation.query[offset];
        const reference = parent.sequence[position - 1];
        if (alternate === reference) continue;
        const key = `${position}:${alternate}`;
        support.set(key, (support.get(key) ?? 0) + 1);
      }
    }
    const passing = new Map<string, PersonalizedSubstitution>();
    for (const [key, count] of support) {
      const [positionText, alternate = ""] = key.split(":");
      const position = Number(positionText);
      const covered = coverage[position];
      if (count < options.minimumNovelSupport || count / Math.max(1, covered) < options.minimumNovelFraction) continue;
      passing.set(key, { position, reference: parent.sequence[position - 1], alternate });
    }
    const patterns = new Map<string, { substitutions: PersonalizedSubstitution[]; observations: number }>();
    for (const observation of items) {
      const substitutions: PersonalizedSubstitution[] = [];
      for (let offset = 0; offset < observation.positions.length; offset += 1) {
        const position = observation.positions[offset];
        const alternate = observation.query[offset];
        const event = passing.get(`${position}:${alternate}`);
        if (event) substitutions.push(event);
      }
      substitutions.sort((left, right) => left.position - right.position || left.alternate.localeCompare(right.alternate));
      // Always propose the interior core independently of stochastic V/N tails.
      // Terminal extensions must subsequently beat that core on their own evidence.
      const core=substitutions.filter(c=>c.position<=parent.sequence.length-TERMINAL_WINDOW);
      const observable=substitutions.filter(c=>c.position<=parent.sequence.length-INITIAL_JUNCTION_BASES);
      // Unidentifiable final-two-base extensions never consume the candidate budget.
      const variants=[...new Map([core,observable].map(cs=>[cs.map(c=>c.position+':'+c.alternate).join('|'),cs])).values()];
      for(const changes of variants){
        if (!changes.length || changes.length > options.maximumNovelSnps) continue;
        const signature = changes.map((item) => `${item.position}:${item.alternate}`).join("|");
        const previous = patterns.get(signature);
        if (previous) previous.observations += 1;
        else patterns.set(signature, { substitutions:changes, observations: 1 });
      }
    }
    for (const pattern of patterns.values()) {
      // A seed need not recur as a pristine exact pattern: further SHM may add
      // other changes to carriers. Count independent linked support below.
      const evidence = proposalEvidence(items, pattern.substitutions, counts, calibration);
      if (evidence.support < options.minimumNovelSupport || evidence.support / Math.max(1, evidence.coverage) < options.minimumNovelFraction) continue;
      const sequence = [...parent.sequence];
      pattern.substitutions.forEach((item) => { sequence[item.position - 1] = item.alternate; });
      const joined = sequence.join("");
      const proposal: CandidateProposal = {
        signature: `${parent.names[0]}\u0000${pattern.substitutions.map((item) => `${item.reference}${item.position}${item.alternate}`).join("+")}`,
        candidate: {
          id: candidateId(parent, pattern.substitutions),
          names: [],
          sequence: joined,
          known: false,
          parent,
          substitutions: pattern.substitutions,
          directSupport: evidence.support,
          directCoverage: evidence.coverage,
          proposalPriority:evidence.priority,
        },
      };
      const existing = proposals.get(joined);
      if (!existing || proposal.candidate.directSupport > existing.candidate.directSupport) proposals.set(joined, proposal);
    }
  }
  const ordered = [...proposals.values()].sort((left, right) => (right.candidate.proposalPriority??0)-(left.candidate.proposalPriority??0)
    || right.candidate.directSupport - left.candidate.directSupport
    || right.candidate.directCoverage - left.candidate.directCoverage
    || left.signature.localeCompare(right.signature));
  const maximum = Math.max(0, Math.floor(options.maximumNovelCandidatesPerGene));
  return { proposals: ordered.slice(0, maximum), truncated: ordered.length > maximum };
}

function exposure(observation: Observation, commonSequence: string, sequencingErrorRate: number): number {
  let mismatches = 0;
  const weights: number[] = [];
  for (let offset = 0; offset < observation.positions.length; offset += 1) {
    const index = observation.positions[offset] - 1;
    if (index>=commonSequence.length-TERMINAL_WINDOW || !BASES.test(commonSequence[index])) continue;
    weights.push(hs5fRate(commonSequence, index));
    if (observation.query[offset] !== commonSequence[index]) mismatches += 1;
  }
  if (!weights.length) return 0.02;
  const corrected = Math.max(0, mismatches - weights.length * sequencingErrorRate);
  const target = Math.min(weights.length * 0.4, corrected + 0.5);
  const meanWeight = weights.reduce((sum, value) => sum + value, 0) / weights.length;
  let tau = Math.max(0, -Math.log(Math.max(1e-9, 1 - target / weights.length)) / Math.max(1e-9, meanWeight));
  for (let iteration = 0; iteration < 5; iteration += 1) {
    let expected = 0;
    let derivative = 0;
    for (const weight of weights) {
      const retained = Math.exp(-tau * weight);
      expected += 1 - retained;
      derivative += weight * retained;
    }
    if (!(derivative > 1e-12)) break;
    tau = Math.max(0, Math.min(2, tau - (expected - target) / derivative));
  }
  return tau;
}

function normalizedEmissions(observations: readonly Observation[], candidates: readonly Candidate[], sequencingErrorRate: number, calibration: AlignedRateCalibration, knownRadius:number, novelRadius:number): Float64Array[] {
  // Every competing hypothesis uses the same calibrated categorical model.
  const rates = candidates.map(candidate => {
    const proxy={...observations[0],gene:candidate.parent.gene,parent:candidate.parent};
    return Array.from(candidate.sequence,(reference,index)=>Float64Array.from("ACGT",alternate=>{
      if(alternate===reference)return 0;
      const prior=calibration.prior(proxy,{position:index+1,reference,alternate});
      return hs5fRate(candidate.sequence,index,alternate)*prior.shape/prior.rate;
    }));
  });
  const context = rates.map(rows=>Float64Array.from(rows,row=>row.reduce((a,b)=>a+b,0)));
  const compatibility = new Map([...new Set(observations.map(o=>o.parent))].map(p=>[p.index,candidates.map(c=>substitutionCompatible(p.sequence,c.sequence,c.known?knownRadius:knownRadius+novelRadius))]));
  const emissions: Float64Array[] = [];
  const error = Math.max(1e-8, Math.min(0.1, sequencingErrorRate));
  for (const observation of observations) {
    const tau = exposure(observation, observation.parent.sequence, error);
    const logLikelihood = new Float64Array(candidates.length);
    let maximum = Number.NEGATIVE_INFINITY;
    for (let candidateIndex = 0; candidateIndex < candidates.length; candidateIndex += 1) {
      const candidate = candidates[candidateIndex];
      if (!compatibility.get(observation.parent.index)![candidateIndex]) { logLikelihood[candidateIndex] = -Infinity; continue; }
      let value = 0;
      for (let offset = 0; offset < observation.positions.length; offset += 1) {
        const index = observation.positions[offset] - 1;
        if (index >= candidate.sequence.length-TERMINAL_WINDOW) continue;
        const germlineBase = candidate.sequence[index];
        if (!BASES.test(germlineBase)) continue;
        const mutation = 1 - Math.exp(-tau * context[candidateIndex][index]);
        const probability = error + (1 - error) * mutation;
        value += observation.query[offset] === germlineBase ? Math.log1p(-error)-tau*context[candidateIndex][index] : Math.log(Math.max(1e-12, error / 3 + (probability-error) * rates[candidateIndex][index]["ACGT".indexOf(observation.query[offset])] / Math.max(1e-12,context[candidateIndex][index])));
      }
      if(observation.terminal)value+=boundaryLogLikelihood(candidate.sequence,observation.terminal.query,observation.terminal.start,(position,base)=>{
        const index=position-1,mutation=1-Math.exp(-tau*context[candidateIndex][index]);
        return base===candidate.sequence[index]?Math.log1p(-error)-tau*context[candidateIndex][index]:Math.log(Math.max(1e-12,error/3+(1-error)*mutation*rates[candidateIndex][index]["ACGT".indexOf(base)]/Math.max(1e-12,context[candidateIndex][index])));
      },calibration.trimming({...observation,gene:candidate.parent.gene}),calibration.junction({...observation,gene:candidate.parent.gene}));
      logLikelihood[candidateIndex] = value;
      maximum = Math.max(maximum, value);
    }
    maximum=Math.max(...logLikelihood);
    emissions.push(Float64Array.from(logLikelihood, (value) => Math.max(1e-300, Math.exp(value - maximum))));
  }
  return emissions;
}

function fitMixture(
  emissions: readonly Float64Array[],
  active: readonly number[],
  options: PersonalizedGermlineOptions,
  initial?: Float64Array,
): MixtureFit {
  const theta = new Float64Array(active.length);
  if (initial?.length === active.length) theta.set(initial);
  else theta.fill(1 / Math.max(1, active.length));
  let totalTheta = theta.reduce((sum, value) => sum + Math.max(0, value), 0);
  if (!(totalTheta > 0)) totalTheta = 1;
  for (let index = 0; index < theta.length; index += 1) theta[index] = Math.max(0, theta[index]) / totalTheta;
  const counts = new Float64Array(active.length);
  const mixture = new Float64Array(emissions.length);
  let logLikelihood = Number.NEGATIVE_INFINITY;
  let converged = active.length <= 1;
  let iterations = 0;
  for (; iterations < Math.max(1, Math.floor(options.maximumIterations)); iterations += 1) {
    counts.fill(0);
    logLikelihood = 0;
    for (let row = 0; row < emissions.length; row += 1) {
      let normalizer = 0;
      for (let index = 0; index < active.length; index += 1) normalizer += theta[index] * emissions[row][active[index]];
      normalizer = Math.max(1e-300, normalizer);
      mixture[row] = normalizer;
      logLikelihood += Math.log(normalizer);
      for (let index = 0; index < active.length; index += 1) counts[index] += theta[index] * emissions[row][active[index]] / normalizer;
    }
    let maximumChange = 0;
    for (let index = 0; index < theta.length; index += 1) {
      const next = counts[index] / Math.max(1, emissions.length);
      maximumChange = Math.max(maximumChange, Math.abs(next - theta[index]));
      theta[index] = next;
    }
    if (maximumChange <= options.convergenceTolerance) {
      converged = true;
      iterations += 1;
      break;
    }
  }
  // The final E-step makes the returned mixture/counts exactly match theta.
  counts.fill(0);
  logLikelihood = 0;
  for (let row = 0; row < emissions.length; row += 1) {
    let normalizer = 0;
    for (let index = 0; index < active.length; index += 1) normalizer += theta[index] * emissions[row][active[index]];
    normalizer = Math.max(1e-300, normalizer);
    mixture[row] = normalizer;
    logLikelihood += Math.log(normalizer);
    for (let index = 0; index < active.length; index += 1) counts[index] += theta[index] * emissions[row][active[index]] / normalizer;
  }
  return { theta, mixture, logLikelihood, counts, iterations, converged };
}

function additionGain(baseline: MixtureFit, candidate: number, emissions: readonly Float64Array[]): { raw: number; frequency: number } {
  const derivative = (frequency: number) => {
    let result = 0;
    for (let row = 0; row < emissions.length; row += 1) {
      const ratio = emissions[row][candidate] / Math.max(1e-300, baseline.mixture[row]);
      result += (ratio - 1) / Math.max(1e-300, (1 - frequency) + frequency * ratio);
    }
    return result;
  };
  if (!(derivative(0) > 0)) return { raw: 0, frequency: 0 };
  let low = 0;
  let high = 0.95;
  if (derivative(high) > 0) low = high;
  else {
    for (let iteration = 0; iteration < 36; iteration += 1) {
      const middle = (low + high) / 2;
      if (derivative(middle) > 0) low = middle;
      else high = middle;
    }
  }
  const frequency = low;
  let raw = 0;
  for (let row = 0; row < emissions.length; row += 1) {
    const ratio = emissions[row][candidate] / Math.max(1e-300, baseline.mixture[row]);
    raw += Math.log(Math.max(1e-300, (1 - frequency) + frequency * ratio));
  }
  return { raw, frequency };
}

function modelPenalty(active: readonly number[], candidates: readonly Candidate[], lineages: number): number {
  const logN = Math.log(Math.max(2, lineages));
  const frequencyParameters = Math.max(0, active.length - 1);
  const searchCost = active.reduce((sum,index)=>sum+(candidates[index].known?0:(candidates[index].searchCost??discoveryCost(candidates[index].sequence.length,candidates[index].substitutions.length,candidates.length))),0);
  return 0.5 * frequencyParameters * logN + searchCost;
}

function fitCandidateSet(
  observations: readonly Observation[],
  candidates: readonly Candidate[],
  emissions: readonly Float64Array[],
  options: PersonalizedGermlineOptions,
): { active: number[]; fit: MixtureFit; gains: Map<number, number> } {
  let first = 0;
  let firstScore = Number.NEGATIVE_INFINITY;
  for (let candidate = 0; candidate < candidates.length; candidate += 1) {
    let score = -modelPenalty([candidate], candidates, observations.length);
    for (const row of emissions) score += Math.log(Math.max(1e-300, row[candidate]));
    if (score > firstScore || (score === firstScore && candidates[candidate].id.localeCompare(candidates[first].id) < 0)) {
      first = candidate;
      firstScore = score;
    }
  }
  let active = [first];
  let fit = fitMixture(emissions, active, options);
  const gains = new Map<number, number>([[first, Number.POSITIVE_INFINITY]]);
  const maximumActive = Math.max(1, Math.floor(options.maximumActiveAllelesPerGene) * new Set(observations.map(o=>o.gene)).size);
  while (active.length < Math.min(maximumActive, candidates.length)) {
    let best = -1;
    let bestAdjusted = Number.NEGATIVE_INFINITY;
    let bestFrequency = 0;
    for (let candidate = 0; candidate < candidates.length; candidate += 1) {
      if (active.includes(candidate)) continue;
      const gain = additionGain(fit, candidate, emissions);
      const adjusted = gain.raw - (modelPenalty([...active,candidate],candidates,observations.length)-modelPenalty(active,candidates,observations.length));
      if (adjusted > bestAdjusted || (adjusted === bestAdjusted && candidates[candidate].id.localeCompare(candidates[best]?.id ?? "") < 0)) {
        best = candidate;
        bestAdjusted = adjusted;
        bestFrequency = gain.frequency;
      }
    }
    if (best < 0 || bestAdjusted < options.minimumLogEvidenceGain || !(bestFrequency > 0)) break;
    const initial = new Float64Array(active.length + 1);
    for (let index = 0; index < active.length; index += 1) initial[index] = fit.theta[index] * (1 - bestFrequency);
    initial[active.length] = bestFrequency;
    active = [...active, best];
    fit = fitMixture(emissions, active, options, initial);
    gains.set(best, bestAdjusted);
  }
  // Correlated additions can make an earlier component redundant. Remove it
  // only when the complete BIC objective improves after an exact sparse refit.
  let currentObjective = fit.logLikelihood - modelPenalty(active, candidates, observations.length);
  let changed = true;
  while (changed && active.length > 1) {
    changed = false;
    let bestActive = active;
    let bestFit = fit;
    let bestObjective = currentObjective;
    for (let remove = 0; remove < active.length; remove += 1) {
      const proposed = active.filter((_, index) => index !== remove);
      const proposedFit = fitMixture(emissions, proposed, options);
      const objective = proposedFit.logLikelihood - modelPenalty(proposed, candidates, observations.length);
      if (objective > bestObjective + 1e-8) {
        bestActive = proposed;
        bestFit = proposedFit;
        bestObjective = objective;
      }
    }
    if (bestActive !== active) {
      active = bestActive;
      fit = bestFit;
      currentObjective = bestObjective;
      changed = true;
    }
  }
  return { active, fit, gains };
}

function inferGene(
  observations: readonly Observation[],
  nodes: readonly ReferenceNode[],
  options: PersonalizedGermlineOptions,
  calibration: AlignedRateCalibration,
  parentCount: number,
): { result: PersonalizedGermlineGeneResult; truncated: boolean } {
  const parents = [...new Set(observations.map((item) => item.parent))];
  const candidatesBySequence = new Map<string, Candidate>();
  for (const node of nodes) {
    if (node.locus !== observations[0].locus || node.sequence.length !== observations[0].parent.sequence.length) continue;
    if (!parents.some((parent) => substitutionCompatible(parent.sequence, node.sequence, options.maximumKnownAlleleSnps))) continue;
    const alias=candidatesBySequence.get(node.sequence);
    if(alias){alias.names.push(...node.names.filter(name=>!alias.names.includes(name)));continue;}
    candidatesBySequence.set(node.sequence, {
      id: node.names.join(","),
      names: [...node.names],
      sequence: node.sequence,
      known: true,
      parent: node,
      substitutions: [],
      directSupport: 0,
      directCoverage: observations.length,
    });
  }
  for (const parent of parents) {
    if (!candidatesBySequence.has(parent.sequence)) candidatesBySequence.set(parent.sequence, {
      id: parent.names.join(","), names: [...parent.names], sequence: parent.sequence, known: true, parent,
      substitutions: [], directSupport: 0, directCoverage: observations.length,
    });
  }
  const novel = proposeNovelCandidates(observations, {...options,maximumNovelCandidatesPerGene: options.maximumNovelCandidatesPerGene * new Set(observations.map(o=>o.gene)).size},calibration);
  const diagnostics: NonNullable<PersonalizedGermlineGeneResult["proposalEvidence"]> = [];
  const accepted: Candidate[] = [];
  for (const proposal of [...novel.proposals].sort((a,b)=>a.candidate.substitutions.length-b.candidate.substitutions.length || b.candidate.directSupport-a.candidate.directSupport)) {
    const candidate=proposal.candidate;
    if(candidatesBySequence.has(candidate.sequence))continue;
    if(candidate.sequence.length!==observations[0].parent.sequence.length)continue;
    const coreChanges=candidate.substitutions.filter(c=>c.position<=candidate.sequence.length-TERMINAL_WINDOW);
    const hasBoundary=coreChanges.length<candidate.substitutions.length;
    const coreSequence=[...candidate.parent.sequence];coreChanges.forEach(c=>coreSequence[c.position-1]=c.alternate);
    if(hasBoundary&&coreChanges.length&&!candidatesBySequence.has(coreSequence.join(''))){
      diagnostics.push({id:candidate.id,gain:-Infinity,support:candidate.directSupport,coverage:candidate.directCoverage,baseline:'Interior haplotype not supported independently',selected:false});continue;
    }
    const subsets=accepted.filter(c=>c.parent===candidate.parent && c.substitutions.length<candidate.substitutions.length && c.substitutions.every(x=>candidate.substitutions.some(y=>x.position===y.position&&x.alternate===y.alternate))).sort((a,b)=>b.substitutions.length-a.substitutions.length||b.directSupport-a.directSupport);
    const baseline=subsets[0];
    const changes=candidate.substitutions.filter(c=>!baseline?.substitutions.some(b=>b.position===c.position&&b.alternate===c.alternate));
    const items=observations.filter(o=>o.parent===candidate.parent && (!baseline||baseline.substitutions.every(c=>baseAt(o,c.position)===c.alternate)));
    const test=testHaplotype(items,changes,candidate.substitutions,calibration,options.sequencingErrorRate,parentCount);
    diagnostics.push({id:candidate.id,gain:test.gain,support:candidate.directSupport,coverage:candidate.directCoverage,baseline:baseline?.id??candidate.parent.names[0],selected:false,reason:test.reason});
    if(!(test.gain>options.minimumLogEvidenceGain))continue;
    candidate.searchCost=discoveryCost(candidate.sequence.length,candidate.substitutions.length,parentCount);
    accepted.push(candidate);candidatesBySequence.set(candidate.sequence,candidate);
  }
  // An accepted allele must also compete with projections built on a different
  // reference backbone. Original-parent-only subset tests cannot establish that
  // the inherited distinguishing bases were actually observed together.
  for(const candidate of [...accepted].sort((a,b)=>b.directSupport-a.directSupport)){
    const isSuperset=(other:Candidate)=>other.parent===candidate.parent&&other.substitutions.length>candidate.substitutions.length&&candidate.substitutions.every(c=>other.substitutions.some(s=>s.position===c.position&&s.alternate===c.alternate));
    const alternatives=accepted.filter(other=>other!==candidate && candidatesBySequence.has(other.sequence) && ((other.parent!==candidate.parent&&other.directSupport>candidate.directSupport)||isSuperset(other)) && substitutionCompatible(other.sequence,candidate.sequence,options.maximumKnownAlleleSnps))
      .sort((a,b)=>hammingDistance(a.sequence,candidate.sequence)-hammingDistance(b.sequence,candidate.sequence)||b.directSupport-a.directSupport);
    const baseline=alternatives[0];if(!baseline)continue;
    const proxyParent={...baseline.parent,sequence:baseline.sequence};
    const changes=candidateSubstitutions(proxyParent,candidate.sequence);
    const items=observations.filter(o=>(o.parent===candidate.parent||o.parent===baseline.parent)&&(!isSuperset(baseline)||candidate.substitutions.every(c=>baseAt(o,c.position)===c.alternate))).map(o=>({...o,gene:proxyParent.gene,parent:proxyParent}));
    const test=testHaplotype(items,changes,changes,calibration,options.sequencingErrorRate,parentCount);
    diagnostics.push({id:candidate.id+" [conditional backbone]",gain:test.gain,support:candidate.directSupport,coverage:items.length,baseline:baseline.id,selected:false,reason:test.reason});
    if(!(test.gain>options.minimumLogEvidenceGain))candidatesBySequence.delete(candidate.sequence);
  }
  const candidates = [...candidatesBySequence.values()].sort((left, right) => Number(right.known) - Number(left.known) || left.id.localeCompare(right.id, undefined, { numeric: true }));
  const emissions = normalizedEmissions(observations, candidates, options.sequencingErrorRate, calibration, options.maximumKnownAlleleSnps, options.maximumNovelSnps);
  for (let candidate = 0; candidate < candidates.length; candidate += 1) {
    let localBest = 0;
    for (const row of emissions) {
      const maximum = Math.max(...row);
      if (Math.abs(row[candidate] - maximum) <= 1e-12) localBest += 1 / Math.max(1, row.reduce((count, value) => count + Number(Math.abs(value - maximum) <= 1e-12), 0));
    }
    if(candidates[candidate].known)candidates[candidate].directSupport = localBest;
  }
  const selected = fitCandidateSet(observations, candidates, emissions, options);
  const activeAlleles = selected.active.map((candidateIndex, activeIndex): PersonalizedGermlineAllele => {
    const candidate = candidates[candidateIndex];
    let localBestLineages = 0;
    for (const row of emissions) {
      const maximum = Math.max(...selected.active.map((index) => row[index]));
      if (Math.abs(row[candidateIndex] - maximum) <= 1e-12) localBestLineages += 1 / Math.max(1, selected.active.reduce((count, index) => count + Number(Math.abs(row[index] - maximum) <= 1e-12), 0));
    }
    return {
      id: candidate.id,
      names: [...candidate.names],
      sequence: candidate.sequence,
      known: candidate.known,
      parentAllele: candidate.parent.names[0],
      substitutions: candidate.known ? candidateSubstitutions(candidate.parent, candidate.sequence) : [...candidate.substitutions],
      frequency: selected.fit.theta[activeIndex],
      expectedLineages: selected.fit.counts[activeIndex],
      localBestLineages,
      directSupport: candidate.directSupport,
      directCoverage: candidate.directCoverage,
      discoveryGain: candidate.known?null:Math.min(...diagnostics.filter(d=>d.id===candidate.id||d.id===candidate.id+" [conditional backbone]").map(d=>d.gain)),
      selectionGain: Number.isFinite(selected.gains.get(candidateIndex) ?? NaN) ? selected.gains.get(candidateIndex)! : null,
    };
  }).sort((left, right) => right.frequency - left.frequency || left.id.localeCompare(right.id, undefined, { numeric: true }));
  return {
    truncated: novel.truncated,
    result: {
      gene: [...new Set(observations.map(o=>o.gene))].sort().join(","),
      representativeLineages: observations.length,
      testedCandidates: candidates.length,
      proposedNovelCandidates: candidates.reduce((sum, candidate) => sum + Number(!candidate.known), 0),
      activeAlleles,
      relativeLogLikelihood: selected.fit.logLikelihood,
      iterations: selected.fit.iterations,
      converged: selected.fit.converged,
      proposalEvidence: diagnostics.map(d=>({...d,selected:activeAlleles.some(a=>a.id===d.id)})),
    },
  };
}

export class PersonalizedGermlineAccumulator {
  private readonly options: PersonalizedGermlineOptions;
  private readonly nodes: ReferenceNode[];
  private readonly byName: Map<string, ReferenceNode>;
  private readonly selected = new Map<string, Observation>();
  private readonly seenLineages = new Set<string>();
  private inputRecords = 0;
  private eligibleRecords = 0;

  constructor(referenceFasta: string, options: Partial<PersonalizedGermlineOptions> = {}) {
    this.options = { ...DEFAULT_PERSONALIZED_GERMLINE_OPTIONS, ...options };
    const references = referenceNodes(referenceFasta);
    this.nodes = references.nodes;
    this.byName = references.byName;
    if (!this.nodes.length) throw new Error("The V reference database contains no usable named nucleotide sequences.");
  }

  add(row: Record<string, string>, ordinal: number, lineageId = 0) {
    this.inputRecords += 1;
    if (!(lineageId > 0)) return;
    const subjectId = (row.subject_id || "unassigned-subject").trim() || "unassigned-subject";
    const unit = `${subjectId}\u0000${Math.floor(lineageId)}`;
    this.seenLineages.add(unit);
    const observation = parseObservation(row, ordinal, Math.floor(lineageId), this.byName, this.options);
    if (!observation) return;
    this.eligibleRecords += 1;
    const existing = this.selected.get(unit);
    if (!existing
      || observation.shmRate < existing.shmRate
      || (observation.shmRate === existing.shmRate && observation.alignedBases > existing.alignedBases)
      || (observation.shmRate === existing.shmRate && observation.alignedBases === existing.alignedBases && observation.ordinal < existing.ordinal)) this.selected.set(unit, observation);
  }

  /** Read-only research input; does not perform inference or alter selection. */
  researchSnapshot() { return { observations: [...this.selected.values()], references: [...this.nodes] }; }

  selectedRepresentativeOrdinals(): number[] {
    return [...this.selected.values()].map((item) => item.ordinal).sort((left, right) => left - right);
  }

  finish(onProgress?: (processedGenes: number, totalGenes: number) => void): PersonalizedGermlineDashboard {
    const observations=[...this.selected.values()];
    const calibration=new AlignedRateCalibration(observations,this.options.sequencingErrorRate);
    const subjectPools=new Map<string,Observation[]>();
    for(const o of observations){const key=`${o.subjectId}\u0000${o.locus}`;const values=subjectPools.get(key)??[];values.push(o);subjectPools.set(key,values);}
    const parentCounts=new Map<string,number>();
    const byGene = new Map<string, Observation[]>();
    for(const [poolKey,pool] of subjectPools){
      const parents=[...new Set(pool.map(o=>o.parent))];parentCounts.set(poolKey,parents.length);
      const roots=new Map(parents.map(p=>[p.index,p.index]));
      const root=(id:number):number=>{let r=id;while(roots.get(r)!==r)r=roots.get(r)!;return r;};
      for(let i=0;i<parents.length;i++)for(let j=0;j<i;j++){
        const a=parents[i],b=parents[j];
        if(substitutionCompatible(a.sequence,b.sequence,this.options.maximumKnownAlleleSnps))roots.set(root(a.index),root(b.index));
      }
      for (const observation of pool) {
        const key = `${poolKey}\u0000${root(observation.parent.index)}\u0000${observation.parent.sequence.length}`;
        const values=byGene.get(key);if(values)values.push(observation);else byGene.set(key,[observation]);
      }
    }
    const entries = [...byGene.entries()].sort(([left], [right]) => left.localeCompare(right, undefined, { numeric: true }));
    const poolByKey = new Map<string, PersonalizedGermlinePool>();
    let testedCandidates = 0;
    let activeKnownAlleles = 0;
    let activeNovelAlleles = 0;
    let proposalTruncations = 0;
    entries.forEach(([key, observations], index) => {
      const [subjectId, locus] = key.split("\u0000");
      onProgress?.(index, entries.length);
      const inferred = inferGene(observations, this.nodes, this.options, calibration, parentCounts.get(`${subjectId}\u0000${locus}`)!);
      if (inferred.truncated) proposalTruncations += 1;
      testedCandidates += inferred.result.testedCandidates;
      activeKnownAlleles += inferred.result.activeAlleles.reduce((sum, allele) => sum + Number(allele.known), 0);
      activeNovelAlleles += inferred.result.activeAlleles.reduce((sum, allele) => sum + Number(!allele.known), 0);
      const poolKey = `${subjectId}\u0000${locus}`;
      let pool = poolByKey.get(poolKey);
      if (!pool) {
        pool = { id: "", subjectId, locus, representativeLineages: 0, genes: [] };
        poolByKey.set(poolKey, pool);
      }
      pool.representativeLineages += observations.length;
      pool.genes.push(inferred.result);
      onProgress?.(index + 1, entries.length);
    });
    const pools = [...poolByKey.values()].sort((left, right) => left.subjectId.localeCompare(right.subjectId, undefined, { numeric: true }) || left.locus.localeCompare(right.locus));
    pools.forEach((pool, index) => {
      pool.id = `personalized_pool_${index + 1}`;
      pool.genes.sort((left, right) => left.gene.localeCompare(right.gene, undefined, { numeric: true }));
    });
    const warnings: string[] = [];
    const skippedLineages = Math.max(0, this.seenLineages.size - this.selected.size);
    if (skippedLineages) warnings.push(`${skippedLineages.toLocaleString()} assigned lineage${skippedLineages === 1 ? " had" : "s had"} no member with a recognized V call and at least ${this.options.minimumAlignedBases} aligned V nucleotides.`);
    if ([...this.selected.values()].some((item) => item.subjectId === "unassigned-subject")) warnings.push("Rows without subject_id were pooled together. Supply subject identifiers before interpreting the result as a per-person genotype.");
    if (proposalTruncations) warnings.push(`${proposalTruncations.toLocaleString()} gene fit${proposalTruncations === 1 ? " reached" : "s reached"} the novel-candidate cap; hypotheses with the largest expected unmutated-lineage support were retained. Increase the candidate cap before interpreting non-recovery in those genes.`);
    if (entries.some(([, observations]) => observations.length < 4)) warnings.push("Some expressed V genes have fewer than four usable lineage representatives; their active sets are weakly identified and should not be treated as genomic absence calls.");
    warnings.push("Novel changes confined to the final two V bases are unresolved under the gene-specific junction-start model. Interior candidate sequences retain their parent terminal bases; those bases are not independently inferred.");
    warnings.push("Discovery scores are model-comparison statistics, not calibrated false-discovery probabilities. Stable full-reference candidates still need independent validation.");
    warnings.push("Only expressed, same-length V-allele hypotheses are identifiable here. Untested genes are retained in downloaded references, and a final full reassignment is still required.");
    return {
      version: 1,
      mode: "lowest-current-v-shm-per-lineage",
      contextModel: "hs5f-hurdle-profiled-junction",
      options: { ...this.options },
      inputRecords: this.inputRecords,
      eligibleRecords: this.eligibleRecords,
      assignedLineages: this.seenLineages.size,
      representativeLineages: this.selected.size,
      skippedLineages,
      pools,
      testedCandidates,
      activeKnownAlleles,
      activeNovelAlleles,
      proposalTruncations,
      warnings,
    };
  }
}

export function personalizedGermlineEvidenceRows(dashboard: PersonalizedGermlineDashboard): Array<Record<string, string | number>> {
  return dashboard.pools.flatMap((pool) => pool.genes.flatMap((gene) => gene.activeAlleles.map((allele) => ({
    pool_id: pool.id,
    subject_id: pool.subjectId,
    locus: pool.locus,
    gene: gene.gene,
    allele: allele.names.join(",") || allele.id,
    known_reference: allele.known ? 1 : 0,
    parent_allele: allele.parentAllele,
    substitutions: allele.substitutions.map((item) => `${item.reference}${item.position}${item.alternate}`).join(";"),
    representative_lineages: gene.representativeLineages,
    inferred_frequency: allele.frequency,
    expected_lineages: allele.expectedLineages,
    local_best_lineages: allele.localBestLineages,
    direct_support: allele.directSupport,
    direct_coverage: allele.directCoverage,
    discovery_log_score: allele.discoveryGain ?? "known",
    bic_adjusted_selection_gain: allele.selectionGain ?? "initial",
  }))));
}

/**
 * Conservative replacement reference for one subject/locus pool: inferred
 * alleles replace only genes that were actually tested. Every untested gene
 * and every other locus remains unchanged.
 */
export function personalizedGermlineFasta(referenceFasta: string, dashboard: PersonalizedGermlineDashboard, poolId: string): string {
  const pool = dashboard.pools.find((item) => item.id === poolId);
  if (!pool) return "";
  const activeNames = new Set(pool.genes.flatMap((gene) => gene.activeAlleles.filter((allele) => allele.known).flatMap((allele) => allele.names)));
  const testedGenes = new Set(pool.genes.flatMap((gene) => gene.gene.split(",")));
  const records = parseReferenceFasta(referenceFasta);
  const retained = records.filter((record) => {
    if (referenceLocus(record.name) !== pool.locus) return true;
    const gene = personalizedAlleleGene(record.name);
    return !testedGenes.has(gene) || activeNames.has(record.name);
  });
  const byName = new Map(records.map((record) => [record.name, record] as const));
  const additions: ReferenceFastaRecord[] = [];
  for (const gene of pool.genes) {
    for (const allele of gene.activeAlleles) {
      if (allele.known) continue;
      const parent = byName.get(allele.parentAllele);
      const metadata = parent?.header.match(/(?:^|\s)(SWIGMETA=[\d,\-]+)/i)?.[1] ?? "";
      const substitutions = allele.substitutions.map((item) => `${item.reference}${item.position}${item.alternate}`).join("+");
      additions.push({
        name: allele.id,
        header: `${allele.id}${metadata ? ` ${metadata}` : ""} SWIG_PERSONALIZED=novel PARENT=${allele.parentAllele} SUBSTITUTIONS=${substitutions} POOL=${pool.id}`,
        sequence: allele.sequence,
      });
    }
  }
  return serializeReferenceFasta([...retained, ...additions]);
}
