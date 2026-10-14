export type ShmMetricKey = "vNtMutations" | "vNtRate" | "vAaReplacements" | "vAaRate" | "synonymous" | "cdrNtRate" | "frameworkNtRate";

export interface ShmRecordMetric {
  ordinal: number;
  lineageId: number;
  sequenceId: string;
  vCall: string;
  jCall: string;
  sampleId: string;
  subjectId: string;
  timepoint: string;
  compartment: string;
  cdr3Nt: string;
  cdr3Aa: string;
  duplicateCount: number;
  comparedNt: number;
  vNtMutations: number;
  vNtRate: number;
  comparedCodons: number;
  vAaReplacements: number;
  vAaRate: number;
  synonymous: number;
  cdrNtCompared: number;
  cdrNtMutations: number;
  cdrNtRate: number;
  frameworkNtCompared: number;
  frameworkNtMutations: number;
  frameworkNtRate: number;
  stratum: string;
}

export interface ShmDistribution {
  label: string;
  values: number[];
  weights: number[];
  records: number;
  abundance: number;
  median: number;
  mean: number;
  maximum: number;
  p95: number;
}

export interface ShmLineageSampleDistribution extends ShmDistribution {
  lineageId: number;
  sampleId: string;
  subjectId: string;
  timepoint: string;
  compartment: string;
}

export interface ShmDashboard {
  analyzedRecords: number;
  analyzedAbundance: number;
  skippedRecords: number;
  sampledRecords: number;
  metric: ShmMetricKey;
  records: ShmRecordMetric[];
  /** Exact lowest V-nucleotide-SHM member seen for each assigned lineage. */
  lowestByLineage?: Array<{ lineageId: number; ordinal: number; vNtRate: number; vNtMutations: number; cdr3Nt: string; cdr3Aa: string }>;
  lineages: ShmDistribution[];
  vGenes: ShmDistribution[];
  strata: ShmDistribution[];
  lineageSamples: ShmLineageSampleDistribution[];
  histogram: Array<{ label: string; count: number; abundance: number }>;
}

export interface ShmAccumulatorOptions {
  metric?: ShmMetricKey;
  maxSamplesPerLineage?: number;
  maxGlobalSamples?: number;
}

const CODON: Record<string, string> = {
  TTT:"F",TTC:"F",TTA:"L",TTG:"L",TCT:"S",TCC:"S",TCA:"S",TCG:"S",TAT:"Y",TAC:"Y",TAA:"*",TAG:"*",TGT:"C",TGC:"C",TGA:"*",TGG:"W",
  CTT:"L",CTC:"L",CTA:"L",CTG:"L",CCT:"P",CCC:"P",CCA:"P",CCG:"P",CAT:"H",CAC:"H",CAA:"Q",CAG:"Q",CGT:"R",CGC:"R",CGA:"R",CGG:"R",
  ATT:"I",ATC:"I",ATA:"I",ATG:"M",ACT:"T",ACC:"T",ACA:"T",ACG:"T",AAT:"N",AAC:"N",AAA:"K",AAG:"K",AGT:"S",AGC:"S",AGA:"R",AGG:"R",
  GTT:"V",GTC:"V",GTA:"V",GTG:"V",GCT:"A",GCC:"A",GCA:"A",GCG:"A",GAT:"D",GAC:"D",GAA:"E",GAG:"E",GGT:"G",GGC:"G",GGA:"G",GGG:"G",
};

function topCall(value: string): string { return value.split(",")[0]?.trim() ?? ""; }
function numeric(value: string): number { const result = Number(value); return Number.isFinite(result) ? result : 0; }
function positiveCount(value: string): number { const result = Math.round(numeric(value)); return result > 0 ? result : 1; }
function valid(base: string): boolean { return /^[ACGT]$/.test(base); }
function rangeValue(row: Record<string,string>, name: string): number { return Math.max(0, Math.floor(numeric(row[name]))); }

function inRegion(queryPosition: number, row: Record<string,string>, prefix: "cdr1" | "cdr2" | "fwr1" | "fwr2" | "fwr3"): boolean {
  const start = rangeValue(row, `${prefix}_start`);
  const end = rangeValue(row, `${prefix}_end`);
  return start > 0 && end >= start && queryPosition >= start && queryPosition <= end;
}

export function computeShmMetric(row: Record<string,string>, ordinal = 0, lineageId = 0, stratum = "all"): ShmRecordMetric | null {
  const query = (row.v_sequence_alignment || row.sequence_alignment || "").toUpperCase().replace(/\./g, "-");
  const germline = (row.v_germline_alignment || row.germline_alignment || "").toUpperCase().replace(/\./g, "-");
  if (!query || !germline) return null;
  const length = Math.min(query.length, germline.length);
  const queryStart = Math.max(1, rangeValue(row, "v_sequence_start") || 1);
  let queryPosition = queryStart - 1;
  let comparedNt = 0, vNtMutations = 0, cdrNtCompared = 0, cdrNtMutations = 0, frameworkNtCompared = 0, frameworkNtMutations = 0;
  const codonQuery: string[] = [], codonGermline: string[] = [];
  let synonymous = 0, vAaReplacements = 0, comparedCodons = 0;

  const flushCodon = () => {
    if (codonQuery.length !== 3 || codonGermline.length !== 3) return;
    const q = codonQuery.join(""), g = codonGermline.join("");
    if (!/^[ACGT]{3}$/.test(q) || !/^[ACGT]{3}$/.test(g)) return;
    comparedCodons += 1;
    if (q !== g) {
      if (CODON[q] === CODON[g]) synonymous += 1;
      else vAaReplacements += 1;
    }
  };

  // Bases of the current codon before the V alignment's first base: codons
  // start at read positions ≡ sequence_frame (mod 3), one-based. Rows
  // without a frame (IgBLAST's) take the V germline's: IMGT V genes start
  // with a whole codon, so v_germline_start − 1 bases precede the alignment.
  const frame = rangeValue(row, "sequence_frame") || rangeValue(row, "v_frame");
  const germlineStart = rangeValue(row, "v_germline_start");
  let codingBases = frame ? (((queryStart - frame) % 3) + 3) % 3 : germlineStart ? (germlineStart - 1) % 3 : 0;
  for (let column = 0; column < length; column += 1) {
    const q = query[column], g = germline[column];
    if (q !== "-") queryPosition += 1;
    if (valid(q) && valid(g)) {
      comparedNt += 1;
      const mismatch = q !== g;
      if (mismatch) vNtMutations += 1;
      const cdr = inRegion(queryPosition, row, "cdr1") || inRegion(queryPosition, row, "cdr2");
      const framework = inRegion(queryPosition, row, "fwr1") || inRegion(queryPosition, row, "fwr2") || inRegion(queryPosition, row, "fwr3");
      if (cdr) { cdrNtCompared += 1; if (mismatch) cdrNtMutations += 1; }
      if (framework) { frameworkNtCompared += 1; if (mismatch) frameworkNtMutations += 1; }
    }
    if (q !== "-" && g !== "-") {
      if (codingBases < 0) codingBases = 0;
      codonQuery.push(q); codonGermline.push(g); codingBases += 1;
      if (codingBases % 3 === 0) { flushCodon(); codonQuery.length = 0; codonGermline.length = 0; }
    } else if (q !== "-" || g !== "-") {
      // An indel disrupts the aligned codon; resume at the next frame boundary.
      codonQuery.length = 0; codonGermline.length = 0; codingBases += q !== "-" ? 1 : 0;
    }
  }
  if (!comparedNt) return null;
  return {
    ordinal, lineageId, sequenceId: row.sequence_id || `record_${ordinal + 1}`, vCall: topCall(row.v_call), jCall: topCall(row.j_call),
    sampleId: row.sample_id || "Unassigned sample", subjectId: row.subject_id || "", timepoint: row.swig_timepoint || "", compartment: row.swig_compartment || "",
    cdr3Nt: row.cdr3 || row.junction || "", cdr3Aa: row.cdr3_aa || row.junction_aa || "",
    duplicateCount: positiveCount(row.duplicate_count || row.consensus_count || "1"), comparedNt, vNtMutations,
    vNtRate: vNtMutations / comparedNt, comparedCodons, vAaReplacements, vAaRate: comparedCodons ? vAaReplacements / comparedCodons : 0,
    synonymous, cdrNtCompared, cdrNtMutations, cdrNtRate: cdrNtCompared ? cdrNtMutations / cdrNtCompared : 0,
    frameworkNtCompared, frameworkNtMutations, frameworkNtRate: frameworkNtCompared ? frameworkNtMutations / frameworkNtCompared : 0,
    stratum,
  };
}

export function metricValue(record: ShmRecordMetric, metric: ShmMetricKey): number { return record[metric]; }

function weightedQuantile(values: Array<{value:number;weight:number}>, quantile: number): number {
  if (!values.length) return 0;
  const sorted = [...values].sort((a,b) => a.value-b.value);
  const threshold = sorted.reduce((sum,item)=>sum+item.weight,0) * Math.max(0, Math.min(1, quantile));
  let total=0;
  for (const item of sorted) { total += item.weight; if (total >= threshold) return item.value; }
  return sorted.at(-1)?.value ?? 0;
}

function distributions(records: ShmRecordMetric[], metric: ShmMetricKey, key: (record:ShmRecordMetric)=>string): ShmDistribution[] {
  const groups = new Map<string,ShmRecordMetric[]>();
  for (const record of records) { const label=key(record); if (!label) continue; const group=groups.get(label)??[]; group.push(record); groups.set(label,group); }
  return [...groups.entries()].map(([label,group]) => {
    const values=group.map((record)=>metricValue(record,metric)); const weights=group.map((record)=>record.duplicateCount);
    const abundance=weights.reduce((a,b)=>a+b,0); const weighted=values.map((value,index)=>({value,weight:weights[index]}));
    return {label,values,weights,records:group.length,abundance,median:weightedQuantile(weighted,.5),mean:abundance?weighted.reduce((sum,item)=>sum+item.value*item.weight,0)/abundance:0,maximum:Math.max(...values),p95:weightedQuantile(weighted,.95)};
  }).sort((a,b)=>b.abundance-a.abundance||a.label.localeCompare(b.label));
}

function lineageSampleDistributions(records: ShmRecordMetric[], metric: ShmMetricKey): ShmLineageSampleDistribution[] {
  const distributionsByPair = distributions(records, metric, (record) => record.lineageId > 0 ? `${record.lineageId}\u0000${record.sampleId}` : "");
  const metadata = new Map<string, ShmRecordMetric>();
  for (const record of records) if (record.lineageId > 0) metadata.set(`${record.lineageId}\u0000${record.sampleId}`, record);
  return distributionsByPair.map((distribution) => {
    const record = metadata.get(distribution.label)!;
    return { ...distribution, label: `Lineage ${record.lineageId} · ${record.sampleId}`, lineageId: record.lineageId, sampleId: record.sampleId, subjectId: record.subjectId, timepoint: record.timepoint, compartment: record.compartment };
  });
}

export class ShmAccumulator {
  private readonly options: Required<ShmAccumulatorOptions>;
  private readonly samples:ShmRecordMetric[]=[];
  private readonly sampledByLineageSample=new Map<string,number>();
  private readonly lowestByLineage=new Map<number,ShmRecordMetric>();
  private skipped = 0;
  private analyzed = 0;
  private abundance = 0;

  constructor(options: ShmAccumulatorOptions = {}) {
    this.options = {metric:options.metric ?? "vNtRate",maxSamplesPerLineage:Math.max(10,options.maxSamplesPerLineage??2000),maxGlobalSamples:Math.max(1000,options.maxGlobalSamples??100000)};
  }

  add(row: Record<string,string>, ordinal: number, lineageId = 0, stratum = "all") {
    const metric=computeShmMetric(row,ordinal,lineageId,stratum);
    if (!metric) { this.skipped += 1; return; }
    this.analyzed += 1; this.abundance += metric.duplicateCount;
    if(lineageId>0){const previous=this.lowestByLineage.get(lineageId);if(!previous||metric.vNtRate<previous.vNtRate||(metric.vNtRate===previous.vNtRate&&(metric.vNtMutations<previous.vNtMutations||(metric.vNtMutations===previous.vNtMutations&&metric.ordinal<previous.ordinal))))this.lowestByLineage.set(lineageId,metric);}
    // Cap each lineage × sample cell independently. A repertoire ordered by
    // sample must not fill a lineage reservoir at the first timepoint and
    // silently erase all later longitudinal observations.
    const pair=`${lineageId}\u0000${metric.sampleId}`;
    const pairSamples=this.sampledByLineageSample.get(pair)??0;
    if(pairSamples>=this.options.maxSamplesPerLineage)return;
    if(this.samples.length<this.options.maxGlobalSamples){this.samples.push(metric);this.sampledByLineageSample.set(pair,pairSamples+1);return;}
    // Deterministic global reservoir: plot/session memory remains bounded even when every row is a singleton lineage.
    const slot=((ordinal*2654435761)>>>0)%this.analyzed;
    if(slot<this.options.maxGlobalSamples){const removed=this.samples[slot];const removedPair=`${removed.lineageId}\u0000${removed.sampleId}`;const removedCount=this.sampledByLineageSample.get(removedPair)??1;if(removedCount<=1)this.sampledByLineageSample.delete(removedPair);else this.sampledByLineageSample.set(removedPair,removedCount-1);this.samples[slot]=metric;this.sampledByLineageSample.set(pair,(this.sampledByLineageSample.get(pair)??0)+1);}
  }

  finish(): ShmDashboard {
    const records=this.samples; const metric=this.options.metric;
    const bins=Array.from({length:20},()=>({count:0,abundance:0}));
    for (const record of records) { const value=metricValue(record,metric); const normalized=metric.toLowerCase().includes("rate")?Math.min(.999999,Math.max(0,value)):Math.min(.999999,Math.max(0,value/50)); const bin=Math.floor(normalized*bins.length); bins[bin].count+=1; bins[bin].abundance+=record.duplicateCount; }
    return {analyzedRecords:this.analyzed,analyzedAbundance:this.abundance,skippedRecords:this.skipped,sampledRecords:records.length,metric,records,
      lowestByLineage:[...this.lowestByLineage.values()].sort((left,right)=>left.lineageId-right.lineageId).map((record)=>({lineageId:record.lineageId,ordinal:record.ordinal,vNtRate:record.vNtRate,vNtMutations:record.vNtMutations,cdr3Nt:record.cdr3Nt,cdr3Aa:record.cdr3Aa})),
      lineages:distributions(records,metric,(record)=>record.lineageId?`Lineage ${record.lineageId}`:"Unassigned"),
      vGenes:distributions(records,metric,(record)=>record.vCall.replace(/\*.*$/,"")),strata:distributions(records,metric,(record)=>record.stratum),
      lineageSamples:lineageSampleDistributions(records,metric),
      histogram:bins.map((bin,index)=>({label:`${index*5}–${(index+1)*5}%`,...bin}))};
  }
}
