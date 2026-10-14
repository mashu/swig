import { useMemo, useRef, useState, type ReactNode } from "react";

import type {
  AirrResultStore,
  FacetValue,
  RepertoirePair,
  RepertoireSnapshot,
} from "./result-store";
import { tableHeader, tableRow } from "./export-formats";
import { contrastingText, sampleColor, type SampleColorMap } from "./sample-colors";

type Metric = "percent" | "count";
type Palette = "teal" | "coral" | "indigo";
type CallSeries = "vCalls" | "dCalls" | "jCalls" | "cCalls" | "isotypes";

const PALETTES: Record<Palette, { main: string; deep: string; pale: string }> = {
  teal: { main: "#08786f", deep: "#064f4b", pale: "#d8eee7" },
  coral: { main: "#df6854", deep: "#943e34", pale: "#f7ddd5" },
  indigo: { main: "#5869a7", deep: "#35416f", pale: "#e0e4f5" },
};

function downloadSvg(svg: SVGSVGElement | null, name: string) {
  if (!svg) return;
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("width", clone.viewBox.baseVal.width.toString());
  clone.setAttribute("height", clone.viewBox.baseVal.height.toString());
  const source = `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}`;
  const url = URL.createObjectURL(new Blob([source], { type: "image/svg+xml;charset=utf-8" }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name.replace(/[^a-z0-9._-]+/gi, "-").toLowerCase();
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function downloadCsv(rows:Array<Record<string,string|number>>,name:string){if(!rows.length)return;const fields=Object.keys(rows[0]);let text=tableHeader(fields,"csv");for(const row of rows)text+=tableRow(fields,row,"csv");const url=URL.createObjectURL(new Blob([text],{type:"text/csv;charset=utf-8"}));const anchor=document.createElement("a");anchor.href=url;anchor.download=name.replace(/\.svg$/i,".csv");anchor.click();window.setTimeout(()=>URL.revokeObjectURL(url),1000);}

function geneName(value: string): string {
  return value.replace(/\*[^,;]+/g, "");
}

function collapse(values: FacetValue[], resolution: "allele" | "gene"): FacetValue[] {
  if (resolution === "allele") return values;
  const combined = new Map<string, number>();
  for (const item of values) {
    const value = geneName(item.value);
    combined.set(value, (combined.get(value) ?? 0) + item.count);
  }
  return [...combined.entries()]
    .map(([value, count]) => ({ value, count }))
    .sort((a, b) => b.count - a.count || a.value.localeCompare(b.value));
}

function valueFor(count: number, total: number, metric: Metric): number {
  return metric === "percent" ? (total ? count / total * 100 : 0) : count;
}

function valueLabel(value: number, metric: Metric): string {
  if (metric === "percent") return `${value < 0.1 && value > 0 ? value.toFixed(2) : value.toFixed(1)}%`;
  return Number.isInteger(value) ? value.toLocaleString() : value.toLocaleString(undefined, { maximumFractionDigits: 2 });
}

function SvgButton({ svg, name }: { svg: React.RefObject<SVGSVGElement | null>; name: string }) {
  return <button className="svg-download" type="button" onClick={() => downloadSvg(svg.current, name)}>Download SVG <span>↓</span></button>;
}

function FrequencyChart({ data, total, metric, palette, title, subtitle, filename }: {
  data: FacetValue[];
  total: number;
  metric: Metric;
  palette: Palette;
  title: string;
  subtitle: string;
  filename: string;
}) {
  const svg = useRef<SVGSVGElement>(null);
  const colors = PALETTES[palette];
  const plotted = data.map((item) => ({ ...item, plotted: valueFor(item.count, total, metric) }));
  const maximum = Math.max(...plotted.map((item) => item.plotted), 1);
  const width = 960;
  const height = Math.max(260, 132 + plotted.length * 32);
  const chartLeft = 250;
  const chartWidth = 570;
  return (
    <article className="figure-card">
      <header><div><span className="section-kicker">Ranked composition</span><h3>{title}</h3><p>{subtitle}</p></div><div className="post-chart-actions"><button type="button" onClick={()=>downloadCsv(plotted.map(item=>({call:item.value,count:item.count,value:item.plotted,measure:metric})),filename)}>Data CSV ↓</button><SvgButton svg={svg} name={filename} /></div></header>
      <div className="svg-frame">
        {plotted.length ? <svg ref={svg} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
          <title>{title}</title><desc>{subtitle}</desc>
          <rect width={width} height={height} fill="#fffdf7" />
          <text x="28" y="42" fill="#132321" fontFamily="Inter,Arial,sans-serif" fontSize="22" fontWeight="750">{title}</text>
          <text x="28" y="66" fill="#677570" fontFamily="Inter,Arial,sans-serif" fontSize="11">{subtitle}</text>
          {[0, .25, .5, .75, 1].map((fraction) => {
            const x = chartLeft + chartWidth * fraction;
            return <g key={fraction}><line x1={x} x2={x} y1="96" y2={height - 30} stroke="#d9d8cf" strokeWidth="1" /><text x={x} y="88" textAnchor={fraction ? "middle" : "start"} fill="#82908b" fontFamily="Inter,Arial,sans-serif" fontSize="9">{valueLabel(maximum * fraction, metric)}</text></g>;
          })}
          {plotted.map((item, index) => {
            const y = 106 + index * 32;
            const barWidth = Math.max(2, item.plotted / maximum * chartWidth);
            return <g key={item.value}>
              <text x={chartLeft - 15} y={y + 14} textAnchor="end" fill="#263936" fontFamily="ui-monospace,SFMono-Regular,Consolas,monospace" fontSize="10">{item.value}</text>
              <rect x={chartLeft} y={y} width={barWidth} height="20" rx="3" fill={colors.main} />
              <text x={Math.min(chartLeft + barWidth + 8, 900)} y={y + 14} fill={colors.deep} fontFamily="Inter,Arial,sans-serif" fontSize="10" fontWeight="700">{valueLabel(item.plotted, metric)}</text>
            </g>;
          })}
          <text x={chartLeft + chartWidth / 2} y={height - 7} textAnchor="middle" fill="#65736e" fontFamily="Inter,Arial,sans-serif" fontSize="10">{metric === "percent" ? "Percent of analyzed records" : "Record count"}</text>
        </svg> : <div className="figure-empty"><strong>No calls to plot</strong><span>Try another locus, population, or segment.</span></div>}
      </div>
    </article>
  );
}

function DistributionChart({ snapshot, series, metric, palette, filename }: {
  snapshot: RepertoireSnapshot;
  series: "cdr3Lengths" | "vIdentityBins";
  metric: Metric;
  palette: Palette;
  filename: string;
}) {
  const svg = useRef<SVGSVGElement>(null);
  const colors = PALETTES[palette];
  const source = snapshot[series];
  const title = series === "cdr3Lengths" ? "CDR3 amino-acid length" : "V-region identity";
  const subtitle = series === "cdr3Lengths" ? "Distribution of called CDR3 lengths." : "One-percentage-point bins from germline V identity.";
  const width = 960;
  const height = 420;
  const left = 72;
  const right = 28;
  const top = 98;
  const bottom = 62;
  const values = source.map((item) => ({ ...item, x: Number(item.value), plotted: valueFor(item.count, snapshot.records, metric) }));
  const maximum = Math.max(...values.map((item) => item.plotted), 1);
  const minX = values.length ? Math.min(...values.map((item) => item.x)) : 0;
  const maxX = values.length ? Math.max(...values.map((item) => item.x)) : 1;
  const range = Math.max(1, maxX - minX + 1);
  const plotWidth = width - left - right;
  const plotHeight = height - top - bottom;
  const barWidth = Math.max(2, plotWidth / range - 1);
  return (
    <article className="figure-card">
      <header><div><span className="section-kicker">Distribution</span><h3>{title}</h3><p>{subtitle}</p></div><div className="post-chart-actions"><button type="button" onClick={()=>downloadCsv(values.map(item=>({bin:item.value,count:item.count,value:item.plotted,measure:metric})),filename)}>Data CSV ↓</button><SvgButton svg={svg} name={filename} /></div></header>
      <div className="svg-frame">
        {values.length ? <svg ref={svg} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
          <title>{title}</title><desc>{subtitle}</desc><rect width={width} height={height} fill="#fffdf7" />
          <text x="28" y="42" fill="#132321" fontFamily="Inter,Arial,sans-serif" fontSize="22" fontWeight="750">{title}</text>
          <text x="28" y="66" fill="#677570" fontFamily="Inter,Arial,sans-serif" fontSize="11">{subtitle}</text>
          {[0, .25, .5, .75, 1].map((fraction) => {
            const y = top + plotHeight * (1 - fraction);
            return <g key={fraction}><line x1={left} x2={width - right} y1={y} y2={y} stroke="#d9d8cf" /><text x={left - 10} y={y + 3} textAnchor="end" fill="#82908b" fontFamily="Inter,Arial,sans-serif" fontSize="9">{valueLabel(maximum * fraction, metric)}</text></g>;
          })}
          {values.map((item) => {
            const x = left + (item.x - minX) / range * plotWidth;
            const barHeight = item.plotted / maximum * plotHeight;
            return <rect key={item.value} x={x} y={top + plotHeight - barHeight} width={barWidth} height={barHeight} fill={colors.main}><title>{item.value}: {valueLabel(item.plotted, metric)}</title></rect>;
          })}
          {Array.from({ length: 6 }, (_, index) => Math.round(minX + (maxX - minX) * index / 5)).map((tick) => {
            const x = left + (tick - minX) / range * plotWidth;
            return <text key={tick} x={x} y={height - bottom + 22} textAnchor="middle" fill="#65736e" fontFamily="Inter,Arial,sans-serif" fontSize="10">{tick}{series === "vIdentityBins" ? "%" : ""}</text>;
          })}
          <text x={left + plotWidth / 2} y={height - 12} textAnchor="middle" fill="#65736e" fontFamily="Inter,Arial,sans-serif" fontSize="10">{series === "cdr3Lengths" ? "CDR3 length (aa)" : "V identity"}</text>
        </svg> : <div className="figure-empty"><strong>No distribution to plot</strong><span>This view requires a called CDR3 or V alignment.</span></div>}
      </div>
    </article>
  );
}

function PairingChart({ pairs, resolution, palette, filename }: {
  pairs: RepertoirePair[];
  resolution: "allele" | "gene";
  palette: Palette;
  filename: string;
}) {
  const svg = useRef<SVGSVGElement>(null);
  const colors = PALETTES[palette];
  const aggregated = useMemo(() => {
    const map = new Map<string, number>();
    for (const pair of pairs) {
      const v = resolution === "gene" ? geneName(pair.v) : pair.v;
      const j = resolution === "gene" ? geneName(pair.j) : pair.j;
      const key = `${v}\u0000${j}`;
      map.set(key, (map.get(key) ?? 0) + pair.count);
    }
    return [...map.entries()].map(([key, count]) => {
      const [v, j] = key.split("\u0000");
      return { v, j, count };
    });
  }, [pairs, resolution]);
  const vTotals = new Map<string, number>();
  const jTotals = new Map<string, number>();
  aggregated.forEach((pair) => {
    vTotals.set(pair.v, (vTotals.get(pair.v) ?? 0) + pair.count);
    jTotals.set(pair.j, (jTotals.get(pair.j) ?? 0) + pair.count);
  });
  const top = (map: Map<string, number>, limit: number) => [...map.entries()].sort((a, b) => b[1] - a[1]).slice(0, limit).map(([value]) => value);
  const vs = top(vTotals, 12);
  const js = top(jTotals, 12);
  const matrix = aggregated.filter((pair) => vs.includes(pair.v) && js.includes(pair.j));
  const maximum = Math.max(...matrix.map((pair) => pair.count), 1);
  const width = 960;
  const left = 240;
  const topOffset = 150;
  const cell = Math.min(48, (width - left - 35) / Math.max(1, js.length));
  const height = Math.max(390, topOffset + vs.length * cell + 70);
  const title = "V–J pairing within rearrangements";
  const subtitle = "Circle area is proportional to record count; axes retain the 12 most frequent calls.";
  return (
    <article className="figure-card figure-card-wide">
      <header><div><span className="section-kicker">Pair structure</span><h3>{title}</h3><p>{subtitle}</p></div><div className="post-chart-actions"><button type="button" onClick={()=>downloadCsv(aggregated.map(pair=>({v_call:pair.v,j_call:pair.j,count:pair.count})),filename)}>Data CSV ↓</button><SvgButton svg={svg} name={filename} /></div></header>
      <div className="svg-frame">
        {matrix.length ? <svg ref={svg} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}>
          <title>{title}</title><desc>{subtitle}</desc><rect width={width} height={height} fill="#fffdf7" />
          <text x="28" y="42" fill="#132321" fontFamily="Inter,Arial,sans-serif" fontSize="22" fontWeight="750">{title}</text>
          <text x="28" y="66" fill="#677570" fontFamily="Inter,Arial,sans-serif" fontSize="11">{subtitle}</text>
          {js.map((j, index) => <text key={j} x={left + index * cell + cell / 2} y={topOffset - 12} transform={`rotate(-42 ${left + index * cell + cell / 2} ${topOffset - 12})`} textAnchor="start" fill="#364844" fontFamily="ui-monospace,SFMono-Regular,Consolas,monospace" fontSize="9">{j}</text>)}
          {vs.map((v, row) => <g key={v}><text x={left - 14} y={topOffset + row * cell + cell / 2 + 3} textAnchor="end" fill="#364844" fontFamily="ui-monospace,SFMono-Regular,Consolas,monospace" fontSize="9">{v}</text>{js.map((j, column) => <rect key={j} x={left + column * cell} y={topOffset + row * cell} width={cell} height={cell} fill={(row + column) % 2 ? "#f4f1e8" : "#faf8f1"} stroke="#e2e0d7" />)}</g>)}
          {matrix.map((pair) => {
            const column = js.indexOf(pair.j);
            const row = vs.indexOf(pair.v);
            const radius = Math.max(2.5, Math.sqrt(pair.count / maximum) * cell * .38);
            return <circle key={`${pair.v}-${pair.j}`} cx={left + column * cell + cell / 2} cy={topOffset + row * cell + cell / 2} r={radius} fill={colors.main} fillOpacity=".82" stroke={colors.deep} strokeWidth="1"><title>{pair.v} × {pair.j}: {pair.count.toLocaleString()}</title></circle>;
          })}
          <text x={left - 170} y={topOffset + vs.length * cell / 2} transform={`rotate(-90 ${left - 170} ${topOffset + vs.length * cell / 2})`} textAnchor="middle" fill="#65736e" fontFamily="Inter,Arial,sans-serif" fontSize="10">V call</text>
          <text x={left + js.length * cell / 2} y={height - 16} textAnchor="middle" fill="#65736e" fontFamily="Inter,Arial,sans-serif" fontSize="10">J call</text>
        </svg> : <div className="figure-empty"><strong>No V–J pairs to plot</strong><span>Both segments must be assigned in the selected population.</span></div>}
      </div>
    </article>
  );
}

function SampleCompositionChart({ samples, colors, filename }: { samples: FacetValue[]; colors: SampleColorMap; filename: string }) {
  const svg = useRef<SVGSVGElement>(null);
  const visible = samples.slice(0, 40);
  const maximum = Math.max(1, ...visible.map((item) => item.count));
  const width = 960;
  const height = Math.max(250, 112 + visible.length * 29);
  const left = 250;
  const bar = 570;
  const title = "Records by biological sample";
  return <article className="figure-card figure-card-wide"><header><div><span className="section-kicker">Study composition</span><h3>{title}</h3><p>Colors follow the study-wide sample palette and remain fixed in downstream trees and sample-stratified figures.</p></div><div className="post-chart-actions"><button type="button" onClick={()=>downloadCsv(samples.map((item)=>({sample_id:item.value,records:item.count,color:sampleColor(item.value,colors)})),filename)}>Data CSV ↓</button><SvgButton svg={svg} name={filename}/></div></header><div className="svg-frame">{visible.length?<svg ref={svg} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={title}><rect width={width} height={height} fill="#fffdf7"/><text x="28" y="42" fill="#132321" fontFamily="Inter,Arial,sans-serif" fontSize="22" fontWeight="750">{title}</text>{visible.map((item,index)=>{const y=76+index*29;const color=sampleColor(item.value,colors);const length=Math.max(3,item.count/maximum*bar);return <g key={item.value}><text x={left-14} y={y+15} textAnchor="end" fill="#263936" fontFamily="ui-monospace,monospace" fontSize="10">{item.value}</text><rect x={left} y={y} width={length} height="20" rx="3" fill={color}/><text x={length>70?left+8:left+length+8} y={y+14} fill={length>70?contrastingText(color):"#263936"} fontFamily="Inter,Arial,sans-serif" fontSize="10" fontWeight="700">{item.count.toLocaleString()}</text></g>;})}</svg>:<div className="figure-empty"><strong>No sample metadata</strong></div>}</div></article>;
}

const SERIES_LABELS: Record<CallSeries, string> = {
  vCalls: "V use",
  dCalls: "D use",
  jCalls: "J use",
  cCalls: "Constant-region calls",
  isotypes: "Isotype / constant class",
};

export function RepertoireDashboard({ store, loci, inputName, samples, sampleColors, sidebarTools, overview }: {
  store: AirrResultStore;
  loci: FacetValue[];
  inputName: string;
  samples: FacetValue[];
  sampleColors: SampleColorMap;
  sidebarTools?: ReactNode;
  overview?: ReactNode;
}) {
  const [locus, setLocus] = useState("");
  const [productiveOnly, setProductiveOnly] = useState(false);
  const [ambiguity, setAmbiguity] = useState<"top" | "fractional">("fractional");
  const [resolution, setResolution] = useState<"allele" | "gene">("allele");
  const [metric, setMetric] = useState<Metric>("percent");
  const [palette, setPalette] = useState<Palette>("teal");
  const [series, setSeries] = useState<CallSeries>("vCalls");
  const [topN, setTopN] = useState(20);
  const [distribution, setDistribution] = useState<"cdr3Lengths" | "vIdentityBins">("cdr3Lengths");
  const [panel, setPanel] = useState<"usage" | "distribution" | "pairing" | "samples">("usage");
  const snapshot = useMemo(
    () => store.repertoire({ locus, productiveOnly, ambiguity }),
    [store, locus, productiveOnly, ambiguity],
  );
  const ranked = useMemo(() => {
    const values = series === "isotypes" ? snapshot.isotypes : collapse(snapshot[series], resolution);
    return [...values].sort((a, b) => b.count - a.count || a.value.localeCompare(b.value)).slice(0, topN);
  }, [resolution, series, snapshot, topN]);
  const fileStem = inputName.replace(/(\.gz)?\.[^.]+$/, "") || "swig";
  const population = productiveOnly ? "productive records" : "all analyzed records";
  const callDetail = ambiguity === "fractional" ? "co-optimal calls split one record's weight equally" : "only the first reported call counted";

  return (
    <section className="repertoire-workspace contextual-workspace">
      <nav className="context-rail" aria-label="Repertoire panels">
        <div className="context-rail-heading"><span>Repertoire</span><small>{snapshot.records.toLocaleString()} records</small></div>
        <button type="button" className={panel==="usage"?"active":""} onClick={()=>setPanel("usage")}><b>01</b><span>Gene use<small>{SERIES_LABELS[series]}</small></span></button>
        <button type="button" className={panel==="distribution"?"active":""} onClick={()=>setPanel("distribution")}><b>02</b><span>Distribution<small>{distribution==="cdr3Lengths"?"CDR3 length":"V identity"}</small></span></button>
        <button type="button" className={panel==="pairing"?"active":""} onClick={()=>setPanel("pairing")}><b>03</b><span>V–J pairs<small>{resolution} labels</small></span></button>
        {samples.length>1&&<button type="button" className={panel==="samples"?"active":""} onClick={()=>setPanel("samples")}><b>04</b><span>Samples<small>{samples.length} observed</small></span></button>}
        {sidebarTools}
      </nav>
      <div className="context-main repertoire-dashboard">
        {overview}
        <header className="repertoire-heading"><div><span className="section-kicker">Repertoire analysis</span><h2>{panel==="usage"?SERIES_LABELS[series]:panel==="distribution"?distribution==="cdr3Lengths"?"CDR3 amino-acid length":"V-region identity":panel==="pairing"?"V–J call pairs":"Sample composition"}</h2><p>Frequencies are calculated per input record. Population and ambiguous-call handling can be changed on this page.</p></div><span className="aggregate-badge">{snapshot.records.toLocaleString()} records in view</span></header>

        <details className="repertoire-global-settings"><summary><span><b>Population and figure settings</b><small>{population} · {locus||"all loci"} · {ambiguity==="fractional"?"fractional ties":"first calls"} · {resolution} labels · {metric}</small></span></summary><div className="figure-controls" aria-label="Repertoire figure controls">
          <label><span>Population</span><select value={productiveOnly ? "productive" : "all"} onChange={(event) => setProductiveOnly(event.target.value === "productive")}><option value="all">All records</option><option value="productive">Productive only</option></select></label>
          <label><span>Locus</span><select value={locus} onChange={(event) => setLocus(event.target.value)}><option value="">All loci</option>{loci.map((item) => <option key={item.value} value={item.value}>{item.value} ({item.count.toLocaleString()})</option>)}</select></label>
          <label><span>Ambiguous calls</span><select value={ambiguity} onChange={(event) => setAmbiguity(event.target.value as "top" | "fractional")}><option value="fractional">Fractional ties</option><option value="top">First call only</option></select></label>
          <label><span>Labels</span><select value={resolution} onChange={(event) => setResolution(event.target.value as "allele" | "gene")}><option value="allele">Alleles</option><option value="gene">Collapse to genes</option></select></label>
          <label><span>Y measure</span><select value={metric} onChange={(event) => setMetric(event.target.value as Metric)}><option value="percent">Percent of records</option><option value="count">Record count</option></select></label>
          <label><span>Palette</span><select value={palette} onChange={(event) => setPalette(event.target.value as Palette)}><option value="teal">Teal</option><option value="coral">Coral</option><option value="indigo">Indigo</option></select></label>
        </div></details>

        {panel==="usage"&&<div className="figure-grid single-panel"><div className="figure-config"><label><span>Ranked series</span><select value={series} onChange={(event) => setSeries(event.target.value as CallSeries)}>{Object.entries(SERIES_LABELS).map(([value, label]) => <option key={value} value={value}>{label}</option>)}</select></label><label><span>Show</span><select value={topN} onChange={(event) => setTopN(Number(event.target.value))}><option value="10">Top 10</option><option value="20">Top 20</option><option value="40">Top 40</option></select></label></div><FrequencyChart data={ranked} total={snapshot.records} metric={metric} palette={palette} title={SERIES_LABELS[series]} subtitle={`${population}; ${series === "isotypes" ? "C anchored after J, ≥12 aligned bases at ≥90% identity" : callDetail}.`} filename={`${fileStem}-${series}.svg`} /></div>}
        {panel==="distribution"&&<div className="figure-grid single-panel"><div className="figure-config single"><label><span>Distribution</span><select value={distribution} onChange={(event) => setDistribution(event.target.value as "cdr3Lengths" | "vIdentityBins")}><option value="cdr3Lengths">CDR3 AA length</option><option value="vIdentityBins">V identity</option></select></label></div><DistributionChart snapshot={snapshot} series={distribution} metric={metric} palette={palette} filename={`${fileStem}-${distribution}.svg`} /></div>}
        {panel==="pairing"&&<div className="figure-grid single-panel"><PairingChart pairs={snapshot.vjPairs} resolution={resolution} palette={palette} filename={`${fileStem}-vj-pairing.svg`} /></div>}
        {panel==="samples"&&samples.length>1&&<div className="figure-grid single-panel"><SampleCompositionChart samples={samples} colors={sampleColors} filename={`${fileStem}-samples.svg`}/></div>}
        <p className="figure-method-note"><span>i</span> Gene-use and proportional pairing views summarize individual rearrangements. V–J pairs do not imply heavy–light pairing unless paired records are explicitly modeled.</p>
      </div>
    </section>
  );
}
