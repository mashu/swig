// app.js: the demonstration page around vdj-assign.js. Files are read and
// assigned in this browser; nothing is uploaded.

import { assignInWorkers, compile, load } from './vdj-assign.js';

const $ = (id) => document.getElementById(id);
const SEGMENTS = ['v', 'd', 'j', 'c'];
const VALUE_OPTIONS = ['profile', 'chain', 'calls', 'd-allele', 'd-presence-rate', 'c-evalue'];
const FLAG_OPTIONS = ['tandem-d', 'no-posterior', 'forward-only', 'no-anchor-window', 'residual'];
const PREVIEW_ROWS = 200;
const PREVIEW_COLUMNS = [
  'sequence_id', 'rev_comp', 'v_call', 'd_call', 'j_call', 'c_call', 'junction_aa',
  'productive', 'v_identity', 'v_posterior', 'd_presence_posterior',
];

let module;
let vdj;
let demo = null;
/** Object URLs of the result links and of the demo's links. */
const urls = { results: [], demo: [] };

/** Show `text` as a notice of `kind` ('warn' or 'bad'). */
function notice(kind, text) {
  const p = document.createElement('p');
  p.className = `notice ${kind}`;
  p.textContent = text;
  $('notices').append(p);
}

/**
 * `file` read by `read` (text or bytes) from a Response over its contents,
 * decompressed when its name ends in .gz.
 */
async function contents(file, read) {
  if (!file.name.endsWith('.gz')) return read(new Response(file));
  try {
    return await read(new Response(file.stream().pipeThrough(new DecompressionStream('gzip'))));
  } catch (error) {
    throw new Error(
      `${file.name}: the browser could not decompress it (${error.message}). Its decompressor ` +
        'reads one gzip member, so files of several (bgzip, concatenated lanes) need ' +
        'decompressing first, for example with zcat.',
    );
  }
}

const fileText = (file) => contents(file, (response) => response.text());
const fileBytes = (file) => contents(file, (response) => response.arrayBuffer());

/**
 * Germline FASTA texts: the chosen files when any is chosen, else the
 * demo's. A chosen database is never mixed with the demo's sets.
 */
async function germlines() {
  const chosen = SEGMENTS.some((segment) => $(segment).files[0]);
  const texts = {};
  for (const segment of SEGMENTS) {
    const file = $(segment).files[0];
    texts[segment] = chosen
      ? file && (await fileText(file))
      : demo?.[`${segment.toUpperCase()}.fasta`];
  }
  return texts;
}

/** The reads: a chosen file, else pasted text, else the demo's. */
async function reads() {
  const file = $('reads').files[0];
  if (file) return fileBytes(file);
  const pasted = $('reads-text').value;
  if (pasted.trim()) return pasted;
  return demo?.['reads.fasta'] ?? '';
}

/** The options as named on the command line. */
function options() {
  const chosen = {};
  for (const id of VALUE_OPTIONS) {
    const value = $(id).value.trim();
    if (value) chosen[id] = value;
  }
  for (const id of FLAG_OPTIONS) {
    if ($(id).checked) chosen[id] = true;
  }
  return chosen;
}

/**
 * A download link for `parts` (a string, or strings that form the file) as
 * `name`, its URL kept in `urls[group]`.
 */
function download(name, parts, group) {
  const blob = new Blob(Array.isArray(parts) ? parts : [parts], { type: 'text/plain' });
  const url = URL.createObjectURL(blob);
  urls[group].push(url);
  const a = document.createElement('a');
  a.href = url;
  a.download = name;
  a.textContent = name;
  return a;
}

/** The first rows of the TSV in `parts`, in the columns worth a glance. */
function preview(parts) {
  // Only as many parts as the preview needs are joined.
  let start = '';
  for (const part of parts) {
    start += part;
    if (start.split('\n').length > PREVIEW_ROWS + 1) break;
  }
  const lines = start.split('\n', PREVIEW_ROWS + 1);
  const header = lines[0].split('\t');
  const shown = PREVIEW_COLUMNS.map((name) => [name, header.indexOf(name)]).filter(([, i]) => i >= 0);
  const table = $('preview');
  table.replaceChildren();
  const head = table.createTHead().insertRow();
  for (const [name] of shown) {
    const th = document.createElement('th');
    th.textContent = name;
    head.append(th);
  }
  const body = table.createTBody();
  for (const line of lines.slice(1)) {
    if (!line) continue;
    const fields = line.split('\t');
    const row = body.insertRow();
    for (const [name, i] of shown) {
      const cell = row.insertCell();
      cell.textContent = fields[i] ?? '';
      if (name === 'junction_aa') cell.className = 'mono';
    }
  }
}

/** Release the object URLs of `group`'s links. */
function revoke(group) {
  for (const url of urls[group]) URL.revokeObjectURL(url);
  urls[group] = [];
}

async function run() {
  $('run').disabled = true;
  $('notices').replaceChildren();
  $('results').hidden = true;
  revoke('results');
  const progress = $('progress');
  progress.value = 0;
  progress.hidden = false;
  try {
    const database = { ...(await germlines()), options: options() };
    if (!database.v || !database.j) throw new Error('Choose a V and a J database, or use the demo data.');
    const input = await reads();
    if (!(typeof input === 'string' ? input.trim() : input.byteLength)) {
      throw new Error('Choose a reads file, paste reads, or use the demo data.');
    }
    const workers = Math.min(64, Math.max(1, Number.parseInt($('workers').value, 10) || 1));
    $('status').textContent = 'preparing the database…';
    const started = performance.now();
    const result = await assignInWorkers({
      module,
      germlines: database,
      reads: input,
      workers,
      onInfo: (info) => {
        for (const warning of info.warnings) notice('warn', warning);
        $('status').textContent = 'assigning…';
      },
      onProgress: (fraction) => {
        progress.value = fraction;
        $('status').textContent = `assigning… ${Math.floor(100 * fraction)} %`;
      },
    });
    const seconds = (performance.now() - started) / 1000;
    const { info } = result;
    $('status').textContent = `done in ${seconds.toFixed(1)} s`;
    $('summary').textContent =
      `${result.reads} reads, ${result.called} with a V(D)J call, in ${seconds.toFixed(1)} s ` +
      `(${Math.round(result.reads / seconds)} reads/s on ${workers} worker${workers > 1 ? 's' : ''}, ` +
      `${info.backend} kernel). Database: ${info.summary}.`;
    const links = $('downloads');
    links.replaceChildren(download('calls.tsv', result.parts, 'results'));
    for (const [name, text] of Object.entries(vdj.prepare(database))) {
      links.append(download(name, text, 'results'));
    }
    preview(result.parts);
    $('preview-note').textContent =
      result.reads > PREVIEW_ROWS ? `The first ${PREVIEW_ROWS} of ${result.reads} rows; calls.tsv has them all.` : '';
    $('results').hidden = false;
  } catch (error) {
    notice('bad', error.message);
    $('status').textContent = '';
  } finally {
    progress.hidden = true;
    $('run').disabled = false;
  }
}

function useDemo() {
  const count = Math.min(200000, Math.max(1, Number.parseInt($('demo-reads').value, 10) || 1));
  const preset = $('demo-preset').value;
  demo = vdj.demo({ seed: 1, reads: count, preset });
  $('demo-note').textContent =
    `Synthetic IGH-like database and ${count} simulated ${preset} reads loaded. Chosen germline files ` +
    'replace the whole demo database, and chosen or pasted reads its reads. It is a test fixture: its ' +
    'V genes are random sequence around conserved motifs, so most reads come out unproductive.';
  revoke('demo');
  const links = document.createElement('span');
  links.className = 'links';
  links.append(
    ' ',
    download('reads.fasta', demo['reads.fasta'], 'demo'),
    download('truth.tsv', demo['truth.tsv'], 'demo'),
  );
  $('demo-note').append(links);
}

$('workers').value = String(Math.min(16, navigator.hardwareConcurrency || 1));
try {
  module = await compile();
  vdj = await load({ module });
  $('engine').textContent = `${vdj.version} · ${vdj.backend}`;
  $('run').disabled = false;
  $('run').addEventListener('click', run);
  $('demo').addEventListener('click', useDemo);
} catch (error) {
  $('engine').textContent = 'unavailable';
  notice('bad', `The WebAssembly module did not load: ${error.message}`);
}
