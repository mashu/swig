// End-to-end tests of the page and of assignInWorkers in headless Chromium:
//
//   npm install --no-save playwright && npx playwright install chromium
//   wasm/build.sh && cargo build --release
//   node web/browser-test.mjs target/release/vdj-assign
//
// The page's and the workers' TSV must equal the native command line's.

import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { createServer } from 'node:http';
import { createRequire } from 'node:module';
import { tmpdir } from 'node:os';
import { extname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { load } from './vdj-assign.js';

const { chromium } = createRequire(import.meta.url)('playwright');
const [cli] = process.argv.slice(2);
assert.ok(cli, 'usage: node web/browser-test.mjs path/to/vdj-assign');

const web = fileURLToPath(new URL('.', import.meta.url));
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.wasm': 'application/wasm' };
const server = createServer((request, response) => {
  const path = new URL(request.url, 'http://localhost').pathname;
  const file = join(web, path === '/' ? 'index.html' : path);
  try {
    const body = readFileSync(file);
    response.writeHead(200, { 'content-type': TYPES[extname(file)] ?? 'application/octet-stream' });
    response.end(body);
  } catch {
    response.writeHead(404).end();
  }
});
await new Promise((resolve) => server.listen(0, '127.0.0.1', resolve));
const origin = `http://127.0.0.1:${server.address().port}`;

// The demo data the page generates, written out for the command line.
const node = await load({ bytes: readFileSync(new URL('vdj_assign_simd.wasm', import.meta.url)) });
const dir = mkdtempSync(join(tmpdir(), 'vdj-assign-browser-'));
const write = (name, text) => {
  writeFileSync(join(dir, name), text);
  return join(dir, name);
};
const demoFiles = (reads) => {
  const demo = node.demo({ seed: 1, reads, preset: 'mid' });
  for (const [name, text] of Object.entries(demo)) write(name, text);
  return demo;
};
const native = (args) =>
  execFileSync(cli, [...args, join(dir, 'reads.fasta')], {
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', 'ignore'],
    maxBuffer: 1 << 28,
  });
const germlineArgs = ['v', 'd', 'j', 'c'].flatMap((k) => [`--${k}`, join(dir, `${k.toUpperCase()}.fasta`)]);

const browser = await chromium.launch();
try {
  const page = await browser.newPage();
  const problems = [];
  page.on('pageerror', (error) => problems.push(error.message));
  page.on('console', (message) => message.type() === 'error' && problems.push(message.text()));
  await page.goto(`${origin}/`);
  await page.waitForFunction(() => document.getElementById('engine').textContent !== 'loading…');
  assert.match(await page.textContent('#engine'), /simd128/);

  // assignInWorkers, with and without callbacks, equals the command line.
  const demo = demoFiles(300);
  const expected = native([...germlineArgs, '--tandem-d']);
  for (const callbacks of [false, true]) {
    const result = await page.evaluate(
      async ({ demo, callbacks }) => {
        const { assignInWorkers } = await import('./vdj-assign.js');
        let infos = 0;
        let last = 0;
        const result = await assignInWorkers({
          germlines: { v: demo['V.fasta'], d: demo['D.fasta'], j: demo['J.fasta'], c: demo['C.fasta'], options: { 'tandem-d': true } },
          reads: demo['reads.fasta'],
          workers: 3,
          block: 32,
          ...(callbacks && { onInfo: () => infos++, onProgress: (p) => (last = p) }),
        });
        return { text: result.parts.join(''), reads: result.reads, called: result.called, infos, last };
      },
      { demo, callbacks },
    );
    assert.equal(result.text, expected, `assignInWorkers (callbacks: ${callbacks}) equals the command line`);
    assert.equal(result.reads, 300);
    assert.ok(result.called > 290);
    if (callbacks) assert.deepEqual([result.infos, result.last], [1, 1]);
  }
  const zero = await page.evaluate(async () => {
    const { assignInWorkers } = await import('./vdj-assign.js');
    return assignInWorkers({ germlines: {}, reads: '', workers: 0 }).catch((error) => error.name);
  });
  assert.equal(zero, 'RangeError');

  // The page: demo data, options, workers; its calls.tsv equals the command line.
  demoFiles(500);
  await page.fill('#demo-reads', '500');
  await page.click('#demo');
  await page.fill('#workers', '2');
  await page.check('#tandem-d');
  await page.click('#run');
  await page.waitForSelector('#results:not([hidden])');
  assert.match(await page.textContent('#summary'), /^500 reads, \d+ with a V\(D\)J call/);
  const links = await page.$$eval('#downloads a', (as) => as.map((a) => a.download));
  assert.deepEqual(links, ['calls.tsv', 'V.fasta', 'D.fasta', 'J.fasta', 'C.fasta', 'annotation.tsv']);
  const calls = await page.evaluate(async () => (await fetch(document.querySelector('#downloads a').href)).text());
  assert.equal(calls, native([...germlineArgs, '--tandem-d']), "the page's calls.tsv equals the command line");
  assert.equal(await page.$$eval('#preview tbody tr', (rows) => rows.length), 200);

  // Chosen germline files replace the whole demo database: no demo D or C.
  await page.setInputFiles('#v', join(dir, 'V.fasta'));
  await page.setInputFiles('#j', join(dir, 'J.fasta'));
  await page.uncheck('#tandem-d');
  await page.click('#run');
  await page.waitForSelector('#results:not([hidden])');
  const vj = await page.$$eval('#downloads a', (as) => as.map((a) => a.download));
  assert.deepEqual(vj, ['calls.tsv', 'V.fasta', 'J.fasta', 'annotation.tsv']);
  const vjCalls = await page.evaluate(async () => (await fetch(document.querySelector('#downloads a').href)).text());
  assert.equal(vjCalls, native(['--v', join(dir, 'V.fasta'), '--j', join(dir, 'J.fasta')]));

  // A bad option is a notice, and the page stays usable.
  await page.fill('#d-presence-rate', '1.5');
  await page.click('#run');
  await page.waitForSelector('#notices p.bad');
  assert.match(await page.textContent('#notices p.bad'), /--d-presence-rate expects a rate in \(0, 1\)/);
  assert.equal(await page.isDisabled('#run'), false);
  assert.deepEqual(problems, []);

  // Opened as a file, the page says it needs a web server.
  const file = await browser.newPage();
  await file.goto(new URL('index.html', import.meta.url).href);
  assert.equal(await file.isVisible('#protocol'), true);
} finally {
  await browser.close();
  server.close();
}
console.log('browser tests passed (assignInWorkers and the page equal the native command line)');
