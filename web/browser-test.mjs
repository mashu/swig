// End-to-end tests of the page and of assignInWorkers in headless browsers:
//
//   npm install --no-save playwright && npx playwright install chromium [firefox webkit]
//   wasm/build.sh && cargo build --release
//   node web/browser-test.mjs target/release/vdj-assign [chromium,firefox,webkit]
//
// index.html is opened from disk, as a file: page with no server, and from a
// static HTTP server. Its calls.tsv, and assignInWorkers' TSV, must equal the
// native command line's; opened from disk, it may request nothing but itself
// and the blobs it makes. Each browser named (default chromium) must be
// installed.

import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { createServer } from 'node:http';
import { createRequire } from 'node:module';
import { tmpdir } from 'node:os';
import { extname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { load } from './vdj-assign.js';

const playwright = createRequire(import.meta.url)('playwright');
const [cli, names = 'chromium'] = process.argv.slice(2);
assert.ok(cli, 'usage: node web/browser-test.mjs path/to/vdj-assign [chromium,firefox,webkit]');

const web = fileURLToPath(new URL('.', import.meta.url));
const page = new URL('index.html', import.meta.url).href;
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

/** A page that records its errors in `problems` and its requests in `requests`. */
async function open(browser, url, problems, requests = []) {
  const tab = await browser.newPage();
  tab.on('pageerror', (error) => problems.push(error.message));
  tab.on('console', (message) => message.type() === 'error' && problems.push(message.text()));
  tab.on('request', (request) => requests.push(request.url()));
  await tab.goto(url);
  await tab.waitForFunction(() => document.getElementById('engine').textContent !== 'loading…');
  assert.match(await tab.textContent('#engine'), /simd128/);
  return tab;
}

/** The text of the file the first download link of the page saves. */
async function firstDownload(tab) {
  const [download] = await Promise.all([tab.waitForEvent('download'), tab.click('#downloads a')]);
  return readFileSync(await download.path(), 'utf8');
}

/** The page assigns 500 demo reads (tandem D, 2 workers) as the command line does. */
async function assignDemo(tab) {
  demoFiles(500);
  await tab.fill('#demo-reads', '500');
  await tab.click('#demo');
  await tab.fill('#workers', '2');
  await tab.check('#tandem-d');
  await tab.click('#run');
  await tab.waitForSelector('#results:not([hidden])');
  assert.match(await tab.textContent('#summary'), /^500 reads, \d+ with a V\(D\)J call/);
  const links = await tab.$$eval('#downloads a', (as) => as.map((a) => a.download));
  assert.deepEqual(links, ['calls.tsv', 'V.fasta', 'D.fasta', 'J.fasta', 'C.fasta', 'annotation.tsv']);
  assert.equal(await firstDownload(tab), native([...germlineArgs, '--tandem-d']), "the page's calls.tsv equals the command line");
  assert.equal(await tab.$$eval('#preview tbody tr', (rows) => rows.length), 200);
}

try {
  for (const name of names.split(',')) {
    const browser = await playwright[name].launch();
    try {
      // Opened from disk: no server, and nothing requested but the page and its blobs.
      const problems = [];
      const requests = [];
      const fromDisk = await open(browser, page, problems, requests);
      await assignDemo(fromDisk);
      const elsewhere = requests.filter((url) => url !== page && !/^(blob|data):/.test(url));
      assert.deepEqual(elsewhere, [], `${name}: the page opened from disk requests only blobs`);

      // From a static server: the page, and assignInWorkers through the module API.
      const tab = await open(browser, `${origin}/`, problems);
      await assignDemo(tab);
      const demo = demoFiles(300);
      const expected = native([...germlineArgs, '--tandem-d']);
      for (const callbacks of [false, true]) {
        const result = await tab.evaluate(
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
        assert.equal(result.text, expected, `${name}: assignInWorkers (callbacks: ${callbacks}) equals the command line`);
        assert.equal(result.reads, 300);
        assert.ok(result.called > 290);
        if (callbacks) assert.deepEqual([result.infos, result.last], [1, 1]);
      }
      const zero = await tab.evaluate(async () => {
        const { assignInWorkers } = await import('./vdj-assign.js');
        return assignInWorkers({ germlines: {}, reads: '', workers: 0 }).catch((error) => error.name);
      });
      assert.equal(zero, 'RangeError');

      // Chosen germline files replace the whole demo database: no demo D or C.
      demoFiles(500);
      await tab.setInputFiles('#v', join(dir, 'V.fasta'));
      await tab.setInputFiles('#j', join(dir, 'J.fasta'));
      await tab.uncheck('#tandem-d');
      await tab.click('#run');
      await tab.waitForSelector('#results:not([hidden])');
      const vj = await tab.$$eval('#downloads a', (as) => as.map((a) => a.download));
      assert.deepEqual(vj, ['calls.tsv', 'V.fasta', 'J.fasta', 'annotation.tsv']);
      assert.equal(await firstDownload(tab), native(['--v', join(dir, 'V.fasta'), '--j', join(dir, 'J.fasta')]));

      // A bad option is a notice, and the page stays usable.
      await tab.fill('#d-presence-rate', '1.5');
      await tab.click('#run');
      await tab.waitForSelector('#notices p.bad');
      assert.match(await tab.textContent('#notices p.bad'), /--d-presence-rate expects a rate in \(0, 1\)/);
      assert.equal(await tab.isDisabled('#run'), false);
      assert.deepEqual(problems, [], `${name}: no page errors`);
      console.log(`${name}: the page from disk and from a server, and assignInWorkers, equal the command line`);
    } finally {
      await browser.close();
    }
  }
} finally {
  server.close();
}
console.log(`browser tests passed (${names})`);
