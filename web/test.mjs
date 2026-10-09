// Tests of the WebAssembly modules through vdj-assign.js, in Node.js:
//
//   wasm/build.sh && node web/test.mjs [path/to/vdj-assign]
//
// Given the native binary, the modules' TSV must equal its output byte for
// byte on the same files.

import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { load, simdSupported } from './vdj-assign.js';

const here = new URL('.', import.meta.url);
const bytes = (file) => readFileSync(new URL(file, here));
const [cli] = process.argv.slice(2);

assert.ok(simdSupported(), 'Node.js runs WebAssembly SIMD');
const simd = await load({ bytes: bytes('vdj_assign_simd.wasm') });
const baseline = await load({ bytes: bytes('vdj_assign.wasm') });
assert.equal(simd.backend, 'simd128');
assert.equal(baseline.backend, 'scalar');

const demo = simd.demo({ seed: 7, reads: 400, preset: 'hard' });
const germlines = { v: demo['V.fasta'], d: demo['D.fasta'], j: demo['J.fasta'], c: demo['C.fasta'] };
assert.equal(demo['reads.fasta'].split('\n').filter((line) => line.startsWith('>')).length, 400);
assert.deepEqual(baseline.demo({ seed: 7, reads: 400, preset: 'hard' }), demo);

const options = { 'tandem-d': true, 'c-evalue': 0.001 };
const tsv = {};
for (const [name, vdj] of [['simd', simd], ['baseline', baseline]]) {
  const run = vdj.createRun({ ...germlines, options });
  assert.match(run.summary, /^V \d+ alleles, D \d+, J \d+, C \d+/);
  assert.equal(run.backend, vdj.backend);
  tsv[name] = run.assign(demo['reads.fasta']);

  // Three workers' blocks, in order, are the whole output.
  const blocks = [];
  for (let worker = 0; worker < 3; worker++) {
    const session = run.session(demo['reads.fasta'], { worker, workers: 3, block: 32 });
    for (let next; (next = session.next()); ) blocks[next.block] = next.rows;
    assert.equal(session.progress, 1);
    session.free();
  }
  assert.equal(`${run.header}\n${blocks.join('')}`, tsv[name]);
  run.free();
}
const lines = tsv.simd.trimEnd().split('\n');
assert.equal(lines.length, 401);
assert.ok(lines[0].startsWith('sequence_id\tsequence\trev_comp\t'));
assert.ok(lines[0].includes('\tc_call\t'));
assert.equal(tsv.baseline, tsv.simd, 'every kernel gives the same rows');

// FASTQ as bytes, and an error that names its line.
const run = simd.createRun(germlines);
const fastq = new TextEncoder().encode('@r1\nACGTACGTAC\n+\nIIIIIIIIII\n');
assert.equal(run.assign(fastq).trimEnd().split('\n').length, 2);
assert.throws(() => run.assign('@r1\nAC GT\n+\nIIIII\n'), /line 2: ' ' in a FASTQ sequence line/);
run.free();

// Errors are exceptions, and the instance stays usable after them.
assert.throws(() => simd.createRun({ ...germlines, options: { profile: 'best' } }), /--profile expects/);
assert.throws(() => simd.createRun({ v: germlines.v }), /a J database is required/);
assert.throws(() => simd.createRun({ ...germlines, v: '@r\nACGT\n+\nIIII\n' }), /V database: .*FASTQ/);
assert.throws(() => simd.demo({ preset: 'impossible' }), /unknown preset/);
const prepared = simd.prepare(germlines);
assert.deepEqual(Object.keys(prepared), ['V.fasta', 'D.fasta', 'J.fasta', 'C.fasta', 'annotation.tsv']);
assert.ok(prepared['annotation.tsv'].startsWith('segment\tallele\t'));

if (cli) {
  const dir = mkdtempSync(join(tmpdir(), 'vdj-assign-'));
  for (const name of ['V.fasta', 'D.fasta', 'J.fasta', 'C.fasta', 'reads.fasta']) {
    writeFileSync(join(dir, name), demo[name]);
  }
  const files = ['v', 'd', 'j', 'c'].flatMap((k) => [`--${k}`, join(dir, `${k.toUpperCase()}.fasta`)]);
  const native = execFileSync(cli, [...files, '--tandem-d', '--c-evalue', '0.001', join(dir, 'reads.fasta')], {
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', 'ignore'],
    maxBuffer: 1 << 28,
  });
  assert.equal(tsv.simd, native, 'the module gives the command line output');
  console.log('the modules equal the native command line byte for byte');
}
console.log(`web tests passed (vdj_assign ${simd.version}; ${lines.length - 1} reads; simd128 and scalar)`);
