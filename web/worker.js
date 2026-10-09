// worker.js: one worker's share of a run (see assignInWorkers in
// vdj-assign.js). It receives { module, germlines, reads, worker, workers,
// block } and posts { type: 'ready', info }, then { type: 'rows', block,
// rows, progress } per block of its share, then { type: 'done' }, or
// { type: 'error', message }.

import { load } from './vdj-assign.js';

self.onmessage = async ({ data }) => {
  const { module, germlines, reads, worker, workers, block } = data;
  let run;
  let session;
  try {
    const vdj = await load({ module });
    run = vdj.createRun(germlines);
    const { summary, warnings, backend, header } = run;
    self.postMessage({ type: 'ready', info: { summary, warnings, backend, header } });
    session = run.session(reads, { worker, workers, block });
    for (let next; (next = session.next()); ) {
      self.postMessage({ type: 'rows', block: next.block, rows: next.rows, progress: session.progress });
    }
    self.postMessage({ type: 'done' });
  } catch (error) {
    self.postMessage({ type: 'error', message: error.message });
  } finally {
    try {
      session?.free();
      run?.free();
    } catch {
      // A trapped instance cannot free; the worker is discarded anyway.
    }
  }
};
