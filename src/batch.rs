//! Parallel batch assignment with per-thread workspaces.
//!
//! **Trick — workspace per thread, dynamic chunking.** Each worker owns one
//! `VdjAssigner::Workspace` for its whole lifetime, so after warm-up the hot
//! path performs no DP-buffer allocation. Work is handed out in chunks from an
//! atomic counter (self-scheduling), which balances load when some queries are
//! much more expensive than others without a work-stealing runtime. Each chunk
//! goes to [`VdjAssigner::assign_many`], which batches every pipeline stage
//! across the chunk so SIMD aligners can fill their lanes. Results keep input
//! order.
//!
//! The calling thread is one of the workers. Helpers are spawned for the
//! others, and a target that cannot spawn threads (WebAssembly without
//! shared memory) runs every chunk on the caller.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use crate::pipeline::{Call, VdjAssigner};

const CHUNK: usize = 64;

/// Assign every query and map each call to an output on the worker thread
/// (so formatting also runs in parallel). `threads == 0` uses all cores.
pub fn assign_map<T, Q, O, F>(assigner: &T, queries: &[Q], threads: usize, map: F) -> Vec<O>
where
    T: VdjAssigner + Sync,
    Q: AsRef<[u8]> + Sync,
    O: Send,
    F: Fn(usize, Call) -> O + Sync,
{
    let threads = match threads {
        0 => thread::available_parallelism().map_or(1, |n| n.get()),
        n => n,
    }
    .min(queries.len().div_ceil(CHUNK).max(1));
    let next = AtomicUsize::new(0);
    let work = || {
        let mut workspace = assigner.workspace();
        let mut calls = Vec::with_capacity(CHUNK);
        let mut done = Vec::new();
        loop {
            let start = next.fetch_add(CHUNK, Ordering::Relaxed);
            if start >= queries.len() {
                break done;
            }
            let end = (start + CHUNK).min(queries.len());
            calls.clear();
            assigner.assign_many(&queries[start..end], &mut workspace, &mut calls);
            let outputs: Vec<O> = calls
                .drain(..)
                .enumerate()
                .map(|(offset, call)| map(start + offset, call))
                .collect();
            done.push((start, outputs));
        }
    };
    let mut chunks: Vec<(usize, Vec<O>)> = thread::scope(|scope| {
        let helpers: Vec<_> = (1..threads)
            .map_while(|_| thread::Builder::new().spawn_scoped(scope, work).ok())
            .collect();
        let mut chunks = work();
        for helper in helpers {
            chunks.extend(helper.join().expect("assignment worker panicked"));
        }
        chunks
    });
    chunks.sort_unstable_by_key(|&(start, _)| start);
    chunks
        .into_iter()
        .flat_map(|(_, outputs)| outputs)
        .collect()
}

/// Assign every query, returning calls in input order.
pub fn assign_all<T, Q>(assigner: &T, queries: &[Q], threads: usize) -> Vec<Call>
where
    T: VdjAssigner + Sync,
    Q: AsRef<[u8]> + Sync,
{
    assign_map(assigner, queries, threads, |_, call| call)
}
