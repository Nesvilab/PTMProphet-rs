// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Minimal data-parallelism helpers over `std::thread::scope` (no external runtime / no rayon).
//!
//! Each helper spawns a **fixed pool** of `worker_count()` OS threads that pull work *dynamically*
//! from a shared queue of small contiguous chunks (a mutex-guarded `chunks`/`chunks_mut` iterator).
//! Dynamic pulling keeps every core busy even when per-item cost is wildly uneven — e.g. PSMs whose
//! placement count (and so match cost) varies by orders of magnitude — unlike a static
//! one-chunk-per-thread split, where a single heavy chunk strands the other cores idle until it
//! finishes. Results are reassembled in **input order**, so the produced numbers/bytes are identical
//! whether run on 1 thread or N: callers stay fully deterministic regardless of the worker count.

use std::sync::Mutex;
use std::sync::OnceLock;

/// Explicit worker count from the `--threads`/`-j` CLI option, set once at startup.
/// The `PTMPROPHET_RS_THREADS` env var still takes precedence over this (an escape hatch).
static THREAD_OVERRIDE: OnceLock<usize> = OnceLock::new();

/// Record the `--threads N` value chosen on the command line. Clamped to ≥ 1; the first call wins.
pub fn set_thread_override(n: usize) {
    let _ = THREAD_OVERRIDE.set(n.max(1));
}

/// Worker-thread count, in precedence order:
/// 1. `PTMPROPHET_RS_THREADS=N` env var (ultimate override for scripts / CI),
/// 2. an explicit `--threads N` from the CLI,
/// 3. all logical CPUs.
///
/// This governs only the internal data parallelism here; the output is identical for any value.
pub fn worker_count() -> usize {
    if let Ok(v) = std::env::var("PTMPROPHET_RS_THREADS") {
        if let Ok(n) = v.trim().parse::<usize>() {
            if n >= 1 {
                return n;
            }
        }
    }
    if let Some(&n) = THREAD_OVERRIDE.get() {
        return n;
    }
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

/// Target number of work chunks **per worker**. More chunks ⇒ finer dynamic balancing (a heavy item
/// strands at most one small chunk's worth of its neighbours before another core picks up the next
/// chunk), at the cost of slightly more queue pulls. Items handled here are coarse (a per-PSM match,
/// a grid-point density), so a fine split is essentially free relative to the work it balances.
const CHUNKS_PER_WORKER: usize = 32;

/// Contiguous block length so `len` items split into ≈ `workers · CHUNKS_PER_WORKER` chunks (≥ 1).
/// Bounds the total queue pulls to ≈ that chunk count regardless of `len`.
fn block_len(len: usize) -> usize {
    let targets = worker_count().saturating_mul(CHUNKS_PER_WORKER).max(1);
    (len / targets).max(1)
}

/// Number of pool threads to spawn: never more than there are chunks of work to pull.
fn pool_size(len: usize, block: usize) -> usize {
    worker_count().min(len.div_ceil(block)).max(1)
}

/// Apply `f` to every element of `items` in place, across a fixed worker pool that pulls disjoint
/// `&mut` sub-slices from a shared queue (so there is no aliasing). Order is irrelevant here.
pub fn par_for_each_mut<T, F>(items: &mut [T], f: F)
where
    T: Send,
    F: Fn(&mut T) + Sync,
{
    if items.is_empty() {
        return;
    }
    let block = block_len(items.len());
    let n_threads = pool_size(items.len(), block);
    let queue = Mutex::new(items.chunks_mut(block));
    let f = &f;
    std::thread::scope(|scope| {
        for _ in 0..n_threads {
            scope.spawn(|| loop {
                // Hold the lock only to pull the next chunk; process it unlocked.
                let chunk = { queue.lock().expect("work queue poisoned").next() };
                match chunk {
                    Some(ch) => {
                        for it in ch.iter_mut() {
                            f(it);
                        }
                    }
                    None => break,
                }
            });
        }
    });
}

/// Map `f` over `items` (mutable), collecting one result per element **in input order**. Workers pull
/// indexed chunks dynamically; each returns its `(chunk_index, results)` and we reassemble by index,
/// so the output is identical regardless of how the chunks were scheduled.
pub fn par_map_mut<T, R, F>(items: &mut [T], f: F) -> Vec<R>
where
    T: Send,
    R: Send,
    F: Fn(&mut T) -> R + Sync,
{
    if items.is_empty() {
        return Vec::new();
    }
    let block = block_len(items.len());
    let n_threads = pool_size(items.len(), block);
    let queue = Mutex::new(items.chunks_mut(block).enumerate());
    let parts: Mutex<Vec<(usize, Vec<R>)>> = Mutex::new(Vec::new());
    let f = &f;
    std::thread::scope(|scope| {
        for _ in 0..n_threads {
            scope.spawn(|| loop {
                let next = { queue.lock().expect("work queue poisoned").next() };
                match next {
                    Some((ci, ch)) => {
                        let mapped: Vec<R> = ch.iter_mut().map(f).collect();
                        parts.lock().expect("results queue poisoned").push((ci, mapped));
                    }
                    None => break,
                }
            });
        }
    });
    reassemble(parts.into_inner().expect("results queue poisoned"))
}

/// Map `f` over `items` (shared), collecting one result per element **in input order**.
pub fn par_map<T, R, F>(items: &[T], f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync,
{
    if items.is_empty() {
        return Vec::new();
    }
    let block = block_len(items.len());
    let n_threads = pool_size(items.len(), block);
    let queue = Mutex::new(items.chunks(block).enumerate());
    let parts: Mutex<Vec<(usize, Vec<R>)>> = Mutex::new(Vec::new());
    let f = &f;
    std::thread::scope(|scope| {
        for _ in 0..n_threads {
            scope.spawn(|| loop {
                let next = { queue.lock().expect("work queue poisoned").next() };
                match next {
                    Some((ci, ch)) => {
                        let mapped: Vec<R> = ch.iter().map(f).collect();
                        parts.lock().expect("results queue poisoned").push((ci, mapped));
                    }
                    None => break,
                }
            });
        }
    });
    reassemble(parts.into_inner().expect("results queue poisoned"))
}

/// Reassemble per-chunk `(chunk_index, results)` into one input-order `Vec`. Chunks tile `[0, len)`
/// contiguously in ascending index, so sorting by index restores input order deterministically.
fn reassemble<R>(mut parts: Vec<(usize, Vec<R>)>) -> Vec<R> {
    parts.sort_by_key(|(ci, _)| *ci);
    let total: usize = parts.iter().map(|(_, v)| v.len()).sum();
    let mut out = Vec::with_capacity(total);
    for (_, v) in parts {
        out.extend(v);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn par_map_preserves_order() {
        let v: Vec<i32> = (0..1000).collect();
        let got = par_map(&v, |&x| x * 2);
        let want: Vec<i32> = v.iter().map(|&x| x * 2).collect();
        assert_eq!(got, want);
    }

    #[test]
    fn par_map_mut_preserves_order_and_mutates() {
        let mut v: Vec<i32> = (0..1000).collect();
        let got = par_map_mut(&mut v, |x| {
            *x += 1;
            *x
        });
        assert_eq!(got, (1..1001).collect::<Vec<i32>>());
        assert_eq!(v, (1..1001).collect::<Vec<i32>>());
    }

    #[test]
    fn par_for_each_mut_handles_empty() {
        let mut v: Vec<i32> = Vec::new();
        par_for_each_mut(&mut v, |x| *x += 1);
        assert!(v.is_empty());
    }

    #[test]
    fn par_for_each_mut_visits_every_element_once() {
        // Uneven sizes around the chunk boundary, to exercise dynamic pulling + the tail chunk.
        for n in [1usize, 7, 12, 13, 100, 1000, 5000] {
            let mut v: Vec<usize> = vec![0; n];
            par_for_each_mut(&mut v, |x| *x += 1);
            assert!(v.iter().all(|&c| c == 1), "n={n}: every element incremented exactly once");
        }
    }

    #[test]
    fn par_map_order_independent_of_thread_count() {
        // Same input, two worker counts → identical output (determinism guarantee).
        let v: Vec<i32> = (0..3333).collect();
        std::env::set_var("PTMPROPHET_RS_THREADS", "1");
        let a = par_map(&v, |&x| x * x);
        std::env::set_var("PTMPROPHET_RS_THREADS", "8");
        let b = par_map(&v, |&x| x * x);
        std::env::remove_var("PTMPROPHET_RS_THREADS");
        assert_eq!(a, b);
        assert_eq!(a, v.iter().map(|&x| x * x).collect::<Vec<i32>>());
    }
}
