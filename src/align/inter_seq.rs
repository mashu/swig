//! Inter-sequence SIMD Gotoh: 16 independent alignments per AVX2 register.
//!
//! **Trick — one alignment per lane.** Striped (Farrar) SIMD puts 16 cells of
//! *one* alignment in a register and needs a lazy correction loop for the
//! horizontal gap state, which makes traceback awkward. Here each 16-bit lane
//! holds the same cell of a *different* job, so every lane runs the ordinary
//! scalar recurrence, traceback codes are the scalar ones, and the scalar
//! traceback walk ([`super::gotoh::traceback_with`]) is reused verbatim. The
//! price is that lanes need work: jobs are grouped by layout and shape and
//! packed 16 at a time; singletons go to the scalar kernel.
//!
//! **Band-relative coordinates.** A diagonal band is stored as rows `i` ×
//! offsets `k = (i − j) − lo`. Every lane walks the same `(row, k)` grid even
//! though bands start at different diagonals, the diagonal neighbour keeps its
//! `k`, the up neighbour is `k − 1` and the left neighbour `k + 1` (so `k` runs
//! downwards). The reference position of cell `(i, k)` depends only on `i − k`,
//! so each batch lays the lanes' profile columns out once, interleaved, and a
//! cell's substitution scores are one contiguous load per base code. Unbanded
//! jobs use a rectangular `(row, column)` layout instead.
//!
//! **Any substitution model.** Scores come from the per-allele
//! [`SubstitutionProfile`]s (5 rows: A C G T N), selected per lane with
//! query-code masks computed once per row — no assumption of match/mismatch
//! structure.
//!
//! **Trick — uniform fast path.** When every lane's profile depends on base
//! identity alone ([`SubstitutionProfile::uniform`], true for every
//! [`crate::scoring::AffineScoring`] profile), the group stores reference
//! *codes* instead of five interleaved profile rows, and a cell's score is
//! one load, one comparison and a blend instead of five masked loads. `N`
//! scores 0 and out-of-matrix cells score `−∞` exactly as in the profile
//! path, so results are identical.
//!
//! **Exact scalar semantics.** Out-of-band and out-of-matrix cells read as
//! `H = 0, I = D = −∞`, exactly like the scalar kernel's never-written
//! columns; lanes narrower than the batch are masked. Cells past a lane's
//! query or reference get a `−∞` substitution, so they can only be reached
//! through gaps from real cells, never feed back into them, and never become
//! a strict new maximum. The best cell is the first strict maximum in
//! row-major order, as in the scalar kernel. Scores, end cells and tracebacks
//! are therefore identical lane by lane (tested exhaustively against
//! [`BandedGotoh`]).
//!
//! **16-bit safety.** Lanes are 16-bit with saturating arithmetic and
//! `−∞ = −16384`. A job enters a lane only if its attainable score
//! (`max profile entry × min(n, m)`) and its costs stay far from that
//! sentinel; others fall back to the scalar kernel.

use super::gotoh::{BandedGotoh, GotohWorkspace, traceback_with};
use super::{AlignJob, Aligner, Alignment, Op, ScoreHit};
use crate::alphabet::N;
use crate::band::Band;
use crate::scoring::{GapCosts, Score, SubstitutionProfile};

/// Jobs per register (AVX2, 16-bit lanes).
pub const LANES: usize = 16;

/// Base codes with a profile row: A C G T N.
const CODES: usize = N as usize + 1;
/// Query code for rows past a lane's query.
const SENTINEL: i16 = N as i16 + 1;
/// Reference code for positions past a lane's reference (uniform layout).
const REF_SENTINEL: i16 = N as i16 + 2;
const NEG: i16 = -16_384;
/// Largest attainable score admitted into a lane.
const MAX_SCORE: i64 = 8_000;
/// Most negative single-step score or gap cost admitted into a lane.
const MIN_STEP: i64 = -1_000;
/// Largest row count / band width admitted into a lane (positions are i16).
const MAX_DIM: usize = 16_000;

/// Batched aligner: AVX2 inter-sequence kernel when available, scalar
/// [`BandedGotoh`] otherwise and for jobs that cannot share a register.
#[derive(Clone, Copy, Debug)]
pub struct InterSequenceGotoh {
    simd: bool,
}

impl InterSequenceGotoh {
    /// Detect AVX2 at run time.
    pub fn new() -> Self {
        Self {
            simd: avx2_available(),
        }
    }

    /// Never use SIMD (reference behaviour, non-x86 targets).
    pub fn scalar_only() -> Self {
        Self { simd: false }
    }

    pub fn uses_simd(&self) -> bool {
        self.simd
    }
}

impl Default for InterSequenceGotoh {
    fn default() -> Self {
        Self::new()
    }
}

fn avx2_available() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Geometry of one job inside a register.
#[derive(Clone, Copy, Debug)]
struct Lane {
    job: usize,
    /// Query rows before lane row 1 (diagonal layout; 0 for rectangular).
    shift: usize,
    /// Lowest diagonal `i − j` of the band (diagonal layout).
    lo: i32,
    /// Band width (diagonal) or reference length (rectangular).
    width: usize,
    /// Lane rows containing in-matrix cells.
    rows: usize,
    n: usize,
    m: usize,
}

impl Lane {
    /// Code of cell `(i, j)` (1-based) in a group's trace, `None` outside the
    /// lane's band or matrix — the contract of `traceback_with`.
    fn trace_code<const DIAG: bool>(
        &self,
        trace: &[[u8; LANES]],
        width: usize,
        lane: usize,
        i: usize,
        j: usize,
    ) -> Option<u8> {
        if i == 0 || j == 0 || i > self.n || j > self.m {
            return None;
        }
        let (r, k) = if DIAG {
            let d = i as i64 - j as i64 - i64::from(self.lo);
            if d < 0 || d >= self.width as i64 || i <= self.shift {
                return None;
            }
            (i - self.shift, d as usize)
        } else {
            (i, j - 1)
        };
        (r <= self.rows).then(|| trace[(r - 1) * width + k][lane])
    }

    /// 1-based `(row, col)` of the lane's best cell from kernel coordinates.
    fn cell_of<const DIAG: bool>(&self, r: usize, k: usize) -> (usize, usize) {
        if DIAG {
            let i = r + self.shift;
            (i, (i as i64 - i64::from(self.lo) - k as i64) as usize)
        } else {
            (r, k + 1)
        }
    }
}

/// Whether a job can share a register, and its geometry.
fn lane_of(job_index: usize, job: &AlignJob<'_>) -> Option<Lane> {
    let (n, m) = (job.query.len(), job.profile.len());
    if n == 0 || m == 0 {
        return None;
    }
    let bound = i64::from(job.profile.max_score().max(0)) * n.min(m) as i64;
    let floor = i64::from(
        job.profile
            .min_score()
            .min(job.gaps.open)
            .min(job.gaps.extend),
    );
    if bound > MAX_SCORE || floor < MIN_STEP {
        return None;
    }
    match job.band {
        Band::Full => (n <= MAX_DIM && m <= MAX_DIM).then_some(Lane {
            job: job_index,
            shift: 0,
            lo: 0,
            width: m,
            rows: n,
            n,
            m,
        }),
        Band::Diagonals { lo, hi } => {
            if hi < lo {
                return None;
            }
            let width = (i64::from(hi) - i64::from(lo) + 1) as usize;
            let shift = lo.max(0) as usize;
            let last = (n as i64).min(m as i64 + i64::from(hi));
            let rows = usize::try_from(last - shift as i64)
                .ok()
                .filter(|&rows| rows > 0)?;
            (width <= MAX_DIM && rows <= MAX_DIM && lo.unsigned_abs() as usize <= MAX_DIM)
                .then_some(Lane {
                    job: job_index,
                    shift,
                    lo,
                    width,
                    rows,
                    n,
                    m,
                })
        }
    }
}

/// Per-group constants handed to the kernel.
struct Shape {
    rows: usize,
    width: usize,
    positions: usize,
    open: [i16; LANES],
    extend: [i16; LANES],
    step: [i16; LANES],
    lane_width: [i16; LANES],
    /// Uniform layout only: per-lane match and mismatch scores.
    matched: [i16; LANES],
    mismatch: [i16; LANES],
}

/// Best cell of each lane in kernel coordinates.
struct LaneBest {
    score: [i16; LANES],
    row: [i16; LANES],
    k: [i16; LANES],
}

/// Reusable buffers for [`InterSequenceGotoh`].
#[derive(Debug, Default)]
pub struct InterSequenceWorkspace {
    scalar: GotohWorkspace,
    scalar_jobs: Vec<usize>,
    diagonal: Vec<Lane>,
    rectangular: Vec<Lane>,
    query: Vec<[i16; LANES]>,
    profile: Vec<[i16; LANES]>,
    h: Vec<[i16; LANES]>,
    ins: Vec<[i16; LANES]>,
    trace: Vec<[u8; LANES]>,
}

/// Best cell and (for traceback passes) the alignment start and columns.
struct LaneResult {
    score: Score,
    row: usize,
    col: usize,
    start: Option<(usize, usize, Vec<Op>)>,
}

/// Lay out the group's queries and, per position, either the five
/// interleaved profile rows or (`UNIFORM`) the reference codes; return its shape.
fn build<const DIAG: bool, const UNIFORM: bool>(
    jobs: &[AlignJob<'_>],
    lanes: &[Lane],
    query: &mut Vec<[i16; LANES]>,
    profile: &mut Vec<[i16; LANES]>,
) -> Shape {
    let rows = lanes.iter().map(|l| l.rows).max().unwrap_or(0);
    let width = lanes.iter().map(|l| l.width).max().unwrap_or(0);
    let positions = if DIAG { rows + width } else { width };
    query.clear();
    query.resize(rows, [SENTINEL; LANES]);
    profile.clear();
    if UNIFORM {
        profile.resize(positions, [REF_SENTINEL; LANES]);
    } else {
        profile.resize(CODES * positions, [NEG; LANES]);
    }
    let mut shape = Shape {
        rows,
        width,
        positions,
        open: [0; LANES],
        extend: [0; LANES],
        step: [0; LANES],
        lane_width: [0; LANES],
        matched: [0; LANES],
        mismatch: [0; LANES],
    };
    for (l, lane) in lanes.iter().enumerate() {
        let job = &jobs[lane.job];
        for r in 1..=lane.rows {
            query[r - 1][l] = i16::from(job.query[r + lane.shift - 1].min(N));
        }
        // Diagonal layout: position t holds reference index j0 = t − width + shift − lo.
        let offset = if DIAG {
            width as i64 - lane.shift as i64 + i64::from(lane.lo)
        } else {
            0
        };
        let positions = 0..positions as i64;
        if UNIFORM {
            let scores = job.profile.uniform().expect("uniform group");
            shape.matched[l] = scores.matched as i16;
            shape.mismatch[l] = scores.mismatch as i16;
            for (j0, &code) in job.profile.codes().iter().enumerate() {
                let t = j0 as i64 + offset;
                if positions.contains(&t) {
                    profile[t as usize][l] = i16::from(code);
                }
            }
        } else {
            for code in 0..CODES {
                let row = job.profile.row(code as u8);
                let column = &mut profile[code * shape.positions..(code + 1) * shape.positions];
                for (j0, &score) in row.iter().enumerate() {
                    let t = j0 as i64 + offset;
                    if positions.contains(&t) {
                        column[t as usize][l] = score as i16;
                    }
                }
            }
        }
        let GapCosts { open, extend } = job.gaps;
        shape.open[l] = open as i16;
        shape.extend[l] = extend as i16;
        shape.step[l] = open.max(extend) as i16;
        shape.lane_width[l] = lane.width as i16;
    }
    shape
}

#[cfg(target_arch = "x86_64")]
mod kernel {
    use super::{CODES, LANES, LaneBest, NEG, SENTINEL, Shape};
    use crate::alphabet::N;
    use std::arch::x86_64::*;

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn load(v: &[i16; LANES]) -> __m256i {
        unsafe { _mm256_loadu_si256(v.as_ptr().cast()) }
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn store(v: &mut [i16; LANES], x: __m256i) {
        unsafe { _mm256_storeu_si256(v.as_mut_ptr().cast(), x) }
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn store_codes(v: &mut [u8; LANES], x: __m256i) {
        unsafe {
            let packed =
                _mm_packus_epi16(_mm256_castsi256_si128(x), _mm256_extracti128_si256::<1>(x));
            _mm_storeu_si128(v.as_mut_ptr().cast(), packed)
        }
    }

    /// Fill one group. `h` and `ins` hold `width + 1` entries; entry 0 is the
    /// virtual out-of-band neighbour (`H = 0`, `I = −∞`) and is never written.
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn fill<const DIAG: bool, const TRACE: bool, const UNIFORM: bool>(
        shape: &Shape,
        query: &[[i16; LANES]],
        profile: &[[i16; LANES]],
        h: &mut Vec<[i16; LANES]>,
        ins: &mut Vec<[i16; LANES]>,
        trace: &mut Vec<[u8; LANES]>,
    ) -> LaneBest {
        unsafe {
            let (rows, width, positions) = (shape.rows, shape.width, shape.positions);
            h.clear();
            h.resize(width + 1, [0; LANES]);
            ins.clear();
            ins.resize(width + 1, [NEG; LANES]);
            if TRACE {
                // Every cell is written before traceback reads it: grow only.
                if trace.len() < rows * width {
                    trace.resize(rows * width, [0; LANES]);
                }
            }
            let zero = _mm256_setzero_si256();
            let neg = _mm256_set1_epi16(NEG);
            let (open, extend, step) = (load(&shape.open), load(&shape.extend), load(&shape.step));
            let lane_width = load(&shape.lane_width);
            let (one, two, three) = (
                _mm256_set1_epi16(1),
                _mm256_set1_epi16(2),
                _mm256_set1_epi16(3),
            );
            let (ins_bit, del_bit) = (_mm256_set1_epi16(4), _mm256_set1_epi16(8));
            let sentinel = _mm256_set1_epi16(SENTINEL);
            let ambiguous_code = _mm256_set1_epi16(i16::from(N));
            let (matched, mismatch) = (load(&shape.matched), load(&shape.mismatch));
            let (mut best, mut best_row, mut best_k) = (zero, zero, zero);

            for r in 1..=rows {
                let q = load(&query[r - 1]);
                let mut code_mask = [zero; CODES];
                if !UNIFORM {
                    for (c, mask) in code_mask.iter_mut().enumerate() {
                        *mask = _mm256_cmpeq_epi16(q, _mm256_set1_epi16(c as i16));
                    }
                }
                let past_query_lanes = _mm256_cmpeq_epi16(q, sentinel);
                let past_query = _mm256_and_si256(past_query_lanes, neg);
                let query_ambiguous = _mm256_cmpeq_epi16(q, ambiguous_code);
                let row_vector = _mm256_set1_epi16(r as i16);
                let mut del = neg;
                let mut del_extends = zero;
                let mut h_diag_carry = zero; // rectangular layout: H(i−1, j−1)

                for step_index in 0..width {
                    // Diagonal layout walks k downwards (left neighbour is k + 1).
                    let k = if DIAG {
                        width - 1 - step_index
                    } else {
                        step_index
                    };
                    let (h_up, i_up, h_diag, t) = if DIAG {
                        (
                            load(&h[k]),
                            load(&ins[k]),
                            load(&h[k + 1]),
                            r + width - 1 - k,
                        )
                    } else {
                        let up = load(&h[k + 1]);
                        let diag = h_diag_carry;
                        h_diag_carry = up;
                        (up, load(&ins[k + 1]), diag, k)
                    };
                    let s = if UNIFORM {
                        // matched / mismatch by identity, 0 for N, −∞ outside the matrix.
                        let reference = load(&profile[t]);
                        let identical = _mm256_cmpeq_epi16(q, reference);
                        let base = _mm256_blendv_epi8(mismatch, matched, identical);
                        let ambiguous = _mm256_or_si256(
                            query_ambiguous,
                            _mm256_cmpeq_epi16(reference, ambiguous_code),
                        );
                        let base = _mm256_andnot_si256(ambiguous, base);
                        let outside = _mm256_or_si256(
                            past_query_lanes,
                            _mm256_cmpgt_epi16(reference, ambiguous_code),
                        );
                        _mm256_blendv_epi8(base, neg, outside)
                    } else {
                        let mut s = past_query;
                        for (c, mask) in code_mask.iter().enumerate() {
                            s = _mm256_or_si256(
                                s,
                                _mm256_and_si256(load(&profile[c * positions + t]), *mask),
                            );
                        }
                        s
                    };

                    let ins_open = _mm256_adds_epi16(h_up, open);
                    let ins_extend = _mm256_adds_epi16(i_up, extend);
                    let mut ins_score = _mm256_max_epi16(ins_open, ins_extend);
                    let diag = _mm256_adds_epi16(h_diag, s);
                    let t0 = _mm256_max_epi16(diag, zero);
                    let t_all = _mm256_max_epi16(t0, ins_score);
                    let mut score = _mm256_max_epi16(t_all, del);

                    if TRACE {
                        let source = _mm256_and_si256(_mm256_cmpgt_epi16(diag, zero), one);
                        let source =
                            _mm256_blendv_epi8(source, two, _mm256_cmpgt_epi16(ins_score, t0));
                        let source =
                            _mm256_blendv_epi8(source, three, _mm256_cmpgt_epi16(del, t_all));
                        let bits = _mm256_or_si256(
                            _mm256_and_si256(_mm256_cmpgt_epi16(ins_extend, ins_open), ins_bit),
                            _mm256_and_si256(del_extends, del_bit),
                        );
                        store_codes(
                            &mut trace[(r - 1) * width + k],
                            _mm256_or_si256(source, bits),
                        );
                    }

                    let del_open = _mm256_adds_epi16(t_all, open);
                    let mut next_extends =
                        _mm256_cmpgt_epi16(_mm256_adds_epi16(del, extend), del_open);
                    let mut next_del = _mm256_max_epi16(del_open, _mm256_adds_epi16(del, step));
                    let k_vector = _mm256_set1_epi16(k as i16);
                    if DIAG {
                        // Lanes whose band is narrower than the group's: k ≥ width_l does not exist.
                        let in_band = _mm256_cmpgt_epi16(lane_width, k_vector);
                        score = _mm256_and_si256(score, in_band);
                        ins_score = _mm256_blendv_epi8(neg, ins_score, in_band);
                        next_del = _mm256_blendv_epi8(neg, next_del, in_band);
                        next_extends = _mm256_and_si256(next_extends, in_band);
                    }
                    store(&mut h[k + 1], score);
                    store(&mut ins[k + 1], ins_score);
                    del = next_del;
                    del_extends = next_extends;

                    let improved = _mm256_cmpgt_epi16(score, best);
                    best = _mm256_max_epi16(best, score);
                    best_row = _mm256_blendv_epi8(best_row, row_vector, improved);
                    best_k = _mm256_blendv_epi8(best_k, k_vector, improved);
                }
            }
            let mut out = LaneBest {
                score: [0; LANES],
                row: [0; LANES],
                k: [0; LANES],
            };
            store(&mut out.score, best);
            store(&mut out.row, best_row);
            store(&mut out.k, best_k);
            out
        }
    }
}

impl InterSequenceGotoh {
    /// Run every SIMD-eligible group; route the rest to `ws.scalar_jobs`.
    fn simd_groups<const TRACE: bool>(
        &self,
        jobs: &[AlignJob<'_>],
        ws: &mut InterSequenceWorkspace,
        mut emit: impl FnMut(usize, Option<LaneResult>),
    ) {
        let InterSequenceWorkspace {
            scalar_jobs,
            diagonal,
            rectangular,
            query,
            profile,
            h,
            ins,
            trace,
            ..
        } = ws;
        scalar_jobs.clear();
        diagonal.clear();
        rectangular.clear();
        for (index, job) in jobs.iter().enumerate() {
            match lane_of(index, job).filter(|_| self.simd) {
                Some(lane) if matches!(job.band, Band::Full) => rectangular.push(lane),
                Some(lane) => diagonal.push(lane),
                None => scalar_jobs.push(index),
            }
        }
        // Similar shapes share registers: less padding per group.
        diagonal.sort_by_key(|l| std::cmp::Reverse((l.rows, l.width)));
        rectangular.sort_by_key(|l| std::cmp::Reverse((l.rows, l.width)));

        #[cfg(target_arch = "x86_64")]
        {
            #[allow(clippy::too_many_arguments)]
            fn run<const DIAG: bool, const TRACE: bool>(
                jobs: &[AlignJob<'_>],
                lanes: &[Lane],
                query: &mut Vec<[i16; LANES]>,
                profile: &mut Vec<[i16; LANES]>,
                h: &mut Vec<[i16; LANES]>,
                ins: &mut Vec<[i16; LANES]>,
                trace: &mut Vec<[u8; LANES]>,
                emit: &mut impl FnMut(usize, Option<LaneResult>),
            ) {
                let uniform = lanes
                    .iter()
                    .all(|l| jobs[l.job].profile.uniform().is_some());
                // SAFETY: only reached when AVX2 was detected at construction.
                let (shape, best) = if uniform {
                    let shape = build::<DIAG, true>(jobs, lanes, query, profile);
                    let best = unsafe {
                        kernel::fill::<DIAG, TRACE, true>(&shape, query, profile, h, ins, trace)
                    };
                    (shape, best)
                } else {
                    let shape = build::<DIAG, false>(jobs, lanes, query, profile);
                    let best = unsafe {
                        kernel::fill::<DIAG, TRACE, false>(&shape, query, profile, h, ins, trace)
                    };
                    (shape, best)
                };
                for (l, lane) in lanes.iter().enumerate() {
                    let score = Score::from(best.score[l]);
                    if score <= 0 {
                        emit(lane.job, None);
                        continue;
                    }
                    let (row, col) = lane.cell_of::<DIAG>(best.row[l] as usize, best.k[l] as usize);
                    let start = TRACE.then(|| {
                        traceback_with(
                            |i, j| lane.trace_code::<DIAG>(trace, shape.width, l, i, j),
                            row,
                            col,
                        )
                    });
                    emit(
                        lane.job,
                        Some(LaneResult {
                            score,
                            row,
                            col,
                            start,
                        }),
                    );
                }
            }
            for chunk in diagonal.chunks(LANES) {
                if chunk.len() < 2 {
                    scalar_jobs.extend(chunk.iter().map(|l| l.job));
                } else {
                    run::<true, TRACE>(jobs, chunk, query, profile, h, ins, trace, &mut emit);
                }
            }
            for chunk in rectangular.chunks(LANES) {
                if chunk.len() < 2 {
                    scalar_jobs.extend(chunk.iter().map(|l| l.job));
                } else {
                    run::<false, TRACE>(jobs, chunk, query, profile, h, ins, trace, &mut emit);
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (&mut emit, query, profile, h, ins, trace);
            scalar_jobs.extend(diagonal.iter().chain(rectangular.iter()).map(|l| l.job));
        }
    }
}

impl Aligner for InterSequenceGotoh {
    type Workspace = InterSequenceWorkspace;

    fn score(
        &self,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
        band: Band,
        workspace: &mut InterSequenceWorkspace,
    ) -> Option<ScoreHit> {
        BandedGotoh.score(query, profile, gaps, band, &mut workspace.scalar)
    }

    fn align(
        &self,
        query: &[u8],
        profile: &SubstitutionProfile,
        gaps: GapCosts,
        band: Band,
        workspace: &mut InterSequenceWorkspace,
    ) -> Option<Alignment> {
        BandedGotoh.align(query, profile, gaps, band, &mut workspace.scalar)
    }

    fn score_batch(
        &self,
        jobs: &[AlignJob<'_>],
        ws: &mut InterSequenceWorkspace,
        out: &mut Vec<Option<ScoreHit>>,
    ) {
        out.clear();
        out.resize(jobs.len(), None);
        self.simd_groups::<false>(jobs, ws, |job, result| {
            out[job] = result.map(|r| ScoreHit {
                score: r.score,
                q_end: r.row as u32,
                r_end: r.col as u32,
            });
        });
        for &index in &ws.scalar_jobs {
            let job = &jobs[index];
            out[index] =
                BandedGotoh.score(job.query, job.profile, job.gaps, job.band, &mut ws.scalar);
        }
    }

    fn align_batch(
        &self,
        jobs: &[AlignJob<'_>],
        ws: &mut InterSequenceWorkspace,
        out: &mut Vec<Option<Alignment>>,
    ) {
        out.clear();
        out.resize(jobs.len(), None);
        self.simd_groups::<true>(jobs, ws, |job, result| {
            out[job] = result.map(|r| {
                let (q_start, r_start, ops) = r.start.expect("traceback pass records the start");
                Alignment {
                    q_start: q_start as u32,
                    q_end: r.row as u32,
                    r_start: r_start as u32,
                    r_end: r.col as u32,
                    score: r.score,
                    ops,
                }
            });
        });
        for &index in &ws.scalar_jobs {
            let job = &jobs[index];
            out[index] =
                BandedGotoh.align(job.query, job.profile, job.gaps, job.band, &mut ws.scalar);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{AffineScoring, HotspotScoring};

    use crate::testing::rng;

    /// A reference plus a mutated, indel-carrying, N-sprinkled query that
    /// contains it at `offset`.
    fn case(next: &mut impl FnMut() -> u64) -> (Vec<u8>, Vec<u8>, i32) {
        let m = 8 + (next() % 120) as usize;
        let reference: Vec<u8> = (0..m).map(|_| (next() % 4) as u8).collect();
        let mut core = reference.clone();
        for _ in 0..(next() % 6) {
            let p = (next() as usize) % core.len();
            core[p] = (next() % 5) as u8; // includes N
        }
        for _ in 0..(next() % 3) {
            let p = (next() as usize) % core.len();
            core.remove(p);
        }
        for _ in 0..(next() % 3) {
            let p = (next() as usize) % core.len();
            core.insert(p, (next() % 4) as u8);
        }
        let lead = (next() % 30) as usize;
        let mut query: Vec<u8> = (0..lead).map(|_| (next() % 4) as u8).collect();
        query.extend(core);
        query.extend((0..(next() % 30)).map(|_| (next() % 4) as u8));
        (query, reference, lead as i32)
    }

    #[test]
    fn identical_to_scalar_for_mixed_batches() {
        let simd = InterSequenceGotoh::new();
        if !simd.uses_simd() {
            eprintln!("AVX2 not available: SIMD path not exercised");
        }
        let scorings = [
            AffineScoring::new(2, -3, -5, -1),
            AffineScoring::new(2, -3, -13, -1),
            AffineScoring::new(2, -3, -17, -2),
            AffineScoring::new(1, -1, -2, -2),
        ];
        let mut next = rng(42);
        let mut ws = InterSequenceWorkspace::default();
        let mut scalar_ws = GotohWorkspace::default();
        let (mut hits, mut alignments) = (Vec::new(), Vec::new());
        for batch in 0..160 {
            let size = 1 + (next() % 40) as usize;
            let cases: Vec<_> = (0..size).map(|_| case(&mut next)).collect();
            // Every other batch mixes in position-specific (hotspot) profiles, so
            // both the uniform fast path and the profile path are exercised.
            let hotspots = batch % 2 == 1;
            let profiles: Vec<_> = cases
                .iter()
                .enumerate()
                .map(|(k, (_, r, _))| {
                    let scoring = scorings[k % 4];
                    let profile = if hotspots && k % 3 == 0 {
                        SubstitutionProfile::build(&HotspotScoring::for_shm(scoring), r)
                    } else {
                        SubstitutionProfile::build(&scoring, r)
                    };
                    (profile, scoring.gaps)
                })
                .collect();
            let jobs: Vec<AlignJob> = cases
                .iter()
                .zip(&profiles)
                .map(|((q, _, lead), (profile, gaps))| {
                    let band = match next() % 4 {
                        0 => Band::Full,
                        1 => Band::around(*lead, 3 + (next() % 30) as i32),
                        2 => Band::around(*lead + (next() % 9) as i32 - 4, 28),
                        _ => Band::Diagonals {
                            lo: *lead - 5,
                            hi: *lead + 60,
                        },
                    };
                    AlignJob {
                        query: q,
                        profile,
                        gaps: *gaps,
                        band,
                    }
                })
                .collect();
            simd.score_batch(&jobs, &mut ws, &mut hits);
            simd.align_batch(&jobs, &mut ws, &mut alignments);
            for (k, job) in jobs.iter().enumerate() {
                let expected_hit =
                    BandedGotoh.score(job.query, job.profile, job.gaps, job.band, &mut scalar_ws);
                let expected_aln =
                    BandedGotoh.align(job.query, job.profile, job.gaps, job.band, &mut scalar_ws);
                assert_eq!(
                    hits[k], expected_hit,
                    "score, batch {batch} job {k} band {:?}",
                    job.band
                );
                assert_eq!(
                    alignments[k], expected_aln,
                    "alignment, batch {batch} job {k} band {:?}",
                    job.band
                );
            }
        }
    }

    /// Degenerate shapes: empty sides, all-N queries, one-row queries,
    /// one-diagonal bands and bands that miss the matrix entirely.
    #[test]
    fn degenerate_jobs_match_the_scalar_kernel() {
        let scoring = AffineScoring::new(2, -3, -5, -1);
        let reference: Vec<u8> = (0..40).map(|i| (i * 7 % 4) as u8).collect();
        let profiles = [
            SubstitutionProfile::build(&scoring, &reference),
            SubstitutionProfile::build(&scoring, &reference[..1]),
            SubstitutionProfile::build(&scoring, &[]),
            SubstitutionProfile::build(&HotspotScoring::for_shm(scoring), &reference),
        ];
        let queries: [Vec<u8>; 5] = [
            reference.clone(),
            Vec::new(),
            vec![N; 30],
            reference[5..6].to_vec(),
            [vec![N; 10], reference[..20].to_vec()].concat(),
        ];
        let bands = [
            Band::Full,
            Band::Diagonals { lo: 0, hi: 0 },
            Band::Diagonals { lo: 10, hi: 10 },
            Band::Diagonals { lo: 500, hi: 520 },
            Band::Diagonals { lo: -520, hi: -500 },
            Band::Diagonals { lo: -3, hi: 3 },
        ];
        let mut jobs = Vec::new();
        for query in &queries {
            for profile in &profiles {
                for &band in &bands {
                    jobs.push(AlignJob {
                        query,
                        profile,
                        gaps: scoring.gaps,
                        band,
                    });
                }
            }
        }
        let (mut hits, mut alignments) = (Vec::new(), Vec::new());
        let mut ws = InterSequenceWorkspace::default();
        InterSequenceGotoh::new().score_batch(&jobs, &mut ws, &mut hits);
        InterSequenceGotoh::new().align_batch(&jobs, &mut ws, &mut alignments);
        let mut scalar = GotohWorkspace::default();
        for (k, job) in jobs.iter().enumerate() {
            let (query, profile, gaps, band) = (job.query, job.profile, job.gaps, job.band);
            assert_eq!(
                hits[k],
                BandedGotoh.score(query, profile, gaps, band, &mut scalar),
                "job {k}"
            );
            assert_eq!(
                alignments[k],
                BandedGotoh.align(query, profile, gaps, band, &mut scalar),
                "job {k}"
            );
        }
    }

    #[test]
    fn oversized_scores_fall_back_to_scalar() {
        let scoring = AffineScoring::new(100, -3, -5, -1);
        let reference: Vec<u8> = (0..200).map(|i| (i % 4) as u8).collect();
        let profile = SubstitutionProfile::build(&scoring, &reference);
        assert!(
            lane_of(
                0,
                &AlignJob {
                    query: &reference,
                    profile: &profile,
                    gaps: scoring.gaps,
                    band: Band::Full
                }
            )
            .is_none()
        );
        let jobs = vec![
            AlignJob {
                query: &reference,
                profile: &profile,
                gaps: scoring.gaps,
                band: Band::Full
            };
            3
        ];
        let mut out = Vec::new();
        InterSequenceGotoh::new().score_batch(
            &jobs,
            &mut InterSequenceWorkspace::default(),
            &mut out,
        );
        assert!(out.iter().all(|hit| hit.map(|h| h.score) == Some(20_000)));
    }
}
