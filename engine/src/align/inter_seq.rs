//! Inter-sequence SIMD Gotoh: many independent alignments per register.
//!
//! **Trick — one alignment per lane.** Striped (Farrar) SIMD puts several
//! cells of *one* alignment in a register and needs a lazy correction loop
//! for the horizontal gap state, which makes traceback awkward. Here each
//! 16-bit lane holds the same cell of a *different* job, so every lane runs
//! the ordinary scalar recurrence, traceback codes are the scalar ones, and
//! the scalar traceback walk ([`super::gotoh::traceback_with`]) is reused
//! verbatim. The price is that lanes need work: jobs are grouped by layout
//! and shape and packed a register at a time; singletons go to the scalar
//! kernel.
//!
//! **Any CPU.** The kernel ([`fill`]) is written once over the [`Lanes`]
//! abstraction and compiled per [`Backend`] — AVX-512BW (32 lanes), AVX2,
//! NEON, WebAssembly `simd128` or portable Rust (16 lanes each) — chosen at
//! run time (at build time for WebAssembly).
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
//! are therefore identical lane by lane on every backend (tested
//! exhaustively against [`BandedGotoh`]).
//!
//! **16-bit safety.** Lanes are 16-bit with saturating arithmetic and
//! `−∞ = −16384`. A job enters a lane only if its attainable score
//! ([`SubstitutionProfile::local_bound`] over `min(n, m)` columns) and its
//! costs stay far from that sentinel; others fall back to the scalar kernel.

use super::gotoh::{BandedGotoh, GotohWorkspace, traceback_with};
#[cfg(target_arch = "aarch64")]
use super::lanes::Neon;
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
use super::lanes::Simd128;
#[cfg(target_arch = "x86_64")]
use super::lanes::{Avx2, Avx512};
use super::lanes::{Kernel, Lanes, Portable};
use super::{AlignJob, Aligner, Alignment, ScoreHit};
use crate::alphabet::N;
use crate::band::Band;
use crate::scoring::{GapCosts, Score, SubstitutionProfile};
use std::ops::Range;

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

/// Instruction set of the inter-sequence kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Backend {
    /// x86-64 AVX-512BW: 32 lanes.
    Avx512,
    /// x86-64 AVX2: 16 lanes.
    Avx2,
    /// AArch64 NEON: 16 lanes in two registers.
    Neon,
    /// WebAssembly `simd128`: 16 lanes in two registers. Present only in
    /// builds with the `simd128` target feature.
    Simd128,
    /// Plain Rust lanes on any target, vectorised by the compiler for the
    /// target's baseline SIMD: 16 lanes.
    Portable,
    /// No batching: every job runs on the scalar kernel.
    Scalar,
}

impl Backend {
    pub const ALL: [Backend; 6] = [
        Self::Avx512,
        Self::Avx2,
        Self::Neon,
        Self::Simd128,
        Self::Portable,
        Self::Scalar,
    ];

    /// Whether this CPU can run the backend (for `Simd128`, whether this
    /// build has it).
    pub fn available(self) -> bool {
        match self {
            #[cfg(target_arch = "x86_64")]
            Self::Avx512 => std::arch::is_x86_feature_detected!("avx512bw"),
            #[cfg(target_arch = "x86_64")]
            Self::Avx2 => std::arch::is_x86_feature_detected!("avx2"),
            #[cfg(target_arch = "aarch64")]
            Self::Neon => std::arch::is_aarch64_feature_detected!("neon"),
            #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
            Self::Simd128 => true,
            Self::Portable | Self::Scalar => true,
            _ => false,
        }
    }

    /// The fastest backend this CPU runs. AVX-512 is used only when the
    /// `avx512` preference is given ([`InterSequenceGotoh::with_backend`]):
    /// on many CPUs 512-bit instructions lower the clock of every core.
    /// Without an explicit backend, portable lanes, which the compiler
    /// vectorises where the target's baseline has SIMD. On WebAssembly
    /// without `simd128`, where 16 emulated lanes were measured several
    /// times slower than the scalar kernel, the scalar kernel is chosen
    /// instead; other targets keep portable lanes.
    pub fn detect() -> Self {
        let emulated = cfg!(all(target_family = "wasm", not(target_feature = "simd128")));
        [Self::Avx2, Self::Neon, Self::Simd128]
            .into_iter()
            .find(|b| b.available())
            .unwrap_or(if emulated {
                Self::Scalar
            } else {
                Self::Portable
            })
    }

    /// Jobs per register.
    pub fn lanes(self) -> usize {
        match self {
            Self::Avx512 => 32,
            Self::Avx2 | Self::Neon | Self::Simd128 | Self::Portable => 16,
            Self::Scalar => 1,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Avx512 => "avx512",
            Self::Avx2 => "avx2",
            Self::Neon => "neon",
            Self::Simd128 => "simd128",
            Self::Portable => "portable",
            Self::Scalar => "scalar",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|b| b.name() == name)
    }
}

/// Batched aligner: the inter-sequence kernel on the best available
/// [`Backend`], and the scalar [`BandedGotoh`] for jobs that cannot share a
/// register.
#[derive(Clone, Copy, Debug)]
pub struct InterSequenceGotoh {
    backend: Backend,
}

impl InterSequenceGotoh {
    /// The fastest backend of this CPU ([`Backend::detect`]).
    pub fn new() -> Self {
        Self {
            backend: Backend::detect(),
        }
    }

    /// Never batch (reference behaviour).
    pub fn scalar_only() -> Self {
        Self {
            backend: Backend::Scalar,
        }
    }

    /// A specific backend, if this CPU runs it.
    pub fn with_backend(backend: Backend) -> Option<Self> {
        backend.available().then_some(Self { backend })
    }

    pub fn backend(&self) -> Backend {
        self.backend
    }
}

impl Default for InterSequenceGotoh {
    fn default() -> Self {
        Self::new()
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
    fn trace_code<const DIAG: bool, const L: usize>(
        &self,
        trace: &[[u8; L]],
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
    let bound = job.profile.local_bound(n.min(m));
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
pub(super) struct Shape<const L: usize> {
    rows: usize,
    width: usize,
    positions: usize,
    open: [i16; L],
    extend: [i16; L],
    step: [i16; L],
    lane_width: [i16; L],
    /// Uniform layout only: per-lane match and mismatch scores.
    matched: [i16; L],
    mismatch: [i16; L],
}

/// Best cell of each lane in kernel coordinates.
pub(super) struct LaneBest<const L: usize> {
    score: [i16; L],
    row: [i16; L],
    k: [i16; L],
}

/// Reusable buffers for [`InterSequenceGotoh`]. Lane buffers are flat and
/// viewed as `[i16; L]` rows of the backend's width.
#[derive(Debug, Default)]
pub struct InterSequenceWorkspace {
    scalar: GotohWorkspace,
    scalar_jobs: Vec<usize>,
    diagonal: Vec<Lane>,
    rectangular: Vec<Lane>,
    query: Vec<i16>,
    profile: Vec<i16>,
    h: Vec<i16>,
    ins: Vec<i16>,
    trace: Vec<u8>,
}

/// `buffer` holding `rows` lane rows filled with `value`, as lane rows.
fn lane_rows<T: Copy, const L: usize>(buffer: &mut Vec<T>, rows: usize, value: T) -> &mut [[T; L]] {
    buffer.clear();
    buffer.resize(rows * L, value);
    buffer.as_chunks_mut::<L>().0
}

/// A lane's best local score and, on a traceback pass, its alignment.
struct LaneResult {
    score: Score,
    alignment: Option<Alignment>,
}

/// Reference indices `j0` whose layout position `t = j0 + offset` lies in
/// `[0, positions)`, and the position of the first.
fn reachable(len: usize, offset: i64, positions: usize) -> (Range<usize>, usize) {
    let lo = (-offset).clamp(0, len as i64);
    let hi = (positions as i64 - offset).clamp(lo, len as i64);
    (lo as usize..hi as usize, (lo + offset).max(0) as usize)
}

/// Lay out the group's queries and, per position, either the five
/// interleaved profile rows or (`UNIFORM`) the reference codes; return its shape.
fn build<const DIAG: bool, const UNIFORM: bool, const L: usize>(
    jobs: &[AlignJob<'_>],
    lanes: &[Lane],
    query: &mut Vec<i16>,
    profile: &mut Vec<i16>,
) -> Shape<L> {
    let rows = lanes.iter().map(|l| l.rows).max().unwrap_or(0);
    let width = lanes.iter().map(|l| l.width).max().unwrap_or(0);
    let positions = if DIAG { rows + width } else { width };
    let query = lane_rows::<i16, L>(query, rows, SENTINEL);
    let profile = if UNIFORM {
        lane_rows::<i16, L>(profile, positions, REF_SENTINEL)
    } else {
        lane_rows::<i16, L>(profile, CODES * positions, NEG)
    };
    let mut shape = Shape {
        rows,
        width,
        positions,
        open: [0; L],
        extend: [0; L],
        step: [0; L],
        lane_width: [0; L],
        matched: [0; L],
        mismatch: [0; L],
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
        // Only the reference indices that land inside the layout: a long
        // reference under a narrow band touches `positions` columns, not `m`.
        let (reach, first) = reachable(job.profile.len(), offset, positions);
        if UNIFORM {
            let scores = job.profile.uniform().expect("uniform group");
            shape.matched[l] = scores.matched as i16;
            shape.mismatch[l] = scores.mismatch as i16;
            for (t, &code) in (first..).zip(&job.profile.codes()[reach]) {
                profile[t][l] = i16::from(code);
            }
        } else {
            for code in 0..CODES {
                let row = job.profile.row(code as u8);
                let column = &mut profile[code * shape.positions..(code + 1) * shape.positions];
                for (t, &score) in (first..).zip(&row[reach.clone()]) {
                    column[t][l] = score as i16;
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

/// Fill one group of `L` lanes. `h` and `ins` hold `width + 1` rows; row 0
/// is the virtual out-of-band neighbour (`H = 0`, `I = −∞`) and is never
/// written. `trace` holds at least `rows · width` rows when `TRACE`.
///
/// # Safety
///
/// The CPU must support `V`'s instruction set; [`Kernel`] entry points
/// compile this under the matching target features.
#[inline(always)]
pub(super) unsafe fn fill<
    V,
    const L: usize,
    const DIAG: bool,
    const TRACE: bool,
    const UNIFORM: bool,
>(
    shape: &Shape<L>,
    query: &[[i16; L]],
    profile: &[[i16; L]],
    h: &mut [[i16; L]],
    ins: &mut [[i16; L]],
    trace: &mut [[u8; L]],
) -> LaneBest<L>
where
    V: Lanes<L>,
{
    // SAFETY: every `Lanes` operation needs only the backend's features.
    unsafe {
        let (rows, width, positions) = (shape.rows, shape.width, shape.positions);
        let zero = V::splat(0);
        let neg = V::splat(NEG);
        let (open, extend, step) = (
            V::load(&shape.open),
            V::load(&shape.extend),
            V::load(&shape.step),
        );
        let lane_width = V::load(&shape.lane_width);
        let (one, two, three) = (V::splat(1), V::splat(2), V::splat(3));
        let (ins_bit, del_bit) = (V::splat(4), V::splat(8));
        let sentinel = V::splat(SENTINEL);
        let ambiguous_code = V::splat(i16::from(N));
        let (matched, mismatch) = (V::load(&shape.matched), V::load(&shape.mismatch));
        let (mut best, mut best_row, mut best_k) = (zero, zero, zero);

        for r in 1..=rows {
            let q = V::load(&query[r - 1]);
            let mut code_mask = [V::none(); CODES];
            if !UNIFORM {
                for (c, mask) in code_mask.iter_mut().enumerate() {
                    *mask = q.eq(V::splat(c as i16));
                }
            }
            let past_query_lanes = q.eq(sentinel);
            let past_query = neg.keep(past_query_lanes);
            let query_ambiguous = q.eq(ambiguous_code);
            let row_vector = V::splat(r as i16);
            let mut del = neg;
            let mut del_extends = V::none();
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
                        V::load(&h[k]),
                        V::load(&ins[k]),
                        V::load(&h[k + 1]),
                        r + width - 1 - k,
                    )
                } else {
                    let up = V::load(&h[k + 1]);
                    let diag = h_diag_carry;
                    h_diag_carry = up;
                    (up, V::load(&ins[k + 1]), diag, k)
                };
                let s = if UNIFORM {
                    // matched / mismatch by identity, 0 for N, −∞ outside the matrix.
                    let reference = V::load(&profile[t]);
                    let base = V::select(q.eq(reference), matched, mismatch);
                    let ambiguous = V::mask_or(query_ambiguous, reference.eq(ambiguous_code));
                    let outside = V::mask_or(past_query_lanes, reference.gt(ambiguous_code));
                    V::select(outside, neg, base.clear(ambiguous))
                } else {
                    let mut s = past_query;
                    for (c, &mask) in code_mask.iter().enumerate() {
                        s = s.or(V::load(&profile[c * positions + t]).keep(mask));
                    }
                    s
                };

                let ins_open = h_up.adds(open);
                let ins_extend = i_up.adds(extend);
                let mut ins_score = ins_open.max(ins_extend);
                let diag = h_diag.adds(s);
                let t0 = diag.max(zero);
                let t_all = t0.max(ins_score);
                let mut score = t_all.max(del);

                if TRACE {
                    let source = one.keep(diag.gt(zero));
                    let source = V::select(ins_score.gt(t0), two, source);
                    let source = V::select(del.gt(t_all), three, source);
                    let bits = ins_bit
                        .keep(ins_extend.gt(ins_open))
                        .or(del_bit.keep(del_extends));
                    source.or(bits).store_bytes(&mut trace[(r - 1) * width + k]);
                }

                let del_open = t_all.adds(open);
                let mut next_extends = del.adds(extend).gt(del_open);
                let mut next_del = del_open.max(del.adds(step));
                let k_vector = V::splat(k as i16);
                if DIAG {
                    // Lanes whose band is narrower than the group's: k ≥ width_l does not exist.
                    let in_band = lane_width.gt(k_vector);
                    score = score.keep(in_band);
                    ins_score = V::select(in_band, ins_score, neg);
                    next_del = V::select(in_band, next_del, neg);
                    next_extends = V::mask_and(next_extends, in_band);
                }
                score.store(&mut h[k + 1]);
                ins_score.store(&mut ins[k + 1]);
                del = next_del;
                del_extends = next_extends;

                if TRACE {
                    // The end cell seeds traceback; the score pass needs only the score.
                    let improved = score.gt(best);
                    best_row = V::select(improved, row_vector, best_row);
                    best_k = V::select(improved, k_vector, best_k);
                }
                best = best.max(score);
            }
        }
        let mut out = LaneBest {
            score: [0; L],
            row: [0; L],
            k: [0; L],
        };
        best.store(&mut out.score);
        best_row.store(&mut out.row);
        best_k.store(&mut out.k);
        out
    }
}

/// Align one group of at most `L` lanes on backend `K` and emit each lane's result.
fn run_group<K, const L: usize, const DIAG: bool, const TRACE: bool>(
    jobs: &[AlignJob<'_>],
    lanes: &[Lane],
    ws: &mut InterSequenceWorkspace,
    emit: &mut impl FnMut(usize, Option<LaneResult>),
) where
    K: Kernel<L>,
{
    let uniform = lanes
        .iter()
        .all(|l| jobs[l.job].profile.uniform().is_some());
    let shape = if uniform {
        build::<DIAG, true, L>(jobs, lanes, &mut ws.query, &mut ws.profile)
    } else {
        build::<DIAG, false, L>(jobs, lanes, &mut ws.query, &mut ws.profile)
    };
    let width = shape.width;
    let query = ws.query.as_chunks::<L>().0;
    let profile = ws.profile.as_chunks::<L>().0;
    let h = lane_rows::<i16, L>(&mut ws.h, width + 1, 0);
    let ins = lane_rows::<i16, L>(&mut ws.ins, width + 1, NEG);
    let cells = if TRACE { shape.rows * width } else { 0 };
    // Every cell is written before traceback reads it: grow only.
    if ws.trace.len() < cells * L {
        ws.trace.resize(cells * L, 0);
    }
    let trace = &mut ws.trace.as_chunks_mut::<L>().0[..cells];
    // SAFETY: `K` was chosen from `Backend::available`.
    let best = unsafe {
        if uniform {
            K::fill::<DIAG, TRACE, true>(&shape, query, profile, h, ins, trace)
        } else {
            K::fill::<DIAG, TRACE, false>(&shape, query, profile, h, ins, trace)
        }
    };
    let trace = &*trace;
    for (l, lane) in lanes.iter().enumerate() {
        let score = Score::from(best.score[l]);
        if score <= 0 {
            emit(lane.job, None);
            continue;
        }
        let alignment = TRACE.then(|| {
            let (row, col) = lane.cell_of::<DIAG>(best.row[l] as usize, best.k[l] as usize);
            let (q_start, r_start, ops) = traceback_with(
                |i, j| lane.trace_code::<DIAG, L>(trace, width, l, i, j),
                row,
                col,
            );
            Alignment {
                q_start: q_start as u32,
                q_end: row as u32,
                r_start: r_start as u32,
                r_end: col as u32,
                score,
                ops,
            }
        });
        emit(lane.job, Some(LaneResult { score, alignment }));
    }
}

/// Run the diagonal and rectangular groups on backend `K`; route
/// singleton groups to `ws.scalar_jobs`.
fn run_groups<K, const L: usize, const TRACE: bool>(
    jobs: &[AlignJob<'_>],
    ws: &mut InterSequenceWorkspace,
    emit: &mut impl FnMut(usize, Option<LaneResult>),
) where
    K: Kernel<L>,
{
    for diagonal in [true, false] {
        let groups = std::mem::take(if diagonal {
            &mut ws.diagonal
        } else {
            &mut ws.rectangular
        });
        for group in groups.chunks(L) {
            match (group.len(), diagonal) {
                (1, _) => ws.scalar_jobs.push(group[0].job),
                (_, true) => run_group::<K, L, true, TRACE>(jobs, group, ws, emit),
                (_, false) => run_group::<K, L, false, TRACE>(jobs, group, ws, emit),
            }
        }
        *if diagonal {
            &mut ws.diagonal
        } else {
            &mut ws.rectangular
        } = groups;
    }
}

impl InterSequenceGotoh {
    /// Run every batchable group; route the rest to `ws.scalar_jobs`.
    fn simd_groups<const TRACE: bool>(
        &self,
        jobs: &[AlignJob<'_>],
        ws: &mut InterSequenceWorkspace,
        mut emit: impl FnMut(usize, Option<LaneResult>),
    ) {
        ws.scalar_jobs.clear();
        ws.diagonal.clear();
        ws.rectangular.clear();
        let batched = self.backend != Backend::Scalar;
        for (index, job) in jobs.iter().enumerate() {
            match lane_of(index, job).filter(|_| batched) {
                Some(lane) if matches!(job.band, Band::Full) => ws.rectangular.push(lane),
                Some(lane) => ws.diagonal.push(lane),
                None => ws.scalar_jobs.push(index),
            }
        }
        // Similar shapes share registers: less padding per group.
        ws.diagonal
            .sort_by_key(|l| std::cmp::Reverse((l.rows, l.width)));
        ws.rectangular
            .sort_by_key(|l| std::cmp::Reverse((l.rows, l.width)));
        let emit = &mut emit;
        match self.backend {
            #[cfg(target_arch = "x86_64")]
            Backend::Avx512 => run_groups::<Avx512, 32, TRACE>(jobs, ws, emit),
            #[cfg(target_arch = "x86_64")]
            Backend::Avx2 => run_groups::<Avx2, 16, TRACE>(jobs, ws, emit),
            #[cfg(target_arch = "aarch64")]
            Backend::Neon => run_groups::<Neon, 16, TRACE>(jobs, ws, emit),
            #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
            Backend::Simd128 => run_groups::<Simd128, 16, TRACE>(jobs, ws, emit),
            Backend::Portable => run_groups::<Portable, 16, TRACE>(jobs, ws, emit),
            Backend::Scalar => {}
            // Construction admits only available backends.
            backend => unreachable!("the {} kernel is not built for this target", backend.name()),
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
            out[job] = result.map(|r| ScoreHit { score: r.score });
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
            out[job] = result.map(|r| r.alignment.expect("a traceback pass aligns"));
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
    use crate::scoring::{AffineScoring, AnchorColumn, HotspotScoring};

    use crate::testing::{random_codes, rng};

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

    /// Every backend this CPU runs.
    fn backends() -> Vec<InterSequenceGotoh> {
        let available: Vec<_> = Backend::ALL
            .into_iter()
            .filter_map(InterSequenceGotoh::with_backend)
            .collect();
        let names: Vec<_> = available.iter().map(|a| a.backend().name()).collect();
        eprintln!("backends tested: {names:?}");
        available
    }

    /// Detection picks a backend this build runs, never the opt-in
    /// AVX-512, and the scalar kernel on WebAssembly without `simd128`.
    #[test]
    fn detection_picks_a_vectorised_backend_or_the_scalar_kernel() {
        let detected = Backend::detect();
        assert!(detected.available());
        assert_ne!(detected, Backend::Avx512);
        if cfg!(target_family = "wasm") {
            let expected = if cfg!(target_feature = "simd128") {
                Backend::Simd128
            } else {
                Backend::Scalar
            };
            assert_eq!(detected, expected);
        }
    }

    #[test]
    fn every_backend_matches_the_scalar_kernel_on_mixed_batches() {
        for simd in backends() {
            mixed_batches_match_the_scalar_kernel(simd);
        }
    }

    fn mixed_batches_match_the_scalar_kernel(simd: InterSequenceGotoh) {
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
            // Up to five registers of the widest backend.
            let size = 1 + (next() % 160) as usize;
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
                let backend = simd.backend();
                assert_eq!(
                    hits[k], expected_hit,
                    "{backend:?} score, batch {batch} job {k} band {:?}",
                    job.band
                );
                assert_eq!(
                    alignments[k], expected_aln,
                    "{backend:?} alignment, batch {batch} job {k} band {:?}",
                    job.band
                );
                assert_eq!(
                    hits[k].map(|hit| hit.score),
                    alignments[k].as_ref().map(|a| a.score),
                    "{backend:?} align vs score, batch {batch} job {k}"
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
        let mut scalar = GotohWorkspace::default();
        for simd in backends() {
            simd.score_batch(&jobs, &mut ws, &mut hits);
            simd.align_batch(&jobs, &mut ws, &mut alignments);
            for (k, job) in jobs.iter().enumerate() {
                let (query, profile, gaps, band) = (job.query, job.profile, job.gaps, job.band);
                assert_eq!(
                    hits[k],
                    BandedGotoh.score(query, profile, gaps, band, &mut scalar),
                    "{:?} job {k}",
                    simd.backend()
                );
                assert_eq!(
                    alignments[k],
                    BandedGotoh.align(query, profile, gaps, band, &mut scalar),
                    "{:?} job {k}",
                    simd.backend()
                );
            }
        }
    }

    /// C-like jobs: anchor profiles (one outsized column) on 1.3 kb
    /// references under 7-diagonal bands, batched with uniform jobs. They
    /// stay in SIMD lanes and match the scalar kernel.
    #[test]
    fn anchor_profiles_on_long_references_match_the_scalar_kernel() {
        let scoring = AffineScoring::new(2, -3, -9, -1);
        let model = AnchorColumn::new(scoring, 1300);
        let mut next = rng(9);
        let mut cases = Vec::new();
        for k in 0..48 {
            let mut reference = vec![N];
            reference.extend(random_codes(&mut next, 1300));
            let x = 4 + (next() % 4) as usize;
            let len = (next() % 1301) as usize;
            let mut query = random_codes(&mut next, x);
            query.extend(&reference[1..=len]);
            for _ in 0..(next() % 4) {
                let p = (next() as usize) % query.len();
                query[p] = (next() % 5) as u8;
            }
            if len > 2 && next().is_multiple_of(2) {
                query.remove(x + 1 + (next() as usize) % (len - 1));
            }
            let tail = (next() % 40) as usize;
            query.extend(random_codes(&mut next, tail));
            let profile = if k % 4 == 3 {
                SubstitutionProfile::build(&scoring, &reference)
            } else {
                SubstitutionProfile::build(&model, &reference)
            };
            let band = Band::Diagonals {
                lo: x as i32 - 4,
                hi: x as i32 + 2,
            };
            cases.push((query, profile, band));
        }
        let jobs: Vec<AlignJob> = cases
            .iter()
            .map(|(query, profile, band)| AlignJob {
                query,
                profile,
                gaps: scoring.gaps,
                band: *band,
            })
            .collect();
        assert!(jobs.iter().all(|job| lane_of(0, job).is_some()));
        let (mut hits, mut alignments) = (Vec::new(), Vec::new());
        let mut ws = InterSequenceWorkspace::default();
        let mut scalar = GotohWorkspace::default();
        for simd in backends() {
            simd.score_batch(&jobs, &mut ws, &mut hits);
            simd.align_batch(&jobs, &mut ws, &mut alignments);
            for (k, job) in jobs.iter().enumerate() {
                let (query, profile, gaps, band) = (job.query, job.profile, job.gaps, job.band);
                let backend = simd.backend();
                assert_eq!(
                    hits[k],
                    BandedGotoh.score(query, profile, gaps, band, &mut scalar),
                    "{backend:?} job {k}"
                );
                assert_eq!(
                    alignments[k],
                    BandedGotoh.align(query, profile, gaps, band, &mut scalar),
                    "{backend:?} job {k}"
                );
            }
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
        let (mut out, mut alignments) = (Vec::new(), Vec::new());
        let mut ws = InterSequenceWorkspace::default();
        InterSequenceGotoh::new().score_batch(&jobs, &mut ws, &mut out);
        assert!(out.iter().all(|hit| hit.map(|h| h.score) == Some(20_000)));
        // The scalar fallback's traceback reports the same score.
        InterSequenceGotoh::new().align_batch(&jobs, &mut ws, &mut alignments);
        assert!(
            alignments
                .iter()
                .all(|a| a.as_ref().map(|a| a.score) == Some(20_000))
        );
    }
}
