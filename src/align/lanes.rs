//! SIMD lanes for the inter-sequence kernel: one trait, one implementation
//! per instruction set.
//!
//! [`Lanes<L>`] is a register of `L` signed 16-bit lanes with exactly the
//! operations the Gotoh recurrence needs. Comparisons yield an associated
//! [`Lanes::Mask`] — a vector of all-ones / all-zeros lanes on AVX2, NEON and
//! the portable backend, a `k` register on AVX-512 — so each instruction set
//! keeps its natural form. The kernel ([`super::inter_seq`]) is written once,
//! generic over `Lanes<L>`, and each backend compiles it inside its own
//! `#[target_feature]` entry point ([`Kernel`]), where every operation inlines
//! to one or two instructions.
//!
//! | Backend | Lanes | Registers per vector |
//! |---|---|---|
//! | `Avx512` (x86-64 AVX-512BW) | 32 | one `zmm`, masks in `k` registers |
//! | `Avx2` (x86-64 AVX2) | 16 | one `ymm` |
//! | `Neon` (AArch64) | 16 | two `q` registers |
//! | [`Portable`] (any target) | 16 | plain arrays; the compiler vectorises them for the target's baseline SIMD (SSE2, NEON, WebAssembly `simd128`, …) |
//!
//! # Safety
//!
//! Every method is `unsafe`: callers must run on a CPU with the backend's
//! target features ([`super::inter_seq::Backend::available`]).

use super::inter_seq::{LaneBest, Shape};

/// `L` signed 16-bit lanes, each holding the same DP cell of a different job.
pub(super) trait Lanes<const L: usize>: Copy {
    /// Per-lane truth values from comparisons.
    type Mask: Copy;

    unsafe fn load(v: &[i16; L]) -> Self;
    unsafe fn store(self, v: &mut [i16; L]);
    /// Store each lane as a byte; lanes must hold `0..=255`.
    unsafe fn store_bytes(self, v: &mut [u8; L]);
    unsafe fn splat(x: i16) -> Self;
    /// Saturating addition.
    unsafe fn adds(self, other: Self) -> Self;
    unsafe fn max(self, other: Self) -> Self;
    unsafe fn or(self, other: Self) -> Self;
    unsafe fn eq(self, other: Self) -> Self::Mask;
    unsafe fn gt(self, other: Self) -> Self::Mask;
    /// `yes` where `mask` holds, else `no`.
    unsafe fn select(mask: Self::Mask, yes: Self, no: Self) -> Self;
    /// `self` where `mask` holds, else 0.
    unsafe fn keep(self, mask: Self::Mask) -> Self;
    /// 0 where `mask` holds, else `self`.
    unsafe fn clear(self, mask: Self::Mask) -> Self;
    unsafe fn mask_or(a: Self::Mask, b: Self::Mask) -> Self::Mask;
    unsafe fn mask_and(a: Self::Mask, b: Self::Mask) -> Self::Mask;
    /// A mask holding in no lane.
    unsafe fn none() -> Self::Mask;
}

/// The inter-sequence kernel compiled for one backend's target features.
pub(super) trait Kernel<const L: usize> {
    /// Fill one group; see [`super::inter_seq`] for the layout.
    ///
    /// # Safety
    ///
    /// The CPU must support the backend.
    unsafe fn fill<const DIAG: bool, const TRACE: bool, const UNIFORM: bool>(
        shape: &Shape<L>,
        query: &[[i16; L]],
        profile: &[[i16; L]],
        h: &mut [[i16; L]],
        ins: &mut [[i16; L]],
        trace: &mut [[u8; L]],
    ) -> LaneBest<L>;
}

/// Implement [`Kernel`] for a backend by compiling the generic kernel under
/// its target features.
macro_rules! kernel {
    ($backend:ty, $lanes:literal $(, $feature:literal)?) => {
        impl Kernel<$lanes> for $backend {
            $(#[target_feature(enable = $feature)])?
            unsafe fn fill<const DIAG: bool, const TRACE: bool, const UNIFORM: bool>(
                shape: &Shape<$lanes>,
                query: &[[i16; $lanes]],
                profile: &[[i16; $lanes]],
                h: &mut [[i16; $lanes]],
                ins: &mut [[i16; $lanes]],
                trace: &mut [[u8; $lanes]],
            ) -> LaneBest<$lanes> {
                // SAFETY: forwarded from the caller.
                unsafe {
                    crate::align::inter_seq::fill::<Self, $lanes, DIAG, TRACE, UNIFORM>(
                        shape, query, profile, h, ins, trace,
                    )
                }
            }
        }
    };
}

/// Plain Rust lanes: correct on every target, vectorised by the compiler.
#[derive(Clone, Copy, Debug)]
pub(super) struct Portable([i16; 16]);

impl Portable {
    #[inline(always)]
    fn map(self, other: Self, f: impl Fn(i16, i16) -> i16) -> Self {
        Self(std::array::from_fn(|i| f(self.0[i], other.0[i])))
    }
}

impl Lanes<16> for Portable {
    /// All-ones (−1) or all-zeros lanes.
    type Mask = Portable;

    #[inline(always)]
    unsafe fn load(v: &[i16; 16]) -> Self {
        Self(*v)
    }
    #[inline(always)]
    unsafe fn store(self, v: &mut [i16; 16]) {
        *v = self.0;
    }
    #[inline(always)]
    unsafe fn store_bytes(self, v: &mut [u8; 16]) {
        *v = self.0.map(|x| x as u8);
    }
    #[inline(always)]
    unsafe fn splat(x: i16) -> Self {
        Self([x; 16])
    }
    #[inline(always)]
    unsafe fn adds(self, other: Self) -> Self {
        self.map(other, i16::saturating_add)
    }
    #[inline(always)]
    unsafe fn max(self, other: Self) -> Self {
        self.map(other, i16::max)
    }
    #[inline(always)]
    unsafe fn or(self, other: Self) -> Self {
        self.map(other, |a, b| a | b)
    }
    #[inline(always)]
    unsafe fn eq(self, other: Self) -> Self {
        self.map(other, |a, b| -i16::from(a == b))
    }
    #[inline(always)]
    unsafe fn gt(self, other: Self) -> Self {
        self.map(other, |a, b| -i16::from(a > b))
    }
    #[inline(always)]
    unsafe fn select(mask: Self, yes: Self, no: Self) -> Self {
        Self(std::array::from_fn(|i| {
            (yes.0[i] & mask.0[i]) | (no.0[i] & !mask.0[i])
        }))
    }
    #[inline(always)]
    unsafe fn keep(self, mask: Self) -> Self {
        self.map(mask, |a, m| a & m)
    }
    #[inline(always)]
    unsafe fn clear(self, mask: Self) -> Self {
        self.map(mask, |a, m| a & !m)
    }
    #[inline(always)]
    unsafe fn mask_or(a: Self, b: Self) -> Self {
        a.map(b, |a, b| a | b)
    }
    #[inline(always)]
    unsafe fn mask_and(a: Self, b: Self) -> Self {
        a.map(b, |a, b| a & b)
    }
    #[inline(always)]
    unsafe fn none() -> Self {
        Self([0; 16])
    }
}

kernel!(Portable, 16);

#[cfg(target_arch = "x86_64")]
pub(super) use x86::{Avx2, Avx512};

#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::{Kernel, LaneBest, Lanes, Shape};
    use std::arch::x86_64::*;

    /// AVX2: 16 lanes in one `ymm` register; masks are vectors.
    #[derive(Clone, Copy, Debug)]
    pub struct Avx2(__m256i);

    impl Lanes<16> for Avx2 {
        type Mask = __m256i;

        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn load(v: &[i16; 16]) -> Self {
            Self(unsafe { _mm256_loadu_si256(v.as_ptr().cast()) })
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn store(self, v: &mut [i16; 16]) {
            unsafe { _mm256_storeu_si256(v.as_mut_ptr().cast(), self.0) }
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn store_bytes(self, v: &mut [u8; 16]) {
            let x = self.0;
            let packed =
                _mm_packus_epi16(_mm256_castsi256_si128(x), _mm256_extracti128_si256::<1>(x));
            unsafe { _mm_storeu_si128(v.as_mut_ptr().cast(), packed) }
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn splat(x: i16) -> Self {
            Self(_mm256_set1_epi16(x))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn adds(self, other: Self) -> Self {
            Self(_mm256_adds_epi16(self.0, other.0))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn max(self, other: Self) -> Self {
            Self(_mm256_max_epi16(self.0, other.0))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn or(self, other: Self) -> Self {
            Self(_mm256_or_si256(self.0, other.0))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn eq(self, other: Self) -> __m256i {
            _mm256_cmpeq_epi16(self.0, other.0)
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn gt(self, other: Self) -> __m256i {
            _mm256_cmpgt_epi16(self.0, other.0)
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn select(mask: __m256i, yes: Self, no: Self) -> Self {
            Self(_mm256_blendv_epi8(no.0, yes.0, mask))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn keep(self, mask: __m256i) -> Self {
            Self(_mm256_and_si256(self.0, mask))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn clear(self, mask: __m256i) -> Self {
            Self(_mm256_andnot_si256(mask, self.0))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn mask_or(a: __m256i, b: __m256i) -> __m256i {
            _mm256_or_si256(a, b)
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn mask_and(a: __m256i, b: __m256i) -> __m256i {
            _mm256_and_si256(a, b)
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn none() -> __m256i {
            _mm256_setzero_si256()
        }
    }

    kernel!(Avx2, 16, "avx2");

    /// AVX-512BW: 32 lanes in one `zmm` register; masks in `k` registers.
    #[derive(Clone, Copy, Debug)]
    pub struct Avx512(__m512i);

    impl Lanes<32> for Avx512 {
        type Mask = __mmask32;

        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn load(v: &[i16; 32]) -> Self {
            Self(unsafe { _mm512_loadu_si512(v.as_ptr().cast()) })
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn store(self, v: &mut [i16; 32]) {
            unsafe { _mm512_storeu_si512(v.as_mut_ptr().cast(), self.0) }
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn store_bytes(self, v: &mut [u8; 32]) {
            unsafe { _mm256_storeu_si256(v.as_mut_ptr().cast(), _mm512_cvtepi16_epi8(self.0)) }
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn splat(x: i16) -> Self {
            Self(_mm512_set1_epi16(x))
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn adds(self, other: Self) -> Self {
            Self(_mm512_adds_epi16(self.0, other.0))
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn max(self, other: Self) -> Self {
            Self(_mm512_max_epi16(self.0, other.0))
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn or(self, other: Self) -> Self {
            Self(_mm512_or_si512(self.0, other.0))
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn eq(self, other: Self) -> __mmask32 {
            _mm512_cmpeq_epi16_mask(self.0, other.0)
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn gt(self, other: Self) -> __mmask32 {
            _mm512_cmpgt_epi16_mask(self.0, other.0)
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn select(mask: __mmask32, yes: Self, no: Self) -> Self {
            Self(_mm512_mask_blend_epi16(mask, no.0, yes.0))
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn keep(self, mask: __mmask32) -> Self {
            Self(_mm512_maskz_mov_epi16(mask, self.0))
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn clear(self, mask: __mmask32) -> Self {
            Self(_mm512_maskz_mov_epi16(!mask, self.0))
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn mask_or(a: __mmask32, b: __mmask32) -> __mmask32 {
            a | b
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn mask_and(a: __mmask32, b: __mmask32) -> __mmask32 {
            a & b
        }
        #[inline]
        #[target_feature(enable = "avx512bw")]
        unsafe fn none() -> __mmask32 {
            0
        }
    }

    kernel!(Avx512, 32, "avx512bw");
}

#[cfg(target_arch = "aarch64")]
pub(super) use arm::Neon;

#[cfg(target_arch = "aarch64")]
mod arm {
    use super::{Kernel, LaneBest, Lanes, Shape};
    use std::arch::aarch64::*;

    /// NEON: 16 lanes in two `q` registers; masks are vectors.
    #[derive(Clone, Copy, Debug)]
    pub struct Neon(int16x8_t, int16x8_t);

    /// Two-register NEON mask (all-ones / all-zeros lanes).
    #[derive(Clone, Copy, Debug)]
    pub struct NeonMask(uint16x8_t, uint16x8_t);

    impl Lanes<16> for Neon {
        type Mask = NeonMask;

        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn load(v: &[i16; 16]) -> Self {
            unsafe { Self(vld1q_s16(v.as_ptr()), vld1q_s16(v.as_ptr().add(8))) }
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn store(self, v: &mut [i16; 16]) {
            unsafe {
                vst1q_s16(v.as_mut_ptr(), self.0);
                vst1q_s16(v.as_mut_ptr().add(8), self.1);
            }
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn store_bytes(self, v: &mut [u8; 16]) {
            let packed = vcombine_u8(vqmovun_s16(self.0), vqmovun_s16(self.1));
            unsafe { vst1q_u8(v.as_mut_ptr(), packed) }
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn splat(x: i16) -> Self {
            Self(vdupq_n_s16(x), vdupq_n_s16(x))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn adds(self, other: Self) -> Self {
            Self(vqaddq_s16(self.0, other.0), vqaddq_s16(self.1, other.1))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn max(self, other: Self) -> Self {
            Self(vmaxq_s16(self.0, other.0), vmaxq_s16(self.1, other.1))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn or(self, other: Self) -> Self {
            Self(vorrq_s16(self.0, other.0), vorrq_s16(self.1, other.1))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn eq(self, other: Self) -> NeonMask {
            NeonMask(vceqq_s16(self.0, other.0), vceqq_s16(self.1, other.1))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn gt(self, other: Self) -> NeonMask {
            NeonMask(vcgtq_s16(self.0, other.0), vcgtq_s16(self.1, other.1))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn select(mask: NeonMask, yes: Self, no: Self) -> Self {
            Self(
                vbslq_s16(mask.0, yes.0, no.0),
                vbslq_s16(mask.1, yes.1, no.1),
            )
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn keep(self, mask: NeonMask) -> Self {
            Self(
                vandq_s16(self.0, vreinterpretq_s16_u16(mask.0)),
                vandq_s16(self.1, vreinterpretq_s16_u16(mask.1)),
            )
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn clear(self, mask: NeonMask) -> Self {
            Self(
                vbicq_s16(self.0, vreinterpretq_s16_u16(mask.0)),
                vbicq_s16(self.1, vreinterpretq_s16_u16(mask.1)),
            )
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn mask_or(a: NeonMask, b: NeonMask) -> NeonMask {
            NeonMask(vorrq_u16(a.0, b.0), vorrq_u16(a.1, b.1))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn mask_and(a: NeonMask, b: NeonMask) -> NeonMask {
            NeonMask(vandq_u16(a.0, b.0), vandq_u16(a.1, b.1))
        }
        #[inline]
        #[target_feature(enable = "neon")]
        unsafe fn none() -> NeonMask {
            NeonMask(vdupq_n_u16(0), vdupq_n_u16(0))
        }
    }

    kernel!(Neon, 16, "neon");
}
