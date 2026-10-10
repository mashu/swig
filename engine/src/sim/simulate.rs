//! Recombination, assembly and noise: one simulated read with exact truth.
//!
//! The generative model is IgSim.jl's, step for step:
//!
//! 1. **Recombination**, in biological order. Sample J and trim its 5′ end.
//!    With probability `include_d` sample a D, trim both ends and add N2; a D
//!    trimmed away entirely means no D. With probability `include_dd` a
//!    second, different D joins that DJ on its 5′ side, separated by N3
//!    (DJ → DDJ). Finally sample V, trim it and add N1. Alleles are sampled
//!    uniformly.
//! 2. **Assembly**: `flank5 + V + N1 + [D + [N3 + D2] + N2] + J + flank3`,
//!    flanks being random DNA. Reads outside `[min_length, max_length]` are
//!    redrawn.
//! 3. **SHM** on body (non-flank) bases: single-base indels at the read's
//!    indel rate (half deletions, half insertions after the base), then
//!    independent substitutions at the read's SHM rate.
//! 4. **Instrument** substitutions on body bases, rate ramping 5′ → 3′.
//! 5. Optionally (not in IgSim) the emitted read is reverse-complemented.
//!
//! With a C database (not in IgSim), [`Simulator::append_constant`] then
//! reads through the J–C splice: the start of the 3′ flank becomes the C
//! allele's germline, from the base that follows the J's last base in the
//! reading frame. It draws from a stream of its own, so the V(D)J part of
//! every read is the same with and without C.
//!
//! Unlike IgSim, which labels bases by region only, every base here also
//! carries its germline **origin** (position in its allele, or none for N,
//! flank and inserted bases). That makes boundaries, junctions and
//! germline-equivalent alleles exact ground truth even across indels.

use std::fmt;
use std::ops::Range;

use super::params::SimParams;
use super::random::{Count, Rng};
use crate::alphabet::decode_base;
use crate::germline::{GermlineSet, SegmentKind};
use crate::translate::open_frames;

/// The C stream's seed is the read stream's seed XOR this
/// ([`Simulator::reads`]): C draws never shift V(D)J draws.
pub const C_STREAM: u64 = 0x0C0C_0C0C_0C0C_0C0C;

/// Where a read base came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Region {
    Flank5,
    V,
    /// V–D (or V–J) non-templated addition.
    N1,
    /// First D in read order.
    D,
    /// Between two tandem Ds.
    N3,
    /// Second D in read order (tandem DD only).
    D2,
    /// D–J non-templated addition.
    N2,
    J,
    /// Constant region, read through the J–C splice
    /// ([`Simulator::append_constant`]).
    C,
    Flank3,
}

impl Region {
    /// The rearranged V(D)J bases, which carry SHM, indels and the
    /// instrument ramp. C is spliced RNA: no SHM, no indels, instrument
    /// errors at the 3′ rate.
    pub fn is_body(self) -> bool {
        !matches!(self, Region::Flank5 | Region::C | Region::Flank3)
    }
}

/// No germline origin (N addition, flank, inserted base).
pub const NO_ORIGIN: u32 = u32::MAX;

/// Exonuclease trims, in nucleotides.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Trims {
    pub v_5p: u32,
    pub v_3p: u32,
    pub d_5p: u32,
    pub d_3p: u32,
    pub d2_5p: u32,
    pub d2_3p: u32,
    pub j_5p: u32,
}

/// One simulated read and its complete ground truth.
#[derive(Clone, Debug)]
pub struct SimulatedRead {
    /// The emitted read (ASCII), reverse-complemented when `reverse_complemented`.
    pub sequence: Vec<u8>,
    /// The read in germline sense; `regions` and `origins` index this.
    pub forward: Vec<u8>,
    pub regions: Vec<Region>,
    /// Germline position of each base in its allele, or [`NO_ORIGIN`].
    pub origins: Vec<u32>,
    /// Allele indices into the simulator's germline sets.
    pub v: usize,
    /// First D in read order.
    pub d: Option<usize>,
    /// Second D in read order (tandem DD).
    pub d2: Option<usize>,
    pub j: usize,
    /// C allele, when the read reaches the constant region.
    pub c: Option<usize>,
    pub trims: Trims,
    pub flank_5p: u32,
    pub flank_3p: u32,
    /// Sampled per-base SHM rate of this read.
    pub shm_rate: f64,
    /// SHM substitutions, indels and instrument substitutions applied (C's
    /// included).
    pub mutations: u32,
    pub indels: u32,
    pub instrument_errors: u32,
    pub reverse_complemented: bool,
}

impl SimulatedRead {
    /// `[first, last + 1)` of the bases labelled `region` in the forward read.
    pub fn span(&self, region: Region) -> Option<Range<usize>> {
        let first = self.regions.iter().position(|&r| r == region)?;
        let last = self.regions.iter().rposition(|&r| r == region)?;
        Some(first..last + 1)
    }

    /// Forward-read index of the base templated from germline position `at`
    /// of `region`'s allele, if it survived trimming and indels.
    pub fn locate(&self, region: Region, at: u32) -> Option<usize> {
        self.regions
            .iter()
            .zip(&self.origins)
            .position(|(&r, &o)| r == region && o == at)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SimError {
    EmptySet(SegmentKind),
    WrongKind {
        expected: SegmentKind,
    },
    /// No read met the length constraints within `max_retries` attempts.
    Retries(usize),
}

impl fmt::Display for SimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySet(kind) => write!(f, "the {} germline set is empty", kind.label()),
            Self::WrongKind { expected } => {
                write!(f, "expected a {} germline set", expected.label())
            }
            Self::Retries(n) => write!(f, "no read met the length limits after {n} attempts"),
        }
    }
}

impl std::error::Error for SimError {}

/// Read generator over V, optional D, J and optional C germline sets.
pub struct Simulator<'a> {
    v: Vec<Vec<u8>>,
    d: Vec<Vec<u8>>,
    j: Vec<Vec<u8>>,
    c: Vec<Vec<u8>>,
    /// Stop-free reading frame of each C allele.
    c_frames: Vec<Option<usize>>,
    sets: [Option<&'a GermlineSet>; 4],
    pub params: SimParams,
}

#[derive(Clone, Copy)]
struct Base {
    base: u8,
    region: Region,
    origin: u32,
}

fn ascii(set: &GermlineSet) -> Vec<Vec<u8>> {
    set.alleles()
        .iter()
        .map(|a| a.codes.iter().map(|&c| decode_base(c)).collect())
        .collect()
}

fn check(set: &GermlineSet, kind: SegmentKind) -> Result<(), SimError> {
    if set.kind() != kind {
        return Err(SimError::WrongKind { expected: kind });
    }
    if set.is_empty() {
        return Err(SimError::EmptySet(kind));
    }
    Ok(())
}

/// `sequence[t5 .. len − t3]` with positions, empty when the trims meet.
fn trimmed(sequence: &[u8], t5: u32, t3: u32) -> Range<usize> {
    let len = sequence.len();
    let end = len.saturating_sub(t3 as usize);
    let start = (t5 as usize).min(end);
    start..end
}

impl<'a> Simulator<'a> {
    pub fn new(
        v: &'a GermlineSet,
        d: Option<&'a GermlineSet>,
        j: &'a GermlineSet,
        params: SimParams,
    ) -> Result<Self, SimError> {
        check(v, SegmentKind::V)?;
        if let Some(d) = d {
            check(d, SegmentKind::D)?;
        }
        check(j, SegmentKind::J)?;
        Ok(Self {
            v: ascii(v),
            d: d.map(ascii).unwrap_or_default(),
            j: ascii(j),
            c: Vec::new(),
            c_frames: Vec::new(),
            sets: [Some(v), d, Some(j), None],
            params,
        })
    }

    /// Read into the constant regions of `c` ([`Self::append_constant`]).
    pub fn with_constant(mut self, c: &'a GermlineSet) -> Result<Self, SimError> {
        check(c, SegmentKind::C)?;
        self.c = ascii(c);
        self.c_frames = c
            .alleles()
            .iter()
            .map(|a| open_frames(&a.codes).next())
            .collect();
        self.sets[SegmentKind::C.index()] = Some(c);
        Ok(self)
    }

    /// `o`, the germline position of C allele `c` that follows the last
    /// base of J allele `j`:
    ///
    /// ```text
    /// o = (f_C + φ_J) mod 3      f_C: C's stop-free frame
    ///                            φ_J = (|J| − anchor_J) mod 3: J bases past its last codon
    /// ```
    ///
    /// so the splice keeps the reading frame. It is 0 for records that
    /// start after the J's last base (IMGT CH1: `φ_J = 1`, `f_C = 2`) and 1
    /// for codon-aligned records that repeat it (`f_C = 0`); 0 without a
    /// frame or a J anchor.
    pub fn constant_offset(&self, c: usize, j: usize) -> usize {
        let j = self.sets[SegmentKind::J.index()].map(|set| set.allele(j as u32));
        let phase = j.and_then(|j| Some((j.codes.len() - j.anchor? as usize) % 3));
        self.c_frames[c].zip(phase).map_or(0, |(f, p)| (f + p) % 3)
    }

    /// The germline set a simulated allele index refers to.
    pub fn set(&self, kind: SegmentKind) -> Option<&'a GermlineSet> {
        self.sets[kind.index()]
    }

    /// The trimmed germline body of segment `region` of `read`, before noise.
    pub fn germline_body(&self, read: &SimulatedRead, region: Region) -> Option<&[u8]> {
        let t = &read.trims;
        let (allele, t5, t3) = match region {
            Region::V => (&self.v[read.v], t.v_5p, t.v_3p),
            Region::D => (&self.d[read.d?], t.d_5p, t.d_3p),
            Region::D2 => (&self.d[read.d2?], t.d2_5p, t.d2_3p),
            Region::J => (&self.j[read.j], t.j_5p, 0),
            _ => return None,
        };
        Some(&allele[trimmed(allele, t5, t3)])
    }

    /// `n` reads from fresh generators seeded with `seed` (V(D)J) and
    /// `seed ^` [`C_STREAM`] (C).
    pub fn reads(&self, n: usize, seed: u64) -> Result<Vec<SimulatedRead>, SimError> {
        let (mut rng, mut c_rng) = (Rng::new(seed), Rng::new(seed ^ C_STREAM));
        (0..n)
            .map(|_| {
                let mut read = self.read(&mut rng)?;
                self.append_constant(&mut read, &mut c_rng);
                Ok(read)
            })
            .collect()
    }

    /// With probability `include_c`, read through the J–C splice: a C
    /// allele drawn uniformly overwrites the start of the 3′ flank from
    /// germline position `o` ([`Self::constant_offset`]), for the whole
    /// flank with probability `c_whole_flank`, else for a length uniform in
    /// `0..=flank`, and never past the allele's end; then instrument errors
    /// at the 3′ rate. Draws only from `rng`, so call it with a stream of
    /// its own; without a C database it does nothing.
    pub fn append_constant(&self, read: &mut SimulatedRead, rng: &mut Rng) {
        let p = &self.params;
        if self.c.is_empty() || !rng.bernoulli(p.include_c) {
            return;
        }
        let allele = rng.below(self.c.len());
        let offset = self.constant_offset(allele, read.j);
        let Some(flank) = read.span(Region::Flank3) else {
            return;
        };
        let len = if rng.bernoulli(p.c_whole_flank) {
            flank.len()
        } else {
            rng.range(0, flank.len() as u32) as usize
        };
        let germline = self.c[allele].get(offset..).unwrap_or_default();
        let len = len.min(germline.len());
        if len == 0 {
            return;
        }
        let rate = p.illumina_error.p3;
        for (k, at) in (flank.start..flank.start + len).enumerate() {
            let mut base = germline[k];
            if rate > 0.0 && rng.uniform() < rate {
                base = rng.substitute(base);
                read.instrument_errors += 1;
            }
            read.forward[at] = base;
            read.regions[at] = Region::C;
            read.origins[at] = (offset + k) as u32;
        }
        read.c = Some(allele);
        read.sequence = if read.reverse_complemented {
            reverse_complement(&read.forward)
        } else {
            read.forward.clone()
        };
    }

    /// One read (IgSim's `gen(rng)`).
    pub fn read(&self, rng: &mut Rng) -> Result<SimulatedRead, SimError> {
        for _ in 0..self.params.max_retries.max(1) {
            if let Some(read) = self.attempt(rng) {
                return Ok(read);
            }
        }
        Err(SimError::Retries(self.params.max_retries))
    }

    fn attempt(&self, rng: &mut Rng) -> Option<SimulatedRead> {
        let p = &self.params;
        let mut trims = Trims::default();

        // 1. DJ, optional DDJ, then V (IgSim's order of draws).
        let j = rng.below(self.j.len());
        trims.j_5p = p.j_trim_5p.sample(rng);
        let j_body = trimmed(&self.j[j], trims.j_5p, 0);
        let (mut d, mut d2) = (None, None);
        let (mut n2, mut n3) = (Vec::new(), Vec::new());
        if !self.d.is_empty() && rng.bernoulli(p.include_d) {
            let first = rng.below(self.d.len());
            let (t5, t3) = (p.d_trim_5p.sample(rng), p.d_trim_3p.sample(rng));
            if !trimmed(&self.d[first], t5, t3).is_empty() {
                n2 = random_n(rng, &p.n2_length);
                d = Some(first);
                (trims.d_5p, trims.d_3p) = (t5, t3);
                if self.d.len() >= 2 && rng.bernoulli(p.include_dd) {
                    let mut other = rng.below(self.d.len() - 1);
                    if other >= first {
                        other += 1;
                    }
                    let (u5, u3) = (p.d_trim_5p.sample(rng), p.d_trim_3p.sample(rng));
                    if !trimmed(&self.d[other], u5, u3).is_empty() {
                        // D + DJ → DDJ: the new D is 5′ of the original one.
                        n3 = random_n(rng, &p.n2_length);
                        d2 = d;
                        (trims.d2_5p, trims.d2_3p) = (t5, t3);
                        d = Some(other);
                        (trims.d_5p, trims.d_3p) = (u5, u3);
                    }
                }
            }
        }
        let v = rng.below(self.v.len());
        trims.v_5p = p.v_trim_5p.sample(rng);
        trims.v_3p = p.v_trim_3p.sample(rng);
        let v_body = trimmed(&self.v[v], trims.v_5p, trims.v_3p);
        if v_body.is_empty() || j_body.is_empty() {
            return None;
        }
        let n1 = random_n(rng, &p.n1_length);

        // 2. Assembly with per-base truth.
        let flank_5p = p.flank_5p.sample(rng);
        let flank_3p = p.flank_3p.sample(rng);
        let mut bases: Vec<Base> = Vec::with_capacity(700);
        let push_random = |bases: &mut Vec<Base>, seq: &[u8], region| {
            bases.extend(seq.iter().map(|&base| Base {
                base,
                region,
                origin: NO_ORIGIN,
            }));
        };
        let push_germline = |bases: &mut Vec<Base>, allele: &[u8], range: Range<usize>, region| {
            bases.extend(range.map(|at| Base {
                base: allele[at],
                region,
                origin: at as u32,
            }));
        };
        push_random(&mut bases, &random_dna(rng, flank_5p), Region::Flank5);
        push_germline(&mut bases, &self.v[v], v_body, Region::V);
        push_random(&mut bases, &n1, Region::N1);
        if let Some(first) = d {
            let range = trimmed(&self.d[first], trims.d_5p, trims.d_3p);
            push_germline(&mut bases, &self.d[first], range, Region::D);
            if let Some(second) = d2 {
                push_random(&mut bases, &n3, Region::N3);
                let range = trimmed(&self.d[second], trims.d2_5p, trims.d2_3p);
                push_germline(&mut bases, &self.d[second], range, Region::D2);
            }
            push_random(&mut bases, &n2, Region::N2);
        }
        push_germline(&mut bases, &self.j[j], j_body, Region::J);
        push_random(&mut bases, &random_dna(rng, flank_3p), Region::Flank3);
        if bases.len() < p.min_length || bases.len() > p.max_length {
            return None;
        }

        // 3. SHM: indels, then substitutions; 4. instrument substitutions.
        let shm_rate = p.body_error_rate.sample(rng);
        let indel_rate = p.indel_rate.sample(rng);
        let mut indels = 0;
        if indel_rate > 0.0 {
            let mut out = Vec::with_capacity(bases.len() + 8);
            for b in bases {
                if b.region.is_body() && rng.uniform() < indel_rate {
                    indels += 1;
                    if rng.uniform() < 0.5 && out.len() > 10 {
                        continue;
                    }
                    out.push(b);
                    out.push(Base {
                        base: rng.base(),
                        region: b.region,
                        origin: NO_ORIGIN,
                    });
                    continue;
                }
                out.push(b);
            }
            bases = out;
        }
        let mut mutations = 0;
        if shm_rate > 0.0 {
            for b in bases.iter_mut().filter(|b| b.region.is_body()) {
                if rng.uniform() < shm_rate {
                    b.base = rng.substitute(b.base);
                    mutations += 1;
                }
            }
        }
        let mut instrument_errors = 0;
        let (p5, p3) = (p.illumina_error.p5, p.illumina_error.p3);
        let body = bases.iter().filter(|b| b.region.is_body()).count();
        if (p5 > 0.0 || p3 > 0.0) && body > 0 {
            for (k, b) in bases.iter_mut().filter(|b| b.region.is_body()).enumerate() {
                let t = if body == 1 {
                    0.0
                } else {
                    k as f64 / (body - 1) as f64
                };
                let rate = p5 + t * (p3 - p5);
                if rate > 0.0 && rng.uniform() < rate {
                    b.base = rng.substitute(b.base);
                    instrument_errors += 1;
                }
            }
        }

        // 5. Optional reverse complement of the emitted read.
        let forward: Vec<u8> = bases.iter().map(|b| b.base).collect();
        let reverse_complemented = rng.bernoulli(p.reverse_complement);
        let sequence = if reverse_complemented {
            reverse_complement(&forward)
        } else {
            forward.clone()
        };
        Some(SimulatedRead {
            sequence,
            forward,
            regions: bases.iter().map(|b| b.region).collect(),
            origins: bases.iter().map(|b| b.origin).collect(),
            v,
            d,
            d2,
            j,
            c: None,
            trims,
            flank_5p,
            flank_3p,
            shm_rate,
            mutations,
            indels,
            instrument_errors,
            reverse_complemented,
        })
    }
}

fn random_dna(rng: &mut Rng, len: u32) -> Vec<u8> {
    (0..len).map(|_| rng.base()).collect()
}

/// A non-templated addition of random length.
fn random_n(rng: &mut Rng, length: &Count) -> Vec<u8> {
    let len = length.sample(rng);
    random_dna(rng, len)
}

/// Reverse complement of an ASCII sequence (`N` and other symbols map to `N`).
pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::random::GatedRate;

    fn sets() -> (GermlineSet, GermlineSet, GermlineSet) {
        let mut v = GermlineSet::new(SegmentKind::V);
        v.push(
            "V1*01",
            b"CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTC",
            None,
        )
        .unwrap();
        v.push(
            "V2*01",
            b"GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTC",
            None,
        )
        .unwrap();
        let mut d = GermlineSet::new(SegmentKind::D);
        d.push("D1*01", b"GTATTACGATATTTTGACTGGTTATTATAAC", None)
            .unwrap();
        d.push("D2*01", b"AGGATATTGTAGTAGTACCAGCTGCTATGCC", None)
            .unwrap();
        let mut j = GermlineSet::new(SegmentKind::J);
        j.push(
            "J1*01",
            b"ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
            None,
        )
        .unwrap();
        (v, d, j)
    }

    #[test]
    fn truth_labels_reconstruct_the_germline_bodies() {
        let (v, d, j) = sets();
        let params = SimParams {
            include_dd: 0.3,
            min_length: 0,
            illumina_error: crate::sim::IlluminaError::NONE,
            ..SimParams::easy()
        };
        let sim = Simulator::new(&v, Some(&d), &j, params).unwrap();
        let mut rng = Rng::new(3);
        let mut tandem = 0;
        for _ in 0..300 {
            let read = sim.read(&mut rng).unwrap();
            assert_eq!(read.forward.len(), read.regions.len());
            for region in [Region::V, Region::D, Region::D2, Region::J] {
                let Some(body) = sim.germline_body(&read, region) else {
                    assert!(read.span(region).is_none());
                    continue;
                };
                // No noise: the labelled bases are the germline body.
                let span = read.span(region).unwrap();
                assert_eq!(&read.forward[span], body, "{region:?}");
            }
            tandem += usize::from(read.d2.is_some());
            if read.d2.is_some() {
                assert_ne!(read.d, read.d2);
            }
        }
        assert!(tandem > 30, "tandem DD should occur at include_dd = 0.3");
    }

    #[test]
    fn noise_is_confined_to_the_body_and_origins_survive_indels() {
        let (v, d, j) = sets();
        let params = SimParams {
            body_error_rate: GatedRate {
                p: 1.0,
                rate: crate::sim::random::Rate::Fixed(0.2),
            },
            indel_rate: GatedRate {
                p: 1.0,
                rate: crate::sim::random::Rate::Fixed(0.05),
            },
            flank_5p: Count::Fixed(20),
            flank_3p: Count::Fixed(20),
            min_length: 0,
            reverse_complement: 0.5,
            ..SimParams::train()
        };
        let sim = Simulator::new(&v, Some(&d), &j, params).unwrap();
        let mut rng = Rng::new(9);
        let (mut rc, mut indels) = (0, 0);
        for _ in 0..200 {
            let read = sim.read(&mut rng).unwrap();
            assert_eq!(read.span(Region::Flank5), Some(0..20));
            let flank = read.span(Region::Flank3).unwrap();
            assert_eq!(flank.len(), 20);
            // Origins of templated V bases strictly increase.
            let v_origins: Vec<u32> = read
                .regions
                .iter()
                .zip(&read.origins)
                .filter(|&(&r, &o)| r == Region::V && o != NO_ORIGIN)
                .map(|(_, &o)| o)
                .collect();
            assert!(v_origins.windows(2).all(|w| w[0] < w[1]));
            if read.reverse_complemented {
                rc += 1;
                assert_eq!(reverse_complement(&read.sequence), read.forward);
            }
            indels += read.indels;
        }
        assert!(rc > 60 && rc < 140);
        assert!(indels > 0);
    }

    #[test]
    fn constant_regions_continue_the_reading_frame_and_leave_vdj_untouched() {
        use crate::alphabet::encode;
        use crate::sim::{SyntheticDatabase, SyntheticSpec};
        use crate::translate::{STOP, translate};
        for (aligned, offset) in [(false, 0), (true, 1)] {
            let spec = SyntheticSpec {
                c_codon_aligned: aligned,
                ..SyntheticSpec::default()
            };
            let db = SyntheticDatabase::generate(&spec, 4);
            let params = SimParams {
                illumina_error: crate::sim::IlluminaError::NONE,
                reverse_complement: 0.5,
                ..SimParams::easy()
            };
            let plain = Simulator::new(&db.v, Some(&db.d), &db.j, params).unwrap();
            let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params)
                .unwrap()
                .with_constant(&db.c)
                .unwrap();
            let (without, with) = (plain.reads(300, 8).unwrap(), sim.reads(300, 8).unwrap());
            let mut constant = 0;
            for (a, b) in without.iter().zip(&with) {
                let flank = a.span(Region::Flank3).unwrap();
                assert_eq!(a.forward[..flank.start], b.forward[..flank.start]);
                assert_eq!((a.v, a.d, a.j, a.trims), (b.v, b.d, b.j, b.trims));
                let oriented = if b.reverse_complemented {
                    reverse_complement(&b.sequence)
                } else {
                    b.sequence.clone()
                };
                assert_eq!(oriented, b.forward);
                let Some(c) = b.c else {
                    assert!(b.span(Region::C).is_none());
                    continue;
                };
                constant += 1;
                let span = b.span(Region::C).unwrap();
                assert_eq!(span.start, flank.start);
                assert_eq!(sim.constant_offset(c, b.j), offset);
                assert_eq!(b.origins[span.start], offset as u32);
                let germline = &db.c.allele(c as u32).codes[offset..offset + span.len()];
                assert_eq!(encode(&b.forward[span.clone()]), germline);
                // The splice keeps the J's reading frame: no stop from the
                // J anchor through C.
                let anchor = b.locate(Region::J, db.j.allele(b.j as u32).anchor.unwrap());
                let spliced = encode(&b.forward[anchor.unwrap()..span.end]);
                let protein = translate(&spliced, 0);
                assert!(!protein[..protein.len() - 1].contains(&STOP), "{aligned}");
            }
            assert!((200..280).contains(&constant), "{constant}");
        }
    }
}
