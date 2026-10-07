//! Read simulation with exact ground truth, and evaluation of calls against it.
//!
//! A dependency-free Rust port of [IgSim.jl](https://github.com/mashu/IgSim.jl)
//! (same generative model and presets), extended with per-base germline
//! origins and optional reverse complementation, plus:
//!
//! * [`synthetic`] — IGH-like germline databases with a family → gene →
//!   allele similarity hierarchy, for regression tests that cannot ship
//!   IMGT-derived data;
//! * [`truth`] — serialisable ground truth per read;
//! * [`evaluate`] — allele, gene, boundary, D-presence and junction accuracy
//!   of any tool's AIRR calls against the truth.
//!
//! ```ignore
//! use vdj_assign::sim::{SimParams, Simulator, Evaluation};
//! let sim = Simulator::new(&v, Some(&d), &j, SimParams::hard().with_reverse_complement(0.5))?;
//! let reads = sim.reads(10_000, 1)?;
//! let mut evaluation = Evaluation::default();
//! for (i, (read, record)) in reads.iter().zip(&records).enumerate() {
//!     let truth = Truth::from_read(&sim, read, &format!("r{i}"));
//!     evaluation.add(&truth, &Assignment::from_record(record));
//! }
//! println!("{evaluation}");
//! ```

pub mod evaluate;
pub mod params;
pub mod random;
pub mod simulate;
pub mod synthetic;
pub mod truth;

pub use evaluate::{Assignment, DPresence, Evaluation, SegmentCall, SegmentStats};
pub use params::{IlluminaError, SimParams};
pub use random::{Count, GatedRate, Rate, Rng};
pub use simulate::{
    NO_ORIGIN, Region, SimError, SimulatedRead, Simulator, Trims, reverse_complement,
};
pub use synthetic::{SyntheticDatabase, SyntheticSpec};
pub use truth::{Slot, TrueSegment, Truth};
