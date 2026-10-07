//! Helpers shared by the unit tests.

use crate::sim::Rng;

/// A deterministic stream of random words.
pub fn rng(seed: u64) -> impl FnMut() -> u64 {
    let mut rng = Rng::new(seed);
    move || rng.next_u64()
}

/// `len` random base codes from `next`.
pub fn random_codes(next: &mut impl FnMut() -> u64, len: usize) -> Vec<u8> {
    (0..len).map(|_| (next() % 4) as u8).collect()
}
