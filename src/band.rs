//! Diagonal bands for the dynamic-programming matrix.
//!
//! Rows are query positions `i = 1..=n`, columns reference positions
//! `j = 1..=m`, and a diagonal is `d = i - j` (the seed offset of the method
//! note). A band keeps only cells whose diagonal lies in `[lo, hi]`; the note's
//! `|j - (i - d̂)| ≤ w` is `Band::around(d̂, w)`.
//!
//! Invariant relied on by [`crate::align::BandedGotoh`]: for consecutive rows the
//! first and last in-band columns never decrease and grow by at most one. Both
//! variants satisfy it, which is what lets the aligner keep a single rolling
//! row of `O(m)` memory and store traceback for band cells only.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Band {
    /// Unrestricted matrix (seedless / exhaustive candidates, short D windows).
    Full,
    /// Inclusive diagonal range `lo ..= hi` of `d = i - j`.
    Diagonals { lo: i32, hi: i32 },
}

impl Band {
    /// Symmetric band of half-width `half_width` around diagonal `center`.
    pub fn around(center: i32, half_width: i32) -> Self {
        Band::Diagonals { lo: center - half_width, hi: center + half_width }
    }

    /// Inclusive 1-based column range of 1-based `row`, or `None` if empty.
    #[inline]
    pub fn columns(&self, row: usize, ref_len: usize) -> Option<(usize, usize)> {
        match *self {
            Band::Full => (ref_len > 0).then_some((1, ref_len)),
            Band::Diagonals { lo, hi } => {
                let i = row as i64;
                let first = (i - i64::from(hi)).max(1);
                let last = (i - i64::from(lo)).min(ref_len as i64);
                (first <= last).then(|| (first as usize, last as usize))
            }
        }
    }

    /// `true` once `row` and every later row lie entirely right of the matrix.
    #[inline]
    pub fn exhausted(&self, row: usize, ref_len: usize) -> bool {
        match *self {
            Band::Full => false,
            Band::Diagonals { hi, .. } => row as i64 - i64::from(hi) > ref_len as i64,
        }
    }

    /// Whether a 1-based cell lies inside the band.
    pub fn contains(&self, row: usize, col: usize) -> bool {
        match *self {
            Band::Full => true,
            Band::Diagonals { lo, hi } => {
                let d = row as i64 - col as i64;
                (i64::from(lo)..=i64::from(hi)).contains(&d)
            }
        }
    }

    /// Number of cells evaluated for an `rows × ref_len` matrix (work estimate).
    pub fn cells(&self, rows: usize, ref_len: usize) -> usize {
        (1..=rows)
            .filter_map(|i| self.columns(i, ref_len))
            .map(|(a, b)| b - a + 1)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_band_rows_are_monotone() {
        let band = Band::around(3, 2);
        let mut previous: Option<(usize, usize)> = None;
        for i in 1..=40 {
            if let Some((a, b)) = band.columns(i, 30) {
                assert!(band.contains(i, a) && band.contains(i, b));
                if let Some((pa, pb)) = previous {
                    assert!(a >= pa && a <= pa + 1 && b >= pb && b <= pb + 1);
                }
                previous = Some((a, b));
            }
        }
        assert!(band.exhausted(40, 30));
        assert!(band.cells(40, 30) < 40 * 30);
    }
}
