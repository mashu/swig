//! Exhaustive retrieval: every allele is a candidate, aligned unbanded.
//!
//! For a small database (J: 10–20 alleles) retrieval cannot save work worth
//! its risk — under heavy SHM a 50-nt J keeps few intact seeds — so every
//! allele is aligned. The query window (right of V, see
//! [`crate::window::DownstreamWindow`]) keeps the matrices small.

use super::{Candidate, RetrievalScratch, Retriever};
use crate::band::Band;
use crate::germline::AlleleId;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExhaustiveRetriever {
    pub alleles: usize,
}

impl Retriever for ExhaustiveRetriever {
    fn retrieve(&self, _: &[u8], _: &mut RetrievalScratch, out: &mut Vec<Candidate>) {
        out.clear();
        out.extend((0..self.alleles).map(|allele| Candidate {
            allele: allele as AlleleId,
            band: Band::Full,
            votes: 0,
        }));
    }
}
