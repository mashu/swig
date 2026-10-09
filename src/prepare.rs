//! Database preparation: deduplicate, annotate, and optionally write out.
//!
//! [`prepare`] runs on every start (it takes milliseconds), so a database never
//! *has* to be prepared. Writing the result ([`PreparedDatabase::write`])
//! produces plain FASTA whose headers carry every anchor and region start,
//! plus `annotation.tsv` with the provenance of each value — for review,
//! manual curation, or sharing a fixed annotation. Header values always win
//! over detection, so edits to the written FASTA are respected on reload.

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::annotate::{AlleleAnnotation, AnnotateOptions, DatabaseAnnotation, annotate};
use crate::germline::{GermlineError, GermlineSet, MergedDuplicate, SegmentKind};
use crate::retrieve::GeneFrames;
use crate::seed::SeedShape;
use crate::translate::open_frames;

/// Longest plausible coding constant region: human C exons total at most
/// 1.6 kb, so a longer record carries introns or UTR. Its anchored jobs
/// would also leave the 16-bit SIMD lanes (anchor bonus `2·|C|` plus up to
/// `2` per column).
pub const MAX_CODING_C: usize = 2000;

fn label(annotations: &[AlleleAnnotation], id: u32) -> Option<String> {
    annotations
        .get(id as usize)
        .map(|a| format!("annotation={}", a.source.label()))
}

pub struct PreparedDatabase {
    pub v: GermlineSet,
    pub d: Option<GermlineSet>,
    pub j: GermlineSet,
    /// Constant regions, when given ([`PreparedDatabase::with_constant`]).
    pub c: Option<GermlineSet>,
    pub merged: Vec<(SegmentKind, MergedDuplicate)>,
    pub annotation: DatabaseAnnotation,
}

/// A database convention that keeps C from being called: the C stage
/// anchors C's base 0 within `s` bases of the J's germline 3′ end
/// ([`crate::constant`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConstantWarning {
    /// A C allele lacks the first `by > s` bases of its gene's longest
    /// allele, so its base 0 never lies next to a J exit.
    Truncated { allele: String, by: u32 },
    /// A C record longer than [`MAX_CODING_C`] or without a stop-free reading
    /// frame: a genomic record (introns, UTR), not the spliced coding region.
    NotCoding { allele: String, len: usize },
    /// A J allele whose 3′ end lies `tail` bases past its anchor where the
    /// database's mode is `mode`, `|tail − mode| > s`: bases past the splice
    /// donor move the J exit off the C splice.
    JEnd {
        allele: String,
        tail: u32,
        mode: u32,
    },
}

impl fmt::Display for ConstantWarning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Truncated { allele, by } => write!(
                f,
                "C allele {allele} lacks the first {by} bases of its gene's longest allele, more than the splice spread: it is never called"
            ),
            Self::NotCoding { allele, len } => write!(
                f,
                "C allele {allele} ({len} nt) is longer than {MAX_CODING_C} nt or has no stop-free reading frame: a genomic record (introns, UTR) rather than the spliced coding region, so its base 0 need not follow the J"
            ),
            Self::JEnd { allele, tail, mode } => write!(
                f,
                "J allele {allele} ends {tail} bases after its anchor where most end {mode} after: bases past the splice donor move the J exit off the C splice, so reads with it get no C call"
            ),
        }
    }
}

pub fn prepare(
    mut v: GermlineSet,
    mut d: Option<GermlineSet>,
    mut j: GermlineSet,
    options: AnnotateOptions,
) -> Result<PreparedDatabase, GermlineError> {
    let mut merged = Vec::new();
    for set in [Some(&mut v), d.as_mut(), Some(&mut j)]
        .into_iter()
        .flatten()
    {
        let kind = set.kind();
        merged.extend(set.deduplicate()?.into_iter().map(|m| (kind, m)));
    }
    let annotation = annotate(&mut v, Some(&mut j), options);
    Ok(PreparedDatabase {
        v,
        d,
        j,
        c: None,
        merged,
        annotation,
    })
}

impl PreparedDatabase {
    /// Add a constant-region database, identical sequences merged as for V,
    /// D and J.
    pub fn with_constant(mut self, mut c: GermlineSet) -> Result<Self, GermlineError> {
        let kind = c.kind();
        self.merged
            .extend(c.deduplicate()?.into_iter().map(|m| (kind, m)));
        self.c = Some(c);
        Ok(self)
    }

    /// The conventions of the database that keep C from being called at
    /// splice spread `spread` (none without C).
    pub fn constant_warnings(&self, spread: u32) -> Vec<ConstantWarning> {
        let Some(c) = &self.c else {
            return Vec::new();
        };
        let mut warnings = Vec::new();
        // Frame offsets: how far into its gene's longest allele each starts.
        let frames = GeneFrames::build(c, &SeedShape::contiguous(9));
        for (id, allele) in (0..).zip(c.alleles()) {
            let by = frames.offset(id);
            if by > spread {
                let allele = allele.name.clone();
                warnings.push(ConstantWarning::Truncated { allele, by });
            }
            if allele.codes.len() > MAX_CODING_C || open_frames(&allele.codes).next().is_none() {
                warnings.push(ConstantWarning::NotCoding {
                    allele: allele.name.clone(),
                    len: allele.codes.len(),
                });
            }
        }
        let tail = |codes: &[u8], anchor: u32| (codes.len() as u32).saturating_sub(anchor);
        let mut counts: HashMap<u32, usize> = HashMap::new();
        for allele in self.j.alleles() {
            if let Some(anchor) = allele.anchor {
                *counts.entry(tail(&allele.codes, anchor)).or_default() += 1;
            }
        }
        // The most frequent tail, the shortest among equals.
        if let Some((&mode, _)) = counts
            .iter()
            .max_by_key(|&(&t, &n)| (n, std::cmp::Reverse(t)))
        {
            for allele in self.j.alleles() {
                let Some(anchor) = allele.anchor else {
                    continue;
                };
                let t = tail(&allele.codes, anchor);
                if t.abs_diff(mode) > spread {
                    warnings.push(ConstantWarning::JEnd {
                        allele: allele.name.clone(),
                        tail: t,
                        mode,
                    });
                }
            }
        }
        warnings
    }

    /// One line for logs.
    pub fn summary(&self) -> String {
        let count = |kind| self.merged.iter().filter(|(k, _)| *k == kind).count();
        let mut text = format!("V {} alleles", self.v.len());
        if let Some(d) = &self.d {
            text += &format!(", D {}", d.len());
        }
        text += &format!(", J {}", self.j.len());
        if let Some(c) = &self.c {
            text += &format!(", C {}", c.len());
        }
        let merged: usize = SegmentKind::ALL.iter().map(|&k| count(k)).sum();
        if merged > 0 {
            let noun = if merged == 1 { "sequence" } else { "sequences" };
            text += &format!(" ({merged} identical {noun} merged as aliases)");
        }
        format!("{text}; {}", self.annotation.summary())
    }

    /// Write `V.fasta`, `D.fasta`, `J.fasta`, `C.fasta` and `annotation.tsv`
    /// into `dir` (D and C when present).
    pub fn write(&self, dir: &Path) -> io::Result<()> {
        std::fs::create_dir_all(dir)?;
        let v_labels = &self.annotation.v;
        let j_labels = &self.annotation.j;
        let mut out = BufWriter::new(File::create(dir.join("V.fasta"))?);
        self.v.write_fasta(&mut out, |id| label(v_labels, id))?;
        out.flush()?;
        if let Some(d) = &self.d {
            let mut out = BufWriter::new(File::create(dir.join("D.fasta"))?);
            d.write_fasta(&mut out, |_| None)?;
            out.flush()?;
        }
        let mut out = BufWriter::new(File::create(dir.join("J.fasta"))?);
        self.j.write_fasta(&mut out, |id| label(j_labels, id))?;
        out.flush()?;
        if let Some(c) = &self.c {
            let mut out = BufWriter::new(File::create(dir.join("C.fasta"))?);
            c.write_fasta(&mut out, |_| None)?;
            out.flush()?;
        }
        let mut out = BufWriter::new(File::create(dir.join("annotation.tsv"))?);
        self.annotation
            .write_tsv(&mut out, &self.v, Some(&self.j))?;
        out.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::decode;
    use crate::testing::{random_codes, rng};

    /// `codons` random sense codons and a stop: open in frame 0.
    fn coding(next: &mut impl FnMut() -> u64, codons: usize) -> Vec<u8> {
        let mut text = String::new();
        while text.len() < 3 * codons {
            let codon = decode(&random_codes(next, 3));
            if !matches!(codon.as_str(), "TAA" | "TAG" | "TGA") {
                text += &codon;
            }
        }
        (text + "TGA").into_bytes()
    }

    #[test]
    fn constant_warnings_name_each_convention_that_blocks_c() {
        let mut next = rng(11);
        let mut v = GermlineSet::new(SegmentKind::V);
        let germline = decode(&random_codes(&mut next, 300));
        v.push("V1*01", germline.as_bytes(), None).unwrap();
        // J tails (bases after the anchor) 11, 11, 13 and 15: mode 11, and
        // only J4 strays by more than s = 3.
        let mut j = GermlineSet::new(SegmentKind::J);
        for (name, tail) in [("J1*01", 11), ("J2*01", 11), ("J3*01", 13), ("J4*01", 15)] {
            let sequence = decode(&random_codes(&mut next, 6 + tail));
            j.push(name, sequence.as_bytes(), Some(6)).unwrap();
        }
        let full = coding(&mut next, 110);
        let mut c = GermlineSet::new(SegmentKind::C);
        c.push("G1*01", &full, None).unwrap();
        c.push("G1*02", &full[3..], None).unwrap(); // within s: fine
        c.push("G1*03", &full[5..], None).unwrap(); // 5′-truncated
        c.push("G2*01", &coding(&mut next, 700), None).unwrap(); // > 2 kb
        let genomic = decode(&random_codes(&mut next, 400)); // stops in every frame
        c.push("G3*01", genomic.as_bytes(), None).unwrap();
        let database = prepare(v, None, j, AnnotateOptions::default())
            .unwrap()
            .with_constant(c)
            .unwrap();
        let allele = |name: &str| name.to_string();
        assert_eq!(
            database.constant_warnings(3),
            vec![
                ConstantWarning::Truncated {
                    allele: allele("G1*03"),
                    by: 5
                },
                ConstantWarning::NotCoding {
                    allele: allele("G2*01"),
                    len: 2103
                },
                ConstantWarning::NotCoding {
                    allele: allele("G3*01"),
                    len: 400
                },
                ConstantWarning::JEnd {
                    allele: allele("J4*01"),
                    tail: 15,
                    mode: 11
                },
            ]
        );
        // At s = 4, J4 is within reach.
        assert_eq!(database.constant_warnings(4).len(), 3);
    }
}
