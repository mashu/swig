//! Database preparation: deduplicate, annotate, and optionally write out.
//!
//! [`prepare`] runs on every start (it takes milliseconds), so a database never
//! *has* to be prepared. Writing the result ([`PreparedDatabase::write`])
//! produces plain FASTA whose headers carry every anchor and region start,
//! plus `annotation.tsv` with the provenance of each value — for review,
//! manual curation, or sharing a fixed annotation. Header values always win
//! over detection, so edits to the written FASTA are respected on reload.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::annotate::{annotate, AlleleAnnotation, AnnotateOptions, DatabaseAnnotation};
use crate::germline::{GermlineError, GermlineSet, MergedDuplicate, SegmentKind};

fn label(annotations: &[AlleleAnnotation], id: u32) -> Option<String> {
    annotations.get(id as usize).map(|a| format!("annotation={}", a.source.label()))
}

pub struct PreparedDatabase {
    pub v: GermlineSet,
    pub d: Option<GermlineSet>,
    pub j: GermlineSet,
    pub merged: Vec<(SegmentKind, MergedDuplicate)>,
    pub annotation: DatabaseAnnotation,
}

pub fn prepare(
    mut v: GermlineSet,
    mut d: Option<GermlineSet>,
    mut j: GermlineSet,
    options: AnnotateOptions,
) -> Result<PreparedDatabase, GermlineError> {
    let mut merged = Vec::new();
    for set in [Some(&mut v), d.as_mut(), Some(&mut j)].into_iter().flatten() {
        let kind = set.kind();
        merged.extend(set.deduplicate()?.into_iter().map(|m| (kind, m)));
    }
    let annotation = annotate(&mut v, Some(&mut j), options);
    Ok(PreparedDatabase { v, d, j, merged, annotation })
}

impl PreparedDatabase {
    /// One line for logs.
    pub fn summary(&self) -> String {
        let count = |kind| self.merged.iter().filter(|(k, _)| *k == kind).count();
        let mut text = format!("V {} alleles", self.v.len());
        if let Some(d) = &self.d {
            text += &format!(", D {}", d.len());
        }
        text += &format!(", J {}", self.j.len());
        let merged: usize = SegmentKind::ALL.iter().map(|&k| count(k)).sum();
        if merged > 0 {
            let noun = if merged == 1 { "sequence" } else { "sequences" };
            text += &format!(" ({merged} identical {noun} merged as aliases)");
        }
        format!("{text}; {}", self.annotation.summary())
    }

    /// Write `V.fasta`, `D.fasta`, `J.fasta` and `annotation.tsv` into `dir`.
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
        let mut out = BufWriter::new(File::create(dir.join("annotation.tsv"))?);
        self.annotation.write_tsv(&mut out, &self.v, Some(&self.j))?;
        out.flush()
    }
}
