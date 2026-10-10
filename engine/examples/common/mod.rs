//! Setup shared by the examples: command-line options and the germline
//! database (`--v/--d/--j` FASTA files, or a synthetic IGH-like database).

use std::error::Error;
use std::str::FromStr;

use vdj_assign::annotate::AnnotateOptions;
use vdj_assign::prelude::*;
use vdj_assign::prepare::{PreparedDatabase, prepare};
use vdj_assign::sim::{SyntheticDatabase, SyntheticSpec};

/// `--name value` options and `--flag`s, in any order.
pub struct Args(Vec<String>);

impl Args {
    pub fn from_env() -> Self {
        Self(std::env::args().skip(1).collect())
    }

    pub fn value(&self, name: &str) -> Option<String> {
        let at = self.0.iter().position(|a| a == name)?;
        self.0.get(at + 1).cloned()
    }

    #[allow(dead_code, reason = "not every example has flags")]
    pub fn flag(&self, name: &str) -> bool {
        self.0.iter().any(|a| a == name)
    }

    /// The parsed value of `name`, or `default` when it is absent.
    pub fn number<T>(&self, name: &str, default: T) -> Result<T, Box<dyn Error>>
    where
        T: FromStr,
        T::Err: Error + 'static,
    {
        Ok(self
            .value(name)
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(default))
    }

    /// The germline database of `--v`, `--d` and `--j`, or without them the
    /// synthetic database of `seed`; deduplicated and annotated.
    pub fn database(&self, seed: u64) -> Result<PreparedDatabase, Box<dyn Error>> {
        let load = |kind: SegmentKind, path: String| -> Result<GermlineSet, Box<dyn Error>> {
            let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
            Ok(GermlineSet::from_fasta(kind, &text).map_err(|e| format!("{path}: {e}"))?)
        };
        let (v, d, j) = match (self.value("--v"), self.value("--j")) {
            (Some(v), Some(j)) => (
                load(SegmentKind::V, v)?,
                self.value("--d")
                    .map(|d| load(SegmentKind::D, d))
                    .transpose()?,
                load(SegmentKind::J, j)?,
            ),
            _ => {
                let db = SyntheticDatabase::generate(&SyntheticSpec::default(), seed);
                (db.v, Some(db.d), db.j)
            }
        };
        let db = prepare(v, d, j, AnnotateOptions::default())?;
        eprintln!("database: {}", db.summary());
        Ok(db)
    }
}
