//! Minimal FASTA reader shared by germline loading and query input.

use std::fmt;

/// One FASTA record. `sequence` keeps the alphabetic symbols; IMGT gap
/// characters (`.`) are recorded in `gaps`, and whitespace, digits and every
/// other symbol are dropped.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastaRecord {
    pub header: String,
    pub sequence: Vec<u8>,
    /// Runs of `.` as `(bases before the run, run length)`, in sequence
    /// order; a run continues across line breaks. Empty (no allocation) for
    /// an ungapped record. Germline loading reads IMGT numbering from them
    /// ([`crate::germline::ImgtNumbering`]).
    pub gaps: Vec<(u32, u32)>,
}

impl FastaRecord {
    /// First whitespace-delimited token of the header.
    pub fn id(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or("")
    }

    /// `key=value` tokens after the id, e.g. `anchor=309`.
    pub fn attribute(&self, key: &str) -> Option<&str> {
        self.header.split_whitespace().skip(1).find_map(|token| {
            let (k, v) = token.split_once('=')?;
            (k == key).then_some(v)
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FastaError {
    SequenceBeforeHeader { line: usize },
    EmptyRecord { header: String },
}

impl fmt::Display for FastaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SequenceBeforeHeader { line } => {
                write!(f, "sequence data before the first header (line {line})")
            }
            Self::EmptyRecord { header } => write!(f, "record '{header}' has no sequence"),
        }
    }
}

impl std::error::Error for FastaError {}

/// Parse FASTA text.
pub fn parse_fasta(text: &str) -> Result<Vec<FastaRecord>, FastaError> {
    let mut records: Vec<FastaRecord> = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if let Some(header) = line.strip_prefix('>') {
            if let Some(last) = records.last()
                && last.sequence.is_empty()
            {
                return Err(FastaError::EmptyRecord {
                    header: last.header.clone(),
                });
            }
            records.push(FastaRecord {
                header: header.trim().to_string(),
                sequence: Vec::new(),
                gaps: Vec::new(),
            });
        } else if !line.is_empty() {
            let record = records
                .last_mut()
                .ok_or(FastaError::SequenceBeforeHeader { line: index + 1 })?;
            for symbol in line.bytes() {
                if symbol.is_ascii_alphabetic() {
                    record.sequence.push(symbol);
                } else if symbol == b'.' {
                    let before = u32::try_from(record.sequence.len()).expect("record fits u32");
                    match record.gaps.last_mut() {
                        Some((at, run)) if *at == before => *run += 1,
                        _ => record.gaps.push((before, 1)),
                    }
                }
            }
        }
    }
    if let Some(last) = records.last()
        && last.sequence.is_empty()
    {
        return Err(FastaError::EmptyRecord {
            header: last.header.clone(),
        });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_records_attributes_and_gaps() {
        let text = ">IGHV1-2*02 anchor=285\n..CAG..GTG.\n..CA-G\n>IGHJ4*02\nACTAC\n";
        let records = parse_fasta(text).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), "IGHV1-2*02");
        assert_eq!(records[0].attribute("anchor"), Some("285"));
        assert_eq!(records[0].sequence, b"CAGGTGCAG");
        // Runs are (bases before, length); one continues across the line break.
        assert_eq!(records[0].gaps, [(0, 2), (3, 2), (6, 3)]);
        assert_eq!(records[1].sequence, b"ACTAC");
        assert!(records[1].gaps.is_empty());
    }

    #[test]
    fn rejects_sequence_before_header() {
        assert!(matches!(
            parse_fasta("ACGT\n>x\nA"),
            Err(FastaError::SequenceBeforeHeader { line: 1 })
        ));
    }
}
