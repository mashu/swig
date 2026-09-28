//! Minimal FASTA reader shared by germline loading and query input.

use std::fmt;

/// One FASTA record. `sequence` keeps the raw symbols (IMGT gap characters,
/// whitespace and digits removed).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastaRecord {
    pub header: String,
    pub sequence: Vec<u8>,
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
            if let Some(last) = records.last() {
                if last.sequence.is_empty() {
                    return Err(FastaError::EmptyRecord {
                        header: last.header.clone(),
                    });
                }
            }
            records.push(FastaRecord {
                header: header.trim().to_string(),
                sequence: Vec::new(),
            });
        } else if !line.is_empty() {
            let record = records
                .last_mut()
                .ok_or(FastaError::SequenceBeforeHeader { line: index + 1 })?;
            record
                .sequence
                .extend(line.bytes().filter(|b| b.is_ascii_alphabetic()));
        }
    }
    if let Some(last) = records.last() {
        if last.sequence.is_empty() {
            return Err(FastaError::EmptyRecord {
                header: last.header.clone(),
            });
        }
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_records_attributes_and_gaps() {
        let text = ">IGHV1-2*02 anchor=285\nCAG..GTG\nCAG\n>IGHJ4*02\nACTAC\n";
        let records = parse_fasta(text).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), "IGHV1-2*02");
        assert_eq!(records[0].attribute("anchor"), Some("285"));
        assert_eq!(records[0].sequence, b"CAGGTGCAG");
        assert_eq!(records[1].sequence, b"ACTAC");
    }

    #[test]
    fn rejects_sequence_before_header() {
        assert!(matches!(
            parse_fasta("ACGT\n>x\nA"),
            Err(FastaError::SequenceBeforeHeader { line: 1 })
        ));
    }
}
