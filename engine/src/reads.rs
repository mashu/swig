//! Streaming FASTA / FASTQ / AIRR-table reader for query reads.
//!
//! The format is detected from the first non-blank character (`>` FASTA,
//! `@` FASTQ, anything else an AIRR Rearrangement table); reads are
//! delivered in batches so arbitrarily large files run in bounded memory.
//! FASTQ records may span several sequence and quality lines (quality lines
//! may start with `@` or `+`; records are delimited by counting quality
//! characters); qualities are validated and kept for the AIRR `quality`
//! column. An AIRR table is tab-separated, or comma-separated when its
//! header line holds no tab, and is read from its `sequence`, `sequence_id`
//! and `quality` columns: Swig re-annotates AIRR input from its sequences.
//! Lines end in `\n` or `\r\n`; a carriage return inside a line (CR-only
//! line endings) is an error rather than one long line. Gzip input is
//! recognised by its magic bytes and rejected with a hint to decompress on
//! the fly (the crate has no dependencies).

use std::fmt;
use std::io::{self, BufRead};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadFormat {
    Fasta,
    Fastq,
    /// An AIRR Rearrangement table (TSV, or CSV without quoting).
    Airr,
}

/// One read: its identifier, raw sequence symbols and FASTQ qualities.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadRecord {
    /// The first whitespace-delimited token of a FASTA/FASTQ header, or an
    /// AIRR row's `sequence_id`; may be empty.
    pub id: String,
    pub sequence: Vec<u8>,
    /// FASTQ qualities or an AIRR row's `quality`; empty for FASTA.
    pub quality: String,
}

/// The identifier of a FASTA/FASTQ header line (without `>`/`@`).
fn header_id(header: &str) -> String {
    header.split_whitespace().next().unwrap_or("").to_string()
}

/// Where an AIRR table keeps what a read needs.
#[derive(Clone, Copy, Debug)]
struct AirrColumns {
    delimiter: char,
    sequence: usize,
    id: Option<usize>,
    quality: Option<usize>,
}

#[derive(Debug)]
pub enum ReadError {
    Io(io::Error),
    Gzip,
    UnknownFormat {
        line: usize,
        found: char,
    },
    SequenceBeforeHeader {
        line: usize,
    },
    ExpectedFastqHeader {
        line: usize,
    },
    MissingSeparator {
        header: String,
    },
    /// A FASTQ sequence line holds whitespace or a control character, which
    /// would split an output column.
    InvalidSymbol {
        line: usize,
        found: char,
    },
    QualityLength {
        header: String,
        sequence: usize,
        quality: usize,
    },
    /// A table (neither FASTA nor FASTQ) without a `sequence` column.
    NoSequenceColumn {
        line: usize,
    },
    /// A carriage return inside a line: CR-only (classic Mac) line endings,
    /// which would otherwise read the whole input as one line.
    CarriageReturn {
        line: usize,
    },
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::Gzip => write!(
                f,
                "input is gzip-compressed; decompress on the fly, e.g. `zcat reads.fastq.gz | swig ... -`"
            ),
            Self::UnknownFormat { line, found } => {
                write!(
                    f,
                    "line {line}: expected FASTA ('>'), FASTQ ('@') or an AIRR table, found {found:?}"
                )
            }
            Self::NoSequenceColumn { line } => write!(
                f,
                "line {line}: expected FASTA ('>'), FASTQ ('@') or an AIRR table with a 'sequence' column"
            ),
            Self::SequenceBeforeHeader { line } => {
                write!(f, "sequence data before the first header (line {line})")
            }
            Self::ExpectedFastqHeader { line } => {
                write!(f, "line {line}: expected a FASTQ header starting with '@'")
            }
            Self::MissingSeparator { header } => {
                write!(f, "FASTQ record '{header}' has no '+' separator line")
            }
            Self::InvalidSymbol { line, found } => {
                write!(f, "line {line}: {found:?} in a FASTQ sequence line")
            }
            Self::QualityLength {
                header,
                sequence,
                quality,
            } => write!(
                f,
                "FASTQ record '{header}': {quality} quality characters for {sequence} bases"
            ),
            Self::CarriageReturn { line } => write!(
                f,
                "line {line}: carriage return inside a line (CR-only line endings?); convert them, e.g. `tr '\\r' '\\n'`"
            ),
        }
    }
}

impl std::error::Error for ReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

pub struct ReadReader<R> {
    input: R,
    format: Option<ReadFormat>,
    started: bool,
    line: usize,
    buffer: String,
    /// FASTA header already read while finishing the previous record.
    pending_header: Option<String>,
    /// The columns of an AIRR table.
    airr: Option<AirrColumns>,
}

impl<R: BufRead> ReadReader<R> {
    pub fn new(input: R) -> Self {
        Self {
            input,
            format: None,
            started: false,
            line: 0,
            buffer: String::new(),
            pending_header: None,
            airr: None,
        }
    }

    /// The underlying input, for instance to see how far a reader over an
    /// in-memory cursor has got.
    pub fn get_ref(&self) -> &R {
        &self.input
    }

    /// Detected format, `None` for empty input.
    pub fn format(&mut self) -> Result<Option<ReadFormat>, ReadError> {
        if !self.started {
            self.started = true;
            self.format = self.detect()?;
        }
        Ok(self.format)
    }

    /// Append up to `max` records to `out`; returns how many were appended
    /// (0 at end of input).
    pub fn read_batch(
        &mut self,
        max: usize,
        out: &mut Vec<ReadRecord>,
    ) -> Result<usize, ReadError> {
        let mut count = 0;
        while count < max {
            match self.next_record()? {
                Some(record) => {
                    out.push(record);
                    count += 1;
                }
                None => break,
            }
        }
        Ok(count)
    }

    pub fn next_record(&mut self) -> Result<Option<ReadRecord>, ReadError> {
        match self.format()? {
            None => Ok(None),
            Some(ReadFormat::Fasta) => self.next_fasta(),
            Some(ReadFormat::Fastq) => self.next_fastq(),
            Some(ReadFormat::Airr) => self.next_airr(),
        }
    }

    /// Skip leading blank space, reject gzip, and classify the first symbol.
    fn detect(&mut self) -> Result<Option<ReadFormat>, ReadError> {
        let mut first_chunk = true;
        loop {
            let chunk = self.input.fill_buf()?;
            if chunk.is_empty() {
                return Ok(None);
            }
            if first_chunk && chunk.starts_with(&[0x1f, 0x8b]) {
                return Err(ReadError::Gzip);
            }
            first_chunk = false;
            let blank = chunk.iter().take_while(|b| b.is_ascii_whitespace()).count();
            let newlines = chunk[..blank].iter().filter(|&&b| b == b'\n').count();
            if blank == chunk.len() {
                self.line += newlines;
                self.input.consume(blank);
                continue;
            }
            let symbol = chunk[blank];
            self.line += newlines;
            self.input.consume(blank);
            return match symbol {
                b'>' => Ok(Some(ReadFormat::Fasta)),
                b'@' => Ok(Some(ReadFormat::Fastq)),
                other if other.is_ascii_graphic() => {
                    self.airr = Some(self.airr_header()?);
                    Ok(Some(ReadFormat::Airr))
                }
                other => Err(ReadError::UnknownFormat {
                    line: self.line + 1,
                    found: char::from(other),
                }),
            };
        }
    }

    /// Read an AIRR table's header line and find its columns.
    fn airr_header(&mut self) -> Result<AirrColumns, ReadError> {
        self.next_line()?;
        let delimiter = if self.buffer.contains('\t') {
            '\t'
        } else {
            ','
        };
        let names: Vec<&str> = self.buffer.split(delimiter).map(str::trim).collect();
        let column = |name: &str| names.iter().position(|&n| n == name);
        Ok(AirrColumns {
            delimiter,
            sequence: column("sequence").ok_or(ReadError::NoSequenceColumn { line: self.line })?,
            id: column("sequence_id"),
            quality: column("quality"),
        })
    }

    /// The next AIRR row with a sequence; rows without one are skipped.
    fn next_airr(&mut self) -> Result<Option<ReadRecord>, ReadError> {
        let columns = self.airr.expect("the header was read on detection");
        while self.next_line()? {
            let values: Vec<&str> = self.buffer.split(columns.delimiter).collect();
            let value = |at: Option<usize>| at.and_then(|i| values.get(i)).map_or("", |v| v.trim());
            let sequence = value(Some(columns.sequence));
            if sequence.is_empty() {
                continue;
            }
            return Ok(Some(ReadRecord {
                id: value(columns.id).to_string(),
                sequence: sequence
                    .bytes()
                    .filter(|b| b.is_ascii_alphabetic())
                    .collect(),
                quality: value(columns.quality).to_string(),
            }));
        }
        Ok(None)
    }

    /// Read one line into `self.buffer` without its line terminator (`\n`
    /// or `\r\n`). A carriage return left inside the line is an error: no
    /// FASTA/FASTQ line or unquoted table field holds one, and with CR-only
    /// line endings the whole input would read as one line.
    fn next_line(&mut self) -> Result<bool, ReadError> {
        self.buffer.clear();
        if self.input.read_line(&mut self.buffer)? == 0 {
            return Ok(false);
        }
        self.line += 1;
        let trimmed = self.buffer.trim_end_matches(['\n', '\r']).len();
        self.buffer.truncate(trimmed);
        if self.buffer.contains('\r') {
            return Err(ReadError::CarriageReturn { line: self.line });
        }
        Ok(true)
    }

    fn next_fasta(&mut self) -> Result<Option<ReadRecord>, ReadError> {
        let header = match self.pending_header.take() {
            Some(header) => header,
            None => loop {
                if !self.next_line()? {
                    return Ok(None);
                }
                if let Some(header) = self.buffer.strip_prefix('>') {
                    break header.trim().to_string();
                }
                if !self.buffer.trim().is_empty() {
                    return Err(ReadError::SequenceBeforeHeader { line: self.line });
                }
            },
        };
        let mut sequence = Vec::new();
        while self.next_line()? {
            if let Some(next) = self.buffer.strip_prefix('>') {
                self.pending_header = Some(next.trim().to_string());
                break;
            }
            sequence.extend(self.buffer.bytes().filter(|b| b.is_ascii_alphabetic()));
        }
        Ok(Some(ReadRecord {
            id: header_id(&header),
            sequence,
            quality: String::new(),
        }))
    }

    fn next_fastq(&mut self) -> Result<Option<ReadRecord>, ReadError> {
        let header = loop {
            if !self.next_line()? {
                return Ok(None);
            }
            if let Some(header) = self.buffer.strip_prefix('@') {
                break header.trim().to_string();
            }
            if !self.buffer.trim().is_empty() {
                return Err(ReadError::ExpectedFastqHeader { line: self.line });
            }
        };
        let mut sequence = Vec::new();
        loop {
            if !self.next_line()? {
                return Err(ReadError::MissingSeparator { header });
            }
            if self.buffer.starts_with('+') {
                break;
            }
            let line = self.buffer.trim();
            // Unlike FASTA, nothing can be dropped: qualities count the
            // bases. No-call symbols such as `.` and `-` read as N.
            if let Some(found) = line.chars().find(|c| c.is_whitespace() || c.is_control()) {
                return Err(ReadError::InvalidSymbol {
                    line: self.line,
                    found,
                });
            }
            sequence.extend(line.bytes());
        }
        let mut quality = String::new();
        while quality.len() < sequence.len() {
            if !self.next_line()? {
                break;
            }
            quality.push_str(self.buffer.trim());
        }
        if quality.len() != sequence.len() {
            return Err(ReadError::QualityLength {
                header,
                sequence: sequence.len(),
                quality: quality.len(),
            });
        }
        Ok(Some(ReadRecord {
            id: header_id(&header),
            sequence,
            quality,
        }))
    }
}

impl<R: BufRead> Iterator for ReadReader<R> {
    type Item = Result<ReadRecord, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_record().transpose()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_all(text: &[u8]) -> Result<Vec<ReadRecord>, ReadError> {
        ReadReader::new(text).collect()
    }

    fn record(id: &str, sequence: &str) -> ReadRecord {
        ReadRecord {
            id: id.into(),
            sequence: sequence.as_bytes().to_vec(),
            quality: String::new(),
        }
    }

    fn fastq(id: &str, sequence: &str, quality: &str) -> ReadRecord {
        ReadRecord {
            quality: quality.into(),
            ..record(id, sequence)
        }
    }

    #[test]
    fn reads_multiline_fasta() {
        let records = read_all(b"\n>r1 desc\nACGT\nac-gt\n>r2\r\nTTTT\r\n").unwrap();
        assert_eq!(
            records,
            vec![record("r1", "ACGTacgt"), record("r2", "TTTT")]
        );
    }

    #[test]
    fn reads_fastq_including_multiline_and_awkward_qualities() {
        let text = b"@r1\nACGT\n+\nIIII\n\n@r2 x\nAC\nGT\n+r2 x\n@+II\n\n@r3\n\n+\n\n";
        let records = read_all(text).unwrap();
        assert_eq!(
            records,
            vec![
                fastq("r1", "ACGT", "IIII"),
                fastq("r2", "ACGT", "@+II"),
                fastq("r3", "", "")
            ]
        );
    }

    #[test]
    fn detects_format_and_reports_problems() {
        let mut reader = ReadReader::new(&b"@r\nAC\n+\nII\n"[..]);
        assert_eq!(reader.format().unwrap(), Some(ReadFormat::Fastq));
        assert!(matches!(
            read_all(b"@r\nACGT\n+\nII\n"),
            Err(ReadError::QualityLength {
                sequence: 4,
                quality: 2,
                ..
            })
        ));
        assert!(matches!(
            read_all(b"@r\nACGT\n"),
            Err(ReadError::MissingSeparator { .. })
        ));
        assert!(matches!(
            read_all(b"@r\nAC\tG\n+\nIIII\n"),
            Err(ReadError::InvalidSymbol {
                line: 2,
                found: '\t'
            })
        ));
        assert_eq!(
            read_all(b"@r\nAC.T-\n+\nIIIII\n").unwrap(),
            [fastq("r", "AC.T-", "IIIII")]
        );
        assert!(matches!(
            read_all(b"ACGT\n"),
            Err(ReadError::NoSequenceColumn { line: 1 })
        ));
        assert!(matches!(
            read_all(b"\x01ACGT\n"),
            Err(ReadError::UnknownFormat { line: 1, .. })
        ));
        assert!(matches!(
            read_all(&[0x1f, 0x8b, 8, 0]),
            Err(ReadError::Gzip)
        ));
        assert!(read_all(b"").unwrap().is_empty());
        assert!(read_all(b"\n  \n").unwrap().is_empty());
    }

    #[test]
    fn rejects_cr_only_line_endings() {
        // Classic Mac line endings would read the whole file as one line:
        // one FASTA read with an empty sequence, or no AIRR row at all.
        for text in [
            &b">r1\rACGT\r>r2\rTTTT\r"[..],
            b"@r1\rACGT\r+\rIIII\r",
            b"sequence_id\tsequence\tv_call\rq1\tACGT\t\r",
            b"\r\r>r1\rACGT\r",
        ] {
            assert!(
                matches!(read_all(text), Err(ReadError::CarriageReturn { line: 1 })),
                "{:?}",
                String::from_utf8_lossy(text)
            );
        }
        assert!(matches!(
            read_all(b">r1\nAC\rGT\n"),
            Err(ReadError::CarriageReturn { line: 2 })
        ));
    }

    #[test]
    fn batches_resume_where_they_stopped() {
        let text = b">a\nA\n>b\nC\n>c\nG\n";
        let mut reader = ReadReader::new(&text[..]);
        let mut out = Vec::new();
        assert_eq!(reader.read_batch(2, &mut out).unwrap(), 2);
        assert_eq!(reader.read_batch(2, &mut out).unwrap(), 1);
        assert_eq!(reader.read_batch(2, &mut out).unwrap(), 0);
        assert_eq!(
            out.iter().map(|r| r.id.as_str()).collect::<Vec<_>>(),
            ["a", "b", "c"]
        );
    }

    #[test]
    fn reads_airr_tables_by_column_name() {
        let tsv = b"sequence_id\tv_call\tsequence\tquality\nq1\tIGHV1\tACGT\tIIII\n\nq2\t\t\t\nq3\t\tAC-gt\t\r\n";
        let records = read_all(tsv).unwrap();
        assert_eq!(
            records,
            vec![fastq("q1", "ACGT", "IIII"), record("q3", "ACgt")]
        );
        let csv = b"sequence,sequence_id\nTTTT,c1\nGG,\n";
        assert_eq!(
            read_all(csv).unwrap(),
            vec![record("c1", "TTTT"), record("", "GG")]
        );
        let mut reader = ReadReader::new(&csv[..]);
        assert_eq!(reader.format().unwrap(), Some(ReadFormat::Airr));
    }
}
