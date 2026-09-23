//! Streaming FASTA / FASTQ reader for query reads.
//!
//! The format is detected from the first non-blank character (`>` FASTA,
//! `@` FASTQ); reads are delivered in batches so arbitrarily large files run
//! in bounded memory. FASTQ records may span several sequence and quality
//! lines (quality lines may start with `@` or `+`; records are delimited by
//! counting quality characters). Qualities are validated and discarded.
//! Gzip input is recognised by its magic bytes and rejected with a hint to
//! decompress on the fly (the crate has no dependencies).

use std::fmt;
use std::io::{self, BufRead};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadFormat {
    Fasta,
    Fastq,
}

/// One read: header line (without `>`/`@`) and raw sequence symbols.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadRecord {
    pub header: String,
    pub sequence: Vec<u8>,
}

impl ReadRecord {
    /// First whitespace-delimited token of the header.
    pub fn id(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or("")
    }
}

#[derive(Debug)]
pub enum ReadError {
    Io(io::Error),
    Gzip,
    UnknownFormat { line: usize, found: char },
    SequenceBeforeHeader { line: usize },
    ExpectedFastqHeader { line: usize },
    MissingSeparator { header: String },
    QualityLength { header: String, sequence: usize, quality: usize },
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::Gzip => write!(
                f,
                "input is gzip-compressed; decompress on the fly, e.g. `zcat reads.fastq.gz | vdj-assign ... -`"
            ),
            Self::UnknownFormat { line, found } => {
                write!(f, "line {line}: expected FASTA ('>') or FASTQ ('@'), found {found:?}")
            }
            Self::SequenceBeforeHeader { line } => write!(f, "sequence data before the first header (line {line})"),
            Self::ExpectedFastqHeader { line } => write!(f, "line {line}: expected a FASTQ header starting with '@'"),
            Self::MissingSeparator { header } => write!(f, "FASTQ record '{header}' has no '+' separator line"),
            Self::QualityLength { header, sequence, quality } => write!(
                f,
                "FASTQ record '{header}': {quality} quality characters for {sequence} bases"
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
}

impl<R: BufRead> ReadReader<R> {
    pub fn new(input: R) -> Self {
        Self { input, format: None, started: false, line: 0, buffer: String::new(), pending_header: None }
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
    pub fn read_batch(&mut self, max: usize, out: &mut Vec<ReadRecord>) -> Result<usize, ReadError> {
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
                other => Err(ReadError::UnknownFormat { line: self.line + 1, found: char::from(other) }),
            };
        }
    }

    /// Read one line into `self.buffer` without its line terminator.
    fn next_line(&mut self) -> Result<bool, ReadError> {
        self.buffer.clear();
        if self.input.read_line(&mut self.buffer)? == 0 {
            return Ok(false);
        }
        self.line += 1;
        let trimmed = self.buffer.trim_end_matches(['\n', '\r']).len();
        self.buffer.truncate(trimmed);
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
        Ok(Some(ReadRecord { header, sequence }))
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
            sequence.extend(self.buffer.trim().bytes());
        }
        let mut quality = 0;
        while quality < sequence.len() {
            if !self.next_line()? {
                break;
            }
            quality += self.buffer.trim().len();
        }
        if quality != sequence.len() {
            return Err(ReadError::QualityLength { header, sequence: sequence.len(), quality });
        }
        Ok(Some(ReadRecord { header, sequence }))
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

    fn record(header: &str, sequence: &str) -> ReadRecord {
        ReadRecord { header: header.into(), sequence: sequence.as_bytes().to_vec() }
    }

    #[test]
    fn reads_multiline_fasta() {
        let records = read_all(b"\n>r1 desc\nACGT\nac-gt\n>r2\r\nTTTT\r\n").unwrap();
        assert_eq!(records, vec![record("r1 desc", "ACGTacgt"), record("r2", "TTTT")]);
        assert_eq!(records[0].id(), "r1");
    }

    #[test]
    fn reads_fastq_including_multiline_and_awkward_qualities() {
        let text = b"@r1\nACGT\n+\nIIII\n\n@r2 x\nAC\nGT\n+r2 x\n@+II\n\n@r3\n\n+\n\n";
        let records = read_all(text).unwrap();
        assert_eq!(records, vec![record("r1", "ACGT"), record("r2 x", "ACGT"), record("r3", "")]);
    }

    #[test]
    fn detects_format_and_reports_problems() {
        let mut reader = ReadReader::new(&b"@r\nAC\n+\nII\n"[..]);
        assert_eq!(reader.format().unwrap(), Some(ReadFormat::Fastq));
        assert!(matches!(read_all(b"@r\nACGT\n+\nII\n"), Err(ReadError::QualityLength { sequence: 4, quality: 2, .. })));
        assert!(matches!(read_all(b"@r\nACGT\n"), Err(ReadError::MissingSeparator { .. })));
        assert!(matches!(read_all(b"ACGT\n"), Err(ReadError::UnknownFormat { line: 1, found: 'A' })));
        assert!(matches!(read_all(&[0x1f, 0x8b, 8, 0]), Err(ReadError::Gzip)));
        assert!(read_all(b"").unwrap().is_empty());
        assert!(read_all(b"\n  \n").unwrap().is_empty());
    }

    #[test]
    fn batches_resume_where_they_stopped() {
        let text = b">a\nA\n>b\nC\n>c\nG\n";
        let mut reader = ReadReader::new(&text[..]);
        let mut out = Vec::new();
        assert_eq!(reader.read_batch(2, &mut out).unwrap(), 2);
        assert_eq!(reader.read_batch(2, &mut out).unwrap(), 1);
        assert_eq!(reader.read_batch(2, &mut out).unwrap(), 0);
        assert_eq!(out.iter().map(|r| r.id()).collect::<Vec<_>>(), ["a", "b", "c"]);
    }
}
