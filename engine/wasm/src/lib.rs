//! The WebAssembly module of swig: a C ABI over linear memory for Swig's
//! Web Workers and its Node.js/Bun command line, with no bindings generator
//! and no imports.
//!
//! Each instance holds one engine. A worker sets the run options by name
//! ([`swig_set_option`], the command line's vocabulary), loads the germline
//! databases ([`swig_init_database`]) and then assigns batches of reads
//! ([`swig_annotate`]): FASTA, FASTQ or AIRR-table text in, AIRR TSV (and
//! the double-D sidecar) out. The module runs the library's [`Run`] exactly
//! as the native `swig` command line does: the same options, the same
//! database preparation and the same rows.
//!
//! # Conventions
//!
//! - Text crosses the boundary as UTF-8 in buffers from [`swig_alloc`]. A
//!   function borrows the buffers it is given; the host frees them with
//!   [`swig_free`].
//! - Results stay readable until the next call: the AIRR TSV (or other
//!   text) through [`swig_result_ptr`]/[`swig_result_len`], the double-D
//!   sidecar through [`swig_double_d_result_ptr`]/[`swig_double_d_result_len`],
//!   an error message through [`swig_error_ptr`]/[`swig_error_len`].
//! - A function returns −1 on failure, with the message in the error
//!   buffer.
//! - A panic traps the instance (this target aborts on panic). Its message
//!   stays readable after [`swig_panic_message`], and the host must then
//!   discard the instance.

use std::cell::RefCell;
use std::sync::Once;

use swig::prelude::*;

/// An instance's engine and buffers.
#[derive(Default)]
struct State {
    options: RunOptions,
    /// The germline FASTA texts, kept to rebuild the run when an option
    /// changes.
    germlines: Option<[Option<String>; 4]>,
    run: Option<Run>,
    result: Vec<u8>,
    double_d: Vec<u8>,
    double_d_count: i32,
    error: String,
    panic: String,
}

thread_local! {
    static STATE: RefCell<State> = RefCell::new(State::default());
}

static PANIC_HOOK: Once = Once::new();

/// Keep a panic's message for [`swig_panic_message`].
fn install_panic_hook() {
    PANIC_HOOK.call_once(|| {
        std::panic::set_hook(Box::new(|info| {
            let message = info.to_string();
            STATE.with(|state| {
                if let Ok(mut state) = state.try_borrow_mut() {
                    state.panic = message;
                }
            });
        }));
    });
}

/// Run `f` on the instance state; on `Err`, keep the message and return −1.
fn with_state(f: impl FnOnce(&mut State) -> Result<i32, String>) -> i32 {
    install_panic_hook();
    STATE.with_borrow_mut(|state| {
        state.error.clear();
        match f(state) {
            Ok(value) => value,
            Err(message) => {
                state.error = message;
                -1
            }
        }
    })
}

/// The UTF-8 text in `len` bytes at `ptr`.
///
/// # Safety
///
/// `ptr` must point to `len` readable bytes (any pointer when `len` is 0).
unsafe fn text<'a>(ptr: *const u8, len: usize) -> Result<&'a str, String> {
    let bytes = match len {
        0 => &[],
        // SAFETY: forwarded from the caller.
        _ => unsafe { std::slice::from_raw_parts(ptr, len) },
    };
    std::str::from_utf8(bytes).map_err(|e| format!("text is not UTF-8: {e}"))
}

/// `text` as a JSON string literal.
fn json_string(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c < ' ' => out.push_str(&format!("\\u{:04x}", u32::from(c))),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

impl State {
    /// The run, built from the options and germlines when it is not yet.
    fn run(&mut self) -> Result<&mut Run, String> {
        if self.run.is_none() {
            let [v, d, j, c] = self
                .germlines
                .as_ref()
                .ok_or("load a reference database before annotation")?;
            let texts = GermlineTexts {
                v: v.as_deref().ok_or("a V germline FASTA is required")?,
                d: d.as_deref(),
                j: j.as_deref().ok_or("a J germline FASTA is required")?,
                c: c.as_deref(),
            };
            let database = PreparedDatabase::from_fasta(texts, self.options.annotate())
                .map_err(|e| format!("{} database: {e}", e.kind().label()))?;
            let summary = format!(
                "{{\"summary\":{},\"warnings\":[{}]}}",
                json_string(&database.summary()),
                database
                    .warnings(CSegmentConfig::default().spread)
                    .iter()
                    .map(|w| json_string(&w.to_string()))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            self.run = Some(Run::new(&self.options, database).map_err(|e| e.to_string())?);
            self.result = summary.into_bytes();
        }
        Ok(self.run.as_mut().expect("just built"))
    }
}

/// `len` bytes for the host to fill; free with [`swig_free`].
#[unsafe(no_mangle)]
pub extern "C" fn swig_alloc(len: usize) -> *mut u8 {
    let mut buffer = Vec::<u8>::with_capacity(len.max(1));
    let ptr = buffer.as_mut_ptr();
    std::mem::forget(buffer);
    ptr
}

/// Free a buffer from [`swig_alloc`].
///
/// # Safety
///
/// `ptr` and `len` must come from one [`swig_alloc`] call, not freed before.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swig_free(ptr: *mut u8, len: usize) {
    // SAFETY: the buffer was allocated as a `Vec<u8>` of capacity `len.max(1)`.
    drop(unsafe { Vec::from_raw_parts(ptr, 0, len.max(1)) });
}

/// Set run option `name` (as on the command line, without `--`) to the
/// value, or set flag `name` (`value_len` 0). The run is rebuilt from the
/// loaded germlines on the next annotation. 0, or −1 for an unknown option
/// or a bad value.
///
/// # Safety
///
/// Each pointer must point to its length in bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swig_set_option(
    name: *const u8,
    name_len: usize,
    value: *const u8,
    value_len: usize,
) -> i32 {
    with_state(|state| {
        // SAFETY: forwarded from the caller.
        let (name, value) = unsafe { (text(name, name_len)?, text(value, value_len)?) };
        let value = (RunOptions::takes_value(name) == Some(true)).then_some(value);
        state.options.set(name, value).map_err(|e| e.to_string())?;
        state.run = None;
        Ok(0)
    })
}

/// Restore every option's default; the run is rebuilt on the next
/// annotation.
#[unsafe(no_mangle)]
pub extern "C" fn swig_reset_options() {
    with_state(|state| {
        state.options = RunOptions::default();
        state.run = None;
        Ok(0)
    });
}

/// Load the germline databases (V and J required, D and C empty when
/// absent) and prepare them under the current options. Returns the number
/// of records read, or −1; the result buffer then holds
/// `{"summary", "warnings"}` as JSON.
///
/// # Safety
///
/// Each pointer must point to its length in bytes.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn swig_init_database(
    v: *const u8,
    v_len: usize,
    d: *const u8,
    d_len: usize,
    j: *const u8,
    j_len: usize,
    c: *const u8,
    c_len: usize,
) -> i32 {
    with_state(|state| {
        let mut texts: [Option<String>; 4] = Default::default();
        let mut records = 0;
        for (slot, (ptr, len)) in
            texts
                .iter_mut()
                .zip([(v, v_len), (d, d_len), (j, j_len), (c, c_len)])
        {
            // SAFETY: forwarded from the caller.
            let fasta = unsafe { text(ptr, len)? };
            if !fasta.trim().is_empty() {
                records += fasta.lines().filter(|l| l.starts_with('>')).count();
                *slot = Some(fasta.to_string());
            }
        }
        state.germlines = Some(texts);
        state.run = None;
        state.run()?;
        Ok(i32::try_from(records).unwrap_or(i32::MAX))
    })
}

/// Assign the reads in `len` bytes at `reads` (FASTA, FASTQ, or an AIRR
/// table with a header line). Returns the number of reads, or −1. The
/// result buffer holds the AIRR TSV, header first and one line per read in
/// input order; with the double-D screen on, the double-D buffer holds the
/// sidecar, header first, a row per screened read.
///
/// # Safety
///
/// `reads` must point to `len` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn swig_annotate(reads: *const u8, len: usize) -> i32 {
    with_state(|state| {
        // SAFETY: forwarded from the caller.
        let input = unsafe { text(reads, len)? };
        let records: Vec<ReadRecord> = ReadReader::new(input.as_bytes())
            .collect::<Result<_, _>>()
            .map_err(|e| e.to_string())?;
        let run = state.run()?;
        let rows = run.assign(&records);
        let mut result = run.airr_header().into_bytes();
        result.push(b'\n');
        result.extend_from_slice(rows.airr.as_bytes());
        let mut double_d = run.double_d_header().into_bytes();
        double_d.push(b'\n');
        double_d.extend_from_slice(rows.double_d.as_bytes());
        state.result = result;
        state.double_d = double_d;
        state.double_d_count = i32::try_from(rows.tally.double_d).unwrap_or(i32::MAX);
        Ok(i32::try_from(rows.tally.reads).unwrap_or(i32::MAX))
    })
}

/// Start of the result buffer.
#[unsafe(no_mangle)]
pub extern "C" fn swig_result_ptr() -> *const u8 {
    STATE.with_borrow(|state| state.result.as_ptr())
}

/// Length of the result buffer in bytes.
#[unsafe(no_mangle)]
pub extern "C" fn swig_result_len() -> usize {
    STATE.with_borrow(|state| state.result.len())
}

/// Start of the double-D sidecar of the last annotation.
#[unsafe(no_mangle)]
pub extern "C" fn swig_double_d_result_ptr() -> *const u8 {
    STATE.with_borrow(|state| state.double_d.as_ptr())
}

/// Length of the double-D sidecar in bytes.
#[unsafe(no_mangle)]
pub extern "C" fn swig_double_d_result_len() -> usize {
    STATE.with_borrow(|state| state.double_d.len())
}

/// Rows of the double-D sidecar of the last annotation.
#[unsafe(no_mangle)]
pub extern "C" fn swig_double_d_count() -> i32 {
    STATE.with_borrow(|state| state.double_d_count)
}

/// Start of the last error message.
#[unsafe(no_mangle)]
pub extern "C" fn swig_error_ptr() -> *const u8 {
    STATE.with_borrow(|state| state.error.as_ptr())
}

/// Length of the last error message in bytes.
#[unsafe(no_mangle)]
pub extern "C" fn swig_error_len() -> usize {
    STATE.with_borrow(|state| state.error.len())
}

/// Leave the last panic's message in the error buffer (empty if none).
#[unsafe(no_mangle)]
pub extern "C" fn swig_panic_message() {
    STATE.with(|state| {
        if let Ok(mut state) = state.try_borrow_mut() {
            state.error = state.panic.clone();
        }
    });
}

/// Leave `{"version", "backend"}` in the result buffer: the engine version
/// and the alignment kernel it runs.
#[unsafe(no_mangle)]
pub extern "C" fn swig_version() {
    with_state(|state| {
        state.result = format!(
            "{{\"version\":{},\"backend\":{}}}",
            json_string(env!("CARGO_PKG_VERSION")),
            json_string(swig::align::Backend::detect().name())
        )
        .into_bytes();
        Ok(0)
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec};

    fn result() -> String {
        STATE.with_borrow(|state| String::from_utf8(state.result.clone()).unwrap())
    }

    fn error() -> String {
        STATE.with_borrow(|state| state.error.clone())
    }

    fn set(name: &str, value: &str) -> i32 {
        // SAFETY: both slices are live for the call.
        unsafe { swig_set_option(name.as_ptr(), name.len(), value.as_ptr(), value.len()) }
    }

    fn init(texts: [&str; 4]) -> i32 {
        let [v, d, j, c] = texts;
        // SAFETY: every slice is live for the call.
        unsafe {
            swig_init_database(
                v.as_ptr(),
                v.len(),
                d.as_ptr(),
                d.len(),
                j.as_ptr(),
                j.len(),
                c.as_ptr(),
                c.len(),
            )
        }
    }

    fn annotate(reads: &str) -> i32 {
        // SAFETY: the slice is live for the call.
        unsafe { swig_annotate(reads.as_ptr(), reads.len()) }
    }

    /// The module's rows are the library run's, for FASTA and AIRR input;
    /// options rebuild the run and bad ones are refused.
    #[test]
    fn the_abi_runs_the_library_run() {
        let db = SyntheticDatabase::generate(&SyntheticSpec::default(), 5);
        let fasta = |set: &swig::germline::GermlineSet| {
            let mut text = Vec::new();
            set.write_fasta(&mut text, |_| None).unwrap();
            String::from_utf8(text).unwrap()
        };
        let (v, d, j, c) = (fasta(&db.v), fasta(&db.d), fasta(&db.j), fasta(&db.c));
        let sim = Simulator::new(&db.v, Some(&db.d), &db.j, SimParams::mid()).unwrap();
        let reads = sim.reads(20, 7).unwrap();
        let text: String = reads
            .iter()
            .enumerate()
            .map(|(i, r)| format!(">r{i}\n{}\n", String::from_utf8_lossy(&r.sequence)))
            .collect();

        assert_eq!(annotate(&text), -1);
        assert!(error().contains("reference database"), "{}", error());
        assert_eq!(
            init([&v, &d, &j, &c]),
            (db.v.len() + db.d.len() + db.j.len() + db.c.len()) as i32
        );
        assert!(result().starts_with("{\"summary\":"), "{}", result());
        assert_eq!(annotate(&text), 20);
        let rows = result();
        let options = RunOptions {
            threads: 1,
            ..RunOptions::default()
        };
        let texts = GermlineTexts {
            v: &v,
            d: Some(&d),
            j: &j,
            c: Some(&c),
        };
        let mut run = Run::new(
            &options,
            PreparedDatabase::from_fasta(texts, options.annotate()).unwrap(),
        )
        .unwrap();
        let mut expected = Vec::new();
        run.write_airr(text.as_bytes(), &mut expected, None)
            .unwrap();
        assert_eq!(rows, String::from_utf8(expected).unwrap());

        // The same reads as an AIRR table give the same rows.
        let table: String = std::iter::once("sequence_id\tsequence\n".to_string())
            .chain(
                reads
                    .iter()
                    .enumerate()
                    .map(|(i, r)| format!("r{i}\t{}\n", String::from_utf8_lossy(&r.sequence))),
            )
            .collect();
        assert_eq!(annotate(&table), 20);
        assert_eq!(result(), rows);

        assert_eq!(set("strand", "sideways"), -1);
        assert!(error().contains("--strand"), "{}", error());
        assert_eq!(set("no-posterior", ""), 0);
        assert_eq!(annotate(&text), 20);
        let header = rows.lines().next().unwrap();
        let posterior = header.split('\t').position(|c| c == "v_posterior").unwrap();
        let row = result().lines().nth(1).unwrap().to_string();
        assert_eq!(row.split('\t').nth(posterior), Some(""));
        swig_reset_options();
        assert_eq!(annotate(&text), 20);
        assert_eq!(result(), rows);

        swig_version();
        assert!(result().contains(env!("CARGO_PKG_VERSION")));
    }
}
