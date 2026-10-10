//! The WebAssembly module of swig: a C ABI over linear memory for
//! JavaScript hosts (browsers, Web Workers, Node.js), with no bindings
//! generator. `web/swig.js` wraps it in a JavaScript API.
//!
//! The module runs the library's [`Run`] exactly as the command line does:
//! the same options by the same names, the same database preparation and
//! the same AIRR rows. A posterior can differ in its last printed digit
//! where it lands within a last-bit difference of `exp` (this target's
//! portable libm against the platform's) of a rounding tie: 2 of 260,000
//! rows in a test, with every call identical.
//!
//! # Conventions
//!
//! - Text crosses the boundary as UTF-8 in buffers from [`vdj_alloc`]. A
//!   function borrows the buffers it is given, and the host frees them with
//!   [`vdj_free`], except [`vdj_session_new`], which takes its reads.
//! - A function's text result (rows, JSON, an error message) is left in the
//!   output buffer, readable through [`vdj_output_ptr`] and
//!   [`vdj_output_len`] until the next call.
//! - Status codes are 0 for success and −1 for failure, with the message in
//!   the output buffer. Constructors return a handle, null on failure.
//! - A panic traps the instance (this target aborts on panic). Its message
//!   stays readable through [`vdj_panic_message`], and the host must then
//!   discard the instance.

use std::cell::RefCell;
use std::fmt::Write as _;
use std::io::Cursor;
use std::rc::Rc;

use swig::align::Backend;
use swig::config::CSegmentConfig;
use swig::prelude::*;
use swig::sim::{SimParams, Simulator, SyntheticDatabase, SyntheticSpec, Truth};

thread_local! {
    static OUTPUT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    static PANIC: RefCell<String> = const { RefCell::new(String::new()) };
}

/// Replace the output buffer's contents.
fn set_output(bytes: &[u8]) {
    OUTPUT.with_borrow_mut(|out| {
        out.clear();
        out.extend_from_slice(bytes);
    });
}

/// Leave `result` in the output buffer: 0 for success, −1 for failure.
fn status(result: Result<String, String>) -> i32 {
    match result {
        Ok(text) => {
            set_output(text.as_bytes());
            0
        }
        Err(message) => {
            set_output(message.as_bytes());
            -1
        }
    }
}

/// `value` as a handle, leaving `text` in the output buffer; on failure
/// null, with the message there.
fn handle<T>(result: Result<(T, String), String>) -> *mut T {
    match result {
        Ok((value, text)) => {
            set_output(text.as_bytes());
            Box::into_raw(Box::new(value))
        }
        Err(message) => {
            set_output(message.as_bytes());
            std::ptr::null_mut()
        }
    }
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
            c if c < ' ' => write!(out, "\\u{:04x}", u32::from(c)).expect("writing to a string"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// A JSON object of `fields`, whose values are JSON already.
fn json_object<'a>(fields: impl IntoIterator<Item = (&'a str, String)>) -> String {
    let fields: Vec<String> = fields
        .into_iter()
        .map(|(key, value)| format!("{}:{value}", json_string(key)))
        .collect();
    format!("{{{}}}", fields.join(","))
}

/// A JSON object mapping file names to their text.
fn json_files<'a>(files: impl IntoIterator<Item = (&'a str, String)>) -> String {
    json_object(
        files
            .into_iter()
            .map(|(name, text)| (name, json_string(&text))),
    )
}

/// Install the panic hook that keeps panic messages for
/// [`vdj_panic_message`]. Call once after instantiation.
#[unsafe(no_mangle)]
pub extern "C" fn vdj_init() {
    std::panic::set_hook(Box::new(|info| {
        let message = info.to_string();
        PANIC.with(|panic| {
            if let Ok(mut panic) = panic.try_borrow_mut() {
                *panic = message;
            }
        });
    }));
}

/// Leave the last panic's message in the output buffer (empty if none).
#[unsafe(no_mangle)]
pub extern "C" fn vdj_panic_message() {
    let message = PANIC.with(|panic| panic.try_borrow().map(|p| p.clone()).unwrap_or_default());
    OUTPUT.with(|out| {
        if let Ok(mut out) = out.try_borrow_mut() {
            out.clear();
            out.extend_from_slice(message.as_bytes());
        }
    });
}

/// `len` bytes for the host to fill; free with [`vdj_free`], or hand to
/// [`vdj_session_new`].
#[unsafe(no_mangle)]
pub extern "C" fn vdj_alloc(len: usize) -> *mut u8 {
    let mut buffer = Vec::<u8>::with_capacity(len);
    let ptr = buffer.as_mut_ptr();
    std::mem::forget(buffer);
    ptr
}

/// Free a buffer from [`vdj_alloc`].
///
/// # Safety
///
/// `ptr` and `len` must come from one [`vdj_alloc`] call, not freed before.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_free(ptr: *mut u8, len: usize) {
    // SAFETY: the buffer was allocated as a `Vec<u8>` of capacity `len`.
    drop(unsafe { Vec::from_raw_parts(ptr, 0, len) });
}

/// Start of the output buffer.
#[unsafe(no_mangle)]
pub extern "C" fn vdj_output_ptr() -> *const u8 {
    OUTPUT.with_borrow(|out| out.as_ptr())
}

/// Length of the output buffer in bytes.
#[unsafe(no_mangle)]
pub extern "C" fn vdj_output_len() -> usize {
    OUTPUT.with_borrow(Vec::len)
}

/// Leave `{"version", "backend"}` in the output buffer: the crate version
/// and the alignment kernel a run uses by default.
#[unsafe(no_mangle)]
pub extern "C" fn vdj_version() {
    let info = json_object([
        ("version", json_string(env!("CARGO_PKG_VERSION"))),
        ("backend", json_string(Backend::detect().name())),
    ]);
    set_output(info.as_bytes());
}

/// Germline texts and options gathered for a run or a preparation.
#[derive(Default)]
pub struct Builder {
    v: Option<String>,
    d: Option<String>,
    j: Option<String>,
    c: Option<String>,
    options: RunOptions,
}

impl Builder {
    /// Set germline set `v`, `d`, `j` or `c` to FASTA text, or any other key
    /// as a [`RunOptions`] name.
    fn set(&mut self, key: &str, value: Option<&str>) -> Result<(), String> {
        let slot = match key {
            "v" => &mut self.v,
            "d" => &mut self.d,
            "j" => &mut self.j,
            "c" => &mut self.c,
            _ => return self.options.set(key, value).map_err(|e| e.to_string()),
        };
        *slot = Some(
            value
                .ok_or_else(|| format!("{key} needs FASTA text"))?
                .to_string(),
        );
        Ok(())
    }

    fn database(&self) -> Result<PreparedDatabase, String> {
        let texts = GermlineTexts {
            v: self.v.as_deref().ok_or("a V database is required")?,
            d: self.d.as_deref(),
            j: self.j.as_deref().ok_or("a J database is required")?,
            c: self.c.as_deref(),
        };
        PreparedDatabase::from_fasta(texts, self.options.annotate())
            .map_err(|e| format!("{} database: {e}", e.kind().label()))
    }
}

/// A new, empty builder.
#[unsafe(no_mangle)]
pub extern "C" fn vdj_builder_new() -> *mut Builder {
    Box::into_raw(Box::default())
}

/// Free a builder that [`vdj_run_new`] did not take.
///
/// # Safety
///
/// `builder` must come from [`vdj_builder_new`] and not be freed or taken.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_builder_free(builder: *mut Builder) {
    // SAFETY: forwarded from the caller.
    drop(unsafe { Box::from_raw(builder) });
}

/// Set `key` (`v`, `d`, `j`, `c` for FASTA text, or an option name as on
/// the command line) to the value, or, with `has_value` 0, set flag `key`.
///
/// # Safety
///
/// `builder` must be live; each pointer must point to its length in bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_builder_set(
    builder: *mut Builder,
    key_ptr: *const u8,
    key_len: usize,
    value_ptr: *const u8,
    value_len: usize,
    has_value: i32,
) -> i32 {
    // SAFETY: forwarded from the caller.
    let (builder, key, value) = unsafe {
        (
            &mut *builder,
            text(key_ptr, key_len),
            text(value_ptr, value_len),
        )
    };
    status(key.and_then(|key| {
        let value = (has_value != 0).then_some(value).transpose()?;
        builder.set(key, value).map(|()| String::new())
    }))
}

/// Prepare the builder's database and leave its files as a JSON object
/// (file name → text) in the output buffer, as `swig prepare`
/// writes them.
///
/// # Safety
///
/// `builder` must be live.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_builder_prepare(builder: *const Builder) -> i32 {
    // SAFETY: forwarded from the caller.
    let builder = unsafe { &*builder };
    status(
        builder
            .database()
            .map(|database| json_files(database.files())),
    )
}

/// A run built from the builder, which it takes. The output buffer holds
/// `{"summary", "warnings", "backend", "header"}`: the database summary,
/// its warnings, the alignment kernel and the AIRR TSV header.
///
/// # Safety
///
/// `builder` must come from [`vdj_builder_new`] and not be freed or taken.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_run_new(builder: *mut Builder) -> *mut Rc<Run> {
    // SAFETY: forwarded from the caller.
    let builder = unsafe { Box::from_raw(builder) };
    handle((|| {
        let database = builder.database()?;
        let warnings: Vec<String> = database
            .warnings(CSegmentConfig::default().spread)
            .iter()
            .map(|w| json_string(&w.to_string()))
            .collect();
        let summary = database.summary();
        let run = Run::new(&builder.options, database).map_err(|e| e.to_string())?;
        let info = json_object([
            ("summary", json_string(&summary)),
            ("warnings", format!("[{}]", warnings.join(","))),
            ("backend", json_string(run.backend().name())),
            ("header", json_string(&run.airr_header())),
        ]);
        Ok((Rc::new(run), info))
    })())
}

/// Free a run. Sessions keep their own reference, so this may precede
/// freeing them.
///
/// # Safety
///
/// `run` must come from [`vdj_run_new`] and not be freed before.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_run_free(run: *mut Rc<Run>) {
    // SAFETY: forwarded from the caller.
    drop(unsafe { Box::from_raw(run) });
}

/// One worker's share of a reads file: every `workers`-th block of
/// `block` consecutive reads, starting at block `worker`.
pub struct Session {
    run: Rc<Run>,
    reader: ReadReader<Cursor<Vec<u8>>>,
    records: Vec<swig::reads::ReadRecord>,
    len: usize,
    next_block: usize,
    worker: usize,
    workers: usize,
    block: usize,
}

impl Session {
    /// Assign the next block of this worker's share, writing its rows to
    /// `out`; `None` at the end of the input.
    fn next(&mut self, out: &mut Vec<u8>) -> Result<Option<usize>, String> {
        loop {
            self.records.clear();
            let read = self
                .reader
                .read_batch(self.block, &mut self.records)
                .map_err(|e| e.to_string())?;
            if read == 0 {
                return Ok(None);
            }
            let block = self.next_block;
            self.next_block += 1;
            if block % self.workers == self.worker {
                self.run
                    .write_rows(&self.records, out)
                    .map_err(|e| e.to_string())?;
                return Ok(Some(block));
            }
        }
    }
}

/// A session over the FASTA or FASTQ reads in the `len` bytes at `reads`,
/// which it takes (the host must not free them). It assigns block `worker`,
/// `worker + workers`, … of `block` reads each, so `workers` sessions over
/// the same reads, one per Web Worker, share the work; their blocks in
/// order are the whole output.
///
/// # Safety
///
/// `run` must be live; `reads` and `len` must come from one [`vdj_alloc`]
/// call, filled by the host.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_session_new(
    run: *const Rc<Run>,
    reads: *mut u8,
    len: usize,
    worker: u32,
    workers: u32,
    block: u32,
) -> *mut Session {
    // SAFETY: forwarded from the caller; the host filled all `len` bytes.
    let (run, bytes) = unsafe { (&*run, Vec::from_raw_parts(reads, len, len)) };
    handle((|| {
        let (worker, workers, block) = (worker as usize, workers as usize, block as usize);
        if workers == 0 || worker >= workers || block == 0 {
            return Err(format!(
                "worker {worker} of {workers} with blocks of {block} reads is not a share of the work"
            ));
        }
        let session = Session {
            run: Rc::clone(run),
            reader: ReadReader::new(Cursor::new(bytes)),
            // Grown by the reader: `block` comes from the host unbounded.
            records: Vec::new(),
            len,
            next_block: 0,
            worker,
            workers,
            block,
        };
        Ok((session, String::new()))
    })())
}

/// Assign the session's next block: its index, with its AIRR rows (one
/// line each, no header) in the output buffer; −1 at the end of the reads,
/// −2 on a read error, with the message in the output buffer.
///
/// # Safety
///
/// `session` must be live.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_session_next(session: *mut Session) -> f64 {
    // SAFETY: forwarded from the caller.
    let session = unsafe { &mut *session };
    let mut rows = Vec::new();
    match session.next(&mut rows) {
        Ok(Some(block)) => {
            set_output(&rows);
            block as f64
        }
        Ok(None) => {
            set_output(&[]);
            -1.0
        }
        Err(message) => {
            set_output(message.as_bytes());
            -2.0
        }
    }
}

/// The fraction of the session's reads consumed so far.
///
/// # Safety
///
/// `session` must be live.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_session_progress(session: *const Session) -> f64 {
    // SAFETY: forwarded from the caller.
    let session = unsafe { &*session };
    match session.len {
        0 => 1.0,
        len => session.reader.get_ref().position() as f64 / len as f64,
    }
}

/// Free a session and its reads.
///
/// # Safety
///
/// `session` must come from [`vdj_session_new`] and not be freed before.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_session_free(session: *mut Session) {
    // SAFETY: forwarded from the caller.
    drop(unsafe { Box::from_raw(session) });
}

/// Demonstration data: the synthetic IGH-like database of `seed` and
/// `reads` reads simulated under IgSim preset `preset` (half of them
/// reverse-complemented), as a JSON object of files: `V.fasta`, `D.fasta`,
/// `J.fasta`, `C.fasta`, `reads.fasta` and `truth.tsv`.
///
/// # Safety
///
/// `preset` must point to `preset_len` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn vdj_demo(
    seed: u32,
    reads: u32,
    preset: *const u8,
    preset_len: usize,
) -> i32 {
    // SAFETY: forwarded from the caller.
    let preset = unsafe { text(preset, preset_len) };
    status(preset.and_then(|preset| demo(u64::from(seed), reads as usize, preset)))
}

fn demo(seed: u64, count: usize, preset: &str) -> Result<String, String> {
    let params = SimParams::preset(preset)
        .ok_or_else(|| format!("unknown preset '{preset}' (easy, mid, train, hard)"))?
        .with_reverse_complement(0.5);
    let db = SyntheticDatabase::generate(&SyntheticSpec::default(), seed);
    let sim = Simulator::new(&db.v, Some(&db.d), &db.j, params)
        .and_then(|sim| sim.with_constant(&db.c))
        .map_err(|e| e.to_string())?;
    let (mut fasta, mut truth) = (String::new(), Truth::tsv_header_with(true));
    // Reads draw from a stream of their own, not the database's.
    let reads = sim
        .reads(count, seed.wrapping_add(1))
        .map_err(|e| e.to_string())?;
    for (index, read) in reads.iter().enumerate() {
        let id = format!("sim{index}");
        let sequence = String::from_utf8_lossy(&read.sequence);
        writeln!(fasta, ">{id}\n{sequence}").expect("writing to a string");
        write!(truth, "\n{}", Truth::from_read(&sim, read, &id).to_tsv())
            .expect("writing to a string");
    }
    truth.push('\n');
    let mut files = Vec::new();
    for (name, set) in [
        ("V.fasta", &db.v),
        ("D.fasta", &db.d),
        ("J.fasta", &db.j),
        ("C.fasta", &db.c),
    ] {
        let mut text = Vec::new();
        set.write_fasta(&mut text, |_| None)
            .expect("writing to memory");
        files.push((
            name,
            String::from_utf8(text).expect("germline FASTA is ASCII"),
        ));
    }
    files.push(("reads.fasta", fasta));
    files.push(("truth.tsv", truth));
    Ok(json_files(files))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The output buffer as text.
    fn output() -> String {
        OUTPUT.with_borrow(|out| String::from_utf8(out.clone()).unwrap())
    }

    /// A buffer from `vdj_alloc` holding `text`.
    fn buffer(text: &str) -> (*mut u8, usize) {
        let ptr = vdj_alloc(text.len());
        // SAFETY: `ptr` has room for `text`.
        unsafe { std::ptr::copy_nonoverlapping(text.as_ptr(), ptr, text.len()) };
        (ptr, text.len())
    }

    /// The value of string field `key` in a flat JSON object of files.
    fn field(json: &str, key: &str) -> String {
        let start = json.find(&format!("\"{key}\":\"")).unwrap() + key.len() + 4;
        let mut value = String::new();
        let mut chars = json[start..].chars();
        while let Some(c) = chars.next() {
            match c {
                '"' => break,
                '\\' => value.push(match chars.next().unwrap() {
                    'n' => '\n',
                    't' => '\t',
                    'r' => '\r',
                    other => other,
                }),
                c => value.push(c),
            }
        }
        value
    }

    fn set(builder: *mut Builder, key: &str, value: Option<&str>) -> i32 {
        let value_text = value.unwrap_or("");
        // SAFETY: the pointers cover the strings.
        unsafe {
            vdj_builder_set(
                builder,
                key.as_ptr(),
                key.len(),
                value_text.as_ptr(),
                value_text.len(),
                i32::from(value.is_some()),
            )
        }
    }

    #[test]
    fn json_strings_escape_what_json_requires() {
        assert_eq!(
            json_string("a\"b\\c\nd\te\u{1}"),
            "\"a\\\"b\\\\c\\nd\\te\\u0001\""
        );
        assert_eq!(
            json_object([("k", json_string("v")), ("n", "[1]".into())]),
            "{\"k\":\"v\",\"n\":[1]}"
        );
    }

    /// Sessions sharing the work give the command line's rows, in order.
    #[test]
    fn sessions_of_several_workers_give_the_rows_of_one_run() {
        let preset = "mid";
        // SAFETY: the pointer covers the preset.
        assert_eq!(
            unsafe { vdj_demo(3, 300, preset.as_ptr(), preset.len()) },
            0
        );
        let files = output();
        let reads = field(&files, "reads.fasta");
        assert_eq!(reads.lines().count(), 600);

        let builder = vdj_builder_new();
        for set_name in ["v", "d", "j", "c"] {
            let fasta = field(&files, &format!("{}.fasta", set_name.to_uppercase()));
            assert_eq!(set(builder, set_name, Some(&fasta)), 0);
        }
        assert_eq!(set(builder, "tandem-d", None), 0);
        assert_eq!(set(builder, "profile", Some("best")), -1);
        assert!(output().contains("--profile expects"));
        // SAFETY: the builder is live.
        assert_eq!(unsafe { vdj_builder_prepare(builder) }, 0);
        assert!(field(&output(), "annotation.tsv").starts_with("segment\tallele\t"));
        // SAFETY: the builder is live, and taken here.
        let run = unsafe { vdj_run_new(builder) };
        assert!(!run.is_null(), "{}", output());
        let header = field(&output(), "header");

        let mut options = RunOptions::default();
        options.set("tandem-d", None).unwrap();
        let texts = |name: &str| field(&files, name);
        let (v, d, j, c) = (
            texts("V.fasta"),
            texts("D.fasta"),
            texts("J.fasta"),
            texts("C.fasta"),
        );
        let database = PreparedDatabase::from_fasta(
            GermlineTexts {
                v: &v,
                d: Some(&d),
                j: &j,
                c: Some(&c),
            },
            options.annotate(),
        )
        .unwrap();
        let mut expected = Vec::new();
        Run::new(&options, database)
            .unwrap()
            .write_airr(reads.as_bytes(), &mut expected)
            .unwrap();
        let expected = String::from_utf8(expected).unwrap();

        let workers = 3;
        let mut blocks: Vec<(usize, String)> = Vec::new();
        for worker in 0..workers {
            let (ptr, len) = buffer(&reads);
            // SAFETY: the run is live; the buffer is filled and taken.
            let session = unsafe { vdj_session_new(run, ptr, len, worker, workers, 64) };
            assert!(!session.is_null(), "{}", output());
            loop {
                // SAFETY: the session is live.
                let block = unsafe { vdj_session_next(session) };
                if block < 0.0 {
                    assert_eq!(block, -1.0, "{}", output());
                    break;
                }
                blocks.push((block as usize, output()));
            }
            // SAFETY: the session is live.
            assert_eq!(unsafe { vdj_session_progress(session) }, 1.0);
            // SAFETY: the session is live and freed once.
            unsafe { vdj_session_free(session) };
        }
        blocks.sort();
        let rows: String = blocks.into_iter().map(|(_, rows)| rows).collect();
        assert_eq!(format!("{header}\n{rows}"), expected);
        // SAFETY: the run is live and freed once.
        unsafe { vdj_run_free(run) };
    }

    #[test]
    fn errors_name_what_failed() {
        let builder = vdj_builder_new();
        assert_eq!(set(builder, "v", Some("@read\nACGT\n+\nIIII\n")), 0);
        assert_eq!(set(builder, "j", Some(">J1*01\nACGTACGT\n")), 0);
        // SAFETY: the builder is live, and taken here.
        let run = unsafe { vdj_run_new(builder) };
        assert!(run.is_null());
        assert!(output().starts_with("V database: "), "{}", output());

        let builder = vdj_builder_new();
        assert_eq!(set(builder, "v", None), -1);
        // SAFETY: the builder is live and freed once.
        unsafe { vdj_builder_free(builder) };
        let (ptr, len) = buffer("");
        // SAFETY: as allocated.
        unsafe { vdj_free(ptr, len) };
    }
}
