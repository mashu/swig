//! The `swig` command line: its subcommands' outputs, refusals and exit
//! codes (0 success, 1 a failure while running, 2 a bad command line).

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A small single-locus (IGH) database with header anchors.
const V: &str = "\
>IGHV1-1*01 anchor=288
TTTCCTCACGCAATTCAAAGCCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGAATACCAAGTTCCGCCTTATTCAGGTCCTAACCGGTAAACCCGGTCTCACCGCCCCCTTAAGGCTGTTCCATCTAGCCAAGTTCATCGGCAGCTGCAATGGAAATAGGCAACTGATATATATTAAAAAGGGTGTTAAGATACATGGCCCGTTGGTGATCCTTCCCCTGAAGCATCGCTTTGTGAAGCGGGACTTGCGCCATATCTGCATACGGGCTCATTCTTCATGTTCGAGTGA
>IGHV1-1*02 anchor=288
TTTCCTCACGCAATTCAAAGCCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGAATACCAAGGTCCGCCTTATTCAGGTCCTAACCGGTAAACCCGCTCTCACCGCCCCCTTAAGGCTGTTCCATCTAGCCAAGTTCATCGGCAGCTGCAATGGAAATAGGCAACTGATATATATTAAAAAGGGTGTTAAGATACATGGCCCGTTGGTGATCCTTCCCCTGAAGCATCGCTTTGTGAAGCGGGACTTGCGCCATATCTGCATACGGGCTCATTCTTCATGTTCGAGTGA
>IGHV1-2*01 anchor=288
TTTCCTCATGCCATTGAAAACCATGTCCGTAATGTAGGCGAAATAGTTAACCATTATACGGACGACTCCAAATTCCTCCTTATTCCGGACCTAACCGGTCAAGCAGGTCTCTCCGCACCCTTAAAGCTGTTGCACCTCGCCAAGTCCAACGGCAGCTGCAATGGAAATAGGCAACGGATATCTATTAAAAAGTATTTTAAGATACATGGCCGGTGCGTGCTCGTCGCCCTGAAGCATTGCTTCGTGAAGAGGGACTTAAGCCAAGCCTGCATACAGGCGCATTCTTCATGTGCCAGAGA
>IGHV1-2*02 anchor=288
TTTCCTCATGCCATTGAAAACCATGTCCGTAATGTAGGCGAAATAGTTAACCATTATACGGACGACTCCAAATTCCTCCTTATTCCGGACCTAACCGGTCAAGCAGGTCTCTCCGCACCCTTAAAGCTTTTGCACCTCGCCAAGTCCAACGGCAGCTGCAATGGAAATAGGCAACGGATATCTATTAAAAAGTATTTTAAGATACATGGCCGGTGCGTGCTCGTCGCCCTGAAGCATTGCTTCGTGAAGAGGGACTTAAGCCAAGCCTGCATACAGGCGCATTCTTCATGTGCCAGAGC
";
const D: &str = "\
>IGHD1-1*01
TATTCCCCCCGCGGCCCACCC
>IGHD2-1*01
GTATTCCTAACGGAGCA
>IGHD3-1*01
AAATCCCACCCGAACTAAGTTTGTCGAACC
";
const J: &str = "\
>IGHJ1*01 anchor=15
GACTCGGTCTCCAGGTGGGGCCAAGGGACCACGGTCATCGTCTCCACAG
>IGHJ2*01 anchor=18
GTTACTAAGGGGTATAATTGGGGCCAAGCGACCACGGTCACCGTCTCCACAG
";
const C: &str = "\
>IGHG1*01
TATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAGTCTTGTACGCCTAGGGGCGCCCCAAAGGTAAACGAACCGTTGCGGTCAATCTTGTCGCGGCTGATGAATTTGAAGCAGTGGCCGGGTGTGTGTGCTCAGGAGTTCGTCCCACACGATAGCGAGAGAACATCCTGTTGGGCTTATAGAATCCCATCGGTTGGAGCCATAACCGCCTCTCGTCGTGTTGATCTATGA
>IGHG2*01
TATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAGTCTTGTGCGCCTAGGGGCGCCCCAAAGGTAAACGAACCGTTGCGGTCAATCTTGTCGCGGCTGATGCATTTGAAGCAGTGGCCGGGAGTGTGTGCTCAGGAGTTCGACCCACACGATAGAGAGAGAACATCCTCTTGGGCTTATAGAATTCCCTCGCTTGGAGCCATAACCGCCCCTCGACGCGTTGATCTATGA
";

/// A fresh directory holding the database as V.fasta, D.fasta, J.fasta and
/// C.fasta.
fn workspace(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("swig-cli-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    for (file, text) in [("V", V), ("D", D), ("J", J), ("C", C)] {
        fs::write(dir.join(format!("{file}.fasta")), text).unwrap();
    }
    dir
}

fn swig(dir: &PathBuf, args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_swig"));
    command.current_dir(dir).args(args);
    command
}

fn run(dir: &PathBuf, args: &[&str]) -> Output {
    swig(dir, args).output().unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// FNV-1a, 64 bits.
fn fnv(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3)
    })
}

const GERMLINES: [&str; 8] = [
    "--v", "V.fasta", "--d", "D.fasta", "--j", "J.fasta", "--c", "C.fasta",
];

/// The reads and truth of a single-locus database are the ones Swig 0.39.1
/// wrote: the simulator keeps IgSim's order of draws.
#[test]
fn simulating_a_single_locus_keeps_its_reads() {
    let dir = workspace("pinned");
    let args = [
        "--out", "sim", "--reads", "300", "--seed", "7", "--preset", "hard", "--dd", "0.3",
    ];
    let output = run(&dir, &[&["simulate"], &GERMLINES[..], &args].concat());
    assert!(output.status.success(), "{}", stderr(&output));
    let hash = |file: &str| fnv(&fs::read(dir.join(file)).unwrap());
    assert_eq!(hash("sim.fasta"), 0x2535_80e9_671d_5e21);
    assert_eq!(hash("sim.truth.tsv"), 0xbc14_6828_0a92_8d33);
    fs::remove_dir_all(dir).unwrap();
}

/// Probabilities and rates outside [0, 1] (a percentage typed for a
/// fraction, NaN) and an empty simulation are refused, naming the option.
#[test]
fn simulate_refuses_out_of_range_parameters() {
    let dir = workspace("ranges");
    let simulate = |extra: &[&str]| {
        let reads: &[&str] = match extra.contains(&"--reads") {
            true => &[],
            false => &["--reads", "20"],
        };
        let base = ["simulate", "--out", "sim"];
        run(&dir, &[&base[..], &GERMLINES[..], reads, extra].concat())
    };
    for (option, value) in [
        ("--rc", "50"),
        ("--rc", "-1"),
        ("--shm", "3"),
        ("--shm", "NaN"),
        ("--shm", "inf"),
        ("--dd", "7"),
        ("--reads", "0"),
    ] {
        let output = simulate(&[option, value]);
        assert_eq!(output.status.code(), Some(1), "{option} {value}");
        assert!(stderr(&output).contains(option), "{}", stderr(&output));
    }
    for (option, value) in [
        ("--rc", "0"),
        ("--rc", "1"),
        ("--shm", "0.05"),
        ("--dd", "1"),
    ] {
        let output = simulate(&[option, value]);
        assert!(
            output.status.success(),
            "{option} {value}: {}",
            stderr(&output)
        );
    }
    fs::remove_dir_all(dir).unwrap();
}

/// evaluate's report goes to a stdout that may be closed or full: that is
/// an error (exit 1), not a panic.
#[test]
fn evaluate_reports_a_failed_write() {
    let dir = workspace("evaluate");
    fs::write(dir.join("truth.tsv"), "sequence_id\tv_call\n").unwrap();
    fs::write(dir.join("calls.tsv"), "sequence_id\tv_call\n").unwrap();
    let args = ["evaluate", "--truth", "truth.tsv", "calls.tsv"];
    let output = run(&dir, &args);
    assert!(output.status.success(), "{}", stderr(&output));
    // A reader that leaves at once.
    let mut child = swig(&dir, &args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    drop(child.stdout.take());
    let output = child.wait_with_output().unwrap();
    assert!(matches!(output.status.code(), Some(0 | 1)), "{output:?}");
    assert!(!stderr(&output).contains("panicked"), "{}", stderr(&output));
    if cfg!(target_os = "linux") {
        let full = fs::File::create("/dev/full").unwrap();
        let output = swig(&dir, &args).stdout(full).output().unwrap();
        assert_eq!(output.status.code(), Some(1), "{}", stderr(&output));
        assert!(
            stderr(&output).contains("No space left"),
            "{}",
            stderr(&output)
        );
    }
    fs::remove_dir_all(dir).unwrap();
}

/// `text` in Latin-1, as some editors save it: `é` is one byte, no UTF-8.
fn latin1(text: &str) -> Vec<u8> {
    text.chars()
        .map(|c| u8::try_from(u32::from(c)).unwrap())
        .collect()
}

/// Germline files are read as the reads are, and as Swig Web reads them: a
/// UTF-8 byte-order mark is skipped, and a byte that is no UTF-8 (a Latin-1
/// header description) reads as U+FFFD. So are evaluate's truth and calls.
#[test]
fn files_are_read_as_utf8_with_replacement() {
    let dir = workspace("encoding");
    let prepare = |out: &str| {
        let args = ["prepare", "--v", "V.fasta", "--j", "J.fasta", "--out", out];
        let output = run(&dir, &args);
        assert!(output.status.success(), "{}", stderr(&output));
        fs::read(dir.join(out).join("V.fasta")).unwrap()
    };
    let plain = prepare("plain");
    let marked = [&b"\xef\xbb\xbf"[..], V.as_bytes()].concat();
    fs::write(dir.join("V.fasta"), marked).unwrap();
    let described = J.replacen("anchor=15", "anchor=15 café", 1);
    fs::write(dir.join("J.fasta"), latin1(&described)).unwrap();
    assert_eq!(prepare("encoded"), plain);
    let args = [
        "simulate", "--v", "V.fasta", "--j", "J.fasta", "--reads", "2",
    ];
    let output = run(&dir, &[&args[..], &["--out", "sim"]].concat());
    assert!(output.status.success(), "{}", stderr(&output));
    let truth = fs::read_to_string(dir.join("sim.truth.tsv")).unwrap();
    let truth = truth.replacen("\nsim0\t", "\nsimé\t", 1);
    fs::write(dir.join("truth.tsv"), latin1(&truth)).unwrap();
    let calls = "sequence_id\tv_call\nsimé\tIGHV1-1*01\n";
    fs::write(dir.join("calls.tsv"), latin1(calls)).unwrap();
    let output = run(&dir, &["evaluate", "--truth", "truth.tsv", "calls.tsv"]);
    assert!(output.status.success(), "{}", stderr(&output));
    // The Latin-1 identifiers match: only sim1 has no row.
    assert!(stderr(&output).contains("1 simulated reads have no row"));
    fs::remove_dir_all(dir).unwrap();
}

/// A command line without an input its mode needs is a usage error (exit 2
/// with the usage), like any other bad command line.
#[test]
fn missing_inputs_are_usage_errors() {
    let dir = workspace("usage");
    for args in [
        &[][..],
        &["--v", "V.fasta", "--j", "J.fasta"],
        &["reads.fasta", "--v", "V.fasta"],
        &["prepare", "--v", "V.fasta", "--j", "J.fasta"],
        &["simulate", "--v", "V.fasta", "--j", "J.fasta"],
        &["evaluate", "--truth", "truth.tsv"],
        &["evaluate", "calls.tsv"],
        &[
            "reads.fasta",
            "--v",
            "V.fasta",
            "--j",
            "J.fasta",
            "--double-d",
            "all",
        ],
    ] {
        let output = run(&dir, args);
        assert_eq!(
            output.status.code(),
            Some(2),
            "{args:?}: {}",
            stderr(&output)
        );
        assert!(stderr(&output).contains("usage:"), "{args:?}");
    }
    // Help on a stdout that cannot take it does not panic either.
    if cfg!(target_os = "linux") {
        let full = fs::File::create("/dev/full").unwrap();
        let output = swig(&dir, &["--help"]).stdout(full).output().unwrap();
        assert!(output.status.success(), "{}", stderr(&output));
    }
    fs::remove_dir_all(dir).unwrap();
}

/// An option that filters or adds to the output without the database it
/// needs is a usage error (exit 2 with the usage), refused before any
/// database is read or the sidecar is created; one that only tunes a stage
/// the run does not have is run, with a warning, and changes nothing.
#[test]
fn options_without_their_database() {
    let dir = workspace("databases");
    let simulate = [
        "--out", "sim", "--reads", "40", "--seed", "3", "--preset", "easy",
    ];
    let output = run(&dir, &[&["simulate"], &GERMLINES[..], &simulate].concat());
    assert!(output.status.success(), "{}", stderr(&output));
    let vj = ["sim.fasta", "--v", "V.fasta", "--j", "J.fasta"];
    let vdj = [&vj[..], &["--d", "D.fasta"]].concat();
    let vjc = [&vj[..], &["--c", "C.fasta"]].concat();
    for (germlines, options, message) in [
        (
            &vdj,
            &["--c-evalue", "0.01"][..],
            "--c-evalue needs a C database (--c)",
        ),
        (
            &vdj,
            &["--c-prefix-identity", "0.97"],
            "--c-prefix-identity needs a C database (--c)",
        ),
        (
            &vjc,
            &["--double-d", "all", "--double-d-out", "dd.tsv"],
            "--double-d needs a D database (--d)",
        ),
        (
            &vjc,
            &["--double-d", "long-span", "--double-d-out", "dd.tsv"],
            "--double-d needs a D database (--d)",
        ),
    ] {
        let output = run(&dir, &[&germlines[..], options].concat());
        let text = stderr(&output);
        assert_eq!(output.status.code(), Some(2), "{options:?}: {text}");
        assert!(text.contains(message), "{options:?}: {text}");
        assert!(text.contains("usage:"), "{options:?}: {text}");
        assert!(!text.contains("database prepared"), "{options:?}: {text}");
        assert!(!dir.join("dd.tsv").exists(), "{options:?}");
    }

    // 0 turns the C-prefix gate off: nothing to refuse or report.
    let output = run(&dir, &[&vdj[..], &["--c-prefix-identity", "0"]].concat());
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(!stderr(&output).contains("c-prefix"), "{}", stderr(&output));

    // D tuning and the screen's rules without a D database or the screen.
    let plain = run(&dir, &vj);
    assert!(plain.status.success(), "{}", stderr(&plain));
    let tuning = [
        "--min-d-match",
        "6",
        "--d-mismatch",
        "-4",
        "--no-anchor-window",
        "--double-d-seed-length",
        "9",
    ];
    let tuned = run(&dir, &[&vj[..], &tuning].concat());
    let text = stderr(&tuned);
    assert!(tuned.status.success(), "{text}");
    assert_eq!(tuned.stdout, plain.stdout);
    assert!(
        text.contains(
            "swig: warning: --min-d-match, --d-mismatch, --no-anchor-window have no \
             effect without a D database (--d)\n"
        ),
        "{text}"
    );
    assert!(
        text.contains(
            "swig: warning: --double-d-seed-length has no effect without --double-d all \
             or long-span\n"
        ),
        "{text}"
    );
    assert!(!stderr(&plain).contains("no effect"), "{}", stderr(&plain));
    fs::remove_dir_all(dir).unwrap();
}

/// Fixture (a) of `tests/open_flank.rs` on this database: each full V
/// allele followed by 100 random bases (a fixed linear congruential stream).
fn v_then_random() -> String {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut base = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        b"ACGT"[(state >> 62) as usize] as char
    };
    let v: Vec<&str> = V.lines().filter(|l| !l.starts_with('>')).collect();
    (0..20)
        .map(|i| {
            let tail: String = (0..100).map(|_| base()).collect();
            format!(">vcut{i}\n{}{tail}\n", v[i % v.len()])
        })
        .collect()
}

/// The `d_call` column of an AIRR TSV.
fn d_column(tsv: &[u8]) -> Vec<String> {
    let text = String::from_utf8_lossy(tsv);
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let d = header.iter().position(|&c| c == "d_call").unwrap();
    lines
        .map(|l| l.split('\t').nth(d).unwrap_or("").to_string())
        .collect()
}

/// FNV-1a of Swig 0.39.3's AIRR rows on [`v_then_random`].
const FNV_0393: u64 = 0x9186_e1b9_be7f_d846;

/// (e) `--d-open-flank off` turns the open-flank cost off and gives Swig
/// 0.39.3's rows on fixture (a), byte for byte: a chance D on all 20
/// reads, where the cost leaves 2.
#[test]
fn the_open_flank_cost_can_be_turned_off() {
    let dir = workspace("open-flank");
    fs::write(dir.join("reads.fa"), v_then_random()).unwrap();
    let call = |extra: &[&str]| {
        let mut args = GERMLINES.to_vec();
        args.extend_from_slice(extra);
        args.push("reads.fa");
        let output = run(&dir, &args);
        assert!(output.status.success(), "{}", stderr(&output));
        output.stdout
    };
    let (on, off) = (call(&[]), call(&["--d-open-flank", "off"]));
    // Swig 0.39.3's output on these reads.
    assert_eq!(fnv(&off), FNV_0393);
    let count = |tsv: &[u8]| d_column(tsv).iter().filter(|d| !d.is_empty()).count();
    assert_eq!((count(&on), count(&off)), (2, 20));
    assert_eq!(call(&["--d-open-flank", "on"]), on);
    let bad = run(&dir, &["--d-open-flank", "maybe", "reads.fa"]);
    assert_eq!(bad.status.code(), Some(2));
    assert!(
        stderr(&bad).contains("--d-open-flank expects on or off"),
        "{}",
        stderr(&bad)
    );
    fs::remove_dir_all(&dir).unwrap();
}
