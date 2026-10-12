//! The `swig` command line: its subcommands' outputs, refusals and exit
//! codes (0 success, 1 a failure while running, 2 a bad command line).

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};

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
