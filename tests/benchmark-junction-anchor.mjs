// Paired before/after check of the clipped-Cys104 recovery (Swig 0.39.6).
//
//   node tests/benchmark-junction-anchor.mjs OLD_SWIG NEW_SWIG [READS] [WORKDIR]
//
// OLD_SWIG and NEW_SWIG are native `swig` binaries (e.g. 0.39.5 and 0.39.6
// release builds). For each germline set and preset the new binary
// simulates READS reads (default 5,000, seed 765127); both binaries annotate
// them with Swig Web's identity floor, and every AIRR column is compared.
// A changed row must be a read that gains a junction and differ only in
// junction-derived columns; recovered junctions are checked against the
// simulator's templated junction, and `swig evaluate` reports junction
// exactness and the spurious-junction rate for both runs.
import { execFileSync } from "node:child_process";
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync, copyFileSync, existsSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { gunzipSync } from "node:zlib";

const [oldBinary, newBinary, readsArgument, workArgument] = process.argv.slice(2);
if (!oldBinary || !newBinary) {
  console.error("usage: node tests/benchmark-junction-anchor.mjs OLD_SWIG NEW_SWIG [READS] [WORKDIR]");
  process.exit(2);
}
const reads = Number(readsArgument ?? 5000);
const work = workArgument ?? mkdtempSync(join(tmpdir(), "swig-junction-anchor-"));
const root = new URL("..", import.meta.url);

// Germline sets: [directory, species in the IMGT pack, locus, chain]; KIMDB
// macaque IGH comes from its bundled FASTA files.
const SETS = [
  ["homo_igh", "Homo sapiens", "IGH", "heavy"],
  ["homo_igk", "Homo sapiens", "IGK", "light"],
  ["homo_igl", "Homo sapiens", "IGL", "light"],
  ["homo_trb", "Homo sapiens", "TRB", "heavy"],
  ["vicugna_igh", "Vicugna pacos", "IGH", "heavy"],
  ["macaca_igh", "Macaca mulatta_AG07107", "IGH", "heavy"],
  ["kimdb_igh", null, "IGH", "heavy"],
  ["mus_igh", "Mus musculus_BALB/c", "IGH", "heavy"],
];
const PRESETS = ["mid", "hard"];
const JUNCTION_FIELDS = new Set([
  "junction", "junction_aa", "junction_length", "junction_aa_length",
  "cdr3", "cdr3_aa", "cdr3_start", "cdr3_end",
  "productive", "vj_in_frame", "complete_vdj",
  "sequence_frame", "sequence_aa", "sequence_alignment_aa", "germline_alignment_aa",
  "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "fwr4_aa", "d_frame",
]);

const pack = JSON.parse(gunzipSync(readFileSync(new URL("public/references/imgt-202632-7-swig-0.7.json.gz", root))));

function writeReferences(directory, species, locus) {
  mkdirSync(directory, { recursive: true });
  if (!species) {
    for (const segment of ["V", "D", "J"]) {
      copyFileSync(new URL(`public/references/kimdb-1.1/Macaca_mulatta/IGH/${segment}.fasta`, root), join(directory, `${segment}.fasta`));
    }
    return;
  }
  const alleles = pack.species.find((entry) => entry.name === species).loci[locus];
  for (const segment of ["V", "D", "J"]) {
    if (!alleles[segment]?.length) continue;
    const fasta = alleles[segment]
      .map(([name, sequence, metadata]) => `>${name}${metadata ? ` SWIGMETA=${metadata.join(",")}` : ""}\n${sequence}\n`)
      .join("");
    writeFileSync(join(directory, `${segment}.fasta`), fasta);
  }
}

function databaseArguments(directory) {
  const args = ["--v", join(directory, "V.fasta"), "--j", join(directory, "J.fasta")];
  if (existsSync(join(directory, "D.fasta"))) args.push("--d", join(directory, "D.fasta"));
  return args;
}

function table(text) {
  const [header, ...lines] = text.trimEnd().split("\n");
  const columns = header.split("\t");
  return lines.map((line) => {
    const values = line.split("\t");
    return Object.fromEntries(columns.map((column, index) => [column, values[index] ?? ""]));
  });
}

function annotate(binary, directory, fasta) {
  return execFileSync(binary, [...databaseArguments(directory), "--min-identity", "0.6", "--threads", "4", fasta], {
    encoding: "utf8", maxBuffer: 1 << 30, stdio: ["ignore", "pipe", "ignore"],
  });
}

function evaluate(truth, calls) {
  const tsv = join(work, "evaluate.tsv");
  writeFileSync(tsv, calls);
  const report = execFileSync(newBinary, ["evaluate", "--truth", truth, tsv], { encoding: "utf8", stdio: ["ignore", "pipe", "ignore"] });
  return report.split("\n", 1)[0].replace(/^.*junction exact/, "exact");
}

const totals = { reads: 0, changed: 0, recovered: 0, templatedEqual: 0, templatedDiffer: 0, rebuilt: 0, productive: 0, otherColumns: 0, previouslyResolved: 0 };
console.log("set          preset  reads  recovered  templated=  templated≠  rebuilt  productive  other-columns  resolved-changed");
for (const [name, species, locus, chain] of SETS) {
  const directory = join(work, "refs", name);
  writeReferences(directory, species, locus);
  for (const preset of PRESETS) {
    const prefix = join(work, `${name}.${preset}`);
    execFileSync(newBinary, ["simulate", ...databaseArguments(directory), "--chain", chain, "--out", prefix,
      "--preset", preset, "--reads", String(reads), "--seed", "765127"], { stdio: "ignore" });
    const oldText = annotate(oldBinary, directory, `${prefix}.fasta`);
    const newText = annotate(newBinary, directory, `${prefix}.fasta`);
    const [before, after] = [table(oldText), table(newText)];
    const truth = new Map(table(readFileSync(`${prefix}.truth.tsv`, "utf8")).map((row) => [row.sequence_id, row]));
    if (before.length !== reads || after.length !== reads) throw new Error(`${name} ${preset}: row count`);
    const count = { changed: 0, recovered: 0, templatedEqual: 0, templatedDiffer: 0, rebuilt: 0, productive: 0, otherColumns: 0, previouslyResolved: 0 };
    before.forEach((old, index) => {
      const row = after[index];
      if (old.sequence_id !== row.sequence_id) throw new Error(`${name} ${preset}: row order`);
      const changed = Object.keys(old).filter((column) => old[column] !== row[column]);
      if (!changed.length) return;
      count.changed += 1;
      const other = changed.filter((column) => !JUNCTION_FIELDS.has(column));
      if (other.length) {
        count.otherColumns += 1;
        console.log(`  ${row.sequence_id}: non-junction columns changed: ${other.join(", ")}`);
      }
      if (old.junction || !row.junction) {
        count.previouslyResolved += 1;
        console.log(`  ${row.sequence_id}: a row with a junction before (or none after) changed: ${changed.join(", ")}`);
        return;
      }
      count.recovered += 1;
      if (row.productive === "T") count.productive += 1;
      const templated = truth.get(row.sequence_id)?.junction;
      if (!templated) count.rebuilt += 1;
      else if (templated === row.junction) count.templatedEqual += 1;
      else {
        count.templatedDiffer += 1;
        console.log(`  ${row.sequence_id}: truth ${templated}, called ${row.junction}`);
      }
    });
    for (const key of Object.keys(count)) totals[key] += count[key];
    totals.reads += reads;
    console.log(`${name.padEnd(12)} ${preset.padEnd(6)} ${String(reads).padStart(6)} ${String(count.recovered).padStart(10)} ${String(count.templatedEqual).padStart(11)} ${String(count.templatedDiffer).padStart(11)} ${String(count.rebuilt).padStart(8)} ${String(count.productive).padStart(11)} ${String(count.otherColumns).padStart(14)} ${String(count.previouslyResolved).padStart(17)}`);
    console.log(`  junctions before: ${evaluate(`${prefix}.truth.tsv`, oldText)}`);
    console.log(`  junctions after:  ${evaluate(`${prefix}.truth.tsv`, newText)}`);
  }
}
console.log(`totals: ${JSON.stringify(totals)}`);
console.log(`work directory: ${work}`);
