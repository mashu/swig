# GenAIRR benchmark: swig, IgBLAST and MiXCR

Simulates reads with [GenAIRR](https://github.com/MuteJester/GenAIRR)
(2.3.1), runs swig, IgBLAST and MiXCR against the same GenAIRR reference,
and scores each tool's V, D and J calls against GenAIRR's truth. Results
and analysis: [`BENCHMARK_GENAIRR_0.39.6.md`](../../BENCHMARK_GENAIRR_0.39.6.md).

## Scoring

GenAIRR's `v_call`/`d_call`/`j_call` is a tie set. Each allele gets one
point per read base that matches it at the true allele's germline index
over the segment's true read region (plus NP bases while they narrow the
tie). The set holds every allele tied at the maximum. A tool's call C
against truth T scores

    |T ∩ C| / |T ∪ C|        (both empty: 1; one empty: 0)

per read and segment, so one correct allele against a two-allele tie
scores 1/2, as does a two-allele call against a one-allele truth.

## Setup

```sh
python3 -m venv venv && venv/bin/pip install GenAIRR==2.3.1 numpy
cargo build --release --manifest-path ../../engine/Cargo.toml --bin swig

# IgBLAST 1.22.0 from bioconda (NCBI's FTP may be unreachable):
mkdir -p igblast && cd igblast
curl -O https://conda.anaconda.org/bioconda/linux-64/igblast-1.22.0-h6a68c12_1.tar.bz2
curl -O https://conda.anaconda.org/bioconda/linux-64/ncbi-vdb-3.1.1-h9948957_3.tar.bz2
mkdir root && for f in *.tar.bz2; do tar -xjf $f -C root; done; cd ..
export IGBLAST_ROOT=$PWD/igblast/root

# MiXCR 3.0.13 from source: the last release without a license key.
git clone --depth 1 --branch v3.0.13 https://github.com/milaboratory/mixcr
(cd mixcr && mvn -q -DskipTests package && mvn -q dependency:copy-dependencies -DincludeScope=runtime -DoutputDirectory=target/deps)
curl -o mixcr/target/deps/jcommander-1.72.jar https://repo1.maven.org/maven2/com/beust/jcommander/1.72/jcommander-1.72.jar
export MIXCR_DIR=$PWD/mixcr
```

## Run

```sh
PYTHON=venv/bin/python ./run_matrix.sh 10000 7          # simulate, run tools -> work/
# Prototype post-processing of swig's output (positional J, D posterior):
PYTHON=venv/bin/python ./improve.sh work/ref_igh work/runs/igh_shm5 parametric:1.0 HUMAN_IGH_OGRDB parametric
python3 merge_calls.py work/runs/igh_shm5.recall.tsv work/runs/igh_shm5.dpost-parametric.tsv work/runs/igh_shm5.final-parametric.tsv
python3 genairr_bench.py score work/runs/igh_shm5.truth.tsv swig=work/runs/igh_shm5.swig.tsv \
  igblast=work/runs/igh_shm5.igblast.tsv final=work/runs/igh_shm5.final-parametric.tsv
python3 loss.py work/runs/igh_shm5.truth.tsv swig=work/runs/igh_shm5.swig.tsv   # where accuracy is lost
```

| File | Role |
|---|---|
| `genairr_bench.py` | `export` GenAIRR's reference for the three tools; `simulate` scenarios; `score` |
| `run_tools.sh`, `run_matrix.sh` | run swig, IgBLAST (`-outfmt 19`) and MiXCR (`align`, `exportAlignments`) |
| `loss.py` | split each tool's loss into missing/spurious calls, under/over-reported ties, wrong allele, wrong gene |
| `recall.py` | prototype: GenAIRR-style positional tie sets at swig's placement, through its CIGAR |
| `d_posterior.py` | prototype: D call from a posterior over allele, offset and trims, reporting the set of best expected score |
| `improve.sh`, `merge_calls.py`, `summary.py` | run the prototypes, combine their calls, tabulate |
| `dposterior/` | the D posterior in Rust, with switches for priors, decision rule and offsets |
| `ablate.sh` | score each `dposterior` variant and swig's `--d-allele max-marginal` |

The Python D posterior is slow, about 0.2 s per read. `dposterior/` runs
the same model at 0.24–0.31 ms per read. Neither is how the engine runs
it. The ablation found that the gain is the trim and N-length priors,
which swig 0.39.7 adds to its placement DAG as terms on its cuts, for
0.6 % more instructions.

```sh
cargo build --release --manifest-path dposterior/Cargo.toml
./ablate.sh work/ref_igh work/runs/igh_shm5
```
