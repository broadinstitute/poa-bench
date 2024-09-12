<h1 align="center">📈 POA-bench</h1>
<h2 align="center">Benchmarking various partial order aligners</h2>

<p>&nbsp;</p>

POA-bench is a tool to benchmark partial order aligners. We constructed multiple datasets from bacterial housekeeping
genes, the human HLA locus, and others, with varying graph sizes and sequence diversity. POA-bench assesses the
runtime, throughput and memory usage of multiple partial order alignment tools. POA-bench uses thin wrappers around
each tool's internal `align` function and thus benchmarks only the time spent aligning sequence, and not anything I/O or
graph modification related.

## Installation

POA-bench only works on Linux.

### Pre-built binaries

TODO

### Conda

TODO

### Building poa-bench from source

#### Rust compiler

POA-bench is written in Rust, and to build and install it, you'll need a recent version of the Rust compiler. The
minimum supported Rust version is 1.70.

1. Download and install `rustup`: https://rustup.rs/
2. Run `rustup update`

#### Building poa-bench

1. Clone the repository.

   ```bash
   git clone https://github.com/broadinstitute/poa-bench
   ```

2. Clone `poasta`, a dependency for `poa-bench` that should reside in the same folder, which is not yet publicly
   available.

   ```bash
   git clone https://github.com/broadinstitute/poasta
   ```

3. Move into the directory.

   ```bash
   cd poa-bench
   ```

4. Build using `cargo`. We enable a flag to ensure the compiler uses all features of your machine's CPU.
   To maximize portability of the binary, however, remove the `RUSTFLAGS="..."` part.

   ```bash
   RUSTFLAGS="-C target-cpu=native" cargo build --release
   ```

5. The built `poa-bench` executable will be available in the directory `target/release/`

### Installing the Python helper tools

`poa-bench` relies on a couple of helper scripts written in Python to
create and manage datasets. These are located in the `python/` folder.
To install them, move to the `python/` directory and run

```bash
pip install .
```

I'd recommend to create a new conda environment first, with at least Python 3.10.

### Additional dependencies

The following additional tools should be available in your `$PATH`:

* `samtools`
* `spoa`

## Datasets

Datasets are defined by a TOML configuration files named `meta.toml`. See the data/ directory for a number of examples. Each
dataset comprises a set of sequences used to construct a graph, and a set of sequences used to benchmark alignment to the
constructed graph. Each FASTA with sequences **must** be compressed with gzip.
To define which sequences to use, specify it as follows in the TOML file:

```toml
[graph_set]
fname = "file_with_graph_seq.fa.gz"

[align_set]
fname = "file_with_aln_seq.fa.gz"
```

The above fields are the only mandatory fields in the TOML file, but you can add other metadata if you want.

## Running benchmarks

To run benchmarks, invoke `poa-bench bench`. This will run the benchmarks for all supported algorithms on all data sets.
To find datasets, it will recursively search a directory `data/` in the current
working directory to find dataset configuration files. Each combination of a dataset and algorithm
is a single job, and `poa-bench` will spawn a new worker process that will perform the benchmark for that job.
To reduce variance in benchmarks, each worker process gets pinned to a specific core. If requested, `poa-bench`
can run multiple jobs in parallel, using the `-j` option. All output files and results are written to a directory
`output/`. The algorithms to run, and the input/output directories are configurable on the command line.

```
Usage: poa-bench bench [OPTIONS]

Options:
  -a, --algorithms [<ALGORITHMS>...]  Specify which algorithms to run, options include poasta and spoa. To specify multiple, separate them by spaces [possible values: poasta, spoa]
  -j, --parallel <PARALLEL>           Number of parallel threads to start. This number will include the main orchestrator process, so the number of actual worker threads will be one less than the number specified [default: 2]
  -d, --datasets-dir <DATASETS_DIR>   [default: data/]
  -o, --output-dir <OUTPUT_DIR>       [default: output/]
  -h, --help                          Print help
  -V, --version                       Print version
```

## Acknowledgements

This program is inspired (and shares some code with) the excellent [pa-bench](https://github.com/pairwise-alignment/pa-bench)
repository.

## Related repositories

* [pa-bench](https://github.com/pairwise-alignment/pa-bench) - Benchmark suite for pairwise aligners
* [poasta](https://github.com/broadinstitute/poasta) - POASTA aligner
* [spoa](https://github.com/rvaser/spoa) - SIMD partial order aligner
* [spoa-rs](https://github.com/broadinstitute/spoa-rs) - Rust bindings to SPOA
* [abPOA](https://github.com/yangao07/abPOA) - adaptive band partial order aligner
* [abpoa-rs](https://github.com/broadinstitute/abpoa-rs) - Rust bindings to abPOA

