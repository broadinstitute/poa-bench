use std::{fs, io, process, fmt};
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::path::{Path, PathBuf};
use std::string::FromUtf8Error;
use clap::{command, Parser, Subcommand, Args};
use walkdir::WalkDir;
use anyhow::{Result, Context};
use rayon::prelude::*;
use core_affinity;
use rayon::ThreadPoolBuilder;

use crate::dataset::DatasetConfig;
use crate::jobs::{Algorithm, Job};

mod dataset;
mod jobs;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct CliParser {
    #[command(subcommand)]
    subcommand: CliSubcommand,
}

#[derive(Subcommand, Debug, Clone)]
enum CliSubcommand {
    /// Main entry point for the benchmark runner
    Bench(BenchArgs),

    /// Runs an individual algorithm on a data set
    RunJob, // TODO: runner args
}

#[derive(Args, Debug, Clone)]
struct BenchArgs {
    /// Specify which algorithms to run, options include poasta and spoa. To specify multiple,
    /// separate them by spaces
    #[clap(value_enum, short, long, num_args=0..)]
    algorithms: Vec<Algorithm>,

    /// Number of parallel threads to start. This number will include the main orchestrator
    /// process, so the number of actual worker threads will be one less than the number specified.
    #[clap(short = 'j', long, default_value="2")]
    parallel: usize,

    #[clap(short, long, default_value="data/")]
    dataset_dir: PathBuf,

    #[clap(short, long, default_value="output/")]
    output_dir: PathBuf,
}

#[derive(Debug)]
enum POABenchError {
    IOError(io::Error),
    Utf8Error(FromUtf8Error),
    BuildGraphError(String),
}

impl Error for POABenchError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        None
    }
}

impl Display for POABenchError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IOError(source) =>
                fmt::Display::fmt(source, f),
            Self::Utf8Error(source) =>
                fmt::Display::fmt(source, f),
            Self::BuildGraphError(stderr) =>
                write!(f, "Could not build graph for data set! Stderr: {}", stderr)
        }
    }
}

impl From<io::Error> for POABenchError {
    fn from(value: io::Error) -> Self {
        Self::IOError(value)
    }
}

impl From<FromUtf8Error> for POABenchError {
    fn from(value: FromUtf8Error) -> Self {
        Self::Utf8Error(value)
    }
}

/// A type that combines the data set name (inferred from its path) and the configuration (read
/// from the TOML file)
#[derive(Clone, Debug)]
pub struct Dataset(String, PathBuf, DatasetConfig);

impl Dataset {
    fn output_dir(&self, base_output_path: &Path) -> PathBuf {
        base_output_path.join(&self.0)
    }

    fn graph_output_fname(&self, base_dir: &Path) -> PathBuf {
        self.output_dir(base_dir).join("graph.poasta")
    }

    fn graph_sequences_fname(&self) -> PathBuf {
        self.1.join(&self.2.graph_set.fname)
    }

    fn align_sequences_fname(&self) -> PathBuf {
        self.1.join(&self.2.align_set.fname)
    }
}


fn build_graph_with_poasta(output_dir: &Path, dataset: &Dataset) -> Result<(), POABenchError> {
    let seq_fname = dataset.graph_sequences_fname();
    let output_fname = dataset.graph_output_fname(output_dir);
    eprintln!("{:?}, {:?}", seq_fname, output_fname);

    if output_fname.exists() {
        let seq_meta = fs::metadata(seq_fname)?;
        let graph_meta = fs::metadata(output_fname)?;

        if seq_meta.modified()? <= graph_meta.modified()? {
            eprintln!("Found up-to-date graph for dataset: {}", dataset.0);
            return Ok(());
        }
    }

    eprintln!("Creating graph for dataset: {}", dataset.0);

    let process = process::Command::new("poasta")
        .arg("align")
        .arg("-o")
        .arg(dataset.graph_output_fname(output_dir))
        .arg(dataset.graph_sequences_fname())
        .output()?;

    if process.status.success() {
        Ok(())
    } else {
        Err(POABenchError::BuildGraphError(String::from_utf8(process.stderr)?))
    }
}


fn build_graphs(output_dir: &Path, datasets: &[Dataset]) -> Result<()> {
    let results: Vec<_> = datasets.par_iter()
        .map(|dataset| {
            build_graph_with_poasta(output_dir, dataset)
        })
        .collect();

    for (dataset, r) in datasets.iter().zip(results) {
        r.with_context(|| format!("Dataset: {:?}", dataset.0))?
    }

    Ok(())
}


fn run_job(job: &Job) {
    eprintln!("JOB algorithm: {:?}, dataset: {:?}", job.algorithm, job.dataset.0);
}


fn bench(bench_args: &BenchArgs) -> Result<()> {
    eprintln!("Finding data sets...");
    let mut datasets = Vec::new();
    for entry in WalkDir::new(&bench_args.dataset_dir)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| !e.file_type().is_dir())
    {
        let fname = entry.file_name().to_string_lossy();

        if fname.ends_with(".toml") {
            let contents = fs::read_to_string(entry.path())
                .with_context(|| format!("Could not read file {:?}!", entry.path()))?;
            let dataset_cfg: DatasetConfig = toml::from_str(&contents)
                .with_context(|| format!("Could not parse data set TOML! Filename: {:?}", entry.path()))?;

            let dataset_dir = entry.path().parent().unwrap();
            let dataset_name = dataset_dir
                .strip_prefix(&bench_args.dataset_dir)?
                .to_string_lossy().to_string();

            eprintln!("Found dataset: {:?}", dataset_name);
            datasets.push(Dataset(dataset_name, dataset_dir.to_owned(), dataset_cfg))
        }
    }

    eprintln!("Building graphs...");
    build_graphs(&bench_args.output_dir, &datasets)?;

    eprintln!("Building job list...");
    let algorithms: &[Algorithm] = if bench_args.algorithms.len() > 0 { &bench_args.algorithms } else { jobs::ALL_ALGORITHMS };

    let mut jobs = Vec::new();
    for algorithm in algorithms {
        for dataset in &datasets {
            jobs.push(Job {
                algorithm: *algorithm,
                dataset,
            })
        }
    }

    let mut cores = core_affinity::get_core_ids()
        .unwrap()
        .into_iter()
        .take(std::cmp::max(2, bench_args.parallel));

    let orchestrator_core = cores.next().unwrap();
    core_affinity::set_for_current(orchestrator_core);

    let worker_cores: Vec<_> = cores.collect();

    std::thread::scope(|scope| {
        let pool = ThreadPoolBuilder::new()
            .num_threads(worker_cores.len())
            .spawn_handler(|worker| {
                // Configure thread to run
                let mut b = std::thread::Builder::new();
                if let Some(name) = worker.name() {
                    b = b.name(name.to_owned());
                }
                if let Some(stack_size) = worker.stack_size() {
                    b = b.stack_size(stack_size);
                }

                // Spawn the thread, and pin it to desired core
                b.spawn_scoped(scope, || {
                    core_affinity::set_for_current(worker_cores[worker.index()]);

                    worker.run()
                })?;

                Ok(())
            })
            .build()?;

        for job in &jobs {
            pool.install(|| run_job(job));
        }

        anyhow::Ok(())
    })?;

    Ok(())
}


fn main() -> Result<()> {
    let args = CliParser::parse();

    match args.subcommand {
        CliSubcommand::Bench(ref bench_args) => bench(bench_args)?,
        CliSubcommand::RunJob => eprintln!("Not implemented yet!")
    }

    Ok(())
}
