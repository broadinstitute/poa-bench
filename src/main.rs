use std::{fs, io, process, fmt};
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::path::{Path, PathBuf};
use std::string::FromUtf8Error;
use clap::{command, Parser, Subcommand, Args};
use walkdir::WalkDir;
use anyhow::{Result, Context};
use rayon::prelude::*;

use crate::dataset::DatasetConfig;

mod dataset;

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
struct Dataset(String, PathBuf, DatasetConfig);

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
