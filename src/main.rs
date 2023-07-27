use std::fs;
use std::path::PathBuf;
use clap::{command, Parser, Subcommand, Args};
use walkdir::WalkDir;
use anyhow::{Result, Context};

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


fn bench(bench_args: &BenchArgs) -> Result<()> {
    for entry in WalkDir::new(&bench_args.dataset_dir)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| !e.file_type().is_dir())
    {
        let fname = entry.file_name().to_string_lossy();

        if fname.ends_with(".toml") {
            let contents = fs::read_to_string(entry.path())
                .with_context(|| format!("Could not read file {:?}!", entry.path()))?;
            let dataset_meta: dataset::Dataset = toml::from_str(&contents)
                .with_context(|| format!("Could not parse data set TOML! Filename: {:?}", entry.path()))?;

            let dataset_name = entry.path().strip_prefix(&bench_args.dataset_dir)?;

            eprintln!("Found dataset: {:?}, graph seq: {:?}", dataset_name, dataset_meta.graph_set.fname);
        }
    }

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
