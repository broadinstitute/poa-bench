use std::{fs, process, thread};
use std::fmt::Debug;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use clap::{Args, command, Parser, Subcommand};
use anyhow::{Context, Result};
use rayon::prelude::*;
use core_affinity;
use core_affinity::CoreId;

use errors::POABenchError;
use crate::dataset::{Dataset, find_datasets};
use crate::jobs::{Algorithm, BenchmarkType, JobResult};

mod errors;
mod worker;
mod dataset;
mod jobs;
mod bench;

#[derive(Parser, Debug, Clone)]
struct CliArgs {
    #[command(subcommand)]
    command: Command
}

#[derive(Subcommand, Debug, Clone)]
enum Command {
    Bench(BenchArgs),
    Worker(worker::WorkerArgs),
}


#[derive(Args, Debug, Clone)]
#[command(author, version, about)]
struct BenchArgs {
    /// Specify which algorithms to run, options include poasta and spoa. To specify multiple,
    /// separate them by spaces.
    #[clap(value_enum, short, long, num_args=0..)]
    algorithms: Vec<Algorithm>,

    /// Specify what benchmark types to run (single sequence vs full MSA). To specify multiple,
    /// separate them by spaces.
    #[clap(value_enum, short, long, num_args=0..)]
    benchmarks: Vec<BenchmarkType>,

    /// Number of parallel threads to start. This number will include the main orchestrator
    /// process, so the number of actual worker threads will be one less than the number specified.
    #[clap(short = 'j', long, default_value="2")]
    parallel: usize,

    #[clap(short, long, default_value="data/")]
    datasets_dir: PathBuf,

    #[clap(short, long, default_value="output/")]
    output_dir: PathBuf,

    /// Filename prefix for final results TSV. poa-bench will output two files:
    /// 1) {prefix}.single_seq.tsv with single sequence benchmark results, and 2) {prefix}.full_msa.tsv
    /// with full MSA benchmark results. Default prefix is "results".
    ///
    /// The filename prefix should *not* contain any directories, use the `output_dir` setting for that.
    #[clap(short='f', long, default_value="results")]
    results_prefix: PathBuf,
}


fn build_graph_with_spoa(output_dir: &Path, dataset: &Dataset) -> Result<()> {
    let seq_fname = dataset.graph_sequences_fname();
    let output_fname = dataset.graph_msa_fname(output_dir);

    if output_fname.exists() {
        let seq_meta = fs::metadata(seq_fname)?;
        let graph_meta = fs::metadata(&output_fname)?;

        if seq_meta.modified()? <= graph_meta.modified()? {
            eprintln!("Found up-to-date graph for dataset: {}", dataset.name());
            return Ok(());
        }
    }

    eprintln!("Creating graph with SPOA for dataset: {}", dataset.name());
    fs::create_dir_all(dataset.output_dir(output_dir))?;

    let process = process::Command::new("spoa")
        .args(&["-l", "1", "-m", "0", "-n", "-4", "-g", "-8", "-e", "-2", "-q", "0", "-c", "0", "-r", "1"])
        .arg(dataset.graph_sequences_fname())
        .output()
        .context("Could not run SPOA, is it installed and available in $PATH?")?;

    if process.status.success() {
        let mut ofile = File::create(&output_fname)
            .map(BufWriter::new)?;
        ofile.write_all(&process.stdout)
            .context("Could not write graph to file!")?;

        Ok(())
    } else {
        Err(POABenchError::BuildGraphError(String::from_utf8(process.stderr)?).into())
    }

}


fn build_graphs(output_dir: &Path, datasets: &[Dataset]) -> Result<()> {
    let results: Vec<_> = datasets.par_iter()
        .map(|dataset| build_graph_with_spoa(output_dir, dataset))
        .collect();

    for (_, r) in datasets.iter().zip(results) {
        r?
    }

    Ok(())
}


fn sort_sequences_by_genetic_distance(output_dir: &Path, dataset: &Dataset) -> Result<()> {
    let graph_seq_fname = dataset.graph_sequences_fname();
    let align_seq_fname = dataset.align_sequences_fname();
    let output_fname = dataset.combined_sorted_fname(output_dir);

    if output_fname.exists() {
        let graph_seq_meta = fs::metadata(graph_seq_fname)
            .context("Can't obtain graph sequences FASTA metadata")?;
        let align_seq_meta = fs::metadata(align_seq_fname)
            .context("Can't obtain align sequences FASTA metadata")?;
        let output_meta = fs::metadata(output_fname)
            .context("Can't obtain combined and sorted FASTA metadata")?;

        let output_modified = output_meta.modified()
            .context("Can't obtain last modified time for combined and sorted FASTA")?;

        if graph_seq_meta.modified()? <= output_modified && align_seq_meta.modified()? <= output_modified {
            eprintln!("Found up-to-date combined and sorted FASTA for dataset: {}", dataset.name());
            return Ok(());
        }
    }

    fs::create_dir_all(dataset.output_dir(output_dir))?;
    let mut fname_without_gz = dataset.combined_sorted_fname(output_dir);
    fname_without_gz.set_extension("");

    let process = process::Command::new("python")
        .args(&["poa-bench-tools", "sort_fasta"])
        .arg(dataset.graph_sequences_fname())
        .arg(dataset.align_sequences_fname())
        .arg("-o")
        .args(&fname_without_gz)
        .output()
        .context("Could not sort FASTA file because poa-bench-tools command is not available!")?;

    if !process.status.success() {
        return Err(POABenchError::BuildGraphError(String::from_utf8(process.stderr)?).into())
    }

    let process = process::Command::new("gzip")
        .arg(&fname_without_gz)
        .output()
        .context("Could not gzip-compress combined and sorted FASTA!")?;

    if !process.status.success() {
        return Err(POABenchError::BuildGraphError(String::from_utf8(process.stderr)?).into())
    }

    Ok(())
}

fn create_sorted_fastas(output_dir: &Path, datasets: &[Dataset]) -> Result<()> {
    let results: Vec<_> = datasets.par_iter()
        .map(|dataset| sort_sequences_by_genetic_distance(output_dir, dataset))
        .collect();

    for (_, r) in datasets.iter().zip(results) {
        r?
    }

    Ok(())
}

fn run_job<'scope, 'env>(
    scope: &'scope thread::Scope<'scope, 'env>,
    proc_exe: &'scope Path,
    datasets_dir: &'scope Path,
    output_dir: &'scope Path,
    algorithm: Algorithm,
    benchmark_type: BenchmarkType,
    dataset: &'scope Dataset,
    core: Option<CoreId>,
    tx: mpsc::Sender<JobResult>
)
where
    'env: 'scope
{
    eprintln!("STARTING JOB algorithm: {:?}, dataset: {:?}", algorithm, dataset.name());

    scope.spawn(move || {
        let mut command = process::Command::new(proc_exe);

        command
            .arg("worker")
            .arg("-d")
            .arg(datasets_dir)
            .arg("-o")
            .arg(output_dir);

        if let Some(core_id) = core {
            command.arg("-c")
                .arg(format!("{}", core_id.id));
        }

        command
            .arg(dataset.name())
            .arg(algorithm.to_str())
            .arg(benchmark_type.to_str())
            .stdout(process::Stdio::piped());

        eprintln!("Running command: {:?}", command);
        let mut child = command.spawn()?;

        let reader = BufReader::new(child.stdout.as_mut().unwrap());

        for line in reader.lines() {
            if let Ok(contents) = line {
                let job_result: serde_json::Result<JobResult> = serde_json::from_str(&contents);
                match job_result {
                    Ok(result) => tx.send(result)?,
                    Err(e) => eprintln!("ERROR Could not parse job result from line: {}\nERROR {}", contents, e)
                }
            }
        }

        if !child.wait().expect("Could not launch worker!").success() {
            tx.send(JobResult::Error)?
        }

        Ok::<(), POABenchError>(())
    });
}


fn bench(bench_args: BenchArgs) -> Result<()> {
    eprintln!("Finding data sets...");
    let datasets = find_datasets(&bench_args.datasets_dir)?;

    eprintln!("Building graphs...");
    build_graphs(&bench_args.output_dir, &datasets)?;

    eprintln!("Creating sorted FASTAs...");
    create_sorted_fastas(&bench_args.output_dir, &datasets)?;

    eprintln!("Building job list...");
    let algorithms: &[Algorithm] = if bench_args.algorithms.len() > 0 { &bench_args.algorithms } else { jobs::ALL_ALGORITHMS };
    let benchmarks: &[BenchmarkType] = if bench_args.algorithms.len() > 0 { &bench_args.benchmarks } else { jobs::ALL_BENCHMARK_TYPES };

    let mut jobs = Vec::new();
    for algorithm in algorithms {
        for benchmark in benchmarks {
            for dataset in &datasets {
                jobs.push((algorithm, benchmark, dataset))
            }
        }
    }

    let proc_exe = std::env::current_exe()
        .expect("Could not obtain executable name!");

    let mut cores = core_affinity::get_core_ids()
        .unwrap()
        .into_iter()
        .take(std::cmp::max(2, bench_args.parallel));

    let orchestrator_core = cores.next().unwrap();
    core_affinity::set_for_current(orchestrator_core);

    let worker_cores: Vec<_> = cores.collect();
    let (tx, rx) = mpsc::channel();
    let mut job_txs = Vec::with_capacity(jobs.len());
    job_txs.push(tx);
    for _ in &jobs[1..] {
        job_txs.push(job_txs[0].clone());
    }

    let results_single_fname = bench_args.results_prefix
        .with_extension(".single_seq.tsv");
    let results_full_msa_fname = bench_args.results_prefix
        .with_extension(".full_msa.tsv");

    let mut tsv_writer_single = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(bench_args.output_dir.join(results_single_fname))?;

    let mut tsv_writer_msa = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(bench_args.output_dir.join(results_full_msa_fname))?;

    thread::scope(|scope| -> Result<()> {
        let mut job_iter = jobs.into_iter().zip(job_txs.into_iter());

        // Start initial jobs, limited to number of worker cores
        for core in &worker_cores {
            if let Some(((algorithm, benchmark, dataset), job_tx)) = job_iter.next() {
                run_job(
                    scope, &proc_exe, &bench_args.datasets_dir, &bench_args.output_dir,
                    *algorithm, *benchmark, dataset, Some(*core), job_tx
                );
            } else {
                break;
            }
        }

        // Receive results, and start new jobs when another finishes
        for result in rx {
            eprintln!("Got result: {:?}", result);

            match result {
                JobResult::SingleSeqMeasurement(
                    algo,
                    dataset,
                    score,
                    graph_edges,
                    seq_name,
                    seq_length,
                    measured
                ) => {
                    tsv_writer_single.write_record(&[
                        &dataset,
                        algo.to_str(),
                        &graph_edges.to_string(),
                        &seq_name,
                        &seq_length.to_string(),
                        &score.to_string(),
                        &measured.runtime.to_string(),
                        &measured.memory_initial.map_or(String::default(), |v| v.to_string()),
                        &measured.memory_total.map_or(String::default(), |v| v.to_string()),
                        &measured.memory.to_string(),
                        &measured.time_start.to_string(),
                        &measured.time_end.to_string()
                    ])?;
                },
                JobResult::FullMSAMeasurement(algo, dataset, measured) => {
                    tsv_writer_msa.write_record(&[
                        &dataset,
                        algo.to_str(),
                        &measured.runtime.to_string(),
                        &measured.memory_initial.map_or(String::default(), |v| v.to_string()),
                        &measured.memory_total.map_or(String::default(), |v| v.to_string()),
                        &measured.memory.to_string(),
                        &measured.time_start.to_string(),
                        &measured.time_end.to_string()
                    ])?;
                },
                JobResult::Finished(core) => {
                    tsv_writer_single.flush()?;
                    tsv_writer_msa.flush()?;

                    if let Some(((algorithm, benchmark, dataset), job_tx)) = job_iter.next() {
                        run_job(
                            scope, &proc_exe, &bench_args.datasets_dir, &bench_args.output_dir,
                            *algorithm, *benchmark, dataset, core.map(|v| CoreId { id: v }),
                            job_tx
                        );
                    }
                },
                JobResult::Error => {
                    return Err(POABenchError::WorkerError.into());
                }
            }
        }

        Ok(())
    })?;

    tsv_writer_single.flush()?;

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = CliArgs::parse();

    match args.command {
        Command::Bench(bench_args) => bench(bench_args)?,
        Command::Worker(worker_args) => worker::main(worker_args)?,
    }

    Ok(())
}