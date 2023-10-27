use std::{fs, process, thread};
use std::fmt::Debug;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use clap::{Args, command, Parser, Subcommand};
use anyhow::{Result};
use rayon::prelude::*;
use core_affinity;
use core_affinity::CoreId;

use errors::POABenchError;
use crate::dataset::{Dataset, find_datasets};
use crate::jobs::{Algorithm, JobResult};

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
    /// separate them by spaces
    #[clap(value_enum, short, long, num_args=0..)]
    algorithms: Vec<Algorithm>,

    /// Number of parallel threads to start. This number will include the main orchestrator
    /// process, so the number of actual worker threads will be one less than the number specified.
    #[clap(short = 'j', long, default_value="2")]
    parallel: usize,

    #[clap(short, long, default_value="data/")]
    datasets_dir: PathBuf,

    #[clap(short, long, default_value="output/")]
    output_dir: PathBuf,

    #[clap(short='f', long, default_value="results.tsv")]
    results_fname: PathBuf,
}


fn build_graph_with_poasta(output_dir: &Path, dataset: &Dataset) -> Result<(), POABenchError> {
    let seq_fname = dataset.graph_sequences_fname();
    let output_fname = dataset.graph_output_fname(output_dir);
    eprintln!("{:?}, {:?}", seq_fname, output_fname);

    if output_fname.exists() {
        let seq_meta = fs::metadata(seq_fname)?;
        let graph_meta = fs::metadata(output_fname)?;

        if seq_meta.modified()? <= graph_meta.modified()? {
            eprintln!("Found up-to-date graph for dataset: {}", dataset.name());
            return Ok(());
        }
    }

    eprintln!("Creating graph for dataset: {}", dataset.name());

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

fn build_graph_with_spoa(output_dir: &Path, dataset: &Dataset) -> Result<(), POABenchError> {
    let seq_fname = dataset.graph_sequences_fname();
    let output_fname = dataset.graph_msa_fname(output_dir);
    eprintln!("{:?}, {:?}", seq_fname, output_fname);

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
        .output()?;

    if process.status.success() {
        let mut ofile = File::create(&output_fname)
            .map(BufWriter::new)?;
        ofile.write_all(&process.stdout)?;

        Ok(())
    } else {
        Err(POABenchError::BuildGraphError(String::from_utf8(process.stderr)?))
    }

}


fn build_graphs(output_dir: &Path, datasets: &[Dataset]) -> Result<(), POABenchError> {
    let results: Vec<_> = datasets.par_iter()
        .map(|dataset| build_graph_with_spoa(output_dir, dataset))
        .collect();

    for (_, r) in datasets.iter().zip(results) {
        r?
    }

    Ok(())
}

fn run_job<'scope, 'env>(scope: &'scope thread::Scope<'scope, 'env>, proc_exe: &'scope Path, datasets_dir: &'scope Path, output_dir: &'scope Path,
           algorithm: Algorithm, dataset: &'scope Dataset, core: Option<CoreId>, tx: mpsc::Sender<JobResult>
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


fn bench(bench_args: BenchArgs) -> Result<(), POABenchError> {
    eprintln!("Finding data sets...");
    let datasets = find_datasets(&bench_args.datasets_dir)?;

    eprintln!("Building graphs...");
    build_graphs(&bench_args.output_dir, &datasets)?;

    eprintln!("Building job list...");
    let algorithms: &[Algorithm] = if bench_args.algorithms.len() > 0 { &bench_args.algorithms } else { jobs::ALL_ALGORITHMS };

    let mut jobs = Vec::new();
    for algorithm in algorithms {
        for dataset in &datasets {
            jobs.push((algorithm, dataset))
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

    let mut tsv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(bench_args.output_dir.join(bench_args.results_fname))?;

    thread::scope(|scope| {
        let mut job_iter = jobs.into_iter().zip(job_txs.into_iter());

        // Start initial jobs, limited to number of worker cores
        for core in &worker_cores {
            if let Some(((algorithm, dataset), job_tx)) = job_iter.next() {
                run_job(
                    scope, &proc_exe, &bench_args.datasets_dir, &bench_args.output_dir,
                    *algorithm, dataset, Some(*core), job_tx
                );
            } else {
                break;
            }
        }

        // Receive results, and start new jobs when another finishes
        for result in rx {
            eprintln!("Got result: {:?}", result);

            match result {
                JobResult::Measurement(algo, dataset, graph_edges, seq_name, seq_length, measured) => {
                    tsv_writer.write_record(&[
                        &dataset,
                        algo.to_str(),
                        &graph_edges.to_string(),
                        &seq_name,
                        &seq_length.to_string(),
                        &measured.score.to_string(),
                        &measured.runtime.to_string(),
                        &measured.memory_initial.map_or(String::default(), |v| v.to_string()),
                        &measured.memory_total.map_or(String::default(), |v| v.to_string()),
                        &measured.memory.to_string(),
                        &measured.time_start.to_string(),
                        &measured.time_end.to_string()
                    ])?;
                },
                JobResult::Finished(core) => {
                    tsv_writer.flush()?;

                    if let Some(((algorithm, dataset), job_tx)) = job_iter.next() {
                        run_job(
                            scope, &proc_exe, &bench_args.datasets_dir, &bench_args.output_dir,
                            *algorithm, dataset, core.map(|v| CoreId { id: v }),
                            job_tx
                        );
                    }
                },
                JobResult::Error => {
                    return Err(POABenchError::WorkerError);
                }
            }
        }

        Ok(())
    })?;

    tsv_writer.flush()?;

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