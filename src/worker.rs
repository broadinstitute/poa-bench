use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::Args;
use core_affinity::CoreId;
use flate2::read::GzDecoder;
use noodles::fasta;

use poasta::aligner::astar::heuristic::Dijkstra;
use poasta::aligner::astar::AlignableGraph;
use poasta::aligner::cost_models::affine::Affine;
use poasta::aligner::{AlignmentMode, GraphAligner, PoastaAligner};
use poasta::graph::io::dot::graph_to_dot;
use poasta::graph::poa::POASeqGraph;

use crate::bench;
use crate::errors::POABenchError;
use crate::dataset::{Dataset, load_dataset, load_align_set_sequences, load_combined_sorted_sequences};
use crate::jobs::{Algorithm, BenchmarkType, JobResult};



#[derive(Args, Debug, Clone)]
pub struct WorkerArgs {
    #[clap(short, long)]
    datasets_dir: PathBuf,

    #[clap(short, long)]
    output_dir: PathBuf,

    #[clap(short, long)]
    core_id: Option<usize>,

    dataset: String,
    algorithm: Algorithm,
    benchmark_type: BenchmarkType,
}

fn make_graph_spoa(dataset: &Dataset) -> Result<spoa_rs::Graph, POABenchError> {
    eprintln!("Preparing SPOA graph for {:?}...", dataset.name());
    let Some(graph_seq_fname) = dataset.graph_sequences_fname() else {
        return Err(POABenchError::BuildGraphError(String::from("No graph sequence set filename")))
    };

    let mut reader = File::open(graph_seq_fname)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    let mut graph = spoa_rs::Graph::new();
    let mut engine = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kNW, 0, -4, -8, -2);

    for record in reader.records() {
        let r = record?;

        let seq = std::str::from_utf8(r.sequence().as_ref())?;
        let (_, aln) = engine.align(seq, &graph);

        graph.add_alignment(aln, seq);
    }

    drop(engine);

    Ok(graph)
}


fn make_graph_abpoa(dataset: &Dataset, output_dir: &Path) -> Result<abpoa_rs::Graph, POABenchError> {
    eprintln!("Loading abPOA graph for {:?}...", dataset.name());
    let graph_seq_fname = dataset.graph_msa_fname(output_dir);

    Ok(abpoa_rs::Graph::from_file(&graph_seq_fname, false))
}

fn perform_alignments_spoa(dataset: &Dataset, graph: &spoa_rs::Graph, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    eprintln!("Performing alignments with SPOA for {:?}...", dataset.name());
    bench::reset_max_rss()?;
    let memory_start = bench::get_maxrss();
    let graph_node_count = graph.node_count();
    let graph_edge_count = graph.edge_count();

    let mut aligner = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kNW, 0, -4, -8, -2);

    for seq in sequences {
        let sequence = std::str::from_utf8(seq.sequence().as_ref())?;
        let (measured, (score, _)) = bench::measure(memory_start, || {
            let (score, alignment) = aligner.align(sequence, graph);

            ((-score) as usize, alignment)
        })?;

        let num_visited = (sequence.len() + 1) * (graph_node_count + 1) * 3;

        let result = JobResult::SingleSeqMeasurement(
            Algorithm::SPOA,
            dataset.name().to_string(),
            score,
            graph_node_count,
            graph_edge_count,
            unsafe { String::from_utf8_unchecked(seq.name().to_vec()) },
            seq.sequence().len(),
            num_visited,
            measured
        );
        let json = match serde_json::to_string(&result) {
            Ok(v) => v,
            Err(e) => return Err(POABenchError::JSONError(e))
        };

        println!("{}", json)
    }

    Ok(())
}

fn perform_alignments_abpoa(dataset: &Dataset, graph: &mut abpoa_rs::Graph, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    eprintln!("Performing alignments with abPOA for {:?}...", dataset.name());
    bench::reset_max_rss()?;

    let memory_start = bench::get_maxrss();
    let graph_node_count = graph.num_nodes() - 2;
    let graph_edge_count = 0;

    let aln_params = abpoa_rs::AlignmentParametersBuilder::new()
        .alignment_mode(abpoa_rs::AlignmentMode::Global)
        .gap_affine_penalties(0, 4, 6, 2)
        .verbosity(abpoa_rs::Verbosity::None)
        .build();

    for seq in sequences {
        let (measured, (score, _)) = bench::measure(memory_start, || {
            let result = graph.align_sequence(
                &aln_params,
                seq.sequence().as_ref(),
            ).unwrap();

            (-result.get_best_score() as usize, result)
        })?;

        let result = JobResult::SingleSeqMeasurement(
            Algorithm::abPOA,
            dataset.name().to_string(),
            score,
            graph_node_count,
            graph_edge_count,
            unsafe { String::from_utf8_unchecked(seq.name().to_vec()) },
            seq.sequence().len(),
            0,
            measured
        );

        let json = match serde_json::to_string(&result) {
            Ok(v) => v,
            Err(e) => return Err(POABenchError::JSONError(e))
        };

        println!("{}", json)
    }

    Ok(())
}

fn bench_single_seq_alignments(dataset: &Dataset, output_dir: &Path, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let mut msa_file = File::open(dataset.graph_msa_fname(output_dir))
        .map(BufReader::new)?;

    let graph = POASeqGraph::<u32>::try_from_fasta_msa(&mut msa_file)?;

    perform_alignments_poasta(dataset, &graph, sequences)
}


fn perform_alignments_poasta<G: AlignableGraph>(dataset: &Dataset, graph: &G, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let cost_model = Affine::new(4, 6, 2);
    let aligner = PoastaAligner::<Dijkstra, i32, u32, _, _>::new(cost_model);

    let memory_start = bench::get_maxrss();
    let graph_node_count = graph.node_count();

    for seq in sequences {
        let (measured, (score, num_visited)) = bench::measure(memory_start, || {
            let result = aligner
                .align(graph, seq.sequence().as_ref(), AlignmentMode::Global).unwrap();

            (result.score.as_usize(), result.num_visited)
        })?;

        let result = JobResult::SingleSeqMeasurement(
            Algorithm::POASTA,
            dataset.name().to_string(),
            score,
            graph_node_count,
            0,
            unsafe { String::from_utf8_unchecked(seq.name().to_vec()) },
            seq.sequence().len(),
            num_visited,
            measured
        );

        let json = match serde_json::to_string(&result) {
            Ok(v) => v,
            Err(e) => return Err(POABenchError::JSONError(e))
        };

        println!("{}", json)
    }

    Ok(())
}

fn bench_full_msa_poasta(
    dataset: &Dataset,
    output_dir: &Path,
    sequences: &[fasta::Record]
) -> Result<(), POABenchError> {
    let aligner = PoastaAligner::<Dijkstra, i32, u32, _, _>::new(Affine::new(4, 6, 2));

    let memory_start = bench::get_maxrss();
    let mut graph = POASeqGraph::<u32>::new();

    let (measured, _) = bench::measure(memory_start, || -> Result<(), POABenchError> {
        for (i, seq) in sequences.iter().enumerate() {
            let weights: Vec<usize> = vec![1; seq.sequence().len()];

            let seq_name = unsafe { std::str::from_utf8_unchecked(seq.name()) };
            if graph.is_empty() {
                graph.add_aligned_sequence(seq_name, seq.sequence(), &weights, None)?;
            } else {
                let result = aligner
                    .align(&graph, seq.sequence(), AlignmentMode::Global)?;

                graph.add_aligned_sequence(seq_name, seq.sequence(), &weights, Some(&result.alignment))?;
                eprintln!("Aligned #{} {} with score {}", i+1, seq_name, result.score.as_usize());
            }
        }

        Ok(())
    })?;

    let result = JobResult::FullMSAMeasurement(
        Algorithm::POASTA,
        dataset.name().to_string(),
        measured
    );

    let json = match serde_json::to_string(&result) {
        Ok(v) => v,
        Err(e) => return Err(POABenchError::JSONError(e))
    };

    println!("{}", json);

    // Save graph to file
    let mut graph_outfile = File::create(dataset.poasta_msa_output(output_dir))?;
    graph_to_dot(&mut graph_outfile, &graph)?;

    Ok(())
}

fn bench_full_msa_spoa(dataset: &Dataset, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let memory_start = bench::get_maxrss();

    let mut graph = spoa_rs::Graph::new();
    let mut engine = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kNW, 0, -4, -8, -2);

    let (measured, _) = bench::measure(memory_start, || -> Result<(), POABenchError> {
        for record in sequences {
            let seq = unsafe { std::str::from_utf8_unchecked(record.sequence().as_ref()) };
            let (_, aln) = engine.align(seq, &graph);

            graph.add_alignment(aln, seq);
        }

        Ok(())
    })?;

    let result = JobResult::FullMSAMeasurement(
        Algorithm::SPOA,
        dataset.name().to_string(),
        measured
    );

    let json = match serde_json::to_string(&result) {
        Ok(v) => v,
        Err(e) => return Err(POABenchError::JSONError(e))
    };

    println!("{}", json);

    Ok(())
}

fn bench_full_msa_abpoa(dataset: &Dataset, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let seqs: Vec<_> = sequences.iter()
        .map(|s| s.sequence().as_ref())
        .collect();
    let weights: Vec<_> = sequences.iter()
        .map(|s| vec![1; s.sequence().len()])
        .collect();
    let names: Vec<_> = sequences.iter()
        .map(|s| s.name())
        .collect();

    let memory_start = bench::get_maxrss();

    let aln_params = abpoa_rs::AlignmentParametersBuilder::new()
        .alignment_mode(abpoa_rs::AlignmentMode::Global)
        .gap_affine_penalties(0, 4, 6, 2)
        .verbosity(abpoa_rs::Verbosity::None)
        .build();

    let mut graph = abpoa_rs::Graph::new(&aln_params);

    let (measured, _) = bench::measure(memory_start, || -> Result<(), POABenchError> {
        let _ = graph.align_and_add_multiple(&aln_params, &seqs, &weights, &names).unwrap();

        Ok(())
    })?;

    let result = JobResult::FullMSAMeasurement(
        Algorithm::abPOA,
        dataset.name().to_string(),
        measured
    );

    let json = match serde_json::to_string(&result) {
        Ok(v) => v,
        Err(e) => return Err(POABenchError::JSONError(e))
    };

    println!("{}", json);

    Ok(())
}

pub fn main(worker_args: WorkerArgs) -> Result<(), POABenchError> {
    if let Some(core_id) = worker_args.core_id {
        core_affinity::set_for_current(CoreId { id: core_id });
    }

    let dataset = load_dataset(&worker_args.datasets_dir, &worker_args.dataset)?;
    let sequences = if worker_args.benchmark_type == BenchmarkType::SingleSequence {
        load_align_set_sequences(&dataset)?
    } else {
        load_combined_sorted_sequences(&dataset, &worker_args.output_dir)?
    };

    match (worker_args.algorithm, worker_args.benchmark_type) {
        (Algorithm::POASTA, BenchmarkType::SingleSequence) =>
            bench_single_seq_alignments(&dataset, &worker_args.output_dir, &sequences)?,
        (Algorithm::POASTA, BenchmarkType::FullMSA) =>
            bench_full_msa_poasta(&dataset, &worker_args.output_dir, &sequences)?,
        (Algorithm::SPOA, BenchmarkType::SingleSequence) => {
            let graph = make_graph_spoa(&dataset)?;
            std::thread::sleep(std::time::Duration::from_secs(1));

            perform_alignments_spoa(&dataset, &graph, &sequences)?;
        },
        (Algorithm::SPOA, BenchmarkType::FullMSA) =>
            bench_full_msa_spoa(&dataset, &sequences)?,
        (Algorithm::abPOA, BenchmarkType::SingleSequence) => {
            let mut graph = make_graph_abpoa(&dataset, &worker_args.output_dir)?;
            std::thread::sleep(std::time::Duration::from_secs(1));

            perform_alignments_abpoa(&dataset, &mut graph, &sequences)?;
        },
        (Algorithm::abPOA, BenchmarkType::FullMSA) =>
            bench_full_msa_abpoa(&dataset, &sequences)?,
    }

    let finished = JobResult::Finished(worker_args.core_id);
    println!("{}", serde_json::to_string(&finished)?);

    Ok(())
}
