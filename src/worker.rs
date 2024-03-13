use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use anyhow::Result;
use clap::Args;
use core_affinity::CoreId;
use flate2::read::GzDecoder;
use poasta::aligner::config::AffineMinGapCost;
use poasta::graphs::AlignableRefGraph;
use poasta::graphs::poa::{POAGraph, POAGraphWithIx};
use poasta::aligner::PoastaAligner;
use poasta::aligner::scoring::{AlignmentType, GapAffine};
use noodles::fasta;
use poasta::bubbles::index::BubbleIndex;
use poasta::io::graph::save_graph;

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
            seq.name().to_string(),
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

fn bench_single_seq_alignments(dataset: &Dataset, output_dir: &Path, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let graph = poasta::io::load_graph_from_fasta_msa(dataset.graph_msa_fname(output_dir))?;

    match graph {
        POAGraphWithIx::U8(ref g) =>
            perform_alignments_poasta(dataset, g, sequences),
        POAGraphWithIx::U16(ref g) =>
            perform_alignments_poasta(dataset, g, sequences),
        POAGraphWithIx::U32(ref g) =>
            perform_alignments_poasta(dataset, g, sequences),
        POAGraphWithIx::USIZE(ref g) =>
            perform_alignments_poasta(dataset, g, sequences),

    }
}


fn perform_alignments_poasta<G: AlignableRefGraph>(dataset: &Dataset, graph: &G, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let bubbles = Rc::new(BubbleIndex::new(graph));

    let scoring = GapAffine::new(4, 2, 6);
    let aligner = PoastaAligner::new(AffineMinGapCost(scoring), AlignmentType::Global);

    let memory_start = bench::get_maxrss();
    let graph_node_count = graph.node_count();
    let graph_edge_count = graph.edge_count();

    for seq in sequences {
        let bubbles_for_aln = bubbles.clone();
        let (measured, (score, _, num_visited)) = bench::measure(memory_start, || {
            let result = aligner
                .align_with_existing_bubbles::<u32, _, _>(graph, seq.sequence(), bubbles_for_aln);

            (u32::from(result.score) as usize, result.alignment, result.num_visited)
        })?;

        let result = JobResult::SingleSeqMeasurement(
            Algorithm::POASTA,
            dataset.name().to_string(),
            score,
            graph_node_count,
            graph_edge_count,
            seq.name().to_string(),
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
    let scoring = GapAffine::new(4, 2, 6);
    let aligner = PoastaAligner::new(AffineMinGapCost(scoring), AlignmentType::Global);

    let memory_start = bench::get_maxrss();
    let mut graph: POAGraph<u32> = POAGraph::new();

    let (measured, _) = bench::measure(memory_start, || -> Result<(), POABenchError> {
        for seq in sequences {
            let weights: Vec<usize> = vec![1; seq.sequence().len()];

            if graph.is_empty() {
                graph.add_alignment_with_weights(seq.name(), seq.sequence(), None, &weights)?;
            } else {
                let result = aligner
                    .align::<u32, _, _>(&graph, seq.sequence());

                graph.add_alignment_with_weights(seq.name(), seq.sequence(), Some(&result.alignment), &weights)?;
                eprintln!("Aligned {} with score {}", seq.name(), result.score);
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
    let graph_with_ix = POAGraphWithIx::U32(graph);
    let mut graph_outfile = File::create(dataset.poasta_msa_output(output_dir))?;
    save_graph(&graph_with_ix, &mut graph_outfile)?;

    Ok(())
}

fn bench_full_msa_spoa(dataset: &Dataset, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let memory_start = bench::get_maxrss();

    let mut graph = spoa_rs::Graph::new();
    let mut engine = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kNW, 0, -4, -8, -2);

    let (measured, _) = bench::measure(memory_start, || -> Result<(), POABenchError> {
        for record in sequences {
            let seq = std::str::from_utf8(record.sequence().as_ref())?;
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
    }

    let finished = JobResult::Finished(worker_args.core_id);
    println!("{}", serde_json::to_string(&finished)?);

    Ok(())
}
