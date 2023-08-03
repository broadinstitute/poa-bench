use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::Args;
use core_affinity::CoreId;
use poasta::graphs::AlignableGraph;
use poasta::graphs::poa::POAGraphWithIx;
use poasta::aligner::PoastaAligner;
use poasta::aligner::scoring::GapAffine;
use noodles::fasta;
use serde::ser::Error;

use crate::bench;
use crate::errors::POABenchError;
use crate::dataset::{Dataset, load_dataset, load_dataset_sequences};
use crate::jobs::{Algorithm, JobResult};



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
}

fn make_graph_spoa(dataset: &Dataset) -> Result<spoa_rs::Graph, POABenchError> {
    eprintln!("Preparing SPOA graph for {:?}...", dataset.name());
    let graph_seq_fname = dataset.graph_sequences_fname();
    let mut reader = File::open(&graph_seq_fname)
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

    Ok(graph)
}

fn perform_alignments_spoa(dataset: &Dataset, graph: &spoa_rs::Graph, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    eprintln!("Performing alignments with SPOA for {:?}...", dataset.name());
    let mut aligner = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kNW, 0, -4, -8, -2);

    let memory_start = bench::get_maxrss();
    let graph_edge_count = graph.edge_count();

    for seq in sequences {
        let sequence = std::str::from_utf8(seq.sequence().as_ref())?;
        let measured = bench::measure(memory_start, || {
            let (score, _) = aligner.align(sequence, graph);

            (-score) as usize
        })?;

        let result = JobResult::Measurement(Algorithm::SPOA, dataset.name().to_string(), graph_edge_count, seq.sequence().len(), measured);
        let json = match serde_json::to_string(&result) {
            Ok(v) => v,
            Err(e) => return Err(POABenchError::JSONError(e))
        };

        println!("{}", json)
    }

    Ok(())
}

fn load_graph_and_align_poasta(dataset: &Dataset, output_dir: &Path, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let graph = poasta::io::load_graph_from_fasta_msa(dataset.graph_msa_fname(output_dir))?;

    match graph {
        POAGraphWithIx::U8(ref g) =>
            perform_alignments_poasta(&dataset, g, &sequences),
        POAGraphWithIx::U16(ref g) =>
            perform_alignments_poasta(&dataset, g, &sequences),
        POAGraphWithIx::U32(ref g) =>
            perform_alignments_poasta(&dataset, g, &sequences),
        POAGraphWithIx::USIZE(ref g) =>
            perform_alignments_poasta(&dataset, g, &sequences),

    }
}


fn perform_alignments_poasta<G: AlignableGraph>(dataset: &Dataset, graph: &G, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let scoring = GapAffine::new(4, 2, 6);
    let mut aligner: PoastaAligner<GapAffine> = PoastaAligner::new(scoring);

    let memory_start = bench::get_maxrss();
    let graph_edge_count = graph.edge_count();

    for seq in sequences {
        let measured = bench::measure(memory_start, || {
            let (score, _) = aligner.align::<u32, usize, _, _, _>(graph, seq.sequence());

            score
        })?;

        let result = JobResult::Measurement(Algorithm::POASTA, dataset.name().to_string(), graph_edge_count, seq.sequence().len(), measured);
        let json = match serde_json::to_string(&result) {
            Ok(v) => v,
            Err(e) => return Err(POABenchError::JSONError(e))
        };

        println!("{}", json)
    }

    Ok(())
}

pub fn main(worker_args: WorkerArgs) -> Result<(), POABenchError> {
    if let Some(core_id) = worker_args.core_id {
        core_affinity::set_for_current(CoreId { id: core_id });
    }

    let dataset = load_dataset(&worker_args.datasets_dir, &worker_args.dataset)?;
    let sequences = load_dataset_sequences(&dataset)?;

    match worker_args.algorithm {
        Algorithm::POASTA => load_graph_and_align_poasta(&dataset, &worker_args.output_dir, &sequences)?,
        Algorithm::SPOA => {
            let graph = make_graph_spoa(&dataset)?;

            perform_alignments_spoa(&dataset, &graph, &sequences)?;
        }
    }

    let finished = JobResult::Finished(worker_args.core_id);
    println!("{}", serde_json::to_string(&finished)?);

    Ok(())
}