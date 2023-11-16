use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::Args;
use core_affinity::CoreId;
use flate2::read::GzDecoder;
use poasta::aligner::config::AffineMinGapCost;
use poasta::graphs::AlignableRefGraph;
use poasta::graphs::poa::POAGraphWithIx;
use poasta::aligner::PoastaAligner;
use poasta::aligner::scoring::{AlignmentType, GapAffine};
use noodles::fasta;
use poasta::bubbles::index::BubbleIndexBuilder;

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
    let graph_edge_count = graph.edge_count();

    let mut aligner = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kNW, 0, -4, -8, -2);

    for seq in sequences {
        let sequence = std::str::from_utf8(seq.sequence().as_ref())?;
        let (measured, alignment) = bench::measure(memory_start, || {
            let (score, alignment) = aligner.align(sequence, graph);

            ((-score) as usize, alignment)
        })?;

        let result = JobResult::Measurement(
            Algorithm::SPOA,
            dataset.name().to_string(),
            graph_edge_count,
            seq.name().to_string(),
            seq.sequence().len(),
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


fn perform_alignments_poasta<G: AlignableRefGraph>(dataset: &Dataset, graph: &G, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let scoring = GapAffine::new(4, 2, 6);
    let mut aligner = PoastaAligner::new(AffineMinGapCost(scoring), AlignmentType::Global);

    let memory_start = bench::get_maxrss();
    let graph_edge_count = graph.edge_count();

    for seq in sequences {
        let (measured, alignment) = bench::measure(memory_start, || {
            let (score, alignment) = aligner.align::<u32, _, _>(graph, seq.sequence());

            (score.into(), alignment)
        })?;

        let result = JobResult::Measurement(
            Algorithm::POASTA,
            dataset.name().to_string(),
            graph_edge_count,
            seq.name().to_string(),
            seq.sequence().len(),
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
            std::thread::sleep(std::time::Duration::from_secs(1));

            perform_alignments_spoa(&dataset, &graph, &sequences)?;
        }
    }

    let finished = JobResult::Finished(worker_args.core_id);
    println!("{}", serde_json::to_string(&finished)?);

    Ok(())
}