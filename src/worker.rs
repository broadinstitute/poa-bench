use std::path::PathBuf;

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
use crate::dataset::load_dataset;
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


fn perform_alignments_poasta<G: AlignableGraph>(dataset: &str, graph: &G, sequences: &[fasta::Record]) -> Result<(), POABenchError> {
    let scoring = GapAffine::new(4, 2, 6);
    let mut aligner: PoastaAligner<GapAffine> = PoastaAligner::new(scoring);

    let memory_start = bench::get_maxrss();
    let graph_edge_count = graph.edge_count();

    for seq in sequences {
        let measured = bench::measure(memory_start, || {
            let (score, _) = aligner.align::<u32, usize, _, _, _>(graph, seq.sequence());

            score
        })?;

        let result = JobResult::Measurement(Algorithm::POASTA, dataset.to_string(), graph_edge_count, seq.sequence().len(), measured);
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

    match worker_args.algorithm {
        Algorithm::POASTA => {
            let loaded_dataset = dataset.load(&worker_args.output_dir)?;

            match loaded_dataset.graph {
                POAGraphWithIx::U8(ref g) =>
                    perform_alignments_poasta(&worker_args.dataset, g, &loaded_dataset.sequences)?,
                POAGraphWithIx::U16(ref g) =>
                    perform_alignments_poasta(&worker_args.dataset, g, &loaded_dataset.sequences)?,
                POAGraphWithIx::U32(ref g) =>
                    perform_alignments_poasta(&worker_args.dataset, g, &loaded_dataset.sequences)?,
                POAGraphWithIx::USIZE(ref g) =>
                    perform_alignments_poasta(&worker_args.dataset, g, &loaded_dataset.sequences)?,

            }
        },
        Algorithm::SPOA => {
            eprintln!("SPOA not implemented yet")
        }
    }

    let finished = JobResult::Finished(worker_args.core_id);
    println!("{}", serde_json::to_string(&finished)?);

    Ok(())
}