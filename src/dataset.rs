use std::io::BufReader;
use std::fs;
use std::path::{Path, PathBuf};
use toml::Table;
use serde::Deserialize;
use crate::errors::POABenchError;
use anyhow::Context;

use poasta::graphs::poa::POAGraphWithIx;
use noodles::fasta;
use flate2::read::GzDecoder;
use walkdir::WalkDir;

/// A type that combines the data set name (inferred from its path) and the configuration (read
/// from the TOML file)
#[derive(Clone, Debug)]
pub struct Dataset(String, PathBuf, DatasetConfig);

impl Dataset {
    pub fn name(&self) -> &str {
        &self.0
    }

    pub fn cfg(&self) -> &DatasetConfig {
        &self.2
    }

    pub fn output_dir(&self, base_output_path: &Path) -> PathBuf {
        base_output_path.join(&self.0)
    }

    pub fn graph_output_fname(&self, base_dir: &Path) -> PathBuf {
        self.output_dir(base_dir).join("graph.poasta")
    }

    pub fn graph_msa_fname(&self, base_dir: &Path) -> PathBuf {
        self.output_dir(base_dir).join("graph.msa.fasta")
    }

    pub fn graph_sequences_fname(&self) -> PathBuf {
        self.1.join(&self.2.graph_set.fname)
    }

    pub fn align_sequences_fname(&self) -> PathBuf {
        self.1.join(&self.2.align_set.fname)
    }

    pub fn load(&self, base_dir: &Path) -> Result<LoadedDataset, POABenchError> {
        // Load the graph represented by the MSA created by SPOA, to keep the graphs
        // consistent across tools
        let graph = poasta::io::load_graph_from_fasta_msa(self.graph_msa_fname(base_dir))?;

        let mut seq_file = fs::File::open(self.align_sequences_fname())
            .map(GzDecoder::new)
            .map(BufReader::new)
            .map(fasta::Reader::new)?;

        let sequences: Vec<_> = seq_file.records().filter_map(Result::ok).collect();

        Ok(LoadedDataset { graph, sequences })

    }
}

pub fn find_datasets(datasets_dir: &Path) -> Result<Vec<Dataset>, POABenchError> {
    let mut datasets = Vec::new();
    for entry in WalkDir::new(datasets_dir)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| !e.file_type().is_dir())
    {
        let fname = entry.file_name().to_string_lossy();

        if fname.ends_with(".toml") {
            let contents = fs::read_to_string(entry.path())?;
            let dataset_cfg: DatasetConfig = toml::from_str(&contents)?;

            let dataset_dir = entry.path().parent().unwrap();
            let dataset_name = dataset_dir
                .strip_prefix(datasets_dir).unwrap()
                .to_string_lossy().to_string();

            datasets.push(Dataset(dataset_name, dataset_dir.to_owned(), dataset_cfg))
        }
    }

    Ok(datasets)
}


pub fn load_dataset(datasets_dir: &Path, dataset_name: &str) -> Result<Dataset, POABenchError> {
    let dataset_dir = datasets_dir.join(dataset_name);
    let meta = dataset_dir.join("meta.toml");

    let contents = fs::read_to_string(&meta)?;
    let dataset_cfg: DatasetConfig = toml::from_str(&contents)?;

    Ok(Dataset(dataset_name.to_string(), dataset_dir, dataset_cfg))
}

pub fn load_dataset_sequences(dataset: &Dataset) -> Result<Vec<fasta::Record>, POABenchError> {
    let mut seq_file = fs::File::open(dataset.align_sequences_fname())
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    Ok(seq_file.records().filter_map(Result::ok).collect())
}


/// Represents a single benchmark data set configuration
///
/// A data set is represented by a set of sequences used for constructing the graph,
/// and a set of sequences used for aligning to the graph
///
/// These structs optionally hold some other metadata, e.g. num sequences, average
/// pairwise distances, etc.
///
/// They are meant to be populated by parsing a data set metadata TOML file.
#[derive(Debug, Clone, Deserialize)]
pub struct DatasetConfig {
    pub clustering_max_dist: Option<f32>,

    pub graph_set: GraphSet,
    pub align_set: AlignSet,
}

#[derive(Debug, Clone, Deserialize)]
pub struct GraphSet {
    pub fname: PathBuf,
    pub num_seqs: Option<usize>,
    pub avg_seq_len: Option<f32>,
    pub avg_pairwise_dist: Option<f32>,
    pub species: Option<Table>
}

#[derive(Debug, Clone, Deserialize)]
pub struct AlignSet {
    pub fname: PathBuf,
    pub num_seqs: Option<usize>,
    pub avg_seq_len: Option<f32>,
    pub avg_pairwise_dist: Option<f32>,
    pub avg_min_dist_to_graph: Option<f32>,
    pub species: Option<Table>
}
