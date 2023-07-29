use std::path::PathBuf;
use toml::Table;
use serde::{Deserialize, Serialize};

/// Represents a single benchmark data set
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

