use std::path::PathBuf;
use toml::Table;
use serde::Deserialize;

/// Represents a single benchmark data set
///
/// A data set is represented by a set of sequences used for constructing the graph,
/// and a set of sequences used for aligning to the graph
///
/// These structs optionally hold some other metadata, e.g. num sequences, average
/// pairwise distances, etc.
///
/// They are meant to be populated by parsing a data set metadata TOML file.
#[derive(Debug, Deserialize)]
struct Dataset {
    clustering_max_dist: Option<f32>,

    graph_set: GraphSet,
    align_set: AlignSet,
}

#[derive(Debug, Deserialize)]
struct GraphSet {
    fname: PathBuf,
    num_seqs: Option<usize>,
    avg_seq_len: Option<f32>,
    avg_pairwise_dist: Option<f32>,
    species: Option<Table>
}

#[derive(Debug, Deserialize)]
struct AlignSet {
    fname: PathBuf,
    num_seqs: Option<usize>,
    avg_seq_len: Option<f32>,
    avg_pairwise_dist: Option<f32>,
    species: Option<Table>
}
