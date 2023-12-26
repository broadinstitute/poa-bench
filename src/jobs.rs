use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use crate::bench::Measured;

#[derive(Copy, Clone, Debug, ValueEnum, Serialize, Deserialize)]
pub enum Algorithm {
    POASTA,
    SPOA
}

impl Algorithm {
    pub fn to_str(&self) -> &str {
        match self {
            Self::POASTA => "poasta",
            Self::SPOA => "spoa"
        }
    }
}

pub static ALL_ALGORITHMS: &[Algorithm] = &[Algorithm::POASTA, Algorithm::SPOA];


#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum, Serialize, Deserialize)]
pub enum BenchmarkType {
    SingleSequence,
    FullMSA,
}

impl BenchmarkType {
    pub fn to_str(&self) -> &str {
        match self {
            Self::SingleSequence => "single-sequence",
            Self::FullMSA => "full-msa",
        }
    }
}

pub static ALL_BENCHMARK_TYPES: &[BenchmarkType] = &[
    BenchmarkType::SingleSequence,
    BenchmarkType::FullMSA
];


#[derive(Debug, Serialize, Deserialize)]
pub struct Job {
    pub algorithm: Algorithm,
    pub dataset: String,
    pub benchmark_type: BenchmarkType
}

/// Used for worker-orchestrator IPC
#[derive(Debug, Serialize, Deserialize)]
pub enum JobResult {
    /// Variant to indicate new measurement results from single sequence alignment
    SingleSeqMeasurement(Algorithm, String, usize, usize, usize, String, usize, usize, Measured),

    /// New measurement from the full MSA benchmark
    FullMSAMeasurement(Algorithm, String, Measured),

    /// Variant to indicate the whole dataset has been processed and optionally indicates which
    /// processor core is now free
    Finished(Option<usize>),

    /// Variant to indicate that a worker did not properly exit correctly
    Error,
}
