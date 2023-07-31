use std::path::{Path, PathBuf};
use crate::Dataset;
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use noodles::fasta;
use crate::bench::Measured;

#[derive(Copy, Clone, Debug, ValueEnum, Serialize, Deserialize)]
pub enum Algorithm {
    POASTA,
    SPOA
}

impl Algorithm {
    pub fn to_str(&self) -> &str {
        match self {
            Self::POASTA => &"poasta",
            Self::SPOA => &"spoa"
        }
    }
}

pub static ALL_ALGORITHMS: &'static [Algorithm] = &[Algorithm::POASTA, Algorithm::SPOA];

#[derive(Debug, Serialize, Deserialize)]
pub struct Job {
    pub algorithm: Algorithm,
    pub dataset: String,
}

/// Used for worker-orchestrator IPC
#[derive(Debug, Serialize, Deserialize)]
pub enum JobResult {
    /// Variant to indicate new measurement results
    Measurement(Algorithm, String, Measured),

    /// Variant to indicate the whole dataset has been processed and optionally indicates which
    /// processor core is now free
    Finished(Option<usize>),

    /// Variant to indicate that a worker did not properly exit correctly
    Error,
}