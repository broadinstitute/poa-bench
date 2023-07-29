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

pub static ALL_ALGORITHMS: &'static [Algorithm] = &[Algorithm::POASTA, Algorithm::SPOA];

#[derive(Debug)]
pub struct Job<'a> {
    pub algorithm: Algorithm,
    pub dataset: &'a Dataset,
    pub output_dir: &'a Path,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct JobResult {
    pub algorithm: Algorithm,
    pub dataset: String,
    pub measured: Measured,
}