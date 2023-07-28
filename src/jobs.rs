use std::path::PathBuf;
use crate::Dataset;
use clap::ValueEnum;

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum Algorithm {
    POASTA,
    SPOA
}

pub static ALL_ALGORITHMS: &'static [Algorithm] = &[Algorithm::POASTA, Algorithm::SPOA];

#[derive(Debug)]
pub struct Job<'a> {
    pub algorithm: Algorithm,
    pub dataset: &'a Dataset,
}