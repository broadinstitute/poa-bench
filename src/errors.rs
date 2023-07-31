use std::sync::mpsc;
use std::error::Error;
use std::fmt::{Display, Formatter, write};
use std::{fmt, io};
use poasta::errors::PoastaError;
use std::string::FromUtf8Error;
use crate::jobs::JobResult;

#[derive(Debug)]
pub enum POABenchError {
    IOError(io::Error),
    SendError(mpsc::SendError<JobResult>),
    Utf8Error(FromUtf8Error),
    BuildGraphError(String),
    PoastaError(PoastaError),
    ParseResultError(serde_json::Error),
    ParseConfigError(toml::de::Error),
    WorkerError,
    MemoryResetError,
}

impl Error for POABenchError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::IOError(source) => Some(source),
            Self::SendError(source) => Some(source),
            Self::Utf8Error(source) => Some(source),
            Self::BuildGraphError(_) => None,
            Self::PoastaError(source) => Some(source),
            Self::ParseResultError(source) => Some(source),
            Self::ParseConfigError(source) => Some(source),
            Self::WorkerError => None,
            Self::MemoryResetError => None,
        }
    }
}

impl Display for POABenchError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::IOError(source) =>
                fmt::Display::fmt(source, f),
            Self::SendError(source) =>
                fmt::Display::fmt(source, f),
            Self::Utf8Error(source) =>
                fmt::Display::fmt(source, f),
            Self::BuildGraphError(stderr) =>
                write!(f, "Could not build graph for data set! Stderr: {}", stderr),
            Self::PoastaError(source) =>
                fmt::Display::fmt(source, f),
            Self::ParseResultError(source) => {
                write!(f, "Could not parse job result: ")?;
                fmt::Display::fmt(source, f)
            }
            Self::ParseConfigError(source) => {
                write!(f, "Could not parse dataset config: ")?;
                fmt::Display::fmt(source, f)
            },
            Self::WorkerError => write!(f, "A worker process did not exit properly!"),
            Self::MemoryResetError => write!(f, "Platform does not support resetting max_rss, memory usage measurements will be incorrect!"),
        }
    }
}

impl From<io::Error> for POABenchError {
    fn from(value: io::Error) -> Self {
        Self::IOError(value)
    }
}

impl From<FromUtf8Error> for POABenchError {
    fn from(value: FromUtf8Error) -> Self {
        Self::Utf8Error(value)
    }
}

impl From<PoastaError> for POABenchError {
    fn from(value: PoastaError) -> Self {
        Self::PoastaError(value)
    }
}

impl From<mpsc::SendError<JobResult>> for POABenchError {
    fn from(value: mpsc::SendError<JobResult>) -> Self {
        Self::SendError(value)
    }
}

impl From<serde_json::Error> for POABenchError {
    fn from(value: serde_json::Error) -> Self {
        Self::ParseResultError(value)
    }
}

impl From<toml::de::Error> for POABenchError {
    fn from(value: toml::de::Error) -> Self {
        Self::ParseConfigError(value)
    }
}
