//! Code to measure runtime and memory for a certain function
//!
//! All code in this module taken from https://github.com/pairwise-alignment/pa-bench
//! and written by Ragnar Groot-Koerkamp and Daniel Liu

use std::{path::Path, process, time::Instant};
use std::fs::File;
use std::io::Write;

use chrono::SubsecRound;
use libc;

use serde::{Serialize, Deserialize};
use crate::errors::POABenchError;

pub type Bytes = u64;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Measured {
    /// Alignment score
    pub score: usize,
    /// Runtime in seconds.
    pub runtime: f32,
    /// max_rss after reading input file.
    pub memory_initial: Option<Bytes>,
    /// max_rss at the end.
    pub memory_total: Option<Bytes>,
    /// Increase in memory usage.
    pub memory: Bytes,
    /// Formatted UTC time when run was started/ended.
    pub time_start: chrono::DateTime<chrono::Utc>,
    pub time_end: chrono::DateTime<chrono::Utc>,
    /// Cpu core running this process at start/end.
    pub cpu_start: Option<i32>,
    pub cpu_end: Option<i32>,
    /// Cpu frequency at start/end.
    pub cpu_freq_start: Option<f32>,
    pub cpu_freq_end: Option<f32>,
}

/// F can return some state that is dropped only after the memory is measured.
pub fn measure<F: FnOnce() -> usize>(memory_start: Bytes, f: F) -> Result<Measured, POABenchError> {
    reset_max_rss()?;

    let cpu_start = get_cpu();
    let cpu_freq_start = cpu_start.and_then(|c| get_cpu_freq(c));
    let time_start = chrono::Utc::now().trunc_subsecs(3);
    let start = Instant::now();

    let score = f();

    let runtime = start.elapsed().as_secs_f32();
    let time_end = chrono::Utc::now().trunc_subsecs(3);
    let memory_total = get_maxrss();
    let memory = memory_total.saturating_sub(memory_start);
    let cpu_end = get_cpu();
    let cpu_freq_end = cpu_end.and_then(|c| get_cpu_freq(c));

    Ok(Measured {
        score,
        runtime,
        memory_initial: Some(memory_start),
        memory_total: Some(memory_total),
        memory,
        time_start,
        time_end,
        cpu_start,
        cpu_end,
        cpu_freq_start,
        cpu_freq_end,
    })
}

/// Returns the maximum resident set size, i.e. the physical memory the thread
/// uses, in bytes.
pub fn get_maxrss() -> Bytes {
    let rusage = unsafe {
        let mut rusage = std::mem::MaybeUninit::uninit();
        libc::getrusage(libc::RUSAGE_SELF, rusage.as_mut_ptr());
        rusage.assume_init()
    };
    let maxrss = rusage.ru_maxrss as _;
    if cfg!(target_os = "macos") {
        // maxrss is in bytes
        maxrss
    } else {
        // maxrss is in kilobytes
        maxrss * 1024
    }
}

fn get_cpu_freq(cur_cpu: i32) -> Option<f32> {
    let path = format!("/sys/devices/system/cpu/cpu{cur_cpu}/cpufreq/scaling_cur_freq");
    let path = &Path::new(&path);
    if !path.exists() {
        return None;
    }

    let val = std::fs::read_to_string(path).ok()?;
    Some(val.trim().parse::<f32>().ok()? / 1000.0)
}

fn get_cpu() -> Option<i32> {
    #[cfg(not(target_os = "macos"))]
    {
        Some(unsafe { libc::sched_getcpu() })
    }
    #[cfg(target_os = "macos")]
    {
        None
    }
}

pub fn reset_max_rss() -> Result<(), POABenchError> {
    #[cfg(target_os = "linux")]
    {
        let pid = process::id();

        // https://stackoverflow.com/questions/45361849/get-memory-high-water-mark-for-a-time-interval
        let mut f = File::create(format!("/proc/{pid}/clear_refs"))?;
        f.write(b"5")?;

        Ok(())
    }
    #[cfg(not(target_os = "linux"))]
    {
        Err(POABenchError::MemoryResetError)
    }
}