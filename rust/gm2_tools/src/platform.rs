//! Host processor and memory detection shared by TipSeek tools.
//!
//! Detection only informs thread placement, memory budgets and log lines; it
//! never changes what a tool computes, so outputs stay identical across
//! machines.

#[cfg(target_os = "linux")]
use std::fs;
use std::sync::OnceLock;
use std::thread;

/// Processor model, core classes and instruction-set features of the host.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct CpuInfo {
    pub model: String,
    pub logical: usize,
    pub physical: usize,
    /// Performance cores on a hybrid processor (Apple P cores).
    pub performance: Option<usize>,
    /// Efficiency cores on a hybrid processor (Apple E cores).
    pub efficiency: Option<usize>,
    pub cache_line_bytes: Option<usize>,
    pub features: Vec<&'static str>,
}

impl CpuInfo {
    /// One-line description, e.g. `Apple M2, 4 performance + 4 efficiency
    /// cores, aarch64 neon crc dotprod`.
    pub fn summary(&self) -> String {
        let cores = match (self.performance, self.efficiency) {
            (Some(performance), Some(efficiency)) if efficiency > 0 => {
                format!("{performance} performance + {efficiency} efficiency cores")
            }
            _ => format!(
                "{} physical core(s), {} logical CPU(s)",
                self.physical, self.logical
            ),
        };
        let model = if self.model.is_empty() {
            "unknown processor"
        } else {
            self.model.as_str()
        };
        format!(
            "{model}, {cores}, {} {}",
            std::env::consts::ARCH,
            if self.features.is_empty() {
                "baseline".to_string()
            } else {
                self.features.join(" ")
            }
        )
    }
}

static CPU_INFO: OnceLock<CpuInfo> = OnceLock::new();

/// Cached host processor description.
pub fn cpu_info() -> &'static CpuInfo {
    CPU_INFO.get_or_init(detect_cpu)
}

fn detect_cpu() -> CpuInfo {
    let logical = thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1);
    let mut info = CpuInfo {
        logical,
        physical: logical,
        features: simd_features(),
        ..CpuInfo::default()
    };
    #[cfg(target_os = "macos")]
    {
        info.model = sysctl_string("machdep.cpu.brand_string").unwrap_or_default();
        if let Some(physical) = sysctl_u64("hw.physicalcpu") {
            info.physical = physical as usize;
        }
        let performance = sysctl_u64("hw.perflevel0.physicalcpu").map(|count| count as usize);
        let efficiency = sysctl_u64("hw.perflevel1.physicalcpu").map(|count| count as usize);
        if performance.is_some() {
            info.performance = performance;
            info.efficiency = Some(efficiency.unwrap_or(0));
        }
        info.cache_line_bytes = sysctl_u64("hw.cachelinesize").map(|bytes| bytes as usize);
    }
    #[cfg(target_os = "linux")]
    {
        info.model = linux_cpu_model().unwrap_or_default();
        info.cache_line_bytes =
            fs::read_to_string("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size")
                .ok()
                .and_then(|value| value.trim().parse().ok());
    }
    info
}

#[cfg(target_os = "linux")]
fn linux_cpu_model() -> Option<String> {
    let contents = fs::read_to_string("/proc/cpuinfo").ok()?;
    contents.lines().find_map(|line| {
        let (key, value) = line.split_once(':')?;
        matches!(key.trim(), "model name" | "Model" | "cpu model").then(|| value.trim().to_owned())
    })
}

/// Instruction-set extensions relevant to TipSeek's integer kernels.
fn simd_features() -> Vec<&'static str> {
    #[allow(unused_mut)]
    let mut features = Vec::new();
    #[cfg(target_arch = "x86_64")]
    {
        macro_rules! probe {
            ($($name:tt),*) => {
                $(if std::arch::is_x86_feature_detected!($name) {
                    features.push($name);
                })*
            };
        }
        probe!("sse4.2", "popcnt", "bmi2", "lzcnt", "avx2", "avx512f", "avx512bw", "avx512vl");
    }
    #[cfg(target_arch = "aarch64")]
    {
        macro_rules! probe {
            ($($name:tt),*) => {
                $(if std::arch::is_aarch64_feature_detected!($name) {
                    features.push($name);
                })*
            };
        }
        probe!("neon", "crc", "lse", "dotprod", "sha3", "i8mm", "sve", "sve2");
    }
    features
}

/// Memory that can be allocated without swapping, in MiB.
///
/// Linux reports `MemAvailable`. macOS counts free, purgeable and
/// file-backed pages, which the kernel reclaims before compressing or
/// swapping; on Apple silicon this memory is also shared with the GPU, so
/// callers should keep their usual safety margin.
pub fn available_memory_mib() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let contents = fs::read_to_string("/proc/meminfo").ok()?;
        contents.lines().find_map(|line| {
            let mut fields = line.split_whitespace();
            (fields.next()? == "MemAvailable:")
                .then(|| fields.next()?.parse::<u64>().ok())??
                .checked_div(1024)
        })
    }
    #[cfg(target_os = "macos")]
    {
        let page = sysctl_u64("hw.pagesize")?;
        let pages = sysctl_u64("vm.page_free_count")?
            + sysctl_u64("vm.page_purgeable_count").unwrap_or(0)
            + sysctl_u64("vm.page_pageable_external_count").unwrap_or(0);
        Some(pages.saturating_mul(page) / (1024 * 1024))
    }
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        None
    }
}

/// Marks the calling thread as latency-sensitive compute work.
///
/// On Apple silicon this asks the scheduler to prefer performance cores
/// (`QOS_CLASS_USER_INITIATED`); elsewhere it does nothing.
pub fn prefer_performance_cores() {
    #[cfg(target_os = "macos")]
    {
        const QOS_CLASS_USER_INITIATED: u32 = 0x19;
        extern "C" {
            fn pthread_set_qos_class_self_np(
                qos_class: u32,
                relative_priority: std::ffi::c_int,
            ) -> std::ffi::c_int;
        }
        // Failure only leaves the default QoS in place.
        unsafe {
            pthread_set_qos_class_self_np(QOS_CLASS_USER_INITIATED, 0);
        }
    }
}

#[cfg(target_os = "macos")]
fn sysctl_raw(name: &str, buffer: &mut [u8]) -> Option<usize> {
    let name = std::ffi::CString::new(name).ok()?;
    let mut length = buffer.len();
    let status = unsafe {
        libc::sysctlbyname(
            name.as_ptr(),
            buffer.as_mut_ptr().cast(),
            &mut length,
            std::ptr::null_mut(),
            0,
        )
    };
    (status == 0).then_some(length)
}

#[cfg(target_os = "macos")]
fn sysctl_u64(name: &str) -> Option<u64> {
    let mut buffer = [0_u8; 8];
    match sysctl_raw(name, &mut buffer)? {
        4 => Some(u64::from(u32::from_ne_bytes(buffer[..4].try_into().ok()?))),
        8 => Some(u64::from_ne_bytes(buffer)),
        _ => None,
    }
}

#[cfg(target_os = "macos")]
fn sysctl_string(name: &str) -> Option<String> {
    let mut buffer = [0_u8; 256];
    let length = sysctl_raw(name, &mut buffer)?;
    let text = &buffer[..length];
    let end = text
        .iter()
        .position(|&byte| byte == 0)
        .unwrap_or(text.len());
    Some(String::from_utf8_lossy(&text[..end]).trim().to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detection_reports_a_consistent_host() {
        let info = cpu_info();
        assert!(info.logical >= 1);
        assert!(info.physical >= 1 && info.physical <= info.logical);
        if let (Some(performance), Some(efficiency)) = (info.performance, info.efficiency) {
            assert_eq!(performance + efficiency, info.physical);
        }
        assert!(info.summary().contains(std::env::consts::ARCH));
        #[cfg(any(target_os = "linux", target_os = "macos"))]
        assert!(available_memory_mib().is_some_and(|mib| mib > 0));
    }
}
