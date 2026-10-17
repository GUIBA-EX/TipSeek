use super::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct MappingBudget {
    pub workers: usize,
    pub mapper_threads: usize,
    pub sort_extra_threads: usize,
    pub sort_memory_mib: u64,
    pub sequential: bool,
}

impl MappingBudget {
    pub fn plan(threads: usize, samples: usize, available_mib: u64) -> Self {
        let threads = threads.max(1);
        // Leave half of available memory for mapper indexes, runtime overhead
        // and other processes. samtools -m is per sorting thread, not per job.
        let sort_budget = (available_mib / 2).max(1);
        let workers = (threads / 2)
            .max(1)
            .min(samples.max(1))
            .min((sort_budget / 256).max(1) as usize);
        let per_worker = threads / workers;
        let sequential = per_worker == 1;
        let mapper_threads = if sequential {
            1
        } else {
            (per_worker / 2).max(1)
        };
        let sort_threads = if sequential {
            1
        } else {
            per_worker - mapper_threads
        };
        Self {
            workers,
            mapper_threads,
            sort_extra_threads: sort_threads - 1,
            sort_memory_mib: (sort_budget / workers as u64 / sort_threads as u64).clamp(1, 512),
            sequential,
        }
    }
}

pub(crate) fn available_memory_mib() -> u64 {
    let mut available = gm2_tools::platform::available_memory_mib().unwrap_or(1024);
    // Container limits may be much lower than host MemAvailable. Walk the
    // process's cgroup-v2 ancestors, including the mounted root.
    if let Ok(groups) = fs::read_to_string("/proc/self/cgroup") {
        if let Some(relative) = groups.lines().find_map(|l| l.strip_prefix("0::")) {
            let root = Path::new("/sys/fs/cgroup");
            let mut group = root.join(relative.trim_start_matches('/'));
            while group.starts_with(root) {
                if let (Ok(limit), Ok(current)) = (
                    fs::read_to_string(group.join("memory.max")),
                    fs::read_to_string(group.join("memory.current")),
                ) {
                    if let (Ok(limit), Ok(current)) =
                        (limit.trim().parse::<u64>(), current.trim().parse::<u64>())
                    {
                        available = available.min(limit.saturating_sub(current) / (1024 * 1024));
                    }
                }
                if group == root || !group.pop() {
                    break;
                }
            }
        }
    }
    available
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn combined_pipeline_budget_and_sort_memory_are_bounded() {
        for threads in 1..65 {
            for memory in [64, 512, 1024, 8192] {
                let b = MappingBudget::plan(threads, 100, memory);
                let concurrent = if b.sequential {
                    1
                } else {
                    b.mapper_threads + b.sort_extra_threads + 1
                };
                assert!(b.workers * concurrent <= threads);
                assert!(
                    b.workers as u64 * (b.sort_extra_threads + 1) as u64 * b.sort_memory_mib
                        <= memory / 2
                );
            }
        }
        assert!(MappingBudget::plan(1, 10, 1024).sequential);
    }
}
