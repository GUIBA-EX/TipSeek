#![no_main]

use gm2_tools::fastx::{FastxFormat, FastxReader};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // FastxReader is path-based; use a single process-local temporary file and
    // accept parser errors as normal outcomes. The CI workflow caps executions
    // and duration, and this target is never run implicitly by normal tests.
    let path = std::env::temp_dir().join(format!("tipseek-fastx-fuzz-{}", std::process::id()));
    if std::fs::write(&path, data).is_err() {
        return;
    }
    for format in [FastxFormat::Fastq, FastxFormat::Fasta] {
        if let Ok(mut reader) = FastxReader::open(&path, format) {
            while let Ok(Some(_)) = reader.next_record() {}
        }
    }
    let _ = std::fs::remove_file(path);
});
