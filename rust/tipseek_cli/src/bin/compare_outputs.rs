//! Byte-level output comparator used while Python and Rust workflows coexist.
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn collect(root: &Path, dir: &Path, paths: &mut Vec<PathBuf>) -> Result<(), String> {
    for entry in fs::read_dir(dir).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let path = entry.path();
        if path.is_dir() {
            collect(root, &path, paths)?;
        } else {
            paths.push(
                path.strip_prefix(root)
                    .map_err(|e| e.to_string())?
                    .to_path_buf(),
            );
        }
    }
    Ok(())
}

fn ignored(path: &Path) -> bool {
    matches!(
        path.file_name().and_then(|name| name.to_str()),
        Some("log.txt") | Some("workflow_profile.tsv") | Some("assembly_profile.tsv")
    )
}

fn main() -> Result<(), String> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.len() != 2 {
        return Err("usage: compare_outputs LEGACY_OUTPUT RUST_OUTPUT".into());
    }
    let (left, right) = (Path::new(&args[0]), Path::new(&args[1]));
    if !left.is_dir() || !right.is_dir() {
        return Err("both arguments must be output directories".into());
    }
    let (mut files, mut other) = (Vec::new(), Vec::new());
    collect(left, left, &mut files)?;
    collect(right, right, &mut other)?;
    files.sort();
    other.sort();
    if files != other {
        return Err("output file sets differ".into());
    }
    for relative in files {
        if !ignored(&relative)
            && fs::read(left.join(&relative)).map_err(|e| e.to_string())?
                != fs::read(right.join(&relative)).map_err(|e| e.to_string())?
        {
            return Err(format!("output differs: {}", relative.display()));
        }
    }
    println!("Compatibility check passed: outputs are identical (runtime logs/profiles ignored).");
    Ok(())
}
