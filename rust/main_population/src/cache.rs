//! Dependency fingerprints for explicit stage restart. Completion receipts are
//! removed before mutation and installed atomically only after successful output.
use super::*;
use sha2::{Digest, Sha256};
use std::io::Read;

pub(crate) fn digest(path: &Path) -> AppResult<String> {
    let mut file = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut hash = Sha256::new();
    let mut buffer = [0_u8; 65536];
    loop {
        let n = file.read(&mut buffer).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }
        hash.update(&buffer[..n]);
    }
    Ok(format!("{:x}", hash.finalize()))
}

// Reads and large alignment files use filesystem identity, avoiding another
// full raw-data scan. Unix ctime also detects edits with restored mtime.
pub(crate) fn identity(path: &Path) -> AppResult<String> {
    let path = path
        .canonicalize()
        .map_err(|e| format!("{}: {e}", path.display()))?;
    let meta = fs::metadata(&path).map_err(|e| e.to_string())?;
    if !meta.is_file() {
        return Err(format!("not a file: {}", path.display()));
    }
    let mut value = format!(
        "{:?}:{}:{:?}",
        path,
        meta.len(),
        meta.modified().map_err(|e| e.to_string())?
    );
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        value.push_str(&format!(
            ":{}:{}:{}:{}",
            meta.dev(),
            meta.ino(),
            meta.ctime(),
            meta.ctime_nsec()
        ));
    }
    Ok(value)
}

pub(crate) fn key(fields: &[String]) -> String {
    let mut hash = Sha256::new();
    hash.update(b"population-cache-v1");
    hash.update(env!("CARGO_PKG_VERSION").as_bytes());
    for field in fields {
        hash.update((field.len() as u64).to_le_bytes());
        hash.update(field.as_bytes());
    }
    format!("{:x}", hash.finalize())
}

fn executable(path: &Path) -> bool {
    let Ok(meta) = fs::metadata(path) else {
        return false;
    };
    if !meta.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        meta.permissions().mode() & 0o111 != 0
    }
    #[cfg(not(unix))]
    {
        true
    }
}

pub(crate) fn tool(program: &str) -> AppResult<String> {
    let path = if Path::new(program).components().count() > 1 {
        PathBuf::from(program)
    } else {
        env::split_paths(&env::var_os("PATH").unwrap_or_default())
            .map(|root| root.join(program))
            .find(|p| executable(p))
            .ok_or_else(|| format!("executable not found: {program}"))?
    };
    if !executable(&path) {
        return Err(format!("not executable: {}", path.display()));
    }
    Ok(format!("{}:{}", identity(&path)?, digest(&path)?))
}

pub(crate) fn samples(samples: &[Sample]) -> AppResult<String> {
    let mut fields = Vec::new();
    for sample in samples {
        fields.push(format!(
            "{:?}|{:?}|{:?}",
            sample.original, sample.internal, sample.vcf_id
        ));
        fields.push(identity(&sample.read1)?);
        fields.push(identity(&sample.read2)?);
    }
    Ok(key(&fields))
}

fn tree(path: &Path, fields: &mut Vec<String>) -> AppResult<()> {
    let mut entries = fs::read_dir(path)
        .map_err(|e| format!("{}: {e}", path.display()))?
        .map(|e| e.map(|e| e.path()))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| e.to_string())?;
    entries.sort();
    for entry in entries {
        // Only direct reference/result files are relevant; never follow directories.
        if entry.is_file() {
            fields.push(format!("{:?}:{}", entry.file_name(), digest(&entry)?));
        }
    }
    Ok(())
}

pub(crate) fn reference_inputs(
    args: &Args,
    samples: &[Sample],
    sample_key: &str,
) -> AppResult<String> {
    let mut fields = vec![
        sample_key.into(),
        format!(
            "{:?}:{:?}:{}",
            args.engine, args.reference_strategy, args.panrefv2_include_low_confidence
        ),
    ];
    if let Some(path) = &args.reference_fasta {
        fields.push(digest(path)?);
    } else if args.engine == Engine::Pseudoref {
        for sample in samples {
            let root = args.output.join(&sample.internal);
            fields.push(digest(&root.join("uce_assembly_summary.csv"))?);
            tree(&root.join("results"), &mut fields)?;
        }
    } else {
        tree(
            args.panref_baits.as_ref().ok_or("missing bait directory")?,
            &mut fields,
        )?;
    }
    Ok(key(&fields))
}

fn receipt(args: &Args, stage: &str) -> PathBuf {
    args.output
        .join("population/cache")
        .join(format!("{stage}.tsv"))
}
fn contents(input: &str, outputs: &[PathBuf]) -> AppResult<String> {
    let mut rows = format!("schema\t1\ninput\t{input}\n");
    for path in outputs {
        rows.push_str(&format!("output\t{}\n", identity(path)?));
    }
    Ok(rows)
}
pub(crate) fn validate(
    args: &Args,
    stage: &str,
    input: &str,
    outputs: &[PathBuf],
) -> AppResult<()> {
    let expected = contents(input, outputs)?;
    let actual = fs::read_to_string(receipt(args, stage))
        .map_err(|_| format!("missing {stage} completion receipt; rerun the {stage} stage"))?;
    if actual != expected {
        return Err(format!(
            "{stage} cache does not match inputs, parameters or outputs; rerun the {stage} stage"
        ));
    }
    Ok(())
}
pub(crate) fn invalidate(args: &Args, stages: &[&str]) -> AppResult<()> {
    for stage in stages {
        match fs::remove_file(receipt(args, stage)) {
            Ok(()) => (),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
            Err(e) => return Err(e.to_string()),
        }
    }
    Ok(())
}
pub(crate) fn save(args: &Args, stage: &str, input: &str, outputs: &[PathBuf]) -> AppResult<()> {
    let value = contents(input, outputs)?;
    let path = receipt(args, stage);
    fs::create_dir_all(path.parent().unwrap()).map_err(|e| e.to_string())?;
    let tmp = path.with_extension("tmp");
    fs::write(&tmp, value).map_err(|e| e.to_string())?;
    fs::rename(tmp, path).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(unix)]
    #[test]
    fn tool_identity_rejects_non_executable_files() {
        use std::os::unix::fs::PermissionsExt;
        let path = env::temp_dir().join(format!("population-tool-{}", process::id()));
        fs::write(&path, "#!/bin/sh\nexit 0\n").unwrap();
        fs::set_permissions(&path, fs::Permissions::from_mode(0o600)).unwrap();
        assert!(tool(path.to_str().unwrap())
            .unwrap_err()
            .contains("not executable"));
        fs::set_permissions(&path, fs::Permissions::from_mode(0o700)).unwrap();
        assert!(tool(path.to_str().unwrap()).is_ok());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn receipts_reject_replaced_outputs_and_interrupted_stages() {
        let root = env::temp_dir().join(format!("population-receipt-{}", process::id()));
        fs::create_dir_all(&root).unwrap();
        let args = Args {
            output: root.clone(),
            ..Args::default()
        };
        let output = root.join("result");
        fs::write(&output, "AAAA").unwrap();
        let outputs = vec![output.clone()];
        assert!(validate(&args, "mapping", "key", &outputs).is_err());
        save(&args, "mapping", "key", &outputs).unwrap();
        validate(&args, "mapping", "key", &outputs).unwrap();
        assert!(validate(&args, "mapping", "different", &outputs).is_err());
        let replacement = root.join("new");
        fs::write(&replacement, "TTTT").unwrap();
        fs::rename(replacement, &output).unwrap();
        assert!(validate(&args, "mapping", "key", &outputs).is_err());
        save(&args, "mapping", "key", &outputs).unwrap();
        invalidate(&args, &["mapping"]).unwrap();
        assert!(validate(&args, "mapping", "key", &outputs).is_err());
        fs::remove_dir_all(root).unwrap();
    }
}
