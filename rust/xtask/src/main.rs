use sha2::{Digest, Sha256};
use std::env;
use std::fs;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

const BINARIES: &[(&str, &str)] = &[
    ("build_consensus", "build_consensus"),
    ("MainFilterNew", "MainFilterNew"),
    ("main_refilter_new", "main_refilter_new"),
    ("uce_filter", "uce_filter"),
    ("main_assembler_original", "main_assembler-original-rust"),
    ("main_assembler", "main_assembler-rust"),
    ("main_population", "main_population"),
    ("fix_alignment", "fix_alignment"),
    ("merge_seq", "merge_seq"),
    ("build_trimed", "build_trimed"),
    ("gm2_stats", "gm2_stats"),
    ("mito_workflow", "mito_workflow"),
    ("gene_workflow", "gene_workflow"),
    ("rad_workflow", "rad_workflow"),
    ("marker_profile", "marker_profile"),
    ("main_repeat", "main_repeat"),
    ("tipseek", "tipseek-rust"),
];

fn usage() {
    eprintln!("Usage: cargo run -p xtask -- <build|clean>");
}

fn workspace_root() -> io::Result<PathBuf> {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .map(Path::to_path_buf)
        .ok_or_else(|| io::Error::other("xtask is not inside a Cargo workspace"))
}

fn run_cargo(root: &Path, args: &[&str]) -> io::Result<()> {
    let status = Command::new("cargo")
        .args(args)
        .current_dir(root)
        .status()?;
    if status.success() {
        Ok(())
    } else {
        Err(io::Error::other("Cargo build failed"))
    }
}

#[cfg(unix)]
fn make_entrypoint(target: &Path, entrypoint: &Path) -> io::Result<()> {
    use std::os::unix::fs::symlink;
    symlink(target, entrypoint)
}

#[cfg(not(unix))]
fn make_entrypoint(target: &Path, entrypoint: &Path) -> io::Result<()> {
    std::os::windows::fs::symlink_file(target, entrypoint)
}

#[cfg(unix)]
fn mark_executable(path: &Path) -> io::Result<()> {
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(path, fs::Permissions::from_mode(0o755))
}

#[cfg(not(unix))]
fn mark_executable(_: &Path) -> io::Result<()> {
    Ok(())
}

fn build(root: &Path) -> io::Result<()> {
    run_cargo(root, &["build", "--release", "--workspace", "--locked"])?;

    let cli = root.join("cli");
    let bin_dir = cli.join("bin");
    if bin_dir.exists() {
        fs::remove_dir_all(&bin_dir)?;
    }
    fs::create_dir_all(&bin_dir)?;

    let release_dir = root.join("target/release");
    for (source, destination) in BINARIES {
        let source = release_dir.join(source);
        let destination = bin_dir.join(destination);
        fs::copy(&source, &destination)?;
        mark_executable(&destination)?;
    }

    let entrypoint = cli.join("tipseek");
    if entrypoint.exists() || entrypoint.is_symlink() {
        fs::remove_file(&entrypoint)?;
    }
    make_entrypoint(Path::new("bin/tipseek-rust"), &entrypoint)?;
    write_release_metadata(root, &cli, &bin_dir)
}

fn sha256(path: &Path) -> io::Result<String> {
    let mut file = fs::File::open(path)?;
    let mut digest = Sha256::new();
    let mut buffer = [0_u8; 65_536];
    loop {
        let count = file.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        digest.update(&buffer[..count]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn release_version(root: &Path) -> io::Result<String> {
    let manifest = fs::read_to_string(root.join("rust/tipseek_cli/Cargo.toml"))?;
    manifest
        .lines()
        .map(str::trim)
        .find_map(|line| {
            line.strip_prefix("version = ")
                .and_then(|value| value.trim_matches('"').split_whitespace().next())
                .map(str::to_owned)
        })
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "CLI package version is missing"))
}

fn write_release_metadata(root: &Path, cli: &Path, bin_dir: &Path) -> io::Result<()> {
    let mut files = fs::read_dir(bin_dir)?.collect::<Result<Vec<_>, _>>()?;
    files.sort_by_key(|entry| entry.file_name());
    let checksums = files
        .iter()
        .map(|entry| {
            let name = entry.file_name().to_string_lossy().into_owned();
            Ok(format!("{}  bin/{name}", sha256(&entry.path())?))
        })
        .collect::<io::Result<Vec<_>>>()?;
    fs::write(
        cli.join("SHA256SUMS"),
        format!("{}\n", checksums.join("\n")),
    )?;
    let version = release_version(root)?;
    let mut namespace_digest = Sha256::new();
    namespace_digest.update(version.as_bytes());
    namespace_digest.update(b"\n");
    namespace_digest.update(checksums.join("\n").as_bytes());
    let namespace_hash = format!("{:x}", namespace_digest.finalize());
    let namespace = format!(
        "https://github.com/GUIBA-EX/TipSeek/releases/binary-sbom/{version}/{}",
        namespace_hash
    );

    // This is a deliberately small SPDX 2.3 document for the distributed
    // binaries. Dependency license policy remains enforced by cargo-deny in
    // CI; package-level dependency inventory can be added without changing
    // the release artifact contract.
    let file_entries = files
        .iter()
        .map(|entry| {
            let name = entry.file_name().to_string_lossy().into_owned();
            let checksum = sha256(&entry.path())?;
            Ok(format!(
                "{{\"SPDXID\":\"SPDXRef-File-{name}\",\"fileName\":\"bin/{name}\",\"checksums\":[{{\"algorithm\":\"SHA256\",\"checksumValue\":\"{checksum}\"}}],\"licenseConcluded\":\"GPL-3.0-or-later\",\"copyrightText\":\"NOASSERTION\"}}"
            ))
        })
        .collect::<io::Result<Vec<_>>>()?;
    let sbom = format!(
        "{{\"spdxVersion\":\"SPDX-2.3\",\"dataLicense\":\"CC0-1.0\",\"SPDXID\":\"SPDXRef-DOCUMENT\",\"name\":\"TipSeek binaries\",\"documentNamespace\":\"{namespace}\",\"creationInfo\":{{\"creators\":[\"Tool: TipSeek xtask\"]}},\"documentDescribes\":[\"SPDXRef-Package-TipSeek\"],\"packages\":[{{\"SPDXID\":\"SPDXRef-Package-TipSeek\",\"name\":\"TipSeek\",\"versionInfo\":\"{version}\",\"downloadLocation\":\"NOASSERTION\",\"licenseConcluded\":\"GPL-3.0-or-later\",\"licenseDeclared\":\"GPL-3.0-or-later\",\"copyrightText\":\"NOASSERTION\"}}],\"files\":[{}]}}",
        file_entries.join(",")
    );
    fs::write(cli.join("SBOM.spdx.json"), sbom)
}

fn clean(root: &Path) -> io::Result<()> {
    run_cargo(root, &["clean"])?;
    let cli = root.join("cli");
    let bin_dir = cli.join("bin");
    if bin_dir.exists() {
        fs::remove_dir_all(bin_dir)?;
    }
    let entrypoint = cli.join("tipseek");
    if entrypoint.exists() || entrypoint.is_symlink() {
        fs::remove_file(entrypoint)?;
    }
    Ok(())
}

fn main() -> ExitCode {
    let root = match workspace_root() {
        Ok(root) => root,
        Err(error) => {
            eprintln!("Unable to locate workspace: {error}");
            return ExitCode::FAILURE;
        }
    };
    let result = match env::args().nth(1).as_deref() {
        Some("build") => build(&root),
        Some("clean") => clean(&root),
        _ => {
            usage();
            return ExitCode::FAILURE;
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("xtask failed: {error}");
            ExitCode::FAILURE
        }
    }
}
