#![cfg(unix)]
use std::{fs, os::unix::fs::PermissionsExt, process::Command};

#[test]
fn standalone_backend_uses_path_and_validates_padding_without_miniprot() {
    let root = std::env::temp_dir().join(format!(
        "tipseek-backend-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let bin = root.join("bin");
    let input = root.join("gene");
    let refs = root.join("refs");
    let proteins = root.join("proteins");
    fs::create_dir_all(&bin).unwrap();
    fs::create_dir_all(input.join("samples/s/candidates")).unwrap();
    fs::create_dir_all(&refs).unwrap();
    fs::create_dir_all(&proteins).unwrap();
    fs::write(
        input.join("samples/s/candidates/f.fasta"),
        ">left\nATGAAAGT\n>right\nAGTTTCCCTAA\n",
    )
    .unwrap();
    fs::write(refs.join("f.fasta"), ">p\nATGAAATTTCCCTAA\n").unwrap();
    fs::write(proteins.join("f.faa"), ">p\nMKFP\n").unwrap();
    let executable = bin.join("exonerate-rs");
    let script = "#!/bin/sh\nif [ \"$1\" = --version ]; then echo exonerate-rs-test; exit 0; fi\nfor arg do target=$arg; done\ncase \"$target\" in\n*/TSK000001/*) echo 'vulgar: p 0 2 . TSK000001 0 6 + 20 M 2 6';;\n*/TSK000002/*) echo 'vulgar: p 2 4 . TSK000002 2 8 + 20 M 2 6';;\n*/TSKPAD001/*) echo 'vulgar: p 0 4 . TSKPAD001 0 116 + 40 M 2 6 5 0 2 I 0 100 3 0 2 M 2 6';;\nesac\n";
    fs::write(&executable, script).unwrap();
    fs::set_permissions(&executable, fs::Permissions::from_mode(0o755)).unwrap();
    let run = |output: &std::path::Path, threads: &str| {
        Command::new(env!("CARGO_BIN_EXE_gene_workflow"))
            .env("PATH", &bin) // Only the fake Exonerate is discoverable.
            .args(["annotate", "--backend", "exonerate-rs", "--input"])
            .arg(&input)
            .arg("--nucleotide-reference")
            .arg(&refs)
            .arg("--protein-reference")
            .arg(&proteins)
            .arg("--out")
            .arg(output)
            .args(["--threads", threads])
            .output()
            .unwrap()
    };
    let out = root.join("ok");
    let result = run(&out, "1");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let parallel = root.join("parallel");
    assert!(run(&parallel, "4").status.success());
    for relative in [
        "models/gene_models.tsv",
        "models/gene_segments.tsv",
        "manifest/fragment_groups.tsv",
        "manifest/fragment_chains.tsv",
        "cds/f.fasta",
    ] {
        assert_eq!(
            fs::read(out.join(relative)).unwrap(),
            fs::read(parallel.join(relative)).unwrap(),
            "{relative}"
        );
    }
    assert!(!out.join("raw_miniprot").exists());
    assert!(fs::read_to_string(out.join("manifest/fragment_groups.tsv"))
        .unwrap()
        .contains("padded_validated"));
    assert!(fs::read_to_string(out.join("cds/f.fasta"))
        .unwrap()
        .contains("ATGAAATTTCCC"));
    assert!(fs::read_to_string(out.join("models/gene_models.tsv"))
        .unwrap()
        .contains("\t1.000000\t.\t"));
    // The source fragments remain clean, but reannotation inserts a stop at
    // the padded path end. It must not supersede them with a complete model.
    fs::write(
        &executable,
        script.replace(
            "TSKPAD001 0 116 + 40 M 2 6 5 0 2 I 0 100 3 0 2 M 2 6",
            "TSKPAD001 0 119 + 40 M 2 6 5 0 2 I 0 100 3 0 2 M 2 6 G 0 3",
        ),
    )
    .unwrap();
    let rejected = root.join("padded_stop");
    assert!(run(&rejected, "4").status.success());
    assert!(
        fs::read_to_string(rejected.join("manifest/fragment_groups.tsv"))
            .unwrap()
            .contains("padded_validation_failed")
    );
    let audit = fs::read_to_string(rejected.join("manifest/padded_models.tsv")).unwrap();
    assert!(audit.contains("internal_stop"));
    assert!(audit.contains("structurally_ineligible;not_complete"));
    assert!(!fs::read_to_string(rejected.join("models/gene_models.tsv"))
        .unwrap()
        .contains("superseded_by_padded_join"));
    for (name, body, expected) in [
        ("nohit", "exit 0", "no_exonerate_model"),
        ("failed", "exit 3", "exonerate_failed"),
    ] {
        fs::write(&executable, format!("#!/bin/sh\nif [ \"$1\" = --version ]; then echo exonerate-rs-test; exit 0; fi\n{body}\n")).unwrap();
        let out = root.join(name);
        assert!(run(&out, "1").status.success());
        assert!(
            fs::read_to_string(out.join("manifest/candidate_manifest.tsv"))
                .unwrap()
                .contains(expected)
        );
    }
    fs::write(&executable, "#!/bin/sh\nif [ \"$1\" = --version ]; then echo exonerate-rs-test; exit 0; fi\necho invalid-output\n").unwrap();
    assert!(!run(&root.join("malformed"), "4").status.success());
    fs::remove_file(&executable).unwrap();
    let preserved = root.join("preserved");
    fs::create_dir_all(&preserved).unwrap();
    fs::write(preserved.join("sentinel"), "keep").unwrap();
    assert!(!run(&preserved, "4").status.success());
    assert_eq!(
        fs::read_to_string(preserved.join("sentinel")).unwrap(),
        "keep"
    );
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn three_contigs_include_short_middle_exon_and_cannot_skip_it() {
    let root = std::env::temp_dir().join(format!("tipseek-three-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    let bin = root.join("bin");
    let input = root.join("gene");
    let refs = root.join("refs");
    let proteins = root.join("proteins");
    for d in [&bin, &refs, &proteins, &input.join("samples/s/candidates")] {
        fs::create_dir_all(d).unwrap();
    }
    fs::write(
        input.join("samples/s/candidates/f.fasta"),
        format!(
            ">left\n{}GT\n>middle\nAG{}GT\n>right\nAG{}\n",
            "GCT".repeat(45),
            "GCT".repeat(10),
            "GCT".repeat(45)
        ),
    )
    .unwrap();
    fs::write(refs.join("f.fasta"), format!(">p\n{}\n", "GCT".repeat(100))).unwrap();
    fs::write(proteins.join("f.faa"), format!(">p\n{}\n", "A".repeat(100))).unwrap();
    let script="#!/bin/sh\nif [ \"$1\" = --version ]; then echo exonerate-rs-test; exit 0; fi\nfor arg do target=$arg; done\ncase \"$target\" in\n*/TSK000001/*) echo 'vulgar: p 0 45 . TSK000001 0 135 + 100 M 45 135';;\n*/TSK000002/*) echo 'vulgar: p 45 55 . TSK000002 2 32 + 50 M 10 30';;\n*/TSK000003/*) echo 'vulgar: p 55 100 . TSK000003 2 137 + 100 M 45 135';;\n*/TSKPAD001/*) echo 'vulgar: p 0 100 . TSKPAD001 0 508 + 250 M 45 135 5 0 2 I 0 100 3 0 2 M 10 30 5 0 2 I 0 100 3 0 2 M 45 135';;\nesac\n";
    let executable = bin.join("exonerate-rs");
    fs::write(&executable, script).unwrap();
    fs::set_permissions(&executable, fs::Permissions::from_mode(0o755)).unwrap();
    let run = |out: &std::path::Path| {
        Command::new(env!("CARGO_BIN_EXE_gene_workflow"))
            .env("PATH", &bin)
            .args(["annotate", "--backend", "exonerate-rs", "--input"])
            .arg(&input)
            .arg("--nucleotide-reference")
            .arg(&refs)
            .arg("--protein-reference")
            .arg(&proteins)
            .arg("--out")
            .arg(out)
            .output()
            .unwrap()
    };
    let out = root.join("ok");
    let result = run(&out);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(fs::read_to_string(out.join("manifest/fragment_groups.tsv"))
        .unwrap()
        .contains("padded_validated"));
    let models = fs::read_to_string(out.join("models/gene_models.tsv")).unwrap();
    assert_eq!(models.matches("superseded_by_padded_join").count(), 3);
    assert!(fs::read_to_string(out.join("cds/f.fasta"))
        .unwrap()
        .contains(&"GCT".repeat(100)));
    assert!(fs::read_to_string(out.join("fragment_cds/f.fasta"))
        .unwrap()
        .contains("reference_ordered_unvalidated"));
    fs::write(
        &executable,
        script.replace(
            "5 0 2 I 0 100 3 0 2 M 10 30 5 0 2 I 0 100 3 0 2",
            "5 0 2 I 0 234 3 0 2 G 10 0",
        ),
    )
    .unwrap();
    let bad = root.join("skipped");
    assert!(run(&bad).status.success());
    assert!(fs::read_to_string(bad.join("manifest/padded_models.tsv"))
        .unwrap()
        .contains("source_fragment_not_recovered"));
    assert!(!fs::read_to_string(bad.join("models/gene_models.tsv"))
        .unwrap()
        .contains("superseded_by_padded_join"));
    fs::remove_dir_all(root).unwrap();
}
