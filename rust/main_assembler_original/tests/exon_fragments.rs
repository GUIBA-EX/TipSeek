use std::{fs, process::Command};

#[test]
fn exhaustive_exon_search_retains_disconnected_short_region() {
    let root = std::env::temp_dir().join(format!("tipseek-seeds-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    let refs = root.join("refs");
    fs::create_dir_all(&refs).unwrap();
    let mut state = 3917u64;
    let mut dna = |length| {
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state % 4) as usize] as char
            })
            .collect::<String>()
    };
    let pieces = [dna(400), dna(60), dna(120)];
    // The middle read-supported region has no exact nucleotide reference seed.
    fs::write(
        refs.join("f.fasta"),
        format!(">ref\n{}{}\n", pieces[0], pieces[2]),
    )
    .unwrap();
    let mut reads = String::new();
    for (region, sequence) in pieces.iter().enumerate() {
        for repeat in 0..5 {
            for start in (0..=sequence.len() - 40).step_by(5) {
                reads.push_str(&format!(
                    ">r{region}_{repeat}_{start}\n{}\n",
                    &sequence[start..start + 40]
                ));
            }
        }
    }
    let mut recovered = Vec::new();
    for exhaust in [false, true] {
        let out = root.join(if exhaust { "all" } else { "legacy" });
        fs::create_dir_all(out.join("filtered")).unwrap();
        fs::write(out.join("filtered/f.fasta"), &reads).unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_main_assembler_original"));
        command
            .arg("-r")
            .arg(&refs)
            .arg("-o")
            .arg(&out)
            .args(["-ka", "21", "-sb", "10000"]);
        if exhaust {
            command.arg("--exhaust-seeds");
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let contigs = fs::read_to_string(out.join("contigs_all/f.fasta")).unwrap();
        recovered.push(
            pieces
                .iter()
                .filter(|p| contigs.contains(p.as_str()))
                .count(),
        );
    }
    assert!(
        recovered[0] < 3,
        "legacy test must demonstrate the lost region"
    );
    assert_eq!(recovered[1], 3);
    fs::remove_dir_all(root).unwrap();
}
