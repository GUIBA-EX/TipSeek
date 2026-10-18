use std::{collections::BTreeSet, fs, process::Command};

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn random_sequence(state: &mut u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(*state >> 32) as usize % 4]
        })
        .collect()
}

#[test]
fn recovers_snp_and_insertion_without_accepting_recombinant() {
    let base = random_sequence(&mut 2718u64, 300);
    let mut snp = base.clone();
    snp[105] = if snp[105] == b'A' { b'C' } else { b'A' };
    let mut insertion = base.clone();
    insertion.splice(190..190, *b"AC");
    run_fixture("snp_indel", vec![base, snp, insertion]);
}

#[test]
fn recovery_is_not_limited_to_sixteen_variants() {
    let base = random_sequence(&mut 2718u64, 300);
    let variants = (0..17)
        .map(|number| {
            let mut sequence = base.clone();
            for bit in 0..5 {
                sequence[100 + bit] = if (number >> bit) & 1 == 0 { b'A' } else { b'C' };
            }
            sequence
        })
        .collect();
    run_fixture("seventeen", variants);
}

fn run_fixture(name: &str, truth: Vec<Vec<u8>>) {
    run_spaced_fixture(name, truth, 1);
}

#[test]
fn recovers_with_sparse_read_start_positions() {
    run_spaced_fixture("sparse", vec![random_sequence(&mut 2718u64, 300)], 13);
}

fn run_spaced_fixture(name: &str, truth: Vec<Vec<u8>>, step: usize) {
    run_panel_fixture(name, truth, step, None);
}

#[test]
fn recovers_a_smaller_disconnected_marker_component() {
    let major = random_sequence(&mut 2718u64, 700);
    let minor = random_sequence(&mut 9817u64, 300);
    run_panel_fixture("disconnected", vec![major, minor.clone()], 3, Some(&minor));
}

fn run_panel_fixture(name: &str, truth: Vec<Vec<u8>>, step: usize, additional: Option<&[u8]>) {
    let root = std::env::temp_dir().join(format!("its2_recovery_{name}_{}", std::process::id()));
    fs::create_dir_all(root.join("ref")).unwrap();
    fs::create_dir_all(root.join("out/filtered")).unwrap();
    let base = &truth[0];
    let mut state = 9371u64;
    let expected_count = truth.len();
    fs::write(
        root.join("ref/marker.fasta"),
        format!(
            ">reference\n{}\n>shorter_panel_entry\n{}\n{}",
            String::from_utf8_lossy(base),
            String::from_utf8_lossy(&base[16..base.len() - 21]),
            additional
                .map(|s| format!(">additional\n{}\n", String::from_utf8_lossy(s)))
                .unwrap_or_default()
        ),
    )
    .unwrap();
    let mut fastq = String::new();
    let flanks: Vec<u8> = (0..80)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 32) as usize % 4]
        })
        .collect();
    let templates: Vec<Vec<u8>> = truth
        .iter()
        .map(|sequence| {
            [
                flanks[..40].to_vec(),
                sequence.clone(),
                flanks[40..].to_vec(),
            ]
            .concat()
        })
        .collect();
    for (variant, sequence) in templates.iter().enumerate() {
        for replicate in 0..3 {
            for start in (0..=sequence.len() - 230).step_by(step) {
                let reads = [
                    sequence[start..start + 150].to_vec(),
                    reverse_complement(&sequence[start + 80..start + 230]),
                ];
                for (mate, read) in reads.iter().enumerate() {
                    fastq.push_str(&format!(
                        "@{variant}_{replicate}_{start}/{}\n{}\n+\n{}\n",
                        mate + 1,
                        String::from_utf8_lossy(read),
                        "I".repeat(read.len())
                    ));
                }
            }
        }
    }
    fs::write(root.join("out/filtered/marker.fq"), fastq).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_main_assembler"))
        .arg("-r")
        .arg(root.join("ref"))
        .arg("-o")
        .arg(root.join("out"))
        .args([
            "--assembly-mode",
            "its2",
            "--assembler-graph-format",
            "gfa",
            "-p",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = fs::read_to_string(root.join("out/results/marker.fasta")).unwrap();
    let recovered: BTreeSet<Vec<u8>> = text
        .split('>')
        .skip(1)
        .map(|record| record.lines().skip(1).flat_map(str::bytes).collect())
        .collect();
    assert_eq!(recovered, truth.into_iter().collect());
    let support = fs::read_to_string(root.join("out/results/marker.its2_support.tsv")).unwrap();
    assert_eq!(
        support
            .lines()
            .filter(|line| line.contains("\tREAD_SUPPORTED\t"))
            .count(),
        expected_count
    );
    assert!(support.lines().all(|line| line.split('\t').count() == 16));
    fs::remove_dir_all(root).unwrap();
}
