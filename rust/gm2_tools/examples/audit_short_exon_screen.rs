//! Diagnostic only: count short observed candidates below the existing 10-aa screen.
use gm2_tools::{
    fasta::read_fasta,
    gene_annotation::{read_dna_fasta, reverse_complement, screen::ProteinScreen, translate_cds},
};
use std::{
    collections::HashSet,
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::PathBuf,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let candidates = PathBuf::from(&args[1]);
    let references = PathBuf::from(&args[2]);
    let mut out = BufWriter::new(File::create(&args[3])?);
    writeln!(out, "family_id\tcandidate\tsequence\texact_informative_8aa")?;
    let mut paths: Vec<_> = fs::read_dir(candidates)?
        .map(|p| p.unwrap().path())
        .collect();
    paths.sort();
    for path in paths {
        let family = path.file_stem().unwrap().to_str().unwrap();
        let reference = references.join(format!("{family}.faa"));
        if !reference.is_file() {
            continue;
        }
        let proteins: Vec<_> = read_fasta(BufReader::new(File::open(reference)?))?
            .into_iter()
            .map(|r| String::from_utf8(r.sequence).unwrap())
            .collect();
        let screen = ProteinScreen::new(proteins.iter().map(String::as_str));
        let mut exact: HashSet<[u8; 8]> = HashSet::new();
        for protein in &proteins {
            for window in protein.as_bytes().windows(8) {
                if window.iter().all(|b| b"ACDEFGHIKLMNPQRSTVWY".contains(b))
                    && window
                        .iter()
                        .fold(0u32, |mask, b| mask | (1 << (b - b'A')))
                        .count_ones()
                        >= 5
                {
                    exact.insert(window.try_into().unwrap());
                }
            }
        }
        for record in read_dna_fasta(&path)? {
            let dna = record.sequence;
            if !(30..60).contains(&dna.len()) || screen.supports(&dna) {
                continue;
            }
            let reverse = reverse_complement(&dna);
            let mut matched = None;
            'oriented: for oriented in [dna.as_str(), reverse.as_str()] {
                for frame in 0..3 {
                    let end = frame + (oriented.len() - frame) / 3 * 3;
                    let Ok(protein) = translate_cds(&oriented[frame..end]) else {
                        continue;
                    };
                    for window in protein.as_bytes().windows(8) {
                        let word: [u8; 8] = window.try_into().unwrap();
                        if exact.contains(&word) {
                            matched = Some(String::from_utf8(word.to_vec())?);
                            break 'oriented;
                        }
                    }
                }
            }
            if let Some(word) = matched {
                writeln!(
                    out,
                    "{family}\t{}\t{dna}\t{word}",
                    record.id.split('|').next_back().unwrap()
                )?;
            }
        }
    }
    out.flush()?;
    Ok(())
}
