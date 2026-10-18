use gm2_tools::fastx::{FastxFormat, FastxReader};
use std::{collections::BTreeMap, fs, path::Path};
use uce_filter_core::{alignment::align_marker, index::UceIndex};

pub const HEADER: &str = "reference_span_status\treference_record\treference_length\treference_start\treference_end\tcandidate_start\tcandidate_end\tstrand\talignment_identity";
const UNMAPPED: &str = "unmapped\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";

// Report reference-relative boundaries independently of reads support and source.
// The exported interval contains sample bases only; missing ends are never filled.
pub fn annotate(
    index: &UceIndex,
    candidates: &Path,
    output: &Path,
) -> Result<BTreeMap<String, String>, String> {
    let mut reader =
        FastxReader::open(candidates, FastxFormat::Fasta).map_err(|e| e.to_string())?;
    let mut rows = BTreeMap::new();
    let mut fasta = String::new();
    while let Some(record) = reader.next_record().map_err(|e| e.to_string())? {
        let header = String::from_utf8_lossy(&record.header);
        let id = header
            .trim_start_matches('>')
            .split_whitespace()
            .next()
            .ok_or("missing candidate ID")?;
        if record.sequence.len() > u16::MAX as usize {
            rows.insert(
                id.into(),
                "candidate_too_long\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA".into(),
            );
            continue;
        }
        let Some(alignment) = align_marker(index, &record.sequence, 0, 32) else {
            rows.insert(id.into(), UNMAPPED.into());
            continue;
        };
        let reference = &index.references[alignment.sequence as usize];
        let length = reference.bases.len();
        let full = alignment.reference_start == 0 && alignment.reference_end as usize == length;
        let status = if full {
            "full_reference_span"
        } else {
            "partial_reference_span"
        };
        let (start, end) = if alignment.strand == 1 {
            (
                alignment.reference_start as usize,
                alignment.reference_end as usize,
            )
        } else {
            (
                length - alignment.reference_end as usize,
                length - alignment.reference_start as usize,
            )
        };
        rows.insert(
            id.into(),
            format!(
                "{status}\t{}\t{length}\t{start}\t{end}\t{}\t{}\t{}\t{:.6}",
                alignment.sequence / 2 + 1,
                alignment.query_start,
                alignment.query_end,
                if alignment.strand == 1 { "+" } else { "-" },
                alignment.identity()
            ),
        );
        let mut sequence =
            record.sequence[alignment.query_start as usize..alignment.query_end as usize].to_vec();
        if alignment.strand == 2 {
            sequence.reverse();
            for base in &mut sequence {
                *base = match *base {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    b'G' => b'C',
                    other => other,
                };
            }
        }
        fasta.push_str(&format!(
            ">{id} {status}\n{}\n",
            String::from_utf8_lossy(&sequence)
        ));
    }
    fs::write(output, fasta).map_err(|e| e.to_string())?;
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spans_preserve_sample_insertions_and_report_missing_ends() {
        let root = std::env::temp_dir().join(format!("profile_spans_{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let reference =
            "ACGTTGCAAGTCGATCGTACCTGAGCTAGTCGACCTAGGCTACGATGCTAGCTTGACGATCGTACGAGTCCAGTGCATGA";
        let refs = root.join("reference.fasta");
        fs::write(&refs, format!(">ref\n{reference}\n")).unwrap();
        let inserted = format!("{}AA{}", &reference[..35], &reference[35..]);
        let terminal_variant = format!("T{}T", &reference[1..reference.len() - 1]);
        let queries = root.join("queries.fasta");
        fs::write(
            &queries,
            format!(
                ">full\nTTTT{reference}CCCC\n>insertion\n{inserted}\n>terminal_variant\nGG{terminal_variant}CC\n>partial\n{}\n",
                &reference[10..]
            ),
        )
        .unwrap();
        let index = UceIndex::build(&refs, 21).unwrap();
        let output = root.join("spans.fasta");
        let spans = annotate(&index, &queries, &output).unwrap();
        assert!(spans["full"].starts_with("full_reference_span\t"));
        assert!(spans["insertion"].starts_with("full_reference_span\t"));
        assert!(spans["partial"].starts_with("partial_reference_span\t"));
        assert!(spans["terminal_variant"].starts_with("full_reference_span\t"));
        let fasta = fs::read_to_string(output).unwrap();
        assert!(fasta.contains(&format!(">full full_reference_span\n{reference}\n")));
        assert!(fasta.contains(&format!(">insertion full_reference_span\n{inserted}\n")));
        assert!(fasta.contains(&format!(
            ">terminal_variant full_reference_span\n{terminal_variant}\n"
        )));
        assert!(fasta.contains(&format!(
            ">partial partial_reference_span\n{}\n",
            &reference[10..]
        )));
        fs::remove_dir_all(root).unwrap();
    }
}
