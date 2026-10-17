//! Audit the production short-candidate gate on frozen reads and assemblies.
use gm2_tools::{
    fasta::read_fasta,
    gene_annotation::{read_dna_fasta, reads::PairAudit, screen::ProteinScreen, Interval},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufReader, Write},
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let gene = PathBuf::from(&args[1]);
    let refs = PathBuf::from(&args[2]);
    let reads = PathBuf::from(&args[3]);
    let mut wanted: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for row in csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args[4])?
        .records()
    {
        let r = row?;
        wanted.entry(r[0].into()).or_default().insert(r[1].into());
    }
    let jobs: Vec<_> = wanted.into_iter().collect();
    let next = AtomicUsize::new(0);
    let rows = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..8 {
            handles.push(scope.spawn(|| -> Result<Vec<String>, String> {
                let mut result = Vec::new();
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some((family, ids)) = jobs.get(i) else { break };
                    let run = || -> Result<Vec<String>, Box<dyn std::error::Error>> {
                        let dna: BTreeMap<_,_> = read_dna_fasta(&gene.join(format!("samples/1_Hx50/candidates/{family}.fasta")))?.into_iter().enumerate().map(|(i,r)| (format!("TSK{:06}",i+1),(r.id.split('|').next_back().unwrap().to_owned(),r.sequence))).collect();
                        let proteins: Vec<_> = read_fasta(BufReader::new(File::open(refs.join(format!("{family}.faa")))?))?.into_iter().map(|r| String::from_utf8(r.sequence)).collect::<Result<_,_>>()?;
                        let screen = ProteinScreen::new(proteins.iter().map(String::as_str));
                        let audit = PairAudit::load(Some(&reads),"1_Hx50",family,&dna)?;
                        let mut rows = Vec::new();
                        for (internal,(candidate,sequence)) in &dna {
                            if !ids.contains(candidate) { continue; }
                            let (count,covered) = audit.segment_support(candidate,Interval::new(0,sequence.len()).unwrap());
                            let accepted = screen.supports_read_backed_short(sequence,count,covered);
                            rows.push(format!("{family}\t{candidate}\t{internal}\t{sequence}\t{count}\t{covered}\t{accepted}\n"));
                        }
                        Ok(rows)
                    };
                    result.extend(run().map_err(|e| format!("{family}: {e}"))?);
                }
                Ok(result)
            }));
        }
        let mut rows = Vec::new();
        for h in handles {
            rows.extend(h.join().map_err(|_| "worker panic".to_string())??);
        }
        Ok::<_, String>(rows)
    })?;
    let mut rows = rows;
    rows.sort();
    let mut out = File::create(&args[5])?;
    writeln!(
        out,
        "family_id\tcandidate\tinternal\tsequence\tunique_reads\tcovered_nt\taccepted"
    )?;
    for row in rows {
        out.write_all(row.as_bytes())?;
    }
    fs::write(format!("{}.complete", args[5]), "completed\n")?;
    Ok(())
}
