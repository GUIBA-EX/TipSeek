#!/usr/bin/env python3
"""Continue real HX50 reads validation from unchanged recruitment, with both backends.
Usage: validate_sco_read_supported_recovery.py PANEL OUTPUT MINIPROT [THREADS] [REFERENCE_ROOT] [SCO_REFERENCE]
"""
import json, os, subprocess, sys, time
from pathlib import Path
from benchmark_exon_shadow import digest


def main():
    panel, out, mini = map(Path, sys.argv[1:4]); threads = sys.argv[4] if len(sys.argv) > 4 else "8"; out.mkdir(parents=True,exist_ok=False)
    reference_root = Path(sys.argv[5]) if len(sys.argv) > 5 else panel / 'references'
    sco_reference = sys.argv[6] if len(sys.argv) > 6 else None
    for reference in (panel / 'references').glob('*.fasta'):
        if digest(reference) != digest(reference_root / reference.name):
            raise ValueError(f'Panel differs from supplied reference: {reference.name}')
    sample = out / '1_Hx50'; sample.mkdir()
    for directory in ['filtered','filtered_pe']:
        (sample / directory).symlink_to(panel / 'after/1_Hx50' / directory, target_is_directory=True)
    workflow = Path('target/release/gene_workflow').resolve()
    assembler = Path('target/release/main_assembler_original').resolve()
    commands = [
        [str(assembler),'-r',str(panel/'references'),'-o',str(sample),'-ka','0','-k_min','21','-k_max','51','-limit_count','2','-iteration','4096','-sb','10000','-cov_min','0','-p','1','--exhaust-seeds'],
        [str(workflow),'classify','--reference',str(panel/'references'),'--contigs',str(sample/'contigs_all'),'--sample','1_Hx50','--out',str(out/'gene'),'--retain-exon-fragments'],
        [str(workflow),'cohort','--reference',str(panel/'references'),'--out',str(out/'gene'),'--sample','1_Hx50'],
    ]
    for backend in ['miniprot','exonerate-rs']:
        command = [str(workflow),'annotate','--backend',backend,'--input',str(out/'gene'),'--nucleotide-reference',str(reference_root),'--out',str(out/backend),'--reads-root',str(out),'--threads',threads]
        if sco_reference: command += ['--sco-reference',sco_reference]
        if backend == 'miniprot': command += ['--miniprot',str(mini)]
        commands.append(command)
    runs = []
    for i, command in enumerate(commands):
        start = time.monotonic()
        with (out/f'run-{i}.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=log, check=True, timeout=1800)
        runs.append({'command':command,'seconds':time.monotonic()-start})
        (out/'runs.json').write_text(json.dumps({'reference_root':str(reference_root),'sco_reference':sco_reference,'assembly_reference_subset_byte_identical':True,'recruitment_reused':str(panel/'after/1_Hx50'), 'binary_sha256':{str(p):digest(p) for p in [workflow,assembler,mini]},'runs':runs},indent=2)+'\n')
        print(i, runs[-1]['seconds'], flush=True)


if __name__ == '__main__': main()
