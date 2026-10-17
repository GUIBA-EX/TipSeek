#!/usr/bin/env python3
"""Diagnostic only: align recruited/refiltered reads directly to the ppect protein.
Reports reference-relative residue union, never exports a CDS or validates linkage.
Usage: audit_sco_read_coverage.py PANEL OUTPUT MINIPROT
"""
import json, re, subprocess, sys, time
from pathlib import Path
from audit_exonerate_reference import fasta
from benchmark_exon_shadow import digest


def main():
    panel, output, mini = map(Path, sys.argv[1:])
    output.mkdir(parents=True, exist_ok=False)
    results = []
    for item in json.loads((panel / 'panel.json').read_text())['families']:
        family = item['family']
        reads = panel / 'after/1_Hx50/filtered' / (family + '.fq')
        target = output / (family + '.reads.fasta')
        count = 0
        with reads.open() as source, target.open('w') as dest:
            while True:
                header = source.readline()
                if not header:
                    break
                seq, plus, quality = [source.readline().strip() for _ in range(3)]
                if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(quality):
                    raise ValueError('invalid FASTQ: ' + str(reads))
                count += 1
                dest.write(f'>r{count}\n{seq}\n')
        protein = fasta(panel / 'verified-miniprot/manifest/derived_proteins' / (family + '.faa'))
        ppect = [(k, v) for k, v in protein.items() if k.startswith('ppect__')]
        if len(ppect) != 1:
            raise ValueError('expected exactly one ppect reference')
        query = output / (family + '.faa')
        query.write_text(f'>{ppect[0][0]}\n{ppect[0][1]}\n')
        raw = output / (family + '.gff')
        command = [str(mini), '-t', '1', '--gff', '-L', '10', '-n', '2', '--outc', '0', '--outs', '0.2', '-N', str(max(32,count)), '--outn', str(max(32,count)), str(target), str(query)]
        start = time.monotonic()
        with raw.open('w') as dest, (output / (family + '.stderr')).open('w') as err:
            subprocess.run(command, stdout=dest, stderr=err, check=True, timeout=300)
        positions, hits = set(), 0
        for line in raw.read_text().splitlines():
            if not line.startswith('##PAF\t'):
                continue
            fields = line.split('\t')[1:]
            tags = {v.split(':',2)[0]:v.split(':',2)[2] for v in fields[12:] if v.count(':') >= 2}
            if int(fields[10]) == 0 or int(fields[9])/int(fields[10]) < .8 or tags.get('fs','0') != '0' or tags.get('st','0') != '0':
                continue
            q = int(fields[2]); hits += 1
            for size, op in re.findall(r'(\d+)([A-Z])', tags.get('cg','')):
                n = int(size)
                if op == 'M':
                    positions.update(range(q,q+n)); q += n
                elif op == 'I': q += n
                elif op in 'UV': q += 1
            if q != int(fields[3]):
                raise ValueError('unsupported CIGAR query advance')
        results.append({'family':family, 'input':str(reads), 'input_sha256':digest(reads), 'refiltered_reads':count, 'passing_alignments':hits, 'reference_aa':len(ppect[0][1]), 'covered_aa':len(positions), 'coverage':len(positions)/len(ppect[0][1]), 'seconds':time.monotonic()-start, 'command':command})
    report = {'scope':'Direct miniprot alignment of refiltered HX50 reads to ppect protein. A diagnostic union, not proof of sequence accuracy, completeness of one copy, or linkage; no CDS reconstruction from read alignments.', 'miniprot_sha256':digest(mini), 'families':results}
    (output / 'report.json').write_text(json.dumps(report,indent=2)+'\n')
    for r in results:
        print(r['family'],r['refiltered_reads'],r['coverage'])


if __name__ == '__main__':
    main()
