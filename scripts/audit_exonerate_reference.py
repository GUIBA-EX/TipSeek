#!/usr/bin/env python3
"""Adjudicate intact alignments against source GFF/CDS, not another aligner.

Usage: audit_exonerate_reference.py DATASET SHADOW_OUTPUT SOURCE_GFF REPORT.tsv
The source annotation is a reference, not experimental biological truth.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys


def fasta(path):
    records = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            records[name] = ''
        elif line.strip():
            records[name] += line.strip().upper()
    return records


def reverse(sequence):
    return sequence.translate(str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN'))[::-1]


def cds(sequence, intervals, strand):
    result = ''.join(sequence[a:b] for a, b in sorted(intervals))
    return reverse(result) if strand == '-' else result


def without_stop(sequence):
    return sequence[:-3] if sequence[-3:] in ('TAA', 'TAG', 'TGA') else sequence


def attributes(text):
    return dict(f.split('=', 1) for f in text.split(';') if '=' in f)


def miniprot(path):
    models = {}
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if len(f) != 9 or line.startswith('#'):
            continue
        attr = attributes(f[8])
        if f[2] == 'mRNA':
            models[attr['ID']] = [f[6], []]
        elif f[2] == 'CDS':
            models[attr['Parent']][1].append((int(f[3])-1, int(f[4])))
    return list(models.values())


def vulgar(path):
    models = []
    for line in path.read_text().splitlines():
        f = line.split()
        if not f or f[0] != 'vulgar:':
            raise ValueError('invalid vulgar output')
        direction = 1 if f[8] == '+' else -1
        position = start = int(f[6])
        intervals = []
        for i in range(10, len(f), 3):
            next_position = position + direction * int(f[i+2])
            if f[i] in ('5', 'I', '3'):
                if start != position:
                    intervals.append(tuple(sorted((start, position))))
                start = next_position
            position = next_position
        if start != position:
            intervals.append(tuple(sorted((start, position))))
        models.append((f[8], intervals))
    return models


def introns(intervals):
    ordered = sorted(intervals)
    return [(a[1], b[0]) for a, b in zip(ordered, ordered[1:])]


def main():
    dataset, output, gff, report = map(Path, sys.argv[1:])
    truth = list(csv.DictReader((dataset/'truth/models.tsv').open(), delimiter='\t'))
    wanted = {r['transcript']: [] for r in truth}
    with gff.open() as handle:
        for line in handle:
            f = line.rstrip().split('\t')
            if len(f) == 9 and f[2] == 'CDS':
                parent = attributes(f[8]).get('Parent')
                if parent in wanted:
                    wanted[parent].append((int(f[3])-1, int(f[4])))
    genomic_targets = {k.split("|")[0]: v for k,v in fasta(dataset/"truth/genomic_targets.fasta").items()}
    rows = []
    source_checks = []
    for item in truth:
        family = item['family']
        dirs = list((output/'raw_exonerate_shadow').glob(f'*/{family}/TSK000001'))
        assert len(dirs) == 1
        directory = dirs[0]
        target = next(iter(fasta(directory/'target.fasta').values()))
        assert target == genomic_targets[family], family
        source = sorted(wanted[item['transcript']])
        origin = source[0][0]
        expected = [(a-origin, b-origin) for a, b in source]
        reference = next(iter(fasta(dataset/'references'/f'{family}.fasta').values()))
        assert cds(target, expected, item['strand']) == reference, family
        source_checks.append(family)
        relative = directory.relative_to(output/'raw_exonerate_shadow')
        alignments = {'miniprot': miniprot((output/'raw_miniprot'/relative).with_suffix('.gff3')),
                      'exonerate_rs': vulgar(directory/'stdout.txt')}
        expected_positions = {p for a,b in expected for p in range(a,b)}
        for engine, models in alignments.items():
            for number, (strand, intervals) in enumerate(models, 1):
                positions = {p for a,b in intervals for p in range(a,b)}
                motifs = []
                for a,b in introns(intervals):
                    sequence = target[a:b] if strand == '+' else reverse(target[a:b])
                    motifs.append(sequence[:2]+'-'+sequence[-2:])
                missed = [b-a for a,b in expected if not positions.intersection(range(a,b))]
                row = {'family': family, 'engine': engine, 'model': number, 'strand': strand,
                       'source_strand': item['strand'], 'same_strand': strand == item['strand'],
                       'source_cds_match_without_stop': strand == item['strand'] and without_stop(cds(target, intervals, strand)) == without_stop(reference),
                       'source_splice_boundaries_match': strand == item['strand'] and introns(intervals) == introns(expected),
                       'source_coding_nt_missing': len(expected_positions-positions),
                       'extra_coding_nt': len(positions-expected_positions),
                       'missed_exon_lengths': ','.join(map(str,missed)),
                       'missed_microexon_lengths': ','.join(str(n) for n in missed if n <= 9),
                       'splice_motifs': ','.join(motifs),
                       'cds_intervals': ','.join(f'{a}-{b}' for a,b in sorted(intervals)),
                       'source_intervals': ','.join(f'{a}-{b}' for a,b in expected)}
                rows.append(row)
    with report.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    summary = {'scope': 'Intact candidates only; concordance with source annotation, not biological ground truth. Both engines used CDS-derived protein queries, so agreement is not independent annotation evidence.',
               'source_cds_genome_checks': len(source_checks), 'source_gff_sha256': hashlib.sha256(gff.read_bytes()).hexdigest(),
               'command': sys.argv, 'engines': {}}
    for engine in ['miniprot', 'exonerate_rs']:
        subset = [r for r in rows if r['engine'] == engine]
        summary['engines'][engine] = {key: len({r['family'] for r in subset if r[key]}) for key in ['source_cds_match_without_stop', 'source_splice_boundaries_match']}
        summary['engines'][engine]['same_strand_models_missing_microexons'] = sum(r['same_strand'] and bool(r['missed_microexon_lengths']) for r in subset)
    report.with_suffix('.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary['engines']))


if __name__ == '__main__':
    main()
