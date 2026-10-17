#!/usr/bin/env python3
"""Explain independent backend acceptance and padding differences.

Usage: audit_exon_backend_differences.py VALIDATION.json REPORT.json [ORIGINAL_EXONERATE]
Optional original replay uses identical saved arguments; no changes to alignments.
Source annotations and either backend are not treated as biological truth.
"""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from benchmark_exon_shadow import digest
from audit_exonerate_reference import fasta, without_stop


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    validation, destination = map(Path, sys.argv[1:3])
    root = Path(json.loads(validation.read_text())['output_root'])
    report = {'validation': str(validation), 'differences': [], 'shared_padded_cds': [], 'original_replay': []}
    replay = Path(tempfile.mkdtemp(prefix='tipseek-original-replay-')) if len(sys.argv) == 4 else None
    if replay:
        report['original'] = {'executable': sys.argv[3], 'sha256': digest(Path(sys.argv[3])), 'outputs': str(replay)}
    for kind in ['intact', 'fragmented']:
        data = {}
        for backend in ['miniprot', 'exonerate-rs']:
            out = root / f'0-{kind}-{backend}'
            data[backend] = {
                'models': rows(out/'models/gene_models.tsv'),
                'segments': rows(out/'models/gene_segments.tsv'),
                'padding': rows(out/'manifest/fragment_groups.tsv'),
                'padded_models': rows(out/'manifest/padded_models.tsv'),
            }
        for family in sorted({r['family_id'] for d in data.values() for r in d['models']}):
            by_backend = {}
            for backend, d in data.items():
                models = [m for m in d['models'] if m['family_id'] == family]
                by_backend[backend] = {
                    'accepted': [m for m in models if m['selected'] == m['eligible_for_resolve'] == '1'],
                    'source_models': [m for m in models if m['candidate'] != 'padded_join_1'],
                    'source_splice_defects': [r for r in d['segments'] if r['family_id'] == family and r['candidate'] != 'padded_join_1' and r['kind'] == 'intron' and r['splice_class'] in ('noncanonical', 'unavailable')],
                    'padding': [r for r in d['padding'] if r['family_id'] == family],
                    'padded_models': [r for r in d['padded_models'] if r['family_id'] == family],
                }
            signature = lambda d: ([(m['candidate'], m['state'], m['query_start_aa'], m['query_end_aa']) for m in d['accepted']], [r['status'] for r in d['padding']])
            if signature(by_backend['miniprot']) != signature(by_backend['exonerate-rs']):
                report['differences'].append({'kind': kind, 'family': family, 'backends': by_backend})
        shared = set.intersection(*[
            {r['family_id'] for r in d['padding'] if r['status'] == 'padded_validated'}
            for d in data.values()
        ])
        for family in sorted(shared):
            sequences = {}
            for backend in data:
                records = fasta(root/f'0-{kind}-{backend}'/'cds'/f'{family}.fasta')
                sequences[backend] = next(seq for name, seq in records.items() if '|padded_join_1|' in name)
            mini, exo = sequences['miniprot'], sequences['exonerate-rs']
            report['shared_padded_cds'].append({'kind': kind, 'family': family,
                'exact_equal': mini == exo, 'terminal_stop_normalized_equal': without_stop(mini) == without_stop(exo)})
        if replay:
            saved = root / f'0-{kind}-exonerate-rs/raw_exonerate'
            for target in sorted(saved.rglob('target.fasta')):
                directory = target.parent
                relative = directory.relative_to(saved)
                out = replay/kind/relative
                out.mkdir(parents=True)
                args = (directory/'arguments.txt').read_text().splitlines()[1:]
                result = subprocess.run([sys.argv[3], *args], capture_output=True, timeout=60)
                (out/'stdout.txt').write_bytes(result.stdout)
                (out/'stderr.txt').write_bytes(result.stderr)
                if result.returncode:
                    raise RuntimeError(f'original failed: {directory}')
                report['original_replay'].append({'kind': kind, 'case': str(relative),
                    'equal_vulgar': result.stdout == (directory/'stdout.txt').read_bytes()})
    destination.write_text(json.dumps(report, indent=2)+'\n')
    print('different family summaries:', len(report['differences']))
    if replay:
        print('original exact:', sum(r['equal_vulgar'] for r in report['original_replay']), '/', len(report['original_replay']))


if __name__ == '__main__':
    main()
