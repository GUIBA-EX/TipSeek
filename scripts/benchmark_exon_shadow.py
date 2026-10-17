#!/usr/bin/env python3
"""Alternate old/new binaries on intact/fragmented real data with exact output checks.
Usage: benchmark_exon_shadow.py OLD NEW DATASET MINIPROT REPORT.json
exonerate-rs is resolved exclusively through the inherited PATH.
"""
import hashlib
import json
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def outputs(root):
    paths = []
    for directory in ['cds', 'proteins', 'exons', 'introns', 'genes', 'genes_flanked',
                      'supercontigs', 'gff3', 'unresolved', 'models']:
        paths.extend(p for p in (root/directory).rglob('*') if p.is_file())
    paths.extend(root/'manifest'/name for name in ['candidate_manifest.tsv', 'fragment_groups.tsv',
                 'reference_translation.tsv', 'annotation_warnings.tsv', 'id_map.tsv'])
    return {str(p.relative_to(root)): digest(p) for p in paths}


def main():
    old, new, data, miniprot, report_path = map(Path, sys.argv[1:])
    exonerate = shutil.which('exonerate-rs')
    if exonerate is None:
        raise RuntimeError('exonerate-rs is not in the current PATH')
    root = Path(tempfile.mkdtemp(prefix='tipseek-exon-benchmark-'))
    binaries = {'before': old.resolve(), 'after': new.resolve()}
    report = {'output_root': str(root), 'runs': [], 'comparisons': [],
              'sha256': {str(p): digest(p) for p in [old, new, miniprot, Path(exonerate)]},
              'resolved_exonerate': exonerate,
              'scope': 'Three alternating repetitions per engine/mode/dataset on this machine; no universal speed or biological-accuracy claim.'}
    expected = {}
    shadow_expected = {}
    for repeat in range(3):
        for kind in ['intact', 'fragmented']:
            for mode in ['baseline', 'shadow']:
                for label in (['before', 'after'] if repeat % 2 == 0 else ['after', 'before']):
                    out = root/f'{repeat}-{kind}-{mode}-{label}'
                    args = [str(binaries[label]), 'annotate', '--input', str(data/kind/'gene'),
                            '--nucleotide-reference', str(data/'references'), '--out', str(out),
                            '--miniprot', str(miniprot), '--threads', '1', '--max-intron', '50000',
                            '--minimum-coverage', '0.20', '--complete-coverage', '0.80', '--flank', '0',
                            '--fragment-padding', '100']
                    if mode == 'shadow': args.append('--exonerate-shadow')
                    start = time.monotonic()
                    with (root/(out.name+'.log')).open('w') as log:
                        result = subprocess.run(args, stdout=log, stderr=log, timeout=600)
                    elapsed = time.monotonic()-start
                    if result.returncode: raise RuntimeError(f'{out}: exit {result.returncode}')
                    snapshot = outputs(out)
                    expected.setdefault(kind, snapshot)
                    if snapshot != expected[kind]: raise AssertionError(f'authoritative output differs: {out}')
                    if mode == 'shadow':
                        shadow = {str(p.relative_to(out)): digest(p) for p in (out/'raw_exonerate_shadow').rglob('*')
                                  if p.is_file() and p.name in ['stdout.txt', 'normalized_models.tsv']}
                        shadow_expected.setdefault(kind, shadow)
                        if shadow != shadow_expected[kind]: raise AssertionError(f'shadow output differs: {out}')
                    report['runs'].append({'repeat': repeat, 'kind': kind, 'mode': mode, 'binary': label,
                                           'seconds': elapsed, 'compared_files': len(snapshot), 'command': args})
                    print(repeat, kind, mode, label, round(elapsed, 3), flush=True)
    for kind in ['intact', 'fragmented']:
        for mode in ['baseline', 'shadow']:
            medians = {label: statistics.median(r['seconds'] for r in report['runs'] if
                       r['kind']==kind and r['mode']==mode and r['binary']==label) for label in binaries}
            report['comparisons'].append({'kind': kind, 'mode': mode, 'median_seconds': medians,
                                          'speedup': medians['before']/medians['after']})
    assert all(digest(Path(path)) == value for path, value in report['sha256'].items()), 'binary changed during benchmark'
    report['passed_output_equivalence'] = True
    report_path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['comparisons']), flush=True)


if __name__ == '__main__':
    main()
