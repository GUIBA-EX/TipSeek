#!/usr/bin/env python3
"""Validate independent exon backends on saved real candidates.
Usage: validate_exon_backends.py GENE_WORKFLOW DATASET MINIPROT REPORT.json
exonerate-rs must be in PATH; no Exonerate path is passed to the workflow.
"""
import csv
import hashlib
import json
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path
from benchmark_exon_shadow import outputs, digest


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    binary, data, miniprot, destination = map(Path, sys.argv[1:])
    exo = shutil.which('exonerate-rs')
    if exo is None: raise RuntimeError('exonerate-rs not in PATH')
    root = Path(tempfile.mkdtemp(prefix='tipseek-independent-backends-'))
    report = {'output_root': str(root), 'runs': [], 'summary': [],
              'sha256': {str(p): digest(p) for p in [binary, miniprot, Path(exo)]},
              'scope': 'Independent annotation of real preassembled candidates, including backend-specific padded joins. Not read recruitment/assembly, biological accuracy, or equal candidate-set validation.'}
    expected = {}
    for repeat in range(3):
        for kind in ['intact', 'fragmented']:
            for backend in (['miniprot', 'exonerate-rs'] if repeat % 2 == 0 else ['exonerate-rs', 'miniprot']):
                out = root/f'{repeat}-{kind}-{backend}'
                args = [str(binary.resolve()), 'annotate', '--backend', backend, '--input', str(data/kind/'gene'),
                        '--nucleotide-reference', str(data/'references'), '--out', str(out), '--threads', '1']
                if backend == 'miniprot': args += ['--miniprot', str(miniprot)]
                start = time.monotonic()
                with (root/(out.name+'.log')).open('w') as log:
                    result = subprocess.run(args, stdout=log, stderr=log, timeout=600)
                elapsed = time.monotonic()-start
                if result.returncode: raise RuntimeError(str(out))
                snapshot = outputs(out)
                snapshot["manifest/padded_models.tsv"] = digest(out/"manifest/padded_models.tsv")
                expected.setdefault((kind, backend), snapshot)
                assert snapshot == expected[(kind, backend)], 'nondeterministic output'
                models = rows(out/'models/gene_models.tsv')
                manifest = rows(out/'manifest/candidate_manifest.tsv')
                fragments = rows(out/'manifest/fragment_groups.tsv')
                assert not any('failed' in r['state'] for r in manifest), 'candidate execution/build failure'
                warnings = rows(out/'manifest/annotation_warnings.tsv')
                assert not warnings, warnings
                if backend == 'exonerate-rs':
                    assert not (out/'raw_miniprot').exists()
                    assert all(r['positive'] in ('.', '') for r in models)
                    assert (out/'manifest/backend.tsv').is_file()
                    for f in (out/'raw_exonerate').rglob('status.txt'):
                        assert f.read_text().startswith('ok\n'), str(f)
                accepted = [m for m in models if m['selected']=='1' and m['eligible_for_resolve']=='1']
                report['runs'].append({'repeat': repeat, 'kind': kind, 'backend': backend, 'seconds': elapsed,
                    'command': args, 'output_files': len(snapshot), 'models': len(models), 'accepted_models': len(accepted),
                    'accepted_families': len({m['family_id'] for m in accepted}),
                    'model_states': dict(Counter(m['state'] for m in models)),
                    'padding_states': dict(Counter(f['status'] for f in fragments))})
                print(repeat, kind, backend, round(elapsed,3), 'accepted',len(accepted), flush=True)
    for kind in ['intact','fragmented']:
        times = {b: statistics.median(r['seconds'] for r in report['runs'] if r['kind']==kind and r['backend']==b)
                 for b in ['miniprot','exonerate-rs']}
        report['summary'].append({'kind':kind,'median_seconds':times,'exonerate_speedup':times['miniprot']/times['exonerate-rs']})
    assert all(digest(Path(p))==v for p,v in report['sha256'].items())
    report['passed_functional_checks'] = True
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['summary']),flush=True)


if __name__ == '__main__': main()
