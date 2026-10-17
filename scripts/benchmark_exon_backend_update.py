#!/usr/bin/env python3
"""Paired standalone Exonerate timing and accepted-sequence regression check.
Usage: benchmark_exon_backend_update.py BEFORE AFTER DATASET REPORT.json
The two binaries must already be built; exonerate-rs is taken from PATH.
"""
import json
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from benchmark_exon_shadow import digest


def accepted_outputs(root):
    # Models/manifests may intentionally gain corrected defect counts.
    return {str(p.relative_to(root)): digest(p) for name in ['cds', 'proteins', 'supercontigs']
            for p in (root/name).rglob('*') if p.is_file()}


def main():
    before, after, data, destination = map(Path, sys.argv[1:])
    exo = shutil.which('exonerate-rs')
    if not exo:
        raise RuntimeError('exonerate-rs must be in PATH')
    binaries = {'before': before.resolve(), 'after': after.resolve()}
    root = Path(tempfile.mkdtemp(prefix='tipseek-backend-update-'))
    report = {'output_root': str(root), 'runs': [], 'summary': [],
              'sha256': {str(p): digest(p) for p in [before, after, Path(exo)]},
              'scope': 'Five paired repetitions of single-thread annotation on this machine; includes filesystem and child-process overhead.'}
    expected = {}
    for repeat in range(5):
        for kind in ['intact', 'fragmented']:
            for label in (['before', 'after'] if repeat % 2 == 0 else ['after', 'before']):
                out = root/f'{repeat}-{kind}-{label}'
                args = [str(binaries[label]), 'annotate', '--backend', 'exonerate-rs',
                        '--input', str(data/kind/'gene'), '--nucleotide-reference', str(data/'references'),
                        '--out', str(out), '--threads', '1']
                usage0 = resource.getrusage(resource.RUSAGE_CHILDREN)
                start = time.monotonic()
                with (root/(out.name+'.log')).open('w') as log:
                    process = subprocess.Popen(args, stdout=log, stderr=log)
                    try:
                        # Blocking wait avoids Python's timed wait polling quantization.
                        code = process.wait()
                    finally:
                        if process.poll() is None:
                            process.kill()
                            process.wait()
                elapsed = time.monotonic()-start
                if code:
                    raise RuntimeError(str(out))
                usage1 = resource.getrusage(resource.RUSAGE_CHILDREN)
                snapshot = accepted_outputs(out)
                expected.setdefault(kind, snapshot)
                assert snapshot == expected[kind], f'accepted outputs changed: {out}'
                report['runs'].append({'repeat': repeat, 'kind': kind, 'binary': label, 'seconds': elapsed,
                    'cpu_seconds': usage1.ru_utime+usage1.ru_stime-usage0.ru_utime-usage0.ru_stime,
                    'command': args})
                print(repeat, kind, label, round(elapsed, 4), flush=True)
    for kind in ['intact', 'fragmented']:
        times = {label: statistics.median(r['seconds'] for r in report['runs'] if r['kind']==kind and r['binary']==label) for label in binaries}
        cpu = {label: statistics.median(r['cpu_seconds'] for r in report['runs'] if r['kind']==kind and r['binary']==label) for label in binaries}
        report['summary'].append({'kind': kind, 'median_seconds': times, 'median_cpu_seconds': cpu,
                                  'speedup': times['before']/times['after']})
    assert all(digest(Path(p)) == h for p, h in report['sha256'].items())
    report['accepted_outputs_unchanged'] = True
    destination.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    main()
