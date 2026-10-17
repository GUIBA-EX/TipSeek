#!/usr/bin/env python3
"""Replay both frozen HX50 backends against one backbone and compare repairs."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
import traceback
from pathlib import Path
from audit_sco_exon_gaps import audit
from report_full_sco_exons import score
from report_sco_projection import rows, fasta, verify, summarize


def save(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('root', type=Path)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--control-binary', type=Path, help='Frozen executable for the before control')
    a = p.parse_args()
    root, out = a.root.resolve(), a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    old = root / 'n-boundary-fix'
    binary = out / 'gene_workflow'
    shutil.copy2(a.binary, binary)
    state = {'status': 'running', 'pid': os.getpid(), 'started_unix': time.time(), 'stages': [],
             'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
             'alignment_reruns': 0, 'reads_recruitment_reruns': 0,
             'scope': 'Reuses read-derived corrected assemblies and all frozen alignments. Not an independent sequencing truth or new reads recruitment run.'}
    save(out / 'status.json', state)
    report = {}
    try:
        for backend in ['miniprot', 'exonerate-rs']:
            annotation = old / (backend + '-merged-cache')
            controls = {}
            for label, executable in [('before', a.control_binary or old / 'components/gene_workflow'), ('after', binary)]:
                dest = out / (backend + '-' + label)
                command = [str(executable), 'recover-sco', '--annotation', str(annotation), '--input', str(old / 'gene-full'), '--out', str(dest), '--threads', '8', '--sco-backbone-manifest', str(root / 'miniprot-backbones.tsv')]
                state['phase'] = backend + '-' + label
                save(out / 'status.json', state)
                start = time.monotonic()
                with (out / (state['phase'] + '.log')).open('w') as log:
                    run = subprocess.Popen(['/usr/bin/time', '-v', *command], stdout=log, stderr=subprocess.STDOUT)
                    state['child_pid'] = run.pid
                    save(out / 'status.json', state)
                    code = run.wait()
                seconds = time.monotonic() - start
                state['stages'].append({'name': state['phase'], 'seconds': seconds, 'returncode': code, 'command': command})
                save(out / 'status.json', state)
                if code:
                    raise RuntimeError(f'{state["phase"]} failed: {code}')
                manifest = rows(dest / 'manifest/sco_recovery.tsv')
                verify(dest, manifest)
                exons = score(root / 'exon_coordinates.json', dest, annotation, root / 'miniprot-backbones.tsv')
                gaps = audit(root / 'exon_coordinates.json', dest, annotation, root / 'miniprot-backbones.tsv')
                assert exons['totals']['reference_exons'] == gaps['bins']['all']['exons'] == 82712
                assert exons['totals']['complete_reference_exon_intervals'] == gaps['bins']['all']['complete']
                save(out / (state['phase'] + '-exons.json'), exons)
                save(out / (state['phase'] + '-gaps.json'), gaps)
                controls[label] = {'seconds': seconds, 'sco': summarize(manifest, 7006), 'exons': exons['totals'], 'gap_bins': gaps['bins'], 'dominated_models': sum(int(r['dominated_models']) for r in manifest)}
            changes = []
            for path in sorted((out / (backend + '-before') / 'sco').glob('*.fasta')):
                before = next(iter(fasta(path).values()))
                after = next(iter(fasta(out / (backend + '-after') / 'sco' / path.name).values()))
                assert len(before) == len(after)
                if before != after:
                    changes.append({'family': path.stem, 'gained_nt': sum(x == 'N' and y != 'N' for x, y in zip(before, after)), 'lost_nt': sum(x != 'N' and y == 'N' for x, y in zip(before, after)), 'substituted_nt': sum(x != 'N' and y != 'N' and x != y for x, y in zip(before, after))})
            controls['changes'] = changes
            report[backend] = controls
            save(out / 'report.json', report)
        state.update(status='completed', phase='completed', finished_unix=time.time())
    except Exception:
        state.update(status='failed', error=traceback.format_exc())
        raise
    finally:
        save(out / 'status.json', state)


if __name__ == '__main__':
    main()
