#!/usr/bin/env python3
"""Audit revised chains against the same HX50 reads without rerunning aligners."""
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
import traceback
from collections import Counter
from pathlib import Path
from report_full_sco_exons import score
from report_sco_projection import rows, verify


def save(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)


def summarize(path):
    totals = {'before': Counter(), 'after': Counter()}
    improvements = []
    for family in sorted(path.glob('*/*.tsv')):
        best = {'before': 0, 'after': 0}
        states = {'before': Counter(), 'after': Counter()}
        available = set()
        family_rows = rows(family)
        for r in family_rows:
            t = totals[r['mode']]
            available.add(r['reads_available'])
            if r['protein_reference'] == '.':
                continue
            t['chains'] += 1
            t['ambiguous_chains'] += r['ambiguous'] == 'true'
            t['search_limited_chains'] += r['search_limited'] == 'true'
            t['discordant_edges_in_reported_paths'] += int(r['discordant_edges'])
            if r['ambiguous'] == 'false':
                t['unambiguous_chains'] += 1
                t['unambiguous_chains_with_some_paired_support'] += int(r['paired_supported_edges']) > 0
                t['unambiguous_chains_with_all_edges_paired_supported'] += int(r['paired_supported_edges']) == int(r['contigs']) - 1
                best[r['mode']] = max(best[r['mode']], int(r['covered_aa']))
                states[r['mode']]['unambiguous'] += 1
        for mode in totals:
            for key in ['total_pairs', 'unique_mapped_pairs', 'distinct_cross_contig_pairs']:
                totals[mode][key] += int(family_rows[0][key])
            totals[mode]['families_audited'] += 1
            totals[mode]['families_with_reads'] += available == {'true'}
            totals[mode]['families_with_unambiguous_chain'] += states[mode]['unambiguous'] > 0
        if best['before'] != best['after']:
            improvements.append({'family': family.stem, 'best_covered_aa': best})
    return {'totals': totals, 'changed_best_coverage': improvements,
            'scope': 'Before/after compare chain source filters using the same repaired read mapper; previous_read_audit contains the original mapper counts. Chains are hypotheses. Reads warnings veto edges before path selection. No reads spanning a long intron does not reject a chain. This is not a measured biological false-merge rate.'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('root', type=Path)
    p.add_argument('--binary', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--chain-control', type=Path)
    a = p.parse_args()
    out, root = a.out.resolve(), a.root.resolve()
    out.mkdir(parents=True, exist_ok=True)
    binary = out / 'gene_workflow'
    shutil.copy2(a.binary, binary)
    state = {'status': 'running', 'pid': os.getpid(), 'started_unix': time.time(), 'stages': [], 'alignment_reruns': 0, 'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest()}
    save(out / 'status.json', state)
    report = {}
    try:
        for backend in ['miniprot', 'exonerate-rs']:
            dest = out / backend
            old = root / 'n-boundary-fix'
            command = [str(binary), 'recover-sco', '--annotation', str(old / (backend + '-merged-cache')), '--input', str(old / 'gene-full'), '--out', str(dest), '--threads', '8', '--sco-backbone-manifest', str(root / 'miniprot-backbones.tsv'), '--reads-root', str(root / 'workflow')]
            state['phase'] = backend
            save(out / 'status.json', state)
            start = time.monotonic()
            with (out / (backend + '.log')).open('w') as log:
                child = subprocess.Popen(['/usr/bin/time', '-v', *command], stdout=log, stderr=subprocess.STDOUT)
                state['child_pid'] = child.pid
                save(out / 'status.json', state)
                code = child.wait()
            state['stages'].append({'backend': backend, 'seconds': time.monotonic() - start, 'returncode': code, 'command': command})
            save(out / 'status.json', state)
            if code:
                raise RuntimeError(f'{backend} failed: {code}')
            manifest = rows(dest / 'manifest/sco_recovery.tsv')
            verify(dest, manifest)
            exons = score(root / 'exon_coordinates.json', dest, old / (backend + '-merged-cache'), root / 'miniprot-backbones.tsv')
            assert exons['totals']['scorable_families'] == 7006
            # Chain hypotheses cannot silently fill reference gaps or overwrite conflicts.
            control = root / 'fragment-audit-20261016' / (backend + '-after')
            for folder in ['sco', 'sco_masks']:
                for f in (dest / folder).iterdir():
                    assert f.read_bytes() == (control / folder / f.name).read_bytes(), f
            if a.chain_control:
                control_chains = a.chain_control / backend / 'chain_audit'
                actual_files = sorted((dest / 'chain_audit').glob('*/*.tsv'))
                expected_files = sorted(control_chains.glob('*/*.tsv'))
                assert len(actual_files) == len(expected_files) == 6935
                for f in actual_files:
                    assert f.read_bytes() == (control_chains / f.relative_to(dest / 'chain_audit')).read_bytes(), f
            result = summarize(dest / 'chain_audit')
            result['chain_audit_byte_identical_to_control'] = bool(a.chain_control)
            assert result['totals']['after']['families_audited'] == len(manifest) == 6935
            result['exons'] = exons['totals']
            old_reads = rows(old / (backend + '-merged-cache') / 'manifest/read_audit.tsv')
            result['previous_read_audit'] = {key: sum(int(r[key]) for r in old_reads) for key in ['total_pairs', 'unique_mapped_pairs', 'distinct_cross_contig_pairs']}
            assert result['totals']['after']['total_pairs'] == result['previous_read_audit']['total_pairs']
            report[backend] = result
            save(out / 'report.json', report)
        state.update(status='completed', phase='completed', finished_unix=time.time())
    except Exception:
        state.update(status='failed', error=traceback.format_exc())
        raise
    finally:
        save(out / 'status.json', state)


if __name__ == '__main__':
    main()
