#!/usr/bin/env python3
"""Count observed contigs once, separately from duplicated reference/model hits."""
import csv
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path


def audit(annotation):
    families = defaultdict(dict)
    with (annotation / 'manifest/id_map.tsv').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['internal_candidate'].startswith('TSKPAD'):
                continue
            families[(r['sample'], r['family_id'])][r['candidate']] = [int(r['sequence_length']), 0]
    names = ['no_manifest_evidence', 'screen_rejected', 'no_alignment_or_failed', 'aligned_below_80pct', 'high_identity_with_structural_flags', 'high_identity_without_reported_structural_flags']
    with (annotation / 'manifest/candidate_manifest.tsv').open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            c = families[(r['sample'], r['family_id'])].get(r['candidate'])
            if c is None:
                continue
            if r['state'] == 'protein_screen_rejected':
                rank = 1
            elif not r['model']:
                rank = 2
            elif float(r['identity']) < .8:
                rank = 3
            else:
                flags = r['qc_flags'].split(';')
                clean = int(r['frameshifts']) == 0 and int(r['internal_stops']) == 0 and r['state'] not in ['frameshifted', 'internal_stop', 'phase_inconsistent', 'ambiguous_cds'] and not any(x in flags for x in ['noncanonical_splice', 'invalid_splice_length', 'missing_paf', 'padded_fragment_join'])
                rank = 5 if clean else 4
            c[1] = max(c[1], rank)
    totals = defaultdict(Counter)
    for candidates in families.values():
        for length, rank in candidates.values():
            label = '<30' if length < 30 else '30-59' if length < 60 else '60-149' if length < 150 else '150-299' if length < 300 else '>=300'
            totals[label][names[rank]] += 1
            totals['all'][names[rank]] += 1
    return {'scope': 'Unique original contigs; high identity without reported structural flags is not proof of exon truth or an unambiguous coordinate projection.', 'bins': dict(totals)}


if __name__ == '__main__':
    result = audit(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
