#!/usr/bin/env python3
"""Classify missing reference exon bases without treating ambiguous bases as recovery."""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from report_sco_projection import fasta, rows


def audit(coordinates, recovery, annotation, backbone):
    index = json.loads(coordinates.read_text())
    refs = {r['family_id']: r['protein_reference'] for r in rows(backbone)}
    translation = {(r['family_id'], r['protein_id']): r for r in rows(annotation / 'manifest/reference_translation.tsv') if r['selected'] == '1'}
    bins = defaultdict(Counter)
    families = []
    near_complete = []
    for family, records in sorted(index['records'].items()):
        ref = refs.get(family, 'ppect' if 'ppect' in records else min(records))
        source = index['records'][family][ref.split('__auto')[0]]
        path = recovery / 'sco' / (family + '.fasta')
        seqs = fasta(path) if path.exists() else {}
        assert len(seqs) <= 1
        seq = next(iter(seqs.values()), '')
        meta = translation.get((family, ref), {})
        frame = int(meta.get('frame', '1')) - 1
        reverse = meta.get('strand') == 'reverse'
        masks = rows(recovery / 'sco_masks' / (family + '.tsv')) if seq else []
        counts = Counter()
        for exon in source['exons']:
            start, end = exon['cds_start'], exon['cds_end']
            if reverse:
                start, end = source['cds_length'] - end, source['cds_length'] - start
            start -= frame
            end -= frame
            n = end - start
            label = '<30' if n < 30 else '30-59' if n < 60 else '60-149' if n < 150 else '150-299' if n < 300 else '>=300'
            segment = ''.join(seq[i] if 0 <= i < len(seq) else 'N' for i in range(start, end))
            missing = segment.count('N')
            reasons = Counter()
            for r in masks:
                reasons[r['reason']] += max(0, min(end, int(r['end_0'])) - max(start, int(r['start_0'])))
            reasons['missing_evidence'] += missing - sum(reasons.values())
            assert reasons['missing_evidence'] >= 0
            category = 'complete' if missing == 0 else 'no_observed_bases' if missing == n else 'only_conflicts' if reasons['missing_evidence'] == 0 else 'only_missing_evidence' if reasons['base_conflict'] + reasons['coordinate_conflict'] == 0 else 'mixed_missing_and_conflicts'
            if 0 < missing <= 6 and missing < n:
                near_complete.append({
                    'family': family, 'reference': ref, 'start': start, 'end': end,
                    'missing_nt': missing, 'category': category,
                    'missing_only_at_ends': segment.strip('N').count('N') == 0,
                    'unknown_positions': [start + i for i, base in enumerate(segment) if base == 'N'],
                    'reasons': dict(reasons),
                })
            for counter in [bins[label], bins['all'], counts]:
                counter['exons'] += 1
                counter[category] += 1
                counter['source_nt'] += n
                counter['observed_nt'] += n - missing
                for reason, number in reasons.items():
                    counter[reason + '_nt'] += number
                if 0 < missing <= 6 and missing < n:
                    counter['incomplete_by_at_most_6nt'] += 1
                    counter['near_complete_' + category] += 1
                    counter['near_complete_at_ends' if segment.strip('N').count('N') == 0 else 'near_complete_internal'] += 1
                if 0 < missing < n and segment.strip('N').count('N') == 0:
                    counter['missing_only_at_exon_ends'] += 1

        families.append({'family': family, **counts})
    return {'scope': 'Exact reference exon intervals, not sample exon boundary truth. Categories are mutually exclusive. Missing outputs remain missing.', 'bins': dict(bins), 'families': families, 'near_complete': near_complete}


def main():
    p = argparse.ArgumentParser()
    for name in ['coordinates', 'recovery', 'annotation', 'backbone', 'out']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    result = audit(a.coordinates, a.recovery, a.annotation, a.backbone)
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['bins'], indent=2))


if __name__ == '__main__':
    main()
