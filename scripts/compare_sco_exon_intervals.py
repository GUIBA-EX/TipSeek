#!/usr/bin/env python3
"""Report each completed/lost reference interval without hiding N conflicts."""
import argparse, json
from collections import Counter
from pathlib import Path
from report_sco_projection import fasta, rows


def compare(coordinates, annotation, backbone, before, after):
    index=json.loads(coordinates.read_text())
    refs={r['family_id']:r['protein_reference'] for r in rows(backbone)}
    translations={(r['family_id'],r['protein_id']):r for r in rows(annotation/'manifest/reference_translation.tsv') if r['selected']=='1'}
    counts=Counter();changes=[]
    for family,records in sorted(index['records'].items()):
        ref=refs.get(family,'ppect' if 'ppect' in records else min(records))
        source=records[ref.split('__auto')[0]];meta=translations.get((family,ref),{})
        sequences=[]
        for directory in [before,after]:
            path=directory/'sco'/(family+'.fasta')
            records_out=fasta(path) if path.exists() else {}
            assert len(records_out)<=1
            sequences.append(next(iter(records_out.values()),''))
        assert len(sequences[0])==len(sequences[1])
        for number,exon in enumerate(source['exons'],1):
            start,end=exon['cds_start'],exon['cds_end']
            if meta.get('strand')=='reverse':start,end=source['cds_length']-end,source['cds_length']-start
            frame=int(meta.get('frame','1'))-1;start-=frame;end-=frame
            segments=[''.join(seq[i] if 0<=i<len(seq) else 'N' for i in range(start,end)) for seq in sequences]
            missing=[s.count('N') for s in segments]
            counts['reference_exons']+=1
            counts['before_complete']+=missing[0]==0;counts['after_complete']+=missing[1]==0
            if (missing[0]==0)!=(missing[1]==0):
                kind='gained_complete' if missing[1]==0 else 'lost_complete';counts[kind]+=1
                if kind=='gained_complete' and missing[0]<=6:counts['gained_from_at_most_6nt']+=1
                changes.append(dict(family=family,reference=ref,source_exon=number,start=start,end=end,change=kind,before_N=missing[0],after_N=missing[1]))
    return dict(totals=dict(counts),changes=changes)


def main():
    p=argparse.ArgumentParser()
    for name in ['coordinates','annotation','backbone','before','after','out']:p.add_argument(name,type=Path)
    a=p.parse_args();result=compare(a.coordinates,a.annotation,a.backbone,a.before,a.after)
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['totals']))

if __name__=='__main__':main()
