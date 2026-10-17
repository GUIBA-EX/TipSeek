#!/usr/bin/env python3
"""Score complete reference exon intervals from exact source CDS/GFF matches.
Usage: report_full_sco_exons.py COORDINATES RECOVERY ANNOTATION OUT.json
Every absent output and clipped source base is missing, never counted complete.
"""
import csv,json,sys
from collections import Counter
from pathlib import Path
from report_sco_projection import fasta,rows

def score(coordinates,recovery,annotation,backbone_manifest=None):
    index=json.loads(coordinates.read_text())
    audits={r['family_id']:r for r in rows(recovery/'manifest/sco_recovery.tsv')}
    translations={(r['family_id'],r['protein_id']):r for r in rows(annotation/'manifest/reference_translation.tsv') if r['selected']=='1'}
    backbones={r["family_id"]:r["protein_reference"] for r in rows(backbone_manifest)} if backbone_manifest else {}
    totals=Counter();details=[];short=Counter()
    for family,records in sorted(index['records'].items()):
        audit=audits.get(family)
        if audit is None:
            reference=backbones.get(family)
            species=reference.split('__auto')[0] if reference else ('ppect' if 'ppect' in records else min(records))
            source=records.get(species)
            if source is None:totals['families_with_unscorable_backbone']+=1;continue
            reference=reference or species;sequence='';frame=0;reverse=False
            totals['families_without_output']+=1
        else:
            reference=audit['protein_reference'];species=reference.split('__auto')[0]
            if family in backbones:assert reference==backbones[family],(family,reference,backbones[family])
            source=records.get(species)
            meta=translations.get((family,reference))
            if source is None or meta is None or meta['source']!='derived_nucleotide':
                totals['families_with_unscorable_backbone']+=1;continue
            frame=int(meta['frame'])-1;reverse=meta['strand']=='reverse'
            assert int(meta['nucleotide_length'])==source['cds_length']
            exported=fasta(recovery/'sco'/(family+'.fasta'));assert len(exported)==1
            sequence=next(iter(exported.values()));assert len(sequence)==int(audit['length_nt'])
            totals['families_with_output']+=1
        totals['scorable_families']+=1
        mask_path=recovery/'sco_masks'/(family+'.tsv')
        mask_rows=rows(mask_path) if audit and mask_path.is_file() else None
        if mask_rows is not None:
            total_masked=0;previous=0
            for r in mask_rows:
                x,y=int(r['start_0']),int(r['end_0'])
                assert 0<=previous<=x<y<=len(sequence)
                assert sequence[x:y]=='N'*(y-x)
                assert r['reason'] in ['missing_evidence','base_conflict','coordinate_conflict']
                total_masked+=y-x;previous=y
            assert total_masked==sequence.count('N')
            totals['families_with_mask_reason_audit']+=1
        family_full=0
        for exon in source['exons']:
            start,end=exon['cds_start'],exon['cds_end']
            if reverse:start,end=source['cds_length']-end,source['cds_length']-start
            start-=frame;end-=frame
            length=end-start
            a=max(0,start);b=min(len(sequence),end)
            segment=sequence[a:b] if b>a else ''
            observed=len(segment)-segment.count('N')
            complete=observed==length
            totals['reference_exons']+=1
            totals['complete_reference_exon_intervals']+=complete
            totals['exons_with_any_observed_bases']+=observed>0
            totals['exons_with_no_observed_bases']+=observed==0
            totals['source_exon_nt']+=length;totals['observed_source_exon_nt']+=observed
            if length<150:
                short['reference_exons']+=1;short['complete_reference_exon_intervals']+=complete;short['exons_with_any_observed_bases']+=observed>0
            family_full+=complete
            if mask_rows is not None or audit is None:
                reasons=Counter()
                if audit is None:reasons['missing_evidence']=length
                else:
                    reasons['missing_evidence']+=length-len(segment)
                    for r in mask_rows:
                        overlap=max(0,min(end,int(r['end_0']))-max(start,int(r['start_0'])))
                        reasons[r['reason']]+=overlap
                assert sum(reasons.values())==length-observed
                totals['exons_audited_for_mask_reasons']+=1
                for reason in ['missing_evidence','base_conflict','coordinate_conflict']:
                    totals['exons_with_'+reason]+=reasons[reason]>0
                    totals[reason+'_nt']+=reasons[reason]
                if reasons['missing_evidence']==0 and reasons['coordinate_conflict']==0:
                    totals['fully_located_exon_intervals_including_base_ambiguity']+=1
        totals['families_with_all_source_exons_observed']+=family_full==len(source['exons'])
        details.append({'family':family,'reference':reference,'reference_exons':len(source['exons']),'complete_reference_exon_intervals':family_full,'all_exons_observed':family_full==len(source['exons'])})
    totals['families_without_any_source_exon_index']=index['reference_families']-len(index['records'])
    return {'scope':'Complete non-N source exon coordinate intervals; exact CDS/GFF reference match. These cross-species intervals are not HX50 boundary truth or same-copy linkage. Missing output and clipped bases count as missing. Unscorable source annotations are excluded and disclosed.','total_reference_families':index['reference_families'],'recovery':str(recovery),'totals':dict(totals),'short_exons_below_150nt':dict(short),'complete_exon_interval_fraction':totals['complete_reference_exon_intervals']/totals['reference_exons'],'families':details}

def main():
    coordinates,recovery,annotation,out=map(Path,sys.argv[1:5]);manifest=Path(sys.argv[5]) if len(sys.argv)>5 else None;result=score(coordinates,recovery,annotation,manifest)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='families'},indent=2))

if __name__=='__main__':main()
