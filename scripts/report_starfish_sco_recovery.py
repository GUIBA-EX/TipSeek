#!/usr/bin/env python3
"""Report comparable read-derived reference coverage, not sample biological truth.
Usage: report_starfish_sco_recovery.py PANEL BEFORE AFTER EXONERATE REPORT.json
Coverage uses M (miniprot) / M,S (Exonerate) query positions; query gaps never count.
"""
import csv,json,re,sys,statistics
from pathlib import Path
from audit_exonerate_reference import fasta
from benchmark_exon_shadow import digest

def rows(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter='\t'))

def covered_intervals(path,backend):
    hits=[]
    for line in path.read_text().splitlines():
        if backend=='miniprot':
            if not line.startswith('##PAF\t'):continue
            f=line.split('\t')[1:]
            if not f[0].startswith('ppect__'):continue
            tags={x.split(':',2)[0]:x.split(':',2)[2] for x in f[12:] if x.count(':')>=2}
            if tags.get('fs','0')!='0' or tags.get('st','0')!='0':continue
            q=int(f[2]);target_start=int(f[7]);positions=set()
            for n,op in re.findall(r'(\d+)([A-Z])',tags.get('cg','')):
                n=int(n)
                if op=='M':positions.update(range(q,q+n));q+=n
                elif op=='I':q+=n
                elif op in 'UV':q+=1  # split codon; exclude it conservatively here
            if q!=int(f[3]):raise ValueError(f'unsupported query advance in {path}: {tags.get("cg")}')
            hits.append((f[0],int(f[1]),positions,target_start,int(f[8]),f[4]))
        else:
            if not line.startswith('vulgar:'):continue
            f=line.split()
            if not f[1].startswith('ppect__'):continue
            q=int(f[2]);positions=set();bad=False
            for i in range(10,len(f),3):
                op,qa,ta=f[i],int(f[i+1]),int(f[i+2])
                if op=='F':bad=True
                if op in ('M','S'):positions.update(range(q,q+qa))
                q+=qa
            if not bad:hits.append((f[1],None,positions,min(int(f[6]),int(f[7])),max(int(f[6]),int(f[7])),f[8]))
    return hits

def summarize(root,panel,backend):
    models=rows(root/'models/gene_models.tsv');manifest=rows(root/'manifest/candidate_manifest.tsv');details=[]
    for item in panel['families']:
        family=item['family'];candidates=[m for m in models if m['family_id']==family and m['candidate']!='padded_join_1']
        # Identity/QC are checked on each normalized model; no requirement to be
        # independently resolve-eligible (which would exclude short exon evidence).
        valid=[m for m in candidates if m['protein_reference'].startswith('ppect__') and float(m['identity'])>=0.8 and m['frameshifts']=='0' and m['internal_stops']=='0' and m['state'] not in ('phase_inconsistent','ambiguous_cds','internal_stop','frameshifted') and not any(f in m['qc_flags'] for f in ['noncanonical_splice','invalid_splice_length','missing_paf'])]
        ids=rows(root/'manifest/id_map.tsv');mapping={x['internal_candidate']:x['candidate'] for x in ids if x['family_id']==family}
        rawname='raw_miniprot' if backend=='miniprot' else 'raw_exonerate'
        paths=(root/rawname).glob(f'*/{family}/*.gff3' if backend=='miniprot' else f'*/{family}/*.vulgar')
        covered=set();longest=0
        for path in paths:
            candidate=mapping.get(path.stem)
            if candidate=='padded_join_1':continue
            for protein,length,positions,start,end,strand in covered_intervals(path,backend):
                matched=[m for m in valid if m['candidate']==candidate and m['protein_reference']==protein and m['strand']==strand and abs(int(m['genomic_start'])-1-start)<=3 and abs(int(m['genomic_end'])-end)<=3]
                if not matched:continue
                covered.update(positions);longest=max(longest,len(positions))
        source_cds=fasta(Path(item['reference']))['ppect'];n=len(source_cds)//3-(source_cds[-3:] in ('TAA','TAG','TGA'));covered={p for p in covered if 0<=p<n}
        exons=sorted(item['exons'],key=lambda x:x[1],reverse=item['exons'][0][3]=='-');offset=0;recovered=0;total=0
        for exon in exons:
            length=exon[2]-exon[1];expected=set(range(offset//3,min(n,(offset+length+2)//3)));offset+=length
            if not expected:continue
            total+=1;recovered+=len(covered & expected)/len(expected)>=0.8
        details.append({'family':family,'source_annotation_max_intron_nt':item['max_intron'],'retained_candidates':len({m['candidate'] for m in manifest if m['family_id']==family and m['candidate']!='padded_join_1'}),'reference_exons':total,'exons_with_at_least_80_percent_aligned_residues':recovered,'aligned_reference_aa':len(covered),'reference_length_aa':n,'cds_union_coverage':len(covered)/n,'longest_single_alignment_coverage':longest/n})
    chains=rows(root/'manifest/fragment_chains.tsv') if (root/'manifest/fragment_chains.tsv').exists() else []
    accepted=[m for m in models if m['selected']=='1' and m['eligible_for_resolve']=='1']
    states={s:sum(c['status']==s for c in chains) for s in sorted({c['status'] for c in chains})}
    return {'output':str(root),'backend':backend,'families':details,'reference_exon_recovery_fraction':sum(r['exons_with_at_least_80_percent_aligned_residues'] for r in details)/sum(r['reference_exons'] for r in details),'mean_cds_union_coverage':statistics.mean(r['cds_union_coverage'] for r in details),'accepted_models':len(accepted),'complete_models':sum(m['state']=='complete' for m in accepted),'chain_states':states,'exported_reference_guided_cds':sum(len(fasta(p)) for p in (root/'fragment_cds').glob('*.fasta')),'biological_wrong_merge_rate':None,'wrong_merge_rate_reason':'No same-individual genome or validated linkage truth for HX50; ambiguity rejection is not a measured false-merge rate.'}

def main():
    panel_path,before,after,exo,destination=map(Path,sys.argv[1:]);panel=json.loads((panel_path/'panel.json').read_text())
    report={'sample':panel['sample'],'species':panel['species'],'input':str(panel_path/'samples.tsv'),'reference_species':panel['reference_annotation_species'],'metric_scope':'Reference-relative proxy on real HX50 reads. 80% query-residue recovery per annotated P. pectinifera exon; CDS completeness excludes query gaps. This is not same-individual biological truth. Multi-contig coverage is a union of evidence, not proof all pieces belong to one copy.','before':summarize(before,panel,'miniprot'),'after':summarize(after,panel,'miniprot'),'exonerate':summarize(exo,panel,'exonerate-rs'),'baseline_runs':json.loads((panel_path/'baseline_runs.json').read_text()),'panel':panel}
    binaries=[Path('/tmp/tipseek-before-multicontig'),Path('/tmp/tipseek-before-exon-seeds'),Path('target/release/gene_workflow'),Path('target/release/main_assembler_original'),Path('target/release/tipseek'),Path('/home/yasuda2/exonerate/target/release/exonerate-rs'),Path('/home/yasuda2/miniforge3/pkgs/miniprot-0.18-h577a1d6_0/bin/miniprot')]
    report['binary_sha256']={str(p.resolve()):digest(p) for p in binaries if p.exists()}
    report['input_files']=[{'path':f,'bytes':Path(f).stat().st_size,'mtime_ns':Path(f).stat().st_mtime_ns} for f in (panel_path/'samples.tsv').read_text().strip().split('\t')[1:]]
    report['annotation_timing']={label:Path(str(root)+'.time').read_text() for label,root in [('miniprot',after),('exonerate-rs',exo)] if Path(str(root)+'.time').exists()}
    report['source_sha256']={str(p):digest(p) for p in sorted(Path('rust').rglob('*.rs'))}
    destination.write_text(json.dumps(report,indent=2)+'\n')
    for label in ['before','after','exonerate']:
        r=report[label];print(label,r['reference_exon_recovery_fraction'],r['mean_cds_union_coverage'],r['chain_states'])
if __name__=='__main__':main()
