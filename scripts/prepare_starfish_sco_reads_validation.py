#!/usr/bin/env python3
"""Prepare a deterministic long-intron SCO panel for the user's real HX50 reads.
Annotation is from P. pectinifera, NOT ground truth for the HX50 individual.
Usage: prepare_starfish_sco_reads_validation.py OUTPUT_DIRECTORY
"""
import csv, json, shutil, sys
from pathlib import Path
from audit_exonerate_reference import fasta, attributes

def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
    sco=Path('/home/yasuda2/starsco/geneminer_sco_refs/fasta')
    cds=fasta(Path('/home/yasuda2/starsco/geneminer_usco_refs/cds/ppect.cds.fa'))
    lookup={s:n for n,s in cds.items()}
    # Some source FASTA IDs include annotations after a pipe.
    features={}
    for line in Path('/home/yasuda2/starsco/annevo_gff/ppect.gff3').open():
        if line.startswith('#'):continue
        f=line.rstrip().split('\t')
        if len(f)==9 and f[2]=='CDS':
            parent=attributes(f[8]).get('Parent')
            features.setdefault(parent,[]).append((f[0],int(f[3])-1,int(f[4]),f[6],f[7]))
    candidates=[]
    for path in sorted(sco.glob('*.fasta')):
        records=fasta(path)
        sequence=records.get('ppect')
        if not sequence:continue
        name=lookup.get(sequence)
        if name is None:continue
        name=name.split('|')[-1]
        exons=sorted(features.get(name,[]),key=lambda x:x[1])
        if len(exons)<3:continue
        introns=[b[1]-a[2] for a,b in zip(exons,exons[1:])]
        if max(introns)<20000:continue
        if not 300<=len(sequence)<=4500:continue
        candidates.append({'family':path.stem,'transcript':name,'max_intron':max(introns),'long_introns':sum(n>=20000 for n in introns),'exons':exons,'cds_nt':len(sequence),'shortest_exon':min(e[2]-e[1] for e in exons),'reference':str(path)})
    # Selection uses reference structure only, never recovery outcomes.
    candidates.sort(key=lambda x:(-x['long_introns'],-x['max_intron'],x['family']))
    selected=candidates[:12]
    if not selected:raise RuntimeError('No source annotation matches; inspect IDs')
    refs=out/'references';refs.mkdir(exist_ok=True)
    for r in selected:shutil.copyfile(r['reference'],refs/(r['family']+'.fasta'))
    reads=Path('/home/yasuda2/staruce/test/01_fastp')
    (out/'samples.tsv').write_text(f'HX50\t{reads}/HX50_1.clean.fq.gz\t{reads}/HX50_2.clean.fq.gz\n')
    (out/'panel.json').write_text(json.dumps({'sample':'HX50','species':'Choriaster granulatus','reads_kind':'real fastp-cleaned paired-end sequencing reads','reference_annotation_species':'Patiria pectinifera','reference_is_sample_truth':False,'selection':'first 12 by long intron count, maximum intron length, family ID; CDS 300–4500 nt and >=3 exons','families':selected},indent=2)+'\n')
    print([(x['family'],x['max_intron'],x['long_introns']) for x in selected])

if __name__=='__main__':main()
