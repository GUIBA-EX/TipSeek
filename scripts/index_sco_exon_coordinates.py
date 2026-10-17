#!/usr/bin/env python3
"""Match each SCO species CDS exactly to source GFF; never infer exon boundaries.
Usage: index_sco_exon_coordinates.py REFS SOURCE_CDS SOURCE_GFF OUTPUT.json
Coordinates are concatenated source-CDS offsets, before selected frame trimming.
"""
import hashlib,json,sys
from collections import defaultdict
from pathlib import Path

def fasta(path):
    header=None;sequence=[]
    with path.open() as stream:
        for line in stream:
            if line.startswith('>'):
                if header is not None:yield header,''.join(sequence).upper()
                header=line[1:].strip();sequence=[]
            else:sequence.append(line.strip())
    if header is not None:yield header,''.join(sequence).upper()

def digest(sequence):return hashlib.sha256(sequence.encode()).hexdigest()

def main():
    refs,cds,gff,out=map(Path,sys.argv[1:])
    wanted=defaultdict(lambda:defaultdict(list))
    targets={}
    for path in sorted(refs.glob('*.fasta')):
        for header,sequence in fasta(path):
            species=header.split()[0]
            key=(path.stem,species)
            assert key not in targets
            targets[key]=sequence
            wanted[species][digest(sequence)].append(path.stem)
    output={'reference_families':len(list(refs.glob('*.fasta'))),'coordinate_scope':'Exact source CDS sequence to GFF exon lengths; cross-species reference proxies, not HX50 truth.','records':{},'species':{}}
    for species,hashes in sorted(wanted.items()):
        source=cds/(species+'.cds.fa');features=gff/(species+'.gff3')
        if not source.exists() or not features.exists():
            output['species'][species]={'status':'missing_source'};continue
        matches=defaultdict(set)
        for header,sequence in fasta(source):
            families=hashes.get(digest(sequence),[])
            for family in families:
                if sequence==targets[(family,species)]:
                    # Source IDs often encode names and the GFF transcript ID after |.
                    matches[family].add(header.split()[0].split('|')[-1])
        parents={p for names in matches.values() for p in names}
        intervals=defaultdict(list)
        with features.open() as stream:
            for line in stream:
                if line.startswith('#'):continue
                f=line.rstrip().split('\t')
                if len(f)!=9 or f[2]!='CDS':continue
                attrs=dict(item.split('=',1) for item in f[8].split(';') if '=' in item)
                for parent in attrs.get('Parent','').split(','):
                    if parent in parents:intervals[parent].append([f[0],int(f[3])-1,int(f[4]),f[6],f[7]])
        usable=0
        for family,names in sorted(matches.items()):
            valid=[]
            for name in sorted(names):
                exons=intervals.get(name,[])
                if not exons or len({(e[0],e[3]) for e in exons})!=1:continue
                exons=sorted(exons,key=lambda e:e[1],reverse=exons[0][3]=='-')
                if sum(e[2]-e[1] for e in exons)!=len(targets[(family,species)]):continue
                lengths=tuple(e[2]-e[1] for e in exons)
                valid.append((name,exons,lengths))
            if not valid or len({v[2] for v in valid})!=1:continue
            name,exons,lengths=valid[0];offset=0;coding=[]
            for e,length in zip(exons,lengths):
                coding.append({'cds_start':offset,'cds_end':offset+length,'genomic':e});offset+=length
            output['records'].setdefault(family,{})[species]={'transcript':name,'equivalent_transcripts':[v[0] for v in valid],'cds_length':offset,'exons':coding,'reference_sequence_sha256':digest(targets[(family,species)])}
            usable+=1
        output['species'][species]={'status':'indexed','exact_CDS_matches':len(matches),'unambiguous_exon_coordinate_records':usable,'cds_file':str(source),'gff_file':str(features)}
        print(species,usable,flush=True)
    out.write_text(json.dumps(output,separators=(',',':'))+'\n')

if __name__=='__main__':main()
