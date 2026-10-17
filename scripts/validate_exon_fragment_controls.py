#!/usr/bin/env python3
"""Controlled fragmentation of real genomic loci (NOT real reads).
Usage: validate_exon_fragment_controls.py PANEL OUTPUT WORKFLOW MINIPROT
Unknown intron sequence is omitted; each remaining edge retains 120 real bases.
"""
import importlib.util,json,subprocess,sys,time,csv
from pathlib import Path
from audit_exonerate_reference import fasta

def main():
    panel,out,workflow,mini=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
    spec=importlib.util.spec_from_file_location('dataset','/home/yasuda2/tipseek_starfish_multigene_test_20260921/build_dataset.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    genome=Path('/home/yasuda2/starsco/tmp_genomes/ppect.genomic.fna');index=module.read_fai(Path(str(genome)+'.fai'))
    refs=out/'references';inputs=out/'gene/samples/ppect/candidates';refs.mkdir(parents=True);inputs.mkdir(parents=True)
    truth=[]
    with genome.open('rb') as handle:
        for item in json.loads((panel/'panel.json').read_text())['families']:
            e=item['exons'];start=e[0][1];end=e[-1][2]
            sequence=module.fetch(handle,index,e[0][0],start+1,end)
            cuts=[(a[2]-start+120,b[1]-start-120) for a,b in zip(e,e[1:]) if b[1]-a[2]>=20000]
            parts=[];last=0
            for left,right in cuts:parts.append(sequence[last:left]);last=right
            parts.append(sequence[last:]);family=item['family']
            (inputs/(family+'.fasta')).write_text(''.join(f'>part_{i+1}\n{s}\n' for i,s in enumerate(parts)))
            cds=fasta(Path(item['reference']))['ppect'];(refs/(family+'.fasta')).write_text(f'>ppect\n{cds}\n')
            truth.append({'family':family,'contigs':len(parts),'expected_cds':cds})
    (out/'truth.json').write_text(json.dumps(truth,indent=2)+'\n')
    report={'scope':'Controlled multi-contig fragmentation of real genomic loci; source annotation is not experimental truth','runs':[]}
    for backend in ['miniprot','exonerate-rs']:
        result=out/backend;args=[str(workflow.resolve()),'annotate','--backend',backend,'--input',str(out/'gene'),'--nucleotide-reference',str(refs),'--out',str(result),'--threads','1']
        if backend=='miniprot':args+=['--miniprot',str(mini)]
        start=time.monotonic()
        with (out/(backend+'.log')).open('w') as log:subprocess.run(args,stdout=log,stderr=log,check=True,timeout=1800)
        def rows(name):return list(csv.DictReader((result/'manifest'/name).open(),delimiter='\t'))
        report['runs'].append({'backend':backend,'seconds':time.monotonic()-start,'command':args,'fragment_groups':rows('fragment_groups.tsv'),'fragment_chains':rows('fragment_chains.tsv')})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
