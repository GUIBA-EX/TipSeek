#!/usr/bin/env python3
"""Run previous assembly/classification/annotation on the same freshly recruited reads.
Usage: run_starfish_sco_baseline.py PANEL OLD_ASSEMBLER OLD_WORKFLOW MINIPROT
Filtering is unchanged and was run from real FASTQ by the new complete workflow.
"""
import json,subprocess,sys,time
from pathlib import Path

def main():
    root,assembler,workflow,mini=map(Path,sys.argv[1:]);refs=root/'references';out=root/'before';sample=out/'1_Hx50';sample.mkdir(parents=True,exist_ok=True)
    (sample/'filtered').symlink_to(root/'after/1_Hx50/filtered',target_is_directory=True)
    commands=[
        [str(assembler),'-r',str(refs),'-o',str(sample),'-ka','0','-k_min','21','-k_max','51','-limit_count','2','-iteration','4096','-sb','10000','-cov_min','0','-p','1'],
        [str(workflow),'classify','--reference',str(refs),'--contigs',str(sample/'contigs_all'),'--sample','Hx50','--out',str(out/'gene')],
        [str(workflow),'cohort','--reference',str(refs),'--out',str(out/'gene'),'--sample','Hx50'],
        [str(workflow),'annotate','--input',str(out/'gene'),'--nucleotide-reference',str(refs),'--out',str(out/'exon'),'--miniprot',str(mini),'--threads','8'],
    ]
    runs=[]
    for number,args in enumerate(commands):
        t=time.monotonic()
        with (root/f'baseline-{number}.log').open('w') as log:subprocess.run(args,check=True,stdout=log,stderr=log,timeout=1800)
        runs.append({'command':args,'seconds':time.monotonic()-t})
    (root/'baseline_runs.json').write_text(json.dumps({'shared_recruitment':str(root/'after/1_Hx50/filtered'),'runs':runs},indent=2)+'\n')
if __name__=='__main__':main()
