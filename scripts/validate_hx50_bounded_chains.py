#!/usr/bin/env python3
"""Verify bounded chain optimization on every HX50 family and time matched subsets."""
import argparse, hashlib, json, os, shutil, subprocess, time, traceback
from pathlib import Path
from validate_hx50_short_boundaries import save
from validate_hx50_chain_replay import summarize


def main():
    p=argparse.ArgumentParser()
    p.add_argument('root',type=Path);p.add_argument('--control',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();control=a.control.resolve();out=a.out.resolve()
    out.mkdir(parents=True,exist_ok=False);binary=out/'gene_workflow';shutil.copy2(a.binary,binary)
    state=dict(status='waiting_for_control',pid=os.getpid(),started_unix=time.time(),stages=[],binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
    save(out/'status.json',state)
    try:
        while True:
            previous=json.loads((control/'status.json').read_text())
            if previous['status']=='failed':raise RuntimeError('control failed: '+str(previous.get('error')))
            if previous['status']=='completed':break
            time.sleep(10)
        report={'scope':'Full-library read audit and cached reconstruction. Exact output/chain comparison to the unbounded-prefix control; paired timing uses the same changed-family subset and eight workers. No reads recruitment or assembly rerun.','backends':{}}
        old=json.loads((control/'report.json').read_text())
        for backend in ['miniprot','exonerate-rs']:
            timings={}
            for label,input_root,extra in [('subset',control/(backend+'-subset'),[]),('full',root/'n-boundary-fix/gene-full',[])]:
                dest=out/(backend+'-'+label)
                cmd=[str(binary),'recover-sco','--annotation',str(control/(backend+'-cache')),'--input',str(input_root),'--out',str(dest),'--threads','8','--sco-backbone-manifest',str(root/'miniprot-backbones.tsv'),'--reads-root',str(root/'workflow')]
                state.update(status='running',phase=backend+'-'+label);save(out/'status.json',state)
                start=time.monotonic()
                with (out/(backend+'-'+label+'.log')).open('w') as log:
                    subprocess.run(['/usr/bin/time','-v',*cmd],stdout=log,stderr=subprocess.STDOUT,check=True)
                timings[label]=time.monotonic()-start
                state['stages'].append(dict(name=state['phase'],seconds=timings[label]));save(out/'status.json',state)
                comparisons={}
                for name in ['sco','sco_masks','chain_audit']:
                    expected=control/(backend+'-all-chains') if name=='chain_audit' else control/backend/name
                    actual=dest/name
                    files=sorted(path for path in actual.rglob('*') if path.is_file())
                    if label=='full':assert {p.relative_to(actual) for p in files}=={p.relative_to(expected) for p in expected.rglob('*') if p.is_file()}
                    for path in files:assert path.read_bytes()==(expected/path.relative_to(actual)).read_bytes(),str(path)
                    comparisons[name]=len(files)
                if label=='full':
                    chains=summarize(dest/'chain_audit')
            before=old['backends'][backend]['changed_family_read_audit_seconds']
            report['backends'][backend]=dict(subset_before_seconds=before,subset_after_seconds=timings['subset'],subset_speedup=before/timings['subset'],full_seconds=timings['full'],byte_identical_files=comparisons,chains=chains)
            save(out/'report.json',report)
        state.update(status='completed',phase='completed',finished_unix=time.time())
    except Exception:
        state.update(status='failed',error=traceback.format_exc());raise
    finally:save(out/'status.json',state)

if __name__=='__main__':main()
