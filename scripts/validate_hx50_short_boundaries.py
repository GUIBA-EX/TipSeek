#!/usr/bin/env python3
"""Incremental HX50 validation: read-backed short candidates and observed stops.
Caches are hard-linked into a new directory; changed files are atomically replaced.
Original alignments are never modified. Only newly screened candidates are aligned.
"""
import argparse, csv, hashlib, json, os, shutil, subprocess, time, traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from audit_sco_exon_gaps import audit
from report_full_sco_exons import score
from report_sco_projection import rows, fasta, verify
from validate_hx50_chain_replay import summarize


def save(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('root', type=Path)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--gate', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    root, out = a.root.resolve(), a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = out / 'gene_workflow'
    shutil.copy2(a.binary, binary)
    state = dict(status='running', pid=os.getpid(), started_unix=time.time(), stages=[],
                 binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
    save(out/'status.json', state)
    old = root/'n-boundary-fix'
    baseline = root/'fragment-audit-20261016/optimized-chains'
    gate = [r for r in rows(a.gate) if r['accepted'] == 'true']
    env = dict(os.environ, PATH=str(root/'components')+os.pathsep+os.environ.get('PATH',''))
    report = {'scope': '7006 SCO; same backbone and reads-derived assemblies. Only read-backed short candidates realigned; unchanged chains reused. Timings are incremental stages, not full reads-to-annotation runtime.', 'gated_candidates': len(gate), 'backends': {}}
    def run(command, label):
        state['phase'] = label; save(out/'status.json',state)
        start = time.monotonic()
        with (out/(label+'.log')).open('w') as log:
            subprocess.run(['/usr/bin/time','-v',*map(str,command)], stdout=log, stderr=subprocess.STDOUT, check=True, env=env)
        seconds = time.monotonic()-start
        state['stages'].append(dict(name=label, seconds=seconds)); save(out/'status.json',state)
        return seconds
    try:
        for backend in ['miniprot','exonerate-rs']:
            cache = out/(backend+'-cache')
            subprocess.run(['cp','-al',str(old/(backend+'-merged-cache')),str(cache)],check=True)
            manifest = cache/'manifest/candidate_manifest.tsv'
            gate_keys={(r['family_id'],r['candidate']) for r in gate}
            with manifest.open() as f:
                reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames
                rejected={(r['family_id'],r['candidate']) for r in reader if r['state']=='protein_screen_rejected' and (r['family_id'],r['candidate']) in gate_keys}
            selected=[r for r in gate if (r['family_id'],r['candidate']) in rejected]
            selected_keys={(r['family_id'],r['candidate']) for r in selected}
            with (cache/'manifest/id_map.tsv').open() as f:
                cached_ids={(r['family_id'],r['candidate']):(r['internal_candidate'],int(r['sequence_length'])) for r in csv.DictReader(f,delimiter='\t') if (r['family_id'],r['candidate']) in selected_keys}
            for r in selected:
                assert cached_ids[(r['family_id'],r['candidate'])] == (r['internal'],len(r['sequence']))
            work=out/(backend+'-new-alignments'); work.mkdir()
            def align(r):
                family,internal=r['family_id'],r['internal']
                target=work/(family+'_'+internal+'.fasta')
                target.write_text('>'+internal+'\n'+r['sequence']+'\n')
                protein=cache/'manifest/derived_proteins'/(family+'.faa')
                if backend=='miniprot':
                    cmd=['miniprot','--gff','-j','1','-G','50000','-L','10','-n','2','--outc','0','--outs','0.20','-N','32','--outn','32','-t','1',str(target),str(protein)]
                    raw=cache/'raw_miniprot/1_Hx50'/family/(internal+'.gff3')
                else:
                    cmd=['exonerate-rs','--model','protein2genome','--querytype','protein','--targettype','dna','--exhaustive','no','--subopt','no','--revcomp','yes','--score','1','--dpmemory','64','--verbose','0','--showalignment','no','--showsugar','no','--showcigar','no','--showvulgar','yes','--maxintron','50000',str(protein),str(target)]
                    raw=cache/'raw_exonerate/1_Hx50'/family/(internal+'.vulgar')
                start=time.monotonic(); result=subprocess.run(cmd,capture_output=True,text=True,env=env,timeout=300)
                if result.returncode: raise RuntimeError(f'{family}/{internal}: {result.stderr}')
                raw.parent.mkdir(parents=True,exist_ok=True)
                tmp=raw.with_suffix(raw.suffix+'.tmp');tmp.write_text(result.stdout);tmp.replace(raw)
                hits=sum(line.startswith('vulgar:') if backend!='miniprot' else '\tmRNA\t' in line for line in result.stdout.splitlines())
                return dict(family=family,candidate=r['candidate'],command=cmd,seconds=time.monotonic()-start,raw_hits=hits)
            state['phase']=backend+'-short-alignments';save(out/'status.json',state)
            start=time.monotonic()
            with ThreadPoolExecutor(max_workers=8) as pool: aligned=list(pool.map(align,selected))
            align_seconds=time.monotonic()-start
            save(out/(backend+'-new-alignments.json'),aligned)
            tmp=manifest.with_suffix('.tmp')
            with tmp.open('w') as f, manifest.open() as source:
                writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader()
                for r in csv.DictReader(source,delimiter='\t'):
                    if (r['family_id'],r['candidate']) in selected_keys and r['state']=='protein_screen_rejected':
                        r['state']='read_backed_short_screen_passed';r['qc_flags']='exact_unique_8aa;two_unique_reads;full_contig_coverage'
                    writer.writerow(r)
            tmp.replace(manifest)
            dest=out/backend
            command=[binary,'recover-sco','--annotation',cache,'--input',old/'gene-full','--out',dest,'--threads','8','--sco-backbone-manifest',root/'miniprot-backbones.tsv']
            seconds=run(command,backend+'-full-recovery')
            verify(dest,rows(dest/'manifest/sco_recovery.tsv'))
            exons=score(root/'exon_coordinates.json',dest,cache,root/'miniprot-backbones.tsv')
            gaps=audit(root/'exon_coordinates.json',dest,cache,root/'miniprot-backbones.tsv')
            save(out/(backend+'-exons.json'),exons);save(out/(backend+'-gaps.json'),gaps)
            changes=[]
            for path in sorted((baseline/backend/'sco').glob('*.fasta')):
                before=next(iter(fasta(path).values())); after=next(iter(fasta(dest/'sco'/path.name).values()))
                assert len(before)==len(after)
                if before!=after:
                    changes.append(dict(family=path.stem,gained_nt=sum(x=='N' and y!='N' for x,y in zip(before,after)),lost_nt=sum(x!='N' and y=='N' for x,y in zip(before,after)),substituted_nt=sum(x!='N' and y!='N' and x!=y for x,y in zip(before,after))))
            # Reaudit changed candidate families; unchanged chain inputs and logic
            # are identical to the frozen full-library read audit.
            subset=out/(backend+'-subset'); candidate_dir=subset/'samples/1_Hx50/candidates';candidate_dir.mkdir(parents=True)
            for family in sorted({r['family_id'] for r in selected}):
                (candidate_dir/(family+'.fasta')).symlink_to(old/'gene-full/samples/1_Hx50/candidates'/(family+'.fasta'))
            chain_dest=out/(backend+'-read-chains')
            chain_seconds=run([binary,'recover-sco','--annotation',cache,'--input',subset,'--out',chain_dest,'--threads','8','--sco-backbone-manifest',root/'miniprot-backbones.tsv','--reads-root',root/'workflow'],backend+'-changed-read-chains')
            combined=out/(backend+'-all-chains')
            subprocess.run(['cp','-al',str(baseline/backend/'chain_audit'),str(combined)],check=True)
            for path in (chain_dest/'chain_audit').glob('*/*.tsv'):
                target=combined/path.parent.name/path.name;tmp=target.with_suffix('.tmp');shutil.copy2(path,tmp);tmp.replace(target)
            for path in (chain_dest/'sco').glob('*.fasta'):
                assert path.read_bytes()==(dest/'sco'/path.name).read_bytes()
            chains=summarize(combined)
            report['backends'][backend]=dict(exons=exons['totals'],fraction=exons['complete_exon_interval_fraction'],gap_bins=gaps['bins'],changes=changes,chains=chains,new_alignments=len(aligned),new_raw_hit_candidates=sum(r['raw_hits']>0 for r in aligned),alignment_seconds=align_seconds,recovery_seconds=seconds,changed_family_read_audit_seconds=chain_seconds)
            save(out/'report.json',report)
        state.update(status='completed',phase='completed',finished_unix=time.time())
    except Exception:
        state.update(status='failed',error=traceback.format_exc());raise
    finally: save(out/'status.json',state)

if __name__=='__main__':main()
