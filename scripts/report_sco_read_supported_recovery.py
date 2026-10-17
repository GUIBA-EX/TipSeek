#!/usr/bin/env python3
"""Summarize HX50 rescue and reads auditing without treating reference order as truth.
Usage: report_sco_read_supported_recovery.py PANEL RUN REPORT.json
"""
import csv, hashlib, json, shutil, sys
from collections import Counter
from pathlib import Path
from report_starfish_sco_recovery import summarize
from benchmark_exon_shadow import digest


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def output_hashes(root):
    result = {}
    for name in ['models','cds','proteins','exons','introns','genes','genes_flanked','supercontigs','gff3','fragment_cds','unresolved']:
        for path in sorted((root/name).rglob('*')):
            if path.is_file(): result[str(path.relative_to(root))] = digest(path)
    for name in ['candidate_manifest.tsv','fragment_chains.tsv','fragment_groups.tsv','read_audit.tsv','fragment_read_support.tsv','completion_audit.tsv','join_sources.tsv']:
        path = root/'manifest'/name
        # read_audit includes the workflow-root path; all other fields must agree.
        if name == 'read_audit.tsv':
            records = rows(path)
            for record in records: record.pop('inputs')
            result[str(path.relative_to(root))] = hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
        else: result[str(path.relative_to(root))] = digest(path)
    return result


def main():
    panel_root, run, destination = map(Path, sys.argv[1:])
    panel = json.loads((panel_root/'panel.json').read_text())
    report = {'scope':'HX50 real-read recovery; unchanged recruitment/refiltering reused. ppect is another species, not same-individual truth. Coverage is a union of coding evidence, not a complete CDS or validated copy. No biological false-merge rate is measurable.', 'runs':json.loads((run/'runs.json').read_text()),'panel':panel,'results':{}}
    for label, root, backend in [
        ('previous_miniprot',panel_root/'verified-miniprot','miniprot'),
        ('previous_exonerate',panel_root/'verified-exonerate','exonerate-rs'),
        ('current_miniprot',run/'miniprot','miniprot'),
        ('current_exonerate',run/'exonerate-rs','exonerate-rs'),
    ]:
        report['results'][label] = summarize(root,panel,backend)
        if label.startswith('current'):
            audits = rows(root/'manifest/read_audit.tsv')
            support = rows(root/'manifest/fragment_read_support.tsv')
            manifest = rows(root/'manifest/candidate_manifest.tsv')
            report['results'][label]['read_audit'] = audits
            report['results'][label]['read_totals'] = {name:sum(int(r[name]) for r in audits) for name in ['total_pairs','unique_mapped_pairs','distinct_cross_contig_pairs']}
            report['results'][label]['edge_status_records'] = dict(Counter(r['status'] for r in support))
            report['results'][label]['candidate_state_records'] = dict(Counter(r['state'] for r in manifest))
            report['results'][label]['completion_audit'] = rows(root/'manifest/completion_audit.tsv')
            report['results'][label]['padded_attempts'] = rows(root/'manifest/fragment_groups.tsv')
    serial = Path('/tmp/tipseek-hx50-recovery-final')
    a, b = output_hashes(serial/'miniprot'), output_hashes(run/'miniprot')
    report['serial_parallel_miniprot'] = {'serial_runs':json.loads((serial/'runs.json').read_text()),'all_compared_outputs_equal':a==b,'compared_files':len(a),'different_files':[name for name in sorted(set(a)|set(b)) if a.get(name)!=b.get(name)],'scope':'8 candidate workers, each single-threaded, versus one worker. Read-root path excluded from audit comparison. Commands/execution metadata intentionally differ.'}
    report['direct_refiltered_read_diagnostic'] = json.loads(Path('/tmp/tipseek-hx50-direct-reads/report.json').read_text())
    sources = ['rust/gm2_tools/src/bin/gene_workflow.rs','rust/gm2_tools/src/gene_annotation/chain.rs','rust/gm2_tools/src/gene_annotation/reads.rs','rust/gm2_tools/src/gene_annotation/screen.rs','rust/main_assembler_original/src/main.rs','rust/tipseek_cli/src/main.rs']
    report['source_sha256'] = {p:digest(Path(p)) for p in sources}
    report['reference_logic_sha256'] = {p:digest(Path(p)) for p in ['Captus/captus/bioformats.py','HybPiper/hybpiper/exonerate_hits.py']}
    exe = shutil.which('exonerate-rs')
    report['exonerate_executable'] = {'path':exe,'sha256':digest(Path(exe))} if exe else None
    report['incomplete_experiments_excluded'] = ['Unscreened all-seed Exonerate','screened greedy orphan Exonerate','unitig preliminary Exonerate','serial final Exonerate: cancelled in favor of final parallel run; no completed timing claimed']
    destination.write_text(json.dumps(report,indent=2)+'\n')
    for label,result in report['results'].items():
        print(label, result['reference_exon_recovery_fraction'],result['mean_cds_union_coverage'],result['complete_models'])
    print('serial_parallel_equal',a==b)


if __name__ == '__main__': main()
