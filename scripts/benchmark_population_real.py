#!/usr/bin/env python3
"""Compare population implementations on existing UCE assemblies and real reads.

Creates isolated output directories and leaves input assemblies untouched. Reads
must be identified by absolute paths in the sample manifest. ADMIXTURE is skipped;
all variant stages, three panels, PCA, individual QC and SAM records are checked.
Requires GNU time, sort, and minibwa/samtools/bcftools/PLINK under --tools.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import time


def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def records(tool, path):
    return subprocess.check_output([str(tool), 'view', '-H', str(path)])


def alignments(tool, path):
    producer = subprocess.Popen([str(tool), 'view', str(path)], stdout=subprocess.PIPE)
    sorter = subprocess.Popen(['sort', '-S', '128M'], stdin=producer.stdout,
                              stdout=subprocess.PIPE, env=dict(os.environ, LC_ALL='C'))
    producer.stdout.close()
    digest = hashlib.sha256()
    count = 0
    for chunk in iter(lambda: sorter.stdout.read(1024 * 1024), b''):
        digest.update(chunk)
        count += chunk.count(b'\n')
    check(sorter.wait() == 0 and producer.wait() == 0, f'SAM comparison failed: {path}')
    return digest.hexdigest(), count


def compare_pca(reference, current):
    for panel in ['all_snps_pca', 'population_pca', 'ld_pruned_pca']:
        for extension in ['eigenval', 'eigenvec']:
            name = f'structure/{panel}.{extension}'
            left = [line.split() for line in (reference / name).read_text().splitlines()]
            right = [line.split() for line in (current / name).read_text().splitlines()]
            check(len(left) == len(right), f'PCA row count differs: {name}')
            if extension == 'eigenvec':
                check([r[:2] for r in left] == [r[:2] for r in right], 'PCA sample order differs')
                left, right = [r[2:] for r in left], [r[2:] for r in right]
            check(left and all(len(a) == len(b) == len(left[0]) for a, b in zip(left, right)),
                  f'PCA dimensions differ: {name}')
            left, right = [[float(v) for v in r] for r in left], [[float(v) for v in r] for r in right]
            check(all(math.isfinite(v) for rows in (left, right) for r in rows for v in r), 'Non-finite PCA output')
            for col in range(len(left[0])):
                plus = max(abs(a[col] - b[col]) for a, b in zip(left, right))
                minus = max(abs(a[col] + b[col]) for a, b in zip(left, right))
                check((min(plus, minus) if extension == 'eigenvec' else plus) < 1e-8,
                      f'PCA differs beyond sign/tolerance: {name}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['baseline', 'updated', 'samples', 'uce-output', 'tools', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--previous', type=Path)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--repeats', type=int, default=2)
    args = parser.parse_args()
    check(args.threads > 0 and args.repeats > 0, 'threads/repeats must be positive')
    root, tools = args.output.resolve(), args.tools.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binaries = {'baseline': args.baseline.resolve(), 'updated': args.updated.resolve()}
    if args.previous:
        binaries['previous'] = args.previous.resolve()
    sources = sorted({p.parent for p in args.uce_output.resolve().glob('*/uce_assembly_summary.csv')})
    check(sources, 'No current UCE assembly summaries found')
    env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ.get('PATH', ''))
    summary = {'binaries': {name: {'path': str(path), 'sha256': sha256(path)} for name, path in binaries.items()},
               'samples': str(args.samples.resolve()), 'samples_sha256': sha256(args.samples),
               'threads': args.threads, 'runs': []}
    for repeat in range(args.repeats):
        order = list(binaries)
        if repeat % 2:
            order.reverse()
        for version in order:
            out = root / f'{version}-{repeat + 1}'
            out.mkdir()
            for source in sources:
                (out / source.name).symlink_to(source, target_is_directory=True)
            cmd = [str(binaries[version]), '--samples', str(args.samples.resolve()), '--output', str(out),
                   '--threads', str(args.threads), '--skip-admixture']
            started = time.monotonic()
            log_path = root / f'{out.name}.log'
            with log_path.open('w') as log:
                result = subprocess.run(['/usr/bin/time', '-v'] + cmd, stdout=log, stderr=log, env=env)
            run = {'version': version, 'repeat': repeat + 1, 'seconds': time.monotonic() - started,
                   'returncode': result.returncode, 'output': str(out)}
            log = log_path.read_text()
            for field, pattern in [('max_child_rss_kib', r'Maximum resident set size \(kbytes\):\s*(\d+)'),
                                   ('filesystem_output_blocks', r'File system outputs:\s*(\d+)')]:
                match = re.search(pattern, log)
                if match:
                    run[field] = int(match.group(1))
            summary['runs'].append(run)
            (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
            print(json.dumps(run), flush=True)
            check(result.returncode == 0, f'{out.name} failed; see {log_path}')
    reference = root / 'baseline-1/population'
    variant_paths = ['variants/cohort.raw.bcf', 'variants/cohort.biallelic.snps.vcf.gz',
                     'variants/cohort.genotype_filtered.vcf.gz', 'variants/cohort.tagged.vcf.gz',
                     'variants/cohort.filtered.vcf.gz', 'structure/all_snps.vcf.gz',
                     'structure/one_snp_per_uce.vcf.gz', 'structure/ld_pruned.vcf.gz']
    expected_vcf = {name: records(tools / 'bcftools', reference / name) for name in variant_paths}
    expected_bam = {p.name: alignments(tools / 'samtools', p) for p in sorted((reference / 'mapping').glob('*.bam'))}
    check(expected_bam and expected_vcf['variants/cohort.filtered.vcf.gz'], 'Empty validation outputs')
    exact_files = ['reference/population_reference.fasta', 'structure/selected_snps.tsv']
    exact_files += [f'structure/{prefix}.{ext}' for prefix in ['all_snps', 'population', 'ld_pruned']
                    for ext in ['bed', 'bim', 'fam']]
    exact_files += [f'structure/qc/individuals.{ext}' for ext in ['imiss', 'het', 'genome']]
    for run in summary['runs']:
        current = Path(run['output']) / 'population'
        for name in exact_files:
            check((reference / name).read_bytes() == (current / name).read_bytes(), f'{current}: {name} differs')
        for name, expected in expected_vcf.items():
            check(records(tools / 'bcftools', current / name) == expected, f'{current}: {name} differs')
        compare_pca(reference, current)
        for name, expected in expected_bam.items():
            check(alignments(tools / 'samtools', current / 'mapping' / name) == expected, f'{current}: {name} differs')
        run['equivalent'] = True
    summary['reference_loci'] = (reference / exact_files[0]).read_text().count('>')
    summary['filtered_snps'] = len(expected_vcf['variants/cohort.filtered.vcf.gz'].splitlines())
    summary['one_snp_loci'] = len(expected_vcf['structure/one_snp_per_uce.vcf.gz'].splitlines())
    summary['alignment_records'] = sum(count for _, count in expected_bam.values())
    summary['median_seconds'] = {v: statistics.median(r['seconds'] for r in summary['runs'] if r['version'] == v)
                                 for v in binaries}
    (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ['runs', 'binaries']}, indent=2))


if __name__ == '__main__':
    main()
