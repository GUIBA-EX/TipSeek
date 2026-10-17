#!/usr/bin/env python3
"""Compare population binaries using deterministic synthetic diploid UCE reads.

Requires real minibwa, samtools, bcftools and PLINK executables in --tools.
This fixture checks result equivalence, not biological performance on real data.
"""
import argparse
import random, subprocess, pathlib, time, json, os
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', type=pathlib.Path, required=True)
parser.add_argument('--updated', type=pathlib.Path, required=True)
parser.add_argument('--tools', type=pathlib.Path, required=True)
parser.add_argument('--output', type=pathlib.Path, required=True)
options = parser.parse_args()
root = options.output.resolve()
root.mkdir(parents=True, exist_ok=False)
tools = options.tools.resolve()
env = dict(os.environ, PATH=str(tools) + ':' + os.environ['PATH'])
rng = random.Random(37)
refs = [''.join(rng.choices('ACGT', k=900)) for _ in range(12)]
records = []
complement = str.maketrans('ACGT', 'TGCA')
for sample in range(12):
    paths = [root / f's{sample}_{mate}.fq' for mate in (1, 2)]
    handles = [p.open('w') for p in paths]
    for locus, ref in enumerate(refs):
        for fragment in range(110):
            dna = list(ref)
            for site in (210, 410, 610):
                if sample % 3 == 2 or (sample % 3 == 1 and fragment % 2):
                    dna[site] = next((b for b in 'ACGT' if b != dna[site]))
            start = rng.randrange(0, 600)
            frag = ''.join(dna[start:start + 300])
            for mate, seq in enumerate((frag[:150], frag[-150:].translate(complement)[::-1])):
                handles[mate].write(f'@s{sample}_l{locus}_p{fragment}/{mate + 1}\n{seq}\n+\n' + 'I' * 150 + '\n')
    for h in handles:
        h.close()
    records.append(f'S{sample}\t{paths[0]}\t{paths[1]}\n')
manifest = root / 'samples.tsv'
manifest.write_text(''.join(records))
header = 'locus,status,accepted,selected_contig_length,slice_supported_bases,slice_support_breadth,max_slice_support_gap,unique_read_count,unique_read_density,low_quality\n'

def setup(name):
    out = root / name
    for sample in range(12):
        d = out / f'{sample + 1}_S{sample}'
        (d / 'results').mkdir(parents=True, exist_ok=True)
        (d / 'uce_assembly_summary.csv').write_text(header + ''.join((f'uce{l},success,true,900,900,1,0,110,0.12,false\n' for l in range(12))))
        for l, seq in enumerate(refs):
            (d / 'results' / f'uce{l}.fasta').write_text(f'>sample\n{seq}\n')
    return out
results = {}

def run(name, binary, threads=4, extra=(), setup_new=True):
    out = setup(name) if setup_new else root / name
    cmd = [str(binary), '--samples', str(manifest), '--output', str(out), '--threads', str(threads), '--skip-admixture'] + list(extra)
    t = time.monotonic()
    with (root / f'{name}-{len(results)}.log').open('w') as log:
        proc = subprocess.run(cmd, stdout=log, stderr=log, env=env)
    elapsed = time.monotonic() - t
    if proc.returncode:
        raise RuntimeError((name, proc.returncode))
    return elapsed
baseline = options.baseline.resolve()
new = options.updated.resolve()
results['baseline_seconds'] = run('baseline', baseline)
results['updated_seconds'] = run('updated', new)
results['one_cpu_seconds'] = run('one_cpu', new, 1)
results['supported_seconds'] = run('supported', new, 4, ['--reference-strategy', 'supported'])

def variants(name):
    p = root / name / 'population/variants/cohort.filtered.vcf.gz'
    return subprocess.check_output([str(tools / 'bcftools'), 'view', '-H', str(p)], env=env)

def calls(name):
    p = root / name / 'population/variants/cohort.filtered.vcf.gz'
    return subprocess.check_output([str(tools / 'bcftools'), 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n', str(p)], env=env)
base = calls('baseline')
assert base, 'empty validation is not informative'
for name in ('updated', 'one_cpu', 'supported'):
    assert calls(name) == base, name
results['called_snps'] = len(base.splitlines())
results['genotypes_identical'] = True
def alignments(run, sample):
    bam = root / run / f'population/mapping/S{sample}.bam'
    return sorted(subprocess.check_output([str(tools / 'samtools'), 'view', str(bam)], env=env).splitlines())

for name in ('updated', 'one_cpu'):
    for sample in range(12):
        assert alignments(name, sample) == alignments('baseline', sample), (name, sample)
results['alignment_records_identical'] = True
# Exercise the final sort writer as well as markdup, including one CPU.
results['no_duplicates_baseline_seconds'] = run(
    'no_duplicates_baseline', baseline, 4, ['--skip-mark-duplicates'])
for name, threads in [('no_duplicates', 4), ('no_duplicates_one_cpu', 1)]:
    results[name + '_seconds'] = run(name, new, threads, ['--skip-mark-duplicates'])
    assert variants(name) == variants('no_duplicates_baseline'), name
    for sample in range(12):
        assert alignments(name, sample) == alignments('no_duplicates_baseline', sample), (name, sample)
# Confirm the BAI is usable for random access, not merely present on disk.
for name in ('updated', 'one_cpu', 'no_duplicates', 'no_duplicates_one_cpu'):
    for sample in range(12):
        bam = root / name / f'population/mapping/S{sample}.bam'
        assert pathlib.Path(str(bam) + '.bai').is_file(), bam
        indexed = subprocess.check_output(
            [str(tools / 'samtools'), 'view', str(bam), 'uce0'], env=env)
        sequential = subprocess.check_output(
            [str(tools / 'samtools'), 'view', str(bam)], env=env)
        expected = [row for row in sequential.splitlines() if row.split(b'\t')[2] == b'uce0']
        assert indexed.splitlines() == expected and expected, bam
results['bai_random_access_verified'] = True
# Compare the complete report, including zero-coverage references and fields
# not currently used by mapping_qc.tsv.
for name in ('updated', 'one_cpu', 'no_duplicates', 'no_duplicates_one_cpu'):
    for sample in range(12):
        bam = root / name / f'population/mapping/S{sample}.bam'
        report = bam.with_suffix('.coverage.tsv')
        if report.exists():
            expected = subprocess.check_output([str(tools / 'samtools'), 'coverage', str(bam)], env=env)
            assert report.read_bytes() == expected, bam
# Exercise flag filtering and CIGAR semantics using real samtools.
edge_sam = root / 'coverage_edge.sam'
rows = ['@HD\tVN:1.6\tSO:unsorted', '@SQ\tSN:chr1\tLN:100', '@SQ\tSN:zero\tLN:100']
for i, (flag, cigar) in enumerate([(0, '20M'), (16, '5S15M'), (0, '10M2I8M'),
                                 (0, '10M3D10M'), (0, '10M20N10M'), (256, '20M'),
                                 (512, '20M'), (1024, '20M'), (2048, '20M'), (4, '*')]):
    rows.append(f'r{i}\t{flag}\tchr1\t10\t60\t{cigar}\t*\t0\t0\t' + 'A' * 20 + '\t' + 'I' * 20)
rows.append('unplaced\t4\t*\t0\t0\t*\t*\t0\t0\t' + 'A' * 20 + '\t' + 'I' * 20)
edge_sam.write_text('\n'.join(rows) + '\n')
edge_bam = root / 'coverage_edge.bam'
subprocess.run([str(tools / 'samtools'), 'sort', '-o', str(edge_bam), str(edge_sam)], check=True, env=env)
subprocess.run([str(tools / 'samtools'), 'index', str(edge_bam)], check=True, env=env)
bed = root / 'coverage_edge.bed'
bed.write_text('chr1\t0\t100\nzero\t0\t100\n')
view = subprocess.Popen([str(tools / 'samtools'), 'view', '-u', '-M', '-L', str(bed), str(edge_bam)],
                        stdout=subprocess.PIPE, env=env)
coverage = subprocess.run([str(tools / 'samtools'), 'coverage', '-'], stdin=view.stdout,
                          stdout=subprocess.PIPE, env=env)
view.stdout.close()
assert view.wait() == 0 and coverage.returncode == 0
assert coverage.stdout == subprocess.check_output([str(tools / 'samtools'), 'coverage', str(edge_bam)], env=env)
results['coverage_reports_identical'] = True
for name in ('updated', 'one_cpu', 'supported'):
    assert (root / name / 'population/structure/selected_snps.tsv').read_bytes() == (root / 'baseline/population/structure/selected_snps.tsv').read_bytes()
raw = root / 'updated/population/variants/cohort.raw.bcf'
stamp = raw.stat().st_mtime_ns
results['refilter_seconds'] = run('updated', new, 4, ['--start-at', 'calling', '--min-dp', '6'], False)
assert raw.stat().st_mtime_ns == stamp
results['raw_calls_reused'] = True
results['selection_seconds'] = run('updated', new, 4, ['--start-at', 'selection', '--min-dp', '6'], False)
(root / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
