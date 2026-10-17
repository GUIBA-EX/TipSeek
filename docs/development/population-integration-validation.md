# Population integration: review and real-read validation

Validated on 2026-10-18 against main-branch commit `d83dea5` (the released population implementation). This check measures computational equivalence and runtime; it is not a biological accuracy or population-structure validation.

## Changes reviewed

- Keep accepted UCE contigs as the default `pseudoref` input and preserve the `sqcl-longest` default.
- Add stage dependency/completion receipts, reference content fingerprints, read-file identities, external-tool fingerprints and raw joint-call reuse. Changed reference inputs must fail before an existing reference is overwritten during restart.
- Count concurrent mapper/sort threads together, constrain sorting memory, and pass the requested thread count to PLINK. Single-thread pipelines spool sequentially.
- Use uncompressed temporary BAMs and pipe fixmate directly into coordinate sorting when the thread budget permits. Final BAMs remain compressed.
- Correct empty R2 fields, collisions in normalized sample IDs, reference copies onto the same inode (including hard links), non-executable PATH candidates, and cleanup of the producer when the pipe consumer cannot start.

## Data and procedure

Real cleaned paired-end UCE reads: `SRR10443675`, `SRR10443764`, `SRR10443767`, `SRR10443774`. Take the first 500,000 pairs per library: 2,000,000 pairs total, with original sequences and qualities, without simulation. Select the first 128 sorted loci shared by existing assemblies and available in the bait directory. This intentionally bounded subset is not a full-library or full-panel benchmark.

Run the current UCE workflow first:

```bash
cli/tipseek filter assemble --assembly-mode uce \
  -f samples.tsv -r references -o uce -p 8
```

The four UCE outputs accepted 41, 38, 30 and 29 loci respectively. Population selected 61 reference loci. Both implementations consume the same accepted contigs, UCE support summaries and original read subset. No synthetic assembly summaries are substituted.

On the same Linux/NVMe host, run each implementation twice with eight requested threads, default filtering and duplicate marking, PLINK enabled and ADMIXTURE disabled. Alternate order baseline–updated–updated–baseline without flushing OS caches. Timings include reference construction, mapping, calling, filtering, three panels, PCA and individual QC; they exclude input preparation, UCE assembly and comparison work. Toolchain: samtools 1.16.1, bcftools 1.16, PLINK 1.9.0-b.7.7; minibwa is identified by its executable SHA-256 in the local provenance record.

## Final-build results

| Implementation | Run 1 | Run 2 | Median |
| --- | ---: | ---: | ---: |
| Released population implementation | 14.38 s | 14.23 s | 14.30 s |
| Reviewed implementation | 10.17 s | 10.04 s | 10.11 s |

Elapsed time decreased by **29.4%** on this subset. An earlier paired trial of the preceding cache/scheduler implementation versus the piped implementation measured medians of 19.13 s and 10.08 s respectively. These measurements do not establish a universal speedup; the tiny synthetic fixture remains dominated by startup and validation overhead.

All four final-build runs agree on:

- The complete reference FASTA.
- All 4,000,012 SAM alignment records, compared as sorted records without BAM headers.
- Every VCF/BCF record across raw calls, normalized SNPs, genotype filtering, tag filling, site filtering and the three panels, including quality, genotype and depth fields. Header command provenance may differ.
- 263 final filtered SNPs and 9 one-SNP-per-UCE loci.
- Selected-SNP tables, all PLINK BED/BIM/FAM panels, individual missingness, heterozygosity and relatedness outputs.
- PCA eigenvalues and eigenvectors within absolute tolerance 1e-8, allowing eigenvector sign reversal and verifying sample order.

A real-data `--start-at calling` replay completed in 0.19 s, left the raw BCF inode/size/mtime unchanged, and reproduced the filtered variant records. Separate synthetic and regression tests cover changed filter thresholds, MAPQ changes, stale references/reads, reordered samples, absent receipts, and single-thread execution.

The speedup trades CPU compression work for disk traffic: GNU time reported 2,872,984 filesystem output blocks for the baseline and 5,322,120 for the updated implementation (about 1.85×). Temporary BAMs are removed after success. This is cumulative reported I/O, not peak temporary disk occupancy. Peak per-process child RSS was about 465 MiB versus 381 MiB; these figures are not summed concurrent memory consumption. Slower storage requires separate benchmarking.

## Reproduction and checks

```bash
python3 scripts/benchmark_population_real.py \
  --baseline /path/to/released/main_population \
  --updated target/release/main_population \
  --samples /path/to/absolute-path-samples.tsv \
  --uce-output /path/to/uce \
  --tools /path/to/tool/bin \
  --output /path/to/new-benchmark-directory \
  --threads 8 --repeats 2
```

The script creates isolated output directories, records executable hashes, preserves logs and timings, and fails on result differences. `--previous` optionally adds a third implementation. The companion `scripts/validate_population_integration.py` provides a deterministic synthetic fixture and verifies raw-call reuse and single-thread equivalence.

All 365 workspace tests, Clippy with warnings denied, Rust formatting, and Python script compilation passed. ADMIXTURE runtime/equivalence, full datasets, non-Linux resource behavior and biological truth were not tested here.

Final tested executable SHA-256: `a0d3c4fbf67a3bb31d44d8331bc7a8bf286b224e5401b7906db8da33ff05c557`.
