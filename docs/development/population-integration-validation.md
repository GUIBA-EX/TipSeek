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


## Follow-up: remove the coordinate-sorted intermediate BAM

Compared against `23ff0fa` on the same read subset, accepted UCE outputs, host and tools, with eight threads. Three fresh runs per binary alternated baseline–updated, updated–baseline, baseline–updated; OS caches were not flushed. The timing boundary and equivalence checks are unchanged.

The retained change streams `fixmate -u` → coordinate `sort -u` → `markdup`. Name-sorted BAMs and compressed final BAMs remain; the complete intermediate coordinate-sorted BAM is no longer written/read. Sort spills, if needed, use a sample-specific prefix inside the mapping directory. A one-thread budget retains the previous sequential path.

Coordinate sorting consumes all input before emitting alignments, so fixmate and markdup do not process alignment records concurrently. Existing sample concurrency and sort memory/thread limits are preserved; the sort phase needs its threads plus one active producer or consumer. All pipeline exit statuses are checked, and failures terminate and reap remaining children.

| Implementation | Run 1 | Run 2 | Run 3 | Median |
| --- | ---: | ---: | ---: | ---: |
| `23ff0fa` | 10.083 s | 10.193 s | 10.240 s | 10.193 s |
| Streaming duplicate marking | 9.649 s | 9.944 s | 9.738 s | 9.738 s |

Median elapsed time decreased **4.46%**; every paired run was faster. GNU time filesystem output blocks fell from **5,322,120 to 2,997,624 (43.68%)** in every run. Maximum child RSS remained approximately 381 MiB. These are cumulative I/O and maximum individual-child RSS, not peak disk usage or aggregate memory. The modest runtime gain applies to this bounded subset and host, not all datasets or storage systems.

Two more elaborate scheduling candidates were rejected: reducing concurrent samples gave medians of 9.958 → 13.977 s; phase-based CPU permits gave 10.139 → 13.537 s. Both preserved results and reduced I/O, but failed the speed requirement. Neither scheduler is included.

All six retained-candidate runs preserve the reference, all **4,000,012 alignment records**, all eight variant record sets, **263 filtered SNPs**, **9 one-SNP-per-UCE loci**, PLINK panels, individual QC and PCA under the checks described above. ADMIXTURE is still excluded. Regression checks additionally exercise every pipeline stage failing, successful streaming, missing executables and the existing single-thread fallback.

Local evidence directory: `tipseek_population_real_validation_20261018/bam-stream-trial` (summary, per-run logs and full outputs). Reproduce with the command above using `23ff0fa` as baseline and `--repeats 3`.

All 366 workspace tests, Clippy with warnings denied and formatting checks passed. The real-tool synthetic fixture also preserved alignments/genotypes, verified one-thread execution and reused raw calls when refiltering. Rebuilding after the added regression test reproduced the benchmarked executable hash.

Tested executable SHA-256: `5530833cc6bee01861b7bf59e2d7c318f0506daf3dc515ed68e15dcd4ab31972`.


## Follow-up: build BAI while writing final BAM

Compared against `023cb41` using the identical real UCE subset, accepted assemblies, eight threads and toolchain above. The remaining redundant pass was standalone `samtools index`: after creating the final BAM, every sample reopened and decoded that complete BAM to build its index. Both final writers (`markdup`, or `sort` when duplicate marking is skipped) now use `--write-index` with an explicit `.bam.bai` destination. Explicit BAI naming preserves the existing restart contract instead of using samtools' default CSI output. The BAM quickcheck, flagstat and coverage checks remain.

Three fresh runs per implementation, in alternating paired order, measured:

| Implementation | Run 1 | Run 2 | Run 3 | Median |
| --- | ---: | ---: | ---: | ---: |
| `023cb41` | 9.696 s | 9.742 s | 9.736 s | 9.736 s |
| Index during final BAM writing | 9.472 s | 9.553 s | 9.431 s | 9.472 s |

Median elapsed time decreased **2.72%**, and each paired run improved. Filesystem output blocks remained **2,997,624**; maximum child RSS remained approximately 381 MiB. The optimization removes a decoding/read pass, not another intermediate file; cached reads are not necessarily physical disk reads. This modest gain is specific to this subset and host.

All six runs preserve all **4,000,012 alignment records**, all eight variant record sets, **263 filtered SNPs**, **9 one-SNP-per-UCE loci**, reference FASTA, PLINK panels, individual QC and PCA. The benchmark now also verifies BAI presence, index counts, complete flagstat reports and mapping QC (excluding output-directory paths).

The real-tool synthetic regression additionally compares the default and skip-duplicate-marking paths with the baseline, verifies one-thread execution in both paths, checks indexed region queries against sequential BAM records, and exercises raw-call reuse and selection restart. Synthetic timings are correctness checks, not speed claims. All 366 workspace tests, Clippy with warnings denied, formatting and Python compilation passed.

Local evidence: `tipseek_population_real_validation_20261018/bam-index-trial` and `bam-index-synthetic`. Reproduce with the same command using `023cb41` as baseline and `--repeats 3`. ADMIXTURE and full-library performance remain outside this validation.

Tested executable SHA-256: `934fe1d1368b4e4e6577f36efb2ca5cae140fc15cb2a312981fcc19f651a3780`.


## Follow-up: indexed coverage QC

Compared against `2d6492e` with the same real UCE reads, accepted assemblies, eight threads and toolchain. The four final BAMs contain 4,000,012 records, including **3,991,220 unplaced unmapped records (99.78%)**. The previous coverage command decoded the entire BAM even though those unplaced records cannot contribute to coverage.

Coverage now uses `samtools view -u -M -L` with full-reference BED intervals and the existing BAI, streaming the selected records into the original `samtools coverage`. It retains the full BAM header, coverage filters, depth limit, rounding and zero-coverage references. Flagstat continues to scan all records. The two single-thread commands fit the existing per-worker budget; one-CPU runs retain sequential full-scan coverage. BED-reserved reference names, empty indexes and zero-length references also use the original path. Invalid FAI lengths fail explicitly. Complete coverage reports are saved for both paths.

The initial candidate measured medians of 9.486 → 9.213 s with all results equal. After adding reference-name fallback and consistent report output, the **final build** measured:

| Implementation | Run 1 | Run 2 | Run 3 | Median |
| --- | ---: | ---: | ---: | ---: |
| `2d6492e` | 9.417 s | 9.480 s | 9.611 s | 9.480 s |
| Indexed coverage QC | 9.256 s | 9.263 s | 9.207 s | 9.256 s |

Median end-to-end elapsed time decreased **2.36%**, with every paired run faster. Filesystem output blocks changed from 2,997,624 to 2,997,664 because the small BED and coverage reports are now retained. Maximum child RSS stayed approximately 381 MiB. This improvement is specific to this dataset and host; a dataset dominated by mapped reads may benefit less or incur extra pipe overhead. Cached BAM reads are not necessarily physical disk reads.

All six final-build runs preserve all alignment and variant records, 263 filtered SNPs, 9 one-SNP-per-UCE loci, reference FASTA, PLINK panels, individual QC and PCA. Comparisons additionally cover complete coverage reports (all fields), mapping QC, flagstat reports and BAI counts. The real-tool regression checks clipped, inserted, deleted and spliced CIGAR operations, reverse/secondary/supplementary/duplicate/QC-failed records, placed/unplaced unmapped reads and a zero-coverage reference. Single-thread, skip-duplicate-marking, BAI random access and restart checks also pass.

All **367** workspace tests passed on a serial rerun, along with Clippy, formatting and Python compilation. The first parallel workspace run encountered an unrelated `gene_workflow` fixture `ExecutableFileBusy` error; it did not recur with serial execution. The population module tests passed separately.

An isolated coverage-only measurement ran the four existing BAMs sequentially per trial, alternating full-scan and indexed modes over three trials and checking report bytes. Median aggregate coverage time was **0.9930 → 0.00894 s (99.1% less)**. This isolates coverage work; it is not the wall-time saving of the concurrently scheduled full workflow. Evidence and executable procedure are in `coverage_qc_profile.json` and `profile_coverage_qc.py` beside the local benchmark directories.

Local evidence: `tipseek_population_real_validation_20261018/bam-qc-final` and `bam-qc-synthetic`. Reproduce using the benchmark command above with `2d6492e` as baseline and `--repeats 3`. ADMIXTURE and full-library performance remain untested.

Final tested executable SHA-256: `fb8d38b9da87cf4824e93bf2fd51abf81dae1e05f5dc8a181ebcfa744b320f02`.


## Follow-up: reduce final BAM compression work

Compared against `9778103` using the same real UCE subset, assemblies, eight threads and external tools. The retained change sets lossless BAM compression to **level 1** in both final writers: duplicate marking and coordinate sorting when duplicate marking is skipped. Index creation, record ordering, filtering, duplicate decisions, full flagstat scanning and coverage semantics are unchanged. Intermediate BAM handling and thread/memory budgets are unchanged.

Flagstat decompression threading was also measured independently on the four existing BAMs: three alternating trials gave aggregate sequential-sample medians of 0.9905 s with zero additional threads and 0.9627 s with one additional thread. One paired trial was slower with the extra thread (1.0044 versus 0.9902 s). This small, inconsistent improvement was not retained. Local measurements are in `flagstat_threads_profile.json`.

Three fresh alternating paired runs of the **retained final build** measured:

| Implementation | Run 1 | Run 2 | Run 3 | Median |
| --- | ---: | ---: | ---: | ---: |
| `9778103` | 9.181 s | 9.212 s | 9.203 s | 9.203 s |
| Final BAM compression level 1 | 7.079 s | 7.108 s | 7.155 s | 7.108 s |

Median end-to-end time decreased **22.77%**, and every paired run improved. The four final BAMs grew from **371,553,820 to 391,165,669 bytes (+5.28%)**, identically across all three runs. Filesystem output blocks increased from **2,997,664 to 3,035,968 (+1.28%)**. Maximum child RSS remained approximately 381 MiB. The storage tradeoff is intentional and measured; this host/subset result does not establish a universal gain on slower storage or other data.

All six runs preserve all **4,000,012 alignment records**, all eight variant record sets, **263 filtered SNPs**, **9 one-SNP-per-UCE loci**, the reference, PLINK panels, individual QC and PCA. Complete flagstat, coverage and mapping QC reports and BAI counts also agree. The real-tool regression covers both final writers, single-thread execution, indexed region queries, coverage edge cases and restart/raw-call reuse.

All **367 workspace tests** passed with serial test execution, along with Clippy with warnings denied and formatting checks. The implementation change is limited to the three output compression options and the tool fixture recognizing those options; no alternative scheduler or flagstat algorithm was introduced.

Local evidence: `tipseek_population_real_validation_20261018/bam-fastcompression-final` (including `summary.json` and `bam_sizes.json`) and `bam-fastcompression-synthetic`. Reproduce with the benchmark command above using `9778103` as baseline and `--repeats 3`. Timings exclude UCE assembly and comparison work, as before; ADMIXTURE remains excluded.

Final tested executable SHA-256: `bc2c7d932f01be78101bc85f536ffb67baa1c17b99b293181818dd847448bca7`.
