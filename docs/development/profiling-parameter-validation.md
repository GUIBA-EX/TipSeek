# Profiling parameter simplification validation

Validated on 2026-10-18 using NC7_3 (SRR36190242), NC2_12 (SRR36190314), and NC1_14 (SRR36190404). Both host and symbiont modes start from the same complete paired raw WGS files, use the existing GenBank/Symb-SHIN reference panels, group maps and opposite-mode decoys, and run with `-p 1`. Independent amplicon controls are not inputs to assembly.

The baseline is the frozen, previously validated ITS2 implementation before CLI simplification, not the v1.6.9 release. The computational component binaries (recruiter, assembler and marker profiler) are byte-identical between versions; only the dispatcher changes. Each run uses a fresh output directory and fresh default index cache. Version order alternates across comparisons. These are single paired elapsed-time measurements, not repeated performance estimates.

| Sample | Mode | Before (s) | After (s) | Biological outputs |
| --- | --- | ---: | ---: | --- |
| NC7_3 | host | 29.41 | 29.33 | Byte-identical |
| NC7_3 | symbiont | 27.33 | 27.32 | Byte-identical |
| NC2_12 | host | 16.71 | 16.59 | Byte-identical |
| NC2_12 | symbiont | 34.44 | 34.33 | Byte-identical |
| NC1_14 | host | 17.46 | 17.42 | Byte-identical |
| NC1_14 | symbiont | 35.57 | 36.57 | Byte-identical |

Total elapsed time: 160.92 → 161.56 seconds (+0.40%). CLI simplification does not change the algorithms; these timings do not establish a speedup.

All top-level biological FASTA/TSV outputs, recruited FASTQ, and GFA graphs match byte-for-byte. QC values also match. Candidate sequences, acceptance/rejection states, phasing diagnostics and abundance estimates are unchanged. Existing low-frequency phasing limitations therefore remain.

The review fixed an output-directory-dependent index cache key: it now uses the original reference rather than its materialized output copy. An additional assembly-free test with recruited real reads verifies that different output directories reuse the same shared index without modifying it and produce identical reference support. Changing reference content at the same path creates a new index.

155 relevant tests and strict Clippy checks pass. The tests cover 600 full-matrix/compact-trace equivalence cases, sequence/quality pairing, conservative phasing and competing models, source attribution, and rejected removed options. The public profiling-specific options are copy mode, group map and decoy; the common cache option remains available. The removed tuning/tool/cache flags are errors, and Themisto is resolved through PATH.

## Reproduction

Build the four release binaries from each revision into separate component directories. For each sample and mode, use separate output directories and run:

```bash
/usr/bin/time -f "elapsed_s=%e user_s=%U system_s=%S max_rss_kib=%M" \
  /path/to/version/tipseek profiling -f sample.tsv \
  -r references/MODE_its2.fasta -o fresh_output \
  --profile-copy-mode MODE --profile-group-map references/MODE_groups.tsv \
  --profile-decoy references/OTHER_MODE_its2.fasta -p 1
```

Keep the same Themisto environment and reference files for both versions. Unset GM2_COMPONENT_DIR so each dispatcher uses its adjacent component binaries. Compare result bytes, not logs or elapsed-time fields.

## Binary identity

| Binary | Before SHA-256 | After SHA-256 |
| --- | --- | --- |
| tipseek | `53b0ccddb32a47ff044fcd3a58bdfe81c3bac208af27fdfccc4d99823a89f567` | `23bf227685d0b223259a614cae3cc7892a018d9c0985aaebd78ce7cee9a9c91c` |
| main_assembler | `cced86b8d624c0f8c3718bdb08bcb2cdf19916699b6bce5b267eb87f1973786a` | `cced86b8d624c0f8c3718bdb08bcb2cdf19916699b6bce5b267eb87f1973786a` |
| uce_filter | `e85adc284b019a6df15490bc2c1535046b8af115b428b0fa52582fcd3f463346` | `e85adc284b019a6df15490bc2c1535046b8af115b428b0fa52582fcd3f463346` |
| marker_profile | `bc596be0d807ab8881933653b23d497022dfb541a25ed6b52c401d7abd765908` | `bc596be0d807ab8881933653b23d497022dfb541a25ed6b52c401d7abd765908` |
