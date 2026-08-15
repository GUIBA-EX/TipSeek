# 3. Profiling

[中文版本](profiling_ZH.md)

Profiling is an assembly-free workflow for recovering reference-level evidence for any amplicon marker from WGS or metagenomic reads. It reports compatibility with reference sequences, not contigs or organismal abundance.

## Workflow

```text
TipSeek k-mer recruitment → Themisto pseudoalignment → reference-level support
```

It runs neither `refilter`, `assemble`, `combine`, nor `tree`. `--profile-kmer-size` sets the same odd k-mer size, from 15 to 31, for both recruitment and Themisto.

## Inputs

Always required:

- One `.fa` or `.fasta` marker library passed to `-r`.
- A sample table with WGS or metagenomic reads.
- `themisto`, either on `PATH` or passed with `--profile-themisto`.

Optional:

- `--profile-group-map`, a two-column TSV: `reference_id<TAB>group`. It is an optional annotation column in the reference-level result.
- `--profile-decoy`, a FASTA of plausible non-target sequences.

The reference ID is the first whitespace-delimited FASTA-header field. A supplied map must cover every reference exactly once; duplicate rows are allowed only when they assign the same group.

## Run and cache

```bash
cli/tipseek profiling \
  -f samples.tsv -r marker_reference.fasta \
  -o output -p 8
```

Add `--profile-group-map marker_groups.tsv` only when the `group` annotation column is useful. Themisto indices use a content-addressed cache. Use `--profile-index-dir` to share it and `--profile-force-rebuild` only when rebuilding is intended. When calling `marker_profile` directly, use a separate cache or `--force-rebuild` after changing the reference, group map, or decoy.

## Outputs and interpretation

Every sample writes `marker_profile/`.

### Primary result: `marker_reference_support.tsv`

One row is written for every hit reference sequence:

- `hit_queries`: number of queries compatible with that reference.
- `fractional_queries`: shared-query support; a query with N candidates contributes `1/N` to each, so it is never counted N times.
- `singleton_queries`: queries compatible with this reference alone.
- `ambiguity_status`: `has_singleton_support` or `shared_only`.

This is evidence of compatibility with a reference sequence, not a claim that the reference is uniquely present, and not an organismal abundance.

`marker_qc.tsv` records pseudoalignment and run parameters. `marker_reference_metadata.tsv` records the Themisto color-to-reference mapping and optional group annotation.

## QC and calibration

Check pseudoaligned queries, each reference’s fractional and singleton support, and the decoy-reference evidence before interpreting a result. Use negative controls, mixtures, and depth-matched downsampling to select reference libraries and pseudoalignment thresholds. Keep the reference library, optional annotation map, k-mer size, thresholds, and decoy strategy fixed for cross-sample comparisons.

See [Filter](filter_EN.md), the [output guide](../manual/EN_US/output.md), and the [command-line guide](../manual/EN_US/command_line.md).
