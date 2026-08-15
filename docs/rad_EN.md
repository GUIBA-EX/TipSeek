# RAD: supplement an existing ipyrad loci matrix with WGS

`rad` supplements new WGS samples missing from an ipyrad `.loci` matrix. It does not invent RAD alleles from absent RAD reads.

| Start with | Get | Boundary |
| --- | --- | --- |
| Completed `.loci` or demultiplexed paired RAD reads + new WGS reads | Independent R1/R2 arms, recovery state, and a validated strict matrix | Does not infer the R1/R2 insert or directly prove allele dropout |

## Recommended route

```bash
# Build a reusable probe from a completed ipyrad matrix
cli/tipseek rad-probe --ipyrad-loci assembly.loci -o rad_probe

# Recover independent arms for new WGS samples
cli/tipseek rad --rad-probe rad_probe/rad_reference \
  -f wgs_samples.tsv -o rad_out -p 8

# Write a separate validated strict matrix
cli/tipseek rad-validate --rad-probe rad_probe/rad_reference \
  --rad-recovery rad_out/rad_recovery -o rad_validate_out
```

## Inputs and probe construction

Prefer `--ipyrad-loci FILE` with a completed `.loci`. Alternatively, use `--ipyrad-params FILE` to run ipyrad (steps `1--7` by default), or `--rad-denovo` with demultiplexed paired RAD reads.

`rad-probe` does not demultiplex reads, identify restriction sites, or trim adapters. Native de novo mode builds a conservative probe; it does not replace full ipyrad clustering and matrix construction.

```text
paired_rad_samples.tsv
sample_id<TAB>R1.fastq.gz<TAB>R2.fastq.gz
```

De novo mode finds candidates with paired-arm multi-seed sketches, then confirms each locus with full-length R1 and R2 distances. Its representative bait is always an observed read. Defaults are `k=31`, depth `3`, support from `2` samples, arm length `60 bp`, and maximum edit distance `3` per arm. Use `--rad-overhang-r2` only when R2 starts at a known second restriction end.

## Recovery and validation

`rad` accepts only new WGS samples absent from the probe. It stops on duplicate samples, mismatched R1/R2 inputs, or malformed arm FASTA.

The default k31 recruitment is the fastest path. Run it first. If recovery is insufficient, try in order:

1. `--rad-linked-recruitment`: offer paired fragments hitting one arm to its sibling arm, capped at 256 fragments by default. Adjust with `--rad-link-max-fragments`.
2. `--rad-fallback-kmers 25`: use a shorter k-mer only for fragments missed by k31. This is slower and less specific.

`rad-validate` requires target breadth `0.80`, identity `0.90`, and an own-locus score at least 5% above the best foreign locus by default. Only samples whose R1 and R2 both pass enter the strict matrix. Complete contigs remain in `rad_recovery/`; matrices contain only the target interval.

## Outputs and interpretation

- `rad_reference/arms/`: multi-allelic R1/R2 bait FASTA per locus.
- `rad_recovery/`: per-WGS-sample recruitment, refiltering, and assembly outputs.
- `rad_matrix/rad_sample_locus.tsv`: arm and joint recovery status per sample × locus.
- `rad_matrix/recovered_arms/`: unaligned FASTA retaining individually supported arms.
- `rad_matrix/paired_arms/`: unaligned FASTA for WGS samples with both recovered arms; not yet validated.
- `rad_validated/rad_validation.tsv`: per sample × locus × arm validation metrics and status.
- `rad_validated/strict_arms/`: original RAD baits plus only WGS samples whose two arms validate.

R1 and R2 are always independent observations. The workflow never bridges their unsequenced insert. All arm FASTA are unaligned and should be aligned only after choosing an appropriate missing-data strategy. `rad-validate` does not recruit or assemble reads again, and it leaves `rad_recovery/` and `rad_matrix/` unchanged.

`rad_missing_wgs_recovered` means that a locus absent from the input RAD matrix was recovered from WGS. It is not, by itself, evidence for restriction-site allele dropout; that interpretation additionally needs restriction-site and cross-locus WGS evidence.
