# TE / repeatome

`tipseek te` is a short-read repeatome workflow for genome-skimming and WGS data. It reports reproducible repeat evidence and abundance. It does not claim complete TE annotation, insertion sites, or genome-wide copy number from capture data.

## Start here

```bash
cli/tipseek te -f te_samples.tsv -o te_out -p 16
```

```text
taxon_id  sample_id  read1  read2
Taxon_A   A01        A01_R1.fq.gz  A01_R2.fq.gz
Taxon_B   B01        B01.fq.gz
```

`read2` is optional. Use `--te-library library.fa` only when a classified, preferably close-relative, library is available. Headers follow `name#Class/Subclass`.

## One workflow, two kinds of repeat

```text
reads → discover → curate → annotate → quantify
                 └→ interspersed (optional recovery lane)
```

The main lane identifies exact repeat units (EQs), preserves paired-read evidence, assembles short supported fragments, and quantifies them against all eligible reads. It is deliberately conservative.

The `interspersed` lane uses shared candidate reads rather than unique EQ assignment: it builds sparse minimizer-overlap components and jointly assembles each component. Run it when the goal is recovery of non-tandem repeat consensuses:

```bash
cli/tipseek te -f te_samples.tsv -o te_out -p 16 --te-stage interspersed
```

## Read the classes as evidence

| Output class | Meaning |
| --- | --- |
| `simple_repeat` | short periodic motif |
| `tandem_repeat_candidate` / `satellite_candidate` | repeated array; short reads do not establish chromosomal location |
| `foldback_like_DNA` | long, read-supported inverted-repeat candidate; not proof of mobility |
| `interspersed_repeat_candidate` | non-periodic component needing structural or homology evidence |
| `unknown_repeat` | insufficient evidence; retain, do not over-classify |

No Dfam or library hit does **not** mean that a candidate is not a repeat. For non-model cnidarians it usually means that a family-level label is not yet justified.

## Outputs worth reading

```text
03_annotate/annotation_evidence.tsv    class, support, period, inverted-repeat score
03_annotate/repeat_families.tsv        conservative similarity groups; EQs stay intact
03_interspersed/clusters.tsv           overlap components and joint consensus structure
03_interspersed/consensus.fasta        consensus candidates for external annotation
04_quantify/repeat_signal.tsv          per-sample EQ abundance and coverage
04_quantify/repeat_landscape.tsv       read-to-consensus divergence proxy
05_compare/repeat_superfamilies.tsv    shared, taxon-shared, or sample-specific families
```

`signal_rpm` is a relative read signal. `estimated_genome_fraction` is only a read-fraction proxy, and only interpretable as genome fraction for comparable random WGS libraries. Never make that interpretation for UCE off-target reads.

## Practical interpretation

Use tandem and foldback calls as structural observations. Use `consensus.fasta` for Dfam, HMM/domain, or curated coral-library annotation. Reserve names such as LINE, LTR, or DNA transposon for candidates with matching structural or protein-domain evidence. A stable unknown category is a valid result.
