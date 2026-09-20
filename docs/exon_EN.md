# Exon/intron annotation

[中文](exon_ZH.md) · [gene workflow](gene_EN.md) · [command-line guide](../manual/EN_US/command_line.md)

`--assembly-mode exon` extends the default gene-family recovery workflow with protein-guided exon/intron annotation. It is intended for nuclear gene families with matching protein references, not for UCEs. The result contains coordinate-preserving gene models, resolve-eligible CDS and proteins, and explicit audit records for rejected or ambiguous models.

## Inputs and command

The nucleotide and protein directories must use matching family stems. Nucleotide baits may contain several species; the corresponding protein file must use the `.faa` extension.

```text
family_reference/
├── geneA.fasta
└── geneB.fa

family_proteins/
├── geneA.faa
└── geneB.faa
```

Run exon mode as a complete workflow; do not add stage subcommands such as `filter` or `assemble`.

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv \
  -r family_reference \
  --gene-protein-reference family_proteins \
  -o exon_output \
  -p 8
```

Required inputs are the sample table, nucleotide family baits, and matching protein references. Annotation requires miniprot 0.18 or later. TipSeek first writes recovered candidates to `exon_output/gene/`, then writes structural results to `exon_output/exon/`.

## Method

1. The default gene workflow recruits, refilters, and assembles reads, then summarizes candidates across samples.
2. Miniprot runs independently on each candidate contig. This prevents one partial candidate from suppressing a complementary candidate through target competition.
3. Embedded PAF supplies protein coordinates, coverage, score, identity, protein CIGAR, frameshifts, and in-frame stops. GFF3 supplies exon/CDS coordinates, strand, and phase.
4. Overlapping models compete deterministically by structural defects, protein coverage, alignment score, identity, splice evidence, and miniprot rank. Non-overlapping models remain independent.
5. Selected models are classified, written to structured audit tables, and routed to resolve-eligible or unresolved output.

`N` and other IUPAC bases are preserved in genomic outputs, so reported coordinates always refer to the candidate sequence used for annotation.

## Model acceptance

Only selected `complete` and `terminal_partial` models can enter `cds/` and later `gene-resolve`. An eligible model must have:

- no frameshift or internal stop;
- explicit, continuous CDS phase;
- CDS length divisible by three;
- translation without `X`;
- no unrecognized noncanonical splice;
- no unresolved model competition.

A single terminal stop codon is retained in the annotated CDS. It and the corresponding terminal `*` are removed together only when protein and codon alignments are prepared. Other states—including `low_coverage`, `frameshifted`, `internal_stop`, `phase_inconsistent`, `ambiguous_cds`, and `ambiguous_model`—remain under `unresolved/` and in the manifests.

## Validated N-padding

N-padding is attempted only when there is one unambiguous pair of selected `terminal_partial` models that covers complementary ends of the same protein and no complete selected model already represents that protein.

TipSeek orients the two source contigs, inserts `--gene-fragment-padding` Ns (100 by default), and reruns miniprot on the derived sequence. The join is accepted only when reannotation produces one complete, resolve-eligible forward model and the entire padded interval lies inside a predicted intron without overlapping any exon.

Consequently:

- accepted CDS and protein sequences contain no synthetic bases;
- padding is retained only in `supercontigs/`, intron output, GFF3, and audit tables;
- source partial models become audit-only with `superseded_by_padded_join`;
- ambiguous pairs and failed reannotations are not joined;
- padding validates a structural bridge but does not estimate the true intron length.

Set `--gene-fragment-padding 0` to disable this step.

## Main options

| Option | Default | Meaning |
| --- | ---: | --- |
| `--gene-miniprot` | `miniprot` | miniprot executable |
| `--gene-max-intron` | `50000` | Maximum intron length accepted by miniprot |
| `--gene-min-model-coverage` | `0.20` | Minimum protein coverage for a partial model |
| `--gene-complete-coverage` | `0.80` | Minimum coverage for a complete model; both protein ends are also required |
| `--gene-fragment-padding` | `100` | Ns inserted for unique two-fragment validation; `0` disables joining |
| `--gene-flank` | `0` | Observed bases added on each side in `genes_flanked/`; never N-padded |

## Outputs

```text
exon_output/
├── gene/                         # recovered candidates and family summaries
└── exon/
    ├── cds/                      # resolve-eligible CDS; no padded N in exons
    ├── proteins/                 # one-to-one with cds/
    ├── exons/                    # one record per supported CDS exon
    ├── introns/                  # one record per predicted intron
    ├── genes/                    # observed single-contig gene spans
    ├── genes_flanked/
    ├── supercontigs/             # selected spans; validated joins may contain intronic N
    ├── gff3/                     # normalized gene/mRNA/exon/CDS/intron hierarchy
    ├── raw_miniprot/             # per-candidate output and padded validation runs
    ├── models/
    │   ├── gene_models.tsv
    │   └── gene_segments.tsv
    ├── manifest/                 # candidate, ID, warning, fragment, and command audits
    └── unresolved/
```

Use `models/gene_models.tsv` for model-level status and eligibility, `models/gene_segments.tsv` for coordinates and splice classes, and `manifest/fragment_groups.tsv` for every attempted padded join. Continue with:

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8
```

See the [gene workflow guide](gene_EN.md) for resolve QC and strict or multicopy species-tree inference.
