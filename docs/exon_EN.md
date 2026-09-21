# Exon/intron annotation

[中文](exon_ZH.md) · [gene workflow](gene_EN.md) · [command-line guide](../manual/EN_US/command_line.md)

`--assembly-mode exon` extends the default gene-family recovery workflow with protein-guided exon/intron annotation. By default it derives protein references directly from coding nucleotide baits, so separate `.faa` files are not required; external proteins remain an optional override. This mode is not for UCEs.

## Inputs and command

Each `.fa`, `.fas`, or `.fasta` file under `-r` defines one family and may contain coding nucleotide baits from several species, such as CDS or concatenated exons that preserve coding frame.

```text
family_reference/
├── geneA.fasta
└── geneB.fa
```

These nucleotide references are translated automatically with the standard nuclear genetic code. If a bait contains introns, is non-coding, uses a non-standard genetic code, or cannot be interpreted as coding sequence, supply an external protein directory instead. Each override must use the matching family stem and `.faa` extension. The directory may be incomplete: families with `.faa` files use them, while the remaining families still use automatic translation.

Run exon mode as a complete workflow; do not add stage subcommands such as `filter` or `assemble`.

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv \
  -r family_reference \
  -o exon_output \
  -p 8
```

To override selected families with external proteins:

```bash
cli/tipseek --assembly-mode exon \
  -f samples.tsv -r family_reference \
  --gene-protein-reference family_proteins \
  -o exon_output -p 8
```

Annotation requires miniprot 0.18 or later. TipSeek first writes recovered candidates to `exon_output/gene/`, then writes structural results to `exon_output/exon/`.

## Method

1. The default gene workflow recruits, refilters, and assembles reads, then summarizes candidates across samples.
2. Without a matching `.faa`, TipSeek examines all six reading frames equally. A terminal stop is allowed and removed from the protein reference; any frame with an internal stop is excluded rather than repaired or masked.
3. A sequence is accepted directly only when it has one stop-free frame or one frame has uniquely strongest CDS endpoint support from an initial methionine and/or terminal stop. The longest such translation becomes the family anchor. Remaining ambiguous sequences are accepted only when one frame has a unique positive amino-acid 3-mer Dice match to that anchor. With no trustworthy anchor, a tied match, or no positive match, TipSeek emits no protein for that sequence and records `ambiguous_reading_frame`; a matching `.faa` is then required if the family has no other usable reference. Every candidate and decision is recorded in `manifest/reference_translation.tsv`.
4. Miniprot runs independently on each candidate contig. This prevents one partial candidate from suppressing a complementary candidate through target competition.
5. Embedded PAF supplies protein coordinates, coverage, score, identity, protein CIGAR, frameshifts, and in-frame stops. GFF3 supplies exon/CDS coordinates, strand, and phase.
6. Overlapping models compete deterministically by structural defects, protein coverage, alignment score, identity, splice evidence, and miniprot rank. Non-overlapping models remain independent.
7. Selected models are classified, written to structured audit tables, and routed to resolve-eligible or unresolved output.

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

N-padding is attempted only when there is one unambiguous pair of selected `terminal_partial` models that covers complementary ends of the same protein and no complete selected model already represents that protein. A candidate pair must also overlap by no more than 10% of the shorter model, add coverage equal to at least 20% of the protein length, reach `--gene-complete-coverage` in total, and reach both protein termini within a 2% (minimum 3-aa) tolerance.

TipSeek orients the two source contigs, inserts `--gene-fragment-padding` Ns (100 by default), and reruns miniprot on the derived sequence. The join is accepted only when reannotation produces one complete, resolve-eligible forward model and the entire padded interval lies inside a predicted intron without overlapping any exon.

This validation depends on each source contig retaining enough genuine intronic and splice-boundary evidence beside the break. If two fragments end exactly at exon boundaries and omit the intronic splice flanks entirely, reannotation will usually be unable to identify the padded interval reliably as an intron, so the join is rejected. N-padding also cannot reconstruct an unrecovered terminal exon or any other coding sequence.

Consequently:

- accepted CDS and protein sequences contain no synthetic bases;
- padding is retained only in `supercontigs/`, intron output, GFF3, and audit tables;
- source partial models become audit-only with `superseded_by_padded_join`;
- ambiguous pairs and failed reannotations are not joined;
- padding validates a structural bridge but does not estimate the true intron length.

In a test using real *Patiria pectinifera* genes, all three intact candidates containing a genuine 2,000-nt intron recovered the correct CDS, strand, exon/intron boundaries, and phase without `.faa` input. Among complementary fragments retaining 120 nt of genuine intronic flank on each side, the pair that met the gates above passed reannotation with 100 Ns and reproduced the true CDS exactly; two other pairs remained unjoined because one overlapped excessively in protein coordinates and the other had a terminal fragment below the coverage gate. This validates conservative acceptance, not guaranteed recovery of complementary terminal fragments from read assembly.

Set `--gene-fragment-padding 0` to disable this step.

## Main options

| Option | Default | Meaning |
| --- | ---: | --- |
| `--gene-protein-reference` | auto-derived | Optional `.faa` directory; matching families override automatic translation |
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
    ├── manifest/
    │   ├── reference_translation.tsv # automatic frame/strand selection audit
    │   ├── derived_proteins/     # automatically derived proteins actually used
    │   └── ...                   # candidate, ID, warning, fragment, and command audits
    └── unresolved/
```

Use `manifest/reference_translation.tsv` for reference source and frame selection, `models/gene_models.tsv` for model-level status and eligibility, `models/gene_segments.tsv` for coordinates and splice classes, and `manifest/fragment_groups.tsv` for every attempted padded join. Continue with:

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8
```

See the [gene workflow guide](gene_EN.md) for resolve QC and strict or multicopy species-tree inference.
