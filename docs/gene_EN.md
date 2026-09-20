# Gene-family recovery and phylogenetic resolution

[中文](gene_ZH.md) · [exon annotation](exon_EN.md) · [command-line guide](../manual/EN_US/command_line.md)

The default `gene` mode recovers candidate contigs for nuclear gene families defined by multi-species nucleotide baits. Candidate count is assembly evidence, **not** an allele count or biological copy-number estimate.

## Candidate recovery

Each `.fa` or `.fasta` file in `family_reference/` defines one family and may contain several species.

```bash
cli/tipseek \
  -f samples.tsv \
  -r family_reference \
  -o gene_output \
  -p 8
```

TipSeek recruits and refilters reads, assembles candidates for every sample and family, and writes the cohort summary to `gene_output/gene/`. It does not assign exon coordinates or claim single-copy status.

For protein-guided CDS, exon, intron, and supercontig output, run `--assembly-mode exon`; see the [exon annotation guide](exon_EN.md). `gene-resolve` consumes that structural `exon/` output, not the unannotated `gene/` directory.

## Workflow boundaries

| Invocation | Input | Role | Main output |
| --- | --- | --- | --- |
| `tipseek` | reads + nucleotide family baits | Candidate recovery and cohort summary | `<output>/gene/` |
| `tipseek --assembly-mode exon` | reads + nucleotide and protein references | Candidate recovery plus structural annotation | `<output>/{gene,exon}/` |
| `tipseek gene-resolve` | `<exon-output>/exon/` | MSA, codon backtranslation, gene trees, and one-to-one clade selection | `gene_resolved/` |
| `tipseek gene-tree` | `gene_resolved/` | Strict or multicopy ASTER2 species tree | tree and provenance |

## Resolve and QC

`gene-resolve` retains only models marked resolve-eligible by exon annotation. It then applies two further QC rounds:

- pre-alignment: translatable candidates of at least `--gene-min-aa-length` (default 30 aa), with `--gene-min-taxa` counted by **distinct sample**;
- post-alignment: occupancy checked again after MAFFT, optional TAPER, and codon backtranslation, with at least `--gene-min-effective-codon-sites` effective codon sites (default 30).

The default is a fast ML gene tree without branch support. `--gene-ufboot` accepts `0` or at least `1000`; only the latter provides branch support for `tree_selection_qc.tsv`.

```bash
cli/tipseek gene-resolve \
  --gene-input exon_output/exon \
  -o gene_resolved \
  -p 8 \
  --gene-outgroup outgroups.tsv \
  --gene-taper /path/to/correction_multi.jl \
  --gene-julia julia \
  --gene-ufboot 1000
```

- `--gene-outgroup`: first TSV/CSV column lists outgroup sample IDs; they must be monophyletic in a gene tree.
- `--gene-taper`: optional masking after protein MSA; malformed, duplicate, or missing headers reject the family to unresolved.
- `occupancy_qc.tsv`: retained candidates, distinct-sample occupancy, median length, thresholds, and rejection reason per family.
- `family_qc.tsv`: alignment QC for families that passed post-alignment checks; it is not an overall success flag.
- `tree_selection_qc.tsv`: occupancy, multicandidate-sample count, and branch support for each strict clade.
- `resolve_manifest.tsv`: final resolved or unresolved status for every family.

## Species trees

```bash
# Strict pseudo-SCO route: one selected leaf per sample
cli/tipseek gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral

# Multicopy route: complete gene trees plus leaf-to-sample mapping
cli/tipseek gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

Both routes require ASTER2 `astral`. The strict route normalizes each selected subtree to one leaf per sample. The multicopy route retains complete gene trees and supplies `leaf_to_species.tsv`. `gene_tree_provenance.tsv` records the command, input paths, and SHA-256 values.

## Main outputs

```text
gene_output/gene/
├── family_summary.tsv
├── family_count_matrix.tsv
├── pseudo_sco/
└── multiple_candidate_families/

gene_resolved/
├── resolved_1to1/                 # CDS and audit tree per strict clade
├── unresolved_multicandidate/     # multicopy, conflicting, or failed families
├── astral_input/resolved_1to1.trees
├── astralpro_input/{multicopy.trees,leaf_to_species.tsv}
├── occupancy_qc.tsv
├── family_qc.tsv
├── tree_selection_qc.tsv
└── resolve_manifest.tsv
```
