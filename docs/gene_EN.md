# Gene subcommands

`gene` handles nuclear gene families defined by multi-species baits. It retains within-sample candidate contigs, then separates reliable one-to-one subtrees from multicopy or ambiguous families with protein annotation and gene trees. Candidate count is assembly evidence, **not** an allele or biological copy-number call.

| Start with | Get | Main boundary |
| --- | --- | --- |
| Reads + multi-species family baits | Candidate contigs, strict one-to-one clades, or multicopy tree inputs | Does not directly claim biological copy number |

## Quick start

Each `family_reference/*.fasta` defines one family and may contain several species. `family_proteins/` contains same-named protein FASTA files.

```bash
# 1. Recover candidates from reads; original-rust is fixed
cli/tipseek gene -f samples.tsv -r family_reference -o gene_output -p 8

# 2. Protein-guided annotation
cli/tipseek gene-annotate --gene-input gene_output/gene \
  --gene-protein-reference family_proteins -o gene_annotation -p 8

# 3. Align, infer gene trees, and select strict one-to-one clades
cli/tipseek gene-resolve --gene-input gene_annotation -o gene_resolved -p 8

# 4a. Strict pseudo-SCO species tree
cli/tipseek gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral

# 4b. Multicopy-family species tree
cli/tipseek gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

`gene` needs `-f/-r/-o`; the other three subcommands need only `--gene-input` and `-o`. Annotation needs miniprot, resolve needs MAFFT and IQ-TREE, and species-tree inference needs ASTER2 `astral`.

## Subcommands

| Subcommand | Input | Role | Main output |
| --- | --- | --- | --- |
| `gene` | reads + family baits | Recruit, refilter, original-rust assembly, and candidate summary | `gene/` |
| `gene-annotate` | `gene/` + protein references | miniprot CDS, exon, intron, and supercontig extraction | `gene_annotation/` |
| `gene-resolve` | `gene_annotation/` | Protein MSA, codon backtranslation, gene tree, and unrooted one-to-one clade selection | `gene_resolved/` |
| `gene-tree` | `gene_resolved/` | ASTER2 strict or multicopy species tree | tree and provenance |

## Resolve and QC

`gene-resolve` runs a fast ML tree by default. It applies two conservative QC rounds: pre-alignment retains only translatable candidates at least `--gene-min-aa-length` long (30 aa by default) and checks `--gene-min-taxa` by **distinct sample**; post-alignment checks occupancy again after MAFFT/TAPER and codon backtranslation, and requires at least `--gene-min-effective-codon-sites` effective codon sites (30 by default). `--gene-ufboot` accepts only `0` (default) or `>=1000`; only the latter supplies usable branch support in `tree_selection_qc.tsv`.

```bash
cli/tipseek gene-resolve --gene-input gene_annotation -o gene_resolved -p 8 \
  --gene-outgroup outgroups.tsv \
  --gene-taper /path/to/correction_multi.jl --gene-julia julia \
  --gene-ufboot 1000
```

- `--gene-outgroup`: first TSV/CSV column lists outgroup sample IDs; they must be monophyletic in a gene tree.
- `--gene-taper`: runs TAPER after AA MSA; malformed, duplicate, or missing headers are rejected to unresolved.
- `occupancy_qc.tsv`: pre/post retained candidates, distinct-sample occupancy, median length, thresholds, and rejection reason per family; multiple candidates from one sample count once.
- `family_qc.tsv`: alignment QC only for families passing post-alignment QC (`alignment_pass`), not an overall resolve-success call.
- `tree_selection_qc.tsv`: candidate occupancy, multicandidate-sample count, and branch support per strict clade.
- `resolve_manifest.tsv`: final resolved/unresolved reason for each family.

## Outputs and interpretation

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

The strict route normalizes every selected subtree to **one leaf per sample** for ASTER2. The multicopy route retains complete gene trees and maps candidate leaves to samples with `leaf_to_species.tsv`. `gene-tree` writes `gene_tree_provenance.tsv` with its command, inputs, and SHA-256 values.
