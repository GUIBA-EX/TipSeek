# Gene and exon workflows

The default `gene` mode handles nuclear gene families defined by multi-species baits. It retains within-sample candidate contigs; `exon`, `gene-resolve`, and `gene-tree` add structural annotation and phylogenetic resolution. Candidate count is assembly evidence, **not** an allele or biological copy-number call.

## Quick start

Each `family_reference/*.fasta` defines one family and may contain several species. `family_proteins/` contains same-named protein FASTA files.

```bash
# Candidate recovery only (the default gene mode)
cli/tipseek -f samples.tsv -r family_reference -o gene_output -p 8

# Candidate recovery followed by protein-guided exon/intron annotation
cli/tipseek --assembly-mode exon \
  -f samples.tsv -r family_reference \
  --gene-protein-reference family_proteins -o exon_output -p 8

# Align, infer gene trees, and select strict one-to-one clades
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8

# Strict pseudo-SCO species tree
cli/tipseek gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral

# Multicopy-family species tree
cli/tipseek gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

The default `gene` mode needs `-f/-r/-o`. `--assembly-mode exon` additionally requires `--gene-protein-reference` and miniprot. Resolve needs MAFFT and IQ-TREE, and species-tree inference needs ASTER2 `astral`. Exon mode is an extension of the default workflow and is not available in UCE mode.

## Modes and post-processing commands

| Invocation | Input | Role | Main output |
| --- | --- | --- | --- |
| `tipseek` | reads + family baits | Recruit, refilter, assemble, and summarize candidates | `<output>/gene/` |
| `tipseek --assembly-mode exon` | reads + family baits + protein references | Default gene recovery followed by miniprot ≥0.18 structural annotation and model QC | `<output>/{gene,exon}/` |
| `gene-resolve` | `<exon-output>/exon/` | Protein MSA, codon backtranslation, gene tree, and unrooted one-to-one clade selection | `gene_resolved/` |
| `gene-tree` | `gene_resolved/` | ASTER2 strict or multicopy species tree | tree and provenance |

## Annotation and structural QC

Exon mode runs miniprot independently for each candidate contig and preserves the GFF3 plus embedded PAF produced by `--gff`. This prevents one partial contig from suppressing a complementary partial contig in miniprot's target competition. PAF supplies protein coverage, alignment score, protein CIGAR, frameshifts, and in-frame stops; GFF3 supplies CDS coordinates and phase. `N` and other IUPAC bases remain in candidates, so coordinates continue to refer to the original contig.

Models on one contig compete only when they directly meet the overlap threshold, in deterministic order by structural defects, protein coverage, alignment score, identity, splice sites, and miniprot rank. A rejected bridge model does not merge two otherwise independent loci. States include `complete`, `terminal_partial`, `low_coverage`, `frameshifted`, `internal_stop`, `phase_inconsistent`, `ambiguous_cds`, and `ambiguous_model`.

By default, only `complete` and `terminal_partial` models enter `cds/` and resolve. They must have no frameshift or internal stop, explicit continuous phase, a CDS length divisible by three, a translation without `X`, no unrecognized noncanonical splice, and no unresolved model competition. A single terminal stop codon is retained in the annotated CDS and removed together with its translated `*` only when preparing protein and codon alignments. Other models remain under `unresolved/` and in the structured manifests.

For one unambiguous pair of selected `terminal_partial` models that map to complementary ends of the same protein, exon mode orients the two complete candidate contigs, inserts 100 `N` bases by default, and reruns miniprot on that single derived sequence. The join is accepted only if the rerun yields one complete, resolve-eligible forward model and the entire padding lies inside a predicted intron without overlapping an exon. The CDS and protein therefore contain no synthetic bases; the N-padded span is retained only in `supercontigs/`, intron output, GFF3, and audit tables. The two source partials become audit-only (`superseded_by_padded_join`). Ambiguous sets or failed reannotation are not joined. This is a structural bridge, not an estimate of the true intron length.

Main exon controls:

- `--gene-max-intron`: maximum intron length; default 50,000 bp.
- `--gene-min-model-coverage`: minimum protein coverage for a partial model; default 0.20.
- `--gene-complete-coverage`: minimum complete-model coverage, also requiring both protein ends; default 0.80.
- `--gene-fragment-padding`: number of `N` bases used for the unique two-fragment validation described above; default 100, and 0 disables joining.
- `--gene-flank`: observed bases added on each side in `genes_flanked/`; default 0 and never N-padded.

## Resolve and QC

`gene-resolve` runs a fast ML tree by default. It applies two conservative QC rounds: pre-alignment retains only translatable candidates at least `--gene-min-aa-length` long (30 aa by default) and checks `--gene-min-taxa` by **distinct sample**; post-alignment checks occupancy again after MAFFT/TAPER and codon backtranslation, and requires at least `--gene-min-effective-codon-sites` effective codon sites (30 by default). `--gene-ufboot` accepts only `0` (default) or `>=1000`; only the latter supplies usable branch support in `tree_selection_qc.tsv`.

```bash
cli/tipseek gene-resolve --gene-input exon_output/exon -o gene_resolved -p 8 \
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
gene_output/gene/               # default gene-mode output
├── family_summary.tsv
├── family_count_matrix.tsv
├── pseudo_sco/
└── multiple_candidate_families/

exon_output/exon/               # structural output from --assembly-mode exon
├── cds/                         # resolve-eligible CDS only; validated joins contain no N in exons
├── proteins/                    # one-to-one with cds/
├── exons/                       # one record per protein-supported CDS exon
├── introns/                     # one record per intron
├── genes/                       # single-contig gene spans only; observed introns
├── genes_flanked/
├── supercontigs/                # selected spans; validated two-contig joins may contain intronic N
├── gff3/                        # normalized gene/mRNA/exon/CDS/intron hierarchy
├── raw_miniprot/                # per-candidate GFF3+PAF/stderr, including padded validation runs
├── models/{gene_models.tsv,gene_segments.tsv}
├── manifest/
└── unresolved/

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
