# Mitochondrial workflow: ordinary circular animal mitochondria

`mito` uses the existing MainFilter, refilter, and Rust UCE assembler to recruit reads from an annotated GenBank reference. It reports a circular sequence only when the sample reads support it.

| Start with | Get | Main boundary |
| --- | --- | --- |
| Annotated GenBank + sample reads | Read-supported circular or partial mitochondrial sequence | Only for ordinary single circular animal mitochondria |

## Scope

This workflow is for **ordinary single circular animal mitochondrial genomes**: gene order should remain broadly comparable with the reference, without multipartite structure, major rearrangement, or heteroplasmy that needs phased reporting. It is not intended for complex plant or fungal mitochondria, fragmented or multichromosomal mitochondria, major rearrangements, or studies requiring specialised heteroplasmy or NUMT handling.

## Minimal command

```bash
cli/tipseek mito \
  -f samples.tsv \
  -o mito_output \
  -p 8 \
  --mito-genbank mitochondrial_reference.gb
```

`--mito-max-reads 320` caps adaptive input; `--no-mito-adaptive-stop` switches to one-pass filtering.
`--cleanup-intermediates` is opt-in: only after final circular confirmation does it remove `filtered/mitochondrion.fq` and adaptive stage copies, while writing `cleanup_manifest.tsv`; the default retains all intermediates.

## Workflow

```text
GenBank gene/rRNA/tRNA + genome/tile baits
→ MainFilter paired-read recruitment
→ one mitochondrial read pool
→ refilter → joint Rust UCE graph assembly
→ finalize and test circularity
→ if not yet circular: contigs_all become sample-specific seeds, recruit once more and reassemble (skipped once already circular)
→ strict overlap, unique GFA paths, and mate-link joins
→ if the primary k fails, unique local k−10 / k / k+10 graph paths only
→ junction-spanning read validation of circularity
```

All baits are written as one mitochondrial locus. The reference is used only for recruitment and seeding: final sequence is neither coordinate-stitched nor reference-filled. Mate links propose adjacency and orientation only; gap bases must be recovered from a unique path in the same filtered read pool. Unresolved gaps remain broken and are never filled with `N`. Multi-k is fallback-only and cannot bypass unique-path or junction-read requirements.

`mito` enables unlimited extension and GFA graph output by default. For distant references or low coverage, explicitly start with sensitive recruitment such as `-kf 17–25 -s 1`. Each adaptive depth is still capped by `--mito-max-reads`; set it at least as high as the input when the full library must be scanned.

## Success criteria

A circular result must have one component, no `N`, every join supported by a strict overlap or unique GFA/read-graph bridge, a closed terminus, and at least `--mito-min-junction-support` reads spanning the final junction. Junction support is measured as the *minimum* depth across a small band of k-mers tiled over the seam — a single lucky k-mer cannot stand in for consistent spanning coverage — and each read is matched allowing one mismatch on either strand, so a lone sequencing error no longer discards a genuine spanning read. A seam k-mer that also recurs inside the contig is not counted as closure evidence. Otherwise partial output is retained but the command fails.

## Common and expert parameters

- `--mito-genbank`: required annotated mitochondrial GenBank reference.
- `--mito-max-reads 320`: at most approximately 1.05M paired-read blocks per adaptive stage; the workflow stops early when two successive stages return an exactly identical circular sequence after cut/strand normalization. Stability is judged per sample: a sample is frozen only after the same verified circle is observed at two consecutive depths. A non-circular result (including no contig yet) is only a checkpoint and continues to deeper read budgets until a circle is confirmed or the maximum is reached.
- `--no-mito-adaptive-stop`: disable staged early stopping and use the normal one-pass `--max-reads` behaviour.

The following hidden expert overrides should be changed only to diagnose a known recruitment, graph-joining, or circularity problem: `--mito-min-overlap`, `--mito-min-overlap-identity`, `--mito-terminal-window`, `--mito-link-kmer`, `--mito-min-link-hits`, `--mito-min-pair-support`, `--mito-bridge-kmer`, `--mito-bridge-min-depth`, `--mito-max-bridge`, and `--mito-min-junction-support`.

## Outputs

- `<sample>/mito/mitochondrial_assembly.fasta`: circular sequence or partial components (the audited raw assembly, never rotated).
- `<sample>/mito/mitochondrial_standardized.fasta`: written only for a verified circle, rotated to a reproducible gene start (tRNA-Phe if present, else a conserved fallback) and oriented to that gene's coding strand for cross-sample comparability. Only existing bases are reordered or reverse-complemented — no reference base is added — and the header records the anchor, strand, offset, and mismatches. Omitted when no anchor can be located confidently.
- `<sample>/mito/mitochondrial_assembly_summary.tsv`: compatible coarse status, `resolution_reason`, joins, junction support, and `ambiguous_bases` / `ambiguous_per_kb` (Ns per 1000 bp of the primary contig).
- `<sample>/mito/mitochondrial_evidence.json`: machine-readable structural, graph, mate-link, and junction evidence.
- `<sample>/mito/mitochondrial_feature_evidence.tsv`: canonical 21-mer similarity to reference features; it reports exact anchor sharing only and cannot establish gene presence or absence in distant samples. `translation_status=not_checked` is not a CDS annotation or translation call.
- `<sample>/mito/mitochondrial_mate_links.tsv`: accepted read-supported links.
- `.gm2_mito_reference/metadata/mitochondrial_genes.tsv`: bait metadata in 0-based half-open coordinates; `segments_0_half_open` preserves every segment of cross-origin or `join(...)` features.
