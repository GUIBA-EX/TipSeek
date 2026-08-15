# Rust CLI full-migration matrix

The release CLI runs through the Rust dispatcher without Python or Biopython at
runtime. The former Python implementation has been removed; release execution has no Python dependency.

| Command family | Native backend(s) | Rust CLI state | Completion evidence |
| --- | --- | --- | --- |
| `filter` / `refilter` / `assemble` original | `MainFilterNew`, `main_refilter_new`, `main_assembler-original-rust` | partial (native standalone + gene route) | paired/single-end recovery and generic original-workflow fixtures |
| UCE recovery | `uce_filter`, `main_assembler-rust`, `gm2_tools` | partial (paired native and legacy-candidate paths) | paired native/legacy fixtures; single-end and shadow failure fixtures remain |
| UCE rescue | UCE recovery backends | partial (whole-contig + terminal-only rounds with per-end read-evidence reconciliation) | broader biological fixtures and end-to-end parity fixtures remain |
| `mito` | `mito_workflow`, UCE tools | partial (native full skeleton) | bait collapse, text refilter, seed rescue and circular adaptive-stop fixture; biological circular/linear fixtures pending |
| `rad`, `rad-probe`, `rad-validate` | `rad_workflow`, optional ipyrad, MainFilter, refilter, original-rust | native route | synthetic `.loci`, inferred ipyrad output, strict/phylogeny arm matrix, partial-arm and CLI-dispatch fixtures |
| `gene` | `gene_workflow`, original-rust tools | partial | single/paired recovery and cohort fixtures |
| gene annotate/resolve/tree | `gene_workflow` + external tools | partial (standalone routes) | annotation dispatch + ASTER/provenance fixtures; resolve fixture pending |
| `profiling` | `marker_profile` | partial (native complete route) | single-marker recruitment/quantification fixture; cache and decoy fixtures pending |
| `population` | `main_population` + external tools | partial (standalone route) | default/panref option forwarding fixtures; staged-output fixture pending |
| `te` | `main_repeat` | partial (standalone route) | default/stage/optional-library forwarding fixture; full stage fixture pending |
| `consensus` | `build_consensus` + minimap2 | partial (standalone route) | filtered FASTQ discovery, mapper/consensus invocation and SAM cleanup fixture |
| `trim` / `combine` / `tree` | `gm2_tools` + external tools | partial (standalone routes) | trim, UCE-aware combine, and coalescent-tree fixtures; full MSA/filter/concat fixtures pending |
| `stats` | `gm2_stats` | complete as a standalone command | native report output fixture |
| cleanup / profiling / cache | Rust CLI filesystem layer | native route | cleanup manifest, reference-cache, failed-profile, and per-stage profile fixtures |

## Migration invariants

1. A completed native command must not invoke `python`, `unix_command.py`, or Biopython.
2. Existing command spelling, defaults, output names, and failure boundaries remain compatible unless explicitly versioned.
3. A cleanup action is allowed only after every native consumer of its input has completed successfully.
4. Each migrated family gets at least one end-to-end fixture and one failure-path fixture before the default engine changes.
5. `compare_outputs` remains available until every matrix row is complete.

## Cutover criterion

`cli/tipseek` points to the Rust binary. Each matrix row remains subject to
its listed fixture coverage; No Python implementation or dispatcher remains in the repository.
