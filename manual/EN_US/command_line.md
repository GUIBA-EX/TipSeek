# TStools Command-Line Guide

TStools (formerly GeneMiner2-UCE) is a command-line workflow for UCEs and other phylogenetic markers. This repository does not include the former GUI, bundled demonstration data, or legacy graphical documentation.

## 1. Building from source

### 1.1 Build dependencies

A complete build requires:

- a C/C++ compiler and [zlib](https://zlib.net/);
- [Rust and Cargo](https://www.rust-lang.org/tools/install/);

Rust/Cargo is required for the complete current build, including the primary read filter, secondary filter, assembler, population workflow, alignment cleanup, sequence merging, reference trimming, and UCE statistics.

On Ubuntu, first install the system dependencies:

```bash
sudo apt install build-essential zlib1g zlib1g-dev
```

### 1.2 Runtime dependencies

The phylogenetic workflow calls a subset of the following external programs, depending on the selected stages:

```bash
conda install -c bioconda \
  aster blast clustalo fasttree iqtree mafft magicblast miniprot minimap2 raxml-ng trimal veryfasttree
```

Optional and population-specific tools:

- `--alignment-filter alifilter` requires a separate AliFilter installation with `AliFilter` in `PATH`.
- `population` requires `minibwa`, samtools, bcftools, and PLINK 1.9.
- ADMIXTURE performs ancestry analysis. If it is absent, population analysis still completes cohort-reference construction, VCF generation, PLINK export, and PCA, while recording its status as `unavailable`.

Population executables are located through `PATH` by default. Override them individually with `--population-minibwa`, `--population-samtools`, `--population-bcftools`, `--population-plink`, and `--population-admixture`.

### 1.3 Download and build

```bash
git clone --depth 1 https://github.com/GUIBA-EX/GeneMiner2-UCE.git
cd GeneMiner2-UCE
cargo run -p xtask -- build
```

The generated entry point is:

```bash
cli/geneminer2 -h
```

Run the same Cargo command again after updating the source.

## 2. Input files

Most read-recovery commands require the following. `gene-annotate`, `gene-resolve`, and `gene-tree` instead use `--gene-input` and `-o`, without `-f/-r`:

- `-f FILE`: tab-delimited sample table;
- `-r DIR`: reference-sequence directory;
- `-o DIR`: output directory.

### 2.1 Sample table

Each row represents one sample. Single-end data use two columns; paired-end data use three:

```text
Sample_A	/data/reads/Sample_A_R1.fq.gz	/data/reads/Sample_A_R2.fq.gz
Sample_B	/data/reads/Sample_B_R1.fq.gz	/data/reads/Sample_B_R2.fq.gz
Sample_C	/data/reads/Sample_C_R1.fq.gz	/data/reads/Sample_C_R2.fq.gz
```

Absolute paths are recommended. Sample names must be unique. For `population`, optional `population` and `batch` fields may follow R2; leave the R2 field empty for single-end samples with metadata. Population analysis writes a manifest linking the original names to internal directory and VCF sample IDs.

### 2.2 Reference directory

Use one FASTA file per locus; the file name becomes the locus name. A file may contain one or several reference sequences:

```text
references/
  uce-0001.fasta
  uce-0002.fasta
```

## 3. Subcommands and default workflows

One or more subcommands can be listed in execution order:

| Subcommand | Function |
| --- | --- |
| `filter` | Recruit reads using reference k-mers |
| `profiling` | Recruit marker reads once and estimate group-level marker signal without assembly |
| `refilter` | Refine per-locus read assignment and filtering |
| `assemble` | Assemble target sequences with the wDBG assembler |
| `gene` | Recover candidate contigs for nuclear gene families |
| `gene-annotate` | Run miniprot annotation from protein references |
| `gene-resolve` | Align candidates, infer gene trees, and select strict one-to-one clades |
| `gene-tree` | Infer a species tree from strict or multicopy gene trees |
| `te` | Discover, curate, annotate, and quantify conservative repeatome units from short reads |
| `population` | Build a cohort UCE reference and generate SNP, PCA, and ADMIXTURE results |
| `consensus` | Generate consensus sequences at heterozygous sites |
| `trim` | Remove flanks relative to reference sequences |
| `combine` | Merge samples, align loci, clean sequences, and filter alignment columns |
| `tree` | Infer coalescent or concatenated trees |
| `stats` | Summarize UCE recovery and optionally generate heatmaps |

When no subcommand is given:

- `--assembly-mode original` (default) is for reference-guided recovery of exons, SCOs, and nuclear or mitochondrial markers; it runs `filter refilter assemble trim combine tree`;
- `--assembly-mode uce` is for UCE recovery from genome skimming or target capture; it runs `filter assemble combine tree`. The fused UCEFilter already includes refilter semantics and omits `trim` so newly recovered UCE flanks are not cut back to the reference interval;
- `profiling` runs one recruitment step followed by Themisto pseudoalignment and reference-level support reporting; it does not assemble or run downstream phylogenetic steps.

Default original-mode example:

```bash
cli/geneminer2 \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references \
  -o /home/user/project/output \
  -p 8
```

## 4. UCE assembly and marker profiling

### 4.1 UCE

UCE mode recovers UCE cores and read-supported flanks from genome-skimming or target-capture reads. The fused `ucefilter` performs recruitment, run-k verification, and adaptive per-locus selection in one scan. Low-depth loci pass through, while only saturated loci have redundant core evidence reduced with bait/contig-edge overhang ladders protected; final per-locus FASTQ is written directly. A complete paired-end fragment is retained atomically, with no GM2 or candidate FASTQ. The default workflow skips `trim` to preserve recovered flanks.

Optional `--uce-alignment-shadow` collects bounded internal alignment evidence without changing adaptive read selection. It is disabled by default.

```bash
cli/geneminer2 \
  -f samples.tsv -r references -o output -p 8 \
  --assembly-mode uce
```

For assembly behavior, backend selection, rescue, cache semantics, and QC, see the [Assembler chapter](../../docs/assembler_EN.md). Use this manual for option definitions and the [output guide](output.md) for file fields.

### 4.2 Marker profiling

Profiling performs one recruitment followed by Themisto pseudoalignment and reference-level support reporting; it does not assemble. It requires one `.fa` or `.fasta` marker library. `--profile-group-map` is optional and adds a group annotation column.

```bash
cli/geneminer2 profiling \
  -f samples.tsv -r marker_reference.fasta \ -o output -p 8
```

Inputs, decoys, cache control, QC, and quantitative interpretation are in the [Profiling chapter](../../docs/profiling_EN.md).

### 4.3 Gene-family recovery and resolution

`gene` is a complete workflow: one bait FASTA defines one family and may contain multiple species. It fixes the backend to `original-rust`, retains candidate contigs, and does not directly claim single-copy status.

```bash
cli/geneminer2 gene -f samples.tsv -r family_reference -o gene_output -p 8
cli/geneminer2 gene-annotate --gene-input gene_output/gene \
  --gene-protein-reference family_proteins -o gene_annotation -p 8
cli/geneminer2 gene-resolve --gene-input gene_annotation -o gene_resolved -p 8
```

`gene-resolve` requires MAFFT and IQ-TREE; `--gene-taper correction_multi.jl` enables optional masking. It applies pre-alignment QC by distinct-sample occupancy and `--gene-min-aa-length` (30 aa by default), then post-alignment QC by occupancy and `--gene-min-effective-codon-sites` (30 by default); see `occupancy_qc.tsv`. `--gene-ufboot` must be `0` (default) or `>=1000`. `family_qc.tsv` is alignment QC for post-QC families, while `tree_selection_qc.tsv` records selected strict clades and occupancy.

```bash
cli/geneminer2 gene-tree --gene-input gene_resolved -o species_strict -p 8 \
  --gene-species-mode strict --gene-aster astral
cli/geneminer2 gene-tree --gene-input gene_resolved -o species_multi -p 8 \
  --gene-species-mode multicopy --gene-aster astral
```

The strict route supplies one leaf per sample to ASTER2. The multicopy route also supplies the candidate-to-sample map. Both require ASTER2 `astral`.

### 4.4 TE / repeatome

`te` is a complete standalone workflow and does not require `-r`. It uses a four-field manifest, `taxon_id sample_id read1 read2` (omit the fourth field for single-end data), and runs `discover → curate → annotate → quantify`.

```bash
cli/geneminer2 te -f te_samples.tsv -o te_output -p 32
```

`--te-library` optionally supplies a classified `name#Class/Subclass` FASTA. Annotation never merges EQs and does not replace complete-TE annotation. See the [TE / repeatome chapter](../../docs/te_EN.md) for thresholds, rerun rules, and output interpretation.

## 5. Population-genetic analysis

### 5.1 Scope and example

`population` uses multiple completed UCE assemblies and their original reads to create a cohort pseudo-reference, joint VCF, PCA, and ADMIXTURE inputs. Each sample must retain `uce_assembly_summary.csv`, accepted contigs in `results/`, and the reads listed in the sample table.

```bash
cli/geneminer2 population \
  -f /home/user/project/samples.tsv \
  -r /home/user/project/references -o output -p 8 \
  --assembly-mode uce --engine panrefv2 \
  --population-admixture-k-min 2 --population-admixture-k-max 6
```

For pseudo-reference strategies, stage restarts, SNP panels, and required QC, see the [Population workflow guide](../../docs/population_EN.md). Individual missingness, heterozygosity, and pairwise PI_HAT are written under `population/structure/qc/`; use `--skip-relatedness-qc` for large cohorts. Use `--population-start-at` only with validated outputs from the requested stage.

## 6. Staged-run examples

Rebuild only a concatenated tree from existing results:

```bash
cli/geneminer2 tree \
  -f samples.tsv -r references -o output \
  -m concatenation --phylo-program iqtree
```

Generate consensus sequences, trim them to references, and combine them:

```bash
cli/geneminer2 consensus trim combine \
  -f samples.tsv -r references -o output \
  -c 0.75 -ts consensus -tm all -tr 0.5 \
  -cs trimmed
```

Generate tables from existing UCE results without plotting heatmaps:

```bash
cli/geneminer2 stats \
  -f samples.tsv -r references -o output \
  --stats-no-heatmap
```

## 7. Parameter reference

The tables below list the main public options and current defaults. Run `cli/geneminer2 -h` for the complete help associated with the checked-out source.

### 7.1 General input and parallelism

| Option | Description |
| --- | --- |
| `-f FILE` | Sample table; required |
| `-r DIR` | Reference-sequence directory; required |
| `-o DIR` | Output directory; required |
| `-p INT\|auto` | Shared CPU budget; default `auto`. Auto counts affinity/cpuset-allowed physical cores and applies cgroup or scheduler limits; an integer overrides detection. |

### 7.2 Read filtering and refiltering

| Option | Description |
| --- | --- |
| `-kf INT` | Filter k-mer size; default `23` in UCE mode and `31` otherwise |
| `-s, --step-size INT` | Read-scanning step; default `4` |
| `--max-reads INT` | Maximum million reads processed per file; `0` means unlimited |
| `--reuse-reference-cache` | Reuse a fingerprinted reference k-mer index; with explicit `original-rust`, also enable its versioned, k-validated binary assembler cache |
| `--reference-cache-dir DIR` | Reference-cache directory; default `output/.gm2_reference_cache`; requires the preceding option |
| `--depth-low-water-mark INT` | Below this depth, attempt relaxed read recruitment; default `50` |
| `--depth-limit INT` | Maximum depth processed during refiltering; default `768` |
| `--file-size-limit INT` | Refiltering file-size limit; default `6` |

### 7.3 Assembly and UCE options

| Option | Description |
| --- | --- |
| `-ka INT` | Assembly k-mer size; default `0` for automatic estimation |
| `--min-ka INT` / `--max-ka INT` | Automatic-estimation range; defaults `21` / `51` |
| `-e, --error-threshold INT` | k-mer error threshold; default `2` |
| `-sb, --soft-boundary VALUE` | Integer, `auto`, or `unlimited`; default `auto` |
| `-i, --search-depth INT` | Search depth; default `4096` |
| `--min-coverage INT` | Minimum contig read depth; default `0` |
| `--assembler-implementation MODE` | `auto` (default) uses `original-rust` in original mode and `uce-rust` in UCE mode; `uce-rust` selects the UCE-oriented Rust assembler; `original` and `original-rust` use the deterministic Rust compatibility implementation; UCE never falls back to another implementation |
| `--assembler-read-chunk-size INT` | Reads loaded per Rust assembler batch; default `8192` |
| `--assembler-kmer-count-threads INT` | K-mer sorting/counting workers per locus; default `0` selects automatically |
| `--assembler-graph-format MODE` | Optional graph output: `none` (default), `gfa`, `dot`, or `both` |
| `--assembly-mode MODE` | `original` or `uce`; default `original` |
| `--assembly-mode uce` | Defaults to k=23/step=4, automatic sensitive recruitment, one bounded rescue round, and fixed safe backbone/QC settings |
| `--uce-recruit-mode fast\|auto` | Recruitment policy; default `auto` in UCE mode and `fast` otherwise |
| `--uce-rescue-reads` | Explicitly enable fixed-k=21 bounded rescue; retained for compatibility because UCE mode enables it by default |
| `--no-uce-rescue-reads` | Disable the default UCE rescue stage |
| `--uce-rescue-rounds 1\|2` | Number of rescue rounds; default `1` |
| `--uce-rescue-inverted-repeat-min-bp INT` | Revert a locus when the current rescue round newly introduces an exact inverted repeat of at least this length; default `150`, `0` disables |
| `--uce-rescue-reverse-reuse-reference-scale FLOAT` | Experimental: scale the reference bonus for a node whose reverse complement is already present in either assembly arm; range `0`--`1`, default `1.0` (disabled), with read depth unchanged |

### 7.4 Population options

| Option | Description |
| --- | --- |
| `--engine MODE` | `pseudoref` (default), legacy `panref`, or `panrefv2`; the latter two use the per-locus bait directory supplied with `-r` |
| `--population-panrefv2-include-low-confidence` | Also write `short` or `low_sample_support` PanRefV2 loci to the mapping FASTA; default keeps only `pass` |
| `--population-reference-strategy MODE` | `pseudoref` only: `sqcl-longest` (default) or `supported` |
| `--population-reference-fasta FILE` | Use a fixed external FASTA as the cohort reference; it is copied into `population/reference/` and has no per-sample contribution table |
| `--population-min-mapq INT` / `--population-min-baseq INT` | Minimum MAPQ / base quality; defaults `20` / `20` |
| `--population-min-dp INT` / `--population-min-gq INT` | Set lower-quality genotypes to missing; defaults `5` / `20` |
| `--population-min-qual FLOAT` | Minimum site QUAL; default `20` |
| `--population-min-call-rate FLOAT` | Minimum non-missing genotype fraction; default `0.8` |
| `--population-min-mac INT` | Minimum minor allele count; default `2` |
| `--population-ld-window INT` / `--population-ld-step INT` | LD-pruning window and step; defaults `50` / `5` SNPs |
| `--population-ld-r2 FLOAT` | LD-pruning r² threshold; default `0.2` |
| `--population-admixture-k-min INT` / `--population-admixture-k-max INT` | ADMIXTURE K range; defaults `2` / `6`; maximum K is capped at sample count |
| `--population-admixture-cv INT` | ADMIXTURE CV folds; default `10`, capped at sample count |
| `--population-start-at STAGE` | Start at `reference` (default), `mapping`, `calling`, or `selection`; later stages reuse validated existing reference, BAM, or filtered VCF outputs |
| `--population-stop-after STAGE` | Stop after `reference`, `mapping`, `calling`, or `selection`; default `selection` |
| `--population-skip-mark-duplicates` | Skip samtools duplicate marking |
| `--population-skip-plink` | Omit PLINK, PCA, LD-pruned, and ADMIXTURE outputs |
| `--population-skip-admixture` | Generate PLINK and PCA without running ADMIXTURE |
| `--population-minibwa PATH` | minibwa executable; default `minibwa` |
| `--population-samtools PATH` | samtools executable; default `samtools` |
| `--population-bcftools PATH` | bcftools executable; default `bcftools` |
| `--population-plink PATH` | PLINK 1.9 executable; default `plink` |
| `--population-admixture PATH` | ADMIXTURE executable; default `admixture` |

### 7.5 Consensus, trimming, and combining

| Option | Description |
| --- | --- |
| `-c, --consensus-threshold FLOAT` | Consensus threshold; default `0.75` |
| `-ts, --trim-source SOURCE` | `assembly` or `consensus`; default is the preceding stage's output |
| `-tm, --trim-mode MODE` | `all`, `longest`, `terminal`, or `isoform`; default `terminal` |
| `-tr, --trim-retention FLOAT` | Reference-trimming retention fraction; default `0` |
| `-cs, --combine-source SOURCE` | `assembly`, `consensus`, or `trimmed`; default is the preceding stage's output |
| `-cd, --clean-difference FLOAT` | Maximum acceptable pairwise alignment difference; default `1.0` |
| `-cn, --clean-sequences INT` | Minimum sequences per alignment; default `0` |
| `--msa-program PROGRAM` | `clustalo` or `mafft`; default `mafft` |
| `--msa-threads INT` | Threads per MSA job; default `1` and cannot exceed `-p` |
| `--alignment-filter PROGRAM` | `trimal`, `alifilter`, or `none`; default `trimal` |
| `--filter-processes INT` | Concurrent alignment-filter jobs; default equals `-p` |
| `--alifilter-model MODEL` | AliFilter model specification or `model.json` path |
| `--strict-combine-errors` | Stop when any locus fails alignment, cleanup, or filtering |
| `--no-alignment` | Skip multiple-sequence alignment |
| `--no-trimal` | Deprecated alias for `--alignment-filter none` |

### 7.6 Tree inference and statistics

| Option | Description |
| --- | --- |
| `-m, --tree-method METHOD` | `coalescent` or `concatenation`; default `coalescent` |
| `-b, --bootstrap INT` | Bootstrap replicates; default `1000` |
| `--phylo-program PROGRAM` | `raxmlng`, `iqtree`, `fasttree`, or `veryfasttree`; default `fasttree` |
| `--stats-no-heatmap` | Do not generate UCE statistics heatmaps |
| `--stats-count-input-reads` | Count input FASTQ reads for `InputReads` and `PctFiltered`; slow on large datasets |

Fractional thresholds use decimals from `0` to `1`. Before changing filters on a full dataset, run a small test with defaults and inspect intermediate VCFs, mapping QC, and locus occupancy.
