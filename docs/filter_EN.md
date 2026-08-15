# MainFilter: read recruitment

[中文](filter_ZH.md) · [output guide](../manual/EN_US/output.md) · [command-line guide](../manual/EN_US/command_line.md)

MainFilter is TipSeek's general read-recruitment layer. It scans reads with reference k-mers and assigns matching reads (or complete paired fragments) to loci for `refilter` and assembly. It is not an aligner and does not itself decide contig validity, marker abundance, or population structure. Default UCE does not use this two-stage route; it uses fused `ucefilter` instead; see [Assembler](assembler_EN.md).

## Place in the workflow

```text
FASTQ/FASTA + one reference FASTA per locus
                ↓
            MainFilter
                ↓
  filtered_pe/ (broad first-pass recruitment)
                ↓
 refilter (disambiguation and depth/size limits)
                ↓
       filtered/ → assembler
```

In `original`, gene, and `--legacy-uce-filter` routes, a hit in either mate retains the complete paired fragment. A core-mapping mate can therefore retain informative flank sequence from its partner. Default UCE applies the paired-fragment rule inside one `ucefilter` scan. Profiling uses the first recruitment only and does not run `refilter`.

## Inputs and references

- The sample table is tab-separated: `sample<TAB>R1<TAB>R2`; omit R2 for single-end data. Each non-comment row must have exactly two or three columns, normalized sample names must be unique, and every listed read file must exist.
- Input may be FASTA, FASTQ, or gzipped FASTQ. Formats must be consistent within one run.
- `original` and `uce` use a reference directory with one `.fa`/`.fasta` per locus. The file stem is the unique locus name.
- Profiling instead accepts one marker-library FASTA file.

`-kf` is the shared reference/read k-mer length and `-s` is the sampling step along reads. The terminal window is always checked so that a read end is not skipped by the step. Windows containing N or another non-ACGTU character are not matched.

## Matching and assignment

MainFilter uses a canonical representation of each k-mer and its reverse complement, so one lookup covers both strands. `-gr` is retained only for compatibility. A read may recruit to more than one locus in the broad first pass; `refilter` handles that ambiguity later.

`ref_reads_count_dict.txt` records reads recruited to each locus in the first pass. It is neither final assembly coverage nor a gene/allele copy-number estimate.

## Cache

Use `--reuse-reference-cache` for an unchanged reference. A dictionary stores the canonical policy, k-mer length, locus names, and a SHA-256 of reference content. A content, k-mer, or format mismatch rejects the old cache and rebuilds it rather than silently reusing it.

The scan step does not alter the reference dictionary; it only changes which read positions are probed. Share a cache directory only among jobs using the same reference and `-kf`.

## Outputs and checks

- `ref_reads_count_dict.txt`: first-pass read count per locus.
- `filtered_pe/`: first-pass reads; normally removable after successful refiltering.
- `filtered/`: refiltered reads supplied to the assembler.
- `large_files/`: appears only when a refilter depth or size limit is reached.

Zero recruited reads do not by themselves demonstrate biological absence. Check reference divergence, `-kf`, `-s`, sequencing depth, and the count table. If many loci are simultaneously zero, first check the reference directory, sample-table paths, and input format.

## Rust MainFilter versus the historical implementation

The historical MainFilter baseline is maintained upstream and is not bundled in this Rust-only repository. The Rust implementation preserves command-line semantics, canonical bidirectional recruitment, and output format. Its purpose is to reduce CPU, allocation, and I/O overhead per read—not to obtain speed by increasing `-s`.

| Area | Typical historical cost | Rust MainFilter approach |
| --- | --- | --- |
| k-mer scan | Interpreter loops and string/object work | 2-bit rolling k-mers: `u64` for `k≤32`, `u128` for `33–64`, byte keys only above 64. |
| Strand matching | Separate forward/reverse work | One canonical key for both strands. |
| Multi-locus hits | Multiple container objects and indirection | Locus IDs packed into contiguous `packed_hits` slices. |
| Reused references | Repeated reference parsing/building | Binary dictionary cache validated by a content SHA-256; safely rebuilt on mismatch. |
| FASTA/FASTQ parsing | `String`, UTF-8 validation, and uppercase copies | Byte-level `read_until` with reusable line/record buffers; the DNA table accepts both cases. |
| gzip input | Small default buffers and more FFI/system calls | 1 MiB outer reader and zlib `gzbuffer`; runtime zlib-ng when available, otherwise system zlib. |
| Per-locus output | Frequent open/close and global-flush jitter | Persistent file handles, per-locus buffers, high/low-water flushing, a bounded buffer pool, and large-buffer release. |
| Encoding and logging | Temporary objects and many small writes | Reused text/GM2 scratch buffers and buffered logging. |

Performance depends on reads, reference size, k, and the filesystem. In a real DK40 UCE test (`k=33`, 3,579 loci, 1 million read pairs), per-locus output decreased from 18.35 s to 3.72 s and peak memory from about 375 MiB to 186 MiB, with byte-identical output files and read counts. See the [performance notes](development/mainfilter-performance.md) for validation details.

## Practical rules

1. Validate reference choice and conservative `-kf`/`-s` values for recovery accuracy before considering speed.
2. Reuse the reference cache across samples; do not manually mix caches built from different references.
3. Treat MainFilter as broad recruitment, not an abundance or absence test.
4. Use normal filtering output for gene/UCE analyses. The internal scan-only mode produces no reads for assembly and is not a recommended analysis path.
