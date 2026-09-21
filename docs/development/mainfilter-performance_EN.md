# MainFilter performance and compatibility

[中文](mainfilter-performance.md)

This note records a real UCE comparison between the current Rust `MainFilterNew` and the historical upstream Haxe/C++ baseline, together with performance design and compatibility boundaries. It addresses engineering performance of read recruitment only; runtime is not a proxy for UCE recovery, assembly accuracy, or population-genetic inference quality.

## Result: real Haxe/C++ comparison

The historical MainFilter baseline is the C++/hxcpp implementation compiled from Haxe, whose source remains in the upstream repository. The historical binary was rebuilt from that source and compared with the current Rust binary using identical real UCE reads, references, options, and output mode.

| Scenario | Historical Haxe/C++ | Current Rust | Improvement |
| --- | ---: | ---: | ---: |
| Build cache + filter | 20.32 s | 3.52 s | **82.7% faster** (5.8×) |
| Reuse cache + filter | 18.35 s | 3.72 s | **79.7% faster** (4.9×) |
| Peak RSS while building cache | 381 MiB | 187 MiB | **51.0% lower** |
| Peak RSS while reusing cache | 375 MiB | 186 MiB | **50.3% lower** |

The filtering stage itself decreased from 18.59 s (first run) / 18.06 s (cache reuse) in Haxe/C++ to 3.42 s / 3.57 s in Rust. Wall time includes gzip decompression, FASTQ parsing, k-mer recruitment, per-locus output, and count-table writing.

## Benchmark protocol

### Implementations

- Historical baseline: upstream `sculab/GeneMiner2` MainFilter source rebuilt with Haxe 4.3.3, hxcpp 4.3.2, `g++ -O2`, and the build options specified by its README.
- Candidate implementation: release `MainFilterNew` from this repository's `rust/main_filter_new` crate.
- Each implementation creates and loads its own dictionary cache; the cache formats are incompatible and are never shared.

### Fixed inputs

- Reads: `validation/reads/fastq/DK40__SRR29729138_{1,2}.fastq.gz`
- Reference: `validation/references/bivalve_uce_2k_v1_loci/`, containing 3,579 loci.
- Options: `-kf 33 -s 3 -m_reads 1 -m 0 -gr`.
- `-m_reads 1` limits each input file to `2^20` records, or 1 Mi read pairs.
- Output: `filtered_pe/`, with one text read file per locus.

Each binary was run twice: first without a cache, then loading the cache produced by its own first run. `/usr/bin/time` recorded wall time, user/system time, and maximum RSS.

## Output equivalence

Recruitment and output were checked as follows:

| Check | Result |
| --- | --- |
| Non-empty locus output files | 2,228 in both implementations |
| Total output bytes | 460,812,556 bytes in both implementations |
| File-by-file comparison of `filtered_pe/` | Identical under `diff -qr` |
| Values in `ref_reads_count_dict.txt` | Identical |

The raw count files cannot be compared directly with `cmp`: the historical Haxe version writes in unordered hash-table order and formats each line as `locus,count,` with a trailing comma, while Rust writes `locus,count` in reference order. After removing the Haxe trailing comma and sorting, all 2,228 locus counts are identical. This is a formatting and row-order difference, not a read-recruitment difference.

## Rust optimizations

| Area | Main cost in the historical Haxe/C++ path | Current Rust approach | Semantic boundary |
| --- | --- | --- | --- |
| k-mer representation | Long-mode and generic-container overhead | 2-bit rolling keys: `u64` for `k≤32`, `u128` for `33–64`, and byte keys only for `k>64`. | Canonical keys, ambiguous bases, and read-end checks are preserved. |
| Strand matching | Separate forward/reverse k-mer handling | Canonical k-mers cover both strands with one lookup. | `-gr` remains only for compatibility. |
| Multi-locus hits | Scattered containers and indirection | `ReferenceHits` plus contiguous `packed_hits`; measured structure size is 12 bytes. | The locus-ID set is unchanged. |
| Index lookup | Generic hash tables and object access | Native integer/byte keys in `AHashMap`. | Hash-table iteration order is not output semantics. |
| Dictionary cache | Older format and many small I/O operations | v4 cache, reference-content SHA-256, a 4 MiB `BufWriter`, and direct loading into `packed_hits`. | Mismatched, damaged, or stale caches are rebuilt safely. |
| FASTA/FASTQ | Haxe string/object path | Byte-level `read_until`, reused line/record buffers, and no whole-read uppercase copy. | Bytes emitted for valid FASTA/FASTQ remain unchanged. |
| gzip | Historical `GzipReader` path | 1 MiB outer-reader and zlib `gzbuffer`; build-time native zlib-ng via `pkg-config`, with runtime probing and system-zlib fallback. | Decompressed bytes are unchanged; zlib-ng remains optional. |
| Output | Buffer flushing and repeated file operations | Persistent handles, per-locus buffers, and batched flushing at 64 MiB high/32 MiB low watermarks. | Record content and ordering in each file remain unchanged. |
| Memory reclamation | Peak capacity for hot loci remains resident | Bounded buffer pool; buffers larger than 1 MiB are released after flushing. | Already written and subsequent records are unaffected. |
| Encoding/logging | Temporary objects and small writes | Reused text/GM2 scratch buffers and buffered logs. | Encoding and log semantics are unchanged. |

## Cache behavior

The dictionary cache depends only on reference content, locus names, k-mer length, and index format. Scan step `-s` changes probe positions on reads but not the reference dictionary. A shared cache is valid only for work using the same reference and `-kf`.

The current Rust cache format is v4. Older Rust caches, historical Haxe caches, damaged caches, and caches for changed reference content are rejected and rebuilt; an incompatible index is never loaded silently.

## Exclusions from workflow-level performance claims

- **Scan-only mode** does not produce assemblable reads. It is an internal diagnostic and not a formal gene/UCE performance result.
- **Increasing `step`** changes sampling density and potential sensitivity. It is an analytical parameter change, not a lossless optimization.
- **Concurrent writes to formal per-locus outputs** introduce merging, ordering, and I/O-contention issues when multiple threads write the same locus set. The production path retains one writer and byte-level validation.
- **Alternative indexes such as SSHash/GGCAT** would change the static-index and cache design. They should be implemented and revalidated separately only after the current index is shown to be the dominant bottleneck.

## Validation requirements after modification

After changing MainFilter, run at least:

```bash
cargo test --manifest-path rust/main_filter_new/Cargo.toml
cargo clippy --manifest-path rust/main_filter_new/Cargo.toml -- -D warnings
cargo build --release --manifest-path rust/main_filter_new/Cargo.toml
```

Changes to recruitment logic, k-mer encoding, caches, parsing, or output also require an old/new binary comparison on the fixed real paired-end sample:

```bash
diff -qr old_run/filtered_pe new_run/filtered_pe
# For count files, normalize trailing commas and sort by locus before comparing values.
```

Release a candidate binary only when outputs remain equivalent, tests pass, and real-data timing shows a stable improvement.
