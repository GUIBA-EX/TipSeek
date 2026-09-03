# TipSeek full deployment

The public full image is the supported zero-setup Linux x86-64 distribution.
It contains the native TipSeek executables plus the redistributable programs
used by the phylogeny, marker-profiling, population, gene-family, RAD, and TE
workflows. ADMIXTURE is the one license-restricted optional extension.
Reference sequences, probes, input reads, custom AliFilter models, and optional
`correction_multi.jl` scripts remain user inputs and are not embedded.

A single musl-linked binary is not sufficient for the full distribution:
musl can make TipSeek's Rust executables more portable, but it cannot bundle
the MAFFT, BLAST, Python/ipyrad, Julia, and other workflow runtimes. The full
container freezes those programs and versions together, leaving only an
Apptainer/Singularity or Docker runtime as a host requirement.

## Included command-line tools

| Area | Programs |
| --- | --- |
| TipSeek core | all executables produced by `cargo run -p xtask -- build`, zlib-ng 2.3.3 |
| Alignment and filtering | MAFFT 7.525, Clustal Omega 1.2.4, trimAl 1.5.1, AliFilter 1.0.1 |
| Similarity and annotation | BLAST 2.17.0, Magic-BLAST 1.7.0, minimap2 2.31, miniprot 0.18 |
| Phylogeny | ASTER 1.25, FastTree 2.2.0, VeryFastTree 4.0.5, IQ-TREE 3.1.3, RAxML-NG 2.0.2 |
| Marker profiling | Themisto 3.2.2 |
| Population | minibwa 0.7, samtools/bcftools 1.23.1, PLINK 1.90b6.21, pigz 2.8; optional ADMIXTURE 1.3.0 extension |
| Gene and RAD optional routes | Julia 1.12.7 and ipyrad 0.9.108 |

Magic-BLAST is intentionally installed in a separate Conda environment because
its zlib 1.2 ABI constraint conflicts with current minimap2's zlib 1.3
constraint. The `/usr/local/bin/magicblast` wrapper selects the isolated binary.
Similarly, samtools and bcftools are pinned to 1.23.1 to share the htslib ABI
required by RAxML-NG 2.0.2. Both combinations are resolved during CI rather than
being forced into an inconsistent environment.

ipyrad 0.9.108 still imports Python's retired `distutils` module, so the image
pins Python 3.11 instead of selecting Python 3.12 or 3.13 and shipping a command
that fails at startup.

## ADMIXTURE license boundary

ADMIXTURE is not redistributed in the public image because its package metadata
states only `Free for Academic Use` and does not provide clear redistribution
terms. TipSeek records ADMIXTURE as unavailable while the rest of the population
workflow remains usable. Academic users who have reviewed and accepted the
applicable terms can create a local extension:

```bash
docker build -f deployment/Dockerfile -t tipseek:full .
docker build -f deployment/Dockerfile.admixture-academic \
  --build-arg TIPSEEK_BASE_IMAGE=tipseek:full \
  -t tipseek:full-academic .
```

That extension is never published by the project workflow. See
[`THIRD_PARTY.md`](THIRD_PARTY.md) for the direct dependency license inventory.

## Apptainer/Singularity

Build a development image from the repository root (reported version `dev`):

```bash
singularity build --fakeroot tipseek-full.sif deployment/apptainer/TipSeek.def
```

For a release build, inject the release version into `tipseek --version`, the
workflow manifest, SPDX SBOM, and image label:

```bash
singularity build --fakeroot \
  --build-arg TIPSEEK_VERSION=vX.Y.Z \
  tipseek-full.sif deployment/apptainer/TipSeek.def
```

After an OCI image is published, a cluster user can instead convert it without
building the software:

```bash
singularity pull tipseek-full.sif docker://ghcr.io/guiba-ex/tipseek:full
```

Run the deployment check:

```bash
singularity exec tipseek-full.sif tipseek-doctor --profile all
```

Run an analysis, binding the host data directory into the container:

```bash
singularity run --bind /data:/data tipseek-full.sif \
  filter assemble \
  -f /data/samples.tsv \
  -r /data/references \
  -o /data/output \
  -p auto \
  --assembly-mode uce
```

## Docker/OCI

Build and check locally:

```bash
docker build -f deployment/Dockerfile -t tipseek:full .
docker run --rm tipseek:full tipseek-doctor --profile all
```

Add `--build-arg TIPSEEK_VERSION=vX.Y.Z` for a local release build. On a `v*`
release tag, the repository workflow injects the actual tag and publishes
`ghcr.io/guiba-ex/tipseek:full`; end users do not need to compile the image.

After the first publication, an administrator must make the package `Public`
in Package settings; new GitHub packages are private by default.

Run with the current user's UID/GID so output files are not owned by root:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v /data:/data \
  tipseek:full \
  filter assemble -f /data/samples.tsv -r /data/references \
  -o /data/output -p auto --assembly-mode uce
```

## Deployment checks

`tipseek-doctor` accepts these profiles:

```text
core  phylogeny  profiling  population  gene  rad  te  all
```

The image build runs `tipseek-doctor --profile all`; a missing required
executable or a broken ipyrad runtime makes the build fail. ADMIXTURE is reported
separately as optional. Published images also retain explicit Conda package
manifests at `/opt/tipseek/conda-*-explicit.txt`, alongside TipSeek's SHA-256,
SPDX, license, and third-party inventory files. License texts supplied by the
resolved Conda packages are preserved under `/opt/tipseek/licenses/conda/`.

## Portability boundary

The image targets Linux x86-64. Containers isolate TipSeek from the host glibc,
zlib, Python, and Conda installations, but they do not hide user data. Bind the
read, reference, and output locations explicitly on clusters that do not bind
them automatically.

The runtime base is pinned to Debian 12 to reduce compatibility risk on older
cluster kernels.

The Conda solve is forced to the baseline `x86_64` archspec so a CI runner does
not silently select packages tuned only for its own CPU. The bundled
VeryFastTree 4.0.5 binary is the exception: its upstream package uses AVX2. On a
pre-AVX2 CPU, select `fasttree`, `iqtree`, or `raxmlng` instead of `veryfasttree`.
