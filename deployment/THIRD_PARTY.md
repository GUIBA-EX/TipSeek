# Third-party software in the full image

TipSeek itself is licensed under GPL-3.0-or-later. The full container is an
aggregate distribution: each bundled program and library remains under its
own license. The table below records the direct command-line dependencies;
the image SBOM and `/opt/tipseek/conda-*-explicit.txt` record the complete
resolved dependency set for each build.

License files supplied inside the resolved Conda packages are preserved under
`/opt/tipseek/licenses/conda/` before the package cache is removed.

| Component | Version | License reported by the package or upstream |
| --- | --- | --- |
| ASTER | 1.25 | AGPL-3.0-or-later |
| NCBI BLAST+ | 2.17.0 | NCBI-PD |
| Clustal Omega | 1.2.4 | GPL-2.0-or-later |
| FastTree | 2.2.0 | GPL-2.0-or-later |
| IQ-TREE | 3.1.3 | GPL-2.0-or-later |
| MAFFT | 7.525 | BSD-3-Clause |
| Magic-BLAST | 1.7.0 package | NCBI-PD |
| minimap2 | 2.31 | MIT |
| miniprot | 0.18 | MIT |
| RAxML-NG | 2.0.2 | AGPL-3.0-or-later |
| trimAl | 1.5.1 | GPL-3.0-or-later |
| VeryFastTree | 4.0.5 | GPL-3.0-only |
| Themisto | 3.2.2 | GPL-2.0 |
| minibwa | 0.7 | GPL-2.0-or-later |
| samtools | 1.23.1 | MIT |
| bcftools | 1.23.1 | GPL (package metadata) |
| PLINK | 1.90b6.21 | GPL (package metadata) |
| pigz | 2.8 | Zlib |
| Julia | 1.12.7 | MIT |
| ipyrad | 0.9.108 | GPL-3.0-only |
| zlib-ng | 2.3.3 | Zlib |
| AliFilter | 1.0.1 | GPL-3.0 |

Conda-installed programs are obtained from the conda-forge and Bioconda
channels. Themisto and AliFilter are downloaded from their official GitHub
release pages and verified against SHA-256 values pinned in the container
recipe.

ADMIXTURE is deliberately not included in the public image. Its Bioconda
metadata states only `Free for Academic Use`, and the redistribution terms are
not explicit. Academic users who have reviewed and accepted the applicable
terms can build the local extension described in the deployment guide. The
extension is not published by the TipSeek container workflow.

This file is an inventory, not legal advice. Consult the upstream license text
for the terms that apply to each component.
