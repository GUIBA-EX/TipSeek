# TipSeek 完整部署

公开的完整镜像是面向 Linux x86-64 用户的免配置发行方式。它包含 TipSeek
自身的全部 Rust 可执行文件，以及系统发育、marker profiling、群体遗传、
gene family、RAD 和 TE 工作流可明确再分发的外部程序；ADMIXTURE 是唯一按
许可证单独处理的可选扩展。参考序列、probe、原始 reads、
自定义 AliFilter 模型和可选的 `correction_multi.jl` 仍由用户提供，不会写进镜像。

这里选择容器而不是只提供一个 musl 静态二进制，是因为 musl 只能解决 TipSeek
自身 Rust 程序的兼容性，不能同时封装 MAFFT、BLAST、Python/ipyrad、Julia 等完整
工作流依赖。完整镜像把这些依赖和版本一起冻结，用户机器只需有
Apptainer/Singularity 或 Docker。

## 镜像包含的软件

| 类别 | 程序 |
| --- | --- |
| TipSeek 核心 | `cargo run -p xtask -- build` 产生的全部程序、zlib-ng 2.3.3 |
| 比对与过滤 | MAFFT 7.525、Clustal Omega 1.2.4、trimAl 1.5.1、AliFilter 1.0.1 |
| 相似性搜索与注释 | BLAST 2.17.0、Magic-BLAST 1.7.0、minimap2 2.31、miniprot 0.18 |
| 系统发育 | ASTER 1.25、FastTree 2.2.0、VeryFastTree 4.0.5、IQ-TREE 3.1.3、RAxML-NG 2.0.2 |
| Marker profiling | Themisto 3.2.2 |
| 群体遗传 | minibwa 0.7、samtools/bcftools 1.23.1、PLINK 1.90b6.21、pigz 2.8；可选 ADMIXTURE 1.3.0 扩展 |
| Gene 和 RAD 可选路线 | Julia 1.12.7、ipyrad 0.9.108 |

Magic-BLAST 1.7.0 要求 zlib 1.2，而当前 minimap2 要求 zlib 1.3，因此镜像
将 Magic-BLAST 放在独立 Conda 环境，并通过包装器调用。RAxML-NG 2.0.2
目前要求 htslib 1.23，因此 samtools 和 bcftools 固定为与其相容的 1.23.1。
这不是随意使用旧版，而是避免把 ABI 不相容的程序强行安装到同一环境。

ipyrad 0.9.108 仍会导入 Python 已移除的 `distutils`，因此镜像固定 Python
3.11；若自动选到 Python 3.12 或 3.13，`ipyrad` 虽然安装成功却无法启动。

## ADMIXTURE 许可边界

公开镜像不再分发 ADMIXTURE：其包元数据只写明 `Free for Academic Use`，没有
提供清楚的再分发条款。缺少它时，TipSeek 会记录 ADMIXTURE 不可用，群体流程的
其余步骤仍可运行。学术用户在自行审阅并接受相应条款后，可以构建仅保存在本地
的扩展镜像：

```bash
docker build -f deployment/Dockerfile -t tipseek:full .
docker build -f deployment/Dockerfile.admixture-academic \
  --build-arg TIPSEEK_BASE_IMAGE=tipseek:full \
  -t tipseek:full-academic .
```

项目工作流不会发布该扩展。直接依赖的软件许可清单见
[`THIRD_PARTY.md`](THIRD_PARTY.md)。

## Apptainer/Singularity（服务器推荐）

在仓库根目录构建开发镜像（程序版本显示为 `dev`）：

```bash
singularity build --fakeroot tipseek-full.sif deployment/apptainer/TipSeek.def
```

构建正式版本时可显式写入版本号；它会同时进入 `tipseek --version`、工作流清单、
SPDX SBOM 和镜像标签：

```bash
singularity build --fakeroot \
  --build-arg TIPSEEK_VERSION=vX.Y.Z \
  tipseek-full.sif deployment/apptainer/TipSeek.def
```

以后 OCI 镜像发布后，普通用户不需要编译，可以直接转换：

```bash
singularity pull tipseek-full.sif docker://ghcr.io/guiba-ex/tipseek:full
```

检查镜像内全部依赖：

```bash
singularity exec tipseek-full.sif tipseek-doctor --profile all
```

运行示例：

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

```bash
docker build -f deployment/Dockerfile -t tipseek:full .
docker run --rm tipseek:full tipseek-doctor --profile all
```

正式版本的本地构建可增加 `--build-arg TIPSEEK_VERSION=vX.Y.Z`。仓库的
GitHub Actions 会在发布 `v*` 标签时自动注入真实标签并发布
`ghcr.io/guiba-ex/tipseek:full`，普通用户不需要自行编译。

首次发布后，维护者还需在组织的 Package settings 中将该包设为 `Public`；这是
GitHub 新包默认私有导致的一次性操作，工作流本身不能替维护者决定公开可见性。

写宿主机结果目录时建议传入当前 UID/GID，避免产生 root 所有的文件：

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v /data:/data \
  tipseek:full \
  filter assemble -f /data/samples.tsv -r /data/references \
  -o /data/output -p auto --assembly-mode uce
```

## 自检和可追溯性

`tipseek-doctor` 支持以下检查范围：

```text
core  phylogeny  profiling  population  gene  rad  te  all
```

构建镜像时会自动运行 `tipseek-doctor --profile all`；缺少必需程序或 ipyrad
无法实际启动都会使构建失败，ADMIXTURE 则单独报告为可选。镜像还会保存完整
Conda 包清单 `/opt/tipseek/conda-*-explicit.txt`，以及 TipSeek 自身的 SHA-256、
SPDX、许可证和第三方软件清单。Conda 包随附的许可证原文保存在
`/opt/tipseek/licenses/conda/`。

该镜像面向 Linux x86-64。容器隔离了宿主机的 glibc、zlib、Python 和 Conda，
但数据文件仍需通过 `--bind` 或 `-v` 映射给容器。
运行层固定为 Debian 12，以降低老集群内核的兼容风险。

Conda 求解被强制为通用 `x86_64` 基线，避免 CI 根据自身 CPU 偷选只适合新机器的
包。唯一例外是上游提供的 VeryFastTree 4.0.5 二进制要求 AVX2；老到不支持 AVX2
的 CPU 应选择 `fasttree`、`iqtree` 或 `raxmlng`，不要选择 `veryfasttree`。
