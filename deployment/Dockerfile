# syntax=docker/dockerfile:1.7

ARG TIPSEEK_BUILDER_IMAGE=rust:1.87.0-bookworm@sha256:6d79f767859a5f025a062505fa9f2c1a041cadafcee71fbcbd226223be462f18
ARG TIPSEEK_RUNTIME_IMAGE=mambaorg/micromamba:2.9.0-debian12-slim@sha256:4f29126d15f60e722538868c7d95f5ed905c9b0f0e9eb833f60f0b966d0b7652
ARG TIPSEEK_VERSION=dev

FROM ${TIPSEEK_BUILDER_IMAGE} AS builder
ARG TIPSEEK_VERSION

ARG THEMISTO_VERSION=3.2.2
ARG THEMISTO_SHA256=1ae148e2eff478ad444da61c5b43f8293ea9c52b72940c7aadd610a0f01217d0
ARG ALIFILTER_VERSION=1.0.1
ARG ALIFILTER_SHA256=db9cb301d25e67edccd4efb644f9619b520062f15c557010d049330b9411a300
ARG ALIFILTER_LICENSE_SHA256=585e25ef8f5946a52bf2aed68d5becfc38be94a8663aa01c1b31d88aa57f1de3

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        binutils \
        build-essential \
        ca-certificates \
        curl \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# The portable container binary deliberately links only the baseline system
# zlib. The runtime image supplies zlib-ng through the existing dlopen path,
# avoiding an absolute RUNPATH to a build-host Conda environment.
ENV TIPSEEK_RELEASE_VERSION=${TIPSEEK_VERSION}
ENV PKG_CONFIG_PATH=/nonexistent
ENV PKG_CONFIG_LIBDIR=/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/share/pkgconfig
RUN cargo run -p xtask -- build \
    && ! readelf -d cli/bin/tipseek-rust | grep -q '/home/' \
    && ! readelf -d cli/bin/tipseek-rust | grep -q 'libz-ng'

RUN mkdir -p /third-party/themisto /third-party/alifilter \
    && curl -fL --retry 3 \
        -o /tmp/themisto.tar.gz \
        "https://github.com/algbio/themisto/releases/download/v${THEMISTO_VERSION}/themisto_linux-v${THEMISTO_VERSION}.tar.gz" \
    && echo "${THEMISTO_SHA256}  /tmp/themisto.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/themisto.tar.gz -C /third-party/themisto --strip-components=1 \
    && curl -fL --retry 3 \
        -o /tmp/alifilter.tar.gz \
        "https://github.com/arklumpus/AliFilter/releases/download/v${ALIFILTER_VERSION}/AliFilter-linux-x64.tar.gz" \
    && echo "${ALIFILTER_SHA256}  /tmp/alifilter.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/alifilter.tar.gz -C /third-party/alifilter \
    && curl -fL --retry 3 \
        -o /third-party/alifilter/LICENSE.md \
        "https://raw.githubusercontent.com/arklumpus/AliFilter/v${ALIFILTER_VERSION}/LICENSE.md" \
    && echo "${ALIFILTER_LICENSE_SHA256}  /third-party/alifilter/LICENSE.md" \
        | sha256sum -c -

FROM ${TIPSEEK_RUNTIME_IMAGE} AS runtime

ARG TIPSEEK_VERSION
LABEL org.opencontainers.image.title="TipSeek full" \
      org.opencontainers.image.description="TipSeek with its redistributable workflow dependencies" \
      org.opencontainers.image.source="https://github.com/GUIBA-EX/TipSeek" \
      org.opencontainers.image.version="${TIPSEEK_VERSION}"

COPY --chown=$MAMBA_USER:$MAMBA_USER deployment/conda/environment-full.yml /tmp/environment-full.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER deployment/conda/environment-magicblast.yml /tmp/environment-magicblast.yml

# Do not let the CI host CPU silently select x86-64-v3 Conda variants.
ENV CONDA_OVERRIDE_ARCHSPEC=x86_64

RUN micromamba install --yes --name base --file /tmp/environment-full.yml \
    && micromamba create --yes --name tipseek-magicblast --file /tmp/environment-magicblast.yml \
    && license_archive=/opt/conda/tipseek-conda-licenses \
    && mkdir -p "$license_archive" \
    && for license_dir in /opt/conda/pkgs/*/info/licenses; do \
        [ -d "$license_dir" ] || continue; \
        package_dir="${license_dir%/info/licenses}"; \
        package_name="${package_dir##*/}"; \
        mkdir -p "$license_archive/$package_name"; \
        cp -aL "$license_dir"/. "$license_archive/$package_name"/; \
    done \
    && micromamba clean --all --yes

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates libgomp1 zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/cli /opt/tipseek/cli
COPY --from=builder /third-party/themisto/themisto /usr/local/bin/themisto
COPY --from=builder /third-party/alifilter/AliFilter /usr/local/bin/AliFilter
COPY LICENSE /opt/tipseek/LICENSE
COPY --from=builder /third-party/themisto/LICENSE.txt /opt/tipseek/licenses/Themisto-LICENSE.txt
COPY --from=builder /third-party/alifilter/LICENSE.md /opt/tipseek/licenses/AliFilter-LICENSE.md
COPY deployment/THIRD_PARTY.md /opt/tipseek/THIRD_PARTY.md
COPY --chmod=0755 deployment/bin/magicblast /usr/local/bin/magicblast
COPY --chmod=0755 deployment/bin/tipseek-doctor /usr/local/bin/tipseek-doctor
COPY --chmod=0755 deployment/bin/tipseek-container-entrypoint /usr/local/bin/tipseek-container-entrypoint

RUN mv /opt/conda/tipseek-conda-licenses /opt/tipseek/licenses/conda \
    && cp -a /opt/conda/lib/libz-ng.so* /usr/local/lib/ \
    && ldconfig \
    && ln -s /opt/tipseek/cli/tipseek /usr/local/bin/tipseek \
    && chmod 0755 /usr/local/bin/themisto /usr/local/bin/AliFilter \
    && mkdir -p /work \
    && chown "$MAMBA_USER:$MAMBA_USER" /work \
    && micromamba list --name base --explicit > /opt/tipseek/conda-base-explicit.txt \
    && micromamba list --name tipseek-magicblast --explicit > /opt/tipseek/conda-magicblast-explicit.txt

ENV TIPSEEK_HOME=/opt/tipseek
ENV PATH=/opt/tipseek/cli:/usr/local/bin:/opt/conda/bin:$PATH
WORKDIR /work
USER $MAMBA_USER

# This catches missing executables and broken packaging during image creation.
RUN tipseek --help >/dev/null \
    && tipseek --version >/dev/null \
    && magicblast -version >/dev/null \
    && tipseek-doctor --profile all

ENTRYPOINT ["/usr/local/bin/tipseek-container-entrypoint"]
CMD ["--help"]
