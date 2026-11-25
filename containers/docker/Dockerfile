# syntax=docker/dockerfile:1
# =============================================================================
# One image per pixi environment (see pixi.toml / pixi.lock).
#
#   docker buildx build -f containers/docker/Dockerfile \
#       --build-arg PIXI_ENV=flecs --platform linux/arm64 -t benchmark/flecs:latest .
#
# CUDA variants (PIXI_ENV=<env>-cuda) are linux/amd64 only; the GPU is not
# needed at build time (CONDA_OVERRIDE_CUDA satisfies the __cuda requirement).
# Method code is NOT copied: the benchmark mounts the repository at run time
# (record the git SHA + image digest with the results for reproducibility).
# =============================================================================

# Pinned by digest (multi-arch: linux/amd64, linux/arm64)
ARG PIXI_IMAGE=ghcr.io/prefix-dev/pixi:0.84.2@sha256:8201f134655a8261f382ee565ddf1498e5f215887bb14e812e1f47c09f4f7196
ARG RUNTIME_IMAGE=ubuntu:24.04@sha256:562111b88caa79c71391c0c4d3c2114a7d5c35915c20561fd893cb1e07929e6b

# ── Build: install exactly the locked environment ────────────────────────────
FROM ${PIXI_IMAGE} AS build
ARG PIXI_ENV
ENV CONDA_OVERRIDE_CUDA=12.0
WORKDIR /app
COPY pixi.toml pixi.lock ./
RUN test -n "${PIXI_ENV}" \
    && pixi install --locked -e "${PIXI_ENV}" \
    && pixi shell-hook -e "${PIXI_ENV}" -s bash > /activate.sh \
    && printf '#!/bin/bash\nsource /activate.sh\nexec "$@"\n' > /entrypoint.sh \
    && chmod 0755 /entrypoint.sh /activate.sh \
    && rm -rf ~/.cache/rattler

# ── Runtime: the environment only ────────────────────────────────────────────
FROM ${RUNTIME_IMAGE}
ARG PIXI_ENV
COPY --from=build /app/.pixi/envs/${PIXI_ENV} /app/.pixi/envs/${PIXI_ENV}
COPY --from=build /activate.sh /entrypoint.sh /
ENV PIXI_ENVIRONMENT_NAME=${PIXI_ENV} \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    XDG_CACHE_HOME=/tmp/.cache
LABEL org.opencontainers.image.title="trajgrn-bench-${PIXI_ENV}" \
      trajgrn.pixi_env="${PIXI_ENV}"
WORKDIR /work
ENTRYPOINT ["/entrypoint.sh"]
CMD ["python", "--version"]
