# Containers

One image per pixi environment used by the methods (`execution.pixi_env` in
`methods_registry.yaml`), built from `pixi.toml` / `pixi.lock` by
`containers/docker/Dockerfile`:

- build stage: `ghcr.io/prefix-dev/pixi` (pinned by digest) runs
  `pixi install --locked -e <env>`;
- runtime stage: `ubuntu:24.04` (pinned by digest) with that environment only,
  activated by `/entrypoint.sh`.

```bash
make -f containers/Makefile list         # targets
make -f containers/Makefile all          # CPU images for the host platform
make -f containers/Makefile all-cuda     # + <env>-cuda images (linux/amd64)
bash containers/build_all.sh --platform linux/amd64 --push flecs   # e.g. for a registry
```

`build_all.sh` works with macOS `/bin/bash` 3.2, builds with `docker buildx`,
prints a summary and exits non-zero if any build failed (logs in a temp dir).

The repository is **not** copied into the images: the benchmark mounts it (and
the data / results folders) at the same paths as on the host and runs the
container as the host user. GPU runs (`execution.gpu: true`) use the
`benchmark/<env>-cuda` images with `docker run --gpus all` (or
`apptainer exec --nv`).

Apptainer: convert an image with
`apptainer build flecs.sif docker-daemon://benchmark/flecs:latest` and set
`execution.container_engine: apptainer` and `container_images` to the `.sif`
paths.
