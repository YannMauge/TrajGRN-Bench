# Containers for method execution

This benchmark supports per-method containers.
Build one image per method and point the config `execution.container_images` to each image.

**Image definitions are in `methods_registry.yaml`.**  Run `python utils/methods_registry.py container-images` to list them.

## Build all images at once

### Using Make (recommended)

```bash
# Build all images
make -f containers/Makefile all

# Build in parallel (e.g. 4 concurrent builds)
make -f containers/Makefile -j4 all

# Build a single image (target name from registry, e.g. flecs, scnode, ...)
make -f containers/Makefile flecs

# List images / show help
make -f containers/Makefile list
make -f containers/Makefile help

# Remove all built images
make -f containers/Makefile clean
```

### Using the shell script

```bash
# Build all images
bash containers/build_all.sh

# Build only specific images (target names from registry)
bash containers/build_all.sh flecs scnode

# Force full rebuild (no layer cache)
bash containers/build_all.sh --no-cache

# Build and push to registry
bash containers/build_all.sh --push

# Parallel builds (4 at a time)
bash containers/build_all.sh -j 4
```

## Docker build examples (one at a time)

Image definitions come from `methods_registry.yaml`. Use the registry to discover the exact build command:

```bash
python utils/methods_registry.py container-images
```

Example output (truncated):
```json
[
  {"target": "flecs", "env_file": "environments/flecs_cpu.yml", "env_name": "flecs", "image": "benchmark/flecs"},
  {"target": "scnode", "env_file": "environments/sc_dynamic.yml", "env_name": "sc_dynamic", "image": "benchmark/scnode"},
  ...
]
```

To build manually from the repository root:
```bash
docker build -f containers/docker/Dockerfile --build-arg ENV_FILE=environments/flecs_cpu.yml --build-arg ENV_NAME=flecs -t benchmark/flecs:latest .
```

## Apptainer

Apptainer can consume the Docker images directly:

```bash
apptainer build scnode.sif docker-daemon://benchmark/scnode:latest
apptainer build flecs.sif docker-daemon://benchmark/flecs:latest
```

Then set the config `execution.method_runners` values to `apptainer` and point
`execution.container_images` at the `.sif` files.
