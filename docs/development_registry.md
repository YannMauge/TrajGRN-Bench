# Methods Registry

The `methods_registry.yaml` file is the **single source of truth** for all
method metadata in TrajGRN-Bench. Every downstream consumer derives its
configuration from this file.

## Why a Registry?

Before the registry, adding a new method required touching 6+ files:

- `benchmark_run.sh` — default method list
- `single_run.sh` — entrypoint parameters
- `ranking_table.py` — capabilities
- `containers/Makefile` — build targets
- `containers/build_all.sh` — build targets
- `configs/benchmark.schema.json` — validation

Now, you only need to edit **one file**: `methods_registry.yaml`.

## Registry Structure

```yaml
schema_version: 1

methods:
  MethodName:
    display_name: Human-Readable Name
    description: >
      Multi-line description of the method.
    github: https://github.com/org/repo
    publication: "Citation string"

    entrypoint:
      type: python               # python | bash
      script: methods/foo/run.py
      pass_output_mode: true     # Pass --output_mode flag
      pass_ko_genes: false       # Pass --ko_genes flag
      supports_perturbation_training: false
      extra_args:                # Always-passed CLI args
        - "--flag"
        - "value"

    execution:
      default_runner: docker     # docker | conda | apptainer
      conda_env: env_name
      container:
        image: benchmark/method:latest
        env_file: environments/method.yml
        env_name: env_name

    capabilities:
      grn_inference: true
      trajectory_reconstruction: true
      perturbation_training: false

iteration_order:
  - MethodOne
  - MethodTwo
  # ... controls display order in tables
```

## Python API

The `utils/methods_registry.py` module provides programmatic access:

```python
from utils.methods_registry import (
    get_registry,
    get_method,
    get_method_keys,
    get_capabilities_table,
    get_execution_config,
)

# Get all registered methods
registry = get_registry()
print(registry["methods"].keys())

# Get a specific method
method = get_method("FLeCS")
print(method["capabilities"])

# Get capabilities as a DataFrame
caps = get_capabilities_table()
print(caps)

# Get execution config (merges registry + user overrides)
exec_cfg = get_execution_config(user_overrides={...})
```

## CLI Reference

For shell scripts (used by `single_run.sh`, `build_all.sh`, etc.):

```bash
# List all method keys
python utils/methods_registry.py keys

# Get entrypoint info for a method
python utils/methods_registry.py entrypoint FLeCS

# Get full execution config as JSON
python utils/methods_registry.py execution-json

# Get capabilities as JSON
python utils/methods_registry.py capabilities-json

# Get default enabled methods JSON
python utils/methods_registry.py default-methods-json

# List container image build specs
python utils/methods_registry.py container-images

# Check if a method supports perturbation training
python utils/methods_registry.py supports-perturbation FLeCS

# Get the conda environment name for a method
python utils/methods_registry.py conda-env FLeCS

# Get the runner mode (docker/conda/apptainer)
python utils/methods_registry.py runner-mode FLeCS

# Get the container image name
python utils/methods_registry.py container-image FLeCS
```

## Auto-Generated Artifacts

When you run `python utils/methods_registry.py` with certain flags, it
generates:

| Consumer | How it uses the registry |
|----------|-------------------------|
| `single_run.sh` | Calls `_registry_query` which invokes registry CLI |
| `benchmark_run.sh` | `python utils/methods_registry.py default-methods-json` |
| `ranking_table.py` | `from utils.methods_registry import get_capabilities_table` |
| `benchmark_run_config.py` | `build_execution_json()` merges registry + overrides |
| `configs/benchmark.schema.json` | `run_methods` propertyNames enum |
| `containers/Makefile` | Targets auto-generated from registry |
| `containers/build_all.sh` | `ALL_IMAGES` auto-generated |
| `docs/_scripts/generate_methods.py` | Auto-generates method pages |
