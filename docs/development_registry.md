# Methods Registry

`methods_registry.yaml` is the single source of truth for method metadata:
entrypoint, CLI flags, execution environment (pixi env / Docker image) and
capabilities. The shell scripts, the config builder, the image build tooling
and the ranking table all read it.

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
      pass_output_mode: true     # pass -u <output_mode>
      pass_ko_genes: false       # pass -k <ko_genes>
      supports_perturbation_training: false   # pass -p; adds a _perturbation_training variant
      extra_args:                # always-passed CLI args
        - "--flag"
        - "value"

    execution:
      default_runner: docker     # docker | apptainer | pixi | conda (legacy)
      pixi_env: foo              # environment in pixi.toml (local runs and image)
      container:
        image: benchmark/foo:latest   # built from pixi env "foo"

    capabilities:
      grn_inference: true
      trajectory_reconstruction: true
      perturbation_training: false

defaults:
  execution:
    runner: pixi
    pixi_env: baselines
    container_engine: docker
    container_image: benchmark/baselines:latest

iteration_order:
  - MethodOne
  - MethodTwo
```

### Runner selection

For each method, the runner is (highest priority first): the config's
`execution.method_runners[method]`, the config's `execution.default_runner`,
the registry's `execution.default_runner`, the registry default.
`execution.container_engine: apptainer` replaces Docker by Apptainer, and
`execution.gpu: true` uses the `<pixi_env>-cuda` environment / image when
`pixi.toml` defines one.

## Python API

```python
from utils.methods_registry import (
    load_registry, list_method_keys, get_iteration_order,
    get_method, get_entrypoint, get_execution, get_capabilities,
    get_capabilities_table, export_execution_json,
    list_container_images, pixi_environments, cuda_pixi_envs,
)

print(get_execution("FLeCS")["pixi_env"])       # "flecs"
print(get_capabilities_table()["GENIE3"])        # {"grn_inference": True, ...}
print(list_container_images()[0])                # {"target": ..., "image": ..., "platforms": [...]}
```

## CLI Reference

```bash
python utils/methods_registry.py keys                      # methods in iteration order
python utils/methods_registry.py entrypoint FLeCS          # type|script|pass_out|pass_ko|pert|extra
python utils/methods_registry.py supports-perturbation FLeCS
python utils/methods_registry.py runner-mode FLeCS         # docker | apptainer | pixi | conda
python utils/methods_registry.py pixi-env FLeCS            # flecs
python utils/methods_registry.py pixi-envs                 # environments declared in pixi.toml
python utils/methods_registry.py container-image FLeCS     # benchmark/flecs:latest
python utils/methods_registry.py container-images          # image build specs (JSON)
python utils/methods_registry.py execution-json            # registry execution settings (JSON)
python utils/methods_registry.py capabilities-json
python utils/methods_registry.py default-methods-json
```

## Consumers

| Consumer | How it uses the registry |
|----------|-------------------------|
| `single_run.sh` | entrypoint, runner, pixi env and image per method (`_registry_query`, `BENCHMARK_EXECUTION_JSON`) |
| `benchmark_run.sh` | default enabled methods |
| `benchmark_run_config.py` | `build_execution_json()` merges registry + config overrides |
| `ranking_table.py` | capabilities (which categories count for each method) |
| `containers/build_all.sh`, `containers/Makefile` | one image per `pixi_env` (+ `-cuda` variants) |
| `docs/_scripts/generate_methods.py` | method pages |
