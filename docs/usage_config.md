# Configuration Guide

TrajGRN-Bench uses a single YAML configuration file to control every aspect of
a benchmark run. This page documents all available options.

## Schema

The configuration is validated against a JSON schema at
`configs/benchmark.schema.json`. Use an editor with YAML schema support (VS
Code, PyCharm) for auto-completion and validation.

## Minimal Example

```yaml
--8<-- "configs/benchmark.example.yaml"
```

## Configuration Sections

### `run_name`

A human-readable name for this benchmark run. Used in output directory names.

```yaml
run_name: my_experiment
```

### `train_data` and `output_mode`

Control which data is used for training and what the methods should output:

```yaml
train_data: future          # full | future | leave-one-out | subsample_full
output_mode: full_test      # full_full | full_train | full_test | no_traj
```

### `simulation`

Controls data generation:

```yaml
simulation:
  mode: simul_replicates    # simul_replicates | simul_ko | false
  simulator_backend: harissa # harissa | boolode | sergio
  replicates_number: 3
  ko_output_genes: none     # none | all | "GENE1,GENE2"
```

### `data`

Input data configuration:

```yaml
data:
  adata_files:              # Paths to input AnnData files
    - benchmark/data/data_1.h5ad
  train_tps: benchmark/data/train_tps.npy
  test_tps: benchmark/data/test_tps.npy
  degradation_rates: null   # Optional path to degradation rates file
```

### `run_methods`

Which methods to run (`"1"` = enabled, `"0"` = disabled):

```yaml
run_methods:
  FLeCS: "1"
  FLeCS-TPs: "0"
  scNODE: "1"
  CardamomOT: "1"
  reference_fitting: "1"
  GENIE3: "1"
  PEARSON: "1"
  WOT: "0"
  TrajectoryNet: "0"
  RENGE: "0"
```

!!! tip "Discover available methods"
    ```bash
    python utils/methods_registry.py keys
    ```

### `execution`

Controls how methods are executed:

```yaml
execution:
  default_runner: docker              # docker | conda | apptainer
  restart_mode: save                  # save (skip completed) | rerun
  conda_envs: {}                      # Only for conda runner
  container_images: {}                # Only for docker/apptainer runner
```

### `metrics`

Which metrics to compute:

```yaml
metrics:
  grn: true           # Compute GRN metrics
  cell: true           # Compute cell-level trajectory metrics
  gene: true           # Compute gene-level metrics
```

### `output`

Output paths:

```yaml
output:
  methods_dir: benchmark/outputs_methods
  metrics_dir: benchmark/outputs_metrics
```

### `future`

Settings for future-prediction mode:

```yaml
future:
  start_tp: 2        # Timepoint index where test period begins
```

### `perturbation`

Perturbation training settings:

```yaml
perturbation:
  training: false     # Enable perturbation-training variants
```

## Advanced: Method-Specific Overrides

You can override execution settings per method:

```yaml
execution:
  default_runner: docker
  method_runners:
    GENIE3: conda          # Run GENIE3 in conda, everything else in Docker
  container_images:
    flecs: benchmark/flecs:latest
  conda_envs:
    GENIE3: genie3_env
```

## Validation

Validate your config against the schema:

```bash
python -c "
import yaml, jsonschema
with open('configs/benchmark.example.yaml') as f:
    cfg = yaml.safe_load(f)
with open('configs/benchmark.schema.json') as f:
    schema = json.load(f)
jsonschema.validate(cfg, schema)
print('Config is valid!')
"
```
