# Running Benchmarks

This guide covers how to run benchmarks, interpret the configuration, and
understand the output.

## Quick Run

The simplest way to run a benchmark is through the shell script:

```bash
bash benchmark_run.sh --config configs/benchmark.example.yaml
```

This single command will:

1. Load the specified dataset(s)
2. Run each enabled method on the training data
3. Compute GRN and trajectory metrics
4. Generate ranking tables and visualizations

## Run Modes

TrajGRN-Bench supports several training/evaluation modes, controlled by
`train_data` and `output_mode`:

### Training Modes (`train_data`)

| Mode | Description |
|------|-------------|
| `full` | Train on all available timepoints |
| `future` | Train on only the first $N$ timepoints, predict the rest |
| `leave-one-out` | Leave one intermediate timepoint out; run all combinations |
| `subsample_full` | All timepoints but only 66% of cells at each |

### Output Modes (`output_mode`)

| Mode | Description |
|------|-------------|
| `full_full` | Simulate all timepoints from $t_0$ |
| `full_train` | Simulate only training timepoints from $t_0$ |
| `full_test` | Simulate test timepoints from last training timepoint |
| `no_traj` | Skip trajectory simulation (GRN-only evaluation) |

### Simulation Modes

| Mode | Description |
|------|-------------|
| `simul_replicates` | Generate replicate simulations |
| `simul_ko` | Simulate gene knockouts |
| `false` | Use existing data (no simulation) |

## Selecting Methods

Methods are selected via the `run_methods` section in your config:

```yaml
run_methods:
  FLeCS: "1"
  scNODE: "1"
  CardamomOT: "1"
  GENIE3: "1"
  PEARSON: "1"
```

Set a method to `"0"` to disable it, or omit it entirely.

To see all available methods:

```bash
python utils/methods_registry.py keys
```

## Output Structure

After a successful run, you'll find:

```
benchmark/
├── outputs_methods/
│   └── <run_name>_<mode>/
│       ├── <method>/
│       │   ├── predicted_grn.csv       # Inferred GRN adjacency matrix
│       │   ├── simulated_adata.h5ad    # Simulated trajectory (AnnData)
│       │   ├── latent_representation/  # Method-specific latent outputs
│       │   └── .done                    # Completion marker
│       └── ...
├── outputs_metrics/
│   └── <run_name>_<mode>/
│       ├── grn_metrics.csv             # GRN inference scores
│       ├── cell_metrics.csv            # Trajectory reconstruction scores
│       ├── ranking_table.png           # Visual ranking chart
│       └── method_capabilities.csv     # Method feature matrix
```

## Running with Python Directly

If you prefer Python over the shell wrapper:

```bash
conda run -n benchmark_runner python benchmark_run_config.py --config configs/benchmark.example.yaml
```

## Restart / Resume

By default, the benchmark skips methods that have a `.done` marker:

```yaml
execution:
  restart_mode: save    # Skip completed methods
  # restart_mode: rerun  # Re-run everything from scratch
```

## Next Steps

- [Customize your configuration](usage_config.md)
- [Interpret your results](usage_results.md)
- [Add a new method](development_add_method.md)
