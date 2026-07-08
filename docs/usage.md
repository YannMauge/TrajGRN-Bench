# Running Benchmarks

This guide covers how to run benchmarks, configure every aspect of a run, and
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

---

## Configuration

TrajGRN-Bench uses a single YAML configuration file to control every aspect of
a benchmark run.

### Schema

The configuration is validated against a JSON schema at
`configs/benchmark.schema.json`. Use an editor with YAML schema support (VS
Code, PyCharm) for auto-completion and validation.

### Minimal Example

```yaml
--8<-- "configs/benchmark.example.yaml"
```

### Configuration Sections

#### `paths`

```yaml
paths:
  input_dir: simulator/custom_network/Data
  adata_dir: benchmark/data
  results_dir: benchmark/outputs_methods
```

#### `simulation`

Controls data generation:

```yaml
simulation:
  mode: simul_replicates    # simul_replicates | simul_ko | false
  replicates_number: 3
  simulator_backend: harissa # harissa | boolode | sergio
  simul_ko_genes: none     # none | all | "GENE1,GENE2"
```

#### `benchmark`

Control which data is used for training and what the methods should output:

```yaml
benchmark:
  train_data: future
  output_mode: full_test
  global_run_name: run_config_test
  future_start_tp: 8
  ko_output_genes: none
  restart_mode: save
  perturbation_training: false
  run_methods:
    FLeCS: false
    FLeCS-TPs: true
    scNODE: true
    reference_fitting: true
    CardamomOT: true
    GENIE3: true
    PEARSON: true
    RENGE: true
    TrajectoryNet: true
    WaddingtonOT: true
```

### Validation

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

---

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
| `full_full` | Simulate all timepoints from t₀ |
| `full_train` | Simulate only training timepoints from t₀ |
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

---

## Output Structure

After a successful run, you'll find:

```
benchmark/
├── outputs_methods/
│   └── <run_name>_<mode>/
│       ├── <method>/
│       │   ├── data_x_GRN.npy       # Inferred GRN adjacency matrix
│       │   ├── data_x_adata.h5ad    # Simulated trajectory (AnnData)
│       │   ├── method_specific/     # Method-specific outputs
│       │   └── .done                # Completion marker
│       └── ...
├── outputs_metrics/
│   └── <run_name>_<mode>/
│       ├── grn_metrics.csv          # GRN inference scores
│       ├── cell_metrics.csv         # Trajectory reconstruction scores
│       ├── ranking_table.png        # Visual ranking chart
│       └── method_capabilities.csv  # Method feature matrix
```

### Output Files (Metrics)

After running a benchmark, you'll find results in
`benchmark/outputs_metrics/<run_name>_<mode>/`:

| File | Description |
|------|-------------|
| `grn_metrics.csv` | GRN inference performance per method |
| `cell_metrics.csv` | Trajectory reconstruction performance per method |
| `ranking_table.png` | Visual ranking of all methods |
| `method_capabilities.csv` | Method feature/capability matrix |

---

## Understanding Metrics

### GRN Metrics

GRN (Gene Regulatory Network) metrics evaluate how well each method recovers
the true regulatory network:

| Metric | Range | Description |
|--------|-------|-------------|
| **AUROC** | [0, 1] | Area Under the ROC Curve — overall edge detection |
| **AUPRC** | [0, 1] | Area Under Precision-Recall Curve — better for imbalanced networks |
| **AUPRC Signed** | [0, 1] | AUPRC considering edge sign (activation vs. inhibition) |
| **Precision@K** | [0, 1] | Precision among top-K predicted edges |

Higher is better for all GRN metrics.

### Cell / Trajectory Metrics

These metrics evaluate how well methods reconstruct cellular trajectories:

| Metric | Range | Description |
|--------|-------|-------------|
| **Wasserstein Distance** | $[0, \infty)$ | Earth Mover's Distance between predicted and true distributions |
| **MMD** | $[0, \infty)$ | Maximum Mean Discrepancy — kernel-based distribution distance |
| **R² (per gene)** | $(-\infty, 1]$ | Gene-wise correlation between predicted and true expression |

Lower is better for distance-based metrics; higher is better for R².

### Gene-Level Metrics

Gene-specific metrics evaluate per-gene prediction quality:

| Metric | Range | Description |
|--------|-------|-------------|
| **Gene-wise R²** | $(-\infty, 1]$ | Per-gene prediction accuracy |
| **Gene-wise MSE** | $[0, \infty)$ | Per-gene mean squared error |

### Ranking Table

The ranking table (`ranking_table.png`) provides a visual summary:

- **Rows**: Methods
- **Columns**: Metric categories (GRN, Cell, Gene, Runtime)
- **Circles**: Proportional to normalized scores — larger/darker = better
- **Number in circle**: Rank position (only top-N shown)
- **Hollow grey circle**: Method does not produce this output type

### Interactive Visualizations

If you ran with visualization enabled, you'll also find:

- `post_analysis/visualizations.py` — Static and interactive plots
- Per-method diagnostic plots in each method's output directory

---

## Exporting for the Docs

To include your results in the documentation (see [Results Gallery](results.md)):

```bash
# Export results to the docs data directory
python docs/_scripts/export_results.py \
  --metrics-dir benchmark/outputs_metrics/my_run_future_full_test \
  --output-dir docs/results_data/my_run
```

Then rebuild the docs. The results will appear automatically in the gallery.

## Comparing Runs

To compare multiple benchmark runs:

```bash
python post_analysis/compute_metrics.py \
  --runs run_1 run_2 run_3 \
  --output comparison.csv
```

This produces a combined CSV with all methods across all runs.

## Next Steps

- [Browse available methods](methods.md)
- [Explore the results gallery](results.md)
- [Add a new method](development_add_method.md)
