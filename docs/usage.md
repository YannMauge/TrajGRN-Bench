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
  simulator_backend: harissa # harissa | boolode | sergio | dyngen
  simul_ko_genes: none     # training knockouts: none | all | "GENE1,GENE2"
  n_cells: 1000            # optional, cells per dataset
  n_time_bins: 10          # optional, timepoints labelled 0, 10, ..., 90
  seed: 0                  # optional, replicate r uses seed*1000 + r
```

All simulators produce `n_time_bins` timepoints with the same labels
(0, 10, …, 90 for 10 bins) and re-simulate knockouts (the knocked-out gene's
production is set to 0, so its targets respond):

| Backend | Model | Notes |
|---------|-------|-------|
| `harissa` | Harissa bursty PDMP | shared 8-gene network (with a mutual-inhibition cycle) |
| `boolode` | BoolODE (BEELINE), `hill` model | shared network; one cell per simulated trajectory |
| `sergio` | upstream SERGIO, differentiation mode | **its own acyclic 8-gene network** (SERGIO cannot simulate cycles or self-loops); a linear chain of states as pseudo-time |
| `dyngen` | dyngen (R) | shared network; time bins relabelled to the shared labels |

Degradation-rate priors: Harissa and BoolODE write their own mRNA / protein
degradation rates (converted to timepoint units) to
`Rates/degradation_rates.txt`, stored in `adata.var['d0' / 'd1']` and used by
CardamomOT. SERGIO (pseudo-time, no protein) and dyngen (progression bins,
aggregated genes) do not define such rates: no priors are given, and
CardamomOT falls back to its half-life table.

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
    dynGENIE3: true
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
| `full` | Train on all available timepoints (no test timepoint: use `full_train` / `full_full`) |
| `future` | Train on the first `future_start_tp` timepoints (index, 2 ≤ value < n), predict the rest |
| `leave-one-out` | Leave one intermediate timepoint out; run all combinations (≥ 3 timepoints) |
| `subsample_full` | Train on 66% of the cells of every timepoint; simulate from the held-out t₀ cells and evaluate against the held-out 34% at every later timepoint (`full_full` or `full_train`) |

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

## Knockouts: training vs prediction

Two settings control knockouts (KOs):

- `simulation.simul_ko_genes`: the **training** KOs. With
  `simulation: simul_ko`, their cells are given to the methods
  (`dataset_id` = gene name). Only `_perturbation_training` variants use them;
  standard variants train on wild-type cells.
- `benchmark.ko_output_genes`: the KOs methods must **predict**
  (`data_x_ko_<gene>_adata.h5ad`).

With `simulation: simul_ko`, both sets are simulated. Methods read their data
from `<run>/inputs/`, which contains wild-type cells plus the training KOs only:
a KO that is only in `ko_output_genes` is held out and used solely as ground
truth, so predicting it is a test on an unseen perturbation. A gene listed in
both settings is allowed, with a warning. With `simulation: false`, the same
rule applies to the `dataset_id` groups of your data.

KO predictions are scored against the simulated KO cells, on the same
evaluated timepoints as wild-type predictions.

## Selecting Methods

Methods are selected via the `run_methods` section in your config:

```yaml
run_methods:
  FLeCS: "1"
  scNODE: "1"
  CardamomOT: "1"
  GENIE3: "1"
  dynGENIE3: "1"
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
pixi run -e runner python benchmark_run_config.py --config configs/benchmark.example.yaml
```

The orchestration steps run in the pixi environment `runner` by default. To
use a conda environment instead (e.g. on a machine without pixi), set
`execution.runner_tool: conda` and `execution.runner_env` in the config, and
export `BENCHMARK_RUNNER_TOOL=conda BENCHMARK_RUNNER_ENV=<env>` so that the
config itself is also parsed with it.

## Restart / Resume

When a method finishes on a data file, `single_run.sh` writes a `.done` marker
containing a fingerprint of everything that determines the result: the data
file content, the train/test split files (`train_tps.npy`, `test_tps.npy`,
`subsample_train_ids.npy`), `output_mode`, `ko_output_genes`, the
perturbation-training flag and the method's registry entrypoint.

```yaml
benchmark:
  restart_mode: save    # Skip methods whose .done fingerprint still matches
  # restart_mode: rerun  # Re-run everything
```

- With `save`, a method is skipped only if its fingerprint is unchanged. If the
  data or settings changed, the differences are logged and the method re-runs.
- Before a method (re-)runs on a data file, its previous outputs for that file
  (`data_X_*`, `project_data_X/`) are removed.
- Changes to a method's source code are not detected: use `rerun` after
  editing a method.
- Markers written by older versions (empty files) are trusted, with a warning.

### Failed methods

A failing method does not stop the run. The remaining methods, data files and
leave-one-out folds still run, metrics and rankings are computed for the methods
that succeeded, and a summary of the failures is printed at the end. Failures are
also listed in `failed_methods.tsv` in each run folder. The run then exits with
code `3`. Re-running with `restart_mode: save` only re-runs the failed methods.
Interrupting with Ctrl-C still stops the whole run.

---

## Output Structure

After a successful run, you'll find:

```
benchmark/
├── outputs_methods/
│   └── <run_name>_<mode>/
│       ├── inputs/                  # What methods read (rebuilt each run)
│       │   ├── data_x.h5ad          # WT + training knockouts only
│       │   ├── train_tps.npy        # Train / test split of this run (or fold)
│       │   └── test_tps.npy
│       ├── train_tps.npy, test_tps.npy   # Same split, read by the metrics
│       ├── GRN_true/                # Ground-truth networks
│       ├── failed_methods.tsv       # Only if some methods failed
│       ├── <method>/
│       │   ├── data_x_GRN.npy       # Inferred GRN adjacency matrix
│       │   ├── data_x_adata.h5ad    # Simulated trajectory (AnnData)
│       │   ├── method_specific/     # Method-specific outputs
│       │   └── data_x.done          # Completion marker (input fingerprint)
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
| `grn_metrics.csv` | GRN metrics per method / replicate (both edge sets) |
| `cell_metrics.csv` | Cell metrics per method / replicate, averaged over evaluated timepoints |
| `cell_metrics_per_timepoint.csv` | Cell metrics for each evaluated timepoint and scale |
| `ranking_table.png` | Visual ranking of all methods (wild-type predictions) |
| `ranking_table_ko_prediction.png` | Ranking of knockout predictions (methods that predict KOs) |
| `method_capabilities.csv` | Method feature/capability matrix |

---

## Understanding Metrics

### GRN Metrics

GRN (Gene Regulatory Network) metrics evaluate how well each method recovers
the true regulatory network. They are **threshold-free**: they score the
ranking of the predicted edge weights, so methods whose scores live on
different scales are compared fairly, and an edge is never lost because it falls
just below a cut-off.

The ground truth `GRN_true/data_X_GRN.npy` holds signed interaction strengths
(`> 0` activation, `< 0` inhibition, `0` no edge). Any real-valued truth works:
an edge exists where the value is non-zero and its sign is the sign of that
value. A truth with no negative entry (e.g. older 0/1 files) is treated as
unsigned, and the sign-aware metrics are then `NaN`.

Each metric is computed on two edge sets (column prefix), always excluding
self-loops:

- `genes_…`: gene → gene edges only;
- `full_…`: gene → gene plus Stimulus → gene edges.

| Metric | Range | Description |
|--------|-------|-------------|
| `edge_auprc` | [0, 1] | Average precision of \|prediction\| for edge existence |
| `edge_auprc_ratio` | $[0, \infty)$ | `edge_auprc` divided by the edge density (1 = random) |
| `edge_auroc` | [0, 1] | ROC AUC of \|prediction\| for edge existence (0.5 = random) |
| `signed_auprc` | [0, 1] | Mean of the AUPRC of activations ranked by +prediction and of inhibitions ranked by −prediction |
| `signed_auprc_ratio` | $[0, \infty)$ | Same, each part divided by its class density |
| `sign_accuracy` | [0, 1] | Fraction of true edges whose predicted sign is correct |
| `spearman` | [−1, 1] | Spearman correlation of signed prediction vs signed truth (uses strengths) |

Only each method's common network (`data_X_GRN.npy`) is scored: the network
learned from wild-type data and, for `_perturbation_training` variants, from the
training perturbations. Knockout-specific network files
(`data_X_ko_<gene>_GRN.npy`) describe a method's model under a knockout, not a
network to recover, and are not scored.

Higher is better for all GRN metrics. Methods that only output non-negative
scores (e.g. importances) can find edges but cannot recover inhibitions, so
they score lower on the signed metrics.

### Cell / Trajectory Metrics

Predicted cells are compared with true cells **separately at each test
timepoint**, and the values are averaged. With `subsample_full`, the true
cells are the held-out ones only (`eval_set` is `heldout`). If a run predicts no test timepoint
(`output_mode: full_train`, or an empty test set), the predicted training
timepoints except the first one are used instead, and `eval_set` is `train`.

Each metric is computed on log1p-transformed expression (`…_log`, used in the
ranking) and on raw expression (`…_raw`):

| Metric | Range | Description |
|--------|-------|-------------|
| `wasserstein` | $[0, \infty)$ | Exact earth mover's distance in a PCA space fitted on the true cells |
| `js_divergence` | $[0, \ln 2]$ | Jensen–Shannon divergence between Gaussian fits in the same PCA space |
| `mean_expr_rmse` | $[0, \infty)$ | RMSE between per-gene mean expressions |
| `gene_corr_similarity` | [−1, 1] | Correlation between the predicted and true gene–gene correlation matrices |

Lower is better, except `gene_corr_similarity` (higher is better).
Knockout predictions (`data_X_ko_<gene>_adata.h5ad`) are compared with the
simulated knockout cells (`obs['dataset_id'] == <gene>` in the true data). A
knockout that was not simulated has no ground truth, so its prediction is not
scored. If a gene is both a training knockout and in `ko_output_genes`, a
warning is printed: perturbation-training variants have seen that knockout.

`n_missing_timepoints` counts evaluated timepoints where the method produced
fewer than 2 cells. Such a timepoint receives the worst value reached by any
method at that timepoint before averaging (`n_imputed_timepoints`; flagged as
`imputed_worst` in `cell_metrics_per_timepoint.csv`).

### Ranking Table

The ranking table (`ranking_table.png`) provides a visual summary:

- **Rows**: Methods
- **Rows used**: the main tables use wild-type results only, so every method is
  compared on the same task; `ranking_table_ko_prediction.png` ranks the
  knockout predictions (cell metrics averaged over knockout genes).
- **Columns**: a non-redundant subset of the metrics: GRN (`genes_edge_auprc_ratio`,
  `genes_signed_auprc_ratio`, `genes_edge_auroc`), Cell (`wasserstein_log`,
  `js_divergence_log`), Gene (`mean_expr_rmse_log`, `gene_corr_similarity_log`)
  and Runtime. Every other metric stays in the CSV files.
- **Circles**: Proportional to normalized scores — larger/darker = better
- **Number in circle**: Rank position (only top-N shown)
- **Hollow grey circle**: the method does not produce this output (per
  `methods_registry.yaml` capabilities)
- **Red circle with a cross**: the method should produce this output but no
  result was found (e.g. it failed)
- **Overall rank**: within each category (GRN, Cell, Gene), methods are ranked and
  the rank is turned into a score (best = 1, last = 1/n). The overall score is
  the mean over the categories the method can produce, according to the
  capabilities in `methods_registry.yaml` (Cell and Gene both need
  `trajectory_reconstruction`). A method that should produce a category but has
  no result there (e.g. it failed) scores 0 in it; categories it cannot produce
  are ignored.

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
