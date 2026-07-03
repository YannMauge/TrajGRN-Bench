# Understanding Results

This page explains how to interpret the output of a benchmark run.

## Output Files

After running a benchmark, you'll find results in
`benchmark/outputs_metrics/<run_name>_<mode>/`:

| File | Description |
|------|-------------|
| `grn_metrics.csv` | GRN inference performance per method |
| `cell_metrics.csv` | Trajectory reconstruction performance per method |
| `ranking_table.png` | Visual ranking of all methods |
| `method_capabilities.csv` | Method feature/capability matrix |

## GRN Metrics

GRN (Gene Regulatory Network) metrics evaluate how well each method recovers
the true regulatory network:

| Metric | Range | Description |
|--------|-------|-------------|
| **AUROC** | [0, 1] | Area Under the ROC Curve — overall edge detection |
| **AUPRC** | [0, 1] | Area Under Precision-Recall Curve — better for imbalanced networks |
| **AUPRC Signed** | [0, 1] | AUPRC considering edge sign (activation vs. inhibition) |
| **Precision@K** | [0, 1] | Precision among top-K predicted edges |

Higher is better for all GRN metrics.

## Cell / Trajectory Metrics

These metrics evaluate how well methods reconstruct cellular trajectories:

| Metric | Range | Description |
|--------|-------|-------------|
| **Wasserstein Distance** | $[0, \infty)$ | Earth Mover's Distance between predicted and true distributions |
| **MMD** | $[0, \infty)$ | Maximum Mean Discrepancy — kernel-based distribution distance |
| **R² (per gene)** | $(-\infty, 1]$ | Gene-wise correlation between predicted and true expression |

Lower is better for distance-based metrics; higher is better for R².

## Gene-Level Metrics

Gene-specific metrics evaluate per-gene prediction quality:

| Metric | Range | Description |
|--------|-------|-------------|
| **Gene-wise R²** | $(-\infty, 1]$ | Per-gene prediction accuracy |
| **Gene-wise MSE** | $[0, \infty)$ | Per-gene mean squared error |

## Ranking Table

The ranking table (`ranking_table.png`) provides a visual summary:

- **Rows**: Methods
- **Columns**: Metric categories (GRN, Cell, Gene, Runtime)
- **Circles**: Proportional to normalized scores — larger/darker = better
- **Number in circle**: Rank position (only top-N shown)
- **Hollow grey circle**: Method does not produce this output type

## Interactive Visualizations

If you ran with visualization enabled, you'll also find:

- `post_analysis/visualizations.py` — Static and interactive plots
- Per-method diagnostic plots in each method's output directory

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
