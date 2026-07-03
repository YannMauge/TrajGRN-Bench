#!/usr/bin/env python3
"""
post_analysis — Benchmark post-analysis package.

Provides modules for computing and visualizing GRN and cell-population metrics
on benchmark run outputs.

Modules:
    _common         Shared constants, regex patterns, and utilities.
    metrics_grn     GRN evaluation metrics (precision, recall, F1, AUPRC, …).
    metrics_cells   Cell-population metrics (Wasserstein, KL, JS, ARI, …).
    visualizations  Plotting functions (GRN heatmaps, cell UMAPs, metrics charts).
    compute_metrics Orchestrator that runs metrics + visualizations for a full run.
"""

from ._common import (
    ADATA_FILE_RE,
    GRN_FILE_RE,
    GRN_TRUE_FOLDER,
    PERTURBATION_TRAINING_SUFFIX,
    TRUE_ADATA_RE,
    discover_replicates,
    flatten_metrics,
    parse_ko_output_genes,
    setup_logging,
)
from .metrics_cells import (
    CellMetricsError,
    compute_cell_metrics,
    compute_wasserstein,
    load_adata,
    read_adata,
)
from .metrics_grn import (
    GRNMetricsError,
    compute_grn_metrics,
    load_grn,
)
from .visualizations import (
    plot_cell_metrics,
    plot_cell_umaps,
    plot_grn_heatmaps,
    plot_grn_metrics,
)

__all__ = [
    # _common
    "ADATA_FILE_RE",
    "GRN_FILE_RE",
    "GRN_TRUE_FOLDER",
    "PERTURBATION_TRAINING_SUFFIX",
    "TRUE_ADATA_RE",
    "discover_replicates",
    "flatten_metrics",
    "parse_ko_output_genes",
    "setup_logging",
    # metrics_cells
    "CellMetricsError",
    "compute_cell_metrics",
    "compute_wasserstein",
    "load_adata",
    "read_adata",
    # metrics_grn
    "GRNMetricsError",
    "compute_grn_metrics",
    "load_grn",
    # visualizations
    "plot_cell_metrics",
    "plot_cell_umaps",
    "plot_grn_heatmaps",
    "plot_grn_metrics",
]
