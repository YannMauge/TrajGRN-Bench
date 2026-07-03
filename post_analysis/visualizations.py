#!/usr/bin/env python3
"""
visualizations.py — Visualization helpers for GRN and cell-population outputs.

Provides:
  - GRN heatmaps (predicted vs true vs difference)
  - Cell UMAP scatter plots (coloured by timepoint)
  - Bar charts and heatmaps summarising metrics CSVs

Standalone usage:
    python visualizations.py grn  <global_run_name> <run_type> [options]
    python visualizations.py cell <global_run_name> <run_type> [options]
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: E402
import warnings  # noqa: E402

for module in ["anndata", "scipy"]:
    warnings.filterwarnings("ignore", module=module)

try:
    from ._common import (  # noqa: E402
        ADATA_FILE_RE,
        GRN_FILE_RE,
        GRN_TRUE_FOLDER,
        PERTURBATION_TRAINING_SUFFIX,
        TRUE_ADATA_RE,
        discover_replicates,
        iter_method_dirs,
        parse_ko_output_genes,
        resolve_benchmark_path,
        setup_logging,
    )
    from .metrics_cells import read_adata  # noqa: E402
except ImportError:
    from _common import (  # noqa: E402  # type: ignore[no-redef]
        ADATA_FILE_RE,
        GRN_FILE_RE,
        GRN_TRUE_FOLDER,
        PERTURBATION_TRAINING_SUFFIX,
        TRUE_ADATA_RE,
        discover_replicates,
        iter_method_dirs,
        parse_ko_output_genes,
        resolve_benchmark_path,
        setup_logging,
    )
    from metrics_cells import read_adata  # noqa: E402  # type: ignore[no-redef]

logger = setup_logging(__name__)

# ============================================================================
# Constants
# ============================================================================

_DEFAULT_PCA_COMPONENTS: int = 50
_UMAP_FIGSIZE: Tuple[int, int] = (12, 5)
_UMAP_DOT_SIZE: int = 25
_UMAP_DOT_ALPHA: float = 0.6
_HEATMAP_FIGSIZE: Tuple[int, int] = (14, 5)

# Metric definitions: (column_name, display_label)
_GRN_METRICS: List[Tuple[str, str]] = [
    ("precision", "Precision"),
    ("recall", "Recall"),
    ("f1_score", "F1 Score"),
    ("auprc", "AUPRC"),
    ("jaccard_index", "Jaccard Index"),
    ("spearman_correlation", "Spearman Corr."),
    ("hamming_distance", "Hamming Distance"),
    ("edge_overlap_edge_overlap_ratio", "Edge Overlap Ratio"),
]

_CELL_METRICS: List[Tuple[str, str]] = [
    ("wasserstein_distance", "Wasserstein Distance"),
    ("pca_wasserstein_distance", "PCA Wasserstein Distance"),
    ("kl_divergence", "KL Divergence"),
    ("jensen_shannon_divergence", "Jensen-Shannon Divergence"),
    ("gene_gene_correlation", "Gene-Gene Correlation"),
    ("top_de_genes_overlap", "Top DE Genes Overlap"),
    ("adjusted_rand_index", "Adjusted Rand Index"),
    ("runtime", "Runtime (s)"),
]


# ============================================================================
# Helpers
# ============================================================================


def signed_log1p(x: np.ndarray) -> np.ndarray:
    """Sign-preserving log1p transform: sign(x) * log(1 + |x|)."""
    return np.sign(x) * np.log1p(np.abs(x))


def _load_grn(path: Path) -> np.ndarray:
    """Load a .npy GRN file and squeeze trailing size-1 dimensions to 2-D."""
    matrix = np.load(str(path)).astype(float)
    while matrix.ndim > 2 and matrix.shape[-1] == 1:
        matrix = matrix[..., 0]
    if matrix.ndim != 2:
        raise ValueError(
            f"GRN '{path}' is not 2-D after squeezing (shape={matrix.shape})"
        )
    return matrix


def _heatmap_panel(ax: plt.Axes, matrix: np.ndarray, title: str, fig: plt.Figure) -> None:
    """Render a single signed-log1p heatmap panel on *ax*."""
    matrix_heatmap = signed_log1p(matrix)
    extremum = float(np.max(np.abs(matrix_heatmap)) or 1.0)
    im = ax.imshow(
        matrix_heatmap, aspect="equal", cmap="RdYlGn",
        vmin=-extremum, vmax=extremum,
    )
    ax.set_box_aspect(1)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Target gene", fontsize=8)
    ax.set_ylabel("Source gene", fontsize=8)
    ax.tick_params(labelsize=8)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.05)
    fig.colorbar(im, cax=cax)


def _resolve_method_name(method: str, perturbation_training: bool) -> str:
    """Append the perturbation-training suffix when needed."""
    if not perturbation_training or method.endswith(PERTURBATION_TRAINING_SUFFIX):
        return method
    return f"{method}{PERTURBATION_TRAINING_SUFFIX}"


# ============================================================================
# GRN heatmaps
# ============================================================================


def plot_grn_heatmaps(run_dir: Path, output_dir: Path) -> List[Path]:
    """Side-by-side heatmaps (predicted | true | difference) for every
    (method, replicate) pair with a matching ground-truth GRN.

    Args:
        run_dir: Path to ``benchmark/outputs_methods/<run_name>/``.
        output_dir: Directory for saved PNG files.

    Returns:
        List of saved figure paths.
    """
    grn_true_dir = run_dir / GRN_TRUE_FOLDER
    if not grn_true_dir.is_dir():
        logger.info("No GRN_true folder in '%s'; skipping heatmaps.", run_dir)
        return []

    true_grns = discover_replicates(grn_true_dir, GRN_FILE_RE)
    if not true_grns:
        logger.info("No ground-truth GRN files in '%s'; skipping heatmaps.", grn_true_dir)
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []

    for method, method_dir, fold in iter_method_dirs(run_dir):

        pred_grns = discover_replicates(method_dir, GRN_FILE_RE)

        if not pred_grns:
            logger.info("Method '%s': no GRN files, skipping heatmaps.", method)
            continue

        for (replicate, perturbation, pt), pred_path in sorted(pred_grns.items()):
            true_key = (
                (replicate, perturbation, False)
                if (replicate, perturbation, False) in true_grns
                else (replicate, None, False)
            )
            if true_key not in true_grns:
                logger.warning(
                    "Method '%s', replicate %d: no matching true GRN, skipping.",
                    method, replicate,
                )
                continue

            try:
                pred_matrix = _load_grn(pred_path)
                true_matrix = _load_grn(true_grns[true_key])
            except Exception as exc:
                logger.error(
                    "Could not load GRN for method='%s', rep=%d: %s",
                    method, replicate, exc,
                )
                continue

            method_name = _resolve_method_name(method, pt)

            fig, axes = plt.subplots(1, 3, figsize=_HEATMAP_FIGSIZE)
            suptitle = f"GRN heatmaps  —  method: {method_name},  replicate: {replicate}"
            if fold is not None:
                suptitle += f",  fold: {fold}"
            fig.suptitle(suptitle, fontsize=11, y=1.01)

            _heatmap_panel(axes[0], pred_matrix, "Predicted GRN", fig)
            _heatmap_panel(axes[1], true_matrix, "True GRN", fig)

            min_n = min(pred_matrix.shape[0], true_matrix.shape[0])
            min_m = min(pred_matrix.shape[1], true_matrix.shape[1])
            _heatmap_panel(
                axes[2],
                pred_matrix[:min_n, :min_m] - true_matrix[:min_n, :min_m],
                "Difference (Pred − True)", fig,
            )

            fig.tight_layout()
            if fold is not None:
                out_path = output_dir / f"grn_heatmap_{method_name}_fold{fold}_rep{replicate}.png"
            else:
                out_path = output_dir / f"grn_heatmap_{method_name}_rep{replicate}.png"
            fig.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            saved.append(out_path)
            logger.info("Saved heatmap: %s", out_path)

    return saved


# ============================================================================
# Precision-Recall curves
# ============================================================================


def plot_grn_pr_curves(
    run_dir: Path,
    output_dir: Path,
    threshold: float = 0.0,
) -> List[Path]:
    """Precision-Recall curves for every (method, replicate) pair with a
    matching ground-truth GRN.

    Each plot shows the PR curve and includes the AUPRC value in the legend.
    The diagonal dashed line represents the baseline (random classifier)
    which equals the fraction of positive edges.

    Args:
        run_dir:   Path to ``benchmark/outputs_methods/<run_name>/``.
        output_dir: Directory for saved PNG files.
        threshold: Binarization threshold for the true GRN.

    Returns:
        List of saved figure paths.
    """
    try:
        from .metrics_grn import compute_precision_recall_curve
    except ImportError:
        from metrics_grn import compute_precision_recall_curve  # type: ignore[no-redef]

    grn_true_dir = run_dir / GRN_TRUE_FOLDER
    if not grn_true_dir.is_dir():
        logger.info("No GRN_true folder in '%s'; skipping PR curves.", run_dir)
        return []

    true_grns = discover_replicates(grn_true_dir, GRN_FILE_RE)
    if not true_grns:
        logger.info("No ground-truth GRN files in '%s'; skipping PR curves.", grn_true_dir)
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []

    for method, method_dir, fold in iter_method_dirs(run_dir):

        pred_grns = discover_replicates(method_dir, GRN_FILE_RE)

        if not pred_grns:
            logger.info("Method '%s': no GRN files, skipping PR curves.", method)
            continue

        for (replicate, perturbation, pt), pred_path in sorted(pred_grns.items()):
            true_key = (
                (replicate, perturbation, False)
                if (replicate, perturbation, False) in true_grns
                else (replicate, None, False)
            )
            if true_key not in true_grns:
                logger.warning(
                    "Method '%s', replicate %d: no matching true GRN, skipping PR curve.",
                    method, replicate,
                )
                continue

            try:
                pred_matrix = _load_grn(pred_path)
                true_matrix = _load_grn(true_grns[true_key])

                # Align shapes
                min_n = min(pred_matrix.shape[0], true_matrix.shape[0])
                pred_matrix = pred_matrix[:min_n, :min_n]
                true_matrix = true_matrix[:min_n, :min_n]

                pr_data = compute_precision_recall_curve(
                    pred_matrix, true_matrix, threshold=threshold
                )
            except Exception as exc:
                logger.error(
                    "PR curve failed for method='%s', rep=%d: %s",
                    method, replicate, exc,
                )
                continue

            method_name = _resolve_method_name(method, pt)
            auprc = pr_data["auprc"]

            # Baseline = fraction of positive edges
            n = true_matrix.shape[0]
            off_diag = ~np.eye(n, dtype=bool)
            n_pos = int(np.sum(np.abs(true_matrix[off_diag]) > threshold))
            n_total = int(np.sum(off_diag))
            baseline = n_pos / n_total if n_total > 0 else 0.0

            fig, ax = plt.subplots(figsize=(7, 6))

            ax.plot(
                pr_data["recall"], pr_data["precision"],
                drawstyle="steps-post", linewidth=2,
                label=f"{method_name} (AUPRC = {auprc:.4f})",
                color="#1f77b4",
            )
            ax.axhline(
                y=baseline, color="grey", linestyle="--", linewidth=1,
                label=f"Baseline ({baseline:.4f})",
            )

            ax.set_xlabel("Recall", fontsize=11)
            ax.set_ylabel("Precision", fontsize=11)
            suptitle = f"Precision-Recall Curve — {method_name}, replicate {replicate}"
            if fold is not None:
                suptitle += f", fold {fold}"
            ax.set_title(suptitle, fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.legend(loc="lower left", fontsize=9)
            ax.grid(alpha=0.3)
            fig.tight_layout()

            if fold is not None:
                out_path = output_dir / f"pr_curve_{method_name}_fold{fold}_rep{replicate}.png"
            else:
                out_path = output_dir / f"pr_curve_{method_name}_rep{replicate}.png"
            fig.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            saved.append(out_path)
            logger.info("Saved PR curve: %s", out_path)

    return saved


# ============================================================================
# Cell UMAP plots
# ============================================================================


def plot_cell_umaps(
    run_dir: Path,
    adata_dir: Path,
    output_dir: Path,
    n_pca_components: int = _DEFAULT_PCA_COMPONENTS,
    ko_output_genes: str = "none",
    simulation: str = "false",
    train_tps: Optional[List[float]] = None,
) -> List[Path]:
    """UMAP scatter plots of predicted vs true cell populations for every
    (method, replicate) pair.

    The UMAP reducer is fitted **once** on the true adata per
    (replicate, perturbation) pair and reused across all methods via
    ``transform()``.  This guarantees that conserved training cells
    appear at the same UMAP coordinates in every panel — both within
    a (pred, true) pair and across different methods.

    Args:
        run_dir: Path to ``benchmark/outputs_methods/<run_name>/``.
        adata_dir: Directory containing true ``data_<i>.h5ad`` files.
        output_dir: Directory for saved PNG files.
        n_pca_components: PCA components before UMAP reduction.
        ko_output_genes: KO filter ("none", "all", or gene list).
        simulation: Simulation mode identifier.
        train_tps: Training-timepoint values.  Cells whose ``timepoint``
            falls in this list are considered *conserved training cells*
            and are drawn with hollow markers at reduced alpha so they
            can be distinguished from predicted cells.

    Returns:
        List of saved figure paths.
    """
    from sklearn.decomposition import PCA

    import scanpy as sc
    import umap as umap_module

    true_adatas = discover_replicates(adata_dir, TRUE_ADATA_RE)
    if not true_adatas:
        logger.info("No true adata files in '%s'; skipping UMAP plots.", adata_dir)
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []

    ko_targets = parse_ko_output_genes(ko_output_genes)
    filter_kos = [g for g in ko_targets if g.lower() != "all"]
    select_all_kos = "all" in [g.lower() for g in ko_targets]

    train_tp_set: Optional[set] = set(train_tps) if train_tps else None

    # ------------------------------------------------------------------
    # Pre-compute PCA + UMAP reducers for each (replicate, perturbation)
    # key so that all methods share the same embedding space.
    # ------------------------------------------------------------------
    reducer_cache: Dict[tuple, Tuple[PCA, umap_module.UMAP, np.ndarray, np.ndarray]] = {}
    # cache key → (pca, reducer, true_X_pca, true_X_umap)

    # Two-pass: first discover all needed keys, then compute reducers once.
    needed_keys: set = set()
    for method, method_dir, _fold in iter_method_dirs(run_dir):
        pred_adatas = discover_replicates(method_dir, ADATA_FILE_RE)
        for (replicate, perturbation, pt), _pred_path in pred_adatas.items():
            target_dataset = perturbation or "WT"
            if filter_kos and target_dataset not in {"WT", *filter_kos}:
                continue
            if not select_all_kos and not filter_kos and target_dataset != "WT":
                continue
            true_key = (
                (replicate, perturbation, False)
                if (replicate, perturbation, False) in true_adatas
                else (replicate, None, False)
            )
            if true_key not in true_adatas:
                continue
            needed_keys.add(true_key)

    for true_key in needed_keys:
        try:
            true_adata = read_adata(
                str(true_adatas[true_key]),
                dataset_id=(true_key[1] or "WT"),
                resimulate_missing_ko=simulation != "false",
                ko_gene=true_key[1] if true_key[1] else None,
            )
        except Exception as exc:
            logger.error(
                "Could not load true adata for key=%s: %s", true_key, exc,
            )
            continue

        # Exclude timepoint 0 cells for simulated data
        if simulation != "false" and "timepoint" in true_adata.obs.columns:
            mask = true_adata.obs["timepoint"] != 0
            true_adata = true_adata[mask, :].copy()
            if true_adata.n_obs == 0:
                logger.warning("All cells at timepoint 0 for key=%s; skipping.", true_key)
                continue

        true_X = np.asarray(true_adata.X).astype(float)

        # Align feature dimensions for simulation mode (drop leading "time" col)
        if simulation != "false" and true_X.shape[1] > 1:
            true_X = true_X[:, 1:]

        n_comps = min(n_pca_components, true_X.shape[1] - 1, true_X.shape[0] - 1)
        if n_comps < 2:
            logger.warning("Too few components (%d) for key=%s; skipping.", n_comps, true_key)
            continue

        try:
            pca = PCA(n_components=n_comps)
            true_pca = pca.fit_transform(true_X)
            reducer = umap_module.UMAP()
            true_umap = reducer.fit_transform(true_pca)
            reducer_cache[true_key] = (pca, reducer, true_pca, true_umap)
        except Exception as exc:
            logger.error("PCA/UMAP failed for true key=%s: %s", true_key, exc)

    if not reducer_cache:
        logger.info("No UMAP reducers computed; skipping UMAP plots.")
        return []

    # ------------------------------------------------------------------
    # Second pass: generate plots for each method using cached reducers.
    # ------------------------------------------------------------------
    for method, method_dir, fold in iter_method_dirs(run_dir):

        pred_adatas = discover_replicates(method_dir, ADATA_FILE_RE)

        if not pred_adatas:
            logger.info("Method '%s': no adata files, skipping UMAP.", method)
            continue

        for (replicate, perturbation, pt), pred_path in sorted(pred_adatas.items()):
            target_dataset = perturbation or "WT"
            if filter_kos and target_dataset not in {"WT", *filter_kos}:
                continue
            if not select_all_kos and not filter_kos and target_dataset != "WT":
                continue

            true_key = (
                (replicate, perturbation, False)
                if (replicate, perturbation, False) in true_adatas
                else (replicate, None, False)
            )
            if true_key not in reducer_cache:
                continue

            method_name = _resolve_method_name(method, pt)

            try:
                pred_adata = read_adata(
                    str(pred_path), dataset_id=target_dataset,
                    resimulate_missing_ko=False,
                    ko_gene=target_dataset if target_dataset != "WT" else None,
                )
                true_adata = read_adata(
                    str(true_adatas[true_key]), dataset_id=target_dataset,
                    resimulate_missing_ko=simulation != "false",
                    ko_gene=target_dataset if target_dataset != "WT" else None,
                )
            except Exception as exc:
                logger.error(
                    "Could not load adata for method='%s', rep=%d, pert='%s': %s",
                    method_name, replicate, target_dataset, exc,
                )
                continue

            # Exclude timepoint 0 cells for simulated data
            if simulation != "false":
                if "timepoint" in pred_adata.obs.columns:
                    pred_mask = pred_adata.obs["timepoint"] != 0
                    pred_adata = pred_adata[pred_mask, :].copy()
                if "timepoint" in true_adata.obs.columns:
                    true_mask = true_adata.obs["timepoint"] != 0
                    true_adata = true_adata[true_mask, :].copy()
                if pred_adata.n_obs == 0 or true_adata.n_obs == 0:
                    logger.warning(
                        "All cells at timepoint 0 for method='%s', rep=%d, pert='%s'; skipping.",
                        method_name, replicate, target_dataset,
                    )
                    continue

            pred_X = np.asarray(pred_adata.X).astype(float)
            true_X = np.asarray(true_adata.X).astype(float)

            # Align feature dimensions for simulation mode
            if simulation != "false":
                if true_X.shape[1] == pred_X.shape[1] + 1:
                    true_adata = true_adata[:, 1:].copy()
                    true_X = true_X[:, 1:]
                elif pred_X.shape[1] == true_X.shape[1] + 1:
                    pred_adata = pred_adata[:, 1:].copy()
                    pred_X = pred_X[:, 1:]

            pca, reducer, true_pca, true_umap = reducer_cache[true_key]

            if pred_X.shape[1] != pca.n_features_in_:
                logger.warning(
                    "Method '%s', rep %d, pert '%s': feature mismatch "
                    "(%d vs PCA's %d), skipping.",
                    method_name, replicate, target_dataset,
                    pred_X.shape[1], pca.n_features_in_,
                )
                continue

            # Project pred through the same PCA + UMAP
            try:
                pred_pca = pca.transform(pred_X)
                pred_umap = reducer.transform(pred_pca)
            except Exception as exc:
                logger.error(
                    "UMAP transform failed for method='%s', rep=%d, pert='%s': %s",
                    method_name, replicate, target_dataset, exc,
                )
                continue

            # Store UMAP in adata for scanpy plotting API
            true_adata.obsm["X_umap"] = true_umap
            pred_adata.obsm["X_umap"] = pred_umap

            # Build training-cell masks for visual distinction
            if train_tp_set is not None:
                true_is_train = true_adata.obs["timepoint"].isin(train_tp_set).values
                pred_is_train = pred_adata.obs["timepoint"].isin(train_tp_set).values
            else:
                true_is_train = np.zeros(true_adata.n_obs, dtype=bool)
                pred_is_train = np.zeros(pred_adata.n_obs, dtype=bool)

            # Plot
            fig, axes = plt.subplots(1, 2, figsize=_UMAP_FIGSIZE)
            suptitle = f"UMAP  —  method: {method_name},  replicate: {replicate}"
            if target_dataset != "WT":
                suptitle += f",  perturbation: {target_dataset}"
            if fold is not None:
                suptitle += f",  fold: {fold}"
            fig.suptitle(suptitle, fontsize=11, y=1.01)

            def _fmt(ax: plt.Axes, title: str) -> None:
                ax.set_title(title, fontsize=9)
                ax.set_xlabel("UMAP 1", fontsize=7)
                ax.set_ylabel("UMAP 2", fontsize=7)
                ax.tick_params(labelsize=6)
                if ax.get_legend():
                    ax.get_legend().set_title("Timepoint")
                    ax.get_legend().set_fontsize(6)

            def _plot_umap_with_train_overlay(
                adata, umap_coords: np.ndarray, is_train: np.ndarray,
                ax: plt.Axes, title: str,
            ) -> None:
                """Plot non-training cells via scanpy, then overlay training
                cells as hollow markers with reduced alpha."""
                # Plot all cells via scanpy for the proper colour legend
                sc.pl.umap(adata, color="timepoint", ax=ax,
                           show=False, size=_UMAP_DOT_SIZE,
                           alpha=_UMAP_DOT_ALPHA)
                # Overlay training cells as hollow markers
                if is_train.any():
                    train_coords = umap_coords[is_train]
                    train_tps_vals = adata.obs["timepoint"].values[is_train]
                    unique_tps = np.unique(train_tps_vals)
                    cmap = plt.get_cmap("tab10")
                    for i, tp in enumerate(unique_tps):
                        mask = train_tps_vals == tp
                        ax.scatter(
                            train_coords[mask, 0], train_coords[mask, 1],
                            s=_UMAP_DOT_SIZE * 1.5,
                            facecolors="none",
                            edgecolors=cmap(i % 10),
                            linewidths=0.5,
                            alpha=0.45,
                            zorder=10,
                        )
                _fmt(ax, title)

            _plot_umap_with_train_overlay(
                pred_adata, pred_umap, pred_is_train,
                axes[0], "Predicted cells",
            )
            _plot_umap_with_train_overlay(
                true_adata, true_umap, true_is_train,
                axes[1], "True cells",
            )

            fig.tight_layout()
            if target_dataset == "WT":
                if fold is not None:
                    out_path = output_dir / f"umap_{method_name}_fold{fold}_rep{replicate}.png"
                else:
                    out_path = output_dir / f"umap_{method_name}_rep{replicate}.png"
            else:
                if fold is not None:
                    out_path = output_dir / f"umap_{method_name}_fold{fold}_rep{replicate}_ko_{target_dataset}.png"
                else:
                    out_path = output_dir / f"umap_{method_name}_rep{replicate}_ko_{target_dataset}.png"
            fig.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            saved.append(out_path)
            logger.info("Saved UMAP plot: %s", out_path)

    return saved


# ============================================================================
# Gene-gene correlation scatter plots
# ============================================================================


def plot_gene_corr_scatters(
    run_dir: Path,
    adata_dir: Path,
    output_dir: Path,
    ko_output_genes: str = "none",
    simulation: str = "false",
) -> List[Path]:
    """Scatter plots of gene-pair correlations: predicted vs true, for every
    (method, replicate) pair.

    For each pair of expression matrices, computes the gene‑gene correlation
    matrix, extracts the upper triangle of unique gene pairs, and plots the
    predicted gene‑pair correlations against the true ones.  The
    correlation‑of‑correlations :math:`R` is reported in the plot title —
    this is the same scalar returned by
    :func:`metrics_cells.compute_gene_corr`.

    Args:
        run_dir: Path to ``benchmark/outputs_methods/<run_name>/``.
        adata_dir: Directory containing true ``data_<i>.h5ad`` files.
        output_dir: Directory for saved PNG files.
        ko_output_genes: KO filter (``"none"``, ``"all"``, or comma‑separated
            gene list).
        simulation: Simulation mode identifier.

    Returns:
        List of saved figure paths.
    """
    true_adatas = discover_replicates(adata_dir, TRUE_ADATA_RE)
    if not true_adatas:
        logger.info("No true adata files in '%s'; skipping gene-corr scatters.", adata_dir)
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []

    ko_targets = parse_ko_output_genes(ko_output_genes)
    filter_kos = [g for g in ko_targets if g.lower() != "all"]
    select_all_kos = "all" in [g.lower() for g in ko_targets]

    for method, method_dir, fold in iter_method_dirs(run_dir):

        pred_adatas = discover_replicates(method_dir, ADATA_FILE_RE)

        if not pred_adatas:
            logger.info("Method '%s': no adata files, skipping gene-corr scatters.", method)
            continue

        for (replicate, perturbation, pt), pred_path in sorted(pred_adatas.items()):
            target_dataset = perturbation or "WT"
            if filter_kos and target_dataset not in {"WT", *filter_kos}:
                continue
            if not select_all_kos and not filter_kos and target_dataset != "WT":
                continue

            true_key = (
                (replicate, perturbation, False)
                if (replicate, perturbation, False) in true_adatas
                else (replicate, None, False)
            )
            if true_key not in true_adatas:
                continue

            method_name = _resolve_method_name(method, pt)

            try:
                pred_adata = read_adata(
                    str(pred_path), dataset_id=target_dataset,
                    resimulate_missing_ko=False,
                    ko_gene=target_dataset if target_dataset != "WT" else None,
                )
                true_adata = read_adata(
                    str(true_adatas[true_key]), dataset_id=target_dataset,
                    resimulate_missing_ko=simulation != "false",
                    ko_gene=target_dataset if target_dataset != "WT" else None,
                )
            except Exception as exc:
                logger.error(
                    "Could not load adata for method='%s', rep=%d, pert='%s': %s",
                    method_name, replicate, target_dataset, exc,
                )
                continue

            pred_X = np.asarray(pred_adata.X).astype(float)
            true_X = np.asarray(true_adata.X).astype(float)

            # Align genes via intersection of var_names (handle order mismatches)
            pred_genes = list(pred_adata.var_names)
            true_genes = list(true_adata.var_names)

            # Align feature dimensions for simulation mode (drop leading "time" col)
            if simulation != "false":
                if true_X.shape[1] == pred_X.shape[1] + 1:
                    true_X = true_X[:, 1:]
                    true_genes = true_genes[1:]
                elif pred_X.shape[1] == true_X.shape[1] + 1:
                    pred_X = pred_X[:, 1:]
                    pred_genes = pred_genes[1:]

            common = sorted(set(pred_genes) & set(true_genes))
            if len(common) < 3:
                logger.warning(
                    "Method '%s', rep %d: too few common genes (%d), skipping.",
                    method_name, replicate, len(common),
                )
                continue

            pred_idx = [pred_genes.index(g) for g in common]
            true_idx = [true_genes.index(g) for g in common]
            pred_X = pred_X[:, pred_idx]
            true_X = true_X[:, true_idx]
            n_genes = len(common)

            # Compute gene-gene correlation matrices (columns = genes)
            try:
                corr_pred = np.corrcoef(pred_X, rowvar=False)
                corr_true = np.corrcoef(true_X, rowvar=False)
            except Exception as exc:
                logger.error(
                    "Correlation computation failed for method='%s', rep=%d: %s",
                    method_name, replicate, exc,
                )
                continue

            # Extract upper triangle (unique gene pairs)
            tri = np.triu_indices_from(corr_pred, k=1)
            cp = corr_pred[tri]
            ct = corr_true[tri]

            # Drop non-finite pairs
            mask = np.isfinite(cp) & np.isfinite(ct)
            cp, ct = cp[mask], ct[mask]
            if len(cp) < 2:
                logger.warning(
                    "Method '%s', rep %d: too few valid gene pairs (%d), skipping.",
                    method_name, replicate, len(cp),
                )
                continue

            # Correlation of correlations (same as compute_gene_corr)
            r_val = float(np.corrcoef(cp, ct)[0, 1]) if len(cp) >= 2 else np.nan

            # --- Plot ---
            fig, ax = plt.subplots(figsize=(7, 6.5))

            # Adjust point size and alpha to the number of data points so that
            # sparse plots have larger, more opaque dots and crowded plots have
            # smaller, more transparent dots.
            n_pairs = len(cp)
            point_size = max(1.0, min(20.0, 200.0 / np.sqrt(max(1.0, n_pairs))))
            point_alpha = max(0.12, min(0.55, 1.5 / np.log10(max(10.0, n_pairs))))

            ax.scatter(
                ct, cp,
                alpha=point_alpha, s=point_size,
                c="grey",
                edgecolors="none",
            )

            # Diagonal reference line
            lim = [min(np.nanmin(ct), np.nanmin(cp)), max(np.nanmax(ct), np.nanmax(cp))]
            ax.plot(lim, lim, "k--", alpha=0.25, linewidth=1, label="Perfect agreement")

            suptitle = f"Gene-pair correlation  —  {method_name},  rep {replicate}"
            if target_dataset != "WT":
                suptitle += f",  KO: {target_dataset}"
            if fold is not None:
                suptitle += f",  fold: {fold}"
            ax.set_title(
                suptitle + f"\nR = {r_val:.4f}  ({n_genes} genes, {len(cp):,} pairs)",
                fontsize=10,
            )
            ax.set_xlabel("True gene-pair correlation", fontsize=9)
            ax.set_ylabel("Predicted gene-pair correlation", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8, loc="upper left")

            fig.tight_layout()

            if target_dataset == "WT":
                if fold is not None:
                    out_path = output_dir / f"gene_corr_scatter_{method_name}_fold{fold}_rep{replicate}.png"
                else:
                    out_path = output_dir / f"gene_corr_scatter_{method_name}_rep{replicate}.png"
            else:
                if fold is not None:
                    out_path = output_dir / f"gene_corr_scatter_{method_name}_fold{fold}_rep{replicate}_ko_{target_dataset}.png"
                else:
                    out_path = output_dir / f"gene_corr_scatter_{method_name}_rep{replicate}_ko_{target_dataset}.png"
            fig.savefig(out_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            saved.append(out_path)
            logger.info("Saved gene-corr scatter: %s (R=%.4f)", out_path, r_val)

    return saved


# ============================================================================
# Shared metrics plotting (bar chart + heatmap)
# ============================================================================


def _plot_metrics_bar(
    metrics_df: pd.DataFrame,
    metric_defs: List[Tuple[str, str]],
    output_dir: Path,
    run_name: str,
    prefix: str,
    cmap_name: str = "tab10",
) -> Optional[Path]:
    """Grouped bar chart: one panel per metric, bars grouped by method,
    coloured by replicate.

    Returns the saved path, or ``None`` if no metrics are present.
    """
    present = [(col, label) for col, label in metric_defs if col in metrics_df.columns]
    if not present:
        logger.info("No recognised %s metric columns; skipping bar chart.", prefix)
        return None

    methods = sorted(metrics_df["method"].unique())
    replicates = sorted(metrics_df["replicate"].unique())
    n_reps = len(replicates)
    n_metrics = len(present)

    cmap = plt.get_cmap(cmap_name)
    rep_colours = {rep: cmap(i % 10) for i, rep in enumerate(replicates)}

    ncols = min(4, n_metrics)
    nrows = (n_metrics + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes_flat = np.array(axes).ravel() if n_metrics > 1 else [axes]

    title_suffix = f" — {run_name}" if run_name else ""
    fig.suptitle(f"{prefix.upper()} Metrics{title_suffix}", fontsize=13)

    bar_width = 0.8 / max(n_reps, 1)
    x = np.arange(len(methods))

    for ax_idx, (col, label) in enumerate(present):
        ax = axes_flat[ax_idx]
        for r_idx, rep in enumerate(replicates):
            subset = metrics_df[metrics_df["replicate"] == rep]
            values = [
                subset.loc[subset["method"] == m, col].values[0]
                if len(subset.loc[subset["method"] == m, col]) > 0
                else np.nan
                for m in methods
            ]
            offset = (r_idx - (n_reps - 1) / 2) * bar_width
            ax.bar(x + offset, values, width=bar_width * 0.9,
                   color=rep_colours[rep], label=f"rep {rep}")

        ax.set_title(label, fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=30, ha="right", fontsize=7)
        ax.tick_params(labelsize=7)
        if ax_idx == 0:
            ax.legend(fontsize=7, loc="upper right")

    for ax in axes_flat[n_metrics:]:
        ax.set_visible(False)

    fig.tight_layout()
    bar_path = output_dir / f"{prefix}_metrics_bar{('_' + run_name) if run_name else ''}.png"
    fig.savefig(bar_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s metrics bar chart: %s", prefix, bar_path)
    return bar_path


def _plot_metrics_heatmap(
    metrics_df: pd.DataFrame,
    metric_defs: List[Tuple[str, str]],
    output_dir: Path,
    run_name: str,
    prefix: str,
    cmap_name: str = "YlOrRd",
) -> Optional[Path]:
    """Summary heatmap: rows = metrics, columns = methods, cells = mean value
    across replicates.  Each row is independently normalised to [0, 1].

    Returns the saved path, or ``None`` if no metrics are present.
    """
    present = [(col, label) for col, label in metric_defs if col in metrics_df.columns]
    if not present:
        return None

    methods = sorted(metrics_df["method"].unique())
    metric_cols = [col for col, _ in present]
    metric_labels = [label for _, label in present]
    n_metrics = len(present)

    avg_df = metrics_df.groupby("method")[metric_cols].mean().reindex(methods)
    data = avg_df.values.T  # shape: (n_metrics, n_methods)

    # Normalize each row independently
    norm = np.zeros_like(data, dtype=float)
    for i in range(data.shape[0]):
        row = data[i, :]
        valid = row[~np.isnan(row)]
        if len(valid) > 0:
            rmin, rmax = np.nanmin(row), np.nanmax(row)
            if rmax > rmin:
                norm[i, :] = (row - rmin) / (rmax - rmin)
            else:
                norm[i, :] = 0.5
        else:
            norm[i, :] = np.nan

    title_suffix = f" — {run_name}" if run_name else ""

    fig, ax = plt.subplots(
        figsize=(max(6, len(methods) * 1.2), max(4, n_metrics * 0.6))
    )
    im = ax.imshow(norm, aspect="auto", cmap=cmap_name, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(methods)))
    ax.set_xticklabels(methods, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(np.arange(n_metrics))
    ax.set_yticklabels(metric_labels, fontsize=8)
    ax.set_title(f"Average {prefix.upper()} Metrics per Method{title_suffix}", fontsize=10)

    for i in range(n_metrics):
        for j in range(len(methods)):
            val = data[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.3f}", ha="center", va="center", fontsize=7)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Normalized value (0-1 per row)", fontsize=7)

    fig.tight_layout()
    heatmap_path = output_dir / f"{prefix}_metrics_heatmap{('_' + run_name) if run_name else ''}.png"
    fig.savefig(heatmap_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s metrics heatmap: %s", prefix, heatmap_path)
    return heatmap_path


# ============================================================================
# Public entry points
# ============================================================================


def plot_grn_metrics(
    metrics_df: pd.DataFrame,
    output_dir: Path,
    run_name: str = "",
) -> Optional[Path]:
    """Bar chart + heatmap for GRN metrics."""
    if metrics_df.empty:
        logger.info("metrics_df is empty; skipping GRN metrics plot.")
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    _plot_metrics_heatmap(metrics_df, _GRN_METRICS, output_dir, run_name, "grn", cmap_name="YlGn")
    return _plot_metrics_bar(metrics_df, _GRN_METRICS, output_dir, run_name, "grn")


def plot_cell_metrics(
    metrics_df: pd.DataFrame,
    output_dir: Path,
    run_name: str = "",
) -> Optional[Path]:
    """Bar chart + heatmap for cell metrics."""
    if metrics_df.empty:
        logger.info("metrics_df is empty; skipping cell metrics plot.")
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    _plot_metrics_heatmap(metrics_df, _CELL_METRICS, output_dir, run_name, "cell")
    return _plot_metrics_bar(metrics_df, _CELL_METRICS, output_dir, run_name, "cell")


# ============================================================================
# Standalone CLI
# ============================================================================


def _build_base_parser(desc: str) -> argparse.ArgumentParser:
    """Shared argument parser for the standalone CLI."""
    p = argparse.ArgumentParser(description=desc)
    p.add_argument("global_run_name", help="Benchmark run name (e.g. test_run_new)")
    p.add_argument("run_type", help="Run type (e.g. future_full_test)")
    p.add_argument("--methods_dir", default=None,
                   help="Root for method outputs (default: benchmark/outputs_methods)")
    p.add_argument("--metrics_dir", default=None,
                   help="Root for metric outputs (default: benchmark/outputs_metrics)")
    p.add_argument("--plots_dir", default=None,
                   help="Output directory for plots (default: <metrics_dir>/<run>/plots)")
    return p


def main_grn(argv: Optional[List[str]] = None) -> None:
    """Standalone CLI: generate GRN heatmaps and metrics plots.

    Usage::

        python visualizations.py grn <global_run_name> <run_type> [options]
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = _build_base_parser("Generate GRN visualizations for a benchmark run.")
    args = parser.parse_args(argv)

    script_dir = Path(__file__).resolve().parent
    methods_root = resolve_benchmark_path(script_dir, args.methods_dir, "benchmark/outputs_methods")
    metrics_root = resolve_benchmark_path(script_dir, args.metrics_dir, "benchmark/outputs_metrics")

    run_name = f"{args.global_run_name}_{args.run_type}"
    run_dir = methods_root / run_name
    metrics_run_dir = metrics_root / run_name
    plots_dir = Path(args.plots_dir) if args.plots_dir else metrics_run_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # --- GRN heatmaps ---
    if run_dir.is_dir():
        logger.info("Generating GRN heatmaps from: %s", run_dir)
        paths = plot_grn_heatmaps(run_dir, plots_dir)
        if paths:
            logger.info("Saved %d GRN heatmap(s) to: %s", len(paths), plots_dir)

        logger.info("Generating PR curves from: %s", run_dir)
        pr_paths = plot_grn_pr_curves(run_dir, plots_dir)
        if pr_paths:
            logger.info("Saved %d PR curve(s) to: %s", len(pr_paths), plots_dir)
    else:
        logger.warning("Methods directory not found: %s", run_dir)

    # --- GRN metrics plots ---
    grn_csv = metrics_run_dir / "grn_metrics.csv"
    if grn_csv.is_file():
        logger.info("Generating GRN metrics plots from: %s", grn_csv)
        plot_grn_metrics(pd.read_csv(grn_csv), plots_dir, run_name=run_name)
    else:
        logger.info("No grn_metrics.csv at '%s'; skipping metrics plots.", grn_csv)

    logger.info("Done. Plots saved to: %s", plots_dir)


def main_cells(argv: Optional[List[str]] = None) -> None:
    """Standalone CLI: generate cell UMAP and metrics plots.

    Usage::

        python visualizations.py cell <global_run_name> <run_type> [options]
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = _build_base_parser("Generate cell visualizations for a benchmark run.")
    parser.add_argument("--adata_dir", default=None,
                        help="Directory with true adata files (default: benchmark/data)")
    parser.add_argument("--n_pca_components", type=int, default=_DEFAULT_PCA_COMPONENTS,
                        help="PCA components before UMAP (default: 50)")
    parser.add_argument("--ko_output_genes", default="none",
                        help='KO targets ("none", "all", or comma-separated genes)')
    parser.add_argument("--simulation", default="false",
                        help='Simulation mode (e.g. "false", "simul_replicates")')
    parser.add_argument("--train_tps_path", default=None,
                        help="Path to .npy file with training timepoints "
                             "(e.g. benchmark/data/train_tps.npy). "
                             "Cells at these timepoints are drawn as "
                             "hollow markers to distinguish conserved "
                             "training cells from predicted cells.")
    args = parser.parse_args(argv)

    script_dir = Path(__file__).resolve().parent
    methods_root = resolve_benchmark_path(script_dir, args.methods_dir, "benchmark/outputs_methods")
    metrics_root = resolve_benchmark_path(script_dir, args.metrics_dir, "benchmark/outputs_metrics")
    adata_dir = resolve_benchmark_path(script_dir, args.adata_dir, "benchmark/data")

    # Resolve training timepoints
    train_tps: Optional[List[float]] = None
    if args.train_tps_path:
        train_tps_path = Path(args.train_tps_path)
        if not train_tps_path.is_absolute():
            train_tps_path = script_dir.parent / train_tps_path
        if train_tps_path.is_file():
            tps_arr = np.load(str(train_tps_path))
            # train_tps.npy is shape (N,2): (index, time_value)
            train_tps = sorted(float(v) for v in (tps_arr[:, 1] if tps_arr.ndim == 2 else tps_arr))
            logger.info("Training timepoints: %s", train_tps)
        else:
            logger.warning("train_tps_path '%s' not found; no training-cell overlay.", train_tps_path)

    run_name = f"{args.global_run_name}_{args.run_type}"
    run_dir = methods_root / run_name
    metrics_run_dir = metrics_root / run_name
    plots_dir = Path(args.plots_dir) if args.plots_dir else metrics_run_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # --- UMAP plots ---
    if run_dir.is_dir():
        logger.info("Generating cell UMAP plots from: %s", run_dir)
        paths = plot_cell_umaps(
            run_dir, adata_dir, plots_dir,
            n_pca_components=args.n_pca_components,
            ko_output_genes=args.ko_output_genes,
            simulation=args.simulation,
            train_tps=train_tps,
        )
        if paths:
            logger.info("Saved %d UMAP plot(s) to: %s", len(paths), plots_dir)

        # --- Gene-gene correlation scatter plots ---
        logger.info("Generating gene-corr scatter plots from: %s", run_dir)
        gc_paths = plot_gene_corr_scatters(
            run_dir, adata_dir, plots_dir,
            ko_output_genes=args.ko_output_genes,
            simulation=args.simulation,
        )
        if gc_paths:
            logger.info("Saved %d gene-corr scatter(s) to: %s", len(gc_paths), plots_dir)
    else:
        logger.warning("Methods directory not found: %s", run_dir)

    # --- Cell metrics plots ---
    cell_csv = metrics_run_dir / "cell_metrics.csv"
    if cell_csv.is_file():
        logger.info("Generating cell metrics plots from: %s", cell_csv)
        plot_cell_metrics(pd.read_csv(cell_csv), plots_dir, run_name=run_name)
    else:
        logger.info("No cell_metrics.csv at '%s'; skipping metrics plots.", cell_csv)

    logger.info("Done. Plots saved to: %s", plots_dir)


def main(argv: Optional[List[str]] = None) -> None:
    """Unified CLI dispatcher.

    Usage::

        python visualizations.py grn  <global_run_name> <run_type> [options]
        python visualizations.py cell <global_run_name> <run_type> [options]
    """
    if argv is None:
        argv = sys.argv[1:]

    if not argv or argv[0] not in ("grn", "cell"):
        print("Usage: python visualizations.py {grn|cell} <global_run_name> <run_type> [options]",
              file=sys.stderr)
        sys.exit(1)

    subcommand = argv[0]
    rest = argv[1:]

    if subcommand == "grn":
        main_grn(rest)
    else:
        main_cells(rest)


if __name__ == "__main__":
    main()
