#!/usr/bin/env python3
"""
compute_metrics.py — Compute GRN and cell metrics for all outputs of a
benchmark run.

For a given benchmark run (identified by global_run_name and run_type), this
script:
  1. Scans ``benchmark/outputs_methods/<global_run_name>_<run_type>/`` for
     method folders and replicates.
  2. Computes GRN metrics (precision, recall, F1, AUPRC, …) by comparing each
     method's ``data_<i>_GRN.npy`` against ``GRN_true/data_<i>_GRN.npy``.
  3. Computes cell metrics (Wasserstein, KL, …) by comparing each method's
     ``data_<i>_adata.h5ad`` against ``benchmark/data/data_<i>.h5ad``.
  4. Stores all results as CSV files under
     ``benchmark/outputs_metrics/<global_run_name>_<run_type>/``, then
     generates visualizations.

Usage:
    python compute_metrics.py <global_run_name> <run_type> <simulation> [options]
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd

try:
    from ._common import (
        ADATA_FILE_RE,
        GRN_FILE_RE,
        GRN_TRUE_FOLDER,
        PERTURBATION_TRAINING_SUFFIX,
        TRUE_ADATA_RE,
        discover_replicates,
        flatten_metrics,
        is_leave_one_out_structure,
        iter_method_dirs,
        parse_ko_output_genes,
        resolve_benchmark_path,
        setup_logging,
    )
    from .metrics_cells import (
        CellMetricsError,
        compute_cell_metrics,
        load_adata,
    )
    from .metrics_grn import (
        GRNMetricsError,
        compute_grn_metrics,
        load_grn,
    )
    from .visualizations import (
        plot_cell_metrics,
        plot_cell_umaps,
        plot_gene_corr_scatters,
        plot_grn_heatmaps,
        plot_grn_metrics,
        plot_grn_pr_curves,
    )
except ImportError:
    from _common import (  # type: ignore[no-redef]
        ADATA_FILE_RE,
        GRN_FILE_RE,
        GRN_TRUE_FOLDER,
        PERTURBATION_TRAINING_SUFFIX,
        TRUE_ADATA_RE,
        discover_replicates,
        flatten_metrics,
        is_leave_one_out_structure,
        iter_method_dirs,
        parse_ko_output_genes,
        resolve_benchmark_path,
        setup_logging,
    )
    from metrics_cells import (  # type: ignore[no-redef]
        CellMetricsError,
        compute_cell_metrics,
        load_adata,
    )
    from metrics_grn import (  # type: ignore[no-redef]
        GRNMetricsError,
        compute_grn_metrics,
        load_grn,
    )
    from visualizations import (  # type: ignore[no-redef]
        plot_cell_metrics,
        plot_cell_umaps,
        plot_gene_corr_scatters,
        plot_grn_heatmaps,
        plot_grn_metrics,
        plot_grn_pr_curves,
    )

for module in ["anndata"]:
    warnings.filterwarnings("ignore", module=module)

logger = setup_logging(__name__)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_all_grn_metrics(
    run_dir: Path,
    simulation: str,
    quantile: float,
) -> pd.DataFrame:
    """
    Compute GRN metrics for every (method, replicate) pair where a true GRN exists.

    Args:
        run_dir:  Path to benchmark/outputs_methods/<run_name>/
        quantile: Quantile threshold used when calling compute_grn_metrics.

    Returns:
        DataFrame with columns [method, replicate, <metric columns>].
        Returns an empty DataFrame if no GRN_true folder or no predictions found.
    """
    grn_true_dir = run_dir / GRN_TRUE_FOLDER
    true_grns = discover_replicates(grn_true_dir, GRN_FILE_RE)

    if not true_grns:
        logger.info(
            f"No ground-truth GRN files found in '{grn_true_dir}'. "
            "Skipping GRN metrics."
        )
        return pd.DataFrame()

    rows = []

    for method, method_dir, fold in iter_method_dirs(run_dir):

        pred_grns = discover_replicates(method_dir, GRN_FILE_RE)

        if not pred_grns:
            logger.info(f"Method '{method}': no GRN files found, skipping GRN metrics.")
            continue

        for (replicate, perturbation, perturbation_training), pred_path in sorted(pred_grns.items()):
            true_key = (
                (replicate, perturbation, False)
                if (replicate, perturbation, False) in true_grns
                else (replicate, None, False)
            )
            if true_key not in true_grns:
                logger.warning(
                    f"Method '{method}', replicate {replicate}: "
                    "no matching ground-truth GRN, skipping."
                )
                continue

            true_path = true_grns[true_key]
            try:
                pred_grn = load_grn(str(pred_path))
                true_grn = load_grn(str(true_path))
                metrics = compute_grn_metrics(pred_grn, true_grn, simulation=simulation, quantile=quantile)
                flat = flatten_metrics(metrics)
                method_name = (
                    method
                    if (not perturbation_training or method.endswith(PERTURBATION_TRAINING_SUFFIX))
                    else f"{method}{PERTURBATION_TRAINING_SUFFIX}"
                )
                flat["method"] = method_name
                flat["replicate"] = replicate
                flat["perturbation"] = perturbation or "WT"
                if fold is not None:
                    flat["fold"] = fold
                rows.append(flat)
                logger.info(
                    f"GRN metrics computed: method={method_name}, replicate={replicate}"
                    + (f", fold={fold}" if fold is not None else "")
                )
            except GRNMetricsError as exc:
                logger.error(
                    f"GRN metrics failed for method='{method}', replicate={replicate}: {exc}"
                )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    # Reorder columns: method, replicate, [fold], then the rest
    priority = ["method", "replicate"]
    if "fold" in df.columns:
        priority.append("fold")
    cols = priority + [c for c in df.columns if c not in priority]
    return df[cols]


def compute_all_cell_metrics(
    run_dir: Path,
    adata_dir: Path,
    simulation: str,
    n_pca_components: int = 50,
    ko_output_genes: str = "none",
    resimulate_missing_ko: bool = True,
) -> pd.DataFrame:
    """
    Compute cell metrics for every (method, replicate) pair where a true adata exists.

    Args:
        run_dir:          Path to benchmark/outputs_methods/<run_name>/
        adata_dir:        Directory containing true adata files (data_<i>.h5ad).
        n_pca_components: Number of PCA components for cell metric computation.

    Returns:
        DataFrame with columns [method, replicate, <metric columns>].
        Returns an empty DataFrame if no true adata files are found.
    """
    # True adata files follow the pattern data_<i>.h5ad (no _adata suffix)
    true_adatas = discover_replicates(adata_dir, TRUE_ADATA_RE)

    if not true_adatas:
        logger.info(
            f"No ground-truth adata files found in '{adata_dir}'. "
            "Skipping cell metrics."
        )
        return pd.DataFrame()

    rows = []
    ko_targets = parse_ko_output_genes(ko_output_genes)
    filter_kos = [g for g in ko_targets if g.lower() != "all"]
    select_all_kos = "all" in [g.lower() for g in ko_targets]

    for method, method_dir, fold in iter_method_dirs(run_dir):

        pred_adatas = discover_replicates(method_dir, ADATA_FILE_RE)

        if not pred_adatas:
            logger.info(f"Method '{method}': no adata files found, skipping cell metrics.")
            continue

        for (replicate, perturbation, perturbation_training), pred_path in sorted(pred_adatas.items()):
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
                logger.warning(
                    f"Method '{method}', replicate {replicate}: "
                    "no matching ground-truth adata, skipping."
                )
                continue

            true_path = true_adatas[true_key]
            try:
                pred_cells, pred_labels, pred_runtime = load_adata(
                    str(pred_path),
                    dataset_id=target_dataset,
                    resimulate_missing_ko=False,
                    ko_gene=target_dataset if target_dataset != "WT" else None,
                )
                true_cells, true_labels, _ = load_adata(
                    str(true_path),
                    dataset_id=target_dataset,
                    resimulate_missing_ko=resimulate_missing_ko and simulation != "false",
                    ko_gene=target_dataset if target_dataset != "WT" else None,
                )
                if simulation != "false":
                    if true_cells.shape[1] == pred_cells.shape[1] + 1:
                        true_cells = true_cells[:, 1:].copy()
                    elif pred_cells.shape[1] == true_cells.shape[1] + 1:
                        pred_cells = pred_cells[:, 1:].copy()
                metrics = compute_cell_metrics(
                    pred_cells,
                    true_cells,
                    pred_labels,
                    true_labels,
                    n_pca_components=n_pca_components,
                    runtime=pred_runtime,
                )
                method_name = (
                    method
                    if (not perturbation_training or method.endswith(PERTURBATION_TRAINING_SUFFIX))
                    else f"{method}{PERTURBATION_TRAINING_SUFFIX}"
                )
                metrics["method"] = method_name
                metrics["replicate"] = replicate
                metrics["perturbation"] = target_dataset
                if fold is not None:
                    metrics["fold"] = fold
                rows.append(metrics)
                logger.info(
                    f"Cell metrics computed: method={method_name}, replicate={replicate}, perturbation={target_dataset}"
                    + (f", fold={fold}" if fold is not None else "")
                )
            except CellMetricsError as exc:
                logger.error(
                    f"Cell metrics failed for method='{method}', replicate={replicate}, perturbation={target_dataset}: {exc}"
                )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    priority = ["method", "replicate"]
    if "fold" in df.columns:
        priority.append("fold")
    cols = priority + [c for c in df.columns if c not in priority]
    return df[cols]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: Optional[list] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description=(
            "Compute GRN and cell metrics for all outputs of a benchmark run "
            "and store results in a structured format."
        )
    )
    parser.add_argument("global_run_name", help="Name of the benchmark run (e.g. test_run_new)")
    parser.add_argument(
        "run_type",
        help=(
            "Combined train_data and output_mode identifier "
            "(e.g. future_full_test, full_no_traj, subsample_full_full_full)"
        ),
    )
    parser.add_argument(
        "simulation",
        help=(
            "Type of simulation"
            "(e.g. Harissa, ...)"
        ),
    )
    parser.add_argument(
        "--quantile",
        type=float,
        default=0.9,
        help="Quantile threshold for GRN binarization (default: 0.9)",
    )
    parser.add_argument(
        "--adata_dir",
        default=None,
        help=(
            "Directory containing the true adata files (data_<i>.h5ad). "
            "Defaults to benchmark/data relative to this script."
        ),
    )
    parser.add_argument(
        "--methods_dir",
        default=None,
        help=(
            "Root directory for method outputs. "
            "Defaults to benchmark/outputs_methods relative to this script."
        ),
    )
    parser.add_argument(
        "--metrics_dir",
        default=None,
        help=(
            "Root directory for metric outputs. "
            "Defaults to benchmark/outputs_metrics relative to this script."
        ),
    )
    parser.add_argument(
        "--no_plots",
        action="store_true",
        default=False,
        help="Skip generating visualization plots (heatmaps and metrics charts).",
    )
    parser.add_argument(
        "--ko_output_genes",
        default="none",
        help='KO targets to score for cell metrics ("none", "all", or comma-separated gene names).',
    )

    args = parser.parse_args(argv)

    script_dir = Path(__file__).resolve().parent

    methods_root = resolve_benchmark_path(script_dir, args.methods_dir, "benchmark/outputs_methods")
    metrics_root = resolve_benchmark_path(script_dir, args.metrics_dir, "benchmark/outputs_metrics")
    adata_dir = resolve_benchmark_path(script_dir, args.adata_dir, "benchmark/data")
    simulation_var = args.simulation

    run_name = f"{args.global_run_name}_{args.run_type}"
    run_dir = methods_root / run_name
    output_dir = metrics_root / run_name

    if not run_dir.is_dir():
        logger.error(f"Methods output directory not found: {run_dir}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Processing run: {run_name}")
    logger.info(f"  Methods dir : {run_dir}")
    logger.info(f"  Output dir  : {output_dir}")
    logger.info(f"  True adata  : {adata_dir}")

    # ------------------------------------------------------------------
    # GRN metrics
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("Computing GRN metrics …")
    grn_df = compute_all_grn_metrics(run_dir, simulation=simulation_var, quantile=args.quantile)

    if grn_df.empty:
        logger.info("No GRN metrics were produced.")
    else:
        grn_path = output_dir / "grn_metrics.csv"
        grn_df.to_csv(grn_path, index=False)
        logger.info(f"GRN metrics saved to: {grn_path}")
        logger.info(f"  Shape: {grn_df.shape}")

    # ------------------------------------------------------------------
    # Cell metrics
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("Computing cell metrics …")
    cell_df = compute_all_cell_metrics(
        run_dir,
        adata_dir,
        simulation=simulation_var,
        ko_output_genes=args.ko_output_genes,
        resimulate_missing_ko=True,
    )

    if cell_df.empty:
        logger.info("No cell metrics were produced.")
    else:
        cell_path = output_dir / "cell_metrics.csv"
        cell_df.to_csv(cell_path, index=False)
        logger.info(f"Cell metrics saved to: {cell_path}")
        logger.info(f"  Shape: {cell_df.shape}")

    # ------------------------------------------------------------------
    # Visualizations
    # ------------------------------------------------------------------
    if not args.no_plots:
        plots_dir = output_dir / "plots"
        logger.info("=" * 60)
        logger.info("Generating GRN visualizations …")

        heatmap_paths = plot_grn_heatmaps(run_dir, plots_dir)
        if heatmap_paths:
            logger.info(f"Saved {len(heatmap_paths)} GRN heatmap(s) to: {plots_dir}")

        pr_curve_paths = plot_grn_pr_curves(run_dir, plots_dir)
        if pr_curve_paths:
            logger.info(f"Saved {len(pr_curve_paths)} PR curve(s) to: {plots_dir}")

        if not grn_df.empty:
            plot_grn_metrics(grn_df, plots_dir, run_name=run_name)
            logger.info(f"GRN metrics plots saved to: {plots_dir}")

        logger.info("=" * 60)
        logger.info("Generating cell visualizations …")

        umap_paths = plot_cell_umaps(
            run_dir,
            adata_dir,
            plots_dir,
            ko_output_genes=args.ko_output_genes,
            simulation=simulation_var,
        )
        if umap_paths:
            logger.info(f"Saved {len(umap_paths)} UMAP plot(s) to: {plots_dir}")

        logger.info("=" * 60)
        logger.info("Generating gene-correlation scatter plots …")
        gc_paths = plot_gene_corr_scatters(
            run_dir,
            adata_dir,
            plots_dir,
            ko_output_genes=args.ko_output_genes,
            simulation=simulation_var,
        )
        if gc_paths:
            logger.info(f"Saved {len(gc_paths)} gene-corr scatter(s) to: {plots_dir}")

        if not cell_df.empty:
            plot_cell_metrics(cell_df, plots_dir, run_name=run_name)
            logger.info(f"Cell metrics plots saved to: {plots_dir}")

    logger.info("=" * 60)
    logger.info("Done.")


if __name__ == "__main__":
    main()
