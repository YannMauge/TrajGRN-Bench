#!/usr/bin/env python3
"""
compute_metrics.py — Compute GRN and cell metrics for all outputs of a
benchmark run.

For a given benchmark run (identified by global_run_name and run_type), this
script:
  1. Scans ``benchmark/outputs_methods/<global_run_name>_<run_type>/`` for
     method folders and replicates.
  2. Computes threshold-free GRN metrics (edge / signed AUPRC, AUROC, …) for
     two edge sets (gene→gene, and with Stimulus→gene) by comparing each
     method's ``data_<i>_GRN.npy`` against ``GRN_true/data_<i>_GRN.npy``.
  3. Computes cell metrics (Wasserstein, JS, mean-expression error, gene-gene
     correlation) per evaluated timepoint (test timepoints, see
     metrics_cells) on log1p and raw scales, comparing each method's
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
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
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
        SCALES,
        MissingTruthError,
        align_genes,
        compute_cell_metrics,
        impute_missing_worst,
        summarize_per_timepoint,
        load_expression,
        select_eval_timepoints,
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
        SCALES,
        MissingTruthError,
        align_genes,
        compute_cell_metrics,
        impute_missing_worst,
        summarize_per_timepoint,
        load_expression,
        select_eval_timepoints,
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

OUTPUT_MODES = ("full_full", "full_train", "full_test", "no_traj")


def _method_label(method: str, perturbation_training: bool) -> str:
    if not perturbation_training or method.endswith(PERTURBATION_TRAINING_SUFFIX):
        return method
    return f"{method}{PERTURBATION_TRAINING_SUFFIX}"


def _order_columns(rows: List[Dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    priority = [c for c in ("method", "replicate", "fold", "perturbation") if c in df.columns]
    return df[priority + [c for c in df.columns if c not in priority]]


def output_mode_from_run_type(run_type: str) -> Optional[str]:
    """``run_type`` is ``<train_data>_<output_mode>``."""
    return next((m for m in OUTPUT_MODES if run_type.endswith(m)), None)


def load_split(fold_dir: Path, adata_dir: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Train/test timepoints saved by benchmark_run.py in the run (or fold)
    folder; falls back to *adata_dir* (last split written there) with a warning."""
    for folder in (fold_dir, adata_dir):
        train_p, test_p = folder / "train_tps.npy", folder / "test_tps.npy"
        if train_p.is_file() and test_p.is_file():
            if folder == adata_dir:
                logger.warning(
                    f"No train/test split saved in {fold_dir}; using {adata_dir} "
                    "(may belong to another run)."
                )
            return np.load(train_p), np.load(test_p)
    return None, None


def compute_all_grn_metrics(run_dir: Path) -> pd.DataFrame:
    """GRN metrics of each method's common network (``data_<i>_GRN.npy``,
    learned from WT and, for ``_perturbation_training`` variants, the training
    perturbations) against the ground truth. KO-specific network files are
    not scored.

    Columns: method, replicate, [fold], perturbation, then
    ``<scope>_<metric>`` (see metrics_grn). A failing item is logged and
    skipped; it never aborts the run.
    """
    true_grns = discover_replicates(run_dir / GRN_TRUE_FOLDER, GRN_FILE_RE)
    if not true_grns:
        logger.info(f"No ground-truth GRN files in '{run_dir / GRN_TRUE_FOLDER}'. Skipping GRN metrics.")
        return pd.DataFrame()

    rows: List[Dict[str, Any]] = []
    for method, method_dir, fold in iter_method_dirs(run_dir):
        pred_grns = discover_replicates(method_dir, GRN_FILE_RE)
        if not pred_grns:
            logger.info(f"Method '{method}': no GRN files found, skipping GRN metrics.")
            continue
        for (replicate, perturbation, perturbation_training), pred_path in pred_grns.items():
            if perturbation is not None:
                # KO-specific networks (data_i_ko_<g>_GRN.npy) are a method's
                # internal state under a KO, not a network to recover: only the
                # common network learned from WT (+ training perturbations) is
                # scored against the ground truth.
                continue
            true_key = (replicate, None, False)
            if true_key not in true_grns:
                logger.warning(f"Method '{method}', replicate {replicate}: no matching ground-truth GRN, skipping.")
                continue
            label = _method_label(method, perturbation_training)
            try:
                metrics = compute_grn_metrics(load_grn(str(pred_path)), load_grn(str(true_grns[true_key])))
            except Exception as exc:  # one bad output must not abort the run
                logger.error(f"GRN metrics failed for method='{label}', replicate={replicate}: {exc}")
                continue
            row: Dict[str, Any] = {"method": label, "replicate": replicate,
                                   "perturbation": perturbation or "WT"}
            if fold is not None:
                row["fold"] = fold
            row.update(flatten_metrics(metrics))
            rows.append(row)
            logger.info(f"GRN metrics computed: method={label}, replicate={replicate}"
                        + (f", fold={fold}" if fold is not None else ""))
    return _order_columns(rows) if rows else pd.DataFrame()


def compute_all_cell_metrics(
    run_dir: Path,
    adata_dir: Path,
    simulation: str,
    output_mode: Optional[str],
    n_pca_components: int = 50,
    ko_output_genes: str = "none",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per-timepoint cell metrics for every (method, replicate[, KO]) pair.

    Returns ``(summary_df, per_timepoint_df)``. Evaluated timepoints come from
    the split saved in each run / fold folder (see metrics_cells). A failing
    item is logged and skipped; it never aborts the run.
    """
    true_adatas = discover_replicates(adata_dir, TRUE_ADATA_RE)
    if not true_adatas:
        logger.info(f"No ground-truth adata files found in '{adata_dir}'. Skipping cell metrics.")
        return pd.DataFrame(), pd.DataFrame()
    if output_mode == "no_traj":
        logger.info("output_mode=no_traj: no trajectories to score.")
        return pd.DataFrame(), pd.DataFrame()

    ko_targets = parse_ko_output_genes(ko_output_genes)
    filter_kos = [g for g in ko_targets if g.lower() != "all"]
    select_all_kos = "all" in [g.lower() for g in ko_targets]
    true_cache: Dict[Tuple[Path, str], Tuple[np.ndarray, List[str], np.ndarray]] = {}
    rows: List[Dict[str, Any]] = []
    tp_rows: List[pd.DataFrame] = []

    for method, method_dir, fold in iter_method_dirs(run_dir):
        pred_adatas = discover_replicates(method_dir, ADATA_FILE_RE)
        if not pred_adatas:
            logger.info(f"Method '{method}': no adata files found, skipping cell metrics.")
            continue
        train_tps, test_tps = load_split(method_dir.parent, adata_dir)
        if train_tps is None:
            logger.error(f"No train/test split found for {method_dir}; skipping cell metrics.")
            continue
        eval_tps, eval_set = select_eval_timepoints(output_mode, train_tps, test_tps)
        # subsample_full: methods trained on a subsample of the cells; score
        # against the held-out cells only.
        sub_path = method_dir.parent / "subsample_train_ids.npy"
        train_subsample = np.load(sub_path, allow_pickle=True) if sub_path.is_file() else None
        if train_subsample is not None:
            eval_set = "heldout"

        for (replicate, perturbation, perturbation_training), pred_path in pred_adatas.items():
            target = perturbation or "WT"
            if filter_kos and target not in {"WT", *filter_kos}:
                continue
            if not select_all_kos and not filter_kos and target != "WT":
                continue
            true_key = (replicate, perturbation, False)
            if true_key not in true_adatas:
                true_key = (replicate, None, False)
            if true_key not in true_adatas:
                logger.warning(f"Method '{method}', replicate {replicate}: no matching ground-truth adata, skipping.")
                continue
            label = _method_label(method, perturbation_training)
            try:
                cache_key = (true_adatas[true_key], target, str(sub_path) if train_subsample is not None else None)
                if cache_key not in true_cache:
                    # KO predictions are only scored against simulated KO cells.
                    tx, tv, tt, _ = load_expression(
                        str(true_adatas[true_key]), dataset_id=target, strict=True,
                        exclude_obs=train_subsample,
                    )
                    true_cache[cache_key] = (tx, tv, tt)
                tx, tv, tt = true_cache[cache_key]
                px, pv, pt, runtime = load_expression(str(pred_path), dataset_id=target)
                px, tx_aligned, _ = align_genes(px, pv, tx, tv)
                summary, per_tp = compute_cell_metrics(
                    px, pt, tx_aligned, tt, eval_tps, n_pca_components=n_pca_components,
                )
            except MissingTruthError as exc:
                logger.warning(f"Not scored: method='{label}', replicate={replicate}, "
                               f"perturbation={target}: {exc}")
                continue
            except Exception as exc:  # one bad output must not abort the run
                logger.error(f"Cell metrics failed for method='{label}', replicate={replicate}, "
                             f"perturbation={target}: {exc}")
                continue
            ident: Dict[str, Any] = {"method": label, "replicate": replicate, "perturbation": target}
            if fold is not None:
                ident["fold"] = fold
            rows.append({**ident, "eval_set": eval_set, **summary,
                         "runtime": runtime if runtime is not None else np.nan})
            if not per_tp.empty:
                tp_rows.append(per_tp.assign(**ident, eval_set=eval_set))
            logger.info(f"Cell metrics computed: method={label}, replicate={replicate}, perturbation={target}, "
                        f"{eval_set} timepoints={eval_tps}" + (f", fold={fold}" if fold is not None else ""))

    if not rows:
        return pd.DataFrame(), pd.DataFrame()
    if not tp_rows:
        return _order_columns(rows), pd.DataFrame()

    # Missing timepoints get the worst value any method reached there, then
    # each method's averages are recomputed from the completed table.
    per_tp_df = impute_missing_worst(
        pd.concat(tp_rows, ignore_index=True), ["replicate", "perturbation", "fold"],
    )
    ident_cols = [c for c in ("method", "replicate", "fold", "perturbation") if c in per_tp_df.columns]
    for row in rows:
        mask = np.ones(len(per_tp_df), dtype=bool)
        for col in ident_cols:
            mask &= (per_tp_df[col] == row.get(col)).to_numpy() if col in row else per_tp_df[col].isna().to_numpy()
        sub = per_tp_df[mask]
        if not sub.empty:
            n_eval = int(len(sub) // len(SCALES))
            row.update(summarize_per_timepoint(sub, n_eval))
            row["n_imputed_timepoints"] = int(sub.loc[sub["scale"] == SCALES[0], "imputed_worst"].sum())
    return _order_columns(rows), _order_columns(per_tp_df.to_dict(orient="records"))


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
    grn_df = compute_all_grn_metrics(run_dir)

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
    output_mode = output_mode_from_run_type(args.run_type)
    cell_df, cell_tp_df = compute_all_cell_metrics(
        run_dir,
        adata_dir,
        simulation=simulation_var,
        output_mode=output_mode,
        ko_output_genes=args.ko_output_genes,
    )

    if cell_df.empty:
        logger.info("No cell metrics were produced.")
    else:
        cell_path = output_dir / "cell_metrics.csv"
        cell_df.to_csv(cell_path, index=False)
        logger.info(f"Cell metrics saved to: {cell_path}")
        if not cell_tp_df.empty:
            tp_path = output_dir / "cell_metrics_per_timepoint.csv"
            cell_tp_df.to_csv(tp_path, index=False)
            logger.info(f"Per-timepoint cell metrics saved to: {tp_path}")
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
