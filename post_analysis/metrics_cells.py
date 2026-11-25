#!/usr/bin/env python3
"""Metrics for evaluating predicted cell populations against true ones.

Cells are compared **per evaluated timepoint** and the per-timepoint values
are averaged. By default the evaluated timepoints are the test timepoints of
the run; when a run predicts no test timepoint (``output_mode=full_train`` or
an empty test set) the predicted training timepoints except the first one are
used instead (reconstruction, ``eval_set = "train"``).

Every metric is computed on two scales and suffixed accordingly:
    ``_log``  log1p of expression (negative predictions clipped to 0)
    ``_raw``  expression as stored in ``.X``

Metrics (lower is better unless stated):
    wasserstein           Exact earth mover's distance (W1, uniform weights,
                          POT ``ot.emd2``) in a PCA space fitted on all true
                          cells, so every timepoint and method share one basis.
    js_divergence         Gaussian Jensen–Shannon divergence in the same PCA
                          space (means *and* covariances are compared).
    mean_expr_rmse        RMSE between per-gene mean expression.
    gene_corr_similarity  Pearson correlation between the gene–gene
                          correlation matrices (higher is better).
Plus ``runtime`` (from ``adata.uns``), ``n_eval_timepoints`` and
``n_missing_timepoints`` (evaluated timepoints with < 2 predicted cells).

Missing timepoints are penalised by :func:`impute_missing_worst` (used by
compute_metrics.py): they receive the worst value any method reached at that
timepoint before the per-timepoint values are averaged.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

try:
    from ._common import setup_logging
except ImportError:
    from _common import setup_logging  # type: ignore[no-redef]

warnings.filterwarnings("ignore", module="anndata")

logger = setup_logging(__name__)

SCALES: Tuple[str, ...] = ("log", "raw")
CELL_METRIC_NAMES: Tuple[str, ...] = (
    "wasserstein", "js_divergence", "mean_expr_rmse", "gene_corr_similarity",
)
MIN_CELLS_PER_TP = 2
# Direction of each metric, used to pick the "worst" value for imputation.
HIGHER_IS_BETTER: Dict[str, bool] = {
    "wasserstein": False, "js_divergence": False,
    "mean_expr_rmse": False, "gene_corr_similarity": True,
}


class CellMetricsError(Exception):
    """Raised when cell metric computation fails."""


class MissingTruthError(CellMetricsError):
    """No ground-truth cells for the requested dataset (e.g. a KO that was
    not simulated): the prediction is simply not scored."""


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_adata(filepath: str, dataset_id: str = "WT", strict: bool = False) -> Any:
    """Read an .h5ad file and subset it to ``obs['dataset_id'] == dataset_id``.

    Files without a ``dataset_id`` column are returned whole. If the column
    exists but *dataset_id* is absent: with ``strict=True`` (ground truths) a
    :class:`CellMetricsError` is raised — a knockout is only scored against a
    *simulated* knockout, never a fabricated one; otherwise (predictions) the
    whole file is returned with a warning.
    """
    import anndata as ad  # type: ignore
    try:
        adata = ad.read_h5ad(filepath)
    except Exception as e:
        raise CellMetricsError(f"Cannot read {filepath}: {e}") from e

    if "dataset_id" not in adata.obs.columns or not dataset_id:
        return adata

    mask = adata.obs["dataset_id"].astype(str) == str(dataset_id)
    if int(mask.sum()) > 0:
        return adata[mask].copy()
    if strict:
        raise MissingTruthError(
            f"No cells with dataset_id='{dataset_id}' in {filepath} (no simulated truth for it)"
        )
    logger.warning("Dataset '%s' not found in %s; using full dataset.", dataset_id, filepath)
    return adata


def cell_times(adata: Any) -> np.ndarray:
    """Per-cell time: ``obs['time']`` where finite, else ``obs['timepoint']``.

    Method outputs usually set ``time`` on predicted cells while concatenated
    training cells carry both columns, so the two are merged per cell.
    """
    times = np.full(adata.n_obs, np.nan)
    for key in ("timepoint", "time"):  # later keys take precedence
        if key in adata.obs.columns:
            values = pd.to_numeric(adata.obs[key], errors="coerce").to_numpy(dtype=float)
            times = np.where(np.isfinite(values), values, times)
    if not np.any(np.isfinite(times)):
        raise CellMetricsError("No usable obs['time'] / obs['timepoint'] column")
    return times


def load_expression(
    filepath: str,
    dataset_id: str = "WT",
    strict: bool = False,
    exclude_obs: Optional[Sequence[str]] = None,
) -> Tuple[np.ndarray, List[str], np.ndarray, Optional[float]]:
    """Return ``(X, var_names, times, runtime)`` with a dense float ``X``.

    *exclude_obs*: cell names to drop (e.g. the training subsample of
    ``subsample_full``, so that only held-out cells remain)."""
    adata = read_adata(filepath, dataset_id, strict=strict)
    if exclude_obs is not None:
        keep = ~adata.obs_names.isin(np.asarray(exclude_obs).astype(str))
        adata = adata[keep].copy()
    if adata.X is None or adata.n_obs == 0:
        raise CellMetricsError(f"No cells / expression matrix in: {filepath}")
    x = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    runtime: Optional[float] = None
    if "runtime" in adata.uns:
        try:
            runtime = float(adata.uns["runtime"])
        except (TypeError, ValueError):
            logger.warning("Bad runtime value in %s; skipping.", filepath)
    return x.astype(float), [str(v) for v in adata.var_names], cell_times(adata), runtime


def _is_stimulus(name: str) -> bool:
    return name.lower().startswith("stimulus")


def align_genes(
    pred_x: np.ndarray, pred_vars: Sequence[str], true_x: np.ndarray, true_vars: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Drop stimulus columns and put predicted genes in the true gene order.

    Matching is by name when every true gene is present in the prediction;
    otherwise positional (with a warning) if the gene counts agree, dropping
    an unnamed leading column of the prediction when it has exactly one extra.
    """
    t_keep = [i for i, v in enumerate(true_vars) if not _is_stimulus(v)]
    genes = [true_vars[i] for i in t_keep]
    true_x = true_x[:, t_keep]
    p_index = {v: i for i, v in enumerate(pred_vars) if not _is_stimulus(v)}
    if all(g in p_index for g in genes):
        return pred_x[:, [p_index[g] for g in genes]], true_x, genes
    p_keep = [i for i, v in enumerate(pred_vars) if not _is_stimulus(v)]
    pred_x = pred_x[:, p_keep]
    if pred_x.shape[1] == len(genes) + 1:
        pred_x = pred_x[:, 1:]
    if pred_x.shape[1] != len(genes):
        raise CellMetricsError(
            f"Cannot align genes: prediction has {pred_x.shape[1]} genes, truth {len(genes)}"
        )
    logger.warning("Predicted gene names do not match the truth; aligning genes by position.")
    return pred_x, true_x, genes


def select_eval_timepoints(
    output_mode: Optional[str],
    train_tps: Optional[np.ndarray],
    test_tps: Optional[np.ndarray],
) -> Tuple[List[float], str]:
    """Timepoint *values* to evaluate and the set they come from.

    ``*_tps`` arrays have rows ``[timepoint_index, timepoint_value]``.
    """
    if output_mode == "no_traj":
        return [], "none"
    test_vals = sorted(float(v) for v in np.asarray(test_tps).reshape(-1, 2)[:, 1]) if test_tps is not None and np.size(test_tps) else []
    train_vals = sorted(float(v) for v in np.asarray(train_tps).reshape(-1, 2)[:, 1]) if train_tps is not None and np.size(train_tps) else []
    if output_mode != "full_train" and test_vals:
        return test_vals, "test"
    return train_vals[1:], "train"


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _transform(x: np.ndarray, scale: str) -> np.ndarray:
    if scale == "log":
        return np.log1p(np.clip(x, 0.0, None))
    if scale == "raw":
        return x
    raise CellMetricsError(f"Unknown scale '{scale}'")


def _subsample(x: np.ndarray, max_cells: int, rng: np.random.Generator) -> np.ndarray:
    if len(x) <= max_cells:
        return x
    return x[rng.choice(len(x), size=max_cells, replace=False)]


def wasserstein_w1(a: np.ndarray, b: np.ndarray) -> float:
    """Exact W1 (earth mover's) distance between two point clouds, uniform weights."""
    import ot  # POT
    cost = ot.dist(a, b, metric="euclidean")
    wa = np.full(len(a), 1.0 / len(a))
    wb = np.full(len(b), 1.0 / len(b))
    return float(ot.emd2(wa, wb, cost, numItermax=1_000_000))


def gaussian_js(a: np.ndarray, b: np.ndarray, ridge: float = 1e-6) -> float:
    """Jensen–Shannon divergence between Gaussian fits of *a* and *b*
    (both in the same coordinates), via the moment-matched mixture."""
    k = a.shape[1]

    def stats(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        cov = np.atleast_2d(np.cov(x, rowvar=False, ddof=0)) + np.eye(k) * ridge
        return x.mean(axis=0), cov

    def kl(mu_p, cov_p, mu_q, cov_q) -> float:
        inv_q = np.linalg.pinv(cov_q)
        diff = mu_q - mu_p
        _, ld_p = np.linalg.slogdet(cov_p)
        _, ld_q = np.linalg.slogdet(cov_q)
        return 0.5 * (np.trace(inv_q @ cov_p) + diff @ inv_q @ diff - k + ld_q - ld_p)

    mu_a, cov_a = stats(a)
    mu_b, cov_b = stats(b)
    mu_m = (mu_a + mu_b) / 2
    cov_m = (cov_a + cov_b) / 2 + np.outer(mu_a - mu_b, mu_a - mu_b) / 4
    return float(0.5 * (kl(mu_a, cov_a, mu_m, cov_m) + kl(mu_b, cov_b, mu_m, cov_m)))


def gene_corr_similarity(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    with np.errstate(invalid="ignore", divide="ignore"):
        ca = np.corrcoef(a, rowvar=False)
        cb = np.corrcoef(b, rowvar=False)
    tri = np.triu_indices(a.shape[1], k=1)
    va, vb = ca[tri], cb[tri]
    ok = np.isfinite(va) & np.isfinite(vb)
    if ok.sum() < 2 or np.ptp(va[ok]) == 0 or np.ptp(vb[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(va[ok], vb[ok])[0, 1])


def compute_cell_metrics(
    pred_x: np.ndarray,
    pred_times: np.ndarray,
    true_x: np.ndarray,
    true_times: np.ndarray,
    eval_tps: Sequence[float],
    n_pca_components: int = 50,
    max_cells: int = 2000,
    seed: int = 0,
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """Per-timepoint metrics averaged over *eval_tps*.

    Returns ``(summary, per_timepoint_df)``; gene columns must already be
    aligned (see :func:`align_genes`).
    """
    from sklearn.decomposition import PCA

    if pred_x.shape[1] != true_x.shape[1]:
        raise CellMetricsError(f"Gene mismatch: {pred_x.shape[1]} vs {true_x.shape[1]}")
    rng = np.random.default_rng(seed)
    rows: List[Dict[str, Any]] = []
    for scale in SCALES:
        p_all, t_all = _transform(pred_x, scale), _transform(true_x, scale)
        n_comp = max(1, min(n_pca_components, t_all.shape[1], t_all.shape[0] - 1))
        pca = PCA(n_components=n_comp, random_state=seed).fit(t_all)
        for tp in eval_tps:
            p_tp = p_all[np.isclose(pred_times, tp)]
            t_tp = t_all[np.isclose(true_times, tp)]
            row: Dict[str, Any] = {"timepoint": tp, "scale": scale,
                                   "n_pred_cells": len(p_tp), "n_true_cells": len(t_tp)}
            if len(p_tp) < MIN_CELLS_PER_TP or len(t_tp) < MIN_CELLS_PER_TP:
                row.update({m: np.nan for m in CELL_METRIC_NAMES})
                rows.append(row)
                continue
            p_s, t_s = _subsample(p_tp, max_cells, rng), _subsample(t_tp, max_cells, rng)
            p_emb, t_emb = pca.transform(p_s), pca.transform(t_s)
            row["wasserstein"] = wasserstein_w1(p_emb, t_emb)
            row["js_divergence"] = gaussian_js(p_emb, t_emb)
            row["mean_expr_rmse"] = float(np.sqrt(np.mean((p_tp.mean(0) - t_tp.mean(0)) ** 2)))
            row["gene_corr_similarity"] = gene_corr_similarity(p_tp, t_tp)
            rows.append(row)

    per_tp = pd.DataFrame(rows)
    return summarize_per_timepoint(per_tp, len(eval_tps)), per_tp


def summarize_per_timepoint(per_tp: pd.DataFrame, n_eval_tps: int) -> Dict[str, Any]:
    """Average per-timepoint metrics into ``<metric>_<scale>`` columns."""
    summary: Dict[str, Any] = {"n_eval_timepoints": n_eval_tps}
    if per_tp.empty:
        summary["n_missing_timepoints"] = 0
        summary.update({f"{m}_{s}": np.nan for m in CELL_METRIC_NAMES for s in SCALES})
        return summary
    first = per_tp[per_tp["scale"] == SCALES[0]]
    summary["n_missing_timepoints"] = int((first["n_pred_cells"] < MIN_CELLS_PER_TP).sum())
    for scale in SCALES:
        sub = per_tp[per_tp["scale"] == scale]
        for metric in CELL_METRIC_NAMES:
            vals = sub[metric].to_numpy(dtype=float)
            summary[f"{metric}_{scale}"] = float(np.nanmean(vals)) if np.any(np.isfinite(vals)) else np.nan
    return summary


def impute_missing_worst(per_tp: pd.DataFrame, group_keys: Sequence[str]) -> pd.DataFrame:
    """Give timepoints where a method produced < MIN_CELLS_PER_TP cells the
    worst value observed for that metric among all methods at the same
    (``group_keys`` + timepoint + scale). Adds a boolean ``imputed_worst``
    column. Values stay NaN if no method has a value there.
    """
    df = per_tp.copy()
    missing = df["n_pred_cells"] < MIN_CELLS_PER_TP
    df["imputed_worst"] = False
    if not missing.any():
        return df
    keys = [k for k in group_keys if k in df.columns] + ["timepoint", "scale"]
    for metric, higher_better in HIGHER_IS_BETTER.items():
        grouped = df.groupby(keys, dropna=False)[metric]
        worst = grouped.transform("min" if higher_better else "max")
        df.loc[missing, metric] = worst[missing]
    df.loc[missing, "imputed_worst"] = True
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Compare predicted and true cells per timepoint")
    parser.add_argument("-p", "--pred_cells", required=True, help="Predicted cells .h5ad")
    parser.add_argument("-t", "--true_cells", required=True, help="True cells .h5ad")
    parser.add_argument("--timepoints", required=True,
                        help="Comma-separated timepoint values to evaluate (e.g. 60,70,80)")
    parser.add_argument("-pc", "--n_pca_components", type=int, default=50)
    parser.add_argument("-o", "--output", default=None, help="Output CSV of per-timepoint metrics")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    try:
        px, pv, pt, runtime = load_expression(args.pred_cells)
        tx, tv, tt, _ = load_expression(args.true_cells)
        px, tx, _ = align_genes(px, pv, tx, tv)
        eval_tps = [float(v) for v in args.timepoints.split(",") if v.strip()]
        summary, per_tp = compute_cell_metrics(px, pt, tx, tt, eval_tps, args.n_pca_components)
    except CellMetricsError as e:
        logger.error(str(e))
        sys.exit(1)
    summary["runtime"] = runtime if runtime is not None else np.nan
    print("\n".join(f"{k:28s} {v}" for k, v in summary.items()))
    if args.output:
        per_tp.to_csv(args.output, index=False)
        logger.info("Saved to %s", args.output)


if __name__ == "__main__":
    main()
