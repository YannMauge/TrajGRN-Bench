#!/usr/bin/env python3
"""Metrics for evaluating cell population similarity.

Provides functions to compute Wasserstein distance, KL/Jensen-Shannon
divergence, ARI, gene-gene correlation, and top-DE-gene overlap between
predicted and true cell populations.
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

try:
    from ._common import (
        ADATA_FILE_RE,
        PERTURBATION_TRAINING_SUFFIX,
        TRUE_ADATA_RE,
        discover_replicates,
        parse_ko_output_genes,
        setup_logging,
    )
except ImportError:
    from _common import (  # type: ignore[no-redef]
        ADATA_FILE_RE,
        PERTURBATION_TRAINING_SUFFIX,
        TRUE_ADATA_RE,
        discover_replicates,
        parse_ko_output_genes,
        setup_logging,
    )

warnings.filterwarnings("ignore", module="anndata")

logger = setup_logging(__name__)


class CellMetricsError(Exception):
    """Raised when cell metric computation fails."""


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _validate_embeddings(e1: np.ndarray, e2: np.ndarray) -> None:
    if e1 is None or e2 is None:
        raise CellMetricsError("Both embeddings must be provided")
    if e1.shape[1] != e2.shape[1]:
        raise CellMetricsError(f"Feature mismatch: {e1.shape[1]} vs {e2.shape[1]}")


def _validate_labels(l1: np.ndarray, l2: np.ndarray, n1: int, n2: int) -> None:
    if l1 is None or l2 is None:
        raise CellMetricsError("Both label arrays must be provided")
    if len(l1) != n1 or len(l2) != n2:
        raise CellMetricsError(f"Label/cell mismatch: {len(l1)}/{n1} vs {len(l2)}/{n2}")


def _is_simulated_adata(adata: Any) -> bool:
    v = adata.uns.get("simulation")
    return v is not None and not (isinstance(v, str) and v.strip().lower() in {"", "none", "false"})


def _resimulate_ko_from_wt(adata: Any, dataset_id: str, ko_gene: Optional[str]) -> Optional[Any]:
    if dataset_id == "WT" or not ko_gene or "dataset_id" not in adata.obs.columns:
        return None
    if not _is_simulated_adata(adata) or ko_gene not in adata.var_names:
        return None
    wt_mask = adata.obs["dataset_id"].astype(str) == "WT"
    if int(wt_mask.sum()) == 0:
        return None
    adata_ko = adata[wt_mask].copy()
    ko_idx = int(list(adata_ko.var_names).index(ko_gene))
    if sparse.issparse(adata_ko.X):
        x = adata_ko.X.tolil(copy=True)
        x[:, ko_idx] = 0
        adata_ko.X = x.tocsr()
    else:
        x = np.asarray(adata_ko.X).copy()
        x[:, ko_idx] = 0
        adata_ko.X = x
    adata_ko.obs["dataset_id"] = dataset_id
    adata_ko.uns["resimulated_from_wt"] = True
    adata_ko.uns["resimulated_ko_gene"] = ko_gene
    return adata_ko


def _safe_inv(cov: np.ndarray) -> np.ndarray:
    try:
        return np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(cov)


def _gaussian_stats(emb: np.ndarray, nc: int) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    from sklearn.decomposition import PCA
    nc = min(nc, emb.shape[1])
    ep = PCA(n_components=nc).fit_transform(emb)
    mu = ep.mean(axis=0)
    cov = np.cov(ep.T, ddof=0)
    if cov.ndim == 1:
        cov = cov.reshape(1, 1)
    cov += np.eye(cov.shape[0]) * 1e-6
    _, ld = np.linalg.slogdet(cov)
    return ep, mu, ld, cov


def _kl_gauss(mu_a, cov_a, ld_a, mu_b, cov_b_inv, ld_b, k: int) -> float:
    return 0.5 * (np.trace(cov_b_inv @ cov_a) + (mu_b - mu_a).T @ cov_b_inv @ (mu_b - mu_a) - k + ld_b - ld_a)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_adata(
    filepath: str,
    dataset_id: str = "WT",
    resimulate_missing_ko: bool = False,
    ko_gene: Optional[str] = None,
) -> Any:
    """Read an .h5ad file and optionally subset to *dataset_id*.

    If *dataset_id* is not found and *resimulate_missing_ko* is True,
    attempt to synthesise a KO dataset from the WT population.
    Falls back to the full dataset with a warning when subsetting fails.

    Returns
    -------
    anndata.AnnData
        The (possibly subsetted) AnnData object.
    """
    import anndata as ad  # type: ignore
    try:
        adata = ad.read_h5ad(filepath)
    except (ad.core.AccessorError, ad.core.ReadError) as e:
        raise CellMetricsError(str(e)) from e

    if "dataset_id" not in adata.obs.columns or not dataset_id:
        return adata

    mask = adata.obs["dataset_id"].astype(str) == str(dataset_id)
    if int(mask.sum()) > 0:
        return adata[mask].copy()

    if resimulate_missing_ko:
        resim = _resimulate_ko_from_wt(adata, str(dataset_id), ko_gene)
        if resim is not None:
            return resim

    logger.warning(
        "Dataset '%s' not found in %s; using full dataset.", dataset_id, filepath
    )
    return adata


def load_adata(
    filepath: str,
    dataset_id: str = "WT",
    resimulate_missing_ko: bool = False,
    ko_gene: Optional[str] = None,
) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[float]]:
    """Load expression matrix, labels and runtime from an .h5ad file.

    Delegates to :func:`read_adata` for I/O and subsetting.
    """
    adata = read_adata(filepath, dataset_id, resimulate_missing_ko, ko_gene)

    if adata.X is None:
        raise CellMetricsError(f"No expression matrix in: {filepath}")
    emb = np.asarray(adata.X).astype(float)
    if emb.shape[0] == 0:
        raise CellMetricsError(f"Empty embeddings in: {filepath}")

    runtime: Optional[float] = None
    if "runtime" in adata.uns:
        try:
            runtime = float(adata.uns["runtime"])
        except (TypeError, ValueError):
            logger.warning("Bad runtime value in %s; skipping.", filepath)

    if "cell_type" not in adata.obs:
        logger.warning("No cell_type column in %s; skipping label metrics.", filepath)
        return emb, None, runtime

    labels = adata.obs["cell_type"].astype(int).values
    if len(labels) == 0:
        logger.warning("Empty labels in %s; skipping label metrics.", filepath)
        return emb, None, runtime
    if len(labels) != emb.shape[0]:
        raise CellMetricsError(
            f"Label/embedding mismatch: {len(labels)} vs {emb.shape[0]}"
        )
    return emb, labels, runtime


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def compute_wasserstein(e1: np.ndarray, e2: np.ndarray, n_components: Optional[int] = None) -> float:
    _validate_embeddings(e1, e2)
    if n_components is not None:
        if n_components < 1:
            raise CellMetricsError("n_components must be >= 1")
        from sklearn.decomposition import PCA
        pca = PCA(n_components=min(n_components, e1.shape[1]))
        e1, e2 = pca.fit_transform(e1), pca.transform(e2)
    from scipy.optimize import linear_sum_assignment  # type: ignore
    from scipy.spatial.distance import cdist  # type: ignore
    cost = cdist(e1, e2, metric="euclidean")
    ri, ci = linear_sum_assignment(cost)
    return float(cost[ri, ci].sum() / min(len(e1), len(e2)))


def compute_kl(e1: np.ndarray, e2: np.ndarray, n_components: int = 50) -> float:
    _validate_embeddings(e1, e2)
    ep1, mu1, ld1, cov1 = _gaussian_stats(e1, n_components)
    _, mu2, ld2, cov2 = _gaussian_stats(e2, n_components)
    return float(_kl_gauss(mu1, cov1, ld1, mu2, _safe_inv(cov2), ld2, ep1.shape[1]))


def compute_js(e1: np.ndarray, e2: np.ndarray, n_components: int = 50) -> float:
    _validate_embeddings(e1, e2)
    ep1, mu1, ld1, cov1 = _gaussian_stats(e1, n_components)
    _, mu2, ld2, cov2 = _gaussian_stats(e2, n_components)
    mu_m = (mu1 + mu2) / 2
    cov_m = (cov1 + cov2) / 2 + np.eye(cov1.shape[0]) * 1e-6
    _, ld_m = np.linalg.slogdet(cov_m)
    inv_m = _safe_inv(cov_m)
    k = ep1.shape[1]
    kl_pm = _kl_gauss(mu1, cov1, ld1, mu_m, inv_m, ld_m, k)
    kl_qm = _kl_gauss(mu2, cov2, ld2, mu_m, inv_m, ld_m, k)
    return float(0.5 * (kl_pm + kl_qm))


def compute_ari(l1: np.ndarray, l2: np.ndarray) -> float:
    _validate_labels(l1, l2, len(l1), len(l2))
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(l1, l2))


def compute_gene_corr(e1: np.ndarray, e2: np.ndarray) -> float:
    _validate_embeddings(e1, e2)
    tri = np.triu_indices_from(np.corrcoef(e1, rowvar=False), k=1)
    c1 = np.corrcoef(e1, rowvar=False)[tri]
    c2 = np.corrcoef(e2, rowvar=False)[tri]
    mask = np.isfinite(c1) & np.isfinite(c2)
    return float(np.corrcoef(c1[mask], c2[mask])[0, 1])


def compute_top_de_overlap(e1: np.ndarray, e2: np.ndarray, top_k: int = 100) -> float:
    _validate_embeddings(e1, e2)
    n = e1.shape[1]
    if n == 0:
        raise CellMetricsError("No features")
    k = max(1, min(int(top_k), n))
    pm = np.nanmean(np.vstack([e1, e2]), axis=0)
    s1 = np.nan_to_num(np.abs(np.nanmean(e1, axis=0) - pm), nan=0.0, posinf=0.0, neginf=0.0)
    s2 = np.nan_to_num(np.abs(np.nanmean(e2, axis=0) - pm), nan=0.0, posinf=0.0, neginf=0.0)
    return float(len(set(np.argpartition(s1, -k)[-k:]) & set(np.argpartition(s2, -k)[-k:])) / k)


# ---------------------------------------------------------------------------
# high-level
# ---------------------------------------------------------------------------

def compute_cell_metrics(
    pred_cells: np.ndarray, true_cells: np.ndarray,
    pred_labels: Optional[np.ndarray] = None, true_labels: Optional[np.ndarray] = None,
    n_pca_components: int = 50, runtime: Optional[float] = None,
) -> Dict[str, Any]:
    _validate_embeddings(pred_cells, true_cells)
    if pred_labels is not None and len(pred_labels) != len(pred_cells):
        raise CellMetricsError(f"Pred label/cell mismatch: {len(pred_labels)} vs {len(pred_cells)}")
    if true_labels is not None and len(true_labels) != len(true_cells):
        raise CellMetricsError(f"True label/cell mismatch: {len(true_labels)} vs {len(true_cells)}")

    ari = compute_ari(pred_labels, true_labels) if (pred_labels is not None and true_labels is not None) else np.nan

    metrics = {
        "wasserstein_distance": compute_wasserstein(pred_cells, true_cells),
        "pca_wasserstein_distance": compute_wasserstein(pred_cells, true_cells, n_components=50),
        "kl_divergence": compute_kl(pred_cells, true_cells, n_components=n_pca_components),
        "jensen_shannon_divergence": compute_js(pred_cells, true_cells, n_components=n_pca_components),
        "gene_gene_correlation": compute_gene_corr(pred_cells, true_cells),
        "top_de_genes_overlap": compute_top_de_overlap(pred_cells, true_cells),
        "adjusted_rand_index": ari,
        "runtime": runtime if runtime is not None else np.nan,
    }
    return metrics


def evaluate_cell_similarity(
    pred_cells: np.ndarray, true_cells: np.ndarray,
    pred_labels: Optional[np.ndarray] = None, true_labels: Optional[np.ndarray] = None,
    n_pca_components: int = 50, runtime: Optional[float] = None,
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    m = compute_cell_metrics(pred_cells, true_cells, pred_labels, true_labels, n_pca_components, runtime)
    return m, pd.DataFrame([m])


# ---------------------------------------------------------------------------
# cli
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(description="Compute cell population similarity metrics")
    parser.add_argument("-p", "--pred_cells", required=True, help="Predicted cells .h5ad")
    parser.add_argument("-t", "--true_cells", required=True, help="True cells .h5ad")
    parser.add_argument("-pc", "--n_pca_components", type=int, default=50)
    parser.add_argument("-o", "--output", default=None, help="Output CSV (optional)")
    args = parser.parse_args(argv)

    try:
        logger.info("Loading predicted: %s", args.pred_cells)
        pc, pl, pr = load_adata(args.pred_cells)
        logger.info("Loading true: %s", args.true_cells)
        tc, tl, _ = load_adata(args.true_cells)

        metrics, df = evaluate_cell_similarity(pc, tc, pl, tl, n_pca_components=args.n_pca_components, runtime=pr)

        print("\n" + "=" * 60)
        print("Cell Population Evaluation Metrics")
        print("=" * 60)
        for key in ["wasserstein_distance", "pca_wasserstein_distance", "kl_divergence",
                     "jensen_shannon_divergence", "adjusted_rand_index",
                     "gene_gene_correlation", "top_de_genes_overlap"]:
            print(f"{key}: {metrics[key]:.4f}")
        if not np.isnan(metrics["runtime"]):
            print(f"runtime: {metrics['runtime']:.4f} s")
        print("=" * 60)

        if args.output:
            df.to_csv(args.output, index=False)
            logger.info("Saved to %s", args.output)
    except CellMetricsError as e:
        logger.error("CellMetricsError: %s", e)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Interrupted.")
        sys.exit(1)
    except Exception as e:
        logger.error("Error: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
