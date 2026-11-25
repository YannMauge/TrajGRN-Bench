"""Shared helpers for method wrappers (TrajGRN-Bench I/O contract).

Kept dependency-light (numpy, pandas, anndata, scipy) and Python 3.8
compatible, because it is imported from every method environment.

Contract (see METHODS_IO_FORMAT.md):
  * input  ``-i data_<i>.h5ad`` in the run's ``inputs/`` folder, next to
    ``train_tps.npy`` / ``test_tps.npy`` (rows ``[index, value]``) and, for
    ``subsample_full``, ``subsample_train_ids.npy``;
  * GRN output ``<data_id>[_ko_<gene>]_GRN.npy``: 2-D, shape
    (n_genes+1, n_genes+1), row/column 0 = Stimulus, ``grn[regulator, target]``.

Usage from a wrapper::

    sys.path.insert(0, str(Path(__file__).resolve().parents[N]))  # methods/
    from _common import load_training_data, save_grn, seed_everything
"""

from __future__ import annotations

import os
import random
from typing import List, Optional, Tuple

import numpy as np

T_KEY = "timepoint"
STIMULUS_PREFIX = "stimulus"


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def seed_everything(seed: int = 0) -> None:
    """Seed python, numpy and (if installed) torch / jax-free libraries."""
    random.seed(seed)
    np.random.seed(seed)
    os.environ.setdefault("PYTHONHASHSEED", str(seed))
    try:
        import torch  # type: ignore
        torch.manual_seed(seed)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def data_id(input_file: str) -> str:
    """``data_<i>`` from the input path."""
    return os.path.splitext(os.path.basename(input_file))[0]


def load_split(input_file: str) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """``(train_tps, test_tps, subsample_ids)`` saved next to *input_file*."""
    folder = os.path.dirname(os.path.abspath(input_file))
    train = np.load(os.path.join(folder, "train_tps.npy"))
    test_path = os.path.join(folder, "test_tps.npy")
    test = np.load(test_path) if os.path.exists(test_path) else np.empty((0, 2))
    sub_path = os.path.join(folder, "subsample_train_ids.npy")
    sub = np.load(sub_path, allow_pickle=True) if os.path.exists(sub_path) else None
    return np.asarray(train).reshape(-1, 2), np.asarray(test).reshape(-1, 2), sub


def is_stimulus(name: str) -> bool:
    return str(name).lower().startswith(STIMULUS_PREFIX)


def read_dense(input_file: str):
    """Read the AnnData with a dense float ``.X``."""
    import anndata as ad
    from scipy import sparse
    adata = ad.read_h5ad(input_file)
    if sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.asarray(adata.X, dtype=float)
    return adata


def load_training_data(
    input_file: str,
    perturbation_training: bool = False,
    use_subsample: bool = True,
):
    """Cells a method may train on.

    * WT cells only unless *perturbation_training* (then WT + the training
      perturbations present in the input file);
    * training timepoints only (``obs['timepoint']`` in ``train_tps``);
    * the ``subsample_full`` subset when ``subsample_train_ids.npy`` exists
      (matched on ``obs_names`` *before* any re-indexing).

    Returns ``(adata_train, train_tps, test_tps)``.
    """
    adata = read_dense(input_file)
    train_tps, test_tps, sub = load_split(input_file)
    mask = np.isin(adata.obs[T_KEY].to_numpy(dtype=float), train_tps[:, 1].astype(float))
    if "dataset_id" in adata.obs.columns and not perturbation_training:
        mask &= (adata.obs["dataset_id"].astype(str) == "WT").to_numpy()
    if use_subsample and sub is not None:
        mask &= np.isin(adata.obs_names.astype(str), np.asarray(sub).astype(str))
    if not mask.any():
        raise ValueError(f"No training cells left in {input_file} after filtering")
    return adata[mask].copy(), train_tps, test_tps


def gene_indices(var_names) -> List[int]:
    """Indices of genes (non-stimulus variables)."""
    return [i for i, v in enumerate(var_names) if not is_stimulus(v)]


def log_normalize_genes(x: np.ndarray, var_names, target_sum: Optional[float] = None) -> np.ndarray:
    """``normalize_total`` + ``log1p`` computed on gene columns only (the
    stimulus column is not counted in the library size); stimulus columns are
    returned as ``log1p`` of their raw value."""
    x = np.asarray(x, dtype=float).copy()
    genes = gene_indices(var_names)
    lib = x[:, genes].sum(axis=1, keepdims=True)
    target = float(np.median(lib[lib > 0])) if target_sum is None and np.any(lib > 0) else (target_sum or 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        x[:, genes] = np.where(lib > 0, x[:, genes] / lib * target, 0.0)
    return np.log1p(np.clip(x, 0.0, None))


def variable_columns(x: np.ndarray, rtol: float = 1e-8) -> np.ndarray:
    """Indices of columns that are not constant (up to float rounding)."""
    x = np.asarray(x, dtype=float)
    scale = 1.0 + np.abs(x).mean(axis=0)
    return np.flatnonzero(np.std(x, axis=0) > rtol * scale)


def run_on_variable_columns(func, x: np.ndarray, *args, **kwargs) -> np.ndarray:
    """Run a GRN method ``func(x_sub, ...) -> (k, k)`` on the non-constant
    columns of *x* only (tree methods such as GENIE3 divide each target by its
    standard deviation) and embed the result in a zero (n, n) matrix."""
    x = np.asarray(x, dtype=float)
    keep = variable_columns(x)
    full = np.zeros((x.shape[1], x.shape[1]))
    if len(keep) >= 2:
        full[np.ix_(keep, keep)] = func(x[:, keep], *args, **kwargs)
    return full


def parse_ko_targets(value: Optional[str], var_names) -> List[str]:
    """KO genes to predict from ``-k``: "none", "all" or comma-separated names."""
    v = (value or "none").strip()
    genes = [str(g) for g in var_names if not is_stimulus(g)]
    if v.lower() in ("", "none"):
        return []
    if v.lower() == "all":
        return genes
    out = [g.strip() for g in v.split(",") if g.strip()]
    unknown = [g for g in out if g not in genes]
    if unknown:
        raise ValueError(f"Unknown KO genes {unknown}; available: {genes}")
    return out


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

def save_grn(
    output_dir: str,
    input_file: str,
    grn: np.ndarray,
    var_names,
    suffix: str = "",
) -> str:
    """Save ``grn[regulator, target]`` following the benchmark contract.

    *grn* is indexed like *var_names* (which must include the stimulus at
    index 0 for a full-size matrix) or like the genes only, in which case a
    zero stimulus row/column is added. Non-finite values become 0. Returns
    the saved path.
    """
    grn = np.asarray(grn, dtype=float)
    grn = np.squeeze(grn)
    names = [str(v) for v in var_names]
    n_all, n_genes = len(names), len(gene_indices(names))
    if grn.ndim != 2 or grn.shape[0] != grn.shape[1]:
        raise ValueError(f"GRN must be a square matrix, got shape {grn.shape}")
    if grn.shape[0] == n_genes and n_all == n_genes + 1:
        full = np.zeros((n_all, n_all))
        g = gene_indices(names)
        full[np.ix_(g, g)] = grn
        grn = full
    has_stimulus = n_genes < n_all
    if grn.shape[0] != n_all or (has_stimulus and not is_stimulus(names[0])):
        raise ValueError(
            f"GRN shape {grn.shape} does not match variables {n_all} with stimulus first ({names[:2]}…)"
        )
    grn = np.nan_to_num(grn, nan=0.0, posinf=0.0, neginf=0.0)
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"{data_id(input_file)}{suffix}_GRN.npy")
    np.save(path, grn)
    return path
