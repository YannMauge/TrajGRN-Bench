#!/usr/bin/env python3
"""
metrics_grn.py — Gene Regulatory Network (GRN) evaluation metrics.

All metrics are threshold-free: they score the *ranking* of the predicted
edge weights, so methods whose scores live on different scales (importances,
correlations, kinetic parameters…) are compared fairly and no edge is lost
because it falls just below an arbitrary cut-off.

Ground truth
    Any real-valued matrix ``true[regulator, target]``: an edge exists where
    ``true != 0``, its sign is ``sign(true)`` (> 0 activation, < 0 inhibition)
    and its magnitude is ``|true|``. Binary 0/1, signed (−1/0/+1) and graded
    (signed strengths) truths are therefore all supported. A truth without any
    negative entry is treated as *unsigned* and the sign-aware metrics are NaN.

Edge sets (``scope``), self-loops always excluded
    ``genes``  gene → gene edges only (rows/cols 1..n).
    ``full``   gene → gene plus Stimulus → gene edges (row 0); the stimulus
               column (edges *into* the stimulus) is never scored.
    Index 0 is the stimulus (see METHODS_IO_FORMAT.md). If the prediction is
    an n×n gene-only matrix while the truth is (n+1)×(n+1), only the ``genes``
    scope is computed.

Metrics (per scope, columns ``<scope>_<metric>``)
    edge_auprc          Average precision of |pred| for edge existence.
    edge_auprc_ratio    edge_auprc / edge density (1 = random ranking).
    edge_auroc          ROC AUC of |pred| for edge existence.
    signed_auprc        Mean of AP(activations ranked by +pred) and
                        AP(inhibitions ranked by −pred): high only if strong
                        positive scores are activations and strong negative
                        scores are inhibitions. Unsigned predictions (e.g.
                        importances ≥ 0) cannot recover inhibitions.
    signed_auprc_ratio  Same, each AP divided by its class density.
    sign_accuracy       Fraction of true edges whose predicted sign is right
                        (a zero prediction counts as wrong).
    spearman            Spearman correlation between signed pred and signed
                        truth over all scored entries (uses graded strengths).
    n_true_edges, n_activations, n_inhibitions, n_scored_pairs.

Usage:
    from metrics_grn import compute_grn_metrics, load_grn
    metrics = compute_grn_metrics(load_grn("pred.npy"), load_grn("true.npy"))
"""

from __future__ import annotations

import argparse
import sys
import warnings
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score

try:
    from ._common import setup_logging
except ImportError:
    from _common import setup_logging  # type: ignore[no-redef]

for module in ["anndata"]:
    warnings.filterwarnings("ignore", module=module)

logger = setup_logging(__name__)

SCOPES: Tuple[str, ...] = ("genes", "full")


class GRNMetricsError(Exception):
    """Custom exception for GRN metrics specific errors."""


# ---------------------------------------------------------------------------
# I/O and alignment
# ---------------------------------------------------------------------------

def load_grn(filepath: str) -> np.ndarray:
    try:
        return np.load(filepath).astype(float)
    except Exception as e:
        raise GRNMetricsError(f"Error loading GRN file {filepath}: {e}") from e


def _squeeze_to_2d(matrix: np.ndarray, name: str) -> np.ndarray:
    matrix = np.squeeze(np.asarray(matrix, dtype=float))
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise GRNMetricsError(
            f"GRN matrix '{name}' must be square 2-D after squeezing, got shape {matrix.shape}"
        )
    return matrix


def align_grns(pred_grn: np.ndarray, true_grn: np.ndarray) -> Tuple[np.ndarray, np.ndarray, bool]:
    """Return ``(pred, true, has_stimulus)`` as aligned square matrices.

    Accepted layouts: identical shapes (index 0 = stimulus when the truth has
    one), or an n×n gene-only prediction against an (n+1)×(n+1) truth, in
    which case the stimulus row/column of the truth is dropped and only the
    ``genes`` scope can be scored. Any other mismatch is an error (no silent
    cropping). Non-finite predicted values are set to 0.
    """
    pred = _squeeze_to_2d(pred_grn, "pred_grn")
    true = _squeeze_to_2d(true_grn, "true_grn")
    has_stimulus = True
    if pred.shape != true.shape:
        if pred.shape[0] == true.shape[0] - 1:
            logger.warning(
                f"Predicted GRN {pred.shape} has no stimulus row/column (truth {true.shape}); "
                "scoring gene→gene edges only."
            )
            true = true[1:, 1:]
            has_stimulus = False
        else:
            raise GRNMetricsError(
                f"GRN shape mismatch: pred={pred.shape}, true={true.shape} "
                "(expected identical shapes, or an n×n prediction for an (n+1)×(n+1) truth)"
            )
    n_bad = int(np.sum(~np.isfinite(pred)))
    if n_bad:
        logger.warning(f"Predicted GRN has {n_bad} non-finite value(s); treating them as 0.")
        pred = np.nan_to_num(pred, nan=0.0, posinf=0.0, neginf=0.0)
    return pred, true, has_stimulus


def scope_mask(n: int, scope: str, has_stimulus: bool = True) -> np.ndarray:
    """Boolean mask of the (regulator, target) pairs scored in *scope*."""
    mask = ~np.eye(n, dtype=bool)  # never score self-loops
    if not has_stimulus:
        if scope != "genes":
            raise GRNMetricsError(f"scope '{scope}' needs a stimulus row (index 0)")
        return mask
    mask[:, 0] = False  # edges into the stimulus are meaningless
    if scope == "genes":
        mask[0, :] = False
    elif scope != "full":
        raise GRNMetricsError(f"Unknown scope '{scope}' (expected one of {SCOPES})")
    return mask


# ---------------------------------------------------------------------------
# Metrics on flattened (pred, true) vectors
# ---------------------------------------------------------------------------

def _average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    if labels.sum() == 0:
        return float("nan")
    return float(average_precision_score(labels, scores))


def edge_metrics(pred: np.ndarray, true: np.ndarray) -> Dict[str, float]:
    """Threshold-free metrics on 1-D arrays of predicted / true edge weights."""
    exists = (true != 0).astype(int)
    pos = (true > 0).astype(int)
    neg = (true < 0).astype(int)
    n_pairs = int(len(true))
    out: Dict[str, float] = {
        "n_scored_pairs": n_pairs,
        "n_true_edges": int(exists.sum()),
        "n_activations": int(pos.sum()),
        "n_inhibitions": int(neg.sum()),
    }

    # Edge existence (sign ignored)
    abs_pred = np.abs(pred)
    density = exists.mean() if n_pairs else float("nan")
    out["edge_auprc"] = _average_precision(exists, abs_pred)
    out["edge_auprc_ratio"] = out["edge_auprc"] / density if density else float("nan")
    out["edge_auroc"] = (
        float(roc_auc_score(exists, abs_pred)) if 0 < exists.sum() < n_pairs else float("nan")
    )

    # Sign-aware metrics (only if the truth is signed)
    if neg.sum() > 0:
        ap_act = _average_precision(pos, pred)
        ap_inh = _average_precision(neg, -pred)
        aps = [a for a in (ap_act, ap_inh) if np.isfinite(a)]
        ratios = [
            a / d for a, d in ((ap_act, pos.mean()), (ap_inh, neg.mean()))
            if np.isfinite(a) and d > 0
        ]
        out["signed_auprc"] = float(np.mean(aps)) if aps else float("nan")
        out["signed_auprc_ratio"] = float(np.mean(ratios)) if ratios else float("nan")
        out["signed_auprc_activation"] = ap_act
        out["signed_auprc_inhibition"] = ap_inh
        true_edges = exists.astype(bool)
        out["sign_accuracy"] = float(np.mean(np.sign(pred[true_edges]) == np.sign(true[true_edges])))
    else:
        for key in ("signed_auprc", "signed_auprc_ratio", "signed_auprc_activation",
                    "signed_auprc_inhibition", "sign_accuracy"):
            out[key] = float("nan")

    # Graded, signed agreement
    if np.ptp(pred) == 0 or np.ptp(true) == 0:
        out["spearman"] = float("nan")
    else:
        from scipy.stats import spearmanr
        out["spearman"] = float(spearmanr(pred, true)[0])
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_grn_metrics(pred_grn: np.ndarray, true_grn: np.ndarray) -> Dict[str, Any]:
    """Compute all GRN metrics for every scope; keys are ``<scope>_<metric>``."""
    pred, true, has_stimulus = align_grns(pred_grn, true_grn)
    if not np.any(true[~np.eye(true.shape[0], dtype=bool)] < 0):
        logger.info("Ground-truth GRN has no negative entries: treated as unsigned (sign metrics = NaN).")
    metrics: Dict[str, Any] = {}
    for scope in SCOPES:
        if scope != "genes" and not has_stimulus:
            continue
        mask = scope_mask(true.shape[0], scope, has_stimulus)
        for key, value in edge_metrics(pred[mask], true[mask]).items():
            metrics[f"{scope}_{key}"] = value
    return metrics


def compute_precision_recall_curve(
    pred_grn: np.ndarray, true_grn: np.ndarray, scope: str = "genes",
) -> Dict[str, Any]:
    """Precision-recall curve of |pred| for edge existence (same edge set and
    scoring as ``<scope>_edge_auprc``)."""
    pred, true, has_stimulus = align_grns(pred_grn, true_grn)
    if not has_stimulus:
        scope = "genes"
    mask = scope_mask(true.shape[0], scope, has_stimulus)
    labels = (true[mask] != 0).astype(int)
    scores = np.abs(pred[mask])
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    return {
        "precision": precision,
        "recall": recall,
        "thresholds": thresholds,
        "auprc": _average_precision(labels, scores),
        "n_positives": int(labels.sum()),
        "n_pairs": int(len(labels)),
    }


def evaluate_grn_similarity(
    pred_grn: np.ndarray, true_grn: np.ndarray
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    metrics = compute_grn_metrics(pred_grn, true_grn)
    return metrics, pd.DataFrame([metrics])


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Compute metrics for evaluating a predicted GRN")
    parser.add_argument("-p", "--pred_grn", required=True, help="Predicted GRN (.npy)")
    parser.add_argument("-t", "--true_grn", required=True, help="True GRN (.npy)")
    parser.add_argument("-o", "--output", default=None, help="Optional output CSV")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    try:
        metrics, metrics_df = evaluate_grn_similarity(load_grn(args.pred_grn), load_grn(args.true_grn))
    except GRNMetricsError as e:
        logger.error(str(e))
        sys.exit(1)
    print("\n".join(f"{k:40s} {v:.4f}" if isinstance(v, float) else f"{k:40s} {v}"
                    for k, v in metrics.items()))
    if args.output:
        metrics_df.to_csv(args.output, index=False)
        logger.info(f"Metrics saved to: {args.output}")


if __name__ == "__main__":
    main()
