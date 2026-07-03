#!/usr/bin/env python3
"""
metrics_grn.py — Gene Regulatory Network (GRN) evaluation metrics.

Provides functions to evaluate the similarity/distance between predicted and
true GRNs: precision, recall, F1-score, AUPRC, Jaccard index, Hamming
distance, Spearman correlation, and edge overlap.

Usage:
    from metrics_grn import compute_grn_metrics, load_grn

    pred_grn = load_grn("path/to/predicted_grn.npy")
    true_grn = load_grn("path/to/true_grn.npy")
    metrics = compute_grn_metrics(pred_grn, true_grn)
"""

from __future__ import annotations

import logging
import sys
import warnings
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

try:
    from ._common import setup_logging
except ImportError:
    from _common import setup_logging  # type: ignore[no-redef]

for module in ["anndata"]:
    warnings.filterwarnings("ignore", module=module)

logger = setup_logging(__name__)


class GRNMetricsError(Exception):
    """Custom exception for GRN metrics specific errors."""

    pass


def load_grn(filepath: str) -> np.ndarray:
    try:
        grn = np.load(filepath).astype(float)
        return grn
    except Exception as e:
        logger.error(f"Error loading GRN file: {e}")
        raise GRNMetricsError(f"Error loading GRN file: {e}")


def _binarize_and_mask_diag(grn: np.ndarray, threshold: float) -> np.ndarray:
    """Convert to boolean via threshold and exclude self-loops (diagonal).

    A strict comparison is used so that the default threshold of 0.0 does not
    classify exact zero-weight edges as positives.
    """
    binary = np.abs(grn) > threshold
    binary[np.diag_indices_from(binary)] = False
    return binary


def compute_precision_recall_f1(
    pred_grn: np.ndarray, true_grn: np.ndarray, threshold: float = 0.0
) -> Dict[str, float]:
    try:
        pred_binary = _binarize_and_mask_diag(pred_grn, threshold)
        true_binary = _binarize_and_mask_diag(true_grn, threshold)

        # Compute confusion matrix elements
        tp = np.sum(pred_binary & true_binary)  # True positives
        fp = np.sum(pred_binary & ~true_binary)  # False positives
        fn = np.sum(~pred_binary & true_binary)  # False negatives
        tn = np.sum(~pred_binary & ~true_binary)  # True negatives

        # Compute metrics
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1_score = (
            2 * precision * recall / (precision + recall)
            if (precision + recall) > 0
            else 0.0
        )

        return {
            "precision": precision,
            "recall": recall,
            "f1_score": f1_score,
            "true_positives": tp,
            "false_positives": fp,
            "false_negatives": fn,
            "true_negatives": tn,
        }

    except Exception as e:
        logger.error(f"Error computing precision/recall/F1: {e}")
        raise GRNMetricsError(f"Error computing precision/recall/F1: {e}")

def compute_jaccard_index(
    pred_grn: np.ndarray, true_grn: np.ndarray, threshold: float = 0.0
) -> float:
    try:
        pred_binary = _binarize_and_mask_diag(pred_grn, threshold)
        true_binary = _binarize_and_mask_diag(true_grn, threshold)

        intersection = np.sum(pred_binary & true_binary)
        union = np.sum(pred_binary | true_binary)
        return float(intersection / union) if union > 0 else 0.0

    except Exception as e:
        logger.error(f"Error computing Jaccard index: {e}")
        raise GRNMetricsError(f"Error computing Jaccard index: {e}")


def compute_hamming_distance(
    pred_grn: np.ndarray, true_grn: np.ndarray, threshold: float = 0.0
) -> float:
    try:
        pred_binary = _binarize_and_mask_diag(pred_grn, threshold)
        true_binary = _binarize_and_mask_diag(true_grn, threshold)
        return float(np.sum(pred_binary != true_binary))

    except Exception as e:
        logger.error(f"Error computing Hamming distance: {e}")
        raise GRNMetricsError(f"Error computing Hamming distance: {e}")


def compute_spearman_correlation(pred_grn: np.ndarray, true_grn: np.ndarray) -> float:
    try:
        from scipy.stats import spearmanr  # type: ignore

        off_diag = ~np.eye(pred_grn.shape[0], dtype=bool)
        pred_flat = pred_grn[off_diag]
        true_flat = true_grn[off_diag]
        correlation, _ = spearmanr(pred_flat, true_flat)
        return float(correlation)

    except ImportError:
        logger.warning("scipy not available, using numpy correlation")
        return float(np.corrcoef(pred_flat, true_flat)[0, 1])  # type: ignore

    except Exception as e:
        logger.error(f"Error computing Spearman correlation: {e}")
        raise GRNMetricsError(f"Error computing Spearman correlation: {e}")


def compute_edge_overlap(
    pred_grn: np.ndarray, true_grn: np.ndarray, threshold: float = 0.0
) -> Dict[str, float]:
    try:
        pred_binary = _binarize_and_mask_diag(pred_grn, threshold)
        true_binary = _binarize_and_mask_diag(true_grn, threshold)

        n_edges_pred = int(np.sum(pred_binary))
        n_edges_true = int(np.sum(true_binary))
        n_edges_overlap = int(np.sum(pred_binary & true_binary))
        max_edges = max(n_edges_pred, n_edges_true)

        return {
            "n_edges_predicted": n_edges_pred,
            "n_edges_true": n_edges_true,
            "n_edges_overlap": n_edges_overlap,
            "edge_overlap_ratio": float(n_edges_overlap / max_edges) if max_edges > 0 else 0.0,
        }

    except Exception as e:
        logger.error(f"Error computing edge overlap: {e}")
        raise GRNMetricsError(f"Error computing edge overlap: {e}")


def compute_auprc(pred_grn: np.ndarray, true_grn: np.ndarray) -> float:
    try:
        n = pred_grn.shape[0]
        off_diag = ~np.eye(n, dtype=bool)
        pred_scores = pred_grn[off_diag]
        true_labels = true_grn[off_diag].astype(int)
        return float(average_precision_score(true_labels, pred_scores))

    except Exception as e:
        logger.error(f"Error computing AUPRC: {e}")
        raise GRNMetricsError(f"Error computing AUPRC: {e}")


def compute_precision_recall_curve(
    pred_grn: np.ndarray,
    true_grn: np.ndarray,
    threshold: float = 0.0,
) -> Dict[str, np.ndarray]:
    """Compute precision-recall curve data for a predicted GRN against a
    ground-truth GRN.

    The ground-truth GRN is binarized using *threshold* (edges with absolute
    weight ≥ threshold are considered present, excluding self-loops).  The
    predicted GRN is used as a continuous score.  Precision and recall are
    computed at all operating points via :func:`sklearn.metrics.precision_recall_curve`.

    Args:
        pred_grn:  Predicted GRN matrix (continuous scores).
        true_grn:  Ground-truth GRN matrix (continuous or binary weights).
        threshold: Binarization threshold for the true GRN.

    Returns:
        Dict with keys ``"precision"``, ``"recall"``, ``"thresholds"`` (arrays),
        and ``"auprc"`` (float).
    """
    from sklearn.metrics import precision_recall_curve

    try:
        n = pred_grn.shape[0]
        off_diag = ~np.eye(n, dtype=bool)
        pred_scores = np.abs(pred_grn[off_diag])
        true_binary = (np.abs(true_grn[off_diag]) > threshold).astype(int)

        precision, recall, thresholds = precision_recall_curve(
            true_binary, pred_scores
        )
        auprc = float(average_precision_score(true_binary, pred_scores))

        return {
            "precision": precision,
            "recall": recall,
            "thresholds": thresholds,
            "auprc": auprc,
        }

    except Exception as e:
        logger.error(f"Error computing precision-recall curve: {e}")
        raise GRNMetricsError(f"Error computing precision-recall curve: {e}")


def _squeeze_to_2d(matrix: np.ndarray, name: str) -> np.ndarray:
    matrix = np.squeeze(matrix)
    if matrix.ndim != 2:
        raise GRNMetricsError(
            f"GRN matrix '{name}' must be 2-D after squeezing, got shape {matrix.shape}"
        )
    return matrix


def compute_grn_metrics(
    pred_grn: np.ndarray, true_grn: np.ndarray, simulation: str , quantile: float = 0.9
) -> Dict[str, Any]:
    try:
        # Validate inputs
        if pred_grn is None or true_grn is None:
            raise GRNMetricsError("Both pred_grn and true_grn must be provided")

        # Squeeze trailing size-1 dimensions (e.g. (9,9,1) -> (9,9))
        pred_grn = _squeeze_to_2d(pred_grn, "pred_grn")
        true_grn = _squeeze_to_2d(true_grn, "true_grn")

        if pred_grn.shape != true_grn.shape:
            logger.warning(
                f"GRN matrices have different shapes: pred={pred_grn.shape}, true={true_grn.shape}. "
                "Using min shape for computation."
            )
            min_n = min(pred_grn.shape[0], true_grn.shape[0])
            min_m = min(pred_grn.shape[1], true_grn.shape[1])
            pred_grn = pred_grn[:min_n, :min_m]
            true_grn = true_grn[:min_n, :min_m]

        # Calculate threshold from quantile if provided
        if not (0 <= quantile <= 1):
            raise GRNMetricsError("quantile must be between 0 and 1")
        # Get absolute values of true GRN (excluding diagonal)
        true_grn_copy = true_grn.copy()
        true_grn_copy[np.diag_indices_from(true_grn_copy)] = 0
        true_grn_abs = np.abs(true_grn_copy)
        threshold = float(np.quantile(true_grn_abs, quantile))
        logger.info(f"Calculated threshold from quantile {quantile}: {threshold}")

        # Compute all metrics
        precision_recall_f1 = compute_precision_recall_f1(pred_grn, true_grn, threshold)
        jaccard = compute_jaccard_index(pred_grn, true_grn, threshold)
        hamming = compute_hamming_distance(pred_grn, true_grn, threshold)
        spearman = compute_spearman_correlation(pred_grn, true_grn)
        auprc = compute_auprc(pred_grn, true_grn)
        edge_overlap = compute_edge_overlap(pred_grn, true_grn, threshold)

        # Compile results
        metrics = {
            **precision_recall_f1,
            "jaccard_index": jaccard,
            "hamming_distance": hamming,
            "spearman_correlation": spearman,
            "auprc": auprc,
            "edge_overlap": edge_overlap,
        }

        return metrics

    except Exception as e:
        logger.error(f"Error computing GRN metrics: {e}")
        raise GRNMetricsError(f"Error computing GRN metrics: {e}")


def evaluate_grn_similarity(
    pred_grn: np.ndarray, true_grn: np.ndarray, quantile: float = 0.9
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    try:
        # Compute metrics
        metrics = compute_grn_metrics(pred_grn, true_grn, quantile=quantile)

        # Create DataFrame for metrics
        metrics_df = pd.DataFrame([metrics])

        return metrics, metrics_df

    except Exception as e:
        logger.error(f"Error evaluating GRN similarity: {e}")
        raise GRNMetricsError(f"Error evaluating GRN similarity: {e}")


def main(argv: Optional[List[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    try:
        # Parse arguments
        import argparse

        parser = argparse.ArgumentParser(
            description="Compute metrics for evaluating GRN predictions"
        )
        parser.add_argument(
            "-p",
            "--pred_grn",
            required=True,
            help="Path to predicted GRN file (h5ad format)",
        )
        parser.add_argument(
            "-t",
            "--true_grn",
            required=True,
            help="Path to true GRN file (h5ad format)",
        )
        parser.add_argument(
            "-q",
            "--quantile",
            type=float,
            default=0.9,
            help="Quantile for threshold calculation from true GRN edges (0-1). If provided, overrides --threshold",
        )
        parser.add_argument(
            "-o", "--output", default=None, help="Output file for metrics (optional)"
        )

        args = parser.parse_args(argv)

        # Load GRNs
        logger.info(f"Loading predicted GRN from: {args.pred_grn}")
        pred_grn = load_grn(args.pred_grn)

        logger.info(f"Loading true GRN from: {args.true_grn}")
        true_grn = load_grn(args.true_grn)

        # Compute metrics
        logger.info("Computing GRN metrics...")
        metrics, metrics_df = evaluate_grn_similarity(
            pred_grn, true_grn, args.quantile
        )

        # Print metrics
        print("\n" + "=" * 60)
        print("GRN Evaluation Metrics")
        print("=" * 60)
        print(f"Precision: {metrics['precision']:.4f}")
        print(f"Recall: {metrics['recall']:.4f}")
        print(f"F1 Score: {metrics['f1_score']:.4f}")
        print(f"AUPRC: {metrics['auprc']:.4f}")
        print(f"Jaccard Index: {metrics['jaccard_index']:.4f}")
        print(f"Hamming Distance: {metrics['hamming_distance']:.4f}")
        print(f"Spearman Correlation: {metrics['spearman_correlation']:.4f}")
        print(
            f"Edge Overlap Ratio: {metrics['edge_overlap']['edge_overlap_ratio']:.4f}"
        )
        print(f"True Positives: {metrics['true_positives']}")
        print(f"False Positives: {metrics['false_positives']}")
        print(f"False Negatives: {metrics['false_negatives']}")
        print(f"True Negatives: {metrics['true_negatives']}")
        print("=" * 60)

        # Save to file if specified
        if args.output:
            metrics_df.to_csv(args.output, index=False)
            logger.info(f"Metrics saved to: {args.output}")

    except GRNMetricsError as e:
        logger.error(f"GRNMetricsError: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("\nInterrupted by user. Cleaning up...")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error during execution: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
