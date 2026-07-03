#!/usr/bin/env python3
"""
ranking_table.py - Visual ranking table for benchmark methods.

Reads grn_metrics.csv and/or cell_metrics.csv produced by compute_metrics.py and
generates a publication-style ranking table PNG where:
  - Each row represents one method.
  - The first column shows the overall rank as a filled horizontal bar
    (darker / wider = better overall rank).
  - Subsequent columns are grouped into "GRN Metrics", "Cell Metrics", "Gene Metrics", and a
    dedicated "Runtime" section. Each column shows a circle whose size and colour
    are proportional to the method's normalised score relative to the best performer.
    The top-N ranked methods also receive a numeric rank label (1, 2, 3 …) inside
    the circle.
  - Overall ranking is computed from category ranks:
    final rank = GRN category rank + Cell category rank + Gene category rank.
  - Methods that produced no output for a metric group receive a hollow grey circle.
  - A method capability table is also exported as method_capabilities.csv.

Usage:
    python ranking_table.py <run_dir> [options]

Arguments:
    run_dir     Directory that contains grn_metrics.csv and/or cell_metrics.csv
                (e.g. benchmark/outputs_metrics/<global_run_name>_<run_type>/)

Options:
    --grn_csv PATH    Path to GRN metrics CSV  (overrides auto-discovery)
    --cell_csv PATH   Path to cell metrics CSV (overrides auto-discovery)
    --output PATH     Output PNG base path for "all methods" table
                      (default: <run_dir>/ranking_table.png). The script also
                      writes "<stem>_native_methods<suffix>" and
                      "<stem>_perturbation_training_methods<suffix>".
    --top_n N         Number of top performers to label with a rank number (default: 3)
    --title TEXT      Optional figure title
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ──────────────────────────────────────────────────────────────────────────────
# Metric catalogue
# (csv_column_name, display_label, higher_is_better)
# ──────────────────────────────────────────────────────────────────────────────

GRN_METRICS: List[Tuple[str, str, bool]] = [
    ("precision",                        "Precision",      True),
    ("recall",                           "Recall",         True),
    ("f1_score",                         "F1 Score",       True),
    ("auprc",                            "AUPRC",          True),
    ("jaccard_index",                    "Jaccard",        True),
    ("spearman_correlation",             "Spearman",       True),
    ("hamming_distance",                 "Hamming",      False),
    ("edge_overlap_edge_overlap_ratio",  "Edge Overlap",   True),
]

CELL_METRICS: List[Tuple[str, str, bool]] = [
    ("wasserstein_distance",      "Wasserstein",  False),
    ("pca_wasserstein_distance",  "PCA Wasserstein",  False),
    ("kl_divergence",             "KL Div.",      False),
    ("jensen_shannon_divergence", "JS Div.",      False),
    ("adjusted_rand_index",       "ARI",             True),
]

GENE_METRICS: List[Tuple[str, str, bool]] = [
    ("gene_gene_correlation", "Gene Corr.", True),
    ("top_de_genes_overlap", "Top DE Overlap", True),
]

RUNTIME_METRIC: List[Tuple[str, str, bool]] = [
    ("runtime", "Runtime", False),
]
PERTURBATION_TRAINING_SUFFIX = "_perturbation_training"

# Import capabilities from the central registry (methods_registry.yaml).
# Falls back to an empty dict if the registry module is unavailable.
try:
    from utils.methods_registry import get_capabilities_table as _get_capabilities_table
    _METHOD_CAPABILITIES: Dict[str, Dict[str, bool]] = _get_capabilities_table()
except Exception:
    _METHOD_CAPABILITIES: Dict[str, Dict[str, bool]] = {}

# ──────────────────────────────────────────────────────────────────────────────
# Layout / style constants
# ──────────────────────────────────────────────────────────────────────────────

# Colour map (Blues: 0 = very light, 1 = dark blue)
_CMAP = plt.cm.Blues

_NO_DATA_FC = "#f0f0f0"   # face colour for "no data" placeholder
_NO_DATA_EC = "#c0c0c0"   # edge colour for "no data" placeholder

_MAX_R = 0.38   # maximum circle radius in data units (cell = 1×1)
_MIN_R = 0.07   # minimum circle radius (score ≈ 0)

_IDENTICAL_VALUE_SCORE = 0.5   # score assigned when all valid values are identical
_TEXT_COLOR_THRESHOLD = 0.40   # normalised score above which rank label uses white text
_MIN_BAR_FILL_RATIO = 0.04     # minimum fill fraction for the overall rank bar

# Inches per data unit (used to compute figure size)
_SCALE = 0.55

# Column / row layout (all in data units)
_NAME_W = 2.4    # width reserved for method-name labels
_HEADER_H = 2.6  # height for rotated column headers + group labels
_LEGEND_H = 2.0  # height for the bottom legend

# Gaps between column groups (data units)
_GAP_OVERALL_GRN = 0.35
_GAP_GRN_CELL = 0.50
_GAP_CELL_GENE = 0.50
_GAP_GENE_RUNTIME = 0.50


# ──────────────────────────────────────────────────────────────────────────────
# Data helpers
# ──────────────────────────────────────────────────────────────────────────────

def _load_csv(path: Optional[Path]) -> Optional[pd.DataFrame]:
    """Load a metrics CSV and average each numeric column per method."""
    if path is None or not path.is_file():
        return None
    df = pd.read_csv(path)
    numeric = [
        c for c in df.select_dtypes(include=np.number).columns
        if c != "replicate"
    ]
    return df.groupby("method")[numeric].mean().reset_index()


def _extract_scores(
    df: Optional[pd.DataFrame],
    metric_defs: List[Tuple[str, str, bool]],
    methods: List[str],
) -> np.ndarray:
    """Return (n_methods, n_metrics) float array; NaN where data is absent."""
    n, m = len(methods), len(metric_defs)
    out = np.full((n, m), np.nan)
    if df is None:
        return out
    idx_map = {meth: i for i, meth in enumerate(methods)}
    for rec in df.to_dict(orient="records"):
        meth = rec.get("method")
        if meth not in idx_map:
            continue
        i = idx_map[meth]
        for j, (col, _, _) in enumerate(metric_defs):
            val = rec.get(col)
            if val is not None and pd.notna(val):
                out[i, j] = float(val)
    return out


def _normalize_col(raw: np.ndarray, higher_is_better: bool) -> np.ndarray:
    """
    Normalise a 1-D metric array to [0, 1] where 1 = best performance.
    NaN values are preserved.
    """
    valid = raw[~np.isnan(raw)]
    if len(valid) == 0:
        return raw.copy()
    mn, mx = valid.min(), valid.max()
    if mx == mn:
        # All values identical → give everyone 0.5
        return np.where(np.isnan(raw), np.nan, _IDENTICAL_VALUE_SCORE)
    if higher_is_better:
        return np.where(np.isnan(raw), np.nan, (raw - mn) / (mx - mn))
    return np.where(np.isnan(raw), np.nan, (mx - raw) / (mx - mn))


def _normalize_matrix(
    raw: np.ndarray,
    metric_defs: List[Tuple[str, str, bool]],
) -> np.ndarray:
    """Normalise each column of (n_methods, n_metrics) independently."""
    norm = np.full_like(raw, np.nan)
    for j, (_, _, hib) in enumerate(metric_defs):
        norm[:, j] = _normalize_col(raw[:, j], hib)
    return norm


def _rank_col(norm_col: np.ndarray) -> np.ndarray:
    """
    Return integer ranks (1 = best) for a normalised column.
    Methods with NaN receive rank 0 (not ranked).
    """
    ranks = np.zeros(len(norm_col), dtype=int)
    valid_idx = np.where(~np.isnan(norm_col))[0]
    if len(valid_idx) == 0:
        return ranks
    order = valid_idx[np.argsort(-norm_col[valid_idx])]
    for rank, idx in enumerate(order, start=1):
        ranks[idx] = rank
    return ranks


def _rank_matrix(norm: np.ndarray) -> np.ndarray:
    """Apply _rank_col to every column of (n_methods, n_cols)."""
    ranks = np.zeros_like(norm, dtype=int)
    for j in range(norm.shape[1]):
        ranks[:, j] = _rank_col(norm[:, j])
    return ranks


def _compute_group_score(norm: np.ndarray) -> np.ndarray:
    """Compute per-method mean score for one metric group."""
    if norm.size == 0:
        return np.full(norm.shape[0], np.nan)
    return np.nanmean(norm, axis=1)


def _category_rank(group_score: np.ndarray, n_methods: int) -> np.ndarray:
    """
    Rank one group and penalize missing values with rank n_methods + 1.
    Lower rank is better.
    """
    ranks = np.full(n_methods, n_methods + 1, dtype=float)
    valid_idx = np.where(~np.isnan(group_score))[0]
    if len(valid_idx) == 0:
        return ranks
    order = valid_idx[np.argsort(-group_score[valid_idx])]
    for rank, idx in enumerate(order, start=1):
        ranks[idx] = float(rank)
    return ranks


def _strip_perturbation_suffix(method: str) -> str:
    if method.endswith(PERTURBATION_TRAINING_SUFFIX):
        return method[: -len(PERTURBATION_TRAINING_SUFFIX)]
    return method


def build_methods_capabilities_table(methods: List[str]) -> pd.DataFrame:
    """Build a method-capability table for the provided methods."""
    rows: List[Dict[str, object]] = []
    for method in sorted(methods):
        base = _strip_perturbation_suffix(method)
        caps = _METHOD_CAPABILITIES.get(
            base,
            {"grn_inference": False, "trajectory_reconstruction": False, "perturbation_training": False},
        )
        rows.append(
            {
                "method": method,
                "base_method": base,
                "grn_inference": bool(caps["grn_inference"]),
                "trajectory_reconstruction": bool(caps["trajectory_reconstruction"]),
                "perturbation_training": bool(caps["perturbation_training"]),
            }
        )
    return pd.DataFrame(rows)


def save_methods_capabilities_plot(capabilities_df: pd.DataFrame, output_path: Path) -> None:
    """Render the method capabilities table as an image."""
    display_df = capabilities_df.copy()
    display_df = display_df.drop(columns=["base_method"], errors="ignore")
    display_df = display_df.rename(
        columns={
            "method": "Method",
            "grn_inference": "GRN inference",
            "trajectory_reconstruction": "Trajectory reconstruction",
            "perturbation_training": "Perturbation training",
        }
    )
    for col in ["GRN inference", "Trajectory reconstruction", "Perturbation training"]:
        if col in display_df.columns:
            display_df[col] = display_df[col].map(lambda v: "✓" if bool(v) else "✗")

    n_rows = max(len(display_df), 1)
    fig_h = max(2.2, 0.45 * n_rows + 1.3)
    fig_w = 8.2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")

    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        loc="center",
        cellLoc="center",
        colLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.2)

    for (row, _col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#dbe7f5")
            cell.set_text_props(weight="bold")
            cell.set_fontsize(7)
        else:
            cell.set_facecolor("#ffffff" if row % 2 else "#f7f9fc")
        cell.set_edgecolor("#d0d7de")
        cell.set_linewidth(0.6)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Method capabilities table plot saved to: {output_path}")


# ──────────────────────────────────────────────────────────────────────────────
# Drawing primitives
# ──────────────────────────────────────────────────────────────────────────────

def _circle_radius(score: float) -> float:
    return _MIN_R + score * (_MAX_R - _MIN_R)


def _score_color(score: float):
    """Map a normalised score in [0, 1] to a Blues colour."""
    return _CMAP(0.18 + 0.78 * score)


def _draw_circle(
    ax: plt.Axes,
    cx: float,
    cy: float,
    score: float,
    rank: int,
    top_n: int,
) -> None:
    """Draw a filled circle; add rank label for top performers."""
    r = _circle_radius(score)
    color = _score_color(score)
    ax.add_patch(
        mpatches.Circle(
            (cx, cy), r,
            facecolor=color,
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
    )
    if 1 <= rank <= top_n:
        text_color = "white" if score > _TEXT_COLOR_THRESHOLD else "#333333"
        ax.text(
            cx, cy, str(rank),
            ha="center", va="center",
            fontsize=6.0, fontweight="bold",
            color=text_color, zorder=4,
        )


def _draw_no_data(ax: plt.Axes, cx: float, cy: float) -> None:
    """Draw an empty (hollow) placeholder circle for missing data."""
    ax.add_patch(
        mpatches.Circle(
            (cx, cy), _MAX_R * 0.48,
            facecolor=_NO_DATA_FC,
            edgecolor=_NO_DATA_EC,
            linewidth=0.6,
            zorder=3,
        )
    )


def _draw_rank_bar(
    ax: plt.Axes,
    cx: float,
    cy: float,
    score: float,
    bar_w: float = 0.82,
) -> None:
    """
    Draw a horizontal rank bar centred at (cx, cy).
    score in [0, 1]; 1 = best overall rank → full dark bar.
    """
    bh = 0.52
    max_w = bar_w
    fill_w = max(max_w * score, max_w * _MIN_BAR_FILL_RATIO)  # ensure at least a thin sliver
    x0 = cx - max_w / 2

    # Light grey background
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x0, cy - bh / 2), max_w, bh,
            boxstyle="round,pad=0",
            facecolor="#e0e0e0", edgecolor="none", zorder=2,
        )
    )
    # Filled foreground
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x0, cy - bh / 2), fill_w, bh,
            boxstyle="round,pad=0",
            facecolor=_score_color(score), edgecolor="none", zorder=3,
        )
    )


# ──────────────────────────────────────────────────────────────────────────────
# Column-position helper
# ──────────────────────────────────────────────────────────────────────────────

class _Layout:
    """Pre-computes x-centres for every column."""

    def __init__(self, n_grn: int, n_cell: int, n_gene: int, n_runtime: int) -> None:
        self.n_grn = n_grn
        self.n_cell = n_cell
        self.n_gene = n_gene
        self.n_runtime = n_runtime

        self.x_overall = 0.5

        grn_start = self.x_overall + 1.0 + _GAP_OVERALL_GRN
        self._grn_centers = [grn_start + i for i in range(n_grn)]

        if n_grn > 0:
            cell_start = grn_start + n_grn + _GAP_GRN_CELL
        else:
            cell_start = grn_start
        self._cell_centers = [cell_start + i for i in range(n_cell)]

        if n_cell > 0:
            gene_start = cell_start + n_cell + _GAP_CELL_GENE
        elif n_grn > 0:
            gene_start = grn_start + n_grn + _GAP_CELL_GENE
        else:
            gene_start = grn_start
        self._gene_centers = [gene_start + i for i in range(n_gene)]

        if n_gene > 0:
            runtime_start = gene_start + n_gene + _GAP_GENE_RUNTIME
        elif n_cell > 0:
            runtime_start = cell_start + n_cell + _GAP_GENE_RUNTIME
        elif n_grn > 0:
            runtime_start = grn_start + n_grn + _GAP_GENE_RUNTIME
        else:
            runtime_start = grn_start
        self._runtime_centers = [runtime_start + i for i in range(n_runtime)]

        # Total data width (from x=0 to right edge)
        if n_runtime > 0:
            self.data_w = self._runtime_centers[-1] + 0.5
        elif n_gene > 0:
            self.data_w = self._gene_centers[-1] + 0.5
        elif n_cell > 0:
            self.data_w = self._cell_centers[-1] + 0.5
        elif n_grn > 0:
            self.data_w = self._grn_centers[-1] + 0.5
        else:
            self.data_w = self.x_overall + 0.5

    def x_grn(self, i: int) -> float:
        return self._grn_centers[i]

    def x_cell(self, i: int) -> float:
        return self._cell_centers[i]

    def x_gene(self, i: int) -> float:
        return self._gene_centers[i]

    def x_runtime(self, i: int) -> float:
        return self._runtime_centers[i]

    def grn_group_span(self) -> Tuple[float, float]:
        if not self._grn_centers:
            return (0.0, 0.0)
        return self._grn_centers[0] - 0.45, self._grn_centers[-1] + 0.45

    def cell_group_span(self) -> Tuple[float, float]:
        if not self._cell_centers:
            return (0.0, 0.0)
        return self._cell_centers[0] - 0.45, self._cell_centers[-1] + 0.45

    def runtime_group_span(self) -> Tuple[float, float]:
        if not self._runtime_centers:
            return (0.0, 0.0)
        return self._runtime_centers[0] - 0.45, self._runtime_centers[-1] + 0.45

    def gene_group_span(self) -> Tuple[float, float]:
        if not self._gene_centers:
            return (0.0, 0.0)
        return self._gene_centers[0] - 0.45, self._gene_centers[-1] + 0.45


# ──────────────────────────────────────────────────────────────────────────────
# Main rendering function
# ──────────────────────────────────────────────────────────────────────────────

def render_ranking_table(
    methods: List[str],
    overall_scores: np.ndarray,
    grn_norm: np.ndarray,
    cell_norm: np.ndarray,
    gene_norm: np.ndarray,
    runtime_norm: np.ndarray,
    grn_ranks: np.ndarray,
    cell_ranks: np.ndarray,
    gene_ranks: np.ndarray,
    runtime_ranks: np.ndarray,
    grn_metric_names: List[str],
    cell_metric_names: List[str],
    gene_metric_names: List[str],
    runtime_metric_name: Optional[str] = None,
    top_n: int = 3,
    title: str = "",
) -> plt.Figure:
    """
    Build and return the ranking-table figure.

    Parameters
    ----------
    methods           : Method names (sorted best → worst by overall rank).
    overall_scores    : Normalised overall score in [0, 1] per method.
    grn_norm          : (n_methods, n_grn) normalised scores; NaN = no data.
    cell_norm         : (n_methods, n_cell) normalised scores; NaN = no data.
    gene_norm         : (n_methods, n_gene) normalised scores; NaN = no data.
    runtime_norm      : (n_methods, n_runtime) normalised scores; NaN = no data.
    grn_ranks         : (n_methods, n_grn) per-metric integer ranks (1 = best).
    cell_ranks        : (n_methods, n_cell) per-metric integer ranks (1 = best).
    gene_ranks        : (n_methods, n_gene) per-metric integer ranks (1 = best).
    runtime_ranks     : (n_methods, n_runtime) per-metric integer ranks (1 = best).
    grn_metric_names  : Display labels for GRN metric columns.
    cell_metric_names : Display labels for Cell metric columns.
    gene_metric_names : Display labels for Gene metric columns.
    runtime_metric_name : Display label for runtime metric column (if present).
    top_n             : How many top performers receive a numeric label.
    title             : Optional figure super-title.

    Returns
    -------
    matplotlib Figure
    """
    n_methods = len(methods)
    n_grn = len(grn_metric_names)
    n_cell = len(cell_metric_names)
    n_gene = len(gene_metric_names)
    n_runtime = 1 if runtime_metric_name else 0

    layout = _Layout(n_grn, n_cell, n_gene, n_runtime)

    # ── figure size ──────────────────────────────────────────────────────────
    PAD = 0.25
    x_min = -_NAME_W - PAD
    x_max = layout.data_w + PAD
    y_min = -_LEGEND_H - PAD
    y_max = n_methods + _HEADER_H + PAD

    x_span = x_max - x_min
    y_span = y_max - y_min

    fig_w = max(x_span * _SCALE, 5.0)
    fig_h = max(y_span * _SCALE, 3.0)

    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.axis("off")

    def y_row(i: int) -> float:
        """Centre y-coordinate for method i (0 = top/best)."""
        return n_methods - i - 0.5

    # ── title ────────────────────────────────────────────────────────────────
    if title:
        ax.text(
            (layout.data_w - _NAME_W) / 2,
            n_methods + _HEADER_H - 0.2,
            title,
            ha="center", va="top",
            fontsize=10, fontweight="bold", color="#111111",
        )

    # ── alternating row backgrounds ──────────────────────────────────────────
    for i in range(n_methods):
        if i % 2 == 0:
            ax.add_patch(
                mpatches.Rectangle(
                    (-_NAME_W, y_row(i) - 0.5),
                    _NAME_W + layout.data_w,
                    1.0,
                    facecolor="#f4f6f8",
                    edgecolor="none",
                    zorder=1,
                )
            )

    # ── horizontal row dividers ──────────────────────────────────────────────
    for i in range(n_methods + 1):
        y_line = n_methods - i
        ax.plot(
            [-_NAME_W, layout.data_w],
            [y_line, y_line],
            color="#d8d8d8", linewidth=0.4, zorder=0,
        )

    # ── group header boxes ───────────────────────────────────────────────────
    grp_y = n_methods + _HEADER_H - 1.0  # group label centre y

    # "Overall" label (above rank bar)
    ax.text(
        layout.x_overall, grp_y,
        "Overall",
        ha="center", va="center",
        fontsize=8.0, fontweight="bold", color="#333333",
    )

    # GRN group box
    if n_grn > 0:
        gx0, gx1 = layout.grn_group_span()
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (gx0, grp_y - 0.35), gx1 - gx0, 0.70,
                boxstyle="round,pad=0.05",
                facecolor="#dde8f7", edgecolor="#99b8e0",
                linewidth=0.8, zorder=2,
            )
        )
        ax.text(
            (gx0 + gx1) / 2, grp_y,
            "GRN Metrics",
            ha="center", va="center",
            fontsize=8.0, fontweight="bold", color="#1a3a6c",
            zorder=3,
        )

    # Cell group box
    if n_cell > 0:
        cx0, cx1 = layout.cell_group_span()
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (cx0, grp_y - 0.35), cx1 - cx0, 0.70,
                boxstyle="round,pad=0.05",
                facecolor="#ddf5e4", edgecolor="#99d4b0",
                linewidth=0.8, zorder=2,
            )
        )
        ax.text(
            (cx0 + cx1) / 2, grp_y,
            "Cell Metrics",
            ha="center", va="center",
            fontsize=8.0, fontweight="bold", color="#1a5c38",
            zorder=3,
        )

    # Runtime group box
    if n_gene > 0:
        gx0, gx1 = layout.gene_group_span()
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (gx0, grp_y - 0.35), gx1 - gx0, 0.70,
                boxstyle="round,pad=0.05",
                facecolor="#f1ddf5", edgecolor="#c99bd4",
                linewidth=0.8, zorder=2,
            )
        )
        ax.text(
            (gx0 + gx1) / 2, grp_y,
            "Gene Metrics",
            ha="center", va="center",
            fontsize=8.0, fontweight="bold", color="#5d1a6c",
            zorder=3,
        )

    # Runtime group box
    if n_runtime > 0:
        rx0, rx1 = layout.runtime_group_span()
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (rx0, grp_y - 0.35), rx1 - rx0, 0.70,
                boxstyle="round,pad=0.05",
                facecolor="#f5e8dd", edgecolor="#d4b399",
                linewidth=0.8, zorder=2,
            )
        )
        ax.text(
            (rx0 + rx1) / 2, grp_y,
            "Runtime",
            ha="center", va="center",
            fontsize=8.0, fontweight="bold", color="#6c3f1a",
            zorder=3,
        )

    # ── column labels (rotated 30°) ──────────────────────────────────────────
    label_y = n_methods + _HEADER_H - 2.3
    rot = 30

    ax.text(
        layout.x_overall, label_y, "Rank",
        ha="left", va="top",
        fontsize=7.0, rotation=rot, rotation_mode="anchor",
        color="#333333",
    )

    for i, name in enumerate(grn_metric_names):
        ax.text(
            layout.x_grn(i), label_y, name,
            ha="left", va="top",
            fontsize=7.0, rotation=rot, rotation_mode="anchor",
            color="#1a3a6c",
        )

    for i, name in enumerate(cell_metric_names):
        ax.text(
            layout.x_cell(i), label_y, name,
            ha="left", va="top",
            fontsize=7.0, rotation=rot, rotation_mode="anchor",
            color="#1a5c38",
        )

    for i, name in enumerate(gene_metric_names):
        ax.text(
            layout.x_gene(i), label_y, name,
            ha="left", va="top",
            fontsize=7.0, rotation=rot, rotation_mode="anchor",
            color="#5d1a6c",
        )

    if n_runtime > 0 and runtime_metric_name is not None:
        ax.text(
            layout.x_runtime(0), label_y, runtime_metric_name,
            ha="left", va="top",
            fontsize=7.0, rotation=rot, rotation_mode="anchor",
            color="#6c3f1a",
        )

    # ── vertical group dividers ──────────────────────────────────────────────
    div_y0 = 0.0
    div_y1 = n_methods

    # Divider after "Overall" column
    xd1 = layout.x_overall + 0.5 + _GAP_OVERALL_GRN * 0.5
    ax.plot([xd1, xd1], [div_y0, div_y1], color="#b0b0b0", linewidth=0.8, zorder=0)

    # Divider between GRN and Cell
    if n_grn > 0 and n_cell > 0:
        xd2 = layout.x_grn(n_grn - 1) + 0.5 + _GAP_GRN_CELL * 0.5
        ax.plot([xd2, xd2], [div_y0, div_y1], color="#b0b0b0", linewidth=0.8, zorder=0)

    # Divider between Cell and Gene
    if n_gene > 0 and n_cell > 0:
        xd3 = layout.x_cell(n_cell - 1) + 0.5 + _GAP_CELL_GENE * 0.5
        ax.plot([xd3, xd3], [div_y0, div_y1], color="#b0b0b0", linewidth=0.8, zorder=0)
    elif n_gene > 0 and n_grn > 0 and n_cell == 0:
        xd3 = layout.x_grn(n_grn - 1) + 0.5 + _GAP_CELL_GENE * 0.5
        ax.plot([xd3, xd3], [div_y0, div_y1], color="#b0b0b0", linewidth=0.8, zorder=0)

    # Divider between Gene and Runtime
    if n_runtime > 0 and n_gene > 0:
        xd4 = layout.x_gene(n_gene - 1) + 0.5 + _GAP_GENE_RUNTIME * 0.5
        ax.plot([xd4, xd4], [div_y0, div_y1], color="#b0b0b0", linewidth=0.8, zorder=0)
    elif n_runtime > 0 and n_cell > 0 and n_gene == 0:
        xd4 = layout.x_cell(n_cell - 1) + 0.5 + _GAP_GENE_RUNTIME * 0.5
        ax.plot([xd4, xd4], [div_y0, div_y1], color="#b0b0b0", linewidth=0.8, zorder=0)
    elif n_runtime > 0 and n_grn > 0 and n_cell == 0 and n_gene == 0:
        xd4 = layout.x_grn(n_grn - 1) + 0.5 + _GAP_GENE_RUNTIME * 0.5
        ax.plot([xd4, xd4], [div_y0, div_y1], color="#b0b0b0", linewidth=0.8, zorder=0)

    # ── per-method rows ──────────────────────────────────────────────────────
    overall_ranks_arr = _rank_col(overall_scores)

    for i, method in enumerate(methods):
        y = y_row(i)

        # Method name
        ax.text(
            -0.12, y, method,
            ha="right", va="center",
            fontsize=8.0, color="#222222",
        )

        # Overall rank bar
        _draw_rank_bar(ax, layout.x_overall, y, float(overall_scores[i]))

        # Overall rank number (small, to the right of the bar)
        ovr = int(overall_ranks_arr[i])
        if ovr > 0:
            ax.text(
                layout.x_overall + 0.44, y,
                f"#{ovr}",
                ha="left", va="center",
                fontsize=5.5, color="#555555",
            )

        # GRN metric circles
        for j in range(n_grn):
            score = grn_norm[i, j]
            rank = int(grn_ranks[i, j])
            if np.isnan(score):
                _draw_no_data(ax, layout.x_grn(j), y)
            else:
                _draw_circle(ax, layout.x_grn(j), y, float(score), rank, top_n)

        # Cell metric circles
        for j in range(n_cell):
            score = cell_norm[i, j]
            rank = int(cell_ranks[i, j])
            if np.isnan(score):
                _draw_no_data(ax, layout.x_cell(j), y)
            else:
                _draw_circle(ax, layout.x_cell(j), y, float(score), rank, top_n)

        # Gene metric circles
        for j in range(n_gene):
            score = gene_norm[i, j]
            rank = int(gene_ranks[i, j])
            if np.isnan(score):
                _draw_no_data(ax, layout.x_gene(j), y)
            else:
                _draw_circle(ax, layout.x_gene(j), y, float(score), rank, top_n)

        # Runtime metric circles
        for j in range(n_runtime):
            score = runtime_norm[i, j]
            rank = int(runtime_ranks[i, j])
            if np.isnan(score):
                _draw_no_data(ax, layout.x_runtime(j), y)
            else:
                _draw_circle(ax, layout.x_runtime(j), y, float(score), rank, top_n)

    # ── legend ───────────────────────────────────────────────────────────────
    leg_y_top = -0.35
    leg_y_graphics = leg_y_top - 0.55

    # Determine x starting positions for the 4 legend blocks
    x1 = -_NAME_W + 0.2            # Overall rank bar
    x2 = x1 + 5                  # Rank (circles)
    x3 = x2 + 5                  # Score color bar
    x4 = x3 + 5                  # No data indicator

    # 1. Overall rank bar legend
    ax.text(
        x1, leg_y_top,
        "Overall rank bar",
        ha="left", va="center",
        fontsize=7.5, fontweight="bold", color="#333333",
    )
    bar_legend_items = [("Bottom", 0.05), ("", 0.33), ("", 0.67), ("Top", 1.0)]
    for k, (lbl, sc) in enumerate(bar_legend_items):
        bx = x1 + 0.4 + k * 0.85
        _draw_rank_bar(ax, bx, leg_y_graphics, sc, bar_w=0.85)
        if lbl:
            ax.text(
                bx, leg_y_graphics - 0.45, lbl,
                ha="center", va="top",
                fontsize=6.0, color="#444444",
            )

    # 2. Circle-size legend (Rank)
    ax.text(
        x2, leg_y_top,
        "Rank (circles)",
        ha="left", va="center",
        fontsize=7.5, fontweight="bold", color="#333333",
    )
    circle_legend_items = [("Bottom", 0.0), ("", 0.33), ("", 0.67), ("Top", 1.0)]
    for k, (lbl, sc) in enumerate(circle_legend_items):
        lx = x2 + 0.3 + k * 0.75
        ax.add_patch(
            mpatches.Circle(
                (lx, leg_y_graphics), _circle_radius(sc),
                facecolor="#8dafd4", edgecolor="white",
                linewidth=0.5, zorder=3,
            )
        )
        if lbl:
            ax.text(
                lx, leg_y_graphics - 0.50, lbl,
                ha="center", va="top",
                fontsize=6.0, color="#444444",
            )

    # 3. Score colour-bar
    cb_w = 2.5
    cb_h = 0.28
    cb_x0 = x3
    cb_y0 = leg_y_graphics - cb_h / 2

    ax.text(
        cb_x0 + cb_w / 2, leg_y_top,
        "Score",
        ha="center", va="center",
        fontsize=7.5, fontweight="bold", color="#333333",
    )

    n_steps = 60
    for k in range(n_steps):
        frac = k / n_steps
        ax.add_patch(
            mpatches.Rectangle(
                (cb_x0 + frac * cb_w, cb_y0),
                cb_w / n_steps, cb_h,
                facecolor=_score_color(frac),
                edgecolor="none",
            )
        )
    ax.add_patch(
        mpatches.Rectangle(
            (cb_x0, cb_y0), cb_w, cb_h,
            facecolor="none", edgecolor="#888888",
            linewidth=0.6,
        )
    )
    ax.text(cb_x0, cb_y0 - 0.22, "0", ha="center", va="top", fontsize=6.0, color="#555555")
    ax.text(cb_x0 + cb_w, cb_y0 - 0.22, "1", ha="center", va="top", fontsize=6.0, color="#555555")

    # 4. No-data indicator
    _draw_no_data(ax, x4, leg_y_graphics)
    ax.text(
        x4 + 0.50, leg_y_graphics, "= No data",
        ha="left", va="center",
        fontsize=7.0, color="#555555",
    )

    return fig


def _build_and_save_ranking_table(
    grn_df: Optional[pd.DataFrame],
    cell_df: Optional[pd.DataFrame],
    methods: List[str],
    output_path: Path,
    top_n: int,
    title: str,
) -> bool:
    """Build and save one ranking table for a selected set of methods."""
    n_methods = len(methods)
    if n_methods == 0:
        return False

    grn_raw = _extract_scores(grn_df, GRN_METRICS, methods)
    cell_raw = _extract_scores(cell_df, CELL_METRICS, methods)
    gene_raw = _extract_scores(cell_df, GENE_METRICS, methods)
    runtime_raw = _extract_scores(cell_df, RUNTIME_METRIC, methods)

    active_grn_idx = [j for j in range(len(GRN_METRICS)) if not np.all(np.isnan(grn_raw[:, j]))]
    active_cell_idx = [j for j in range(len(CELL_METRICS)) if not np.all(np.isnan(cell_raw[:, j]))]
    active_gene_idx = [j for j in range(len(GENE_METRICS)) if not np.all(np.isnan(gene_raw[:, j]))]
    active_runtime_idx = [j for j in range(len(RUNTIME_METRIC)) if not np.all(np.isnan(runtime_raw[:, j]))]

    grn_raw_a = grn_raw[:, active_grn_idx] if active_grn_idx else np.empty((n_methods, 0))
    cell_raw_a = cell_raw[:, active_cell_idx] if active_cell_idx else np.empty((n_methods, 0))
    gene_raw_a = gene_raw[:, active_gene_idx] if active_gene_idx else np.empty((n_methods, 0))
    runtime_raw_a = runtime_raw[:, active_runtime_idx] if active_runtime_idx else np.empty((n_methods, 0))

    active_grn_defs = [GRN_METRICS[j] for j in active_grn_idx]
    active_cell_defs = [CELL_METRICS[j] for j in active_cell_idx]
    active_gene_defs = [GENE_METRICS[j] for j in active_gene_idx]
    active_runtime_defs = [RUNTIME_METRIC[j] for j in active_runtime_idx]

    grn_norm = (
        _normalize_matrix(grn_raw_a, active_grn_defs)
        if grn_raw_a.size > 0
        else grn_raw_a.copy()
    )
    cell_norm = (
        _normalize_matrix(cell_raw_a, active_cell_defs)
        if cell_raw_a.size > 0
        else cell_raw_a.copy()
    )
    gene_norm = (
        _normalize_matrix(gene_raw_a, active_gene_defs)
        if gene_raw_a.size > 0
        else gene_raw_a.copy()
    )
    runtime_norm = (
        _normalize_matrix(runtime_raw_a, active_runtime_defs)
        if runtime_raw_a.size > 0
        else runtime_raw_a.copy()
    )

    grn_group = _compute_group_score(grn_norm)
    cell_group = _compute_group_score(cell_norm)
    gene_group = _compute_group_score(gene_norm)
    grn_rank = _category_rank(grn_group, n_methods)
    cell_rank = _category_rank(cell_group, n_methods)
    gene_rank = _category_rank(gene_group, n_methods)

    overall_rank_sum = grn_rank + cell_rank + gene_rank
    sort_idx = np.argsort(overall_rank_sum)
    overall_order_rank = np.empty(n_methods, dtype=float)
    overall_order_rank[sort_idx] = np.arange(1, n_methods + 1, dtype=float)
    if n_methods == 1:
        overall_scores = np.array([1.0], dtype=float)
    else:
        overall_scores = 1.0 - (overall_order_rank - 1.0) / float(n_methods - 1)

    methods_s = [methods[i] for i in sort_idx]
    overall_s = overall_scores[sort_idx]
    grn_norm_s = grn_norm[sort_idx] if grn_norm.size > 0 else grn_norm
    cell_norm_s = cell_norm[sort_idx] if cell_norm.size > 0 else cell_norm
    gene_norm_s = gene_norm[sort_idx] if gene_norm.size > 0 else gene_norm
    runtime_norm_s = runtime_norm[sort_idx] if runtime_norm.size > 0 else runtime_norm

    grn_ranks_s = _rank_matrix(grn_norm_s) if grn_norm_s.size > 0 else np.zeros_like(grn_norm_s, dtype=int)
    cell_ranks_s = _rank_matrix(cell_norm_s) if cell_norm_s.size > 0 else np.zeros_like(cell_norm_s, dtype=int)
    gene_ranks_s = _rank_matrix(gene_norm_s) if gene_norm_s.size > 0 else np.zeros_like(gene_norm_s, dtype=int)
    runtime_ranks_s = _rank_matrix(runtime_norm_s) if runtime_norm_s.size > 0 else np.zeros_like(runtime_norm_s, dtype=int)

    grn_names = [name for _, name, _ in active_grn_defs]
    cell_names = [name for _, name, _ in active_cell_defs]
    gene_names = [name for _, name, _ in active_gene_defs]
    runtime_name = active_runtime_defs[0][1] if active_runtime_defs else None

    fig = render_ranking_table(
        methods=methods_s,
        overall_scores=overall_s,
        grn_norm=grn_norm_s,
        cell_norm=cell_norm_s,
        gene_norm=gene_norm_s,
        runtime_norm=runtime_norm_s,
        grn_ranks=grn_ranks_s,
        cell_ranks=cell_ranks_s,
        gene_ranks=gene_ranks_s,
        runtime_ranks=runtime_ranks_s,
        grn_metric_names=grn_names,
        cell_metric_names=cell_names,
        gene_metric_names=gene_names,
        runtime_metric_name=runtime_name,
        top_n=top_n,
        title=title,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Ranking table saved to: {output_path}")
    return True


# ──────────────────────────────────────────────────────────────────────────────
# CLI entry-point
# ──────────────────────────────────────────────────────────────────────────────

def main(argv: Optional[List[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Generate a visual ranking table from benchmark metrics CSVs."
    )
    parser.add_argument(
        "run_dir",
        help=(
            "Directory containing grn_metrics.csv and/or cell_metrics.csv "
            "(e.g. benchmark/outputs_metrics/<global_run_name>_<run_type>/)"
        ),
    )
    parser.add_argument(
        "--grn_csv", default=None,
        help="Explicit path to GRN metrics CSV (overrides auto-discovery).",
    )
    parser.add_argument(
        "--cell_csv", default=None,
        help="Explicit path to cell metrics CSV (overrides auto-discovery).",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output PNG path (default: <run_dir>/ranking_table.png).",
    )
    parser.add_argument(
        "--top_n", type=int, default=3,
        help="Number of top performers to label with a rank number (default: 3).",
    )
    parser.add_argument(
        "--title", default="",
        help="Optional figure title.",
    )

    args = parser.parse_args(argv)

    run_dir = Path(args.run_dir)

    grn_csv_path = Path(args.grn_csv) if args.grn_csv else run_dir / "grn_metrics.csv"
    cell_csv_path = Path(args.cell_csv) if args.cell_csv else run_dir / "cell_metrics.csv"
    output_path = Path(args.output) if args.output else run_dir / "ranking_table.png"
    output_native_path = output_path.with_name(f"{output_path.stem}_native_methods{output_path.suffix}")
    output_perturb_path = output_path.with_name(
        f"{output_path.stem}_perturbation_training_methods{output_path.suffix}"
    )

    # ── load data ────────────────────────────────────────────────────────────
    grn_df = _load_csv(grn_csv_path)
    cell_df = _load_csv(cell_csv_path)

    if grn_df is None and cell_df is None:
        print(
            "ERROR: No metrics CSVs found. "
            "Provide at least grn_metrics.csv or cell_metrics.csv in the run_dir.",
            file=sys.stderr,
        )
        sys.exit(1)

    if grn_df is None:
        print(f"INFO: No GRN metrics CSV found at '{grn_csv_path}'. GRN columns will be empty.")
    if cell_df is None:
        print(f"INFO: No cell metrics CSV found at '{cell_csv_path}'. Cell columns will be empty.")

    # ── collect all methods ──────────────────────────────────────────────────
    methods_set: set[str] = set()
    if grn_df is not None:
        methods_set.update(grn_df["method"].tolist())
    if cell_df is not None:
        methods_set.update(cell_df["method"].tolist())
    methods = sorted(methods_set)
    n_methods = len(methods)

    if n_methods == 0:
        print("ERROR: No methods found in the CSVs.", file=sys.stderr)
        sys.exit(1)

    native_methods = sorted([m for m in methods if not m.endswith(PERTURBATION_TRAINING_SUFFIX)])
    perturb_methods = sorted([m for m in methods if m.endswith(PERTURBATION_TRAINING_SUFFIX)])

    capabilities_df = build_methods_capabilities_table(native_methods)
    capabilities_path = run_dir / "method_capabilities.csv"
    capabilities_df.to_csv(capabilities_path, index=False)
    print(f"Method capabilities table saved to: {capabilities_path}")
    capabilities_plot_path = run_dir / "method_capabilities.png"
    save_methods_capabilities_plot(capabilities_df, capabilities_plot_path)

    if args.title:
        title_both = f"{args.title} (All methods)"
        title_native = f"{args.title} (Native methods)"
        title_perturb = f"{args.title} (Perturbation training methods)"
    else:
        title_both = f"Method Ranking — {run_dir.name} (All methods)"
        title_native = f"Method Ranking — {run_dir.name} (Native methods)"
        title_perturb = f"Method Ranking — {run_dir.name} (Perturbation training methods)"

    _build_and_save_ranking_table(
        grn_df=grn_df,
        cell_df=cell_df,
        methods=methods,
        output_path=output_path,
        top_n=args.top_n,
        title=title_both,
    )

    if native_methods:
        _build_and_save_ranking_table(
            grn_df=grn_df,
            cell_df=cell_df,
            methods=native_methods,
            output_path=output_native_path,
            top_n=args.top_n,
            title=title_native,
        )
    else:
        print("INFO: No native methods found; skipping native-only ranking table.")

    if perturb_methods:
        _build_and_save_ranking_table(
            grn_df=grn_df,
            cell_df=cell_df,
            methods=perturb_methods,
            output_path=output_perturb_path,
            top_n=args.top_n,
            title=title_perturb,
        )
    else:
        print(
            "INFO: No perturbation-training methods found; "
            "skipping perturbation-only ranking table."
        )


if __name__ == "__main__":
    main()
