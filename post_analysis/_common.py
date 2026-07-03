#!/usr/bin/env python3
"""
_common.py — Shared constants, regex patterns, and utilities for the
post_analysis package.

Imported by metrics_cells, metrics_grn, visualizations, and compute_metrics.
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PERTURBATION_TRAINING_SUFFIX: str = "_perturbation_training"
GRN_TRUE_FOLDER: str = "GRN_true"

# ---------------------------------------------------------------------------
# Compiled regex patterns
# ---------------------------------------------------------------------------

ADATA_FILE_RE: re.Pattern = re.compile(
    r"^data_(\d+)(?:_ko_([A-Za-z0-9_-]+))?(_perturbation_training)?_adata\.h5ad$"
)
TRUE_ADATA_RE: re.Pattern = re.compile(
    r"^data_(\d+)(?:_ko_([A-Za-z0-9_-]+))?(_perturbation_training)?\.h5ad$"
)
GRN_FILE_RE: re.Pattern = re.compile(
    r"^data_(\d+)(?:_ko_([A-Za-z0-9_-]+))?(_perturbation_training)?_GRN\.npy$"
)
RUN_FOLD_RE: re.Pattern = re.compile(r"^run_(\d+)$")

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------


def setup_logging(name: Optional[str] = None) -> logging.Logger:
    """Configure and return a logger with a standard format."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(name or __name__)


# ---------------------------------------------------------------------------
# File discovery utilities
# ---------------------------------------------------------------------------


def discover_replicates(folder: Path, file_re: re.Pattern) -> Dict[tuple, Path]:
    """Return ``{(replicate, perturbation, perturbation_training): path}``
    for every file in *folder* whose name matches *file_re*."""
    result: Dict[tuple, Path] = {}
    if not folder.is_dir():
        return result
    for f in sorted(folder.iterdir()):
        m = file_re.match(f.name)
        if m:
            result[(int(m.group(1)), m.group(2), bool(m.group(3)))] = f
    return result


def is_leave_one_out_structure(run_dir: Path) -> bool:
    """Return ``True`` if *run_dir* contains ``run_N`` fold subdirectories
    instead of flat method directories."""
    if not run_dir.is_dir():
        return False
    children = list(run_dir.iterdir())
    # Must have at least one run_N dir and no flat method dirs
    # (aside from GRN_true).
    has_run_dirs = any(
        d.is_dir() and RUN_FOLD_RE.match(d.name) for d in children
    )
    if not has_run_dirs:
        return False
    # Confirm there are no "method-like" dirs (i.e. dirs that are not
    # run_N and not GRN_true) at the top level.
    non_fold_dirs = [
        d for d in children
        if d.is_dir() and not RUN_FOLD_RE.match(d.name) and d.name != GRN_TRUE_FOLDER
    ]
    return len(non_fold_dirs) == 0


def iter_method_dirs(
    run_dir: Path,
) -> Iterator[Tuple[str, Path, Optional[int]]]:
    """Yield ``(method_name, method_path, fold)`` for every method directory
    in *run_dir*, handling both flat and leave-one-out nested structures.

    In a **flat** structure::

        run_dir/
          MethodA/
          MethodB/
          GRN_true/

    the fold is ``None``.  In a **leave-one-out** structure::

        run_dir/
          run_1/
            MethodA/
            MethodB/
          run_2/
            MethodA/
          GRN_true/

    the fold is the integer extracted from ``run_N``.
    """
    if not run_dir.is_dir():
        return

    if is_leave_one_out_structure(run_dir):
        for fold_dir in sorted(run_dir.iterdir()):
            m = RUN_FOLD_RE.match(fold_dir.name)
            if not m or not fold_dir.is_dir():
                continue
            fold = int(m.group(1))
            for method_dir in sorted(fold_dir.iterdir()):
                if method_dir.is_dir() and method_dir.name != GRN_TRUE_FOLDER:
                    yield (method_dir.name, method_dir, fold)
    else:
        for method_dir in sorted(run_dir.iterdir()):
            if method_dir.is_dir() and method_dir.name != GRN_TRUE_FOLDER:
                yield (method_dir.name, method_dir, None)


def parse_ko_output_genes(value: Optional[str]) -> List[str]:
    """Parse a KO-output-genes specification string.

    Accepted values:
      - ``None``, ``""``, ``"none"``  → empty list (WT only)
      - ``"all"``                      → ``["all"]`` (all KO datasets)
      - comma-separated gene names    → list of gene names
    """
    if value is None:
        return []
    normalized = str(value).strip()
    if normalized.lower() in {"", "none"}:
        return []
    if normalized.lower() == "all":
        return ["all"]
    return [gene.strip() for gene in normalized.split(",") if gene.strip()]


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def flatten_metrics(metrics: Dict[str, object]) -> Dict[str, object]:
    """Flatten nested dicts one level deep (e.g. edge_overlap sub-dict)."""
    flat: Dict[str, object] = {}
    for k, v in metrics.items():
        if isinstance(v, dict):
            for sub_k, sub_v in v.items():
                flat[f"{k}_{sub_k}"] = sub_v
        else:
            flat[k] = v
    return flat


def resolve_benchmark_path(
    script_dir: Path,
    cli_arg: Optional[str],
    default_subpath: str,
) -> Path:
    """Resolve a path from a CLI argument, falling back to a default relative
    to the benchmark project root.

    Args:
        script_dir: ``Path(__file__).resolve().parent``
        cli_arg: Value from argparse (may be ``None``).
        default_subpath: Path suffix relative to the project root, e.g.
                         ``"benchmark/outputs_methods"``.

    Returns:
        Resolved ``Path``.
    """
    if cli_arg:
        return Path(cli_arg)
    return script_dir.parent / default_subpath
