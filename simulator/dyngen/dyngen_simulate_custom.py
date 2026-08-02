#!/usr/bin/env python3
"""
dyngen_simulate_custom.py — Unified Python wrapper for dyngen simulation.

Calls dyngen (R package) via Rscript subprocess.  R and dyngen must be
available on PATH (provided by the benchmark_runner conda environment).

The R script uses dyngen's SSA cell simulation with experiment_synchronised()
to produce time-stratified expression data.  Ground-truth GRN is written
directly from the benchmark config in common_config.py.

Supports both regular (WT-only) and KO simulation modes via --ko_genes.

Usage (WT only):
    python dyngen_simulate_custom.py -o <output_folder> [-a <adata_folder>]
        [-n <n_runs>] [--n_genes 8] [--n_cells 1000] [--n_time_bins 20]
        [--seed 42]

Usage (KO mode):
    python dyngen_simulate_custom.py -o <output_folder> -k <all|Gene0,Gene1,...>
        [-n <n_runs>] [--n_genes 8] [--n_cells 1000] [--n_time_bins 20]
        [--seed 42]
"""

from __future__ import annotations

import getopt
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from simulator.common_config import (
    SimulatorCommonConfig,
    build_dyngen_backbone_spec,
    get_benchmark_grn,
)

# ── Constants ────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
R_SCRIPT = SCRIPT_DIR / "dyngen_simulate.R"


class DyngenRunConfig:
    """dyngen-specific parameters layered on top of SimulatorCommonConfig."""

    def __init__(self, common: SimulatorCommonConfig, seed: int = 42):
        self.common = common
        self.seed = seed

    @property
    def output_folder(self) -> Path:
        return self.common.output_folder

    @property
    def n_runs(self) -> int:
        return self.common.n_runs

    @property
    def n_genes(self) -> int:
        return self.common.n_genes

    @property
    def n_cells(self) -> int:
        return self.common.n_cells

    @property
    def n_time_bins(self) -> int:
        return self.common.n_time_bins

    @property
    def ko_genes(self) -> str:
        return self.common.ko_genes


def parse_cli(argv: Sequence[str]) -> DyngenRunConfig:
    """Parse command-line arguments, combining common + dyngen-specific params."""
    output_folder: Optional[str] = None
    adata_folder: Optional[str] = None
    n_runs = 2
    ko_genes = "none"
    n_genes = 8
    n_cells = 1000
    n_time_bins = 20
    seed = 42

    opts, _ = getopt.getopt(
        list(argv),
        "ho:a:n:k:",
        [
            "output_folder=", "adata_folder=", "n_runs=", "ko_genes=",
            "n_genes=", "n_cells=", "n_time_bins=", "seed=",
        ],
    )

    for opt, arg in opts:
        if opt in ("-o", "--output_folder"):
            output_folder = arg
        elif opt in ("-a", "--adata_folder"):
            adata_folder = arg
        elif opt in ("-n", "--n_runs"):
            n_runs = int(arg)
        elif opt in ("-k", "--ko_genes"):
            ko_genes = arg
        elif opt == "--n_genes":
            n_genes = int(arg)
        elif opt == "--n_cells":
            n_cells = int(arg)
        elif opt == "--n_time_bins":
            n_time_bins = int(arg)
        elif opt == "--seed":
            seed = int(arg)

    if output_folder is None:
        raise ValueError("--output_folder is required")

    common = SimulatorCommonConfig(
        output_folder=Path(output_folder),
        adata_folder=Path(adata_folder) if adata_folder else None,
        n_runs=n_runs,
        ko_genes=ko_genes,
        n_cells=n_cells,
        n_genes=n_genes,
        n_time_bins=n_time_bins,
    )

    return DyngenRunConfig(common=common, seed=seed)


def _write_ground_truth(cfg: DyngenRunConfig) -> None:
    """Write the benchmark GRN as ground-truth .npy files in the True/ directory.

    Uses the GRN defined in common_config.py (same as Harissa/SERGIO/BoolODE).
    Saves a binary (n_genes+1)×(n_genes+1) adjacency matrix per run, where
    index 0 is the stimulus.
    """
    try:
        grn = get_benchmark_grn(cfg.n_genes)
    except ValueError:
        print("[dyngen] No benchmark GRN defined for n_genes={}; skipping ground truth.".format(cfg.n_genes))
        return

    inter = grn.build_true_interaction_matrix()
    true_dir = cfg.output_folder.parent / "True"
    true_dir.mkdir(parents=True, exist_ok=True)

    for run_idx in range(1, cfg.n_runs + 1):
        np.save(true_dir / f"inter_{run_idx}.npy", inter)
    print(f"[dyngen] Wrote ground truth GRN ({cfg.n_genes} genes) to {true_dir}")


def _run_rscript(cfg: DyngenRunConfig) -> None:
    """Execute the R simulation script via Rscript."""
    output_folder_abs = cfg.output_folder.resolve()
    args: List[str] = [
        "Rscript", str(R_SCRIPT),
        "--output_folder", str(output_folder_abs),
        "--n_runs", str(cfg.n_runs),
        "--n_genes", str(cfg.n_genes),
        "--n_cells", str(cfg.n_cells),
        "--n_time_bins", str(cfg.n_time_bins),
        "--seed", str(cfg.seed),
        "--ko_genes", cfg.ko_genes,
    ]

    # Pass benchmark GRN spec if defined for this gene count
    grn_arg = build_dyngen_backbone_spec(cfg.n_genes)
    if grn_arg is not None:
        args.extend(["--benchmark_grn", grn_arg])

    print(f"[dyngen] Running: {' '.join(args)}")
    result = subprocess.run(args, cwd=SCRIPT_DIR, check=False)

    if result.returncode != 0:
        print(
            f"[dyngen] ERROR: R simulation failed with code {result.returncode}",
            file=sys.stderr,
        )
        sys.exit(result.returncode)


def main(argv: Sequence[str]) -> None:
    cfg = parse_cli(argv)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)

    # Write ground truth GRN from benchmark config (shared across all simulators)
    _write_ground_truth(cfg)

    # Run the dyngen R pipeline (handles both WT and KO internally)
    _run_rscript(cfg)

    mode = "KO" if cfg.ko_genes.lower() not in ("", "none") else "WT-only"
    print(f"[dyngen] {mode} simulation completed successfully.")


if __name__ == "__main__":
    main(sys.argv[1:])
