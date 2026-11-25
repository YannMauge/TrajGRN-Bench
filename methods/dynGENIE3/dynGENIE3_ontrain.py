#!/usr/bin/env python3
"""dynGENIE3_ontrain.py — benchmark entrypoint for dynGENIE3
(Huynh-Thu & Geurts, Scientific Reports 2018; vendored upstream code).

dynGENIE3 needs time series of the same samples, which snapshot scRNA-seq
does not provide. We build ``N_SERIES`` bootstrapped pseudo-bulk series: for
each series and each training timepoint, the mean (log-normalised) expression
of a random subset (``SUBSET_FRACTION``, with replacement) of the WT cells at
that timepoint. dynGENIE3 natively accepts several series, so each tree sees
N_SERIES × (T − 1) transitions instead of T − 1 with a single averaged series.

dynGENIE3 returns VIM[i, j] = importance of gene i for target j
(regulator × target), the benchmark orientation. It has no seed argument:
numpy's global RNG (used by its random forests) is seeded before the run.

Usage:
    python dynGENIE3_ontrain.py -i <input_h5ad> -o <output_folder>
"""
import getopt
import os
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dynGENIE3 import dynGENIE3  # noqa: E402
from _common import (  # noqa: E402
    T_KEY, load_training_data, log_normalize_genes, save_grn, seed_everything, variable_columns,
)

warnings.filterwarnings("ignore")

SEED = 0
N_SERIES = 20
SUBSET_FRACTION = 0.5


def bootstrap_series(x, times, n_series=N_SERIES, fraction=SUBSET_FRACTION, seed=SEED):
    """Return (TS_data, time_points) for dynGENIE3 from snapshot cells."""
    rng = np.random.default_rng(seed)
    tps = np.sort(np.unique(times))
    if len(tps) < 2:
        raise ValueError(f"dynGENIE3 needs at least 2 training timepoints, got {tps.tolist()}")
    idx_by_tp = [np.flatnonzero(times == tp) for tp in tps]
    series = []
    for _ in range(n_series):
        rows = []
        for idx in idx_by_tp:
            k = max(1, int(round(fraction * len(idx))))
            rows.append(x[rng.choice(idx, size=k, replace=True)].mean(axis=0))
        series.append(np.vstack(rows))
    return series, [tps.astype(np.float64) for _ in range(n_series)]


def main(argv):
    inputfile, outputfolder = "", ""
    opts, _ = getopt.getopt(argv, "hi:o:", ["input=", "outputfolder="])
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-o", "--outputfolder"):
            outputfolder = arg
    if not os.path.exists(inputfile):
        raise FileNotFoundError(f"Error reading h5ad file: {inputfile}")

    seed_everything(SEED)
    adata, _, _ = load_training_data(inputfile)
    x = log_normalize_genes(adata.X, adata.var_names)
    # Constant variables (e.g. a knocked-out gene) carry no signal and break
    # the per-target normalisation of the random forests: score them 0.
    keep = variable_columns(x)
    ts_data, time_points = bootstrap_series(x[:, keep], adata.obs[T_KEY].to_numpy(dtype=float))

    VIM_sub, _, _, _, _ = dynGENIE3(ts_data, time_points, alpha="from_data")
    VIM = np.zeros((x.shape[1], x.shape[1]))
    VIM[np.ix_(keep, keep)] = VIM_sub
    path = save_grn(outputfolder, inputfile, VIM, adata.var_names)
    print(f"dynGENIE3 completed. GRN saved to {path}")


if __name__ == "__main__":
    main(sys.argv[1:])
