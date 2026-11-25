#!/usr/bin/env python3
"""
dynGENIE3_ontrain.py — Benchmark entrypoint for dynGENIE3.

dynGENIE3 is a semi-parametric approach for GRN inference from time-series
expression data (Huynh-Thu & Geurts, Scientific Reports 2018).

Usage:
    python dynGENIE3_ontrain.py -i <input_h5ad> -o <output_folder>

The script aggregates scRNA-seq cells by timepoint (mean expression) to
produce pseudo-bulk time series suitable for dynGENIE3's finite-difference
formulation, then saves the VIM matrix as `<data_id>_GRN.npy`.
"""

import numpy as np
import scipy
import scanpy as sc
import anndata as ad
import sys
import getopt
import os

from dynGENIE3 import dynGENIE3

import warnings
for module in ["anndata", "scipy", "scanpy"]:
    warnings.filterwarnings("ignore", module=module)


def main(argv):
    t_key = "timepoint"
    lognorm = 1
    inputfile = ""
    outputfolder = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["input=", "outputfolder="])
    except getopt.GetoptError as err:
        print(err, file=sys.stderr)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = arg
        if opt in ("-o", "--outputfolder"):
            outputfolder = arg

    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    # ── Load data ───────────────────────────────────────────────────────
    data_path = os.path.join(inputfile)
    if os.path.exists(data_path):
        adata = ad.read_h5ad(data_path)
    else:
        raise FileNotFoundError(f"Error reading h5ad file: {data_path}")

    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    # ── Filter to training timepoints ───────────────────────────────────
    path_data_folder = os.path.dirname(inputfile)
    train_tps = np.load(os.path.join(path_data_folder, "train_tps.npy"))[:, 1]
    adata = adata[adata.obs[t_key].isin(train_tps)].copy()

    # ── Normalize (optional, matching GENIE3 baseline) ───────────────────
    if lognorm:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    # ── Aggregate cells per timepoint (mean) for pseudo-bulk time series ─
    # dynGENIE3 expects TS_data as a list of arrays, each of shape
    # (n_timepoints, n_genes), with one row per time point.
    unique_tps = sorted(adata.obs[t_key].unique())
    mean_expr = np.zeros((len(unique_tps), adata.n_vars))
    for i, tp in enumerate(unique_tps):
        mask = adata.obs[t_key] == tp
        mean_expr[i, :] = adata.X[mask].mean(axis=0)

    TS_data = [mean_expr]
    time_points = [np.array(unique_tps, dtype=np.float64)]

    # ── Run dynGENIE3 ───────────────────────────────────────────────────
    VIM, alphas, _, _, _ = dynGENIE3(TS_data, time_points, alpha="from_data")

    # ── Save GRN ────────────────────────────────────────────────────────
    data_id = os.path.splitext(os.path.basename(inputfile))[0]
    np.save(os.path.join(outputfolder, f"{data_id}_GRN"), VIM)

    # ── Write .done marker ──────────────────────────────────────────────
    done_path = os.path.join(outputfolder, f"{data_id}.done")
    open(done_path, "w").close()

    print(f"dynGENIE3 completed. GRN saved to {outputfolder}")


if __name__ == "__main__":
    main(sys.argv[1:])
