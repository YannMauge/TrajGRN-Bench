#!/usr/bin/env python3
"""
inference_renge.py - RENGE GRN inference and expression prediction wrapper.

RENGE (REgulatory Network inference using GEne perturbation data) infers
gene regulatory networks from time-series single-cell data, optionally
leveraging CRISPR perturbation (KO) information when available.

Usage:
    python inference_renge.py -i <input_h5ad> -o <output_folder>
                              [-u <output_mode>] [-p <perturbation_training>]
                              [-k <ko_genes>]

Input:
    -i, --input                  Path to input H5AD file
    -o, --outputfolder           Path to output directory

Optional:
    -u, --output_mode            Output simulation mode: full_test, full_train,
                                 full_full, no_traj (default: full_test)
    -p, --perturbation_training  Use KO perturbation groups during training
                                 (default: false)
    -k, --ko_genes               KO target genes for output prediction
                                 ("none", "all", or comma-separated; default: none)

Output:
    <data_id>_GRN.npy            G x G regulatory coefficient matrix
    <data_id>_adata.h5ad         Predicted expression (if output_mode != no_traj)
"""

import sys
import getopt
import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy

import warnings
for module in ["anndata", "scipy", "renge", "scanpy"]:
    warnings.filterwarnings("ignore", module=module)

from renge import Renge


def _parse_bool(value: str) -> bool:
    return str(value).strip().lower() in ("true", "1", "yes", "y")


def _prepare_renge_input(
    adata_train,
    t_key: str = "timepoint",
):
    """
    Prepare X (KO indicators + time) and E (expression) DataFrames for RENGE.

    RENGE expects:
      - X: C x (G+1) DataFrame. First G columns are one-hot KO indicators
           (all zeros for WT). Last column is the sampling time.
      - E: C x G DataFrame of normalized, log-transformed expression.

    Returns (X, E, gene_names).
    """
    if scipy.sparse.issparse(adata_train.X):
        expr = adata_train.X.toarray()
    else:
        expr = adata_train.X

    # Normalize and log-transform (RENGE expects this)
    adata_norm = adata_train.copy()
    sc.pp.normalize_total(adata_norm)
    sc.pp.log1p(adata_norm)
    expr_log = adata_norm.X
    if scipy.sparse.issparse(expr_log):
        expr_log = expr_log.toarray()

    gene_names = list(adata_train.var_names)
    n_cells = adata_train.n_obs
    n_genes = len(gene_names)

    # Build X: KO indicator columns (all zeros for WT) + time column
    X_data = np.zeros((n_cells, n_genes + 1))

    # Fill KO columns from dataset_id if available
    if "dataset_id" in adata_train.obs.columns:
        for i, ds_id in enumerate(adata_train.obs["dataset_id"]):
            ds_str = str(ds_id)
            if ds_str != "WT" and ds_str in gene_names:
                gene_idx = gene_names.index(ds_str)
                X_data[i, gene_idx] = 1

    # Fill time column (last column)
    if t_key in adata_train.obs.columns:
        X_data[:, -1] = adata_train.obs[t_key].values
    else:
        raise KeyError(f"Time key '{t_key}' not found in adata.obs")

    X_cols = gene_names + ["time"]
    X = pd.DataFrame(X_data, index=adata_train.obs_names, columns=X_cols)
    E = pd.DataFrame(expr_log, index=adata_train.obs_names, columns=gene_names)

    return X, E, gene_names


def _build_prediction_input(
    gene_names,
    timepoints,
    n_cells_per_tp: int = 1,
    ko_gene=None,
):
    """
    Build X_pred DataFrame for predicting expression at given timepoints.

    RENGE's predict() expects X with real-valued KO columns (negative for KO)
    and a time column. For WT prediction, KO columns are all zeros.
    """
    n_genes = len(gene_names)
    ko_idx = gene_names.index(ko_gene) if ko_gene is not None else None
    rows = []
    for t in timepoints:
        for _ in range(n_cells_per_tp):
            row = np.zeros(n_genes + 1)
            if ko_idx is not None:
                row[ko_idx] = -1.0
            row[-1] = float(t)
            rows.append(row)

    X_cols = gene_names + ["time"]
    X_pred = pd.DataFrame(rows, columns=X_cols)
    return X_pred


def _parse_ko_targets(ko_output_genes: str, gene_names):
    if ko_output_genes.lower() == "all":
        ko_targets = [gene for gene in gene_names if gene != "Stimulus"]
    elif ko_output_genes.lower() not in ["none", ""]:
        ko_targets = [gene.strip() for gene in ko_output_genes.split(",") if gene.strip()]
    else:
        ko_targets = []

    invalid_kos = [gene for gene in ko_targets if gene not in gene_names]
    if invalid_kos:
        raise ValueError(f"Unknown KO gene(s): {invalid_kos}")
    return ko_targets


def main(argv):
    t_key = "timepoint"
    output_mode = "full_test"
    perturbation_training = False
    ko_output_genes = "none"

    inputfile = ""
    outputfolder = ""

    try:
        opts, args = getopt.getopt(
            argv,
            "hi:o:u:p:k:",
            ["input=", "outputfolder=", "output_mode=",
             "perturbation_training=", "ko_genes="],
        )
    except getopt.GetoptError as err:
        print(err, file=sys.stderr)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = str(arg)
        if opt in ("-o", "--outputfolder"):
            outputfolder = str(arg)
        if opt in ("-u", "--output_mode"):
            output_mode = str(arg)
        if opt in ("-p", "--perturbation_training"):
            perturbation_training = _parse_bool(arg)
        if opt in ("-k", "--ko_genes"):
            ko_output_genes = str(arg)

    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    # --- Load data ---
    data_path = os.path.join(inputfile)
    if os.path.exists(data_path):
        adata = ad.read_h5ad(data_path)
    else:
        raise FileNotFoundError(f"Input file not found: {data_path}")

    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    # --- Filter to training timepoints ---
    path_data_folder = os.path.dirname(inputfile)
    train_tps_path = os.path.join(path_data_folder, "train_tps.npy")
    if not os.path.exists(train_tps_path):
        raise FileNotFoundError(f"train_tps.npy not found in {path_data_folder}")
    train_tps = np.load(train_tps_path)
    train_tp_values = train_tps[:, 1] if train_tps.ndim > 1 else train_tps

    adata_train = adata[adata.obs[t_key].isin(train_tp_values), :].copy()

    # --- Handle perturbation training ---
    has_ko = "dataset_id" in adata_train.obs.columns and \
             len(np.unique(adata_train.obs["dataset_id"])) > 1

    if not perturbation_training and has_ko:
        adata_train = adata_train[
            adata_train.obs["dataset_id"] == "WT"
        ].copy()

    if perturbation_training and "dataset_id" in adata_train.obs.columns:
        print("Training with perturbation data")
    elif has_ko:
        print("Training on WT only")

    # Apply cell subsampling if subsample_train_ids.npy is present
    subsample_ids_path = os.path.join(path_data_folder, "subsample_train_ids.npy")
    if os.path.exists(subsample_ids_path):
        subsample_ids = np.load(subsample_ids_path, allow_pickle=True)
        adata_train = adata_train[
            adata_train.obs_names.isin(subsample_ids)
        ].copy()

    print(f"Training on {adata_train.n_obs} cells, {adata_train.n_vars} genes")
    print(f"Perturbation training: {perturbation_training}")

    # --- Prepare RENGE input ---
    X, E, gene_names = _prepare_renge_input(adata_train, t_key=t_key)
    ko_targets = _parse_ko_targets(ko_output_genes, gene_names)

    # --- Run RENGE ---
    start_runtime = time.time()
    print("Starting RENGE model training (hyperparameter optimization + fitting)...")

    reg = Renge()
    A = reg.estimate_hyperparams_and_fit(X, E, n_trials=3)
    print(f"RENGE training completed in {time.time() - start_runtime:.1f}s")

    # --- Save GRN ---
    save_path_prefix = os.path.splitext(os.path.basename(inputfile))[0]
    # A is G x G DataFrame: A[i, j] = regulatory coefficient from gene j to gene i
    grn = A.values  # numpy array
    np.save(os.path.join(outputfolder, f"{save_path_prefix}_GRN"), grn)
    print(f"Saved GRN: {save_path_prefix}_GRN.npy")

    if ko_targets:
        for ko_gene in ko_targets:
            ko_idx = gene_names.index(ko_gene)
            ko_grn = grn.copy()
            ko_grn[:, ko_idx] = 0
            np.save(
                os.path.join(outputfolder, f"{save_path_prefix}_ko_{ko_gene}_GRN"),
                ko_grn,
            )
            print(f"Saved KO GRN: {save_path_prefix}_ko_{ko_gene}_GRN.npy")

    # --- Predict expression (if requested) ---
    if output_mode != "no_traj":
        test_tps_path = os.path.join(path_data_folder, "test_tps.npy")
        test_tp_values = None
        if os.path.exists(test_tps_path):
            test_tps = np.load(test_tps_path)
            test_tp_values = test_tps[:, 1] if test_tps.ndim > 1 else test_tps

        # Build prediction input for relevant timepoints
        if output_mode == "full_test" and test_tp_values is not None:
            pred_timepoints = test_tp_values
        elif output_mode == "full_train":
            pred_timepoints = train_tp_values
        elif output_mode == "full_full":
            all_tps = np.unique(adata.obs[t_key].values)
            pred_timepoints = all_tps
        else:
            pred_timepoints = np.array([]) if test_tp_values is None else test_tp_values

        n_cells_per_tp = max(1, adata_train.n_obs // max(1, len(train_tp_values)))

        def simulate_and_save(X_pred, output_suffix: str, dataset_id: str):
            E_pred = reg.predict(X_pred)
            print(f"Predicted expression for timepoints : {pred_timepoints} ")

            pred_obs = pd.DataFrame({
                "timepoint": np.repeat(pred_timepoints, n_cells_per_tp),
                "time": np.repeat(pred_timepoints, n_cells_per_tp),
                "dataset_id": dataset_id,
            })
            pred_adata = ad.AnnData(
                X=np.array(E_pred),
                obs=pred_obs,
                var=pd.DataFrame(index=gene_names),
            )

            if output_mode == "full_test":
                pred_adata = ad.concat((pred_adata, adata_train))

            pred_adata.uns["runtime"] = time.time() - start_runtime
            if adata.uns.get("simulation") is not None:
                pred_adata = pred_adata[:, 1:].copy()

            pred_adata.write_h5ad(
                os.path.join(outputfolder, f"{save_path_prefix}{output_suffix}_adata.h5ad")
            )
            print(f"Saved predicted expression: {save_path_prefix}{output_suffix}_adata.h5ad")

        if len(pred_timepoints) > 0:
            try:
                X_pred = _build_prediction_input(
                    gene_names,
                    pred_timepoints,
                    n_cells_per_tp=n_cells_per_tp,
                )
                simulate_and_save(X_pred, "", "WT")

                for ko_gene in ko_targets:
                    X_pred_ko = _build_prediction_input(
                        gene_names,
                        pred_timepoints,
                        n_cells_per_tp=n_cells_per_tp,
                        ko_gene=ko_gene,
                    )
                    simulate_and_save(X_pred_ko, f"_ko_{ko_gene}", ko_gene)
            except Exception as exc:
                print(f"Warning: Expression prediction failed: {exc}", file=sys.stderr)
        else:
            print("No timepoints to predict; skipping expression prediction.")

    print("RENGE inference completed successfully.")


if __name__ == "__main__":
    main(sys.argv[1:])
