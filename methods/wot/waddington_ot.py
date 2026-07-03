#!/usr/bin/env python3
from __future__ import annotations

"""
waddington_ot.py - Benchmark wrapper for Waddington-OT.

This wrapper uses Waddington-OT transport maps and interpolation to generate
predicted gene-expression clouds at the requested timepoints.

Usage:
    python waddington_ot.py -i <adata.h5ad> -o <output_dir> -u <output_mode>

output_mode:
    full_test  - simulate test timepoints and concatenate the real training cells
    full_train - simulate all training timepoints from the first training interval
    full_full  - simulate all training and test timepoints
    no_traj    - skip prediction
"""

import getopt
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import scipy

import wot

warnings.filterwarnings("ignore")

T_KEY = "timepoint"
N_SIM_CELLS = 100


def _load_adata(inputfile: str) -> ad.AnnData:
    adata = ad.read_h5ad(inputfile)
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    if "dataset_id" in adata.obs.columns and len(np.unique(adata.obs["dataset_id"])) > 1:
        adata = adata[adata.obs["dataset_id"] == "WT"].copy()
    return adata


def _apply_subsampling(adata_train: ad.AnnData, subsample_ids_path: str) -> ad.AnnData:
    if not os.path.exists(subsample_ids_path):
        return adata_train
    subsample_ids = set(np.load(subsample_ids_path, allow_pickle=True).tolist())
    keep_mask = np.array([cell_id in subsample_ids for cell_id in adata_train.obs_names])
    return adata_train[keep_mask].copy()


def _training_timepoints(train_tps_arr: np.ndarray) -> np.ndarray:
    return np.unique(train_tps_arr[:, 1].astype(float))


def _anchor_pair(target_tp: float, train_tps: np.ndarray) -> tuple[float, float]:
    if len(train_tps) < 2:
        raise ValueError("WaddingtonOT requires at least two training timepoints.")
    if target_tp <= train_tps[0]:
        return float(train_tps[0]), float(train_tps[1])
    if target_tp >= train_tps[-1]:
        return float(train_tps[-2]), float(train_tps[-1])
    idx = int(np.searchsorted(train_tps, target_tp, side="right") - 1)
    idx = max(0, min(idx, len(train_tps) - 2))
    return float(train_tps[idx]), float(train_tps[idx + 1])


def _build_ot_model(adata_train: ad.AnnData):
    if "g2" in adata_train.obs.columns:
        return wot.ot.OTModel(adata_train, day_field=T_KEY, growth_rate_field="g2")
    if "cell_growth_rate" in adata_train.obs.columns:
        return wot.ot.OTModel(adata_train, day_field=T_KEY, growth_rate_field="cell_growth_rate")
    return wot.ot.OTModel(adata_train, day_field=T_KEY)


def _predict_timepoint(
    target_tp: float,
    train_tps: np.ndarray,
    pair_cache: dict[tuple[float, float], ad.AnnData],
    adata_train: ad.AnnData,
    n_sim_cells: int,
) -> np.ndarray:
    t0, t1 = _anchor_pair(target_tp, train_tps)
    tmap = pair_cache[(t0, t1)]

    p0 = adata_train[adata_train.obs[T_KEY].astype(float) == t0]
    p1 = adata_train[adata_train.obs[T_KEY].astype(float) == t1]
    if scipy.sparse.issparse(p0.X):
        p0_x = p0.X.toarray()
    else:
        p0_x = np.asarray(p0.X, dtype=np.float64)
    if scipy.sparse.issparse(p1.X):
        p1_x = p1.X.toarray()
    else:
        p1_x = np.asarray(p1.X, dtype=np.float64)

    interp_frac = (target_tp - t0) / (t1 - t0)
    predicted = wot.ot.interpolate_with_ot(p0_x, p1_x, tmap.X, interp_frac, n_sim_cells)
    return np.maximum(predicted, 0)


def main(argv):
    inputfile = ""
    outputfolder = ""
    output_mode = "full_test"

    try:
        opts, _ = getopt.getopt(argv, "hi:o:u:", ["input=", "outputfolder=", "output_mode="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-o", "--outputfolder"):
            outputfolder = arg
        elif opt in ("-u", "--output_mode"):
            output_mode = arg

    if output_mode == "no_traj":
        print("No pred, did not compute anything")
        return

    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    if not os.path.exists(inputfile):
        raise FileNotFoundError(f"Error reading h5ad file: {inputfile}")

    adata = _load_adata(inputfile)
    data_path = os.path.dirname(inputfile)
    train_tps_arr = np.load(os.path.join(data_path, "train_tps.npy"))
    test_tps_arr = np.load(os.path.join(data_path, "test_tps.npy"))

    start_runtime = time.time()

    train_tps = _training_timepoints(train_tps_arr)
    if len(train_tps) < 2:
        raise ValueError("WaddingtonOT requires at least two training timepoints.")

    train_mask = adata.obs[T_KEY].isin(train_tps_arr[:, 1])
    adata_train = adata[train_mask].copy()

    subsample_ids_path = os.path.join(data_path, "subsample_train_ids.npy")
    adata_train = _apply_subsampling(adata_train, subsample_ids_path)

    ot_model = _build_ot_model(adata_train)

    pair_cache: dict[tuple[float, float], ad.AnnData] = {}
    for t0, t1 in zip(train_tps[:-1], train_tps[1:]):
        pair_cache[(float(t0), float(t1))] = ot_model.compute_transport_map(float(t0), float(t1))

    if output_mode == "full_test":
        target_tps = test_tps_arr[:, 1].astype(float)
    elif output_mode == "full_train":
        target_tps = train_tps
    else:
        target_tps = np.unique(np.concatenate([train_tps, test_tps_arr[:, 1].astype(float)]))

    all_pred = []
    all_pred_tps = []
    for tp in target_tps:
        predicted = _predict_timepoint(float(tp), train_tps, pair_cache, adata_train, N_SIM_CELLS)
        all_pred.append(predicted)
        all_pred_tps.extend([float(tp)] * N_SIM_CELLS)

    pred_array = np.vstack(all_pred)
    pred_adata = ad.AnnData(pred_array)
    pred_adata.obs["time"] = all_pred_tps
    pred_adata.obs[T_KEY] = all_pred_tps
    pred_adata.var_names = adata_train.var_names

    if output_mode == "full_test":
        pred_train_data = ad.concat((pred_adata, adata_train))
    else:
        pred_train_data = pred_adata

    pred_train_data.uns["runtime"] = time.time() - start_runtime
    if adata.uns.get("simulation") is not None:
        pred_train_data = pred_train_data[:, 1:].copy()

    save_path_prefix = os.path.splitext(os.path.basename(inputfile))[0]
    pred_train_data.write_h5ad(os.path.join(outputfolder, save_path_prefix + "_adata.h5ad"))


if __name__ == "__main__":
    main(sys.argv[1:])