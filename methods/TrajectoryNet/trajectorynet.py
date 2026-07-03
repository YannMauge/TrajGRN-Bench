#!/usr/bin/env python3
"""
trajectorynet.py - Benchmark wrapper for TrajectoryNet.

TrajectoryNet learns a Continuous Normalizing Flow (CNF) from scRNA-seq
trajectory data. This script:
  1. Reduces gene expression data to PCA space using training cells only.
  2. Trains TrajectoryNet on that embedding.
  3. Samples predicted cells at the requested timepoints by integrating the
     learned ODE from the base (Gaussian) distribution.
  4. Projects the PCA predictions back to gene space and writes the output.

Usage:
    python trajectorynet.py -i <adata.h5ad> -o <output_dir> -u <output_mode>

output_mode:
    full_test  - predict test timepoints; concatenate with actual train cells
    full_train - predict all training timepoints from the base distribution
    full_full  - predict all timepoints (train + test)
    no_traj    - skip prediction
"""

from __future__ import annotations

import getopt
import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

import anndata as ad
import numpy as np
import scipy
import torch
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

T_KEY = "timepoint"
N_SIM_CELLS = 100
N_PCA_DIMS = 50
N_ITER = 50


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _load_adata(inputfile: str) -> ad.AnnData:
    adata = ad.read_h5ad(inputfile)
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    if len(np.unique(adata.obs["dataset_id"])) > 1:
        adata = adata[adata.obs["dataset_id"] == "WT"].copy()
    return adata


def _apply_subsampling(
    adata_train: ad.AnnData, subsample_ids_path: str
) -> ad.AnnData:
    """Keep only cells that appear in subsample_train_ids.npy."""
    if not os.path.exists(subsample_ids_path):
        return adata_train
    subsample_ids = set(
        np.load(subsample_ids_path, allow_pickle=True).tolist()
    )
    keep_mask = np.array([c in subsample_ids for c in adata_train.obs_names])
    return adata_train[keep_mask].copy()


def _build_tn_adata(adata_train: ad.AnnData, pca: PCA) -> ad.AnnData:
    """
    Build the AnnData fed to TrajectoryNet.

    Requirements:
      - adata.obsm contains the PCA embedding under key 'X_pca'.
      - adata.obs['sample_labels'] holds integer time labels (the raw
        timepoint indices, which are already 0-indexed integers in this
        benchmark).
    """
    adata_tn = adata_train.copy()
    adata_tn.obsm["X_pca"] = pca.transform(adata_tn.X)
    adata_tn.obs["sample_labels"] = adata_tn.obs[T_KEY].astype(int)
    return adata_tn


def _train_trajectorynet(
    tmp_adata_path: str, save_dir: str, n_pca_dims: int, n_iter: int
) -> None:
    """Train TrajectoryNet via subprocess (avoids argparse global-state issues)."""
    cmd = [
        sys.executable, "-m", "TrajectoryNet.main",
        "--dataset", tmp_adata_path,
        "--embedding_name", "pca",
        "--max_dim", str(n_pca_dims),
        "--niter", str(n_iter),
        "--save", save_dir,
        "--use_cpu",
    ]
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _load_model(tmp_adata_path: str, save_dir: str, n_pca_dims: int, n_iter: int):
    """Reconstruct the TrajectoryNet model and load the saved checkpoint."""
    from TrajectoryNet import dataset
    from TrajectoryNet.parse import parser
    from TrajectoryNet.train_misc import (
        build_model_tabular,
        create_regularization_fns,
    )

    args_list = [
        "--dataset", tmp_adata_path,
        "--embedding_name", "pca",
        "--max_dim", str(n_pca_dims),
        "--niter", str(n_iter),
        "--save", save_dir,
        "--use_cpu",
    ]
    args = parser.parse_args(args_list)
    args.data = dataset.SCData.factory(args.dataset, args)
    args.timepoints = args.data.get_unique_times()
    # int_tps covers all integer indices from 0 to max(timepoints), inclusive,
    # mapped to continuous time using time_scale (default 0.5).
    args.int_tps = (
        np.arange(max(args.timepoints) + 1) + 1.0
    ) * args.time_scale

    regularization_fns, _ = create_regularization_fns(args)
    model = build_model_tabular(
        args, args.data.get_shape()[0], regularization_fns
    )
    device = torch.device("cpu")
    model.to(device)

    chkpt_path = os.path.join(save_dir, "checkpt.pth")
    chkpt = torch.load(chkpt_path, map_location=device)
    model.load_state_dict(chkpt["state_dict"])
    model.eval()

    return model, args, device


def _sample_at_timepoint(
    tp_val: int,
    model,
    args,
    device: torch.device,
    n_sim_cells: int,
    n_pca_dims: int,
) -> np.ndarray:
    """
    Sample `n_sim_cells` predicted cells at integer timepoint `tp_val`.

    Integration proceeds stepwise from the base Gaussian distribution (t=0)
    forward to the target time, extrapolating beyond training range when needed.
    """
    from TrajectoryNet.main import get_transforms

    # Build integration-time sequence from the first timepoint up to tp_val.
    # Each step is one time_scale unit; we integrate from (it - time_scale) → it.
    target_int_tps = [
        (i + 1.0) * args.time_scale for i in range(tp_val + 1)
    ]

    sample_fn, _ = get_transforms(device, args, model, target_int_tps)
    z_samples = torch.randn(n_sim_cells, n_pca_dims).to(device)

    with torch.no_grad():
        pred_pca = sample_fn(z_samples).cpu().numpy()

    return pred_pca


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv):
    inputfile = ""
    outputfolder = ""
    output_mode = "full_test"

    try:
        opts, _ = getopt.getopt(
            argv, "hi:o:u:", ["input=", "outputfolder=", "output_mode="]
        )
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

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    adata = _load_adata(inputfile)

    path_data_folder = os.path.dirname(inputfile)
    train_tps_arr = np.load(os.path.join(path_data_folder, "train_tps.npy"))
    test_tps_arr = np.load(os.path.join(path_data_folder, "test_tps.npy"))

    start_runtime = time.time()

    n_genes = adata.shape[1]
    n_pca_dims = min(N_PCA_DIMS, n_genes)

    # ------------------------------------------------------------------
    # Prepare training data
    # ------------------------------------------------------------------
    train_mask = adata.obs[T_KEY].isin(train_tps_arr[:, 1])
    adata_train = adata[train_mask].copy()

    subsample_ids_path = os.path.join(path_data_folder, "subsample_train_ids.npy")
    adata_train = _apply_subsampling(adata_train, subsample_ids_path)

    # Fit PCA on training cells only.
    pca = PCA(n_components=n_pca_dims)
    pca.fit(adata_train.X)

    adata_tn = _build_tn_adata(adata_train, pca)

    # ------------------------------------------------------------------
    # Train TrajectoryNet
    # ------------------------------------------------------------------
    tmp_dir = tempfile.mkdtemp(prefix="trajectorynet_")
    try:
        tmp_adata_path = os.path.join(tmp_dir, "input.h5ad")
        save_dir = os.path.join(tmp_dir, "model")
        adata_tn.write_h5ad(tmp_adata_path)

        _train_trajectorynet(tmp_adata_path, save_dir, n_pca_dims, N_ITER)

        # ------------------------------------------------------------------
        # Load model and sample predictions
        # ------------------------------------------------------------------
        model, args, device = _load_model(
            tmp_adata_path, save_dir, n_pca_dims, N_ITER
        )

        # Determine which timepoints to predict.
        if output_mode == "full_test":
            target_tps = test_tps_arr[:, 1]
        elif output_mode == "full_train":
            target_tps = train_tps_arr[:, 1]
        else:  # full_full
            target_tps = np.unique(
                np.concatenate([train_tps_arr[:, 1], test_tps_arr[:, 1]])
            )

        all_pred_pca = []
        all_pred_tps = []

        for tp in target_tps:
            pred_pca = _sample_at_timepoint(
                int(tp), model, args, device, N_SIM_CELLS, n_pca_dims
            )
            all_pred_pca.append(pred_pca)
            all_pred_tps.extend([tp] * N_SIM_CELLS)

        # Project PCA predictions back to gene space.
        all_pred_pca_arr = np.vstack(all_pred_pca)
        all_pred_gene = pca.inverse_transform(all_pred_pca_arr)
        all_pred_gene = np.maximum(all_pred_gene, 0)  # avoid negative expression
    except Exception as exc:
        print(f"TrajectoryNet failed: {exc}", file=sys.stderr)
        return 1
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Build and save output AnnData
    # ------------------------------------------------------------------
    pred_adata = ad.AnnData(all_pred_gene)
    pred_adata.obs["time"] = all_pred_tps
    pred_adata.obs[T_KEY] = all_pred_tps
    pred_adata.var_names = adata_train.var_names

    if output_mode == "full_test":
        pred_train_data = ad.concat((pred_adata, adata_train))
    else:
        pred_train_data = pred_adata

    pred_train_data.uns["runtime"] = time.time() - start_runtime

    if adata.uns.get("simulation") is not None:
        # Remove the first gene (Stimulus) if data is from simulation.
        pred_train_data = pred_train_data[:, 1:].copy()

    save_path_prefix = os.path.splitext(os.path.basename(inputfile))[0]
    pred_train_data.write_h5ad(
        os.path.join(outputfolder, save_path_prefix + "_adata.h5ad")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
