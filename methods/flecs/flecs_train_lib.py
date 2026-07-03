import os
import time
import warnings

import anndata as ad
import numpy as np
import scanpy as sc
import scipy
import torch
from joblib import Parallel, delayed
from torch.utils.data import DataLoader
from tqdm import tqdm

from flecs.sc.dataset import Paul15Dataset
from flecs.sc.model import GRNCellPop
from flecs.sc.utils import compute_cell_knn_paths, train_epoch
from flecs.trajectory import simulate_deterministic_trajectory
from flecs.utils import set_seed

for module in ["anndata", "scipy", "torchdiffeq"]:
    warnings.filterwarnings("ignore", module=module)


def deterministic_pred(mycellpop, adata, device, t_start, sim_path_length, n_timepoints, batch_divisor):
    t_key = "timepoint"
    start_adata = adata[adata.obs[t_key] == t_start]
    batch_size = max(1, start_adata.shape[0] // batch_divisor)
    mycellpop.change_batch_size(batch_size)

    all_traj = []
    with torch.no_grad():
        for i in tqdm(range(start_adata.shape[0] // batch_size)):
            gene_expr = torch.tensor(
                start_adata.X[i * batch_size : (i + 1) * batch_size, :, None], dtype=torch.float32
            ).to(device)
            mycellpop.set_visible_state(gene_expr)
            traj = simulate_deterministic_trajectory(
                mycellpop, torch.linspace(t_start, t_start + sim_path_length, n_timepoints)
            )
            traj = traj[:, :, : mycellpop.n_genes, 0]
            all_traj.append(traj)

    all_traj_cat = torch.cat(all_traj)
    obs_times = torch.cat(
        [torch.linspace(t_start, t_start + sim_path_length, n_timepoints)[None]] * len(all_traj_cat)
    ).reshape(-1)
    traj_adata = ad.AnnData(all_traj_cat.reshape(-1, adata.shape[1]).cpu().numpy(), obs=obs_times.numpy())
    traj_adata.obs = traj_adata.obs.rename(columns={0: "time"})
    return traj_adata


def parse_ko_targets(ko_output_genes, var_names):
    ko_targets = []
    if ko_output_genes.lower() == "all":
        ko_targets = [gene for gene in var_names if gene != "Stimulus"]
    elif ko_output_genes.lower() not in ["none", ""]:
        ko_targets = [gene.strip() for gene in ko_output_genes.split(",") if gene.strip()]

    invalid_kos = [gene for gene in ko_targets if gene not in var_names]
    if invalid_kos:
        raise ValueError(f"Unknown KO gene(s): {invalid_kos}")

    return ko_targets


def compute_grn_adj_mat(adata_train, t_key):
    unique_t_values = np.unique(adata_train.obs[t_key])
    n_genes = adata_train.shape[1]

    def compute_pearson_matrix(t):
        subset = adata_train[adata_train.obs[t_key] == t].X
        p_matrix_t = np.zeros((n_genes, n_genes))
        for i in range(n_genes):
            for j in range(n_genes):
                x = subset[:, i]
                y = subset[:, j]
                if np.all(x == x[0]) or np.all(y == y[0]):
                    p_matrix_t[i, j] = 0
                else:
                    p_matrix_t[i, j] = scipy.stats.pearsonr(x, y).statistic
        return p_matrix_t

    list_p_matrices = Parallel(n_jobs=-1)(delayed(compute_pearson_matrix)(t) for t in unique_t_values)
    p_matrices_3d = np.stack(list_p_matrices)
    grn_adj_mat = np.max(np.abs(p_matrices_3d), axis=0)
    grn_adj_mat = (
        np.sign(
            np.take_along_axis(
                p_matrices_3d,
                np.expand_dims(np.argmax(np.abs(p_matrices_3d), axis=0), axis=0),
                axis=0,
            )[0]
        )
        * grn_adj_mat
    )
    return grn_adj_mat + 1e-6


def build_knn_paths(adata_train, t_key, train_tps):
    start_cells = list(
        adata_train.obs.index.get_indexer_for(
            adata_train.obs[adata_train.obs[t_key] == train_tps[:, 1][0]].index
        )
    )
    end_cells = list(
        adata_train.obs.index.get_indexer_for(
            adata_train.obs[adata_train.obs[t_key] == train_tps[:, 1][-1]].index
        )
    )
    unsorted_shortest_paths = compute_cell_knn_paths(
        adata_train, start_cells, end_cells, option="PCA", n_bins=None
    )
    adata_train.uns["unsorted_shortest_paths"] = unsorted_shortest_paths
    return unsorted_shortest_paths


def build_force_timepoints_paths(adata_train, t_key, train_tps):
    import ot

    end_cells = list(
        adata_train.obs.index.get_indexer_for(
            adata_train.obs[adata_train.obs[t_key] == train_tps[:, 1][-1]].index
        )
    )

    ot_cell_paths = [[cell] for cell in end_cells]
    for t in tqdm(range(len(train_tps) - 1)):
        t_curr = train_tps[:, 1][-1 - t]
        t_prev = train_tps[:, 1][-2 - t]
        pop_size1 = len(adata_train[adata_train.obs[t_key] == t_curr])
        pop_size2 = len(adata_train[adata_train.obs[t_key] == t_prev])
        dist_mat = scipy.spatial.distance_matrix(
            adata_train[adata_train.obs[t_key] == t_curr].X,
            adata_train[adata_train.obs[t_key] == t_prev].X,
            p=2,
        )
        transport = ot.emd([1 / pop_size1] * pop_size1, [1 / pop_size2] * pop_size2, dist_mat)

        for i in range(len(ot_cell_paths)):
            current_cell = ot_cell_paths[i][0]
            cell_iloc = adata_train[adata_train.obs[t_key] == t_curr].obs.index.get_loc(int(current_cell))
            prev_cell_iloc = np.argmax(transport[cell_iloc])
            prev_cell_idx = int(adata_train[adata_train.obs[t_key] == t_prev].obs.iloc[prev_cell_iloc].name)
            ot_cell_paths[i].insert(0, prev_cell_idx)

    ot_cell_paths_dict = {path[-1]: path for path in ot_cell_paths}
    unsorted_shortest_paths = {str(k): v for k, v in ot_cell_paths_dict.items()}
    adata_train.uns["unsorted_shortest_paths"] = unsorted_shortest_paths
    return unsorted_shortest_paths


def choose_knn_path_length(adata_train, unsorted_shortest_paths, train_tps):
    path_length = 2
    unsorted_dataset = Paul15Dataset(adata_train, unsorted_shortest_paths, path_length=path_length)
    n_max_traj = len(unsorted_dataset)

    path_length = int(len(train_tps[:, 1]) / 2)
    unsorted_dataset = Paul15Dataset(adata_train, unsorted_shortest_paths, path_length=path_length)
    n_traj = len(unsorted_dataset)
    while n_traj / n_max_traj < 0.5:
        path_length -= 1
        unsorted_dataset = Paul15Dataset(adata_train, unsorted_shortest_paths, path_length=path_length)
        n_traj = len(unsorted_dataset)

    return path_length, unsorted_dataset


def train_flecs_model(adata_train, unsorted_dataset, path_length):
    train_len = len(unsorted_dataset) - 1
    valid_len = 1
    learning_rate = 0.005
    batch_size = 5

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if torch.cuda.is_available():
        torch.cuda.set_device(device)

    test_len = len(unsorted_dataset) - train_len - valid_len
    if test_len < 0:
        raise ValueError("Not enough samples for these train / test len")

    train_dataset, _, _ = torch.utils.data.random_split(
        unsorted_dataset,
        [train_len, valid_len, test_len],
        generator=torch.Generator().manual_seed(0),
    )

    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    mycellpop = GRNCellPop(adata=adata_train, batch_size=batch_size, n_latent_var=0, use_2nd_order_interactions=False).to(device)

    optimizer = torch.optim.Adam(mycellpop.parameters(), lr=learning_rate)
    loss = torch.nn.MSELoss()

    for epoch in range(50):
        print("New epoch", epoch)
        train_epoch(mycellpop, train_dataloader, optimizer, path_length=path_length, loss=loss)

    return mycellpop, device


def compute_simulation_parameters(output_mode, train_tps, test_tps):
    if output_mode == "full_test":
        t_start = train_tps[:, 1][-1]
        target_tps = test_tps[:, 1]
        n_sim = len(target_tps) + 1
        sim_path_length = target_tps[-1] - t_start
    elif output_mode == "full_train":
        t_start = train_tps[:, 1][0]
        target_tps = train_tps[:, 1]
        n_sim = len(target_tps)
        sim_path_length = target_tps[-1] - t_start
    else:
        all_sim_tps = np.unique(np.concatenate([train_tps[:, 1], test_tps[:, 1]]))
        t_start = all_sim_tps[0]
        target_tps = all_sim_tps
        n_sim = len(target_tps)
        sim_path_length = target_tps[-1] - t_start

    return t_start, target_tps, n_sim, sim_path_length


def simulate_and_save(
    mycellpop,
    adata,
    adata_train,
    device,
    outputfolder,
    save_path_prefix,
    output_suffix,
    grn_weights_flat,
    n_genes,
    start_runtime,
    output_mode,
    t_start,
    sim_path_length,
    n_sim,
    sim_adata,
    batch_divisor,
):
    t_key = "timepoint"

    with torch.no_grad():
        mycellpop["gene", "regulates", "gene"].simple_conv_weights.data[0, :, 0] = torch.tensor(
            grn_weights_flat,
            dtype=mycellpop["gene", "regulates", "gene"].simple_conv_weights.data.dtype,
        ).to(device)

    if output_mode != "no_traj":
        pred_data = deterministic_pred(
            mycellpop=mycellpop,
            adata=sim_adata,
            device=device,
            t_start=t_start,
            sim_path_length=sim_path_length,
            n_timepoints=n_sim,
            batch_divisor=batch_divisor,
        )
        pred_data.var_names = adata_train.var_names
        pred_data.obs[t_key] = pred_data.obs["time"]

        if output_mode == "full_test":
            pred_data = pred_data[pred_data.obs["time"] != t_start].copy()
            pred_train_data = ad.concat((pred_data, adata_train))
        else:
            pred_train_data = pred_data

        pred_train_data.uns["runtime"] = time.time() - start_runtime
        if adata.uns.get("simulation") is not None:
            pred_train_data = pred_train_data[:, 1:].copy()
        ad.settings.allow_write_nullable_strings = True
        pred_train_data.write_h5ad(os.path.join(outputfolder, save_path_prefix + output_suffix + "_adata.h5ad"))

    inter_flecs = np.reshape(grn_weights_flat, (n_genes, n_genes, 1))
    np.save(os.path.join(outputfolder, save_path_prefix + output_suffix + "_GRN"), inter_flecs)


def run_training(
    inputfile,
    outputfolder,
    output_mode,
    ko_output_genes,
    path_strategy,
    batch_divisor,
):
    set_seed(0)
    t_key = "timepoint"

    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    if not os.path.exists(inputfile):
        raise FileNotFoundError("Error reading h5ad file")

    adata = ad.read_h5ad(inputfile)
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    path_data_folder = os.path.dirname(inputfile)
    train_tps = np.load(os.path.join(path_data_folder, "train_tps.npy"))
    test_tps = np.load(os.path.join(path_data_folder, "test_tps.npy"))

    adata_train = adata[adata.obs[t_key].isin(train_tps[:, 1]), :].copy()
    if len(np.unique(adata_train.obs["dataset_id"])) > 1:
        adata_train = adata_train[adata_train.obs["dataset_id"] == "WT"].copy()

    if path_strategy == "force_timepoints":
        adata_train.obs.reset_index(inplace=True)

    subsample_ids_path = os.path.join(path_data_folder, "subsample_train_ids.npy")
    subsample_ids = None
    if os.path.exists(subsample_ids_path):
        subsample_ids = np.load(subsample_ids_path, allow_pickle=True)
        adata_train = adata_train[adata_train.obs_names.isin(subsample_ids)].copy()

    start_runtime = time.time()

    adata_train.varp["grn_adj_mat"] = compute_grn_adj_mat(adata_train, t_key)

    sc.tl.pca(adata_train)
    sc.pp.neighbors(adata_train)
    sc.tl.umap(adata_train)

    if path_strategy == "force_timepoints":
        unsorted_shortest_paths = build_force_timepoints_paths(adata_train, t_key, train_tps)
        path_length = len(np.unique(adata_train.obs[t_key]))
        unsorted_dataset = Paul15Dataset(adata_train, unsorted_shortest_paths, path_length=path_length)
    else:
        unsorted_shortest_paths = build_knn_paths(adata_train, t_key, train_tps)
        path_length, unsorted_dataset = choose_knn_path_length(adata_train, unsorted_shortest_paths, train_tps)

    mycellpop, device = train_flecs_model(adata_train, unsorted_dataset, path_length)

    save_path_prefix = os.path.splitext(os.path.basename(inputfile))[0]
    n_genes = adata_train.shape[1]

    ko_targets = parse_ko_targets(ko_output_genes, adata_train.var_names)

    sim_adata = adata_train
    t_start = None
    n_sim = None
    sim_path_length = None
    if output_mode != "no_traj":
        t_start, _, n_sim, sim_path_length = compute_simulation_parameters(output_mode, train_tps, test_tps)
        if subsample_ids is not None and t_start == train_tps[:, 1][0]:
            sim_adata = adata[~adata.obs_names.isin(subsample_ids)].copy()

    original_weights = mycellpop["gene", "regulates", "gene"].simple_conv_weights.data[0, :, 0].cpu().numpy().copy()

    simulate_and_save(
        mycellpop=mycellpop,
        adata=adata,
        adata_train=adata_train,
        device=device,
        outputfolder=outputfolder,
        save_path_prefix=save_path_prefix,
        output_suffix="",
        grn_weights_flat=original_weights,
        n_genes=n_genes,
        start_runtime=start_runtime,
        output_mode=output_mode,
        t_start=t_start,
        sim_path_length=sim_path_length,
        n_sim=n_sim,
        sim_adata=sim_adata,
        batch_divisor=batch_divisor,
    )

    for ko_gene in ko_targets:
        ko_idx = list(adata_train.var_names).index(ko_gene)
        ko_weights = original_weights.copy()
        ko_weights[ko_idx::n_genes] = 0
        simulate_and_save(
            mycellpop=mycellpop,
            adata=adata,
            adata_train=adata_train,
            device=device,
            outputfolder=outputfolder,
            save_path_prefix=save_path_prefix,
            output_suffix=f"_ko_{ko_gene}",
            grn_weights_flat=ko_weights,
            n_genes=n_genes,
            start_runtime=start_runtime,
            output_mode=output_mode,
            t_start=t_start,
            sim_path_length=sim_path_length,
            n_sim=n_sim,
            sim_adata=sim_adata,
            batch_divisor=batch_divisor,
        )

    with torch.no_grad():
        mycellpop["gene", "regulates", "gene"].simple_conv_weights.data[0, :, 0] = torch.tensor(
            original_weights,
            dtype=mycellpop["gene", "regulates", "gene"].simple_conv_weights.data.dtype,
        ).to(device)
