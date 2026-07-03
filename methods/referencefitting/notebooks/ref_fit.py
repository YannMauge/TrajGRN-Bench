import os
num_threads = "8"
os.environ["OMP_NUM_THREADS"] = num_threads
os.environ["OPENBLAS_NUM_THREADS"] = num_threads
os.environ["MKL_NUM_THREADS"] = num_threads
os.environ["VECLIB_MAXIMUM_THREADS"] = num_threads
os.environ["NUMEXPR_NUM_THREADS"] = num_threads

import numpy as np
import scanpy as sc
import anndata as ad
import scipy
import os
import sys, getopt
import torch
import time
sys.path.append("./methods/referencefitting/src/")
import rf
import util
import ot
import warnings
for module in ["anndata"]:
    warnings.filterwarnings("ignore", module=module)

def main(argv):
    #options
    t_key = 'timepoint'
    output_mode = 'full_test'
    perturbation_training = False
    ko_output_genes = "none"

    def future_pred(estimator, adata, T_start, n_timepoints):
        Xs_copy = estimator.Xs[adata].copy()[T_start]
        Xs_pred = [Xs_copy]
        A = estimator.A * estimator.Ms[adata]
        b = estimator.b * estimator.Ms[adata][0, :]
        t = 1 /estimator.T
        P = torch.linalg.matrix_exp(t*A)
        for n in range(n_timepoints):
            future_pred = torch.relu(((Xs_pred[n] / estimator.std) @ P + t * b) * estimator.std)
            Xs_pred.append(future_pred)
        return Xs_pred
    
    inputfile = ''
    outputfolder = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:u:p:k:", ["input=", "outputfolder=", "output_mode=", "perturbation_training=", "ko_genes="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = '{}'.format(arg)
        if opt in ("-o", "--outputfolder"):
            outputfolder = '{}'.format(arg)
        if opt in ("-u", "--output_mode"):
            output_mode = str(arg)
        if opt in ("-p", "--perturbation_training"):
            perturbation_training = arg.lower() == 'true'
        if opt in ("-k", "--ko_genes"):
            ko_output_genes = str(arg)
    
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    data_path = os.path.join(inputfile)
    if os.path.exists(data_path):
        adata = ad.read_h5ad(data_path)
    else:
        raise FileNotFoundError(
            "There is no data available. Create a subfolder 'Data' in your main folder "
            "and put inside a count table named 'data.h5ad'."
        )

    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    
    path_data_folder = os.path.dirname(inputfile)
    train_tps = np.load(path_data_folder+"/train_tps.npy")
    test_tps = np.load(path_data_folder+"/test_tps.npy")

    adata_train = adata[adata.obs[t_key].isin(train_tps[:, 1]),:]
    if perturbation_training == False and "dataset_id" in adata_train.obs:
        if len(np.unique(adata_train.obs["dataset_id"])) > 1:
            adata_train = adata_train[adata_train.obs["dataset_id"] == "WT"].copy()

    # Apply cell subsampling if subsample_train_ids.npy is present
    subsample_ids_path = os.path.join(path_data_folder, "subsample_train_ids.npy")
    if os.path.exists(subsample_ids_path):
        subsample_ids = np.load(subsample_ids_path, allow_pickle=True)
        adata_train = adata_train[adata_train.obs_names.isin(subsample_ids)].copy()

    start_runtime = time.time()
    print("Starting reference fitting model training")

    options = {
        "lr" : 0.05, 
        "reg_sinkhorn" : 0.1,
        "reg_A" : 1, 
        "reg_A_elastic" : 0.5, 
        "iter" : 1000, #1000
        "ot_coupling" : True,
        "optimizer" : torch.optim.Adam,
        "n_pca_components" : -1
    }

    if perturbation_training:
        print("Training with perturbation data")
        # Keep WT first, then KO groups in alphabetical order for stable estimator inputs.
        dataset_groups = sorted(np.unique(adata_train.obs['dataset_id'].astype(str)), key=lambda x: (x != 'WT', x))
        adatas = [adata_train[adata_train.obs['dataset_id'].astype(str)==i].copy() for i in dataset_groups]
        kos = [None if g == 'WT' else g for g in dataset_groups]
        estim_alt_wt = rf.Estimator(adatas, kos = kos,
                lr = options["lr"],
                reg_sinkhorn = options["reg_sinkhorn"], 
                reg_A = options["reg_A"], 
                reg_A_elastic = options["reg_A_elastic"], 
                iter = options["iter"], 
                ot_coupling = options["ot_coupling"],
                optimizer = options["optimizer"],
                norm = False,
                pca = "common",
                t_key = t_key)
    else:
        estim_alt_wt = rf.Estimator([adata_train], kos = [None],
                lr = options["lr"],
                reg_sinkhorn = options["reg_sinkhorn"], 
                reg_A = options["reg_A"], 
                reg_A_elastic = options["reg_A_elastic"], 
                iter = options["iter"], 
                ot_coupling = options["ot_coupling"],
                optimizer = options["optimizer"],
                norm = False,
                t_key = t_key)
    
    estim_alt_wt.fit(print_iter=10, alg = "alternating", update_couplings_iter=250)

    #t = 1/estim_alt_wt.T
    #P = torch.linalg.matrix_exp(t*estim_alt_wt.A)
    #grn = P.cpu().numpy()
    grn = estim_alt_wt.A.cpu().numpy()
    save_path_prefix = os.path.splitext(os.path.basename(inputfile))[0]
    np.save(outputfolder+save_path_prefix+'_GRN', grn)

    ko_targets = []
    if ko_output_genes.lower() == "all":
        ko_targets = [gene for gene in adata_train.var_names if gene != "Stimulus"]
    elif ko_output_genes.lower() not in ["none", ""]:
        ko_targets = [gene.strip() for gene in ko_output_genes.split(",") if gene.strip()]
    invalid_kos = [gene for gene in ko_targets if gene not in adata_train.var_names]
    if invalid_kos:
        raise ValueError(f"Unknown KO gene(s): {invalid_kos}")

    if output_mode != 'no_traj':
        # ref_fit uses squeezed time indices (0-based, 1 unit between samples)
        if output_mode == 'full_test':
            # Simulate from the last training timepoint forward to test timepoints
            T_start_idx = len(train_tps[:, 1]) - 1
            n_steps = len(test_tps[:, 1])
            sim_time_values = np.concatenate([[train_tps[:, 1][-1]], test_tps[:, 1]])
        elif output_mode == 'full_train':
            # Simulate all training timepoints starting from the first timepoint
            T_start_idx = 0
            n_steps = len(train_tps[:, 1]) - 1
            sim_time_values = train_tps[:, 1]
        else:  # full_full
            # Simulate from first timepoint through all train and test timepoints
            T_start_idx = 0
            n_steps = len(train_tps[:, 1]) + len(test_tps[:, 1]) - 1
            sim_time_values = np.unique(np.concatenate([train_tps[:, 1], test_tps[:, 1]]))

        # When subsampling is active, replace the estimator's internal data
        # with cells that are not in subsample_ids for all timepoints
        is_subsampling = os.path.exists(subsample_ids_path)
        if is_subsampling:
            subsample_ids = np.load(subsample_ids_path, allow_pickle=True)
            # Replace Xs with cells not in subsample_ids for all timepoints
            for tp_idx, tp_val in enumerate(train_tps[:, 1]):
                tp_cells = adata[adata.obs[t_key] == tp_val]
                not_in_subsample_mask = ~tp_cells.obs_names.isin(subsample_ids)
                X_not_subsampled = tp_cells[not_in_subsample_mask].X
                if scipy.sparse.issparse(X_not_subsampled):
                    X_not_subsampled = X_not_subsampled.toarray()
                estim_alt_wt.Xs[0][tp_idx] = torch.tensor(X_not_subsampled, dtype=torch.float64)

        def simulate_and_save(estimator, output_suffix, dataset_id):
            n_cells_at_T_start = int(estimator.Xs[0][T_start_idx].shape[0])
            pred_data = future_pred(estimator, 0, T_start_idx, n_steps)
            pred_array = np.array([pred_data[i].cpu().numpy() for i in range(len(pred_data))])
            pred_array = np.reshape(pred_array, newshape=(len(pred_data) * n_cells_at_T_start, adata_train.shape[1]))
            list_times_pred = []
            for i in range(len(pred_data)):
                list_times_pred += [sim_time_values[i]] * n_cells_at_T_start
            pred_adata = ad.AnnData(pred_array)
            pred_adata.obs['time'] = list_times_pred
            pred_adata.obs[t_key] = list_times_pred
            pred_adata.obs['dataset_id'] = dataset_id
            pred_adata.var_names = adata_train.var_names

            if output_mode == 'full_test':
                pred_adata = pred_adata[pred_adata.obs[t_key] != sim_time_values[0]].copy()
                pred_train_data = ad.concat((pred_adata, adata_train))
            else:
                pred_train_data = pred_adata

            pred_train_data.uns['runtime'] = time.time() - start_runtime
            if adata.uns['simulation'] != None:
                pred_train_data = pred_train_data[:, 1:].copy()
            pred_train_data.write_h5ad(outputfolder+save_path_prefix+output_suffix+'_adata.h5ad')

        simulate_and_save(estim_alt_wt, '', 'WT')

    if ko_targets:
        base_A = estim_alt_wt.A.clone()
        for ko_gene in ko_targets:
            ko_idx = list(adata_train.var_names).index(ko_gene)
            estim_alt_wt.A = base_A.clone()
            # KO simulation: zero all incoming edges toward ko_gene.
            estim_alt_wt.A[:, ko_idx] = 0
            np.save(outputfolder+save_path_prefix+f'_ko_{ko_gene}_GRN', estim_alt_wt.A.cpu().numpy())
            if output_mode != 'no_traj':
                simulate_and_save(estim_alt_wt, f'_ko_{ko_gene}', ko_gene)
        estim_alt_wt.A = base_A
    

if __name__ == "__main__":
   main(sys.argv[1:])
