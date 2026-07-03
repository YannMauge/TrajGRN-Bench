import matplotlib.pyplot as plt
import torch
import numpy as np
import scanpy as sc
import anndata as ad
import scipy

from optim.running import constructscNODEModel, scNODETrainWithPreTrain, scNODEPredict
import os
import sys, getopt
import time

import warnings
warnings.filterwarnings("ignore", module="anndata")

def main(argv):  
    #options
    t_key = 'timepoint'
    inputfile = ''
    outputfolder = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:u:", ["input=", "outputfolder=", "output_mode="])
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
    
    if output_mode != 'no_traj':

        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
    
        data_path = os.path.join(inputfile)
        if os.path.exists(data_path):
            adata = ad.read_h5ad(data_path)
        else:
            raise FileNotFoundError(
                "Error reading h5ad file"
            )

        if scipy.sparse.issparse(adata.X):
            adata.X = adata.X.toarray()
        if len(np.unique(adata.obs["dataset_id"])) > 1:
            adata = adata[adata.obs["dataset_id"] == "WT"].copy()

        path_data_folder = os.path.dirname(inputfile)
        train_tps_arr = np.load(path_data_folder+"/train_tps.npy")
        test_tps_arr = np.load(path_data_folder+"/test_tps.npy")
        start_runtime = time.time()

        cell_tps = adata.obs[t_key] + 1
        n_tps = len(np.unique(adata.obs[t_key]))
        n_genes = adata.shape[1]
        cell_types = None
        data = adata.X
    
        # Convert to torch project
        traj_data = [torch.FloatTensor(data[np.where(cell_tps == t)[0]]) for t in np.unique(cell_tps)]
        if cell_types is not None:
                traj_cell_types = [cell_types[np.where(cell_tps == t)[0]] for t in np.unique(cell_tps)]
        
        all_tps = list(range(n_tps))
        train_data = [traj_data[np.where(np.unique(adata.obs[t_key])==t)[0][0]] for t in train_tps_arr[:, 1]]

        # Apply cell subsampling if subsample_train_ids.npy is present
        subsample_ids_path = os.path.join(path_data_folder, "subsample_train_ids.npy")
        if os.path.exists(subsample_ids_path):
            subsample_ids = set(np.load(subsample_ids_path, allow_pickle=True).tolist())
            subsampled_train_data = []
            for idx, t_val in enumerate(train_tps_arr[:, 1]):
                tp_mask = adata.obs[t_key].values == t_val
                tp_cells = adata.obs_names[tp_mask]
                subsample_mask = np.array([c in subsample_ids for c in tp_cells])
                subsampled_train_data.append(train_data[idx][subsample_mask])
            train_data = subsampled_train_data

        test_data = [traj_data[np.where(np.unique(adata.obs[t_key])==t)[0][0]] for t in test_tps_arr[:, 1]]
        tps = torch.FloatTensor(all_tps)
        train_tps_torch = torch.FloatTensor(train_tps_arr[:, 1])
        test_tps_torch = torch.FloatTensor(test_tps_arr[:, 1])
        n_cells = [each.shape[0] for each in traj_data]
        
        # ======================================================
        # Model training
        pretrain_iters = 200
        pretrain_lr = 1e-3
        latent_coeff = 1.0 # regularization coefficient: beta
        epochs = 10 #10
        iters = 100 #100
        batch_size = 32
        lr = 1e-3
        act_name = "relu"
        n_sim_cells = 100
        latent_dim = 50
        drift_latent_size = [50, 50]
        enc_latent_list = None
        dec_latent_list = None
    
        latent_ode_model = constructscNODEModel(
            n_genes, latent_dim=latent_dim, enc_latent_list=enc_latent_list, 
            dec_latent_list=dec_latent_list, drift_latent_size=drift_latent_size,
            latent_enc_act="none", latent_dec_act=act_name, drift_act=act_name,
            ode_method="euler"
        )
        latent_ode_model, loss_list, recon_obs, first_latent_dist, latent_seq = scNODETrainWithPreTrain(
            train_data, train_tps_torch, latent_ode_model, latent_coeff=latent_coeff, epochs=epochs, iters=iters,
            batch_size=batch_size, lr=lr, pretrain_iters=pretrain_iters, pretrain_lr=pretrain_lr
        )
    
        train_times = train_tps_arr[:, 1].tolist()
        save_path_prefix = os.path.splitext(os.path.basename(inputfile))[0]
        adata_train = adata[adata.obs[t_key].isin(train_tps_arr[:, 1]), :].copy()

        # Helper: find the traj_data index for a given timepoint value.
        # traj_data is ordered by np.unique(adata.obs[t_key]), so the index of a
        # value in that sorted array equals the index into traj_data.
        def tp_to_traj_idx(tp_value):
            return int(np.where(np.unique(adata.obs[t_key]) == tp_value)[0][0])

        # Determine which timepoints to simulate based on output_mode
        if output_mode == 'full_test':
            # Simulate test timepoints starting from the last training timepoint
            pred_tps = test_tps_torch
            start_tp_val = train_tps_arr[:, 1][-1]
            start_cells = traj_data[tp_to_traj_idx(start_tp_val)]
        elif output_mode == 'full_train':
            # Simulate all training timepoints starting from the first timepoint
            pred_tps = train_tps_torch
            start_tp_val = train_tps_arr[:, 1][0]
            start_cells = traj_data[tp_to_traj_idx(start_tp_val)]
        else:  # full_full
            # Simulate all timepoints (train + test) starting from the first timepoint
            all_sim_tps = np.unique(np.concatenate([train_tps_arr[:, 1], test_tps_arr[:, 1]]))
            pred_tps = torch.FloatTensor(all_sim_tps)
            start_tp_val = all_sim_tps[0]
            start_cells = traj_data[tp_to_traj_idx(start_tp_val)]

        # When subsampling exists, use cells not in subsample_ids as starting cells
        subsample_ids_path = os.path.join(path_data_folder, "subsample_train_ids.npy")
        if os.path.exists(subsample_ids_path):
            subsample_ids = set(np.load(subsample_ids_path, allow_pickle=True).tolist())
            tp_mask = adata.obs[t_key].values == start_tp_val
            tp_cells = adata.obs_names[tp_mask]
            not_subsampled_mask = np.array([c not in subsample_ids for c in tp_cells])
            start_cells = torch.FloatTensor(adata[tp_mask].X[not_subsampled_mask])

        all_recon_obs = scNODEPredict(latent_ode_model, start_cells, pred_tps, n_cells=n_sim_cells)
        all_recon_obs = np.maximum(all_recon_obs, 0) #Band-aid to avoid lognorm crash in case of bad model
        list_times_pred = []
        for i in pred_tps.cpu().numpy():
            list_times_pred += [i] * n_sim_cells
        pred_adata = ad.AnnData(np.reshape(all_recon_obs, newshape=(n_sim_cells*len(pred_tps), n_genes)))
        pred_adata.obs['time'] = list_times_pred
        pred_adata.obs[t_key] = list_times_pred
        pred_adata.var_names = adata_train.var_names

        if output_mode == 'full_test':
            pred_train_data = ad.concat((pred_adata, adata_train))
        else:
            # full_train / full_full: output is the full simulated trajectory
            pred_train_data = pred_adata

        pred_train_data.uns['runtime'] = time.time() - start_runtime
        if adata.uns['simulation'] != None:
        # Remove the 1st gene (Stimulus)
            pred_train_data = pred_train_data[:, 1:].copy()
        pred_train_data.write_h5ad(outputfolder+save_path_prefix+'_adata.h5ad')
    else:
        print('No pred, did not compute anything')

if __name__ == "__main__":
   main(sys.argv[1:])
    
    
    
    
    
    
    
    
