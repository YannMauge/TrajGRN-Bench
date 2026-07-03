import numpy as np
import scanpy as sc
import sys
import getopt
import anndata as ad
import os
import pandas as pd
import time
from CardamomOT import extract_degradation_rates

def main(argv):
    t_key = 'timepoint'
    inputfile = ''
    ko_output_genes = 'none'
    try:
        opts, args = getopt.getopt(argv, "hi:a:k:p:", ["ifile=","adata_file=","ko_genes=","perturbation_training="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inputfile = arg
        if opt in ("-a", "--adata_file"):
            adata_file = arg
        if opt in ("-k", "--ko_genes"):
            ko_output_genes = str(arg)
        if opt in ("-p", "--perturbation_training"):
            perturbation_training = arg.lower() == 'true'
    
    p = '{}/'.format(inputfile)  # Name of the file where are the data

    data_path = os.path.join(p, 'Data', 'data_full.h5ad')
    if os.path.exists(data_path):
        adata = ad.read_h5ad(data_path)
    else:
        raise FileNotFoundError(
            "There is no data available."
        )

    path_data_folder = os.path.dirname(adata_file)
    train_tps = np.load(path_data_folder+"/train_tps.npy")
    test_tps = np.load(path_data_folder+"/test_tps.npy")

    adata_train = adata[adata.obs[t_key].isin(train_tps[:, 1])].copy()
    adata_train.obs['time'] = adata_train.obs[t_key]
    start_runtime = adata.uns['runtime']

    adata_sim = ad.read_h5ad(p+'CardamomOT/adata_sim_stim1.0_prior1.0.h5ad')
    adata_sim.obs[t_key] = adata_sim.obs['time']
    adata_sim = adata_sim[adata_sim.obs[t_key].isin(test_tps[:, 1])].copy()
    if perturbation_training:
        if "dataset_id" in adata_sim.obs.columns:
            if len(np.unique(adata_sim.obs["dataset_id"])) > 1:
                adata_sim = adata_sim[adata_sim.obs["dataset_id"] == "WT"].copy()
        else:
            # Resample adata_sim to match number of cells in adata_train per timepoint
            adata_sim_resampled = []
            n_cells_train = len(adata_train[adata_train.obs[t_key] == train_tps[0, 1]])
            for timepoint in adata_sim.obs[t_key].unique():
                adata_sim_tp = adata_sim[adata_sim.obs[t_key] == timepoint]
                if len(adata_sim_tp) > 0:
                    # Sample n_cells_train cells from adata_sim_tp (without replacement)
                    indices = np.random.choice(len(adata_sim_tp), min(n_cells_train, len(adata_sim_tp)), replace=False)
                    adata_sim_resampled.append(adata_sim_tp[indices].copy())
            if adata_sim_resampled:
                adata_sim = ad.concat(adata_sim_resampled)

    sim_train_data = ad.concat((adata_sim,adata_train))
    sim_train_data.uns['runtime'] = time.time() - start_runtime

    ad.settings.allow_write_nullable_strings = True
    sim_train_data.write_h5ad(p+'CardamomOT/adata_sim_final.h5ad')

if __name__ == "__main__":
   main(sys.argv[1:]) 