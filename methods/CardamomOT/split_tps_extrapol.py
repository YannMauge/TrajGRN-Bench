import numpy as np
import scanpy as sc
import sys
import getopt
import anndata as ad
import os
import pandas as pd
import time
from CardamomOT import extract_degradation_rates

import warnings
for module in ["anndata", "scipy", "torchdiffeq"]:
    warnings.filterwarnings("ignore", module=module)

def main(argv):
    t_key = 'timepoint'
    inputfile = ''
    output_mode = 'full_test'
    perturbation_training = False
    try:
        opts, args = getopt.getopt(argv, "hi:a:u:p:", ["ifile=","adata_file=","output_mode=","perturbation_training="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inputfile = arg
        if opt in ("-a", "--adata_file"):
            adata_file = arg
        if opt in ("-u", "--output_mode"):
            output_mode = arg
        if opt in ("-p", "--perturbation_training"):
            perturbation_training = arg.lower() == 'true'
    
    p = '{}/'.format(inputfile)  # Name of the file where the data is stored.

    data_path = os.path.join(p, 'Data', 'data_orig.h5ad')
    if os.path.exists(data_path):
        adata = ad.read_h5ad(data_path)
    else:
        raise FileNotFoundError(
            "There is no data available. Create a subfolder 'Data' in your main folder "
            "and put inside a count table named 'data_orig.h5ad'."
        )

    path_data_folder = os.path.dirname(adata_file)
    train_tps = np.load(path_data_folder+"/train_tps.npy")[:,1]
    test_tps = np.load(path_data_folder+"/test_tps.npy")[:,1] # in hours
    start_runtime = time.time()

    if perturbation_training == False:
            if len(np.unique(adata.obs["dataset_id"])) > 1:
                adata = adata[adata.obs["dataset_id"] == "WT"].copy()
    adata = adata[adata.obs[t_key].isin(train_tps)].copy()
    print(f"[split_tps_extrapol] Loaded data with {adata.n_obs} cells and {adata.n_vars} genes.")

    if adata.uns['simulation'] != None:
    # Remove the 1st gene (Stimulus)
        adata = adata[:, 1:].copy()

    # Apply cell subsampling if subsample_train_ids.npy is present
    subsample_ids_path = os.path.join(path_data_folder, "subsample_train_ids.npy")
    adata_test = None
    if os.path.exists(subsample_ids_path):
        subsample_ids = np.load(subsample_ids_path, allow_pickle=True)
        adata_test = adata[~adata.obs_names.isin(subsample_ids)].copy()  # cells NOT in subsample_ids
        adata = adata[adata.obs_names.isin(subsample_ids)].copy()

    csv_path = path_data_folder+ "/halflife/table_halflife_mamalian.csv"
    if 'd0' not in adata.var.columns:
        df = pd.read_csv(csv_path, sep=',')
        deg = extract_degradation_rates(df, adata.var_names)
        adata.var['d0'] = deg[0]
        adata.var['d1'] = deg[1]
        
        # Apply the same degradation rates to test data if it exists
        if adata_test is not None:
            adata_test.var['d0'] = deg[0]
            adata_test.var['d1'] = deg[1]
    
    adata.uns['runtime'] = start_runtime
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(p+'Data/data_full.h5ad')
    
    # Save test data if it exists
    if adata_test is not None:
        adata_test.uns['runtime'] = start_runtime
        adata_test.write_h5ad(p+'Data/data_test.h5ad')

    # Determine which timepoints to simulate based on output_mode
    if output_mode == 'full_test':
        # Simulate test timepoints (current default behaviour)
        times_to_simulate = test_tps
    elif output_mode == 'full_train':
        # Simulate all training timepoints
        times_to_simulate = train_tps
    elif output_mode == 'full_full':
        # Simulate all train and test timepoints
        times_to_simulate = np.unique(np.concatenate([train_tps, test_tps]))
    else:  # no_traj
        # No simulation requested; write an empty file so downstream scripts don't fail
        times_to_simulate = np.array([], dtype=int)

    np.savetxt(p+'Data/times_to_simulate.txt', times_to_simulate, fmt='%d')

if __name__ == "__main__":
   main(sys.argv[1:])
