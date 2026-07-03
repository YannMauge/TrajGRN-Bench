import numpy as np
import scipy
import scanpy as sc
import anndata as ad
import pandas as pd
from GENIE3 import *
import sys, getopt
import os

import warnings
for module in ["anndata", "scipy", "torchdiffeq", "scanpy"]:
    warnings.filterwarnings("ignore", module=module)

def main(argv):
    t_key = 'timepoint'
    lognorm = 1
    inputfile = ''
    outputfolder = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["input=", "outputfolder="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = '{}'.format(arg)
        if opt in ("-o", "--outputfolder"):
            outputfolder = '{}'.format(arg)

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

    path_data_folder = os.path.dirname(inputfile)
    train_tps = np.load(path_data_folder+"/train_tps.npy")[:,1]
    adata = adata[adata.obs[t_key].isin(train_tps)].copy()

    if lognorm:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    
    pred = GENIE3(adata.X)
    save_path_prefix = os.path.splitext(os.path.basename(inputfile))[0]
    np.save(outputfolder+save_path_prefix+'_GRN',pred)

if __name__ == "__main__":
   main(sys.argv[1:])