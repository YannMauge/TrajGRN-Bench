"""Pearson co-expression baseline: grn[i, j] = Pearson correlation of genes
i and j over the WT training cells (log-normalised; symmetric, so it carries
no direction). The stimulus column is kept as a variable (log1p of its raw
value) but not counted in the library size."""
import getopt
import os
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import load_training_data, log_normalize_genes, save_grn  # noqa: E402

warnings.filterwarnings("ignore")


def main(argv):
    inputfile, outputfolder = "", ""
    opts, _ = getopt.getopt(argv, "hi:o:", ["input=", "outputfolder="])
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-o", "--outputfolder"):
            outputfolder = arg
    if not os.path.exists(inputfile):
        raise FileNotFoundError(f"Error reading h5ad file: {inputfile}")

    adata, _, _ = load_training_data(inputfile)
    x = log_normalize_genes(adata.X, adata.var_names)
    with np.errstate(invalid="ignore", divide="ignore"):
        pred = np.corrcoef(x, rowvar=False)
    pred = np.nan_to_num(pred, nan=0.0)  # constant variables have no correlation
    path = save_grn(outputfolder, inputfile, pred, adata.var_names)
    print(f"Pearson GRN saved to {path}")


if __name__ == "__main__":
    main(sys.argv[1:])
