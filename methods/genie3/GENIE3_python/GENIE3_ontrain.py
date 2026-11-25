"""GENIE3 baseline (Huynh-Thu et al. 2010, vendored Python implementation).

Runs on the WT training cells (log-normalised; the stimulus column is a
variable but not counted in the library size). GENIE3 returns
VIM[i, j] = importance of gene i for target j (regulator × target), which is
the benchmark orientation. GENIE3 has no seed argument: numpy's global RNG is
seeded before the run.
"""
import getopt
import os
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from GENIE3 import GENIE3  # noqa: E402
from _common import (  # noqa: E402
    load_training_data, log_normalize_genes, run_on_variable_columns, save_grn, seed_everything,
)

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

    seed_everything(0)
    adata, _, _ = load_training_data(inputfile)
    x = log_normalize_genes(adata.X, adata.var_names)
    pred = run_on_variable_columns(GENIE3, x)  # constant variables -> 0
    path = save_grn(outputfolder, inputfile, pred, adata.var_names)
    print(f"GENIE3 GRN saved to {path}")


if __name__ == "__main__":
    main(sys.argv[1:])
