"""BoolODE WT + knockout replicates for the benchmark GRN (see boolode_driver.py).

Each KO is re-simulated with the KO gene's production rate and initial value
set to 0 (same seeds as WT). Writes data_<r>_ko_WT.txt and data_<r>_ko_<gene>.txt.

Usage: python boolode_simulate_custom_ko.py -o <output_folder> -a <adata_folder>
       -n <n_runs> -k <all|none|Gene1,Gene2> [--n_cells N] [--n_time_bins B]
       [--seed S] [--sim_time T] [--step H] [--modeltype hill|heaviside]
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from simulator.common_config import parse_common_cli, parse_ko_labels  # noqa: E402
from boolode_driver import EXTRA_OPTIONS, simulate_runs  # noqa: E402


def main(argv):
    cfg = parse_common_cli(argv, include_ko=True, default_n_cells=1000, default_n_genes=8,
                           extra_options=EXTRA_OPTIONS)
    if cfg.output_folder.exists():
        print(f"Output folder {cfg.output_folder} already exists, deleting...")
        shutil.rmtree(cfg.output_folder)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)
    simulate_runs(cfg, ko_labels=parse_ko_labels(cfg.ko_genes, cfg.gene_names), write_ko_files=True)


if __name__ == "__main__":
    main(sys.argv[1:])
