"""BoolODE WT replicates for the benchmark GRN (see boolode_driver.py).

Usage: python boolode_simulate_custom.py -o <output_folder> -a <adata_folder>
       -n <n_runs> [--n_cells N] [--n_time_bins B] [--seed S]
       [--sim_time T] [--step H] [--modeltype hill|heaviside]
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from simulator.common_config import parse_common_cli  # noqa: E402
from boolode_driver import EXTRA_OPTIONS, simulate_runs  # noqa: E402


def main(argv):
    cfg = parse_common_cli(argv, include_ko=False, default_n_cells=1000, default_n_genes=8,
                           extra_options=EXTRA_OPTIONS)
    if cfg.output_folder.exists():
        print(f"Output folder {cfg.output_folder} already exists, deleting...")
        shutil.rmtree(cfg.output_folder)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)
    simulate_runs(cfg, ko_labels=[], write_ko_files=False)


if __name__ == "__main__":
    main(sys.argv[1:])
