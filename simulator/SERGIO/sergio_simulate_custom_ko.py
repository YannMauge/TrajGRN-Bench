"""SERGIO WT + knockout replicates (KO gene production set to 0 via the input files; see sergio_driver.py).

Usage: python sergio_simulate_custom_ko.py -o <output_folder> -a <adata_folder> -n <n_runs> -k <all|none|Gene1,...> [--n_cells N] [--n_time_bins B] [--seed S]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from simulator.common_config import parse_common_cli, parse_ko_labels  # noqa: E402
from sergio_driver import EXTRA_OPTIONS, simulate_runs  # noqa: E402


def main(argv):
    cfg = parse_common_cli(argv, include_ko=True, default_n_cells=1000, default_n_genes=8,
                           extra_options=EXTRA_OPTIONS)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)
    simulate_runs(cfg, ko_labels=parse_ko_labels(cfg.ko_genes, cfg.gene_names), write_ko_files=True)


if __name__ == "__main__":
    main(sys.argv[1:])
