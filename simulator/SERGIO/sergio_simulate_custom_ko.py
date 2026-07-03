import sys
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from simulator.common_config import parse_common_cli, parse_ko_labels
from sergio_simulate_custom import _simulate_sergio_matrix, _write_shared_files, _write_true_network


def _write_ko_data_file(
    output_folder: Path,
    run_index: int,
    dataset_id: str,
    counts: np.ndarray,
    n_time_bins: int,
):
    n_genes, n_cells = counts.shape
    edges = np.linspace(0, n_cells, n_time_bins + 1, dtype=int)
    timepoints = np.zeros(n_cells, dtype=int)
    for idx in range(n_time_bins):
        timepoints[edges[idx] : edges[idx + 1]] = idx

    data = np.zeros((n_cells + 1, n_genes + 2), dtype=int)
    data[0, 1:] = np.arange(n_genes + 1)
    data[1:, 0] = timepoints
    data[1:, 1] = 100 * (timepoints > 0)
    data[1:, 2:] = counts.T
    np.savetxt(output_folder / f"data_{run_index}_ko_{dataset_id}.txt", data.T, fmt="%d", delimiter="\t")


def main(argv):
    cfg = parse_common_cli(argv, include_ko=True, default_n_cells=1000, default_n_genes=8)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)
    ko_labels = parse_ko_labels(cfg.ko_genes, cfg.gene_names)

    for run_index in range(1, cfg.n_runs + 1):
        print(f"Run {run_index}...")
        wt_counts = _simulate_sergio_matrix(
            n_cells=cfg.n_cells,
            n_genes=cfg.n_genes,
            n_bins=cfg.n_time_bins,
            run_seed=run_index,
        )
        _write_shared_files(
            cfg.output_folder,
            n_genes=wt_counts.shape[0],
            degradation_rates=cfg.degradation_rates,
        )
        _write_true_network(cfg.output_folder, run_index, n_genes=wt_counts.shape[0])
        _write_ko_data_file(cfg.output_folder, run_index, "WT", wt_counts, cfg.n_time_bins)

        for ko_label in ko_labels:
            ko_counts = wt_counts.copy()
            ko_counts[cfg.gene_names.index(ko_label), :] = 0
            _write_ko_data_file(cfg.output_folder, run_index, ko_label, ko_counts, cfg.n_time_bins)


if __name__ == "__main__":
    main(sys.argv[1:])
