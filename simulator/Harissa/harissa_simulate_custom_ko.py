import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from harissa import NetworkModel
from simulator.common_config import (
    get_benchmark_grn,
    parse_common_cli,
    parse_ko_labels,
    write_shared_files,
)

def _build_base_model(gene_count: int) -> NetworkModel:
    grn = get_benchmark_grn(gene_count)
    model = NetworkModel(gene_count)
    model.d[0] = 0.25
    model.d[1] = 0.05

    model.basal[1:] = grn.basal
    for regulator, target, strength in grn.harissa_edges:
        model.inter[regulator, target] = strength
    return model


def _simulate_dataset(model: NetworkModel, timepoints: np.ndarray, n_cells: int, n_genes: int) -> np.ndarray:
    data = np.zeros((n_cells + 1, n_genes + 2), dtype="int")
    data[0][1:] = np.arange(n_genes + 1)
    data[1:, 0] = timepoints
    data[1:, 1] = 100 * (timepoints > 0)
    for cell_idx in range(n_cells):
        sim = model.simulate(timepoints[cell_idx], burnin=5)
        data[cell_idx + 1, 2:] = np.random.poisson(sim.m[-1])
    return data


def main(argv):
    cfg = parse_common_cli(argv, include_ko=True, default_n_cells=1000, default_n_genes=8)
    np.random.seed(0)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)
    time = cfg.build_time_vector(total_time=100)
    ko_labels = parse_ko_labels(cfg.ko_genes, cfg.gene_names)
    print( f"Simulating with KO genes: {ko_labels}...")

    true_dir = cfg.output_folder.parent / "True"
    true_dir.mkdir(parents=True, exist_ok=True)

    base_model = _build_base_model(cfg.n_genes)
    write_shared_files(cfg.output_folder, cfg.gene_names, cfg.degradation_rates)
    np.savetxt(
        cfg.output_folder / "Rates" / "degradation_rates.txt",
        base_model.d.T,
        fmt="%1.3f",
        delimiter="\t",
    )

    for run_idx in range(cfg.n_runs):
        print(f"Run {run_idx + 1}...")
        run_model = _build_base_model(cfg.n_genes)
        inter = 1 * (abs(run_model.inter) > 0)
        np.save(true_dir / f"inter_{run_idx + 1}", inter)

        wt_data = _simulate_dataset(run_model, time, cfg.n_cells, cfg.n_genes)
        np.savetxt(
            cfg.output_folder / f"data_{run_idx + 1}_ko_WT.txt",
            wt_data.T,
            fmt="%d",
            delimiter="\t",
        )

        for ko_label in ko_labels:
            ko_model = _build_base_model(cfg.n_genes)
            ko_column = cfg.gene_names.index(ko_label) + 1
            ko_model.inter[:, ko_column] = 0
            ko_model.basal[ko_column] = -100
            ko_data = _simulate_dataset(ko_model, time, cfg.n_cells, cfg.n_genes)
            np.savetxt(
                cfg.output_folder / f"data_{run_idx + 1}_ko_{ko_label}.txt",
                ko_data.T,
                fmt="%d",
                delimiter="\t",
            )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
