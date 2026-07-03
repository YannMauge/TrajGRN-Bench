import csv
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import numpy as np
from simulator.common_config import (
    get_benchmark_grn,
    parse_common_cli,
    write_shared_files as write_common_shared_files,
)


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_SERGIO_PATH = SCRIPT_PATH.parent / "SERGIO"


def _load_sergio_class():
    if DEFAULT_SERGIO_PATH.exists():
        parent = str(DEFAULT_SERGIO_PATH.parent)
        if parent not in sys.path:
            sys.path.insert(0, parent)

    sergio_path = os.environ.get("SERGIO_PATH")
    if sergio_path and sergio_path not in sys.path:
        sys.path.insert(0, sergio_path)

    if not hasattr(np, "int"):
        np.int = int  # type: ignore[attr-defined]
    if not hasattr(np, "float"):
        np.float = float  # type: ignore[attr-defined]

    try:
        from SERGIO.sergio import sergio as sergio_class
    except Exception as exc:  # pragma: no cover - defensive import guard
        raise ImportError(
            "SERGIO could not be imported. Clone https://github.com/PayamDiba/SERGIO "
            "and set SERGIO_PATH to the cloned repository root, or place it at "
            "simulator/SERGIO."
        ) from exc
    return sergio_class


def _write_csv_rows(path: Path, rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


def _simulate_sergio_matrix(n_cells: int, n_genes: int, n_bins: int, run_seed: int) -> np.ndarray:
    grn = get_benchmark_grn(n_genes)
    cells_per_bin = max(1, n_cells // n_bins)
    total_cells = cells_per_bin * n_bins

    targets, regs = grn.build_sergio_inputs(n_bins)

    output_dir = SCRIPT_PATH.parent / "tmp_sergio_inputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    targets_file = output_dir / "targets.csv"
    regs_file = output_dir / "regs.csv"
    _write_csv_rows(targets_file, targets)
    _write_csv_rows(regs_file, regs)

    sergio_class = _load_sergio_class()
    np.random.seed(run_seed)
    simulator = sergio_class(
        number_genes=n_genes,
        number_bins=n_bins,
        number_sc=cells_per_bin,
        noise_params=1.0,
        decays=0.8,
        sampling_state=12,
        noise_type="dpd",
    )
    simulator.build_graph(
        input_file_taregts=str(targets_file),
        input_file_regs=str(regs_file),
        shared_coop_state=2,
    )
    simulator.simulate()
    expr = simulator.getExpressions()
    expr_matrix = np.concatenate(expr, axis=1)
    expr_matrix = np.clip(expr_matrix[:, :total_cells], 0.0, None)

    rng = np.random.default_rng(run_seed)
    return rng.poisson(expr_matrix).astype(int)


def _write_shared_files(output_folder: Path, n_genes: int, degradation_rates: tuple[float, float] = (0.75, 0.15)):
    gene_names = [f"Gene{i}" for i in range(n_genes)]
    write_common_shared_files(output_folder, gene_names, degradation_rates)


def _write_true_network(output_folder: Path, run_index: int, n_genes: int):
    true_dir = output_folder.parent / "True"
    true_dir.mkdir(parents=True, exist_ok=True)
    np.save(true_dir / f"inter_{run_index}", get_benchmark_grn(n_genes).build_true_interaction_matrix())


def _write_data_file(output_folder: Path, run_index: int, counts: np.ndarray, n_time_bins: int):
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
    np.savetxt(output_folder / f"data_{run_index}.txt", data.T, fmt="%d", delimiter="\t")


def main(argv):
    cfg = parse_common_cli(argv, include_ko=False, default_n_cells=1000, default_n_genes=8)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)
    for run_index in range(1, cfg.n_runs + 1):
        print(f"Run {run_index}...")
        counts = _simulate_sergio_matrix(
            n_cells=cfg.n_cells,
            n_genes=cfg.n_genes,
            n_bins=cfg.n_time_bins,
            run_seed=run_index,
        )
        _write_shared_files(cfg.output_folder, n_genes=counts.shape[0], degradation_rates=cfg.degradation_rates)
        _write_true_network(cfg.output_folder, run_index, n_genes=counts.shape[0])
        _write_data_file(cfg.output_folder, run_index, counts, cfg.n_time_bins)


if __name__ == "__main__":
    main(sys.argv[1:])
