import sys
from pathlib import Path
import shutil
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import numpy as np
import pandas as pd
from simulator.common_config import (
    get_benchmark_grn,
    parse_common_cli,
    write_boolode_model_files,
    write_shared_files,
)


SCRIPT_PATH = Path(__file__).resolve()
BOOL_ODE_ROOT = SCRIPT_PATH.parent
if str(BOOL_ODE_ROOT) not in sys.path:
    sys.path.insert(0, str(BOOL_ODE_ROOT))

import BoolODE as bo  # noqa: E402


def _run_boolode(output_folder: Path, run_index: int, n_cells: int, n_genes: int, max_time: float, n_ratio_time: int) -> Path:
    model_dir = output_folder / "boolode_model"
    run_name = f"benchmark_run_{run_index}"
    output_dir = output_folder / "boolode_outputs"
    model_definition, model_initial_conditions = write_boolode_model_files(model_dir, n_genes)

    global_settings = bo.GlobalSettings(
        model_dir=str(model_dir),
        output_dir=str(output_dir),
        do_simulations=True,
        do_post_processing=False,
        modeltype="heaviside",
    )
    job_settings = bo.JobSettings(
        [
            {
                "name": run_name,
                "model_definition": model_definition,
                "model_initial_conditions": model_initial_conditions,
                "simulation_time": float(max_time),
                "num_cells": int(n_cells),
                "do_parallel": True,
                "sample_cells": True,
                "integration_step_size": 1/n_ratio_time,
            }
        ]
    )
    post_settings = bo.PostProcSettings(None, None, None, None, None)
    bo.BoolODE(job_settings, global_settings, post_settings).execute_jobs()
    return output_dir / run_name


def _load_boolode_outputs(run_dir: Path):
    exp_df = pd.read_csv(run_dir / "ExpressionData.csv", index_col=0)
    pseudo_df = pd.read_csv(run_dir / "PseudoTime.csv")

    if "Cell ID" not in pseudo_df.columns or "Time" not in pseudo_df.columns:
        raise ValueError("PseudoTime.csv must contain 'Cell ID' and 'Time' columns")

    time_by_cell = {
        str(cell_id): float(time)
        for cell_id, time in zip(pseudo_df["Cell ID"], pseudo_df["Time"])
    }
    ordered_cells = [str(c) for c in exp_df.columns]
    if not all(c in time_by_cell for c in ordered_cells):
        missing = [c for c in ordered_cells if c not in time_by_cell][:5]
        raise ValueError(f"Missing time annotations for cells: {missing}")

    raw_times = np.array([time_by_cell[c] for c in ordered_cells], dtype=float)
    genes = [str(g) for g in exp_df.index.tolist()]

    expr_values = np.clip(exp_df.T.to_numpy(dtype=float), 0.0, None)
    order = np.argsort(raw_times, kind="stable")
    return genes, raw_times[order], expr_values[order, :]


def _discretize_timepoints(raw_times: np.ndarray, n_time_bins: int, n_ratio_time: int) -> np.ndarray:
    n_cells = raw_times.shape[0]
    timepoints = np.zeros(n_cells, dtype=int)
    for i in range(n_cells):
        timepoints[i] = raw_times[i] // n_ratio_time
    return timepoints



def _write_data_file(output_folder: Path, run_index: int, times, counts, stimulus_value: int):
    n_cells, n_genes = counts.shape
    data = np.zeros((n_cells + 1, n_genes + 2), dtype=int)
    data[0, 1:] = np.arange(n_genes + 1)
    data[1:, 0] = times
    data[1:, 1] = int(stimulus_value) * (times > 0)
    data[1:, 2:] = counts

    output_file = output_folder / f"data_{run_index}.txt"
    np.savetxt(output_file, data.T, fmt="%d", delimiter="\t")


def _write_true_network(output_folder: Path, run_index: int, n_genes: int):
    true_dir = output_folder.parent / "True"
    true_dir.mkdir(parents=True, exist_ok=True)
    np.save(true_dir / f"inter_{run_index}", get_benchmark_grn(n_genes).build_true_interaction_matrix())


def main(argv):
    n_ratio_time = 0.2
    cfg = parse_common_cli(argv, include_ko=False)
    if cfg.output_folder.exists():
        print(f"Output folder {cfg.output_folder} already exists, deleting...")
        shutil.rmtree(cfg.output_folder)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)

    for run_index in range(1, cfg.n_runs + 1):
        print(f"Run {run_index}...")
        run_dir = _run_boolode(
            cfg.output_folder, run_index, n_cells=cfg.n_cells, n_genes=cfg.n_genes, max_time=cfg.n_time_bins,
        n_ratio_time=n_ratio_time)
        genes, raw_times, expr_values = _load_boolode_outputs(run_dir)

        # Rescale expression values so the global maximum is 100 (before discretization)
        if expr_values.size > 0:
            max_expr = float(np.max(expr_values))
            if max_expr > 0.0:
                expr_values = expr_values * (100.0 / max_expr)

        times = _discretize_timepoints(raw_times, cfg.n_time_bins, n_ratio_time=n_ratio_time)
        write_shared_files(cfg.output_folder, cfg.gene_names, cfg.degradation_rates)
        _write_true_network(cfg.output_folder, run_index, cfg.n_genes)

        rng = np.random.default_rng(run_index)
        counts = rng.poisson(expr_values).astype(int)
        _write_data_file(cfg.output_folder, run_index, times, counts, cfg.stimulus_value)


if __name__ == "__main__":
    main(sys.argv[1:])
