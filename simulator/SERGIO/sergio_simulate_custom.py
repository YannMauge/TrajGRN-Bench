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


def _compute_grn_layers(grn, n_genes: int):
    """
    Compute topological layers of the GRN by BFS from master regulators.

    Returns a list of lists, where layers[0] contains the master regulators
    (no incoming edges, ignoring self-loops) and each subsequent layer
    contains genes whose regulators are all in previous layers.

    Mutual-regulation cycles (e.g. A↔B) are handled by promoting both
    genes to the same layer.
    """
    sergio_edges = grn.sergio_edges

    # Build incoming edge sets (ignoring self-loops for layering)
    incoming = {i: set() for i in range(n_genes)}
    for regulator, target, _strength in sergio_edges:
        if regulator != target:
            incoming[target].add(regulator)

    # Master regulators: genes with no incoming edges
    masters = sorted(i for i in range(n_genes) if not incoming[i])

    layers = []
    visited = set()
    frontier = set(masters)

    while frontier:
        layers.append(sorted(frontier))
        visited.update(frontier)

        # Find next set of genes: a gene is ready if all its regulators
        # are either already visited OR are mutually regulated (part of
        # a cycle where each depends on the other).
        remaining = set(range(n_genes)) - visited
        next_frontier = set()
        for target in remaining:
            regs = incoming[target]
            if not regs:
                next_frontier.add(target)
                continue
            unvisited_regs = regs - visited
            if not unvisited_regs:
                next_frontier.add(target)
            else:
                # Allow if all unvisited regulators mutually depend on target
                all_mutual = all(
                    target in incoming.get(r, set()) for r in unvisited_regs
                )
                if all_mutual:
                    next_frontier.add(target)

        frontier = next_frontier

        # Safety: if stuck but genes remain, promote all remaining at once
        if not frontier and remaining:
            layers.append(sorted(remaining))
            break

    return layers


def _build_full_graph_inputs(grn, active_genes, inactive_genes, n_time_bins, total_time=100):
    """
    Build SERGIO inputs for the FULL set of n_genes.

    Active genes keep their normal regulatory edges (but only edges where
    both regulator AND target are active). Inactive genes are added as
    master regulators with near-zero production rates so they appear in
    all timepoints with very low expression, preventing all-zero columns
    that break downstream methods.
    """
    active_set = set(active_genes)

    # Filter edges: both regulator and target must be active
    incoming = {g: [] for g in active_genes}
    for regulator, target, strength in grn.sergio_edges:
        if regulator in active_set and target in active_set:
            incoming[target].append((regulator, strength))

    # Build targets_rows for active target genes
    targets_rows = []
    for target in active_genes:
        regs = incoming[target]
        if not regs:
            continue  # master regulator within active set
        row = [float(target), float(len(regs))]
        for reg, _strength in regs:
            row.append(float(reg))
        for _reg, strength in regs:
            row.append(float(strength))
        row += [2.0] * len(regs)  # shared_coop_state
        targets_rows.append(row)

    # Master regulators: active genes with no incoming edges from actives,
    # PLUS all inactive genes (with near-zero production rates)
    step = total_time // n_time_bins
    regulator_rate = [float(i * step) for i in range(n_time_bins)]

    active_masters = sorted(g for g in active_genes if not incoming[g])
    regs_rows = [[float(m)] + regulator_rate for m in active_masters]

    # Inactive genes: master regulators with very low rate (0.01 per bin)
    low_rate = [0.01] * n_time_bins
    for g in inactive_genes:
        regs_rows.append([float(g)] + low_rate)

    return targets_rows, regs_rows


def _simulate_sergio_matrix(
    n_cells: int, n_genes: int, n_bins: int, run_seed: int
) -> np.ndarray:
    """
    Simulate gene expression with a cascade effect: each successive group of
    time bins corresponds to a deeper layer of the GRN being active.

    For each GRN layer, a separate SERGIO simulation is run with a subgraph
    containing only genes up to that layer. Inactive genes (deeper layers)
    are simulated as master regulators with near-zero production rates so
    they appear at all timepoints with very low expression, avoiding
    all-zero columns in training data.
    """
    grn = get_benchmark_grn(n_genes)
    layers = _compute_grn_layers(grn, n_genes)
    n_layers = len(layers)
    cells_per_bin = max(1, n_cells // n_bins)
    total_cells = cells_per_bin * n_bins

    sergio_class = _load_sergio_class()
    rng = np.random.default_rng(run_seed)

    # Build gene-to-layer mapping
    gene_to_layer = {}
    for layer_idx, layer_genes in enumerate(layers):
        for g in layer_genes:
            gene_to_layer[g] = layer_idx

    # Determine bins_per_layer: spread layers across all bins, ensuring
    # each layer gets at least 1 bin.
    bins_per_layer = max(1, n_bins // n_layers)

    # Initialize full expression matrix
    expr_full = np.zeros((n_genes, total_cells), dtype=np.float32)

    for layer_idx in range(n_layers):
        # Active genes: all genes in layers 0..layer_idx
        active_genes = sorted(
            g for g in range(n_genes) if gene_to_layer[g] <= layer_idx
        )
        inactive_genes = sorted(set(range(n_genes)) - set(active_genes))

        # Build subgraph inputs with ALL n_genes present:
        # - active genes keep their normal regulatory structure
        # - inactive genes become master regulators with rate 0
        targets_rows, regs_rows = _build_full_graph_inputs(
            grn, active_genes, inactive_genes, bins_per_layer
        )

        np.random.seed(run_seed + layer_idx)
        simulator = sergio_class(
            number_genes=n_genes,
            number_bins=bins_per_layer,
            number_sc=cells_per_bin,
            noise_params=1.0,
            decays=0.8,
            sampling_state=12,
            noise_type="dpd",
        )
        simulator.build_graph_from_rows(
            targets_rows=targets_rows,
            regs_rows=regs_rows,
            shared_coop_state=2,
        )
        simulator.simulate()
        expr_layer = simulator.getExpressions()
        # expr_layer shape: (bins_per_layer, n_genes, cells_per_bin)
        expr_flat = np.concatenate(expr_layer, axis=1)

        # Place into full expression matrix
        start_col = layer_idx * bins_per_layer * cells_per_bin
        end_col = start_col + bins_per_layer * cells_per_bin
        expr_full[:, start_col:end_col] = expr_flat[:, :end_col - start_col]

    expr_full = np.clip(expr_full[:, :total_cells], 0.0, None)
    return rng.poisson(expr_full).astype(int)


def _write_shared_files(output_folder: Path, n_genes: int, degradation_rates: tuple[float, float] = (0.75, 0.15)):
    gene_names = [f"Gene{i}" for i in range(n_genes)]
    write_common_shared_files(output_folder, gene_names, degradation_rates)


def _write_true_network(output_folder: Path, run_index: int, n_genes: int):
    true_dir = output_folder.parent / "True"
    true_dir.mkdir(parents=True, exist_ok=True)
    np.save(true_dir / f"inter_{run_index}", get_benchmark_grn(n_genes).build_true_interaction_matrix())


def _write_data_file(output_folder: Path, run_index: int, counts: np.ndarray, time_vector: np.ndarray):
    n_genes, n_cells = counts.shape
    timepoints = time_vector[:n_cells]

    data = np.zeros((n_cells + 1, n_genes + 2), dtype=int)
    data[0, 1:] = np.arange(n_genes + 1)
    data[1:, 0] = timepoints
    data[1:, 1] = 100 * (timepoints > 0)
    data[1:, 2:] = counts.T
    np.savetxt(output_folder / f"data_{run_index}.txt", data.T, fmt="%d", delimiter="\t")


def main(argv):
    cfg = parse_common_cli(argv, include_ko=False, default_n_cells=1000, default_n_genes=8)
    cfg.output_folder.mkdir(parents=True, exist_ok=True)
    time = cfg.build_time_vector(total_time=100)
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
        _write_data_file(cfg.output_folder, run_index, counts, time)


if __name__ == "__main__":
    main(sys.argv[1:])
