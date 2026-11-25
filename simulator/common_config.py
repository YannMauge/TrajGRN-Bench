from __future__ import annotations

import getopt
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np


@dataclass
class SimulatorCommonConfig:
    output_folder: Path
    adata_folder: Optional[Path]
    n_runs: int = 2
    ko_genes: str = "all"
    n_cells: int = 100
    n_genes: int = 8
    n_time_bins: int = 20
    stimulus_value: int = 100
    degradation_rates: tuple[float, float] = (0.75, 0.15)
    seed: int = 0
    # Simulator-specific options (e.g. BoolODE --sim_time), as strings.
    extra: Dict[str, str] = field(default_factory=dict)

    @property
    def gene_names(self) -> List[str]:
        return [f"Gene{i}" for i in range(self.n_genes)]

    def time_labels(self, total_time: int = 100) -> np.ndarray:
        """Timepoint labels shared by all simulators: n_time_bins evenly
        spaced integers starting at 0 (0, 10, …, 90 for 10 bins)."""
        return np.arange(self.n_time_bins) * (total_time // self.n_time_bins)

    def run_seed(self, run_index: int) -> int:
        """Reproducible, distinct seed for replicate *run_index* (1-based)."""
        return int(self.seed) * 1000 + int(run_index)

    def build_time_vector(self, total_time: int = 100) -> np.ndarray:
        bins = self.time_labels(total_time)
        edges = np.linspace(0, self.n_cells, self.n_time_bins + 1, dtype=int)
        time = np.zeros(self.n_cells, dtype=int)
        for idx, value in enumerate(bins):
            time[edges[idx] : edges[idx + 1]] = value
        return time


@dataclass(frozen=True)
class BenchmarkGRNConfig:
    n_genes: int
    basal: tuple[float, ...]
    harissa_edges: tuple[tuple[int, int, float], ...]

    @property
    def sergio_edges(self) -> tuple[tuple[int, int, float], ...]:
        return tuple(
            (regulator - 1, target - 1, strength)
            for regulator, target, strength in self.harissa_edges
            if regulator > 0
        )

    def build_true_interaction_matrix(self) -> np.ndarray:
        """Signed ground-truth GRN, shape (n_genes+1, n_genes+1).

        ``inter[regulator, target]`` holds the signed interaction strength
        (> 0 activation, < 0 inhibition, 0 no edge); index 0 is the stimulus.
        """
        inter = np.zeros((self.n_genes + 1, self.n_genes + 1), dtype=float)
        for regulator, target, strength in self.harissa_edges:
            inter[regulator, target] = float(strength)
        return inter

    def build_sergio_inputs(self, n_time_bins: int) -> tuple[List[List[float]], List[List[float]]]:
        incoming = {gene_idx: [] for gene_idx in range(self.n_genes)}
        for regulator, target, strength in self.sergio_edges:
            incoming[target].append((regulator, strength))

        targets: List[List[float]] = []
        for target_idx in range(self.n_genes):
            regulators = incoming[target_idx]
            if not regulators:
                continue
            row: List[float] = [float(target_idx), float(len(regulators))]
            row += [float(reg_idx) for reg_idx, _ in regulators]
            row += [float(strength) for _, strength in regulators]
            row += [2.0] * len(regulators)
            targets.append(row)

        master_regulators = [idx for idx in range(self.n_genes) if not incoming[idx]]
        total_time = 100
        step = total_time // n_time_bins
        regulator_rate = [float(i * step) for i in range(n_time_bins)]
        regs = [[float(idx)] + regulator_rate for idx in master_regulators]
        return targets, regs

    def build_boolode_rules(self) -> tuple[tuple[str, str], ...]:
        incoming_positive = {gene_idx: [] for gene_idx in range(1, self.n_genes + 1)}
        incoming_negative = {gene_idx: [] for gene_idx in range(1, self.n_genes + 1)}
        stimulus_targets = set()

        for regulator, target, strength in self.harissa_edges:
            if target < 1 or target > self.n_genes:
                continue
            if regulator == 0:
                if strength > 0:
                    stimulus_targets.add(target)
                continue
            if strength >= 0:
                incoming_positive[target].append(regulator)
            else:
                incoming_negative[target].append(regulator)

        rules: list[tuple[str, str]] = []
        for target in range(1, self.n_genes + 1):
            activators = sorted(set(incoming_positive[target]))
            inhibitors = sorted(set(incoming_negative[target]))

            gene_name = f"g{target}"
            if not activators and not inhibitors:
                if target in stimulus_targets:
                    rule = f"( {gene_name} )"
                else:
                    rule = f"( {gene_name} )"
            elif activators and not inhibitors:
                activator_expr = " or ".join(f"g{idx}" for idx in activators)
                rule = f"( {activator_expr} )"
            elif inhibitors and not activators:
                inhibitor_expr = " or ".join(f"g{idx}" for idx in inhibitors)
                rule = f"not ( {inhibitor_expr} )"
            else:
                activator_expr = " or ".join(f"g{idx}" for idx in activators)
                inhibitor_expr = " or ".join(f"g{idx}" for idx in inhibitors)
                rule = f"(( {activator_expr} ) and not( {inhibitor_expr} ))"
            rules.append((gene_name, rule))

        return tuple(rules)

    def build_boolode_initial_conditions(self) -> tuple[tuple[str, str], ...]:
        initial_genes = sorted(
            {
                f"g{target}"
                for regulator, target, strength in self.harissa_edges
                if regulator == 0 and strength > 0 and 1 <= target <= self.n_genes
            },
            key=lambda name: int(name[1:]),
        )
        if not initial_genes:
            initial_genes = ["g1"]
        return ((str(initial_genes), str([1] * len(initial_genes))),)


_BENCHMARK_GRN_BY_GENES = {
    8: BenchmarkGRNConfig(
        n_genes=8,
        basal=(-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0),
        harissa_edges=(
            (0, 1, 10.0),
            (1, 2, 10.0),
            (1, 3, 10.0),
            (3, 2, -10.0),
            (2, 3, -10.0),
            (2, 2, 5.0),
            (3, 3, 5.0),
            (2, 4, 10.0),
            (3, 5, 10.0),
            (2, 5, -10.0),
            (3, 4, -10.0),
            (4, 7, -10.0),
            (5, 6, -10.0),
            (4, 6, 10.0),
            (5, 7, 10.0),
            (7, 8, 10.0),
            (6, 8, -10.0),
        ),
    ),
}


# SERGIO cannot simulate cycles or self-loops (Dibaeinia & Sinha 2020), so it
# gets its own acyclic network with the same entry point (stimulus -> Gene0, a
# master regulator whose production rate ramps along the differentiation
# chain), the same number of genes and a mix of activations / inhibitions,
# built from motifs SERGIO handles: incoherent FFL (Gene0 -> Gene1 -| Gene2,
# Gene0 -> Gene2), convergent regulation (Gene1 -> Gene3 |- Gene2), coherent FFL
# (Gene1 -> Gene3 -> Gene4, Gene1 -> Gene4), mixed convergence on Gene6 and a
# cascade to Gene7. Strengths are SERGIO interaction strengths K (1-5 as in the
# SERGIO demos; sign = activation / repression). Nodes: 0 = stimulus,
# k = Gene{k-1}.
_SERGIO_GRN_BY_GENES = {
    8: BenchmarkGRNConfig(
        n_genes=8,
        basal=(0.0,) * 8,  # unused by SERGIO (master-regulator rates instead)
        harissa_edges=(
            (0, 1, 3.0),
            (1, 2, 4.0),
            (1, 3, 3.0),
            (2, 3, -3.0),
            (2, 4, 4.0),
            (3, 4, -2.0),
            (2, 5, 2.0),
            (4, 5, 3.0),
            (3, 6, 4.0),
            (5, 7, -3.0),
            (6, 7, 3.0),
            (7, 8, 4.0),
            (5, 8, -2.0),
        ),
    ),
}


def get_benchmark_grn(n_genes: int, simulator: str = "default") -> BenchmarkGRNConfig:
    """Benchmark GRN for *n_genes*; ``simulator="sergio"`` returns the
    acyclic SERGIO-specific network, any other value the shared network."""
    table = _SERGIO_GRN_BY_GENES if simulator.lower() == "sergio" else _BENCHMARK_GRN_BY_GENES
    grn = table.get(n_genes)
    if grn is None:
        raise ValueError(f"Benchmark GRN is not defined for n_genes={n_genes} (simulator={simulator})")
    return grn


def write_boolode_model_files(model_dir: Path, n_genes: int) -> tuple[str, str]:
    grn = get_benchmark_grn(n_genes)
    model_dir.mkdir(parents=True, exist_ok=True)

    definition_file = model_dir / f"benchmark_grn_{n_genes}.txt"
    initial_conditions_file = model_dir / f"benchmark_grn_{n_genes}_ics.txt"

    model_rows = np.array([("Gene", "Rule"), *grn.build_boolode_rules()], dtype=object)
    ics_rows = np.array([("Genes", "Values"), *grn.build_boolode_initial_conditions()], dtype=object)

    np.savetxt(definition_file, model_rows, fmt="%s", delimiter="\t")
    np.savetxt(initial_conditions_file, ics_rows, fmt="%s", delimiter="\t")
    return definition_file.name, initial_conditions_file.name


def parse_common_cli(
    argv: Sequence[str],
    *,
    include_ko: bool = False,
    default_n_cells: int = 3000,
    default_n_genes: int = 8,
    default_n_time_bins: int = 10,
    default_seed: int = 0,
    extra_options: Sequence[str] = (),
) -> SimulatorCommonConfig:
    """Common simulator CLI.

    -o/--output_folder (required), -a/--adata_folder, -n/--n_runs,
    [-k/--ko_genes], --n_cells, --n_time_bins, --seed, and any
    simulator-specific ``--<name>`` listed in *extra_options*.
    """
    output_folder: Optional[Path] = None
    adata_folder: Optional[Path] = None
    n_runs = 2
    ko_genes = "all"
    n_cells, n_time_bins, seed = default_n_cells, default_n_time_bins, default_seed
    extra: Dict[str, str] = {}

    short_opts = "ho:a:n:k:" if include_ko else "ho:a:n:"
    long_opts = ["output_folder=", "adata_folder=", "n_runs=", "n_cells=", "n_time_bins=", "seed="]
    long_opts += [f"{name}=" for name in extra_options]
    if include_ko:
        long_opts.append("ko_genes=")

    opts, _ = getopt.getopt(list(argv), short_opts, long_opts)
    for opt, arg in opts:
        if opt in ("-o", "--output_folder"):
            output_folder = Path(arg)
        elif opt in ("-a", "--adata_folder"):
            adata_folder = Path(arg)
        elif opt in ("-n", "--n_runs"):
            n_runs = int(arg)
        elif include_ko and opt in ("-k", "--ko_genes"):
            ko_genes = str(arg)
        elif opt == "--n_cells":
            n_cells = int(arg)
        elif opt == "--n_time_bins":
            n_time_bins = int(arg)
        elif opt == "--seed":
            seed = int(arg)
        elif opt.startswith("--") and opt[2:] in extra_options:
            extra[opt[2:]] = str(arg)

    if output_folder is None:
        raise ValueError("--output_folder is required")
    if n_time_bins < 2 or n_cells < n_time_bins:
        raise ValueError(f"Need n_time_bins >= 2 and n_cells >= n_time_bins (got {n_time_bins}, {n_cells})")

    return SimulatorCommonConfig(
        output_folder=output_folder,
        adata_folder=adata_folder,
        n_runs=n_runs,
        ko_genes=ko_genes,
        n_cells=n_cells,
        n_genes=default_n_genes,
        n_time_bins=n_time_bins,
        seed=seed,
        extra=extra,
    )


def parse_ko_labels(ko_genes: str, gene_names: Sequence[str]) -> List[str]:
    value = ko_genes.strip()
    if value.lower() in {"", "none"}:
        return []
    if value.lower() == "all":
        return list(gene_names)

    labels = [label.strip() for label in value.split(",") if label.strip()]
    invalid = [label for label in labels if label not in gene_names]
    if invalid:
        raise ValueError(
            f"Invalid ko_genes value(s): {invalid}. Allowed names: {list(gene_names)}, all, none"
        )
    return labels


def build_dyngen_backbone_spec(n_genes: int) -> Optional[str]:
    """Convert the common benchmark GRN into a compact string for dyngen's R script.

    Format: ``reg>tar:eff:str;reg>tar:eff:str;…``
    - *reg* = 0 for stimulus edges, otherwise 1‑based gene index
    - *eff* = +1 (activation) or −1 (inhibition)
    - *str* = absolute interaction strength

    Returns None if no benchmark GRN is defined for *n_genes*.
    """
    try:
        grn = get_benchmark_grn(n_genes)
    except ValueError:
        return None

    parts: List[str] = []
    for regulator, target, strength in grn.harissa_edges:
        effect = 1 if strength >= 0 else -1
        abs_str = int(abs(strength))
        parts.append(f"{regulator}>{target}:{effect}:{abs_str}")

    return ";".join(parts)


def clear_previous_outputs(output_folder: Path) -> None:
    """Remove data files and ground truths of a previous simulation, so that
    leftovers (e.g. more replicates or KO files of an earlier run) are never
    converted together with the new data."""
    for pattern, folder in (("data_*.txt", output_folder), ("inter_*.npy", output_folder.parent / "True")):
        if folder.is_dir():
            for f in folder.glob(pattern):
                f.unlink()


def write_shared_files(output_folder: Path, gene_names: Sequence[str], degradation_rates) -> None:
    """Write ``panel_genes.txt`` and ``Rates/degradation_rates.txt``.

    *degradation_rates*: the simulator's own mRNA / protein degradation rates
    **in timepoint-label units** (``(d0, d1)`` for all variables, or an array
    of shape (n_genes + 1, 2) with the stimulus first), copied to
    ``adata.var['d0' / 'd1']`` and used as priors by CardamomOT. ``None`` when
    the simulator does not define such rates (no file is written, so methods
    fall back to their own priors, e.g. CardamomOT's half-life table).
    """
    rates_dir = output_folder / "Rates"
    rates_file = rates_dir / "degradation_rates.txt"

    panel_gene = np.array(
        [["0"] + [str(i + 1) for i in range(len(gene_names))], ["Stimulus"] + list(gene_names)],
        dtype=object,
    ).T
    np.savetxt(output_folder / "panel_genes.txt", panel_gene, fmt="%s", delimiter="\t")

    if degradation_rates is None:
        if rates_file.exists():
            rates_file.unlink()
        return
    rates_dir.mkdir(parents=True, exist_ok=True)
    rates = np.asarray(degradation_rates, dtype=float)
    if rates.ndim == 1:
        rates = np.tile(rates[None, :], (len(gene_names) + 1, 1))
    if rates.shape != (len(gene_names) + 1, 2):
        raise ValueError(f"degradation rates must have shape ({len(gene_names) + 1}, 2), got {rates.shape}")
    np.savetxt(rates_file, rates, fmt="%.6g", delimiter="\t")
