"""SERGIO driver shared by the WT and KO entry scripts.

Runs **unmodified upstream SERGIO** (``SERGIO/sergio.py``, PayamDiba/SERGIO
master, Dibaeinia & Sinha 2020) in its published *differentiation* (dynamics)
mode. Only orchestration lives here:

* Network: the SERGIO-specific acyclic benchmark GRN
  (``get_benchmark_grn(n, "sergio")``) — upstream SERGIO cannot simulate cycles
  or self-loops. Gene0 is the master regulator (stimulus entry point).
* Time: a linear differentiation chain of ``n_time_bins`` states
  (bifurcation matrix with 0.2 on the super-diagonal, as in the SERGIO demo);
  state b is labelled with the shared timepoint label of bin b (0, 10, …, 90).
  Gene0's production rate ramps from 0.2 (state 0, stimulus off) to 4.0.
* Parameters: the SERGIO differentiation demo values — noise 0.3 ('dpd'),
  splicing noise 0.1, decay 0.8, splice ratio 1.5, sampling_state 1, Hill
  coefficient 2. Expression = unspliced + spliced, converted to UMI counts
  with SERGIO's ``convert_to_UMIcounts_dynamics`` (Poisson).
* Knockouts (no native support upstream): through the input files only — a
  KO master regulator gets production rate 0 in every state, a KO target gets
  interaction strength 0 on all its incoming edges. Its production is then 0,
  and with 'dpd' noise (∝ sqrt(production), sqrt(decay·x)) it stays at 0;
  its outgoing edges are unchanged, so its targets respond (upstream's Hill
  function returns 0 / 1 for a regulator at exactly 0). WT and KOs use the
  same seed (common random numbers).
* Seeding: upstream uses numpy's global RNG only, so ``np.random.seed`` with a
  per-replicate seed makes runs reproducible and replicates distinct.
* Compatibility: upstream uses ``np.int`` / ``np.float`` (removed in
  numpy 1.24); the aliases are restored before importing SERGIO.
* KO convergence: upstream's dynamics convergence test is
  ``|mean - ss| < 0.2 * ss``, which can never hold for a gene whose steady
  state is exactly 0 (a KO gene), so the simulation would never end.
  ``_SergioKO`` (a subclass; upstream file untouched) also accepts
  ``mean == ss == 0``. WT genes all have ss > 0, so WT runs are unaffected.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent.parent))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

for _alias, _type in (("int", int), ("float", float), ("bool", bool)):
    if not hasattr(np, _alias) or _alias in getattr(np, "__deprecated_attrs__", {}):
        setattr(np, _alias, _type)

from simulator.common_config import (  # noqa: E402
    SimulatorCommonConfig,
    clear_previous_outputs,
    get_benchmark_grn,
    write_shared_files,
)
from SERGIO.sergio import sergio  # noqa: E402


class _SergioKO(sergio):
    """Upstream SERGIO; only the dynamics convergence test also accepts genes
    held at a steady state of exactly 0 (knocked-out genes)."""

    def check_convergence_dynamics_(self, binID, num_init_cells):
        for g in self.binDict[binID]:
            if not g.converged_ and g.ss_U_ == 0:
                if all(np.all(np.asarray(g.Conc[i][-10:]) == 0) for i in range(num_init_cells)):
                    g.setConverged()
            if g.converged_ and not g.converged_S_ and g.ss_S_ == 0:
                if all(np.all(np.asarray(g.Conc_S[i][-10:]) == 0) for i in range(num_init_cells)):
                    g.setConverged_S()
        return super().check_convergence_dynamics_(binID, num_init_cells)

# SERGIO differentiation-demo parameters (Demo/demo_run.ipynb)
NOISE_U = 0.3
NOISE_S = 0.1
DECAY = 0.8
SPLICE_RATIO = 1.5
SAMPLING_STATE = 1
COOP = 2
TRANSITION = 0.2
MR_RATE_OFF, MR_RATE_MAX = 0.2, 4.0
EXTRA_OPTIONS = ()


def _write_inputs(folder: Path, n_genes: int, n_bins: int, ko_gene: Optional[int]) -> tuple:
    """SERGIO input files (targets with shared coop state, master regulators)."""
    grn = get_benchmark_grn(n_genes, "sergio")
    incoming: Dict[int, List[tuple]] = {g: [] for g in range(n_genes)}
    masters = []
    for reg, tgt, k in grn.harissa_edges:
        if reg == 0:
            masters.append(tgt - 1)
        else:
            incoming[tgt - 1].append((reg - 1, k))
    targets_file, regs_file = folder / "targets.txt", folder / "regs.txt"
    with open(targets_file, "w") as fh:
        for tgt in range(n_genes):
            regs = incoming[tgt]
            if not regs:
                continue
            ks = [0.0 if tgt == ko_gene else float(k) for _, k in regs]
            row = [float(tgt), float(len(regs))] + [float(r) for r, _ in regs] + ks
            fh.write(",".join(repr(v) for v in row) + "\n")
    ramp = [MR_RATE_OFF] + list(np.linspace(1.0, MR_RATE_MAX, n_bins - 1))
    with open(regs_file, "w") as fh:
        for mr in sorted(set(masters)):
            rates = [0.0] * n_bins if mr == ko_gene else ramp
            fh.write(",".join(repr(float(v)) for v in [mr] + list(rates)) + "\n")
    return str(targets_file), str(regs_file)


def _simulate(n_genes: int, n_bins: int, cells_per_bin: int, seed: int, ko_gene: Optional[int]) -> np.ndarray:
    """Return UMI counts (bins, genes, cells_per_bin) of U + S."""
    bif = np.zeros((n_bins, n_bins))
    for b in range(n_bins - 1):
        bif[b, b + 1] = TRANSITION
    np.random.seed(seed)
    with tempfile.TemporaryDirectory() as tmp:
        targets, regs = _write_inputs(Path(tmp), n_genes, n_bins, ko_gene)
        sim = _SergioKO(
            number_genes=n_genes, number_bins=n_bins, number_sc=cells_per_bin,
            noise_params=NOISE_U, decays=DECAY, sampling_state=SAMPLING_STATE,
            noise_params_splice=NOISE_S, noise_type="dpd", dynamics=True,
            bifurcation_matrix=bif, splice_ratio=SPLICE_RATIO,
        )
        sim.build_graph(input_file_taregts=targets, input_file_regs=regs, shared_coop_state=COOP)
        sim.simulate_dynamics()
        u, s = sim.getExpressions_dynamics()
        u_counts, s_counts = sim.convert_to_UMIcounts_dynamics(u, s)
    return u_counts + s_counts


def simulate_runs(cfg: SimulatorCommonConfig, ko_labels: Sequence[str], write_ko_files: bool) -> None:
    """Simulate ``cfg.n_runs`` replicates (WT and, if requested, each KO)."""
    n_bins = cfg.n_time_bins
    cells_per_bin = cfg.n_cells // n_bins
    if cells_per_bin * n_bins != cfg.n_cells:
        print(f"[SERGIO] n_cells={cfg.n_cells} not divisible by {n_bins} states: "
              f"using {cells_per_bin * n_bins} cells ({cells_per_bin} per state)")
    labels = cfg.time_labels()
    times = np.repeat(labels, cells_per_bin)
    clear_previous_outputs(cfg.output_folder)

    # SERGIO's decay is per SERGIO time step along a pseudo-time chain (states
    # are run to convergence, so they have no fixed duration) and SERGIO has no
    # protein: its rates are not defined in timepoint units -> no priors.
    write_shared_files(cfg.output_folder, cfg.gene_names, None)
    true_dir = cfg.output_folder.parent / "True"
    true_dir.mkdir(parents=True, exist_ok=True)
    truth = get_benchmark_grn(cfg.n_genes, "sergio").build_true_interaction_matrix()

    for run_index in range(1, cfg.n_runs + 1):
        seed = cfg.run_seed(run_index)
        print(f"Run {run_index} (seed {seed})...")
        np.save(true_dir / f"inter_{run_index}", truth)
        conditions = [("WT", None)] + [(lab, cfg.gene_names.index(lab)) for lab in ko_labels]
        for label, ko_gene in conditions:
            counts = _simulate(cfg.n_genes, n_bins, cells_per_bin, seed, ko_gene)
            # (bins, genes, cells) -> (cells, genes), bins in time order
            cells = np.concatenate([counts[b].T for b in range(n_bins)], axis=0).astype(int)
            name = f"data_{run_index}.txt" if not write_ko_files else f"data_{run_index}_ko_{label}.txt"
            _write_data_file(cfg.output_folder / name, times, cells, cfg.stimulus_value)


def _write_data_file(path: Path, times: np.ndarray, counts: np.ndarray, stimulus_value: int) -> None:
    n_cells, n_genes = counts.shape
    data = np.zeros((n_cells + 1, n_genes + 2), dtype=int)
    data[0, 1:] = np.arange(n_genes + 1)
    data[1:, 0] = times
    data[1:, 1] = int(stimulus_value) * (times > 0)
    data[1:, 2:] = counts
    np.savetxt(path, data.T, fmt="%d", delimiter="\t")
