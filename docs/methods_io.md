# Method Input/Output Specification

This page documents the input/output contract every method entrypoint follows.
`methods/_common.py` implements it (loading, filtering, GRN saving); see
[Adding a Method](development_add_method.md) for a template.

## CLI Convention

`single_run.sh` calls Python entrypoints as:

```
python <script> -i <inputs>/data_<i>.h5ad -o <results>/<Method>/ \
    [-u <output_mode>]            # full_full | full_train | full_test | no_traj
    [-k <ko_genes>]               # none | all | Gene1,Gene2 (KOs to predict)
    [-p <true|false>]             # perturbation training
    [extra_args ...]              # from the registry
```

Flags `-u`, `-k`, `-p` are passed when the registry entry declares
`pass_output_mode`, `pass_ko_genes`, `supports_perturbation_training`. Bash
entrypoints receive positional arguments:
`<input> <output_dir> [output_mode] [perturbation_training] [ko_genes]`.
The script must exit with a non-zero code on failure.

## Inputs (the run's `inputs/` folder)

```
inputs/
├── data_<i>.h5ad              # WT cells + training perturbations only
├── train_tps.npy              # [[index, timepoint_value], ...]
├── test_tps.npy
├── subsample_train_ids.npy    # subsample_full only: obs_names of training cells
└── halflife/                  # half-life table (CardamomOT priors for real genes)
```

```
adata
├── .X          # counts [n_cells × n_vars]; var 0 is "Stimulus" for simulated data
├── .obs        # 'timepoint' / 'time' (timepoint value), 'dataset_id' ('WT' or KO gene)
├── .var        # index = variable names; 'd0' / 'd1' = mRNA / protein degradation rates in
│               #   timepoint units, only when the simulator defines them (Harissa, BoolODE)
└── .uns        # 'simulation' (simulator name), …
```

Training cells: `obs['timepoint']` in `train_tps[:, 1]`, `dataset_id == "WT"`
unless perturbation training, and `obs_names` in `subsample_train_ids` when
present (`_common.load_training_data`).

## Outputs (in `-o`)

### `data_<i>_GRN.npy` (GRN-capable methods)

- 2-D array, shape `(n_vars, n_vars)` — with simulated data
  `(n_genes + 1, n_genes + 1)`, index 0 = Stimulus, then genes in
  `adata.var_names` order (an `n_genes × n_genes` matrix is accepted and scored
  on gene→gene edges only).
- **`grn[regulator, target]`**; positive = activation, negative = inhibition
  (unsigned scores are fine but cannot score on the sign-aware metrics).
- Any scale: metrics are threshold-free (rankings of the weights).
- `_common.save_grn` checks shape/stimulus and pads a gene-only matrix.

### `data_<i>_adata.h5ad` (trajectory-capable methods)

Predicted cells with `obs['time']` = timepoint value, genes as in the input
(stimulus optional, matched by name), on the count scale. Timepoints by
`output_mode`:

| `output_mode` | Predicted timepoints | Also contains |
|---------------|----------------------|---------------|
| `full_test` | test timepoints, from the last training timepoint | the real training cells |
| `full_train` | training timepoints, from the first one | — |
| `full_full` | all timepoints, from the first one | — |
| `no_traj` | none | — |

### `data_<i>_ko_<gene>_adata.h5ad` (KO predictions, with `-k`)

Same format, cells predicted under the knockout of `<gene>`. Scored only if
that knockout was simulated (held-out KOs are never in the inputs).

KO-specific networks (`data_<i>_ko_<gene>_GRN.npy`) may be written but are not
scored: only the common network is evaluated.

### Completion marker

Written by `single_run.sh` (`data_<i>.done`, containing a fingerprint of the
inputs and settings) when the script exits with code 0 — methods do not write
it.

## Examples

- `methods/pearson/inference_pearson.py` — minimal baseline on `_common`
- `methods/dynGENIE3/dynGENIE3_ontrain.py` — GRN method on `_common`
- `methods/referencefitting/notebooks/ref_fit.py` — GRN + trajectories + KO predictions
- `methods/CardamomOT/run_carda_benchmark.sh` — bash entrypoint
