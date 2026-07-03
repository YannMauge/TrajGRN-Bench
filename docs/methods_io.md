# Method Input/Output Specification

This page documents the standardized input/output formats for all methods in
TrajGRN-Bench. Every method must follow these conventions.

## CLI Convention

All method entrypoints must accept these arguments:

```
python <script> \
    --adata <path>          # Input AnnData (.h5ad)
    --output_dir <path>      # Where to write results
    --train_tps <path>       # Training timepoints (.npy)
    --test_tps <path>        # Test timepoints (.npy)
    --output_mode <mode>     # full_full | full_train | full_test | no_traj
    [--ko_genes <genes>]     # Optional: comma-separated KO target genes
```

## Input: AnnData

```
adata
├── .X          # Gene expression [n_cells × n_genes] — dense or sparse
├── .obs        # Cell metadata — must have 'timepoint' or 'time'
├── .var        # Gene metadata — index = gene names
└── .uns        # Optional unstructured metadata
```

## Input: Timepoints

NumPy array, shape `(n_timepoints, 2)`:

```python
# [[timepoint_id, timepoint_value], ...]
[[0, 0.0], [1, 1.0], [2, 2.0]]
```

Filter training data: `adata[adata.obs['timepoint'].isin(train_tps[:, 1])]`

## Required Outputs

### `predicted_grn.csv` (if GRN-capable)

$N \times N$ adjacency matrix with gene names as both row and column labels:

```csv
, GeneA, GeneB, GeneC
GeneA, 0.0, 0.8, -0.3
GeneB, 0.1, 0.0, 0.5
GeneC, 0.0, 0.2, 0.0
```

- Values represent edge weights (positive = activation, negative = inhibition)
- Must be symmetric? No — directed networks are allowed

### `simulated_adata.h5ad` (if trajectory-capable)

AnnData object with simulated expression in `.X`:

```python
import scanpy as sc
adata_sim = sc.read_h5ad("simulated_adata.h5ad")
# adata_sim.X shape: [n_simulated_cells × n_genes]
# adata_sim.obs should contain 'timepoint'
```

### `.done` marker (always required)

An empty file signaling successful completion:

```bash
touch <output_dir>/.done
```

## Optional Outputs

| File | Description |
|------|-------------|
| `latent_representation/` | Method-specific latent embeddings |
| `kinetic_parameters.csv` | Inferred kinetic rates per gene |
| `*.png` | Diagnostic plots |

## Examples

See working implementations in:

- `methods/CardamomOT/infer_test.py` — Python entrypoint
- `methods/flecs/flecs_train.py` — Python entrypoint with extra args
- `methods/scnode/scnode.py` — Python entrypoint
- `methods/pearson/inference_pearson.py` — Simple baseline example
