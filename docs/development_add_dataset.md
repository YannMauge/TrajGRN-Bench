# Adding a New Dataset

This guide explains how to prepare and integrate a new benchmarking dataset.

## Dataset Format

Datasets must be in **AnnData (HDF5)** format (`.h5ad`). See the full [I/O
specification](methods_io.md) for the complete schema.

### Minimal Requirements

```
adata
├── .X          # Gene expression matrix [n_cells × n_genes] (dense or sparse)
├── .obs        # Cell metadata (must contain 'timepoint' or 'time' column)
├── .var        # Gene metadata (index = gene names)
└── .uns        # Optional: 'simulation' flag, ground-truth network, etc.
```

### Preparing Your Data

#### From a CSV / TSV

```python
import pandas as pd
import scanpy as sc
import anndata

# Load expression matrix
expr = pd.read_csv("expression.csv", index_col=0)  # genes × cells

# Load cell metadata
meta = pd.read_csv("metadata.csv", index_col=0)     # cells × attributes

# Create AnnData
adata = anndata.AnnData(X=expr.T.values, obs=meta, var=pd.DataFrame(index=expr.index))

# Ensure timepoint column exists
assert "timepoint" in adata.obs.columns, "Metadata must have 'timepoint' column"

# Save
adata.write_h5ad("benchmark/data/my_dataset.h5ad")
```

#### From a Seurat Object (R)

```r
library(Seurat)
library(SeuratDisk)

# Load Seurat object
obj <- readRDS("my_data.rds")

# Convert and save
SaveH5Seurat(obj, filename = "my_data.h5seurat")
Convert("my_data.h5seurat", dest = "h5ad")
```

### Timepoint Files

Create NumPy arrays for train/test timepoint splits:

```python
import numpy as np

# Define which timepoints are training vs test
# Format: [[timepoint_id, timepoint_value], ...]
train_tps = np.array([[0, 0.0], [1, 1.0], [2, 2.0]])
test_tps = np.array([[3, 3.0], [4, 4.0]])

np.save("benchmark/data/my_train_tps.npy", train_tps)
np.save("benchmark/data/my_test_tps.npy", test_tps)
```

### Ground Truth (Optional)

If you have a ground-truth GRN:

```python
# N×N adjacency matrix with gene names as index/columns
grn = pd.read_csv("ground_truth_grn.csv", index_col=0)
# Place in benchmark/data/True/ or alongside your dataset
```

## Dataset Placement

Place your files in the standard locations:

```
benchmark/data/
├── my_dataset.h5ad          # Your AnnData file
├── my_train_tps.npy         # Training timepoints
├── my_test_tps.npy          # Test timepoints
└── True/                    # Ground truth (optional)
    └── my_dataset_grn.csv   # Ground-truth GRN
```

## Configuration

Add your dataset to the benchmark config:

```yaml
data:
  adata_files:
    - benchmark/data/my_dataset.h5ad
  train_tps: benchmark/data/my_train_tps.npy
  test_tps: benchmark/data/my_test_tps.npy
```

For multiple datasets:

```yaml
data:
  adata_files:
    - benchmark/data/data_1.h5ad
    - benchmark/data/my_dataset.h5ad
```

## Testing

Verify your dataset loads correctly:

```python
import scanpy as sc
import numpy as np

adata = sc.read_h5ad("benchmark/data/my_dataset.h5ad")
print(f"Shape: {adata.shape}")
print(f"Timepoints: {adata.obs['timepoint'].unique()}")
print(f"Genes: {adata.var_names[:5].tolist()}")
assert "timepoint" in adata.obs.columns or "time" in adata.obs.columns
```
