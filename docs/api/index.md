# API Reference

TrajGRN-Bench's Python API is organized into several modules:

| Module | Description |
|--------|-------------|
| [`utils`](utils.md) | Utility functions and the methods registry API |
| `simulator` | Simulation drivers (Harissa, BoolODE, SERGIO) — see [Development Guide](../development.md) |
| [`post_analysis`](post_analysis.md) | Metrics computation and visualization |

## Key Classes

### Methods Registry

The central registry API — the canonical reference for all method metadata.

```python
from utils.methods_registry import get_registry, get_capabilities_table

# Load the full registry
registry = get_registry()

# Get capabilities as a DataFrame
caps = get_capabilities_table()
print(caps)
```

### Benchmark Runner

The config-driven benchmark orchestrator:

```python
from benchmark_run_config import BenchmarkRunner

runner = BenchmarkRunner("configs/benchmark.example.yaml")
runner.run()
```

### Metrics

Compute GRN and cell-level evaluation metrics:

```python
from post_analysis.compute_metrics import compute_all_metrics

metrics = compute_all_metrics(
    methods_dir="benchmark/outputs_methods/my_run/",
    ground_truth_grn="benchmark/data/True/grn.csv",
    ground_truth_adata="benchmark/data/data_1.h5ad",
)
```

## Module Index

- [utils](utils.md) — `methods_registry`, `utils`
- `simulator` — `common_config`, `Harissa`, `BoolODE`, `SERGIO` (see [Development Guide](../development.md))
- [post_analysis](post_analysis.md) — `compute_metrics`, `metrics_grn`, `metrics_cells`, `visualizations`
