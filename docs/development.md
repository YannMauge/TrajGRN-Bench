# Development Overview

This section covers how to extend TrajGRN-Bench with new methods, datasets,
metrics, and simulators.

## Architecture

```mermaid
graph TD
    benchmark_run.sh --> Method Entrypoints
    Method Entrypoints --> FLeCS & scNODE & CardamomOT & ...
    FLeCS & scNODE & CardamomOT & ... --> Method Outputs
    Method Outputs --> compute_metrics.py
    compute_metrics.py --> GRN Metrics & Cell Metrics
    GRN Metrics & Cell Metrics --> ranking_table.py
    ranking_table.py --> Ranking PNG + CSV
```

## Key Concepts

### Methods Registry

The canonical reference for all method metadata is
`methods_registry.yaml` (in the repo root). Every downstream tool
derives its configuration from this file:

- `benchmark_run.sh` — default enabled methods list
- `single_run.sh` — per-method execution parameters
- `ranking_table.py` — capabilities table
- `containers/Makefile` — Docker image build targets
- `configs/benchmark.schema.json` — config validation

[:fontawesome-solid-code: Registry API Reference →](development_registry.md)

### Method Entrypoints

Each method has a standardized entrypoint script that receives:

1. Input AnnData file path
2. Output directory path
3. Training timepoints
4. Output mode specification
5. Optional parameters (KO genes, etc.)

All methods follow the [I/O specification](methods_io.md) for
consistency.

### Simulators

Simulators generate synthetic scRNA-seq data for benchmarking. Located in
`simulator/`, they support:

- **Harissa** — default mechanistic simulator
- **BoolODE** — Boolean network-based simulator
- **SERGIO** — (in progress) gene regulation simulator

## Guides

<div class="grid cards" markdown>

-   :material-plus-circle:{ .lg .middle } **Adding a Method**

    ---

    Step-by-step guide to integrate a new inference method into the benchmark.

    [:fontawesome-solid-arrow-right: Read guide](development_add_method.md)

-   :material-database-plus:{ .lg .middle } **Adding a Dataset**

    ---

    How to prepare and add new benchmarking datasets.

    [:fontawesome-solid-arrow-right: Read guide](development_add_dataset.md)

-   :material-file-cog:{ .lg .middle } **Methods Registry**

    ---

    Deep dive into the registry system and its Python API.

    [:fontawesome-solid-arrow-right: Read guide](development_registry.md)

</div>

## Contributing

We welcome contributions! Here's how you can help improve TrajGRN-Bench:

- **Add a Method** — Integrate a new GRN or trajectory inference method. See [Adding a Method](development_add_method.md).
- **Add a Dataset** — Contribute a new benchmarking dataset. See [Adding a Dataset](development_add_dataset.md).
- **Add a Metric** — Propose a new evaluation metric for GRN or trajectory quality.
- **Report a Bug** — Open an issue with detailed reproduction steps.
- **Improve Docs** — Fix typos, add examples, or clarify explanations.
