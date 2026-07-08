# TrajGRN-Bench: Trajectory and Gene Regulatory Network Inference Benchmark

**TrajGRN-Bench** is a comprehensive benchmark of mechanistic models for **joint gene regulatory network (GRN) and single-cell RNA-seq trajectory inference**. We evaluate and compare multiple joint-inference approaches alongside complementary GRN and trajectory methods on simulated datasets.

## Key Features

- **Standardized metrics** — Unified metrics for GRN inference (AUROC, AUPRC, etc.) and trajectory reconstruction (Wasserstein, EMD, etc.).
- **Flexible configuration** — A single YAML config drives methods, datasets, metrics, and output layout.
- **Docker-based execution** — Each method runs in its own container for reproducible, conflict-free runs.
- **Interactive rankings** — Auto-generated ranking tables and plots for exploring method performance.

## Documentation

Full documentation is available at the [TrajGRN-Bench docs site](https://yannmauge.github.io/TrajGRN-Bench/) (or build locally with `mkdocs serve`).

| Section | Description |
|---|---|
| [Home](https://yannmauge.github.io/TrajGRN-Bench/) | Project overview and key features |
| [Installation](https://yannmauge.github.io/TrajGRN-Bench/install/) | Prerequisites, runner env setup, and container builds |
| [User Guide](https://yannmauge.github.io/TrajGRN-Bench/usage/) | Running benchmarks |
| [Usage](https://yannmauge.github.io/TrajGRN-Bench/usage/) | Running benchmarks, configuration, and interpreting results |
| [Methods Overview](https://yannmauge.github.io/TrajGRN-Bench/methods/) | All included methods and their capabilities |
| [Method I/O](https://yannmauge.github.io/TrajGRN-Bench/methods_io/) | Input/output specification for methods |
| [Results Gallery](https://yannmauge.github.io/TrajGRN-Bench/results/) | Interactive rankings and plots |
| [Developer Guide](https://yannmauge.github.io/TrajGRN-Bench/development/) | Overview for contributors |
| [Adding a Method](https://yannmauge.github.io/TrajGRN-Bench/development_add_method/) | Step-by-step method integration |
| [Adding a Dataset](https://yannmauge.github.io/TrajGRN-Bench/development_add_dataset/) | Step-by-step dataset integration |
| [Registry](https://yannmauge.github.io/TrajGRN-Bench/development_registry/) | Methods registry reference |
| [API Reference](https://yannmauge.github.io/TrajGRN-Bench/api/) | Python API documentation |
| [FAQ](https://yannmauge.github.io/TrajGRN-Bench/faq/) | Troubleshooting and common questions |
| [Contributing](https://yannmauge.github.io/TrajGRN-Bench/contributing/) | How to contribute |

## Quick Start

```bash
# 1. Create the runner environment
conda env create -f environments/benchmark_runner.yml

# 2. Build method Docker images
make -f containers/Makefile -j4 all

# 3. Run the benchmark
bash benchmark_run.sh --config configs/benchmark.example.yaml
```

Outputs are written to `benchmark/outputs_methods/` and `benchmark/outputs_metrics/`.

> See the [Installation guide](https://yannmauge.github.io/TrajGRN-Bench/install/) for detailed setup instructions, including optional conda-based execution.

## Methods Included

| Method | Type | Reference |
|---|---|---|
| CardamotOT | GRN + Trajectory | [bioRxiv 2026](https://doi.org/10.64898/2026.03.31.715390) |
| FLeCS | GRN + Trajectory | [arXiv:2503.20027](https://arxiv.org/abs/2503.20027) |
| Reference Fitting | GRN + Trajectory | [arXiv:2409.06879](https://arxiv.org/abs/2409.06879) |
| RENGE | GRN (perturbation-aware) | [Comm Biol 2023](https://doi.org/10.1038/s42003-023-05594-4) |
| scNODE | Trajectory (VAE+NeuralODE) | [bioRxiv 2023](https://doi.org/10.1101/2023.11.22.568346) |
| TrajectoryNet | Trajectory (CNF/OT) | [ICML 2020](https://arxiv.org/abs/2002.04461) |
| GENIE3 | GRN (baseline) | [PLoS ONE 2010](https://doi.org/10.1371/journal.pone.0012776) |
| Pearson | GRN (baseline) | Coexpression baseline |
| Waddington-OT | Trajectory (OT) | [Cell 2019](https://doi.org/10.1016/j.cell.2019.01.006) |

See the [Methods Overview](https://yannmauge.github.io/TrajGRN-Bench/methods/) on the docs site for full details, DOIs, and GitHub links.

## Repository Structure

```
.
├── benchmark_run.sh / .py        # Main benchmark entrypoints
├── benchmark_run_config.py       # Config-driven runner
├── single_run.sh                 # Single replicate / modality runner
├── ranking_table.py              # Method ranking generation
├── configs/                      # YAML config files and JSON schema
├── containers/                   # Docker/Apptainer build recipes
├── environments/                 # Conda environment YAMLs
├── methods/                      # Per-method inference scripts
├── simulator/                    # Simulation drivers (SERGIO, Harissa, BoolODE)
├── post_analysis/                # Metrics computation and visualizations
├── utils/                        # Shared utilities and registry
├── docs/                         # Documentation (MkDocs Material site)
└── benchmark/
    ├── data/                     # Input datasets (.h5ad)
    ├── outputs_methods/          # Method inference outputs
    └── outputs_metrics/          # Computed metrics
```

