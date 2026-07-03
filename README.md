# TrajGRN-Bench: Trajectory and Gene Regulatory Network Inference Benchmark

**TrajGRN-Bench** is a comprehensive benchmark of mechanistic models for joint gene regulatory network (GRN) and single-cell RNA-seq trajectory inference. We evaluate and compare multiple joint-inference approaches alongside complementary GRN and trajectory methods on simulated datasets. 

We use :
- Standardized metrics for GRN and trajectory inference.
- A flexible configuration system to run different methods and metrics on different datasets.
- Docker-based execution for repeatable runs across different computational environments.

See [METHODS.md](METHODS.md) for a concise list of included methods and their DOI, GitHub links, and short descriptions.

## Installation

### Default setup (Docker method execution)

This repository now defaults to **Docker-based method execution** to avoid installing
multiple method-specific conda environments. You only need a single runner environment
to orchestrate simulations and metrics.

1. Install Docker.
2. Create the runner environment:
```
conda env create -f ./environments/benchmark_runner.yml
```
3. Build the method images (see `containers/README.md` for full examples):
```
docker build -f containers/docker/Dockerfile --build-arg ENV_FILE=environments/flecs_cpu.yml --build-arg ENV_NAME=flecs -t benchmark/flecs:latest .
docker build -f containers/docker/Dockerfile --build-arg ENV_FILE=environments/sc_dynamic.yml --build-arg ENV_NAME=sc_dynamic -t benchmark/scnode:latest .
docker build -f containers/docker/Dockerfile --build-arg ENV_FILE=environments/reference_fitting.yml --build-arg ENV_NAME=reference_fitting -t benchmark/referencefitting:latest .
docker build -f containers/docker/Dockerfile --build-arg ENV_FILE=environments/cardamom_ot.yml --build-arg ENV_NAME=cardamom_env -t benchmark/cardamomot:latest .
```

Run the benchmark (default):
```
bash benchmark_run.sh --config configs/benchmark.example.yaml
```

### Optional: method conda environments (legacy)

If you prefer local conda environments instead of containers, create the method envs:

If you have a GPU and want to run on GPU:
```
conda env create -f ./environments/cardamom_ot.yml
conda env create -f ./environments/sc_dynamic.yml
conda env create -f ./environments/reference_fitting.yml
conda env create -f ./environments/flecs_gpu.yml
```

Outputs will be found in './benchmark/outputs_methods/' and './benchmark/outputs_metrics/'

The example config defaults to Docker runners (`execution.default_runner: docker`).
If you want to use conda environments, set `execution.default_runner: conda` and
fill `execution.conda_envs` in your config.

You can also invoke the config runner directly:

```
conda run -n benchmark_runner python benchmark_run_config.py --config configs/benchmark.example.yaml
```

## Repository structure

```
.
├── benchmark_run.sh              # Main entrypoint (supports --config)
├── benchmark_run.py              # Benchmark orchestration logic
├── benchmark_run_config.py       # Config-driven runner
├── compute_metrics.py            # Metrics orchestration
├── single_run.sh                 # Runs a single inference run
├── visualize_cells.py            # Cell-level metric plots
├── visualize_grn.py              # GRN-level metric plots
├── metrics_cells.py              # Cell-level metrics
├── metrics_grn.py                # GRN-level metrics
├── ranking_table.py              # Ranking table based on metrics
├── README.md
├── benchmark/
│   ├── data/                     # Benchmark inputs
│   ├── outputs_methods/          # Method outputs (GRNs, adatas, etc.)
│   └── outputs_metrics/          # Metrics outputs
├── configs/                      # Benchmark configs and schema
├── containers/                   # Docker/Apptainer build recipes
├── environments/                 # Conda environment files
├── methods/                      # Methods used in the benchmark
├── simulator/                    # Simulation drivers and utilities
└── utils/
    └── utils.py
```
