# Installation

This guide walks you through setting up TrajGRN-Bench on your machine.

## Prerequisites

- **Linux** or **macOS** (Windows via WSL2)
- **Docker** (for containerized method execution)
- **Conda** (Miniconda or Anaconda, for the runner environment)

## Step 1: Clone the Repository

```bash
git clone https://github.com/YannMauge/TrajGRN-Bench.git
cd TrajGRN-Bench
```

## Step 2: Create the Runner Environment

The runner environment manages simulations, orchestrates method calls, and
computes metrics. It does **not** need method-specific dependencies — those live
inside Docker containers.

```bash
conda env create -f environments/benchmark_runner.yml
conda activate benchmark_runner
```

??? tip "Verify the installation"
    ```bash
    python -c "import anndata, scanpy, pandas; print('Runner env OK')"
    ```

## Step 3: Build Method Containers

Each method runs in its own Docker container. Build all images at once:

=== "Using Make (recommended)"
    ```bash
    # Build all images (4 parallel builds)
    make -f containers/Makefile -j4 all

    # Or build a single method
    make -f containers/Makefile flecs
    ```

=== "Using the shell script"
    ```bash
    bash containers/build_all.sh -j 4
    ```

=== "Manual (one at a time)"
    ```bash
    docker build -f containers/docker/Dockerfile \
      --build-arg ENV_FILE=environments/flecs_cpu.yml \
      --build-arg ENV_NAME=flecs \
      -t benchmark/flecs:latest .
    ```

??? info "List available images"
    ```bash
    python utils/methods_registry.py container-images
    ```

## Step 4: Verify the Setup

Run a minimal test to verify everything works:

```bash
bash benchmark_run.sh --config configs/benchmark.example.yaml
```

The first run will download datasets and may take a while. Check
`benchmark/outputs_metrics/` for result files.

## Optional: Conda-Only Execution (Legacy)

If you don't want to use Docker, you can create method-specific conda
environments directly:

```bash
# Create method environments (without Docker)
conda env create -f environments/flecs_cpu.yml
conda env create -f environments/sc_dynamic.yml
conda env create -f environments/reference_fitting.yml
conda env create -f environments/cardamom_ot.yml

# Run with conda runner
python benchmark_run_config.py --config configs/benchmark.example.yaml
```

Then set `execution.default_runner: conda` in your config file.

!!! warning
    Conda environments may conflict with each other. Docker is the recommended
    approach for reproducible runs.

## GPU Support

For GPU-accelerated methods, use the GPU environment files:

```bash
docker build -f containers/docker/Dockerfile \
  --build-arg ENV_FILE=environments/flecs_gpu.yml \
  --build-arg ENV_NAME=flecs \
  -t benchmark/flecs:latest .
```

Ensure you have the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) installed.

## Next Steps

- [Run your first benchmark](usage.md)
- [Understand the configuration](usage_config.md)
- [Explore available methods](methods.md)
