# Installation

This guide walks you through setting up TrajGRN-Bench on your machine.

## Prerequisites

- **Linux** (x86-64 or ARM64) or **macOS on Apple Silicon** (Intel Macs are not
  supported; Windows via WSL2)
- **[pixi](https://pixi.sh)** ≥ 0.84 — environment manager (all environments
  are declared in `pixi.toml` and locked in `pixi.lock`)
- **Docker** with `buildx` (default way to run the methods)
- For GPU runs: an NVIDIA GPU, CUDA 12 drivers and the
  [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)

```bash
# macOS
brew install pixi
# Linux
curl -fsSL https://pixi.sh/install.sh | bash
```

## Step 1: Clone the Repository

```bash
git clone https://github.com/YannMauge/TrajGRN-Bench.git
cd TrajGRN-Bench
```

## Step 2: Create the Runner Environment

The runner environment (pixi env `runner`) runs the simulators, prepares the
train/test splits, orchestrates method calls and computes the metrics. Method
dependencies live in their own environments / images.

```bash
pixi install -e runner
# dyngen (R) is not on conda-forge: install it once from a dated CRAN snapshot
pixi run -e runner install-dyngen
```

??? tip "Verify the installation"
    ```bash
    pixi run -e runner python -c "import anndata, scanpy, harissa; print('Runner env OK')"
    ```

## Step 3: Build Method Containers

Each method runs in a Docker image built from its pixi environment
(`benchmark/<pixi_env>:latest`, see `execution.pixi_env` in
`methods_registry.yaml`):

| Image | Methods |
|-------|---------|
| `benchmark/baselines` | GENIE3, dynGENIE3, PEARSON |
| `benchmark/cardamom` | CardamomOT |
| `benchmark/flecs` | FLeCS, FLeCS-TPs |
| `benchmark/scnode` | scNODE |
| `benchmark/trajectorynet` | TrajectoryNet, WaddingtonOT |
| `benchmark/reference-fitting` | reference_fitting |
| `benchmark/renge` | RENGE |

=== "Using Make"
    ```bash
    make -f containers/Makefile all            # CPU images, native platform
    make -f containers/Makefile all-cuda       # + CUDA images (linux/amd64)
    make -f containers/Makefile flecs scnode   # selected images
    make -f containers/Makefile list
    ```

=== "Using the shell script"
    ```bash
    bash containers/build_all.sh -j 2                  # CPU images
    bash containers/build_all.sh --cuda -j 2           # + CUDA images
    bash containers/build_all.sh --platform linux/amd64 flecs
    ```

=== "Manual (one image)"
    ```bash
    docker buildx build -f containers/docker/Dockerfile \
      --build-arg PIXI_ENV=flecs --platform linux/arm64 --load \
      -t benchmark/flecs:latest .
    ```

Images contain only the locked environment; the repository is mounted at run
time. Record the git commit and image digests with your results.

## Step 4: Verify the Setup

```bash
bash benchmark_run.sh --config configs/benchmark.example.yaml
```

Check `benchmark/outputs_metrics/` for result files.

## Optional: Run Methods Without Docker (pixi)

Every method environment can be created locally:

```bash
pixi install -e flecs          # or baselines, cardamom, scnode, trajectorynet,
                               #    reference-fitting, renge
```

Then set `execution.default_runner: pixi` (or per method in
`execution.method_runners`) in your config.

## GPU Support

PyTorch-based methods have CUDA 12 variants (`cardamom-cuda`, `flecs-cuda`,
`scnode-cuda`, `trajectorynet-cuda`, `reference-fitting-cuda`; linux/amd64).
Build their images with `make -f containers/Makefile all-cuda` and set

```yaml
execution:
  gpu: true
```

Methods with a CUDA variant then use `benchmark/<env>-cuda` images (or the
`<env>-cuda` pixi env with `default_runner: pixi`); Docker runs get
`--gpus all`, Apptainer runs `--nv`. Methods without a CUDA variant (GENIE3,
dynGENIE3, PEARSON, RENGE) run on CPU.

## Conda (unsupported)

The conda environment files were removed in favour of `pixi.toml`. The
`conda` runner modes remain for users who maintain their own environments:
`execution.runner_tool: conda` + `runner_env`, and `default_runner: conda`
with the per-method `conda_env` names of the registry.

## Next Steps

- [Run benchmarks & understand results](usage.md)
- [Explore available methods](methods.md)
