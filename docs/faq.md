# FAQ & Troubleshooting

## Installation Issues

### "conda: command not found"

Install Miniconda:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### "Docker daemon not running"

Start Docker:

```bash
sudo systemctl start docker        # Linux
# Or use Docker Desktop on macOS/Windows
```

Verify:

```bash
docker run hello-world
```

### "permission denied" when running Docker

Add yourself to the `docker` group:

```bash
sudo usermod -aG docker $USER
newgrp docker
```

### Conda environment creation fails

Try using Mamba (faster solver):

```bash
conda install -c conda-forge mamba
mamba env create -f environments/benchmark_runner.yml
```

### Container build fails with "no space left on device"

Docker images can be large. Clean up:

```bash
docker system prune -a
```

## Runtime Issues

### Method fails with "Killed" or OOM

The method ran out of memory. Options:

1. Reduce the dataset size
2. Subsample cells (use `train_data: subsample_full`)
3. Increase Docker memory limit in Docker Desktop settings

### ".done file not found" error

This means the method crashed or was interrupted. Check the method's log output
in the terminal. Re-run with `restart_mode: rerun` to force re-execution.

### Container can't access GPU

Ensure NVIDIA Container Toolkit is installed:

```bash
nvidia-container-cli info
```

If not installed:

```bash
# Ubuntu/Debian
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
  sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
```

### "AnnData file not found" or "KeyError: timepoint"

Your dataset is missing the required structure. Verify:

```python
import scanpy as sc
adata = sc.read_h5ad("benchmark/data/data_1.h5ad")
print(adata.obs.columns)  # Should include 'timepoint' or 'time'
print(adata.shape)
```

## Results Questions

### Why are some metrics NaN?

A method may not support certain output types. For example, GENIE3 doesn't
produce trajectory outputs, so cell-level metrics will be NaN.

### How do I compare results across runs?

Use `post_analysis/compute_metrics.py` with multiple run directories:

```bash
python post_analysis/compute_metrics.py \
  --runs run_1 run_2 run_3 \
  --output comparison.csv
```

### Ranking table is empty

Ensure both `grn_metrics.csv` and `cell_metrics.csv` exist in the metrics output
directory. Run `compute_metrics.py` first if needed.

## Development Questions

### How do I add a new method?

See [Adding a Method](development_add_method.md). The short version: add an
entry to `methods_registry.yaml` and write an entrypoint script.

### Can I use a language other than Python?

Yes. Set `entrypoint.type: bash` in the registry and write a bash script that
wraps your method. The script must accept the same CLI arguments.

### How do I test my method without Docker?

Set `execution.default_runner: conda` in your config and create the conda
environment. Run the entrypoint script directly for debugging.

## Still Stuck?

- Check existing [GitHub Issues](https://github.com/YOUR_ORG/TrajGRN-Bench/issues)
- Open a new issue with:
    - Full error message
    - Your config file (sanitized)
    - OS and Docker version
    - Steps to reproduce
