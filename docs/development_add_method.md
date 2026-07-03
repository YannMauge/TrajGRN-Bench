# Adding a New Method

This guide walks you through integrating a new inference method into
TrajGRN-Bench. In most cases, you only need to touch **two files**: the method
entrypoint script and `methods_registry.yaml`.

## Step-by-Step

### 1. Write the Method Entrypoint Script

Create a Python (or bash) script in `methods/<your_method>/` that follows the
standard CLI convention:

```bash
python methods/<your_method>/infer.py \
    --adata <input.h5ad> \
    --output_dir <output_dir> \
    --train_tps <train_tps.npy> \
    --test_tps <test_tps.npy> \
    --output_mode full_test \
    [--ko_genes GENE1,GENE2]
```

Your script must produce:

| File | Required | Description |
|------|:--------:|-------------|
| `predicted_grn.csv` | If GRN-capable | $N \times N$ adjacency matrix, CSV, genes as row/col names |
| `simulated_adata.h5ad` | If trajectory-capable | AnnData with `.X` containing simulated expression |
| `.done` | Always | Empty marker file indicating successful completion |

See the full [Method I/O Specification](methods_io.md) for details.

??? example "Minimal Python entrypoint template"
    ```python
    #!/usr/bin/env python3
    """Minimal method entrypoint for TrajGRN-Bench."""
    import argparse, os, numpy as np, pandas as pd, scanpy as sc

    def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--adata", required=True)
        parser.add_argument("--output_dir", required=True)
        parser.add_argument("--train_tps", required=True)
        parser.add_argument("--test_tps", required=True)
        parser.add_argument("--output_mode", default="full_test")
        parser.add_argument("--ko_genes", default=None)
        args = parser.parse_args()

        os.makedirs(args.output_dir, exist_ok=True)

        # Load data
        adata = sc.read_h5ad(args.adata)
        train_tps = np.load(args.train_tps)

        # --- Your method logic here ---
        grn_matrix = np.random.randn(adata.n_vars, adata.n_vars)

        # Save GRN
        pd.DataFrame(grn_matrix, index=adata.var_names,
                     columns=adata.var_names).to_csv(
                     os.path.join(args.output_dir, "predicted_grn.csv"))

        # Mark completion
        open(os.path.join(args.output_dir, ".done"), "w").close()

    if __name__ == "__main__":
        main()
    ```

### 2. Create a Conda Environment File

Create `environments/<your_method>.yml` with all dependencies:

```yaml
name: your_method
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - numpy
  - pandas
  - scanpy
  - pip
  - pip:
      - your-method-package
```

### 3. Register in `methods_registry.yaml`

Add your method under the `methods:` key:

```yaml
methods:
  # ... existing methods ...

  YourMethod:
    display_name: Your Method
    description: >
      Short description of your method and what it does.
    github: https://github.com/yourorg/your-method
    publication: "Your Paper (Year) DOI:..."
    entrypoint:
      type: python
      script: methods/your_method/infer.py
      pass_output_mode: true
      pass_ko_genes: false
      supports_perturbation_training: false
    execution:
      default_runner: docker
      conda_env: your_method
      container:
        image: benchmark/yourmethod:latest
        env_file: environments/your_method.yml
        env_name: your_method
    capabilities:
      grn_inference: true
      trajectory_reconstruction: true
      perturbation_training: false
```

Then add `YourMethod` to the `iteration_order` list at the bottom of the file.

### 4. Build the Docker Image

```bash
docker build -f containers/docker/Dockerfile \
  --build-arg ENV_FILE=environments/your_method.yml \
  --build-arg ENV_NAME=your_method \
  -t benchmark/yourmethod:latest .
```

### 5. Test Locally

```bash
# Test your entrypoint directly
python methods/your_method/infer.py \
  --adata benchmark/data/data_1.h5ad \
  --output_dir /tmp/test_output \
  --train_tps benchmark/data/train_tps.npy \
  --test_tps benchmark/data/test_tps.npy \
  --output_mode full_test

# Test via the benchmark runner (single method)
python benchmark_run_config.py --config configs/benchmark.example.yaml
```

### 6. Verify Metrics

```bash
python post_analysis/compute_metrics.py \
  --methods-dir benchmark/outputs_methods/<run_name>/YourMethod \
  --output-dir /tmp/test_metrics
```

## Checklist

- [ ] Entrypoint script follows CLI convention
- [ ] Produces `predicted_grn.csv` (if GRN-capable)
- [ ] Produces `simulated_adata.h5ad` (if trajectory-capable)
- [ ] Creates `.done` marker on completion
- [ ] Conda environment file in `environments/`
- [ ] Entry added to `methods_registry.yaml`
- [ ] Entry added to `iteration_order` in registry
- [ ] Docker image builds successfully
- [ ] Metrics compute correctly

## Need Help?

- Check the [Method I/O Specification](methods_io.md)
- See examples in `methods/CardamomOT/`, `methods/flecs/`, etc.
- Open an issue on GitHub
