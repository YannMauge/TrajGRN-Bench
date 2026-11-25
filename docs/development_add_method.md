# Adding a New Method

Integrating a method touches four places: the **entrypoint script**, a **pixi
environment** (`pixi.toml`), an entry in **`methods_registry.yaml`**, and its
**Docker image** (built automatically from the pixi environment).

## 1. Write the Method Entrypoint Script

Create `methods/<your_method>/infer.py` (or a bash script). `single_run.sh`
calls it as:

```bash
python methods/<your_method>/infer.py -i <inputs>/data_1.h5ad -o <results>/<Method>/ \
    [-u <output_mode>] [-k <ko_genes>] [-p <perturbation_training>] [extra_args...]
```

| Flag | Passed when (registry) | Meaning |
|------|------------------------|---------|
| `-i` | always | input AnnData in the run's `inputs/` folder |
| `-o` | always | output folder (ends with `/`) |
| `-u` | `pass_output_mode: true` | `full_full`, `full_train`, `full_test` or `no_traj` |
| `-k` | `pass_ko_genes: true` | KO genes to predict: `none`, `all` or `Gene1,Gene2` |
| `-p` | `supports_perturbation_training: true` | `true` / `false` |

Bash entrypoints get positional arguments instead:
`<input> <output_dir> [output_mode] [perturbation_training] [ko_genes]`.

Next to the input file are `train_tps.npy` / `test_tps.npy` (rows
`[index, timepoint_value]`) and, for `subsample_full`,
`subsample_train_ids.npy`. The input contains WT cells plus the *training*
perturbations only (`obs['dataset_id']`); held-out knockouts are never visible
to methods.

### Outputs

| File | Required | Description |
|------|:--------:|-------------|
| `data_<i>_GRN.npy` | if GRN-capable | 2-D, shape (n_genes+1, n_genes+1), row/col 0 = Stimulus, `grn[regulator, target]`, signed if the method infers signs |
| `data_<i>_adata.h5ad` | if trajectory-capable | predicted cells with `obs['time']` (timepoint values); for `full_test` also the training cells |
| `data_<i>_ko_<gene>_adata.h5ad` | if `-k` and KO-capable | predicted cells under the knockout |

The `.done` marker is written by `single_run.sh` (with an input fingerprint);
the script only has to exit with a non-zero code on failure.

### Use the shared helpers

`methods/_common.py` implements the contract; using it avoids the classic bugs
(wrong orientation, KO cells in the training set, subsampling after
re-indexing, stimulus counted in library sizes):

```python
#!/usr/bin/env python3
"""Minimal method entrypoint for TrajGRN-Bench."""
import getopt, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # methods/
from _common import load_training_data, log_normalize_genes, save_grn, seed_everything


def main(argv):
    opts, _ = getopt.getopt(argv, "i:o:u:k:p:")
    args = dict(opts)
    inputfile, outdir = args["-i"], args["-o"]
    perturbation_training = args.get("-p", "false").lower() == "true"

    seed_everything(0)
    # WT (+ training KOs if perturbation_training) cells at training timepoints
    adata, train_tps, test_tps = load_training_data(inputfile, perturbation_training)
    x = log_normalize_genes(adata.X, adata.var_names)

    grn = np.corrcoef(x, rowvar=False)            # your inference here: regulator x target
    save_grn(outdir, inputfile, grn, adata.var_names)


if __name__ == "__main__":
    main(sys.argv[1:])
```

See [Method I/O Specification](methods_io.md) for the full contract.

## 2. Add a pixi Environment

In `pixi.toml`, add a feature and an environment (conda-forge first; PyPI only
for packages missing there):

```toml
[feature.yourmethod.dependencies]
python = "3.12.*"
pytorch = "*"

[feature.yourmethod.pypi-dependencies]
your-method-package = "==1.2.3"

[environments]
yourmethod = { features = ["sci", "yourmethod"], no-default-feature = true }
# optional GPU variant (linux-64): PyTorch with CUDA 12
yourmethod-cuda = { features = ["sci", "yourmethod", "cuda"], no-default-feature = true }
```

Then lock and test locally:

```bash
pixi lock
pixi install -e yourmethod
pixi run -e yourmethod python methods/your_method/infer.py -i ... -o /tmp/out/
```

The lock must cover linux-64, linux-aarch64 and osx-arm64. Pin exact versions
when the method's publication does.

## 3. Register in `methods_registry.yaml`

```yaml
methods:
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
      extra_args: []
    execution:
      default_runner: docker
      pixi_env: yourmethod
      container:
        image: benchmark/yourmethod:latest
    capabilities:
      grn_inference: true
      trajectory_reconstruction: true
      perturbation_training: false
```

Add `YourMethod` to `iteration_order` at the bottom of the file and to
`benchmark.run_methods` in your config.

## 4. Build the Docker Image

```bash
make -f containers/Makefile yourmethod          # benchmark/yourmethod:latest
make -f containers/Makefile yourmethod-cuda     # GPU variant (linux/amd64)
```

## 5. Run and Check the Metrics

```yaml
benchmark:
  run_methods:
    YourMethod: true
```

```bash
bash benchmark_run.sh --config configs/your_config.yaml
```

Results: `benchmark/outputs_metrics/<run>/grn_metrics.csv`,
`cell_metrics.csv` and the ranking tables. A failing method is listed in
`failed_methods.tsv` and the run continues with the other methods.

## Checklist

- [ ] Entrypoint parses `-i -o` (+ `-u -k -p` as declared) and exits non-zero on failure
- [ ] Trains on WT cells (or WT + training perturbations) at training timepoints only
- [ ] `data_<i>_GRN.npy`: 2-D, stimulus at index 0, `grn[regulator, target]`
- [ ] `data_<i>_adata.h5ad` with `obs['time']` for the requested `output_mode`
- [ ] Random seeds fixed
- [ ] pixi environment added and `pixi.lock` updated
- [ ] Registry entry (+ `iteration_order`) with `pixi_env` and image
- [ ] Docker image builds and the method runs in `benchmark_run.sh`
