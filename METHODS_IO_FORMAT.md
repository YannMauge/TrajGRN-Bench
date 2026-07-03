# Methods Input/Output (I/O) Format Specification

This document specifies the standardized input/output formats for inference methods and simulators in the Benchmark project. All methods should follow these conventions for consistency and interoperability.

## Table of Contents

1. [Input Data Formats](#input-data-formats)
2. [Output Data Formats](#output-data-formats)
3. [File Organization](#file-organization)
4. [Command-Line Interface Convention](#command-line-interface-convention)
5. [Step-by-Step Guide for New Method Authors](#step-by-step-guide-for-new-method-authors)
   - [Overview of Files You'll Touch](#overview-of-files-youll-touch)
   - [Step 1: Write Your Method Entrypoint Script](#step-1-write-your-method-entrypoint-script)
   - [Step 2: Create a Conda Environment File](#step-2-create-a-conda-environment-file)
   - [Step 3: Register Your Method in methods_registry.yaml](#step-3-register-your-method-in-methods_registryyaml)
   - [Step 4: Local Testing (Conda)](#step-4-local-testing-conda)
   - [Step 5: Build and Test the Docker Image](#step-5-build-and-test-the-docker-image)
   - [Step 6: Verify Metrics Are Computed](#step-6-verify-metrics-are-computed)
   - [Quick Checklist](#quick-checklist)
6. [Common Data Structures](#common-data-structures)
7. [Examples](#examples)

---

## Input Data Formats

### Primary Gene Expression Data

**Format:** H5AD (HDF5-based AnnData format)

**File Location:** `benchmark/data/data_1.h5ad`, `benchmark/data/data_2.h5ad`, etc.

**Structure (AnnData object):**
```
adata
├── .X                          # Gene expression count matrix [n_cells × n_genes]
│                               # Dense array or sparse matrix
├── .obs (DataFrame)            # Cell metadata
│   ├── 'timepoint' or 'time'   # (Required) Timepoint index for each cell [0, 1, 2, ...]
│   ├── 'dataset_id'            # (Optional) Perturbation condition: 'WT', 'KO_GENE1', etc.
│   └── [other metadata]        # Any additional cell-level annotations
├── .var (DataFrame)            # Gene metadata
│   ├── index                   # Gene names (matching network matrix column/row names)
│   ├── 'd0'                    # (Optional) Basal degradation rate coefficient
│   ├── 'd1'                    # (Optional) Stimulus-dependent degradation rate
│   └── [other metadata]        # Gene-level annotations
├── .uns (dict)                 # Unstructured metadata
│   ├── 'simulation'            # (Optional) Indicates if data is simulated (True/False/None)
│   └── [other metadata]        # Method-specific parameters
└── .varm, .obsm, .obsp         # (Optional) Embeddings and neighbors
```

**Key Requirements:**
- `adata.X` must be convertible to dense array: `adata.X.toarray()` or already dense
- Time key must be present in `.obs` (name varies: 'timepoint', 'time', etc.)
- Gene names in `.var.index` must match network file gene lists
- Sparse matrices should be explicitly converted to dense arrays

### Timepoint Information

**Files:** `benchmark/data/train_tps.npy`, `benchmark/data/test_tps.npy`

**Format:** NumPy array, shape `(n_timepoints, 2)`

**Content:**
```python
# Each row: [timepoint_id, timepoint_value]
# Example:
# [[0, 0.0],   # Training timepoint 0 at time 0.0
#  [1, 1.0],   # Training timepoint 1 at time 1.0
#  [2, 2.0]]   # Training timepoint 2 at time 2.0
```

**Usage:**
- Filter training data: `adata[adata.obs['timepoint'].isin(train_tps[:, 1])]`
- Defines which timepoints are used for model fitting vs. evaluation

### Degradation Rates (Optional)

**Files:** `benchmark/data/halflife/degradation_rates.txt` (or similar)

**Format:** Tab-delimited text file

**Content:**
```
# Gene degradation rates, tab-separated
# Typically shape: [2 × n_genes] where:
# Row 0: Basal degradation rate (d0)
# Row 1: Stimulus-dependent rate (d1)

0.25	0.30	0.28	...	0.22
0.05	0.08	0.06	...	0.04
```

### Subsample Information (Optional)

**File:** `benchmark/data/subsample_train_ids.npy`

**Format:** NumPy array of cell identifiers

**Content:**
```python
# Cell IDs to use for training (subset of total cells)
# Example: ['CELL_0001', 'CELL_0003', 'CELL_0015', ...]
subsample_ids = np.load('benchmark/data/subsample_train_ids.npy', allow_pickle=True)
```

---

## Output Data Formats

### Learned Gene Regulatory Network (GRN)

**Format:** NumPy array (`.npy` binary or `.txt` text)

**File Location:** `benchmark/outputs_methods/<method_name>/`

**Naming Convention:**
```
<data_id>_GRN.npy                    # WT (wild-type) GRN
<data_id>_ko_<GENE_NAME>_GRN.npy     # KO (knockout) GRN
```

**Structure:**
```python
# Shape: (n_genes + n_stimuli, n_genes + n_stimuli)
# Typical: (n_genes + 1, n_genes + 1) with 1 stimulus

# First n_stimuli rows/columns: stimulus genes (typically 'Stimulus' or 'Stimulus_0')
# Remaining rows/columns: target genes in same order as adata.var_names

# Content: Interaction strengths or adjacency
# - Continuous values: regulatory strength (e.g., -0.5 to 1.0)
# - Binary: 0 (no interaction) or 1 (interaction exists)
# - Gene regulatory network matrix[regulator, target] = strength
```

**Example (8 genes + 1 stimulus):**
```python
# Rows/cols: [Stimulus, Gene1, Gene2, ..., Gene8]
grn = np.load('data_1_GRN.npy')  # Shape: (9, 9)
# grn[0, :] = outgoing edges FROM Stimulus
# grn[:, 1] = incoming edges TO Gene1
# grn[2, 3] = Stimulus on Gene2 affecting Gene3
```

### Predicted Gene Expression Data

**Format:** H5AD (AnnData)

**File Location:** `benchmark/outputs_methods/<method_name>/`

**Naming Convention:**
```
<data_id>_adata.h5ad                        # WT predictions
<data_id>_ko_<GENE_NAME>_adata.h5ad         # KO predictions
```

**Structure:**
```
pred_adata
├── .X                          # Predicted expression counts [n_cells × n_genes]
├── .obs
│   ├── 'time' or 'timepoint'   # Timepoint index
│   ├── 'dataset_id'            # 'WT' or KO gene name
│   └── [other metadata]
├── .var                        # Gene names (matches input)
├── .uns
│   ├── 'runtime'               # Execution time in seconds
│   └── [other metadata]
└── [embeddings, etc.]
```

**Key Requirements:**
- Predictions typically from last training timepoint forward
- Include original training data cells with predictions
- Remove stimulus gene column (if present in input) from output
- Store runtime as `pred_adata.uns['runtime']`


## File Organization

### Benchmark Data Directory Structure

```
benchmark/
├── data/
│   ├── data_1.h5ad                          # Primary expression data
│   ├── data_2.h5ad
│   ├── train_tps.npy                        # Timepoints for training
│   ├── test_tps.npy                         # Timepoints for evaluation
│   ├── subsample_train_ids.npy              # (Optional) subset of cells
│   ├── halflife/
│   │   └── degradation_rates.txt            # Degradation rate info
│   └── True/                                # Ground truth (if synthetic data)
│       ├── inter_1.npy
│       ├── inter_2.npy
│       └── ...
├── outputs_methods/
│   ├── test_run_harissa_future_full_test/
│   │   ├── data_1_GRN.npy
│   │   ├── data_1_adata.h5ad
│   │   ├── data_1_ko_GENE1_GRN.npy
│   │   ├── data_1_ko_GENE1_adata.h5ad
│   │   └── ...
│   ├── test_run_boolode_1_future_full_test/
│   └── ...
└── outputs_metrics/
    ├── test_run_harissa_future_full_test/
    │   ├── data_1_grn_mse.npy
    │   ├── data_1_expression_mae.npy
    │   └── ...
    └── ...
```

---

## Command-Line Interface Convention

All methods must accept input/output via command-line arguments. Two invocation styles are supported:

### Python Entrypoints (`getopt` or `argparse`)

```bash
python method_script.py -i <input_h5ad> -o <output_folder> [optional_args]
```

**Required Arguments:**
- `-i, --input`: Path to input H5AD file
- `-o, --outputfolder`: Path to output directory (will be created if needed)

**Optional Arguments (passed automatically by the framework when applicable):**
- `-u, --output_mode`: Output type (`full_test`, `full_train`, `full_full`, `no_traj`)
- `-p, --perturbation_training`: String `"true"` or `"false"` indicating whether KO perturbation groups are available for training
- `-k, --ko_genes`: Comma-separated list of genes for knockout prediction (`"none"`, `"all"`, or e.g. `"GENE1,GENE2"`)

**Method-specific extra args** (defined in `methods_registry.yaml` → `entrypoint.extra_args`):
- Example: `--path_strategy knn` (FLeCS), `--n_pca_dims 50` (custom)
- These are passed directly after the standard args

You may use either `getopt` (e.g., `methods/pearson/inference_pearson.py`) or `argparse` (e.g., `methods/flecs/flecs_train.py`).

### Bash Entrypoints (positional arguments)

```bash
bash method_script.sh <input_h5ad> <output_folder> [output_mode] [perturbation_training] [ko_genes]
```

**Positional Arguments (in order):**
1. `<input_h5ad>` — Path to input H5AD file
2. `<output_folder>` — Path to output directory
3. `[output_mode]` — (Optional) `full_test`, `full_train`, `full_full`, or `no_traj`
4. `[perturbation_training]` — (Optional) `"true"` or `"false"`
5. `[ko_genes]` — (Optional) `"none"`, `"all"`, or comma-separated gene list

Example: `methods/CardamomOT/run_carda_benchmark.sh`

---

## Step-by-Step Guide for New Method Authors

This guide walks you through integrating a new inference method into the benchmark. You will create a method script, a conda environment, register it in the registry, test locally, and build a Docker image.

### Overview of Files You'll Touch

| File | Purpose |
|------|---------|
| `methods/<your_method>/` | Your method code (create a new directory) |
| `environments/<your_env>.yml` | Conda environment specification |
| `methods_registry.yaml` | Register your method's metadata, entrypoint, and capabilities |
| `configs/benchmark.schema.json` | Add your method name to the allowed list (optional, for config validation) |

No other files need to be modified. The build system, shell scripts, and Python tools all derive their configuration from `methods_registry.yaml`.

---

### Step 1: Write Your Method Entrypoint Script

Your script is the **single entrypoint** the benchmark framework calls. It must implement the I/O contract described above.

#### 1a. Choose Python or Bash

- **Python** (`type: python` in registry): Receives named arguments `-i`, `-o`, `-u`, `-p`, `-k` plus any `extra_args`. Use `argparse` or `getopt`.
- **Bash** (`type: bash` in registry): Receives positional arguments: `<input> <outputfolder> [output_mode] [perturbation_training] [ko_genes]`.

#### 1b. Create the script

Create `methods/<your_method>/run_<your_method>.py` (or `.sh`). A minimal Python skeleton:

```python
#!/usr/bin/env python3
"""
run_<your_method>.py — Benchmark entrypoint for <YourMethod>.

Usage:
    python run_<your_method>.py -i <input_h5ad> -o <output_folder>
                                [-u <output_mode>] [-p <perturbation_training>]
                                [-k <ko_genes>]
"""
import sys, os, time
import argparse
import numpy as np
import scanpy as sc
import anndata as ad
import scipy

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputfile", required=True)
    parser.add_argument("-o", "--outputfolder", required=True)
    parser.add_argument("-u", "--output_mode", default="full_test")
    parser.add_argument("-p", "--perturbation_training", default="false")
    parser.add_argument("-k", "--ko_genes", default="none")
    return parser.parse_args(argv)

def main(argv):
    t_start = time.time()
    args = parse_args(argv)

    # 1. Load input data
    adata = ad.read_h5ad(args.inputfile)
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()

    # 2. Filter to training timepoints
    data_dir = os.path.dirname(args.inputfile)
    train_tps = np.load(os.path.join(data_dir, "train_tps.npy"))[:, 1]
    t_key = "timepoint" if "timepoint" in adata.obs.columns else "time"
    adata_train = adata[adata.obs[t_key].isin(train_tps)].copy()

    # 3. Run your method (infer GRN, predict expression, etc.)
    grn = your_method_infer_grn(adata_train)

    # 4. Save outputs
    os.makedirs(args.outputfolder, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.inputfile))[0]
    np.save(os.path.join(args.outputfolder, f"{base}_GRN"), grn)

    # 5. (Optional) Save predicted expression as AnnData
    if args.output_mode != "no_traj":
        pred_adata = your_method_predict(adata_train, grn)
        pred_adata.uns["runtime"] = time.time() - t_start
        pred_adata.write(os.path.join(args.outputfolder, f"{base}_adata.h5ad"))

if __name__ == "__main__":
    main(sys.argv[1:])
```

#### 1c. Mandatory Outputs

| Your method supports… | Must produce per input file |
|------------------------|-----------------------------|
| GRN inference (`grn_inference: true`) | `<data_id>_GRN.npy` — shape `(G, G)` NumPy array |
| Trajectory prediction (`trajectory_reconstruction: true`) | `<data_id>_adata.h5ad` — AnnData with `.X`, `.obs['time']`, `.uns['runtime']` |
| Perturbation training (`perturbation_training: true`) | When `-p true`, also produce: `<data_id>_ko_<GENE>_GRN.npy` and `<data_id>_ko_<GENE>_adata.h5ad` |

**GRN matrix conventions:**
- Shape: `(n_genes, n_genes)` (no stimulus row/column in output)
- Gene order must match `adata.var_names`
- Values: continuous (signed interaction strength) or binary (0/1 adjacency)

**Predicted AnnData conventions:**
- `.X`: dense array, shape `(n_predicted_cells, n_genes)`
- `.obs`: must contain a `'time'` column with timepoint values
- `.var`: gene names matching input
- `.uns['runtime']`: total wall-clock time in seconds (float)
- Typically includes predictions from the last training timepoint forward

---

### Step 2: Create a Conda Environment File

Create `environments/<your_env>.yml` with all dependencies your method needs:

```yaml
name: your_method_env
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - numpy
  - scipy
  - pandas
  - anndata
  - scanpy
  - pip
  - pip:
    - your-method-package
```

**Rules:**
- The `name:` line will be stripped during Docker builds; the actual env name comes from `methods_registry.yaml`.
- Install and test locally: `conda env create -f environments/<your_env>.yml`
- Place the file in `environments/` alongside existing env files.

---

### Step 3: Register Your Method in `methods_registry.yaml`

Add a new entry under `methods:` in `methods_registry.yaml`. Follow the existing pattern:

```yaml
  YourMethod:
    display_name: Your Method
    description: >
      Brief description of what your method does (1-2 sentences).
    github: https://github.com/yourname/yourmethod
    publication: "Venue (Year) DOI:..."
    # ── Entrypoint ──────────────────────────────────────────────────
    entrypoint:
      type: python                              # python | bash
      script: methods/your_method/run_your_method.py
      pass_output_mode: true                    # does it accept -u?
      pass_ko_genes: true                       # does it accept -k?
      supports_perturbation_training: false     # does it accept -p and produce KO outputs?
      extra_args: []                            # extra CLI flags, e.g. ["--myflag", "value"]
    # ── Execution environment ───────────────────────────────────────
    execution:
      default_runner: conda                     # conda | docker | apptainer
      conda_env: your_method_env                # must match your env file's intent
      container:
        image: benchmark/yourmethod:latest
        env_file: environments/your_env.yml
        env_name: your_method_env
    # ── Capabilities ────────────────────────────────────────────────
    capabilities:
      grn_inference: true
      trajectory_reconstruction: true
      perturbation_training: false
```

**Field reference:**

| Field | Description |
|-------|-------------|
| `entrypoint.type` | `python` (named flags) or `bash` (positional args) |
| `entrypoint.script` | Path from repo root to your script |
| `entrypoint.pass_output_mode` | If `true`, `-u <output_mode>` is appended (python) or passed as 3rd positional arg (bash) |
| `entrypoint.pass_ko_genes` | If `true`, `-k <ko_genes>` is appended (python) or passed as 5th positional arg (bash) |
| `entrypoint.supports_perturbation_training` | If `true`, the framework calls your script twice: once with `-p false`, once with `-p true` |
| `entrypoint.extra_args` | Additional CLI flags passed verbatim after standard args |
| `execution.default_runner` | `conda` (local), `docker` (container), or `apptainer` |
| `execution.conda_env` | Conda environment name used when `default_runner: conda` |
| `execution.container.image` | Docker image tag (used by build system and `single_run.sh`) |
| `execution.container.env_file` | Path to your conda environment YAML (used by Docker build) |
| `execution.container.env_name` | Environment name inside the container (often same as `conda_env`) |
| `capabilities.grn_inference` | Does your method produce a `_GRN.npy`? |
| `capabilities.trajectory_reconstruction` | Does your method produce a `_adata.h5ad`? |
| `capabilities.perturbation_training` | Can your method leverage KO perturbation data? |

Also append your method key to `iteration_order:` at the bottom of the file:

```yaml
iteration_order:
  - FLeCS
  # ... existing methods ...
  - WaddingtonOT
  - YourMethod                           # <-- add here
```

**(Optional)** To enable config-file validation, add your method name to `configs/benchmark.schema.json` under `properties.benchmark.properties.run_methods.propertyNames.enum`.

---

### Step 4: Local Testing (Conda)

Before building a container, verify your method runs correctly with conda.

#### 4a. Validate your registry entry

```bash
# List all registered methods — yours should appear
python utils/methods_registry.py keys

# Inspect your method's entrypoint spec
python utils/methods_registry.py entrypoint YourMethod
# Output: python|methods/your_method/run_your_method.py|true|true|false|

# Check execution config
python utils/methods_registry.py execution-json | python -m json.tool
```

#### 4b. Test your method directly

```bash
# Activate your conda environment
conda activate your_method_env

# Run on a single data file
python methods/your_method/run_your_method.py \
    -i benchmark/data/data_1.h5ad \
    -o /tmp/test_output/ \
    -u full_test \
    -k none \
    -p false

# Verify outputs
ls /tmp/test_output/
# Should contain: data_1_GRN.npy  (and data_1_adata.h5ad if applicable)
```

#### 4c. Run via the benchmark framework (single data file)

```bash
# Export execution config so single_run.sh knows about your method
export BENCHMARK_EXECUTION_JSON=$(python utils/methods_registry.py execution-json)

# Create a run_methods JSON that only enables your method
run_methods='{"YourMethod":1}'

# Dry-run on one data file
bash single_run.sh \
    benchmark/data \
    /tmp/test_results/ \
    "$run_methods" \
    full_test \
    none \
    false
```

#### 4d. Run the full benchmark with your method

```bash
# Enable your method in configs/benchmark.example.yaml (or your own config)
# Then run:
bash benchmark_run.sh --config configs/benchmark.example.yaml
```

---

### Step 5: Build and Test the Docker Image

Once your method works locally with conda, build its container image.

#### 5a. Build the image

```bash
# Option A: Using Make (recommended)
make -f containers/Makefile yourmethod

# Option B: Using the shell script
bash containers/build_all.sh YourMethod

# Option C: Manually
docker build \
    -f containers/docker/Dockerfile \
    --build-arg ENV_FILE=environments/your_env.yml \
    --build-arg ENV_NAME=your_method_env \
    -t benchmark/yourmethod:latest \
    .
```

#### 5b. Test the container locally

```bash
# Run your method inside the container
docker run --rm \
    -v "$(pwd):/work" \
    -w /work \
    benchmark/yourmethod:latest \
    python methods/your_method/run_your_method.py \
        -i benchmark/data/data_1.h5ad \
        -o /tmp/test_output/ \
        -u full_test

# Check the output was written
ls /tmp/test_output/
```

#### 5c. Run the benchmark with Docker runner

Set your method's `execution.default_runner` to `docker` in `methods_registry.yaml`, or override it in your benchmark config:

```yaml
execution:
  method_runners:
    YourMethod: docker
```

Then run the benchmark as usual. The framework will pull your image and execute your method inside the container.

#### 5d. (Optional) Apptainer / Singularity

```bash
# Convert the Docker image to a SIF file
apptainer build yourmethod.sif docker-daemon://benchmark/yourmethod:latest

# Test
apptainer exec --bind "$(pwd):/work" yourmethod.sif \
    python methods/your_method/run_your_method.py \
        -i benchmark/data/data_1.h5ad \
        -o /tmp/test_output/
```

Then set `method_runners.YourMethod: apptainer` and point `container_images.YourMethod` at your `.sif` file.

---

### Step 6: Verify Metrics Are Computed

After a successful benchmark run, check that metrics were generated for your method:

```bash
ls benchmark/outputs_metrics/<run_name>/
# Look for: data_1_grn_mse.npy, data_1_expression_mae.npy, etc.

# Run post-analysis to generate visualizations and ranking tables
python post_analysis/compute_metrics.py  # (check actual script name)
```

---

### Quick Checklist

- [ ] Method script created in `methods/<your_method>/`
- [ ] Script accepts `-i`, `-o` (and optionally `-u`, `-p`, `-k`)
- [ ] Script produces `<data_id>_GRN.npy` and/or `<data_id>_adata.h5ad`
- [ ] Conda environment file created in `environments/`
- [ ] Entry added to `methods_registry.yaml` with correct capabilities
- [ ] Method name added to `iteration_order`
- [ ] Method name added to `configs/benchmark.schema.json` (optional)
- [ ] Local conda test passes (`single_run.sh` with only your method)
- [ ] Docker image builds successfully
- [ ] Docker-based test passes
- [ ] Full benchmark run produces metrics for your method

---

## Version History

- **v1.1** (2026-07-03): Added step-by-step guide for new method authors; clarified CLI conventions for Python vs Bash entrypoints.
- **v1.0** (2026-06-29): Initial standardized I/O format specification
