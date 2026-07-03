# 🧬 TrajGRN-Bench

**Trajectory and Gene Regulatory Network Inference Benchmark**

---

TrajGRN-Bench is a comprehensive benchmark of mechanistic models for **joint gene
regulatory network (GRN) and single-cell RNA-seq trajectory inference**. We
evaluate and compare multiple joint-inference approaches alongside complementary
GRN and trajectory methods on simulated datasets.

## What We Provide

<div class="grid cards" markdown>

-   :material-chart-bell-curve:{ .lg .middle } **Standardized Metrics**

    ---

    Unified metrics for GRN inference (AUROC, AUPRC, etc.) and trajectory
    reconstruction (Wasserstein distance, EMD, etc.) — making fair comparisons
    across methods.

-   :material-cog-sync:{ .lg .middle } **Flexible Configuration**

    ---

    A single YAML config drives everything: which methods run, which datasets
    are used, what metrics to compute, and how results are stored.

-   :material-docker:{ .lg .middle } **Docker-Based Execution**

    ---

    Each method runs in its own container. No dependency conflicts. Reproducible
    across any machine with Docker installed.

-   :material-table-star:{ .lg .middle } **Interactive Rankings**

    ---

    Auto-generated ranking tables and interactive visualizations let you explore
    method performance across all metrics at a glance.

</div>

## Quick Start

```bash
# 1. Install the runner environment
conda env create -f environments/benchmark_runner.yml

# 2. Build method containers
make -f containers/Makefile -j4 all

# 3. Run the benchmark
bash benchmark_run.sh --config configs/benchmark.example.yaml
```

## Methods at a Glance

--8<-- "methods_glance_table.md"

[:material-book-open-page-variant: See full method details →](methods.md)

## Project Status

This benchmark is under active development. See the [GitHub
repository](https://github.com/YannMauge/TrajGRN-Bench) for the latest updates
and to contribute.

---

<div class="grid" markdown>

[:fontawesome-solid-rocket: Get Started](install.md){ .md-button .md-button--primary }
[:fontawesome-solid-book: User Guide](usage.md){ .md-button }

</div>
