# Contributing to TrajGRN-Bench

We welcome contributions! This guide covers how to contribute new methods,
datasets, metrics, and improvements.

## Code of Conduct

Please be respectful and constructive. We follow the
[Contributor Covenant](https://www.contributor-covenant.org/).

## Ways to Contribute

<div class="grid cards" markdown>

-   :material-puzzle-plus:{ .lg .middle } **Add a Method**

    ---

    Integrate a new GRN/trajectory inference method. See [Adding a Method](development_add_method.md).

-   :material-database-plus:{ .lg .middle } **Add a Dataset**

    ---

    Contribute a new benchmarking dataset. See [Adding a Dataset](development_add_dataset.md).

-   :material-chart-box-plus:{ .lg .middle } **Add a Metric**

    ---

    Propose a new evaluation metric for GRN or trajectory quality.

-   :material-bug:{ .lg .middle } **Report a Bug**

    ---

    Open an issue with detailed reproduction steps.

-   :material-book-edit:{ .lg .middle } **Improve Docs**

    ---

    Fix typos, add examples, or clarify explanations.

</div>

## Development Workflow

### 1. Fork and Clone

```bash
git clone https://github.com/YOUR_USERNAME/TrajGRN-Bench.git
cd TrajGRN-Bench
git remote add upstream https://github.com/YOUR_ORG/TrajGRN-Bench.git
```

### 2. Create a Branch

```bash
git checkout -b feature/my-new-method
```

### 3. Make Your Changes

Follow the relevant guide:

- [Adding a Method](development_add_method.md)
- [Adding a Dataset](development_add_dataset.md)
- [Methods Registry](development_registry.md)

### 4. Test

```bash
# Run the full benchmark with your changes
bash benchmark_run.sh --config configs/benchmark.example.yaml

# Build and serve documentation locally
pip install mkdocs-material mkdocstrings
mkdocs serve
```

### 5. Submit a Pull Request

1. Push your branch
2. Open a PR against `main`
3. Describe your changes and why they're needed
4. Link any related issues

## Pull Request Checklist

- [ ] Method/dataset follows I/O specification
- [ ] Entry added to `methods_registry.yaml` (if applicable)
- [ ] Conda environment file provided (if new method)
- [ ] Docker image builds successfully
- [ ] Benchmark runs complete without errors
- [ ] Documentation updated (if applicable)
- [ ] Tests pass (if applicable)

## Coding Conventions

### Python

- Follow [PEP 8](https://pep8.org/)
- Use type hints where practical
- Document functions with Google-style docstrings
- Use `pathlib.Path` for file paths

### Shell Scripts

- Use `#!/bin/bash` shebang
- Quote all variable expansions: `"$var"`
- Use `[[ ]]` for conditionals
- Run through [shellcheck](https://www.shellcheck.net/)

### YAML

- Use 2-space indentation
- Quote strings that contain special characters
- Include comments for non-obvious settings

## Documentation

We use **MkDocs Material**. To contribute to docs:

```bash
pip install mkdocs-material "mkdocstrings[python]"
cd TrajGRN-Bench
mkdocs serve
```

Visit `http://localhost:8000` to preview.

## Adding a Metric

To add a new evaluation metric:

1. Implement the metric function in `post_analysis/metrics_grn.py` or
   `post_analysis/metrics_cells.py`
2. Register it in `post_analysis/compute_metrics.py`
3. Update `ranking_table.py` to display the new metric
4. Add documentation in this guide

## Questions?

Open a [GitHub Discussion](https://github.com/YOUR_ORG/TrajGRN-Bench/discussions)
or issue.
