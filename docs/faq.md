# FAQ & Troubleshooting

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

## Troubleshooting

- Check existing [GitHub Issues](https://github.com/YannMauge/TrajGRN-Bench/issues)
- Open a new issue with:
    - Full error message
    - Your config file (sanitized)
    - OS and Docker version
    - Steps to reproduce
