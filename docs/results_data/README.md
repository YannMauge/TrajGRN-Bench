# Results Data

Place benchmark result exports here. Each subdirectory represents one benchmark run.

## Structure

```
results_data/
├── my_run/
│   ├── metadata.yaml       # Run label, description, date, config
│   ├── grn_metrics.csv     # GRN metric scores
│   ├── cell_metrics.csv    # Cell/trajectory scores
│   ├── gene_metrics.csv    # (optional) Gene-level scores
│   ├── ranking_table.png   # (optional) Visual ranking
│   └── notes.md            # Free-form notes
```

## Adding Results

```bash
python docs/_scripts/export_results.py \
  --metrics-dir benchmark/outputs_metrics/<your_run> \
  --output-dir docs/results_data/<your_run> \
  --label "My Run" \
  --description "Description of this run"
```

Results are auto-rendered in the docs when you rebuild.
