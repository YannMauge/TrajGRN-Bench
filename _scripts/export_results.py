#!/usr/bin/env python3
"""
export_results.py — Export benchmark run results to the docs data directory.

Usage:
    python docs/_scripts/export_results.py \\
        --metrics-dir benchmark/outputs_metrics/my_run_future_full_test \\
        --output-dir docs/results_data/my_run \\
        --label "My Benchmark Run" \\
        --description "Full benchmark on dataset 1"

This creates:
    docs/results_data/my_run/
    ├── metadata.yaml       # Run metadata
    ├── grn_metrics.csv     # Copy of GRN metrics
    ├── cell_metrics.csv    # Copy of cell metrics
    ├── gene_metrics.csv    # Copy of gene metrics (if exists)
    ├── ranking_table.png   # Copy of ranking chart (if exists)
    └── notes.md            # User-editable notes placeholder
"""

import argparse
import shutil
from datetime import date
from pathlib import Path
from typing import Optional


def export_results(
    metrics_dir: Path,
    output_dir: Path,
    label: str = "",
    description: str = "",
    methods: Optional[list] = None,
    config_summary: Optional[dict] = None,
) -> None:
    """Export benchmark metric files to the docs results data directory."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy metric CSV files
    for csv_name in ["grn_metrics.csv", "cell_metrics.csv", "gene_metrics.csv"]:
        src = metrics_dir / csv_name
        if src.exists():
            shutil.copy2(src, output_dir / csv_name)
            print(f"  ✓ Copied {csv_name}")
        else:
            print(f"  - {csv_name} not found (skipping)")

    # Copy ranking table image
    ranking_src = metrics_dir / "ranking_table.png"
    if ranking_src.exists():
        shutil.copy2(ranking_src, output_dir / "ranking_table.png")
        print(f"  ✓ Copied ranking_table.png")

    # Write metadata.yaml
    meta = {
        "label": label or metrics_dir.parent.name,
        "description": description,
        "date": str(date.today()),
    }
    if config_summary:
        meta["config"] = config_summary
    if methods:
        meta["methods"] = methods

    import yaml
    with open(output_dir / "metadata.yaml", "w") as f:
        yaml.dump(meta, f, default_flow_style=False, sort_keys=False)
    print(f"  ✓ Wrote metadata.yaml")

    # Create notes.md placeholder if it doesn't exist
    notes_path = output_dir / "notes.md"
    if not notes_path.exists():
        notes_path.write_text(
            f"# {label or output_dir.name}\n\n"
            f"*Add your notes and observations about this run here.*\n"
        )
        print(f"  ✓ Created notes.md placeholder")

    print(f"\nDone! Results exported to: {output_dir}")
    print(f"Rebuild docs with: mkdocs build")


def main():
    parser = argparse.ArgumentParser(
        description="Export benchmark results to docs data directory"
    )
    parser.add_argument(
        "--metrics-dir", required=True, type=Path,
        help="Path to metrics output directory (contains grn_metrics.csv, etc.)"
    )
    parser.add_argument(
        "--output-dir", required=True, type=Path,
        help="Output directory under docs/results_data/"
    )
    parser.add_argument(
        "--label", default="",
        help="Human-readable label for this run"
    )
    parser.add_argument(
        "--description", default="",
        help="Description of this run"
    )
    parser.add_argument(
        "--methods", nargs="*", default=None,
        help="List of methods used in this run"
    )
    parser.add_argument(
        "--config", type=Path, default=None,
        help="Path to benchmark config YAML (for extracting summary)"
    )

    args = parser.parse_args()

    config_summary = None
    if args.config and args.config.exists():
        import yaml
        with open(args.config) as f:
            cfg = yaml.safe_load(f)
        config_summary = {
            "train_data": cfg.get("train_data", "?"),
            "output_mode": cfg.get("output_mode", "?"),
            "simulator": cfg.get("simulation", {}).get("simulator_backend", "?"),
        }

    export_results(
        metrics_dir=args.metrics_dir,
        output_dir=args.output_dir,
        label=args.label,
        description=args.description,
        methods=args.methods,
        config_summary=config_summary,
    )


if __name__ == "__main__":
    main()
