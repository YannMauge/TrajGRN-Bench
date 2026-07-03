#!/usr/bin/env python3
"""
generate_results.py — MkDocs gen-files plugin script.

Scans docs/results_data/ for benchmark result directories and generates
Markdown tables and Plotly interactive charts embedded in results.md.

Each result directory should contain:
  - metadata.yaml   (run label, description, date, config summary, methods list)
  - grn_metrics.csv (GRN metric scores)
  - cell_metrics.csv (cell/trajectory metric scores)
  - gene_metrics.csv (optional gene-level scores)
  - ranking_table.png (optional visual ranking)
  - notes.md (optional free-form notes)
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import mkdocs_gen_files


RESULTS_DATA_DIR = Path("docs/results_data")
RESULTS_OUTPUT = Path("results_generated.md")


def find_result_dirs(base: Path) -> List[Path]:
    """Find all result directories containing metadata.yaml."""
    if not base.exists():
        return []
    result_dirs = []
    for item in sorted(base.iterdir()):
        if item.is_dir() and (item / "metadata.yaml").exists():
            result_dirs.append(item)
    return result_dirs


def read_metadata(result_dir: Path) -> Optional[Dict[str, Any]]:
    """Read metadata.yaml from a result directory."""
    meta_file = result_dir / "metadata.yaml"
    if not meta_file.exists():
        return None
    try:
        import yaml
        with open(meta_file) as f:
            return yaml.safe_load(f)
    except ImportError:
        return {"label": result_dir.name, "description": "No metadata available"}


def read_csv_safe(csv_path: Path) -> Optional[str]:
    """Read a CSV file and return it as a Markdown table, or None."""
    if not csv_path.exists():
        return None
    try:
        import pandas as pd
        df = pd.read_csv(csv_path)
        return df.to_markdown(index=False)
    except ImportError:
        with open(csv_path) as f:
            return f"```csv\n{f.read()}\n```"


def generate_plotly_chart(csv_path: Path, title: str, metric_cols: List[str]) -> str:
    """Generate a Plotly bar chart as an HTML div from a CSV file."""
    if not csv_path.exists():
        return ""
    try:
        import pandas as pd
    except ImportError:
        return f"[CSV data available]({csv_path})"

    df = pd.read_csv(csv_path)

    # Find method column
    method_col = None
    for col in ["method", "Method", "name", "Name"]:
        if col in df.columns:
            method_col = col
            break

    if method_col is None:
        return f"```csv\n{df.head().to_csv(index=False)}\n```"

    # Find metric columns that exist
    available_metrics = [c for c in metric_cols if c in df.columns]
    if not available_metrics:
        return f"```csv\n{df.head().to_csv(index=False)}\n```"

    methods = df[method_col].tolist()

    # Build Plotly traces
    traces = []
    for metric in available_metrics:
        traces.append({
            "x": methods,
            "y": df[metric].tolist(),
            "name": metric,
            "type": "bar",
        })

    chart_data = json.dumps({
        "data": traces,
        "layout": {
            "title": title,
            "barmode": "group",
            "xaxis": {"title": "Method"},
            "yaxis": {"title": "Score"},
            "height": 400,
        },
    })

    return f"""<div id="chart_{title.replace(' ', '_').replace('/', '_')}"></div>
<script>
(function() {{
    var data = {chart_data};
    if (typeof Plotly !== 'undefined') {{
        Plotly.newPlot(
            document.getElementById('chart_{title.replace(' ', '_').replace('/', '_')}'),
            data.data,
            data.layout,
            {{responsive: true}}
        );
    }} else {{
        document.getElementById('chart_{title.replace(' ', '_').replace('/', '_')}').innerHTML =
            '<p><em>Plotly not loaded. Install via CDN or view as table below.</em></p>';
    }}
}})();
</script>"""


def generate_results_page(result_dirs: List[Path]) -> str:
    """Generate the full results Markdown page content."""
    lines = []

    if not result_dirs:
        lines.append(":fontawesome-solid-circle-info: **No results data yet.**")
        lines.append("")
        lines.append("Run a benchmark and export results to populate this page:")
        lines.append("```bash")
        lines.append("python docs/_scripts/export_results.py \\")
        lines.append("  --metrics-dir benchmark/outputs_metrics/<your_run> \\")
        lines.append("  --output-dir docs/results_data/<your_run>")
        lines.append("```")
        return "\n".join(lines)

    lines.append(f"## Available Results ({len(result_dirs)} runs)")
    lines.append("")

    for i, result_dir in enumerate(result_dirs):
        meta = read_metadata(result_dir)
        label = meta.get("label", result_dir.name) if meta else result_dir.name
        desc = meta.get("description", "") if meta else ""

        lines.append(f"### {label}")
        lines.append("")

        if desc:
            lines.append(f"{desc}")
            lines.append("")

        # Metadata table
        if meta:
            lines.append("| Field | Value |")
            lines.append("|-------|-------|")
            if "date" in meta:
                lines.append(f"| Date | {meta['date']} |")
            if "config" in meta:
                for key, val in meta["config"].items():
                    lines.append(f"| Config: `{key}` | `{val}` |")
            if "methods" in meta:
                lines.append(f"| Methods | {', '.join(meta['methods'])} |")
            lines.append("")

        # GRN metrics
        grn_csv = result_dir / "grn_metrics.csv"
        if grn_csv.exists():
            lines.append("#### GRN Metrics")
            lines.append("")
            grn_table = read_csv_safe(grn_csv)
            if grn_table:
                lines.append(grn_table)
                lines.append("")

            # Plotly chart
            chart = generate_plotly_chart(
                grn_csv,
                f"{label} — GRN Metrics",
                ["AUROC", "AUPRC", "AUPRC_signed", "Precision@K"],
            )
            if chart:
                lines.append(chart)
                lines.append("")

        # Cell metrics
        cell_csv = result_dir / "cell_metrics.csv"
        if cell_csv.exists():
            lines.append("#### Trajectory Metrics")
            lines.append("")
            cell_table = read_csv_safe(cell_csv)
            if cell_table:
                lines.append(cell_table)
                lines.append("")

            chart = generate_plotly_chart(
                cell_csv,
                f"{label} — Trajectory Metrics",
                ["Wasserstein_dist", "MMD", "R2"],
            )
            if chart:
                lines.append(chart)
                lines.append("")

        # Gene metrics
        gene_csv = result_dir / "gene_metrics.csv"
        if gene_csv.exists():
            lines.append("#### Gene-Level Metrics")
            lines.append("")
            gene_table = read_csv_safe(gene_csv)
            if gene_table:
                lines.append(gene_table)
                lines.append("")

        # Ranking image
        ranking_png = result_dir / "ranking_table.png"
        if ranking_png.exists():
            lines.append(f"![Ranking Table]({ranking_png})")
            lines.append("")

        # Notes
        notes_md = result_dir / "notes.md"
        if notes_md.exists():
            with open(notes_md) as f:
                lines.append(f.read())
            lines.append("")

        if i < len(result_dirs) - 1:
            lines.append("---")
            lines.append("")

    return "\n".join(lines)


def main():
    result_dirs = find_result_dirs(RESULTS_DATA_DIR)
    content = generate_results_page(result_dirs)

    with mkdocs_gen_files.open(RESULTS_OUTPUT, "w") as f:
        f.write(content)


if __name__ == "__main__":
    main()
