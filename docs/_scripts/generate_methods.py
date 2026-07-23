#!/usr/bin/env python3
"""
generate_methods.py — MkDocs gen-files plugin script.

Auto-generates method documentation from methods_registry.yaml:
  - methods_capabilities.md  — capability matrix (embedded in methods.md)
  - methods_descriptions.md  — full method descriptions by category
  - methods_glance_table.md  — compact table for index.md
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path


REGISTRY_PATH = Path("methods_registry.yaml")
OUTPUT_CAPABILITIES = Path("methods_capabilities.md")
OUTPUT_DESCRIPTIONS = Path("methods_descriptions.md")
OUTPUT_GLANCE = Path("methods_glance_table.md")

# ── helpers ──────────────────────────────────────────────────────────────────

def _badge(value: bool) -> str:
    """Material icon badge for capability matrix."""
    return ":material-check:{ .green }" if value else ":material-close:{ .red }"


def _emoji(value: bool) -> str:
    """Unicode emoji for glance table."""
    return "\u2705" if value else "\u2014"


def _load_registry() -> dict:
    import yaml
    with open(REGISTRY_PATH) as f:
        return yaml.safe_load(f)


# ── categorisation ───────────────────────────────────────────────────────────

def _categorise(methods: dict, iteration_order: list[str]) -> dict[str, list[tuple[str, dict]]]:
    """Group methods by capability profile.

    Returns
    -------
    dict with keys ``"joint"``, ``"grn_only"``, ``"traj_only"``.
    """
    cats: dict[str, list[tuple[str, dict]]] = {
        "joint": [],
        "grn_only": [],
        "traj_only": [],
    }
    for key in iteration_order:
        m = methods.get(key)
        if m is None:
            continue
        caps = m.get("capabilities", {})
        grn = caps.get("grn_inference", False)
        traj = caps.get("trajectory_reconstruction", False)
        if grn and traj:
            cats["joint"].append((key, m))
        elif grn:
            cats["grn_only"].append((key, m))
        elif traj:
            cats["traj_only"].append((key, m))
    return cats


# ── variant detection (methods sharing the same GitHub URL) ──────────────────

def _build_variant_map(methods: dict, iteration_order: list[str]) -> dict[str, list[str]]:
    """Return ``{primary_key: [variant_keys...]}`` for methods sharing a GitHub URL.

    The first key in *iteration_order* for a given URL becomes the primary; all
    later keys are considered variants.
    """
    by_url: dict[str, list[str]] = defaultdict(list)
    for key in iteration_order:
        m = methods.get(key)
        if m is None:
            continue
        url = (m.get("github") or "").strip()
        if url:
            by_url[url].append(key)

    variant_map: dict[str, list[str]] = {}
    for keys in by_url.values():
        if len(keys) > 1:
            variant_map[keys[0]] = keys[1:]
    return variant_map


# ── per-method markdown ──────────────────────────────────────────────────────

def _method_entry(
    key: str,
    m: dict,
    variant_map: dict[str, list[str]],
    methods: dict,
) -> str:
    """Render a single method's description block."""
    display = m.get("display_name", key)
    description = m.get("description", "").strip()
    github = (m.get("github") or "").strip()
    publication = (m.get("publication") or "").strip()
    caps = m.get("capabilities", {})

    lines: list[str] = [f"### {display}", ""]

    # GitHub button or local-script path
    if github:
        lines.append(
            f"[:fontawesome-brands-github: Repository]({github})"
            "{ .md-button .md-button--primary }"
        )
    else:
        script = m.get("entrypoint", {}).get("script", "")
        if script:
            lines.append(f"Local script: `{script}`")
    lines.append("")

    # Description
    lines.append(description)
    lines.append("")

    # Publication
    if publication:
        lines.append(f"- **Publication:** {publication}")

    # Capability badges
    grn = _badge(caps.get("grn_inference", False))
    traj = _badge(caps.get("trajectory_reconstruction", False))
    pert = _badge(caps.get("perturbation_training", False))
    lines.append(f"- **GRN:** {grn} **Trajectory:** {traj} **Perturbation:** {pert}")

    # Variant note
    variants = variant_map.get(key, [])
    if variants:
        variant_displays = [
            f"`{methods.get(v, {}).get('display_name', v)}`"
            for v in variants
        ]
        lines.append(f"- **Variant{'s' if len(variants) > 1 else ''}:** {', '.join(variant_displays)}")

    lines.append("")
    return "\n".join(lines)


# ── generators ───────────────────────────────────────────────────────────────

def generate_descriptions(registry: dict) -> str:
    """Generate full method descriptions organised by capability category."""
    methods = registry.get("methods", {})
    iteration_order: list[str] = registry.get("iteration_order", list(methods.keys()))
    categories = _categorise(methods, iteration_order)
    variant_map = _build_variant_map(methods, iteration_order)

    sections: list[str] = []

    def _append_category(title: str, intro: str, items: list[tuple[str, dict]]) -> None:
        if not items:
            return
        sections.append(f"## {title}")
        sections.append("")
        sections.append(intro)
        sections.append("")
        for key, m in items:
            sections.append(_method_entry(key, m, variant_map, methods))
            sections.append("---")
            sections.append("")

    _append_category(
        "Joint GRN + Trajectory Methods",
        "These methods jointly infer both the gene regulatory network and "
        "cellular trajectories from time-course scRNA-seq data.",
        categories["joint"],
    )
    _append_category(
        "GRN-Only Methods",
        "",
        categories["grn_only"],
    )
    _append_category(
        "Trajectory-Only Methods",
        "",
        categories["traj_only"],
    )

    # Strip trailing "---\n\n"
    result = "\n".join(sections)
    while result.endswith("---\n\n") or result.endswith("---\n"):
        result = result.rstrip()
        if result.endswith("---"):
            result = result[:-3].rstrip()

    return result


def generate_capabilities_table(registry: dict) -> str:
    """Generate the detailed capability-matrix table."""
    methods = registry.get("methods", {})
    iteration_order: list[str] = registry.get("iteration_order", list(methods.keys()))

    lines = [
        "## Capability Matrix",
        "",
        "| Method | GRN Inference | Trajectory | Perturbation Training |",
        "|--------|:------------:|:----------:|:---------------------:|",
    ]

    for key in iteration_order:
        if key not in methods:
            continue
        m = methods[key]
        caps = m.get("capabilities", {})
        display = m.get("display_name", key)
        grn = _badge(caps.get("grn_inference", False))
        traj = _badge(caps.get("trajectory_reconstruction", False))
        pert = _badge(caps.get("perturbation_training", False))
        lines.append(f"| **{display}** | {grn} | {traj} | {pert} |")

    lines.append("")
    return "\n".join(lines)


def generate_glance_table(registry: dict) -> str:
    """Generate the compact 'Methods at a Glance' table for index.md."""
    methods = registry.get("methods", {})
    iteration_order: list[str] = registry.get("iteration_order", list(methods.keys()))

    lines = [
        "| Method | GRN | Trajectory | Perturbations |",
        "|--------|:---:|:----------:|:-------------:|",
    ]

    for key in iteration_order:
        if key not in methods:
            continue
        m = methods[key]
        caps = m.get("capabilities", {})
        display = m.get("display_name", key)
        grn = _emoji(caps.get("grn_inference", False))
        traj = _emoji(caps.get("trajectory_reconstruction", False))
        pert = _emoji(caps.get("perturbation_training", False))
        lines.append(f"| **{display}** | {grn} | {traj} | {pert} |")

    lines.append("")
    return "\n".join(lines)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    """Called by the mkdocs-gen-files plugin during ``mkdocs build``.

    Writes generated markdown to the real ``docs/`` directory so that
    ``pymdownx.snippets`` (``--8<--`` syntax) can resolve includes.
    """
    _write_generated_files()


def _write_generated_files() -> None:
    """Write all generated markdown files to the docs/ directory."""
    registry = _load_registry()
    docs_dir = Path("docs")
    if not docs_dir.is_dir():
        docs_dir = Path(".")

    content_map = {
        OUTPUT_CAPABILITIES: generate_capabilities_table(registry),
        OUTPUT_DESCRIPTIONS: generate_descriptions(registry),
        OUTPUT_GLANCE: generate_glance_table(registry),
    }
    for path, content in content_map.items():
        (docs_dir / path).write_text(content)


if __name__ == "__main__":
    import sys

    if "--standalone" in sys.argv:
        _write_generated_files()
        print(f"Wrote {OUTPUT_CAPABILITIES}, {OUTPUT_DESCRIPTIONS}, {OUTPUT_GLANCE}")
    else:
        _write_generated_files()
