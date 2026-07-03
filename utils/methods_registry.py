"""
methods_registry.py — Read-only API for methods_registry.yaml.

Provides the canonical source of method metadata for all Python consumers
(benchmark_run_config.py, ranking_table.py, compute_metrics.py, etc.).

Usage:
    from utils.methods_registry import (
        load_registry,
        get_method,
        list_method_keys,
        get_entrypoint,
        get_execution,
        get_capabilities,
        get_capabilities_table,
        get_iteration_order,
        get_default_methods_json,
        export_execution_json,
        export_registry_json,
    )

    registry = load_registry()

    # Look up a single method
    method = get_method("FLeCS")
    print(method["entrypoint"]['script'])   # methods/flecs/flecs_train.py

    # Get capabilities as a dict (for ranking_table.py)
    caps = get_capabilities_table()         # {name: {grn_inference: True, ...}}

    # Export execution config as JSON (for BENCHMARK_EXECUTION_JSON)
    exec_json = export_execution_json()
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml as _yaml  # type: ignore


# ── Cache ────────────────────────────────────────────────────────────────────
_registry: Optional[Dict[str, Any]] = None
_registry_path: Optional[Path] = None


def _repo_root() -> Path:
    """Return the repository root (parent of utils/)."""
    return Path(__file__).resolve().parent.parent


def _default_registry_path() -> Path:
    return _repo_root() / "methods_registry.yaml"


def load_registry(path: Optional[str | Path] = None) -> Dict[str, Any]:
    """Load (and cache) the methods registry from YAML."""
    global _registry, _registry_path
    p = Path(path) if path else _default_registry_path()
    if _registry is not None and _registry_path == p:
        return _registry
    with open(p, "r", encoding="utf-8") as fh:
        _registry = _yaml.safe_load(fh)
    _registry_path = p
    if _registry is None:
        _registry = {}
    return _registry


def _methods_section() -> Dict[str, Any]:
    return load_registry().get("methods", {})


def _defaults_section() -> Dict[str, Any]:
    return load_registry().get("defaults", {})


# ── Lookup helpers ───────────────────────────────────────────────────────────

def list_method_keys() -> List[str]:
    """Return all registered method keys."""
    return list(_methods_section().keys())


def get_iteration_order() -> List[str]:
    """Return the ordered list of method keys for iteration."""
    order = load_registry().get("iteration_order", [])
    return order if order else list_method_keys()


def has_method(key: str) -> bool:
    """Check whether a method key is registered."""
    return key in _methods_section()


def get_method(key: str) -> Dict[str, Any]:
    """Return the full metadata dict for a method, or raise KeyError."""
    methods = _methods_section()
    if key not in methods:
        raise KeyError(f"Method '{key}' not found in registry. Known: {list(methods)}")
    return methods[key]


# ── Entrypoint ───────────────────────────────────────────────────────────────

def get_entrypoint(key: str) -> Dict[str, Any]:
    """Return {type, script, pass_output_mode, pass_ko_genes,
                supports_perturbation_training, extra_args}."""
    return get_method(key).get("entrypoint", {})


# ── Execution ─────────────────────────────────────────────────────────────────

def get_execution(key: str) -> Dict[str, Any]:
    """Return {default_runner, conda_env, container: {image, env_file, env_name}}."""
    return get_method(key).get("execution", {})


def export_execution_json() -> str:
    """
    Export execution config as a JSON string suitable for
    BENCHMARK_EXECUTION_JSON (consumed by single_run.sh).
    """
    defaults = _defaults_section().get("execution", {})
    methods = _methods_section()

    result: Dict[str, Any] = {
        "default_runner": defaults.get("runner", "conda"),
        "container_engine": defaults.get("container_engine", "docker"),
        "conda_envs": {"default": defaults.get("conda_env", "cardamom_env")},
        "container_images": {"default": defaults.get("container_image", "benchmark/cardamom:latest")},
        "method_runners": {},
    }

    for key, m in methods.items():
        ex = m.get("execution", {})
        result["conda_envs"][key] = ex.get("conda_env", defaults.get("conda_env", "cardamom_env"))
        container = ex.get("container", {})
        result["container_images"][key] = container.get("image", defaults.get("container_image", ""))
        runner = ex.get("default_runner")
        if runner and runner != defaults.get("runner", "conda"):
            result["method_runners"][key] = runner

    return json.dumps(result, separators=(",", ":"))


# ── Capabilities ─────────────────────────────────────────────────────────────

def get_capabilities(key: str) -> Dict[str, bool]:
    """Return {grn_inference, trajectory_reconstruction, perturbation_training}."""
    return get_method(key).get("capabilities", {})


def get_capabilities_table() -> Dict[str, Dict[str, bool]]:
    """
    Return the capabilities dict in the format expected by ranking_table.py:
        {method_key: {grn_inference: bool, trajectory_reconstruction: bool,
                       perturbation_training: bool}}
    """
    table: Dict[str, Dict[str, bool]] = {}
    for key in list_method_keys():
        table[key] = get_capabilities(key)
    return table


# ── Default methods JSON (for benchmark_run.sh fallback) ─────────────────────

def get_default_methods_json() -> str:
    """
    Return a JSON string mapping every registered method to 1,
    suitable as the default --run_methods value.
    """
    enabled = {key: 1 for key in list_method_keys()}
    return json.dumps(enabled, separators=(",", ":"))


# ── Full registry as JSON (for shell scripts that want the whole thing) ──────

def export_registry_json() -> str:
    """Export the entire registry as a compact JSON string."""
    return json.dumps(load_registry(), separators=(",", ":"), default=str)


# ── Container image list (for build systems) ─────────────────────────────────

def list_container_images() -> List[Dict[str, str]]:
    """
    Return [{target, env_file, env_name, image}, ...] for every method
    that has a distinct container image, deduplicated by image tag.
    """
    seen: set = set()
    images: List[Dict[str, str]] = []
    for key in get_iteration_order():
        if not has_method(key):
            continue
        ex = get_execution(key)
        container = ex.get("container", {})
        image = container.get("image", "")
        if not image or image in seen:
            continue
        seen.add(image)
        images.append({
            "target": key.lower().replace("-", "").replace("_", ""),
            "env_file": container.get("env_file", ""),
            "env_name": container.get("env_name", ""),
            "image": image.rstrip(":latest"),
        })
    return images


# ── CLI (for shell scripts to call) ──────────────────────────────────────────

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python utils/methods_registry.py <command> [args...]", file=sys.stderr)
        print("Commands: keys, execution-json, capabilities-json, default-methods-json,", file=sys.stderr)
        print("          registry-json, container-images, entrypoint <key>", file=sys.stderr)
        sys.exit(1)

    cmd = sys.argv[1]

    if cmd == "keys":
        print("\n".join(get_iteration_order()))
    elif cmd == "execution-json":
        print(export_execution_json())
    elif cmd == "capabilities-json":
        print(json.dumps(get_capabilities_table(), separators=(",", ":")))
    elif cmd == "default-methods-json":
        print(get_default_methods_json())
    elif cmd == "registry-json":
        print(export_registry_json())
    elif cmd == "container-images":
        print(json.dumps(list_container_images(), separators=(",", ":")))
    elif cmd == "entrypoint":
        if len(sys.argv) < 3:
            print("Usage: python utils/methods_registry.py entrypoint <method_key>", file=sys.stderr)
            sys.exit(1)
        key = sys.argv[2]
        ep = get_entrypoint(key)
        # Output as pipe-delimited fields matching method_spec() format:
        # type|script|pass_out|pass_ko|supports_perturbation|extra_arg1 extra_arg2...
        extra = " ".join(ep.get("extra_args", []))
        print("|".join([
            ep.get("type", ""),
            ep.get("script", ""),
            str(ep.get("pass_output_mode", False)).lower(),
            str(ep.get("pass_ko_genes", False)).lower(),
            str(ep.get("supports_perturbation_training", False)).lower(),
            extra,
        ]))
    elif cmd == "supports-perturbation":
        if len(sys.argv) < 3:
            sys.exit(1)
        key = sys.argv[2]
        ep = get_entrypoint(key)
        print("1" if ep.get("supports_perturbation_training", False) else "0")
    elif cmd == "conda-env":
        if len(sys.argv) < 3:
            sys.exit(1)
        key = sys.argv[2]
        print(get_execution(key).get("conda_env", "cardamom_env"))
    elif cmd == "runner-mode":
        if len(sys.argv) < 3:
            sys.exit(1)
        key = sys.argv[2]
        defaults = _defaults_section().get("execution", {})
        print(get_execution(key).get("default_runner", defaults.get("runner", "conda")))
    elif cmd == "container-image":
        if len(sys.argv) < 3:
            sys.exit(1)
        key = sys.argv[2]
        container = get_execution(key).get("container", {})
        print(container.get("image", ""))
    else:
        print(f"Unknown command: {cmd}", file=sys.stderr)
        sys.exit(1)
