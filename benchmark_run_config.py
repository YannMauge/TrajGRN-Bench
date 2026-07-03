#!/usr/bin/env python3
"""
benchmark_run_config.py - Config-driven benchmark runner.

Loads a YAML/JSON configuration file and invokes benchmark_run.sh
with the derived arguments, exporting execution settings for methods.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Dict, Optional


DEFAULTS = {
    "train_data": "future",
    "output_mode": "full_test",
    "simulation": "simul_ko",
    "global_run_name": "test_run_3",
    "future_start_tp": 10,
    "replicates_number": 2,
    "simul_ko_genes": "all",
    "ko_output_genes": "none",
    "restart_mode": "save",
    "perturbation_training": True,
    "simulator_backend": "harissa",
}


class ConfigError(RuntimeError):
    pass


def _load_yaml(path: Path) -> Dict[str, Any]:
    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise ConfigError(
            "PyYAML is required to load YAML configs. Install it in the runner env."
        ) from exc
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError("Config root must be a mapping.")
    return data


def _load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ConfigError("Config root must be a mapping.")
    return data


def load_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise ConfigError(f"Config file not found: {path}")
    suffix = path.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        return _load_yaml(path)
    if suffix == ".json":
        return _load_json(path)
    raise ConfigError("Config must be a .yaml, .yml, or .json file.")


# Canonical boolean-like values mapped to int (1/0 for bash JSON compatibility)
_BOOL_MAP: Dict[str, int] = {
    "1": 1, "true": 1, "yes": 1, "y": 1,
    "0": 0, "false": 0, "no": 0, "n": 0, "": 0,
}


def _normalize_bool(value: Any) -> int:
    if isinstance(value, bool):
        return 1 if value else 0
    if isinstance(value, (int, float)):
        return 1 if value else 0
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in _BOOL_MAP:
            return _BOOL_MAP[normalized]
    raise ConfigError(f"Invalid boolean-like value: {value!r}")


def _require(mapping: Dict[str, Any], key: str) -> Any:
    if key not in mapping:
        raise ConfigError(f"Missing required config key: {key}")
    return mapping[key]


def _get_section(config: Dict[str, Any], name: str) -> Dict[str, Any]:
    section = config.get(name, {})
    if section is None:
        return {}
    if not isinstance(section, dict):
        raise ConfigError(f"Config section '{name}' must be a mapping.")
    return section


def build_execution_json(config: Dict[str, Any]) -> str:
    """
    Build the BENCHMARK_EXECUTION_JSON string by merging the registry defaults
    with any user-provided overrides from the config's ``execution`` section.
    """
    try:
        from utils.methods_registry import export_execution_json as _registry_exec_json
        import json as _json
        base = _json.loads(_registry_exec_json())
    except Exception:
        base: Dict[str, Any] = {}

    execution = _get_section(config, "execution")

    # Merge user overrides (shallow merge for each sub-dict)
    for key in ("conda_envs", "container_images", "method_runners"):
        if key in execution and isinstance(execution[key], dict):
            base.setdefault(key, {})
            base[key].update(execution[key])

    for key in ("default_runner", "container_engine"):
        if key in execution and execution[key] is not None:
            base[key] = execution[key]

    # Also handle top-level conda_envs.default override
    if "runner_env" in execution:
        base.setdefault("conda_envs", {})
        base["conda_envs"]["default"] = execution["runner_env"]

    import json as _json
    return _json.dumps(base, separators=(",", ":"))


def build_args(config: Dict[str, Any]) -> Dict[str, Any]:
    benchmark = _get_section(config, "benchmark")
    simulation = _get_section(config, "simulation")
    paths = _get_section(config, "paths")
    execution = _get_section(config, "execution")
    simulation = _get_section(config, "simulation")
    paths = _get_section(config, "paths")
    execution = _get_section(config, "execution")

    train_data = benchmark.get("train_data", DEFAULTS["train_data"])
    output_mode = benchmark.get("output_mode", DEFAULTS["output_mode"])
    global_run_name = benchmark.get("global_run_name", DEFAULTS["global_run_name"])
    run_methods = _require(benchmark, "run_methods")

    if not isinstance(run_methods, dict):
        raise ConfigError("benchmark.run_methods must be a mapping.")

    run_methods_json = {
        name: _normalize_bool(value) for name, value in run_methods.items()
    }

    future_start_tp = benchmark.get("future_start_tp", DEFAULTS["future_start_tp"])
    ko_output_genes = benchmark.get("ko_output_genes", DEFAULTS["ko_output_genes"])
    restart_mode = benchmark.get("restart_mode", DEFAULTS["restart_mode"])
    if restart_mode not in ("save", "rerun"):
        raise ConfigError(f"restart_mode must be 'save' or 'rerun', got {restart_mode!r}")

    perturbation_training = _normalize_bool(
        benchmark.get("perturbation_training", DEFAULTS["perturbation_training"])
    )

    simulation_mode = simulation.get("mode", DEFAULTS["simulation"])
    replicates_number = simulation.get("replicates_number", DEFAULTS["replicates_number"])
    simul_ko_genes = simulation.get("simul_ko_genes", DEFAULTS["simul_ko_genes"])
    simulator_backend = simulation.get("simulator_backend", DEFAULTS["simulator_backend"])

    input_dir = paths.get("input_dir")
    adata_dir = paths.get("adata_dir")
    results_dir = paths.get("results_dir")

    runner_env = execution.get("runner_env")

    return {
        "train_data": train_data,
        "output_mode": output_mode,
        "simulation": simulation_mode,
        "global_run_name": global_run_name,
        "run_methods_json": json.dumps(run_methods_json, separators=(",", ":")),
        "future_start_tp": future_start_tp,
        "replicates_number": replicates_number,
        "simul_ko_genes": simul_ko_genes,
        "ko_output_genes": ko_output_genes,
        "restart_mode": restart_mode,
        "perturbation_training": perturbation_training,
        "simulator_backend": simulator_backend,
        "input_dir": input_dir,
        "adata_dir": adata_dir,
        "results_dir": results_dir,
        "execution_json": build_execution_json(config),
        "runner_env": runner_env,
    }


def run_benchmark(config_path: Path) -> None:
    config = load_config(config_path)

    # Validate against JSON Schema early to catch misconfigurations
    repo_root = Path(__file__).resolve().parent
    schema_path = repo_root / "configs" / "benchmark.schema.json"
    if schema_path.exists():
        import jsonschema
        with schema_path.open("r", encoding="utf-8") as schema_fh:
            schema = json.load(schema_fh)
        jsonschema.validate(config, schema)
    else:
        print("Warning: Schema file not found, skipping validation.", file=sys.stderr)

    args = build_args(config)

    command = [
        "bash",
        "benchmark_run.sh",
        args["train_data"],
        args["output_mode"],
        args["simulation"],
        args["global_run_name"],
        args["run_methods_json"],
        str(args["future_start_tp"]),
        str(args["replicates_number"]),
        str(args["simul_ko_genes"]),
        str(args["simulator_backend"]),
        "true" if args["perturbation_training"] else "false",
        str(args["ko_output_genes"]),
        str(args["restart_mode"]),
    ]

    env = os.environ.copy()
    env["BENCHMARK_EXECUTION_JSON"] = args["execution_json"]
    if args["runner_env"]:
        env["BENCHMARK_RUNNER_ENV"] = str(args["runner_env"])
    if args["input_dir"]:
        env["BENCHMARK_INPUT_DIR"] = str(args["input_dir"])
    if args["adata_dir"]:
        env["BENCHMARK_ADATA_DIR"] = str(args["adata_dir"])
    if args["results_dir"]:
        env["BENCHMARK_RESULTS_DIR"] = str(args["results_dir"])

    subprocess.run(command, check=True, env=env, cwd=repo_root)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run benchmark from a YAML/JSON config file.")
    parser.add_argument("-c", "--config", required=True, help="Path to YAML/JSON config file.")
    args = parser.parse_args(argv)

    try:
        run_benchmark(Path(args.config))
    except ConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 2
    except subprocess.CalledProcessError as exc:
        return exc.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
