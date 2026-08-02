#!/usr/bin/env python3
"""
dyngen_simulate_custom_ko.py — Thin wrapper for dyngen KO simulation.

Delegates to dyngen_simulate_custom.py, which handles both WT-only and KO
modes via the --ko_genes flag.  benchmark_run.sh already passes -k, so we
simply forward all arguments.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the dyngen module is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from simulator.dyngen.dyngen_simulate_custom import main

if __name__ == "__main__":
    main(sys.argv[1:])
