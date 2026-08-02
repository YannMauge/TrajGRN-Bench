"""dyngen simulator backend for TrajGRN-Bench.

Uses dyngen's gold standard (deterministic ODE trajectory) to produce
time-varying expression data matching the benchmark GRN from common_config.py.

Python wrappers:
  - dyngen_simulate_custom.py    : WT simulation
  - dyngen_simulate_custom_ko.py : KO simulation (forwards to custom.py)

Both call dyngen_simulate.R via Rscript subprocess.
"""
