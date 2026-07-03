#!/usr/bin/env python3
"""
benchmark_run.py - Combined benchmark runner for all analysis modes.

Supports four training data modes (--train_data):
  - full:          train on all timepoints
  - future:        train on only the first x timepoints
  - leave-one-out: leave one intermediate timepoint out as test; run all combinations
  - subsample_full: all timepoints but only 66% of cells at each timepoint

Supports four output simulation modes (--output_mode):
  - full_full:  simulate all train and test timepoints from the first timepoint
                (when train_data==subsample_full, t0 is included in the test set)
  - full_train: simulate all training timepoints from the first timepoint
  - full_test:  simulate all test timepoints from the last timepoint before test
  - no_traj:    do not output simulation

Usage:
    python benchmark_run.py --train_data <train_data> --output_mode <output_mode>
                            -a <adata_dir> -n <name> -r <results_dir>
                            -m <run_methods_json> [-t <future_start_tp>]

Arguments:
    --train_data          Training data mode: full, future, leave-one-out, subsample_full
    --output_mode         Output simulation mode: full_full, full_train, full_test, no_traj
    -a, --adata_dir       Directory containing the adata files
    -n, --name            Name for the benchmark run
    -r, --results_dir     Directory to store results
    -m, --run_methods     JSON object specifying which methods to run
    -t, --future_start_tp Timepoint index to start future prediction (for future train_data only)
"""

import argparse
import json
import logging
import os
import subprocess
import sys

import anndata as ad
import numpy as np
from typing import List, Optional, Tuple
import warnings
for module in ["anndata", "scipy"]:
    warnings.filterwarnings("ignore", module=module)


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def parse_bool_arg(value: str) -> bool:
    """Parse CLI boolean-like values."""
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(
        f"Invalid boolean value: {value!r}. Use true/false."
    )


class BenchmarkRunnerError(Exception):
    """Custom exception for BenchmarkRunner errors."""

    pass


class BenchmarkRunner:
    """
    Combined benchmark runner supporting all training data and output simulation modes.

    Handles argument parsing, data loading, train/test split
    preparation, and execution of single_run.sh for each mode.
    """

    VALID_TRAIN_DATA = ["full", "future", "leave-one-out", "subsample_full"]
    VALID_OUTPUT_MODES = ["full_full", "full_train", "full_test", "no_traj"]

    def __init__(self) -> None:
        self.train_data: Optional[str] = None
        self.output_mode: Optional[str] = None
        self.adata_dir: Optional[str] = None
        self.global_run_name: Optional[str] = None
        self.results_dir: Optional[str] = None
        self.run_methods_json: Optional[str] = None
        self.future_start_tp: Optional[int] = None
        self.perturbation_training: bool = True
        self.restart_mode: str = "save"
        self.ko_output_genes: str = "none"
        self.adata: Optional[ad.AnnData] = None
        self.unique_tps: Optional[np.ndarray] = None
        self.n_timepoints: int = 0

    # -------------------------------------------------------------------------
    # Argument parsing
    # -------------------------------------------------------------------------

    def parse_arguments(self, argv: List[str]) -> None:
        """Parse command-line arguments."""
        parser = argparse.ArgumentParser(
            description="Combined benchmark runner supporting multiple training data "
                        "and output simulation modes."
        )
        parser.add_argument(
            "--train_data",
            required=True,
            choices=self.VALID_TRAIN_DATA,
            help="Training data mode: full, future, leave-one-out, subsample_full",
        )
        parser.add_argument(
            "--output_mode",
            required=True,
            choices=self.VALID_OUTPUT_MODES,
            help="Output simulation mode: full_full, full_train, full_test, no_traj",
        )
        parser.add_argument(
            "-a", "--adata_dir",
            required=True,
            help="Directory containing the adata files",
        )
        parser.add_argument(
            "-n", "--name",
            required=True,
            dest="global_run_name",
            help="Name for the benchmark run",
        )
        parser.add_argument(
            "-r", "--results_dir",
            required=True,
            help="Directory to store results",
        )
        parser.add_argument(
            "-m", "--run_methods",
            required=True,
            dest="run_methods_json",
            help="JSON object specifying which methods to run",
        )
        parser.add_argument(
            "-t", "--future_start_tp",
            type=int,
            default=None,
            help="Timepoint index to start future prediction (required for future train_data)",
        )
        parser.add_argument(
            "--ko_output_genes",
            default="none",
            help='KO targets for output simulation ("none", "all", or comma-separated gene names)',
        )
        parser.add_argument(
            "--perturbation_training",
            type=parse_bool_arg,
            default=True,
            help="Whether to run perturbation-training variants of supported methods (default: true)",
        )
        parser.add_argument(
            "--restart_mode",
            default="save",
            choices=["save", "rerun"],
            help="save: skip completed methods; rerun: re-run all methods (default: save)",
        )
        parser.add_argument(
            "--ko_genes",
            default=None,
            help=argparse.SUPPRESS,
        )

        args = parser.parse_args(argv)
        self.train_data = args.train_data
        self.output_mode = args.output_mode
        self.adata_dir = args.adata_dir
        self.global_run_name = args.global_run_name
        self.results_dir = args.results_dir
        self.run_methods_json = args.run_methods_json
        self.future_start_tp = args.future_start_tp
        self.restart_mode = args.restart_mode
        self.perturbation_training = args.perturbation_training
        self.ko_output_genes = str(args.ko_output_genes)
        if args.ko_genes is not None:
            if self.ko_output_genes == "none":
                logger.warning("--ko_genes is deprecated; use --ko_output_genes instead")
                self.ko_output_genes = str(args.ko_genes)
            elif self.ko_output_genes != str(args.ko_genes):
                logger.warning(
                    "--ko_genes is deprecated and ignored because --ko_output_genes was also provided."
                )

        logger.info(
            f"Arguments parsed: train_data={self.train_data}, output_mode={self.output_mode}, "
            f"adata_dir={self.adata_dir}, name={self.global_run_name}, "
            f"perturbation_training={self.perturbation_training}, "
            f"ko_output_genes={self.ko_output_genes}, "
            f"restart_mode={self.restart_mode}"
        )

    # -------------------------------------------------------------------------
    # Data loading
    # -------------------------------------------------------------------------

    @property
    def adata_file(self) -> str:
        """Path to the primary adata file."""
        return os.path.join(self.adata_dir, "data_1.h5ad")

    def load_adata(self) -> None:
        """Load and validate the AnnData object."""
        try:
            self.adata = ad.read_h5ad(self.adata_file)
        except Exception as e:
            raise BenchmarkRunnerError(f"Failed to load AnnData file: {e}")

    def get_timepoints(self) -> None:
        """Extract unique timepoints from the AnnData object."""
        if "timepoint" not in self.adata.obs.columns:
            raise BenchmarkRunnerError(
                "'timepoint' column not found in adata.obs"
            )
        self.unique_tps = np.unique(self.adata.obs["timepoint"])
        self.n_timepoints = len(self.unique_tps)
        logger.info(
            f"Extracted {self.n_timepoints} unique timepoints: {self.unique_tps}"
        )

    # -------------------------------------------------------------------------
    # Train/test split preparation (train_data-specific)
    # -------------------------------------------------------------------------

    def _splits_full(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Full mode: train on all timepoints, no test set.
        """
        if self.n_timepoints <= 0:
            raise BenchmarkRunnerError(
                f"n_timepoints must be positive, got {self.n_timepoints}"
            )
        train_idx = np.arange(self.n_timepoints)
        train_tps = np.column_stack((train_idx, self.unique_tps[train_idx]))
        test_tps = np.empty((0, 2), dtype=int)
        logger.info(
            f"Full splits: train={train_tps.shape}, test={test_tps.shape}"
        )
        return train_tps, test_tps

    def _splits_future(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Future mode: train on timepoints before future_start_tp,
        test on timepoints from future_start_tp onward.
        """
        if self.future_start_tp >= self.n_timepoints:
            raise BenchmarkRunnerError(
                f"future_start_tp ({self.future_start_tp}) must be less than "
                f"n_timepoints ({self.n_timepoints})"
            )
        train_idx = np.arange(self.future_start_tp)
        test_idx = np.arange(self.future_start_tp, self.n_timepoints)
        train_tps = np.column_stack((train_idx, self.unique_tps[train_idx]))
        test_tps = np.column_stack((test_idx, self.unique_tps[test_idx]))
        logger.info(
            f"Future splits: train={train_tps.shape}, test={test_tps.shape}"
        )
        return train_tps, test_tps

    def _splits_leave_one_out(self, test_tp_idx: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Leave-one-out mode: train on all timepoints except test_tp_idx,
        test on test_tp_idx. Only intermediate timepoints (not first or last)
        are valid test candidates.
        """
        if test_tp_idx <= 0 or test_tp_idx >= self.n_timepoints - 1:
            raise BenchmarkRunnerError(
                f"test_tp_idx must be in the range [1, {self.n_timepoints - 2}] "
                f"(got {test_tp_idx}); the first (index 0) and last "
                f"(index {self.n_timepoints - 1}) timepoints cannot be held out"
            )
        train_idx = np.concatenate(
            [np.arange(test_tp_idx), np.arange(test_tp_idx + 1, self.n_timepoints)]
        )
        train_tps = np.column_stack((train_idx, self.unique_tps[train_idx]))
        test_tps = np.column_stack(([test_tp_idx], [self.unique_tps[test_tp_idx]]))
        return train_tps, test_tps

    def _splits_subsample_full(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Subsample-full mode: train on all timepoints using only 66% of cells
        at each timepoint. The first timepoint (t0) is included in the test set
        so that full-trajectory simulations can be evaluated against the complete
        t0 population.
        """
        if self.n_timepoints <= 0:
            raise BenchmarkRunnerError(
                f"n_timepoints must be positive, got {self.n_timepoints}"
            )
        train_idx = np.arange(self.n_timepoints)
        train_tps = np.column_stack((train_idx, self.unique_tps[train_idx]))
        # t0 is kept as the test reference for simulation evaluation
        test_tps = np.column_stack(([0], [self.unique_tps[0]]))
        logger.info(
            f"Subsample-full splits: train={train_tps.shape}, test={test_tps.shape}"
        )
        return train_tps, test_tps

    # -------------------------------------------------------------------------
    # I/O helpers
    # -------------------------------------------------------------------------

    def save_timepoints(
        self,
        train_tps: np.ndarray,
        test_tps: np.ndarray,
        output_dir: Optional[str] = None,
    ) -> None:
        """Save train/test timepoint arrays to *output_dir* (default: adata_dir)."""
        output_dir = output_dir or self.adata_dir
        try:
            np.save(os.path.join(output_dir, "train_tps.npy"), train_tps)
            np.save(os.path.join(output_dir, "test_tps.npy"), test_tps)
            logger.info(f"Saved timepoints to {output_dir}")
        except Exception as e:
            raise BenchmarkRunnerError(f"Failed to save timepoints: {e}")

    def save_subsample_ids(self, output_dir: Optional[str] = None) -> None:
        """
        For subsample_full mode: save the IDs of the 66% of cells selected
        at each timepoint for training.  Methods that support subsampling will
        load this file and restrict their training set accordingly.
        """
        output_dir = output_dir or self.adata_dir
        rng = np.random.default_rng(42)
        selected_ids: List[str] = []
        for tp in self.unique_tps:
            cell_mask = self.adata.obs["timepoint"] == tp
            cells_at_tp = self.adata.obs_names[cell_mask]
            n_subsample = max(1, int(len(cells_at_tp) * 0.66))
            chosen = rng.choice(cells_at_tp, n_subsample, replace=False)
            selected_ids.extend(chosen.tolist())
        try:
            np.save(
                os.path.join(output_dir, "subsample_train_ids.npy"),
                np.array(selected_ids),
            )
            logger.info(
                f"Saved subsample_train_ids ({len(selected_ids)} cells) to {output_dir}"
            )
        except Exception as e:
            raise BenchmarkRunnerError(f"Failed to save subsample IDs: {e}")

    def create_results_directory(self, suffix: str) -> str:
        """Create and return a results sub-directory named <global_run_name><suffix>."""
        results_dir_run = os.path.join(
            self.results_dir, f"{self.global_run_name}{suffix}"
        )
        try:
            os.makedirs(results_dir_run, exist_ok=True)
            logger.info(f"Created results directory: {results_dir_run}")
        except Exception as e:
            raise BenchmarkRunnerError(
                f"Failed to create results directory: {e}"
            )
        return results_dir_run

    # -------------------------------------------------------------------------
    # Subprocess execution
    # -------------------------------------------------------------------------

    def execute_single_run(
        self,
        results_dir_run: str,
        ko_output_genes: str = "none",
    ) -> None:
        """
        Execute single_run.sh with the given results directory.

        Args:
            results_dir_run:       Path to the results directory for this run.
            ko_output_genes:       Passed to single_run.sh for KO output simulation.
        """
        command = [
            "bash",
            "./single_run.sh",
            self.adata_dir,
            results_dir_run,
            self.run_methods_json,
            self.output_mode,
            ko_output_genes,
            str(self.perturbation_training).lower(),
            self.restart_mode,
        ]

        logger.info(f"Running: {' '.join(command)}")

        if not os.environ.get("BENCHMARK_EXECUTION_JSON"):
            logger.warning(
                "BENCHMARK_EXECUTION_JSON is not set; methods will use default conda env fallbacks"
            )

        try:
            subprocess.run(command, check=True, env=os.environ.copy())
        except KeyboardInterrupt:
            raise BenchmarkRunnerError("Interrupted by user")
        except BenchmarkRunnerError:
            raise
        except Exception as e:
            raise BenchmarkRunnerError(f"Failed to execute single_run.sh: {e}")

    # -------------------------------------------------------------------------
    # Train-data runners
    # -------------------------------------------------------------------------

    def run_full(self) -> None:
        """Execute benchmark with all timepoints as training data."""
        train_tps, test_tps = self._splits_full()
        self.save_timepoints(train_tps, test_tps)
        suffix = f"_full_{self.output_mode}"
        results_dir_run = self.create_results_directory(suffix)
        self.execute_single_run(
            results_dir_run,
            ko_output_genes=self.ko_output_genes,
        )
        self.save_timepoints(train_tps, test_tps, results_dir_run)
        logger.info("Full training analysis completed successfully")

    def run_future(self) -> None:
        """Execute benchmark with future prediction training split."""
        if self.future_start_tp is None:
            raise BenchmarkRunnerError(
                "future_start_tp is required for train_data=future"
            )
        train_tps, test_tps = self._splits_future()
        self.save_timepoints(train_tps, test_tps)
        suffix = f"_future_{self.output_mode}"
        results_dir_run = self.create_results_directory(suffix)
        self.execute_single_run(
            results_dir_run,
            ko_output_genes=self.ko_output_genes,
        )
        self.save_timepoints(train_tps, test_tps, results_dir_run)
        logger.info("Future prediction analysis completed successfully")

    def run_leave_one_out(self) -> None:
        """
        Execute leave-one-out analysis (one run per held-out intermediate timepoint).

        Iterates over all timepoints strictly between the first and last, leaving
        each one out as the test set in turn.
        """
        for i in range(1, self.n_timepoints - 1):
            train_tps, test_tps = self._splits_leave_one_out(i)
            self.save_timepoints(train_tps, test_tps)
            suffix = f"_leave-one-out_{self.output_mode}/run_{i}"
            results_dir_run = self.create_results_directory(suffix)
            self.execute_single_run(
                results_dir_run,
                ko_output_genes=self.ko_output_genes,
            )
        logger.info("Leave-one-out analysis completed successfully")

    def run_subsample_full(self) -> None:
        """
        Execute benchmark with all timepoints but only 66% of cells per timepoint
        as training data.
        """
        train_tps, test_tps = self._splits_subsample_full()
        self.save_timepoints(train_tps, test_tps)
        self.save_subsample_ids()
        suffix = f"_subsample_full_{self.output_mode}"
        results_dir_run = self.create_results_directory(suffix)
        self.execute_single_run(
            results_dir_run,
            ko_output_genes=self.ko_output_genes,
        )
        self.save_timepoints(train_tps, test_tps, results_dir_run)
        self.save_subsample_ids(results_dir_run)
        logger.info("Subsample-full analysis completed successfully")

    # -------------------------------------------------------------------------
    # Main entry point
    # -------------------------------------------------------------------------

    def run(self) -> None:
        """Dispatch to the appropriate train-data runner after shared setup."""
        self.load_adata()
        self.get_timepoints()
        if self.perturbation_training:
            if "dataset_id" not in self.adata.obs.columns:
                raise BenchmarkRunnerError(
                    "perturbation_training=true requires adata.obs['dataset_id']"
                )
            n_groups = self.adata.obs["dataset_id"].nunique(dropna=True)
            if n_groups < 2:
                groups = sorted(self.adata.obs["dataset_id"].dropna().unique().tolist())
                raise BenchmarkRunnerError(
                    "perturbation_training=true requires at least two dataset_id groups; "
                    f"found {n_groups}: {groups}"
                )

        # Clean up any stale subsample_train_ids.npy from a previous run so
        # that methods only apply subsampling when it was explicitly requested.
        # The file is recreated below by run_subsample_full() when needed.
        subsample_ids_path = os.path.join(self.adata_dir, "subsample_train_ids.npy")
        if os.path.exists(subsample_ids_path):
            os.remove(subsample_ids_path)
            logger.info(f"Removed stale subsample_train_ids.npy from {self.adata_dir}")

        dispatch = {
            "full": self.run_full,
            "future": self.run_future,
            "leave-one-out": self.run_leave_one_out,
            "subsample_full": self.run_subsample_full,
        }
        dispatch[self.train_data]()


def main(argv: List[str]) -> None:
    """Main entry point for benchmark_run.py."""
    try:
        runner = BenchmarkRunner()
        runner.parse_arguments(argv)
        runner.run()
    except BenchmarkRunnerError as e:
        logger.error(str(e))
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])
