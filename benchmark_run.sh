#!/bin/bash
# =============================================================================
# benchmark_run.sh - Main benchmark runner script
# =============================================================================
# Usage: ./benchmark_run.sh --config <config.yaml>          (recommended)
#        ./benchmark_run.sh <train_data> <output_mode> <simulation> <global_run_name>
#            <run_methods_json> <future_start_tp> <replicates_number> <simul_ko_genes>
#            <simulator_backend> <perturbation_training> <ko_output_genes> <restart_mode>
#
# train_data: Defines the training data given to methods. One of:
#   - full          : train on all timepoints
#   - future        : train on only the first x timepoints
#   - leave-one-out : leave one intermediate timepoint out as test; run all combinations
#   - subsample_full: all timepoints but only 66% of cells at each timepoint
#
# output_mode: Defines the output simulated set. One of:
#   - full_full  : simulate all train and test timepoints from the first timepoint
#                  (when train_data==subsample_full, t0 is also in the test set)
#   - full_train : simulate all training timepoints from the first timepoint
#   - full_test  : simulate all test timepoints from the last timepoint before test
#   - no_traj    : do not output simulation
#
# simulation: "simul_replicates" for replicates, "simul_ko" for knockout, "false" for existing data
# simulator_backend: "harissa" (default), "boolode", "sergio", or "dyngen"
# global_run_name: Name for the benchmark run
# run_methods_json: JSON object specifying which methods to run
# future_start_tp: Timepoint to start future prediction (for future train_data only)
# replicates_number: Number of simulation runs/replicates when simulating data
# ko_output_genes: KO target list for output simulation ("none", "all", or comma-separated genes; default: none)
# simul_ko_genes: training KOs to simulate and give to methods ("none", "all" or comma-separated; default: all)
# perturbation_training: Whether to run perturbation-training variants of supported methods ("true" or "false"; default: false)
# restart_mode: "save" to skip methods with .done markers, "rerun" to re-run all (default: save)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Orchestration tool for the runner environment (simulators, splits, metrics):
#   pixi (default when available): pixi env "${BENCHMARK_RUNNER_ENV:-runner}"
#   conda (legacy):                 conda env "${BENCHMARK_RUNNER_ENV:-benchmark_runner}"
runner_tool="${BENCHMARK_RUNNER_TOOL:-}"
if [ -z "$runner_tool" ]; then
    if command -v pixi >/dev/null 2>&1 && [ -f "$SCRIPT_DIR/pixi.toml" ]; then runner_tool=pixi; else runner_tool=conda; fi
fi
if [ "$runner_tool" = "pixi" ]; then
    runner_env="${BENCHMARK_RUNNER_ENV:-runner}"
else
    runner_env="${BENCHMARK_RUNNER_ENV:-benchmark_runner}"
fi

# Run a command inside the runner environment.
run_in_runner() {
    if [ "$runner_tool" = "pixi" ]; then
        (
            eval "$(pixi shell-hook --manifest-path "$SCRIPT_DIR/pixi.toml" -e "$runner_env" -s bash)" || exit 1
            exec "$@"
        )
    else
        conda run -n "$runner_env" --live-stream "$@"
    fi
}

# Config shortcut (the config parser needs PyYAML / jsonschema: runner env).
# The config's execution.runner_tool only applies after parsing: to parse the
# config with conda, export BENCHMARK_RUNNER_TOOL=conda (and BENCHMARK_RUNNER_ENV).
if [[ "${1:-}" == "--config" || "${1:-}" == "-c" ]]; then
    if [ -z "${2:-}" ]; then
        echo "Usage: $0 --config <config.yaml>" >&2
        exit 1
    fi
    run_in_runner python benchmark_run_config.py --config "$2"
    exit $?
fi

# Set default paths
input_dir="${BENCHMARK_INPUT_DIR:-${SCRIPT_DIR}/simulator/custom_network/Data}"
adata_dir="${BENCHMARK_ADATA_DIR:-${SCRIPT_DIR}/benchmark/data}"
results_dir="${BENCHMARK_RESULTS_DIR:-${SCRIPT_DIR}/benchmark/outputs_methods/}"
metrics_dir="${BENCHMARK_METRICS_DIR:-${SCRIPT_DIR}/benchmark/outputs_metrics}"

# =============================================================================
# Configuration variables
# =============================================================================
train_data="${1:-future}"
output_mode="${2:-full_test}"
simulation="${3:-simul_replicates}"
global_run_name="${4:-test_run_3}"

DEFAULT_METHODS_JSON=$(run_in_runner python "$SCRIPT_DIR/utils/methods_registry.py" default-methods-json 2>/dev/null || cat <<EOF
{"FLeCS":1,"FLeCS-TPs":1,"scNODE":1,"reference_fitting":1,"CardamomOT":1,"GENIE3":1,"dynGENIE3":1,"PEARSON":1,"RENGE":1,"TrajectoryNet":1,"WaddingtonOT":1}
EOF
)
run_methods_json="${5:-$DEFAULT_METHODS_JSON}"

future_start_tp="${6:-10}"
replicates_number="${7:-2}"
simul_ko_genes="${8:-all}"
simulator_backend="${9:-harissa}"
perturbation_training="${10:-false}"
ko_output_genes="${11:-none}"
restart_mode="${12:-save}"

# =============================================================================
# Logging and error handling
# =============================================================================
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2; }
log_warning() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $1" >&2; }

set -e
trap 'log_error "Script failed at line $LINENO"' ERR

# Never delete data or results on failure: inputs may be user-provided, and
# completed method outputs (.done markers) allow restart_mode=save to resume.
# Set to 1 when benchmark_run.py reports that some method runs failed
# (exit code 3); the run continues with metrics and a summary is printed.
methods_failed=0

# Print every failed method run of this benchmark (all folds), as listed by
# single_run.sh in <run folder>/failed_methods.tsv.
print_failure_summary() {
    local run_root="$results_dir/${global_run_name}_${train_data}_${output_mode}"
    log_error "=========================================="
    log_error "Failed method runs (method / data file / exit code):"
    local f fold
    while IFS= read -r f; do
        fold="${f#"$run_root"}"; fold="${fold%/failed_methods.tsv}"; fold="${fold#/}"
        tail -n +2 "$f" | while IFS=$'\t' read -r m d rc; do
            log_error "  ${fold:+[$fold] }$m / $d / exit $rc"
        done
    done < <(find "$run_root" -name failed_methods.tsv 2>/dev/null | sort)
    log_error "Fix the errors and re-run: restart_mode=save only re-runs the failed methods."
    log_error "=========================================="
}

cleanup() {
    local exit_code=$?
    [ "$methods_failed" = "1" ] && print_failure_summary
    if [ $exit_code -ne 0 ] && [ $exit_code -ne 3 ]; then
        log_error "Benchmark run failed (exit code $exit_code). No files were removed; fix the error and re-run (restart_mode=save resumes completed methods)."
    fi
    exit $exit_code
}

trap cleanup EXIT

# =============================================================================
# Input validation
# =============================================================================
validate_inputs() {
    [[ "$train_data" =~ ^(full|future|leave-one-out|subsample_full)$ ]] || {
        log_error "Invalid train_data: $train_data. Must be one of: full, future, leave-one-out, subsample_full"
        exit 1
    }

    [[ "$output_mode" =~ ^(full_full|full_train|full_test|no_traj)$ ]] || {
        log_error "Invalid output_mode: $output_mode. Must be one of: full_full, full_train, full_test, no_traj"
        exit 1
    }
    
    [[ "$simulation" =~ ^(simul_replicates|simul_ko|false)$ ]] || {
        log_error "Invalid simulation: $simulation"
        exit 1
    }

    [[ "$simulator_backend" =~ ^(harissa|boolode|sergio|dyngen)$ ]] || {
        log_error "Invalid simulator_backend: $simulator_backend. Must be one of: harissa, boolode, sergio, dyngen"
        exit 1
    }

    case "$perturbation_training" in
        true|false) ;;
        *)
        log_error "Invalid perturbation_training: $perturbation_training. Must be 'true' or 'false'"
        exit 1
        ;;
    esac

    case "$restart_mode" in
        save|rerun) ;;
        *)
        log_error "Invalid restart_mode: $restart_mode. Must be 'save' or 'rerun'"
        exit 1
        ;;
    esac
    
    [[ "$train_data" == "future" ]] && {
        [[ "$future_start_tp" =~ ^[0-9]+$ ]] || {
            log_error "future_start_tp must be a non-negative integer"
            exit 1
        }
    }

    [[ "$replicates_number" =~ ^[0-9]+$ ]] || {
        log_error "replicates_number must be a positive integer"
        exit 1
    }

    if ! echo "$run_methods_json" | jq -e . >/dev/null 2>&1; then
        log_error "Invalid JSON in run_methods_json. Evaluated value was:"
        echo "$run_methods_json" >&2
        exit 1
    fi
}

# =============================================================================
# Simulation helpers
# =============================================================================

_simulator_dir() {
    case "$simulator_backend" in
        boolode) echo "boolode" ;;
        sergio)   echo "SERGIO" ;;
        dyngen)   echo "dyngen" ;;
        *)        echo "Harissa" ;;
    esac
}

_simulator_label() {
    case "$simulator_backend" in
        boolode) echo "BoolODE" ;;
        sergio)   echo "SERGIO" ;;
        dyngen)   echo "dyngen" ;;
        *)        echo "Harissa" ;;
    esac
}

# Run a simulator script with common arguments
# Usage: _run_simulator <script_name> [extra_args...]
_run_simulator() {
    local script_name="$1"; shift
    local sim_dir
    sim_dir=$(_simulator_dir)
    local label
    label=$(_simulator_label)

    # Optional shared simulator settings (from the config's simulation section)
    local sim_opts=()
    [ -n "${BENCHMARK_SIM_N_CELLS:-}" ] && sim_opts+=(--n_cells "$BENCHMARK_SIM_N_CELLS")
    [ -n "${BENCHMARK_SIM_N_TIME_BINS:-}" ] && sim_opts+=(--n_time_bins "$BENCHMARK_SIM_N_TIME_BINS")
    [ -n "${BENCHMARK_SIM_SEED:-}" ] && sim_opts+=(--seed "$BENCHMARK_SIM_SEED")

    log "Running $label/$script_name ${sim_opts[*]:-}"
    run_in_runner python "simulator/$sim_dir/$script_name" \
        -o "$input_dir" -a "$adata_dir" -n "$replicates_number" "${sim_opts[@]}" "$@" || {
        log_error "$label simulation failed"
        exit 1
    }
}

# =============================================================================
# Simulation functions
# =============================================================================

run_simulation_replicates() {
    _run_simulator "${simulator_backend}_simulate_custom.py"
    log "Replicates simulation completed successfully"
}

# KO genes to simulate: training KOs (simul_ko_genes) plus the KOs to predict
# (ko_output_genes), so held-out KOs have a ground truth. Methods only see the
# training KOs (benchmark_run.py --train_ko_genes builds their inputs).
ko_genes_to_simulate() {
    local train="$simul_ko_genes" test="$ko_output_genes"
    if [ "$train" = "all" ] || [ "$test" = "all" ]; then echo "all"; return; fi
    if [ "$test" = "none" ] || [ -z "$test" ]; then echo "$train"; return; fi
    if [ "$train" = "none" ] || [ -z "$train" ]; then echo "$test"; return; fi
    printf '%s,%s\n' "$train" "$test" | tr ',' '\n' | sed '/^$/d' | awk '!seen[$0]++' | paste -sd, -
}

run_simulation_ko() {
    local genes
    genes=$(ko_genes_to_simulate)
    log "Simulating KOs: $genes (training: $simul_ko_genes, to predict: $ko_output_genes)"
    _run_simulator "${simulator_backend}_simulate_custom_ko.py" -k "$genes"
    log "KO simulation completed successfully"
}

convert_simulated_datasets() {
    log "Converting simulated datasets to h5ad format"
    local gene_file_name="$input_dir/panel_genes.txt"
    # Degradation-rate priors are optional (only simulators that define them
    # in timepoint units write this file: Harissa, BoolODE)
    local degradation_file_name="$input_dir/Rates/degradation_rates.txt"

    [ -e "$gene_file_name" ] || {
        log_error "Required simulation file not found: $gene_file_name"
        exit 1
    }
    if [ -e "$degradation_file_name" ]; then
        log "Using simulator degradation rates as priors: $degradation_file_name"
    else
        log "No simulator degradation rates (methods use their own priors)"
    fi
    
    local failed_files=()
    local simulation_name
    simulation_name=$(_simulator_label)

    for i in "$input_dir"/data*; do
        [ -e "$i" ] || continue
        local filename_no_path=$(basename "${i%.txt}")
        local dataset_id="WT"
        if [[ "$simulation" == "simul_ko" && "$filename_no_path" == *_ko_* ]]; then
            dataset_id="${filename_no_path##*_ko_}"
        fi
        log "Converting file: $filename_no_path"
        
        run_in_runner python simulator/utils/convert_sim_data_to_ad.py \
            -i "$i" -g "$gene_file_name" -d "$degradation_file_name" -o "$adata_dir" \
            -n "$filename_no_path" -s "$dataset_id" -p "$simulation_name" || failed_files+=("$filename_no_path")
    done
    
    [ ${#failed_files[@]} -eq 0 ] || {
        log_error "Failed to convert ${#failed_files[@]} file(s): ${failed_files[*]}"
        exit 1
    }
    log "All simulated datasets converted successfully"
}

merge_ko_datasets() {
    log "Merging all ko datasets into a single anndata object"
    run_in_runner python simulator/utils/merge_ko_datasets.py \
        -o "$adata_dir" || {
        log_error "Merging failed"
        exit 1
    }
    log "KO datasets merged successfully"
}

# =============================================================================
# Main execution
# =============================================================================
main() {
    log "=========================================="
    log "Starting benchmark run"
    log "=========================================="
    log "Configuration: train_data=$train_data, output_mode=$output_mode, simulation=$simulation, simulator_backend=$simulator_backend, run_name=$global_run_name, ko_output_genes=$ko_output_genes, perturbation_training=$perturbation_training"
    log "=========================================="
    
    validate_inputs
    
    # Handle simulation if needed
    if [ "$simulation" != "false" ]; then
        log "Simulation set to: $simulation. Generating new data."
        [ "$simulation" = "simul_replicates" ] && run_simulation_replicates || run_simulation_ko
        convert_simulated_datasets
        [ "$simulation" = "simul_ko" ] && merge_ko_datasets
        
        # Save ground truth network files
        log "Saving ground truth network files"
        grn_true_dir="$results_dir/${global_run_name}_${train_data}_${output_mode}/GRN_true"
        mkdir -p "$grn_true_dir"
        
        # Simulators write the true network next to their output folder
        # (<input_dir>/../True), so follow BENCHMARK_INPUT_DIR / paths.input_dir.
        true_dir="$(cd "$input_dir/.." && pwd)/True"
        for ((i=1; i<=replicates_number; i++)); do
            source_file="$true_dir/inter_${i}.npy"
            target_file="$grn_true_dir/data_${i}_GRN.npy"
            
            if [ -e "$source_file" ]; then
                cp "$source_file" "$target_file"
                log "Saved $source_file to $target_file"
            else
                log_warning "Source file not found: $source_file"
            fi
        done
    else
        log "Using existing data"
    fi

    log "Running benchmark analysis: train_data=$train_data, output_mode=$output_mode"
    
    if [ "$simulation" = "simul_replicates" ] && [ "$ko_output_genes" != "none" ]; then
        log_warning "ko_output_genes=$ko_output_genes but simulation=simul_replicates simulates no knockout: KO predictions will not be scored."
    fi
    local args=(-a "$adata_dir" -n "$global_run_name" -r "$results_dir" -m "$run_methods_json" --train_data "$train_data" --output_mode "$output_mode" --ko_output_genes "$ko_output_genes" --train_ko_genes "$simul_ko_genes" --perturbation_training "$perturbation_training" --restart_mode "$restart_mode")
    [[ "$train_data" == "future" ]] && args+=(-t "$future_start_tp")
    
    local bench_rc=0
    BENCHMARK_EXECUTION_JSON="${BENCHMARK_EXECUTION_JSON:-}" run_in_runner python benchmark_run.py "${args[@]}" || bench_rc=$?
    if [ "$bench_rc" -eq 3 ]; then
        methods_failed=1
        log_warning "Some method runs failed; computing metrics for the methods that succeeded."
    elif [ "$bench_rc" -ne 0 ]; then
        log_error "benchmark_run.py failed"
        exit 1
    fi

    log "=========================================="
    log "Method runs finished"
    log "=========================================="
    log "Calculating metrics and generating plots"
    log "=========================================="
    run_in_runner python post_analysis/compute_metrics.py "$global_run_name" "${train_data}_${output_mode}" "$simulation" \
        --ko_output_genes "$ko_output_genes" --adata_dir "$adata_dir" --methods_dir "$results_dir" --metrics_dir "$metrics_dir" || {
        log_error "compute_metrics.py failed"
        exit 1
    }
    run_in_runner python ranking_table.py "${metrics_dir}/${global_run_name}_${train_data}_${output_mode}/" || {
        log_error "ranking_table.py failed"
        exit 1
    }
    log "Metrics calculation and plotting completed successfully"

    if [ "$methods_failed" = "1" ]; then
        exit 3   # summary printed by the EXIT trap
    fi
}

main "$@"
