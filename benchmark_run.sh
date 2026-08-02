#!/bin/bash
# =============================================================================
# benchmark_run.sh - Main benchmark runner script
# =============================================================================
# Usage: ./benchmark_run.sh <train_data> <output_mode> <simulation> <global_run_name>
#        <run_methods_json> <future_start_tp> <replicates_number> [ko_output_genes] [simulator_backend] [perturbation_training]
#        ./benchmark_run.sh --config <config.yaml>
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
# perturbation_training: Whether to run perturbation-training variants of supported methods ("true" or "false"; default: true)
# restart_mode: "save" to skip methods with .done markers, "rerun" to re-run all (default: save)
# =============================================================================

# Config shortcut
if [[ "${1:-}" == "--config" || "${1:-}" == "-c" ]]; then
    if [ -z "${2:-}" ]; then
        echo "Usage: $0 --config <config.yaml>" >&2
        exit 1
    fi
    exec python benchmark_run_config.py --config "$2"
fi

# Set default paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
input_dir="${BENCHMARK_INPUT_DIR:-${SCRIPT_DIR}/simulator/custom_network/Data}"
adata_dir="${BENCHMARK_ADATA_DIR:-${SCRIPT_DIR}/benchmark/data}"
results_dir="${BENCHMARK_RESULTS_DIR:-${SCRIPT_DIR}/benchmark/outputs_methods/}"
runner_env="${BENCHMARK_RUNNER_ENV:-cardamom_env}"

# =============================================================================
# Configuration variables
# =============================================================================
train_data="${1:-future}"
output_mode="${2:-full_test}"
simulation="${3:-simul_replicates}"
global_run_name="${4:-test_run_3}"

DEFAULT_METHODS_JSON=$(python "$SCRIPT_DIR/utils/methods_registry.py" default-methods-json 2>/dev/null || cat <<EOF
{"FLeCS":"1","FLeCS-TPs":"1","scNODE":"1","reference_fitting":"1","CardamomOT":"1","GENIE3":"1","PEARSON":"1"}
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

cleanup() {
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        log_error "Cleaning up after failure..."
        if [ -n "${input_dir:-}" ] && [ -d "$input_dir" ]; then
            rm -rf "$input_dir"/*
        fi
        if [ -n "${results_dir:-}" ] && [ -d "$results_dir" ]; then
            rm -rf "$results_dir/${global_run_name}_${train_data}_${output_mode}" 2>/dev/null
        fi
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
        sergio)   echo "sergio" ;;
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

    log "Running $label/$script_name"
    conda run -n "$runner_env" --live-stream python "simulator/$sim_dir/$script_name" \
        -o "$input_dir" -a "$adata_dir" -n "$replicates_number" "$@" || {
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

run_simulation_ko() {
    _run_simulator "${simulator_backend}_simulate_custom_ko.py" -k "$simul_ko_genes"
    log "KO simulation completed successfully"
}

convert_simulated_datasets() {
    log "Converting simulated datasets to h5ad format"
    local gene_file_name=("$input_dir"/panel_gene*)
    local degradation_file_name=("$input_dir"/Rates/degradation_rates.txt*)
    
    [ -e "$gene_file_name" ] && [ -e "$degradation_file_name" ] || {
        log_error "Required simulation files not found in $input_dir"
        exit 1
    }
    
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
        
        conda run -n "$runner_env" --live-stream python simulator/utils/convert_sim_data_to_ad.py \
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
    conda run -n "$runner_env" --live-stream python simulator/utils/merge_ko_datasets.py \
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
        
        for ((i=1; i<=replicates_number; i++)); do
            source_file="$SCRIPT_DIR/simulator/custom_network/True/inter_${i}.npy"
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
    
    local args=(-a "$adata_dir" -n "$global_run_name" -r "$results_dir" -m "$run_methods_json" --train_data "$train_data" --output_mode "$output_mode" --ko_output_genes "$ko_output_genes" --perturbation_training "$perturbation_training" --restart_mode "$restart_mode")
    [[ "$train_data" == "future" ]] && args+=(-t "$future_start_tp")
    
    BENCHMARK_EXECUTION_JSON="${BENCHMARK_EXECUTION_JSON:-}" conda run -n "$runner_env" --live-stream python benchmark_run.py "${args[@]}" || {
        log_error "benchmark_run.py failed"
        exit 1
    }

    log "=========================================="
    log "Benchmark run completed successfully"
    log "=========================================="
    log "Calculating metrics and generating plots"
    log "=========================================="
    conda run -n "$runner_env" --live-stream python post_analysis/compute_metrics.py "$global_run_name" "${train_data}_${output_mode}" "$simulation" "--ko_output_genes" "$ko_output_genes" || {
        log_error "compute_metrics.py failed"
        exit 1
    }
    conda run -n "$runner_env" --live-stream python ranking_table.py "benchmark/outputs_metrics/${global_run_name}_${train_data}_${output_mode}/" || {
        log_error "ranking_table.py failed"
        exit 1
    }
    log "Metrics calculation and plotting completed successfully"
}

main "$@"
