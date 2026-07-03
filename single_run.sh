#!/bin/bash
# Single run script for benchmarking methods
# Usage: single_run.sh <adata_dir> <results_dir> <run_methods_json> <output_mode> [ko_output_genes] [perturbation_training] [restart_mode]
#
# restart_mode: "save" to skip methods with .done markers, "rerun" to re-run all (default: save)
#
# output_mode: Controls what simulation the methods produce. One of:
#   full_full  - simulate all train and test timepoints from the first timepoint
#   full_train - simulate all training timepoints from the first timepoint
#   full_test  - simulate all test timepoints from the last timepoint before test
#   no_traj    - do not output simulation

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# Configuration and Argument Parsing
# =============================================================================

adata_dir="$1"
results_dir="$2"
run_methods_json="$3"
output_mode="${4:-full_test}"
ko_output_genes="${5:-none}"
perturbation_training="${6:-true}"
restart_mode="${7:-save}"

# =============================================================================
# Logging Functions
# =============================================================================

log()         { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }
log_error()   { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2; }
log_warning() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $1" >&2; }

# =============================================================================
# Execution Configuration (conda/docker/apptainer)
# =============================================================================

EXECUTION_JSON="${BENCHMARK_EXECUTION_JSON:-}"

# ── Registry query helper ───────────────────────────────────────────
# Calls utils/methods_registry.py for any CLI command and caches results
# in a temp file keyed by method+command to avoid repeated Python startup.
# Usage: _registry_query <command> <method_key>
_registry_query() {
    local cmd="$1" method="$2"
    python "$SCRIPT_DIR/utils/methods_registry.py" "$cmd" "$method"
}

method_runner_mode() {
    local method="$1"
    if [ -n "$EXECUTION_JSON" ]; then
        local mode
        mode=$(echo "$EXECUTION_JSON" | jq -r --arg method "$method" '.method_runners[$method] // .default_runner // empty' 2>/dev/null)
        if [ -n "$mode" ] && [ "$mode" != "null" ]; then
            echo "$mode"
            return
        fi
        local engine
        engine=$(echo "$EXECUTION_JSON" | jq -r '.container_engine // empty' 2>/dev/null)
        if [ -n "$engine" ] && [ "$engine" != "null" ]; then
            echo "$engine"
            return
        fi
    fi
    # Fall back to registry default
    _registry_query runner-mode "$method"
}

method_conda_env() {
    local method="$1"
    if [ -n "$EXECUTION_JSON" ]; then
        local env
        env=$(echo "$EXECUTION_JSON" | jq -r --arg method "$method" '.conda_envs[$method] // .conda_envs.default // empty' 2>/dev/null)
        if [ -n "$env" ] && [ "$env" != "null" ]; then
            echo "$env"
            return
        fi
    fi
    _registry_query conda-env "$method"
}

method_container_image() {
    local method="$1"
    if [ -n "$EXECUTION_JSON" ]; then
        local image
        image=$(echo "$EXECUTION_JSON" | jq -r --arg method "$method" '.container_images[$method] // .container_images.default // empty' 2>/dev/null)
        if [ -n "$image" ] && [ "$image" != "null" ]; then
            echo "$image"
            return
        fi
    fi
    _registry_query container-image "$method"
}

method_exec() {
    local method="$1"
    shift
    local wrap_conda="true"
    if [ "${1:-}" = "--no-conda-wrap" ]; then
        wrap_conda="false"
        shift
    fi

    local mode
    mode=$(method_runner_mode "$method")

    case "$mode" in
        conda)
            if [ "$wrap_conda" = "true" ]; then
                local env
                env=$(method_conda_env "$method")
                conda run -n "$env" --live-stream "$@"
            else
                if [ "$method" = "CardamomOT" ]; then
                    BENCHMARK_CARDAMOMOT_CONDA_ENV="$(method_conda_env "$method")" "$@"
                else
                    "$@"
                fi
            fi
            ;;
        docker)
            local image
            image=$(method_container_image "$method")
            if [ -z "$image" ]; then
                log_error "No container image configured for method '$method'."
                return 1
            fi
            if ! command -v docker &> /dev/null 2>&1; then
                log_error "docker not found. Install docker or switch runner mode to conda."
                return 1
            fi
            docker run --rm -v "$SCRIPT_DIR":/work -w /work -e BENCHMARK_RUNNER_MODE=container "$image" "$@"
            ;;
        apptainer)
            local image
            image=$(method_container_image "$method")
            if [ -z "$image" ]; then
                log_error "No container image configured for method '$method'."
                return 1
            fi
            if ! command -v apptainer &> /dev/null 2>&1; then
                log_error "apptainer not found. Install apptainer or switch runner mode to conda."
                return 1
            fi
            apptainer exec --bind "$SCRIPT_DIR":/work --env BENCHMARK_RUNNER_MODE=container "$image" "$@"
            ;;
        *)
            log_error "Unsupported runner mode '$mode' for method '$method'."
            return 1
            ;;
    esac
}

# =============================================================================
# Error Handling
# =============================================================================

# Set trap for cleanup on error
cleanup() {
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        log_error "Script failed at line $LINENO. Exit code: $exit_code"
        # Remove partial results if results_dir is set and exists
        if [ -n "${results_dir:-}" ] && [ -d "$results_dir" ]; then
            log_warning "Removing partial results from $results_dir"
            rm -rf "$results_dir"/*
        fi
        exit $exit_code
    fi
}

trap cleanup EXIT

# =============================================================================
# Input Validation (Optimized)
# =============================================================================

# Validate required arguments
if [ -z "$adata_dir" ] || [ -z "$results_dir" ] || [ -z "$run_methods_json" ]; then
    log_error "Missing required arguments"
    echo "Usage: single_run.sh <adata_dir> <results_dir> <run_methods_json> [output_mode] [ko_output_genes] [perturbation_training]"
    exit 1
fi

case "$perturbation_training" in
    true|false) ;;
    *)
        log_error "Invalid perturbation_training: $perturbation_training. Must be 'true' or 'false'."
        exit 1
        ;;
esac

case "$restart_mode" in
    save|rerun) ;;
    *)
        log_error "Invalid restart_mode: $restart_mode. Must be 'save' or 'rerun'."
        exit 1
        ;;
esac

# Check if jq is available (cached check)
if ! command -v jq &> /dev/null 2>&1; then
    log_error "jq not found. Please install jq for JSON parsing."
    exit 1
fi

# Check if run_methods_json is valid JSON (cached check)
if ! echo "$run_methods_json" | jq . &> /dev/null 2>&1; then
    log_error "Invalid JSON in run_methods_json"
    exit 1
fi

if [ -n "$EXECUTION_JSON" ]; then
    if ! echo "$EXECUTION_JSON" | jq . &> /dev/null 2>&1; then
        log_error "Invalid JSON in BENCHMARK_EXECUTION_JSON"
        exit 1
    fi
fi

# Cache enabled methods once and use a compact membership test.
enabled_methods=$(echo "$run_methods_json" | jq -r 'to_entries[] | select(.value == 1) | .key' | tr '\n' '|')
enabled_methods="${enabled_methods%|}"

method_enabled() {
    case "|$enabled_methods|" in
        *"|$1|"*) echo 1 ;;
        *) echo 0 ;;
    esac
}

# Validate conda environments only for enabled methods that need conda.
needs_conda=0
REQUIRED_CONDA_ENVS=""
for method in $(_registry_query keys); do
    [ "$(method_enabled "$method")" = "1" ] || continue
    [ "$(method_runner_mode "$method")" = "conda" ] || continue

    needs_conda=1
    env_name=$(method_conda_env "$method")
    case " $REQUIRED_CONDA_ENVS " in
        *" $env_name "*) ;;
        *) REQUIRED_CONDA_ENVS="$REQUIRED_CONDA_ENVS $env_name" ;;
    esac
done

if [ "$needs_conda" = "1" ]; then
    if ! command -v conda &> /dev/null 2>&1; then
        log_error "conda not found. Please install conda first."
        exit 1
    fi
    log "Validating conda environments..."
    for env in $REQUIRED_CONDA_ENVS; do
        if ! conda env list 2>/dev/null | grep -q "$env"; then
            log_error "Conda environment '$env' not found. Please create it first."
            exit 1
        fi
    done
    log "All conda environments validated successfully"
else
    log "Skipping conda environment validation (no conda-based methods selected)"
fi

# Validate directories exist and are writable (optimized single pass)
for dir in "$adata_dir" "$results_dir"; do
    if [ ! -d "$dir" ]; then
        log_error "Directory does not exist: $dir"
        exit 1
    fi
    if [ ! -w "$dir" ]; then
        log_error "Directory is not writable: $dir"
        exit 1
    fi
done

# Create results directory if it doesn't exist
results_dir_run="$results_dir"
if [ ! -d "$results_dir_run" ]; then
    mkdir -p "$results_dir_run"
    log "Created results directory: $results_dir_run"
fi

# =============================================================================
# Method registry and execution helpers
# =============================================================================
# All method metadata now lives in methods_registry.yaml.
# The _registry_query helper (defined above) calls utils/methods_registry.py.

method_spec() {
    # Returns: type|script|pass_out|pass_ko|supports_perturbation|extra_args...
    _registry_query entrypoint "$1"
}

method_supports_perturbation() {
    _registry_query supports-perturbation "$1"
}

# Usage: run_method_variant <method_name> <perturbation_flag> <output_subdir>
run_method_variant() {
    local method_name="$1"
    local pert_flag="${2:-false}"
    local out_subdir="${3:-$method_name}"

    local label="$method_name"
    [ "$pert_flag" = "true" ] && label="$method_name (perturbation training)"

    local cmd_type script pass_out pass_ko supports_perturbation method_extra
    IFS='|' read -r cmd_type script pass_out pass_ko supports_perturbation method_extra <<EOF
$(method_spec "$method_name")
EOF

    # ── Save / restart logic ──────────────────────────────────────────
    local data_basename
    data_basename=$(basename "$i" .h5ad)
    local done_marker="$results_dir_run/$out_subdir/${data_basename}.done"

    if [ "$restart_mode" = "save" ] && [ -f "$done_marker" ]; then
        log "Skipping $label on $i (.done marker exists)"
        return 0
    fi

    if [ "$restart_mode" = "rerun" ] && [ -f "$done_marker" ]; then
        log "Re-running $label on $i (restart_mode=rerun, removing old marker)"
        rm -f "$done_marker"
    fi

    log "Running $label on $i"

    case "$cmd_type" in
        python)
            local extra_args=()
            [ "$pass_out" = "true" ] && extra_args+=(-u "$output_mode")
            [ "$pass_ko" = "true" ] && extra_args+=(-k "$ko_output_genes")
            [ "$supports_perturbation" = "true" ] && extra_args+=(-p "$pert_flag")
            [ -n "$method_extra" ] && extra_args+=($method_extra)

            if ! method_exec "$method_name" python -u "$script" -i "$i" -o "$results_dir_run/$out_subdir/" "${extra_args[@]}" 2>&1; then
                log_error "$label failed on $i. Check the error output above for details."
                return 1
            fi
            ;;
        bash)
            local extra_args=()
            [ "$pass_out" = "true" ] && extra_args+=("$output_mode")
            [ "$supports_perturbation" = "true" ] && extra_args+=("$pert_flag")
            [ "$pass_ko" = "true" ] && extra_args+=("$ko_output_genes")

            if ! method_exec "$method_name" --no-conda-wrap bash "$script" "$i" "$results_dir_run/$out_subdir/" "${extra_args[@]}" 2>&1; then
                log_error "$label failed on $i. Check the error output above for details."
                return 1
            fi
            ;;
        *)
            log_error "Unknown command type '$cmd_type' for method '$method_name'."
            return 1
            ;;
    esac

    # Write .done marker on success
    mkdir -p "$(dirname "$done_marker")"
    touch "$done_marker"
    log "$label completed successfully on $i"
    return 0
}

# Run one or both variants of a method (standard + optional perturbation-training)
run_method() {
    local method_name="$1"
    [ "$(method_enabled "$method_name")" != "1" ] && { log "Skipping $method_name"; return 0; }

    run_method_variant "$method_name" "false" || return 1

    if [ "$(method_supports_perturbation "$method_name")" = "1" ] && [ "$perturbation_training" = "true" ]; then
        run_method_variant "$method_name" "true" "${method_name}_perturbation_training" || return 1
    fi
    return 0
}

# =============================================================================
# Process Each Data File
# =============================================================================

log "Starting processing of data files in: $adata_dir"

if [ ! -d "$adata_dir" ]; then
    log_error "Adata directory does not exist: $adata_dir"
    exit 1
fi

shopt -s nullglob
data_files=("$adata_dir"/data*)
shopt -u nullglob

if [ ${#data_files[@]} -eq 0 ]; then
    log "No data files found in $adata_dir"
    log "All data files processed"
    log "Benchmark run completed successfully"
    exit 0
fi

for i in "${data_files[@]}"; do
    [ -e "$i" ] || continue
    log "Processing data file: $i"

    for method in $(_registry_query keys); do
        run_method "$method"
    done

    log "Completed processing: $i"
done

log "All data files processed"
log "Benchmark run completed successfully"
