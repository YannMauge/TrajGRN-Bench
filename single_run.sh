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
perturbation_training="${6:-false}"
restart_mode="${7:-save}"

# =============================================================================
# Logging Functions
# =============================================================================

# SHA-256 of stdin, portable across Linux (sha256sum) and macOS (shasum).
sha256_stdin() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum | awk '{print $1}'
    else
        shasum -a 256 | awk '{print $1}'
    fi
}

log()         { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }
log_error()   { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2; }
log_warning() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $1" >&2; }

# =============================================================================
# Execution Configuration (conda/docker/apptainer)
# =============================================================================

EXECUTION_JSON="${BENCHMARK_EXECUTION_JSON:-}"

# ── Registry query helper ───────────────────────────────────────────
# Calls utils/methods_registry.py for any CLI command and caches the answer
# in a temp file keyed by command+method, so each query starts Python once
# per run (the registry does not change during a run).
# Usage: _registry_query <command> [method_key]
_REGISTRY_CACHE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/trajgrn-registry.XXXXXX")"
_registry_query() {
    local cmd="$1" method="${2:-}"
    local cache="$_REGISTRY_CACHE_DIR/${cmd}__${method}"
    if [ ! -f "$cache" ]; then
        python "$SCRIPT_DIR/utils/methods_registry.py" "$cmd" ${method:+"$method"} > "$cache.tmp" \
            && mv "$cache.tmp" "$cache" || { rm -f "$cache.tmp"; return 1; }
    fi
    cat "$cache"
}

method_runner_mode() {
    local method="$1"
    if [ -n "$EXECUTION_JSON" ]; then
        local mode
        # method_runners already encodes the precedence resolved by
        # benchmark_run_config.py; container_engine picks the engine used for
        # container runs (docker or apptainer).
        mode=$(echo "$EXECUTION_JSON" | jq -r --arg method "$method" '
            (.method_runners[$method] // .default_runner // empty) as $m
            | if $m == "docker" and (.container_engine // "") != "" then .container_engine else $m end' 2>/dev/null)
        if [ -n "$mode" ] && [ "$mode" != "null" ]; then
            echo "$mode"
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

method_pixi_env() {
    local method="$1" env=""
    if [ -n "$EXECUTION_JSON" ]; then
        env=$(echo "$EXECUTION_JSON" | jq -r --arg method "$method" '.pixi_envs[$method] // .pixi_envs.default // empty' 2>/dev/null)
    fi
    [ -n "$env" ] && [ "$env" != "null" ] || env=$(_registry_query pixi-env "$method")
    echo "$env"
}

# GPU runs (execution.gpu: true): use the "<env>-cuda" pixi environment /
# "benchmark/<env>-cuda" image when the method's environment has one.
gpu_enabled() {
    [ -n "$EXECUTION_JSON" ] && [ "$(echo "$EXECUTION_JSON" | jq -r '.gpu // false' 2>/dev/null)" = "true" ]
}

has_cuda_variant() {  # <pixi_env>
    [ -n "$EXECUTION_JSON" ] && echo "$EXECUTION_JSON" | jq -e --arg env "$1" '(.cuda_envs // []) | index($env) != null' >/dev/null 2>&1
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

# Directories that must be visible inside containers, as unique absolute
# paths. They are mounted at the *same* path as on the host, so any path
# passed to a method (absolute or relative to the cwd) resolves identically.
container_mount_dirs() {
    local d abs seen="|"
    for d in "$SCRIPT_DIR" "$PWD" "$adata_dir" "$results_dir"; do
        [ -n "$d" ] && [ -d "$d" ] || continue
        abs="$(cd "$d" && pwd)"
        case "$seen" in *"|$abs|"*) continue ;; esac
        seen="$seen$abs|"
        printf '%s\n' "$abs"
    done
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
        pixi)
            # Activate the method's pixi environment in a subshell and exec the
            # command directly (arguments are passed through unchanged).
            local env
            env=$(method_pixi_env "$method")
            if gpu_enabled && has_cuda_variant "$env"; then
                env="${env}-cuda"
            fi
            if ! command -v pixi >/dev/null 2>&1; then
                log_error "pixi not found. Install pixi (https://pixi.sh) or switch runner mode."
                return 1
            fi
            (
                eval "$(pixi shell-hook --manifest-path "$SCRIPT_DIR/pixi.toml" -e "$env" -s bash)" || exit 1
                export BENCHMARK_RUNNER_MODE=pixi
                exec "$@"
            )
            ;;
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
            local docker_args=(--rm -w "$PWD" -e BENCHMARK_RUNNER_MODE=container)
            if gpu_enabled; then
                local cuda_image="${image%%:*}-cuda:${image##*:}"
                if has_cuda_variant "$(method_pixi_env "$method")"; then
                    image="$cuda_image"
                fi
                docker_args+=(--gpus all)
            fi
            local mount_dir
            while IFS= read -r mount_dir; do
                docker_args+=(-v "$mount_dir:$mount_dir")
            done < <(container_mount_dirs)
            # Run as the host user so outputs are owned by them (matters on
            # Linux); HOME and tool caches go to /tmp since that uid has no
            # home directory inside the image.
            docker_args+=(--user "$(id -u):$(id -g)" -e HOME=/tmp
                          -e MPLCONFIGDIR=/tmp/matplotlib -e NUMBA_CACHE_DIR=/tmp/numba
                          -e XDG_CACHE_HOME=/tmp/.cache)
            docker run "${docker_args[@]}" "$image" "$@"
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
            local apptainer_args=(--pwd "$PWD" --env BENCHMARK_RUNNER_MODE=container)
            if gpu_enabled; then
                if has_cuda_variant "$(method_pixi_env "$method")"; then
                    image="${image%%:*}-cuda:${image##*:}"
                fi
                apptainer_args+=(--nv)
            fi
            local mount_dir
            while IFS= read -r mount_dir; do
                apptainer_args+=(--bind "$mount_dir:$mount_dir")
            done < <(container_mount_dirs)
            apptainer exec "${apptainer_args[@]}" "$image" "$@"
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

# Report failures without deleting anything: completed methods keep their
# outputs and .done markers, and unfinished ones (no .done marker) are simply
# re-run and overwritten next time (restart_mode=save).
cleanup() {
    local exit_code=$?
    [ -n "${_REGISTRY_CACHE_DIR:-}" ] && rm -rf "${_REGISTRY_CACHE_DIR:?}"
    # Exit code 3 = some methods failed (already summarised above).
    if [ $exit_code -ne 0 ] && [ $exit_code -ne 3 ]; then
        log_error "single_run.sh failed (exit code $exit_code). Completed method outputs in ${results_dir:-<unset>} were kept; methods without a .done marker will be re-run."
    fi
    exit $exit_code
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
# Accepts 1, "1", true and "true" (the encodings allowed by the config schema).
enabled_methods=$(echo "$run_methods_json" | jq -r 'to_entries[] | select(.value == 1 or .value == true or .value == "1" or .value == "true") | .key' | tr '\n' '|')
enabled_methods="${enabled_methods%|}"
if [ -z "$enabled_methods" ]; then
    log_warning "No method is enabled in run_methods_json: $run_methods_json"
fi

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
        if ! conda env list 2>/dev/null | awk '!/^#/ {print $1}' | grep -qxF "$env"; then
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
    local out_dir="$results_dir_run/$out_subdir"
    local done_marker="$out_dir/${data_basename}.done"

    # Everything that determines this result. Stored in the .done marker;
    # in save mode a result is reused only if the fingerprint still matches.
    local fp_text fp
    fp_text=$(printf '%s\n' \
        "fingerprint_version=1" \
        "data=$data_basename $current_data_hash" \
        "splits=$split_fingerprint" \
        "method=$method_name perturbation_training=$pert_flag" \
        "entrypoint=$cmd_type|$script|$pass_out|$pass_ko|$supports_perturbation|$method_extra" \
        "output_mode=$output_mode" \
        "ko_output_genes=$ko_output_genes")
    fp=$(printf '%s\n' "$fp_text" | sha256_stdin)

    if [ -f "$done_marker" ]; then
        local stored_fp
        stored_fp=$(sed -n 's/^fingerprint=//p' "$done_marker" | head -n 1)
        if [ "$restart_mode" = "rerun" ]; then
            log "Re-running $label on $i (restart_mode=rerun)"
        elif [ -z "$stored_fp" ]; then
            log_warning "Skipping $label on $i: legacy .done marker without fingerprint (inputs cannot be checked; use restart_mode=rerun to recompute)"
            return 0
        elif [ "$stored_fp" = "$fp" ]; then
            log "Skipping $label on $i (.done marker matches current inputs and settings)"
            return 0
        else
            log "Inputs or settings changed for $label on $i since the last run; re-running. Changes (< old, > new):"
            diff <(sed '1d' "$done_marker") <(printf '%s\n' "$fp_text") | grep '^[<>]' || true
        fi
        rm -f "$done_marker"
    fi

    # Clear this method's previous outputs for this data file only, so that
    # stale files (or CardamomOT's cached project folder, which would reuse
    # the old data) cannot survive into the new result.
    if [ -d "$out_dir" ]; then
        local stale
        for stale in "$out_dir/${data_basename}_"* "$out_dir/project_${data_basename}"; do
            [ -e "$stale" ] || continue
            rm -rf -- "$stale"
        done
    fi

    mkdir -p "$out_dir"
    log "Running $label on $i"

    local method_rc=0
    case "$cmd_type" in
        python)
            local extra_args=()
            [ "$pass_out" = "true" ] && extra_args+=(-u "$output_mode")
            [ "$pass_ko" = "true" ] && extra_args+=(-k "$ko_output_genes")
            [ "$supports_perturbation" = "true" ] && extra_args+=(-p "$pert_flag")
            [ -n "$method_extra" ] && extra_args+=($method_extra)

            method_exec "$method_name" python -u "$script" -i "$i" -o "$results_dir_run/$out_subdir/" "${extra_args[@]}" 2>&1 || method_rc=$?
            ;;
        bash)
            local extra_args=()
            [ "$pass_out" = "true" ] && extra_args+=("$output_mode")
            [ "$supports_perturbation" = "true" ] && extra_args+=("$pert_flag")
            [ "$pass_ko" = "true" ] && extra_args+=("$ko_output_genes")
            [ -n "$method_extra" ] && extra_args+=($method_extra)

            method_exec "$method_name" --no-conda-wrap bash "$script" "$i" "$results_dir_run/$out_subdir/" "${extra_args[@]}" 2>&1 || method_rc=$?
            ;;
        *)
            log_error "Unknown command type '$cmd_type' for method '$method_name'."
            method_rc=1
            ;;
    esac

    if [ "$method_rc" -ne 0 ]; then
        # 130 = interrupted (Ctrl-C): stop the whole run instead of continuing.
        if [ "$method_rc" -eq 130 ]; then
            log_error "$label interrupted on $i; stopping."
            exit 130
        fi
        log_error "$label failed on $i (exit code $method_rc). Continuing with the remaining methods; see the error output above."
        record_failure "$out_subdir" "$i" "$method_rc"
        return 1
    fi

    # Write .done marker (with input fingerprint) on success
    mkdir -p "$out_dir"
    printf 'fingerprint=%s\n%s\n' "$fp" "$fp_text" > "$done_marker"
    log "$label completed successfully on $i"
    return 0
}

# Run one or both variants of a method (standard + optional perturbation-training)
run_method() {
    local method_name="$1"
    [ "$(method_enabled "$method_name")" != "1" ] && { log "Skipping $method_name"; return 0; }

    # Variants are independent: a failure is recorded and the next one runs.
    run_method_variant "$method_name" "false" || true

    if [ "$(method_supports_perturbation "$method_name")" = "1" ] && [ "$perturbation_training" = "true" ]; then
        run_method_variant "$method_name" "true" "${method_name}_perturbation_training" || true
    fi
    return 0
}

# ── Failure bookkeeping (continue-on-error) ─────────────────────────────
# Failed method runs are listed in $failures_file (TSV) so benchmark_run.sh
# can aggregate them across folds. Exit code 3 = some method runs failed.
failures_file="$results_dir_run/failed_methods.tsv"
rm -f "$failures_file"
n_failures=0

record_failure() {
    [ -f "$failures_file" ] || printf 'method\tdata_file\texit_code\n' > "$failures_file"
    printf '%s\t%s\t%s\n' "$1" "$(basename "$2")" "$3" >> "$failures_file"
    n_failures=$((n_failures + 1))
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

# Fingerprint of the train/test split shared by all data files of this run
# (encodes train_data, future_start_tp, leave-one-out fold, subsampling).
split_fingerprint=""
for split_file in train_tps.npy test_tps.npy subsample_train_ids.npy; do
    if [ -f "$adata_dir/$split_file" ]; then
        split_fingerprint="$split_fingerprint$split_file:$(sha256_stdin < "$adata_dir/$split_file") "
    else
        split_fingerprint="$split_fingerprint$split_file:absent "
    fi
done

for i in "${data_files[@]}"; do
    [ -e "$i" ] || continue
    log "Processing data file: $i"
    current_data_hash=$(sha256_stdin < "$i")

    for method in $(_registry_query keys); do
        run_method "$method"
    done

    log "Completed processing: $i"
done

log "All data files processed"
if [ "$n_failures" -gt 0 ]; then
    log_error "$n_failures method run(s) failed in $results_dir_run (listed in $failures_file):"
    column -t -s $'\t' "$failures_file" >&2 2>/dev/null || cat "$failures_file" >&2
    exit 3
fi
log "Benchmark run completed successfully"
