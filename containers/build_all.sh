#!/usr/bin/env bash
# -------------------------------------------------------------------
# build_all.sh — Build all benchmark Docker images in one command.
#
# Usage:
#   bash containers/build_all.sh              # build all images
#   bash containers/build_all.sh flecs scnode # build only listed images
#   bash containers/build_all.sh --no-cache   # force full rebuild
#   bash containers/build_all.sh --push       # build then push to registry
#   bash containers/build_all.sh -j 4         # parallel builds (GNU parallel or xargs)
#
# Images are defined in methods_registry.yaml.
# Run `python utils/methods_registry.py container-images` to list them.
# -------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Declare all known images from the registry (methods_registry.yaml).
# Image list is built dynamically:  name -> "env_file|env_name|image_tag"
declare -A ALL_IMAGES
while IFS='|' read -r target env_file env_name image; do
    ALL_IMAGES["$target"]="${env_file}|${env_name}|${image}"
done < <(python "$REPO_ROOT/utils/methods_registry.py" container-images 2>/dev/null | python -c "
import json, sys
for d in json.load(sys.stdin):
    print('|'.join([d['target'], d['env_file'], d['env_name'], d['image']]))
")

DOCKER="${DOCKER:-docker}"
DOCKERFILE="containers/docker/Dockerfile"
NO_CACHE=""
PUSH=""
PARALLEL=0
TARGETS=()

# ── Parse arguments ─────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-cache) NO_CACHE="--no-cache" ;;
    --push)     PUSH="1" ;;
    -j)
      PARALLEL="${2:-1}"
      shift
      ;;
    *)
      if [[ -n "${ALL_IMAGES[$1]:-}" ]]; then
        TARGETS+=("$1")
      else
        echo "Warning: unknown image '$1', skipping" >&2
      fi
      ;;
  esac
  shift
done

# Default to all images if none specified
if [[ ${#TARGETS[@]} -eq 0 ]]; then
  TARGETS=("${!ALL_IMAGES[@]}")
fi

# ── Build a single image ────────────────────────────────────────────
build_one() {
  local name="$1"
  local IFS='|' parts
  read -ra parts <<< "${ALL_IMAGES[$name]}"
  local env_file="${parts[0]}"
  local env_name="${parts[1]}"
  local tag="${parts[2]}:latest"

  echo "========================================"
  echo "Building $name → $tag"
  echo "  env:  $env_file"
  echo "  name: $env_name"
  echo "========================================"

  cd "$REPO_ROOT"
  $DOCKER build \
    $NO_CACHE \
    -f "$DOCKERFILE" \
    --build-arg "ENV_FILE=$env_file" \
    --build-arg "ENV_NAME=$env_name" \
    -t "$tag" \
    .

  if [[ -n "${PUSH:-}" ]]; then
    echo "Pushing $tag ..."
    $DOCKER push "$tag"
  fi

  echo "Done: $name"
  echo ""
}

# ── Main ────────────────────────────────────────────────────────────
echo "Will build ${#TARGETS[@]} image(s): ${TARGETS[*]}"
echo "Repository root: $REPO_ROOT"
echo ""

if [[ "$PARALLEL" -gt 0 ]]; then
  # Parallel builds via background jobs
  running=0
  for name in "${TARGETS[@]}"; do
    build_one "$name" &
    running=$((running + 1))
    if [[ $running -ge $PARALLEL ]]; then
      wait -n 2>/dev/null || true
      running=$((running - 1))
    fi
  done
  wait
else
  for name in "${TARGETS[@]}"; do
    build_one "$name"
  done
fi

echo "========================================"
echo "All builds completed."
echo "========================================"
