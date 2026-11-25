#!/usr/bin/env bash
# =============================================================================
# Build the benchmark images: one per pixi environment used by a method
# (methods_registry.yaml -> execution.pixi_env), from pixi.toml / pixi.lock.
#
# Usage: containers/build_all.sh [options] [target ...]
#   targets            pixi env names (default: all), e.g. flecs scnode flecs-cuda
#   --platform P       linux/arm64 | linux/amd64 | native (default: native)
#   --cuda             also build the <env>-cuda images (linux/amd64 only)
#   --no-cache         build without layer cache
#   --push             push after building (multi-platform builds need --push)
#   -j N               N builds in parallel (default 1)
#   --list             list targets and exit
# Works with macOS /bin/bash 3.2. Exits non-zero if any build failed.
# =============================================================================
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOCKERFILE="containers/docker/Dockerfile"
DOCKER="${DOCKER:-docker}"
PLATFORM="native"; WITH_CUDA=0; NO_CACHE=""; PUSH=0; PARALLEL=1; LIST=0
REQUESTED=()

while [ $# -gt 0 ]; do
  case "$1" in
    --platform) PLATFORM="${2:?}"; shift ;;
    --cuda) WITH_CUDA=1 ;;
    --no-cache) NO_CACHE="--no-cache" ;;
    --push) PUSH=1 ;;
    -j) PARALLEL="${2:?}"; shift ;;
    --list) LIST=1 ;;
    -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
    *) REQUESTED+=("$1") ;;
  esac
  shift
done

cd "$REPO_ROOT"

# Python with PyYAML: the pixi "runner" env, else a conda "benchmark_runner"
# env, else any python3 that has PyYAML.
registry() {
  local out
  if command -v pixi >/dev/null 2>&1; then
    out=$( ( eval "$(pixi shell-hook -e runner -s bash 2>/dev/null)" && python utils/methods_registry.py "$@" ) 2>/dev/null ) \
      && [ -n "$out" ] && { echo "$out"; return 0; }
  fi
  if command -v conda >/dev/null 2>&1; then
    out=$(conda run -n benchmark_runner python utils/methods_registry.py "$@" 2>/dev/null) \
      && [ -n "$out" ] && { echo "$out"; return 0; }
  fi
  python3 utils/methods_registry.py "$@"
}

# Lines: target|image|platforms(comma-separated)
TARGET_LINES=$(registry container-images | python3 -c '
import json, sys
for d in json.load(sys.stdin):
    print("|".join([d["target"], d["image"], ",".join(d["platforms"])]))
') || { echo "ERROR: cannot read targets from methods_registry.yaml" >&2; exit 1; }

if [ "$LIST" = 1 ]; then
  echo "$TARGET_LINES" | awk -F'|' '{printf "  %-24s %s:latest  [%s]\n", $1, $2, $3}'
  exit 0
fi

if [ "$PLATFORM" = "native" ]; then
  case "$($DOCKER info --format '{{.Architecture}}' 2>/dev/null)" in
    aarch64|arm64) PLATFORM="linux/arm64" ;;
    *) PLATFORM="linux/amd64" ;;
  esac
fi

# Select targets
SELECTED=()
while IFS='|' read -r target image platforms; do
  [ -n "$target" ] || continue
  if [ ${#REQUESTED[@]} -gt 0 ]; then
    case " ${REQUESTED[*]} " in *" $target "*) ;; *) continue ;; esac
  else
    case "$target" in *-cuda) [ "$WITH_CUDA" = 1 ] || continue ;; esac
  fi
  case "$target" in *-cuda) plat="linux/amd64" ;; *) plat="$PLATFORM" ;; esac
  SELECTED+=("$target|$image|$plat")
done <<< "$TARGET_LINES"

if [ ${#SELECTED[@]} -eq 0 ]; then
  echo "No target selected (see --list)." >&2; exit 1
fi

LOG_DIR="$(mktemp -d "${TMPDIR:-/tmp}/trajgrn-build.XXXXXX")"
build_one() {  # target|image|platform
  local target image plat
  IFS='|' read -r target image plat <<< "$1"
  local log="$LOG_DIR/$target.log" out=(--load)
  [ "$PUSH" = 1 ] && out=(--push)
  echo "[build] $target -> $image:latest ($plat)  log: $log"
  if $DOCKER buildx build $NO_CACHE -f "$DOCKERFILE" \
       --build-arg "PIXI_ENV=$target" --platform "$plat" \
       "${out[@]}" -t "$image:latest" . > "$log" 2>&1; then
    echo "ok" > "$LOG_DIR/$target.status"; echo "[done]  $target"
  else
    echo "failed" > "$LOG_DIR/$target.status"; echo "[FAIL]  $target (see $log)"
  fi
}

running=0
for spec in "${SELECTED[@]}"; do
  build_one "$spec" &
  running=$((running + 1))
  if [ "$running" -ge "$PARALLEL" ]; then wait; running=0; fi
done
wait

failed=0
echo "========================================"
for spec in "${SELECTED[@]}"; do
  t="${spec%%|*}"; st=$(cat "$LOG_DIR/$t.status" 2>/dev/null || echo failed)
  printf "  %-24s %s\n" "$t" "$st"
  [ "$st" = ok ] || failed=$((failed + 1))
done
echo "========================================"
[ "$failed" -eq 0 ] || { echo "$failed build(s) failed; logs in $LOG_DIR" >&2; exit 1; }
