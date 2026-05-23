#!/usr/bin/env bash
# sinGRAAL v12 — cloud-ready entrypoint
#
# Environment variables (set in cloud provider UI):
#
#   TARGET_X       64-hex-char x-coordinate of target EC point (required)
#   TARGET_Y       64-hex-char y-coordinate of target EC point (required)
#   RANGE_BITS     search range in bits (default: 135)
#   NUM_ANIMALS    kangaroos per GPU (default: 262144)
#   DP_BITS        distinguished-point difficulty override (default: auto)
#
#   COORDINATOR    HOST:PORT of pool coordinator — enables worker mode
#                  example: COORDINATOR=10.0.0.1:5135
#
#   SERVE          set to "1" to run as pool coordinator (binds :5135)
#   BIND           coordinator bind address (default: 0.0.0.0:5135)
#
#   ALL_GPUS       set to "1" to use all GPUs (default: 1 GPU, device 0)
#   DEVICE         specific CUDA device index (default: 0, ignored if ALL_GPUS=1)
#
#   NO_CHECKPOINT  set to "1" to disable checkpoint save/load
#   CHECKPOINT     checkpoint file path (default: /data/kangaroo.ckpt)
#
# If no env vars are set and arguments are provided, they are passed directly.
#
# Examples (RunPod / vast.ai / Lambda Labs):
#
#   # Worker connecting to coordinator at 10.0.0.1:5135:
#   docker run --gpus all \
#     -e TARGET_X=... -e TARGET_Y=... \
#     -e COORDINATOR=10.0.0.1:5135 \
#     -e ALL_GPUS=1 \
#     singraal:v12-4090
#
#   # Coordinator (expose port 5135):
#   docker run --gpus all -p 5135:5135 \
#     -e TARGET_X=... -e TARGET_Y=... \
#     -e SERVE=1 \
#     singraal:v12-4090
#
#   # Solo run (no coordinator):
#   docker run --gpus all \
#     -e TARGET_X=... -e TARGET_Y=... \
#     -e ALL_GPUS=1 \
#     singraal:v12-4090

set -euo pipefail

KANGAROO=/usr/local/bin/kangaroo

# If arguments are passed directly, run with them (bypass env-var logic)
if [[ $# -gt 0 ]]; then
    exec "$KANGAROO" "$@"
fi

# ─── Build argument list from environment ────────────────────────────────────

ARGS=()

# Target (required unless --research / --benchmark-c)
if [[ -n "${TARGET_X:-}" ]] && [[ -n "${TARGET_Y:-}" ]]; then
    ARGS+=("--target-x" "$TARGET_X" "--target-y" "$TARGET_Y")
elif [[ -z "${SERVE:-}" ]] && [[ -z "${RESEARCH:-}" ]]; then
    echo "[sinGRAAL] ERROR: TARGET_X and TARGET_Y must be set (or pass args directly)"
    echo "           Set SERVE=1 for coordinator mode without immediate target."
    echo "           Run with --help to see all flags."
    exec "$KANGAROO" --help
fi

# Range
ARGS+=("--range-bits" "${RANGE_BITS:-135}")

# Animals per GPU
ARGS+=("--num-animals" "${NUM_ANIMALS:-262144}")

# DP bits (optional override)
if [[ -n "${DP_BITS:-}" ]]; then
    ARGS+=("--dp-bits" "$DP_BITS")
fi

# GPU selection
if [[ "${ALL_GPUS:-0}" == "1" ]]; then
    ARGS+=("--all-gpus")
else
    ARGS+=("--device" "${DEVICE:-0}")
fi

# Checkpoint
if [[ "${NO_CHECKPOINT:-0}" == "1" ]]; then
    ARGS+=("--no-checkpoint")
else
    ARGS+=("--checkpoint" "${CHECKPOINT:-/data/kangaroo.ckpt}")
fi

# Mode selection
if [[ "${SERVE:-0}" == "1" ]]; then
    ARGS+=("--serve" "--bind" "${BIND:-0.0.0.0:5135}")
elif [[ -n "${COORDINATOR:-}" ]]; then
    ARGS+=("--coordinator" "$COORDINATOR")
fi

# Research modes (for debugging/benchmarking)
if [[ "${RESEARCH:-0}" == "1" ]]; then
    ARGS+=("--research")
fi
if [[ "${BENCHMARK_C:-0}" == "1" ]]; then
    ARGS+=("--benchmark-c")
    ARGS+=("--trials" "${TRIALS:-200}")
fi

echo "[sinGRAAL] v12 starting — $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo 'CPU mode')"
echo "[sinGRAAL] args: ${ARGS[*]}"

exec "$KANGAROO" "${ARGS[@]}"
