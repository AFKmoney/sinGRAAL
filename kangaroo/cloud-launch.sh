#!/usr/bin/env bash
# sinGRAAL — Cloud GPU Launch Script
# ====================================
# Usage: ./cloud-launch.sh [OPTIONS]
#
# Quick start:
#   ./cloud-launch.sh --target-x <64hex> --target-y <64hex>
#
# Distributed (run on each worker machine):
#   COORDINATOR=192.168.1.100:5135 ./cloud-launch.sh

set -euo pipefail

# ─── Puzzle #135 target (public) ─────────────────────────────────────────────
# Source: https://privatekeys.pw/puzzles/bitcoin-puzzle-tx
TARGET_X="${TARGET_X:-}"
TARGET_Y="${TARGET_Y:-}"
RANGE_BITS="${RANGE_BITS:-135}"
COORDINATOR="${COORDINATOR:-}"
CHECKPOINT="${CHECKPOINT:-/data/singraal-checkpoint.bin}"
BIND_ADDR="${BIND_ADDR:-0.0.0.0:5135}"

BINARY="./target/release/kangaroo"

# ─── Detect GPU architecture ─────────────────────────────────────────────────
detect_arch() {
    if command -v nvidia-smi &>/dev/null; then
        local cc
        cc=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')
        case "$cc" in
            80|86|87) echo "sm_${cc}" ;;   # A100, RTX 3090, H100
            89|90)    echo "sm_${cc}" ;;   # RTX 4090, H100
            75)       echo "sm_75"    ;;   # RTX 2080
            *)        echo "sm_80"    ;;   # safe default
        esac
    else
        echo "sm_80"
    fi
}

# ─── Build ────────────────────────────────────────────────────────────────────
build() {
    local arch
    arch=$(detect_arch)
    echo "Building sinGRAAL v13 for CUDA arch ${arch}..."
    CUDA_ARCH="${arch}" cargo build --release --features cuda
    echo "Build complete: $BINARY"
}

# ─── Verify research modules ──────────────────────────────────────────────────
verify() {
    echo "=== Running research verification ==="
    $BINARY --research4d --range-bits 135
    echo ""
    echo "=== Empirical C measurement (48-bit, 100 trials) ==="
    $BINARY --benchmark-c --range-bits 48 --trials 100
    echo ""
    echo "=== GLS 4D foundation ==="
    $BINARY --gls4d --range-bits 64
}

# ─── Coordinator mode ─────────────────────────────────────────────────────────
run_coordinator() {
    echo "Starting sinGRAAL COORDINATOR on ${BIND_ADDR}"
    echo "Workers connect with: --coordinator ${BIND_ADDR}"
    exec $BINARY \
        --serve \
        --bind "${BIND_ADDR}" \
        --target-x "${TARGET_X}" \
        --target-y "${TARGET_Y}" \
        --range-bits "${RANGE_BITS}"
}

# ─── Worker mode ──────────────────────────────────────────────────────────────
run_worker() {
    echo "Starting sinGRAAL WORKER → coordinator ${COORDINATOR}"
    exec $BINARY \
        --coordinator "${COORDINATOR}" \
        --all-gpus \
        --range-bits "${RANGE_BITS}" \
        --checkpoint "${CHECKPOINT}"
}

# ─── Standalone mode ─────────────────────────────────────────────────────────
run_standalone() {
    [[ -z "$TARGET_X" ]] && { echo "ERROR: set TARGET_X"; exit 1; }
    [[ -z "$TARGET_Y" ]] && { echo "ERROR: set TARGET_Y"; exit 1; }
    echo "Starting sinGRAAL STANDALONE — all local GPUs"
    exec $BINARY \
        --target-x "${TARGET_X}" \
        --target-y "${TARGET_Y}" \
        --range-bits "${RANGE_BITS}" \
        --all-gpus \
        --checkpoint "${CHECKPOINT}"
}

# ─── Main ─────────────────────────────────────────────────────────────────────
case "${1:-launch}" in
    build)       build ;;
    verify)      verify ;;
    coordinator) run_coordinator ;;
    worker)      run_worker ;;
    launch|*)
        build
        if [[ -n "$COORDINATOR" ]]; then
            run_worker
        elif [[ "${SERVE:-}" == "1" ]]; then
            run_coordinator
        else
            run_standalone
        fi
        ;;
esac
