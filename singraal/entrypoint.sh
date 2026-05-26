#!/usr/bin/env bash
# sinGRAAL entrypoint — puzzle #135 cloud GPU
# Toutes les options via variables d'environnement.
set -euo pipefail

BINARY=/usr/local/bin/singraal
ARGS=()

# ── Puzzle cible ──────────────────────────────────────────────────────────────
ARGS+=(--target-x "${TARGET_X:?TARGET_X required}")
ARGS+=(--target-y "${TARGET_Y:?TARGET_Y required}")
ARGS+=(--range-bits "${RANGE_BITS:-135}")

# ── Mode : coordinateur, worker, ou standalone ────────────────────────────────
if [[ "${SERVE:-0}" == "1" ]]; then
    ARGS+=(--serve)
    ARGS+=(--bind "${BIND:-0.0.0.0:5135}")
    [[ -n "${CHECKPOINT:-}" ]] && ARGS+=(--checkpoint "$CHECKPOINT")
    echo "[singraal] mode coordinateur — bind ${BIND:-0.0.0.0:5135}"
elif [[ -n "${COORDINATOR:-}" ]]; then
    ARGS+=(--coordinator "$COORDINATOR")
    ARGS+=(--all-gpus)
    ARGS+=(--num-animals "${NUM_ANIMALS:-262144}")
    echo "[singraal] mode worker → coordinateur $COORDINATOR"
else
    ARGS+=(--all-gpus)
    ARGS+=(--num-animals "${NUM_ANIMALS:-262144}")
    [[ -n "${CHECKPOINT:-}" ]] && ARGS+=(--checkpoint "$CHECKPOINT")
    echo "[singraal] mode standalone — puzzle #135"
fi

echo "[singraal] cmd: $BINARY ${ARGS[*]}"
exec "$BINARY" "${ARGS[@]}"
