#!/usr/bin/env bash
# sinGRAAL entrypoint — puzzle #135 cloud GPU
# Toutes les options via variables d'environnement.
#
# Détection GPU automatique :
#   - Liste les GPU via nvidia-smi (nom, VRAM, compute capability).
#   - Si NUM_ANIMALS n'est pas fixé, le calibre selon la VRAM de la 1ʳᵉ carte.
#   Le binaire ayant été compilé en fatbinary multi-architecture, il s'exécute
#   sur n'importe quel modèle sans recompilation.
set -euo pipefail

BINARY=/usr/local/bin/singraal
ARGS=()

# ── Détection GPU ───────────────────────────────────────────────────────────────
detect_gpus() {
    if ! command -v nvidia-smi >/dev/null 2>&1; then
        echo "[singraal] nvidia-smi introuvable — exécution CPU/fallback"
        return
    fi
    echo "[singraal] GPU détectés :"
    nvidia-smi --query-gpu=index,name,memory.total,compute_cap \
        --format=csv,noheader 2>/dev/null \
        | while IFS=',' read -r idx name mem cc; do
            printf "  [%s]%s | %s | compute_cap=%s\n" \
                "$(echo "$idx" | xargs)" "$name" \
                "$(echo "$mem" | xargs)" "$(echo "$cc" | xargs)"
        done
}

# Calibre NUM_ANIMALS selon la VRAM (≈ MiB → animaux), si non fixé.
# Empirique : ~3 KiB d'état GPU par animal → on vise ~25% de la VRAM.
auto_num_animals() {
    if [[ -n "${NUM_ANIMALS:-}" ]]; then
        echo "$NUM_ANIMALS"
        return
    fi
    local mib
    mib=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null \
            | head -n1 | xargs || echo "")
    if [[ -z "$mib" ]]; then
        echo 262144
        return
    fi
    # animaux ≈ (VRAM_MiB * 0.25 * 1024) / 3 KiB, arrondi à la puissance de 2 inférieure.
    local cand=$(( mib * 1024 / 4 / 3 ))
    local p=1024
    while (( p * 2 <= cand )); do p=$(( p * 2 )); done
    # bornes raisonnables : 2^16 .. 2^20
    (( p < 65536 )) && p=65536
    (( p > 1048576 )) && p=1048576
    echo "$p"
}

detect_gpus
NUM_ANIMALS_RESOLVED=$(auto_num_animals)

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
    ARGS+=(--num-animals "$NUM_ANIMALS_RESOLVED")
    echo "[singraal] mode worker → coordinateur $COORDINATOR  (animals=$NUM_ANIMALS_RESOLVED)"
else
    ARGS+=(--all-gpus)
    ARGS+=(--num-animals "$NUM_ANIMALS_RESOLVED")
    [[ -n "${CHECKPOINT:-}" ]] && ARGS+=(--checkpoint "$CHECKPOINT")
    echo "[singraal] mode standalone — puzzle #135  (animals=$NUM_ANIMALS_RESOLVED)"
fi

echo "[singraal] cmd: $BINARY ${ARGS[*]}"
exec "$BINARY" "${ARGS[@]}"
