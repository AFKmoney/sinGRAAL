#!/usr/bin/env bash
set -euo pipefail

# ─── Bannière ────────────────────────────────────────────────────────────────
cat <<'BANNER'
╔══════════════════════════════════════════════════════════════════════════════╗
║          sinGRAAL — Solveur Bitcoin Puzzle #135  (ECDLP secp256k1)         ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Innovations actives :                                                       ║
║   #1  GLV décomposition 4D : k = k₁ + λk₂ + (1+λ)k₃ + (1-λ)k₄ (mod n)   ║
║   #2  6-automorphismes secp256k1 : φ(P)=(β·x,y), ×6 couverture par tile   ║
║   #3  Semaev S₃ : polynôme trivarié S₃(xA,xB,xP)=0, matrice Macaulay 28×28║
║   #4  LLL Kill-Switch certifié m=3 : λ₁²≥p⁶/28 → tile tuée en O(µs)      ║
║   #5  GSDD : Galois Symmetry + Nested Field Decomposition (Frobenius+CRT)  ║
║   #6  AnchorTable L2 : points pré-calculés en cache CPU L2                 ║
║   #7  Rayon parallel : tous les cœurs CPU (work-stealing)                  ║
║   #8  Marche bidirectionnelle tame↑/wild↓                                  ║
║   #9  Distribution Halton LDS (bases 2,3,5,7) — équidistribution optimale  ║
║   #10 Détection DP hiérarchique (hard/easy) avec table locale GPU          ║
╚══════════════════════════════════════════════════════════════════════════════╝
BANNER

# ─── Mode selftest ────────────────────────────────────────────────────────────
if [[ "${TARGET_X:-}" == "selftest" ]]; then
    echo "[INFO] Mode selftest activé (range_bits=40)"
    exec bsgs2d --selftest --range-bits 40
fi

# ─── Mode gsdd-selftest ───────────────────────────────────────────────────────
if [[ "${TARGET_X:-}" == "gsdd-selftest" ]]; then
    BLOCK="${BLOCK_BITS:-10}"
    echo "[INFO] Mode gsdd-selftest activé (range_bits=40, block_bits=${BLOCK})"
    exec bsgs2d --gsdd-selftest --range-bits 40 --block-bits "${BLOCK}"
fi

# ─── Validation des paramètres obligatoires ────────────────────────────────────
if [[ -z "${TARGET_X:-}" || "${TARGET_X}" == "PUZZLE135_X" ]]; then
    echo "[ERREUR] TARGET_X n'est pas défini ou contient la valeur placeholder." >&2
    echo "         Définissez TARGET_X avec la coordonnée x du point cible du puzzle #135." >&2
    echo "         Exemple : docker run -e TARGET_X=<hex64> -e TARGET_Y=<hex64> ..." >&2
    exit 1
fi

if [[ -z "${TARGET_Y:-}" || "${TARGET_Y}" == "PUZZLE135_Y" ]]; then
    echo "[ERREUR] TARGET_Y n'est pas défini ou contient la valeur placeholder." >&2
    echo "         Définissez TARGET_Y avec la coordonnée y du point cible du puzzle #135." >&2
    echo "         Exemple : docker run -e TARGET_X=<hex64> -e TARGET_Y=<hex64> ..." >&2
    exit 1
fi

# ─── Paramètres ──────────────────────────────────────────────────────────────
RANGE="${RANGE_BITS:-135}"
BLOCK="${BLOCK_BITS:-20}"
THREADS="${THREADS:-0}"

echo "[INFO] Démarrage du solveur"
echo "       TARGET_X    = ${TARGET_X}"
echo "       TARGET_Y    = ${TARGET_Y}"
echo "       RANGE_BITS  = ${RANGE}"
echo "       BLOCK_BITS  = ${BLOCK}"
echo "       THREADS     = ${THREADS} (0 = tous les cœurs)"
echo ""

# ─── Lancement du solveur ────────────────────────────────────────────────────
RESULT_FILE="/data/solution.txt"

if bsgs2d \
    --solve \
    --target-x "${TARGET_X}" \
    --target-y "${TARGET_Y}" \
    --range-bits "${RANGE}" \
    --block-bits "${BLOCK}" \
    --threads "${THREADS}" \
    | tee /tmp/solver_output.txt; then

    # Chercher la clé dans la sortie
    KEY=$(grep -oP '(?<=k = )[0-9a-fA-F]+' /tmp/solver_output.txt | head -1 || true)
    if [[ -n "${KEY}" ]]; then
        echo ""
        echo "╔══════════════════════════════════════════════════════╗"
        echo "║  SOLUTION TROUVÉE                                    ║"
        echo "║  k = ${KEY}"
        echo "╚══════════════════════════════════════════════════════╝"
        {
            echo "# sinGRAAL — Bitcoin Puzzle #135 Solution"
            echo "# Date: $(date -u)"
            echo "TARGET_X=${TARGET_X}"
            echo "TARGET_Y=${TARGET_Y}"
            echo "RANGE_BITS=${RANGE}"
            echo "k=${KEY}"
        } > "${RESULT_FILE}"
        echo "[INFO] Solution sauvegardée dans ${RESULT_FILE}"
        exit 0
    fi
fi

echo "[INFO] Solveur terminé sans solution trouvée dans cette session."
exit 0
