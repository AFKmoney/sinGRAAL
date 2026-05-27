#!/usr/bin/env bash
set -euo pipefail

cat <<'BANNER'
╔══════════════════════════════════════════════════════════════════════════════╗
║          sinGRAAL — Solveur Bitcoin Puzzle #135  (ECDLP secp256k1)         ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Pipeline complet — toutes innovations actives :                             ║
║   #1  GLV 4D           k = k₁+λk₂+(1+λ)k₃+(1-λ)k₄ (mod n)                ║
║   #2  Semaev S₃        polynôme trivarié → matrice Macaulay 28×28 (m=3)    ║
║   #3  LLL Kill-Switch  λ₁²≥p⁶/28 → tile tuée en O(µs), >99.99% kill-rate  ║
║   #4  6-automorphismes φ(P)=(β·x,y), ×6 couverture k-space par tile        ║
║   #5  GSDD             Frobenius + CRT sur facteurs de (n-1)                ║
║   #6  AnchorTable L2   points pré-calculés en cache CPU L1/L2               ║
║   #7  Rayon parallel   work-stealing sur tous les cœurs CPU                 ║
║   #8  Tame/Wild bidi   convergence dirigée tame↑ wild↓                      ║
║   #9  Halton LDS       équidistribution des sauts (bases 2,3,5,7)           ║
║   #10 DP hiérarchique  hard/easy distinguished points                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
BANNER

# ─── Mode selftest ────────────────────────────────────────────────────────────
if [[ "${TARGET_X:-}" == "selftest" ]]; then
    echo "[INFO] Selftest automatique range_bits=40"
    exec bsgs2d --selftest --range-bits 40
fi

if [[ "${TARGET_X:-}" == "gsdd-selftest" ]]; then
    BLOCK="${BLOCK_BITS:-10}"
    echo "[INFO] GSDD selftest range_bits=40 block_bits=${BLOCK}"
    exec bsgs2d --gsdd-selftest --range-bits 40 --block-bits "${BLOCK}"
fi

# ─── Validation TARGET_X / TARGET_Y ──────────────────────────────────────────
if [[ -z "${TARGET_X:-}" || "${TARGET_X}" == "PUZZLE135_X" ]]; then
    echo "[ERREUR] TARGET_X non défini. Passez -e TARGET_X=<hex64> au conteneur." >&2
    exit 1
fi
if [[ -z "${TARGET_Y:-}" || "${TARGET_Y}" == "PUZZLE135_Y" ]]; then
    echo "[ERREUR] TARGET_Y non défini. Passez -e TARGET_Y=<hex64> au conteneur." >&2
    exit 1
fi

# ─── Paramètres (avec valeurs par défaut) ────────────────────────────────────
RANGE="${RANGE_BITS:-135}"
BLOCK="${BLOCK_BITS:-20}"
HALF="${HALF_BITS:-0}"          # 0 = auto (ceil(range/2)+2)
MLEVEL="${M_LEVEL:-3}"          # niveau Jochemsz-May : 3 = dim=28 (optimal)
THREADS="${THREADS:-0}"         # 0 = tous les cœurs
MODE="${SOLVER_MODE:-solve}"    # solve | optimized | parallel | gsdd

# Construire --half-bits seulement si spécifié
HALF_ARG=""
if [[ "${HALF}" != "0" ]]; then
    HALF_ARG="--half-bits ${HALF}"
fi

echo "[INFO] Paramètres :"
echo "       TARGET_X     = ${TARGET_X}"
echo "       TARGET_Y     = ${TARGET_Y}"
echo "       RANGE_BITS   = ${RANGE}"
echo "       BLOCK_BITS   = ${BLOCK}"
echo "       HALF_BITS    = ${HALF} (0=auto)"
echo "       M_LEVEL      = ${MLEVEL}  (dim=$(( MLEVEL==3 ? 28 : 15 )))"
echo "       THREADS      = ${THREADS} (0=tous les cœurs)"
echo "       SOLVER_MODE  = ${MODE}"
echo ""

# ─── Sélection du mode ────────────────────────────────────────────────────────
case "${MODE}" in
    solve|gsdd)
        # Full stack : GLV-4D + Semaev-S₃ + LLL-m=3 + 6-aut + AnchorL2 + Rayon + GSDD
        MODE_FLAG="--solve"
        ;;
    optimized)
        # Dispatcher optimisé : AnchorTable L2 + 6-aut + Rayon (sans GSDD complet)
        MODE_FLAG="--optimized"
        ;;
    parallel)
        # Dispatcher parallèle Rayon uniquement
        MODE_FLAG="--parallel"
        ;;
    dispatch)
        MODE_FLAG="--dispatch"
        ;;
    *)
        echo "[ERREUR] SOLVER_MODE invalide : ${MODE}. Valeurs : solve|optimized|parallel|dispatch" >&2
        exit 1
        ;;
esac

# ─── Lancement ────────────────────────────────────────────────────────────────
RESULT_FILE="/data/solution.txt"
mkdir -p /data

CMD="bsgs2d \
    ${MODE_FLAG} \
    --target-x ${TARGET_X} \
    --target-y ${TARGET_Y} \
    --range-bits ${RANGE} \
    --block-bits ${BLOCK} \
    --m-level ${MLEVEL} \
    --threads ${THREADS} \
    ${HALF_ARG}"

echo "[INFO] Commande : ${CMD}"
echo ""

if ${CMD} | tee /tmp/solver_out.txt; then
    KEY=$(grep -oP '(?<=k = )[0-9a-fA-F]+' /tmp/solver_out.txt | head -1 || true)
    if [[ -n "${KEY}" ]]; then
        echo ""
        echo "╔══════════════════════════════════════════════════════════════════╗"
        echo "║  SOLUTION TROUVÉE                                                ║"
        printf "║  k = %-60s║\n" "${KEY}"
        echo "╚══════════════════════════════════════════════════════════════════╝"
        {
            echo "# sinGRAAL — Bitcoin Puzzle #135 Solution"
            echo "# Date: $(date -u)"
            echo "TARGET_X=${TARGET_X}"
            echo "TARGET_Y=${TARGET_Y}"
            echo "RANGE_BITS=${RANGE}"
            echo "SOLVER_MODE=${MODE}"
            echo "k=${KEY}"
        } > "${RESULT_FILE}"
        echo "[INFO] Solution sauvegardée : ${RESULT_FILE}"
        exit 0
    fi
fi

echo "[INFO] Session terminée sans solution. Relancer pour continuer la recherche."
exit 0
