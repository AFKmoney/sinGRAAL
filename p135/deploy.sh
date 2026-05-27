#!/usr/bin/env bash
# sinGRAAL — Déploiement cloud Bitcoin Puzzle #135
# Usage : ./deploy.sh [local|runpod|vast|lambda] [N_INSTANCES]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

# ─── Coordonnées puzzle #135 ──────────────────────────────────────────────────
# Définir ces variables AVANT de lancer le script :
#   export TARGET_X="<hex64>"
#   export TARGET_Y="<hex64>"
TARGET_X="${TARGET_X:-PUZZLE135_X}"
TARGET_Y="${TARGET_Y:-PUZZLE135_Y}"

# ─── Paramètres solveur ───────────────────────────────────────────────────────
RANGE_BITS="${RANGE_BITS:-135}"
BLOCK_BITS="${BLOCK_BITS:-20}"
HALF_BITS="${HALF_BITS:-0}"
M_LEVEL="${M_LEVEL:-3}"
THREADS="${THREADS:-0}"
SOLVER_MODE="${SOLVER_MODE:-solve}"

IMAGE_NAME="singraal-p135"
IMAGE_TAG="latest"
FULL_IMAGE="${IMAGE_NAME}:${IMAGE_TAG}"

# ─── Aide ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage : $0 [local|runpod|vast|lambda] [N_INSTANCES]

  local   N   — Lance N instances locales (docker-compose)
  runpod  N   — Génère les commandes RunPod CLI pour N pods GPU
  vast    N   — Génère les commandes vast.ai CLI pour N instances GPU
  lambda  N   — Génère les commandes Lambda Labs SSH+Docker pour N instances

Variables d'environnement à définir avant le lancement :
  TARGET_X      — Coordonnée x du point cible puzzle #135 (hex 64 chars)
  TARGET_Y      — Coordonnée y du point cible puzzle #135 (hex 64 chars)
  RANGE_BITS    — Taille espace recherche (défaut: 135)
  BLOCK_BITS    — Taille des tiles : 2^BLOCK_BITS (défaut: 20, RAM ~16 GB)
  HALF_BITS     — Dimension grille GLV (défaut: 0=auto)
  M_LEVEL       — Niveau Jochemsz-May (défaut: 3 = matrice 28×28)
  THREADS       — Threads CPU (défaut: 0=tous)
  SOLVER_MODE   — Mode : solve|optimized|parallel (défaut: solve)

Exemples :
  export TARGET_X="abc123..." TARGET_Y="def456..."
  $0 local 4          # 4 instances locales
  $0 runpod 8         # commandes pour 8 pods RunPod
  $0 vast 16          # commandes pour 16 instances vast.ai
EOF
    exit 1
}

# ─── Vérification coordonnées ─────────────────────────────────────────────────
check_coords() {
    if [[ "${TARGET_X}" == "PUZZLE135_X" || "${TARGET_Y}" == "PUZZLE135_Y" ]]; then
        echo "[AVERTISSEMENT] TARGET_X/Y non définis — les commandes générées contiennent des placeholders."
        echo "  Définissez les coordonnées avec : export TARGET_X=... TARGET_Y=..."
        echo ""
    fi
}

# ─── Build image locale ───────────────────────────────────────────────────────
build_image() {
    echo "[BUILD] Construction de l'image ${FULL_IMAGE}..."
    docker build \
        -f "${SCRIPT_DIR}/Dockerfile" \
        -t "${FULL_IMAGE}" \
        "${REPO_ROOT}"
    echo "[OK] Image prête : ${FULL_IMAGE}"
}

# ─── LOCAL ────────────────────────────────────────────────────────────────────
deploy_local() {
    N="${1:-1}"
    build_image
    echo "[INFO] Lancement de ${N} instance(s)..."
    TARGET_X="${TARGET_X}" \
    TARGET_Y="${TARGET_Y}" \
    RANGE_BITS="${RANGE_BITS}" \
    BLOCK_BITS="${BLOCK_BITS}" \
    HALF_BITS="${HALF_BITS}" \
    M_LEVEL="${M_LEVEL}" \
    THREADS="${THREADS}" \
    SOLVER_MODE="${SOLVER_MODE}" \
    docker-compose -f "${SCRIPT_DIR}/docker-compose.yml" \
        up --build --scale solver="${N}" -d
    echo "[OK] ${N} instance(s) démarrée(s)."
    echo "     Logs    : docker-compose -f ${SCRIPT_DIR}/docker-compose.yml logs -f"
    echo "     Résultat: docker-compose -f ${SCRIPT_DIR}/docker-compose.yml exec solver cat /data/solution.txt"
}

# ─── RUNPOD ───────────────────────────────────────────────────────────────────
deploy_runpod() {
    N="${1:-1}"
    cat <<EOF
# ════════════════════════════════════════════════════════════════
#  sinGRAAL — RunPod : ${N} pod(s) GPU
#  Pré-requis : runpodctl configuré (runpodctl config --apiKey <KEY>)
# ════════════════════════════════════════════════════════════════

# 1. Pusher l'image vers Docker Hub
docker tag ${FULL_IMAGE} <DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}
docker push <DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}

# 2. Créer ${N} pod(s) (RTX 4090 ou A100)
EOF
    for i in $(seq 1 "${N}"); do
        cat <<EOF
runpodctl create pods \\
    --name "singraal-p135-${i}" \\
    --imageName "<DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}" \\
    --gpuType "NVIDIA GeForce RTX 4090" \\
    --containerDiskSize 30 \\
    --volumeSize 10 \\
    --env "TARGET_X=${TARGET_X}" \\
    --env "TARGET_Y=${TARGET_Y}" \\
    --env "RANGE_BITS=${RANGE_BITS}" \\
    --env "BLOCK_BITS=${BLOCK_BITS}" \\
    --env "HALF_BITS=${HALF_BITS}" \\
    --env "M_LEVEL=${M_LEVEL}" \\
    --env "THREADS=0" \\
    --env "SOLVER_MODE=${SOLVER_MODE}"

EOF
    done
    cat <<EOF
# 3. Suivre les pods
runpodctl get pods

# 4. Logs d'un pod
runpodctl logs <POD_ID>
EOF
}

# ─── VAST.AI ──────────────────────────────────────────────────────────────────
deploy_vast() {
    N="${1:-1}"
    cat <<EOF
# ════════════════════════════════════════════════════════════════
#  sinGRAAL — vast.ai : ${N} instance(s)
#  Pré-requis : pip install vastai && vastai set api-key <KEY>
# ════════════════════════════════════════════════════════════════

# 1. Pusher l'image
docker tag ${FULL_IMAGE} <DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}
docker push <DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}

# 2. Chercher des offres GPU
vastai search offers 'reliability > 0.99 num_gpus=1 gpu_ram >= 20' --order dph_total

# 3. Louer ${N} instance(s) (remplacer OFFER_ID)
EOF
    for i in $(seq 1 "${N}"); do
        cat <<EOF
vastai create instance <OFFER_ID_${i}> \\
    --image "<DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}" \\
    --disk 30 \\
    --env "-e TARGET_X=${TARGET_X} -e TARGET_Y=${TARGET_Y} -e RANGE_BITS=${RANGE_BITS} -e BLOCK_BITS=${BLOCK_BITS} -e HALF_BITS=${HALF_BITS} -e M_LEVEL=${M_LEVEL} -e SOLVER_MODE=${SOLVER_MODE}"

EOF
    done
    cat <<EOF
# 4. Voir les instances
vastai show instances

# 5. Logs
vastai logs <INSTANCE_ID>
EOF
}

# ─── LAMBDA LABS ──────────────────────────────────────────────────────────────
deploy_lambda() {
    N="${1:-1}"
    cat <<EOF
# ════════════════════════════════════════════════════════════════
#  sinGRAAL — Lambda Labs : ${N} instance(s) A100/H100
#  Lambda Labs = instances bare-metal SSH, Docker manuel
# ════════════════════════════════════════════════════════════════

# 1. Pusher l'image
docker tag ${FULL_IMAGE} <DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}
docker push <DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}

# 2. Lancer les instances via Lambda dashboard ou CLI
#    puis pour CHAQUE instance, se connecter par SSH et exécuter :

EOF
    for i in $(seq 1 "${N}"); do
        cat <<EOF
# ── Instance ${i}/${N} ──
ssh ubuntu@<IP_INSTANCE_${i}> << 'ENDSSH'
docker run -d \\
    --name singraal-p135 \\
    -e TARGET_X="${TARGET_X}" \\
    -e TARGET_Y="${TARGET_Y}" \\
    -e RANGE_BITS="${RANGE_BITS}" \\
    -e BLOCK_BITS="${BLOCK_BITS}" \\
    -e HALF_BITS="${HALF_BITS}" \\
    -e M_LEVEL="${M_LEVEL}" \\
    -e SOLVER_MODE="${SOLVER_MODE}" \\
    -v /home/ubuntu/results:/data \\
    <DOCKERHUB_USER>/${IMAGE_NAME}:${IMAGE_TAG}
docker logs -f singraal-p135
ENDSSH

EOF
    done
}

# ─── MAIN ─────────────────────────────────────────────────────────────────────
PROVIDER="${1:-help}"
N="${2:-1}"

check_coords

case "${PROVIDER}" in
    local)   deploy_local "${N}" ;;
    runpod)  deploy_runpod "${N}" ;;
    vast)    deploy_vast "${N}" ;;
    lambda)  deploy_lambda "${N}" ;;
    help|-h|--help) usage ;;
    *)
        echo "[ERREUR] Provider inconnu : ${PROVIDER}" >&2
        usage
        ;;
esac
