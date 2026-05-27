#!/usr/bin/env bash
# ─── sinGRAAL — Script de déploiement cloud pour Bitcoin Puzzle #135 ──────────
#
# Usage : ./deploy.sh [runpod|vast|lambda|local] [N_INSTANCES]
#
# NOTE : Remplacer TARGET_X et TARGET_Y par les vraies coordonnées du puzzle #135
#        avant le déploiement.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

# ─── Coordonnées du puzzle #135 ───────────────────────────────────────────────
# TODO: Remplacer ces valeurs par les vraies coordonnées du point cible puzzle #135
TARGET_X="${PUZZLE135_TARGET_X:-PUZZLE135_X}"
TARGET_Y="${PUZZLE135_TARGET_Y:-PUZZLE135_Y}"
RANGE_BITS=135
BLOCK_BITS=20
THREADS=0

# Image Docker
IMAGE_NAME="singraal-p135"
IMAGE_TAG="latest"
FULL_IMAGE="${IMAGE_NAME}:${IMAGE_TAG}"

# ─── Aide ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage : $0 [runpod|vast|lambda|local] [N_INSTANCES]

  runpod  N   — Affiche les commandes RunPod CLI pour déployer N pods GPU
  vast    N   — Affiche les commandes vast.ai CLI pour louer N instances GPU
  lambda  N   — Affiche les commandes Lambda Labs CLI pour N instances GPU
  local   N   — Lance N instances locales avec docker-compose

Exemples :
  $0 local 4             # 4 instances locales
  $0 runpod 8            # commandes pour 8 pods RunPod
  $0 vast 16             # commandes pour 16 instances vast.ai

Variables d'environnement :
  PUZZLE135_TARGET_X   — Coordonnée x du point cible (hex 64 chars)
  PUZZLE135_TARGET_Y   — Coordonnée y du point cible (hex 64 chars)
EOF
    exit 1
}

# ─── Validation des paramètres ────────────────────────────────────────────────
PROVIDER="${1:-local}"
N="${2:-1}"

if [[ "${TARGET_X}" == "PUZZLE135_X" || "${TARGET_Y}" == "PUZZLE135_Y" ]]; then
    echo "[AVERTISSEMENT] TARGET_X/TARGET_Y non configurés."
    echo "  Définissez PUZZLE135_TARGET_X et PUZZLE135_TARGET_Y avant de déployer."
    echo "  Le déploiement continuera avec les placeholders pour génération des commandes."
    echo ""
fi

# ─── Build de l'image Docker locale ──────────────────────────────────────────
build_image() {
    echo "[INFO] Construction de l'image Docker ${FULL_IMAGE}..."
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
    docker build \
        -f "${SCRIPT_DIR}/Dockerfile" \
        -t "${FULL_IMAGE}" \
        "${REPO_ROOT}"
    echo "[OK] Image construite : ${FULL_IMAGE}"
}

# ─── Déploiement local ────────────────────────────────────────────────────────
deploy_local() {
    N_INST="${1:-1}"
    echo "[INFO] Déploiement local de ${N_INST} instance(s) avec docker-compose..."
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    TARGET_X="${TARGET_X}" \
    TARGET_Y="${TARGET_Y}" \
    RANGE_BITS="${RANGE_BITS}" \
    BLOCK_BITS="${BLOCK_BITS}" \
    THREADS="${THREADS}" \
    docker-compose -f "${SCRIPT_DIR}/docker-compose.yml" \
        up --build --scale solver="${N_INST}" -d

    echo "[OK] ${N_INST} instance(s) démarrée(s)."
    echo "     Logs : docker-compose -f ${SCRIPT_DIR}/docker-compose.yml logs -f"
    echo "     Résultats : docker volume inspect p135_solver_data"
}

# ─── RunPod ───────────────────────────────────────────────────────────────────
deploy_runpod() {
    N_PODS="${1:-1}"
    cat <<EOF
# ─── Commandes RunPod CLI pour déployer ${N_PODS} pod(s) ──────────────────────
# Pré-requis : runpodctl installé et configuré (runpodctl config --apiKey <KEY>)
# Documentation : https://docs.runpod.io/runpodctl/

# 1. Pusher l'image vers Docker Hub (ou registry privé)
docker tag ${FULL_IMAGE} <VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}
docker push <VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}

# 2. Créer ${N_PODS} pod(s) GPU (exemple avec RTX 4090, 24 GB VRAM)
EOF

    for i in $(seq 1 "${N_PODS}"); do
        cat <<EOF
runpodctl create pods \\
    --name "singraal-p135-${i}" \\
    --imageName "<VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}" \\
    --gpuType "NVIDIA GeForce RTX 4090" \\
    --containerDiskSize 20 \\
    --volumeSize 10 \\
    --env "TARGET_X=${TARGET_X}" \\
    --env "TARGET_Y=${TARGET_Y}" \\
    --env "RANGE_BITS=${RANGE_BITS}" \\
    --env "BLOCK_BITS=${BLOCK_BITS}" \\
    --env "THREADS=0"

EOF
    done

    cat <<EOF
# 3. Surveiller les pods
runpodctl get pods

# 4. Voir les logs d'un pod
runpodctl logs <POD_ID>

# 5. Arrêter tous les pods
runpodctl remove pods <POD_ID_1> <POD_ID_2> ...
EOF
}

# ─── vast.ai ──────────────────────────────────────────────────────────────────
deploy_vast() {
    N_INST="${1:-1}"
    cat <<EOF
# ─── Commandes vast.ai pour déployer ${N_INST} instance(s) ────────────────────
# Pré-requis : vastai CLI installé (pip install vastai) et clé API configurée
# Documentation : https://vast.ai/docs/

# 1. Pusher l'image vers Docker Hub
docker tag ${FULL_IMAGE} <VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}
docker push <VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}

# 2. Chercher les offres GPU disponibles (RTX 4090, ~24 GB)
vastai search offers 'reliability > 0.99 num_gpus=1 gpu_name=RTX_4090' --order dph_total

# 3. Louer ${N_INST} instance(s) (remplacer OFFER_ID par l'ID trouvé ci-dessus)
EOF

    for i in $(seq 1 "${N_INST}"); do
        cat <<EOF
vastai create instance OFFER_ID \\
    --image "<VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}" \\
    --disk 20 \\
    --env "-e TARGET_X=${TARGET_X} -e TARGET_Y=${TARGET_Y} -e RANGE_BITS=${RANGE_BITS} -e BLOCK_BITS=${BLOCK_BITS} -e THREADS=0" \\
    # Instance ${i}/${N_INST}

EOF
    done

    cat <<EOF
# 4. Lister les instances actives
vastai show instances

# 5. Voir les logs d'une instance
vastai logs <INSTANCE_ID>

# 6. Détruire toutes les instances
vastai destroy instance <INSTANCE_ID>
EOF
}

# ─── Lambda Labs ───────────────────────────────────────────────────────────────
deploy_lambda() {
    N_INST="${1:-1}"
    cat <<EOF
# ─── Commandes Lambda Labs pour déployer ${N_INST} instance(s) ────────────────
# Pré-requis : lambdacli installé et clé API configurée
# Documentation : https://docs.lambdalabs.com/cloud/

# 1. Pusher l'image vers Docker Hub ou Lambda Container Registry
docker tag ${FULL_IMAGE} <VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}
docker push <VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}

# 2. Lister les types d'instances GPU disponibles
lambda instances list-types

# 3. Lancer ${N_INST} instance(s) (exemple avec gpu_1x_a100_sxm4)
EOF

    for i in $(seq 1 "${N_INST}"); do
        cat <<EOF
lambda instances launch \\
    --instance-type gpu_1x_a100_sxm4 \\
    --name "singraal-p135-${i}" \\
    --ssh-key-names <VOTRE_SSH_KEY_NAME>
    # Note : Lambda Labs utilise des instances bare-metal + SSH
    # Après le lancement, se connecter par SSH et lancer Docker manuellement :
    # ssh ubuntu@<IP>
    # docker run -d \\
    #   -e TARGET_X=${TARGET_X} \\
    #   -e TARGET_Y=${TARGET_Y} \\
    #   -e RANGE_BITS=${RANGE_BITS} \\
    #   -e BLOCK_BITS=${BLOCK_BITS} \\
    #   -v /data:/data \\
    #   <VOTRE_DOCKERHUB_USERNAME>/${IMAGE_NAME}:${IMAGE_TAG}

EOF
    done
}

# ─── Main ─────────────────────────────────────────────────────────────────────
case "${PROVIDER}" in
    local)
        build_image
        deploy_local "${N}"
        ;;
    runpod)
        deploy_runpod "${N}"
        ;;
    vast)
        deploy_vast "${N}"
        ;;
    lambda)
        deploy_lambda "${N}"
        ;;
    help|-h|--help)
        usage
        ;;
    *)
        echo "[ERREUR] Provider inconnu : ${PROVIDER}" >&2
        usage
        ;;
esac
