# ─── Stage 1 : Build ──────────────────────────────────────────────────────────
FROM rust:latest AS builder

# Dépendances système pour Rayon (libc, pthread) et compilation générale
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copier le code source bsgs2d
COPY bsgs2d/ ./bsgs2d/

# Compiler en mode release avec optimisations maximales
WORKDIR /build/bsgs2d
RUN cargo build --release

# ─── Stage 2 : Runtime ────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Créer répertoire de données pour les résultats
RUN mkdir -p /data

# Copier le binaire compilé
COPY --from=builder /build/bsgs2d/target/release/bsgs2d /usr/local/bin/bsgs2d

# Copier le script d'entrée
COPY p135/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# ─── Variables d'environnement par défaut ─────────────────────────────────────
# NOTE: Remplacer TARGET_X et TARGET_Y par les vraies coordonnées du puzzle #135
ENV TARGET_X=PUZZLE135_X
ENV TARGET_Y=PUZZLE135_Y
ENV RANGE_BITS=135
ENV BLOCK_BITS=20
ENV THREADS=0

VOLUME ["/data"]

ENTRYPOINT ["/entrypoint.sh"]
