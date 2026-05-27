# ─── Stage 1 : Build ──────────────────────────────────────────────────────────
FROM rust:latest AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build/bsgs2d

# Cache des dépendances Cargo
COPY bsgs2d/Cargo.toml bsgs2d/Cargo.lock* ./
RUN mkdir src && echo 'fn main(){}' > src/main.rs \
    && cargo build --release 2>/dev/null || true \
    && rm -rf src

# Compilation complète avec optimisations maximales
COPY bsgs2d/src ./src
COPY bsgs2d/build.rs ./build.rs 2>/dev/null || true
RUN cargo build --release

# ─── Stage 2 : Runtime ────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /data

COPY --from=builder /build/bsgs2d/target/release/bsgs2d /usr/local/bin/bsgs2d
COPY p135/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# ─── Variables d'environnement ────────────────────────────────────────────────
# Obligatoires (remplacer par les vraies coordonnées du puzzle #135)
ENV TARGET_X=PUZZLE135_X
ENV TARGET_Y=PUZZLE135_Y

# Paramètres du solveur
ENV RANGE_BITS=135
ENV BLOCK_BITS=20
ENV HALF_BITS=0
ENV M_LEVEL=3
ENV THREADS=0

# Mode : solve (full stack) | optimized | parallel | dispatch
ENV SOLVER_MODE=solve

VOLUME ["/data"]
ENTRYPOINT ["/entrypoint.sh"]
