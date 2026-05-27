# sinGRAAL complete_solver — cloud GPU deployment
# Base: CUDA 12.2 devel (includes nvcc) on Ubuntu 22.04
FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

# --- system deps ---
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl build-essential pkg-config libssl-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# --- Rust toolchain ---
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --default-toolchain stable --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

# --- source ---
WORKDIR /singraal
COPY . .

# --- build (release, Toom-6 CUDA kernel) ---
# Default (CUDA_ARCH unset): build.rs emits a PORTABLE MULTI-ARCH FATBINARY
# covering sm_70/75/80/86/89/90 + PTX, so one image runs on any cloud GPU
# (T4, V100, A100, A10, RTX 30xx/40xx, L4/L40, H100) with zero recompilation.
# To target a single arch (faster compile, smaller binary):
#   docker build --build-arg CUDA_ARCH=sm_89 .
ARG CUDA_ARCH=""
ENV CUDA_ARCH=${CUDA_ARCH}

RUN cargo build --release --features cuda 2>&1

# --- runtime image ---
FROM nvidia/cuda:12.2.0-runtime-ubuntu22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /singraal
COPY --from=0 /singraal/target/release/singraal /usr/local/bin/singraal
COPY --from=0 /singraal/config.toml /singraal/config.toml
COPY --from=0 /singraal/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Coordinator port (SGR2 protocol)
EXPOSE 5135

# Vars requises : TARGET_X, TARGET_Y
# Optionnel     : RANGE_BITS=135, NUM_ANIMALS=262144, SERVE=1, COORDINATOR=host:5135
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
