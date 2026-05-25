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
# Override CUDA_ARCH at build time: docker build --build-arg CUDA_ARCH=sm_86 .
ARG CUDA_ARCH=sm_80
ENV CUDA_ARCH=${CUDA_ARCH}

RUN cargo build --release 2>&1

# --- runtime image ---
FROM nvidia/cuda:12.2.0-runtime-ubuntu22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /singraal
COPY --from=0 /singraal/target/release/singraal /usr/local/bin/singraal
COPY --from=0 /singraal/config.toml /singraal/config.toml

# Coordinator port (SGR2 protocol)
EXPOSE 5135

# Default: standalone solver on puzzle #135
# Override CMD for coordinator/worker modes
CMD ["singraal", "--config", "/singraal/config.toml"]
