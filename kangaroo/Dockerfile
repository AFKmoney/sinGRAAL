# sinGRAAL Kangaroo — Cloud GPU Docker image
# Base: CUDA 12.4 + Ubuntu 22.04
# Build: docker build -t singraal-kangaroo .
# Run:   docker run --gpus all singraal-kangaroo \
#          --target-x <hex64> --target-y <hex64> --range-bits 135

FROM nvidia/cuda:12.4.1-devel-ubuntu22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# System deps
RUN apt-get update && apt-get install -y \
    curl build-essential pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Rust toolchain
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /build
COPY . .

# Build with CUDA support (sm_80 = A100; override with --build-arg CUDA_ARCH=sm_86 for RTX 30xx)
ARG CUDA_ARCH=sm_80
ENV CUDA_ARCH=${CUDA_ARCH}

RUN cargo build --release 2>&1

# ─── Runtime image ────────────────────────────────────────────────────────────
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

RUN apt-get update && apt-get install -y libstdc++6 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/kangaroo /usr/local/bin/kangaroo

ENTRYPOINT ["/usr/local/bin/kangaroo"]
CMD ["--help"]
