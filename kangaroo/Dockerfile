# sinGRAAL Kangaroo v12 — Cloud GPU Docker image
#
# CUDA_ARCH target selection (pass with --build-arg CUDA_ARCH=<arch>):
#   sm_70  — Tesla V100    (Lambda Labs, old AWS p3)
#   sm_75  — RTX 2080 Ti
#   sm_80  — A100, RTX 30xx (default — best for cloud: vast.ai, RunPod A100)
#   sm_86  — RTX 3090 / 3080 Ti
#   sm_89  — RTX 4090 / 4080 (Ada Lovelace — best for RunPod H4, vast.ai RTX 4090)
#   sm_90  — H100 (Lambda Labs, CoreWeave, RunPod H100)
#
# Build examples:
#   docker build -t singraal:v12-a100   --build-arg CUDA_ARCH=sm_80 .
#   docker build -t singraal:v12-4090   --build-arg CUDA_ARCH=sm_89 .
#   docker build -t singraal:v12-h100   --build-arg CUDA_ARCH=sm_90 .
#
# Run standalone:
#   docker run --gpus all singraal:v12-4090 \
#     --target-x <hex64> --target-y <hex64> --range-bits 135 --all-gpus
#
# Run as pool worker (connect to coordinator):
#   docker run --gpus all \
#     -e COORDINATOR=10.0.0.1:5135 \
#     -e TARGET_X=<hex64> \
#     -e TARGET_Y=<hex64> \
#     singraal:v12-4090
#
# Run as coordinator:
#   docker run --gpus all -p 5135:5135 \
#     -e SERVE=1 \
#     -e TARGET_X=<hex64> \
#     -e TARGET_Y=<hex64> \
#     singraal:v12-4090

FROM nvidia/cuda:12.4.1-devel-ubuntu22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    curl build-essential pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Rust toolchain
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --default-toolchain stable --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /build
COPY . .

ARG CUDA_ARCH=sm_80
ENV CUDA_ARCH=${CUDA_ARCH}

# Build with CUDA (build.rs auto-detects nvcc and enables the cuda feature)
RUN cargo build --release --features cuda 2>&1 | tail -5 \
    && echo "Build complete: $(ls -lh target/release/kangaroo)"

# ─── Runtime image ────────────────────────────────────────────────────────────
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

RUN apt-get update && apt-get install -y libstdc++6 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/kangaroo /usr/local/bin/kangaroo
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Default checkpoint directory (mount a volume here to persist across restarts)
VOLUME ["/data"]
WORKDIR /data

ENTRYPOINT ["/entrypoint.sh"]
