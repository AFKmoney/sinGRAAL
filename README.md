<div align="center">

# sinGRAAL v13 — GPU-Accelerated ECDLP Solver

### Bitcoin Puzzle #135 · Kangaroo + Coppersmith · Cloud-Ready

[![Rust](https://img.shields.io/badge/Rust-1.70+-CE422B?logo=rust)](https://www.rust-lang.org/)
[![CUDA](https://img.shields.io/badge/CUDA-12.2+-76B900?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![Docker](https://img.shields.io/badge/Docker-ready-2496ED?logo=docker)](cloud/CLOUD_GPU_GUIDE.md)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

</div>

---

## What is sinGRAAL?

sinGRAAL is a high-performance solver for the **Elliptic Curve Discrete Logarithm Problem (ECDLP)** on secp256k1, targeting **Bitcoin Puzzle #135** (135-bit private key).

It combines three independent attack layers:

| Layer | Algorithm | Reduces |
|---|---|---|
| **GPU walk engine** | Pollard Kangaroo + 6-automorphisms + GLV 4D | C constant → 0.55 |
| **Algebraic filter** | Coppersmith bivariate (Semaev S₃, m=2/3) | Pruning of (x₁,x₂) blocks |
| **Distributed pool** | TCP coordinator + N workers | Linear GPU scaling |

---

## Target

```
Puzzle #135 public key (compressed):
02145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16

TARGET_X = 145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16
TARGET_Y = (decompressed y-coordinate)

Range: k ∈ [2^134, 2^135)   →   135-bit private key
Reward: 135 BTC
```

---

## Architecture

### 1. GPU Walk Engine (`singraal/cuda/kangaroo_toom6.cu`)

**Pollard Kangaroo** with 5 stacked optimizations:

#### 6-Automorphisms (÷√6 ≈ 2.45×)
secp256k1 has an efficiently computable endomorphism ψ(x,y) = (β·x, y) where β is a cube root of 1 mod p. Points P, ψ(P), ψ²(P) and their negatives share the same **canonical_x**:
```
canonical_x(x) = min(x, β·x mod p, β²·x mod p)
```
One distinguished point detection covers 6 equivalent positions → 6× fewer DPs needed.

#### Bidirectional Walk (÷√2 ≈ 1.41×)
- **Tame animals** walk forward (+)
- **Wild animals** walk backward (−)
- Directed convergence reduces expected steps by √2

#### GLV 4D Halton Jump Table (×0.95)
The GLV decomposition k = k₁ + λk₂ defines a 2D lattice. The jump table uses 4 orthogonal directions:
```
dk₁: direction G
dk₂: direction φG  (Frobenius endomorphism)
dk₃: direction G+φG
dk₄: direction G−φG
```
Jump sizes are drawn from a **Halton low-discrepancy sequence** (bases 2,3,5,7) across 29 geometric bands → minimal variance, no clustering.

#### 5-Level Decorrelation
| Level | Mechanism |
|---|---|
| L1 | Full-state hash (canonical_x + scalar + k₁..k₄) |
| L2 | Per-animal stable scramble (Knuth hash on thread ID) |
| L3 | Anti-cycle ring buffer (4 entries, detects short cycles) |
| L4 | Stagnation escape (32M steps without DP → reset) |
| L5 | Periodic LCG scramble evolution (every 2²⁰ steps) |

#### Toom-Cook-6 Field Arithmetic (−31% multiplications)
Each field element is split into 6 × 43-bit limbs, evaluated at {0,1,...,9,∞}, multiplied, then interpolated via Newton differences:
- Schoolbook 256×256: **16 MAD** operations
- Toom-Cook-6: **11 MAD** operations (−31%)
- Impact on fp_inv (256S + 15M): 240 MAD → 165 MAD per inversion

#### Performance Formula

```
E[steps] = C × √(range / 6)    with C ≈ 0.55

Puzzle #135 (range = 2^135):
  E[steps] = 0.55 × 2^67.5 / √6 ≈ 2^65.3 total steps
```

| Config | Gstep/s | ETA |
|---|---|---|
| 1× RTX 4090 | ~1.5 | ~600 years |
| 8× RTX 4090 | ~12 | ~75 years |
| 100× RTX 4090 | ~150 | ~6 years |
| 1 000× RTX 4090 | ~1 500 | ~220 days |
| 10 000× RTX 4090 | ~15 000 | ~22 days |

---

### 2. Coppersmith Algebraic Filter (`bsgs2d/src/coppersmith.rs`)

Uses **Semaev's summation polynomial S₃** to prune (x₁, x₂) coordinate pairs before running the full kangaroo walk on a block.

**S₃ polynomial for secp256k1** (b=7):
```
S₃(x₁,x₂,x₃) = (x₁−x₂)²·x₃² − 2[x₁x₂(x₁+x₂)+14]·x₃ + (x₁x₂)²−28(x₁+x₂)
```

Fix x₃ = TARGET_X, get bivariate f(x₁,x₂) = S₃(x₁, x₂, target).

**Jochemsz-May Macaulay matrix** with Howgrave-Graham bound:

| Parameter m | Matrix size | HG bound | block_bits |
|---|---|---|---|
| m = 2 | 15 × 15 | p⁴/15 ≈ 2^1021 | ~5 bits |
| m = 3 | 28 × 28 | p⁶/28 ≈ 2^1526 | ~8 bits |

**LLL reduction**: incremental BigRational Gram-Schmidt — O(n) per swap, guaranteed convergence (no Cohen integral oscillation).

If LLL finds a short vector with norm² < HG bound, f has a small root (x₁,x₂) → block is viable and fed to the GPU engine.

---

### 3. Distributed Coordinator (`singraal/src/coordinator.rs`)

Binary TCP protocol (SGR2):
```
Handshake magic: "SGR2"
Worker → Coord: [n_dps: u32][n × 68 bytes (canonical_x + scalar + is_wild)]
Coord → Worker: [0x00] = OK | [0x01][32 bytes] = FOUND
```

Scales linearly with GPU count. Coordinator saves checkpoint every 60s.

---

## Quick Start

### Prerequisites
- CUDA 12.2+ with `nvcc`
- Rust stable (1.70+)
- GPU: RTX 3090 / 4090 / A100 / H100

### Build

```bash
git clone https://github.com/afkmoney/singraal
cd singraal/singraal

# With CUDA (production)
CUDA_ARCH=sm_89 cargo build --release --features cuda

# CPU-only (testing)
cargo build --release
```

### Run Standalone (single machine)

```bash
./target/release/singraal \
  --target-x 145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16 \
  --target-y <64-hex-chars> \
  --range-bits 135 \
  --all-gpus \
  --num-animals 262144
```

### Run Distributed (coordinator + workers)

**Coordinator** (one machine):
```bash
./target/release/singraal \
  --serve \
  --target-x <hex> --target-y <hex> \
  --range-bits 135 \
  --bind 0.0.0.0:5135 \
  --checkpoint /data/coordinator.ckpt
```

**Workers** (N machines):
```bash
./target/release/singraal \
  --coordinator <coordinator_ip>:5135 \
  --all-gpus \
  --num-animals 262144
```

---

## Cloud Deployment

### Docker (single GPU)

```bash
cd singraal
docker build -t singraal:v13 --build-arg CUDA_ARCH=sm_89 .

docker run --gpus all \
  -e TARGET_X=145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16 \
  -e TARGET_Y=<hex> \
  -e RANGE_BITS=135 \
  -v /data/singraal:/data \
  singraal:v13
```

### Docker Compose (4 GPUs, same machine)

```bash
cd cloud
export TARGET_X=145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16
export TARGET_Y=<hex>
export CUDA_ARCH=sm_89

docker compose build
docker compose up -d
docker compose logs -f coordinator
```

### CUDA_ARCH by GPU

| GPU | CUDA_ARCH |
|---|---|
| H100 | `sm_90` |
| A100 | `sm_80` |
| RTX 4090 | `sm_89` |
| RTX 3090 / 3080 Ti | `sm_86` |
| RTX 2080 Ti | `sm_75` |
| V100 | `sm_70` |

See [cloud/CLOUD_GPU_GUIDE.md](cloud/CLOUD_GPU_GUIDE.md) for RunPod, vast.ai, and Lambda Labs step-by-step setup.

---

## Repository Structure

```
singraal/                ← Main solver (use this)
  src/
    main.rs              — Orchestration, DP table, CLI
    secp.rs              — secp256k1 CPU arithmetic (scalar_mul, canonical_x)
    glv.rs               — 6-automorphism recovery (recover_k_6aut)
    coordinator.rs       — Distributed TCP protocol (SGR2)
  cuda/
    kangaroo_toom6.cu    — Main CUDA kernel (walk, DP, Jacobian init)
    secp256k1_toom6.cuh  — secp256k1 GPU arithmetic + Toom-Cook-6
  config.toml            — Puzzle #135 parameters
  Dockerfile             — Production CUDA build
  entrypoint.sh          — Docker entry (SERVE / COORDINATOR / standalone)

bsgs2d/                  ← Coppersmith algebraic filter
  src/
    coppersmith.rs       — Semaev S₃, Macaulay m=2/3, incremental LLL

cloud/
  docker-compose.yml     — Coordinator + 4 workers (multi-GPU)
  CLOUD_GPU_GUIDE.md     — RunPod / vast.ai / Lambda deployment

solver/                  ← WASM web interface (separate, experimental)
src/                     ← TypeScript frontend visualizer
```

---

## Monitoring

Coordinator output (every 10s):
```
[coord] 12.3M DPs total | 45.2k DP/s | 32 workers | table=8192
```

Worker output:
```
[GPU 0] 4.20B steps (18.3%) | 12450 DPs | 1.03 Gstep/s | 41.2 DP/s | table=8192 | ETA~1850.2d
```

---

## License

MIT — see [LICENSE](LICENSE).
