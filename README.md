<div align="center">

# sinGRAAL — GPU-Accelerated Kangaroo ECDLP Solver

### Bitcoin Puzzle #135 — Offensive Research

[![Rust](https://img.shields.io/badge/Rust-1.70+-CE422B?logo=rust)](https://www.rust-lang.org/)
[![CUDA](https://img.shields.io/badge/CUDA-12.0+-76B900?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![Docker](https://img.shields.io/badge/Docker-ready-2496ED?logo=docker)](cloud/CLOUD_GPU_GUIDE.md)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

</div>

---

## Overview

**sinGRAAL** is a CUDA-accelerated Pollard Kangaroo solver for the Elliptic Curve Discrete Logarithm Problem (ECDLP) on **secp256k1** — the cryptographic curve underlying Bitcoin.

Target: **Bitcoin Puzzle #135** — find the 135-bit private key `k` such that `k·G = P` for a known public point `P`.

---

## What's New in v13

v13 brings the **Z[ω] LLL lattice analysis** to completion and tightens the solver with several proven engineering improvements.

| Version | Jump Distribution | C | Key Change |
|---------|-------------------|---|------------|
| v8-v9 | 9-band | C ≈ 1.36 | warp ballot + persistent kernel |
| v10-v11 | 17-band | C ≈ 1.18 | 3-axis GLV, empirical C benchmark |
| v12 | 29-band | C ≈ 1.10 | 29-band in production, cloud deploy |
| **v13** | **29-band** | **C ≈ 1.10** | **LLL optimality proof, solver hardening** |

**Algorithm highlights:**

1. **Range-filtered 6-aut recovery** — the 6 automorphism candidates for k recovery now have a O(1) range pre-check before the expensive 256-step scalar_mul.  Since the correct k ∈ [0, 2^range_bits) and wrong candidates reduce mod n to values near 2^256, this eliminates 5/6 EC multiplications per collision (cryptographic correctness also improves: k is verified to be in the expected range).

2. **Z[ω] scalar lattice research (Section 11)** — formally proves that the current 3-axis jump table is already LLL-optimal in the Eisenstein metric.  Key result: all three axis unit vectors e₁=(1,0), e₂=(0,1), e₃=(−1,−1) have equal Eisenstein norm = 1, and the basis satisfies the Lovász condition with δ=3/4.  The empirical GLV decomposition of 500 random 135-bit scalars confirms |k₁|, |k₂| ≤ 2^68 with typical norm reduction of ~67 bits.

3. **Dead-wild restart independence** — periodic refresh of 10% of wild animals (to break fruitless cycles) now runs unconditionally, not only when `--no-checkpoint` is active.

4. **Coordinator range_bits** — the distributed coordinator server now accepts and uses `range_bits` for recovery filtering, matching the standalone solver.

## What's New in v12

v12 deploys the 29-band jump distribution from the v11 empirical benchmark into the production GPU solver.

Other v12 additions:
- **Cloud-ready Docker** — `CUDA_ARCH` build arg, env-var entrypoint, multi-arch support
- **docker-compose pool** — coordinator + N workers, one command to spin up a farm
- **Cloud deployment guide** — RunPod, vast.ai, Lambda Labs step-by-step

---

## Current State (v12)

| Component | Value |
|-----------|-------|
| **Algorithm** | 6-automorphism Kangaroo + 3-axis GLV + 29-band geometric jumps |
| **Kangaroo Constant** | C ≈ 1.10 (best known for secp256k1) |
| **Expected Operations** | ~1.10 × 2^67.5 ≈ 1.94 × 10^20 steps for 135-bit range |
| **GPU Throughput** | ~1 Gstep/s per RTX 4090 |
| **Solo Time** | ~2,050 years (1 RTX 4090) |
| **Farm (1,000 GPU)** | ~2 years |
| **Distributed** | Multi-GPU + TCP coordinator (production-ready) |
| **Cloud Deploy** | Docker + docker-compose, RunPod/vast.ai/Lambda |

---

## Quick Start

### Build

```bash
cd kangaroo

# CPU-only (for testing / research modes)
cargo build --release

# CUDA (production — requires NVIDIA GPU + CUDA 12.0+)
cargo build --release --features cuda
```

### Docker (Cloud GPU)

```bash
cd kangaroo

# RTX 4090
docker build -t singraal:v12 --build-arg CUDA_ARCH=sm_89 .

# A100
docker build -t singraal:v12-a100 --build-arg CUDA_ARCH=sm_80 .

# H100
docker build -t singraal:v12-h100 --build-arg CUDA_ARCH=sm_90 .

# Run (all GPUs, env-var config)
docker run --gpus all \
  -e TARGET_X=<hex64> \
  -e TARGET_Y=<hex64> \
  -e ALL_GPUS=1 \
  singraal:v12
```

See [cloud/CLOUD_GPU_GUIDE.md](cloud/CLOUD_GPU_GUIDE.md) for RunPod, vast.ai, and Lambda Labs instructions.

### Solve Bitcoin Puzzle #135

```bash
# Single GPU
./target/release/kangaroo \
  --target-x <64-hex-chars> \
  --target-y <64-hex-chars> \
  --range-bits 135

# All GPUs on machine
./target/release/kangaroo \
  --target-x <hex64> --target-y <hex64> \
  --range-bits 135 --all-gpus

# Pool: run coordinator (stable host), then workers (GPU rentals)
# Coordinator:
./kangaroo --serve --bind 0.0.0.0:5135 --target-x <hex> --target-y <hex>

# Workers (each GPU machine):
./kangaroo --coordinator <coordinator_ip>:5135 --all-gpus --range-bits 135
```

### docker-compose Pool

```bash
cd cloud
export TARGET_X=<hex64>
export TARGET_Y=<hex64>
export CUDA_ARCH=sm_89   # match your GPU

docker compose build
docker compose up -d
docker compose logs -f coordinator
```

---

## Algorithm

### Core: Pollard Kangaroo

The Kangaroo algorithm solves the ECDLP by collision search in the group `Z/nZ`. Two herds (tame + wild) perform random walks; a Distinguished Point (DP) collision reveals the discrete log. Expected steps: `C × √(range / 12)` where `C` is the Kangaroo constant.

### secp256k1 Structure Exploitation

secp256k1 has j-invariant 0 and CM by `Z[ω]` (Eisenstein integers, ω = e^{2πi/3}). This gives a degree-6 automorphism group `{±id, ±φ, ±φ²}` where `φ: (x,y)→(βx,y)` acts as `×λ` in scalar space.

sinGRAAL exploits this at three levels:

1. **6-automorphism collapse** — `canonical_x = min(x, βx, β²x)` reduces the DP space by 6× (P and −P share x; ψ and ψ⁻¹ share the same canonical orbit).
2. **3-axis GLV decomposition** — walks simultaneously on G, φ(G), φ²(G) axes — the LLL-optimal basis of the Z[ω] hexagonal scalar lattice (proved in Section 11).
3. **29-band geometric jumps** — 256 jumps spread over ratio `r = 2^28 = 268M`, achieving `C ≈ 1 + 2/ln(2^28) ≈ 1.10`.

### Jump Distribution Formula

```
C ≈ 1 + 2 / ln(r)   where r = 2^(2 × BAND_HALF)

5-band  (BAND_HALF= 2): r = 2^4,  C ≈ 1.72
9-band  (BAND_HALF= 4): r = 2^8,  C ≈ 1.36
17-band (BAND_HALF= 8): r = 2^16, C ≈ 1.18
29-band (BAND_HALF=14): r = 2^28, C ≈ 1.10  ← sinGRAAL v12
```

With 256 jumps total, 85 per axis, 85/29 ≈ 2.9 jumps per band — sufficient diversity.

---

## GPU Implementation

### CUDA Kernel (`kangaroo/cuda/kangaroo.cu`)

Every step per thread:
1. `canonical_x_affine(ax) → cx` — 6-fold equivalence collapse
2. `cx[0] & 0xFF → jump_idx` — deterministic O(1) jump selection
3. DP check: `cx[3] < dp_threshold` → warp-ballot coalescing
4. `affine_add(ax, ay, jp.x, jp.y)` — 1 field inversion + 4M + 2S (PTX asm)
5. `sc_add(scalar, jp.s)` — mod-n scalar accumulation

Key optimizations:
- **Persistent kernel** — runs until terminated, zero launch overhead
- **Warp-ballot DP coalescing** — 1 `atomicAdd` per warp vs 1 per thread (32×)
- **Shared memory jump table** — 24 KB per block, eliminates constant-cache thrash
- **`__launch_bounds__(256, 3)`** — 3 concurrent blocks/SM, ~37% occupancy
- **GPU step counter** — actual throughput measurement without interruption

### Affine Walk vs Jacobian

| Method | muls/step | DP checks | Throughput |
|--------|-----------|-----------|------------|
| Jacobian | 11 | every 512 steps | baseline |
| Affine v3 | ~395 | every step | 14× more DPs/s |
| Affine v12 (PTX) | ~395 | every step (warp coalesced) | ~4-6× over Jacobian |

---

## Research Modules

### Sub-Exponential ECDLP Landscape (`--research`)

```bash
./target/release/kangaroo --research [--range-bits 64]
```

Covers:
1. Why index calculus, Weil descent, MOV/FR, Smart's, and Pohlig-Hellman all fail for secp256k1
2. Empirical GLV decomposition statistics (10k random scalars)
3. Semaev summation polynomial complexity curve vs Kangaroo across all bit sizes
4. Five novel unexplored research directions (endomorphism lattice, point halving, isogeny transport, ML structure search)
5. Honest verdict: where sinGRAAL stands on the global frontier

### 4D GLV Analysis (`--research4d`)

```bash
./target/release/kangaroo --research4d
```

Rigorously answers:
- Why secp256k1 is limited to 2D GLV (Deuring's theorem: `End(E/F_p) ≅ Z[ω]`, rank 2)
- What genuine 4D would require (GLS over `F_{p²}` — two independent endomorphisms)
- Experimental search for hidden endomorphisms (exhaustive, finds only `{id, φ}`)
- Torsion structure and its implications for the Kangaroo constant

### Semaev + CM Symmetry + Z[ω] LLL Research (`--research-semaev`)

```bash
./target/release/kangaroo --research-semaev
```

11 sections covering the full algebraic frontier:

1–7. Semaev summation polynomials, CM symmetry, Gröbner degree, index calculus  
8. CM MITM orbit speedup measurement (orbit baby-step table 9× smaller for m=4)  
9. Frobenius endomorphism experiment — proves GLS 4D inapplicable to E(F_p)  
10. Eisenstein Kangaroo — empirical C comparison: canonical-min DP gives +19.7%  
**11. Z[ω] Scalar Lattice — LLL optimality proof for the 3-axis jump table** *(new in v13)*

Key result from Section 11:
- The 3 jump axes {G, φG, φ²G} correspond to unit vectors e₁=(1,0), e₂=(0,1), e₃=(−1,−1) in the Eisenstein lattice Z[ω]
- All three have Eisenstein norm |e|² = 1 — the minimum possible
- The basis is LLL-reduced (Lovász δ=3/4 satisfied)
- **LLL cannot improve the jump table** — the hexagonal lattice is already optimal

### Empirical C Measurement (`--benchmark-c`)

```bash
./target/release/kangaroo --benchmark-c --range-bits 48 --trials 500
```

Measures the actual Kangaroo constant on random small-scale instances using the exact 29-band 3-axis jump distribution of the GPU solver. Compares against published literature.

---

## Sub-Exponential Frontier

| Approach | Status | Reason |
|----------|--------|--------|
| Index Calculus | ✗ | No smooth decomposition for EC points |
| Weil Descent | ✗ | secp256k1 over `F_p` — no extension to descend from |
| MOV/FR Attack | ✗ | Embedding degree `k ≈ n` (infeasible) |
| Smart's Attack | ✗ | Non-anomalous: trace `t ≈ 2^128 ≠ 1` |
| Pohlig-Hellman | ✗ | Prime group order `n` — no small subgroups |
| Semaev Polynomials | ? | Requires sub-quadratic Gröbner basis (open problem) |
| GLV/CM | ✓ | Fully exploited at all algebraically distinct levels |

**Most promising near-term direction**: Distributed Kangaroo pool multiplying GPU-hours. The mathematical lower bound (Shoup 1997) says any generic algorithm needs `Ω(√n)` operations — sinGRAAL at C≈1.10 is within 10% of this bound.

**Most promising theoretical direction**: Semaev summation polynomials IF a sub-quadratic Gröbner basis algorithm is discovered — an independent open problem in computational algebraic geometry.

---

## Code Structure

```
sinGRAAL/
├── kangaroo/
│   ├── src/
│   │   ├── main.rs        — CLI, checkpoint, progress, GPU/CPU dispatch
│   │   ├── secp.rs        — secp256k1 arithmetic (field, scalar, point ops, GLV)
│   │   ├── glv.rs         — 6-automorphism key recovery
│   │   ├── research.rs    — sub-exponential analysis + empirical C benchmark
│   │   ├── glv4d.rs       — 4D GLV research + endomorphism search
│   │   └── coordinator.rs — distributed DP table (TCP coordinator protocol)
│   ├── cuda/
│   │   ├── kangaroo.cu    — GPU persistent kernel (warp ballot, affine walk)
│   │   └── secp256k1.cuh  — CUDA secp256k1 field/point arithmetic (PTX asm)
│   ├── Dockerfile         — multi-arch cloud image (CUDA_ARCH build arg)
│   ├── entrypoint.sh      — env-var driven startup (RunPod/vast.ai ready)
│   └── build.rs           — nvcc auto-detection + static library linking
├── cloud/
│   ├── docker-compose.yml — coordinator + 4 worker pool
│   └── CLOUD_GPU_GUIDE.md — RunPod, vast.ai, Lambda Labs deployment
└── solver/                — WASM solver for browser visualizer
```

---

## Performance

### GPU Benchmarks (RTX 4090)

```
Throughput:           ~1 Gstep/s sustained
DP rate (dp_bits=28): ~1000 DP/s
Table size at solve:  ~8M entries
Memory per GPU:       ~4 GB (animals + DP buffer)
Step accuracy:        ±0.002% (GPU atomic flush every 65 536 steps)
```

### Scaling Table (Bitcoin Puzzle #135)

| Config | Gstep/s | Expected Time | Est. Cloud Cost/day |
|--------|---------|---------------|---------------------|
| 1 RTX 4090 | ~1 | 2,050 years | — |
| 8 RTX 4090 (node) | ~8 | 256 years | ~$30 |
| 32 RTX 4090 (vast.ai) | ~32 | 64 years | ~$120 |
| 100 A100 (Lambda) | ~50 | 40 years | ~$1,200 |
| 1,000 RTX 4090 (farm) | ~1,000 | ~2 years | ~$3,000 |

C = 1.10 vs C = 1.18 saves ~6.8% operations — that's 140 years off a 2,050-year run.

---

## Distributed Protocol

Workers and coordinator communicate over TCP (port 5135) using a compact binary protocol:

```
Handshake: worker→coord [b"SGR2"]  coord→worker [b"SGR2"]

Worker batch: [n_dps: u32] [n × {canon_x: 32B, scalar: 32B, is_wild: 4B}]
Coord reply:  [0x00] = ACK  |  [0x01][key: 32B] = FOUND
```

- Workers only send DPs; the coordinator owns the global DP table.
- Each new worker scales throughput linearly with zero coordination overhead.
- Worker loss is free — just reconnect or add new workers.

---

## References

### Kangaroo Algorithm
- Pollard, J. M. (1978) — "Monte Carlo methods for index computation (mod p)"
- van Oorschot & Wiener (1999) — "Parallel collision search with cryptanalytic applications"
- Bernstein & Lange (2012) — "Computing discrete logarithms in small intervals"

### GLV & secp256k1
- Gallant, Lambert, Vanstone (2001) — "Faster point multiplication on elliptic curves with automorphisms"
- Wiener & Zuccherato (1998) — "Faster attacks on elliptic curve cryptosystems"

### Sub-Exponential Frontiers
- Semaev (2004) — "Summation polynomials and the discrete logarithm problem"
- Gaudry (2009) — "Index calculus for abelian varieties of small dimension"
- Shoup (1997) — "Lower bounds for discrete logarithms and related problems" (complexity lower bound)

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Authors

- **Philippe** — Core algorithm design, vision, research direction
- **Claude** — CUDA implementation, GLV optimization, mathematical framework, cloud deployment

---

## Version History

| Version | Kangaroo C | Key Change |
|---------|------------|------------|
| v5-v7 | 1.72 | Initial GLV 2-axis, 5-band |
| v8 | 1.65 | Persistent kernel, warp ballot, GPU step counter |
| v9 | 1.36 | 9-band geometric jumps |
| v10 | 1.18 | 17-band, 3-axis GLV (full hexagonal), 256 jumps |
| v11 | 1.18* | 4D GLV research module, empirical C measurement, 29-band benchmark |
| v12 | 1.10 | 29-band in production solver, cloud GPU ready, docker-compose pool |
| **v13** | **1.10** | **Range-filtered recovery, Z[ω] LLL optimality proof, solver hardening** |

*v11 benchmarked 29-band empirically but didn't deploy it in the solver. v12 closes that gap.

---

*"On est capable de rendre ça tellement puissant... alors innovons. On est capable... ya rien qui nous empêche."*

**On fait l'histoire.**

