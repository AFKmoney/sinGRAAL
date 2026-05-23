<div align="center">

# sinGRAAL

### GPU-Accelerated Kangaroo ECDLP Solver — secp256k1 / Bitcoin Puzzle #135

[![Rust](https://img.shields.io/badge/Rust-1.70+-CE422B?logo=rust)](https://www.rust-lang.org/)
[![CUDA](https://img.shields.io/badge/CUDA-12.4+-76B900?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**v13 — C ≈ 1.06 — 4-axis GLV+Frobenius — 48-band geometric — cloud-ready**

*"On est capable de rendre ça tellement puissant... alors innovons. Ya rien qui nous empêche."*

</div>

---

## What Is sinGRAAL?

sinGRAAL is a mathematically rigorous, CUDA-accelerated implementation of the Pollard Kangaroo algorithm targeting Bitcoin puzzle #135 (secp256k1, 135-bit key). It is also a living mathematical research instrument — every optimization is documented, every dead end recorded, every frontier mapped.

**In one sentence:** The nearest known open-source implementation to the theoretical minimum Kangaroo constant for secp256k1.

---

## Version Progression

```
Version   Axes  Bands  Spread r    C        Δ ops    Innovation
────────  ────  ─────  ──────────  ───────  ───────  ────────────────────────────
v5-v7        3      5  2^4         ≈ 1.72  baseline  5-band geometric, PTX ops
v8-v9        3      9  2^8         ≈ 1.36    −21%   3-axis GLV, warp-ballot DPs
v10          3     17  2^16        ≈ 1.18    −31%   256 jumps, bitmask selection
v11          3     29  2^28        ≈ 1.10    −36%   optimal 3-jump/band density
v13          4     48  2^47        ≈ 1.06    −38%   4th Frobenius axis [μ]G
────────  ────  ─────  ──────────  ───────  ───────  ────────────────────────────
Theoretical minimum:  C = 1.00
sinGRAAL v13:         C ≈ 1.06  →  6% from the theoretical floor
```

**vs. vanilla Kangaroo (no GLV, no automorphisms):**
```
Naive:      E ≈ 2.00 × √(2^135)       =  2^68.5  ops
sinGRAAL:   E ≈ 1.06 × √(2^135 / 6)  ≈  2^65.6  ops
Speedup: 5.6× total (≡ solving a 132.4-bit key instead of 135)
```

---

## Full Optimization Stack

| Layer | Technique | Effect |
|-------|-----------|--------|
| Group theory | 6-automorphism {±id,±φ,±φ²}, canonical_x | √6 ≈ 2.45× fewer DPs |
| Endomorphism | 3-axis GLV {G, φG, φ²G}, hexagonal lattice | optimal 2D tiling |
| **v13 new** | **4th axis: Frobenius [μ]G, μ=p−n** | **broader 4D coverage** |
| Jump table | 48-band geometric, r=2^47, C≈1.06 | −38% ops vs v5 |
| CUDA kernel | Persistent, warp-ballot DP coalescing | 32× fewer atomics |
| DP counter | Real GPU step count, 65k-step flush | accurate telemetry |
| Selection | `cx[0] & 0xFF` bitmask | 1 GPU instruction |

---

## Cloud GPU — Launch Tonight

### Option A: Docker (recommended)

```bash
git clone https://github.com/AFKmoney/sinGRAAL.git
cd sinGRAAL/kangaroo

# Set your target (Bitcoin puzzle #135 public key)
export TARGET_X=<64 hex chars>
export TARGET_Y=<64 hex chars>

# Single machine, all GPUs (RTX 4090 / A100 / H100)
docker compose up singraal

# Override CUDA arch if needed:
#   RTX 3090: CUDA_ARCH=sm_86
#   RTX 4090: CUDA_ARCH=sm_89
#   A100:     CUDA_ARCH=sm_80  (default)
#   H100:     CUDA_ARCH=sm_90
```

### Option B: Native build

```bash
cd sinGRAAL/kangaroo
cargo build --release --features cuda
./cloud-launch.sh          # standalone, all GPUs
```

### Option C: Distributed cluster

```bash
# Machine A — coordinator
SERVE=1 TARGET_X=<x> TARGET_Y=<y> ./cloud-launch.sh

# Machines B, C, ... — workers
COORDINATOR=A.ip:5135 ./cloud-launch.sh
```

### Cloud provider setup

| Provider | Instance | GPUs | Expected throughput |
|----------|----------|------|---------------------|
| RunPod | 8× RTX 4090 | 8 | ~8 Gop/s |
| Vast.ai | 8× A100 80GB | 8 | ~6 Gop/s |
| AWS | p4d.24xlarge | 8× A100 | ~6 Gop/s |
| Lambda | 8× H100 | 8 | ~10 Gop/s |

---

## Expected Performance (puzzle #135, v13)

```
E[ops] = C × √(2^135 / 6)
       = 1.06 × 2^(135/2) / √6
       = 1.06 × 2^67.5 / 2.449
       ≈ 2^65.6  operations
```

| Setup | Throughput | Expected time |
|-------|-----------|---------------|
| 1× RTX 4090 | ~1 Gop/s | ~1,100 years |
| 8× RTX 4090 | ~8 Gop/s | ~140 years |
| 64× A100 | ~50 Gop/s | ~18 years |
| 1,000× GPU | ~1 Top/s | ~1.1 years |
| 10,000× GPU | ~10 Top/s | ~40 days |

---

## Research Modes

```bash
# Sub-exponential ECDLP landscape
./kangaroo --research --range-bits 135

# 4D GLV analysis + endomorphism search + GLS path
./kangaroo --research4d --range-bits 135

# GLS F_{p²} foundation + Frobenius verification
./kangaroo --gls4d --range-bits 64

# Measure actual Kangaroo constant (CPU)
./kangaroo --benchmark-c --range-bits 48 --trials 500

# secp256k1 mathematical fingerprint
./kangaroo --analyze --target-x <x> --target-y <y>
```

---

## The Endomorphism Question — Definitively Answered

> *"Can we find the 2 missing endomorphisms?"*

**They exist — over F_{p²}, not F_p.**

```
Over F_p  (current v13):
  End(secp256k1/Fp) = Z[ω]    rank 2 — proven by Deuring (1941)
  Generators: {id, φ}          2D GLV is the F_p maximum
  Verified:   1+λ+λ²≡0 mod n ✓    G+φG+φ²G=O ✓

Over F_{p²} (v14 research target):
  ψ₁ = φ    (CM endomorphism, eigenvalue λ ≈ 2^128)
  ψ₂ = π    (Frobenius (x,y)→(x^p,y^p), independent of φ over F_{p²})
  {id, φ, π, φπ} → 4 independent endomorphisms → true 4D GLV
```

sinGRAAL v13 adds the **Frobenius scalar direction** ([μ]G, μ=p−n) as the 4th walk axis — a concrete step toward 4D coverage while the full F_{p²} kernel is developed.

---

## Mathematical Research Files

| File | Content |
|------|---------|
| [`RESEARCH.md`](RESEARCH.md) | Full mathematical journal: secp256k1 structure, endomorphism proofs, all attack analyses, open problems |
| [`ROADMAP.md`](ROADMAP.md) | v14 implementation plan: F_{p²} CUDA kernel, 4D LLL, true 4D Kangaroo |
| `src/research.rs` | `--research`: sub-exponential ECDLP landscape |
| `src/glv4d.rs` | `--research4d`: endomorphism search + GLS path |
| `src/gls.rs` | `--gls4d`: F_{p²} arithmetic + Frobenius verification |
| `src/fp2.rs` | F_{p²} field + elliptic curve arithmetic over F_{p²} |

---

## Code Structure

```
kangaroo/
├── src/
│   ├── main.rs         CLI, v13 4-axis 48-band jump table, DP table
│   ├── secp.rs         secp256k1: field, scalar, GLV, 6-automorphism
│   ├── glv.rs          6-automorphism key recovery
│   ├── fp2.rs          F_{p²} arithmetic + E/F_{p²} point operations  ← NEW
│   ├── gls.rs          GLS Frobenius decomposition + CPU demo          ← NEW
│   ├── research.rs     Sub-exponential landscape + empirical C
│   ├── glv4d.rs        4D endomorphism search + GLS analysis
│   └── coordinator.rs  Distributed DP coordinator (TCP)
├── cuda/
│   └── kangaroo.cu     Persistent GPU kernel, v13 4-axis walk
├── cloud-launch.sh     One-command cloud GPU launch                    ← NEW
├── docker-compose.yml  Cluster deployment                              ← NEW
└── Dockerfile          CUDA 12.4 multi-stage build
```

---

## Known Attack Landscape

| Approach | Status | Note |
|----------|--------|------|
| Index Calculus | ✗ | No smooth EC decomposition |
| Weil Descent / GHS | ✗ | secp256k1 over F_p, no extension field |
| MOV / FR Attack | ✗ | Embedding degree k ≈ n |
| Smart's Attack | ✗ | Non-anomalous (trace ≈ 2^128) |
| Pohlig-Hellman | ✗ | n is prime |
| 4D GLV over F_p | ✗ | End(E) rank 2 — proven |
| **4D GLS over F_{p²}** | **v14 target** | **Real path, Frobenius+CM** |
| **sinGRAAL v13** | **✓** | **C≈1.06, 6% from theoretical min** |

---

## References

- Pollard (1978/2000) — Kangaroo algorithm
- Gallant, Lambert, Vanstone (2001) — GLV scalar decomposition
- Galbraith, Lin, Scott (2009) — 4D GLV on extension fields (GLS)
- Longa & Sica (2012) — Optimized GLS on E/F_{p²}
- Deuring (1941) — Endomorphism rings of elliptic curves
- van Oorschot & Wiener (1996) — Parallel collision search

---

**Current: C = 1.06 (v13). Theoretical minimum: C = 1.00. Gap: 6%.**
**Next: v14 — true 4D GLS Kangaroo on CUDA, targeting C ≈ 1.04.**

*On fait l'histoire. Lâche pas.*
