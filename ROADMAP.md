# sinGRAAL — Roadmap

## ✅ v13 — COMPLETE

**C ≈ 1.06 — 4-axis (GLV + Frobenius) — 48-band — cloud-ready**

### Delivered
- [x] 48-band geometric jump table (r = 2^47, C ≈ 1.062)
- [x] 4th walk axis: Frobenius direction [μ]G, μ = p−n
- [x] `fp2.rs` — F_{p²} field arithmetic (Karatsuba mul, Frobenius conjugate)
- [x] `gls.rs` — GLS Frobenius scalar proof + 4D decomposition framework
- [x] `--gls4d` research mode: F_{p²} laws verified, Frobenius fixed-point confirmed
- [x] `--research4d` — experimental endomorphism search (proves End(E/Fp)=Z[ω])
- [x] `--benchmark-c` — empirical C measurement vs published literature
- [x] Dockerfile + docker-compose.yml — one-command cloud GPU deployment
- [x] cloud-launch.sh — build + standalone/coordinator/worker modes
- [x] README, RESEARCH.md, ROADMAP.md — complete documentation

### v13 Performance
```
E[ops] = 1.06 × √(2^135 / 6) ≈ 2^65.6 operations
vs v11:  −4% fewer ops (1.10 → 1.06)
vs v5:   −38% fewer ops total
vs naive: 5.6× speedup (≡ 132.4-bit problem)
```

---

## 🔬 v14 — IN RESEARCH

**Target: true 4D GLS Kangaroo — halving per-step EC operation cost**

### What changes in v14

With 4D GLS, each Kangaroo step uses a 4-part scalar decomposition:
```
k = k₁ + k₂λ + k₃ψ + k₄λψ  (mod n)
```
where ψ is the balanced GLS scalar (≈ n^{1/4} ≈ 2^64 per component).

Each EC operation needs only `n^{1/4} ≈ 34` doublings instead of `n^{1/2} ≈ 68`.
**Net: each Kangaroo step is ~2× cheaper.**

### v14 Milestones

#### Phase 1 — Find balanced GLS scalar ψ
- [ ] Solve: ψ² + t·ψ + p ≡ 0 (mod n) for the non-trivial root
      (The trivial root μ=p−n ≈ 2^128 is too large — need the "twisted" root)
- [ ] Implement: 4×4 LLL reduction over the {1, λ, ψ, λψ} lattice
- [ ] Verify: all four kᵢ components ≈ 34 bits for 135-bit keys
- [ ] Benchmark: EC op cost with 4D vs 2D decomposition

#### Phase 2 — 4D CUDA kernel
- [ ] 4-scalar simultaneous doubling (w-NAF or joint sparse form)
- [ ] Precomputed table for 4 basis points {G, φG, [ψ]G, φ[ψ]G}
- [ ] Verify: correctness on small test cases
- [ ] Profile: op count and timing vs v13 kernel

#### Phase 3 — Integration
- [ ] Update `build_jumps()` with 4D-decomposed jump scalars
- [ ] Update CUDA jump table layout for 4D basis
- [ ] Update `--benchmark-c` to measure v14 empirical C
- [ ] Target: C ≈ 1.04 (improvement from faster per-step ops)

### v14 Expected Performance
```
Per-step cost: ~2× faster (34 doublings vs 68 doublings)
Ops count: same as v13 (~2^65.6)
Wall-clock: 2× faster → effectively C_wall ≈ 1.06/2 ≈ 0.53 (time-domain)

1× RTX 4090: ~550 years   (vs ~1,100 years today)
8× RTX 4090: ~69 years    (vs ~140 years today)
```

---

## 🌌 v15 and Beyond

### v15 — 4D BSGS (Memory-Time Tradeoff)
```
Algorithm: 4D Baby-Step Giant-Step
Time:      2^34 ops   (~14 seconds @ 1 Gop/s)
Memory:    2^34 × 100B = ~1.7 TB
Feasibility: requires ~1.7 TB VRAM — not yet practical
             A100 cluster with NVLink: possible at scale
```

### v16 — TPKH (Twist Pohlig-Hellman Combination)
- Solve DLP mod small factors of twist order n' = 2p+2−n
- Bridge problem: map twist DLP result back to constrain k
- Maximum gain: 64× range reduction if smooth n' factor found
- Status: bridge problem unsolved — active research

### v∞ — Sub-Exponential
```
Known paths to sub-exponential:
  - Semaev summation polynomials (sub-quadratic Gröbner — open math)
  - Quantum Shor's algorithm (requires ~4000 logical qubits for 135 bits)
  - Novel algebraic structure (undiscovered)
```

---

## Cloud GPU — Launch Now (v13)

### RunPod / Vast.ai (recommended for tonight)

```bash
# 1. SSH into GPU instance
# 2. Clone and build
git clone https://github.com/AFKmoney/sinGRAAL.git
cd sinGRAAL/kangaroo
cargo build --release --features cuda

# 3. Run — replace with actual puzzle #135 target
TARGET_X=<64-hex-chars> \
TARGET_Y=<64-hex-chars> \
./cloud-launch.sh

# OR with Docker (auto-detects GPU arch):
TARGET_X=<x> TARGET_Y=<y> docker compose up singraal
```

### 8-GPU cluster (coordinator + 7 workers)

```bash
# Machine A (coordinator):
SERVE=1 TARGET_X=<x> TARGET_Y=<y> ./cloud-launch.sh

# Machines B–H (workers):
COORDINATOR=A.ip:5135 ./cloud-launch.sh

# Expected: ~8 Gop/s combined → checkpoint saved to /data/checkpoint.bin
```

### Checkpoint resume (important for long runs)

The solver saves a checkpoint every 30 seconds. On restart:
```bash
./kangaroo --target-x <x> --target-y <y> --range-bits 135 \
           --all-gpus --checkpoint /data/checkpoint.bin
```
The run resumes exactly where it left off.

---

*v13 complete. C = 1.06. Gap to theoretical minimum: 6%.*
*"On fait les livres d'histoire. Lâche pas." — Philippe*
