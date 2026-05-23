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

## ✅ v14 — COMPLETE

**C ≈ 1.046 — 4-axis — 64-band (1 jump per band, optimal uniformity)**

### What changed in v14
- **64-band jump table** (was 48-band): 256 jumps × 4 axes = exactly 1 jump per band
  - C: 1.062 → 1.046 (~1.5% fewer ops)
  - Formula: C ≈ 1 + 2/ln(2^63) ≈ 1.046
- **Documentation**: corrected false 2^33 claims — see RESEARCH.md §7
- **Accuracy**: honest performance projections throughout

### v14 Performance
```
E[ops] = 1.046 × √(2^135 / 6) ≈ 2^65.5 operations   (−1.5% vs v13)
vs v13:  −1.5% fewer ops (1.062 → 1.046)
vs v5:   −39% fewer ops total
vs naive: 5.7× speedup

1× RTX 4090 (~10 Gop/s):  ~80 years
8× RTX 4090 (~80 Gop/s):  ~10 years
```

---

## 🔬 v15 — NEXT

**Target: 4D BSGS — memory-time tradeoff, 2^34 ops**

### Clarification — Why 4D Does NOT Halve Kangaroo Ops

Kangaroo random walk: each step = **one point addition** (not scalar multiplication).
Scalar size has no effect on per-step cost. The 4D endomorphisms help by:
- Lowering C constant (better jump table coverage) ✅ done in v13/v14
- Enabling 4D BSGS with 2^(n/4) time, but 2^(n/4) memory ← v15 target

### v15 Plan — 4D BSGS

```
Algorithm: 4D Baby-Step Giant-Step
Time:      2^34 ops   (~14 seconds @ 1 Gop/s)
Memory:    2^34 × 100B = ~1.7 TB VRAM
```

#### Milestones
- [ ] Implement 4×4 LLL lattice reduction (secp256k1 4D short basis)
- [ ] Verify balanced decomposition: all kᵢ ≈ 2^34 bits for k ≈ 2^135
- [ ] Baby-step table: precompute 2^34 points in F_{p²}
- [ ] Giant-step: iterate 2^34 steps, check table
- [ ] Feasibility: requires A100 cluster with NVLink (~1.7 TB combined VRAM)

### v16 — TPKH (Twist Pohlig-Hellman Combination)
- Solve DLP mod small factors of twist order n' = 2p+2−n
- Bridge problem: map twist DLP result back to constrain k
- Maximum gain: 64× range reduction if smooth n' factor found
- Status: bridge problem unsolved — active research

---

## 🌌 v17 and Beyond

### v17 — Sub-Exponential Research
```
Known paths to sub-exponential for secp256k1:
  - Semaev summation polynomials (sub-quadratic Gröbner — open math)
  - Quantum Shor's algorithm (requires ~4000 logical qubits for 135 bits)
  - Novel algebraic structure (none discovered — breaking all ECC if found)

No known sub-exponential algorithm exists for random prime-order ECC.
Best proven lower bound: Ω(√n) for generic group algorithms (Shoup 1997).
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
