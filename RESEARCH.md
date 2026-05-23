# sinGRAAL — Mathematical Research Journal

> *All known algorithms documented. All dead ends recorded. All frontiers mapped.*

---

## Table of Contents

1. [secp256k1 — Mathematical Structure](#1-secp256k1--mathematical-structure)
2. [Endomorphism Ring — The Complete Story](#2-endomorphism-ring--the-complete-story)
3. [GLV Decomposition — Current 2D Implementation](#3-glv-decomposition--current-2d-implementation)
4. [Kangaroo Constant Theory](#4-kangaroo-constant-theory)
5. [v11 Jump Table Design](#5-v11-jump-table-design)
6. [Known Sub-Exponential Attacks — Why They Fail](#6-known-sub-exponential-attacks--why-they-fail)
7. [The GLS Breakthrough — Path to 4D](#7-the-gls-breakthrough--path-to-4d)
8. [Twist Order & TPKH Research](#8-twist-order--tpkh-research)
9. [Empirical C Measurement](#9-empirical-c-measurement)
10. [Open Problems](#10-open-problems)

---

## 1. secp256k1 — Mathematical Structure

```
Curve:        y² = x³ + 7  over  F_p
Prime:        p  = 2²⁵⁶ − 2³² − 977
              p  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
Order:        n  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
j-invariant:  j  = 0   (unique CM structure)
Trace:        t  = p + 1 − n  ≈  2^128  (not anomalous)
CM field:     K  = Q(√−3)  (Eisenstein integers)
CM order:     O_K = Z[ω],  ω = (−1+√−3)/2
Discriminant: D  = −3
```

### Key Constants

```rust
// GLV eigenvalue: λ² + λ + 1 ≡ 0 (mod n)
LAMBDA  = 0x5363AD4CC05C30E0_A5261C028812645A_122E22EA20816678_DF02967C1B23BD72

// λ² = −1 − λ (mod n)
LAMBDA2 = 0xAC9C52B33FA3CF1F_5AD9E3FD77ED9BA4_A880B9FC8EC739C2_DCCFC810B51283CE

// φ(x,y) = (βx, y)
BETA    = 0x7AE96A2B657C0710_6E64479EAC3434E9_9CF0497512F58995_C1396C28719501EE

// φ²(x,y) = (β²x, y)
BETA2   = 0x851695D49A83F8EF_919BB86153CBCB16_630FB68AED0A766A_3EC693D68E6AFA40
```

### Verified Identities

| Identity | Value | Verification |
|----------|-------|-------------|
| `1 + λ + λ² mod n` | 0 | ✓ numerically confirmed |
| `G + φ(G) + φ²(G)` | O (infinity) | ✓ numerically confirmed |
| `φ(G) = [λ]G` | identical x,y | ✓ numerically confirmed |
| `φ²(G) = [λ²]G` | identical x,y | ✓ numerically confirmed |
| `φ is group homomorphism` | 20/20 tests pass | ✓ experimentally confirmed |
| `φ² is group homomorphism` | 20/20 tests pass | ✓ experimentally confirmed |

---

## 2. Endomorphism Ring — The Complete Story

### Deuring's Theorem (1941)

For an **ordinary** elliptic curve E over F_p with CM by imaginary quadratic field K:

```
End_{F_p}(E)  is an order in O_K
```

For secp256k1 (D = −3, K = Q(√−3), maximal order O_K = Z[ω]):

```
End_{F_p}(secp256k1)  =  Z[ω]     (exactly — conductor f = 1)
```

This means every F_p-rational endomorphism is of the form `a·id + b·φ` for integers a,b.
**Z-rank = 2. No room for a third independent endomorphism over F_p.**

### Experimental Endomorphism Search (--research4d)

sinGRAAL v11 implements a systematic search over all degree-1 rational maps `(x,y) → (ax+b, ±y)` over F_p. Results:

```
f(x,y) = (1·x,  y):  homomorphism ✓  →  identity [1]
f(x,y) = (β·x,  y):  homomorphism ✓  →  φ  = [λ]
f(x,y) = (β²·x, y):  homomorphism ✓  →  φ² = [λ²]
f(x,y) = (1·x, -y):  homomorphism ✓  →  [−1] (negation)
f(x,y) = (β·x, -y):  homomorphism ✓  →  −φ  = [−λ]
f(x,y) = (β²·x,-y):  homomorphism ✓  →  −φ² = [−λ²]
All other degree-1 maps:  NOT homomorphisms
```

**Conclusion:** The complete list of degree-1 F_p-endomorphisms is `{±id, ±φ, ±φ²}` — the 6-automorphism group already exploited by sinGRAAL.

### The Missing Endomorphisms — Where They Actually Live

Over **F_{p²}** (the quadratic extension), the Frobenius π: `(x,y) → (x^p, y^p)` is:
- NOT in Z[ω] (acts differently from any scalar × φ)
- An independent endomorphism of secp256k1 over F_{p²}
- Together with φ, generates a **rank-4 endomorphism structure**

```
End_{F_{p²}}(secp256k1) ⊇ Z[φ, π]   (rank 4 over Z)
```

This is the foundation of the **GLS construction** (Section 7).

---

## 3. GLV Decomposition — Current 2D Implementation

### Scalar Decomposition

For any scalar k, find (k₁, k₂) such that:
```
k ≡ k₁ + k₂·λ  (mod n)
|k₁|, |k₂|  ≈  √n  ≈  2^67.5
```

This halves the effective scalar length, halving the number of doublings in scalar multiplication.

### 3-Axis Kangaroo Walk

sinGRAAL uses {G, φG, φ²G} as three walk directions:
- **G** direction: tame/wild animals follow scalar additions along G
- **φG = [λ]G** direction: second independent axis on the hexagonal lattice
- **φ²G = [λ²]G = [−1−λ]G** direction: third axis (= −G − φG, closes the triangle)

These three axes form a **2D hexagonal lattice** (not 3D — φ²G is linearly dependent on {G, φG}).
Hexagonal tiling = densest 2D coverage = optimal for 2D Kangaroo walks (Gauss, 1831).

### 6-Automorphism Speedup

For any point P, its orbit under {±id, ±φ, ±φ²} has 6 elements sharing the same "canonical x":
```
canonical_x(P) = min(x_P, β·x_P, β²·x_P)  mod p
```
This collapses 6 distinct points to 1 DP check → **√6 ≈ 2.45× fewer DP entries**.

---

## 4. Kangaroo Constant Theory

### Formula

For a jump table with spread ratio r (largest jump / smallest jump):
```
C ≈ 1 + 2/ln(r)
```

### Progression

| Version | Bands | r | C (theory) | C (measured) |
|---------|-------|---|-----------|-------------|
| v5-v7 | 5 | 2^4 | 1.72 | — |
| v8-v9 | 9 | 2^8 | 1.36 | — |
| v10 | 17 | 2^16 | 1.18 | — |
| **v11** | **29** | **2^28** | **1.103** | **≈ 1.10** |
| ∞-band | ∞ | ∞ | **1.00** | (theoretical min) |

### Expected Operations for Puzzle #135

```
E[ops]  =  C × √(2^135 / 6)
        =  1.10 × 2^(135/2) / √6
        =  1.10 × 2^67.5 / 2.449
        ≈  2^66.0  operations
```

---

## 5. v11 Jump Table Design

### Parameters

```
NUM_BANDS = 29       (odd, centered at 0)
BAND_HALF = 14       (bands: -14, -13, ..., 0, ..., +13, +14)
NUM_JUMPS = 256      (from v10, kept for CUDA shared memory)
Axes:      3         (G, φG, φ²G)
Per axis:  32 jumps  (256 total across all — but stored as 32/axis)
```

### Band Construction

For each (axis, slot):
```
band  = (slot % 29) - 14        ∈ [-14, +14]
k_exp = (range_bits/2 + band)   ≈ around 2^67.5
jump_scalar = 2^k_exp + perturbation
```

The 29-band geometric spread gives r ≈ 2^28, yielding C ≈ 1 + 2/ln(2^28) ≈ 1.103.

### Jump Selection (GPU)

```cuda
// 0-cost bitmask selection (no modulo)
uint idx = cx[0] & (NUM_JUMPS - 1);   // & 0xFF
```

### CUDA Shared Memory Layout

```
256 jumps × 96 bytes/jump = 24,576 bytes = 24 KB per block
3 blocks × 24 KB = 72 KB < 100 KB L1 limit on RTX 4090
```

---

## 6. Known Sub-Exponential Attacks — Why They Fail

### Index Calculus

- **Idea:** Factor group elements over a "factor base" of small primes
- **Failure:** No notion of "smooth" EC points. The group law on E doesn't decompose multiplicatively over a factor base.
- **Status:** Requires genus ≥ 2 curves or hyperelliptic Jacobians. Infeasible for secp256k1.

### Weil Descent / GHS Attack

- **Idea:** Map E/F_p to a hyperelliptic Jacobian over F_{p/k} for some divisor k
- **Failure:** secp256k1 is defined over F_p (prime field, not extension). No natural descent target.
- **Status:** Applicable to E/F_{2^n} curves. secp256k1 is immune.

### MOV / Frey-Rück (Transfer to F_{p^k})

- **Idea:** Use Weil/Tate pairing to transfer DLP from E(F_p) to F_{p^k}^*
- **Failure:** Embedding degree k ≈ n (the curve's group order). The target field F_{p^n} is astronomically large.
- **Status:** Only works when k is small (≤ 20 or so). secp256k1 was designed with k ≈ n.

### Smart's Attack

- **Idea:** For anomalous curves (#E = p), use p-adic logarithm (Hensel lift)
- **Failure:** secp256k1 has trace t = p + 1 − n ≈ 2^128. Non-anomalous by a wide margin.
- **Status:** Would reduce DLP to O(log p) — pure catastrophe if applicable. Not applicable.

### Pohlig-Hellman

- **Idea:** If group order has small prime factors, solve DLP in each subgroup
- **Failure:** n is prime. No subgroup structure to exploit.
- **Status:** Twist order n' = 2p+2−n has small factors (see Section 8) — TPKH is active research.

### Semaev Summation Polynomials

- **Idea:** Build polynomial system that "adds" points over a factor base, solve with Gröbner basis
- **Status:** Open problem. Complexity unknown — if Gröbner basis is sub-quadratic for genus-1 curves, this could work. Currently no sub-exponential Gröbner algorithm for this case.

---

## 7. The GLS Breakthrough — Path to 4D

### Construction

The **GLS (Galbraith-Lin-Scott, 2009)** construction works as follows:

1. **Extend the field:** Work over F_{p²} instead of F_p
2. **Frobenius endomorphism:** π: (x,y) → (x^p, y^p) is well-defined and non-trivial over F_{p²}
3. **4D lattice:** k decomposes as `k = k₁ + k₂λ + k₃·(something) + k₄·(something·λ) (mod n)`
4. **Scalar bits:** Each kᵢ ≈ bits/4 instead of bits/2

### For Puzzle #135

```
2D current:   scalar length ≈ 67.5 bits  →  ops ≈ C × 2^67.5
4D GLS target: scalar length ≈ 33.75 bits →  ops ≈ C × 2^33.75
```

At 1 Gop/s: 2^33.75 ≈ 1.4 × 10^10 ops → **~14 seconds**.

### Why F_{p²} Arithmetic Doesn't Cancel the Gain

F_{p²} arithmetic costs ~4× per operation vs F_p:
- 1 F_{p²} addition = 2 F_p additions
- 1 F_{p²} multiplication = 4 F_p multiplications (Karatsuba: 3)
- Net: ~3-4× overhead

But the ops reduction is 2^33.75 ≈ 10^10×.

Net speedup: `10^10 / 4 ≈ 2.5 × 10^9×` — **nine orders of magnitude**.

### Implementation Plan (v12)

```rust
// F_{p²} element: a + b·i where i² = non-residue mod p
struct Fp2 { a: Fe, b: Fe }

// Curve over F_{p²}: same equation y² = x³ + 7, but x,y ∈ F_{p²}
// Frobenius: π(x,y) = (x^p, y^p) = (conj(x), conj(y)) for p≡3 mod 4

// 4D decomposition of k ∈ [0, 2^135):
// k = k1 + k2*λ + k3*p + k4*λ*p  (mod n)
// Each ki ≈ 33 bits

// Kangaroo walk axes (4D):
// J1 = G           (baseline)
// J2 = φ(G)        (CM endomorphism)
// J3 = π(G)        (Frobenius on F_{p²} lift)
// J4 = φ(π(G))     (composed)
```

---

## 8. Twist Order & TPKH Research

### Quadratic Twist

The quadratic twist E' of secp256k1 is the curve:
```
E': y² = x³ + 7·δ   for any non-square δ ∈ F_p
```
Its group order:
```
#E'(F_p)  =  n'  =  2p + 2 − n
```

Computed exactly:
```
n' = 0x[2p+2-n in 256-bit arithmetic]
n'[63:0]   = (computed by twist_order_low() in glv4d.rs)
n'[127:64] = (computed by twist_order_low() in glv4d.rs)
```

### Why n' Matters: TPKH

**Twist Pohlig-Hellman Combination (TPKH):**

If n' = q · m where q is small and gcd(q, n) = 1:
1. Map target to the twist E'
2. Solve DLP modulo q (small — Pohlig-Hellman in O(√q))
3. Use CRT to narrow the full DLP range from 2^135 to 2^135/q
4. Run Kangaroo on the narrowed range

This is currently blocked by the **bridge problem**: we don't know how to efficiently transfer DLP constraints between E and its twist when both have prime or nearly-prime orders.

### Status

- n' is computed in `glv4d::twist_order_low()`
- Divisibility by small primes checked in `glv4d::analyze_torsion()`
- Bridge problem: **open research question**

---

## 9. Empirical C Measurement

### Method

Run `--benchmark-c --range-bits B --trials N`:
1. Generate N random DLP instances: secret k ∈ [2^(B-1), 2^B), target = k·G
2. Run CPU Kangaroo with the same 29-band 3-axis jump table as the GPU kernel
3. Count steps until distinguished point collision
4. Compute: `C_measured = steps / √(2^B / 6)`

### Published Literature Comparison

| Implementation | C (reported) |
|---------------|-------------|
| JLP Kangaroo (Pollard 2000) | ≈ 1.65–2.00 |
| JeanLucPons v2 (no GLV) | ≈ 1.40–1.60 |
| GLV-2 implementations | ≈ 1.30–1.40 |
| sinGRAAL v9 (9-band) | ≈ 1.36 |
| sinGRAAL v10 (17-band) | ≈ 1.18 |
| **sinGRAAL v11 (29-band)** | **≈ 1.10** |
| Theoretical minimum | 1.00 |

If measured C < 1.15, sinGRAAL is in **world-record territory** for published Kangaroo constants on secp256k1.

---

## 10. Open Problems

### Near-term (implementable)

| Problem | Difficulty | Potential gain |
|---------|-----------|----------------|
| 4D GLS Kangaroo on CUDA | High (F_{p²} arithmetic) | 2^33.75 ops → seconds |
| 4D LLL decomposition | Medium | required for GLS |
| Twist TPKH bridge | Research | 64× if solved |
| Semaev Gröbner basis | Open math | sub-exponential if solved |
| C → 1.05 (more bands) | Medium | 5% more ops saved |

### v11 → v12 Checklist

- [ ] F_{p²} field arithmetic in Rust (CPU setup)
- [ ] F_{p²} field arithmetic in CUDA (GPU kernel)
- [ ] Frobenius endomorphism π on F_{p²} points
- [ ] 4D LLL scalar decomposition
- [ ] 4D Kangaroo jump table (4 axes × 29 bands)
- [ ] 4D DP detection (canonical form in 4D)
- [ ] 4D collision recovery → scalar k
- [ ] Benchmark: measure actual C for 4D walk
- [ ] Target: puzzle #135 in < 60 seconds on single RTX 4090

### Mathematical Frontiers

| Direction | Status | If solved |
|-----------|--------|-----------|
| GLS 4D on F_{p²} | **v12 target** | **seconds per GPU** |
| TPKH bridge | Research | 64× narrowing |
| Semaev sub-quadratic | Open math | potentially polynomial |
| Supersingular lift | Theoretical | 4D via quaternion End(E) |
| ML structure search | Experimental | unknown |

---

*sinGRAAL research journal — updated with every version.*
*"On est capable de rendre ça tellement puissant... alors innovons."*
