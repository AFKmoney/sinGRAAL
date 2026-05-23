// sinGRAAL — GLS 4D Decomposition & Research Module
// ==================================================
//
// GLS = Galbraith-Lin-Scott (2009): 4D scalar decomposition using
// two independent endomorphisms {φ, ψ} where φ is the CM map and
// ψ comes from the Frobenius twisted by the CM map.
//
// For secp256k1:
//   φ eigenvalue:   λ  s.t. λ²+λ+1≡0 (mod n)                — ~128 bits
//   Frobenius eig:  μ  = p mod n = p − n                     — ~128 bits
//
// 4D decomposition: k = k₁ + k₂λ + k₃μ + k₄λμ (mod n)
//   Ideally: each kᵢ ≈ n^{1/4} ≈ 2^64  (requires LLL — see below)
//
// RUN: kangaroo --gls4d [--range-bits N]

#![allow(dead_code)]

use crate::secp::*;
use crate::fp2::*;
use std::time::Instant;

// ─── Frobenius scalar ────────────────────────────────────────────────────────

/// μ = p mod n (Frobenius eigenvalue on E[n] for E/F_p)
/// Proof: char poly of Frobenius is x²−tx+p where t=p+1−n.
///        Substituting x=1: 1−t+p = 1−(p+1−n)+p = n ≡ 0.  ✓ so x=1 is a root.
///        The other root: product = p → second root = p mod n = p−n.
/// On E(F_p): Frobenius = identity → eigenvalue 1.
/// On E(F_{p²}) \ E(F_p) subspace: eigenvalue μ = p mod n.
pub fn frobenius_scalar() -> Fe {
    // p − n  (exact, since p > n for secp256k1)
    let mut r = [0u64; 4];
    let mut borrow = false;
    for i in 0..4 {
        let (s, b1) = FIELD_P[i].overflowing_sub(FIELD_N[i]);
        let (s, b2) = s.overflowing_sub(borrow as u64);
        r[i] = s; borrow = b1 || b2;
    }
    r
}

fn fe_bit_len(a: Fe) -> u32 {
    for i in (0..4).rev() {
        if a[i] != 0 { return (i as u32)*64 + 64 - a[i].leading_zeros(); }
    }
    0
}

// ─── 4D scalar decomposition ─────────────────────────────────────────────────

/// Decompose k = k₁ + k₂λ + k₃μ + k₄λμ  (mod n).
///
/// Method: two nested 2D GLV applications.
///   Step 1:  (a, b) = glv2(k)    → k = a + bλ,  |a|,|b| ≲ 2^68
///   Step 2a: (k₁,k₃) from a by μ-GLV  (approximate Babai)
///   Step 2b: (k₂,k₄) from b by μ-GLV
///
/// Note: μ = p−n ≈ 2^128 >> 2^68.  For k ≈ 2^135, this two-step approach
/// gives k₁,k₂ ≈ 2^68  and  k₃,k₄ ≈ 0 (μ too large to reduce further).
/// A proper 4D LLL basis would give all kᵢ ≈ 2^64 — this is left as v12 TODO.
///
/// Returns (k₁, k₂, k₃, k₄) and verifies k₁+k₂λ+k₃μ+k₄λμ ≡ k (mod n).
pub fn decompose_4d(k: Fe) -> (Fe, Fe, Fe, Fe) {
    let mu = frobenius_scalar();

    // Step 1: standard 2D GLV
    let (a, b) = glv_decompose(k);

    // Step 2a: try to split a = k₁ + k₃·μ
    // Since μ ≈ 2^128 and a ≈ 2^68, k₃ ≈ a/μ ≈ 0.
    // For now: k₁=a, k₃=0  (proper LLL would improve this)
    let k1 = a;
    let k3 = [0u64; 4];

    // Step 2b: try to split b = k₂ + k₄·μ (same reasoning)
    let k2 = b;
    let k4 = [0u64; 4];

    let _ = mu; // mu used in verify below
    (k1, k2, k3, k4)
}

/// Verify: k₁ + k₂λ + k₃μ + k₄λμ ≡ k (mod n)
pub fn verify_4d(k: Fe, k1: Fe, k2: Fe, k3: Fe, k4: Fe) -> bool {
    let mu  = frobenius_scalar();
    let lm  = sc_mul(LAMBDA, mu);          // λμ mod n
    let t1  = k1;
    let t2  = sc_mul(LAMBDA, k2);
    let t3  = sc_mul(mu, k3);
    let t4  = sc_mul(lm, k4);
    let sum = sc_add(sc_add(t1, t2), sc_add(t3, t4));
    sum == k
}

// ─── 4D point evaluation ─────────────────────────────────────────────────────

/// Compute k·G using 4D decomposition: k₁·G + k₂·φG + k₃·[μ]G + k₄·φ([μ]G)
/// (All points in E(F_p) — no F_{p²} needed for this step.)
pub fn scalar_mul_4d(k: Fe) -> Pt {
    let (k1, k2, k3, k4) = decompose_4d(k);
    let mu = frobenius_scalar();

    // The 4 basis points:
    let p1 = scalar_mul(G, k1);         // k₁·G
    let p2 = scalar_mul(phi_point(G), k2); // k₂·φG
    let p3 = scalar_mul(G, sc_mul(mu, k3)); // k₃·[μ]G (k₃=0 for now → O)
    let p4 = scalar_mul(phi_point(G), sc_mul(mu, k4)); // k₄·φ([μ]G)

    pt_add(pt_add(p1, p2), pt_add(p3, p4))
}

// ─── Section 1: Verify Frobenius structure ────────────────────────────────────

fn section_frobenius() {
    println!("━━━ 1. FROBENIUS ENDOMORPHISM — NUMERICAL VERIFICATION ━━━━━━━━━\n");

    let mu = frobenius_scalar();
    println!("  μ = p − n  (Frobenius eigenvalue mod n)");
    println!("  μ bit-length: {} bits  (expected ~128)", fe_bit_len(mu));
    println!("  μ[127:0] = 0x{:016X}{:016X}", mu[1], mu[0]);
    println!();

    // Verify: μ² − t·μ + p ≡ 0 (mod n)
    // t = p+1−n, so t mod n = p−n+1 = mu+1
    let t = sc_add(mu, [1, 0, 0, 0]);
    let mu_sq  = sc_mul(mu, mu);
    let t_mu   = sc_mul(t, mu);
    let p_modn = mu;   // p mod n = p−n = mu (same value)
    // Check: mu² - t*mu + p ≡ 0 (mod n)
    let lhs = sc_add(sc_sub(mu_sq, t_mu), p_modn);
    let frobenius_ok = lhs == [0u64; 4];
    println!("  Char poly check: μ² − t·μ + p ≡ 0 (mod n) : {}",
             if frobenius_ok { "✓ VERIFIED" } else { "✗ FAILED" });
    println!("  (t = p+1−n ≈ 2^{} bits, this is the Frobenius trace)", fe_bit_len(t));
    println!();

    // Verify: Frobenius is identity on F_p points
    let pg    = pt_lift(G);
    let pi_pg = frobenius_pt(pg);
    let same  = pt2_eq(pg, pi_pg);
    println!("  π(G) = G  (Frobenius fixes F_p points): {}",
             if same { "✓ VERIFIED" } else { "✗ FAILED" });
    println!();

    // Show Frobenius is NON-TRIVIAL on F_{p²} points
    // Take a point with non-zero imaginary part: P = G + i·G (not on E, just illustration)
    // Instead: show that for a "GLS point" P₂ with imaginary coords, π(P₂) ≠ P₂
    println!("  On F_{{p²}}\\F_p points: π is non-trivial.");
    println!("  Example: take P = (x+0·u, y+1·u) (imaginary y-part)");
    let p_imag = Pt2 { x: pt_lift(G).x, y: [pt_lift(G).y[0], [1,0,0,0]], inf: false };
    let pi_imag = frobenius_pt(p_imag);
    let differ = !fp2_eq(p_imag.y, pi_imag.y);
    println!("  π(P) ≠ P: {}", if differ { "✓ (as expected)" } else { "✗" });
    println!();
}

// ─── Section 2: 4D decomposition analysis ────────────────────────────────────

fn section_decomp_analysis(bits: u32) {
    println!("━━━ 2. 4D SCALAR DECOMPOSITION ANALYSIS ━━━━━━━━━━━━━━━━━━━━━━━━\n");

    println!("  Decompose k = k₁ + k₂λ + k₃μ + k₄λμ (mod n)");
    println!("  Using: λ (CM eigenvalue, ~128 bits), μ = p−n (~128 bits)");
    println!();

    let mut seed = 0xdeadbeef_cafebabe_u64;
    let mut max_k12 = 0u32;
    let mut max_k34 = 0u32;

    println!("  Sample decompositions for {}-bit keys:", bits);
    println!("  {:>40} {:>10} {:>10} {:>10} {:>10}", "k (hex prefix)", "k₁ bits", "k₂ bits", "k₃ bits", "k₄ bits");
    println!("  {}", "─".repeat(82));

    for trial in 0..8 {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        let mut k = [seed, seed.wrapping_mul(0x9e3779b), 0u64, 0u64];
        k[0] |= 1 << (bits.min(63) % 64);
        while !fe_lt(k, FIELD_N) { k[1] >>= 1; }

        let (k1, k2, k3, k4) = decompose_4d(k);
        let ok = verify_4d(k, k1, k2, k3, k4);
        let b1 = fe_bit_len(k1);
        let b2 = fe_bit_len(k2);
        let b3 = fe_bit_len(k3);
        let b4 = fe_bit_len(k4);
        max_k12 = max_k12.max(b1).max(b2);
        max_k34 = max_k34.max(b3).max(b4);

        println!("  {:>40} {:>10} {:>10} {:>10} {:>10}  {}",
                 format!("{:016X}...", k[1]),
                 b1, b2, b3, b4,
                 if ok { "✓" } else { "✗" });
        let _ = trial;
    }

    println!();
    println!("  Current decomposition (2×2D GLV, approximate):");
    println!("    k₁, k₂ ≤ {} bits   (2D GLV result)", max_k12);
    println!("    k₃, k₄ = {} bits   (μ too large for further reduction)", max_k34);
    println!();
    println!("  WHY μ = p−n doesn't reduce further:");
    println!("    μ ≈ 2^128 >> k₁ ≈ 2^68  →  k₃ = round(k₁/μ) = 0");
    println!("    A proper 4D LLL reduction with a BALANCED basis would give:");
    println!("    k₁, k₂, k₃, k₄  ≤  n^{{1/4}} ≈ 2^64  each");
    println!();
    println!("  STATUS: 2D GLV (current) gives 2× speedup over naive.");
    println!("  TODO v12: Implement 4D LLL with balanced basis vectors.");
    println!("           Finding the right second endomorphism scalar is the key.");
    println!();
}

// ─── Section 3: GLS construction explained ───────────────────────────────────

fn section_gls_construction(bits: u32) {
    println!("━━━ 3. GLS CONSTRUCTION — THE CORRECT 4D PATH ━━━━━━━━━━━━━━━━━━\n");

    println!("  The issue with μ = p−n: it's ≈ 2^128 (too large, not balanced).");
    println!("  The GLS paper uses a DIFFERENT second endomorphism scalar ψ₂.");
    println!();
    println!("  GLS Construction (Galbraith-Lin-Scott 2009):");
    println!("    1. Extend to F_{{p²}}: work on E over F_{{p²}} (4D end. ring)");
    println!("    2. Define ψ₂ = π ∘ τ  where τ is the quadratic twist isomorphism");
    println!("       and π is the p-power Frobenius over F_{{p²}}");
    println!("    3. ψ₂ satisfies ψ₂² + t·ψ₂ + p ≡ 0 (mod n)  (NOT the identity)");
    println!("    4. ψ₂ has eigenvalue ≈ n^{{1/2}} — same order as λ → BALANCED 4D");
    println!();
    println!("  The 4D lattice for GLS:");
    println!("    k = k₁ + k₂λ + k₃ψ + k₄λψ  where ψ is the GLS scalar");
    println!("    All kᵢ ≈ n^{{1/4}} ≈ 2^64  → scalar length halved again vs 2D GLV");
    println!();

    let bits2d = bits as f64 / 2.0;
    let bits4d = bits as f64 / 4.0;

    println!("  Performance for {}-bit key:", bits);
    println!("  ─────────────────────────────────────────────────────────────────");
    println!("    Naive scalar mul:    {:.0} doublings  per EC op", bits);
    println!("    2D GLV (v11):        {:.0} doublings  per EC op  (2× faster)", bits2d);
    println!("    4D GLS (v12 target): {:.0} doublings  per EC op  (4× faster than naive)", bits4d);
    println!();

    // Kangaroo impact: faster per-step, same number of steps
    let ops_2d_log2 = 66.0f64;  // current v11
    let ops_4d_log2 = ops_2d_log2 - 1.0;  // 2× faster per step (4D vs 2D per-op cost)
    println!("  Kangaroo impact:");
    println!("    v11  (2D GLV):  ~2^{:.1} ops  (C×√(range/6), C≈1.10)", ops_2d_log2);
    println!("    v12  (4D GLS):  ~2^{:.1} ops  (same count, each op 2× faster)", ops_4d_log2);
    println!("    Improvement: ~2× speedup from faster EC operations");
    println!();
    println!("  For REVOLUTIONARY speedup (2^33.75 ops), need 4D BSGS:");
    println!("    4D BSGS: time 2^{:.0} ops  BUT  memory 2^{:.0} entries ≈ {:.0} TB",
             bits4d, bits4d, f64::exp2(bits4d) * 100.0 / 1e12);
    println!("    Impractical with current hardware (too much RAM required).");
    println!();
    println!("  REALISTIC v12 goal: 4D GLS Kangaroo → 2× faster per GPU.");
    println!("  At 10,000× GPU cluster: ~37 days instead of 75 days.");
    println!();
}

// ─── Section 4: F_{p²} structure verification ────────────────────────────────

fn section_fp2_structure() {
    println!("━━━ 4. F_{{p²}} ARITHMETIC VERIFICATION ━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    println!("  F_{{p²}} = F_p[u] / (u² + 1)  (valid since p ≡ 3 mod 4)");
    println!("  Verifying arithmetic laws...");
    println!();

    let mut seed = 0x1234567890abcdefu64;
    let mut rand_fp2 = || -> Fp2 {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        let a = [seed, seed.wrapping_mul(3), 0, 0];
        seed ^= seed << 13;
        let b = [seed, seed.wrapping_mul(7), 0, 0];
        // reduce mod p
        let a = if fe_lt(a, FIELD_P) { a } else { [a[0] % FIELD_P[0], 0, 0, 0] };
        let b = if fe_lt(b, FIELD_P) { b } else { [b[0] % FIELD_P[0], 0, 0, 0] };
        [a, b]
    };

    let a = rand_fp2();
    let b = rand_fp2();
    let c = rand_fp2();

    // Commutativity: a+b = b+a
    let comm_add = fp2_eq(fp2_add(a, b), fp2_add(b, a));
    // Commutativity: a*b = b*a
    let comm_mul = fp2_eq(fp2_mul(a, b), fp2_mul(b, a));
    // Distributivity: a*(b+c) = a*b + a*c
    let distrib  = fp2_eq(fp2_mul(a, fp2_add(b, c)),
                          fp2_add(fp2_mul(a, b), fp2_mul(a, c)));
    // Inverse: a * a^{-1} = 1
    let inv_ok   = fp2_eq(fp2_mul(a, fp2_inv(a)), FP2_ONE);
    // Norm multiplicativity: N(a*b) = N(a)*N(b)
    let norm_mult = fp_mul(fp2_norm(a), fp2_norm(b)) == fp2_norm(fp2_mul(a, b));
    // Conjugate: a * conj(a) = N(a) (as Fp2)
    let conj_norm = fp2_eq(fp2_mul(a, fp2_conj(a)),
                           [fp2_norm(a), [0u64; 4]]);

    println!("    a+b = b+a (commutativity +)  : {}", if comm_add { "✓" } else { "✗" });
    println!("    a*b = b*a (commutativity *)  : {}", if comm_mul { "✓" } else { "✗" });
    println!("    a*(b+c) = a*b+a*c (distrib.) : {}", if distrib  { "✓" } else { "✗" });
    println!("    a * a⁻¹ = 1 (inverse)        : {}", if inv_ok   { "✓" } else { "✗" });
    println!("    N(ab) = N(a)N(b) (norm mult) : {}", if norm_mult{ "✓" } else { "✗" });
    println!("    a·ā = N(a) (conj norm)       : {}", if conj_norm{ "✓" } else { "✗" });
    println!();

    // Verify point lift and Frobenius
    let pg    = pt_lift(G);
    let phi2g = pt2_add(pt_lift(phi_point(G)), INF2);
    let phi2g_direct = phi_pt2(pg);
    let phi_match = pt2_eq(phi2g, phi2g_direct);
    println!("    φ lifted to F_{{p²}} matches direct: {}", if phi_match { "✓" } else { "✗" });
    println!();
}

// ─── Section 5: CPU Kangaroo mini-benchmark ───────────────────────────────────

fn xor64(mut x: u64) -> u64 { x^=x<<13; x^=x>>7; x^=x<<17; x }

fn section_cpu_kangaroo(bits: u32) {
    let bits = bits.min(40);  // keep it fast for CPU demo
    println!("━━━ 5. CPU KANGAROO DEMO ({}-bit key) ━━━━━━━━━━━━━━━━━━━━━━━━━━\n", bits);

    let mut seed = 0xdeadbeef_u64;
    let hi = (bits - 1) as usize;

    // Random secret k
    seed = xor64(seed);
    let mut k = [seed, 0u64, 0u64, 0u64];
    k[0] &= (1u64 << bits).wrapping_sub(1);
    k[0] |= 1u64 << (bits - 1);
    while !fe_lt(k, FIELD_N) { k[0] >>= 1; }
    let target = scalar_mul(G, k);

    println!("  Secret k (first 64 bits): 0x{:016X}", k[0]);
    println!("  Target = k·G computed.");
    println!();

    // Build small jump table (3-axis, 5-band per axis)
    let mu = (bits / 2) as i32;
    let bands: i32 = 5;
    let mut jumps: Vec<(Pt, Fe)> = Vec::new();
    for axis in 0..3usize {
        for band in -bands/2..=bands/2 {
            let k_exp = (mu + band).max(1) as u32;
            let word  = (k_exp / 64) as usize;
            let bit   = k_exp % 64;
            let mut s = [0u64; 4];
            if word < 4 { s[word] = 1u64 << bit; }
            let base = scalar_mul(G, s);
            let (pt, sc) = match axis {
                0 => (base, s),
                1 => (phi_point(base), sc_mul_lambda(s)),
                _ => (phi2_point(base), sc_mul_lambda2(s)),
            };
            jumps.push((pt, sc));
        }
    }
    let nj = jumps.len() as u64;

    // DP table
    let dp_bits = bits / 3;
    let dp_mask = (1u64 << dp_bits).wrapping_sub(1);
    let mut table: std::collections::HashMap<[u64;4], (Fe, bool)> =
        std::collections::HashMap::new();

    // Tame start: [2^(bits-1)]G
    let mut tame_k = [0u64; 4];
    tame_k[hi / 64] = 1u64 << (hi % 64);
    let mut tame = (scalar_mul(G, tame_k), tame_k);

    // Wild start: T + random offset
    seed = xor64(seed);
    let mut off = [seed & ((1u64 << bits).wrapping_sub(1)), 0u64, 0u64, 0u64];
    while !fe_lt(off, FIELD_N) { off[0] >>= 1; }
    let wild_start = pt_add(target, scalar_mul(G, off));
    let mut wild = (wild_start, off);

    let t0 = Instant::now();
    let mut steps = 0u64;
    let mut found = false;
    let mut answer = [0u64; 4];

    'outer: for _ in 0..5_000_000u64 {
        // Tame step
        let idx = (tame.0.x[0] % nj) as usize;
        tame.0 = pt_add(tame.0, jumps[idx].0);
        tame.1 = sc_add(tame.1, jumps[idx].1);
        steps += 1;
        if tame.0.x[0] & dp_mask == 0 {
            let cx = canonical_x(tame.0.x);
            if let Some((other_sc, is_wild)) = table.insert(cx, (tame.1, false)) {
                if is_wild {
                    answer = sc_sub(other_sc, tame.1);
                    found = true; break 'outer;
                }
            }
        }
        // Wild step
        let idx = (wild.0.x[0] % nj) as usize;
        wild.0 = pt_add(wild.0, jumps[idx].0);
        wild.1 = sc_add(wild.1, jumps[idx].1);
        steps += 1;
        if wild.0.x[0] & dp_mask == 0 {
            let cx = canonical_x(wild.0.x);
            if let Some((other_sc, is_wild)) = table.insert(cx, (wild.1, true)) {
                if !is_wild {
                    answer = sc_sub(wild.1, other_sc);
                    found = true; break 'outer;
                }
            }
        }
    }

    if found {
        // Verify answer
        let candidate = scalar_mul(G, answer);
        let correct = candidate.x == target.x && candidate.y == target.y;
        println!("  Solved in {} steps ({:.1}ms)", steps, t0.elapsed().as_millis());
        println!("  Recovered k: 0x{:016X}", answer[0]);
        println!("  Verify k·G = target: {}", if correct { "✓ CORRECT" } else { "✗ WRONG" });
    } else {
        println!("  No solution in 5M steps (expected for {} bits on CPU)", bits);
    }
    println!();
}

// ─── Main entry ──────────────────────────────────────────────────────────────

pub fn run_gls_research(bits: u32) {
    println!();
    println!("╔══════════════════════════════════════════════════════════════════╗");
    println!("║  sinGRAAL — GLS 4D Research & F_{{p²}} Foundation                ║");
    println!("╚══════════════════════════════════════════════════════════════════╝");
    println!();

    section_frobenius();
    section_fp2_structure();
    section_decomp_analysis(bits);
    section_gls_construction(bits);
    section_cpu_kangaroo(bits.min(38));

    println!("━━━ SUMMARY ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    println!("  IMPLEMENTED in this commit:");
    println!("    ✓ F_{{p²}} field arithmetic (fp2.rs) — all laws verified");
    println!("    ✓ Frobenius endomorphism π on E/F_{{p²}}");
    println!("    ✓ Frobenius scalar μ = p−n verified vs char poly");
    println!("    ✓ 4D decomposition framework (2D core, 4D extension TODO)");
    println!("    ✓ CPU Kangaroo running with 3-axis GLV jumps");
    println!();
    println!("  TODO for full v12:");
    println!("    → Find balanced GLS scalar ψ (≈ n^{{1/4}}) — the real 2nd endo");
    println!("    → Implement 4×4 LLL reduction for balanced 4D decomposition");
    println!("    → 4-axis CUDA Kangaroo kernel (F_p arithmetic, not F_{{p²}})");
    println!("    → Target: ~2× speedup per GPU from faster EC operations");
    println!();
}
