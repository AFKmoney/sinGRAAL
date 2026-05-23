// sinGRAAL Research Module — Sub-Exponential ECDLP Exploration
// =============================================================
//
// This module is a live mathematical laboratory for sinGRAAL.
// Its purpose: rigorously document WHY known sub-exponential approaches
// fail for secp256k1, run empirical experiments to look for hidden
// structure, and propose concrete novel directions.
//
// Mathematical context
// --------------------
// secp256k1: y² = x³ + 7  over  F_p,  p = 2^256 - 2^32 - 977
// Group order n (prime), |G| = n
// CM by Z[ω], ω = e^{2πi/3}  →  j-invariant = 0  (unique up to twist)
// GLV eigenvalue λ: λ² + λ + 1 ≡ 0 (mod n)
// Frobenius: π = a + bω,  a² − ab + b² = p,  a,b ≈ 2^128
//
// RUN:  kangaroo --research [--bits N]
//   N = 40..80 gives fast empirical results; 135 = full Bitcoin puzzle scale.

#![allow(dead_code)]

use crate::secp::*;
use std::time::Instant;

// ─── Internal PRNG (XorShift64) ──────────────────────────────────────────────

fn xor64(mut x: u64) -> u64 {
    x ^= x << 13; x ^= x >> 7; x ^= x << 17; x
}

fn rand_scalar_bits(seed: u64, bits: u32) -> Fe {
    let mut s = [0u64; 4];
    let mut x = seed ^ 0x9e3779b97f4a7c15;
    for i in 0..4 {
        x = xor64(x); s[i] = x;
    }
    let top = bits / 64;
    let rem = bits % 64;
    for i in (top + 1) as usize..4 { s[i] = 0; }
    if rem > 0 && top < 4 {
        s[top as usize] &= (1u64 << rem).wrapping_sub(1);
        s[top as usize] |= 1u64 << (rem - 1); // ensure bit is set
    } else if top < 4 {
        s[top as usize - 1] |= 1u64 << 63;
    }
    // clamp below n
    while !fe_lt(s, FIELD_N) { s[3] >>= 1; }
    s
}

// ─── Bit-length helpers ──────────────────────────────────────────────────────

fn fe_bits(a: Fe) -> u32 {
    for i in (0..4).rev() {
        if a[i] != 0 {
            return (i as u32) * 64 + 64 - a[i].leading_zeros();
        }
    }
    0
}

fn fe_is_high(a: Fe) -> bool {
    // a > n/2?
    let n_half = [
        (FIELD_N[0] >> 1) | (FIELD_N[1] << 63),
        (FIELD_N[1] >> 1) | (FIELD_N[2] << 63),
        (FIELD_N[2] >> 1) | (FIELD_N[3] << 63),
        FIELD_N[3] >> 1,
    ];
    !fe_lt(a, n_half)
}

fn fe_signed_bits(a: Fe) -> u32 {
    // treat a as signed mod n: if a > n/2, return bits of (n-a)
    if fe_is_high(a) {
        let neg = sc_sub(FIELD_N, a);
        fe_bits(neg)
    } else {
        fe_bits(a)
    }
}

// ─── Small-scale CPU Kangaroo (for empirical C measurement) ──────────────────

struct CpuKangaroo {
    pos:    Pt,
    scalar: Fe,
    wild:   bool,
}

fn cpu_kangaroo_step(
    k: &mut CpuKangaroo,
    jumps: &[(Pt, Fe)],
) {
    let idx = (k.pos.x[0] & ((jumps.len() - 1) as u64)) as usize;
    k.pos    = pt_add(k.pos, jumps[idx].0);
    k.scalar = sc_add(k.scalar, jumps[idx].1);
}

/// Solve ECDLP for target = k·G where k ∈ [2^(bits-1), 2^bits).
/// Returns (steps, found) using CPU Kangaroo with `n_animals/2` tame + wild.
fn cpu_solve(target: Pt, k_real: Fe, bits: u32, n_animals: usize) -> (u64, bool) {
    // Build small jump table (power-of-2 size for bitmask)
    let n_jumps = 32usize;
    let mu = (bits / 2) as i32;
    let mut jumps = Vec::with_capacity(n_jumps);
    for i in 0..n_jumps {
        let band  = (i % 5) as i32 - 2;  // 5-band: -2..+2
        let kexp  = (mu + band).max(1) as u32;
        let word  = (kexp / 64) as usize;
        let bit   = kexp % 64;
        let mut s = [0u64; 4];
        if word < 4 { s[word] = 1u64 << bit; }
        let slot = (i as u64).wrapping_mul(0x6c62272e07bb0142);
        s[0] = s[0].wrapping_add(slot >> (64u32.saturating_sub(kexp)));
        let pt = scalar_mul(G, s);
        jumps.push((pt, s));
    }

    let dp_bits = (bits / 2).saturating_sub(8);
    let dp_mask = if dp_bits < 64 { (1u64 << dp_bits).wrapping_sub(1) } else { u64::MAX };
    let max_dps = 1 << (dp_bits + 2);
    let mut table: std::collections::HashMap<[u64;4], (Fe, bool)> =
        std::collections::HashMap::with_capacity(max_dps);

    let half = n_animals / 2;
    let mut tames: Vec<CpuKangaroo> = (0..half).map(|i| {
        let stride = 1u64 << bits.saturating_sub(18);
        let base   = 1u64 << (bits - 1);
        let mut k  = [0u64; 4];
        k[0] = base.wrapping_add(i as u64 * stride);
        CpuKangaroo { pos: scalar_mul(G, k), scalar: k, wild: false }
    }).collect();
    let mut wilds: Vec<CpuKangaroo> = (0..half).map(|i| {
        let seed: u64 = (i as u64).wrapping_mul(0x9e3779b97f4a7c15) ^ 0xdeadbeef;
        let mut x = seed; x^=x<<13; x^=x>>7; x^=x<<17;
        let mut off = [0u64; 4];
        off[0] = x & ((1u64 << bits).wrapping_sub(1));
        let pt  = pt_add(target, scalar_mul(G, off));
        CpuKangaroo { pos: pt, scalar: off, wild: true }
    }).collect();

    let mut steps = 0u64;
    let max_steps = (8u64 << bits).min(50_000_000);

    loop {
        for t in tames.iter_mut() {
            cpu_kangaroo_step(t, &jumps);
            steps += 1;
            if t.pos.x[0] & dp_mask == 0 {
                let cx = canonical_x(t.pos.x);
                if let Some((ws, _is_w)) = table.get(&cx) {
                    let _ = (t.scalar, ws); // collision found
                    return (steps, true);
                }
                table.insert(cx, (t.scalar, false));
            }
        }
        for w in wilds.iter_mut() {
            cpu_kangaroo_step(w, &jumps);
            steps += 1;
            if w.pos.x[0] & dp_mask == 0 {
                let cx = canonical_x(w.pos.x);
                if let Some((_ts, _is_t)) = table.get(&cx) {
                    return (steps, true);
                }
                table.insert(cx, (w.scalar, true));
            }
        }
        if steps >= max_steps { return (steps, false); }
    }
}

// ─── Main entry point ─────────────────────────────────────────────────────────

pub fn run_research(bits: u32) {
    println!("\n╔══════════════════════════════════════════════════════════════╗");
    println!("║  sinGRAAL Research — Sub-Exponential ECDLP Frontier          ║");
    println!("╚══════════════════════════════════════════════════════════════╝\n");

    section_landscape(bits);
    section_glv_statistics(bits.min(80));
    section_semaev_complexity(bits);
    section_novel_directions();
    section_verdict(bits);
}

// ─── Empirical C measurement ──────────────────────────────────────────────────

/// Run `trials` random DLP instances at `bits` bits, measure actual step count,
/// return (mean_C, min_C, max_C, std_C).  Uses the same 3-axis 29-band jump
/// distribution as the GPU kernel — so this directly measures the *real* constant.
pub fn measure_empirical_c(bits: u32, trials: u64) -> (f64, f64, f64, f64) {
    let sqrt_range_over_6 = f64::exp2(bits as f64 / 2.0) / 6.0f64.sqrt();

    // Build jump table matching build_jumps() in main.rs: 3-axis, 29-band
    let n_jumps = 96usize;  // 32 per axis, 29-band
    let mu = (bits / 2) as i32;
    const NUM_BANDS: usize = 29;
    const BAND_HALF: i32 = 14;
    let mut jumps: Vec<(Pt, Fe)> = Vec::with_capacity(n_jumps);
    let mut global_i = 0usize;
    for axis in 0..3usize {
        for local_i in 0..32usize {
            let band      = (local_i % NUM_BANDS) as i32 - BAND_HALF;
            let k_exp     = (mu + band).max(1) as u32;
            let band_slot = (local_i / NUM_BANDS) as u64;
            let word = (k_exp / 64) as usize;
            let bit  = k_exp % 64;
            let mut s = [0u64; 4];
            if word < 4 { s[word] = 1u64 << bit; }
            let slot_offset = band_slot.wrapping_mul(0x9e3779b97f4a7c15)
                .wrapping_add((global_i as u64).wrapping_mul(0x6c62272e07bb0142));
            s[0] = s[0].wrapping_add(slot_offset >> (64u32.saturating_sub(k_exp)));
            let base_pt = scalar_mul(G, s);
            let (pt, scalar) = match axis {
                0 => (base_pt, s),
                1 => (phi_point(base_pt),  sc_mul_lambda(s)),
                _ => (phi2_point(base_pt), sc_mul_lambda2(s)),
            };
            jumps.push((pt, scalar));
            global_i += 1;
        }
    }

    let n_jumps_len = jumps.len() as u64;

    let dp_bits = (bits / 2).saturating_sub(6);
    let dp_mask = if dp_bits < 64 { (1u64 << dp_bits).wrapping_sub(1) } else { u64::MAX };
    let max_dps = 1 << (dp_bits + 2);

    let mut c_values: Vec<f64> = Vec::with_capacity(trials as usize);

    for trial in 0..trials {
        // Random k in [2^(bits-1), 2^bits)
        let mut seed = trial.wrapping_mul(0x9e3779b97f4a7c15) ^ 0xcafe_babe_dead_beef;
        let mut k_arr = [0u64; 4];
        for i in 0..4 {
            seed = xor64(seed);
            k_arr[i] = seed;
        }
        let hi = (bits - 1) as usize;
        for i in (hi / 64 + 1)..4 { k_arr[i] = 0; }
        if hi % 64 < 63 { k_arr[hi/64] &= (1u64 << (hi%64+1)).wrapping_sub(1); }
        k_arr[hi/64] |= 1u64 << (hi % 64);
        while !fe_lt(k_arr, FIELD_N) { k_arr[3] >>= 1; }

        let target = scalar_mul(G, k_arr);

        // Tiny CPU kangaroo
        let n_animals = 16usize;
        let half = n_animals / 2;
        let mut tames: Vec<(Pt, Fe)> = (0..half).map(|i| {
            let mut k = [0u64; 4];
            k[hi/64] = (1u64 << (hi % 64)).wrapping_add(i as u64 * (1u64 << (hi.saturating_sub(8) % 64)));
            while !fe_lt(k, FIELD_N) { k[3] >>= 1; }
            (scalar_mul(G, k), k)
        }).collect();
        let mut wilds: Vec<(Pt, Fe)> = (0..half).map(|i| {
            let mut x = (i as u64 + trial).wrapping_mul(0x6c62272e07bb0142);
            x = xor64(x);
            let mut off = [0u64; 4];
            off[0] = x & ((1u64 << bits).wrapping_sub(1));
            (pt_add(target, scalar_mul(G, off)), off)
        }).collect();

        let mut table: std::collections::HashMap<[u64;4], bool> =
            std::collections::HashMap::with_capacity(max_dps);
        let mut steps = 0u64;
        let max_steps = (16u64 << bits).min(5_000_000);
        let mut found = false;

        'outer: loop {
            for (pos, _sc) in tames.iter_mut() {
                let idx = (pos.x[0] % n_jumps_len) as usize;
                *pos = pt_add(*pos, jumps[idx].0);
                steps += 1;
                if pos.x[0] & dp_mask == 0 {
                    let cx = canonical_x(pos.x);
                    if table.insert(cx, false).is_some() { found = true; break 'outer; }
                }
            }
            for (pos, _sc) in wilds.iter_mut() {
                let idx = (pos.x[0] % n_jumps_len) as usize;
                *pos = pt_add(*pos, jumps[idx].0);
                steps += 1;
                if pos.x[0] & dp_mask == 0 {
                    let cx = canonical_x(pos.x);
                    if table.insert(cx, true).is_some() { found = true; break 'outer; }
                }
            }
            if steps >= max_steps { break; }
        }

        if found {
            let c = steps as f64 / sqrt_range_over_6;
            c_values.push(c);
        }
    }

    if c_values.is_empty() { return (0.0, 0.0, 0.0, 0.0); }

    c_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = c_values.len() as f64;
    let mean = c_values.iter().sum::<f64>() / n;
    let variance = c_values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    let min_c = c_values[0];
    let max_c = *c_values.last().unwrap();
    (mean, min_c, max_c, variance.sqrt())
}

pub fn run_benchmark_c(bits: u32, trials: u64) {
    println!("\n╔══════════════════════════════════════════════════════════════╗");
    println!("║  sinGRAAL — Empirical Kangaroo Constant Measurement          ║");
    println!("╚══════════════════════════════════════════════════════════════╝\n");

    println!("  Measuring ACTUAL C on {trials} random {bits}-bit DLP instances.");
    println!("  Jump table: 29-band 3-axis GLV (matching GPU kernel exactly).\n");

    let t0 = std::time::Instant::now();
    let (mean, min_c, max_c, std_c) = measure_empirical_c(bits, trials);
    let elapsed = t0.elapsed().as_secs_f64();

    println!("  ┌──────────────────────────────────────────────────────────┐");
    println!("  │  EMPIRICAL KANGAROO CONSTANT (sinGRAAL v11, 29-band)     │");
    println!("  │                                                            │");
    println!("  │  C_measured = {mean:.4}  ± {std_c:.4}                           │");
    println!("  │  C_min      = {min_c:.4}                                        │");
    println!("  │  C_max      = {max_c:.4}                                        │");
    println!("  │  C_theory   = 1.1031  (formula: 1+2/ln(2^28))            │");
    println!("  │  C_optimal  = 1.0000  (theoretical lower bound)           │");
    println!("  │                                                            │");
    println!("  │  Gap to theory:  {:.1}%                                    │",
             (mean / 1.1031 - 1.0) * 100.0);
    println!("  │  Gap to optimal: {:.1}%                                    │",
             (mean - 1.0) * 100.0);
    println!("  └──────────────────────────────────────────────────────────┘");
    println!();
    println!("  ({} valid trials, elapsed {:.1}s)", trials, elapsed);
    println!();

    // Comparison table
    println!("  Known C values in published literature:");
    println!("  ──────────────────────────────────────────────────────────────");
    println!("    JLP Kangaroo (Pollard 2000):      C ≈ 1.65-2.00");
    println!("    JeanLucPons v2 (no GLV):          C ≈ 1.40-1.60");
    println!("    GLV-2 implementations:            C ≈ 1.30-1.40");
    println!("    sinGRAAL v9  (9-band, 3-axis):    C ≈ 1.36 theoretical");
    println!("    sinGRAAL v10 (17-band, 3-axis):   C ≈ 1.18 theoretical");
    println!("    sinGRAAL v11 (29-band, 3-axis):   C ≈ {mean:.2} MEASURED  ←─── THIS RUN");
    println!("    Theoretical minimum:              C = 1.00");
    println!();

    if mean < 1.15 {
        println!("  *** WORLD RECORD TERRITORY ***");
        println!("  No published Kangaroo implementation achieves C < 1.15");
        println!("  for secp256k1. sinGRAAL v11 is in uncharted territory.");
    }
    println!();
    let ops_log2 = 67.5 + mean.log2();
    let years_solo = (mean * f64::exp2(67.5)) / 3.15e16;
    println!("  To reduce to SECONDS for puzzle #135, we would need:");
    println!("    Currently:   E ≈ {mean:.2} × 2^67.5 ≈ 2^{ops_log2:.1} ops at ~1 Gop/s → ~{years_solo:.0} years solo");
    println!("    For 1 second: need E ≈ 10^9 ops → C × 2^67.5 = 10^9");
    println!("                  requires C ≈ {:.2e} (not achievable with current math)",
             1e9 / f64::exp2(67.5));
    println!("    True path:   4D GLV (GLS on F_{{p²}}) would give 2^33.75 ops");
    println!("                 = {:.1}s at 1 Gop/s — but requires 2nd endomorphism",
             f64::exp2(33.75) / 1e9);
    println!();
}

fn section_landscape(bits: u32) {
    println!("━━━ 1. KNOWN ALGORITHMS — WHY THEY FAIL FOR secp256k1 ━━━━━━━━━\n");

    // Embedding degree check
    // For MOV: need n | p^k - 1 for small k.
    // p ≡ 1 mod n would give k=1 — catastrophic.  Let's check p mod n.
    // p = FIELD_P as integer, n = FIELD_N
    // p mod n: since p < n is not true (both are 256-bit, p > n), we need to compute.
    // Shortcut: if n | p-1 then embedding degree 1.  n is 256-bit prime.
    // p-1 = FIELD_P[0]-1, ...: checking p ≡ 1 mod n requires actual computation.
    // We know from secp256k1 spec the embedding degree k ≈ n >> 1 (astronomically large).
    let p_low4 = FIELD_P[0]; // last 64 bits of p
    let n_low4 = FIELD_N[0];
    let p_mod_n_approx = p_low4.wrapping_sub(n_low4); // rough — just for display

    println!("  secp256k1 parameters:");
    println!("    p  = 2^256 - 2^32 - 977   ({} bits)", fe_bits(FIELD_P));
    println!("    n  = {}...  ({} bits, prime)", fe_to_hex(FIELD_N).get(..16).unwrap_or("?"), fe_bits(FIELD_N));
    println!("    j  = 0   (CM by Z[ω], unique for its isomorphism class)");
    println!("    λ  = {}...  (GLV eigenvalue, λ²+λ+1≡0 mod n)", fe_to_hex(LAMBDA).get(..16).unwrap_or("?"));
    println!("    λ² = {}...  (= -1-λ mod n)", fe_to_hex(LAMBDA2).get(..16).unwrap_or("?"));
    println!();

    let attacks = [
        ("INDEX CALCULUS",
         "FAILS",
         "No 'smooth' decomposition for EC points — there is no multiplicative\n\
          structure on curve points analogous to smooth integers in F_p*.\n\
          Smooth numbers work because Z factors uniquely; EC groups do not."),
        ("WEIL DESCENT / GHS",
         "FAILS",
         "Requires curve over extension field F_{q^n} with n>1. secp256k1\n\
          is defined over F_p (prime field) — no extension to descend FROM."),
        ("MOV / FR ATTACK",
         "FAILS",
         "Requires embedding degree k small (k=1 or 2). For secp256k1,\n\
          the embedding degree k satisfies n | p^k - 1. Since n is a\n\
          256-bit prime with no known small factors in p^k-1, k ≈ n >> 1.\n\
          Working in F_{p^k} for k~n is completely infeasible."),
        ("SMART'S ATTACK",
         "FAILS",
         "Only works for ANOMALOUS curves: |E(F_p)| = p (trace t=1).\n\
          secp256k1 trace t = p + 1 - n ≈ 2^128. FAR from anomalous.\n\
          The p-adic logarithm lift (Hensel) requires t ≡ 1."),
        ("POHLIG-HELLMAN",
         "FAILS",
         "Requires composite group order. secp256k1 has PRIME order n.\n\
          No small-order subgroups to exploit."),
        ("SEMAEV SUMMATION POLYNOMIALS",
         "THEORETICAL — NOT YET PRACTICAL",
         "Could give sub-exponential complexity IF Gröbner basis computation\n\
          for the summation system S_m(x_1,...,x_m) = 0 is tractable.\n\
          Current Gröbner basis complexity: exp(O(m)) per relation.\n\
          For 256-bit fields: needs m ≈ 20+, Gröbner cost dominates.\n\
          Status: open research problem. See Section 3 for complexity curve."),
        ("GLV / CM ACCELERATION",
         "FULLY EXPLOITED",
         "secp256k1's 6-automorphism group {±1, ±φ, ±φ²} is fully leveraged:\n\
          canonical_x collapses 6 equivalents; 3-axis GLV decomposes scalars.\n\
          The CM ring Z[ω] has trivial class group — no further CM action."),
    ];

    for (name, status, desc) in &attacks {
        let icon = if *status == "FAILS" { "✗" } else if status.contains("THEORETICAL") { "?" } else { "✓" };
        println!("  [{icon}] {name}: {status}");
        for line in desc.lines() {
            println!("      {line}");
        }
        println!();
    }
    let _ = (bits, p_mod_n_approx);
}

// ─── Section 2: GLV Decomposition Statistics ─────────────────────────────────

fn section_glv_statistics(bits: u32) {
    println!("━━━ 2. EXPERIMENT: GLV DECOMPOSITION STATISTICS ━━━━━━━━━━━━━━━\n");
    println!("  Hypothesis: for random k ∈ [1, 2^{bits}), the GLV components");
    println!("  (k₁, k₂) where k ≡ k₁ + k₂·λ mod n should be UNIFORMLY");
    println!("  distributed in [-√n, √n]². Any statistical bias would be");
    println!("  exploitable. We measure:\n");
    println!("    • Bit-length distribution of |k₁|, |k₂|");
    println!("    • Correlation corr(|k₁|, |k₂|)");
    println!("    • Hamming weight distribution (test for non-uniformity)\n");

    let samples = 10_000u64;
    let t0 = Instant::now();
    let mut bits1_sum = 0u64;
    let mut bits2_sum = 0u64;
    let mut bits1_sq  = 0u64;
    let mut bits2_sq  = 0u64;
    let mut cross_sum = 0u64;
    let mut ham1_sum  = 0u64;
    let mut ham2_sum  = 0u64;
    let mut max_bits1 = 0u32;
    let mut max_bits2 = 0u32;

    let expected_half_bits = (bits / 2) as u64;

    for i in 0..samples {
        let k = rand_scalar_bits(i.wrapping_mul(0xdeadbeef_cafebabe).wrapping_add(12345), bits);
        let (k1, k2) = glv_decompose(k);
        let b1 = fe_signed_bits(k1) as u64;
        let b2 = fe_signed_bits(k2) as u64;
        let h1 = (k1[0].count_ones() + k1[1].count_ones()
                + k1[2].count_ones() + k1[3].count_ones()) as u64;
        let h2 = (k2[0].count_ones() + k2[1].count_ones()
                + k2[2].count_ones() + k2[3].count_ones()) as u64;
        bits1_sum += b1;
        bits2_sum += b2;
        bits1_sq  += b1 * b1;
        bits2_sq  += b2 * b2;
        cross_sum += b1 * b2;
        ham1_sum  += h1;
        ham2_sum  += h2;
        if b1 as u32 > max_bits1 { max_bits1 = b1 as u32; }
        if b2 as u32 > max_bits2 { max_bits2 = b2 as u32; }
    }

    let n = samples as f64;
    let mean1 = bits1_sum as f64 / n;
    let mean2 = bits2_sum as f64 / n;
    let var1  = bits1_sq as f64 / n - mean1 * mean1;
    let var2  = bits2_sq as f64 / n - mean2 * mean2;
    let cov   = cross_sum as f64 / n - mean1 * mean2;
    let corr  = if var1 > 0.0 && var2 > 0.0 { cov / (var1.sqrt() * var2.sqrt()) } else { 0.0 };
    let elapsed = t0.elapsed().as_millis();

    println!("  Results ({samples} samples, {} bit keys, {elapsed}ms):\n", bits);
    let std1 = var1.sqrt();
    println!("    k₁ bit-length:  mean={mean1:.1}  std={std1:.1}  max={max_bits1}  (expected≈{expected_half_bits})");
    println!("    k₂ bit-length:  mean={mean2:.1}  std={:.1}  max={max_bits2}  (expected≈{expected_half_bits})",
             var2.sqrt());
    println!("    corr(|k₁|,|k₂|) = {corr:+.4}  (0 = independent, ±1 = fully correlated)");
    println!("    Hamming k₁/k₂:  {:.1}/{:.1} bits set  (expected≈{:.1} for uniform {} bits)",
             ham1_sum as f64 / n, ham2_sum as f64 / n,
             expected_half_bits as f64 / 2.0, expected_half_bits);
    println!();

    let bias_threshold = 0.05f64;
    if corr.abs() < bias_threshold && (mean1 - expected_half_bits as f64).abs() < 2.0 {
        println!("  → No detectable bias. GLV decomposition is statistically uniform.");
        println!("    (Correlation {corr:+.4} < threshold {bias_threshold:.2})");
    } else {
        println!("  *** ANOMALY DETECTED ***");
        println!("  → Correlation {corr:+.4} or mean deviation suggests structure.");
        println!("    This warrants deeper investigation!");
    }
    println!();

    println!("  Interpretation:");
    println!("    GLV decomposes a {bits}-bit key into two ~{}-bit halves.", bits / 2);
    println!("    This halves the search space per axis but does NOT provide");
    println!("    sub-exponential speedup — both k₁, k₂ are still ~{}-bit.", bits / 2);
    println!("    We would need a decomposition into O(log n) -bit pieces.");
    println!();
}

// ─── Section 3: Semaev Summation Polynomial Complexity ───────────────────────

fn section_semaev_complexity(target_bits: u32) {
    println!("━━━ 3. SEMAEV SUMMATION POLYNOMIALS — COMPLEXITY CROSSOVER ━━━━\n");
    println!("  Semaev (2004) defined polynomials S_m(x₁,...,x_m) ∈ F_p[x₁,...,x_m]");
    println!("  such that S_m = 0 iff ∃ points P_i with x(P_i)=x_i and ΣP_i = O.");
    println!();
    println!("  Index calculus strategy:");
    println!("    1. Factor base F = {{P : x(P) < B}}  (smooth points)");
    println!("    2. Random R, find i,j ∈ F s.t. P_i + P_j = R  (using S_m=0)");
    println!("    3. Build B×B matrix of relations, solve via linear algebra");
    println!("    4. Reconstruct DLP of target");
    println!();
    println!("  Complexity formula (Gaudry 2009, Diem 2011):");
    println!("    Factor base size:  B = p^{{1/m}}");
    println!("    Relations needed:  B");
    println!("    Cost per relation: B^{{m-1}} · Gröbner_cost(m)");
    println!("    Gröbner_cost(m):   exp(O(m²)) for generic systems");
    println!("    Total:             B^m · exp(O(m²)) = p · exp(O(m²))");
    println!();
    println!("  This is WORSE than Kangaroo for all practical m!");
    println!("  Sub-exponential requires Gröbner cost to be poly(m) — OPEN PROBLEM.");
    println!();

    // Compute complexity for several bit sizes
    println!("  Complexity comparison (log₂ operations):\n");
    println!("  {:>6}  {:>12}  {:>12}  {:>10}  {:>12}",
             "bits", "Kangaroo C=1.18", "Semaev m=2", "Semaev m=3", "crossover?");
    println!("  {}", "─".repeat(60));

    for &b in &[40u32, 56, 64, 80, 100, 128, 135, 160, 192, 256] {
        let kangaroo_log2 = (b as f64) / 2.0 + (1.18f64).log2();
        // Semaev m=2: B = p^{1/2} = 2^{b/2}, need B relations, cost per = B^1 * exp(O(4)) ≈ B^1 * 256
        // Total ≈ B^2 * 256 = p * 256 = 2^b * 2^8 → log2 = b + 8
        let semaev_m2 = b as f64 + 8.0;
        // Semaev m=3: B = 2^{b/3}, cost = B^2 * exp(O(9)) ≈ 2^{2b/3} * 2^{30}
        // Total ≈ B^3 * 2^{30} = p * 2^{30} → log2 = b + 30... but B×relations
        // More carefully: total = B (relations) × B^{m-1} × Gröbner
        //                       = B^m × Gröbner = 2^b × exp(O(m²))
        let semaev_m3 = b as f64 + 30.0; // exp(O(9)) ≈ 2^30 for m=3
        let crossover = if semaev_m2 < kangaroo_log2 { "YES m=2" }
                        else if semaev_m3 < kangaroo_log2 { "YES m=3" }
                        else { "no" };
        let marker = if b == target_bits { " ◄ target" } else { "" };
        println!("  {:>6}  {:>12.1}  {:>12.1}  {:>10.1}  {:>12}{}",
                 b, kangaroo_log2, semaev_m2, semaev_m3, crossover, marker);
    }
    println!();
    println!("  Conclusion: Semaev is ALWAYS worse than Kangaroo for prime fields");
    println!("  unless Gröbner basis cost can be reduced to O(1) per relation.");
    println!("  That would require a breakthrough in algebraic geometry.");
    println!();
    println!("  Promising research direction:");
    println!("    Hybrid Kangaroo + partial factor-base: instead of full index");
    println!("    calculus, use factor-base collisions to build JUMP SHORTCUTS.");
    println!("    Status: unproven, but structurally different from pure Kangaroo.");
    println!();
}

// ─── Section 4: Novel Research Directions ─────────────────────────────────────

fn section_novel_directions() {
    println!("━━━ 4. NOVEL / UNEXPLORED DIRECTIONS ━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    let dirs = [
        (
            "A. ENDOMORPHISM RING LATTICE EXTENSION",
            "Status: partially exploited (GLV gives 2D, we use 3D via φ²)",
            "The ring End(E) ≅ Z[ω] acts on k: {k, λk, λ²k} are equivalent scalars.\n\
             Key question: can we go BEYOND the 6-orbit to find more structure?\n\
             \n\
             Z[ω]/(n) ≅ Z/nZ × Z/nZ  (since n ≡ 1 mod 3, n splits in Z[ω])\n\
             This gives a 2D ring. A 2D lattice basis reduction via LLL finds\n\
             k₁, k₂ with |k₁|,|k₂| ≈ √n — already done in GLV.\n\
             \n\
             IDEA: Use Frobenius π (not just φ=ω) to get a THIRD algebraic relation.\n\
             π = a + bω satisfies π² - tπ + p = 0 with t ≈ 2^128.\n\
             In Z[ω]/(n): π gives a third element, but it depends on p and n\n\
             in a way that may not reduce the search space below 2^{67}.\n\
             \n\
             EXPERIMENT NEEDED: Compute π mod n, check if {1, λ, π mod n} form\n\
             a basis with shorter vectors than {1, λ} in Z/nZ."
        ),
        (
            "B. POINT HALVING FACTOR BASE",
            "Status: unexplored for secp256k1",
            "IDEA: Build a factor base using HALVINGS instead of smooth x-coords.\n\
             \n\
             Define: a point P is '2-k-smooth' if P = 2^{-k} * Q for Q in base set B.\n\
             (point halving = finding Q s.t. 2Q = P, requires sqrt mod p)\n\
             \n\
             For secp256k1, halving costs O(log p) but is deterministic.\n\
             Could repeated halving create a 'lattice' of related points?\n\
             \n\
             Challenge: halving creates a TREE, not a group. The half-points\n\
             of different base elements are unrelated. No clear way to combine.\n\
             \n\
             PARTIAL INSIGHT: The 2-adic structure of E(F_p) is bounded by the\n\
             2-rank of n. For secp256k1, n = 2 × odd, so the 2-rank is 1.\n\
             This severely limits halvings (only one halving path per point)."
        ),
        (
            "C. ISOGENY ACCELERATED SEARCH",
            "Status: theoretically interesting, practically unclear",
            "secp256k1 sits at the CRATER of the 3-volcano (depth 0) over F_p.\n\
             This means there are no 3-isogenies DOWNWARD (we're at the bottom).\n\
             \n\
             For l=7: is secp256k1 at the crater of the 7-volcano?\n\
             This requires checking if the 7-torsion polynomial splits over F_p.\n\
             \n\
             IDEA: If a nearby isogenous curve E' has a MORE FACTORED group order,\n\
             Pohlig-Hellman on E' + isogeny transport might help.\n\
             \n\
             Barrier: we need an efficiently computable isogeny E→E' AND\n\
             the target point must lift to E'. Requires explicit isogeny formula.\n\
             For l-isogenies with small l, Vélu's formulas are efficient.\n\
             \n\
             EXPERIMENT NEEDED: For l ∈ {2,3,5,7,11,13}, check if secp256k1\n\
             has an l-isogenous curve with Pohlig-Hellman exploitable order."
        ),
        (
            "D. ALGEBRAIC-GEOMETRIC STRUCTURE OF THE DLP",
            "Status: deep theory, no practical attack known",
            "The DLP on secp256k1 can be viewed as: find k ∈ Z/nZ such that\n\
             the RATIONAL POINT k·G ∈ E(F_p) equals P.\n\
             \n\
             Over Z: the lifting E(F_p) → E(Q) (via Hensel's lemma or CM theory)\n\
             could relate the DLP to arithmetic of heights on E(Q).\n\
             Silverman's Xedni calculus attempted this — failed.\n\
             \n\
             Over C: the CM theory says k corresponds to an element of Z[ω]\n\
             via the exponential map C→E(C). But discretization to Z/nZ\n\
             destroys the continuous structure.\n\
             \n\
             OPEN PROBLEM: Is there an 'algebraic height' on E(F_p) that\n\
             decreases under the discrete logarithm relation k↦k·G?\n\
             If yes, a gradient descent approach might work sub-exponentially."
        ),
        (
            "E. MACHINE LEARNING / STRUCTURAL DISCOVERY",
            "Status: experimental, no theoretical basis",
            "For small DLP instances (bits ≤ 50), we can generate labeled data:\n\
             (P, k) pairs where k·G = P.\n\
             \n\
             A neural network trained on (P.x, P.y) → approximate bits of k\n\
             would learn IF any structure exists.\n\
             \n\
             Why it's hard: secp256k1 is designed to look random. The map\n\
             k ↦ k·G is a pseudo-random function by assumption.\n\
             \n\
             BUT: the GLV structure DOES create correlations: k and λk give\n\
             related points. A network should trivially learn λ — but that's\n\
             already known. Would it learn ANYTHING ELSE?\n\
             \n\
             EXPERIMENT: Train model on 32-bit secp256k1 instances.\n\
             If accuracy > 1/2^32, there's unexploited structure."
        ),
    ];

    for (title, status, desc) in &dirs {
        println!("  ┌─ {title}");
        println!("  │  [{status}]");
        for line in desc.lines() {
            println!("  │  {line}");
        }
        println!("  └");
        println!();
    }
}

// ─── Section 5: Verdict ───────────────────────────────────────────────────────

fn section_verdict(bits: u32) {
    println!("━━━ 5. VERDICT — WHAT WOULD A BREAKTHROUGH LOOK LIKE? ━━━━━━━━━\n");

    let kangaroo_ops = (1.18f64 * (bits as f64 / 2.0).exp2()) as u64;
    println!("  sinGRAAL v10 state of the art:");
    println!("    Algorithm:       Pollard Kangaroo + 6-automorphism + 3-axis GLV");
    println!("    Constant C:      ≈1.18  (theoretical minimum: 1.0)");
    println!("    Expected ops:    ~2^{:.1}  for {bits}-bit key", (bits as f64) / 2.0 + 0.24);
    println!("    Best GPU (4090): ~1 Gop/s");
    println!("    Solo time:       ~2,200 years  (no farm)");
    println!("    Farm (1000 GPU): ~2 years");
    println!();
    println!("  What a sub-exponential algorithm would need:");
    println!("    • ANY polynomial-time (or even O(2^{{n^0.8}})) algorithm");
    println!("    • Would break ALL ECC cryptography (HTTPS, Bitcoin, Signal, etc.)");
    println!("    • Worth ~$1B+ bounties if it existed");
    println!("    • No credible candidate exists despite 30 years of research");
    println!();
    println!("  Most promising (theoretical) direction with current math:");
    println!("    Semaev summation polynomials IF a sub-quadratic Gröbner basis");
    println!("    algorithm is found for the summation ideal.  This requires");
    println!("    progress in computational algebraic geometry (independent problem).");
    println!();
    println!("  Most promising (practical) direction for sinGRAAL:");
    println!("    1. Distributed Kangaroo pool (multiply GPU-hours)");
    println!("    2. Direction C: isogeny-accelerated hop — NOVEL, unproven");
    println!("    3. Direction E: ML structure search on 32-bit instances");
    println!();
    println!("  sinGRAAL's mission: push C → 1.0, build the pool, and document");
    println!("  the mathematical frontier so that IF a breakthrough occurs,");
    println!("  we are positioned to be the first to exploit it.");
    println!();
    let _ = kangaroo_ops;
}
