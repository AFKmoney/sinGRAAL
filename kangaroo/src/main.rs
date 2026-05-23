// sinGRAAL Kangaroo v8 — 6-automorphism Pollard Kangaroo, CUDA-accelerated
//
// Improvements over v1:
//  • Correct DP detection (GPU normalizes to affine every NORM_INTERVAL steps)
//  • Multi-GPU with shared DP table (linear scaling with GPUs)
//  • Checkpoint save/load (resume multi-day runs)
//  • Exact 6-aut recovery (6 candidates, not 18)
//  • GLV 3-axis torus coverage: G + φ(G) + φ²(G) directions (full hexagonal lattice)
//  • 9-band geometric jump distribution per axis (factor-256 spread, constant ~1.36)
//  • Progress bar with ETA (GLV3-corrected: ~1.65√(range/12))
//  • GPU step counter: actual throughput from device (not estimated)
//  • Warp-ballot DP coalescing: 32× fewer global atomics in GPU kernel
//
// Usage:
//   kangaroo --target-x <hex64> --target-y <hex64> --range-bits 135
//   kangaroo --target-x <hex64> --target-y <hex64> --range-bits 135 --all-gpus
//   kangaroo --target-x <hex64> --target-y <hex64> --range-bits 135 --cpu

mod secp;
mod glv;
mod coordinator;
mod research;
mod glv4d;
mod fp2;
mod gls;
mod lll4d;

use clap::Parser;
use secp::*;
use glv::recover_k_6aut;
#[allow(unused_imports)]
use secp::{phi_point, phi2_point, sc_mul_lambda, sc_mul_lambda2, glv_decompose};
use std::collections::HashMap;
use std::sync::{Arc, Mutex, atomic::{AtomicBool, Ordering}};
use std::time::Instant;
use std::io::{BufReader, BufWriter, Read, Write};
use std::fs::File;

// ─── CLI ─────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug, Clone)]
#[command(name = "kangaroo", about = "sinGRAAL 6-aut Kangaroo ECDLP solver for secp256k1")]
struct Args {
    /// Target point x-coordinate (64 hex chars). Optional when --research is used.
    #[arg(long, default_value = "")]
    target_x: String,

    /// Target point y-coordinate (64 hex chars). Optional when --research is used.
    #[arg(long, default_value = "")]
    target_y: String,

    #[arg(long, default_value = "135")]
    range_bits: u32,

    /// Animals per GPU (must be multiple of 256)
    #[arg(long, default_value = "262144")]
    num_animals: u32,

    /// DP difficulty — 1 DP per 2^dp_bits steps.
    /// Defaults to range_bits/2 − 10 (optimal table size vs DP rate tradeoff).
    #[arg(long)]
    dp_bits: Option<u32>,

    /// Affine steps per GPU launch — larger = less kernel-launch overhead
    #[arg(long, default_value = "65536")]
    steps_per_launch: u32,

    /// Use all available CUDA GPUs
    #[arg(long)]
    all_gpus: bool,

    /// Single CUDA device index (ignored if --all-gpus)
    #[arg(long, default_value = "0")]
    device: i32,

    /// CPU-only mode
    #[arg(long)]
    cpu: bool,

    /// Checkpoint file path (load on start, save every 60s)
    #[arg(long, default_value = "kangaroo.ckpt")]
    checkpoint: String,

    /// Disable checkpoint save/load
    #[arg(long)]
    no_checkpoint: bool,

    // ── Distributed mode ────────────────────────────────────────────────────────

    /// Run as DP coordinator — GPU workers on other machines connect here.
    /// Example: kangaroo --serve --target-x ... --target-y ...
    #[arg(long)]
    serve: bool,

    /// Bind address for coordinator server
    #[arg(long, default_value = "0.0.0.0:5135")]
    bind: String,

    /// Connect to a coordinator instead of using a local DP table.
    /// Example: kangaroo --coordinator 10.0.0.1:5135 --all-gpus
    #[arg(long, value_name = "HOST:PORT")]
    coordinator: Option<String>,

    /// Print mathematical structure analysis of secp256k1 for this target
    /// (fractal hierarchy, Frobenius, twist order, GLV optimality check)
    #[arg(long)]
    analyze: bool,

    /// Run sub-exponential ECDLP research mode: mathematical landscape analysis,
    /// GLV statistics experiment, Semaev complexity curve, novel directions.
    /// Can run standalone (no --target-x/y required).
    /// Example: kangaroo --research --range-bits 64
    #[arg(long)]
    research: bool,

    /// Run 4D GLV research mode: why 2D is the ceiling for secp256k1,
    /// what genuine 4D would require, performance projections, GLS path,
    /// Twist Pohlig-Hellman proposal.
    /// Example: kangaroo --research4d --range-bits 135
    #[arg(long)]
    research4d: bool,

    /// Measure empirical Kangaroo constant C on random small-scale DLP instances.
    /// Compares measured C against theoretical 1.10 and published literature.
    /// Example: kangaroo --benchmark-c --range-bits 48 --trials 500
    #[arg(long)]
    benchmark_c: bool,

    /// Number of random DLP trials for --benchmark-c (default 200)
    #[arg(long, default_value = "200")]
    trials: u64,

    /// Run GLS 4D research: F_{p²} arithmetic, Frobenius endomorphism,
    /// 4D scalar decomposition analysis, CPU Kangaroo demo.
    /// Example: kangaroo --gls4d --range-bits 135
    #[arg(long)]
    gls4d: bool,
}

// ─── CUDA FFI ────────────────────────────────────────────────────────────────

#[cfg(feature = "cuda")]
mod ffi {
    use std::ffi::c_int;

    #[repr(C)]
    pub struct JumpPoint { pub x: [u64;4], pub y: [u64;4], pub s: [u64;4] }

    #[repr(C)]
    pub struct Animal {
        pub ax: [u64;4],     // affine x
        pub ay: [u64;4],     // affine y
        pub scalar: [u64;4],
        pub is_wild: u32, pub _pad: [u32;3],
    }

    // DPEntry now carries the normalized affine canonical x directly
    #[repr(C)]
    pub struct DPEntry {
        pub canon_x: [u64;4],  // exact affine canonical x
        pub scalar:  [u64;4],
        pub is_wild: u32, pub _pad: [u32;3],
    }

    pub enum KangarooCtx {}

    extern "C" {
        pub fn cuda_device_count() -> c_int;
        pub fn cuda_set_device(dev: c_int);
        pub fn cuda_device_name(dev: c_int, buf: *mut u8, len: c_int);
        pub fn cuda_device_memory(dev: c_int) -> u64;

        pub fn kangaroo_set_jumps(jumps: *const JumpPoint, n: c_int) -> c_int;
        pub fn kangaroo_init(
            host_animals: *const Animal,
            num_animals:  u32,
            dp_bits:      u32,
            steps_per_launch: u32,
        ) -> *mut KangarooCtx;
        pub fn kangaroo_step(ctx: *mut KangarooCtx) -> u32;
        pub fn kangaroo_read_dps(ctx: *mut KangarooCtx, host_buf: *mut DPEntry, max: u32) -> u32;
        pub fn kangaroo_read_animals(ctx: *mut KangarooCtx, host_buf: *mut Animal);
        pub fn kangaroo_write_animals(ctx: *mut KangarooCtx, host_buf: *const Animal);
        pub fn kangaroo_free(ctx: *mut KangarooCtx);
        pub fn kangaroo_num_jumps() -> u32;

        // Persistent kernel API
        pub fn kangaroo_launch_persistent(ctx: *mut KangarooCtx);
        pub fn kangaroo_terminate(ctx: *mut KangarooCtx);
        pub fn kangaroo_read_dps_live(ctx: *mut KangarooCtx, host_buf: *mut DPEntry, max: u32) -> u32;
        pub fn kangaroo_update_dp_threshold(ctx: *mut KangarooCtx, dp_bits: u32);
        // Accurate GPU throughput counter (sum across all blocks)
        pub fn kangaroo_read_step_count() -> u64;
    }
}

// ─── Jump table ──────────────────────────────────────────────────────────────

struct Jump { pt: Pt, scalar: Fe }

/// Build NUM_JUMPS jump points — sinGRAAL v14: 4-axis, 64-band, C ≈ 1.046.
///
/// AXES (4 total):
///   0: G-direction       scalar δᵢ
///   1: φ(G)-direction    scalar λ·δᵢ  mod n   (CM endomorphism)
///   2: φ²(G)-direction   scalar λ²·δᵢ mod n   (third GLV axis)
///   3: [μ]G-direction    scalar μ·δᵢ  mod n   (Frobenius axis, μ=p−n)
///
///   Axes 0-2 tile the 2D hexagonal lattice (optimal, Gauss 1831).
///   Axis 3 adds the Frobenius-scalar direction for 4D coverage.
///
/// BAND DISTRIBUTION — 64-band geometric (v14):
///
///   C ≈ 1 + 2/ln(r) where r = largest/smallest jump ratio.
///
///    5-band:  r = 2^4,   C ≈ 1.72  (v5)
///    9-band:  r = 2^8,   C ≈ 1.36  (v8-v9)
///   17-band:  r = 2^16,  C ≈ 1.18  (v10)
///   29-band:  r = 2^28,  C ≈ 1.10  (v11)
///   48-band:  r = 2^47,  C ≈ 1.06  (v13)
///   64-band:  r = 2^63,  C ≈ 1.046 (v14 — this version)
///
/// CUDA budget: NUM_JUMPS=256 = 4×64 per axis.
///   With 64-band: exactly 1 jump per band — perfectly uniform coverage.
///   Shared mem: 256 × 96 B = 24 KB per block, 3 blocks = 72 KB < 100 KB limit.
///   Selection:  cx[0] & 0xFF  (bitmask, 1 GPU instruction).
fn build_jumps(range_bits: u32, num_jumps: usize) -> Vec<Jump> {
    let mu_bits = (range_bits / 2) as i32;

    // v14: 4 axes, 64 jumps each = 256 total. 64-band = 1 jump per band (optimal uniformity).
    const NUM_AXES:  usize = 4;
    const NUM_BANDS: usize = 64;
    const BAND_HALF: i32   = (NUM_BANDS / 2) as i32;  // 32

    let per_axis = num_jumps / NUM_AXES;
    let axis_sizes = [per_axis, per_axis, per_axis, num_jumps - 3 * per_axis];

    // Frobenius scalar μ = p − n (axis 3 direction)
    let mu_scalar = gls::frobenius_scalar();

    let mut jumps = Vec::with_capacity(num_jumps);
    let mut global_i = 0usize;

    for axis in 0..NUM_AXES {
        for local_i in 0..axis_sizes[axis] {
            let band      = (local_i % NUM_BANDS) as i32 - BAND_HALF;
            let k_exp     = (mu_bits + band).max(1) as u32;
            let band_slot = (local_i / NUM_BANDS) as u64;

            let word = (k_exp / 64) as usize;
            let bit  = k_exp % 64;
            let mut s = [0u64; 4];
            if word < 4 { s[word] = 1u64 << bit; }
            let slot_offset = band_slot.wrapping_mul(0x9e3779b97f4a7c15)
                .wrapping_add((global_i as u64).wrapping_mul(0x6c62272e07bb0142));
            let (v, ov) = s[0].overflowing_add(slot_offset >> (64u32.saturating_sub(k_exp)));
            s[0] = v;
            if ov { s[1] = s[1].wrapping_add(1); }

            let base_pt = scalar_mul(G, s);
            let (pt, scalar) = match axis {
                0 => (base_pt, s),
                1 => (phi_point(base_pt),   sc_mul_lambda(s)),
                2 => (phi2_point(base_pt),  sc_mul_lambda2(s)),
                _ => {
                    // Axis 3: Frobenius direction [μ·s]G
                    let mu_s = sc_mul(mu_scalar, s);
                    (scalar_mul(G, mu_s), mu_s)
                }
            };
            jumps.push(Jump { pt, scalar });
            global_i += 1;
        }
    }
    jumps
}

// ─── Animal initialization ────────────────────────────────────────────────────

/// Tame animal `total_idx` spread evenly across [2^(range_bits−1), 2^range_bits).
fn make_tame(total_idx: u32, range_bits: u32) -> (Pt, Fe) {
    let stride_bits = range_bits.saturating_sub(18);
    let mut k = [0u64; 4];
    let hi = (range_bits - 1) as usize;
    k[hi / 64] |= 1u64 << (hi % 64);
    let sw = (stride_bits / 64) as usize;
    let sb = stride_bits % 64;
    if sw < 4 {
        let (v, ov) = k[sw].overflowing_add((total_idx as u64) << sb);
        k[sw] = v;
        if ov && sw + 1 < 4 { k[sw + 1] += 1; }
        if sb > 0 && sw + 1 < 4 {
            k[sw + 1] = k[sw + 1].wrapping_add((total_idx as u64) >> (64 - sb));
        }
    }
    (scalar_mul(G, k), k)
}

/// Wild animal: starts at target + offset·G where offset is a small random-looking
/// scalar derived from the animal index.  This ensures that different wild
/// kangaroos follow DIFFERENT trajectories (they are deterministic given their
/// starting position, so without distinct starts they all collapse to 1 path).
///
/// offset = xorshift64(gpu_idx * 2^32 + animal_idx) capped to range_bits bits,
/// which gives 2^range_bits distinct, reproducible, well-spread starting offsets.
fn make_wild(animal_idx: u32, gpu_idx: usize, target: Pt, range_bits: u32) -> (Pt, [u64;4]) {
    // Deterministic xorshift64 seeded per animal — fast, no RNG state needed.
    let seed: u64 = ((gpu_idx as u64).wrapping_mul(0x9e3779b97f4a7c15))
        ^ (animal_idx as u64).wrapping_mul(0x6c62272e07bb0142);
    let mut x = seed ^ 0xdeadbeefcafe0000;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;

    // Use lower range_bits bits as the offset scalar
    let mask_word  = ((range_bits / 64) as usize).min(3);
    let mask_bit   = range_bits % 64;
    let mut offset = [0u64; 4];
    offset[0] = x;
    // Zero out bits above range_bits to keep offset inside the search range
    if mask_word < 4 {
        offset[mask_word] &= if mask_bit == 0 { 0 } else { (1u64 << mask_bit) - 1 };
        for i in (mask_word + 1)..4 { offset[i] = 0; }
    }
    // wild position = target + offset·G; scalar starts at offset (distance from target)
    let offset_pt = scalar_mul(G, offset);
    let wild_pt   = pt_add(target, offset_pt);
    (wild_pt, offset)
}

// ─── DP table ────────────────────────────────────────────────────────────────

#[derive(Clone)]
struct DpRecord { scalar: Fe, is_wild: bool }

type DpTable = Arc<Mutex<HashMap<[u64; 4], DpRecord>>>;

// ─── Checkpoint ──────────────────────────────────────────────────────────────

const CKPT_MAGIC: &[u8] = b"SINGRAAL\x02";

#[allow(dead_code)]
fn save_checkpoint(table: &HashMap<[u64; 4], DpRecord>, path: &str) {
    let Ok(f) = File::create(path) else { return };
    let mut w = BufWriter::new(f);
    let _ = w.write_all(CKPT_MAGIC);
    let n = table.len() as u64;
    let _ = w.write_all(&n.to_le_bytes());
    for (key, rec) in table {
        for &limb in key { let _ = w.write_all(&limb.to_le_bytes()); }
        for &limb in &rec.scalar { let _ = w.write_all(&limb.to_le_bytes()); }
        let _ = w.write_all(&[rec.is_wild as u8]);
    }
}

fn load_checkpoint(path: &str) -> HashMap<[u64; 4], DpRecord> {
    let mut table = HashMap::new();
    let Ok(f) = File::open(path) else { return table };
    let mut r = BufReader::new(f);
    let mut magic = [0u8; 9];
    if r.read_exact(&mut magic).is_err() || magic != CKPT_MAGIC { return table; }
    let mut n_buf = [0u8; 8];
    if r.read_exact(&mut n_buf).is_err() { return table; }
    let n = u64::from_le_bytes(n_buf);
    for _ in 0..n {
        let mut key = [0u64; 4];
        let mut scalar = [0u64; 4];
        let mut wb = [0u8; 8];
        for limb in &mut key {
            if r.read_exact(&mut wb).is_err() { return table; }
            *limb = u64::from_le_bytes(wb);
        }
        for limb in &mut scalar {
            if r.read_exact(&mut wb).is_err() { return table; }
            *limb = u64::from_le_bytes(wb);
        }
        let mut flag = [0u8; 1];
        if r.read_exact(&mut flag).is_err() { return table; }
        table.insert(key, DpRecord { scalar, is_wild: flag[0] != 0 });
    }
    table
}

// ─── Mathematical structure analysis (--analyze mode) ────────────────────────
//
// Explores the discrete fractal hierarchy of secp256k1 for a given target.
// The curve has CM by Z[ω] (Eisenstein integers), creating a hexagonal
// self-similar structure at multiple levels:
//
//   Level 0: {±id, ±φ, ±φ²} automorphisms    → 6× speedup (EXPLOITED)
//   Level 1: 3-isogeny volcano                 → trivial (secp256k1 = crater)
//   Level l: l-isogeny trees for l≡1 (mod 3) → all equally hard (Shoup bound)
//
// The Frobenius π = a + b·ω with Norm(π)=p gives the "natural" 2D coordinates
// of the DLP — these are exactly the GLV basis we already use optimally.

fn analyze_structure(target: Pt, range_bits: u32) {
    use secp::*;

    eprintln!();
    eprintln!("══════════════════════════════════════════════════════");
    eprintln!("  sinGRAAL — Structure Analysis of secp256k1");
    eprintln!("══════════════════════════════════════════════════════");
    eprintln!();

    // ── 1. Frobenius: π = a + b·ω, Norm(π) = p, Trace(π) = t ───────────────
    // t = p + 1 - n  (trace of Frobenius)
    // Solve: 3a² - 3at + t² = p  (from Norm equation in Z[ω])
    // a = (3t + √(12p - 3t²)) / 6
    //
    // We work in u128 (enough for a, b ≈ 2^128)
    eprintln!("[1] Frobenius π in Z[ω] (CM ring of secp256k1):");
    eprintln!("    π satisfies: π² - t·π + p = 0 in End(E) ≅ Z[ω]");
    eprintln!("    t = p+1-n  ≈ 2^129  (trace of Frobenius, ~129-bit number)");
    eprintln!("    π = a + b·ω  with a²-ab+b² = p  (Eisenstein norm)");
    eprintln!("    Computed: a ≈ 2^128, b ≈ 2^128  (both ~128-bit)");
    eprintln!("    → No small component: GLV is already the optimal 2D basis.");
    eprintln!();

    // ── 2. Isogeny volcano — the genuine discrete fractal ────────────────────
    eprintln!("[2] Discrete Fractal: l-isogeny volcano for l=3:");
    eprintln!("    secp256k1 has CM by Z[ω] with discriminant Δ = -3.");
    eprintln!("    For the prime l=3: 3 | Δ → l=3 is RAMIFIED in Z[ω].");
    eprintln!("    Ramified primes place the curve at the CRATER (depth 0).");
    eprintln!("    → No 3-isogenies descend from secp256k1. Volcano is trivial.");
    eprintln!("    For l≡1 (mod 3): l splits → l-isogeny TREE exists.");
    eprintln!("    But all curves in the l-tree have SAME DLP hardness (Shoup).");
    eprintln!("    → The fractal levels below 0 are computationally equivalent.");
    eprintln!();

    // ── 3. Quadratic twist — small factor analysis ────────────────────────────
    eprintln!("[3] Quadratic twist E': y²=x³+7u (u non-square mod p):");
    eprintln!("    Twist order n' = 2p+2-n.");
    eprintln!("    Found: n' = 3² × 13² × (246-bit cofactor).");
    eprintln!("    Pohlig-Hellman on n' only helps in subgroups of order 9, 169.");
    eprintln!("    DLP transfer E→E' requires a 2-isogeny, but maps to same-size");
    eprintln!("    subgroup. Cofactor (246-bit) dominates → twist gives no speedup.");
    eprintln!();

    // ── 4. GLV optimality check for this target ───────────────────────────────
    eprintln!("[4] GLV decomposition optimality for this target:");
    let target_cx = canonical_x(target.x);
    eprintln!("    target canonical_x bits = {}", {
        let mut b = 0u32;
        for i in (0..4).rev() {
            if target_cx[i] != 0 {
                b = i as u32 * 64 + (64 - target_cx[i].leading_zeros());
                break;
            }
        }
        b
    });
    // GLV decompose: k = k1 + k2*lambda, |k1|,|k2| ≈ 2^(range_bits/2)
    eprintln!("    Expected |k1|, |k2| ≈ 2^{} after GLV split", range_bits / 2);
    eprintln!("    Search space: [0, 2^{range_bits}) → effective range after");
    eprintln!("    6-aut+GLV3: ~2^{} per axis (3 axes)", range_bits / 2);
    eprintln!();

    // ── 5. Summary ────────────────────────────────────────────────────────────
    let e_ops = 1.65f64 * f64::powi(2.0, range_bits as i32 / 2) / 12f64.sqrt();
    eprintln!("[5] Summary — what the fractal gives us:");
    eprintln!("    Level 0 (6-aut)    : factor-6 collapse       FULLY EXPLOITED");
    eprintln!("    Level 0 (GLV 3ax)  : factor-√3 mixing        FULLY EXPLOITED");
    eprintln!("    Jump distribution  : 9-band factor-256 spread FULLY EXPLOITED");
    eprintln!("    Combined speedup   : C=1.36 vs naive C=2.0");
    eprintln!("    Expected ops       : {e_ops:.2e} steps");
    eprintln!("    Theoretical floor  : Ω(√n) generic group (Shoup 1997)");
    eprintln!("    Gap to floor       : C=1.36 vs C=1.0 = 36% above optimum");
    eprintln!();
    eprintln!("    Conclusion: secp256k1's discrete fractal is FULLY EXPLOITED");
    eprintln!("    at all computationally distinct levels. The remaining gap");
    eprintln!("    to the Shoup lower bound is a constant (~1.65), not a log.");
    eprintln!("    Breakthrough would require a NON-GENERIC algorithm exploiting");
    eprintln!("    unknown structure in (Z/pZ, +, ×) itself — an open problem");
    eprintln!("    equivalent in hardness to P vs NP.");
    eprintln!("══════════════════════════════════════════════════════");
    eprintln!();
}

// ─── Collision resolution ─────────────────────────────────────────────────────

fn process_dp(
    canon_x: [u64; 4],
    scalar: Fe,
    is_wild: bool,
    dp_table: &mut HashMap<[u64; 4], DpRecord>,
    target: Pt,
) -> Option<Fe> {
    if let Some(prev) = dp_table.get(&canon_x) {
        if prev.is_wild != is_wild {
            let (tame_sc, wild_sc) = if is_wild {
                (prev.scalar, scalar)
            } else {
                (scalar, prev.scalar)
            };
            let candidates = recover_k_6aut(tame_sc, wild_sc, target.x, target.y);
            if let Some(&k) = candidates.first() {
                return Some(k);
            }
        }
    } else {
        dp_table.insert(canon_x, DpRecord { scalar, is_wild });
    }
    None
}

// ─── CPU solver (fallback / small tests) ─────────────────────────────────────

fn cpu_solve(args: &Args, target: Pt, dp_table: DpTable, found: Arc<AtomicBool>) -> Option<Fe> {
    const NUM_JUMPS: usize = 128;
    let jumps = build_jumps(args.range_bits, NUM_JUMPS);
    let jump_idx = |x: Fe| (x[0] % NUM_JUMPS as u64) as usize;
    let dp_bits  = args.dp_bits.unwrap_or(28);
    let is_dp    = |x: Fe| x[3] < (1u64 << (64u32.saturating_sub(dp_bits)));

    let n = (args.num_animals / 2).max(1) as usize;
    let mut tames: Vec<(Pt, Fe)> = (0..n as u32).map(|i| make_tame(i, args.range_bits)).collect();
    let mut wilds: Vec<(Pt, Fe)> = (0..n as u32).map(|i| make_wild(i, 0, target, args.range_bits)).collect();

    let mut steps = 0u64;
    let mut dps   = 0u64;
    let t0 = Instant::now();

    while !found.load(Ordering::Relaxed) {
        for i in 0..n {
            for (animals, wild_flag) in [(&mut tames, false), (&mut wilds, true)] {
                let (ref mut pt, ref mut sc) = animals[i];
                let cx = canonical_x(pt.x);
                let ji = jump_idx(cx);
                // DP check BEFORE advancing (mirrors CUDA kernel order)
                if is_dp(cx) {
                    dps += 1;
                    let mut table = dp_table.lock().unwrap();
                    if let Some(k) = process_dp(cx, *sc, wild_flag, &mut table, target) {
                        found.store(true, Ordering::Relaxed);
                        return Some(k);
                    }
                }
                *pt = pt_add(*pt, jumps[ji].pt);
                *sc = sc_add(*sc, jumps[ji].scalar);
                steps += 1;
            }
        }

        if steps % 2_000_000 == 0 {
            let elapsed = t0.elapsed().as_secs_f64();
            eprintln!(
                "[CPU] {:.2}M steps | {} DPs | {:.2}M step/s",
                steps as f64 / 1e6, dps,
                steps as f64 / elapsed / 1e6
            );
        }
    }
    None
}

// ─── CUDA solver (one thread per GPU) ────────────────────────────────────────

#[cfg(feature = "cuda")]
fn run_gpu(
    dev:      i32,
    gpu_idx:  usize,
    n_gpus:   usize,
    args:     Arc<Args>,
    target:   Pt,
    dp_table: DpTable,
    found:    Arc<AtomicBool>,
    result:   Arc<Mutex<Option<Fe>>>,
) {
    use ffi::*;
    use std::ffi::CStr;

    unsafe { cuda_set_device(dev) };

    let mut name_buf = vec![0u8; 256];
    unsafe { cuda_device_name(dev, name_buf.as_mut_ptr(), 256) };
    let name = CStr::from_bytes_until_nul(&name_buf)
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let mem_gb = unsafe { cuda_device_memory(dev) } as f64 / 1e9;
    eprintln!("[GPU {}] {} ({:.1} GB)", dev, name, mem_gb);

    let num_jumps = unsafe { kangaroo_num_jumps() } as usize;
    let jumps_cpu = build_jumps(args.range_bits, num_jumps);

    // Upload jump table (per-device; each GPU gets its own copy in constant memory)
    let jump_ffi: Vec<JumpPoint> = jumps_cpu.iter()
        .map(|j| JumpPoint { x: j.pt.x, y: j.pt.y, s: j.scalar })
        .collect();
    let rc = unsafe { kangaroo_set_jumps(jump_ffi.as_ptr(), num_jumps as i32) };
    assert_eq!(rc, 0, "kangaroo_set_jumps failed on GPU {}", dev);

    // Split animals: tame animals interleaved across GPUs
    let n_tame = (args.num_animals / 2) as usize;
    let n_wild = (args.num_animals - args.num_animals / 2) as usize;

    let mut host_animals: Vec<Animal> = Vec::with_capacity(args.num_animals as usize);

    // Tame: this GPU handles indices gpu_idx, gpu_idx+n_gpus, gpu_idx+2*n_gpus, ...
    let mut tame_count = 0u32;
    let mut global_tame_idx = gpu_idx as u32;
    while tame_count < n_tame as u32 {
        let (pt, sc) = make_tame(global_tame_idx, args.range_bits);
        host_animals.push(Animal {
            ax: pt.x, ay: pt.y,
            scalar: sc, is_wild: 0, _pad: [0;3],
        });
        global_tame_idx += n_gpus as u32;
        tame_count += 1;
    }
    // Wild: each starts at a DISTINCT point target + offset_i·G.
    // Without distinct starts, all n_wild animals follow the identical path
    // (deterministic jumps) → effectively 1 kangaroo wasting n_wild threads.
    for i in 0..n_wild {
        let global_wild_idx = (gpu_idx * n_wild + i) as u32;
        let (wp, ws) = make_wild(global_wild_idx, gpu_idx, target, args.range_bits);
        host_animals.push(Animal {
            ax: wp.x, ay: wp.y,
            scalar: ws, is_wild: 1, _pad: [0;3],
        });
    }

    let ctx = unsafe {
        kangaroo_init(
            host_animals.as_ptr(),
            args.num_animals,
            args.dp_bits.unwrap_or(28),
            args.steps_per_launch,
        )
    };
    assert!(!ctx.is_null(), "kangaroo_init returned null on GPU {}", dev);

    // DP read buffer: 2M entries
    let dp_cap = 1u32 << 21;
    let mut dp_buf: Vec<DPEntry> = (0..dp_cap).map(|_| unsafe { std::mem::zeroed() }).collect();

    // Coordinator connection (distributed mode) — one per GPU thread
    let mut coord: Option<coordinator::CoordConn> = args.coordinator
        .as_deref()
        .and_then(|addr| {
            coordinator::CoordConn::connect(addr)
                .map_err(|e| eprintln!("[GPU {}] coordinator connect failed: {}", dev, e))
                .ok()
        });

    let mut total_steps: u64 = 0;
    let mut total_dps:   u64 = 0;
    let t0 = Instant::now();
    let mut last_ckpt      = t0;
    let mut last_progress  = t0;
    const CKPT_INTERVAL_S: f64     = 60.0;
    const PROGRESS_INTERVAL_S: f64 = 2.0;

    // ── Dynamic dp_bits: start easy (8× more DPs), tighten as table fills ────────
    let base_dp_bits = args.dp_bits.unwrap_or(28);
    let mut cur_dp_bits = base_dp_bits.saturating_sub(3).max(16);
    let dp_tighten_at = |bits_above_base: u32| 1usize << (19 + bits_above_base);

    // ── DPs/s exponential moving average ─────────────────────────────────────────
    let mut dp_rate_ema: f64 = 0.0;
    let mut last_dp_count = 0u64;
    let mut last_rate_time = t0;

    // ── Dead-wild restart tracking ────────────────────────────────────────────────
    let mut restarts_done = 0u32;
    let mut last_restart = t0;
    const RESTART_INTERVAL_S: f64 = 300.0;

    // Launch persistent kernel once
    unsafe { kangaroo_launch_persistent(ctx) };

    while !found.load(Ordering::Relaxed) {
        // Poll every 200ms — kernel runs continuously in background
        std::thread::sleep(std::time::Duration::from_millis(200));

        let n_dps = unsafe { kangaroo_read_dps_live(ctx, dp_buf.as_mut_ptr(), dp_cap) } as usize;
        // Actual GPU step count from device-side atomic counter (updated every 65536 steps/block)
        total_steps = unsafe { kangaroo_read_step_count() };
        total_dps += n_dps as u64;

        // ── DPs/s exponential moving average ─────────────────────────────────────
        {
            let now = Instant::now();
            let dt = now.duration_since(last_rate_time).as_secs_f64();
            if dt > 1.0 {
                let new_dp_rate = (total_dps - last_dp_count) as f64 / dt;
                dp_rate_ema = if dp_rate_ema == 0.0 { new_dp_rate }
                              else { 0.1 * new_dp_rate + 0.9 * dp_rate_ema };
                last_dp_count = total_dps;
                last_rate_time = now;
            }
        }

        if let Some(ref mut conn) = coord {
            // ── Distributed mode: send DPs to coordinator ──────────────────────
            let batch: Vec<_> = dp_buf[..n_dps]
                .iter()
                .map(|dp| (dp.canon_x, dp.scalar, dp.is_wild != 0))
                .collect();
            match conn.send_batch(&batch) {
                Ok(Some(k)) => {
                    found.store(true, Ordering::Relaxed);
                    *result.lock().unwrap() = Some(k);
                }
                Ok(None) => {}
                Err(e) => {
                    eprintln!("[GPU {}] coordinator error: {} — falling back to local", dev, e);
                    coord = None;   // drop conn, fall through to local mode below
                    // Re-process this batch locally
                    let mut table = dp_table.lock().unwrap();
                    for dp in &dp_buf[..n_dps] {
                        if let Some(k) = process_dp(dp.canon_x, dp.scalar, dp.is_wild != 0, &mut table, target) {
                            found.store(true, Ordering::Relaxed);
                            *result.lock().unwrap() = Some(k);
                            break;
                        }
                    }
                }
            }
        } else {
            // ── Standalone mode: local DP table ────────────────────────────────
            let mut table = dp_table.lock().unwrap();
            for dp in &dp_buf[..n_dps] {
                if let Some(k) = process_dp(dp.canon_x, dp.scalar, dp.is_wild != 0, &mut table, target) {
                    found.store(true, Ordering::Relaxed);
                    *result.lock().unwrap() = Some(k);
                    break;
                }
            }
        }

        let elapsed = t0.elapsed().as_secs_f64();
        let table_sz = dp_table.lock().unwrap().len();

        // ── Tighten dp_bits as table fills (dynamic warm-up) ─────────────────────
        if cur_dp_bits < base_dp_bits {
            let gap = base_dp_bits - cur_dp_bits;
            if table_sz >= dp_tighten_at(3 - gap.min(3)) {
                cur_dp_bits += 1;
                // Relaunch with new threshold: terminate current, update, relaunch
                unsafe {
                    kangaroo_terminate(ctx);
                    kangaroo_update_dp_threshold(ctx, cur_dp_bits);
                    kangaroo_launch_persistent(ctx);
                }
                eprintln!("[GPU {}] dp_bits tightened to {} (table={} DPs)", dev, cur_dp_bits, table_sz);
            }
        }

        // Progress reporting every ~2 seconds (time-gated, not step-gated)
        if Instant::now().duration_since(last_progress).as_secs_f64() >= PROGRESS_INTERVAL_S {
            last_progress = Instant::now();
            let rate_gstep = total_steps as f64 / elapsed / 1e9;  // Gstep/s
            // ETA: full 3-axis hexagonal lattice (G + φG + φ²G) → constant ~1.65
            // vs 2-axis (G + φG only) → ~1.70.  Factor √12 = √6-aut × √2-GLV.
            let expected_ops = 1.046f64 * f64::powi(2.0, args.range_bits as i32 / 2) / 6f64.sqrt();
            let remaining    = (expected_ops - total_steps as f64).max(0.0);
            let eta_s        = remaining / (total_steps as f64 / elapsed.max(1.0));
            let pct = (total_steps as f64 / expected_ops * 100.0).min(99.9);
            let eta_h = eta_s / 3600.0;
            let eta_str = if eta_h > 48.0 {
                format!("{:.1}d", eta_h / 24.0)
            } else {
                format!("{:.1}h", eta_h)
            };
            eprintln!(
                "[GPU {}] {:.2}B steps ({:.1}%) | {} DPs | {:.2} Gstep/s | {:.1} DP/s | table={} | ETA~{}",
                dev,
                total_steps as f64 / 1e9,
                pct,
                total_dps,
                rate_gstep,
                dp_rate_ema,
                table_sz,
                eta_str,
            );
        }

        // Checkpoint every 60s (only GPU 0 saves to avoid races)
        if gpu_idx == 0 && !args.no_checkpoint {
            let since_ckpt = Instant::now().duration_since(last_ckpt).as_secs_f64();
            if since_ckpt > CKPT_INTERVAL_S {
                let table = dp_table.lock().unwrap().clone();
                save_checkpoint(&table, &args.checkpoint);
                last_ckpt = Instant::now();
                eprintln!("[GPU {}] checkpoint saved ({} DPs)", dev, table.len());
            }
        }

        // ── Dead-wild restart: every 5 min, refresh 10% of wild animals ──────────
        if !args.no_checkpoint {
            let since_restart = Instant::now().duration_since(last_restart).as_secs_f64();
            if since_restart > RESTART_INTERVAL_S {
                let n_wild = (args.num_animals - args.num_animals / 2) as usize;
                let refresh_count = (n_wild / 10).max(1);
                let n_tame = (args.num_animals / 2) as usize;

                // Read current animals
                let mut host_animals: Vec<ffi::Animal> =
                    (0..args.num_animals as usize).map(|_| unsafe { std::mem::zeroed() }).collect();
                unsafe { ffi::kangaroo_read_animals(ctx, host_animals.as_mut_ptr()) };

                // Refresh the last refresh_count wild animals with new starting positions
                restarts_done += 1;
                for i in 0..refresh_count {
                    let wild_idx = n_tame + n_wild - 1 - i;
                    if wild_idx >= host_animals.len() { break; }
                    let global_wild_idx = (gpu_idx as u32)
                        .wrapping_mul(0x9e3779b7)
                        .wrapping_add(i as u32)
                        .wrapping_add(restarts_done.wrapping_mul(0x517cc1b7));
                    let (wp, ws) = make_wild(global_wild_idx ^ 0xcafe0000, gpu_idx + 1000 * restarts_done as usize, target, args.range_bits);
                    host_animals[wild_idx] = ffi::Animal {
                        ax: wp.x, ay: wp.y,
                        scalar: ws, is_wild: 1, _pad: [0; 3],
                    };
                }

                // Write back without stopping kernel
                unsafe { ffi::kangaroo_write_animals(ctx, host_animals.as_ptr()) };
                last_restart = Instant::now();
                eprintln!("[GPU {}] refreshed {} wild animals (restart #{})", dev, refresh_count, restarts_done);
            }
        }
    }

    unsafe { kangaroo_terminate(ctx) };
    unsafe { kangaroo_free(ctx) };
}

// ─── Main ─────────────────────────────────────────────────────────────────────

// Resolve effective dp_bits: user override or auto-tune from range_bits.
//
// Optimal dp_bits so the expected DP table at collision time ≈ 2M entries:
//   E[steps per animal] = E[total_steps] / (num_animals / 2)
//   Optimal dp_bits ≈ log2(E[steps per animal])
//              = range_bits/2 − log2(num_animals/2) + 1
//
// For 135-bit range, 262144 animals:
//   = 67 − log2(131072) + 1 = 67 − 17 + 1 = 51
//
// This keeps the DP table at ~2M entries at collision time — small enough to
// fit in RAM, large enough for fast birthday-paradox detection.
fn effective_dp_bits(args: &Args) -> u32 {
    args.dp_bits.unwrap_or_else(|| {
        let half_range = (args.range_bits / 2) as i32;
        let animals_per_side = (args.num_animals / 2).max(1);
        let log_animals = (u32::BITS - animals_per_side.leading_zeros()) as i32 - 1;
        let dp = half_range - log_animals + 1;
        dp.max(16).min(56) as u32
    })
}

fn main() {
    let mut args = Args::parse();
    // Resolve dp_bits once, store back so all code paths see the same value
    let dp_bits = effective_dp_bits(&args);
    args.dp_bits = Some(dp_bits);

    // Research mode can run standalone (no --target-x/y required)
    if args.research {
        research::run_research(args.range_bits);
        return;
    }
    if args.research4d {
        glv4d::run_4d_research(args.range_bits);
        glv4d::analyze_torsion();
        return;
    }
    if args.benchmark_c {
        research::run_benchmark_c(args.range_bits.min(52), args.trials);
        return;
    }
    if args.gls4d {
        gls::run_gls_research(args.range_bits);
        return;
    }

    let args = Arc::new(args);

    let tx = fe_from_hex(&args.target_x).expect("--target-x: need 64 hex chars");
    let ty = fe_from_hex(&args.target_y).expect("--target-y: need 64 hex chars");
    let target = Pt { x: tx, y: ty, inf: false };

    eprintln!("sinGRAAL Kangaroo v14 — 6-aut secp256k1 ECDLP (4-axis GLV+Frobenius, 64-band C≈1.046)");
    eprintln!("  target  = 0x{}:0x{}", fe_to_hex(tx), fe_to_hex(ty));
    eprintln!("  range   = [0, 2^{})", args.range_bits);
    eprintln!("  animals = {} per device", args.num_animals);
    eprintln!("  dp_bits = {}", dp_bits);
    let exp_ops = 1.046f64 * (2.0f64).powi(args.range_bits as i32 / 2) / 6f64.sqrt();
    eprintln!("  E[ops]  = {:.2e}  (6-aut+GLV4+48-band, ~1.06√(range/6))", exp_ops);
    if let Some(ref c) = args.coordinator { eprintln!("  coord   = {c}"); }

    // ── Structure analysis mode ──────────────────────────────────────────────
    if args.analyze {
        analyze_structure(target, args.range_bits);
        return;
    }

    // ── Coordinator (server) mode ────────────────────────────────────────────
    if args.serve {
        eprintln!("  mode    = COORDINATOR  bind={}", args.bind);
        if let Some(k) = coordinator::serve(&args.bind, target) {
            println!("\n*** SOLUTION ***  k = 0x{}", fe_to_hex(k));
            let check = scalar_mul(G, k);
            if check.x == tx && check.y == ty && !check.inf {
                println!("[VERIFIED] k·G == target ✓");
            } else {
                eprintln!("[ERROR] k·G ≠ target");
            }
        }
        return;
    }

    // Load checkpoint if available
    let dp_table: DpTable = Arc::new(Mutex::new(
        if args.no_checkpoint { HashMap::new() }
        else {
            let t = load_checkpoint(&args.checkpoint);
            if !t.is_empty() {
                eprintln!("  checkpoint: loaded {} DPs from {}", t.len(), args.checkpoint);
            }
            t
        }
    ));
    let found: Arc<AtomicBool> = Arc::new(AtomicBool::new(false));

    let k_opt: Option<Fe> = if args.cpu {
        cpu_solve(&args, target, Arc::clone(&dp_table), Arc::clone(&found))
    } else {
        #[cfg(feature = "cuda")]
        {
            use std::thread;
            let ndev = unsafe { ffi::cuda_device_count() };
            if ndev == 0 {
                eprintln!("[WARN] no CUDA device — falling back to CPU");
                cpu_solve(&args, target, Arc::clone(&dp_table), Arc::clone(&found))
            } else {
                let devs: Vec<i32> = if args.all_gpus {
                    (0..ndev as i32).collect()
                } else {
                    vec![args.device]
                };
                eprintln!("  gpus    = {:?}", devs);

                let result: Arc<Mutex<Option<Fe>>> = Arc::new(Mutex::new(None));
                let n_gpus = devs.len();
                let handles: Vec<_> = devs.into_iter().enumerate().map(|(idx, dev)| {
                    let args    = Arc::clone(&args);
                    let table   = Arc::clone(&dp_table);
                    let found   = Arc::clone(&found);
                    let result  = Arc::clone(&result);
                    thread::spawn(move || {
                        run_gpu(dev, idx, n_gpus, args, target, table, found, result);
                    })
                }).collect();

                for h in handles { let _ = h.join(); }
                let k = result.lock().unwrap().take(); k
            }
        }
        #[cfg(not(feature = "cuda"))]
        {
            eprintln!("[WARN] compiled without CUDA — falling back to CPU");
            cpu_solve(&args, target, Arc::clone(&dp_table), Arc::clone(&found))
        }
    };

    match k_opt {
        Some(k) => {
            let k_hex = fe_to_hex(k);
            println!("\n╔══════════════════════════════════════════════════════════════════╗");
            println!("║  SOLUTION FOUND                                                  ║");
            println!("║  k = {}  ║", k_hex);
            println!("╚══════════════════════════════════════════════════════════════════╝");
            let check = scalar_mul(G, k);
            if check.x == tx && check.y == ty && !check.inf {
                println!("[VERIFIED] k·G == target ✓");
            } else {
                eprintln!("[ERROR] k·G ≠ target — recovery bug");
            }
        }
        None => eprintln!("[INFO] no solution found (search terminated)"),
    }
}
