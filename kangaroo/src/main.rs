// sinGRAAL Kangaroo — 6-automorphism Pollard Kangaroo, CUDA-accelerated
// Solves secp256k1 ECDLP: find k ∈ [0, 2^range_bits) such that k·G = target
//
// Usage:
//   kangaroo --target-x <hex> --target-y <hex> --range-bits 135 [options]

mod secp;

use clap::Parser;
use secp::*;
use std::collections::HashMap;
use std::time::Instant;

// ─── CLI ─────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(name = "kangaroo", about = "6-aut Kangaroo ECDLP solver for secp256k1")]
struct Args {
    /// Target point X coordinate (hex, 64 chars)
    #[arg(long)]
    target_x: String,

    /// Target point Y coordinate (hex, 64 chars)
    #[arg(long)]
    target_y: String,

    /// Key space: search in [0, 2^range_bits)
    #[arg(long, default_value = "135")]
    range_bits: u32,

    /// Number of Kangaroo animals (GPU threads)
    #[arg(long, default_value = "131072")]
    num_animals: u32,

    /// Distinguished-point difficulty (28 = ~1 per 2^28 steps)
    #[arg(long, default_value = "28")]
    dp_bits: u32,

    /// Steps per GPU launch (tune for latency vs. throughput)
    #[arg(long, default_value = "4096")]
    steps_per_launch: u32,

    /// CUDA device index
    #[arg(long, default_value = "0")]
    device: i32,

    /// Enable CPU-only mode (no CUDA)
    #[arg(long)]
    cpu: bool,
}

// ─── CUDA FFI ────────────────────────────────────────────────────────────────

#[cfg(feature = "cuda")]
mod cuda_ffi {
    use std::ffi::c_int;

    // Mirror of C structs (must match secp256k1.cuh layout)
    #[repr(C)]
    pub struct JumpPoint {
        pub x: [u64; 4],
        pub y: [u64; 4],
        pub s: [u64; 4],
    }

    #[repr(C)]
    pub struct Animal {
        pub x:       [u64; 4],
        pub y:       [u64; 4],
        pub z:       [u64; 4],
        pub scalar:  [u64; 4],
        pub is_wild: u32,
        pub _pad:    [u32; 3],
    }

    #[repr(C)]
    pub struct DPEntry {
        pub x_jac:  [u64; 4],
        pub z_jac:  [u64; 4],
        pub scalar: [u64; 4],
        pub is_wild: u32,
        pub _pad:   [u32; 3],
    }

    pub enum KangarooCtx {}

    extern "C" {
        pub fn cuda_device_count() -> c_int;
        pub fn cuda_device_name(dev: c_int, buf: *mut u8, len: c_int);
        pub fn cuda_device_memory(dev: c_int) -> u64;

        pub fn kangaroo_set_jumps(jumps: *const JumpPoint, n: c_int) -> c_int;
        pub fn kangaroo_init(
            host_animals: *const Animal,
            num_animals:  u32,
            dp_bits:      u32,
        ) -> *mut KangarooCtx;
        pub fn kangaroo_step(ctx: *mut KangarooCtx, steps_per_launch: u32) -> u32;
        pub fn kangaroo_read_dps(
            ctx:      *mut KangarooCtx,
            host_buf: *mut DPEntry,
            max:      u32,
        ) -> u32;
        pub fn kangaroo_free(ctx: *mut KangarooCtx);
    }
}

// ─── Jump table ──────────────────────────────────────────────────────────────

struct Jump { pt: Pt, scalar: Fe }

fn build_jumps(range_bits: u32, num_jumps: usize) -> Vec<Jump> {
    // Jump distances are powers of 2: 2^(range_bits/2 - num_jumps/2 + i)
    // Provides good mixing across the range
    let base_exp = (range_bits / 2).saturating_sub(num_jumps as u32 / 2);
    let mut jumps = Vec::with_capacity(num_jumps);
    for i in 0..num_jumps {
        let exp = base_exp + i as u32;
        // scalar = 2^exp
        let mut s = [0u64; 4];
        let word = (exp / 64) as usize;
        let bit  = exp % 64;
        if word < 4 { s[word] = 1u64 << bit; }
        let pt = scalar_mul(G, s);
        jumps.push(Jump { pt, scalar: s });
    }
    jumps
}

// ─── Animal initialization ────────────────────────────────────────────────────

fn make_tame(idx: u32, range_bits: u32) -> (Pt, Fe) {
    // Tame start: k_t = 2^(range_bits-1) + idx * stride
    // Spread across [2^(range_bits-1), 2^range_bits)
    let stride_bits = range_bits.saturating_sub(17); // ~2^17 animals max
    let mut k = [0u64; 4];
    let bit = (range_bits - 1) as usize;
    k[bit / 64] |= 1u64 << (bit % 64);
    // Add idx * 2^stride_bits
    let sbit = stride_bits as usize;
    let add_word = sbit / 64;
    if add_word < 4 {
        let (v, overflow) = k[add_word].overflowing_add((idx as u64) << (sbit % 64));
        k[add_word] = v;
        if overflow && add_word + 1 < 4 { k[add_word + 1] += 1; }
    }
    (scalar_mul(G, k), k)
}

fn make_wild(target: Pt) -> (Pt, Fe) {
    // Wild starts at target; scalar accumulator starts at 0
    (target, [0u64; 4])
}

// ─── DP table & collision detection ──────────────────────────────────────────

struct DpRecord {
    scalar:  Fe,
    is_wild: bool,
}

// Try to recover k from a tame/wild collision
// tame: scalar_t such that scalar_t*G = DP_point
// wild: scalar_w such that target + scalar_w*G = DP_point
// → target = (scalar_t - scalar_w)*G → k = scalar_t - scalar_w (mod n)
fn recover_k(tame_sc: Fe, wild_sc: Fe) -> Fe {
    sc_sub(tame_sc, wild_sc)
}

// ─── CPU-only solver ──────────────────────────────────────────────────────────

fn cpu_solve(args: &Args, target: Pt) -> Option<Fe> {
    const NUM_JUMPS: usize = 32;
    let jumps = build_jumps(args.range_bits, NUM_JUMPS);

    let jump_idx = |x: Fe| -> usize { (x[0] % NUM_JUMPS as u64) as usize };
    let is_dp    = |x: Fe| -> bool  { x[3] < (1u64 << (64 - args.dp_bits.min(63))) };

    // Initialize animals: half tame, half wild
    let n_tame = (args.num_animals / 2).max(1) as usize;
    let n_wild = n_tame;

    let mut tames: Vec<(Pt, Fe)> = (0..n_tame as u32).map(|i| make_tame(i, args.range_bits)).collect();
    let mut wilds: Vec<(Pt, Fe)> = (0..n_wild as u32).map(|_| make_wild(target)).collect();

    let mut dp_table: HashMap<[u64; 4], DpRecord> = HashMap::new();
    let mut total_steps = 0u64;
    let mut total_dps = 0u64;
    let t0 = Instant::now();

    loop {
        for i in 0..n_tame {
            let (ref mut pt, ref mut sc) = tames[i];
            let ji = jump_idx(pt.x);
            *pt = pt_add(*pt, jumps[ji].pt);
            *sc = sc_add(*sc, jumps[ji].scalar);
            total_steps += 1;

            let cx = canonical_x(pt.x);
            if is_dp(cx) {
                total_dps += 1;
                if let Some(prev) = dp_table.get(&cx) {
                    if prev.is_wild != false {
                        let k = recover_k(*sc, prev.scalar);
                        if !verify_k(k, target) { continue; }
                        print_stats(total_steps, total_dps, t0.elapsed().as_secs_f64());
                        return Some(k);
                    }
                } else {
                    dp_table.insert(cx, DpRecord { scalar: *sc, is_wild: false });
                }
            }
        }

        for i in 0..n_wild {
            let (ref mut pt, ref mut sc) = wilds[i];
            let ji = jump_idx(pt.x);
            *pt = pt_add(*pt, jumps[ji].pt);
            *sc = sc_add(*sc, jumps[ji].scalar);
            total_steps += 1;

            let cx = canonical_x(pt.x);
            if is_dp(cx) {
                total_dps += 1;
                if let Some(prev) = dp_table.get(&cx) {
                    if prev.is_wild == false {
                        // Try all 6-automorphism combinations for wild scalar adjustment
                        for k in try_recover_6aut(prev.scalar, *sc, cx, target) {
                            print_stats(total_steps, total_dps, t0.elapsed().as_secs_f64());
                            return Some(k);
                        }
                    }
                } else {
                    dp_table.insert(cx, DpRecord { scalar: *sc, is_wild: true });
                }
            }
        }

        if total_steps % 1_000_000 == 0 {
            let elapsed = t0.elapsed().as_secs_f64();
            eprintln!(
                "[CPU] steps={:.2}M  DPs={}  rate={:.2}Mstep/s",
                total_steps as f64 / 1e6,
                total_dps,
                total_steps as f64 / elapsed / 1e6
            );
        }
    }
}

// Try 6-automorphism k recovery combinations:
// The canonical form maps 6 orbit points to same canonical_x.
// When tame walks in canonical space and wild starts at target,
// the collision at canonical_x means:
//   tame_scalar*G = ψ^a(target) ± ψ^b(scalar_w * G) for some a,b in {0,1,2}
//   → k = ±λ^a * (tame_sc ∓ λ^b * wild_sc) (mod n)
fn try_recover_6aut(tame_sc: Fe, wild_sc: Fe, _canonical: Fe, target: Pt) -> Vec<Fe> {
    let lambdas = [
        [1u64, 0, 0, 0],  // λ^0 = 1
        LAMBDA,            // λ^1
        sc_mul(LAMBDA, LAMBDA), // λ^2
    ];
    let mut candidates = Vec::new();
    for &lam_a in &lambdas {
        for &lam_b in &lambdas {
            for sign in [true, false] {
                let adjusted_wild = sc_mul(lam_b, wild_sc);
                let k = if sign {
                    sc_mul(lam_a, sc_sub(tame_sc, adjusted_wild))
                } else {
                    sc_mul(lam_a, sc_add(tame_sc, adjusted_wild))
                };
                if verify_k(k, target) {
                    candidates.push(k);
                }
            }
        }
    }
    candidates
}

fn verify_k(k: Fe, target: Pt) -> bool {
    let computed = scalar_mul(G, k);
    computed.x == target.x && computed.y == target.y && !computed.inf
}

fn print_stats(steps: u64, dps: u64, elapsed: f64) {
    eprintln!(
        "[SOLVED] steps={:.2}M  DPs={}  elapsed={:.1}s  rate={:.2}Mstep/s",
        steps as f64 / 1e6,
        dps,
        elapsed,
        steps as f64 / elapsed / 1e6
    );
}

// ─── CUDA solver ─────────────────────────────────────────────────────────────

#[cfg(feature = "cuda")]
fn cuda_solve(args: &Args, target: Pt) -> Option<Fe> {
    use cuda_ffi::*;
    use std::ffi::CStr;

    // Print device info
    let ndev = unsafe { cuda_device_count() };
    if ndev == 0 {
        eprintln!("No CUDA devices found — falling back to CPU");
        return cpu_solve(args, target);
    }
    let dev = args.device;
    let mut name_buf = vec![0u8; 256];
    unsafe { cuda_device_name(dev, name_buf.as_mut_ptr(), 256) };
    let name = CStr::from_bytes_until_nul(&name_buf).map(|c| c.to_string_lossy().into_owned()).unwrap_or_default();
    let mem_gb = unsafe { cuda_device_memory(dev) } as f64 / 1e9;
    eprintln!("[CUDA] device {}: {} ({:.1} GB)", dev, name, mem_gb);

    const NUM_JUMPS: usize = 32;
    let jumps_cpu = build_jumps(args.range_bits, NUM_JUMPS);

    // Upload jump table
    let mut jump_ffi: Vec<JumpPoint> = jumps_cpu.iter().map(|j| JumpPoint {
        x: j.pt.x,
        y: j.pt.y,
        s: j.scalar,
    }).collect();
    let rc = unsafe { kangaroo_set_jumps(jump_ffi.as_mut_ptr(), NUM_JUMPS as i32) };
    assert_eq!(rc, 0, "kangaroo_set_jumps failed");

    // Build animals
    let n_tame = args.num_animals / 2;
    let n_wild = args.num_animals - n_tame;
    let mut host_animals: Vec<Animal> = Vec::with_capacity(args.num_animals as usize);

    for i in 0..n_tame {
        let (pt, sc) = make_tame(i, args.range_bits);
        host_animals.push(Animal {
            x: pt.x, y: pt.y, z: [1,0,0,0],
            scalar: sc, is_wild: 0, _pad: [0;3],
        });
    }
    for _ in 0..n_wild {
        let (pt, sc) = make_wild(target);
        host_animals.push(Animal {
            x: pt.x, y: pt.y, z: [1,0,0,0],
            scalar: sc, is_wild: 1, _pad: [0;3],
        });
    }

    let ctx = unsafe {
        kangaroo_init(host_animals.as_ptr(), args.num_animals, args.dp_bits)
    };
    assert!(!ctx.is_null(), "kangaroo_init returned null");

    let mut dp_table: HashMap<[u64; 4], DpRecord> = HashMap::new();
    let mut dp_host_buf: Vec<DPEntry> = vec![unsafe { std::mem::zeroed() }; 1 << 20];
    let mut total_steps = 0u64;
    let mut total_dps   = 0u64;
    let t0 = Instant::now();

    loop {
        let _count = unsafe { kangaroo_step(ctx, args.steps_per_launch) };
        total_steps += args.num_animals as u64 * args.steps_per_launch as u64;

        let n_dps = unsafe {
            kangaroo_read_dps(ctx, dp_host_buf.as_mut_ptr(), (1u32 << 20))
        } as usize;
        total_dps += n_dps as u64;

        for dp in &dp_host_buf[..n_dps] {
            // Convert Jacobian DP to affine canonical x
            let z_inv = fp_inv(dp.z_jac);
            let z2    = fp_sqr(z_inv);
            let ax    = fp_mul(dp.x_jac, z2);
            let cx    = canonical_x(ax);

            let is_wild = dp.is_wild != 0;
            if let Some(prev) = dp_table.get(&cx) {
                if prev.is_wild != is_wild {
                    let (tame_sc, wild_sc) = if is_wild {
                        (prev.scalar, dp.scalar)
                    } else {
                        (dp.scalar, prev.scalar)
                    };
                    for k in try_recover_6aut(tame_sc, wild_sc, cx, target) {
                        unsafe { kangaroo_free(ctx) };
                        print_stats(total_steps, total_dps, t0.elapsed().as_secs_f64());
                        return Some(k);
                    }
                }
            } else {
                dp_table.insert(cx, DpRecord { scalar: dp.scalar, is_wild });
            }
        }

        let elapsed = t0.elapsed().as_secs_f64();
        if (total_steps / (args.num_animals as u64 * args.steps_per_launch as u64)) % 20 == 0 {
            eprintln!(
                "[GPU]  steps={:.2}M  DPs={}  rate={:.2}Mstep/s  dp_table={}",
                total_steps as f64 / 1e6,
                total_dps,
                total_steps as f64 / elapsed / 1e6,
                dp_table.len(),
            );
        }
    }
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let args = Args::parse();

    let tx = fe_from_hex(&args.target_x).expect("invalid --target-x (need 64 hex chars)");
    let ty = fe_from_hex(&args.target_y).expect("invalid --target-y (need 64 hex chars)");
    let target = Pt { x: tx, y: ty, inf: false };

    eprintln!("sinGRAAL Kangaroo v0.1 — 6-automorphism secp256k1 ECDLP solver");
    eprintln!("  target  = {},{}", fe_to_hex(tx), fe_to_hex(ty));
    eprintln!("  range   = 2^{}", args.range_bits);
    eprintln!("  animals = {}", args.num_animals);
    eprintln!("  dp_bits = {}", args.dp_bits);
    eprintln!("  mode    = {}", if args.cpu { "CPU" } else { "CUDA" });

    let result = if args.cpu {
        cpu_solve(&args, target)
    } else {
        #[cfg(feature = "cuda")]
        { cuda_solve(&args, target) }
        #[cfg(not(feature = "cuda"))]
        {
            eprintln!("[WARN] compiled without CUDA support, falling back to CPU");
            cpu_solve(&args, target)
        }
    };

    match result {
        Some(k) => {
            println!("\n[SUCCESS] k = {}", fe_to_hex(k));
            // Also verify
            let check = scalar_mul(G, k);
            if check.x == tx && check.y == ty {
                println!("[VERIFIED] k*G == target ✓");
            } else {
                eprintln!("[ERROR] verification failed — k*G != target");
            }
        }
        None => {
            eprintln!("[INFO] search ended without result");
        }
    }
}
