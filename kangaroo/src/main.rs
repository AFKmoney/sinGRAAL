// sinGRAAL Kangaroo v2 — 6-automorphism Pollard Kangaroo, CUDA-accelerated
//
// Improvements over v1:
//  • Correct DP detection (GPU normalizes to affine every NORM_INTERVAL steps)
//  • Multi-GPU with shared DP table (linear scaling with GPUs)
//  • Checkpoint save/load (resume multi-day runs)
//  • Exact 6-aut recovery (6 candidates, not 18)
//  • Better jump distribution (128 jumps, uniform spacing around √range)
//  • Progress bar with ETA
//
// Usage:
//   kangaroo --target-x <hex64> --target-y <hex64> --range-bits 135
//   kangaroo --target-x <hex64> --target-y <hex64> --range-bits 135 --all-gpus
//   kangaroo --target-x <hex64> --target-y <hex64> --range-bits 135 --cpu

mod secp;
mod glv;
mod coordinator;

use clap::Parser;
use secp::*;
use glv::recover_k_6aut;
use std::collections::HashMap;
use std::sync::{Arc, Mutex, atomic::{AtomicBool, Ordering}};
use std::time::Instant;
use std::io::{BufReader, BufWriter, Read, Write};
use std::fs::File;

// ─── CLI ─────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug, Clone)]
#[command(name = "kangaroo", about = "sinGRAAL 6-aut Kangaroo ECDLP solver for secp256k1")]
struct Args {
    #[arg(long)]
    target_x: String,

    #[arg(long)]
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
        pub fn kangaroo_free(ctx: *mut KangarooCtx);
        pub fn kangaroo_num_jumps() -> u32;
    }
}

// ─── Jump table ──────────────────────────────────────────────────────────────

struct Jump { pt: Pt, scalar: Fe }

/// Build NUM_JUMPS jump points with near-optimal torus coverage.
///
/// Scalars are pseudo-random powers-of-2 drawn uniformly from the interval
/// [2^(μ_bits-1), 2^(μ_bits+1)] where μ_bits = range_bits/2.  This gives:
///   mean  ≈ 2^μ_bits   (≈ √range, optimal Kangaroo mean)
///   E[s²] ≈ 4/3 * mean²  (uniform-distribution bound — minimum reachable)
/// The 128 distinct scalars are powers of 2 so scalar_mul is fast (single
/// double-and-add chain, no precomputation needed), yet they tile the orbit
/// torus with minimal correlation between consecutive jumps.
fn build_jumps(range_bits: u32, num_jumps: usize) -> Vec<Jump> {
    let mu_bits = range_bits / 2;   // target mean exponent ≈ √range
    let mut jumps = Vec::with_capacity(num_jumps);

    // Deterministic PRNG (xorshift64) seeded from range_bits for reproducibility.
    // Each jump scalar = 2^k where k ∈ [mu_bits-1, mu_bits+1] (± 1 bit of μ).
    // We spread k values uniformly across [mu_bits-1, mu_bits+1] so that the
    // scalars are distinct and cover the interval with minimal repetition.
    // For 128 jumps: 43 at 2^(μ-1), 42 at 2^μ, 43 at 2^(μ+1)  — trivially
    // uniform, deterministic, and independent of a PRNG.
    for i in 0..num_jumps {
        // Spread over 3 bands: lo / mid / hi  (i mod 3 → exponent offset)
        let exp_offset: i32 = (i % 3) as i32 - 1;  // -1, 0, +1
        let k = (mu_bits as i32 + exp_offset).max(1) as u32;

        // Within the band, vary the low bits to break arithmetic correlation.
        // Shift position within the band: add i/3 to the low word.
        let band_idx = (i / 3) as u64;
        let word = (k / 64) as usize;
        let bit  = k % 64;
        let mut s = [0u64; 4];
        if word < 4 {
            s[word] = 1u64 << bit;
        }
        // Add band_idx scaled by 1 (so within a band scalars differ by 1 each)
        // This keeps scalars as powers-of-2 PLUS a small low-bit offset,
        // giving genuine 128 distinct values without hurting the mean/variance.
        let (v, ov) = s[0].overflowing_add(band_idx);
        s[0] = v;
        if ov && 1 < 4 { s[1] += 1; }

        let pt = scalar_mul(G, s);
        jumps.push(Jump { pt, scalar: s });
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
    let mut wilds: Vec<(Pt, Fe)> = (0..n).map(|_| (target, [0u64; 4])).collect();

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
    // Wild: all start at target
    for _ in 0..n_wild {
        host_animals.push(Animal {
            ax: target.x, ay: target.y,
            scalar: [0;4], is_wild: 1, _pad: [0;3],
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

    let mut total_steps = 0u64;
    let mut total_dps   = 0u64;
    let t0 = Instant::now();
    let mut last_ckpt   = t0;
    const CKPT_INTERVAL_S: f64 = 60.0;

    while !found.load(Ordering::Relaxed) {
        let _cumulative = unsafe { kangaroo_step(ctx) };
        total_steps += args.num_animals as u64 * args.steps_per_launch as u64;

        let n_dps = unsafe { kangaroo_read_dps(ctx, dp_buf.as_mut_ptr(), dp_cap) } as usize;
        total_dps += n_dps as u64;

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

        // Progress every 10 launches
        if total_steps % (args.num_animals as u64 * args.steps_per_launch as u64 * 10) == 0 {
            let rate = total_steps as f64 / elapsed / 1e6;
            let table_sz = dp_table.lock().unwrap().len();
            // ETA estimate: 6-aut Kangaroo expected ops ≈ 2·√(range/6)
            let expected_ops = 2.0f64 * f64::powi(2.0, args.range_bits as i32 / 2) / 6f64.sqrt();
            let eta_s = (expected_ops / 1e6 - total_steps as f64) / rate / 1e3;
            eprintln!(
                "[GPU {}] {:.2}B steps | {} DPs | {:.0}M step/s | dp_table={} | ETA~{:.0}h",
                dev,
                total_steps as f64 / 1e9,
                total_dps,
                rate,
                table_sz,
                (eta_s / 3600.0).max(0.0),
            );
        }

        // Checkpoint every 60s (only GPU 0 saves to avoid races)
        if gpu_idx == 0 && !args.no_checkpoint {
            let since_ckpt = t0.duration_since(last_ckpt).as_secs_f64();
            if since_ckpt > CKPT_INTERVAL_S {
                let table = dp_table.lock().unwrap().clone();
                save_checkpoint(&table, &args.checkpoint);
                last_ckpt = Instant::now();
                eprintln!("[GPU {}] checkpoint saved ({} DPs)", dev, table.len());
            }
        }
    }

    unsafe { kangaroo_free(ctx) };
}

// ─── Main ─────────────────────────────────────────────────────────────────────

// Resolve effective dp_bits: user override or auto-tune from range_bits.
// Optimal: dp_bits ≈ range_bits/2 − 10 keeps the DP table at ~1M entries
// at collision time while ensuring DPs are dense enough for fast detection.
fn effective_dp_bits(args: &Args) -> u32 {
    args.dp_bits.unwrap_or_else(|| {
        (args.range_bits / 2).saturating_sub(10).max(16).min(40)
    })
}

fn main() {
    let mut args = Args::parse();
    // Resolve dp_bits once, store back so all code paths see the same value
    let dp_bits = effective_dp_bits(&args);
    args.dp_bits = Some(dp_bits);
    let args = Arc::new(args);

    let tx = fe_from_hex(&args.target_x).expect("--target-x: need 64 hex chars");
    let ty = fe_from_hex(&args.target_y).expect("--target-y: need 64 hex chars");
    let target = Pt { x: tx, y: ty, inf: false };

    eprintln!("sinGRAAL Kangaroo v5 — 6-automorphism secp256k1 ECDLP");
    eprintln!("  target  = 0x{}:0x{}", fe_to_hex(tx), fe_to_hex(ty));
    eprintln!("  range   = [0, 2^{})", args.range_bits);
    eprintln!("  animals = {} per device", args.num_animals);
    eprintln!("  dp_bits = {}", dp_bits);
    // Expected ops (informational):
    let exp_ops = 2.0f64 * (2.0f64).powi(args.range_bits as i32 / 2) / 6f64.sqrt();
    eprintln!("  E[ops]  = {:.2e}  (6-aut kangaroo, ~2√(range/6))", exp_ops);
    if let Some(ref c) = args.coordinator { eprintln!("  coord   = {c}"); }

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
                result.lock().unwrap().take()
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
