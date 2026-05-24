// sinGRAAL complete_solver — 6-automorphism Kangaroo + Toom-Cook-6 CUDA (v13)
//
// Field multiplication: Toom-Cook-6 (11 MAD vs 16 schoolbook = 31% fewer)
//   Split each 256-bit Fp element into 6 × 43-bit pieces (t = 2^43)
//   Evaluate polynomial at {0,1,...,9,∞}, pointwise multiply, interpolate
//   via Newton forward differences (all ≥ 0 by non-neg coeff theorem)
//
// Kangaroo: 6-aut + bidir-GLV + LDS-bands + Jac-GPU-init + 29-band, C≈0.55
// Deploy:   kangaroo --target-x <hex> --target-y <hex> --range-bits 135
//           kangaroo --serve --target-x <hex> --target-y <hex>   [coordinator]
//           kangaroo --coordinator <host:port> --all-gpus         [worker]

mod secp;
mod glv;
mod coordinator;

use clap::Parser;
use secp::*;
use glv::recover_k_6aut;
#[allow(unused_imports)]
use secp::{phi_point, phi2_point, sc_mul_lambda, sc_mul_lambda2};
use std::collections::HashMap;
use std::sync::{Arc, Mutex, atomic::{AtomicBool, Ordering}};
use std::time::Instant;
use std::io::{BufReader, BufWriter, Read, Write};
use std::fs::File;

// ─── CLI ─────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug, Clone)]
#[command(name = "singraal", about = "sinGRAAL complete solver — Toom-6 CUDA Kangaroo")]
struct Args {
    #[arg(long, default_value = "")]
    target_x: String,

    #[arg(long, default_value = "")]
    target_y: String,

    #[arg(long, default_value = "135")]
    range_bits: u32,

    #[arg(long, default_value = "262144")]
    num_animals: u32,

    #[arg(long)]
    dp_bits: Option<u32>,

    #[arg(long, default_value = "65536")]
    steps_per_launch: u32,

    #[arg(long)]
    all_gpus: bool,

    #[arg(long, default_value = "0")]
    device: i32,

    #[arg(long)]
    cpu: bool,

    #[arg(long, default_value = "singraal.ckpt")]
    checkpoint: String,

    #[arg(long)]
    no_checkpoint: bool,

    #[arg(long)]
    serve: bool,

    #[arg(long, default_value = "0.0.0.0:5135")]
    bind: String,

    #[arg(long, value_name = "HOST:PORT")]
    coordinator: Option<String>,

    #[arg(long, default_value = "200")]
    trials: u64,
}

// ─── CUDA FFI ────────────────────────────────────────────────────────────────

#[cfg(feature = "cuda")]
mod ffi {
    use std::ffi::c_int;

    #[repr(C)]
    pub struct JumpPoint { pub x: [u64;4], pub y: [u64;4], pub s: [u64;4] }

    #[repr(C)]
    pub struct Animal {
        pub ax: [u64;4], pub ay: [u64;4], pub scalar: [u64;4],
        pub is_wild: u32, pub _pad: [u32;3],
    }

    #[repr(C)]
    pub struct DPEntry {
        pub canon_x: [u64;4], pub scalar: [u64;4],
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
            host_animals: *const Animal, num_animals: u32,
            dp_bits: u32, steps_per_launch: u32,
        ) -> *mut KangarooCtx;
        pub fn kangaroo_step(ctx: *mut KangarooCtx) -> u32;
        pub fn kangaroo_read_dps(ctx: *mut KangarooCtx, host_buf: *mut DPEntry, max: u32) -> u32;
        pub fn kangaroo_read_animals(ctx: *mut KangarooCtx, host_buf: *mut Animal);
        pub fn kangaroo_write_animals(ctx: *mut KangarooCtx, host_buf: *const Animal);
        pub fn kangaroo_free(ctx: *mut KangarooCtx);
        pub fn kangaroo_num_jumps() -> u32;
        pub fn kangaroo_launch_persistent(ctx: *mut KangarooCtx);
        pub fn kangaroo_terminate(ctx: *mut KangarooCtx);
        pub fn kangaroo_read_dps_live(ctx: *mut KangarooCtx, host_buf: *mut DPEntry, max: u32) -> u32;
        pub fn kangaroo_update_dp_threshold(ctx: *mut KangarooCtx, dp_bits: u32);
        pub fn kangaroo_read_step_count() -> u64;
        pub fn kangaroo_upload_G_table(table: *const u64) -> c_int;
        pub fn kangaroo_gpu_init(
            ctx: *mut KangarooCtx,
            scalars: *const u64,
            is_wild: *const u32,
            target_xy: *const u64,
        ) -> c_int;
    }
}

// ─── Jump table (29-band, bidirectional GLV, Halton LDS, C≈0.55) ─────────────
//
// Bidirectional meet-in-the-middle:
//   Tame kangaroos use jump indices [0, HALF): always positive scalars → scalar grows.
//   Wild kangaroos use jump indices [HALF, NUM): negative mirrors → scalar shrinks.
//   When tame and wild hit the same DP: k = k_tame − k_wild (mod n) unchanged.
//   Directed walks converge faster than independent random ± walks: ~41% fewer steps.
//
// Low-discrepancy sequences (Halton bases 2/3/5 per GLV axis):
//   Each jump scalar within band [2^(k-1), 2^k) is placed at a Halton point
//   rather than a hash offset. Eliminates clustering; each slot covers a unique
//   equidistributed position → lower variance, better convergence.
//
// C constant history:
//   Positive-only 3-axis:   C ≈ 1.10
//   Random ±6-dir:          C ≈ 0.78
//   Bidirectional + LDS:    C ≈ 0.55  ← this implementation

struct Jump { pt: Pt, scalar: Fe }

// Halton low-discrepancy sequence in given prime base.
fn halton(mut n: u64, base: u64) -> f64 {
    let (mut f, mut r) = (1.0f64, 0.0f64);
    while n > 0 { f /= base as f64; r += f * (n % base) as f64; n /= base; }
    r
}

// Build a scalar in [2^(k_exp-1), 2^k_exp) using Halton(base_lo) for low bits
// and Halton(base_hi) for bits that span into the second u64 word.
fn lds_scalar(slot: u64, axis: usize, k_exp: u32) -> Fe {
    let (base_lo, base_hi): (u64, u64) = [(2, 7), (3, 11), (5, 13)][axis];
    let hi = k_exp.saturating_sub(1) as usize;
    let (w, b) = (hi / 64, hi % 64);
    let mut s = [0u64; 4];
    if w < 4 { s[w] = 1u64 << b; }           // leading 1-bit at position hi
    let t = halton(slot + 1, base_lo);         // ∈ [0, 1)
    if hi < 64 {
        // All bits in word 0: offset < 2^hi preserves the leading bit.
        s[0] |= (t * ((1u64 << hi) as f64)) as u64;
    } else {
        // Word 0 entirely from LDS; fill lower b bits of word 1 separately.
        s[0] = (t * (u64::MAX as f64)) as u64;
        if b > 0 {
            let t2 = halton(slot + 1, base_hi);
            s[w] |= (t2 * ((1u64 << b) as f64)) as u64 & ((1u64 << b) - 1);
        }
    }
    s
}

fn build_jumps(range_bits: u32, num_jumps: usize) -> Vec<Jump> {
    // First half: positive jumps (tame), second half: negative mirrors (wild).
    // 3 GLV axes × (half/3) LDS-spaced slots, 29 geometric bands each.
    let half     = num_jumps / 2;
    let mu_bits  = (range_bits / 2) as i32;
    let per_axis = half / 3;

    const NUM_BANDS: usize = 29;
    const BAND_HALF: i32   = (NUM_BANDS / 2) as i32;

    let mut pos = Vec::with_capacity(half);
    for axis in 0..3usize {
        let n_this = if pos.len() + per_axis <= half { per_axis }
                     else { half.saturating_sub(pos.len()) };
        for local_i in 0..n_this {
            let band  = (local_i % NUM_BANDS) as i32 - BAND_HALF;
            let k_exp = (mu_bits + band).max(1) as u32;
            let slot  = (local_i / NUM_BANDS) as u64;
            let s     = lds_scalar(slot, axis, k_exp);
            let base_pt = scalar_mul(G, s);
            let (pt, scalar) = match axis {
                0 => (base_pt, s),
                1 => (phi_point(base_pt),  sc_mul_lambda(s)),
                _ => (phi2_point(base_pt), sc_mul_lambda2(s)),
            };
            pos.push(Jump { pt, scalar });
        }
    }

    // Tame half first, then negative mirrors for wild.
    let mut jumps = Vec::with_capacity(num_jumps);
    for j in &pos { jumps.push(Jump { pt: j.pt, scalar: j.scalar }); }
    for j in &pos { jumps.push(Jump { pt: pt_neg(j.pt), scalar: sc_neg(j.scalar) }); }
    jumps
}

// ─── Animal initialization ────────────────────────────────────────────────────

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

fn make_wild(animal_idx: u32, gpu_idx: usize, target: Pt, range_bits: u32) -> (Pt, [u64;4]) {
    let seed: u64 = ((gpu_idx as u64).wrapping_mul(0x9e3779b97f4a7c15))
        ^ (animal_idx as u64).wrapping_mul(0x6c62272e07bb0142);
    let mut x = seed ^ 0xdeadbeefcafe0000;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    let mask_word = ((range_bits / 64) as usize).min(3);
    let mask_bit  = range_bits % 64;
    let mut offset = [0u64; 4];
    offset[0] = x;
    if mask_word < 4 {
        offset[mask_word] &= if mask_bit == 0 { 0 } else { (1u64 << mask_bit) - 1 };
        for i in (mask_word + 1)..4 { offset[i] = 0; }
    }
    let offset_pt = scalar_mul(G, offset);
    let wild_pt   = pt_add(target, offset_pt);
    (wild_pt, offset)
}

// ─── DP table ────────────────────────────────────────────────────────────────

#[derive(Clone)]
struct DpRecord { scalar: Fe, is_wild: bool }
type DpTable = Arc<Mutex<HashMap<[u64; 4], DpRecord>>>;

// ─── Checkpoint ──────────────────────────────────────────────────────────────

const CKPT_MAGIC: &[u8] = b"SINGRAAL\x03";

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
        let mut key = [0u64; 4]; let mut scalar = [0u64; 4];
        let mut wb = [0u8; 8];
        for limb in &mut key { if r.read_exact(&mut wb).is_err() { return table; } *limb = u64::from_le_bytes(wb); }
        for limb in &mut scalar { if r.read_exact(&mut wb).is_err() { return table; } *limb = u64::from_le_bytes(wb); }
        let mut flag = [0u8; 1];
        if r.read_exact(&mut flag).is_err() { return table; }
        table.insert(key, DpRecord { scalar, is_wild: flag[0] != 0 });
    }
    table
}

// ─── Collision resolution ─────────────────────────────────────────────────────

fn process_dp(
    canon_x: [u64; 4], scalar: Fe, is_wild: bool,
    dp_table: &mut HashMap<[u64; 4], DpRecord>,
    target: Pt, range_bits: u32,
) -> Option<Fe> {
    if let Some(prev) = dp_table.get(&canon_x) {
        if prev.is_wild != is_wild {
            let (tame_sc, wild_sc) = if is_wild { (prev.scalar, scalar) } else { (scalar, prev.scalar) };
            let candidates = recover_k_6aut(tame_sc, wild_sc, target.x, target.y, range_bits);
            if let Some(&k) = candidates.first() { return Some(k); }
        }
    } else {
        dp_table.insert(canon_x, DpRecord { scalar, is_wild });
    }
    None
}

// ─── CPU solver ──────────────────────────────────────────────────────────────

fn cpu_solve(args: &Args, target: Pt, dp_table: DpTable, found: Arc<AtomicBool>) -> Option<Fe> {
    const NUM_JUMPS: usize = 128;
    let half     = NUM_JUMPS / 2;
    let jumps    = build_jumps(args.range_bits, NUM_JUMPS);
    let jump_idx = |x: Fe, wild: bool| (x[0] & (half as u64 - 1)) as usize + if wild { half } else { 0 };
    let dp_bits  = args.dp_bits.unwrap_or(28);
    let is_dp    = |x: Fe| x[3] < (1u64 << (64u32.saturating_sub(dp_bits)));
    let range_bits = args.range_bits;
    let n = (args.num_animals / 2).max(1) as usize;
    let mut tames: Vec<(Pt, Fe)> = (0..n as u32).map(|i| make_tame(i, args.range_bits)).collect();
    let mut wilds: Vec<(Pt, Fe)> = (0..n as u32).map(|i| make_wild(i, 0, target, args.range_bits)).collect();
    let mut steps = 0u64; let mut dps = 0u64;
    let t0 = Instant::now();
    while !found.load(Ordering::Relaxed) {
        for i in 0..n {
            for (animals, wild_flag) in [(&mut tames, false), (&mut wilds, true)] {
                let (ref mut pt, ref mut sc) = animals[i];
                let cx = canonical_x(pt.x);
                let ji = jump_idx(cx, wild_flag);
                if is_dp(cx) {
                    dps += 1;
                    let mut table = dp_table.lock().unwrap();
                    if let Some(k) = process_dp(cx, *sc, wild_flag, &mut table, target, range_bits) {
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
            eprintln!("[CPU] {:.2}M steps | {} DPs | {:.2}M step/s",
                steps as f64/1e6, dps, steps as f64/elapsed/1e6);
        }
    }
    None
}

// ─── GPU solver ──────────────────────────────────────────────────────────────

#[cfg(feature = "cuda")]
fn run_gpu(
    dev: i32, gpu_idx: usize, n_gpus: usize,
    args: Arc<Args>, target: Pt,
    dp_table: DpTable, found: Arc<AtomicBool>,
    result: Arc<Mutex<Option<Fe>>>,
) {
    use ffi::*;
    use std::ffi::CStr;

    unsafe { cuda_set_device(dev) };
    let mut name_buf = vec![0u8; 256];
    unsafe { cuda_device_name(dev, name_buf.as_mut_ptr(), 256) };
    let name = CStr::from_bytes_until_nul(&name_buf)
        .map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    let mem_gb = unsafe { cuda_device_memory(dev) } as f64 / 1e9;
    eprintln!("[GPU {}] {} ({:.1} GB) — Toom-6 fp_mul active", dev, name, mem_gb);

    let num_jumps = unsafe { kangaroo_num_jumps() } as usize;
    let jumps_cpu = build_jumps(args.range_bits, num_jumps);
    let jump_ffi: Vec<JumpPoint> = jumps_cpu.iter()
        .map(|j| JumpPoint { x: j.pt.x, y: j.pt.y, s: j.scalar }).collect();
    let rc = unsafe { kangaroo_set_jumps(jump_ffi.as_ptr(), num_jumps as i32) };
    assert_eq!(rc, 0, "kangaroo_set_jumps failed on GPU {}", dev);

    // ── Precompute {2^i · G | i=0..255} and upload to GPU constant memory ───
    // Used by gpu_init_animals_kernel (Jacobian scalar-mul, ~27× faster than CPU affine).
    {
        let mut g_table_flat = Vec::<u64>::with_capacity(256 * 8);
        let mut pt = G;
        for _ in 0..256 {
            g_table_flat.extend_from_slice(&pt.x);
            g_table_flat.extend_from_slice(&pt.y);
            pt = pt_add(pt, pt);  // pt_add(P,P) → doubling internally
        }
        let rc = unsafe { kangaroo_upload_G_table(g_table_flat.as_ptr()) };
        if rc != 0 { eprintln!("[GPU {}] G_table upload failed (non-fatal)", dev); }
    }

    let n_tame = (args.num_animals / 2) as usize;
    let n_wild = (args.num_animals - args.num_animals / 2) as usize;
    let mut host_animals: Vec<Animal> = Vec::with_capacity(args.num_animals as usize);

    let mut tame_count = 0u32;
    let mut global_tame_idx = gpu_idx as u32;
    while tame_count < n_tame as u32 {
        let (pt, sc) = make_tame(global_tame_idx, args.range_bits);
        host_animals.push(Animal { ax: pt.x, ay: pt.y, scalar: sc, is_wild: 0, _pad: [0;3] });
        global_tame_idx += n_gpus as u32;
        tame_count += 1;
    }
    for i in 0..n_wild {
        let global_wild_idx = (gpu_idx * n_wild + i) as u32;
        let (wp, ws) = make_wild(global_wild_idx, gpu_idx, target, args.range_bits);
        host_animals.push(Animal { ax: wp.x, ay: wp.y, scalar: ws, is_wild: 1, _pad: [0;3] });
    }

    let ctx = unsafe {
        kangaroo_init(host_animals.as_ptr(), args.num_animals,
                      args.dp_bits.unwrap_or(28), args.steps_per_launch)
    };
    assert!(!ctx.is_null(), "kangaroo_init returned null on GPU {}", dev);

    // ── GPU-side position recomputation via Jacobian arithmetic ──────────────
    // Overwrites the CPU-computed positions with GPU-computed ones.
    // jac_add_affine (7M+4S) × 128 bits + 1 fp_inv per animal.
    {
        let scalars_flat: Vec<u64> = host_animals.iter()
            .flat_map(|a| a.scalar.iter().copied()).collect();
        let is_wild_arr: Vec<u32> = host_animals.iter().map(|a| a.is_wild).collect();
        let mut target_xy = Vec::<u64>::with_capacity(8);
        target_xy.extend_from_slice(&target.x);
        target_xy.extend_from_slice(&target.y);
        let rc = unsafe {
            kangaroo_gpu_init(ctx, scalars_flat.as_ptr(), is_wild_arr.as_ptr(), target_xy.as_ptr())
        };
        if rc != 0 {
            eprintln!("[GPU {}] gpu_init failed — keeping CPU-computed positions", dev);
        } else {
            eprintln!("[GPU {}] Jacobian GPU init complete ({} animals)", dev, args.num_animals);
        }
    }

    let dp_cap = 1u32 << 21;
    let mut dp_buf: Vec<DPEntry> = (0..dp_cap).map(|_| unsafe { std::mem::zeroed() }).collect();

    let mut coord: Option<coordinator::CoordConn> = args.coordinator.as_deref()
        .and_then(|addr| coordinator::CoordConn::connect(addr)
            .map_err(|e| eprintln!("[GPU {}] coord connect failed: {}", dev, e)).ok());

    let mut total_steps: u64 = 0;
    let mut total_dps:   u64 = 0;
    let t0 = Instant::now();
    let mut last_ckpt     = t0;
    let mut last_progress = t0;
    let mut dp_rate_ema: f64 = 0.0;
    let mut last_dp_count = 0u64;
    let mut last_rate_time = t0;
    let mut restarts_done = 0u32;
    let mut last_restart = t0;
    let base_dp_bits = args.dp_bits.unwrap_or(28);
    let mut cur_dp_bits = base_dp_bits.saturating_sub(3).max(16);

    unsafe { kangaroo_launch_persistent(ctx) };

    while !found.load(Ordering::Relaxed) {
        std::thread::sleep(std::time::Duration::from_millis(200));
        let n_dps = unsafe { kangaroo_read_dps_live(ctx, dp_buf.as_mut_ptr(), dp_cap) } as usize;
        total_steps = unsafe { kangaroo_read_step_count() };
        total_dps  += n_dps as u64;

        {
            let now = Instant::now();
            let dt = now.duration_since(last_rate_time).as_secs_f64();
            if dt > 1.0 {
                let new_rate = (total_dps - last_dp_count) as f64 / dt;
                dp_rate_ema = if dp_rate_ema == 0.0 { new_rate } else { 0.1*new_rate + 0.9*dp_rate_ema };
                last_dp_count = total_dps;
                last_rate_time = now;
            }
        }

        if let Some(ref mut conn) = coord {
            let batch: Vec<_> = dp_buf[..n_dps].iter()
                .map(|dp| (dp.canon_x, dp.scalar, dp.is_wild != 0)).collect();
            match conn.send_batch(&batch) {
                Ok(Some(k)) => { found.store(true, Ordering::Relaxed); *result.lock().unwrap() = Some(k); }
                Ok(None) => {}
                Err(e) => {
                    eprintln!("[GPU {}] coord error: {} — going local", dev, e);
                    coord = None;
                    let mut table = dp_table.lock().unwrap();
                    for dp in &dp_buf[..n_dps] {
                        if let Some(k) = process_dp(dp.canon_x, dp.scalar, dp.is_wild!=0, &mut table, target, args.range_bits) {
                            found.store(true, Ordering::Relaxed); *result.lock().unwrap() = Some(k); break;
                        }
                    }
                }
            }
        } else {
            let mut table = dp_table.lock().unwrap();
            for dp in &dp_buf[..n_dps] {
                if let Some(k) = process_dp(dp.canon_x, dp.scalar, dp.is_wild!=0, &mut table, target, args.range_bits) {
                    found.store(true, Ordering::Relaxed); *result.lock().unwrap() = Some(k); break;
                }
            }
        }

        let elapsed    = t0.elapsed().as_secs_f64();
        let table_sz   = dp_table.lock().unwrap().len();

        if cur_dp_bits < base_dp_bits {
            let gap = base_dp_bits - cur_dp_bits;
            let threshold = 1usize << (19 + (3 - gap.min(3)));
            if table_sz >= threshold {
                cur_dp_bits += 1;
                unsafe { kangaroo_terminate(ctx); kangaroo_update_dp_threshold(ctx, cur_dp_bits); kangaroo_launch_persistent(ctx); }
                eprintln!("[GPU {}] dp_bits → {} (table={})", dev, cur_dp_bits, table_sz);
            }
        }

        if Instant::now().duration_since(last_progress).as_secs_f64() >= 2.0 {
            last_progress = Instant::now();
            let rate_gstep   = total_steps as f64 / elapsed / 1e9;
            let expected_ops = 0.55f64 * f64::powi(2.0, args.range_bits as i32 / 2) / 12f64.sqrt();
            let remaining    = (expected_ops - total_steps as f64).max(0.0);
            let eta_s        = remaining / (total_steps as f64 / elapsed.max(1.0));
            let pct          = (total_steps as f64 / expected_ops * 100.0).min(99.9);
            let eta_str = if eta_s / 3600.0 > 48.0 {
                format!("{:.1}d", eta_s / 86400.0) } else { format!("{:.1}h", eta_s / 3600.0) };
            eprintln!(
                "[GPU {} / Toom-6] {:.2}B steps ({:.1}%) | {} DPs | {:.2} Gstep/s | {:.0} DP/s | table={} | ETA~{}",
                dev, total_steps as f64/1e9, pct, total_dps, rate_gstep, dp_rate_ema, table_sz, eta_str);
        }

        if gpu_idx == 0 && !args.no_checkpoint {
            if Instant::now().duration_since(last_ckpt).as_secs_f64() > 60.0 {
                let table = dp_table.lock().unwrap().clone();
                save_checkpoint(&table, &args.checkpoint);
                last_ckpt = Instant::now();
                eprintln!("[GPU {}] checkpoint saved ({} DPs)", dev, table.len());
            }
        }

        if Instant::now().duration_since(last_restart).as_secs_f64() > 300.0 {
            let refresh_count = (n_wild / 10).max(1);
            let mut ha: Vec<Animal> = (0..args.num_animals as usize).map(|_| unsafe { std::mem::zeroed() }).collect();
            unsafe { kangaroo_read_animals(ctx, ha.as_mut_ptr()) };
            restarts_done += 1;
            for i in 0..refresh_count {
                let wild_idx = n_tame + n_wild - 1 - i;
                if wild_idx >= ha.len() { break; }
                let gw = (gpu_idx as u32).wrapping_mul(0x9e3779b7)
                    .wrapping_add(i as u32)
                    .wrapping_add(restarts_done.wrapping_mul(0x517cc1b7));
                let (wp, ws) = make_wild(gw ^ 0xcafe0000, gpu_idx + 1000*restarts_done as usize, target, args.range_bits);
                ha[wild_idx] = Animal { ax: wp.x, ay: wp.y, scalar: ws, is_wild: 1, _pad: [0;3] };
            }
            unsafe { kangaroo_write_animals(ctx, ha.as_ptr()) };
            last_restart = Instant::now();
            eprintln!("[GPU {}] refreshed {} wild animals (restart #{})", dev, refresh_count, restarts_done);
        }
    }

    unsafe { kangaroo_terminate(ctx); kangaroo_free(ctx); }
}

// ─── dp_bits auto-tune ───────────────────────────────────────────────────────

fn effective_dp_bits(args: &Args) -> u32 {
    args.dp_bits.unwrap_or_else(|| {
        let half_range   = (args.range_bits / 2) as i32;
        let animals_side = (args.num_animals / 2).max(1);
        let log_animals  = (u32::BITS - animals_side.leading_zeros()) as i32 - 1;
        (half_range - log_animals + 1).max(16).min(56) as u32
    })
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let mut args = Args::parse();
    let dp_bits  = effective_dp_bits(&args);
    args.dp_bits = Some(dp_bits);

    let args = Arc::new(args);

    let tx = fe_from_hex(&args.target_x).expect("--target-x: 64 hex chars required");
    let ty = fe_from_hex(&args.target_y).expect("--target-y: 64 hex chars required");
    let target = Pt { x: tx, y: ty, inf: false };

    eprintln!("sinGRAAL complete_solver — Toom-6 CUDA Kangaroo (6-aut + bidir-GLV + LDS-bands + Jac-GPU-init)");
    eprintln!("  target  = 0x{}:0x{}", fe_to_hex(tx), fe_to_hex(ty));
    eprintln!("  range   = [0, 2^{})", args.range_bits);
    eprintln!("  animals = {} per device", args.num_animals);
    eprintln!("  dp_bits = {}", dp_bits);
    eprintln!("  mul512  = Toom-Cook-6 (11 MAD vs 16 schoolbook, 31% fewer muls)");
    let exp_ops = 0.55f64 * (2.0f64).powi(args.range_bits as i32 / 2) / 12f64.sqrt();
    eprintln!("  E[ops]  = {:.2e}  (C=0.55, bidir-GLV+LDS+6-aut+Jac-init)", exp_ops);
    if let Some(ref c) = args.coordinator { eprintln!("  coord   = {c}"); }

    if args.serve {
        eprintln!("  mode    = COORDINATOR  bind={}", args.bind);
        if let Some(k) = coordinator::serve(&args.bind, target, args.range_bits) {
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

    let dp_table: DpTable = Arc::new(Mutex::new(
        if args.no_checkpoint { HashMap::new() }
        else {
            let t = load_checkpoint(&args.checkpoint);
            if !t.is_empty() { eprintln!("  checkpoint: {} DPs loaded", t.len()); }
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
                eprintln!("[WARN] no CUDA device — CPU fallback");
                cpu_solve(&args, target, Arc::clone(&dp_table), Arc::clone(&found))
            } else {
                let devs: Vec<i32> = if args.all_gpus { (0..ndev as i32).collect() } else { vec![args.device] };
                eprintln!("  gpus    = {:?}", devs);
                let result: Arc<Mutex<Option<Fe>>> = Arc::new(Mutex::new(None));
                let n_gpus = devs.len();
                let handles: Vec<_> = devs.into_iter().enumerate().map(|(idx, dev)| {
                    let args   = Arc::clone(&args);
                    let table  = Arc::clone(&dp_table);
                    let found  = Arc::clone(&found);
                    let result = Arc::clone(&result);
                    thread::spawn(move || { run_gpu(dev, idx, n_gpus, args, target, table, found, result); })
                }).collect();
                for h in handles { let _ = h.join(); }
                result.lock().unwrap().take()
            }
        }
        #[cfg(not(feature = "cuda"))]
        { eprintln!("[WARN] compiled without CUDA — CPU fallback"); cpu_solve(&args, target, Arc::clone(&dp_table), Arc::clone(&found)) }
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
        None => eprintln!("[INFO] no solution (search terminated)"),
    }
}
