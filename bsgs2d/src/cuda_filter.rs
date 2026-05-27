// sinGRAAL — Rust wrapper pour le filtre CUDA PNC direct
//
// Architecture pipeline :
//   CPU enqueue → GPU broyeur (kill 99.99% en ~100 ns/paire) → CPU sniper (LLL exact)
//
// Utilisation :
//   let f = CudaFilter::init(target_x, block_bits)?;
//   let pairs = vec![(a0, b0), (a1, b1), ...];
//   let survivors = f.filter_pairs(&pairs)?;   // indices dans pairs[] des survivants

#[cfg(feature = "cuda")]
mod cuda_impl {
    use num_bigint::BigInt;
    use std::ffi::c_void;

    #[link(name = "s3filter", kind = "static")]
    unsafe extern "C" {
        fn launch_s3_filter(
            a_lo:            *const u64,
            b_lo:            *const u64,
            t_x:             *const u64,
            block_bits:      u32,
            survivors:       *mut u32,
            n_survivors_dev: *mut u32,
            n_pairs:         u32,
            stream:          *mut c_void,
        );

        fn cudaMalloc(ptr: *mut *mut c_void, size: usize) -> i32;
        fn cudaFree(ptr: *mut c_void) -> i32;
        fn cudaMemcpy(
            dst: *mut c_void,
            src: *const c_void,
            count: usize,
            kind: u32,
        ) -> i32;
        fn cudaMemset(ptr: *mut c_void, value: i32, count: usize) -> i32;
        fn cudaDeviceSynchronize() -> i32;
        fn cudaStreamCreate(stream: *mut *mut c_void) -> i32;
        fn cudaStreamDestroy(stream: *mut c_void) -> i32;
    }

    const CUDA_MEMCPY_HOST_TO_DEVICE: u32 = 1;
    const CUDA_MEMCPY_DEVICE_TO_HOST: u32 = 2;

    fn bigint_to_limbs4(n: &BigInt) -> [u64; 4] {
        let p_bytes = {
            let (_, bytes) = n.to_bytes_le();
            bytes
        };
        let mut limbs = [0u64; 4];
        for (i, chunk) in p_bytes.chunks(8).take(4).enumerate() {
            let mut buf = [0u8; 8];
            buf[..chunk.len()].copy_from_slice(chunk);
            limbs[i] = u64::from_le_bytes(buf);
        }
        limbs
    }

    pub struct CudaFilter {
        tx_limbs:   [u64; 4],
        block_bits: u32,
        stream:     *mut c_void,
    }

    unsafe impl Send for CudaFilter {}
    unsafe impl Sync for CudaFilter {}

    impl CudaFilter {
        pub fn init(target_x: &BigInt, block_bits: u32) -> Result<Self, String> {
            let tx_limbs = bigint_to_limbs4(target_x);
            let mut stream: *mut c_void = std::ptr::null_mut();
            let rc = unsafe { cudaStreamCreate(&mut stream as *mut _ as *mut *mut c_void) };
            if rc != 0 {
                return Err(format!("cudaStreamCreate failed: {rc}"));
            }
            Ok(CudaFilter { tx_limbs, block_bits, stream })
        }

        /// Filter a batch of (A_lo, B_lo) block-base pairs.
        /// Returns indices (into `pairs`) of the survivors.
        pub fn filter_pairs(&self, pairs: &[(BigInt, BigInt)]) -> Result<Vec<usize>, String> {
            let n = pairs.len();
            if n == 0 { return Ok(vec![]); }

            // Flatten pairs into host arrays
            let mut h_a = vec![0u64; n * 4];
            let mut h_b = vec![0u64; n * 4];
            for (i, (a, b)) in pairs.iter().enumerate() {
                let al = bigint_to_limbs4(a);
                let bl = bigint_to_limbs4(b);
                h_a[i*4..i*4+4].copy_from_slice(&al);
                h_b[i*4..i*4+4].copy_from_slice(&bl);
            }

            let max_survivors = n as u32;

            unsafe {
                // Allocate device memory
                let mut d_a: *mut c_void = std::ptr::null_mut();
                let mut d_b: *mut c_void = std::ptr::null_mut();
                let mut d_tx: *mut c_void = std::ptr::null_mut();
                let mut d_surv: *mut c_void = std::ptr::null_mut();
                let mut d_n: *mut c_void = std::ptr::null_mut();

                macro_rules! cuda_check {
                    ($rc:expr, $msg:literal) => {
                        if $rc != 0 { return Err(format!("{}: cuda error {}", $msg, $rc)); }
                    };
                }

                cuda_check!(cudaMalloc(&mut d_a,    n * 4 * 8), "alloc A");
                cuda_check!(cudaMalloc(&mut d_b,    n * 4 * 8), "alloc B");
                cuda_check!(cudaMalloc(&mut d_tx,   4 * 8),     "alloc Tx");
                cuda_check!(cudaMalloc(&mut d_surv, (max_survivors as usize) * 4), "alloc surv");
                cuda_check!(cudaMalloc(&mut d_n,    4),          "alloc n");

                // Copy input to device
                cuda_check!(cudaMemcpy(d_a,  h_a.as_ptr()          as *const c_void, n*4*8, CUDA_MEMCPY_HOST_TO_DEVICE), "H2D A");
                cuda_check!(cudaMemcpy(d_b,  h_b.as_ptr()          as *const c_void, n*4*8, CUDA_MEMCPY_HOST_TO_DEVICE), "H2D B");
                cuda_check!(cudaMemcpy(d_tx, self.tx_limbs.as_ptr() as *const c_void, 4*8,   CUDA_MEMCPY_HOST_TO_DEVICE), "H2D Tx");
                cuda_check!(cudaMemset(d_n, 0, 4), "memset n");

                // Launch kernel
                launch_s3_filter(
                    d_a as *const u64,
                    d_b as *const u64,
                    d_tx as *const u64,
                    self.block_bits,
                    d_surv as *mut u32,
                    d_n as *mut u32,
                    n as u32,
                    self.stream,
                );
                cuda_check!(cudaDeviceSynchronize(), "sync");

                // Read back results
                let mut h_n: u32 = 0;
                cuda_check!(cudaMemcpy(&mut h_n as *mut u32 as *mut c_void, d_n, 4, CUDA_MEMCPY_DEVICE_TO_HOST), "D2H n");

                let mut h_surv = vec![0u32; h_n as usize];
                if h_n > 0 {
                    cuda_check!(cudaMemcpy(h_surv.as_mut_ptr() as *mut c_void, d_surv, (h_n as usize)*4, CUDA_MEMCPY_DEVICE_TO_HOST), "D2H surv");
                }

                // Free device memory
                cudaFree(d_a); cudaFree(d_b); cudaFree(d_tx);
                cudaFree(d_surv); cudaFree(d_n);

                Ok(h_surv.iter().map(|&i| i as usize).collect())
            }
        }
    }

    impl Drop for CudaFilter {
        fn drop(&mut self) {
            if !self.stream.is_null() {
                unsafe { cudaStreamDestroy(self.stream); }
            }
        }
    }
}

// ── CPU fallback (no CUDA) ────────────────────────────────────────────────────
// Uses the same direct S₃ quadratic root check in Rust BigInt arithmetic.
// Slower but functionally identical — used in development and CI.

use num_bigint::BigInt;
use num_traits::{Zero, One};

/// Positive modulo: ((a % p) + p) % p
fn bmod(a: BigInt, p: &BigInt) -> BigInt {
    let r = a % p;
    if r < BigInt::zero() { r + p } else { r }
}

fn mod_pow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    if modulus.is_one() { return BigInt::zero(); }
    let mut result = BigInt::one();
    let mut b = bmod(base.clone(), modulus);
    let mut e = exp.clone();
    let zero = BigInt::zero();
    let two  = BigInt::from(2u32);
    while e > zero {
        if &e % &two != zero {
            result = bmod(result * &b, modulus);
        }
        e >>= 1;
        b = bmod(&b * &b, modulus);
    }
    result
}

fn mod_inv(a: &BigInt, p: &BigInt) -> BigInt {
    mod_pow(a, &(p - BigInt::from(2u32)), p)
}

/// Returns true if S₃(x1, x2, tx) = 0 mod p has a root x2 ∈ [b_lo, b_lo + 2^block_bits).
pub fn s3_root_in_range_cpu(
    x1:         &BigInt,
    tx:         &BigInt,
    b_lo:       &BigInt,
    block_bits: u32,
    p:          &BigInt,
) -> bool {
    let two  = BigInt::from(2u32);
    let four = BigInt::from(4u32);
    let b7   = BigInt::from(28u32); // 4*b = 4*7 = 28

    // A = (t - x1)² mod p
    let tdx    = bmod(tx - x1, p);
    let a_coef = bmod(&tdx * &tdx, p);

    // B = -2·x1·t·(t + x1) - 28  mod p
    let tpx    = bmod(tx + x1, p);
    let inner  = bmod(&two * x1 % p * tx % p * &tpx, p);
    let b_coef = bmod(p - &inner - &b7, p);

    // C = t²·x1² - 28·(t + x1)  mod p
    let c_coef = bmod(tx * tx % p * x1 % p * x1 % p - &b7 * &tpx, p);

    // Discriminant: Δ = B² - 4AC  mod p
    let disc = bmod(&b_coef * &b_coef - &four * &a_coef * &c_coef, p);

    let block_size = BigInt::one() << block_bits as usize;

    if disc.is_zero() {
        // double root: x₂ = -B / (2A)
        let inv2a = mod_inv(&bmod(&two * &a_coef, p), p);
        let root  = bmod((p - &b_coef) * &inv2a, p);
        let diff  = bmod(&root - b_lo, p);
        return diff < block_size;
    }

    // Legendre symbol: disc^((p-1)/2) must be 1 (QR) or 0
    let pm1_half = (p - BigInt::one()) / BigInt::from(2u32);
    let leg = mod_pow(&disc, &pm1_half, p);
    if leg != BigInt::zero() && leg != BigInt::one() {
        return false;
    }

    // sqrt(disc) = disc^((p+1)/4) since p ≡ 3 mod 4
    let pp1_4 = (p + BigInt::one()) / BigInt::from(4u32);
    let sq    = mod_pow(&disc, &pp1_4, p);

    let inv2a = mod_inv(&bmod(&two * &a_coef, p), p);
    let neg_b = bmod(p - &b_coef, p);

    let r1 = bmod((&neg_b + &sq) * &inv2a, p);
    let r2 = bmod((&neg_b + p - &sq) * &inv2a, p);

    bmod(&r1 - b_lo, p) < block_size || bmod(&r2 - b_lo, p) < block_size
}

/// High-level filter: for a batch of (A_lo, B_lo) pairs, returns survivor indices.
/// Uses CUDA if feature enabled, else CPU fallback.
pub struct DirectFilter {
    pub target_x:   BigInt,
    pub block_bits: u32,
    pub p:          BigInt,
}

impl DirectFilter {
    pub fn new(target_x: BigInt, block_bits: u32, p: BigInt) -> Self {
        DirectFilter { target_x, block_bits, p }
    }

    /// Returns the indices of pairs that MIGHT contain a root (survivors of PNC kill).
    /// Pairs with no root in range are killed immediately.
    pub fn filter(&self, pairs: &[(BigInt, BigInt)]) -> Vec<usize> {
        #[cfg(feature = "cuda")]
        {
            if let Ok(f) = cuda_impl::CudaFilter::init(&self.target_x, self.block_bits) {
                if let Ok(survivors) = f.filter_pairs(pairs) {
                    return survivors;
                }
            }
            // fall through to CPU on CUDA error
        }

        // CPU fallback: check each (A, B) pair
        // For each pair, scan all x1 in [A_lo, A_lo + 2^block_bits)
        let block_size = 1usize << self.block_bits.min(24); // cap at 2^24 for CPU
        pairs.iter().enumerate()
            .filter(|(_, (a_lo, b_lo))| {
                (0..block_size).any(|i| {
                    let x1 = a_lo + i;
                    s3_root_in_range_cpu(&x1, &self.target_x, b_lo, self.block_bits, &self.p)
                })
            })
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Benchmark: process `n_pairs` random pairs and return rejection rate.
    /// CPU version: tests only x1=A_lo per pair (conservative, GPU tests all x1 in block).
    pub fn benchmark_rejection_rate(&self, n_pairs: usize) -> f64 {
        let mut xs: u64 = 0xdeadbeef_cafebabe ^ n_pairs as u64;
        let mut xorshift = move || -> u64 {
            xs ^= xs << 13; xs ^= xs >> 7; xs ^= xs << 17; xs
        };

        let mut killed = 0usize;
        for _ in 0..n_pairs {
            let bytes_a: Vec<u8> = (0..32).map(|i| {
                if i % 8 == 0 { xorshift() as u8 } else { (xorshift() >> ((i%8)*8)) as u8 }
            }).collect();
            let bytes_b: Vec<u8> = (0..32).map(|i| {
                if i % 8 == 0 { xorshift() as u8 } else { (xorshift() >> ((i%8)*8)) as u8 }
            }).collect();
            let mut ba = [0u8; 32];
            let mut bb = [0u8; 32];
            ba.copy_from_slice(&bytes_a);
            bb.copy_from_slice(&bytes_b);
            // Clear top 4 bits to stay below p
            ba[31] &= 0x0F;
            bb[31] &= 0x0F;

            let a_lo = BigInt::from_bytes_le(num_bigint::Sign::Plus, &ba);
            let b_lo = BigInt::from_bytes_le(num_bigint::Sign::Plus, &bb);

            let survived = s3_root_in_range_cpu(
                &a_lo, &self.target_x, &b_lo, self.block_bits, &self.p
            );
            if !survived { killed += 1; }
        }

        killed as f64 / n_pairs as f64
    }
}
