/*
 * sinGRAAL — PNC Direct Filter CUDA Kernel
 *
 * Principe : S₃(x₁, x₂, t) est quadratique en x₂ (coefficients en x₁).
 * Pour chaque x₁ candidat on résout directement les racines x₂ mod p
 * et on vérifie si l'une tombe dans [B_lo, B_lo + 2^block_bits).
 *
 * Pas de LLL ici — la structure quadratique de S₃ permet un kill-switch
 * exactement en O(fp_mul) par (x₁, B-bloc).
 *
 * Architecture :
 *   Un thread CUDA = un x₁ dans [A_lo, A_lo + 2^block_bits)
 *   Warp-ballot : un seul survivor par (A-bloc, B-bloc) si racine trouvée
 *
 * secp256k1 : p = 2^256 − 2^32 − 977,  p ≡ 3 mod 4
 *   → sqrt(a) = a^((p+1)/4) mod p  (pas de Tonelli-Shanks)
 *   → (p+1)/4 = (2^256 − 2^32 − 976) / 4  (voir constante EXP_SQRT)
 *
 * Mémoire par thread :
 *   ~20 Fe (registres 256-bit) = ~160 × 64-bit regs
 *   Bien en-dessous du max (255 regs/thread)
 */

#include <stdint.h>
#include <cuda_runtime.h>

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef unsigned __int128 u128;

/* ── 256-bit field element (little-endian 64-bit limbs) ──────────────────── */
struct Fe { u64 d[4]; };

/* ── secp256k1 prime p ────────────────────────────────────────────────────── */
__device__ __constant__ Fe FE_P = {{
    0xFFFFFFFEFFFFFC2FULL,
    0xFFFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFFULL
}};

/* (p+1)/4 for sqrt exponent — p ≡ 3 mod 4 */
__device__ __constant__ Fe EXP_SQRT = {{
    0xBFFFFFFFBFFFFF0BULL,  /* (p+1)/4 limb 0 */
    0xFFFFFFFFFFFFFFFFULL,
    0xFFFFFFFFFFFFFFFFULL,
    0x3FFFFFFFFFFFFFFFULL
}};

__device__ __constant__ Fe FE_ZERO = {{0,0,0,0}};
__device__ __constant__ Fe FE_ONE  = {{1,0,0,0}};
__device__ __constant__ u64 FE_B   = 7ULL; /* secp256k1 b=7 */

/* ── Basic carry-propagating add / sub ────────────────────────────────────── */

__device__ inline int fe_is_zero(Fe a) {
    return (a.d[0]|a.d[1]|a.d[2]|a.d[3]) == 0;
}

__device__ inline int fe_lt(Fe a, Fe b) {
    for (int i = 3; i >= 0; i--) {
        if (a.d[i] < b.d[i]) return 1;
        if (a.d[i] > b.d[i]) return 0;
    }
    return 0; /* equal */
}

/* Lazy reduce: subtract p if a >= p */
__device__ Fe fe_reduce(Fe a, u64 carry) {
    /* carry is 0 or 1 (overflow word above d[3]) */
    if (carry == 0 && !fe_lt(FE_P, a) && !(a.d[0]==FE_P.d[0] && a.d[1]==FE_P.d[1] &&
        a.d[2]==FE_P.d[2] && a.d[3]==FE_P.d[3])) {
        /* a < p already (or equal — but p is prime so rarely == p) */
        if (fe_lt(a, FE_P)) return a;
    }
    /* subtract p */
    Fe r;
    u128 t = (u128)a.d[0] + 0x10000003D1ULL + ((u128)carry << 64);
    /* p = 2^256 - 2^32 - 977  →  -p mod 2^256 = 2^32 + 977 = 0x10000003D1 */
    /* a - p = a + (2^256 - p) mod 2^256 = a + (2^32+977) mod 2^256 */
    /* reuse simpler path: a >= p → a - p */
    t = (u128)a.d[0] - FE_P.d[0] + (carry ? (u128)1 << 64 : 0);
    r.d[0] = (u64)t;
    t = (u128)a.d[1] - FE_P.d[1] + (u64)(t >> 64);  /* borrow */
    /* This is signed borrow — use __int128 */
    /* Cleaner approach: */
    __int128 s = (__int128)a.d[0] - (__int128)FE_P.d[0];
    r.d[0] = (u64)s; s >>= 64;
    s += (__int128)a.d[1] - (__int128)FE_P.d[1]; r.d[1] = (u64)s; s >>= 64;
    s += (__int128)a.d[2] - (__int128)FE_P.d[2]; r.d[2] = (u64)s; s >>= 64;
    s += (__int128)a.d[3] - (__int128)FE_P.d[3] + (__int128)carry;
    r.d[3] = (u64)s;
    return r;
}

__device__ Fe fe_add(Fe a, Fe b) {
    Fe r;
    u128 t = (u128)a.d[0] + b.d[0];  r.d[0] = (u64)t;
    t = (t>>64) + a.d[1] + b.d[1];   r.d[1] = (u64)t;
    t = (t>>64) + a.d[2] + b.d[2];   r.d[2] = (u64)t;
    t = (t>>64) + a.d[3] + b.d[3];   r.d[3] = (u64)t;
    u64 carry = (u64)(t >> 64);
    /* reduce if carry or a+b >= p */
    if (carry || !fe_lt(r, FE_P)) {
        __int128 s = (__int128)r.d[0] - (__int128)FE_P.d[0];
        r.d[0] = (u64)s; s >>= 64;
        s += (__int128)r.d[1] - (__int128)FE_P.d[1]; r.d[1] = (u64)s; s >>= 64;
        s += (__int128)r.d[2] - (__int128)FE_P.d[2]; r.d[2] = (u64)s; s >>= 64;
        s += (__int128)r.d[3] - (__int128)FE_P.d[3] + (__int128)carry;
        r.d[3] = (u64)s;
    }
    return r;
}

__device__ Fe fe_sub(Fe a, Fe b) {
    Fe r;
    __int128 s = (__int128)a.d[0] - b.d[0];
    r.d[0] = (u64)s; s >>= 64;
    s += (__int128)a.d[1] - b.d[1]; r.d[1] = (u64)s; s >>= 64;
    s += (__int128)a.d[2] - b.d[2]; r.d[2] = (u64)s; s >>= 64;
    s += (__int128)a.d[3] - b.d[3]; r.d[3] = (u64)s;
    /* if borrow (s negative → top bit set), add p */
    if ((int64_t)r.d[3] < 0 || (s >> 64) < 0) {
        u128 t = (u128)r.d[0] + FE_P.d[0]; r.d[0] = (u64)t;
        t = (t>>64) + r.d[1] + FE_P.d[1];  r.d[1] = (u64)t;
        t = (t>>64) + r.d[2] + FE_P.d[2];  r.d[2] = (u64)t;
        t = (t>>64) + r.d[3] + FE_P.d[3];  r.d[3] = (u64)t;
    }
    return r;
}

__device__ Fe fe_neg(Fe a) { return fe_sub(FE_ZERO, a); }

/* 256×256 → 256 multiplication mod p (schoolbook + Mersenne reduction) */
__device__ Fe fe_mul(Fe a, Fe b) {
    /* Step 1: compute 512-bit product r[0..7] */
    u64 r[8] = {0};
    u128 t;
    #pragma unroll
    for (int i = 0; i < 4; i++) {
        u128 carry = 0;
        for (int j = 0; j < 4; j++) {
            t = (u128)a.d[i] * b.d[j] + r[i+j] + carry;
            r[i+j] = (u64)t;
            carry = t >> 64;
        }
        r[i+4] += (u64)carry;
    }
    /*
     * Step 2: reduce mod p = 2^256 - 2^32 - 977
     * r[4..7] × 2^256 ≡ r[4..7] × (2^32 + 977) mod p
     *
     * hi = r[4..7] (256-bit), lo = r[0..3] (256-bit)
     * result = lo + hi × (2^32 + 977)
     *        = lo + hi × 2^32 + hi × 977
     */
    u64 hi[4] = {r[4], r[5], r[6], r[7]};
    Fe lo = {{r[0], r[1], r[2], r[3]}};

    /* hi × 977 */
    u64 c977[4];
    {
        u128 carry = 0;
        for (int i = 0; i < 4; i++) {
            t = (u128)hi[i] * 977ULL + carry;
            c977[i] = (u64)t;
            carry = t >> 64;
        }
        /* carry is overflow — handle below via second reduction */
    }

    /* hi × 2^32 (left-shift by 32 bits) */
    u64 c32[4];
    c32[0] = hi[0] << 32;
    c32[1] = (hi[0] >> 32) | (hi[1] << 32);
    c32[2] = (hi[1] >> 32) | (hi[2] << 32);
    c32[3] = (hi[2] >> 32) | (hi[3] << 32);
    u64 overflow_32 = hi[3] >> 32;

    /* Add all contributions to lo */
    Fe mid;
    {
        u128 carry = 0;
        t = (u128)lo.d[0] + c977[0] + c32[0]; mid.d[0] = (u64)t; carry = t >> 64;
        t = (u128)lo.d[1] + c977[1] + c32[1] + carry; mid.d[1] = (u64)t; carry = t >> 64;
        t = (u128)lo.d[2] + c977[2] + c32[2] + carry; mid.d[2] = (u64)t; carry = t >> 64;
        t = (u128)lo.d[3] + c977[3] + c32[3] + carry; mid.d[3] = (u64)t; carry = t >> 64;
        /* carry + overflow_32 × 2^32 may still overflow — second pass */
        u64 extra = (u64)carry + overflow_32;
        if (extra) {
            /* extra × (2^32 + 977) is tiny (< 2^33), add to low limbs */
            t = (u128)mid.d[0] + extra * 977ULL; mid.d[0] = (u64)t; carry = t >> 64;
            t = (u128)mid.d[1] + (extra << 32) + carry; mid.d[1] = (u64)t; carry = t >> 64;
            t = (u128)mid.d[2] + carry; mid.d[2] = (u64)t; carry = t >> 64;
            t = (u128)mid.d[3] + carry; mid.d[3] = (u64)t;
        }
    }

    /* Final conditional subtraction */
    if (!fe_lt(mid, FE_P)) {
        Fe res;
        __int128 s = (__int128)mid.d[0] - (__int128)FE_P.d[0];
        res.d[0] = (u64)s; s >>= 64;
        s += (__int128)mid.d[1] - FE_P.d[1]; res.d[1] = (u64)s; s >>= 64;
        s += (__int128)mid.d[2] - FE_P.d[2]; res.d[2] = (u64)s; s >>= 64;
        s += (__int128)mid.d[3] - FE_P.d[3]; res.d[3] = (u64)s;
        return res;
    }
    return mid;
}

__device__ Fe fe_sqr(Fe a) { return fe_mul(a, a); }

/* Scalar multiply: a × c (64-bit scalar) */
__device__ Fe fe_mul_small(Fe a, u64 c) {
    Fe r;
    u128 carry = 0;
    for (int i = 0; i < 4; i++) {
        u128 t = (u128)a.d[i] * c + carry;
        r.d[i] = (u64)t;
        carry = t >> 64;
    }
    /* carry × 2^256 ≡ carry × (2^32 + 977) mod p */
    u64 ex = (u64)carry;
    if (ex) {
        u128 t = (u128)r.d[0] + ex * 977ULL; r.d[0] = (u64)t; u64 c2 = (u64)(t>>64);
        t = (u128)r.d[1] + (ex << 32) + c2; r.d[1] = (u64)t; c2 = (u64)(t>>64);
        t = (u128)r.d[2] + c2; r.d[2] = (u64)t; c2 = (u64)(t>>64);
        t = (u128)r.d[3] + c2; r.d[3] = (u64)t;
    }
    if (!fe_lt(r, FE_P)) {
        __int128 s = (__int128)r.d[0] - FE_P.d[0];
        r.d[0] = (u64)s; s >>= 64;
        s += (__int128)r.d[1] - FE_P.d[1]; r.d[1] = (u64)s; s >>= 64;
        s += (__int128)r.d[2] - FE_P.d[2]; r.d[2] = (u64)s; s >>= 64;
        s += (__int128)r.d[3] - FE_P.d[3]; r.d[3] = (u64)s;
    }
    return r;
}

/*
 * Square root mod p — fast because p ≡ 3 (mod 4):
 *   sqrt(a) = a^((p+1)/4) mod p
 *
 * Returns sqrt if a is a QR, else returns garbage.
 * Caller must verify: sqr(result) == a.
 */
__device__ Fe fe_sqrt(Fe a) {
    /* binary exponentiation by EXP_SQRT = (p+1)/4 */
    Fe result = FE_ONE;
    Fe base   = a;
    u64 exp[4]; /* copy from constant */
    exp[0] = EXP_SQRT.d[0];
    exp[1] = EXP_SQRT.d[1];
    exp[2] = EXP_SQRT.d[2];
    exp[3] = EXP_SQRT.d[3];
    for (int w = 0; w < 4; w++) {
        u64 e = exp[w];
        for (int b = 0; b < 64; b++) {
            if (e & 1ULL) result = fe_mul(result, base);
            base = fe_sqr(base);
            e >>= 1;
        }
    }
    return result;
}

/* Modular inverse via Fermat: a^(p-2) mod p */
__device__ Fe fe_inv(Fe a) {
    /* p - 2 = FE_P - 2 */
    Fe exp = FE_P;
    /* exp = p - 2 */
    exp.d[0] -= 2ULL;  /* p.d[0] >= 2 so no borrow needed */
    Fe result = FE_ONE;
    Fe base = a;
    for (int w = 0; w < 4; w++) {
        u64 e = exp.d[w];
        for (int b = 0; b < 64; b++) {
            if (e & 1ULL) result = fe_mul(result, base);
            base = fe_sqr(base);
            e >>= 1;
        }
    }
    return result;
}

/*
 * ── S₃ Semaev Polynomial (secp256k1, b=7) ────────────────────────────────
 *
 * S₃(x₁, x₂, t) as quadratic in x₂:
 *   A_coef = (t - x₁)²
 *   B_coef = -2·x₁·t·(t + x₁) - 28
 *   C_coef = t²·x₁² - 28·(t + x₁)
 *
 *   S₃ = A_coef·x₂² + B_coef·x₂ + C_coef  (all mod p)
 *
 * Roots: x₂ = (-B ± sqrt(B²-4AC)) / (2A)
 *
 * Returns: 1 if root x₂ ∈ [B_lo, B_lo + 2^block_bits), 0 otherwise.
 * Sets *root_out to the root if found (optional, may be NULL).
 */
__device__ int s3_root_in_range(
    Fe x1,        /* x₁ candidate */
    Fe t,         /* target x-coordinate T.x */
    Fe B_lo,      /* lower bound of x₂ range */
    u32 block_bits, /* size = 2^block_bits */
    Fe* root_out  /* optional output */
) {
    Fe two = {{2,0,0,0}};
    Fe four= {{4,0,0,0}};

    /* A = (t - x₁)² */
    Fe tdx = fe_sub(t, x1);
    Fe A = fe_sqr(tdx);

    /* B = -2·x₁·t·(t + x₁) - 28 */
    Fe tpx = fe_add(t, x1);
    Fe xt  = fe_mul(x1, t);
    Fe B   = fe_mul(fe_mul_small(xt, 2), tpx);
    B = fe_neg(B);
    B = fe_sub(B, (Fe){{28,0,0,0}});

    /* C = t²·x₁² - 28·(t + x₁) */
    Fe t2 = fe_sqr(t);
    Fe x2 = fe_sqr(x1);
    Fe C  = fe_sub(fe_mul(t2, x2), fe_mul_small(tpx, 28));

    /* Discriminant: Δ = B² - 4AC */
    Fe disc = fe_sub(fe_sqr(B), fe_mul(fe_mul_small(A, 4), C));

    /* Check if Δ is a quadratic residue: Δ^((p-1)/2) must be 1 or 0 */
    if (fe_is_zero(disc)) {
        /* double root: x₂ = -B / (2A) */
        Fe inv2A = fe_inv(fe_mul_small(A, 2));
        Fe root  = fe_mul(fe_neg(B), inv2A);
        /* range check: root - B_lo < 2^block_bits (unsigned comparison) */
        Fe diff = fe_sub(root, B_lo);
        int in_range = (block_bits >= 64)
            ? (diff.d[1]==0 && diff.d[2]==0 && diff.d[3]==0 && diff.d[0] < (1ULL<<(block_bits-0)))
            : (diff.d[1]==0 && diff.d[2]==0 && diff.d[3]==0 && diff.d[0] < (1ULL<<block_bits));
        if (in_range && root_out) *root_out = root;
        return in_range;
    }

    /* Legendre test: disc^((p-1)/2) */
    {
        /* (p-1)/2 = (FE_P - 1) / 2  →  shift FE_P right by 1 */
        u64 exp_leg[4];
        exp_leg[0] = (FE_P.d[0]>>1) | (FE_P.d[1]<<63);
        exp_leg[1] = (FE_P.d[1]>>1) | (FE_P.d[2]<<63);
        exp_leg[2] = (FE_P.d[2]>>1) | (FE_P.d[3]<<63);
        exp_leg[3] = FE_P.d[3]>>1;
        /* Since p-1 is even, this gives (p-1)/2 */
        Fe leg = FE_ONE;
        Fe base = disc;
        for (int w = 0; w < 4; w++) {
            u64 e = exp_leg[w];
            for (int b = 0; b < 64; b++) {
                if (e & 1ULL) leg = fe_mul(leg, base);
                base = fe_sqr(base);
                e >>= 1;
            }
        }
        /* leg == 1: QR; leg == p-1: non-residue; leg == 0: disc==0 (handled above) */
        if (!( (leg.d[0]==1 && leg.d[1]==0 && leg.d[2]==0 && leg.d[3]==0) ))
            return 0; /* not a QR → no real root */
    }

    Fe sq  = fe_sqrt(disc);
    Fe inv2A = fe_inv(fe_mul_small(A, 2));

    /* Two roots: r1 = (-B + sq) / 2A,  r2 = (-B - sq) / 2A */
    Fe negB = fe_neg(B);
    Fe r1 = fe_mul(fe_add(negB, sq), inv2A);
    Fe r2 = fe_mul(fe_sub(negB, sq), inv2A);

    /* Range check helper */
    auto in_range_check = [&](Fe root) -> int {
        Fe diff = fe_sub(root, B_lo);
        if (diff.d[2] || diff.d[3]) return 0;
        if (block_bits >= 128) return 1;
        if (block_bits >= 64) {
            if (diff.d[1] != 0) return 0;
            return 1; /* diff.d[0] is any value, block covers 2^64+ */
        }
        return (diff.d[1]==0 && diff.d[0] < (1ULL<<block_bits));
    };

    if (in_range_check(r1)) {
        if (root_out) *root_out = r1;
        return 1;
    }
    if (in_range_check(r2)) {
        if (root_out) *root_out = r2;
        return 1;
    }
    return 0;
}

/* ── Main CUDA Kernel ─────────────────────────────────────────────────────── */

/*
 * s3_direct_filter_kernel
 *
 * Each CUDA thread handles ONE x₁ value = A_lo[pair_idx] + thread_idx_in_block
 * (thread_idx_in_block ∈ [0, block_size) where block_size = 2^block_bits).
 *
 * If any thread in the warp/block finds a root in [B_lo, B_lo+block_size),
 * the (A-bloc, B-bloc) pair is a SURVIVOR.
 *
 * Grid: (n_pairs, ceil(block_size / THREADS_PER_PAIR)) blocks
 * Block: THREADS_PER_PAIR threads (must be ≤ 1024)
 */
#define THREADS_PER_PAIR 256

__global__ void s3_direct_filter_kernel(
    const u64* __restrict__ A_lo,       /* [n_pairs × 4] left block bases */
    const u64* __restrict__ B_lo,       /* [n_pairs × 4] right block bases */
    const u64* __restrict__ T_x,        /* target x-coordinate [4] */
    u32  block_bits,                    /* 2^block_bits = block size */
    u32* __restrict__ survivors,        /* output: survivor pair indices */
    u32* __restrict__ n_survivors,      /* atomic counter */
    u32  n_pairs
) {
    u32 pair_idx   = blockIdx.x;
    u32 thread_off = blockIdx.y * THREADS_PER_PAIR + threadIdx.x;
    u32 block_size = 1u << block_bits;

    if (pair_idx >= n_pairs || thread_off >= block_size) return;

    /* Load block bases into registers */
    Fe A, B, Tx;
    A.d[0] = A_lo[pair_idx*4+0];  A.d[1] = A_lo[pair_idx*4+1];
    A.d[2] = A_lo[pair_idx*4+2];  A.d[3] = A_lo[pair_idx*4+3];
    B.d[0] = B_lo[pair_idx*4+0];  B.d[1] = B_lo[pair_idx*4+1];
    B.d[2] = B_lo[pair_idx*4+2];  B.d[3] = B_lo[pair_idx*4+3];
    Tx.d[0] = T_x[0]; Tx.d[1] = T_x[1];
    Tx.d[2] = T_x[2]; Tx.d[3] = T_x[3];

    /* x₁ = A + thread_off */
    Fe x1 = A;
    u128 t = (u128)x1.d[0] + thread_off;
    x1.d[0] = (u64)t; t >>= 64;
    x1.d[1] += (u64)t;  /* no further carry needed for small thread_off */

    int found = s3_root_in_range(x1, Tx, B, block_bits, NULL);

    /* Warp-level OR reduction */
    u32 ballot = __ballot_sync(0xFFFFFFFF, found);

    /* First thread in warp with a hit records the pair */
    if (ballot && threadIdx.x % 32 == __ffs(ballot) - 1) {
        u32 pos = atomicAdd(n_survivors, 1u);
        survivors[pos] = pair_idx;
    }
}

/*
 * Wrapper for coarse block_bits > THREADS_PER_PAIR:
 * launches multiple Y-grid blocks per pair.
 */
extern "C" void launch_s3_filter(
    const u64* A_lo,
    const u64* B_lo,
    const u64* T_x,
    u32  block_bits,
    u32* survivors,
    u32* n_survivors_dev,
    u32  n_pairs,
    cudaStream_t stream
) {
    u32 block_size   = 1u << block_bits;
    u32 y_blocks     = (block_size + THREADS_PER_PAIR - 1) / THREADS_PER_PAIR;
    dim3 grid(n_pairs, y_blocks);
    dim3 block(THREADS_PER_PAIR);
    s3_direct_filter_kernel<<<grid, block, 0, stream>>>(
        A_lo, B_lo, T_x, block_bits, survivors, n_survivors_dev, n_pairs
    );
}
