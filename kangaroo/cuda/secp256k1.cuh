// secp256k1 device arithmetic for CUDA Kangaroo
// Field elements: uint64_t[4], little-endian (limb 0 = bits 0..63)
// Uses __uint128_t for intermediate products (supported in CUDA device code)

#pragma once
#include <stdint.h>

typedef unsigned long long  u64;
typedef unsigned __int128   u128;
typedef unsigned int        u32;

// ─── secp256k1 constants ──────────────────────────────────────────────────────

// p = 2^256 − 2^32 − 977
#define P0  0xFFFFFFFEFFFFFC2FULL
#define P1  0xFFFFFFFFFFFFFFFFULL
#define P2  0xFFFFFFFFFFFFFFFFULL
#define P3  0xFFFFFFFFFFFFFFFFULL

// β: primitive cube root of 1 mod p  →  ψ(x,y) = (β·x, y)
#define BETA0  0xC1396C28719501EEULL
#define BETA1  0x9CF0497512F58995ULL
#define BETA2_  0x6E64479EAC3434E9ULL   // trailing _ to avoid collision with BETA2_ prefix macros
#define BETA3  0x7AE96A2B657C0710ULL

// β² = p − 1 − β
#define BETA2_0  0x3EC693D68E6AFA40ULL
#define BETA2_1  0x630FB68AED0A766AULL
#define BETA2_2  0x919BB86153CBCB16ULL
#define BETA2_3  0x851695D49A83F8EFULL

// n = group order
#define N0  0xBFD25E8CD0364141ULL
#define N1  0xBAAEDCE6AF48A03BULL
#define N2  0xFFFFFFFFFFFFFFFEULL
#define N3  0xFFFFFFFFFFFFFFFFULL

// λ = GLV eigenvalue (ψ(P) = λ·P in scalar field)
#define LAM0  0xDF02967C1B23BD72ULL
#define LAM1  0x122E22EA20816678ULL
#define LAM2  0xA5261C028812645AULL
#define LAM3  0x5363AD4CC05C30E0ULL

// ─── Comparison ───────────────────────────────────────────────────────────────

__device__ __forceinline__
bool fe_lt(const u64 a[4], const u64 b[4]) {
    for (int i = 3; i >= 0; i--) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}

__device__ __forceinline__
bool fe_eq(const u64 a[4], const u64 b[4]) {
    return a[0]==b[0] && a[1]==b[1] && a[2]==b[2] && a[3]==b[3];
}

// ─── Field mod p ──────────────────────────────────────────────────────────────

__device__ __forceinline__
void fp_sub_p_inplace(u64 r[4]) {
    const u64 p[4] = {P0, P1, P2, P3};
    u128 s; u64 borrow = 0;
    for (int i = 0; i < 4; i++) {
        s = (u128)r[i] - p[i] - borrow;
        r[i] = (u64)s;
        borrow = (s >> 127) & 1;
    }
}

__device__ __forceinline__
void fp_add(const u64 a[4], const u64 b[4], u64 r[4]) {
    u128 s; u64 carry = 0;
    for (int i = 0; i < 4; i++) {
        s = (u128)a[i] + b[i] + carry;
        r[i] = (u64)s; carry = (u64)(s >> 64);
    }
    const u64 p[4] = {P0, P1, P2, P3};
    if (carry || !fe_lt(r, p)) fp_sub_p_inplace(r);
}

__device__ __forceinline__
void fp_sub(const u64 a[4], const u64 b[4], u64 r[4]) {
    if (!fe_lt(a, b)) {
        u128 s; u64 borrow = 0;
        for (int i = 0; i < 4; i++) {
            s = (u128)a[i] - b[i] - borrow;
            r[i] = (u64)s; borrow = (s >> 127) & 1;
        }
    } else {
        u64 tmp[4];
        u128 s; u64 borrow = 0;
        for (int i = 0; i < 4; i++) {
            s = (u128)b[i] - a[i] - borrow;
            tmp[i] = (u64)s; borrow = (s >> 127) & 1;
        }
        const u64 p[4] = {P0, P1, P2, P3};
        borrow = 0;
        for (int i = 0; i < 4; i++) {
            s = (u128)p[i] - tmp[i] - borrow;
            r[i] = (u64)s; borrow = (s >> 127) & 1;
        }
    }
}

__device__ __forceinline__
void fp_neg(const u64 a[4], u64 r[4]) {
    if ((a[0]|a[1]|a[2]|a[3]) == 0) { r[0]=r[1]=r[2]=r[3]=0; return; }
    const u64 p[4] = {P0, P1, P2, P3};
    fp_sub(p, a, r);
}

// 256×256 → 512-bit product
__device__ __forceinline__
void mul512(const u64 a[4], const u64 b[4], u64 t[8]) {
    for (int i = 0; i < 8; i++) t[i] = 0;
    for (int i = 0; i < 4; i++) {
        u64 carry = 0;
        for (int j = 0; j < 4; j++) {
            u128 prod = (u128)a[i] * b[j] + t[i+j] + carry;
            t[i+j] = (u64)prod; carry = (u64)(prod >> 64);
        }
        t[i+4] += carry;
    }
}

// Montgomery-style reduction mod p = 2^256 − 2^32 − 977
__device__ __forceinline__
void reduce512(u64 t[8], u64 r[4]) {
    const u64 lo[4] = {t[0], t[1], t[2], t[3]};
    const u64 hi[4] = {t[4], t[5], t[6], t[7]};

    u64 c[5] = {};
    u128 s; u64 carry = 0;

    for (int i = 0; i < 4; i++) {
        s = (u128)hi[i] * 977ULL + c[i] + carry;
        c[i] = (u64)s; carry = (u64)(s >> 64);
    }
    c[4] = carry; carry = 0;

    s = (u128)c[0] + (hi[0] << 32);
    c[0] = (u64)s; carry = (u64)(s >> 64);
    for (int i = 1; i < 4; i++) {
        s = (u128)c[i] + (hi[i] << 32) + (hi[i-1] >> 32) + carry;
        c[i] = (u64)s; carry = (u64)(s >> 64);
    }
    c[4] += (hi[3] >> 32) + carry;

    u64 r5[5] = {}; carry = 0;
    for (int i = 0; i < 4; i++) {
        s = (u128)lo[i] + c[i] + carry;
        r5[i] = (u64)s; carry = (u64)(s >> 64);
    }
    r5[4] = c[4] + carry;

    if (r5[4] > 0) {
        u64 e = r5[4];
        u128 ex = (u128)e * (977ULL + ((u128)1ULL << 32));
        s = (u128)r5[0] + (u64)ex;
        r5[0] = (u64)s; carry = (u64)(s >> 64);
        s = (u128)r5[1] + (u64)(ex >> 64) + carry;
        r5[1] = (u64)s; carry = (u64)(s >> 64);
        s = (u128)r5[2] + carry;
        r5[2] = (u64)s; carry = (u64)(s >> 64);
        r5[3] += carry; r5[4] = 0;
    }

    for (int i = 0; i < 4; i++) r[i] = r5[i];
    const u64 p[4] = {P0, P1, P2, P3};
    if (!fe_lt(r, p)) fp_sub_p_inplace(r);
}

__device__ __forceinline__
void fp_mul(const u64 a[4], const u64 b[4], u64 r[4]) {
    u64 t[8]; mul512(a, b, t); reduce512(t, r);
}

__device__ __forceinline__
void fp_sqr(const u64 a[4], u64 r[4]) { fp_mul(a, a, r); }

// a^e mod p — constant-time binary exponentiation
__device__
void fp_pow(const u64 a[4], const u64 e[4], u64 r[4]) {
    u64 base[4] = {a[0],a[1],a[2],a[3]};
    r[0]=1; r[1]=r[2]=r[3]=0;
    u64 exp[4] = {e[0],e[1],e[2],e[3]};
    while (exp[0]|exp[1]|exp[2]|exp[3]) {
        if (exp[0] & 1) { u64 tmp[4]; fp_mul(r, base, tmp); for(int i=0;i<4;i++) r[i]=tmp[i]; }
        for (int i=0;i<3;i++) exp[i]=(exp[i]>>1)|(exp[i+1]<<63); exp[3]>>=1;
        u64 tmp[4]; fp_sqr(base, tmp); for(int i=0;i<4;i++) base[i]=tmp[i];
    }
}

__device__ __forceinline__
void fp_inv(const u64 a[4], u64 r[4]) {
    // p − 2 = 0xFFFFFFFEFFFFFC2D || 0xFF...FF × 3
    const u64 pm2[4] = {0xFFFFFFFEFFFFFC2DULL, P1, P2, P3};
    fp_pow(a, pm2, r);
}

// ─── GLV endomorphisms (affine, no inversion) ────────────────────────────────

__device__ __forceinline__
void psi_x(const u64 ax[4], u64 r[4]) {
    const u64 beta[4] = {BETA0, BETA1, BETA2_, BETA3};
    fp_mul(beta, ax, r);
}

__device__ __forceinline__
void psi2_x(const u64 ax[4], u64 r[4]) {
    const u64 beta2[4] = {BETA2_0, BETA2_1, BETA2_2, BETA2_3};
    fp_mul(beta2, ax, r);
}

// canonical affine x = min(ax, β·ax, β²·ax)
// Takes AFFINE x only — no Jacobian Z parameter
__device__ __forceinline__
void canonical_x_affine(const u64 ax[4], u64 r[4]) {
    u64 x1[4], x2[4];
    psi_x(ax, x1);
    psi2_x(ax, x2);
    for (int i = 0; i < 4; i++) r[i] = ax[i];
    if (fe_lt(x1, r)) for (int i = 0; i < 4; i++) r[i] = x1[i];
    if (fe_lt(x2, r)) for (int i = 0; i < 4; i++) r[i] = x2[i];
}

// ─── Jacobian → affine normalization ─────────────────────────────────────────

// Converts (X:Y:Z) Jacobian to affine (ax, ay).
// Cost: 1 fp_inv + 4 fp_mul + 1 fp_sqr  (amortize with NORM_INTERVAL)
__device__
void pt_normalize(const u64 X[4], const u64 Y[4], const u64 Z[4],
                  u64 ax[4], u64 ay[4]) {
    u64 zinv[4], z2[4], z3[4];
    fp_inv(Z, zinv);
    fp_sqr(zinv, z2);
    fp_mul(zinv, z2, z3);
    fp_mul(X, z2, ax);
    fp_mul(Y, z3, ay);
}

// ─── Jacobian mixed add: Jacobian P + Affine Q → Jacobian R ──────────────────
// madd-2007-bl (EFD), cost 7M + 4S, no inversion

__device__
void pt_add_mixed(
    const u64 px[4], const u64 py[4], const u64 pz[4],
    const u64 qx[4], const u64 qy[4],
    u64 rx[4], u64 ry[4], u64 rz[4]
) {
    u64 z1z1[4], u2[4], s2[4], h[4], r_[4], hh[4], hhh[4], v[4], tmp[4], tmp2[4];

    fp_sqr(pz,  z1z1);          // Z1Z1 = Z1²
    fp_mul(qx,  z1z1, u2);      // U2   = X2·Z1Z1
    fp_mul(pz,  z1z1, tmp);     // Z1³
    fp_mul(qy,  tmp,  s2);      // S2   = Y2·Z1³
    fp_sub(u2,  px,   h);       // H    = U2 − X1
    fp_sub(s2,  py,   r_);      // R    = S2 − Y1
    fp_sqr(h,   hh);            // HH   = H²
    fp_mul(h,   hh,   hhh);     // HHH  = H·HH
    fp_mul(px,  hh,   v);       // V    = X1·HH
    fp_sqr(r_,  tmp);           // R²
    fp_sub(tmp, hhh,  tmp2);    // R² − HHH
    fp_add(v,   v,    tmp);     // 2V
    fp_sub(tmp2,tmp,  rx);      // X3   = R² − HHH − 2V
    fp_sub(v,   rx,   tmp);     // V − X3
    fp_mul(r_,  tmp,  tmp2);    // R·(V−X3)
    fp_mul(py,  hhh,  tmp);     // Y1·HHH
    fp_sub(tmp2,tmp,  ry);      // Y3   = R·(V−X3) − Y1·HHH
    fp_mul(h,   pz,   rz);      // Z3   = H·Z1
}

// ─── Scalar arithmetic mod n ──────────────────────────────────────────────────

__device__ __forceinline__
void sc_add(const u64 a[4], const u64 b[4], u64 r[4]) {
    u128 s; u64 carry = 0;
    for (int i = 0; i < 4; i++) {
        s = (u128)a[i] + b[i] + carry;
        r[i] = (u64)s; carry = (u64)(s >> 64);
    }
    const u64 n[4] = {N0, N1, N2, N3};
    if (carry || !fe_lt(r, n)) {
        u64 borrow = 0;
        for (int i = 0; i < 4; i++) {
            s = (u128)r[i] - n[i] - borrow;
            r[i] = (u64)s; borrow = (s >> 127) & 1;
        }
    }
}

// ─── Structures ───────────────────────────────────────────────────────────────

struct JumpPoint {
    u64 x[4];   // affine
    u64 y[4];   // affine
    u64 s[4];   // scalar mod n
};

struct Animal {
    u64 x[4], y[4], z[4];  // Jacobian
    u64 scalar[4];           // accumulated jump scalar mod n
    u32 is_wild;
    u32 pad[3];
};

// DPEntry carries the NORMALIZED canonical x (no Z conversion needed on host)
struct DPEntry {
    u64 canon_x[4];  // exact affine canonical x = min(x, β·x, β²·x)
    u64 scalar[4];   // accumulated scalar mod n at the DP
    u32 is_wild;
    u32 pad[3];
};
