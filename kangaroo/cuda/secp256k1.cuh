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
    u64 r0, r1, r2, r3, carry;
    asm("add.cc.u64   %0, %4, %8;\n\t"
        "addc.cc.u64  %1, %5, %9;\n\t"
        "addc.cc.u64  %2, %6, %10;\n\t"
        "addc.cc.u64  %3, %7, %11;\n\t"
        "addc.u64     %4, 0, 0;"
        : "=l"(r0),"=l"(r1),"=l"(r2),"=l"(r3),"=l"(carry)
        : "l"(a[0]),"l"(a[1]),"l"(a[2]),"l"(a[3]),
          "l"(b[0]),"l"(b[1]),"l"(b[2]),"l"(b[3]));
    r[0]=r0; r[1]=r1; r[2]=r2; r[3]=r3;
    const u64 p[4] = {P0, P1, P2, P3};
    if (carry || !fe_lt(r, p)) fp_sub_p_inplace(r);
}

__device__ __forceinline__
void fp_sub(const u64 a[4], const u64 b[4], u64 r[4]) {
    u64 r0, r1, r2, r3, borrow;
    asm("sub.cc.u64   %0, %4, %8;\n\t"
        "subc.cc.u64  %1, %5, %9;\n\t"
        "subc.cc.u64  %2, %6, %10;\n\t"
        "subc.cc.u64  %3, %7, %11;\n\t"
        "subc.u64     %4, 0, 0;"
        : "=l"(r0),"=l"(r1),"=l"(r2),"=l"(r3),"=l"(borrow)
        : "l"(a[0]),"l"(a[1]),"l"(a[2]),"l"(a[3]),
          "l"(b[0]),"l"(b[1]),"l"(b[2]),"l"(b[3]));
    r[0]=r0; r[1]=r1; r[2]=r2; r[3]=r3;
    if (borrow & 1) {
        // a < b: add p back
        asm("add.cc.u64   %0, %0, %4;\n\t"
            "addc.cc.u64  %1, %1, %5;\n\t"
            "addc.cc.u64  %2, %2, %6;\n\t"
            "addc.u64     %3, %3, %7;"
            : "+l"(r[0]),"+l"(r[1]),"+l"(r[2]),"+l"(r[3])
            : "l"((u64)P0),"l"((u64)P1),"l"((u64)P2),"l"((u64)P3));
    }
}

__device__ __forceinline__
void fp_neg(const u64 a[4], u64 r[4]) {
    if ((a[0]|a[1]|a[2]|a[3]) == 0) { r[0]=r[1]=r[2]=r[3]=0; return; }
    const u64 p[4] = {P0, P1, P2, P3};
    fp_sub(p, a, r);
}

// 256×256 → 512-bit product — PTX carry-chain (mad.lo.cc / madc.hi)
__device__ __forceinline__
void mul512(const u64 a[4], const u64 b[4], u64 t[8]) {
    u64 t0,t1,t2,t3,t4,t5,t6,t7;
    asm(
        // ── Row 0: a0 × b[0..3] ─────────────────────────────────────────
        "mul.lo.u64     %0,  %8, %12;\n\t"
        "mul.hi.u64     %1,  %8, %12;\n\t"
        "mad.lo.cc.u64  %1,  %8, %13, %1;\n\t"
        "madc.hi.u64    %2,  %8, %13,  0;\n\t"
        "mad.lo.cc.u64  %2,  %8, %14, %2;\n\t"
        "madc.hi.u64    %3,  %8, %14,  0;\n\t"
        "mad.lo.cc.u64  %3,  %8, %15, %3;\n\t"
        "madc.hi.u64    %4,  %8, %15,  0;\n\t"
        // ── Row 1: a1 × b[0..3] ─────────────────────────────────────────
        "mad.lo.cc.u64  %1,  %9, %12, %1;\n\t"
        "madc.hi.cc.u64 %2,  %9, %12, %2;\n\t"
        "addc.u64       %3, %3,  0;\n\t"
        "mad.lo.cc.u64  %2,  %9, %13, %2;\n\t"
        "madc.hi.cc.u64 %3,  %9, %13, %3;\n\t"
        "addc.u64       %4, %4,  0;\n\t"
        "mad.lo.cc.u64  %3,  %9, %14, %3;\n\t"
        "madc.hi.cc.u64 %4,  %9, %14, %4;\n\t"
        "addc.u64       %5,  0,  0;\n\t"
        "mad.lo.cc.u64  %4,  %9, %15, %4;\n\t"
        "madc.hi.u64    %5,  %9, %15, %5;\n\t"
        // ── Row 2: a2 × b[0..3] ─────────────────────────────────────────
        "mad.lo.cc.u64  %2, %10, %12, %2;\n\t"
        "madc.hi.cc.u64 %3, %10, %12, %3;\n\t"
        "addc.u64       %4, %4,  0;\n\t"
        "mad.lo.cc.u64  %3, %10, %13, %3;\n\t"
        "madc.hi.cc.u64 %4, %10, %13, %4;\n\t"
        "addc.u64       %5, %5,  0;\n\t"
        "mad.lo.cc.u64  %4, %10, %14, %4;\n\t"
        "madc.hi.cc.u64 %5, %10, %14, %5;\n\t"
        "addc.u64       %6,  0,  0;\n\t"
        "mad.lo.cc.u64  %5, %10, %15, %5;\n\t"
        "madc.hi.u64    %6, %10, %15, %6;\n\t"
        // ── Row 3: a3 × b[0..3] ─────────────────────────────────────────
        "mad.lo.cc.u64  %3, %11, %12, %3;\n\t"
        "madc.hi.cc.u64 %4, %11, %12, %4;\n\t"
        "addc.u64       %5, %5,  0;\n\t"
        "mad.lo.cc.u64  %4, %11, %13, %4;\n\t"
        "madc.hi.cc.u64 %5, %11, %13, %5;\n\t"
        "addc.u64       %6, %6,  0;\n\t"
        "mad.lo.cc.u64  %5, %11, %14, %5;\n\t"
        "madc.hi.cc.u64 %6, %11, %14, %6;\n\t"
        "addc.u64       %7,  0,  0;\n\t"
        "mad.lo.cc.u64  %6, %11, %15, %6;\n\t"
        "madc.hi.u64    %7, %11, %15, %7;\n\t"
        : "=l"(t0),"=l"(t1),"=l"(t2),"=l"(t3),
          "=l"(t4),"=l"(t5),"=l"(t6),"=l"(t7)
        : "l"(a[0]),"l"(a[1]),"l"(a[2]),"l"(a[3]),
          "l"(b[0]),"l"(b[1]),"l"(b[2]),"l"(b[3])
    );
    t[0]=t0; t[1]=t1; t[2]=t2; t[3]=t3;
    t[4]=t4; t[5]=t5; t[6]=t6; t[7]=t7;
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

// 256×256 → 512-bit squaring — PTX, 10 muls (6 cross×2 + 4 diag)
__device__ __forceinline__
void sqr512(const u64 a[4], u64 t[8]) {
    u64 t0,t1,t2,t3,t4,t5,t6,t7;
    // ── Phase 1: cross products (each appears twice → computed once, doubled) ──
    asm(
        // (0,1) → slots 1,2
        "mul.lo.u64     %1, %8, %9;\n\t"
        "mul.hi.u64     %2, %8, %9;\n\t"
        // (0,2) → slots 2,3
        "mad.lo.cc.u64  %2, %8, %10, %2;\n\t"
        "madc.hi.u64    %3, %8, %10,  0;\n\t"
        // (1,2) → slots 3,4
        "mad.lo.cc.u64  %3, %9, %10, %3;\n\t"
        "madc.hi.u64    %4, %9, %10,  0;\n\t"
        // (0,3) → slots 3,4  (add into existing; carry may overflow)
        "mad.lo.cc.u64  %3, %8, %11, %3;\n\t"
        "madc.hi.cc.u64 %4, %8, %11, %4;\n\t"
        "addc.u64       %5,  0,  0;\n\t"
        // (1,3) → slots 4,5
        "mad.lo.cc.u64  %4, %9, %11, %4;\n\t"
        "madc.hi.cc.u64 %5, %9, %11, %5;\n\t"
        "addc.u64       %6,  0,  0;\n\t"
        // (2,3) → slots 5,6
        "mad.lo.cc.u64  %5, %10, %11, %5;\n\t"
        "madc.hi.u64    %6, %10, %11, %6;\n\t"
        "mov.u64        %0, 0;\n\t"
        "mov.u64        %7, 0;\n\t"
        : "=l"(t0),"=l"(t1),"=l"(t2),"=l"(t3),
          "=l"(t4),"=l"(t5),"=l"(t6),"=l"(t7)
        : "l"(a[0]),"l"(a[1]),"l"(a[2]),"l"(a[3])
    );
    // ── Phase 2: double the cross-product accumulator ─────────────────────────
    asm(
        "add.cc.u64   %0, %0, %0;\n\t"
        "addc.cc.u64  %1, %1, %1;\n\t"
        "addc.cc.u64  %2, %2, %2;\n\t"
        "addc.cc.u64  %3, %3, %3;\n\t"
        "addc.cc.u64  %4, %4, %4;\n\t"
        "addc.cc.u64  %5, %5, %5;\n\t"
        "addc.cc.u64  %6, %6, %6;\n\t"
        "addc.u64     %7, %7, %7;\n\t"
        : "+l"(t0),"+l"(t1),"+l"(t2),"+l"(t3),
          "+l"(t4),"+l"(t5),"+l"(t6),"+l"(t7)
    );
    // ── Phase 3: add diagonal a[i]² at positions 2i, 2i+1 ────────────────────
    asm(
        "mad.lo.cc.u64  %0, %8, %8, %0;\n\t"
        "madc.hi.cc.u64 %1, %8, %8, %1;\n\t"
        "madc.lo.cc.u64 %2, %9, %9, %2;\n\t"
        "madc.hi.cc.u64 %3, %9, %9, %3;\n\t"
        "madc.lo.cc.u64 %4, %10, %10, %4;\n\t"
        "madc.hi.cc.u64 %5, %10, %10, %5;\n\t"
        "madc.lo.cc.u64 %6, %11, %11, %6;\n\t"
        "madc.hi.u64    %7, %11, %11, %7;\n\t"
        : "+l"(t0),"+l"(t1),"+l"(t2),"+l"(t3),
          "+l"(t4),"+l"(t5),"+l"(t6),"+l"(t7)
        : "l"(a[0]),"l"(a[1]),"l"(a[2]),"l"(a[3])
    );
    t[0]=t0; t[1]=t1; t[2]=t2; t[3]=t3;
    t[4]=t4; t[5]=t5; t[6]=t6; t[7]=t7;
}

__device__ __forceinline__
void fp_sqr(const u64 a[4], u64 r[4]) {
    u64 t[8]; sqr512(a, t); reduce512(t, r);
}

// a^e mod p — binary exponentiation (kept for generic use)
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

// Optimized a^(p-2) mod p — secp256k1 addition chain, register-minimal form
//
// Peak simultaneous live arrays: 6 local (x2,x3,x22,x44,x88,t) + param a = 7
// = 28 u64s = 56 32-bit registers (vs 48 u64s / 96 regs in the naive version).
// Shorter intermediates (x6, x9, x11, x176, x220, x223) are all routed through t.
// Cost: 256S + 15M (same as before).
__device__
void fp_inv(const u64 a[4], u64 r[4]) {
    u64 x2[4], x3[4], x22[4], x44[4], x88[4], t[4];
    int j;
#define SQR(d)     fp_sqr(d, d)
#define MUL(d, s)  fp_mul(d, s, d)
#define CPY(d, s)  for(int _i=0;_i<4;_i++) (d)[_i]=(s)[_i]

    fp_sqr(a,  x2); fp_mul(x2, a,  x2);  // x2  = a^3        1S+1M
    fp_sqr(x2, x3); fp_mul(x3, a,  x3);  // x3  = a^7        1S+1M

    // x6, x9, x11 all routed through t (released before the next named value)
    CPY(t, x3);
    for(j=0;j<3;j++) SQR(t); MUL(t, x3); // t   = a^(2^6-1) 3S+1M
    for(j=0;j<3;j++) SQR(t); MUL(t, x3); // t   = a^(2^9-1) 3S+1M
    for(j=0;j<2;j++) SQR(t); MUL(t, x2); // t   = a^(2^11-1) 2S+1M

    // x22: need t (=x11) and x22 simultaneously for the final mul only
    CPY(x22, t);
    for(j=0;j<11;j++) SQR(x22); MUL(x22, t); // x22=a^(2^22-1) 11S+1M
    // t (=x11) released here

    // x44 through t, then saved — both live briefly during x88
    CPY(t, x22);
    for(j=0;j<22;j++) SQR(t); MUL(t, x22);   // t  = a^(2^44-1) 22S+1M
    CPY(x44, t);   // save x44; t still used for x88

    for(j=0;j<44;j++) SQR(t); MUL(t, x44);   // t  = a^(2^88-1) 44S+1M
    CPY(x88, t);   // save x88; t reused for x176 (peak: a,x2,x3,x22,x44,x88,t=7)

    for(j=0;j<88;j++) SQR(t); MUL(t, x88);   // t  = a^(2^176-1) 88S+1M
    // x88 released — 6 live: a,x2,x3,x22,x44,t

    for(j=0;j<44;j++) SQR(t); MUL(t, x44);   // t  = a^(2^220-1) 44S+1M
    // x44 released — 5 live: a,x2,x3,x22,t

    for(j=0;j< 3;j++) SQR(t); MUL(t,  x3);   // t  = a^(2^223-1)  3S+1M
    // x3 released — 4 live: a,x2,x22,t

    // p-2 tail
    for(j=0;j<23;j++) SQR(t); MUL(t, x22);   // 23S+1M; x22 released
    for(j=0;j< 5;j++) SQR(t); MUL(t,   a);   //  5S+1M
    for(j=0;j< 3;j++) SQR(t); MUL(t,  x2);   //  3S+1M; x2 released
    for(j=0;j< 2;j++) SQR(t);
    fp_mul(t, a, r);                           //  2S+1M  → a^(p-2)

#undef SQR
#undef MUL
#undef CPY
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

// ─── Pure affine + affine point addition ─────────────────────────────────────
// Cost: 1 fp_inv + 4M + 2S  (no Z coordinate, result is affine)
// Degenerate cases (same x) are astronomically unlikely in a Kangaroo walk;
// handled gracefully (result = 0,0) but practically never triggered.

__device__
void affine_add(
    const u64 ax[4], const u64 ay[4],   // animal (affine)
    const u64 bx[4], const u64 by[4],   // jump point (affine)
    u64 rx[4], u64 ry[4]                // result (affine)
) {
    u64 dx[4], dy[4], dinv[4], lam[4], lam2[4], tmp[4];
    fp_sub(bx, ax, dx);         // dx = bx − ax
    fp_sub(by, ay, dy);         // dy = by − ay

    // Degenerate: dx == 0  (same x-coordinate)
    if ((dx[0]|dx[1]|dx[2]|dx[3]) == 0) {
        if ((dy[0]|dy[1]|dy[2]|dy[3]) == 0) {
            // Point doubling: λ = 3x²/(2y)
            u64 x2[4], num[4], den[4];
            fp_sqr(ax, x2);
            const u64 three[4] = {3,0,0,0};
            fp_mul(three, x2, num);        // 3x²
            fp_add(ay, ay, den);           // 2y
            fp_inv(den, dinv);
            fp_mul(num, dinv, lam);
        } else {
            // Opposite points: result = point at infinity (set to 0)
            for (int i=0;i<4;i++) { rx[i]=0; ry[i]=0; }
            return;
        }
    } else {
        fp_inv(dx, dinv);
        fp_mul(dy, dinv, lam);             // λ = dy/dx
    }

    fp_sqr(lam, lam2);                     // λ²
    fp_sub(lam2, ax, tmp);                 // λ² − ax
    fp_sub(tmp,  bx, rx);                  // x3 = λ² − ax − bx
    fp_sub(ax,   rx, tmp);                 // ax − x3
    fp_mul(lam, tmp, lam2);                // λ(ax−x3)
    fp_sub(lam2, ay, ry);                  // y3 = λ(ax−x3) − ay
}

// ─── Scalar arithmetic mod n ──────────────────────────────────────────────────

__device__ __forceinline__
void sc_add(const u64 a[4], const u64 b[4], u64 r[4]) {
    u64 r0, r1, r2, r3, carry;
    asm("add.cc.u64   %0, %4, %8;\n\t"
        "addc.cc.u64  %1, %5, %9;\n\t"
        "addc.cc.u64  %2, %6, %10;\n\t"
        "addc.cc.u64  %3, %7, %11;\n\t"
        "addc.u64     %4, 0, 0;"
        : "=l"(r0),"=l"(r1),"=l"(r2),"=l"(r3),"=l"(carry)
        : "l"(a[0]),"l"(a[1]),"l"(a[2]),"l"(a[3]),
          "l"(b[0]),"l"(b[1]),"l"(b[2]),"l"(b[3]));
    r[0]=r0; r[1]=r1; r[2]=r2; r[3]=r3;
    const u64 n[4] = {N0, N1, N2, N3};
    if (carry || !fe_lt(r, n)) {
        u64 bw;
        asm("sub.cc.u64   %0, %0, %4;\n\t"
            "subc.cc.u64  %1, %1, %5;\n\t"
            "subc.cc.u64  %2, %2, %6;\n\t"
            "subc.cc.u64  %3, %3, %7;\n\t"
            "subc.u64     %4, 0, 0;"
            : "+l"(r[0]),"+l"(r[1]),"+l"(r[2]),"+l"(r[3]),"=l"(bw)
            : "l"((u64)N0),"l"((u64)N1),"l"((u64)N2),"l"((u64)N3));
    }
}

// ─── Structures ───────────────────────────────────────────────────────────────

struct JumpPoint {
    u64 x[4];   // affine
    u64 y[4];   // affine
    u64 s[4];   // scalar mod n
};

// Affine animal — no Z coordinate; jump_idx and DP check both use canonical affine x
struct Animal {
    u64 ax[4];       // affine x
    u64 ay[4];       // affine y
    u64 scalar[4];   // accumulated jump scalar mod n
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
