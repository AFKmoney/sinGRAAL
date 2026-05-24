// sinGRAAL — secp256k1 CUDA arithmetic with Toom-Cook-6 field multiplication
//
// Toom-Cook-6 for 256×256 → 512-bit integer multiplication:
//
//   Split each 256-bit input A into 6 pieces of 43 bits (t = 2^43):
//     A = a[0] + a[1]*t + a[2]*t^2 + a[3]*t^3 + a[4]*t^4 + a[5]*t^5
//     Pieces: a[0..4] < 2^43,  a[5] < 2^41  (41 bits)
//
//   Evaluate A(x) at 11 points {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ∞}:
//     A(k) = a[0] + k*(a[1] + k*(a[2] + k*(a[3] + k*(a[4] + k*a[5]))))  [Horner]
//     Max at k=9: A(9) < 2^61  → fits in u64
//     A(∞) = a[5] < 2^41
//
//   Compute 11 pointwise products C(k) = A(k) * B(k):
//     Each C(k) < (2^61)^2 = 2^122  → fits in u128
//
//   Interpolate via Newton forward differences:
//     Key property: polynomial with non-negative coefficients evaluated at
//     non-negative integers has ALL forward differences ≥ 0.
//     → pure unsigned arithmetic throughout, no sign tracking needed.
//     D[k][j] = D[k-1][j+1] − D[k-1][j]   (non-negative u128 subtraction)
//     c_k = Δ^k g(0) / k!                   (exact integer division)
//     All intermediate values < 2^122 → fit in u128.
//
//   Reconstruct:
//     t[0..7] = Σ c_k * 2^(43k)   for k=0..10
//
//   Multiplication count vs schoolbook (4×4):
//     Schoolbook: 16 MAD instructions
//     Toom-6:     11 MAD instructions  (31% reduction)
//     Plus forward-difference table: 45 u128 subtractions
//     Plus 9 exact divisions by k! (emulated on GPU via NVCC)
//
// Default: Toom-6. Compile with -DUSE_SCHOOLBOOK to compare.

#pragma once
#include <stdint.h>

typedef unsigned long long  u64;
typedef unsigned __int128   u128;
typedef unsigned int        u32;

// ─── secp256k1 constants ──────────────────────────────────────────────────────

#define P0  0xFFFFFFFEFFFFFC2FULL
#define P1  0xFFFFFFFFFFFFFFFFULL
#define P2  0xFFFFFFFFFFFFFFFFULL
#define P3  0xFFFFFFFFFFFFFFFFULL

#define BETA0  0xC1396C28719501EEULL
#define BETA1  0x9CF0497512F58995ULL
#define BETA2_ 0x6E64479EAC3434E9ULL
#define BETA3  0x7AE96A2B657C0710ULL

#define BETA2_0  0x3EC693D68E6AFA40ULL
#define BETA2_1  0x630FB68AED0A766AULL
#define BETA2_2  0x919BB86153CBCB16ULL
#define BETA2_3  0x851695D49A83F8EFULL

#define N0  0xBFD25E8CD0364141ULL
#define N1  0xBAAEDCE6AF48A03BULL
#define N2  0xFFFFFFFFFFFFFFFEULL
#define N3  0xFFFFFFFFFFFFFFFFULL

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

// ─── Schoolbook 256×256 → 512-bit (reference, -DUSE_SCHOOLBOOK) ───────────────

__device__ __forceinline__
void mul512_schoolbook(const u64 a[4], const u64 b[4], u64 t[8]) {
    u64 t0,t1,t2,t3,t4,t5,t6,t7;
    asm(
        "mul.lo.u64     %0,  %8, %12;\n\t"
        "mul.hi.u64     %1,  %8, %12;\n\t"
        "mad.lo.cc.u64  %1,  %8, %13, %1;\n\t"
        "madc.hi.u64    %2,  %8, %13,  0;\n\t"
        "mad.lo.cc.u64  %2,  %8, %14, %2;\n\t"
        "madc.hi.u64    %3,  %8, %14,  0;\n\t"
        "mad.lo.cc.u64  %3,  %8, %15, %3;\n\t"
        "madc.hi.u64    %4,  %8, %15,  0;\n\t"
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

// ─── Toom-Cook-6 helpers ──────────────────────────────────────────────────────

// Split 256-bit A (4 u64 limbs, little-endian) into 6 pieces of 43 bits:
//   A = p[0] + p[1]*2^43 + p[2]*2^86 + p[3]*2^129 + p[4]*2^172 + p[5]*2^215
//   p[0..4] < 2^43,  p[5] < 2^41   (41 bits since 5*43+41 = 256)
__device__ __forceinline__
void toom6_split(const u64 a[4], u64 p[6]) {
    const u64 M43 = (1ULL << 43) - 1;
    p[0] =  a[0] & M43;
    p[1] = (a[0] >> 43) | ((a[1] & 0x3FFFFFULL) << 21);   // bits 43..85
    p[2] = (a[1] >> 22) | ((a[2] & 1ULL) << 42);          // bits 86..128
    p[3] = (a[2] >> 1)  & M43;                             // bits 129..171
    p[4] = (a[2] >> 44) | ((a[3] & 0x7FFFFFULL) << 20);   // bits 172..214
    p[5] =  a[3] >> 23;                                    // bits 215..255 (41 bits)
}

// Evaluate a degree-5 polynomial at point k using Horner's method.
// p[0..5] < 2^43.  For k ≤ 9: result < 2^61 (fits in u64).
__device__ __forceinline__
u64 toom6_eval(const u64 p[6], u64 k) {
    u64 v = p[5];
    v = p[4] + k * v;
    v = p[3] + k * v;
    v = p[2] + k * v;
    v = p[1] + k * v;
    v = p[0] + k * v;
    return v;
}

// Add a u128 value at bit position bit_pos into a 9-limb u64 accumulator.
// Handles arbitrary bit offsets without UB (shift=0 case separated).
__device__ __forceinline__
void toom6_add_at(u64 out[9], u128 val, int bit_pos) {
    int limb0 = bit_pos / 64;
    int shift  = bit_pos % 64;
    u64 lo = (u64)val;
    u64 hi = (u64)(val >> 64);

    u64 part0, part1, part2;
    if (shift == 0) {
        part0 = lo;
        part1 = hi;
        part2 = 0;
    } else {
        part0 =  lo << shift;
        part1 = (lo >> (64 - shift)) | (hi << shift);
        part2 =  hi >> (64 - shift);
    }

    u128 s = (u128)out[limb0] + part0;
    out[limb0] = (u64)s;
    s = (u128)out[limb0 + 1] + part1 + (s >> 64);
    out[limb0 + 1] = (u64)s;
    if (limb0 + 2 < 9) {
        s = (u128)out[limb0 + 2] + part2 + (s >> 64);
        out[limb0 + 2] = (u64)s;
        u64 carry = (u64)(s >> 64);
        for (int i = limb0 + 3; carry && i < 9; i++) {
            s = (u128)out[i] + carry;
            out[i]  = (u64)s;
            carry   = (u64)(s >> 64);
        }
    }
}

// ─── Toom-Cook-6: 256×256 → 512-bit ──────────────────────────────────────────
//
// Mathematical guarantee:
//   All 11 products C(k) are non-negative u128.
//   All Newton forward differences D[k][j] are non-negative (theorem: polynomial
//   with non-negative coefficients, positive integer evaluation points).
//   → unsigned arithmetic throughout, no sign bits.
//   → All intermediates < 2^122 (bounded by max evaluation C(9) < 2^122).
//
// Exact division by k! (k=2..9):
//   Divisors: 2, 6, 24, 120, 720, 5040, 40320, 362880.
//   NVCC emulates u128/u64 division via software (~10-20 instructions each).
//   Result guaranteed exact (divisibility from Newton forward difference theory).
__device__
void mul512_toom6(const u64 a[4], const u64 b[4], u64 t[8]) {
    // ── Step 1: Split into 6 pieces each ─────────────────────────────────────
    u64 ap[6], bp[6];
    toom6_split(a, ap);
    toom6_split(b, bp);

    // ── Step 2: Evaluate at {0,1,...,9} using Horner (11 evaluations × 6 ops)
    u64 av[10], bv[10];
    for (int k = 0; k < 10; k++) {
        av[k] = toom6_eval(ap, (u64)k);
        bv[k] = toom6_eval(bp, (u64)k);
    }

    // ── Step 3: 11 pointwise multiplications (the core "11 MADs" ─────────────
    // cv[0..9] = products at evaluation points {0..9}
    // cv[10]   = a[5]*b[5] = product at ∞ (leading coefficient c10 directly)
    u128 cv[11];
    for (int k = 0; k < 10; k++) cv[k] = (u128)av[k] * bv[k];
    cv[10] = (u128)ap[5] * bp[5];   // c10 = A(∞)*B(∞)

    // ── Step 4: Subtract c10 contribution from each g(j) = C(j) − c10·j^10 ──
    // j^10 for j=0..9 (precomputed): 0^10=0, 1^10=1, ..., 9^10=3486784401
    // These are exact: g(j) ≥ 0 since C(j) ≥ c10·j^10 for non-neg coefficients.
    const u64 j10[10] = {
        0ULL, 1ULL, 1024ULL, 59049ULL, 1048576ULL,
        9765625ULL, 60466176ULL, 282475249ULL, 1073741824ULL, 3486784401ULL
    };
    for (int j = 0; j < 10; j++) cv[j] -= cv[10] * j10[j];
    // cv[0..9] now hold g(0)..g(9) where g has degree ≤ 9, coefficients c0..c9.

    // ── Step 5: Newton forward differences to interpolate c0..c9 ─────────────
    // D[j] = D^(k)[j] updated in place.  At level k: D[j] = Δ^k g(j).
    // All differences non-negative (polynomial with non-neg coefficients).
    u128 D[10];
    for (int j = 0; j < 10; j++) D[j] = cv[j];

    u128 coeff[11];
    coeff[0]  = D[0];         // c0 = g(0)  (no division, k!=0!=1)
    coeff[10] = cv[10];       // c10 from ∞ evaluation

    for (int k = 1; k < 10; k++) {
        for (int j = 0; j < 10 - k; j++) D[j] = D[j + 1] - D[j];
        coeff[k] = D[0];      // Δ^k g(0) before division
    }

    // ── Step 6: Divide coeff[k] by k! (exact integer division) ───────────────
    // coeff[k] = Δ^k g(0) / k!  gives the polynomial coefficient c_k.
    // Result bounds: c_k ≤ min(k+1,11-k) * 2^86 < 2^89.  Fits in u128.
    // coeff[1] /= 1  (no-op)
    coeff[2] /= 2ULL;
    coeff[3] /= 6ULL;
    coeff[4] /= 24ULL;
    coeff[5] /= 120ULL;
    coeff[6] /= 720ULL;
    coeff[7] /= 5040ULL;
    coeff[8] /= 40320ULL;
    coeff[9] /= 362880ULL;

    // ── Step 7: Reconstruct 512-bit product t = Σ c_k * 2^(43k) ─────────────
    // Bit positions: k*43 for k=0..10
    //   k=0:  bit   0  (limb 0, shift  0)
    //   k=1:  bit  43  (limb 0, shift 43)
    //   k=2:  bit  86  (limb 1, shift 22)
    //   k=3:  bit 129  (limb 2, shift  1)
    //   k=4:  bit 172  (limb 2, shift 44)
    //   k=5:  bit 215  (limb 3, shift 23)
    //   k=6:  bit 258  (limb 4, shift  2)
    //   k=7:  bit 301  (limb 4, shift 45)
    //   k=8:  bit 344  (limb 5, shift 24)
    //   k=9:  bit 387  (limb 6, shift  3)
    //   k=10: bit 430  (limb 6, shift 46)
    // c_k < 2^89 → spans at most 3 u64 limbs from any alignment.
    u64 out[9] = {0,0,0,0,0,0,0,0,0};
    toom6_add_at(out, coeff[ 0],   0);
    toom6_add_at(out, coeff[ 1],  43);
    toom6_add_at(out, coeff[ 2],  86);
    toom6_add_at(out, coeff[ 3], 129);
    toom6_add_at(out, coeff[ 4], 172);
    toom6_add_at(out, coeff[ 5], 215);
    toom6_add_at(out, coeff[ 6], 258);
    toom6_add_at(out, coeff[ 7], 301);
    toom6_add_at(out, coeff[ 8], 344);
    toom6_add_at(out, coeff[ 9], 387);
    toom6_add_at(out, coeff[10], 430);

    for (int i = 0; i < 8; i++) t[i] = out[i];
    // out[8] must be 0 (A*B < 2^512, no overflow beyond 8 limbs).
}

// ─── Select mul512 implementation ────────────────────────────────────────────

__device__ __forceinline__
void mul512(const u64 a[4], const u64 b[4], u64 t[8]) {
#ifdef USE_SCHOOLBOOK
    mul512_schoolbook(a, b, t);
#else
    mul512_toom6(a, b, t);   // default: Toom-Cook-6
#endif
}

// ─── Fast reduction mod p = 2^256 − 2^32 − 977 ───────────────────────────────

__device__ __forceinline__
void reduce512(const u64 t[8], u64 r[4]) {
    u64 h0=t[4], h1=t[5], h2=t[6], h3=t[7];
    u64 r0,r1,r2,r3,r4;
    asm(
        "mul.lo.u64     %0, %5, 977;\n\t"
        "mul.hi.u64     %1, %5, 977;\n\t"
        "mad.lo.cc.u64  %1, %6, 977, %1;\n\t"
        "madc.hi.cc.u64 %2, %6, 977,  0;\n\t"
        "mad.lo.cc.u64  %2, %7, 977, %2;\n\t"
        "madc.hi.cc.u64 %3, %7, 977,  0;\n\t"
        "mad.lo.cc.u64  %3, %8, 977, %3;\n\t"
        "madc.hi.u64    %4, %8, 977,  0;\n\t"
        : "=l"(r0),"=l"(r1),"=l"(r2),"=l"(r3),"=l"(r4)
        : "l"(h0),"l"(h1),"l"(h2),"l"(h3)
    );
    asm(
        "add.cc.u64    %0, %0, %5;\n\t"
        "addc.cc.u64   %1, %1, %6;\n\t"
        "addc.cc.u64   %2, %2, %7;\n\t"
        "addc.cc.u64   %3, %3, %8;\n\t"
        "addc.u64      %4, %4, %9;\n\t"
        : "+l"(r0),"+l"(r1),"+l"(r2),"+l"(r3),"+l"(r4)
        : "l"(h0 << 32), "l"(h0 >> 32), "l"(h1 >> 32), "l"(h2 >> 32), "l"(h3 >> 32)
    );
    asm(
        "add.cc.u64    %0, %0, %4;\n\t"
        "addc.cc.u64   %1, %1, %5;\n\t"
        "addc.cc.u64   %2, %2, %6;\n\t"
        "addc.u64      %3, %3,  0;\n\t"
        : "+l"(r1),"+l"(r2),"+l"(r3),"+l"(r4)
        : "l"(h1 << 32), "l"(h2 << 32), "l"(h3 << 32)
    );
    asm(
        "add.cc.u64    %0, %0, %5;\n\t"
        "addc.cc.u64   %1, %1, %6;\n\t"
        "addc.cc.u64   %2, %2, %7;\n\t"
        "addc.cc.u64   %3, %3, %8;\n\t"
        "addc.u64      %4, %4,  0;\n\t"
        : "+l"(r0),"+l"(r1),"+l"(r2),"+l"(r3),"+l"(r4)
        : "l"(t[0]),"l"(t[1]),"l"(t[2]),"l"(t[3])
    );
    if (r4) {
        u64 e = r4;
        asm(
            "mad.lo.cc.u64   %0, %4, 977, %0;\n\t"
            "madc.hi.cc.u64  %1, %4, 977, %1;\n\t"
            "addc.cc.u64     %2, %2, 0;\n\t"
            "addc.u64        %3, %3, 0;\n\t"
            : "+l"(r0),"+l"(r1),"+l"(r2),"+l"(r3)
            : "l"(e)
        );
        asm(
            "add.cc.u64    %0, %0, %2;\n\t"
            "addc.cc.u64   %1, %1, %3;\n\t"
            : "+l"(r0),"+l"(r1)
            : "l"(e << 32), "l"(e >> 32)
        );
    }
    r[0]=r0; r[1]=r1; r[2]=r2; r[3]=r3;
    const u64 p[4] = {P0, P1, P2, P3};
    if (!fe_lt(r, p)) fp_sub_p_inplace(r);
}

__device__ __forceinline__
void fp_mul(const u64 a[4], const u64 b[4], u64 r[4]) {
    u64 t[8]; mul512(a, b, t); reduce512(t, r);
}

// 256×256 → 512-bit squaring (10 cross-products + 4 diagonals, PTX-optimized)
__device__ __forceinline__
void sqr512(const u64 a[4], u64 t[8]) {
    u64 t0,t1,t2,t3,t4,t5,t6,t7;
    asm(
        "mul.lo.u64     %1, %8, %9;\n\t"
        "mul.hi.u64     %2, %8, %9;\n\t"
        "mad.lo.cc.u64  %2, %8, %10, %2;\n\t"
        "madc.hi.u64    %3, %8, %10,  0;\n\t"
        "mad.lo.cc.u64  %3, %9, %10, %3;\n\t"
        "madc.hi.u64    %4, %9, %10,  0;\n\t"
        "mad.lo.cc.u64  %3, %8, %11, %3;\n\t"
        "madc.hi.cc.u64 %4, %8, %11, %4;\n\t"
        "addc.u64       %5,  0,  0;\n\t"
        "mad.lo.cc.u64  %4, %9, %11, %4;\n\t"
        "madc.hi.cc.u64 %5, %9, %11, %5;\n\t"
        "addc.u64       %6,  0,  0;\n\t"
        "mad.lo.cc.u64  %5, %10, %11, %5;\n\t"
        "madc.hi.u64    %6, %10, %11, %6;\n\t"
        "mov.u64        %0, 0;\n\t"
        "mov.u64        %7, 0;\n\t"
        : "=l"(t0),"=l"(t1),"=l"(t2),"=l"(t3),
          "=l"(t4),"=l"(t5),"=l"(t6),"=l"(t7)
        : "l"(a[0]),"l"(a[1]),"l"(a[2]),"l"(a[3])
    );
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

// Field inversion: a^(p-2) mod p — addition chain, 256S + 15M
__device__
void fp_inv(const u64 a[4], u64 r[4]) {
    u64 x2[4], x3[4], x22[4], x44[4], x88[4], t[4];
    int j;
#define SQR(d)     fp_sqr(d, d)
#define MUL(d, s)  fp_mul(d, s, d)
#define CPY(d, s)  for(int _i=0;_i<4;_i++) (d)[_i]=(s)[_i]

    fp_sqr(a,  x2); fp_mul(x2, a,  x2);
    fp_sqr(x2, x3); fp_mul(x3, a,  x3);
    CPY(t, x3);
    for(j=0;j<3;j++) SQR(t); MUL(t, x3);
    for(j=0;j<3;j++) SQR(t); MUL(t, x3);
    for(j=0;j<2;j++) SQR(t); MUL(t, x2);
    CPY(x22, t);
    for(j=0;j<11;j++) SQR(x22); MUL(x22, t);
    CPY(t, x22);
    for(j=0;j<22;j++) SQR(t); MUL(t, x22);
    CPY(x44, t);
    for(j=0;j<44;j++) SQR(t); MUL(t, x44);
    CPY(x88, t);
    for(j=0;j<88;j++) SQR(t); MUL(t, x88);
    for(j=0;j<44;j++) SQR(t); MUL(t, x44);
    for(j=0;j< 3;j++) SQR(t); MUL(t,  x3);
    for(j=0;j<23;j++) SQR(t); MUL(t, x22);
    for(j=0;j< 5;j++) SQR(t); MUL(t,   a);
    for(j=0;j< 3;j++) SQR(t); MUL(t,  x2);
    for(j=0;j< 2;j++) SQR(t);
    fp_mul(t, a, r);
#undef SQR
#undef MUL
#undef CPY
}

// ─── GLV endomorphisms ────────────────────────────────────────────────────────

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

__device__ __forceinline__
void canonical_x_affine(const u64 ax[4], u64 r[4]) {
    u64 x1[4], x2[4];
    psi_x(ax, x1);
    psi2_x(ax, x2);
    for (int i = 0; i < 4; i++) r[i] = ax[i];
    if (fe_lt(x1, r)) for (int i = 0; i < 4; i++) r[i] = x1[i];
    if (fe_lt(x2, r)) for (int i = 0; i < 4; i++) r[i] = x2[i];
}

// ─── Affine point addition ────────────────────────────────────────────────────

__device__
void affine_add(
    const u64 ax[4], const u64 ay[4],
    const u64 bx[4], const u64 by[4],
    u64 rx[4], u64 ry[4]
) {
    u64 dx[4], dy[4], dinv[4], lam[4], lam2[4], tmp[4];
    fp_sub(bx, ax, dx);
    fp_sub(by, ay, dy);

    if ((dx[0]|dx[1]|dx[2]|dx[3]) == 0) {
        if ((dy[0]|dy[1]|dy[2]|dy[3]) == 0) {
            u64 x2[4], num[4], den[4];
            fp_sqr(ax, x2);
            const u64 three[4] = {3,0,0,0};
            fp_mul(three, x2, num);
            fp_add(ay, ay, den);
            fp_inv(den, dinv);
            fp_mul(num, dinv, lam);
        } else {
            for (int i=0;i<4;i++) { rx[i]=0; ry[i]=0; }
            return;
        }
    } else {
        fp_inv(dx, dinv);
        fp_mul(dy, dinv, lam);
    }

    fp_sqr(lam, lam2);
    fp_sub(lam2, ax, tmp);
    fp_sub(tmp,  bx, rx);
    fp_sub(ax,   rx, tmp);
    fp_mul(lam, tmp, lam2);
    fp_sub(lam2, ay, ry);
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

struct JumpPoint { u64 x[4]; u64 y[4]; u64 s[4]; };

struct Animal {
    u64 ax[4]; u64 ay[4]; u64 scalar[4];
    u32 is_wild; u32 pad[3];
};

struct DPEntry {
    u64 canon_x[4]; u64 scalar[4];
    u32 is_wild; u32 pad[3];
};
