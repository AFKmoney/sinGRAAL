// F_{p²} = F_p[u] / (u² + 1)
//
// Valid since p ≡ 3 (mod 4)  →  −1 is a non-residue mod p
// (Verified: secp256k1 prime p = 2²⁵⁶−2³²−977, p ≡ 3 mod 4)
//
// Element representation: a + b·u  stored as  [a, b]  (Fp2 = [Fe; 2])
// Frobenius (Galois conjugate over F_p): π(a+bu) = a−bu
// since u^p = (u²)^{(p-1)/2}·u = (−1)^{(p-1)/2}·u = −u for p≡3 mod 4

#![allow(dead_code)]

use crate::secp::*;

pub type Fp2 = [Fe; 2];  // [real, imag]

pub const FP2_ZERO: Fp2 = [[0; 4], [0; 4]];
pub const FP2_ONE:  Fp2 = [[1, 0, 0, 0], [0; 4]];

// ─── Fp2 arithmetic ──────────────────────────────────────────────────────────

pub fn fp2_add(a: Fp2, b: Fp2) -> Fp2 {
    [fp_add(a[0], b[0]), fp_add(a[1], b[1])]
}

pub fn fp2_sub(a: Fp2, b: Fp2) -> Fp2 {
    [fp_sub(a[0], b[0]), fp_sub(a[1], b[1])]
}

pub fn fp2_neg(a: Fp2) -> Fp2 {
    [fp_neg(a[0]), fp_neg(a[1])]
}

// (a+bu)(c+du) = (ac−bd) + (ad+bc)u
// Karatsuba: t₀=ac, t₁=bd, t₂=(a+b)(c+d)  → real=t₀−t₁, imag=t₂−t₀−t₁
pub fn fp2_mul(a: Fp2, b: Fp2) -> Fp2 {
    let t0 = fp_mul(a[0], b[0]);
    let t1 = fp_mul(a[1], b[1]);
    let t2 = fp_mul(fp_add(a[0], a[1]), fp_add(b[0], b[1]));
    [fp_sub(t0, t1), fp_sub(fp_sub(t2, t0), t1)]
}

pub fn fp2_sqr(a: Fp2) -> Fp2 { fp2_mul(a, a) }

// Galois conjugate (Frobenius over F_{p²}/F_p):  a+bu → a−bu
pub fn fp2_conj(a: Fp2) -> Fp2 { [a[0], fp_neg(a[1])] }

// Norm: N(a+bu) = a² + b²  ∈ F_p
pub fn fp2_norm(a: Fp2) -> Fe { fp_add(fp_sqr(a[0]), fp_sqr(a[1])) }

// (a+bu)^{−1} = (a−bu) / (a²+b²)
pub fn fp2_inv(a: Fp2) -> Fp2 {
    let n_inv = fp_inv(fp2_norm(a));
    [fp_mul(a[0], n_inv), fp_mul(fp_neg(a[1]), n_inv)]
}

// Scale by a field element:  α·(a+bu) = (α·a) + (α·b)·u
pub fn fp2_scale(a: Fp2, s: Fe) -> Fp2 { [fp_mul(s, a[0]), fp_mul(s, a[1])] }

pub fn fp2_eq(a: Fp2, b: Fp2) -> bool { a[0] == b[0] && a[1] == b[1] }
pub fn fp2_is_zero(a: Fp2) -> bool { a == FP2_ZERO }

// ─── Elliptic curve E over F_{p²}: y² = x³ + 7 ──────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Pt2 { pub x: Fp2, pub y: Fp2, pub inf: bool }

pub const INF2: Pt2 = Pt2 { x: FP2_ZERO, y: FP2_ZERO, inf: true };

/// Lift an F_p point into E(F_{p²}) via the natural embedding a → (a, 0).
pub fn pt_lift(p: Pt) -> Pt2 {
    if p.inf { return INF2; }
    Pt2 { x: [p.x, [0; 4]], y: [p.y, [0; 4]], inf: false }
}

pub fn pt2_neg(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 { x: p.x, y: fp2_neg(p.y), inf: false }
}

pub fn pt2_eq(a: Pt2, b: Pt2) -> bool {
    if a.inf && b.inf { return true; }
    if a.inf || b.inf { return false; }
    fp2_eq(a.x, b.x) && fp2_eq(a.y, b.y)
}

pub fn pt2_add(a: Pt2, b: Pt2) -> Pt2 {
    if a.inf { return b; }
    if b.inf { return a; }
    if fp2_eq(a.x, b.x) {
        if !fp2_eq(a.y, b.y) { return INF2; }
        return pt2_dbl(a);
    }
    let dx  = fp2_sub(b.x, a.x);
    let dy  = fp2_sub(b.y, a.y);
    let m   = fp2_mul(dy, fp2_inv(dx));
    let x3  = fp2_sub(fp2_sub(fp2_sqr(m), a.x), b.x);
    let y3  = fp2_sub(fp2_mul(m, fp2_sub(a.x, x3)), a.y);
    Pt2 { x: x3, y: y3, inf: false }
}

fn pt2_dbl(a: Pt2) -> Pt2 {
    if a.inf { return a; }
    let three: Fe = [3, 0, 0, 0];
    let x2      = fp2_sqr(a.x);
    let three_x2 = fp2_scale(x2, three);
    let two_y   = fp2_add(a.y, a.y);
    let m       = fp2_mul(three_x2, fp2_inv(two_y));
    let x3      = fp2_sub(fp2_sub(fp2_sqr(m), a.x), a.x);
    let y3      = fp2_sub(fp2_mul(m, fp2_sub(a.x, x3)), a.y);
    Pt2 { x: x3, y: y3, inf: false }
}

pub fn pt2_scalar_mul(mut p: Pt2, mut k: Fe) -> Pt2 {
    let mut r = INF2;
    for _ in 0..256 {
        if k[0] & 1 == 1 { r = pt2_add(r, p); }
        p = pt2_dbl(p);
        let mut carry = 0u64;
        for j in (0..4).rev() {
            let nk = (k[j] >> 1) | (carry << 63);
            carry = k[j] & 1;
            k[j] = nk;
        }
    }
    r
}

// ─── Frobenius endomorphism on E/F_{p²} ──────────────────────────────────────

/// π: (x, y) → (x^p, y^p) = (conj(x), conj(y))
/// On lifted F_p points: π(P) = P (identity — Frobenius fixes F_p points).
/// On general F_{p²} points: π is a non-trivial group automorphism.
pub fn frobenius_pt(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 { x: fp2_conj(p.x), y: fp2_conj(p.y), inf: false }
}

/// CM endomorphism φ lifted to F_{p²}: (x,y) → (β·x, y)
pub fn phi_pt2(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 { x: fp2_scale(p.x, BETA), y: p.y, inf: false }
}
