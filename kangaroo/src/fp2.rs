// sinGRAAL — F_{p²} arithmetic for GLS 4D endomorphism research
//
// F_{p²} = F_p[i] / (i² + 1)
//
// Choice: i² = -1. Valid because p ≡ 3 (mod 4) → -1 is non-square in F_p.
//
// Frobenius on F_{p²}/F_p:
//   (a + b·i)^p = a^p + b^p · i^p
//               = a   + b   · i^p     [Fermat: a^p=a, b^p=b for a,b ∈ F_p]
//   i^p = i^(p mod 4). Since p ≡ 3 (mod 4): i^p = i^3 = -i.
//   → Frobenius = complex conjugation: (a+bi)^p = a - b·i
//
// GLS twist endomorphism ψ on E: y²=x³+7 over F_{p²}:
//   Quadratic twist (d=-1): E^{-1}: y² = x³ - 7
//   Twist iso ι: E^{-1} → E: (x,y) → (-x, i·y)       [verified: (-y)²=(−x)³+7 ✓]
//   ψ = ι ∘ π ∘ ι^{-1}  where π = Frobenius
//   ι^{-1}(X,Y) = (-X, Y/i) = (-X, -iY)
//   π(-X, -iY) = (conj(-X), conj(-iY)) = (-conj(X), i·conj(Y))
//   ι(-conj(X), i·conj(Y)) = (conj(X), -conj(Y))
//   ∴ ψ(x,y) = (conj(x), -conj(y))
//
// Reference F_{p²} test point: ι(2,1) = (p-2, i)
//   Verify: i² = -1 ≡ (p-2)³+7 = (-2)³+7 = -8+7 = -1  ✓

#![allow(dead_code)]

use crate::secp::{Fe, fp_add, fp_sub, fp_mul, fp_neg, fp_inv, FIELD_P, Pt, BETA, BETA2};

// ─── F_{p²} element ──────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Fp2 {
    pub a: Fe,  // real part
    pub b: Fe,  // coefficient of i
}

impl Fp2 {
    pub const ZERO: Self = Fp2 { a: [0; 4], b: [0; 4] };
    pub const ONE:  Self = Fp2 { a: [1, 0, 0, 0], b: [0; 4] };
    pub const I:    Self = Fp2 { a: [0; 4], b: [1, 0, 0, 0] };

    pub fn from_fp(a: Fe) -> Self { Fp2 { a, b: [0; 4] } }
    pub fn is_fp(self) -> bool { self.b == [0; 4] }

    /// Frobenius / complex conjugation: (a+bi) → (a-bi)
    pub fn conj(self) -> Self { Fp2 { a: self.a, b: fp_neg(self.b) } }
}

pub fn fp2_add(x: Fp2, y: Fp2) -> Fp2 {
    Fp2 { a: fp_add(x.a, y.a), b: fp_add(x.b, y.b) }
}

pub fn fp2_sub(x: Fp2, y: Fp2) -> Fp2 {
    Fp2 { a: fp_sub(x.a, y.a), b: fp_sub(x.b, y.b) }
}

pub fn fp2_neg(x: Fp2) -> Fp2 {
    Fp2 { a: fp_neg(x.a), b: fp_neg(x.b) }
}

/// (a+bi)(c+di) = (ac−bd) + (ad+bc)i   [i²=−1]
pub fn fp2_mul(x: Fp2, y: Fp2) -> Fp2 {
    Fp2 {
        a: fp_sub(fp_mul(x.a, y.a), fp_mul(x.b, y.b)),
        b: fp_add(fp_mul(x.a, y.b), fp_mul(x.b, y.a)),
    }
}

pub fn fp2_sqr(x: Fp2) -> Fp2 {
    fp2_mul(x, x)
}

/// (a+bi)^{-1} = (a−bi) / (a²+b²)
pub fn fp2_inv(x: Fp2) -> Fp2 {
    let norm = fp_add(fp_mul(x.a, x.a), fp_mul(x.b, x.b));
    let n = fp_inv(norm);
    Fp2 { a: fp_mul(x.a, n), b: fp_neg(fp_mul(x.b, n)) }
}

pub fn fp2_from_hex(a_hi: u64, a_lo: u64, b_hi: u64, b_lo: u64) -> Fp2 {
    Fp2 {
        a: [a_lo, a_hi, 0, 0],
        b: [b_lo, b_hi, 0, 0],
    }
}

// ─── E(F_{p²}) point ─────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
pub struct Pt2 { pub x: Fp2, pub y: Fp2, pub inf: bool }

impl Pt2 {
    pub const INF: Self = Pt2 { x: Fp2::ZERO, y: Fp2::ZERO, inf: true };

    pub fn from_pt(p: Pt) -> Self {
        if p.inf { return Pt2::INF; }
        Pt2 { x: Fp2::from_fp(p.x), y: Fp2::from_fp(p.y), inf: false }
    }

    pub fn eq(self, b: Self) -> bool {
        if self.inf && b.inf { return true; }
        if self.inf || b.inf { return false; }
        self.x == b.x && self.y == b.y
    }

    pub fn is_fp_rational(self) -> bool {
        self.inf || (self.x.is_fp() && self.y.is_fp())
    }
}

pub fn pt2_neg(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 { x: p.x, y: fp2_neg(p.y), inf: false }
}

fn pt2_dbl(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    // lam = 3x² / 2y
    let x2  = fp2_sqr(p.x);
    let num = fp2_add(fp2_add(x2, x2), x2);   // 3x²
    let den = fp2_add(p.y, p.y);               // 2y
    let lam = fp2_mul(num, fp2_inv(den));
    let x3  = fp2_sub(fp2_sqr(lam), fp2_add(p.x, p.x));
    let y3  = fp2_sub(fp2_mul(lam, fp2_sub(p.x, x3)), p.y);
    Pt2 { x: x3, y: y3, inf: false }
}

pub fn pt2_add(a: Pt2, b: Pt2) -> Pt2 {
    if a.inf { return b; }
    if b.inf { return a; }
    let dx = fp2_sub(b.x, a.x);
    let dy = fp2_sub(b.y, a.y);
    if dx == Fp2::ZERO {
        return if dy == Fp2::ZERO { pt2_dbl(a) } else { Pt2::INF };
    }
    let lam = fp2_mul(dy, fp2_inv(dx));
    let x3  = fp2_sub(fp2_sub(fp2_sqr(lam), a.x), b.x);
    let y3  = fp2_sub(fp2_mul(lam, fp2_sub(a.x, x3)), a.y);
    Pt2 { x: x3, y: y3, inf: false }
}

pub fn pt2_scalar_mul(mut p: Pt2, mut k: Fe) -> Pt2 {
    let mut r = Pt2::INF;
    let mut i = 0u32;
    while i < 256 {
        if k[0] & 1 == 1 { r = pt2_add(r, p); }
        p = pt2_dbl(p);
        k = shr1(k);
        i += 1;
    }
    r
}

fn shr1(mut k: Fe) -> Fe {
    let mut carry = 0u64;
    for i in (0..4).rev() {
        let next_carry = k[i] & 1;
        k[i] = (k[i] >> 1) | (carry << 63);
        carry = next_carry;
    }
    k
}

// ─── GLV endomorphism on E(F_{p²}) ──────────────────────────────────────────

/// φ on E(F_{p²}): (x,y) → (β·x, y).  Acts as [λ] on E(F_p).
pub fn phi2(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 {
        x: Fp2 { a: fp_mul(BETA, p.x.a), b: fp_mul(BETA, p.x.b) },
        y: p.y,
        inf: false,
    }
}

pub fn phi2_sq(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 {
        x: Fp2 { a: fp_mul(BETA2, p.x.a), b: fp_mul(BETA2, p.x.b) },
        y: p.y,
        inf: false,
    }
}

/// ψ: GLS twist endomorphism on E(F_{p²}).
/// ψ(x,y) = (conj(x), -conj(y))
///
/// On E(F_p): conj(x)=x, conj(y)=y → ψ(P) = (x,-y) = -P  [scalar n-1]
/// On twist subgroup: acts as +id                           [scalar 1]
pub fn gls_psi(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 { x: p.x.conj(), y: fp2_neg(p.y.conj()), inf: false }
}

/// Frobenius π: (x,y) → (x^p, y^p) = (conj(x), conj(y))
pub fn frobenius(p: Pt2) -> Pt2 {
    if p.inf { return p; }
    Pt2 { x: p.x.conj(), y: p.y.conj(), inf: false }
}

// ─── The canonical twist-image test point ────────────────────────────────────
//
// ι(2, 1) = (-2 mod p, i) ∈ E(F_{p²}) \ E(F_p)
// Verify: i² = -1 = (-2)³+7 = -8+7 = -1  ✓
// x-coord is in F_p but y-coord (= i) is NOT in F_p.

pub fn twist_test_point() -> Pt2 {
    // x = -2 mod p  (F_p-rational)
    // y = i = 0 + 1·i  (NOT F_p-rational)
    let x_val: Fe = [FIELD_P[0] - 2, FIELD_P[1], FIELD_P[2], FIELD_P[3]];
    Pt2 {
        x: Fp2::from_fp(x_val),
        y: Fp2::I,
        inf: false,
    }
}

/// Verify P lies on E/F_{p²}: y² = x³ + 7
pub fn on_curve2(p: Pt2) -> bool {
    if p.inf { return true; }
    let y2  = fp2_sqr(p.y);
    let x3  = fp2_mul(fp2_sqr(p.x), p.x);
    let rhs = fp2_add(x3, Fp2::from_fp([7, 0, 0, 0]));
    y2 == rhs
}

// ─── 4D LLL decomposition for puzzle scalars ─────────────────────────────────
//
// The 4D GLS endomorphisms on E(F_{p²}) are {id, φ, ψ, φ∘ψ} with
// "scalar actions" that are well-defined only on each eigenspace separately.
//
// For a scalar k ∈ [0, n):
//   On E(F_p) sub:       {+1, +λ, -1, -λ} — only 2D independent (pairs cancel)
//   On twist E^d(F_p):   {+1, +λ, +1, +λ} — same 2 values
//
// True 4D requires k ≈ n (full 256-bit range), not k < 2^135.
//
// We compute: p̃ = p mod n  (the Frobenius scalar eigenvalue on E(F_p))
//             and check how much a 135-bit k can benefit.

pub fn p_mod_n() -> Fe {
    // p mod n: since p > n for secp256k1, compute p - n
    use crate::secp::FIELD_N;
    let mut r = [0u64; 4];
    let mut borrow = false;
    for i in 0..4 {
        let (s, b1) = FIELD_P[i].overflowing_sub(FIELD_N[i]);
        let (s, b2) = s.overflowing_sub(borrow as u64);
        r[i] = s;
        borrow = b1 || b2;
    }
    r
}
