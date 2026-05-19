// secp256k1 field arithmetic and point operations
// All field elements: little-endian [u64; 4] (limb 0 = bits 0..63)

// p = 2^256 - 2^32 - 977
pub const P: [u64; 4] = [
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
];

// n = group order
pub const N: [u64; 4] = [
    0xBFD25E8CD0364141,
    0xBAAEDCE6AF48A03B,
    0xFFFFFFFFFFFFFFFE,
    0xFFFFFFFFFFFFFFFF,
];

// Generator
pub const GX: [u64; 4] = [
    0x59F2815B16F81798,
    0x029BFCDB2DCE28D9,
    0x55A06295CE870B07,
    0x79BE667EF9DCBBAC,
];
pub const GY: [u64; 4] = [
    0x9C47D08FFB10D4B8,
    0xFD17B448A6855419,
    0x5DA4FBFC0E1108A8,
    0x483ADA7726A3C465,
];

// β: cube root of unity mod p  (ψ(x,y) = (β·x, y))
pub const BETA: [u64; 4] = [
    0xC1396C28719501EE,
    0x9CF0497512F58995,
    0x6E64479EAC3434E9,
    0x7AE96A2B657C0710,
];

// β² = p - 1 - β  (since β² + β + 1 ≡ 0 mod p)
pub const BETA2: [u64; 4] = [
    0x3EC693D68E6AFA40,
    0x630FB68AED0A766A,
    0x919BB86153CBCB16,
    0x851695D49A83F8EF,
];

// λ: eigenvalue of ψ in scalar field  (ψ(P) = λ·P)
pub const LAMBDA: [u64; 4] = [
    0xDF02967C1B23BD72,
    0x122E22EA20816678,
    0xA5261C028812645A,
    0x5363AD4CC05C30E0,
];

pub type Fe = [u64; 4];

// ─── Comparison ───────────────────────────────────────────────────────────────

pub fn fe_lt(a: Fe, b: Fe) -> bool {
    for i in (0..4).rev() {
        if a[i] < b[i] { return true; }
        if a[i] > b[i] { return false; }
    }
    false
}

// ─── Field arithmetic mod p ───────────────────────────────────────────────────

fn sub_p(a: Fe) -> Fe {
    let mut r = [0u64; 4];
    let mut borrow = false;
    for i in 0..4 {
        let (s, b1) = a[i].overflowing_sub(P[i]);
        let (s, b2) = s.overflowing_sub(borrow as u64);
        r[i] = s;
        borrow = b1 || b2;
    }
    r
}

pub fn fp_add(a: Fe, b: Fe) -> Fe {
    let mut r = [0u64; 4];
    let mut carry = false;
    for i in 0..4 {
        let (s, c1) = a[i].overflowing_add(b[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        r[i] = s;
        carry = c1 || c2;
    }
    if carry || !fe_lt(r, P) { sub_p(r) } else { r }
}

pub fn fp_sub(a: Fe, b: Fe) -> Fe {
    let mut r = [0u64; 4];
    let mut borrow = false;
    for i in 0..4 {
        let (s, b1) = a[i].overflowing_sub(b[i]);
        let (s, b2) = s.overflowing_sub(borrow as u64);
        r[i] = s;
        borrow = b1 || b2;
    }
    if borrow {
        let mut carry = false;
        for i in 0..4 {
            let (s, c1) = r[i].overflowing_add(P[i]);
            let (s, c2) = s.overflowing_add(carry as u64);
            r[i] = s;
            carry = c1 || c2;
        }
    }
    r
}

pub fn fp_neg(a: Fe) -> Fe {
    if a == [0u64; 4] { return a; }
    fp_sub(P, a)
}

// 256×256 → 512 schoolbook
fn mul512(a: Fe, b: Fe) -> [u64; 8] {
    let mut r = [0u64; 8];
    for i in 0..4 {
        let mut carry = 0u64;
        for j in 0..4 {
            let uv = (a[i] as u128) * (b[j] as u128)
                   + r[i + j] as u128
                   + carry as u128;
            r[i + j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        r[i + 4] = carry;
    }
    r
}

// Reduce 512-bit number mod p using p = 2^256 − (2^32 + 977)
fn reduce512(t: [u64; 8]) -> Fe {
    let lo = [t[0], t[1], t[2], t[3]];
    let hi = [t[4], t[5], t[6], t[7]];

    // hi * (2^32 + 977) as 320-bit (5 limbs)
    // hi * 977
    let mut acc977 = [0u64; 5];
    let mut carry = 0u128;
    for i in 0..4 {
        let p = hi[i] as u128 * 977 + carry;
        acc977[i] = p as u64;
        carry = p >> 64;
    }
    acc977[4] = carry as u64;

    // hi * 2^32 (shift left 32 bits within 320-bit)
    let hi32 = [
        hi[0] << 32,
        (hi[0] >> 32) | (hi[1] << 32),
        (hi[1] >> 32) | (hi[2] << 32),
        (hi[2] >> 32) | (hi[3] << 32),
        hi[3] >> 32,
    ];

    // sum = hi*(2^32+977)
    let mut sum5 = [0u64; 5];
    let mut c = false;
    for i in 0..5 {
        let (s, c1) = acc977[i].overflowing_add(hi32[i]);
        let (s, c2) = s.overflowing_add(c as u64);
        sum5[i] = s;
        c = c1 || c2;
    }

    // r = lo + sum5[0..4], plus overflow into sum5[4]
    let mut r = [0u64; 5];
    let mut c = false;
    for i in 0..4 {
        let (s, c1) = lo[i].overflowing_add(sum5[i]);
        let (s, c2) = s.overflowing_add(c as u64);
        r[i] = s;
        c = c1 || c2;
    }
    r[4] = sum5[4].wrapping_add(c as u64);

    // If r[4] > 0, one more reduction pass
    if r[4] > 0 {
        let extra = r[4] as u128 * (977 + (1u128 << 32));
        let mut c = false;
        let (s, c1) = r[0].overflowing_add(extra as u64);
        let (s, c2) = s.overflowing_add(c as u64);
        r[0] = s; c = c1 || c2;
        let high = (extra >> 64) as u64;
        let (s, c1) = r[1].overflowing_add(high);
        let (s, c2) = s.overflowing_add(c as u64);
        r[1] = s; c = c1 || c2;
        for i in 2..4 { let (s, c1) = r[i].overflowing_add(c as u64); r[i] = s; c = c1; }
        r[4] = 0;
        let _ = c;
    }

    let mut res = [r[0], r[1], r[2], r[3]];
    if !fe_lt(res, P) { res = sub_p(res); }
    res
}

pub fn fp_mul(a: Fe, b: Fe) -> Fe {
    reduce512(mul512(a, b))
}

pub fn fp_sqr(a: Fe) -> Fe {
    fp_mul(a, a)
}

// a^e mod p via binary exp
pub fn fp_pow(mut base: Fe, mut exp: Fe) -> Fe {
    let mut r = [0u64; 4]; r[0] = 1;
    loop {
        if exp[0] & 1 == 1 { r = fp_mul(r, base); }
        // exp >>= 1
        let mut ne = [0u64; 4];
        for i in 0..3 { ne[i] = (exp[i] >> 1) | (exp[i+1] << 63); }
        ne[3] = exp[3] >> 1;
        exp = ne;
        if exp == [0u64; 4] { break; }
        base = fp_sqr(base);
    }
    r
}

pub fn fp_inv(a: Fe) -> Fe {
    // p - 2
    let pm2: Fe = [0xFFFFFFFEFFFFFC2D, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF];
    fp_pow(a, pm2)
}

// ─── Scalar arithmetic mod n ──────────────────────────────────────────────────

pub fn sc_add(a: Fe, b: Fe) -> Fe {
    let mut r = [0u64; 4];
    let mut carry = false;
    for i in 0..4 {
        let (s, c1) = a[i].overflowing_add(b[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        r[i] = s; carry = c1 || c2;
    }
    // reduce mod n
    let mut borrow = false;
    if carry || !fe_lt(r, N) {
        let mut s = [0u64; 4];
        for i in 0..4 {
            let (v, b1) = r[i].overflowing_sub(N[i]);
            let (v, b2) = v.overflowing_sub(borrow as u64);
            s[i] = v; borrow = b1 || b2;
        }
        s
    } else { r }
}

// ─── EC Point ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
pub struct Pt { pub x: Fe, pub y: Fe, pub inf: bool }

pub const INF: Pt = Pt { x: [0;4], y: [0;4], inf: true };

pub fn pt_neg(a: Pt) -> Pt {
    if a.inf { a } else { Pt { x: a.x, y: fp_neg(a.y), inf: false } }
}

pub fn pt_add(a: Pt, b: Pt) -> Pt {
    if a.inf { return b; }
    if b.inf { return a; }
    if a.x == b.x {
        if a.y != b.y { return INF; }
        return pt_dbl(a);
    }
    let lam = fp_mul(fp_sub(b.y, a.y), fp_inv(fp_sub(b.x, a.x)));
    let x3 = fp_sub(fp_sub(fp_sqr(lam), a.x), b.x);
    let y3 = fp_sub(fp_mul(lam, fp_sub(a.x, x3)), a.y);
    Pt { x: x3, y: y3, inf: false }
}

fn pt_dbl(a: Pt) -> Pt {
    if a.inf || a.y == [0u64;4] { return INF; }
    let lam = fp_mul(fp_mul([3,0,0,0], fp_sqr(a.x)), fp_inv(fp_add(a.y, a.y)));
    let x3 = fp_sub(fp_sub(fp_sqr(lam), a.x), a.x);
    let y3 = fp_sub(fp_mul(lam, fp_sub(a.x, x3)), a.y);
    Pt { x: x3, y: y3, inf: false }
}

// Double-and-add scalar multiplication
pub fn scalar_mul(mut p: Pt, mut k: Fe) -> Pt {
    let mut r = INF;
    loop {
        if k[0] & 1 == 1 { r = pt_add(r, p); }
        let mut nk = [0u64; 4];
        for i in 0..3 { nk[i] = (k[i] >> 1) | (k[i+1] << 63); }
        nk[3] = k[3] >> 1;
        k = nk;
        if k == [0u64;4] { break; }
        p = pt_dbl(p);
    }
    r
}

pub fn gen() -> Pt { Pt { x: GX, y: GY, inf: false } }

// ─── GLV endomorphism ─────────────────────────────────────────────────────────

/// ψ(x, y) = (β·x mod p, y)
pub fn psi(p: Pt) -> Pt {
    if p.inf { p } else { Pt { x: fp_mul(BETA, p.x), y: p.y, inf: false } }
}

/// ψ²(x, y) = (β²·x mod p, y)
pub fn psi2(p: Pt) -> Pt {
    if p.inf { p } else { Pt { x: fp_mul(BETA2, p.x), y: p.y, inf: false } }
}

// ─── 6-orbit canonical form ───────────────────────────────────────────────────
//
// The 6 automorphisms of secp256k1 (units of End(E) = Z[ζ₃]):
//   { id, ψ, ψ², -id, -ψ, -ψ² }
//
// For each point P, canonical(P) = min x-coordinate over {x, βx, β²x},
// with tie-broken by y parity convention.
// orbit_tag encodes which of the 6 was the canonical representative.

pub struct Canonical {
    pub x: Fe,       // the canonical (minimal) x
    pub tag: u8,     // 0..5: which automorphism maps canonical → P
}

pub fn canonical(p: Pt) -> Canonical {
    if p.inf { return Canonical { x: [0;4], tag: 0 }; }

    let x0 = p.x;
    let x1 = fp_mul(BETA, p.x);
    let x2 = fp_mul(BETA2, p.x);

    // Find min x
    let (min_x, endo) = if fe_lt(x0, x1) {
        if fe_lt(x0, x2) { (x0, 0u8) } else { (x2, 2u8) }
    } else {
        if fe_lt(x1, x2) { (x1, 1u8) } else { (x2, 2u8) }
    };

    // y parity bit: use lowest bit of y
    let y_bit = (p.y[0] & 1) as u8;
    let tag = endo * 2 + y_bit;

    Canonical { x: min_x, tag }
}

// ─── Distinguished point check ────────────────────────────────────────────────

/// Returns true if the canonical x has `dp_bits` trailing zero bits in its most significant limb
pub fn is_dp(canonical_x: Fe, dp_bits: u32) -> bool {
    if dp_bits == 0 { return true; }
    if dp_bits > 64 { return canonical_x[3] == 0 && is_dp([canonical_x[0], canonical_x[1], canonical_x[2], 0], dp_bits - 64); }
    canonical_x[3].trailing_zeros() >= dp_bits
}

// ─── Hex conversion ───────────────────────────────────────────────────────────

pub fn fe_to_hex(a: Fe) -> String {
    format!("{:016x}{:016x}{:016x}{:016x}", a[3], a[2], a[1], a[0])
}

pub fn fe_from_hex(s: &str) -> Option<Fe> {
    let s = s.trim_start_matches("0x");
    if s.len() > 64 { return None; }
    let padded = format!("{:0>64}", s);
    let mut r = [0u64; 4];
    for i in 0..4 {
        let chunk = &padded[(3-i)*16..(3-i)*16+16];
        r[i] = u64::from_str_radix(chunk, 16).ok()?;
    }
    Some(r)
}
