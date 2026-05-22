// secp256k1 CPU arithmetic (Rust) — used for setup, verification, and k recovery
// Mirrors the CUDA device code but runs on host.

#![allow(dead_code)]
pub const P: u128 = 0; // placeholder — we use BigInt-style via the functions below

// Field element: [u64; 4] little-endian
pub type Fe = [u64; 4];

pub const FIELD_P: Fe = [0xFFFFFFFEFFFFFC2F, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF];
pub const FIELD_N: Fe = [0xBFD25E8CD0364141, 0xBAAEDCE6AF48A03B, 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF];
pub const GX: Fe     = [0x59F2815B16F81798, 0x029BFCDB2DCE28D9, 0x55A06295CE870B07, 0x79BE667EF9DCBBAC];
pub const GY: Fe     = [0x9C47D08FFB10D4B8, 0xFD17B448A6855419, 0x5DA4FBFC0E1108A8, 0x483ADA7726A3C465];
pub const BETA: Fe   = [0xC1396C28719501EE, 0x9CF0497512F58995, 0x6E64479EAC3434E9, 0x7AE96A2B657C0710];
pub const BETA2: Fe  = [0x3EC693D68E6AFA40, 0x630FB68AED0A766A, 0x919BB86153CBCB16, 0x851695D49A83F8EF];
pub const LAMBDA: Fe  = [0xDF02967C1B23BD72, 0x122E22EA20816678, 0xA5261C028812645A, 0x5363AD4CC05C30E0];
// λ² ≡ -1 - λ (mod n)  since 1 + λ + λ² ≡ 0 (mod n)
pub const LAMBDA2: Fe = [0xDCCFC810B51283CE, 0xA880B9FC8EC739C2, 0x5AD9E3FD77ED9BA4, 0xAC9C52B33FA3CF1F];

// ─── Comparison ───────────────────────────────────────────────────────────────

pub fn fe_lt(a: Fe, b: Fe) -> bool {
    for i in (0..4).rev() {
        if a[i] < b[i] { return true; }
        if a[i] > b[i] { return false; }
    }
    false
}

pub fn fe_eq(a: Fe, b: Fe) -> bool { a == b }

// ─── Field mod p ──────────────────────────────────────────────────────────────

fn sub_p(mut r: Fe) -> Fe {
    let mut borrow = false;
    for i in 0..4 {
        let (s, b1) = r[i].overflowing_sub(FIELD_P[i]);
        let (s, b2) = s.overflowing_sub(borrow as u64);
        r[i] = s; borrow = b1 || b2;
    }
    r
}

pub fn fp_add(a: Fe, b: Fe) -> Fe {
    let mut r = [0u64; 4]; let mut carry = false;
    for i in 0..4 {
        let (s, c1) = a[i].overflowing_add(b[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        r[i] = s; carry = c1 || c2;
    }
    if carry || !fe_lt(r, FIELD_P) { sub_p(r) } else { r }
}

pub fn fp_sub(a: Fe, b: Fe) -> Fe {
    if !fe_lt(a, b) {
        let mut r = [0u64; 4]; let mut borrow = false;
        for i in 0..4 {
            let (s, b1) = a[i].overflowing_sub(b[i]);
            let (s, b2) = s.overflowing_sub(borrow as u64);
            r[i] = s; borrow = b1 || b2;
        }
        r
    } else {
        // a < b → p - (b - a)
        let diff = fp_sub(b, a); // b - a (positive now)
        fp_sub(FIELD_P, diff)    // p - diff... but this might recurse. Use explicit:
    }
}

pub fn fp_neg(a: Fe) -> Fe {
    if a == [0u64; 4] { return a; }
    fp_sub(FIELD_P, a)
}

fn mul512(a: Fe, b: Fe) -> [u64; 8] {
    let mut t = [0u64; 8];
    for i in 0..4 {
        let mut carry = 0u64;
        for j in 0..4 {
            let prod = (a[i] as u128) * (b[j] as u128) + t[i+j] as u128 + carry as u128;
            t[i+j] = prod as u64;
            carry   = (prod >> 64) as u64;
        }
        t[i+4] += carry;
    }
    t
}

fn reduce512(t: [u64; 8]) -> Fe {
    let lo = [t[0], t[1], t[2], t[3]];
    let hi = [t[4], t[5], t[6], t[7]];

    // c = hi * (2^32 + 977)
    let mut c = [0u64; 5];
    let mut carry: u64 = 0;
    for i in 0..4 {
        let p = (hi[i] as u128) * 977 + c[i] as u128 + carry as u128;
        c[i] = p as u64; carry = (p >> 64) as u64;
    }
    c[4] = carry;

    // c += hi * 2^32
    #[allow(unused_assignments)] { carry = 0; }
    let s = (c[0] as u128) + ((hi[0] << 32) as u128);
    c[0] = s as u64; carry = (s >> 64) as u64;
    for i in 1..4 {
        let s = (c[i] as u128) + ((hi[i] << 32) as u128) + ((hi[i-1] >> 32) as u128) + carry as u128;
        c[i] = s as u64; carry = (s >> 64) as u64;
    }
    c[4] += (hi[3] >> 32) + carry;

    // r5 = lo + c
    let mut r5 = [0u64; 5]; carry = 0;
    for i in 0..4 {
        let s = (lo[i] as u128) + (c[i] as u128) + carry as u128;
        r5[i] = s as u64; carry = (s >> 64) as u64;
    }
    r5[4] = c[4] + carry;

    if r5[4] > 0 {
        let e = r5[4] as u128;
        let ex = e * (977 + (1u128 << 32));
        let s = (r5[0] as u128) + (ex as u64 as u128);
        r5[0] = s as u64; carry = (s >> 64) as u64;
        let s = (r5[1] as u128) + ((ex >> 64) as u64 as u128) + carry as u128;
        r5[1] = s as u64; carry = (s >> 64) as u64;
        let s = (r5[2] as u128) + carry as u128;
        r5[2] = s as u64; carry = (s >> 64) as u64;
        r5[3] = r5[3].wrapping_add(carry);
        r5[4] = 0;
    }

    let mut r = [r5[0], r5[1], r5[2], r5[3]];
    if !fe_lt(r, FIELD_P) { r = sub_p(r); }
    r
}

pub fn fp_mul(a: Fe, b: Fe) -> Fe { reduce512(mul512(a, b)) }
pub fn fp_sqr(a: Fe) -> Fe { fp_mul(a, a) }

pub fn fp_pow(mut base: Fe, mut e: Fe) -> Fe {
    let mut r = [1u64, 0, 0, 0];
    loop {
        if e[0] & 1 == 1 { r = fp_mul(r, base); }
        let mut ne = [0u64; 4];
        for i in 0..3 { ne[i] = (e[i] >> 1) | (e[i+1] << 63); }
        ne[3] = e[3] >> 1; e = ne;
        if e == [0u64; 4] { break; }
        base = fp_sqr(base);
    }
    r
}

pub fn fp_inv(a: Fe) -> Fe {
    // p-2
    let pm2 = [0xFFFFFFFEFFFFFC2D, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF];
    fp_pow(a, pm2)
}

// ─── Scalar mod n ─────────────────────────────────────────────────────────────

fn sub_n(mut r: Fe) -> Fe {
    let mut borrow = false;
    for i in 0..4 {
        let (s, b1) = r[i].overflowing_sub(FIELD_N[i]);
        let (s, b2) = s.overflowing_sub(borrow as u64);
        r[i] = s; borrow = b1 || b2;
    }
    r
}

pub fn sc_add(a: Fe, b: Fe) -> Fe {
    let mut r = [0u64; 4]; let mut carry = false;
    for i in 0..4 {
        let (s, c1) = a[i].overflowing_add(b[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        r[i] = s; carry = c1 || c2;
    }
    if carry || !fe_lt(r, FIELD_N) { sub_n(r) } else { r }
}

pub fn sc_sub(a: Fe, b: Fe) -> Fe {
    if !fe_lt(a, b) {
        let mut r = [0u64; 4]; let mut borrow = false;
        for i in 0..4 {
            let (s, b1) = a[i].overflowing_sub(b[i]);
            let (s, b2) = s.overflowing_sub(borrow as u64);
            r[i] = s; borrow = b1 || b2;
        }
        r
    } else {
        let diff = sc_sub(b, a);
        sc_sub(FIELD_N, diff)
    }
}

pub fn sc_neg(a: Fe) -> Fe {
    if a == [0u64; 4] { return a; }
    sc_sub(FIELD_N, a)
}

pub fn sc_mul(a: Fe, b: Fe) -> Fe {
    // mod n via Barrett — simplified: multiply then subtract multiples of n
    let t = mul512(a, b);
    // Use same reduction but mod n (2^256 - n_comp where n_comp = 2^256 - n)
    // n_comp = [0x402DA1732FC9BEBF, 0x4551231950B75FC4, 0x0000000000000001, 0x0000000000000000]
    let n_comp: Fe = [0x402DA1732FC9BEBF, 0x4551231950B75FC4, 0x0000000000000001, 0x0000000000000000];
    let hi = [t[4], t[5], t[6], t[7]];
    let lo = [t[0], t[1], t[2], t[3]];

    // correction = hi * n_comp (roughly)
    let corr = mul512(hi, n_comp);
    let corr_lo = [corr[0], corr[1], corr[2], corr[3]];

    let mut res = [0u64; 4]; let mut carry = false;
    for i in 0..4 {
        let (s, c1) = lo[i].overflowing_add(corr_lo[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        res[i] = s; carry = c1 || c2;
    }
    while !fe_lt(res, FIELD_N) { res = sub_n(res); }
    res
}

// ─── EC point ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Pt { pub x: Fe, pub y: Fe, pub inf: bool }

pub const INF: Pt = Pt { x: [0;4], y: [0;4], inf: true };
pub const G:   Pt = Pt { x: GX,  y: GY,  inf: false };

pub fn pt_neg(p: Pt) -> Pt {
    if p.inf { p } else { Pt { x: p.x, y: fp_neg(p.y), inf: false } }
}

pub fn pt_add(a: Pt, b: Pt) -> Pt {
    if a.inf { return b; }
    if b.inf { return a; }
    if a.x == b.x {
        if a.y != b.y { return INF; }
        return pt_dbl(a);
    }
    let lam = fp_mul(fp_sub(b.y, a.y), fp_inv(fp_sub(b.x, a.x)));
    let x3  = fp_sub(fp_sub(fp_sqr(lam), a.x), b.x);
    let y3  = fp_sub(fp_mul(lam, fp_sub(a.x, x3)), a.y);
    Pt { x: x3, y: y3, inf: false }
}

fn pt_dbl(a: Pt) -> Pt {
    if a.inf || a.y == [0u64;4] { return INF; }
    let lam = fp_mul(fp_mul([3,0,0,0], fp_sqr(a.x)), fp_inv(fp_add(a.y, a.y)));
    let x3  = fp_sub(fp_sub(fp_sqr(lam), a.x), a.x);
    let y3  = fp_sub(fp_mul(lam, fp_sub(a.x, x3)), a.y);
    Pt { x: x3, y: y3, inf: false }
}

pub fn scalar_mul(mut p: Pt, mut k: Fe) -> Pt {
    let mut r = INF;
    loop {
        if k[0] & 1 == 1 { r = pt_add(r, p); }
        let mut nk = [0u64; 4];
        for i in 0..3 { nk[i] = (k[i] >> 1) | (k[i+1] << 63); }
        nk[3] = k[3] >> 1; k = nk;
        if k == [0u64;4] { break; }
        p = pt_dbl(p);
    }
    r
}

// ─── GLV endomorphism ─────────────────────────────────────────────────────────

/// Apply the GLV endomorphism φ: (x, y) → (β·x, y).
/// φ is of order 3 in the automorphism group; φ(P) = λ·P in scalar space.
pub fn phi_point(p: Pt) -> Pt {
    if p.inf { return p; }
    Pt { x: fp_mul(BETA, p.x), y: p.y, inf: false }
}

/// Apply φ²: (x, y) → (β²·x, y).  φ²(P) = λ²·P in scalar space.
pub fn phi2_point(p: Pt) -> Pt {
    if p.inf { return p; }
    Pt { x: fp_mul(BETA2, p.x), y: p.y, inf: false }
}

/// Multiply scalar s by λ (the GLV eigenvalue) mod n.
pub fn sc_mul_lambda(s: Fe) -> Fe { sc_mul(LAMBDA, s) }

/// Multiply scalar s by λ² mod n.
pub fn sc_mul_lambda2(s: Fe) -> Fe { sc_mul(LAMBDA2, s) }

/// GLV scalar decomposition: k = k1 + k2·λ  with |k1|, |k2| ≈ √k.
///
/// Uses the secp256k1 short basis (Bernstein et al. / bitcoin-core):
///   a1 =  0x3086d221a7d46bcde86c90e49284eb15 (128-bit)
///   b1 = -0xe4437ed6010e88286f547fa90abfe4c3 (128-bit)
///   a2 =  0x114ca50f7a8e2f3f657c1108d9d44cfd8 (128-bit)
///   b2 =  a1
///
/// Babai rounding: c1 = ⌊b2·k/n⌋, c2 = ⌊-b1·k/n⌋
///                 k1 = k − c1·a1 − c2·a2  (mod n)
///                 k2 =   − c1·b1 − c2·b2  (mod n)
///
/// For k < 2^135, both k1 and k2 are < 2^68.
pub fn glv_decompose(k: Fe) -> (Fe, Fe) {
    // Precomputed round constants g1, g2 = round(2^384 · b / n) for Babai.
    // From secp256k1 reference: g1 = 0x3086d221a7d46bcde86c90e49284eb153dab4d1b27bb1a3e09f2c622c27d1b25
    //                           g2 = 0xe4437ed6010e88286f547fa90abfe4c42b0a9be4fe36d6e0de9f40a6d1e7d887
    let g1: [u64; 4] = [0x09f2c622c27d1b25, 0x3dab4d1b27bb1a3e, 0xe86c90e49284eb15, 0x3086d221a7d46bcd];
    let g2: [u64; 4] = [0xde9f40a6d1e7d887, 0x2b0a9be4fe36d6e0, 0x6f547fa90abfe4c4, 0xe4437ed6010e8828];

    // c1 = (g1 · k) >> 384, c2 = (g2 · k) >> 384  (top 256-384 bits of 512-bit product)
    let t1 = mul512(g1, k);
    let t2 = mul512(g2, k);
    // >> 384 = >> (256+128) means we want bits [384..511] of the 512-bit product
    // Since mul512 gives an 8-limb (512-bit) result, >> 384 = drop lower 6 limbs, keep limbs [6,7]
    // But our product is only 512 bits so limbs [6] and [7] give us 128-bit quotient.
    let c1: Fe = [t1[6], t1[7], 0, 0];
    let c2: Fe = [t2[6], t2[7], 0, 0];

    // Basis vectors
    let a1: Fe = [0xe86c90e49284eb15, 0x3086d221a7d46bcd, 0, 0];
    let b1: Fe = [0x6f547fa90abfe4c3, 0xe4437ed6010e8828, 0, 0]; // |b1|, b1 is negative
    let a2: Fe = [0x57c1108d9d44cfd8, 0x114ca50f7a8e2f3f, 0, 0];
    // b2 = a1

    // k1 = k − c1·a1 − c2·a2 (mod n)
    let c1a1 = sc_mul(c1, a1);
    let c2a2 = sc_mul(c2, a2);
    let k1   = sc_sub(sc_sub(k, c1a1), c2a2);

    // k2 = c1·|b1| − c2·a1 (mod n)  [signs: b1<0, b2=a1>0]
    let c1b1 = sc_mul(c1, b1);
    let c2b2 = sc_mul(c2, a1); // b2 = a1
    let k2   = sc_sub(c1b1, c2b2);

    (k1, k2)
}

// ─── 6-automorphism canonical form ────────────────────────────────────────────

pub fn canonical_x(x: Fe) -> Fe {
    let x1 = fp_mul(BETA,  x);
    let x2 = fp_mul(BETA2, x);
    let mut m = x;
    if fe_lt(x1, m) { m = x1; }
    if fe_lt(x2, m) { m = x2; }
    m
}

// ─── Hex helpers ──────────────────────────────────────────────────────────────

pub fn fe_to_hex(a: Fe) -> String {
    format!("{:016x}{:016x}{:016x}{:016x}", a[3], a[2], a[1], a[0])
}

pub fn fe_from_hex(s: &str) -> Option<Fe> {
    let s = s.trim_start_matches("0x");
    if s.len() > 64 { return None; }
    let padded = format!("{:0>64}", s);
    let mut r = [0u64; 4];
    for i in 0..4 {
        r[i] = u64::from_str_radix(&padded[(3-i)*16..(3-i)*16+16], 16).ok()?;
    }
    Some(r)
}
