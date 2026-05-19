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
pub const LAMBDA: Fe = [0xDF02967C1B23BD72, 0x122E22EA20816678, 0xA5261C028812645A, 0x5363AD4CC05C30E0];

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
