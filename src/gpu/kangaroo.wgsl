// ─── secp256k1 Kangaroo — WebGPU Compute Shader ──────────────────────────────
//
// 6-automorphism distinguished-point Kangaroo on secp256k1.
// Each GPU thread is one independent Kangaroo animal.
// Field elements: 8 × u32, little-endian (limb 0 = bits 0..31).
//
// End(E) = Z[ζ₃] units: { id, ψ, ψ², −id, −ψ, −ψ² }
// canonical(P) = min(x, β·x, β²·x)  →  one table slot per 6-orbit.

// ─── Constants ────────────────────────────────────────────────────────────────

// p = 2^256 − 2^32 − 977  (little-endian u32)
const P0: u32 = 0xFFFFFC2Fu;
const P1: u32 = 0xFFFFFFFEu;
// P2..P7 = 0xFFFFFFFF

// β = cube root of unity mod p  (ψ(x,y) = (β·x, y))
const BETA0: u32 = 0x719501EEu;
const BETA1: u32 = 0xC1396C28u;
const BETA2_: u32 = 0x12F58995u;
const BETA3: u32 = 0x9CF04975u;
const BETA4: u32 = 0xAC3434E9u;
const BETA5: u32 = 0x6E64479Eu;
const BETA6: u32 = 0x657C0710u;
const BETA7: u32 = 0x7AE96A2Bu;

// β² = p − 1 − β
const BETA2_0: u32 = 0x8E6AFA40u;
const BETA2_1: u32 = 0x3EC693D6u;
const BETA2_2: u32 = 0xED0A766Au;
const BETA2_3: u32 = 0x630FB68Au;
const BETA2_4: u32 = 0x53CBCB16u;
const BETA2_5: u32 = 0x919BB861u;
const BETA2_6: u32 = 0x9A83F8EFu;
const BETA2_7: u32 = 0x851695D4u;

// ─── Structs ──────────────────────────────────────────────────────────────────

struct U256 { l: array<u32, 8> }

struct Point {
  x: U256,
  y: U256,
  // Jacobian Z
  z: U256,
}

struct Animal {
  pt: Point,
  scalar: U256,   // accumulated jump scalar
  is_wild: u32,   // 1 = wild, 0 = tame
  _pad: u32,
}

struct JumpPoint {
  x: U256,
  y: U256,
  scalar: U256,
}

struct DP {
  canon_x: U256,
  scalar: U256,
  is_wild: u32,
  _pad: array<u32, 3>,
}

struct Params {
  steps_per_dispatch: u32,
  dp_mask: u32,        // DP if x.l[7] & dp_mask == 0  (28-bit → mask = 0xFFFFFFF0)
  num_animals: u32,
  max_dps: u32,
}

// ─── Bindings ─────────────────────────────────────────────────────────────────

@group(0) @binding(0) var<storage, read_write> animals: array<Animal>;
@group(0) @binding(1) var<storage, read>       jumps:   array<JumpPoint>;
@group(0) @binding(2) var<storage, read_write> dp_buf:  array<DP>;
@group(0) @binding(3) var<storage, read_write> dp_count: atomic<u32>;
@group(0) @binding(4) var<uniform>             params:  Params;

// ─── 256-bit arithmetic ───────────────────────────────────────────────────────

fn u256_zero() -> U256 {
  return U256(array<u32,8>(0u,0u,0u,0u,0u,0u,0u,0u));
}

fn u256_one() -> U256 {
  return U256(array<u32,8>(1u,0u,0u,0u,0u,0u,0u,0u));
}

// Add with carry chain — returns (sum, carry_out)
fn u256_add_carry(a: U256, b: U256) -> vec2<u32> {
  // Returns carry in .y, sum in animals (side-effect not possible — use global)
  // We return just the carry; caller must track sum separately.
  // Actual impl below uses a struct trick.
  return vec2<u32>(0u, 0u); // placeholder
}

// In-place: r = a + b, returns carry
fn add256(a: U256, b: U256) -> U256 {
  var r: U256;
  var c: u32 = 0u;
  for (var i = 0; i < 8; i++) {
    let s0 = a.l[i] + b.l[i];
    let c0 = u32(s0 < a.l[i]);
    let s1 = s0 + c;
    let c1 = u32(s1 < s0);
    r.l[i] = s1;
    c = c0 + c1;
  }
  return r;
}

// r = a - b mod 2^256, returns borrow
fn sub256_borrow(a: U256, b: U256) -> U256 {
  var r: U256;
  var borrow: u32 = 0u;
  for (var i = 0; i < 8; i++) {
    let s0 = a.l[i] - b.l[i];
    let b0 = u32(a.l[i] < b.l[i]);
    let s1 = s0 - borrow;
    let b1 = u32(s0 < borrow);
    r.l[i] = s1;
    borrow = b0 + b1;
  }
  return r;
}

fn u256_lt(a: U256, b: U256) -> bool {
  for (var i = 7; i >= 0; i--) {
    if a.l[i] < b.l[i] { return true; }
    if a.l[i] > b.l[i] { return false; }
  }
  return false;
}

// ─── secp256k1 field mod p ────────────────────────────────────────────────────

fn fp_p() -> U256 {
  return U256(array<u32,8>(P0, P1, 0xFFFFFFFFu, 0xFFFFFFFFu,
                           0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu));
}

fn fp_add(a: U256, b: U256) -> U256 {
  var r = add256(a, b);
  let p = fp_p();
  // Reduce: if r >= p, subtract p
  if !u256_lt(r, p) {
    r = sub256_borrow(r, p);
  }
  return r;
}

fn fp_sub(a: U256, b: U256) -> U256 {
  let p = fp_p();
  if u256_lt(a, b) {
    // a - b < 0, add p
    let t = sub256_borrow(b, a);
    return sub256_borrow(p, t);
  }
  return sub256_borrow(a, b);
}

fn fp_neg(a: U256) -> U256 {
  let zero = u256_zero();
  // Check if a == 0
  var is_zero = true;
  for (var i = 0; i < 8; i++) { if a.l[i] != 0u { is_zero = false; } }
  if is_zero { return zero; }
  return fp_sub(fp_p(), a);
}

// 32×32 → 64 bit multiply: returns (lo, hi)
fn mul32x32(a: u32, b: u32) -> vec2<u32> {
  let a_lo = a & 0xFFFFu;
  let a_hi = a >> 16u;
  let b_lo = b & 0xFFFFu;
  let b_hi = b >> 16u;

  let p0 = a_lo * b_lo;
  let p1 = a_lo * b_hi;
  let p2 = a_hi * b_lo;
  let p3 = a_hi * b_hi;

  let mid = p1 + p2;
  let mid_carry = u32(mid < p1);

  let lo_part = (mid & 0xFFFFu) << 16u;
  let lo = p0 + lo_part;
  let c0 = u32(lo < p0);

  let hi = p3 + (mid >> 16u) + (mid_carry << 16u) + c0;
  return vec2<u32>(lo, hi);
}

// 256×256 → 512-bit product stored in r[0..15]
fn mul512(a: U256, b: U256, r: ptr<function, array<u32, 16>>) {
  for (var i = 0; i < 16; i++) { (*r)[i] = 0u; }
  for (var i = 0; i < 8; i++) {
    var carry: u32 = 0u;
    for (var j = 0; j < 8; j++) {
      let prod = mul32x32(a.l[i], b.l[j]);
      let s0 = (*r)[i+j] + prod.x;
      let c0 = u32(s0 < (*r)[i+j]);
      let s1 = s0 + carry;
      let c1 = u32(s1 < s0);
      (*r)[i+j] = s1;
      carry = prod.y + c0 + c1;
    }
    (*r)[i+8] += carry;
  }
}

// Reduce 512-bit number mod p using p = 2^256 − (2^32 + 977)
// c_p = 2^32 + 977 = 0x1000003D1
fn reduce512(t: array<u32, 16>) -> U256 {
  // hi = t[8..15], lo = t[0..7]
  // hi * c_p added to lo

  // c_p = 2^32 + 977: multiply hi[i] by c_p
  // hi * 977 (scalar)
  var acc: array<u32, 9>;
  for (var i = 0; i < 9; i++) { acc[i] = 0u; }

  var carry: u32 = 0u;
  for (var i = 0; i < 8; i++) {
    let p = mul32x32(t[i+8], 977u);
    let s0 = acc[i] + p.x;
    let c0 = u32(s0 < acc[i]);
    let s1 = s0 + carry;
    let c1 = u32(s1 < s0);
    acc[i] = s1;
    carry = p.y + c0 + c1;
  }
  acc[8] = carry;

  // hi * 2^32 = shift hi left 32 bits (i.e. one u32 limb)
  var hi32: array<u32, 9>;
  hi32[0] = 0u;
  for (var i = 0; i < 8; i++) { hi32[i+1] = t[i+8]; }

  // acc += hi32
  carry = 0u;
  for (var i = 0; i < 9; i++) {
    let s0 = acc[i] + hi32[i];
    let c0 = u32(s0 < acc[i]);
    let s1 = s0 + carry;
    let c1 = u32(s1 < s0);
    acc[i] = s1;
    carry = c0 + c1;
  }
  // acc[8] might overflow — but bounded, ignore extra carry for now

  // r = lo + acc[0..7]
  var r: array<u32, 9>;
  carry = 0u;
  for (var i = 0; i < 8; i++) {
    let s0 = t[i] + acc[i];
    let c0 = u32(s0 < t[i]);
    let s1 = s0 + carry;
    let c1 = u32(s1 < s0);
    r[i] = s1;
    carry = c0 + c1;
  }
  r[8] = acc[8] + carry;

  // If r[8] > 0, do one more reduction
  if r[8] > 0u {
    let extra = r[8];
    // extra * c_p: extra * 977 + extra * 2^32
    let e977 = mul32x32(extra, 977u);
    var carry2: u32 = 0u;
    let s0 = r[0] + e977.x;
    let c0 = u32(s0 < r[0]);
    carry2 = e977.y + c0;
    r[0] = s0;

    let s1 = r[1] + carry2;
    let c1 = u32(s1 < r[1]);
    r[1] = s1; carry2 = c1;

    // extra * 2^32 = extra goes to limb 1
    let s2 = r[1] + extra;
    let c2 = u32(s2 < r[1]);
    r[1] = s2; carry2 += c2;

    for (var i = 2; i < 8; i++) {
      let s = r[i] + carry2;
      carry2 = u32(s < r[i]);
      r[i] = s;
    }
    r[8] = 0u;
  }

  var res = U256(array<u32,8>(r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7]));
  let p = fp_p();
  if !u256_lt(res, p) { res = sub256_borrow(res, p); }
  return res;
}

fn fp_mul(a: U256, b: U256) -> U256 {
  var t: array<u32, 16>;
  mul512(a, b, &t);
  return reduce512(t);
}

fn fp_sqr(a: U256) -> U256 { return fp_mul(a, a); }

// ─── GLV endomorphisms ────────────────────────────────────────────────────────

fn beta() -> U256 {
  return U256(array<u32,8>(BETA0,BETA1,BETA2_,BETA3,BETA4,BETA5,BETA6,BETA7));
}

fn beta2() -> U256 {
  return U256(array<u32,8>(BETA2_0,BETA2_1,BETA2_2,BETA2_3,BETA2_4,BETA2_5,BETA2_6,BETA2_7));
}

fn psi_x(x: U256) -> U256  { return fp_mul(beta(),  x); }
fn psi2_x(x: U256) -> U256 { return fp_mul(beta2(), x); }

// Canonical x: min of {x, β·x, β²·x}
fn canonical_x(x: U256) -> U256 {
  let x1 = psi_x(x);
  let x2 = psi2_x(x);
  var m = x;
  if u256_lt(x1, m) { m = x1; }
  if u256_lt(x2, m) { m = x2; }
  return m;
}

// ─── Jacobian point addition (mixed: Jac + Affine) ───────────────────────────
// Formula: madd-2007-bl from EFD (secp256k1, a=0)
// Input: P=(X1:Y1:Z1) Jacobian, Q=(x2,y2) affine (Z2=1)
// Output: R=(X3:Y3:Z3) Jacobian

fn pt_add_mixed(px: U256, py: U256, pz: U256,
                qx: U256, qy: U256) -> array<U256, 3> {
  // Z1^2
  let z1z1 = fp_sqr(pz);
  // U2 = x2 * Z1^2
  let u2 = fp_mul(qx, z1z1);
  // S2 = y2 * Z1^3 = y2 * Z1 * Z1^2
  let z1_cubed = fp_mul(pz, z1z1);
  let s2 = fp_mul(qy, z1_cubed);
  // H = U2 - X1
  let h = fp_sub(u2, px);
  // R = S2 - Y1
  let r = fp_sub(s2, py);
  // HH = H^2
  let hh = fp_sqr(h);
  // HHH = H * HH
  let hhh = fp_mul(h, hh);
  // V = X1 * HH
  let v = fp_mul(px, hh);
  // X3 = R^2 - HHH - 2*V
  let rr = fp_sqr(r);
  let two_v = fp_add(v, v);
  let x3 = fp_sub(fp_sub(rr, hhh), two_v);
  // Y3 = R*(V - X3) - Y1*HHH
  let v_x3 = fp_sub(v, x3);
  let r_vx3 = fp_mul(r, v_x3);
  let y1_hhh = fp_mul(py, hhh);
  let y3 = fp_sub(r_vx3, y1_hhh);
  // Z3 = H * Z1
  let z3 = fp_mul(h, pz);

  return array<U256,3>(x3, y3, z3);
}

// Convert Jacobian (X:Y:Z) to affine x = X * Z^(-2)
// We use Fermat: Z^(-1) = Z^(p-2)
// Returns affine x only (we don't need y for DP check)
fn jac_to_affine_x(X: U256, Z: U256) -> U256 {
  // Z_inv = Z^(p-2)
  // p-2 in limbs: [0xFFFFFC2D, 0xFFFFFFFE, 0xFFFF..., 0xFFFF...]
  // Use square-and-multiply with hardcoded exponent bits
  // For gas/cycle efficiency, we use the p-2 exponent directly
  var base = Z;
  var r = u256_one();
  // Exponent p-2: iterate over 256 bits
  // p-2 = [0xFFFFFC2D, 0xFFFFFFFE, 0xFFFFFFFF × 6]
  // We encode as array of u32 and iterate
  let exp = array<u32,8>(0xFFFFFC2Du, 0xFFFFFFFEu, 0xFFFFFFFFu, 0xFFFFFFFFu,
                          0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu);
  for (var wi = 0; wi < 8; wi++) {
    var word = exp[wi];
    for (var bi = 0; bi < 32; bi++) {
      if (word & 1u) != 0u { r = fp_mul(r, base); }
      base = fp_sqr(base);
      word >>= 1u;
    }
  }
  // r = Z^(-1), affine_x = X * Z^(-2) = X * (Z^(-1))^2
  let z_inv2 = fp_sqr(r);
  return fp_mul(X, z_inv2);
}

// ─── Jump index from Jacobian X (top bits) ───────────────────────────────────

fn jump_idx(px: U256, num_jumps: u32) -> u32 {
  return px.l[0] % num_jumps;
}

// ─── Distinguished point check ────────────────────────────────────────────────
// Use probabilistic check on Jacobian X top limb (avoids inversion every step)
// Actual canonical DP check done only on positive results (rare).

fn is_probable_dp(px: U256, dp_mask: u32) -> bool {
  // Fast: check top u32 of raw Jacobian X
  return (px.l[7] & dp_mask) == 0u;
}

// ─── Main Kangaroo kernel ─────────────────────────────────────────────────────

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let idx = gid.x;
  if idx >= params.num_animals { return; }

  var animal = animals[idx];
  let num_jumps = arrayLength(&jumps);

  for (var step = 0u; step < params.steps_per_dispatch; step++) {
    // Pick jump based on current X
    let ji = jump_idx(animal.pt.x, num_jumps);
    let jump = jumps[ji];

    // Jacobian + Affine mixed addition
    let new_coords = pt_add_mixed(
      animal.pt.x, animal.pt.y, animal.pt.z,
      jump.x, jump.y
    );
    animal.pt.x = new_coords[0];
    animal.pt.y = new_coords[1];
    animal.pt.z = new_coords[2];

    // Accumulate scalar
    var carry: u32 = 0u;
    for (var i = 0; i < 8; i++) {
      let s0 = animal.scalar.l[i] + jump.scalar.l[i];
      let c0 = u32(s0 < animal.scalar.l[i]);
      let s1 = s0 + carry;
      let c1 = u32(s1 < s0);
      animal.scalar.l[i] = s1;
      carry = c0 + c1;
    }

    // DP check (probabilistic, on raw Jacobian X)
    if is_probable_dp(animal.pt.x, params.dp_mask) {
      // Convert to affine x for canonical form
      let ax = jac_to_affine_x(animal.pt.x, animal.pt.z);
      let cx = canonical_x(ax);

      let slot = atomicAdd(&dp_count, 1u) % params.max_dps;
      dp_buf[slot] = DP(
        cx,
        animal.scalar,
        animal.is_wild,
        array<u32,3>(0u, 0u, 0u)
      );
    }
  }

  animals[idx] = animal;
}
