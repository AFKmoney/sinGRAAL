// secp256k1 GLV — exact 6-automorphism k recovery
//
// The 6-automorphism group of secp256k1 is {±id, ±ψ, ±ψ²} where
// ψ(x,y) = (β·x, y) acts as scalar multiplication by λ.
//
// In scalar space the 6 automorphisms are: {1, n−1, λ, n−λ, λ², n−λ²}.
//
// When a tame animal with accumulated scalar T and a wild animal with
// accumulated scalar W land on the same canonical x, some automorphism
// α satisfies:  T·G = α(W·G + P)
//               T·G = α_scalar·W·G + α_scalar·P
//             (T − α_scalar·W)·G = α_scalar·P
//
// Since P = k·G: k = (T − α_scalar·W) · α_scalar⁻¹  ... but more simply:
// rearranging gives k = (T·α_scalar⁻¹ − W)  or  k = α_scalar·T − W.
// We just try all 6 automorphisms and verify.

#![allow(dead_code)]

use crate::secp::{Fe, FIELD_N, sc_mul, sc_sub, scalar_mul, G};

// λ = GLV eigenvalue: ψ(P) = λ·P in the scalar field
pub const LAMBDA: Fe = [
    0xDF02967C1B23BD72, 0x122E22EA20816678,
    0xA5261C028812645A, 0x5363AD4CC05C30E0,
];

// ─── Range check ──────────────────────────────────────────────────────────────
//
// Returns true iff k < 2^range_bits.
// For puzzle #135, valid keys satisfy k < 2^135, so range_bits = 135.
// Of the 6 automorphism candidates, exactly 1 (the correct one) lands in
// [0, 2^range_bits). The others reduce mod n to values ≈ n ≈ 2^256, which
// fail this check — saving 5 of 6 expensive scalar_mul calls per collision.
fn in_range(k: Fe, range_bits: u32) -> bool {
    if range_bits == 0 || range_bits >= 256 { return true; }
    let word = (range_bits / 64) as usize;
    let bit  = range_bits % 64;
    if word < 4 {
        // Top word: all bits ≥ `bit` must be zero
        let hi_mask = if bit == 0 { !0u64 } else { !((1u64 << bit) - 1) };
        if k[word] & hi_mask != 0 { return false; }
        // All words above `word` must be zero
        for i in (word + 1)..4 {
            if k[i] != 0 { return false; }
        }
    }
    true
}

// ─── Exact 6-automorphism k recovery ─────────────────────────────────────────
//
// Given a tame/wild DP collision on canonical x:
//   tame_sc·G  =  some canonical representative  Q
//   P + wild_sc·G  =  aut(Q)  for some aut in the 6-group
//
// Trying every automorphism α: k = α·tame_sc − wild_sc
// (6 candidates; we verify each by computing k·G and comparing with target)
//
// range_bits: expected bit-length of the answer (0 = no filter).
// When range_bits > 0, candidates outside [0, 2^range_bits) are skipped
// before the expensive scalar_mul — typically eliminates 5 of 6 candidates.

pub fn recover_k_6aut(
    tame_sc: Fe,
    wild_sc: Fe,
    target_x: Fe,
    target_y: Fe,
    range_bits: u32,
) -> Vec<Fe> {
    let lam2  = sc_mul(LAMBDA, LAMBDA);
    let neg1  = sc_sub(FIELD_N, [1, 0, 0, 0]);
    let negl  = sc_sub(FIELD_N, LAMBDA);
    let negl2 = sc_sub(FIELD_N, lam2);

    let automorphisms: [Fe; 6] = [
        [1, 0, 0, 0],   // +id
        neg1,           // −id
        LAMBDA,         // +ψ
        negl,           // −ψ
        lam2,           // +ψ²
        negl2,          // −ψ²
    ];

    let mut results = Vec::new();
    for &alpha in &automorphisms {
        let k = sc_sub(sc_mul(alpha, tame_sc), wild_sc);
        // Skip the expensive EC multiply if k is clearly out of range.
        // The correct k ∈ [0, 2^range_bits); wrong candidates reduce mod n
        // to values near n ≈ 2^256 — this filter eliminates them in O(1).
        if !in_range(k, range_bits) { continue; }
        let computed = scalar_mul(G, k);
        if !computed.inf
            && computed.x == target_x
            && computed.y == target_y
        {
            results.push(k);
        }
    }
    results
}
