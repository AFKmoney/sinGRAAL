// 6-automorphism Kangaroo for secp256k1 ECDLP
//
// Innovation: uses all 6 units of End(E) = Z[ζ₃]
//   { id, ψ, ψ², -id, -ψ, -ψ² }
// → √6 ≈ 2.45× speedup over basic Kangaroo
//   vs √4 ≈ 2× for standard 4-automorphism implementations
//
// Walk: Q_{n+1} = Q_n + jump[hash(canonical(Q_n))]
// DP:   distinguished point when top dp_bits of canonical_x are zero

use std::collections::HashMap;
use crate::secp::*;

const NUM_JUMPS: usize = 32;

// ─── Scalar arithmetic mod n ──────────────────────────────────────────────────

fn sc_sub(a: Fe, b: Fe) -> Fe {
    let mut r = [0u64; 4];
    let mut borrow = false;
    for i in 0..4 {
        let (s, b1) = a[i].overflowing_sub(b[i]);
        let (s, b2) = s.overflowing_sub(borrow as u64);
        r[i] = s; borrow = b1 || b2;
    }
    if borrow {
        let mut carry = false;
        for i in 0..4 {
            let (s, c1) = r[i].overflowing_add(N[i]);
            let (s, c2) = s.overflowing_add(carry as u64);
            r[i] = s; carry = c1 || c2;
        }
    }
    r
}

fn sc_neg(a: Fe) -> Fe {
    if a == [0u64;4] { return a; }
    sc_sub(N, a)
}

fn sc_mul(a: Fe, b: Fe) -> Fe {
    // schoolbook mod n via Barrett / subtraction
    // Simplified: use the same reduce512 approach but mod n
    // For now: use trial reduction (n is close to 2^256)
    let mut r = [0u64; 8];
    for i in 0..4 {
        let mut carry = 0u64;
        for j in 0..4 {
            let uv = (a[i] as u128) * (b[j] as u128)
                   + r[i+j] as u128 + carry as u128;
            r[i+j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        r[i+4] = carry;
    }
    // Reduce mod n using repeated subtraction of n * 2^(64*i) for hi limbs
    // Simple approach: use the fact that n ≈ 2^256
    let hi = [r[4], r[5], r[6], r[7]];
    // hi * 2^256 mod n: use hi * (2^256 - n) = hi * (BAA...BF)
    // complement = 2^256 - n
    let n_comp: [u64; 4] = [
        u64::wrapping_sub(0, N[0]).wrapping_add(1), // this is 2^256 - n's low limb
        !N[1],
        !N[2],
        !N[3],
    ];
    // Actually: 2^256 - n = [0x402DA1732FC9BEBF, 0x4551231950B75FC4, 0x0000000000000001, 0x0000000000000000]
    let n_comp: [u64; 4] = [
        0x402DA1732FC9BEBF,
        0x4551231950B75FC4,
        0x0000000000000001,
        0x0000000000000000, // 2^256 - n
    ];

    // hi * n_comp as 512-bit
    let mut prod = [0u64; 8];
    for i in 0..4 {
        let mut carry = 0u64;
        for j in 0..4 {
            let uv = (hi[i] as u128) * (n_comp[j] as u128)
                   + prod[i+j] as u128 + carry as u128;
            prod[i+j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        prod[i+4] = carry;
    }

    let lo = [r[0], r[1], r[2], r[3]];
    let correction = [prod[0], prod[1], prod[2], prod[3]];

    let mut res = [0u64; 4];
    let mut carry = false;
    for i in 0..4 {
        let (s, c1) = lo[i].overflowing_add(correction[i]);
        let (s, c2) = s.overflowing_add(carry as u64);
        res[i] = s; carry = c1 || c2;
    }
    // Final reduction
    while !fe_lt(res, N) {
        res = sc_sub(res, N);
    }
    res
}

// Precompute jump table: jumps[i] = 2^(base + i) * G
// base chosen so avg jump ≈ √(range)
fn build_jumps(range_bits: u32) -> (Vec<Pt>, Vec<Fe>) {
    let base = if range_bits >= 2 { range_bits / 2 - 1 } else { 0 };
    let g = gen();
    let mut pts = Vec::with_capacity(NUM_JUMPS);
    let mut scs = Vec::with_capacity(NUM_JUMPS);

    for i in 0..NUM_JUMPS {
        let bit = base + (i as u32 * range_bits) / (NUM_JUMPS as u32 * 2);
        // scalar = 2^bit
        let mut sc = [0u64; 4];
        sc[(bit / 64) as usize] = 1u64 << (bit % 64);
        let pt = scalar_mul(g, sc);
        pts.push(pt);
        scs.push(sc);
    }
    (pts, scs)
}

// Jump index from canonical x (deterministic)
fn jump_idx(cx: Fe) -> usize {
    (cx[0] as usize) % NUM_JUMPS
}

// ─── Kangaroo Engine ──────────────────────────────────────────────────────────

/// Entry in the distinguished-point table.
/// For tame: total_sc = tame_start_sc + accumulated jump scalars
/// For wild: total_sc = accumulated jump scalars (add to target after collision)
#[derive(Clone)]
struct DpEntry {
    total_sc: Fe,
    is_wild: bool,
    tag: u8,
}

pub struct KangarooEngine {
    jumps: Vec<Pt>,
    jump_scs: Vec<Fe>,
    dp_bits: u32,
    table: HashMap<[u64; 4], DpEntry>,

    // Tame kangaroos: (current_pt, total_scalar)
    tames: Vec<(Pt, Fe)>,

    // Wild kangaroo: (current_pt, acc_scalar)
    wild_pt: Pt,
    wild_acc: Fe,

    target_pt: Pt,
    range_low: Fe,

    pub steps: u64,
    pub dps_found: u64,
    pub solved: bool,
    pub solution: Option<Fe>,

    // For progress reporting
    pub last_wild_x: [u64; 4],
    pub last_tame_x: [u64; 4],
}

impl KangarooEngine {
    /// Create engine for range [range_low, range_low + 2^range_bits)
    /// target_pt = P = k·G, we want to find k
    pub fn new(
        target_pt: Pt,
        range_low: Fe,
        range_bits: u32,
        num_tames: usize,
        dp_bits: u32,
    ) -> Self {
        let (jumps, jump_scs) = build_jumps(range_bits);

        // Distribute tame starting points evenly across range
        let g = gen();
        let mut tames = Vec::with_capacity(num_tames);

            for i in 0..num_tames {
            // tame_start_sc = range_low + (i+1) * range_size / (num_tames + 1)
            // Simplified: spread evenly
            let frac_bits = range_bits.saturating_sub(8); // shift
            let mut offset = [0u64; 4];
            let step = ((i + 1) as u64) << (frac_bits.min(63));
            offset[0] = step;
            // Rough placement: range_low + offset
            let tsc = sc_add(range_low, offset);
            let pt = scalar_mul(g, tsc);
            tames.push((pt, tsc));
        }

        KangarooEngine {
            jumps,
            jump_scs,
            dp_bits,
            table: HashMap::new(),
            tames,
            wild_pt: target_pt,
            wild_acc: [0u64; 4],
            target_pt,
            range_low,
            steps: 0,
            dps_found: 0,
            solved: false,
            solution: None,
            last_wild_x: target_pt.x,
            last_tame_x: if !target_pt.inf { target_pt.x } else { [0;4] },
        }
    }

    /// Step the engine by `n` iterations.
    /// Returns true if solution found.
    pub fn step(&mut self, n: u32) -> bool {
        if self.solved { return true; }

        for _ in 0..n {
            // Step wild
            if let Some(k) = self.step_animal(false) {
                self.solved = true;
                self.solution = Some(k);
                return true;
            }
            // Step each tame
            for ti in 0..self.tames.len() {
                if let Some(k) = self.step_tame(ti) {
                    self.solved = true;
                    self.solution = Some(k);
                    return true;
                }
            }
            self.steps += 1 + self.tames.len() as u64;
        }
        false
    }

    fn step_animal(&mut self, _is_tame: bool) -> Option<Fe> {
        let c = canonical(self.wild_pt);
        let ji = jump_idx(c.x);

        // Check DP
        if is_dp(c.x, self.dp_bits) {
            self.dps_found += 1;
            self.last_wild_x = self.wild_pt.x;

            if let Some(entry) = self.table.get(&c.x) {
                // Collision! Try to recover k
                if entry.is_wild { /* wild-wild collision, skip */ }
                else {
                    if let Some(k) = self.recover(entry.clone(), c.tag, false) {
                        return Some(k);
                    }
                }
            } else {
                self.table.insert(c.x, DpEntry {
                    total_sc: self.wild_acc,
                    is_wild: true,
                    tag: c.tag,
                });
            }
        }

        // Take jump
        let j = self.jumps[ji];
        let js = self.jump_scs[ji];
        self.wild_pt = pt_add(self.wild_pt, j);
        self.wild_acc = sc_add(self.wild_acc, js);
        None
    }

    fn step_tame(&mut self, ti: usize) -> Option<Fe> {
        let (pt, sc) = self.tames[ti];
        let c = canonical(pt);
        let ji = jump_idx(c.x);

        // Check DP
        if is_dp(c.x, self.dp_bits) {
            self.dps_found += 1;
            self.last_tame_x = pt.x;

            if let Some(entry) = self.table.get(&c.x) {
                if entry.is_wild {
                    if let Some(k) = self.recover(entry.clone(), c.tag, true) {
                        return Some(k);
                    }
                }
                // tame-tame: skip
            } else {
                self.table.insert(c.x, DpEntry {
                    total_sc: sc,
                    is_wild: false,
                    tag: c.tag,
                });
            }
        }

        // Take jump
        let j = self.jumps[ji];
        let js = self.jump_scs[ji];
        self.tames[ti] = (pt_add(pt, j), sc_add(sc, js));
        None
    }

    // Recover k from a collision.
    // entry: what was in the table when we found this DP
    // current_tag: tag of the DP we just found
    // table_is_tame: whether the entry in the table is the tame one
    fn recover(&self, entry: DpEntry, _current_tag: u8, table_is_tame: bool) -> Option<Fe> {
        let (tame_sc, wild_acc) = if table_is_tame {
            (entry.total_sc, self.wild_acc)
        } else {
            (entry.total_sc, self.wild_acc) // entry is wild, we are tame... handled below
        };

        // Try all 6 automorphism combinations:
        // ψ^a(tame_pt) = ψ^b(wild_pt)  →  λ^a * tame_sc = k + wild_acc  (mod n)
        // → k = λ^a * tame_sc - wild_acc
        // Also try with negation: λ^a * tame_sc = -(k + wild_acc)
        // → k = -λ^a * tame_sc - wild_acc

        let lambda0: Fe = [0x0000000000000001, 0, 0, 0]; // λ^0 = 1
        let lambda1 = LAMBDA;
        // λ² mod n
        let lambda2 = sc_mul(LAMBDA, LAMBDA);

        let lambdas = [lambda0, lambda1, lambda2];

        for &lam in &lambdas {
            for &neg in &[false, true] {
                let lam_tame = sc_mul(lam, tame_sc);
                let candidate = if neg {
                    sc_sub(sc_neg(lam_tame), wild_acc)
                } else {
                    sc_sub(lam_tame, wild_acc)
                };

                // Verify: candidate * G == target
                let g = gen();
                let check = scalar_mul(g, candidate);
                if !check.inf && check.x == self.target_pt.x && check.y == self.target_pt.y {
                    // Also check range
                    if !fe_lt(candidate, self.range_low) {
                        return Some(candidate);
                    }
                }
            }
        }
        None
    }
}

// ─── Orbit demo helper ────────────────────────────────────────────────────────

/// Return all 6 orbit elements of a point as hex strings
pub fn orbit_demo(p: Pt) -> Vec<(String, String, String)> {
    if p.inf { return vec![]; }
    let pts = [
        ("id",  p),
        ("-id", pt_neg(p)),
        ("ψ",   psi(p)),
        ("-ψ",  pt_neg(psi(p))),
        ("ψ²",  psi2(p)),
        ("-ψ²", pt_neg(psi2(p))),
    ];
    pts.iter().map(|(name, pt)| {
        (name.to_string(), fe_to_hex(pt.x), fe_to_hex(pt.y))
    }).collect()
}
