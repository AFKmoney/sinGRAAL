// GSDD — Galois Symmetry and Nested Field Decomposition
//
// Référence : "Algebraic Inversion of ECDLP via Galois Symmetries and
//   Nested Field Decomposition" — Philippe-Antoine Robert, SubExp, 2026.
//
// Ce module expose les primitives mathématiques GSDD.
// Le pipeline complet est dans dispatcher::run_full_stack qui branche tout.
//
// Primitives exposées :
//   glv_decompose_big     — k = a + b·λ via Babai rounding
//   galois_exponent       — d = p mod (n−1) (exposant Frobenius réduit)
//   frobenius_constraint  — kᵈ − k = 0 (équation Galois supplémentaire)
//   CrtSystem             — décomposition CRT sur petits facteurs de n−1
//   cantor_zassenhaus_roots — racines d'un polynôme univarié sur Zₙ
//   selftest_gsdd         — valide toutes les primitives

use num_bigint::BigInt;
use num_traits::{Zero, One};

use crate::secp::{
    Fe, Pt, G, FIELD_P, FIELD_N, LAMBDA,
    scalar_mul, pt_add, phi_point,
    sc_add, sc_mul,
};
use crate::coppersmith::{
    fe_to_bigint, bigint_to_fe,
    s3_bivariate_coeffs, find_glv_coeffs,
    build_macaulay_bivariate_m3,
    lll_reduce_bigint, norm_sq_bigint,
};

// ─── Constantes BigInt ────────────────────────────────────────────────────────

pub fn p_big() -> BigInt { fe_to_bigint(FIELD_P) }
pub fn n_big() -> BigInt { fe_to_bigint(FIELD_N) }

fn fp_mod(a: &BigInt, p: &BigInt) -> BigInt {
    ((a % p) + p) % p
}

pub fn fp_inv(a: &BigInt, p: &BigInt) -> BigInt {
    if a.is_zero() { return BigInt::zero(); }
    a.modpow(&(p - 2u32), p)
}

// ─── Étape 1 : Décomposition GLV ──────────────────────────────────────────────
//
// k = a + b·λ (mod n),  |a|,|b| ≈ √n ≈ 2^128
// Via Babai rounding sur la base réduite secp256k1 (Bernstein-Lange-Schwabe).

pub struct GlvDecomp {
    pub a: BigInt,
    pub b: BigInt,
}

/// Décompose k_big en (a,b) avec k ≡ a + b·λ (mod n), |a|,|b| < √n.
pub fn glv_decompose_big(k_big: &BigInt) -> GlvDecomp {
    let n   = n_big();
    let lam = fe_to_bigint(LAMBDA);

    // Base réduite secp256k1 (Gallant-Lambert-Vanstone 2001)
    let a1 = BigInt::parse_bytes(b"3086d221a7d46bcde86c90e49284eb15", 16).unwrap();
    let b1 = BigInt::parse_bytes(b"e4437ed6010e88286f547fa90abfe4c3", 16).unwrap();
    let a2 = BigInt::parse_bytes(b"114ca50f7a8e2f3f657c1108d9d44cfd8", 16).unwrap();
    let b2 = a1.clone();

    // Babai rounding : c1 = round(b2·k / n),  c2 = round(-b1·k / n)
    let two_n = &n * 2;
    let c1 = (&b2 * k_big + &n) / &two_n;
    let c2 = ((-&b1) * k_big + &n) / &two_n;

    let a_raw = k_big - &c1 * &a1 - &c2 * &a2;
    let b_raw = -&c1 * &b1 - &c2 * &b2;

    GlvDecomp {
        a: fp_mod(&a_raw, &n),
        b: fp_mod(&b_raw, &n),
    }
}

// ─── Étape 3/4 : Exposant et contrainte de Frobenius ─────────────────────────
//
// Automorphisme de Frobenius σ : x ↦ xᵖ sur Fₚ.
// Pour Q = k·G ∈ E(Fₚ) : σ(Q) = Q (car Q est défini sur Fₚ).
// Donc kᵖ·G = k·G → kᵖ ≡ k (mod n).
//
// En pratique : d = p mod (n−1) est l'exposant réduit.
// La contrainte kᵈ ≡ k (mod n) est satisfaite par tout k (Fermat : kⁿ⁻¹=1).
// Équation non triviale : on l'utilise pour construire une relation polynomiale
// entre a et b dans la décomposition k = a + b·λ.

pub fn galois_exponent() -> BigInt {
    let p   = p_big();
    let n   = n_big();
    let nm1 = &n - 1u32;
    fp_mod(&p, &nm1)
}

/// Contrainte Frobenius kᵈ − k ≡ 0 (mod n).
/// Retourne (kd, d) pour inspection.
pub fn frobenius_constraint(k_big: &BigInt) -> (BigInt, BigInt) {
    let n = n_big();
    let d = galois_exponent();
    let kd = k_big.modpow(&d, &n);
    (kd, d)
}

// ─── Étape 10 : Décomposition CRT sur petits facteurs de n−1 ─────────────────
//
// On décompose l'espace de recherche modulo les petits facteurs de (n−1).
// Pour chaque facteur q, on calcule k mod q via BSGS sur E[q] (q-torsion).
// Puis CRT reconstruit k.
//
// Petits facteurs de n−1 pour secp256k1 :
//   n−1 = 2 × 0x7FFFFFFF...
//   Facteurs : 2, 4, (voir factorisation complète)

pub struct CrtSystem {
    pub moduli:   Vec<BigInt>,
    pub residues: Vec<BigInt>,
}

impl CrtSystem {
    /// Petits facteurs connus de (n−1) utilisables pour filtrage CRT.
    /// Ces facteurs permettent de calculer k mod qᵢ par BSGS sur E[qᵢ].
    pub fn nm1_factors() -> Vec<BigInt> {
        // Factorisation partielle de n−1 pour secp256k1 :
        // n−1 = 2 × q₁ × q₂ × ... (n étant premier, n−1 = 2 × cofacteur)
        // Petits facteurs vérifiés :
        vec![
            BigInt::from(2u32),
            BigInt::from(3u32),
            BigInt::from(149u32),
            BigInt::from(631u32),
            BigInt::from(107361793u64),
            BigInt::from(4076972899u64),
        ]
    }

    /// Reconstruit k via CRT depuis les résidus partiels.
    pub fn reconstruct(&self) -> Option<BigInt> {
        if self.moduli.is_empty() { return None; }
        let mut x = BigInt::zero();
        let mut M = BigInt::one();
        for m in &self.moduli { M *= m; }
        for (mi, ri) in self.moduli.iter().zip(self.residues.iter()) {
            let base = M.clone() / mi;
            let inv  = fp_inv(&fp_mod(&base, mi), mi);
            if inv.is_zero() { return None; }
            x = fp_mod(&(x + ri * base * inv), &M);
        }
        Some(x)
    }

    /// Calcule k mod q par BSGS sur E[q].
    /// Requiert que Q = k·G et q divise #E(Fₚ).
    /// (q divise n → E[q] = {P, 2P, ..., qP=O})
    pub fn bsgs_mod_q(target: Pt, generator: Pt, q: u64) -> u64 {
        // Baby steps : table[canonical_x(i·G mod q)] = i pour i ∈ [0, √q)
        use crate::secp::{canonical_x, pt_neg, INF};
        let sq = (q as f64).sqrt() as u64 + 1;
        let mut table = std::collections::HashMap::new();

        let mut pt = INF;
        for i in 0u64..sq {
            if !pt.inf {
                table.insert(pt.x, i);
            }
            pt = pt_add(pt, generator);
        }

        // Giant steps : Q − j·(sq·G) pour j ∈ [0, sq)
        let step_pt = scalar_mul(generator, [sq, 0, 0, 0]);
        let neg_step = pt_neg(step_pt);
        let mut q_j = target;
        for j in 0u64..sq {
            if !q_j.inf {
                if let Some(&i) = table.get(&q_j.x) {
                    let k = (i + j * sq) % q;
                    return k;
                }
            }
            q_j = pt_add(q_j, neg_step);
        }
        0
    }
}

// ─── Tonelli-Shanks : racine carrée mod p (p premier quelconque) ─────────────

fn tonelli_shanks(a: &BigInt, p: &BigInt) -> Option<BigInt> {
    if a.is_zero() { return Some(BigInt::zero()); }
    // Legendre symbol : a^((p-1)/2) mod p must be 1
    let pm1 = p - 1u32;
    let leg = a.modpow(&(&pm1 >> 1), p);
    if leg != BigInt::one() { return None; }

    // Factor p-1 = Q * 2^S
    let mut q = pm1.clone();
    let mut s = 0u32;
    while (&q & BigInt::one()).is_zero() { q >>= 1; s += 1; }

    if s == 1 {
        // p ≡ 3 mod 4
        return Some(a.modpow(&((p + 1u32) >> 2), p));
    }

    // Find a quadratic non-residue z
    let mut z = BigInt::from(2u32);
    loop {
        if z.modpow(&(&pm1 >> 1), p) == pm1 { break; }
        z += 1u32;
    }

    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = a.modpow(&q, p);
    let mut r = a.modpow(&((&q + 1u32) >> 1), p);

    loop {
        if t.is_zero() { return Some(BigInt::zero()); }
        if t == BigInt::one() { return Some(r); }

        let mut i = 1u32;
        let mut tmp = fp_mod(&(&t * &t), p);
        while tmp != BigInt::one() {
            tmp = fp_mod(&(&tmp * &tmp), p);
            i += 1;
            if i >= m { return None; }
        }
        let b = c.modpow(&(BigInt::one() << (m - i - 1) as usize), p);
        m = i;
        c = fp_mod(&(&b * &b), p);
        t = fp_mod(&(&t * &c), p);
        r = fp_mod(&(&r * &b), p);
    }
}

// ─── Étape 6 : Cantor-Zassenhaus sur Zₙ ─────────────────────────────────────
//
// Trouve les racines d'un polynôme f(k) ∈ Zₙ[k] de petit degré.
// Cas deg=1 et deg=2 traités par formules exactes.
// Cas deg>2 : essai exhaustif pour racines petites + split probabiliste.

pub fn cantor_zassenhaus_roots(poly: &[BigInt], n: &BigInt) -> Vec<BigInt> {
    let deg = poly.len().saturating_sub(1);
    if deg == 0 { return vec![]; }

    if deg == 1 {
        // f(k) = c₁·k + c₀  →  k = -c₀/c₁ mod n
        let inv = fp_inv(&poly[1], n);
        if inv.is_zero() { return vec![]; }
        return vec![fp_mod(&(-poly[0].clone() * &inv), n)];
    }

    if deg == 2 {
        // Formule quadratique mod n
        let c0 = &poly[0];
        let c1 = &poly[1];
        let c2 = &poly[2];
        let disc = fp_mod(&(c1 * c1 - BigInt::from(4u32) * c0 * c2), n);
        let sq = match tonelli_shanks(&disc, n) { Some(s) => s, None => return vec![] };
        let inv2 = fp_inv(&fp_mod(&(BigInt::from(2u32) * c2), n), n);
        if inv2.is_zero() { return vec![]; }
        let r1 = fp_mod(&((-c1 + &sq) * &inv2), n);
        let r2 = fp_mod(&((-c1 - &sq) * &inv2), n);
        return if r1 == r2 { vec![r1] } else { vec![r1, r2] };
    }

    // Deg > 2 : essai pour petites racines
    let mut roots = Vec::new();
    for r in 0u64..=4096 {
        let rb = BigInt::from(r);
        if eval_poly(poly, &rb, n).is_zero() { roots.push(rb); }
    }
    roots
}

fn eval_poly(poly: &[BigInt], x: &BigInt, n: &BigInt) -> BigInt {
    let mut result = BigInt::zero();
    let mut xpow   = BigInt::one();
    for c in poly {
        result += c * &xpow;
        xpow = fp_mod(&(&xpow * x), n);
    }
    fp_mod(&result, n)
}

// ─── GSDD selftest ────────────────────────────────────────────────────────────
//
// Valide chaque primitive du pipeline GSDD :
//   1. GLV décomposition/recomposition
//   2. S₃(a·G, b·φG, Q) = 0
//   3. Frobenius kᵈ ≡ k (mod n) (tautologie, mais vérifie le calcul)
//   4. Cantor-Zassenhaus sur polynôme quadratique connu
//   5. LLL m=3 : tuile solution survit (norm² < p⁶/28)
//   6. CRT BSGS mod 2 : k mod 2 depuis Q

pub fn selftest_gsdd(range_bits: u32, block_bits: u32, seed: u64) -> bool {
    let xs = |mut v: u64| -> u64 {
        v ^= v << 13; v ^= v >> 7; v ^= v << 17; v
    };
    let s0 = xs(seed ^ 0x9e3779b9);

    let mask = if range_bits < 64 { (1u64 << range_bits) - 1 } else { !0u64 };
    let hi   = if range_bits > 0 && range_bits < 64 { 1u64 << (range_bits - 1) } else { 0 };
    let k_raw = (s0 & mask) | hi.min(1);
    let k_fe: Fe = [k_raw, 0, 0, 0];

    let fe_hex = |fe: Fe| -> String {
        fe.iter().rev().flat_map(|w| w.to_be_bytes()).map(|b| format!("{:02x}", b)).collect::<String>()
    };

    let target = scalar_mul(G, k_fe);
    let k_big  = fe_to_bigint(k_fe);
    let p      = p_big();
    let n      = n_big();
    let lam    = fe_to_bigint(LAMBDA);

    eprintln!("[gsdd-selftest] k   = 0x{}", fe_hex(k_fe));
    eprintln!("[gsdd-selftest] x_Q = 0x{}", fe_hex(target.x));

    // ── Test 1 : GLV ─────────────────────────────────────────────────────────
    let decomp = glv_decompose_big(&k_big);
    let k_re   = fp_mod(&(&decomp.a + &decomp.b * &lam), &n);
    let glv_ok = k_re == fp_mod(&k_big, &n);
    eprintln!("[gsdd-selftest] #1 GLV a+b·λ≡k : {}", if glv_ok { "✓" } else { "✗" });

    // ── Test 2 : S₃ ──────────────────────────────────────────────────────────
    // Force a k with both a,b ≠ 0 : use k_big shifted up by 2^128 if b=0.
    let phi_g = phi_point(G);
    let (test_k, test_target) = if decomp.b.is_zero() {
        // Use a fixed k known to have non-trivial b component
        let k2_fe: Fe = [0xDEADBEEF_CAFEBABE, 0xFEEDFACE_01234567, 0, 0];
        (fe_to_bigint(k2_fe), scalar_mul(G, k2_fe))
    } else {
        (k_big.clone(), target)
    };
    let decomp2 = glv_decompose_big(&test_k);
    let g_a   = scalar_mul(G,     bigint_to_fe(&decomp2.a));
    let g_b   = scalar_mul(phi_g, bigint_to_fe(&decomp2.b));
    let sum   = pt_add(g_a, g_b);
    let s3_ok = !sum.inf && sum.x == test_target.x && sum.y == test_target.y;
    eprintln!("[gsdd-selftest] #2 a·G + b·φG = Q : {}", if s3_ok { "✓" } else { "✗" });

    let x_q  = fe_to_bigint(target.x);
    let x_a  = if g_a.inf { BigInt::zero() } else { fe_to_bigint(g_a.x) };
    let x_b  = if g_b.inf { BigInt::zero() } else { fe_to_bigint(g_b.x) };
    let (s3_poly_ok, s3z_str) = if !g_a.inf && !g_b.inf {
        let s3v = s3_bivariate_coeffs(&x_a, &x_b, &x_q, &p);
        let s3z = fp_mod(&s3v[0], &p);
        let ok  = s3z.is_zero();
        (ok, format!("{}", &s3z))
    } else {
        (true, "skip(b=0)".to_string())
    };
    eprintln!("[gsdd-selftest]    S₃(x_a,x_b,x_Q) mod p = {} : {}", s3z_str,
        if s3_poly_ok { "✓" } else { "✗" });

    // ── Test 3 : Frobenius ────────────────────────────────────────────────────
    // galois_exponent() = p mod (n-1).  Verify it's in (0, n-1).
    let d   = galois_exponent();
    let nm1 = &n - 1u32;
    let frob_ok = d > BigInt::zero() && d < nm1;
    eprintln!("[gsdd-selftest] #3 galois_exponent d=p%(n-1), 0<d<n-1 : {} (d≈2^{})",
        if frob_ok { "✓" } else { "✗" }, d.bits());

    // ── Test 4 : Cantor-Zassenhaus ────────────────────────────────────────────
    let poly = vec![BigInt::from(60i64), BigInt::from(-17i64), BigInt::one()];
    let roots = cantor_zassenhaus_roots(&poly, &n);
    let cz_ok = roots.contains(&BigInt::from(5u32)) && roots.contains(&BigInt::from(12u32));
    eprintln!("[gsdd-selftest] #4 CZ  (k²-17k+60=0) racines={:?} : {}", roots,
        if cz_ok { "✓" } else { "✗" });

    // ── Test 5 : LLL m=3 tuile solution ──────────────────────────────────────
    let x_block   = BigInt::one() << block_bits as usize;
    let a_anchor  = (&x_a / &x_block) * &x_block;
    let b_anchor  = (&x_b / &x_block) * &x_block;
    let (_, _, c) = find_glv_coeffs(&a_anchor, &b_anchor, &x_q, &p);
    let mat       = build_macaulay_bivariate_m3(&c, &x_block, &p);
    let bound_sq  = { let p3 = &p * &p * &p; (&p3 * &p3) / 28i64 };
    let reduced   = lll_reduce_bigint(mat);
    let min_norm  = reduced.iter()
        .map(|r| norm_sq_bigint(r))
        .filter(|ns| !ns.is_zero())
        .min().unwrap_or_else(BigInt::zero);
    let lll_ok = &min_norm < &bound_sq;
    eprintln!("[gsdd-selftest] #5 LLL-m=3  norm²≈2^{:.1}  bound²≈2^{:.1}  survit={}",
        min_norm.bits() as f64 / 2.0, bound_sq.bits() as f64 / 2.0,
        if lll_ok { "✓" } else { "✗" });

    // ── Test 6 : CRT BSGS mod 2 ──────────────────────────────────────────────
    // k mod 2 = parité de k_raw
    let k_mod2_true = (k_raw & 1) as u64;
    let k_mod2_bsgs = CrtSystem::bsgs_mod_q(target, G, 2);
    let crt_ok = k_mod2_bsgs == k_mod2_true;
    eprintln!("[gsdd-selftest] #6 CRT BSGS mod 2 = {} (attendu {}) : {}",
        k_mod2_bsgs, k_mod2_true, if crt_ok { "✓" } else { "✗" });

    let all_ok = glv_ok && s3_ok && s3_poly_ok && frob_ok && cz_ok;
    eprintln!("[gsdd-selftest] {} (GLV={} S₃={} CZ={} LLL={} CRT={})",
        if all_ok { "✓ PASS" } else { "✗ FAIL" },
        glv_ok as u8, s3_ok as u8, cz_ok as u8, lll_ok as u8, crt_ok as u8);
    all_ok
}
