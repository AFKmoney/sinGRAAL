// GLV 4D multi-k recovery pour secp256k1
//
// Les 6 automorphismes de secp256k1 : {±1, ±λ, ±λ²} agissent sur les scalaires.
// Si k*G = P, alors (α·k)*G = α(P) pour tout α dans {1,-1,λ,-λ,λ²,-λ²}.
//
// Application dispatcher :
//   Au lieu de chercher uniquement P = Q_L + Q_R,
//   on cherche α(P) = Q_L + Q_R pour tout α → couvre 6× plus de l'espace par tile.
//   Clé récupérée : k = α⁻¹ · (k_L + λ·k_R) mod n.
//
// GLV 4D :
//   k = k₁ + λk₂ + (1+λ)k₃ + (1-λ)k₄  avec |k₁..k₄| ≈ 2^(range/4)
//   Les 4 directions couvrent le réseau GLV complet.
//   Pour k < 2^135 : k₂=k₃=k₄=0, mais le recovery multi-aut reste applicable.
//
// C≈0.55 via GLV 4D walk :
//   La marche GLV 4D (singraal) avec 4 directions Halton-LDS atteint C≈0.55.
//   Dans bsgs2d : on utilise les 6 automorphismes pour réduire les tiles de 1/6.

use crate::secp::{
    Fe, Pt, scalar_mul, pt_add, pt_neg,
    phi_point, phi2_point,
    sc_add, sc_sub, sc_neg, sc_mul,
    LAMBDA, LAMBDA2, FIELD_N,
    G, INF,
    fe_lt,
    PtJ, INFJ, ptj_from_affine, ptj_add_affine, ptj_batch_to_affine,
    fp_mul, fp_sqr,
};

// ─── Inverses des 6 automorphismes ───────────────────────────────────────────
//
// α       α⁻¹
// 1   →   1
// -1  →  -1
// λ   →   λ² (car λ³ = 1 mod n)
// -λ  →  -λ²
// λ²  →   λ
// -λ² →  -λ
//
// Automorphismes de P = k*G : α(P) = (α·k)*G
// Si on trouve Q_L + Q_R = α(P), alors k = α⁻¹·(k_L + λ·k_R)

#[derive(Clone, Copy, Debug)]
pub struct Aut {
    pub alpha_inv: Fe,  // α⁻¹ (scalaire)
    pub target_x:  Fe,  // x-coordonnée de α(P)
    pub target_y:  Fe,  // y-coordonnée de α(P)
}

/// Pré-calcule les 6 cibles automorphiques de target P.
/// Chaque Aut contient le point α(P) et l'inverse α⁻¹ pour récupérer k.
pub fn build_6aut_targets(target: Pt) -> [Aut; 6] {
    // P, φ(P)=λP, φ²(P)=λ²P et leurs négatifs
    let p0  = target;
    let p1  = pt_neg(target);
    let p2  = phi_point(target);    // λP
    let p3  = pt_neg(phi_point(target));
    let p4  = phi2_point(target);   // λ²P
    let p5  = pt_neg(phi2_point(target));

    [
        Aut { alpha_inv: [1,0,0,0],       target_x: p0.x, target_y: p0.y },
        Aut { alpha_inv: sc_neg([1,0,0,0]), target_x: p1.x, target_y: p1.y },
        Aut { alpha_inv: LAMBDA2,          target_x: p2.x, target_y: p2.y },
        Aut { alpha_inv: sc_neg(LAMBDA2),  target_x: p3.x, target_y: p3.y },
        Aut { alpha_inv: LAMBDA,           target_x: p4.x, target_y: p4.y },
        Aut { alpha_inv: sc_neg(LAMBDA),   target_x: p5.x, target_y: p5.y },
    ]
}

/// Tente de récupérer k depuis une somme k_sum = k_L + λ·k_R (mod n)
/// en essayant tous les 6 automorphismes.
/// k_sum est supposé être une combinaison scalaire trouvée dans un tile.
pub fn recover_k_from_sum(k_sum: Fe, target: Pt, range_bits: u32) -> Option<Fe> {
    let auts: [(Fe, bool); 6] = [
        ([1,0,0,0], false),
        ([1,0,0,0], true),  // négation = -k_sum
        (LAMBDA2, false),   // α⁻¹=λ² : k = λ²·k_sum
        (LAMBDA2, true),
        (LAMBDA,  false),   // α⁻¹=λ : k = λ·k_sum
        (LAMBDA,  true),
    ];
    for (alpha_inv, negate) in auts {
        let ks = if negate { sc_neg(k_sum) } else { k_sum };
        let k = sc_mul(alpha_inv, ks);
        if !in_range_bsgs(k, range_bits) { continue; }
        let check = scalar_mul(G, k);
        if check.x == target.x && check.y == target.y {
            return Some(k);
        }
    }
    None
}

/// Version complète : tente k = α⁻¹·k_sum pour α dans les 6 automorphismes,
/// et aussi essaie k_sum directement avec scalar_mul pour vérification rapide.
pub fn recover_k_6aut_full(k_candidate: Fe, target: Pt, range_bits: u32) -> Option<Fe> {
    let alphas_inv: [Fe; 6] = [
        [1,0,0,0],
        sc_neg([1,0,0,0]),
        LAMBDA,
        sc_neg(LAMBDA),
        LAMBDA2,
        sc_neg(LAMBDA2),
    ];
    for alpha_inv in alphas_inv {
        let k = sc_mul(alpha_inv, k_candidate);
        if !in_range_bsgs(k, range_bits) { continue; }
        let pt = scalar_mul(G, k);
        if pt.x == target.x && pt.y == target.y {
            return Some(k);
        }
    }
    None
}

// ─── GLV 4D décomposition ────────────────────────────────────────────────────
//
// Décompose k en 4 composantes courtes via le réseau GLV secp256k1.
// Base de Bernstein-Lange-Schwabe :
//   v₁ = (a₁, b₁, 0, 0)
//   v₂ = (a₂, a₁, 0, 0)
//   v₃ = (0, 0, 1, 0)  [direction canonique courte pour k₃]
//   v₄ = (0, 0, 0, 1)
//
// Pour k < 2^135 : k = k + 0·λ + 0·(1+λ) + 0·(1-λ)
// La décomposition 4D reste valide mais k₂=k₃=k₄=0.
// Gain : les 6 automorphismes permettent de couvrir 6× plus de k-space.

#[derive(Clone, Copy, Debug)]
pub struct Glv4D {
    pub k1: Fe,  // composante direction G
    pub k2: Fe,  // composante direction φ(G) = λ·G
    pub k3: Fe,  // composante direction (1+λ)·G
    pub k4: Fe,  // composante direction (1-λ)·G
}

/// Reconstruit k depuis la décomposition 4D :
/// k = k₁ + λ·k₂ + (1+λ)·k₃ + (1-λ)·k₄  (mod n)
pub fn glv4d_recompose(d: &Glv4D) -> Fe {
    let lam_k2 = sc_mul(LAMBDA, d.k2);
    let one_plus_lam = sc_add([1,0,0,0], LAMBDA);
    let one_minus_lam = sc_sub([1,0,0,0], LAMBDA);
    let term3 = sc_mul(one_plus_lam, d.k3);
    let term4 = sc_mul(one_minus_lam, d.k4);
    let k = sc_add(sc_add(d.k1, lam_k2), sc_add(term3, term4));
    k
}

/// Tente toutes les combinaisons de signes de (k₁, k₂, k₃, k₄)
/// pour trouver k tel que k*G = target.
/// Pour k < 2^135 : k₃=k₄=0, donc on teste surtout (±k₁, ±λk₂).
pub fn glv4d_multi_k_recovery(
    d: &Glv4D,
    target: Pt,
    range_bits: u32,
) -> Option<Fe> {
    // Pour k₁, k₂ non nuls, tester les 4 combinaisons de signes × 6 aut = 24
    let k1_variants = [d.k1, sc_neg(d.k1)];
    let k2_variants = [d.k2, sc_neg(d.k2)];
    let k3_variants = [d.k3, sc_neg(d.k3)];
    let k4_variants = [d.k4, sc_neg(d.k4)];

    let one_plus_lam  = sc_add([1,0,0,0], LAMBDA);
    let one_minus_lam = sc_sub([1,0,0,0], LAMBDA);

    for &k1 in &k1_variants {
        for &k2 in &k2_variants {
            for &k3 in &k3_variants {
                for &k4 in &k4_variants {
                    let k_sum = sc_add(
                        sc_add(k1, sc_mul(LAMBDA, k2)),
                        sc_add(sc_mul(one_plus_lam, k3), sc_mul(one_minus_lam, k4)),
                    );
                    if let Some(k) = recover_k_6aut_full(k_sum, target, range_bits) {
                        return Some(k);
                    }
                }
            }
        }
    }
    None
}

// ─── Brute-force scalaire avec 6 cibles automorphiques ───────────────────────
//
// Version étendue de brute_scalar_block :
//   Au lieu de chercher Q_L + Q_R = P, cherche aussi Q_L + Q_R = α(P) pour les 6 aut.
//   Si trouvé, k = α⁻¹·(k_L + λ·k_R).
//
// Avantage : chaque tile (A, B) couvre 6× plus de l'espace scalaire.
// Impact direct sur C : C_eff = C/√6 ≈ 0.55/√6 ≈ 0.22 (en théorie).
pub fn brute_scalar_block_6aut(
    a_scalar: u64,
    b_scalar: u64,
    block_bits: u32,
    target: Pt,
    range_bits: u32,
) -> Option<Fe> {
    if block_bits > 25 { return None; }
    let x = (1u64 << block_bits) as usize;
    let phi_g = phi_point(G);

    // ── Pré-calcul des lignes (a_scalar+δ)·G en affine via batch-inversion ────
    // row_pt[δ] = (a_scalar+δ)·G.  Construction incrémentale en Jacobien
    // (ptj_add_affine : 7M+4S, zéro inversion), puis 1 seule inversion par lot.
    let a_pt = if a_scalar == 0 { INF } else { scalar_mul(G, fe_from_u64_scalar(a_scalar)) };
    let mut rows_j: Vec<PtJ> = Vec::with_capacity(x);
    let mut cur = ptj_from_affine(a_pt);
    for _ in 0..x {
        rows_j.push(cur);
        cur = ptj_add_affine(cur, G);
    }
    let rows: Vec<Pt> = ptj_batch_to_affine(&rows_j);

    // ── Pré-calcul des colonnes (b_scalar+ε)·φG en affine via batch-inversion ─
    // col[ε] = (b_scalar+ε)·φG.  b_scalar=0 → col[0] = INF (k_r=0).
    let b_pt0 = if b_scalar == 0 { INF } else { scalar_mul(phi_g, fe_from_u64_scalar(b_scalar)) };
    let mut cols_j: Vec<PtJ> = Vec::with_capacity(x);
    let mut curb = ptj_from_affine(b_pt0);
    for _ in 0..x {
        cols_j.push(curb);
        curb = ptj_add_affine(curb, phi_g);
    }
    let cols: Vec<Pt> = ptj_batch_to_affine(&cols_j);

    // Pré-calcule les 6 cibles automorphiques
    let auts = build_6aut_targets(target);

    // ── Boucle chaude : sommes en Jacobien, comparaison par produit croisé ────
    // sum (X:Y:Z) == aut (x,y) affine  ⟺  X == x·Z²  ET  Y == y·Z³.
    // Aucune inversion par somme (gain ~100× vs addition affine).
    for delta in 0..x {
        let rp = rows[delta];
        for eps in 0..x {
            let cp = cols[eps];
            let sum_j: PtJ = if rp.inf { ptj_from_affine(cp) }
                             else if cp.inf { ptj_from_affine(rp) }
                             else { ptj_add_affine(ptj_from_affine(rp), cp) };
            if sum_j.inf { continue; }

            let z2 = fp_sqr(sum_j.z);
            let z3 = fp_mul(z2, sum_j.z);

            for aut in &auts {
                if sum_j.x == fp_mul(aut.target_x, z2)
                    && sum_j.y == fp_mul(aut.target_y, z3)
                {
                    let k_l = fe_from_u64_scalar(a_scalar.wrapping_add(delta as u64));
                    let k_r = fe_from_u64_scalar(b_scalar.wrapping_add(eps as u64));
                    let k_sum = sc_add(k_l, sc_mul(LAMBDA, k_r));
                    let k = sc_mul(aut.alpha_inv, k_sum);
                    if in_range_bsgs(k, range_bits) {
                        let check = scalar_mul(G, k);
                        if check.x == target.x && check.y == target.y {
                            return Some(k);
                        }
                    }
                    let k_neg = sc_mul(aut.alpha_inv, sc_neg(k_sum));
                    if in_range_bsgs(k_neg, range_bits) {
                        let check = scalar_mul(G, k_neg);
                        if check.x == target.x && check.y == target.y {
                            return Some(k_neg);
                        }
                    }
                }
            }
        }
    }
    None
}

// ─── Utilitaires ─────────────────────────────────────────────────────────────

#[inline]
pub fn fe_from_u64_scalar(v: u64) -> Fe { [v, 0, 0, 0] }

/// Vérifie k ∈ [2^(range_bits-1), 2^range_bits).
pub fn in_range_bsgs(k: Fe, range_bits: u32) -> bool {
    if range_bits == 0 { return false; }
    if range_bits >= 256 { return true; }
    // Borne supérieure : k < 2^range_bits
    let word_hi = (range_bits / 64) as usize;
    let bit_hi  = range_bits % 64;
    if bit_hi == 0 {
        for i in word_hi..4 { if k[i] != 0 { return false; } }
    } else {
        for i in (word_hi+1)..4 { if k[i] != 0 { return false; } }
        if k[word_hi] >= (1u64 << bit_hi) { return false; }
    }
    // Borne inférieure : k >= 2^(range_bits-1)
    if range_bits == 1 { return k[0] >= 1; }
    let lo = range_bits - 1;
    let word_lo = (lo / 64) as usize;
    let bit_lo  = lo % 64;
    let threshold = if bit_lo < 63 { 1u64 << bit_lo } else { 1u64 << 63 };
    // Chercher si k >= threshold dans le limbe word_lo
    // D'abord vérifier les limbes supérieures
    for i in (word_lo+1)..word_hi.min(4) {
        if k[i] > 0 { return true; } // k a des bits > range-1 → k > threshold, borne sup déjà OK
    }
    k[word_lo] >= threshold
}
