// Coppersmith LLL Block Pruner
//
// Objectif : pour un bloc [start_k, start_k + 2^block_bits),
// décider en O(dim³) si la clé secrète k peut y être,
// SANS calculer k·G pour chaque k du bloc.
//
// Fondement mathématique :
//   S₃(x₁, x₂, x_P) = (x₁-x₂)²x_P² − 2(x₁+x₂)(x₁x₂+7)x_P + (x₁x₂-7)² = 0
//   ⟺ les points d'abscisses x₁, x₂, x_P se somment à O sur la courbe.
//
//   k = A + δ,  δ ∈ [0, X),  X = 2^block_bits
//   On cherche δ petit tel que (A+δ)·G = P.
//
//   Polynôme en δ (mod p) — linéarisé via les coordonnées affines :
//     f(δ) = S₃(x((A+δ)·G), x(k_R·G), x_P) ≡ 0 (mod p)
//
//   Coppersmith (Howgrave-Graham 1997) :
//     Si ||LLL(Macaulay_matrix · diag(X^i))||_min < p^m / sqrt(dim),
//     il n'existe aucun δ ∈ [0,X) solution → BLOC REJETÉ.
//
// Limitation connue (documentée honnêtement) :
//   La formulation ci-dessus fixe x(k_R·G) comme paramètre connu.
//   Pour un filtre O(1) sur (δ, k_R) simultanément, il faudrait
//   Coppersmith bivarié (Jochemsz-May 2006), dont les bornes pratiques
//   nécessitent block_bits < ~64 pour être utiles sur secp256k1.
//   L'implémentation actuelle itère sur k_R et filtre sur δ.

use num_bigint::{BigInt, ToBigInt};
use num_traits::{Zero, One, Signed, ToPrimitive};
use crate::secp::{Fe, fp_mul, fp_sub, fp_add, fp_neg, fp_inv, FIELD_P, fe_lt, BETA, BETA2};

// ─── Conversion Fe ↔ BigInt ──────────────────────────────────────────────────

pub fn fe_to_bigint(a: Fe) -> BigInt {
    let mut bytes = [0u8; 32];
    for i in 0..4 {
        let b = a[3 - i].to_be_bytes();
        bytes[i*8..(i+1)*8].copy_from_slice(&b);
    }
    BigInt::from_bytes_be(num_bigint::Sign::Plus, &bytes)
}

pub fn bigint_to_fe(n: &BigInt) -> Fe {
    let p = fe_to_bigint(FIELD_P);
    let n = ((n % &p) + &p) % &p;
    let bytes = n.to_bytes_be().1;
    let mut padded = [0u8; 32];
    let start = 32usize.saturating_sub(bytes.len());
    padded[start..].copy_from_slice(&bytes[..bytes.len().min(32)]);
    let mut r = [0u64; 4];
    for i in 0..4 {
        r[3-i] = u64::from_be_bytes(padded[i*8..(i+1)*8].try_into().unwrap());
    }
    r
}

fn p_bigint() -> BigInt { fe_to_bigint(FIELD_P) }

/// Applique β^k (k ∈ {0,1,2}) à une coordonnée x (BigInt mod p).
/// β^0·x = x,  β^1·x = β·x mod p,  β^2·x = β²·x mod p.
pub fn beta_pow_bigint(x: &BigInt, k: u8) -> BigInt {
    match k {
        0 => x.clone(),
        1 => {
            let beta = fe_to_bigint(BETA);
            fp_mod(&(x * &beta))
        }
        2 => {
            let beta2 = fe_to_bigint(BETA2);
            fp_mod(&(x * &beta2))
        }
        _ => x.clone(),
    }
}

/// Cherche parmi les 9 combinaisons (β^i·x_L, β^j·x_R) celle où c₀₀=0.
/// Retourne (i, j, coeffs) pour la combinaison gagnante, ou (0,0,coeffs_0) si aucune.
pub fn find_glv_coeffs(
    x_l: &BigInt,
    x_r: &BigInt,
    x_p: &BigInt,
    p:   &BigInt,
) -> (u8, u8, [BigInt; 6]) {
    for i in 0u8..3 {
        for j in 0u8..3 {
            let xl_k = beta_pow_bigint(x_l, i);
            let xr_k = beta_pow_bigint(x_r, j);
            let c = s3_bivariate_coeffs(&xl_k, &xr_k, x_p, p);
            if c[0].is_zero() {
                return (i, j, c);
            }
        }
    }
    // Aucune combinaison exacte — retourner (0,0) avec coeffs bruts
    let c = s3_bivariate_coeffs(x_l, x_r, x_p, p);
    (0, 0, c)
}

fn fp_mod(a: &BigInt) -> BigInt {
    let p = p_bigint();
    ((a % &p) + &p) % &p
}

// ─── S₃ coefficients en δ ─────────────────────────────────────────────────────
//
// S₃(A + δ, x_R, x_P) = c₀ + c₁·δ + c₂·δ²  (mod p)
//
// avec :
//   c₂ = (x_P - x_R)²                         (mod p)
//   c₁ = 2·(A - x_R)·x_P² − 2·(A·x_R + x_R² + A + x_R)·x_P
//          + 2·x_R·(A·x_R − 7)                 (mod p)  [dérivée de S₃ en x₁=A]
//   c₀ = S₃(A, x_R, x_P)                       (mod p)
pub fn s3_poly_coeffs(a: &BigInt, x_r: &BigInt, x_p: &BigInt, p: &BigInt) -> [BigInt; 3] {
    // c₂ = (x_P - x_R)² mod p
    let diff = fp_mod(&(x_p - x_r));
    let c2 = fp_mod(&(&diff * &diff));

    // c₀ = S₃(A, x_R, x_P) mod p
    // S₃(x1,x2,x3) = (x1-x2)²·x3² - 2·(x1+x2)·(x1·x2+7)·x3 + (x1·x2-7)²
    let d01 = fp_mod(&(a - x_r));
    let s01 = fp_mod(&(a + x_r));
    let pr  = fp_mod(&(a * x_r));
    let p7  = fp_mod(&(&pr + 7));
    let m7  = fp_mod(&(&pr - 7));
    let xp2 = fp_mod(&(x_p * x_p));
    let c0  = fp_mod(&(fp_mod(&(&d01 * &d01)) * &xp2
              - 2 * fp_mod(&(&s01 * &p7)) * x_p
              + fp_mod(&(&m7 * &m7))));

    // c₁ = ∂S₃/∂x₁ at x₁=A
    //    = 2·(A−x_R)·x_P² − 2·(A·x_R+7+x_R²+x_R)·x_P + 2·x_R·(A·x_R−7)
    // Dérivée exacte :
    //  ∂/∂x₁ [(x₁-x₂)²x₃²] = 2(x₁-x₂)x₃²
    //  ∂/∂x₁ [-2(x₁+x₂)(x₁x₂+7)x₃] = -2(x₁x₂+7+x₂(x₁+x₂))x₃ = -2(2x₁x₂+x₂²+7)x₃
    //  ∂/∂x₁ [(x₁x₂-7)²] = 2x₂(x₁x₂-7)
    let t1 = fp_mod(&(2 * &d01 * &xp2));
    let t2 = fp_mod(&(2 * (2 * fp_mod(&(a * x_r)) + fp_mod(&(x_r * x_r)) + 7) * x_p));
    let t3 = fp_mod(&(2 * x_r * &m7));
    let c1 = fp_mod(&(&t1 - &t2 + &t3));

    // Réduction finale mod p
    [fp_mod(&c0), fp_mod(&c1), fp_mod(&c2)]
}

// ─── Matrice de Coppersmith (Howgrave-Graham) pour f(δ) = c₀+c₁δ+c₂δ² ──────
//
// Pour m niveaux de shift et X = 2^block_bits :
//   Lignes : {x^i · f(xX)^j · p^(m-j)} pour i=0.., j=0..m
//
// Ici : m=1, dim=4 (minimal pour degré 2)
//   v₀ = [c₀,     c₁·X,     c₂·X²,  0      ]  (f)
//   v₁ = [0,      c₀·X,     c₁·X²,  c₂·X³  ]  (x·f)
//   v₂ = [p,      0,        0,       0      ]  (p)
//   v₃ = [0,      p·X,      0,       0      ]  (p·x)
pub fn build_macaulay_matrix(
    coeffs: &[BigInt; 3],
    x_big: &BigInt,  // X = 2^block_bits
    p: &BigInt,
    dim: usize,      // 4 ou 8
) -> Vec<Vec<BigInt>> {
    let mut mat = vec![vec![BigInt::zero(); dim]; dim];
    let [c0, c1, c2] = coeffs;
    let x1 = x_big;
    let x2 = fp_mod(&(x1 * x1));
    let x3 = fp_mod(&(&x2 * x1));

    if dim >= 4 {
        // Ligne 0 : f(δ/X · X) = c0 + c1·X·(δ/X) + c2·X²·(δ/X)² → coefficients à δ^i
        mat[0][0] = c0.clone();
        mat[0][1] = fp_mod(&(c1 * x1));
        mat[0][2] = fp_mod(&(c2 * &x2));
        // Ligne 1 : δ · f
        mat[1][1] = fp_mod(&(c0 * x1));
        mat[1][2] = fp_mod(&(c1 * &x2));
        mat[1][3] = fp_mod(&(c2 * &x3));
        // Ligne 2 : p
        mat[2][0] = p.clone();
        // Ligne 3 : p·δ (→ p·X au niveau mis à l'échelle)
        mat[3][1] = fp_mod(&(p * x1));
    }

    // Extension à dim=8 : ajouter f², x·f², p·f, p·x·f
    if dim >= 8 {
        let x4 = fp_mod(&(&x3 * x1));
        let x5 = fp_mod(&(&x4 * x1));
        let x6 = fp_mod(&(&x5 * x1));
        let x7 = fp_mod(&(&x6 * x1));

        // f² = c0² + 2c0c1·X·δ + (c1²+2c0c2)·X²·δ² + 2c1c2·X³·δ³ + c2²·X⁴·δ⁴
        let f2_0 = fp_mod(&(c0 * c0));
        let f2_1 = fp_mod(&(2 * c0 * c1 * x1));
        let f2_2 = fp_mod(&((c1 * c1 + 2 * c0 * c2) * &x2));
        let f2_3 = fp_mod(&(2 * c1 * c2 * &x3));
        let f2_4 = fp_mod(&(c2 * c2 * &x4));

        mat[4][0] = f2_0.clone();
        mat[4][1] = f2_1.clone();
        mat[4][2] = f2_2.clone();
        mat[4][3] = f2_3.clone();
        mat[4][4] = f2_4.clone();
        // δ·f²
        mat[5][1] = fp_mod(&(&f2_0 * x1));
        mat[5][2] = fp_mod(&(&f2_1 * x1));
        mat[5][3] = fp_mod(&(&f2_2 * x1));
        mat[5][4] = fp_mod(&(&f2_3 * x1));
        mat[5][5] = fp_mod(&(&f2_4 * x1));
        // p·f
        mat[6][0] = fp_mod(&(p * c0));
        mat[6][1] = fp_mod(&(p * c1 * x1));
        mat[6][2] = fp_mod(&(p * c2 * &x2));
        // p·δ·f
        mat[7][1] = fp_mod(&(p * c0 * x1));
        mat[7][2] = fp_mod(&(p * c1 * &x2));
        mat[7][3] = fp_mod(&(p * c2 * &x3));
    }

    mat
}

// ─── LLL sur matrice BigInt (δ=3/4) ─────────────────────────────────────────
//
// LLL classique, condition de Lovász vérifiée sur les normes exactes.
// Arithmétique BigInt pour les vecteurs, f64 pour les tests de Lovász.

pub fn dot_bigint(a: &[BigInt], b: &[BigInt]) -> BigInt {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm_sq_f64(v: &[BigInt]) -> f64 {
    v.iter().map(|x| { let f = x.to_f64().unwrap_or(f64::MAX / 2.0); f * f }).sum()
}

pub fn lll_reduce_bigint(mut b: Vec<Vec<BigInt>>) -> Vec<Vec<BigInt>> {
    let n = b.len();
    if n <= 1 { return b; }

    let mut k = 1usize;
    while k < n {
        // ── Size reduce b[k] against b[k-1..0] ──────────────────────────────
        for j in (0..k).rev() {
            let n_jj = dot_bigint(&b[j], &b[j]);
            if n_jj.is_zero() { continue; }
            let n_kj = dot_bigint(&b[k], &b[j]);
            // µ = round(n_kj / n_jj)
            let double: BigInt = &n_kj * 2;
            let mu = if double.abs() > n_jj.abs() {
                // |µ| > 1/2 → subtract
                let q = &n_kj / &n_jj;
                let r = &n_kj - &q * &n_jj;
                if (2 * r.abs()) > n_jj.abs() {
                    if n_kj.is_positive() { q + 1 } else { q - 1 }
                } else { q }
            } else {
                BigInt::zero()
            };
            if !mu.is_zero() {
                let bj = b[j].clone();
                for l in 0..b[k].len() {
                    b[k][l] -= &mu * &bj[l];
                }
            }
        }

        // ── Lovász condition : 4·||b_k||² >= 3·||b_{k-1}||² ────────────────
        let nk  = norm_sq_f64(&b[k]);
        let nk1 = norm_sq_f64(&b[k-1]);
        if nk1 == 0.0 || 4.0 * nk >= 3.0 * nk1 {
            k += 1;
        } else {
            b.swap(k, k-1);
            if k > 1 { k -= 1; }
        }
    }
    b
}

// ─── Norme euclidienne d'un vecteur BigInt ────────────────────────────────────

pub fn norm_sq_bigint(v: &[BigInt]) -> BigInt {
    v.iter().map(|x| x * x).sum()
}

// ─── Filtre de Coppersmith principal ─────────────────────────────────────────

pub struct LatticePruner {
    pub p:        BigInt,
    pub target_x: BigInt,
    pub dim:      usize,   // 4 ou 8
}

impl LatticePruner {
    pub fn new(target_x: Fe, dim: usize) -> Self {
        LatticePruner {
            p:        fe_to_bigint(FIELD_P),
            target_x: fe_to_bigint(target_x),
            dim,
        }
    }

    /// Teste si le bloc [start_a, start_a + 2^block_bits) peut contenir k_left
    /// pour un x_right donné.
    ///
    /// Retourne false (REJET) si Howgrave-Graham prouve qu'aucune racine δ ∈ [0,X)
    /// ne satisfait S₃(start_a + δ, x_right, target_x) ≡ 0 (mod p).
    /// Retourne true (SURVIE) si LLL ne peut pas exclure l'existence d'une racine.
    pub fn is_block_viable(&self, start_a: &BigInt, x_right: &BigInt, block_bits: u32) -> bool {
        let x_big = BigInt::one() << block_bits as usize;  // X = 2^block_bits
        let coeffs = s3_poly_coeffs(start_a, x_right, &self.target_x, &self.p);

        // Cas trivial : c₀ = 0 → δ=0 est racine → survie immédiate
        if coeffs[0].is_zero() { return true; }

        let mat = build_macaulay_matrix(&coeffs, &x_big, &self.p, self.dim);
        let reduced = lll_reduce_bigint(mat);

        // Vecteur le plus court après LLL
        let shortest_norm_sq = reduced.iter()
            .map(|row| norm_sq_bigint(row))
            .filter(|n| !n.is_zero())
            .min()
            .unwrap_or_else(BigInt::zero);

        // HG : survie (true) si et seulement si norm_min < p/√dim
        // Contraposé : norm_min ≥ p/√dim → aucune racine → REJET (false)
        let bound_sq = (&self.p * &self.p) / (self.dim as i64);
        shortest_norm_sq < bound_sq
    }

    /// Compte le taux de rejet sur N blocs aléatoires (benchmark)
    pub fn benchmark_rejection_rate(
        &self,
        block_bits: u32,
        n_blocks: u64,
        range_bits: u32,
    ) -> f64 {
        let mut rejected = 0u64;
        let mut rng_state = 0x13579bdf2468ace0u64;
        let xs = |v: u64| -> u64 { let v=v^(v<<13); let v=v^(v>>7); v^(v<<17) };

        for i in 0..n_blocks {
            rng_state = xs(rng_state ^ i.wrapping_mul(0x9e3779b97f4a7c15));
            let start_a = BigInt::from(rng_state & ((1u64 << range_bits.min(63)) - 1));
            rng_state = xs(rng_state);
            let x_right = BigInt::from(rng_state & ((1u64 << range_bits.min(63)) - 1));

            if !self.is_block_viable(&start_a, &x_right, block_bits) {
                rejected += 1;
            }
        }
        rejected as f64 / n_blocks as f64
    }

    /// Bivarié m=2 (Jochemsz-May) : teste si la paire (δ∈[0,X), ε∈[0,X)) peut contenir
    /// une solution de S₃(A+δ, B+ε, x_P) ≡ 0 (mod p).
    /// Essaie les 9 combinaisons GLV (β^i·A, β^j·B) — De-GLV fix.
    /// Matrice 15×15, bound p⁴/15. false = REJET prouvé.
    pub fn is_block_pair_viable(&self, start_a: &BigInt, start_b: &BigInt, block_bits: u32) -> bool {
        let x_big = BigInt::one() << block_bits as usize;
        // Chercher le twist β^i, β^j qui annule c₀₀ sur les centres de blocs.
        // Les blocs doivent être twistés AVANT d'entrer dans Macaulay — sinon c₀₀≠0
        // et le filtre LLL tue le bloc solution même quand la clé est dedans.
        let center_a = start_a + (&x_big >> 1);
        let center_b = start_b + (&x_big >> 1);
        let (ti, tj, _) = find_glv_coeffs(&center_a, &center_b, &self.target_x, &self.p);
        let a_tw = beta_pow_bigint(start_a, ti);
        let b_tw = beta_pow_bigint(start_b, tj);
        let (_i, _j, coeffs) = find_glv_coeffs(&a_tw, &b_tw, &self.target_x, &self.p);
        if coeffs[0].is_zero() { return true; }

        let mat = build_macaulay_bivariate_m2(&coeffs, &x_big, &self.p);
        let reduced = lll_reduce_bigint(mat);

        let shortest_norm_sq = reduced.iter()
            .map(|row| norm_sq_bigint(row))
            .filter(|n| !n.is_zero())
            .min()
            .unwrap_or_else(BigInt::zero);

        // HG m=2, dim=15 : survie si norm < p²/√15
        let p2 = &self.p * &self.p;
        let bound_sq = (&p2 * &p2) / 15i64;
        shortest_norm_sq < bound_sq
    }

    /// Benchmark du taux de rejet bivarié m=2 sur N paires aléatoires
    pub fn benchmark_bivariate_rejection_rate(
        &self,
        block_bits: u32,
        n_blocks: u64,
        range_bits: u32,
    ) -> f64 {
        let mut rejected = 0u64;
        let mut rng_state = 0xfedcba9876543210u64;
        let xs = |v: u64| -> u64 { let v=v^(v<<13); let v=v^(v>>7); v^(v<<17) };

        for i in 0..n_blocks {
            rng_state = xs(rng_state ^ i.wrapping_mul(0x9e3779b97f4a7c15));
            let start_a = BigInt::from(rng_state & ((1u64 << range_bits.min(63)) - 1));
            rng_state = xs(rng_state);
            let start_b = BigInt::from(rng_state & ((1u64 << range_bits.min(63)) - 1));

            if !self.is_block_pair_viable(&start_a, &start_b, block_bits) {
                rejected += 1;
            }
        }
        rejected as f64 / n_blocks as f64
    }
}

// ─── Coefficients bivariés de S₃(A+x, B+y, w) ───────────────────────────────
//
// Vraie formule Semaev pour y²=x³+b (a=0, b=7) :
//   S₃(u,v,w) = (u-v)²w² - 2[uv(u+v)+2b]w + u²v²-4b(u+v)
//
// Dérivé de la condition P₁+P₂+P₃=O :
//   (x₃+x₁+x₂)(x₁-x₂)² = (y₁-y₂)²
//   Élimination de y₁y₂ via y_i²=x_i³+b → polynôme en x₁,x₂,x₃ seul.
//
// Coefficients de S₃(A+x, B+y, w) = Σ c_{ij} x^i y^j :
//   c00 = (A-B)²w² - 2(AB(A+B)+14)w + A²B²-28(A+B)
//   c10 = 2(A-B)w² - 2B(2A+B)w + 2AB²-28
//   c01 = -2(A-B)w² - 2A(A+2B)w + 2A²B-28
//   c20 = (w-B)²
//   c02 = (w-A)²
//   c11 = -2w² - 4(A+B)w + 4AB
pub fn s3_bivariate_coeffs(
    a: &BigInt,
    b: &BigInt,
    x_p: &BigInt,
    _p: &BigInt,
) -> [BigInt; 6] {
    let p   = p_bigint();
    let w   = x_p;
    let w2  = fp_mod(&(w * w));
    let ab  = fp_mod(&(a * b));       // A·B
    let apb = fp_mod(&(a + b));       // A+B
    let amb = fp_mod(&(a - b));       // A-B
    let ab2 = fp_mod(&(&ab * b));     // A·B²
    let a2b = fp_mod(&(&ab * a));     // A²·B
    let ab_sq = fp_mod(&(&ab * &ab)); // (A·B)²

    // c00 = (A-B)²w² - 2(AB(A+B)+14)w + A²B²-28(A+B)
    let d2  = fp_mod(&(&amb * &amb));
    let t1  = fp_mod(&(&ab * &apb + 14));  // AB(A+B)+14
    let c00 = fp_mod(&(&d2 * &w2 + &p + &p
        - 2 * fp_mod(&(&t1 * w))
        + &ab_sq + &p
        - 28 * &apb));

    // c10 = 2(A-B)w² - 2B(2A+B)w + 2AB²-28
    let b_2apb = fp_mod(&(b * fp_mod(&(2 * a + b)))); // B(2A+B)
    let c10 = fp_mod(&(2 * &amb * &w2 + &p + &p
        - 2 * fp_mod(&(&b_2apb * w))
        + 2 * &ab2 + &p - 28));

    // c01 = -2(A-B)w² - 2A(A+2B)w + 2A²B-28
    let a_ap2b = fp_mod(&(a * fp_mod(&(a + 2 * b)))); // A(A+2B)
    let c01 = fp_mod(&(&p - 2 * &amb * &w2 + &p + &p
        - 2 * fp_mod(&(&a_ap2b * w))
        + 2 * &a2b + &p - 28));

    // c20 = (w-B)² = w² - 2Bw + B²
    let c20 = fp_mod(&(&w2 + &p - 2 * b * w + b * b));

    // c02 = (w-A)² = w² - 2Aw + A²
    let c02 = fp_mod(&(&w2 + &p - 2 * a * w + a * a));

    // c11 = -2w² - 4(A+B)w + 4AB
    let c11 = fp_mod(&(&p - 2 * &w2 + &p
        - 4 * &apb * w
        + 4 * &ab));

    [fp_mod(&c00), fp_mod(&c10), fp_mod(&c01), fp_mod(&c20), fp_mod(&c11), fp_mod(&c02)]
}

// ─── Matrice de Macaulay bivariée (Jochemsz-May m=1) ─────────────────────────
//
// Monomômes (colonnes, scalés) : {1, x·X, y·Y, x²·X², xy·XY, y²·Y²}
//
// Ligne 0 : f(x,y) scalé            → [c00, c10·X, c01·Y, c20·X², c11·XY, c02·Y²]
// Ligne 1 : p                       → [p,   0,     0,    0,      0,      0      ]
// Ligne 2 : p·x scalé               → [0,   p·X,   0,    0,      0,      0      ]
// Ligne 3 : p·y scalé               → [0,   0,     p·Y,  0,      0,      0      ]
// Ligne 4 : p·x² scalé              → [0,   0,     0,    p·X²,   0,      0      ]
// Ligne 5 : p·xy scalé              → [0,   0,     0,    0,      p·XY,   0      ]
//
// La colonne y² n'a pas de ligne p dédiée — c02=w²≠0 assure le rang plein.
pub fn build_macaulay_bivariate(
    coeffs: &[BigInt; 6],
    x_big: &BigInt,
    y_big: &BigInt,
    p: &BigInt,
) -> Vec<Vec<BigInt>> {
    let dim = 6usize;
    let mut mat = vec![vec![BigInt::zero(); dim]; dim];
    let [c00, c10, c01, c20, c11, c02] = coeffs;

    let x2 = x_big * x_big;
    let y2 = y_big * y_big;
    let xy = x_big * y_big;

    // Ligne 0 : f scalé
    mat[0][0] = c00.clone();
    mat[0][1] = c10 * x_big;
    mat[0][2] = c01 * y_big;
    mat[0][3] = c20 * &x2;
    mat[0][4] = c11 * &xy;
    mat[0][5] = c02 * &y2;

    // Lignes p-shift
    mat[1][0] = p.clone();
    mat[2][1] = p * x_big;
    mat[3][2] = p * y_big;
    mat[4][3] = p * &x2;
    mat[5][4] = p * &xy;

    mat
}

// ─── Matrice Jochemsz-May m=2 bivariée (15×15) ───────────────────────────────
//
// Monomômes (cols) ordonnés par degré total croissant :
//   0:(0,0)  1:(1,0)  2:(0,1)  3:(2,0)  4:(1,1)  5:(0,2)
//   6:(3,0)  7:(2,1)  8:(1,2)  9:(0,3)
//   10:(4,0) 11:(3,1) 12:(2,2) 13:(1,3) 14:(0,4)
//   Colonne (a,b) scalée par X^a·Y^b  (X=Y=2^block_bits).
//
// 15 lignes :
//   Ligne 0      : f²(xX,yY)
//   Lignes 1-6   : p·{1,x,y,x²,xy,y²}·f(xX,yY)
//   Lignes 7-14  : p²·{1,x,y,x²,xy,y²,x³,x²y}  (diagonale)
//
// Borne HG m=2 : survie ⟺ norm_min² < p⁴/15
pub fn build_macaulay_bivariate_m2(coeffs: &[BigInt; 6], x_big: &BigInt, p: &BigInt) -> Vec<Vec<BigInt>> {
    let [c00, c10, c01, c20, c11, c02] = coeffs;
    let dim = 15usize;
    let mut mat = vec![vec![BigInt::zero(); dim]; dim];

    let x  = x_big;
    let x2 = x * x;
    let x3 = &x2 * x;
    let x4 = &x3 * x;
    // X = Y (blocs carrés)
    let y2  = &x2;
    let y3  = &x3;
    let y4  = &x4;
    let xy   = x * x;      // X·Y = X²
    let x2y  = &x2 * x;    // X²·Y = X³
    let xy2  = x * &x2;    // X·Y² = X³
    let x2y2 = &x2 * &x2;  // X²·Y² = X⁴
    let x3y  = &x3 * x;    // X³·Y = X⁴
    let xy3  = x * &x3;    // X·Y³ = X⁴
    let p2   = p * p;

    // ── Ligne 0 : f²(xX,yY) ──────────────────────────────────────────────────
    mat[0][ 0] = c00 * c00;
    mat[0][ 1] = 2 * c00 * c10 * x;
    mat[0][ 2] = 2 * c00 * c01 * x;           // Y=X
    mat[0][ 3] = (c10 * c10 + 2 * c00 * c20) * &x2;
    mat[0][ 4] = (2 * c10 * c01 + 2 * c00 * c11) * &xy;
    mat[0][ 5] = (c01 * c01 + 2 * c00 * c02) * y2;
    mat[0][ 6] = 2 * c10 * c20 * &x3;
    mat[0][ 7] = (2 * c10 * c11 + 2 * c01 * c20) * &x2y;
    mat[0][ 8] = (2 * c10 * c02 + 2 * c01 * c11) * &xy2;
    mat[0][ 9] = 2 * c01 * c02 * y3;
    mat[0][10] = c20 * c20 * &x4;
    mat[0][11] = 2 * c20 * c11 * &x3y;
    mat[0][12] = (2 * c20 * c02 + c11 * c11) * &x2y2;
    mat[0][13] = 2 * c11 * c02 * &xy3;
    mat[0][14] = c02 * c02 * y4;

    // ── Lignes 1-6 : p·x^a·y^b·f scalé, (a,b) ∈ {(0,0)…(0,2)} ─────────────
    // Ligne 1: (0,0)
    mat[1][0] = p * c00;
    mat[1][1] = p * c10 * x;
    mat[1][2] = p * c01 * x;
    mat[1][3] = p * c20 * &x2;
    mat[1][4] = p * c11 * &xy;
    mat[1][5] = p * c02 * y2;
    // Ligne 2: (1,0) — shift x
    mat[2][1] = p * c00 * x;
    mat[2][3] = p * c10 * &x2;
    mat[2][4] = p * c01 * &xy;
    mat[2][6] = p * c20 * &x3;
    mat[2][7] = p * c11 * &x2y;
    mat[2][8] = p * c02 * &xy2;
    // Ligne 3: (0,1) — shift y=x
    mat[3][2] = p * c00 * x;
    mat[3][4] = p * c10 * &xy;
    mat[3][5] = p * c01 * y2;
    mat[3][7] = p * c20 * &x2y;
    mat[3][8] = p * c11 * &xy2;
    mat[3][9] = p * c02 * y3;
    // Ligne 4: (2,0) — shift x²
    mat[4][ 3] = p * c00 * &x2;
    mat[4][ 6] = p * c10 * &x3;
    mat[4][ 7] = p * c01 * &x2y;
    mat[4][10] = p * c20 * &x4;
    mat[4][11] = p * c11 * &x3y;
    mat[4][12] = p * c02 * &x2y2;
    // Ligne 5: (1,1) — shift xy
    mat[5][ 4] = p * c00 * &xy;
    mat[5][ 7] = p * c10 * &x2y;
    mat[5][ 8] = p * c01 * &xy2;
    mat[5][11] = p * c20 * &x3y;
    mat[5][12] = p * c11 * &x2y2;
    mat[5][13] = p * c02 * &xy3;
    // Ligne 6: (0,2) — shift y²
    mat[6][ 5] = p * c00 * y2;
    mat[6][ 8] = p * c10 * &xy2;
    mat[6][ 9] = p * c01 * y3;
    mat[6][12] = p * c20 * &x2y2;
    mat[6][13] = p * c11 * &xy3;
    mat[6][14] = p * c02 * y4;

    // ── Lignes 7-14 : p²·monomôme (diagonale) ────────────────────────────────
    mat[ 7][ 0] = p2.clone();
    mat[ 8][ 1] = &p2 * x;
    mat[ 9][ 2] = &p2 * x;      // Y=X
    mat[10][ 3] = &p2 * &x2;
    mat[11][ 4] = &p2 * &xy;
    mat[12][ 5] = &p2 * y2;
    mat[13][ 6] = &p2 * &x3;
    mat[14][ 7] = &p2 * &x2y;

    mat
}
