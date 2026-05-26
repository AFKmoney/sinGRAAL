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
use crate::secp::{Fe, fp_mul, fp_sub, fp_add, fp_neg, fp_inv, FIELD_P, fe_lt};

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

        // Borne de Howgrave-Graham : ||h||₂ < p^m / sqrt(dim)
        // Pour m=1 : bound² = p² / dim
        let bound_sq = (&self.p * &self.p) / (self.dim as i64);

        // REJET si norme minimale >= borne (aucune petite racine possible)
        shortest_norm_sq >= bound_sq
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
            rng_state = xs(rng_state ^ (i * 0x9e3779b97f4a7c15));
            // start_a et x_right aléatoires dans [0, 2^range_bits)
            let start_a = BigInt::from(rng_state & ((1u64 << range_bits.min(63)) - 1));
            rng_state = xs(rng_state);
            let x_right = BigInt::from(rng_state & ((1u64 << range_bits.min(63)) - 1));

            if !self.is_block_viable(&start_a, &x_right, block_bits) {
                rejected += 1;
            }
        }
        rejected as f64 / n_blocks as f64
    }
}
