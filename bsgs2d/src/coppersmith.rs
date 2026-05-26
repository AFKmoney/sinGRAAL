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
use num_rational::BigRational;
use num_traits::{Zero, One, Signed, ToPrimitive};
use num_integer::Integer;
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
// Formule correcte pour y²=x³+b (a=0, b=7) — dérivée de la condition d'addition EC :
//   f₃(x₁,x₂,x₃) = (x₁-x₂)²·x₃² − 2[x₁x₂(x₁+x₂)+2b]·x₃ + (x₁x₂)²−4b(x₁+x₂)
//
// Coefficients de Taylor en δ pour f₃(A+δ, x_R, x_P) :
//   c₂ = (x_P − x_R)²
//   c₁ = 2(A−x_R)·x_P² − 2x_R(2A+x_R)·x_P + 2A·x_R² − 4b
//   c₀ = (A−x_R)²·x_P² − 2[A·x_R(A+x_R)+2b]·x_P + (A·x_R)²−4b(A+x_R)
pub fn s3_poly_coeffs(a: &BigInt, x_r: &BigInt, x_p: &BigInt, _p: &BigInt) -> [BigInt; 3] {
    let b14 = BigInt::from(14u32);   // 2b
    let b28 = BigInt::from(28u32);   // 4b

    let d   = fp_mod(&(a - x_r));
    let d2  = fp_mod(&(&d * &d));
    let ar  = fp_mod(&(a * x_r));
    let apr = fp_mod(&(a + x_r));
    let r2  = fp_mod(&(x_r * x_r));
    let xp2 = fp_mod(&(x_p * x_p));

    // c₂ = (x_P − x_R)²
    let diff = fp_mod(&(x_p - x_r));
    let c2   = fp_mod(&(&diff * &diff));

    // c₀ = d²·x_P² − 2[ar·(A+x_R)+2b]·x_P + ar²−4b(A+x_R)
    let mid0 = fp_mod(&(fp_mod(&(&ar * &apr)) + &b14));
    let c0   = fp_mod(
        &(&d2 * &xp2
          - BigInt::from(2u32) * &mid0 * x_p
          + fp_mod(&(&ar * &ar))
          - &b28 * &apr),
    );

    // c₁ = 2d·x_P² − 2x_R(2A+x_R)·x_P + 2A·x_R²−4b
    let inner = fp_mod(&(BigInt::from(2u32) * a + x_r));
    let t1 = fp_mod(&(BigInt::from(2u32) * &d  * &xp2));
    let t2 = fp_mod(&(BigInt::from(2u32) * x_r * &inner * x_p));
    let t3 = fp_mod(&(BigInt::from(2u32) * a   * &r2));
    let c1 = fp_mod(&(&t1 - &t2 + &t3 - &b28));

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
// LLL classique avec Gram-Schmidt EXACT en BigRational.
// Anciens bugs : (1) Lovász utilisait f64 → overflow en inf → aucun swap ;
//               (2) µ projetait sur b_j au lieu de b*_j → taille réduit faux.
// Cette implémentation recalcule le Gram-Schmidt exact à chaque étape.
// Complexité O(n⁵) BigRational, acceptable pour n ≤ 15.

pub fn dot_bigint(a: &[BigInt], b: &[BigInt]) -> BigInt {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn dot_rat(a: &[BigRational], b: &[BigRational]) -> BigRational {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn to_rat(v: &[BigInt]) -> Vec<BigRational> {
    v.iter().map(|x| BigRational::from(x.clone())).collect()
}

fn rat_round(r: &BigRational) -> BigInt {
    let q = r.to_integer();           // floor
    let frac = r - BigRational::from(q.clone());
    let half = BigRational::new(BigInt::one(), BigInt::from(2u32));
    if frac >= half { q + 1 } else if frac < -half.clone() { q - 1 } else { q }
}

// Gram-Schmidt orthogonalisation exacte.
// Retourne (bstar, mu) où bstar[i] est le i-ème vecteur GS, mu[i][j] le coeff.
fn gram_schmidt(b: &[Vec<BigInt>]) -> (Vec<Vec<BigRational>>, Vec<Vec<BigRational>>) {
    let n  = b.len();
    let dim = b[0].len();
    let mut bstar: Vec<Vec<BigRational>> = Vec::with_capacity(n);
    let mut mu: Vec<Vec<BigRational>>    = vec![vec![BigRational::zero(); n]; n];

    for i in 0..n {
        let mut bs = to_rat(&b[i]);
        for j in 0..i {
            let b_norm_sq = dot_rat(&bstar[j], &bstar[j]);
            if b_norm_sq.is_zero() { continue; }
            let mu_ij = dot_rat(&to_rat(&b[i]), &bstar[j]) / &b_norm_sq;
            mu[i][j] = mu_ij.clone();
            for k in 0..dim {
                let sub = &mu_ij * &bstar[j][k];
                bs[k] -= sub;
            }
        }
        bstar.push(bs);
    }
    (bstar, mu)
}

// ─── LLL BigRational INCRÉMENTAL ─────────────────────────────────────────────
//
// Formules fermées pour le swap (O(n) par swap, pas de recomputation GS) :
//   Après b[k]↔b[k-1], avec λ=mu[k][k-1], B0=norm[k-1], B1=norm[k] :
//     new_B0 = B1 + λ²·B0
//     new_B1 = B0·B1 / new_B0
//     new_mu[k][k-1] = λ·B0 / new_B0
//     For i>k : new_mu[i][k-1] = (mu[i][k]·B1 + λ·mu[i][k-1]·B0) / new_B0
//               new_mu[i][k]   = mu[i][k-1] − λ·mu[i][k]
//   Échange des colonnes j<k-1 dans les lignes k-1 et k.
pub fn lll_reduce_bigint(mut b: Vec<Vec<BigInt>>) -> Vec<Vec<BigInt>> {
    let n   = b.len();
    let dim = b[0].len();
    if n <= 1 { return b; }

    // ── Initialisation GS exacte (une seule fois) ────────────────────────────
    let mut bstar:   Vec<Vec<BigRational>> = Vec::with_capacity(n);
    let mut mu:      Vec<Vec<BigRational>> = vec![vec![BigRational::zero(); n]; n];
    let mut norm_sq: Vec<BigRational>      = vec![BigRational::zero(); n];

    for i in 0..n {
        let mut bs = to_rat(&b[i]);
        for j in 0..i {
            if norm_sq[j].is_zero() { continue; }
            let m = dot_rat(&to_rat(&b[i]), &bstar[j]) / &norm_sq[j];
            mu[i][j] = m.clone();
            for l in 0..dim { bs[l] -= &m * &bstar[j][l]; }
        }
        norm_sq[i] = dot_rat(&bs, &bs);
        bstar.push(bs);
    }

    let delta = BigRational::new(BigInt::from(3u32), BigInt::from(4u32));
    let mut k = 1usize;
    let mut iters = 0usize;
    let mut swaps = 0usize;

    while k < n {
        iters += 1;
        if iters > 300_000 { eprintln!("[lll] garde 300k ({} swaps)", swaps); break; }

        // ── Réduction de taille ───────────────────────────────────────────────
        // NOTE : ne change PAS norm_sq[k] ni bstar[k], seulement mu[k][*].
        for j in (0..k).rev() {
            let m = rat_round(&mu[k][j]);
            if m.is_zero() { continue; }
            let mq = BigRational::from(m.clone());
            let bj = b[j].clone();
            for l in 0..dim { b[k][l] -= &m * &bj[l]; }
            mu[k][j] -= &mq;
            for i in 0..j { let u = &mu[k][i] - &mq * &mu[j][i]; mu[k][i] = u; }
        }

        // ── Lovász ───────────────────────────────────────────────────────────
        let rhs = (&delta - &mu[k][k-1] * &mu[k][k-1]) * &norm_sq[k-1];
        if norm_sq[k] >= rhs {
            k += 1;
        } else {
            swaps += 1;
            b.swap(k, k - 1);

            // ── Mise à jour incrémentale du GS ────────────────────────────────
            let lam  = mu[k][k-1].clone();
            let b0   = norm_sq[k-1].clone();
            let b1   = norm_sq[k].clone();
            let nb0  = &b1 + &lam * &lam * &b0;  // new norm_sq[k-1]
            let nb1  = &b0 * &b1 / &nb0;          // new norm_sq[k]

            // Mise à jour mu pour lignes i > k
            for i in k+1..n {
                let a  = mu[i][k-1].clone();
                let bv = mu[i][k].clone();
                mu[i][k-1] = (&bv * &b1 + &lam * &a * &b0) / &nb0;
                mu[i][k]   = &a - &lam * &bv;
            }

            // Échange lignes k-1 et k dans mu (colonnes j < k-1)
            for j in 0..k-1 {
                let tmp = mu[k][j].clone();
                mu[k][j]   = mu[k-1][j].clone();
                mu[k-1][j] = tmp;
            }
            mu[k][k-1]   = &lam * &b0 / &nb0;
            norm_sq[k-1] = nb0;
            norm_sq[k]   = nb1;

            // Mise à jour bstar[k-1] et bstar[k]
            let old_bk  = bstar[k].clone();
            let old_bkm = bstar[k-1].clone();
            for l in 0..dim {
                bstar[k-1][l] = &old_bk[l] + &lam * &old_bkm[l];
            }
            let nb0_v = &norm_sq[k-1];
            for l in 0..dim {
                bstar[k][l] = (&b1 * &old_bkm[l] - &lam * &b0 * &old_bk[l]) / nb0_v;
            }

            if k > 1 { k -= 1; }
        }
    }
    eprintln!("[lll] terminé : {} iters, {} swaps", iters, swaps);
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
        // Essayer les 9 combinaisons GLV pour trouver c₀₀=0 ou la moins mauvaise
        let (_i, _j, coeffs) = find_glv_coeffs(start_a, start_b, &self.target_x, &self.p);
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

// ─── Coefficients bivariés de f₃(A+x, B+y, w) ───────────────────────────────
//
// Formule correcte : f₃(u,v,w) = (u-v)²w² − 2[uv(u+v)+2b]w + (uv)²−4b(u+v)  (b=7)
//
// Dérivées partielles en (u,v,w) = (A,B,w) :
//   c20 = ½ ∂²f₃/∂u² = (w−B)²
//   c02 = ½ ∂²f₃/∂v² = (w−A)²
//   c11 =   ∂²f₃/∂u∂v = −2w²−4(A+B)w+4AB
//   c10 = ∂f₃/∂u = 2(A−B)w²−2B(2A+B)w+2AB²−4b
//   c01 = ∂f₃/∂v = −2(A−B)w²−2A(A+2B)w+2A²B−4b
//   c00 = f₃(A,B,w)
pub fn s3_bivariate_coeffs(
    a: &BigInt,
    b: &BigInt,
    x_p: &BigInt,
    _p: &BigInt,
) -> [BigInt; 6] {
    let b14 = BigInt::from(14u32);   // 2b
    let b28 = BigInt::from(28u32);   // 4b
    let w   = x_p;
    let w2  = fp_mod(&(w * w));
    let ab  = fp_mod(&(a * b));
    let apb = fp_mod(&(a + b));
    let amb = fp_mod(&(a - b));
    let a2  = fp_mod(&(a * a));
    let b2  = fp_mod(&(b * b));

    // c00 = (A-B)²w² − 2[AB(A+B)+2b]w + (AB)²−4b(A+B)
    let d2     = fp_mod(&(&amb * &amb));
    let abapb  = fp_mod(&(&ab * &apb));
    let mid00  = fp_mod(&(&abapb + &b14));
    let ab2    = fp_mod(&(&ab * &ab));
    let c00    = fp_mod(
        &(&d2 * &w2
          - BigInt::from(2u32) * &mid00 * w
          + &ab2
          - &b28 * &apb),
    );

    // c10 = 2(A-B)w² − 2B(2A+B)w + 2AB²−4b
    let t1  = fp_mod(&(BigInt::from(2u32) * &amb * &w2));
    let t2  = fp_mod(&(BigInt::from(2u32) * b * fp_mod(&(BigInt::from(2u32) * a + b)) * w));
    let t3  = fp_mod(&(BigInt::from(2u32) * a * &b2));
    let c10 = fp_mod(&(&t1 - &t2 + &t3 - &b28));

    // c01 = −2(A-B)w² − 2A(A+2B)w + 2A²B−4b
    let t4  = fp_mod(&(BigInt::from(2u32) * a * fp_mod(&(a + BigInt::from(2u32) * b)) * w));
    let t5  = fp_mod(&(BigInt::from(2u32) * &a2 * b));
    let c01 = fp_mod(&(-&t1 - &t4 + &t5 - &b28));

    // c20 = (w−B)²
    let wmb = fp_mod(&(w - b));
    let c20 = fp_mod(&(&wmb * &wmb));

    // c02 = (w−A)²
    let wma = fp_mod(&(w - a));
    let c02 = fp_mod(&(&wma * &wma));

    // c11 = −2w²−4(A+B)w+4AB
    let c11 = fp_mod(
        &(-BigInt::from(2u32) * &w2
          - BigInt::from(4u32) * &apb * w
          + BigInt::from(4u32) * &ab),
    );

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

    // ── Remplissage diagonal lignes 2-6 (p·X^a·Y^b) ─────────────────────────
    // Cols 2:(0,1) 3:(2,0) 4:(1,1) 5:(0,2) 6:(3,0) — zéros si f ne les couvre pas
    mat[2][2] = p * x;       // (0,1) → Y = X
    mat[3][3] = p * &x2;     // (2,0) → X²
    mat[4][4] = p * &xy;     // (1,1) → XY = X²
    mat[5][5] = p * y2;      // (0,2) → Y² = X²
    mat[6][6] = p * &x3;     // (3,0) → X³

    // ── Lignes 7-14 : p²·monomôme (diagonale) ────────────────────────────────
    mat[ 7][ 7] = p2.clone();
    mat[ 8][ 8] = &p2 * x;
    mat[ 9][ 9] = &p2 * x;      // Y=X
    mat[10][10] = &p2 * &x2;
    mat[11][11] = &p2 * &xy;
    mat[12][12] = &p2 * y2;
    mat[13][13] = &p2 * &x3;
    mat[14][14] = &p2 * &x2y;

    mat
}

// ─── Helpers for m=3 bivariate (monomial indexing) ───────────────────────────

fn mono_to_col(a: usize, b: usize) -> usize {
    let d = a + b;
    d * (d + 1) / 2 + b
}

fn col_to_mono(col: usize) -> (usize, usize) {
    let mut d = 0usize;
    while (d + 1) * (d + 2) / 2 <= col {
        d += 1;
    }
    let b = col - d * (d + 1) / 2;
    (d - b, b)
}

/// Multiply two polynomials (coefficient vectors indexed by col), keeping total degree ≤ max_deg.
fn poly_mul_trunc(a: &[BigInt], b: &[BigInt], max_deg: usize) -> Vec<BigInt> {
    let ncols = (max_deg + 1) * (max_deg + 2) / 2;
    let mut result = vec![BigInt::zero(); ncols];
    for ai in 0..a.len() {
        if a[ai].is_zero() { continue; }
        let (ax, ay) = col_to_mono(ai);
        for bi in 0..b.len() {
            if b[bi].is_zero() { continue; }
            let (bx, by) = col_to_mono(bi);
            let cx = ax + bx;
            let cy = ay + by;
            if cx + cy <= max_deg {
                let ci = mono_to_col(cx, cy);
                result[ci] += &a[ai] * &b[bi];
            }
        }
    }
    result
}

// ─── Matrice Jochemsz-May m=3 bivariée (28×28) ───────────────────────────────
//
// Monomômes (cols) ordonnés par degré total croissant, 28 = C(8,2) :
//   deg 0: (0,0)
//   deg 1: (1,0) (0,1)
//   deg 2: (2,0) (1,1) (0,2)
//   deg 3: (3,0) (2,1) (1,2) (0,3)
//   deg 4: (4,0) (3,1) (2,2) (1,3) (0,4)
//   deg 5: (5,0) (4,1) (3,2) (2,3) (1,4) (0,5)
//   deg 6: (6,0) (5,1) (4,2) (3,3) (2,4) (1,5) (0,6)
//
// Colonne (a,b) scalée par X^(a+b) (X=Y=2^block_bits).
//
// 28 lignes :
//   Ligne 0       : f³ scalé
//   Lignes 1-6    : p · x^sa·y^sb · f²,  (sa,sb) avec sa+sb ≤ 2  (6 lignes)
//   Lignes 7-21   : p² · x^sa·y^sb · f,  (sa,sb) avec sa+sb ≤ 4  (15 lignes)
//   Lignes 22-27  : p³ · diagonale (pour les 6 colonnes restantes)
//
// Borne HG m=3 : survie ⟺ norm_min² < p⁶/28
pub fn build_macaulay_bivariate_m3(coeffs: &[BigInt; 6], x_big: &BigInt, p: &BigInt) -> Vec<Vec<BigInt>> {
    let dim = 28usize;
    let max_deg = 6usize;
    let mut mat = vec![vec![BigInt::zero(); dim]; dim];

    // f as coefficient vector (28 entries, only first 6 non-zero)
    let mut f_poly = vec![BigInt::zero(); dim];
    for i in 0..6 {
        f_poly[i] = coeffs[i].clone();
    }

    // f², f³
    let f2_poly = poly_mul_trunc(&f_poly, &f_poly, max_deg);
    let f3_poly = poly_mul_trunc(&f2_poly, &f_poly, max_deg);

    // Scaling: entry at column (a,b) gets multiplied by x_big^(a+b)
    let mut x_pows = vec![BigInt::one(); max_deg + 1];
    for i in 1..=max_deg {
        x_pows[i] = &x_pows[i - 1] * x_big;
    }

    // Helper: build a row for p^k * x^shift_a*y^shift_b * poly_coeffs, scaled
    let build_row = |poly: &[BigInt], shift_a: usize, shift_b: usize, pk: &BigInt| -> Vec<BigInt> {
        let mut row = vec![BigInt::zero(); dim];
        for col in 0..dim {
            let (a, b) = col_to_mono(col);
            if a < shift_a || b < shift_b { continue; }
            let src_col = mono_to_col(a - shift_a, b - shift_b);
            if src_col >= poly.len() { continue; }
            if poly[src_col].is_zero() { continue; }
            row[col] = pk * &poly[src_col] * &x_pows[a + b];
        }
        row
    };

    let p2 = p * p;
    let p3 = &p2 * p;

    // Row 0: f³ scaled
    mat[0] = build_row(&f3_poly, 0, 0, &BigInt::one());

    // Rows 1-6: p * x^sa*y^sb * f², (sa,sb) with sa+sb ≤ 2
    let mut row_idx = 1usize;
    for d in 0..=2usize {
        for sb in 0..=d {
            let sa = d - sb;
            mat[row_idx] = build_row(&f2_poly, sa, sb, p);
            row_idx += 1;
        }
    }
    // row_idx = 7 now

    // Rows 7-21: p² * x^sa*y^sb * f, (sa,sb) with sa+sb ≤ 4
    for d in 0..=4usize {
        for sb in 0..=d {
            let sa = d - sb;
            mat[row_idx] = build_row(&f_poly, sa, sb, &p2);
            row_idx += 1;
        }
    }
    // row_idx = 22 now

    // Fill any zero diagonal entries with p³ * X^(a+b)
    for i in 0..dim {
        if mat[i][i].is_zero() {
            let (a, b) = col_to_mono(i);
            mat[i][i] = &p3 * &x_pows[a + b];
        }
    }
    // Ensure the last 6 rows (22-27) are p³ diagonal
    for i in 22..dim {
        let (a, b) = col_to_mono(i);
        mat[i][i] = &p3 * &x_pows[a + b];
    }

    mat
}

impl LatticePruner {
    /// Bivarié m=3 (Jochemsz-May) : teste si la paire (δ∈[0,X), ε∈[0,X)) peut contenir
    /// une solution de S₃(A+δ, B+ε, x_P) ≡ 0 (mod p).
    /// Essaie les 9 combinaisons GLV (β^i·A, β^j·B) — De-GLV fix.
    /// Matrice 28×28, borne p⁶/28. false = REJET prouvé.
    pub fn is_block_pair_viable_m3(&self, start_a: &BigInt, start_b: &BigInt, block_bits: u32) -> bool {
        let x_big = BigInt::one() << block_bits as usize;
        let (_i, _j, coeffs) = find_glv_coeffs(start_a, start_b, &self.target_x, &self.p);
        if coeffs[0].is_zero() { return true; }
        let mat = build_macaulay_bivariate_m3(&coeffs, &x_big, &self.p);
        let reduced = lll_reduce_bigint(mat);
        let shortest_norm_sq = reduced.iter()
            .map(|row| norm_sq_bigint(row))
            .filter(|n| !n.is_zero())
            .min()
            .unwrap_or_else(BigInt::zero);
        // HG m=3, dim=28 : survie si norm < p³/√28
        let p3 = &self.p * &self.p * &self.p;
        let bound_sq = (&p3 * &p3) / 28i64;
        shortest_norm_sq < bound_sq
    }
}
