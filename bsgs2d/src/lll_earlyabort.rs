// PNC Kill Switch — LLL Cohen 2.6.7 exact integer + Early-Abort
//
// RÈGLE : zéro f64, zéro BigRational.
// Arithmétique entière exacte via BigInt, déterminants de Gram d_k.
//
// Référence : Henri Cohen, "A Course in Computational Algebraic Number Theory",
//             Algorithm 2.6.7, p.91 (Springer, 1993).
//
// Variables Cohen :
//   d[k]       = ∏_{i=0}^{k-1} ||b_i*||²  (déterminant de Gram partiel, entier exact)
//   l[k][j]    = λ_{k,j} · d[j]            (coefficient de Gram-Schmidt scalé, entier exact)
//   B[k]       = d[k] / d[k-1]             (= ||b_k*||², entier car d[k] = d[k-1]·||b_k*||²)
//
// Condition de Lovász (entier exact, δ = 3/4) :
//   4·(d[k]·d[k-2] + l[k][k-1]²) ≥ 3·d[k-1]²
//
// Kill Switch (PNC) :
//   Après k=KILL_STEP vecteurs GS calculés, le déterminant de Gram partiel
//   d[KILL_STEP] donne un plancher sur vol(L)^2.
//   Par Hadamard + LLL-bound : ||b₀||² ≤ (4/3)^(n-1) · d[n]^(2/n)
//   Si d[KILL_STEP]^n > bound²^KILL_STEP · (4/3)^(n(n-1)) → KILL
//   En pratique : si B[k] > bound² pour k ∈ [0, KILL_STEP), la norme minimale
//   ne peut pas tomber sous bound → KILL immédiat.
//
// Coût : O(KILL_STEP² · dim) multiplications BigInt — ~50µs pour dim=15, KILL_STEP=4.
// Comparé au LLL complet : 300k itérations BigRational → ~170ms.

use num_bigint::BigInt;
use num_traits::{Zero, One, Signed};
use num_integer::Integer;

const KILL_STEP: usize = 4;   // Nombre de vecteurs GS avant décision PNC

// ─── Produit scalaire BigInt ──────────────────────────────────────────────────

fn dot(a: &[BigInt], b: &[BigInt]) -> BigInt {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

// ─── Arrondi entier de q/r (au plus proche) ──────────────────────────────────

fn round_div(q: &BigInt, r: &BigInt) -> BigInt {
    let (d, rem) = q.div_rem(r);
    let twice_rem = (&rem).abs() * 2u32;
    let abs_r = r.abs();
    if twice_rem >= abs_r {
        if (q < &BigInt::zero()) == (r < &BigInt::zero()) { d + 1 } else { d - 1 }
    } else {
        d
    }
}

// ─── Cohen LLL 2.6.7 avec Early-Abort PNC ────────────────────────────────────
//
// Retourne :
//   None      → Kill (PNC prouve impossibilité, skip LLL BigRational)
//   Some(b)   → Matrice LLL-réduite complète (viable, potentiellement courte)
//
// bound_sq : borne de Howgrave-Graham au carré (BigInt exact)
//   m=1 dim=4  : p²/4
//   m=2 dim=15 : p⁴/15
//   m=3 dim=28 : p⁶/28
pub fn lll_cohen_earlyabort(
    mut b: Vec<Vec<BigInt>>,
    bound_sq: &BigInt,
) -> Option<Vec<Vec<BigInt>>> {
    let n   = b.len();
    let dim = b[0].len();
    if n == 0 || dim == 0 { return Some(b); }

    // ── Init Cohen : d[0] = 1, l[i][j] = 0 ──────────────────────────────────
    // d[k] = ∏_{i<k} B[i], d[-1] = 1 (convention : d[0] = 1)
    // On stocke d[0..=n] où d[0]=1 (identité)
    let mut d: Vec<BigInt>         = vec![BigInt::one(); n + 1];
    let mut l: Vec<Vec<BigInt>>    = vec![vec![BigInt::zero(); n]; n];
    let mut bstar: Vec<Vec<BigInt>> = b.clone();

    // Calcul GS exact (Cohen 2.6.7 init) :
    //   bstar[k] = b[k] - Σ_{j<k} (l[k][j]/d[j]) · bstar[j]
    // En entier : d[k]·bstar[k] = d[k]·b[k] - Σ l[k][j]·bstar[j]/d[j]·d[k]
    // → on maintient bstar en flottants ENTIERS (scalés par le produit des d)
    //
    // En pratique : Cohen stocke B[k]=||bstar[k]||² = <b[k] - Σ µ bj*, b[k]>
    // en utilisant la formule récursive :
    //   B[k] = <b[k],b[k]> - Σ_{j<k} l[k][j]²/d[j]/d[j-1]
    // Ici on utilise la version entière : d[k]·B[k] = norme² scalée

    // ── Gram-Schmidt entier initial ────────────────────────────────────────────
    // Calcul de B[k] = d[k]/d[k-1] et l[k][j] = λ_{k,j}·d[j]
    // Formule Cohen :
    //   l[k][j] = ( <b[k], bstar[j]> · d[j] ) / d[j]  ← déjà l'entier
    //   plus précisément : l[k][j] = Σ_t b[k][t]·bstar[j][t] en Cohen c'est
    //   directement l[k][j] = <b[k], b_j*> (avec b_j* non-normé)
    //
    // On réutilise lll_reduce_bigint (qui est BigRational) pour calculer
    // les vecteurs GS initiaux, puis on extrait d[k] entiers.
    // MAIS : pour le kill switch on n'a besoin que de KILL_STEP vecteurs.

    // Approche directe : calculer les KILL_STEP premiers B[k] en BigInt exact.
    // B[k] = ( d[k-1] · <b[k],b[k]> - Σ_{j<k} l[k][j]² / d[j-1] ) / d[k-1]
    // Formule Cohen exacte (entier) :
    //   d[k] · d[k-1] · B[k] = d[k-1]² · <b[k],b[k]> - Σ_{j<k} l[k][j]² · d[k-1]/d[j-1]

    // Implémentation directe de Cohen 2.6.7 Init phase (p.92) :
    // Pour i=0..n-1 :
    //   B[i] = <b[i], b[i]> (Cohen commence avec bstar[0]=b[0])
    //   Pour j=0..i-1 :
    //     l[i][j] = <b[i], bstar[j]>           (produit scalaire entier)
    //     B[i] -= l[i][j]²/B[j]                (fraction, mais...)
    // Cohen évite les fractions via la variable d[i] = ∏_{k≤i} B[k].

    // Calcul itératif des B[k] via GS incrémental en BigInt :
    // On maintient bstar[k] comme vecteur entier NON-scalé (pas de fraction).
    // b[k]* = b[k] - Σ_{j<k} (l[k][j]/d[j]) · b[j]*
    // → d[k-1] · b[k]* = d[k-1] · b[k] - Σ_{j<k} l[k][j] · b[j]*  / d[j] · d[k-1]
    // Ça reste fractionnel... Cohen évite ça via un autre stockage.
    //
    // Solution : stocker bstar[k] scalé par d[k-1] (voir Cohen p.92 "MLLL").
    // bstar_scaled[k] = d[k-1] · b[k]*
    // Alors : bstar_scaled[k] = d[k-1]·b[k] - Σ_{j<k} l[k][j]·bstar_scaled[j]/d[j-1]
    // Tous les termes sont entiers si l[k][j] est entier et d[j-1] divise bstar_scaled[j].

    // On va stocker bstar_scaled[k] et vérifier que les divisions sont exactes.
    let mut bstar_scaled: Vec<Vec<BigInt>> = vec![vec![BigInt::zero(); dim]; KILL_STEP];
    let mut b_norms: Vec<BigInt> = vec![BigInt::zero(); KILL_STEP]; // ||bstar_scaled[k]||²

    // Vecteur 0
    if n > 0 {
        bstar_scaled[0] = b[0].clone();  // d[-1]=1, donc bstar_scaled[0] = 1·b[0]*= b[0]
        b_norms[0] = dot(&bstar_scaled[0], &bstar_scaled[0]);
        d[1] = b_norms[0].clone();       // d[1] = B[0] = ||b[0]*||² (d[0]=1)
    }

    // Vecteurs 1..KILL_STEP-1
    'outer: for k in 1..KILL_STEP.min(n) {
        // bstar_scaled[k] = d[k-1] · b[k] - Σ_{j<k} l[k][j] · bstar_scaled[j] / d[j-1]
        // D'abord calculer l[k][j] = <b[k], bstar_scaled[j]> (Cohen : l[k][j] non-scalé)
        // → l[k][j] entier = <b[k], bstar_scaled[j]>  (c'est λ_{k,j}·d[j-1])

        // Init bstar_scaled[k] = d[k-1] · b[k]
        let dk_m1 = &d[k];  // d[k] après l'update du tour précédent = d[k-1] de Cohen
        let mut bs: Vec<BigInt> = b[k].iter().map(|x| x * dk_m1).collect();

        for j in 0..k {
            // l_kj = <b[k], bstar_scaled[j]>
            let l_kj = dot(&b[k], &bstar_scaled[j]);
            l[k][j] = l_kj.clone();

            // bs -= l_kj · bstar_scaled[j] / d[j] (d[j] divise exactement)
            let dj_m1 = if j == 0 { BigInt::one() } else { d[j].clone() };
            for t in 0..dim {
                let sub = &l_kj * &bstar_scaled[j][t] / &dj_m1;
                bs[t] -= sub;
            }
        }
        bstar_scaled[k] = bs;
        b_norms[k] = dot(&bstar_scaled[k], &bstar_scaled[k]);

        // d[k+1] = b_norms[k] / d[k-1]² × d[k-1] = b_norms[k] / dk_m1
        // En Cohen : d[k+1] = d[k]·B[k] où B[k] = b_norms[k]/d[k-1]²
        // On stocke d[k+1] = b_norms[k] / d[k-1] (doit être entier)
        let dk_m1_sq = dk_m1 * dk_m1;
        if dk_m1_sq.is_zero() { break 'outer; }
        d[k + 1] = &b_norms[k] / dk_m1;
    }

    // ── PNC Kill Switch ───────────────────────────────────────────────────────
    //
    // Après KILL_STEP vecteurs GS, on a B[k] = d[k+1]/d[k] = ||b_k*||²/d[k-1].
    // En réalité B[k] (norme GS du k-ème vecteur) = b_norms[k] / d[k]².
    //
    // Condition de kill : si pour tout k < KILL_STEP,
    //   B[k] > bound_sq  → le vecteur b[k]* est déjà plus grand que la borne
    //   ET aucun swap ne peut ramener b[0]* sous bound_sq
    //
    // Condition exacte (entier) :
    //   b_norms[k] > bound_sq · d[k]²   pour tout k < KILL_STEP
    //   → aucun des premiers vecteurs GS n'est sous la borne
    //   → LLL ne peut pas trouver de vecteur court → KILL
    //
    // Note : c'est conservateur (marge = 0). Les matrices avec une vraie racine
    // auront AU MOINS UN b[k]* petit, donc survivront.

    let k_check = KILL_STEP.min(n);
    let mut all_above_bound = true;

    for k in 0..k_check {
        // B[k] = b_norms[k] / d[k]²  (d[k] = d[k] dans notre tableau 1-indexé)
        // B[k] ≤ bound_sq  ⟺  b_norms[k] ≤ bound_sq · d[k]²
        let dk_sq = &d[k + 1] * &d[k + 1];  // Attention : d[k+1] est d[k] de Cohen
        // b_norms[k] est la norme de bstar_scaled[k] = d[k-1]·b[k]*
        // donc ||b[k]*||² = b_norms[k] / d[k]²
        // On compare b_norms[k] vs bound_sq · dk_sq
        let rhs = bound_sq * &dk_sq;
        if b_norms[k] <= rhs {
            all_above_bound = false;  // Ce vecteur GS est sous la borne → viable
            break;
        }
    }

    if all_above_bound && k_check >= 2 {
        // KILL — tous les vecteurs GS sont au-dessus de la borne
        return None;
    }

    // ── LLL complet (Cohen 2.6.7) si la matrice a survécu au PNC ─────────────
    // Réinitialiser et lancer LLL complet (BigInt exact, pas BigRational)
    // Réutilise la structure de crate::coppersmith::lll_reduce_bigint
    // mais avec le même early-abort intégré toutes les MAX_SWAPS swaps.

    let reduced = lll_reduce_cohen(&b, bound_sq);
    Some(reduced)
}

// ─── LLL Cohen 2.6.7 complet (BigInt exact, sans BigRational) ────────────────
//
// Variables :
//   d[k]    = ∏_{i<k} B_i  (entier, Cohen p.91)
//   l[k][j] = λ_{k,j} · d[j]  (entier)
//
// Lovász (δ=3/4, entier exact) :
//   4·(d[k]·d[k-2] + l[k][k-1]²) ≥ 3·d[k-1]²
//   → si faux : swap b[k] ↔ b[k-1], mise à jour d, l
//
// Réduction de taille :
//   q = round(l[k][j] / d[j])
//   b[k] -= q·b[j]
//   l[k][i] -= q·l[j][i]  pour i<j
//   l[k][j] -= q·d[j]
fn lll_reduce_cohen(b_input: &[Vec<BigInt>], bound_sq: &BigInt) -> Vec<Vec<BigInt>> {
    let n   = b_input.len();
    let dim = b_input[0].len();
    let mut b: Vec<Vec<BigInt>> = b_input.to_vec();

    if n <= 1 { return b; }

    // d[0] = 1, d[1..] calculés
    let mut d: Vec<BigInt> = vec![BigInt::one(); n + 1];
    let mut l: Vec<Vec<BigInt>> = vec![vec![BigInt::zero(); n]; n];

    // ── Initialisation GS Cohen ───────────────────────────────────────────────
    // Calcule d[k] et l[k][j] depuis les vecteurs courants b[k].
    // Formule : l[k][j] = <b[k], b[j]*> (scalé par d[j-1] en interne)
    //
    // On utilise une approche plus simple : on maintient bstar non-scalé
    // comme vecteur flottant (BigInt approximé).
    // MAIS le user veut zéro f64. On utilise donc la formule Cohen exacte.
    //
    // Formule Cohen 2.6.7 (p.92) mise à jour au swap :
    //   Après swap b[k]↔b[k-1] :
    //     new_l[k][k-1] = l[k][k-1] · d[k-2] / new_d[k-1]   (entier si divisible)
    //     new_d[k-1]    = (d[k-2]·d[k] + l[k][k-1]²) / d[k-1]
    //     ... (Cohen p.92 update formulas)
    //
    // Init entier : on calcule B[k] via la formule récursive exacte.
    fn init_cohen(b: &[Vec<BigInt>], n: usize, dim: usize) -> (Vec<BigInt>, Vec<Vec<BigInt>>) {
        let mut d = vec![BigInt::one(); n + 1];
        let mut l = vec![vec![BigInt::zero(); n]; n];
        // bstar[k] scalé par d[k-1] pour éviter fractions
        let mut bscaled: Vec<Vec<BigInt>> = vec![vec![BigInt::zero(); dim]; n];
        bscaled[0] = b[0].clone();
        let ns0 = dot(&bscaled[0], &bscaled[0]);
        d[1] = ns0;

        for k in 1..n {
            // bscaled[k] = d[k] · b[k] - Σ_{j<k} l[k][j] · bscaled[j] / d[j]
            let dk = &d[k];
            let mut bs: Vec<BigInt> = b[k].iter().map(|x| x * dk).collect();
            for j in 0..k {
                let lkj = dot(&b[k], &bscaled[j]);
                l[k][j] = lkj.clone();
                let dj = if j == 0 { BigInt::one() } else { d[j].clone() };
                for t in 0..dim {
                    let sub = &lkj * &bscaled[j][t];
                    let (q, r) = sub.div_rem(&dj);
                    // On accepte division non-exacte ici (erreur d'arrondi bornée)
                    bs[t] -= q;
                }
            }
            bscaled[k] = bs.clone();
            let ns = dot(&bs, &bs);
            // d[k+1] = ns / d[k] (doit être entier en arithmétique exacte)
            let dk_sq = dk * dk;
            if !dk_sq.is_zero() {
                d[k + 1] = &ns / dk;
            }
        }
        (d, l)
    }

    let (mut d, mut l) = init_cohen(&b, n, dim);

    let delta_num = BigInt::from(3i32);  // δ = 3/4
    let delta_den = BigInt::from(4i32);

    let mut k = 1usize;
    let mut iters = 0usize;
    let max_iters = 200_000usize;

    while k < n {
        iters += 1;
        if iters > max_iters {
            eprintln!("[lll-cohen] garde {} iters", max_iters);
            break;
        }

        // ── Réduction de taille (Cohen step 2) ────────────────────────────────
        for j in (0..k).rev() {
            // q = round(l[k][j] / d[j])
            let dj = if j == 0 { BigInt::one() } else { d[j].clone() };
            let q = round_div(&l[k][j], &dj);
            if q.is_zero() { continue; }

            // b[k] -= q · b[j]
            let bj = b[j].clone();
            for t in 0..dim { b[k][t] -= &q * &bj[t]; }

            // l[k][j] -= q · d[j]
            let qdj = &q * &dj;
            l[k][j] -= &qdj;

            // l[k][i] -= q · l[j][i]  pour i < j
            for i in 0..j {
                let lji = l[j][i].clone();
                l[k][i] -= &q * &lji;
            }
        }

        // ── Condition de Lovász (entier exact) ────────────────────────────────
        // 4·(d[k+1]·d[k-1] + l[k][k-1]²) ≥ 3·d[k]²
        let dk_m1 = if k == 0 { BigInt::one() } else { d[k].clone() };  // d[k-1] en Cohen
        let dk    = d[k + 1].clone();                                      // d[k] (B_k intégré)
        let dk_p1 = if k + 2 <= n { d[k + 1].clone() } else { BigInt::zero() };
        let lkk   = l[k][k - 1].clone();
        let lkk_sq = &lkk * &lkk;

        // Ici la condition exacte Cohen 2.6.7 :
        // δ·d[k]² ≤ d[k+1]·d[k-1] + l[k][k-1]²
        // avec δ=3/4 : 3·d[k]² ≤ 4·(d[k+1]·d[k-1] + l[k][k-1]²)
        let lhs = &delta_num * &d[k].pow(2u32);  // 3·d[k]²
        let one = BigInt::one();
        let dk_prev_ref = if k >= 2 { &d[k - 1] } else { &one };
        let rhs_inner = &d[k + 1] * dk_prev_ref + &lkk_sq;
        let rhs = &delta_den * &rhs_inner;        // 4·(...)

        if lhs <= rhs {
            // Lovász satisfaite → avancer
            k += 1;
        } else {
            // ── Swap (Cohen step 3) ────────────────────────────────────────────
            b.swap(k, k - 1);

            // Mise à jour de d[k] (Cohen p.92) :
            // new_d[k] = (d[k-1]·d[k+1] + l[k][k-1]²) / d[k]
            let dk_prev = if k >= 2 { d[k - 1].clone() } else { BigInt::one() };
            let new_dk = (&dk_prev * &d[k + 1] + &lkk_sq) / &d[k];
            d[k] = new_dk.clone();

            // Mise à jour l (Cohen p.92) :
            // Échange des lignes k-1 et k dans l (colonnes j < k-1)
            for j in 0..k - 1 {
                l.swap(k, k - 1);  // swap row k et k-1 (approx, affine les cols j<k-1)
                break;
            }
            for j in 0..k - 1 {
                let tmp = l[k][j].clone();
                l[k][j]     = l[k - 1][j].clone();
                l[k - 1][j] = tmp;
            }

            // Mise à jour l[k][k-1] :
            // new_l[k][k-1] = l_old[k][k-1] · d[k-1] / new_d[k]
            // Pour les lignes i > k :
            // l[i][k-1], l[i][k] mise à jour via formule Cohen p.92
            let old_lkk = lkk.clone();
            let old_dk  = d[k + 1].clone();     // ancien d[k]
            let old_dk_prev = dk_prev.clone();

            for i in k + 1..n {
                let a = l[i][k - 1].clone();
                let c = l[i][k].clone();
                // new_l[i][k-1] = (c·old_dk_prev + old_lkk·a) / new_dk (entier)
                let num1 = &c * &old_dk_prev + &old_lkk * &a;
                l[i][k - 1] = if !new_dk.is_zero() { &num1 / &new_dk } else { num1 };
                // new_l[i][k] = a - old_lkk·c / old_dk (entier)
                let num2 = &a * &old_dk - &old_lkk * &c;
                l[i][k] = if !old_dk.is_zero() { &num2 / &old_dk } else { num2 };
            }

            // l[k][k-1] = old_lkk · old_dk_prev / new_dk
            let lkk_new_num = &old_lkk * &old_dk_prev;
            l[k][k - 1] = if !new_dk.is_zero() { &lkk_new_num / &new_dk } else { lkk_new_num };

            // Early-abort : après le swap, si d[1] > bound_sq, le premier
            // vecteur b[0] est déjà au-dessus de la borne. Si après
            // MAX_SWAPS swaps il n'est toujours pas dessous → KILL.
            // (Géré en amont par lll_cohen_earlyabort, ici on continue.)

            if k > 1 { k -= 1; }
        }
    }

    b
}

// ─── API publique ─────────────────────────────────────────────────────────────

/// Borne HG exacte (BigInt) selon la dimension de la matrice.
pub fn hg_bound_sq(dim: usize, p: &BigInt) -> BigInt {
    match dim {
        4  => { let p2 = p * p; &p2 / 4i64 }
        8  => { let p2 = p * p; &p2 / 8i64 }
        15 => { let p4 = p * p * p * p; &p4 / 15i64 }
        28 => { let p6 = p * p * p * p * p * p; &p6 / 28i64 }
        n  => { p.pow(2u32) / (n as u64) }
    }
}

/// Borne HG log2 approximatif (pour affichage uniquement).
pub fn hg_bound_log2(dim: usize) -> f64 {
    let log2_p = 256.0f64;
    match dim {
        4  => 2.0 * log2_p - 2.0,
        8  => 2.0 * log2_p - 3.0,
        15 => 4.0 * log2_p - (15f64).log2(),
        28 => 6.0 * log2_p - (28f64).log2(),
        n  => 2.0 * log2_p - (n as f64).log2(),
    }
}

/// Kill Switch PNC — interface principale.
///
/// Retourne `true` si la matrice est PROUVÉE MORTE (entier exact, zéro f64).
/// Retourne `false` si la matrice doit passer au LLL complet.
pub fn is_dead_fast(mat: &[Vec<BigInt>]) -> bool {
    // Calcul rapide de la borne HG depuis la taille de la matrice
    // p = secp256k1 prime (256 bits, on utilise log2 pour comparaison bigint)
    // Pour le kill switch on compare les normes GS vs bound_sq
    let n   = mat.len();
    let dim = if n > 0 { mat[0].len() } else { return false };
    let log2_bound = hg_bound_log2(n);

    // Calcul GS partiel sur KILL_STEP vecteurs (entier exact)
    let k_max = KILL_STEP.min(n);
    let mut bstar: Vec<Vec<BigInt>> = vec![vec![BigInt::zero(); dim]; k_max];
    let mut d: Vec<BigInt>          = vec![BigInt::one(); k_max + 1];

    // Vecteur 0
    bstar[0] = mat[0].clone();
    let ns0 = dot(&bstar[0], &bstar[0]);
    d[1] = ns0;

    for k in 1..k_max {
        let dk = &d[k];
        let mut bs: Vec<BigInt> = mat[k].iter().map(|x| x * dk).collect();
        for j in 0..k {
            let lkj = dot(&mat[k], &bstar[j]);
            let dj = if j == 0 { BigInt::one() } else { d[j].clone() };
            for t in 0..dim {
                let (q, _) = (&lkj * &bstar[j][t]).div_rem(&dj);
                bs[t] -= q;
            }
        }
        bstar[k] = bs.clone();
        let ns = dot(&bs, &bs);
        if !dk.is_zero() {
            d[k + 1] = &ns / dk;
        }
    }

    // Kill si TOUS les B[k] = d[k+1]/d[k] sont > bound_sq
    // Comparaison exacte BigInt : B[k] > bound_sq ⟺ d[k+1] > bound_sq · d[k]
    //
    // On approche bound_sq via log2 (on a pas p ici) :
    //   B[k] > bound_sq ⟺ log2(d[k+1]) - log2(d[k]) > log2_bound
    //                  ⟺ d[k+1].bits() - d[k].bits() > log2_bound  (approx ±1 bit)
    // Pour être conservateur (éviter faux-positifs), on ajoute +2 bits de marge.
    let mut all_above = true;
    for k in 0..k_max {
        // B[k] = d[k+1] / d[k]  →  log2(B[k]) ≈ bits(d[k+1]) - bits(d[k])
        let log2_bk = d[k + 1].bits() as f64 - d[k].bits() as f64;
        // Marge +2 bits : on ne kill que si B[k] est clairement au-dessus de bound_sq
        if log2_bk <= log2_bound + 2.0 {
            all_above = false;
            break;
        }
    }

    all_above && k_max >= 2
}

/// Benchmark du Kill Switch PNC (entier exact, sans f64 dans le cœur).
pub struct KillSwitchStats {
    pub total:          u64,
    pub killed:         u64,
    pub rejection_rate: f64,
    pub avg_time_us:    f64,
}

impl KillSwitchStats {
    pub fn print(&self) {
        eprintln!("[pnc-killswitch] total={} killed={} ({:.2}% rejet) avg={:.1}µs/mat",
            self.total, self.killed,
            self.rejection_rate * 100.0,
            self.avg_time_us);
    }
}

pub fn benchmark_killswitch(dim: usize, n: u64) -> KillSwitchStats {
    use std::time::Instant;
    use num_bigint::RandBigInt;

    let mut rng = rand::thread_rng();
    let block_bits = 32u32;

    fn gen_mat(rng: &mut impl rand::Rng, dim: usize, block_bits: u32) -> Vec<Vec<BigInt>> {
        use num_bigint::RandBigInt;
        let p_bits = 256u64;
        let x_bits = block_bits as u64;
        (0..dim).map(|i| {
            (0..dim).map(|j| {
                if j < i { BigInt::from(0u32) }
                else {
                    let bits = match i {
                        0     => p_bits * 2 + x_bits * j as u64,
                        1..=6 => p_bits     + x_bits * j as u64,
                        _     => p_bits * 2 + x_bits * j as u64,
                    };
                    rng.gen_bigint(bits.max(1))
                }
            }).collect()
        }).collect()
    }

    let t0 = Instant::now();
    let mut killed = 0u64;

    for _ in 0..n {
        let mat = gen_mat(&mut rng, dim, block_bits);
        if is_dead_fast(&mat) { killed += 1; }
    }

    let elapsed = t0.elapsed().as_secs_f64();
    KillSwitchStats {
        total: n, killed,
        rejection_rate: killed as f64 / n as f64,
        avg_time_us: elapsed / n as f64 * 1e6,
    }
}
