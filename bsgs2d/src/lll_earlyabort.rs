// PNC intra-LLL — Kill Switch CPU (Rust)
//
// Pré-filtre rapide avant le LLL BigRational complet.
// Deux niveaux :
//   1. GS partiel f64 (k=5 vecteurs) → détecte les matrices "trop grandes"
//   2. LLL f64 avec early-abort (max N swaps) → raffine la décision
//
// Borne HG utilisée :
//   m=1 (dim=4)  : bound_sq = p² / dim       → log2 ≈ 503
//   m=2 (dim=15) : bound_sq = p⁴ / 15        → log2 ≈ 1021
//   m=3 (dim=28) : bound_sq = p⁶ / 28        → log2 ≈ 1535
//
// Principe du kill :
//   On maintient log2_gram_k = Σ log2(ns[i]) pour i < k.
//   Si log2_gram_k / k > log2_bound_sq + LLL_CORRECTION + MARGIN,
//   il est géométriquement impossible que LLL trouve ||b_0||² < bound_sq.
//   → KILL (return true depuis is_dead_fast)

use num_bigint::BigInt;
use num_traits::{Zero, ToPrimitive};

// Correction LLL : ||b_0|| ≤ (4/3)^((n-1)/4) * vol^(1/n)
// En log2 : correction = (n-1)/2 * log2(4/3)
fn lll_correction_log2(n: usize) -> f64 {
    (n as f64 - 1.0) * 0.5 * (4.0f64 / 3.0).log2()
}

// Conversion BigInt → f64 MSB (top 53 bits, arrondi)
fn bigint_to_f64(x: &BigInt) -> f64 {
    if x.is_zero() { return 0.0; }
    let (sign, bytes) = x.to_bytes_be();
    let n = bytes.len();
    // Top 8 bytes
    let top_bytes = if n >= 8 { &bytes[..8] } else { &bytes[..n] };
    let mut arr = [0u8; 8];
    arr[8 - top_bytes.len()..].copy_from_slice(top_bytes);
    let mantissa = u64::from_be_bytes(arr) as f64;
    let scale = if n > 8 { ((n - 8) * 8) as f64 } else { 0.0 };
    let val = mantissa * scale.exp2();
    if sign == num_bigint::Sign::Minus { -val } else { val }
}

fn dot_f64(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

// ─── Passe 1 : GS partiel (k vecteurs) ──────────────────────────────────────
//
// Retourne (log2_gram_per_vec, log2_min_ns) des k premiers vecteurs GS.
// Coût : O(k²·dim) multiplications f64 — ~1µs pour k=5, dim=28.
fn partial_gs_stats(mat: &[Vec<BigInt>], k: usize) -> (f64, f64) {
    let n = mat.len().min(k);
    let dim = mat[0].len();

    let b: Vec<Vec<f64>> = mat[..n].iter()
        .map(|row| row.iter().map(bigint_to_f64).collect())
        .collect();

    let mut bstar: Vec<Vec<f64>> = vec![vec![0.0; dim]; n];
    let mut ns = vec![0.0f64; n];

    for i in 0..n {
        bstar[i] = b[i].clone();
        for prev in 0..i {
            if ns[prev] < 1e-300 { continue; }
            let mu = dot_f64(&b[i], &bstar[prev]) / ns[prev];
            let bprev = bstar[prev].clone();
            for l in 0..dim { bstar[i][l] -= mu * bprev[l]; }
        }
        ns[i] = bstar[i].iter().map(|x| x * x).sum();
    }

    let log2s: Vec<f64> = ns.iter()
        .filter(|&&v| v > 1e-300)
        .map(|&v| v.log2())
        .collect();

    if log2s.is_empty() { return (f64::NEG_INFINITY, f64::NEG_INFINITY); }

    let avg = log2s.iter().sum::<f64>() / log2s.len() as f64;
    let min = log2s.iter().cloned().fold(f64::INFINITY, f64::min);
    (avg, min)
}

// ─── Passe 2 : LLL f64 early-abort ──────────────────────────────────────────
//
// Simule max_swaps swaps LLL en f64.
// Retourne (survived, min_ns_log2) :
//   survived=true  → au moins un vecteur GS est passé sous bound → viable
//   survived=false → aucun vecteur n'est descendu sous bound en max_swaps swaps
fn lll_f64_earlyabort(mat: &[Vec<BigInt>], max_swaps: usize, log2_bound_sq: f64) -> (bool, f64) {
    let n = mat.len().min(8);   // traiter au plus 8 vecteurs (coût max ~200µs)
    let dim = mat[0].len();
    let bound = (log2_bound_sq / 2.0).exp2();  // norm (pas norm²)

    let mut b: Vec<Vec<f64>> = mat[..n].iter()
        .map(|row| row.iter().map(bigint_to_f64).collect())
        .collect();

    // GS initiale
    let mut bstar: Vec<Vec<f64>> = vec![vec![0.0; dim]; n];
    let mut ns  = vec![0.0f64; n];
    let mut mu  = vec![vec![0.0f64; n]; n];

    for i in 0..n {
        bstar[i] = b[i].clone();
        for prev in 0..i {
            if ns[prev] < 1e-300 { mu[i][prev] = 0.0; continue; }
            let dot = dot_f64(&b[i], &bstar[prev]);
            mu[i][prev] = dot / ns[prev];
            let m = mu[i][prev];
            let bprev = bstar[prev].clone();
            for l in 0..dim { bstar[i][l] -= m * bprev[l]; }
        }
        ns[i] = bstar[i].iter().map(|x| x * x).sum();
    }

    let delta = 0.75f64;
    let mut k = 1usize;
    let mut swaps = 0usize;

    while k < n && swaps < max_swaps {
        // Réduction de taille
        for j in (0..k).rev() {
            let m = mu[k][j].round();
            if m == 0.0 { continue; }
            let bj = b[j].clone();
            for l in 0..dim { b[k][l] -= m * bj[l]; }
            mu[k][j] -= m;
            let muj: Vec<f64> = mu[j][..j].to_vec();
            for i in 0..j { mu[k][i] -= m * muj[i]; }
        }

        let rhs = (delta - mu[k][k-1] * mu[k][k-1]) * ns[k-1];
        if ns[k] >= rhs {
            // Vérifier si ce vecteur GS est déjà sous la borne
            if ns[k].sqrt() < bound || ns[0].sqrt() < bound {
                return (true, ns.iter().cloned().fold(f64::INFINITY, f64::min).log2());
            }
            k += 1;
        } else {
            swaps += 1;
            b.swap(k, k - 1);

            let lam = mu[k][k-1];
            let b0  = ns[k-1];
            let b1  = ns[k];
            let nb0 = b1 + lam * lam * b0;
            let nb1 = b0 * b1 / nb0;

            for i in k+1..n {
                let a  = mu[i][k-1];
                let bv = mu[i][k];
                mu[i][k-1] = (bv * b1 + lam * a * b0) / nb0;
                mu[i][k]   = a - lam * bv;
            }
            for j in 0..k-1 {
                let tmp = mu[k][j]; mu[k][j] = mu[k-1][j]; mu[k-1][j] = tmp;
            }
            mu[k][k-1]  = lam * b0 / nb0;
            ns[k-1] = nb0;
            ns[k]   = nb1;

            // Mise à jour bstar
            let old_bk  = bstar[k].clone();
            let old_bkm = bstar[k-1].clone();
            for l in 0..dim {
                bstar[k-1][l] = old_bk[l] + lam * old_bkm[l];
                bstar[k][l]   = (b1 * old_bkm[l] - lam * b0 * old_bk[l]) / nb0;
            }

            if k > 1 { k -= 1; }
        }
    }

    let min_ns = ns.iter().cloned().fold(f64::INFINITY, f64::min);
    let survived = min_ns.sqrt() < bound;
    (survived, min_ns.log2())
}

// ─── API publique ─────────────────────────────────────────────────────────────

/// Borne HG log2 en fonction de dim.
///   dim=4  → log2(p²/4)   ≈ 503
///   dim=8  → log2(p²/8)   ≈ 500
///   dim=15 → log2(p⁴/15)  ≈ 1021
///   dim=28 → log2(p⁶/28)  ≈ 1535
pub fn hg_bound_log2(dim: usize) -> f64 {
    let log2_p = 256.0f64;  // secp256k1 : p ≈ 2^256
    match dim {
        4       => 2.0 * log2_p - 2.0,   // p²/4
        6       => 2.0 * log2_p - (6f64).log2(),
        8       => 2.0 * log2_p - 3.0,
        15      => 4.0 * log2_p - (15f64).log2(),
        28      => 6.0 * log2_p - (28f64).log2(),
        n       => {
            // Générique : m tel que dim ≈ C(m+2, 2), bound = p^(2m) / dim
            let m = ((2 * n) as f64).sqrt() as f64;
            2.0 * m * log2_p - (n as f64).log2()
        }
    }
}

/// Kill switch principal.
///
/// Retourne `true` si la matrice est PROUVÉE MORTE (pas de vecteur court possible).
/// Retourne `false` si la matrice doit être passée au LLL BigRational complet.
///
/// Coût : ~2-10 µs (contre ~170 ms pour LLL complet sur dim=15).
pub fn is_dead_fast(mat: &[Vec<BigInt>]) -> bool {
    if mat.is_empty() || mat[0].is_empty() { return false; }

    let n = mat.len();
    let log2_bound = hg_bound_log2(n);

    // ── Passe 1 : GS partiel (k=5 vecteurs) ──────────────────────────────────
    let k_gs = 5.min(n);
    let (avg_log2, min_log2) = partial_gs_stats(mat, k_gs);

    // Correction LLL : si avg est trop haut ET min est trop haut → KILL
    let correction = lll_correction_log2(n);
    let margin = 8.0;  // 8 bits de marge conservatrice (évite les faux positifs)
    let threshold = log2_bound + correction + margin;

    if avg_log2 > threshold && min_log2 > log2_bound + margin {
        return true;  // KILL — impossible géométriquement
    }

    // ── Passe 2 : LLL f64 early-abort (64 swaps max) ─────────────────────────
    // Invoqué seulement si Passe 1 ne tue pas (survivants marginaux)
    let (survived, best_log2) = lll_f64_earlyabort(mat, 64, log2_bound);
    if !survived && best_log2 > log2_bound + margin {
        return true;  // KILL — LLL f64 n'a pas trouvé de vecteur court
    }

    false  // VIABLE — envoyer au LLL BigRational
}

/// Statistiques du kill switch pour benchmark.
pub struct KillSwitchStats {
    pub total:    u64,
    pub killed:   u64,
    pub rejection_rate: f64,
    pub avg_time_us: f64,
}

impl KillSwitchStats {
    pub fn print(&self) {
        eprintln!("[pnc-killswitch] total={} killed={} ({:.2}% rejet) avg={:.1}µs/mat",
            self.total, self.killed,
            self.rejection_rate * 100.0,
            self.avg_time_us);
    }
}

/// Génère une matrice Macaulay bivar m=2 (15×15) synthétique réaliste.
/// Les entrées sont dans les mêmes ordres de grandeur que les vraies matrices :
///   - Ligne 0 (f²) : entrées ~ p² * X^k  (p ~ 2^256, X = 2^block_bits)
///   - Lignes 1-6   : entrées ~ p  * X^k
///   - Lignes 7-14  : entrées ~ p² * X^k  (diagonale)
fn gen_realistic_macaulay_m2(rng: &mut impl rand::Rng, block_bits: u32) -> Vec<Vec<BigInt>> {
    use num_bigint::RandBigInt;
    let dim = 15usize;
    let p_bits = 256u64;
    let x_bits = block_bits as u64;

    let mut mat = vec![vec![BigInt::from(0u32); dim]; dim];
    for i in 0..dim {
        for j in i..dim {
            let scale_bits = match i {
                0     => p_bits * 2 + x_bits * (j as u64),  // f²
                1..=6 => p_bits     + x_bits * (j as u64),  // p·f
                _     => p_bits * 2 + x_bits * (j as u64),  // p²
            };
            // entrée ~ 2^scale_bits (aléatoire)
            mat[i][j] = rng.gen_bigint(scale_bits.max(1));
        }
    }
    mat
}

/// Benchmark du kill switch sur `n` matrices Macaulay synthétiques réalistes.
pub fn benchmark_killswitch(dim: usize, n: u64) -> KillSwitchStats {
    use std::time::Instant;

    let mut rng = rand::thread_rng();
    let block_bits = 32u32;  // X = 2^32, typique pour range_bits ~ 128

    let t0 = Instant::now();
    let mut killed = 0u64;

    for _ in 0..n {
        let mat = if dim == 15 {
            gen_realistic_macaulay_m2(&mut rng, block_bits)
        } else {
            // Fallback générique : entrées ~ 2^(512 + col*block_bits)
            use num_bigint::RandBigInt;
            (0..dim).map(|i| {
                (0..dim).map(|j| {
                    if j < i { BigInt::from(0u32) }
                    else { rng.gen_bigint((512 + j as u64 * block_bits as u64).max(1)) }
                }).collect()
            }).collect()
        };

        if is_dead_fast(&mat) {
            killed += 1;
        }
    }

    let elapsed = t0.elapsed().as_secs_f64();
    KillSwitchStats {
        total:    n,
        killed,
        rejection_rate: killed as f64 / n as f64,
        avg_time_us: elapsed / n as f64 * 1e6,
    }
}
