// PNC intra-LLL — Kill Switch Gram-Schmidt partiel
//
// Principe : si après k=4 vecteurs GS la norme partielle prouve que
// le réseau ne peut PAS contenir de vecteur < borne HG, le thread
// s'arrête immédiatement sans faire tourner le LLL complet.
//
// Borne HG m=2 (matrice 15×15) :
//   survie ⟺ ||b_0||² < p⁴/15  →  log2(bound_sq) ≈ 1021
//
// Inégalité de Hadamard (kill switch) :
//   vol(L) = ∏ ||b_i*||  ≤  ∏ ||b_i||  (bornes de Hadamard)
//   λ_1 ≥ vol(L)^(1/n) / γ_n          (Minkowski)
//   Si vol^(2/n) > bound_sq * (4/3)^(n-1) → aucun vecteur court → KILL
//
// En pratique (après k=4 GS) :
//   Gram_k = ∏_{i=0}^{k-1} ||b_i*||²   (déterminant de Gram partiel)
//   Si Gram_k^(1/k) > bound_sq * MARGIN → KILL
//
// Architecture CUDA :
//   - 1 thread = 1 matrice
//   - Matrice en f64 MSB (128 entrées max, 15×15 ou 28×28)
//   - GS partiel k=5 en registres (pas de shared mem pour les vecteurs)
//   - Si kill : atomicAdd sur counter_killed
//   - Survie : résultat en output_viables[tid]

#include <cuda_runtime.h>
#include <math.h>
#include <stdint.h>

#define MAX_DIM 28
#define GS_STEPS 5          // Nombre de vecteurs GS à calculer avant décision
#define MARGIN_BITS 8       // Marge conservatrice : bound * 2^8 avant de killer

// ─── Kernel principal ────────────────────────────────────────────────────────
//
// matrices   : [n_mats × dim × dim] en row-major f64 (MSB des BigInt)
// dim        : dimension de la matrice (15 ou 28)
// log2_bound : log2(bound_sq) de la condition HG (≈1021 pour m=2)
// results    : [n_mats] uint8, 1=viable(passer au LLL CPU), 0=KILL
// n_mats     : nombre de matrices dans le batch

__global__ void lll_killswitch_kernel(
    const double* __restrict__ matrices,
    int           dim,
    double        log2_bound_sq,
    uint8_t*      results,
    int           n_mats
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_mats) return;

    const double* mat = matrices + (long long)tid * dim * dim;

    // ── Chargement en registres (MAX_DIM×MAX_DIM, mais on lit dim×dim) ───────
    // On travaille sur GS_STEPS vecteurs max, stockés inline.
    // bstar[i][j] : composante j du i-ème vecteur GS (i < GS_STEPS)
    // On les calcule un par un pour économiser les registres.

    double bstar[GS_STEPS][MAX_DIM];   // vecteurs GS courants
    double ns[GS_STEPS];               // norm_sq des vecteurs GS

    int k_max = (dim < GS_STEPS) ? dim : GS_STEPS;

    for (int i = 0; i < k_max; i++) {
        // Charger la ligne i de la matrice
        for (int j = 0; j < dim; j++) {
            bstar[i][j] = mat[i * dim + j];
        }

        // Soustraire les projections sur les vecteurs GS précédents
        for (int prev = 0; prev < i; prev++) {
            if (ns[prev] < 1e-300) continue;

            // mu = <b[i], bstar[prev]> / ns[prev]
            double dot = 0.0;
            for (int j = 0; j < dim; j++) {
                dot += mat[i * dim + j] * bstar[prev][j];
            }
            double mu = dot / ns[prev];

            for (int j = 0; j < dim; j++) {
                bstar[i][j] -= mu * bstar[prev][j];
            }
        }

        // Calculer la norme
        double sq = 0.0;
        for (int j = 0; j < dim; j++) {
            sq += bstar[i][j] * bstar[i][j];
        }
        ns[i] = sq;
    }

    // ── Calcul du déterminant de Gram partiel en log2 ────────────────────────
    // log2_gram_k = Σ log2(ns[i]) pour i=0..k_max-1
    double log2_gram = 0.0;
    int valid_vecs = 0;

    for (int i = 0; i < k_max; i++) {
        if (ns[i] < 1e-300) continue;
        log2_gram += log2(ns[i]);
        valid_vecs++;
    }

    if (valid_vecs == 0) {
        results[tid] = 0;  // Matrice nulle → KILL
        return;
    }

    // ── Kill Switch : inégalité de Hadamard + Minkowski ───────────────────────
    //
    // λ_1 ≥ (Gram_k^(1/k) / C_n)  avec C_n ≈ (4/3)^((n-1)/4)
    // En log2 : log2(λ_1²) ≥ log2_gram/k - (n-1)/2 * log2(4/3)
    //
    // Kill si log2(λ_1²) > log2_bound_sq (bound HG)
    // i.e. log2_gram/k > log2_bound_sq + (n-1)/2 * log2(4/3) - MARGIN_BITS
    //                      ↑ si λ_1 > bound même après la facteur LLL, KILL

    double lll_correction = (dim - 1) * 0.5 * log2(4.0 / 3.0);
    double threshold = log2_bound_sq + lll_correction + MARGIN_BITS;
    double log2_gram_per_vec = log2_gram / valid_vecs;

    // Vecteur GS min (meilleur candidat pour b_0 après LLL)
    double min_ns_log2 = log2(ns[0] + 1e-300);
    for (int i = 1; i < k_max; i++) {
        if (ns[i] > 1e-300) {
            double v = log2(ns[i]);
            if (v < min_ns_log2) min_ns_log2 = v;
        }
    }

    // KILL si la moyenne et le min des normes GS sont tous les deux >> bound
    // (double critère pour éviter les faux-positifs)
    int kill = (log2_gram_per_vec > threshold) && (min_ns_log2 > log2_bound_sq + MARGIN_BITS);

    results[tid] = kill ? 0 : 1;
}

// ─── Kernel batch : Lovász early-abort pendant LLL ───────────────────────────
//
// Version plus fine : simule k_lll_steps de l'algo LLL (swaps + réduction)
// et s'arrête dès que norm_sq[0] descend sous bound_sq → VIABLE
// ou dès que le progrès est trop lent → KILL
//
// Utilisé pour les matrices qui ont survécu au kill switch GS partiel ci-dessus.

#define MAX_LLL_SWAPS 64    // Nombre max de swaps avant d'arrêter

__global__ void lll_earlyabort_kernel(
    const double* __restrict__ matrices,
    int           dim,
    double        log2_bound_sq,
    uint8_t*      results,
    int           n_mats
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_mats) return;

    if (results[tid] == 0) return;  // Déjà killed par le GS partiel

    const double* src = matrices + (long long)tid * dim * dim;

    // Charger la matrice en shared/registres (n ≤ 15 vecteurs de dim ≤ 28)
    // Pour économiser les registres, on limite à 8×28
    const int K = 8;
    double b[K][MAX_DIM];
    double bstar[K][MAX_DIM];
    double ns[K];
    double mu[K][K];   // mu[i][j] = <b[i], bstar[j]> / ns[j]

    int n = (dim < K) ? dim : K;

    // Init
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
            b[i][j] = src[i * dim + j];
        }
    }

    // GS initiale
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) bstar[i][j] = b[i][j];
        for (int prev = 0; prev < i; prev++) {
            if (ns[prev] < 1e-300) { mu[i][prev] = 0.0; continue; }
            double dot = 0.0;
            for (int j = 0; j < dim; j++) dot += b[i][j] * bstar[prev][j];
            mu[i][prev] = dot / ns[prev];
            for (int j = 0; j < dim; j++) bstar[i][j] -= mu[i][prev] * bstar[prev][j];
        }
        double sq = 0.0;
        for (int j = 0; j < dim; j++) sq += bstar[i][j] * bstar[i][j];
        ns[i] = sq;
    }

    // LLL avec early-abort
    double delta = 0.75;
    int k = 1;
    int swaps = 0;

    while (k < n && swaps < MAX_LLL_SWAPS) {
        // Réduction de taille
        for (int j = k - 1; j >= 0; j--) {
            double m = round(mu[k][j]);
            if (m == 0.0) continue;
            for (int l = 0; l < dim; l++) b[k][l] -= m * b[j][l];
            mu[k][j] -= m;
            for (int i = 0; i < j; i++) mu[k][i] -= m * mu[j][i];
        }

        // Lovász
        double rhs = (delta - mu[k][k-1] * mu[k][k-1]) * ns[k-1];
        if (ns[k] >= rhs) {
            // Check : est-ce que ns[0] est déjà sous la borne ?
            if (ns[0] < pow(2.0, log2_bound_sq)) {
                results[tid] = 1;   // VIABLE — envoyer au LLL CPU
                return;
            }
            k++;
        } else {
            swaps++;
            // Swap b[k] ↔ b[k-1]
            for (int j = 0; j < dim; j++) {
                double tmp = b[k][j]; b[k][j] = b[k-1][j]; b[k-1][j] = tmp;
            }
            // Mise à jour GS (formule de Cohen)
            double lam  = mu[k][k-1];
            double b0   = ns[k-1];
            double b1   = ns[k];
            double nb0  = b1 + lam * lam * b0;
            double nb1  = b0 * b1 / nb0;

            for (int i = k + 1; i < n; i++) {
                double a = mu[i][k-1], bv = mu[i][k];
                mu[i][k-1] = (bv * b1 + lam * a * b0) / nb0;
                mu[i][k]   = a - lam * bv;
            }
            for (int j = 0; j < k-1; j++) {
                double tmp = mu[k][j]; mu[k][j] = mu[k-1][j]; mu[k-1][j] = tmp;
            }
            mu[k][k-1] = lam * b0 / nb0;
            ns[k-1] = nb0;
            ns[k]   = nb1;

            // Mise à jour bstar[k-1] et bstar[k]
            double old_bk[MAX_DIM], old_bkm[MAX_DIM];
            for (int j = 0; j < dim; j++) { old_bk[j] = bstar[k][j]; old_bkm[j] = bstar[k-1][j]; }
            for (int j = 0; j < dim; j++) bstar[k-1][j] = old_bk[j] + lam * old_bkm[j];
            for (int j = 0; j < dim; j++) bstar[k][j] = (b1 * old_bkm[j] - lam * b0 * old_bk[j]) / nb0;

            if (k > 1) k--;
        }
    }

    // Après MAX_LLL_SWAPS : vérifier le meilleur vecteur courant
    double min_ns = ns[0];
    for (int i = 1; i < n; i++) if (ns[i] < min_ns) min_ns = ns[i];

    results[tid] = (min_ns < pow(2.0, log2_bound_sq)) ? 1 : 0;
}

// ─── Interface C exportée pour Rust FFI ──────────────────────────────────────

extern "C" {

// Filtre un batch de matrices.
// matrices_f64 : [n_mats × dim × dim] doubles (MSB des BigInt)
// dim          : 15 (m=2) ou 28 (m=3)
// log2_bound_sq: log2 de la borne HG² (≈1021.0 pour m=2, ≈1535.0 pour m=3)
// out_viable   : [n_mats] uint8 (1=envoyer au LLL CPU, 0=skip)
// Retourne le nombre de matrices viables.
int lll_killswitch_batch(
    const double* matrices_f64,
    int           n_mats,
    int           dim,
    double        log2_bound_sq,
    uint8_t*      out_viable
) {
    double*  d_mats    = nullptr;
    uint8_t* d_results = nullptr;
    size_t   mat_bytes = (size_t)n_mats * dim * dim * sizeof(double);
    size_t   res_bytes = (size_t)n_mats * sizeof(uint8_t);

    cudaMalloc(&d_mats,    mat_bytes);
    cudaMalloc(&d_results, res_bytes);

    cudaMemcpy(d_mats, matrices_f64, mat_bytes, cudaMemcpyHostToDevice);
    cudaMemset(d_results, 1, res_bytes);  // init à viable

    int threads = 256;
    int blocks  = (n_mats + threads - 1) / threads;

    // Passe 1 : GS partiel (kill rapide)
    lll_killswitch_kernel<<<blocks, threads>>>(
        d_mats, dim, log2_bound_sq, d_results, n_mats
    );

    // Passe 2 : LLL early-abort sur les survivants
    lll_earlyabort_kernel<<<blocks, threads>>>(
        d_mats, dim, log2_bound_sq, d_results, n_mats
    );

    cudaDeviceSynchronize();
    cudaMemcpy(out_viable, d_results, res_bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_mats);
    cudaFree(d_results);

    int count = 0;
    for (int i = 0; i < n_mats; i++) count += out_viable[i];
    return count;
}

} // extern "C"
