// sinGRAAL — BSGS secp256k1 avec 6-automorphismes (1D) et GLV 2D
//
// ─── Algorithme 1D (défaut) ───────────────────────────────────────────────────
//
//   k = j·M + b,  baby table T[canonical_x(b·G)] = b pour b ∈ [1,M]
//   Giant : Q_j = P − j·M·G, chercher canonical_x(Q_j) ∈ T (N/M steps)
//   Optimal : M = √N,  total ≈ 2·√N  (N = 2^range_bits)
//
// ─── Algorithme 2D GLV (--glv) ────────────────────────────────────────────────
//
//   Décomposition GLV : k = k₁ + λ·k₂ (mod n)
//   Endomorphisme : φ(P) = (β·x, y) = λ·P  →  φ(G) = λ·G
//
//   Baby table 2D :
//     T[canonical_x((a₁·G + a₂·φG).x)] = (a₁, a₂)
//     pour (a₁, a₂) ∈ [0,M₁) × [0,M₂)   —  M₁·M₂ entrées
//
//   Giant 2D :
//     Q_{j₁,j₂} = P − j₁·M₁·G − j₂·M₂·φG
//     j₁ ∈ [0, N₁/M₁),  j₂ ∈ [0, N₂/M₂)
//
//   Récupération (6-aut) :
//     Q = α·(a₁·G + a₂·φG)  pour α ∈ {±1, ±λ, ±λ²}
//     k = j₁·M₁ + j₂·M₂·λ + α·(a₁ + a₂·λ)  (mod n)
//
// ─── Analyse des bornes N₁, N₂ (secp256k1) ──────────────────────────────────
//
//   Babai rounding : c₁ = ⌊g₁·k / 2^384⌋,  c₂ = ⌊g₂·k / 2^384⌋
//   g₁ ≈ 2^189,  g₂ ≈ 2^191
//
//   Pour k < 2^189 : c₁ = c₂ = 0  →  k₁ = k,  k₂ = 0  (décomp triviale)
//   Pour k ∈ [2^189, 2^256) : c₁, c₂ > 0  →  |k₁|, |k₂| ≤ √n ≈ 2^128
//
//   Conséquence pour puzzle #135 (k ∈ [2^134, 2^135)) :
//     k < 2^189  →  k₁ = k (135 bits),  k₂ = 0
//     → GLV 2D dégénère en BSGS 1D sur k₁
//     → Pas d'accélération sur puzzle #135 via GLV 2D
//
//   Pour résoudre puzzle #135 en "quelques GPU" :
//     Kangaroo singraal/ (C≈0.55, ~2^65 ops) reste la seule voie praticable.
//
// ─── Faisabilité BSGS ────────────────────────────────────────────────────────
//   range_bits=40 : M=2^20 → table 40 MB,  giant 2^20  → secondes  ✓
//   range_bits=50 : M=2^25 → table 1.3 GB, giant 2^25  → minutes   ✓
//   puzzle #135   : M=2^67 → 590 exaoctets              → impossible ✗

mod secp;
mod lll;
mod coppersmith;
mod lll_earlyabort;
mod dispatcher;
mod glv4d;
mod gsdd;

use clap::Parser;
use secp::*;
use std::collections::HashMap;
use std::time::Instant;
use rayon;

// ─── CLI ─────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(name = "bsgs2d", about = "BSGS 1D/2D-GLV/Semaev-MitM 6-aut secp256k1")]
struct Args {
    #[arg(long, default_value = "")]
    target_x: String,
    #[arg(long, default_value = "")]
    target_y: String,

    /// k ∈ [0, 2^range_bits).
    #[arg(long, default_value = "40")]
    range_bits: u32,

    /// baby_bits : M = 2^baby_bits par dimension (BSGS 1D/2D).
    #[arg(long)]
    baby_bits: Option<u32>,

    /// Activer le mode 2D GLV : baby table (a₁,a₂), giant (j₁,j₂).
    #[arg(long)]
    glv: bool,

    /// Activer le mode Semaev Meet-in-the-Middle (arbre de Semaev).
    #[arg(long)]
    semaev: bool,

    /// Bits par bloc pour Semaev MitM (défaut : range_bits/4, arrondi en haut).
    #[arg(long)]
    block_bits: Option<u32>,

    /// Afficher l'analyse LLL du réseau GLV avant la recherche.
    #[arg(long)]
    lll: bool,

    /// Activer le filtre Coppersmith LLL en amont du MitM Semaev.
    #[arg(long)]
    prune: bool,

    /// Benchmarker le taux de rejet du filtre Coppersmith univarié (N blocs aléatoires).
    #[arg(long)]
    prune_bench: bool,

    /// Benchmarker le taux de rejet du filtre Coppersmith bivarié (δ et ε simultanés).
    #[arg(long)]
    prune_bivar: bool,

    /// Golden Block Test : valider que le bloc solution survit au filtre bivarié m=2.
    #[arg(long)]
    golden_test: bool,

    /// Benchmarker le Kill Switch PNC (GS partiel + LLL f64 early-abort).
    #[arg(long)]
    pnc_bench: bool,

    /// Test du juge de paix : évaluer S₃(x_L, x_R, x_P) sur la vraie racine → doit donner 0.
    #[arg(long)]
    s3test: bool,

    /// Afficher les estimations sans lancer la recherche.
    #[arg(long)]
    estimate_only: bool,

    /// Auto-test : générer k aléatoire, chercher, vérifier.
    #[arg(long)]
    selftest: bool,

    /// Graine hex u64 pour selftest.
    #[arg(long, default_value = "0x135")]
    seed: String,

    /// Lancer le Dispatcher GLV toroïdal (PNC+LLL sur grille scalaire 2D).
    #[arg(long)]
    dispatch: bool,

    /// Bits de la demi-dimension de la grille GLV (défaut : ceil(range_bits/2)+2).
    #[arg(long)]
    half_bits: Option<u32>,

    /// Niveau Jochemsz-May pour la matrice de Macaulay (2=dim15, 3=dim28).
    #[arg(long, default_value = "3")]
    m_level: u32,

    /// Benchmark du Dispatcher : N paires aléatoires, mesure PNC+LLL rejection rate.
    #[arg(long)]
    dispatch_bench: bool,

    /// Nombre de paires pour --dispatch-bench.
    #[arg(long, default_value = "20")]
    dispatch_bench_n: u64,

    /// Test Golden Block Dispatcher : vérifie que la tuile solution survit au filtre.
    #[arg(long)]
    golden_dispatch: bool,

    /// Dispatcher parallèle Rayon (innovations #24 GLV 6-aut + #28 batch + #30 multi-k LLL).
    #[arg(long)]
    parallel: bool,

    /// Dispatcher optimisé : AnchorTable L2 + auto-tune + 6-aut + Rayon (#21+#22+#24+#28+#30).
    #[arg(long)]
    optimized: bool,

    /// Auto-calibration de block_bits (mesure kill_rate, recommande le meilleur block_bits).
    #[arg(long)]
    auto_tune: bool,

    /// Nombre de threads Rayon (défaut : tous les cœurs disponibles).
    #[arg(long, default_value = "0")]
    threads: usize,

    /// Pipeline complet : GLV-4D + Semaev-S₃ + Frobenius + LLL-m3 + 6-aut + AnchorL2 + Rayon.
    #[arg(long)]
    solve: bool,

    /// Lancer le pipeline GSDD complet (Galois Symmetry + Nested Field Decomposition).
    #[arg(long)]
    gsdd: bool,

    /// Selftest GSDD : valide GLV, S₃, Frobenius, Cantor-Zassenhaus, LLL m=3.
    #[arg(long)]
    gsdd_selftest: bool,
}

// ─── Utilitaires scalaires ───────────────────────────────────────────────────

#[inline]
pub fn fe_from_u64(v: u64) -> Fe { [v, 0, 0, 0] }

pub fn in_range(k: Fe, range_bits: u32) -> bool {
    if range_bits >= 256 { return true; }
    let word = (range_bits / 64) as usize;
    let bit  = range_bits % 64;
    let mask = if bit == 0 { !0u64 } else { !((1u64 << bit) - 1) };
    if word < 4 && k[word] & mask != 0 { return false; }
    for i in (word + 1)..4 { if k[i] != 0 { return false; } }
    true
}

// ─── Récupération 6-aut (1D) ─────────────────────────────────────────────────
//
// Collision canonical_x(Q_j) = canonical_x(b·G)
//   ⟹  Q_j = α·(b·G)  →  k = k_giant + α·b
fn recover_k_1d(k_giant: Fe, baby_b: u64, target: Pt, range_bits: u32) -> Option<Fe> {
    let kb = fe_from_u64(baby_b);
    let alphas: [Fe; 6] = [
        [1,0,0,0], sc_neg([1,0,0,0]),
        LAMBDA,    sc_neg(LAMBDA),
        LAMBDA2,   sc_neg(LAMBDA2),
    ];
    for alpha in alphas {
        let k = sc_add(k_giant, sc_mul(alpha, kb));
        if !in_range(k, range_bits) { continue; }
        let pt = scalar_mul(G, k);
        if !pt.inf && pt.x == target.x && pt.y == target.y { return Some(k); }
    }
    None
}

// ─── Récupération 6-aut (2D GLV) ─────────────────────────────────────────────
//
// Q_{j₁,j₂} = P − j₁·M₁·G − j₂·M₂·φG = α·(a₁·G + a₂·φG)
//
// k·G = P
// k = j₁·M₁ + j₂·M₂·λ + α·(a₁ + a₂·λ)  (mod n)
fn recover_k_2d(
    j1_m1: Fe,       // j₁·M₁ (mod n)
    j2_m2_lam: Fe,   // j₂·M₂·λ (mod n)
    a1: u64, a2: u64,
    target: Pt, range_bits: u32,
) -> Option<Fe> {
    let k_giant = sc_add(j1_m1, j2_m2_lam);
    // baby scalar = a₁ + a₂·λ
    let baby = sc_add(fe_from_u64(a1), sc_mul(LAMBDA, fe_from_u64(a2)));
    let alphas: [Fe; 6] = [
        [1,0,0,0], sc_neg([1,0,0,0]),
        LAMBDA,    sc_neg(LAMBDA),
        LAMBDA2,   sc_neg(LAMBDA2),
    ];
    for alpha in alphas {
        let k = sc_add(k_giant, sc_mul(alpha, baby));
        if !in_range(k, range_bits) { continue; }
        let pt = scalar_mul(G, k);
        if !pt.inf && pt.x == target.x && pt.y == target.y { return Some(k); }
    }
    None
}

// ─── Baby table 1D ───────────────────────────────────────────────────────────
//
// Version rapide : Jacobien + batch-inversion (Montgomery trick).
// Réduit le nombre d'inversions de M à M/W (W = 1024).
// ~24× plus rapide que la version naïve (affine + 1 inv/step).
fn build_baby_table_1d(m: u64) -> HashMap<[u64; 4], u64> {
    use secp::{ptj_from_affine, ptj_add_affine, ptj_batch_to_affine, PtJ};
    const W: u64 = 1024;
    let mut table = HashMap::with_capacity(m as usize);
    // Precompute step_W = W*G (used to advance batch starts)
    let step_g = if m < W { G } else { scalar_mul(G, [W, 0, 0, 0]) };
    let mut batch_start_pt = INF; // 0*G = INF (will be advanced to 1*G at first step)
    // batch_start_pt begins at (batch_idx * W) * G, but we start the first batch at 1*G
    // so initial batch_start = 0*G = INF; first point in batch = 1*G = INF + G
    let mut batch_scalar: u64 = 1; // scalar of first point in current batch

    loop {
        if batch_scalar > m { break; }
        let batch_end = (batch_scalar + W - 1).min(m);
        let count = (batch_end - batch_scalar + 1) as usize;

        // Fill batch in Jacobian using mixed J+A additions from batch_start_pt
        let mut jpts: Vec<PtJ> = Vec::with_capacity(count);
        // batch_start_pt is (batch_scalar - 1)*G. First point = batch_scalar*G.
        // Use Jacobian+Affine addition (no inversion needed).
        let mut cur = ptj_add_affine(ptj_from_affine(batch_start_pt), G);
        jpts.push(cur);
        for _ in 1..count {
            cur = ptj_add_affine(cur, G);
            jpts.push(cur);
        }

        // Batch convert to affine (1 inversion per batch)
        let affine_pts = ptj_batch_to_affine(&jpts);

        // Insert into table
        for (i, aff) in affine_pts.iter().enumerate() {
            let b = batch_scalar + i as u64;
            if !aff.inf {
                let cx = canonical_x(aff.x);
                table.entry(cx).or_insert(b);
            }
        }

        // Advance batch_start_pt by W steps
        batch_scalar += W;
        batch_start_pt = pt_add(batch_start_pt, step_g);
    }
    table
}

// ─── Baby table 2D GLV ───────────────────────────────────────────────────────
//
// T[canonical_x((a₁·G + a₂·φG).x)] = (a₁, a₂)
// Construction incrémentale : O(M₁·M₂) additions EC, zéro scalar_mul.
fn build_baby_table_2d(m1: u64, m2: u64) -> HashMap<[u64; 4], (u64, u64)> {
    let phi_g = phi_point(G);
    let capacity = (m1 * m2) as usize;
    let mut table = HashMap::with_capacity(capacity);

    let mut row = INF; // a₁·G (commence à 0·G = INF)
    for a1 in 0..m1 {
        let mut pt = row; // a₁·G + 0·φG
        for a2 in 0..m2 {
            if !pt.inf {
                let cx = canonical_x(pt.x);
                table.entry(cx).or_insert((a1, a2));
            }
            pt = pt_add(pt, phi_g); // a₁·G + (a₂+1)·φG
        }
        row = pt_add(row, G); // (a₁+1)·G
    }
    table
}

// ─── Giant search 1D ─────────────────────────────────────────────────────────

fn giant_search_1d(
    target:     Pt,
    table:      &HashMap<[u64; 4], u64>,
    m:          u64,
    giant_max:  u64,
    range_bits: u32,
) -> Option<Fe> {
    let mg       = scalar_mul(G, fe_from_u64(m));
    let m_scalar = fe_from_u64(m);
    let t0       = Instant::now();
    let mut q        = target;
    let mut k_giant  = [0u64; 4];
    let mut steps    = 0u64;

    for _ in 0..giant_max {
        if !q.inf {
            let cx = canonical_x(q.x);
            if let Some(&b) = table.get(&cx) {
                if let Some(k) = recover_k_1d(k_giant, b, target, range_bits) {
                    eprintln!("\r[giant1D] ✓ step={steps} ({:.2}s)          ",
                        t0.elapsed().as_secs_f64());
                    return Some(k);
                }
            }
        }
        q       = pt_add(q, pt_neg(mg));
        k_giant = sc_add(k_giant, m_scalar);
        steps  += 1;
        if steps & 0xFFFFF == 0 {
            eprint!("\r[giant1D] {steps}  ({:.1}%)  {:.1}M step/s   ",
                steps as f64 / giant_max as f64 * 100.0,
                steps as f64 / t0.elapsed().as_secs_f64() / 1e6);
        }
    }
    eprintln!();
    None
}

// ─── Giant search 2D GLV ──────────────────────────────────────────────────────
//
// Q_{j₁,j₂} = P − j₁·M₁·G − j₂·M₂·φG
// Boucle externe : j₂ (direction φG), interne : j₁ (direction G).
// Avance par soustraction incrémentale : zéro scalar_mul par step.
fn giant_search_2d(
    target:     Pt,
    table:      &HashMap<[u64; 4], (u64, u64)>,
    m1: u64, m2: u64,
    n1: u64,  n2: u64,
    range_bits: u32,
) -> Option<Fe> {
    let m1_g     = scalar_mul(G,          fe_from_u64(m1));
    let m2_phi_g = scalar_mul(phi_point(G), fe_from_u64(m2));
    let m1_sc    = fe_from_u64(m1);
    let m2_lam   = sc_mul(LAMBDA, fe_from_u64(m2)); // M₂·λ

    let max_j1 = n1 / m1 + 2;
    let max_j2 = n2 / m2 + 2;

    let t0    = Instant::now();
    let total = max_j1.saturating_mul(max_j2);
    let mut steps = 0u64;

    let mut q_row     = target;              // P − j₂·M₂·φG  (j₂=0 initially)
    let mut j2_lam_sc = [0u64; 4];          // j₂·M₂·λ mod n

    for _j2 in 0..max_j2 {
        let mut q        = q_row;
        let mut j1_m1_sc = [0u64; 4]; // j₁·M₁ mod n

        for _j1 in 0..max_j1 {
            if !q.inf {
                let cx = canonical_x(q.x);
                if let Some(&(a1, a2)) = table.get(&cx) {
                    if let Some(k) = recover_k_2d(j1_m1_sc, j2_lam_sc, a1, a2, target, range_bits) {
                        eprintln!("\r[giant2D] ✓ step={steps} ({:.2}s)          ",
                            t0.elapsed().as_secs_f64());
                        return Some(k);
                    }
                }
            }
            q        = pt_add(q, pt_neg(m1_g));
            j1_m1_sc = sc_add(j1_m1_sc, m1_sc);
            steps   += 1;
            if steps & 0xFFFFF == 0 {
                eprint!("\r[giant2D] {steps}  ({:.1}%)  {:.1}M step/s   ",
                    steps as f64 / total as f64 * 100.0,
                    steps as f64 / t0.elapsed().as_secs_f64() / 1e6);
            }
        }
        q_row     = pt_add(q_row, pt_neg(m2_phi_g));
        j2_lam_sc = sc_add(j2_lam_sc, m2_lam);
    }
    eprintln!();
    None
}

// ─── Analyse de faisabilité ───────────────────────────────────────────────────

fn print_feasibility(range_bits: u32, baby_bits: u32, glv: bool) {
    let fmt_b = |b: u128| -> String {
        if b >= 1u128<<50 { format!("~2^{}", (b as f64).log2() as u32) }
        else if b >= 1u128<<40 { format!("{:.0} TB", b as f64/(1u128<<40) as f64) }
        else if b >= 1u128<<30 { format!("{:.0} GB", b as f64/(1u128<<30) as f64) }
        else if b >= 1u128<<20 { format!("{:.0} MB", b as f64/(1u128<<20) as f64) }
        else { format!("{} B", b) }
    };
    let fmt_n = |n: u128| -> String {
        if n >= 1u128<<100 { format!("~2^{}", (n as f64).log2() as u32) }
        else if n >= 1_000_000_000 { format!("{:.1}G", n as f64/1e9) }
        else if n >= 1_000_000 { format!("{:.1}M", n as f64/1e6) }
        else if n >= 1_000 { format!("{:.1}K", n as f64/1e3) }
        else { format!("{}", n) }
    };

    println!("╔════════════════════════════════════════════════════════════╗");
    if glv {
        // 2D mode: baby = M₁×M₂ = M², giant = (N₁/M) × (N₂/M)
        // N₁ = range_bits, N₂ depends on GLV bounds
        let k2_bits = glv_k2_bits(range_bits);
        let m = 1u128 << baby_bits;
        let m2 = if k2_bits == 0 { 1u128 }
                 else { m.min(1u128 << k2_bits) };
        let n1_exp = range_bits as i64;
        let n2_exp = k2_bits as i64;
        let g1_exp = (n1_exp - baby_bits as i64).max(0);
        let g2_exp = (n2_exp - baby_bits as i64).max(0);
        let baby_entries = m.saturating_mul(m2);
        let baby_ram = baby_entries.saturating_mul(48);
        println!("║  BSGS-2D GLV secp256k1 — Faisabilité                      ║");
        println!("╠════════════════════════════════════════════════════════════╣");
        println!("║  range_bits = {range_bits:<5}  k ∈ [0, 2^{range_bits})");
        println!("║  baby_bits  = {baby_bits:<5}  M = 2^{baby_bits} par dimension");
        println!("║  Baby table = {} × {} = {}  entrées", fmt_n(m), fmt_n(m2), fmt_n(baby_entries));
        println!("║  RAM baby   = {}", fmt_b(baby_ram));
        println!("║  N₁ = 2^{n1_exp}  (borne k₁)");
        println!("║  N₂ = 2^{n2_exp}  (borne k₂, Babai GLV)");
        println!("║  Giant j₁   ≈ 2^{g1_exp}  steps");
        println!("║  Giant j₂   ≈ 2^{g2_exp}  steps");
        println!("║  Total giant≈ 2^{}  (j₁×j₂)", g1_exp + g2_exp);
        println!("║");
        if k2_bits == 0 {
            println!("║  ℹ k < 2^189 → k₂=0 (GLV trivial) → 2D = 1D sur k₁");
        }
        let ok = baby_ram < 1u128<<40 && (g1_exp + g2_exp) < 35;
        if ok { println!("║  ✓ Faisable"); }
        else {
            if baby_ram >= 1u128<<40 { println!("║  ✗ RAM infaisable : {}", fmt_b(baby_ram)); }
            if g1_exp + g2_exp >= 35 { println!("║  ✗ Giant infaisable : ~2^{}", g1_exp+g2_exp); }
            println!("║     → singraal/ (kangaroo C≈0.55)");
        }
    } else {
        // 1D mode
        let m = 1u128 << baby_bits;
        let giant_exp = range_bits as i64 - baby_bits as i64;
        let baby_ram  = m.saturating_mul(40);
        let giant     = if giant_exp <= 0 { 1u128 } else { 1u128 << giant_exp.min(100) };
        println!("║  BSGS-1D 6-aut secp256k1 — Faisabilité                    ║");
        println!("╠════════════════════════════════════════════════════════════╣");
        println!("║  range_bits = {range_bits:<5}  k ∈ [0, 2^{range_bits})");
        println!("║  baby_bits  = {baby_bits:<5}  M = {}  entrées", fmt_n(m));
        println!("║  RAM baby   = {}", fmt_b(baby_ram));
        println!("║  Giant      ≈ {}  (2^{range_bits}/M)", fmt_n(giant));
        println!("║");
        let ok_ram   = baby_ram < 1u128<<40;
        let ok_giant = giant_exp < 35;
        if ok_ram && ok_giant {
            let total = m + giant;
            let secs  = total as f64 / 1e7;
            if secs < 60.0 { println!("║  ✓ Faisable (~{secs:.1}s CPU)"); }
            else if secs < 3600.0 { println!("║  ✓ Faisable (~{:.1} min CPU)", secs/60.0); }
            else { println!("║  ~ Lent  (~{:.1} h CPU)", secs/3600.0); }
        } else {
            if !ok_ram   { println!("║  ✗ RAM infaisable : {}", fmt_b(baby_ram)); }
            if !ok_giant { println!("║  ✗ Giant infaisable : ~2^{giant_exp}"); }
            println!("║     → singraal/ (kangaroo C≈0.55)");
        }
    }
    println!("╚════════════════════════════════════════════════════════════╝");
}

// ─── Borne N₂ (k₂) selon range_bits ──────────────────────────────────────────
//
// Babai rounding sur secp256k1 : c₁ = ⌊g₁·k/2^384⌋ avec g₁ ≈ 2^189.
// Pour k < 2^189 : c₁ = c₂ = 0 → k₂ = 0.
// Pour k ≥ 2^189 : k₂ ∈ [0, 2^128) environ.
fn glv_k2_bits(range_bits: u32) -> u32 {
    if range_bits <= 189 { 0 }
    else { (range_bits as i32 - 128).max(0) as u32 }
}

// ─── Générateur de clé test ───────────────────────────────────────────────────

fn random_key(seed: u64, range_bits: u32) -> Fe {
    let xs = |mut v: u64| -> u64 { v^=v<<13; v^=v>>7; v^=v<<17; v };
    let w0 = xs(seed ^ 0x9e3779b97f4a7c15);
    let w1 = xs(w0  ^ 0x1234567890abcdef);
    let w2 = xs(w1  ^ 0xfedcba9876543210);
    let w3 = xs(w2  ^ 0x0f1e2d3c4b5a6978);
    let mut k = [w0, w1, w2, w3];
    let mw = (range_bits / 64) as usize;
    let mb = range_bits % 64;
    if mw < 4 {
        k[mw] &= if mb == 0 { 0 } else { (1u64 << mb) - 1 };
        for i in (mw + 1)..4 { k[i] = 0; }
    }
    if range_bits > 0 {
        let hw = ((range_bits - 1) / 64) as usize;
        let hb = (range_bits - 1) % 64;
        if hw < 4 { k[hw] |= 1u64 << hb; }
    }
    k
}

// ─── Semaev Tree Meet-in-the-Middle ──────────────────────────────────────────
//
// k = Σ_{i=0}^{B-1} v_i · 2^(i·block_bits)   v_i ∈ [0, 2^block_bits)
// G_i = 2^(i·block_bits) · G
//
// Gauche (blocs 0..L) :
//   T[canonical_x(Σ v_i·G_i)] = k_gauche    pour tous (v_0,...,v_{L-1})
//
// Droite (blocs L..B) :
//   query = P − Σ v_j·G_j
//   si canonical_x(query) ∈ T → récupération 6-aut : k = α·k_G + k_D

// 2^exp comme scalaire Fe (exp < 256)
fn pow2_fe(exp: u32) -> Fe {
    let mut s = [0u64; 4];
    if exp < 256 {
        s[(exp / 64) as usize] |= 1u64 << (exp % 64);
    }
    s
}

// Récupération 6-aut après collision canonical_x :
//   α(left_sum) + right_sum = P  →  k = α·k_gauche + k_droite
fn semaev_recover_6aut(sc_left: Fe, sc_right: Fe, target: Pt, range_bits: u32) -> Option<Fe> {
    let alphas: [Fe; 6] = [
        [1,0,0,0], sc_neg([1,0,0,0]),
        LAMBDA,    sc_neg(LAMBDA),
        LAMBDA2,   sc_neg(LAMBDA2),
    ];
    for alpha in alphas {
        let k = sc_add(sc_mul(alpha, sc_left), sc_right);
        if !in_range(k, range_bits) { continue; }
        let pt = scalar_mul(G, k);
        if !pt.inf && pt.x == target.x && pt.y == target.y { return Some(k); }
    }
    None
}

// Construction récursive de la table gauche (incrémentale, 0 scalar_mul)
fn semaev_left_recurse(
    table:        &mut HashMap<[u64;4], Fe>,
    base_pts:     &[Pt],
    base_scalars: &[Fe],
    block_size:   u64,
    depth:        usize,
    pt:           Pt,
    sc:           Fe,
) {
    if depth == base_pts.len() {
        if !pt.inf {
            let cx = canonical_x(pt.x);
            table.entry(cx).or_insert(sc);
        }
        return;
    }
    let gi = base_pts[depth];
    let si = base_scalars[depth];
    let mut cur_pt = pt;
    let mut cur_sc = sc;
    for _ in 0..block_size {
        semaev_left_recurse(table, base_pts, base_scalars, block_size, depth + 1, cur_pt, cur_sc);
        cur_pt = pt_add(cur_pt, gi);
        cur_sc = sc_add(cur_sc, si);
    }
}

// Recherche droite récursive (incrémentale, early exit dès solution trouvée)
#[allow(clippy::too_many_arguments)]
fn semaev_right_recurse(
    table:        &HashMap<[u64;4], Fe>,
    neg_pts:      &[Pt],   // -G_{L+i} pour chaque dimension droite
    base_scalars: &[Fe],   // pas scalaire par dimension droite
    block_size:   u64,
    depth:        usize,
    query:        Pt,      // P − partial_right_sum courant
    sc_right:     Fe,
    target:       Pt,
    range_bits:   u32,
    found:        &mut Option<Fe>,
    steps:        &mut u64,
    t0:           &Instant,
    total_steps:  u64,
) {
    if found.is_some() { return; }

    if depth == neg_pts.len() {
        *steps += 1;
        if *steps & 0x3FFFF == 0 {
            eprint!("\r[semaev-R] {}/{total_steps}  ({:.1}%)  {:.2}s   ",
                steps, *steps as f64 / total_steps as f64 * 100.0,
                t0.elapsed().as_secs_f64());
        }
        if !query.inf {
            let cx = canonical_x(query.x);
            if let Some(&sc_left) = table.get(&cx) {
                if let Some(k) = semaev_recover_6aut(sc_left, sc_right, target, range_bits) {
                    *found = Some(k);
                }
            }
        }
        return;
    }
    let gi_neg = neg_pts[depth];
    let si     = base_scalars[depth];
    let mut cur_query = query;
    let mut cur_sc    = sc_right;
    for _ in 0..block_size {
        if found.is_some() { return; }
        semaev_right_recurse(
            table, neg_pts, base_scalars, block_size,
            depth + 1, cur_query, cur_sc, target, range_bits,
            found, steps, t0, total_steps,
        );
        cur_query = pt_add(cur_query, gi_neg);
        cur_sc    = sc_add(cur_sc, si);
    }
}

fn run_semaev(target: Pt, range_bits: u32, block_bits: u32) -> Option<Fe> {
    // Nombre de blocs (toujours pair pour équilibrer gauche/droite)
    let n_blocks = {
        let nb = ((range_bits + block_bits - 1) / block_bits).max(2);
        if nb % 2 == 1 { nb + 1 } else { nb }
    };
    let left_count  = (n_blocks / 2) as usize;
    let right_count = left_count;
    let block_size  = 1u64 << block_bits;

    let left_entries  = (block_size as u128).pow(left_count as u32);
    let right_entries = (block_size as u128).pow(right_count as u32);

    eprintln!("[semaev] range_bits={range_bits}  block_bits={block_bits}  n_blocks={n_blocks}  block_size={block_size}");
    eprintln!("[semaev] gauche : {left_count} blocs → {left_entries} entrées");
    eprintln!("[semaev] droite : {right_count} blocs → {right_entries} requêtes");

    // G_i = 2^(i·block_bits) · G
    let base_pts: Vec<Pt> = (0..n_blocks)
        .map(|i| scalar_mul(G, pow2_fe(i * block_bits)))
        .collect();

    // Scalaires de pas par dimension : 2^(i·block_bits)
    let base_scalars: Vec<Fe> = (0..n_blocks)
        .map(|i| pow2_fe(i * block_bits))
        .collect();

    // ── Table gauche ──────────────────────────────────────────────────────────
    eprintln!("[semaev-L] Construction...");
    let t_left = Instant::now();
    let mut table: HashMap<[u64;4], Fe> = HashMap::with_capacity(left_entries as usize);
    semaev_left_recurse(
        &mut table,
        &base_pts[..left_count],
        &base_scalars[..left_count],
        block_size, 0, INF, [0u64;4],
    );
    eprintln!("[semaev-L] {} entrées en {:.3}s", table.len(), t_left.elapsed().as_secs_f64());

    // ── Recherche droite ──────────────────────────────────────────────────────
    eprintln!("[semaev-R] Recherche...");
    let neg_right_pts: Vec<Pt> = base_pts[left_count..].iter().map(|&p| pt_neg(p)).collect();
    let total_steps = right_entries.min(u64::MAX as u128) as u64;
    let mut found: Option<Fe> = None;
    let mut steps = 0u64;
    let t0 = Instant::now();

    semaev_right_recurse(
        &table, &neg_right_pts, &base_scalars[left_count..],
        block_size, 0, target, [0u64;4], target, range_bits,
        &mut found, &mut steps, &t0, total_steps,
    );
    eprintln!("\r[semaev-R] terminé en {:.3}s  ({steps} requêtes)          ",
        t0.elapsed().as_secs_f64());
    found
}

// ─── Golden Block Test ────────────────────────────────────────────────────────
//
// Valide que le filtre bivarié m=2 laisse SURVIVRE le bloc contenant la vraie clé.
//
// Protocole :
//   k = random_key(seed, range_bits)   → clé connue
//   target = k · G
//   Split 2 blocs : k = v_L · G_L + v_R · G_R
//     G_L = G,  G_R = 2^(range_bits/2) · G
//     v_L = k & mask,  v_R = k >> half
//   x_L = x(v_L · G_L),  x_R = x(v_R · G_R)
//   Vérif : S₃(x_L, x_R, x_P) = 0  ✓
//
//   Block base :  A = x_L - (x_L mod X),  B = x_R - (x_R mod X)
//   Roots dans le bloc :  δ₀ = x_L mod X,  ε₀ = x_R mod X
//   Appel filtre :  is_block_pair_viable(A, B, block_bits)  → doit retourner TRUE
fn run_golden_block_test(seed: u64, range_bits: u32, block_bits: u32) {
    use coppersmith::{fe_to_bigint, s3_bivariate_coeffs, find_glv_coeffs,
                     build_macaulay_bivariate_m2,
                     lll_reduce_bigint, norm_sq_bigint, LatticePruner};
    use num_bigint::BigInt;
    use num_traits::Zero;

    let k = random_key(seed, range_bits);
    let target = scalar_mul(G, k);
    eprintln!("[golden] k     = 0x{}", fe_to_hex(k));
    eprintln!("[golden] x_P   = 0x{}", fe_to_hex(target.x));

    // Split : v_L = lower (range_bits/2) bits, v_R = upper bits
    let half = range_bits / 2;
    let g_r  = scalar_mul(G, pow2_fe(half));

    // v_L : k & ((1<<half)-1)
    let mut v_l = k;
    let wl = (half / 64) as usize;
    let bl = half % 64;
    if wl < 4 { v_l[wl] &= if bl == 0 { 0 } else { (1u64<<bl)-1 }; }
    for i in (wl+1)..4 { v_l[i] = 0; }

    // v_R : k >> half  (multi-word right shift)
    let mut v_r = [0u64; 4];
    let shift_words = (half / 64) as usize;
    let shift_bits  = (half % 64) as u32;
    for i in 0..(4-shift_words) {
        v_r[i] = k[i + shift_words] >> shift_bits;
        if shift_bits > 0 && i + shift_words + 1 < 4 {
            v_r[i] |= k[i + shift_words + 1] << (64 - shift_bits);
        }
    }

    let pt_l = scalar_mul(G,   v_l);
    let pt_r = scalar_mul(g_r, v_r);
    let sum  = pt_add(pt_l, pt_r);

    let ok = !sum.inf && sum.x == target.x && sum.y == target.y;
    eprintln!("[golden] split check (pt_L + pt_R == target) : {ok}");
    if !ok {
        eprintln!("[golden] ERREUR split — abort");
        return;
    }

    let x_l = fe_to_bigint(pt_l.x);
    let x_r = fe_to_bigint(pt_r.x);
    let x_p = fe_to_bigint(target.x);
    let p   = fe_to_bigint(FIELD_P);

    // Vérif S₃(x_L, x_R, x_P) ≡ 0 (mod p) — tester les 9 combinaisons GLV
    let (ei, ej, coeffs_exact) = find_glv_coeffs(&x_l, &x_r, &x_p, &p);
    eprintln!("[golden] S₃(β^{}·x_L, β^{}·x_R, x_P) mod p = {}  (doit être 0)",
        ei, ej, &coeffs_exact[0]);
    if !coeffs_exact[0].is_zero() {
        eprintln!("[golden] AVERTISSEMENT : aucune combinaison GLV ne donne c₀₀=0 — vérifier split");
    }

    // Block boundaries (sur x_L et x_R bruts — les β^k sont appliqués dans find_glv_coeffs)
    let xblk = BigInt::from(1u64) << block_bits as usize;
    let a_base = &x_l - (&x_l % &xblk);
    let b_base = &x_r - (&x_r % &xblk);
    let delta0 = &x_l - &a_base;
    let eps0   = &x_r - &b_base;

    eprintln!("[golden] block_bits={block_bits}  X={}", &xblk);
    eprintln!("[golden] δ₀ = {delta0}  (doit être < {xblk})");
    eprintln!("[golden] ε₀ = {eps0}  (doit être < {xblk})");

    // Calculer manuellement norme LLL vs borne pour le bloc solution
    // Utiliser find_glv_coeffs sur les bases de blocs (β^i·A, β^j·B)
    let (gi, gj, coeffs) = find_glv_coeffs(&a_base, &b_base, &x_p, &p);
    eprintln!("[golden] GLV combo (i={gi}, j={gj}) : c₀₀ = S₃(β^{gi}·A, β^{gj}·B, x_P)");
    eprintln!("[golden] c₀₀ ≠ 0 : {}", !coeffs[0].is_zero());

    let mat = build_macaulay_bivariate_m2(&coeffs, &xblk, &p);

    // ── Sondes diagnostiques avant LLL ───────────────────────────────────────
    let mut max_bits: u64 = 0;
    let mut zero_rows: Vec<usize> = Vec::new();
    for (i, row) in mat.iter().enumerate() {
        let all_zero = row.iter().all(|v| v.is_zero());
        if all_zero { zero_rows.push(i); }
        for val in row {
            let b = val.bits();
            if b > max_bits { max_bits = b; }
        }
    }
    eprintln!("[debug-lll] dim={}×{}  max_coeff_bits={}  zero_rows={:?}",
        mat.len(), mat[0].len(), max_bits, zero_rows);

    // Vérifier la dépendance linéaire basique : lignes proportionnelles ?
    let dim = mat.len();
    let mut colinear_pairs: Vec<(usize,usize)> = Vec::new();
    'outer: for i in 0..dim {
        for j in (i+1)..dim {
            // Vérifie si row[i] et row[j] sont proportionnelles
            let nz_i: Vec<_> = mat[i].iter().enumerate().filter(|(_,v)| !v.is_zero()).collect();
            let nz_j: Vec<_> = mat[j].iter().enumerate().filter(|(_,v)| !v.is_zero()).collect();
            if nz_i.is_empty() || nz_j.is_empty() { continue; }
            if nz_i.len() != nz_j.len() { continue; }
            // ratio du premier terme non-nul
            let (ki, vi) = nz_i[0]; let (kj, vj) = nz_j[0];
            if ki != kj { continue; }
            // ratio = vi/vj ; vérifie tous les termes
            let prop = nz_i.iter().zip(nz_j.iter()).all(|((ai, av), (aj, bv))| {
                ai == aj && vi * (*bv) == vj * (*av)
            });
            if prop { colinear_pairs.push((i, j));
                      if colinear_pairs.len() >= 3 { break 'outer; } }
        }
    }
    if !colinear_pairs.is_empty() {
        eprintln!("[debug-lll] LIGNES COLINÉAIRES : {:?}", colinear_pairs);
    } else {
        eprintln!("[debug-lll] Aucune colinéarité triviale détectée");
    }
    // ── Test A : diagonale ────────────────────────────────────────────────────
    eprintln!("[debug-lll] Diagonale (bits de mat[i][i]) :");
    for i in 0..mat.len() {
        eprintln!("  ligne {:2} : {:4} bits  (val mod 2^32 = {})",
            i, mat[i][i].bits(),
            &mat[i][i] % BigInt::from(1u64 << 32));
    }

    // ── Test B : dump matrice pour SageMath / vérif déterminant ──────────────
    {
        use std::io::Write;
        let mut f = std::fs::File::create("/tmp/macaulay_matrix.txt").unwrap();
        writeln!(f, "# 15x15 Macaulay m=2 pour S3 Semaev secp256k1").unwrap();
        writeln!(f, "# block_bits={block_bits}  seed=0x{seed:x}  range_bits={range_bits}",
            block_bits=block_bits, seed=seed, range_bits=range_bits).unwrap();
        writeln!(f, "M = Matrix(ZZ, [").unwrap();
        for (i, row) in mat.iter().enumerate() {
            let row_str: Vec<String> = row.iter().map(|v| v.to_string()).collect();
            let comma = if i < mat.len()-1 { "," } else { "" };
            writeln!(f, "  [{}]{}", row_str.join(", "), comma).unwrap();
        }
        writeln!(f, "])").unwrap();
        writeln!(f, "print('det =', M.det())").unwrap();
        writeln!(f, "print('rank =', M.rank())").unwrap();
        writeln!(f, "L = M.LLL()").unwrap();
        writeln!(f, "norms = sorted([v.norm() for v in L.rows() if v != 0])").unwrap();
        writeln!(f, "print('shortest norm =', norms[0] if norms else 'N/A')").unwrap();
    }
    eprintln!("[debug-lll] Matrice dumpée dans /tmp/macaulay_matrix.txt");
    eprintln!("[debug-lll] Pour vérifier : sage /tmp/macaulay_matrix.txt");
    eprintln!("[debug-lll]   ou : python3 -c \"exec(open('/tmp/macaulay_matrix.txt').read())\" (avec fpylll)");
    // ─────────────────────────────────────────────────────────────────────────

    let reduced = lll_reduce_bigint(mat);

    let shortest = reduced.iter()
        .map(|row| norm_sq_bigint(row))
        .filter(|n| !n.is_zero())
        .min()
        .unwrap_or_else(BigInt::zero);

    let p2      = &p * &p;
    let bound   = (&p2 * &p2) / 15i64;

    let survives = &shortest < &bound;

    eprintln!("[golden] norm²  ≈ 2^{:.1}", shortest.bits() as f64);
    eprintln!("[golden] bound² ≈ 2^{:.1}  (p⁴/15)", bound.bits() as f64);
    eprintln!("[golden] Golden Block {} → filtre {}",
        if survives { "SURVIT ✓" } else { "TUÉ ✗" },
        if survives { "GO (vert)" } else { "NO-GO (rouge)" });

    // Tester aussi via l'API haut-niveau
    let pruner  = LatticePruner::new(target.x, 8);
    let api_ok  = pruner.is_block_pair_viable(&a_base, &b_base, block_bits);
    eprintln!("[golden] is_block_pair_viable API : {api_ok}");
}

// ─── Test du juge de paix : S₃(x_L, x_R, x_P) == 0 ─────────────────────────
//
// Protocole :
//   k  = random_key(seed, range_bits)   — clé connue
//   P  = k·G                            — cible
//   half = range_bits / 2
//   P_L = v_L · G          où v_L = k & ((1<<half)-1)
//   P_R = v_R · 2^half·G   où v_R = k >> half
//   → P_L + P_R = P  (split exact)
//   S₃(x_L, x_R, x_P) = (x_L-x_R)²·x_P² − 2(x_L+x_R)(x_L·x_R+7)·x_P + (x_L·x_R−7)²
//
// Si S₃ = 0 : arithmétique modulaire OK.
// Si S₃ ≠ 0 : overflow dans fp_mod (ex : (x*y-7)² calculé avant réduction).
fn run_s3_direct_test(seed: u64, range_bits: u32) {
    use coppersmith::{fe_to_bigint, s3_bivariate_coeffs, find_glv_coeffs};
    use num_bigint::BigInt;
    use num_traits::Zero;

    // ── Sanité : G + 2G = 3G  →  S₃(x(G), x(2G), x(3G)) doit être 0 ─────────
    {
        let g1  = G;
        let g2  = pt_dbl_pub(G);
        let g3  = pt_add(g1, g2);
        let p   = fe_to_bigint(FIELD_P);
        let x1  = fe_to_bigint(g1.x);
        let x2  = fe_to_bigint(g2.x);
        let x3  = fe_to_bigint(g3.x);
        let c   = s3_bivariate_coeffs(&x1, &x2, &x3, &p);
        println!("[s3-sanity] G + 2G = 3G : S₃(x(G), x(2G), x(3G)) = {}", &c[0]);
        if c[0].is_zero() {
            println!("[s3-sanity] ✓ formule S₃ correcte");
        } else {
            println!("[s3-sanity] ✗ BUG dans la formule S₃ elle-même (nb bits={})", c[0].bits());
        }
    }
    println!();

    let k      = random_key(seed, range_bits);
    let target = scalar_mul(G, k);
    let half   = range_bits / 2;
    let g_r    = scalar_mul(G, pow2_fe(half));

    // v_L = k & ((1 << half) - 1)
    let mut v_l = k;
    let wl = (half / 64) as usize;
    let bl = half % 64;
    if wl < 4 { v_l[wl] &= if bl == 0 { 0 } else { (1u64 << bl) - 1 }; }
    for i in (wl + 1)..4 { v_l[i] = 0; }

    // v_R = k >> half
    let mut v_r = [0u64; 4];
    let sw = (half / 64) as usize;
    let sb = (half % 64) as u32;
    for i in 0..(4 - sw) {
        v_r[i] = k[i + sw] >> sb;
        if sb > 0 && i + sw + 1 < 4 {
            v_r[i] |= k[i + sw + 1] << (64 - sb);
        }
    }

    let pt_l = scalar_mul(G,   v_l);
    let pt_r = scalar_mul(g_r, v_r);
    let sum  = pt_add(pt_l, pt_r);

    println!("[s3-test] k      = 0x{}", fe_to_hex(k));
    println!("[s3-test] x_P    = 0x{}", fe_to_hex(target.x));
    println!("[s3-test] x_L    = 0x{}", fe_to_hex(pt_l.x));
    println!("[s3-test] x_R    = 0x{}", fe_to_hex(pt_r.x));

    let split_ok = !sum.inf && sum.x == target.x && sum.y == target.y;
    println!("[s3-test] split P_L + P_R == P : {}", if split_ok { "✓" } else { "✗ ERREUR split" });
    if !split_ok {
        eprintln!("[s3-test] ABORT — split incorrect, impossible de tester S₃");
        return;
    }

    let p   = fe_to_bigint(FIELD_P);
    let x_l = fe_to_bigint(pt_l.x);
    let x_r = fe_to_bigint(pt_r.x);
    let x_p = fe_to_bigint(target.x);

    // Évaluation directe de S₃(x_L, x_R, x_P) — c00 de s3_bivariate_coeffs(A=x_L, B=x_R)
    let coeffs = s3_bivariate_coeffs(&x_l, &x_r, &x_p, &p);
    let c00 = &coeffs[0];

    println!("[s3-test] S₃(x_L, x_R, x_P) mod p = {}", c00);
    if c00.is_zero() {
        println!("[s3-test] ✓ CORRECT — S₃ = 0 (arithmétique modulaire OK)");
    } else {
        println!("[s3-test] ✗ BUG — S₃ ≠ 0");
        println!("[s3-test]   bits de la valeur = {}", c00.bits());
        println!("[s3-test]   → vérifier fp_mod sur (x*y-7)² ou (x-y)²");
    }

    // Essayer les 9 combinaisons GLV (β^i·x_L, β^j·x_R) — l'une doit donner 0
    println!();
    println!("[s3-test] Recherche combinaison GLV (β^i·x_L, β^j·x_R) avec S₃=0 :");
    let (gi, gj, coeffs_glv) = find_glv_coeffs(&x_l, &x_r, &x_p, &p);
    println!("[s3-test]   β^{}·x_L, β^{}·x_R → S₃ = {}", gi, gj, &coeffs_glv[0]);
    if coeffs_glv[0].is_zero() {
        println!("[s3-test] ✓ Combinaison GLV correcte trouvée (i={gi}, j={gj})");
    } else {
        println!("[s3-test] ✗ Aucune combinaison GLV ne donne S₃=0 — bug dans split ou β");
    }
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let args = Args::parse();

    let baby_bits = args.baby_bits.unwrap_or_else(|| {
        let bb = (args.range_bits / 2).max(4).min(24);
        eprintln!("[auto] baby_bits = {bb}  (M = 2^{bb})");
        bb
    });

    // ── GSDD Selftest ─────────────────────────────────────────────────────────
    if args.gsdd_selftest {
        let seed_s = args.seed.trim_start_matches("0x");
        let seed   = u64::from_str_radix(seed_s, 16).expect("--seed: hex u64");
        let block_bits = args.block_bits.unwrap_or(8);
        gsdd::selftest_gsdd(args.range_bits, block_bits, seed);
        return;
    }

    // ── Analyse LLL (optionnelle) ─────────────────────────────────────────────
    if args.lll {
        lll::print_lll_report(args.range_bits);
    }

    // ── Test juge de paix S₃ ────────────────────────────────────────────────
    if args.s3test {
        let seed_s = args.seed.trim_start_matches("0x");
        let seed   = u64::from_str_radix(seed_s, 16).expect("--seed: hex u64");
        run_s3_direct_test(seed, args.range_bits);
        return;
    }

    // ── Golden Block Test ────────────────────────────────────────────────────
    if args.golden_test {
        let block_bits = args.block_bits.unwrap_or(5);
        let seed_s = args.seed.trim_start_matches("0x");
        let seed   = u64::from_str_radix(seed_s, 16).expect("--seed: hex u64");
        run_golden_block_test(seed, args.range_bits, block_bits);
        return;
    }

    // ── Benchmark filtre PNC (Kill Switch Cohen 2.6.7 exact) ───────────────
    if args.pnc_bench {
        let block_bits = args.block_bits.unwrap_or(5);
        eprintln!("[pnc-bench] Kill Switch PNC  block_bits={block_bits}");
        let stats = lll_earlyabort::benchmark_killswitch(15, 1000);
        stats.print();
    }

    // ── Benchmark filtre Coppersmith univarié ────────────────────────────────
    if args.prune_bench {
        let block_bits = args.block_bits.unwrap_or(5);
        let dummy_x = GX;
        let pruner = coppersmith::LatticePruner::new(dummy_x, 8);
        let n_blocks = 1000u64;
        eprintln!("[prune-bench] univarié  block_bits={block_bits}  dim=8  n_blocks={n_blocks}");
        let t0 = std::time::Instant::now();
        let rate = pruner.benchmark_rejection_rate(block_bits, n_blocks, args.range_bits);
        eprintln!("[prune-bench] Taux de rejet : {:.1}%  ({:.2}s)",
            rate * 100.0, t0.elapsed().as_secs_f64());
        eprintln!("  > 99% → filtre efficace");
        eprintln!("  ~20%  → baseline univarié observé");
        if !args.estimate_only { println!(); }
    }

    // ── Benchmark filtre Coppersmith bivarié ─────────────────────────────────
    if args.prune_bivar {
        let block_bits = args.block_bits.unwrap_or(5);
        let dummy_x = GX;
        let pruner = coppersmith::LatticePruner::new(dummy_x, 8);
        let n_blocks = 1000u64;
        eprintln!("[prune-bivar] bivarié  block_bits={block_bits}  dim=15 (m=2)  n_blocks={n_blocks}");
        eprintln!("[prune-bivar] lattice Jochemsz-May m=2, S₃(A+δ, B+ε, x_P), δ,ε ∈ [0,2^{block_bits})");
        let t0 = std::time::Instant::now();
        let rate = pruner.benchmark_bivariate_rejection_rate(block_bits, n_blocks, args.range_bits);
        eprintln!("[prune-bivar] Taux de rejet : {:.1}%  ({:.2}s)",
            rate * 100.0, t0.elapsed().as_secs_f64());
        eprintln!("[prune-bivar] Interprétation :");
        eprintln!("  > 99% → chaque paire (bloc_G, bloc_D) prouvée vide en O(dim³) LLL");
        eprintln!("  ~50%  → gain ×2 sur univarié");
        eprintln!("  < 20% → bivarié moins efficace qu'univarié (bornes HG trop larges)");
        if !args.estimate_only { println!(); }
    }

    if !args.semaev {
        print_feasibility(args.range_bits, baby_bits, args.glv);
        println!();
    }

    // ── Benchmark Kill Switch PNC ─────────────────────────────────────────────
    if args.pnc_bench {
        let block_bits = args.block_bits.unwrap_or(5);
        eprintln!("[pnc-bench] Kill Switch PNC — GS partiel f64 + LLL f64 early-abort");
        eprintln!("[pnc-bench] dim=15 (m=2 bivarié)  block_bits={block_bits}");
        eprintln!("[pnc-bench] Borne HG : log2(p⁴/15) = {:.1}", lll_earlyabort::hg_bound_log2(15));
        let stats = lll_earlyabort::benchmark_killswitch(15, 200);
        stats.print();
        eprintln!("[pnc-bench] Interprétation :");
        eprintln!("  > 80% → Kill Switch efficace (économie majeure de LLL BigRational)");
        eprintln!("  ~0%   → matrices toutes viables (Kill Switch conservateur)");
        if !args.estimate_only { println!(); }
    }

    // ── Benchmark Dispatcher ─────────────────────────────────────────────────
    if args.dispatch_bench {
        let block_bits = args.block_bits.unwrap_or(5);
        dispatcher::benchmark_dispatcher(
            args.range_bits, block_bits, args.m_level, args.dispatch_bench_n
        );
        if args.estimate_only { return; }
    }

    if args.estimate_only { return; }

    // ── Cible ─────────────────────────────────────────────────────────────────
    let (target, expected_k) = if args.selftest {
        let seed_s = args.seed.trim_start_matches("0x");
        let seed   = u64::from_str_radix(seed_s, 16)
            .expect("--seed: entier hex u64");
        let k  = random_key(seed, args.range_bits);
        let pt = scalar_mul(G, k);
        eprintln!("[selftest] k     = 0x{}", fe_to_hex(k));
        eprintln!("[selftest] k·G.x = 0x{}", fe_to_hex(pt.x));

        // Afficher la décomposition GLV pour diagnostic
        let (k1, k2) = glv_decompose(k);
        eprintln!("[glv-decomp] k₁ = 0x{}", fe_to_hex(k1));
        eprintln!("[glv-decomp] k₂ = 0x{}", fe_to_hex(k2));
        let k_check = sc_add(k1, sc_mul(LAMBDA, k2));
        eprintln!("[glv-check]  k₁+λk₂ == k : {}", k_check == k);

        (pt, Some(k))
    } else {
        let tx = fe_from_hex(&args.target_x).expect("--target-x: 64 hex chars");
        let ty = fe_from_hex(&args.target_y).expect("--target-y: 64 hex chars");
        (Pt { x: tx, y: ty, inf: false }, None)
    };

    // ── Paramètres ────────────────────────────────────────────────────────────
    let m    = 1u64 << baby_bits;
    let n_val = if args.range_bits >= 64 {
        u64::MAX
    } else {
        (1u64 << args.range_bits).saturating_sub(1)
    };

    // ── Dispatch 1D / 2D / Semaev ────────────────────────────────────────────
    // ── Golden Dispatch Test ──────────────────────────────────────────────────
    if args.golden_dispatch {
        let block_bits = args.block_bits.unwrap_or(5);
        let seed_s = args.seed.trim_start_matches("0x");
        let seed   = u64::from_str_radix(seed_s, 16).expect("--seed: hex u64");
        let k  = random_key(seed, args.range_bits);
        eprintln!("[golden-dispatch] k = 0x{}", fe_to_hex(k));
        let passed = dispatcher::golden_block_test(k, args.range_bits, block_bits, args.m_level);
        if passed { eprintln!("[golden-dispatch] ✓ PASS"); }
        else       { eprintln!("[golden-dispatch] ✗ FAIL"); }
        return;
    }

    // Configurer Rayon avant tout dispatch parallèle
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap_or(());
    }

    // ── Auto-tune block_bits ──────────────────────────────────────────────────
    if args.auto_tune {
        let n_samples = args.dispatch_bench_n;
        let result = dispatcher::auto_tune_block_bits(
            args.range_bits, args.m_level, n_samples, 0.999
        );
        result.print();
        return;
    }

    let result = if args.solve {
        // ── FULL STACK : toutes innovations branchées ─────────────────────────
        let block_bits = args.block_bits.unwrap_or_else(|| {
            let bb = ((args.range_bits + 3) / 4).max(3).min(20);
            eprintln!("[auto] solve block_bits = {bb}");
            bb
        });
        let half_bits = args.half_bits.unwrap_or_else(|| {
            let hb = (args.range_bits / 2 + 2).min(64);
            eprintln!("[auto] solve half_bits = {hb}");
            hb
        });
        let mut cfg = dispatcher::DispatcherConfig::new(target, args.range_bits, block_bits);
        cfg.half_bits = half_bits;
        cfg.m_level   = 3;  // LLL m=3 dim=28 (Jochemsz-May calibré)
        cfg.verbose   = true;
        dispatcher::run_full_stack(&cfg)
    } else if args.gsdd {
        // GSDD = Full Stack avec m=3 + verbose (alias pédagogique de --solve)
        let block_bits = args.block_bits.unwrap_or_else(|| {
            let bb = ((args.range_bits + 3) / 4).max(3).min(20);
            eprintln!("[auto] gsdd block_bits = {bb}");
            bb
        });
        let half_bits = args.half_bits.unwrap_or_else(|| {
            (args.range_bits / 2 + 2).min(64)
        });
        let mut cfg = dispatcher::DispatcherConfig::new(target, args.range_bits, block_bits);
        cfg.half_bits = half_bits;
        cfg.m_level   = 3;
        cfg.verbose   = true;
        eprintln!("[main] Mode GSDD — Galois Symmetry + Nested Field Decomposition (m=3)");
        dispatcher::run_full_stack(&cfg)
    } else if args.optimized || args.parallel || args.dispatch {
        let block_bits = args.block_bits.unwrap_or_else(|| {
            let bb = ((args.range_bits + 3) / 4).max(3).min(22);
            eprintln!("[auto] dispatcher block_bits = {bb}");
            bb
        });
        let half_bits = args.half_bits.unwrap_or_else(|| {
            let hb = (args.range_bits / 2 + 2).min(64);
            eprintln!("[auto] dispatcher half_bits = {hb}");
            hb
        });
        let mut cfg = dispatcher::DispatcherConfig::new(target, args.range_bits, block_bits);
        cfg.half_bits = half_bits;
        cfg.m_level   = args.m_level;
        cfg.verbose   = true;
        if args.optimized {
            eprintln!("[main] Mode dispatcher OPTIMISÉ (#21 L2 + #22 auto-tune + #24+#28+#30)");
            dispatcher::run_dispatcher_optimized(&cfg)
        } else if args.parallel {
            eprintln!("[main] Mode dispatcher PARALLÈLE Rayon (innovations #24+#28+#30)");
            dispatcher::run_dispatcher_parallel(&cfg)
        } else {
            dispatcher::run_dispatcher(&cfg)
        }
    } else if args.semaev {
        let block_bits = args.block_bits.unwrap_or_else(|| {
            let bb = ((args.range_bits + 3) / 4).max(3).min(20);
            eprintln!("[auto] block_bits = {bb}");
            bb
        });
        run_semaev(target, args.range_bits, block_bits)
    } else if args.glv {
        let k2_bits = glv_k2_bits(args.range_bits);
        let n2 = if k2_bits == 0 { 1u64 }
                 else if k2_bits >= 64 { u64::MAX }
                 else { 1u64 << k2_bits };
        let m2 = m.min(n2.max(1));

        eprintln!("[bsgs2d-glv] M₁={m}  M₂={m2}  N₁≈2^{}  N₂≈2^{k2_bits}",
            args.range_bits);
        eprintln!("[baby2D] Construction ({} × {} = {} entrées)...", m, m2, m*m2);
        let t_baby = Instant::now();
        let table  = build_baby_table_2d(m, m2);
        eprintln!("[baby2D] {} entrées en {:.2}s",
            table.len(), t_baby.elapsed().as_secs_f64());
        println!();

        eprintln!("[giant2D] Recherche...");
        giant_search_2d(target, &table, m, m2, n_val, n2, args.range_bits)
    } else {
        let giant_max = n_val / m + 2;
        eprintln!("[bsgs1D] M={m}  giant_max≈{giant_max}");
        eprintln!("[baby1D] Construction ({m} entrées)...");
        let t_baby = Instant::now();
        let table  = build_baby_table_1d(m);
        eprintln!("[baby1D] {} entrées en {:.2}s",
            table.len(), t_baby.elapsed().as_secs_f64());
        println!();

        eprintln!("[giant1D] Recherche (max {giant_max} steps)...");
        giant_search_1d(target, &table, m, giant_max, args.range_bits)
    };  // end result

    // ── Résultat ──────────────────────────────────────────────────────────────
    match result {
        Some(k) => {
            println!();
            println!("╔══════════════════════════════════════════════════════════════════╗");
            println!("║  SOLUTION TROUVÉE                                                ║");
            println!("║  k = {}  ║", fe_to_hex(k));
            println!("╚══════════════════════════════════════════════════════════════════╝");
            let check = scalar_mul(G, k);
            if !check.inf && check.x == target.x && check.y == target.y {
                println!("[VÉRIFIÉ] k·G == cible ✓");
            } else {
                eprintln!("[ERREUR] k·G ≠ cible — bug dans recover_k");
            }
            if let Some(ek) = expected_k {
                if ek == k { println!("[SELFTEST ✓] k trouvé == k attendu"); }
                else {
                    eprintln!("[SELFTEST ✗]");
                    eprintln!("  attendu : 0x{}", fe_to_hex(ek));
                    eprintln!("  trouvé  : 0x{}", fe_to_hex(k));
                }
            }
        }
        None => {
            eprintln!("[INFO] Aucune solution dans l'espace parcouru.");
            if let Some(ek) = expected_k {
                eprintln!("[SELFTEST ✗] k attendu = 0x{}", fe_to_hex(ek));
            }
        }
    }
}
