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

use clap::Parser;
use secp::*;
use std::collections::HashMap;
use std::time::Instant;

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

    /// Afficher les estimations sans lancer la recherche.
    #[arg(long)]
    estimate_only: bool,

    /// Auto-test : générer k aléatoire, chercher, vérifier.
    #[arg(long)]
    selftest: bool,

    /// Graine hex u64 pour selftest.
    #[arg(long, default_value = "0x135")]
    seed: String,
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

fn build_baby_table_1d(m: u64) -> HashMap<[u64; 4], u64> {
    let mut table = HashMap::with_capacity(m as usize);
    let mut pt = G;
    for b in 1..=m {
        if !pt.inf {
            let cx = canonical_x(pt.x);
            table.entry(cx).or_insert(b);
        }
        pt = pt_add(pt, G);
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

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let args = Args::parse();

    let baby_bits = args.baby_bits.unwrap_or_else(|| {
        let bb = (args.range_bits / 2).max(4).min(24);
        eprintln!("[auto] baby_bits = {bb}  (M = 2^{bb})");
        bb
    });

    if !args.semaev {
        print_feasibility(args.range_bits, baby_bits, args.glv);
        println!();
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
    let result = if args.semaev {
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
    };

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
