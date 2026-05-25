// sinGRAAL — BSGS avec 6-automorphismes secp256k1
//
// ─── Algorithme ───────────────────────────────────────────────────────────────
//
// Objectif : trouver k ∈ [0, 2^range_bits) tel que k·G = P.
//
// Décomposition : k = j·M + b,  b ∈ [0,M),  j ∈ [0, N/M)  (N = 2^range_bits)
//   P − j·M·G = b·G  →  chercher b dans la table baby.
//
// Baby table T :
//   T[canonical_x(b·G)] = b   pour b ∈ [1, M]
//   canonical_x(Q) = min(x(Q), β·x(Q) mod p, β²·x(Q) mod p)
//   → 6 automorphismes {±id, ±ψ, ±ψ²} donnent la même canonical_x
//   → la clé de table détecte aussi les faux-positifs via recover_k_6aut
//
// Giant step :
//   Q₀ = P,  Q_{j+1} = Q_j − M·G,  j = 0,1,...,N/M−1
//   Pour chaque j : vérifier canonical_x(Q_j) ∈ T
//   Collision : Q_j = α·b·G  →  k = j·M + α·b (mod n),  tester k·G == P
//
// Complexité :
//   Baby  : M  additions EC
//   Giant : N/M  additions EC  (N = 2^range_bits)
//   Optimal : M = √N  →  total ≈ 2·√N  ops
//   Mémoire : M × 40 octets
//
// ─── Rôle des 6-automorphismes ───────────────────────────────────────────────
// canonical_x relie les 6 scalaires {b, n−b, λb, n−λb, λ²b, n−λ²b} (mod n)
// au même representant x-coordinate.
//
// Pour range_bits ≪ 128 : λ ≈ 2^128, donc λb et λ²b sont >> N pour b < M=√N.
// → seul α=±1 peut donner une collision dans [0,N). Benefit: vérification
//   rapide de 6 candidats dans recover_k_6aut, sans coût de stockage extra.
//
// Pour range_bits > 128 : les orbites λ et λ² tombent aussi dans [0,N), le
// baby table peut alors être réduit à M = √(N/6) entrées (vrai gain 6×).
//
// ─── Faisabilité ─────────────────────────────────────────────────────────────
//   range_bits=40 : M=2^20  → table 40 MB,   giant 2^20  → ~secondes  ✓
//   range_bits=50 : M=2^25  → table 1.3 GB,  giant 2^25  → ~minutes   ✓
//   range_bits=60 : M=2^30  → table 43 GB,   giant 2^30  → ~heures    ⚠
//
//   puzzle #135 (N=2^135) : M=2^67 → 590 exaoctets → impossible.
//   → Utiliser singraal/ (kangaroo, C≈0.55, ~2^65 ops).

mod secp;

use clap::Parser;
use secp::*;
use std::collections::HashMap;
use std::time::Instant;

// ─── CLI ─────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(name = "bsgs2d", about = "BSGS 6-aut GLV secp256k1 — utile jusqu'à ~60 bits")]
struct Args {
    /// x-coordinate de la cible (64 hex chars). Omis en mode --selftest.
    #[arg(long, default_value = "")]
    target_x: String,

    /// y-coordinate de la cible (64 hex chars). Omis en mode --selftest.
    #[arg(long, default_value = "")]
    target_y: String,

    /// k ∈ [0, 2^range_bits).
    #[arg(long, default_value = "40")]
    range_bits: u32,

    /// M = 2^baby_bits baby steps. Table ≈ M × 40 octets.
    /// Défaut optimal : ceil(range_bits/2 - 1) pour minimiser le total.
    #[arg(long)]
    baby_bits: Option<u32>,

    /// Afficher les estimations sans lancer la recherche.
    #[arg(long)]
    estimate_only: bool,

    /// Auto-test : générer k, chercher, vérifier.
    #[arg(long)]
    selftest: bool,

    /// Graine hex u64 pour selftest.
    #[arg(long, default_value = "0x135")]
    seed: String,
}

// ─── Utilitaires ─────────────────────────────────────────────────────────────

#[inline]
fn fe_from_u64(v: u64) -> Fe { [v, 0, 0, 0] }

fn in_range(k: Fe, range_bits: u32) -> bool {
    if range_bits >= 256 { return true; }
    let word = (range_bits / 64) as usize;
    let bit  = range_bits % 64;
    let mask = if bit == 0 { !0u64 } else { !((1u64 << bit) - 1) };
    if word < 4 && k[word] & mask != 0 { return false; }
    for i in (word + 1)..4 { if k[i] != 0 { return false; } }
    true
}

// ─── Récupération 6-automorphismes ───────────────────────────────────────────
//
// Collision : canonical_x(Q_j) == canonical_x(b·G)
//   ⟹ Q_j = α·(b·G)  pour un α ∈ {1,−1,λ,−λ,λ²,−λ²}
//
// P = j·M·G + Q_j = j·M·G + α·b·G = (j·M + α·b)·G  (mod n)
//
// On essaie les 6 α, filtre in_range (élimine 5/6 en O(1)), vérifie k·G == P.
fn recover_k_6aut(k_giant: Fe, baby_b: u64, target: Pt, range_bits: u32) -> Option<Fe> {
    let k_baby = fe_from_u64(baby_b);
    let candidates: [Fe; 6] = [
        [1, 0, 0, 0],
        sc_neg([1, 0, 0, 0]),
        LAMBDA,
        sc_neg(LAMBDA),
        LAMBDA2,
        sc_neg(LAMBDA2),
    ];
    for alpha in candidates {
        let k = sc_add(k_giant, sc_mul(alpha, k_baby));
        if !in_range(k, range_bits) { continue; }
        let pt = scalar_mul(G, k);
        if !pt.inf && pt.x == target.x && pt.y == target.y {
            return Some(k);
        }
    }
    None
}

// ─── Phase baby ───────────────────────────────────────────────────────────────
//
// T[canonical_x(b·G)] = b  pour b = 0,...,M-1.
// Avance par additions successives (pas de scalar_mul par step).
fn build_baby_table(m: u64) -> HashMap<[u64; 4], u64> {
    let mut table = HashMap::with_capacity(m as usize);
    let mut pt = G; // commencer à 1·G (b=0 → point à l'infini, pas de DP utile)
    // b=0 : 0·G = INF → on démarre à b=1 pour éviter le point à l'infini
    for b in 1..=m {
        if !pt.inf {
            let cx = canonical_x(pt.x);
            table.entry(cx).or_insert(b);
        }
        pt = pt_add(pt, G);
    }
    table
}

// ─── Phase giant ──────────────────────────────────────────────────────────────
//
// Q₀ = P,  Q_{j+1} = Q_j − M·G
// Scalaire k_giant = j·M  maintenu incrémentalement (sc_add, pas de scalar_mul).
fn giant_search(
    target:     Pt,
    table:      &HashMap<[u64; 4], u64>,
    m:          u64,
    giant_max:  u64,  // nombre de giant steps
    range_bits: u32,
) -> Option<Fe> {
    // M·G : un seul scalar_mul.
    let m_g   = scalar_mul(G, fe_from_u64(m));
    let m_mod = fe_from_u64(m); // M mod n pour incrémenter k_giant

    let t0  = Instant::now();
    let mut q       = target;           // Q_j = P − j·M·G
    let mut k_giant = [0u64; 4];       // j·M  (mod n)
    let mut steps   = 0u64;

    for _j in 0..giant_max {
        if !q.inf {
            let cx = canonical_x(q.x);
            if let Some(&b) = table.get(&cx) {
                if let Some(k) = recover_k_6aut(k_giant, b, target, range_bits) {
                    eprintln!(
                        "\r[giant] ✓ collision step={steps} ({:.2}s)          ",
                        t0.elapsed().as_secs_f64()
                    );
                    return Some(k);
                }
            }
        }

        // Q_{j+1} = Q_j − M·G
        q       = pt_add(q, pt_neg(m_g));
        k_giant = sc_add(k_giant, m_mod);
        steps  += 1;

        if steps & 0xFFFFF == 0 {
            let secs = t0.elapsed().as_secs_f64();
            let pct  = steps as f64 / giant_max as f64 * 100.0;
            eprint!("\r[giant] {steps}  ({pct:.2}%)  {:.1}M step/s   ",
                steps as f64 / secs / 1e6);
        }
    }
    eprintln!();
    None
}

// ─── Analyse de faisabilité ───────────────────────────────────────────────────

fn print_feasibility(range_bits: u32, baby_bits: u32) {
    // M = 2^baby_bits baby steps, N = 2^range_bits
    // Giant = N / M = 2^(range_bits - baby_bits)
    let m = 1u128 << baby_bits;
    let giant_exp = range_bits as i64 - baby_bits as i64;
    let giant = if giant_exp <= 0 { 1u128 } else { 1u128 << giant_exp.min(100) };
    let baby_ram = m.saturating_mul(40); // ~40 octets par entrée (clé 32B + valeur 8B)

    let fmt_b = |b: u128| -> String {
        if b >= 1u128<<40 { format!("{:.0} TB", b as f64/(1u128<<40) as f64) }
        else if b >= 1u128<<30 { format!("{:.0} GB", b as f64/(1u128<<30) as f64) }
        else if b >= 1u128<<20 { format!("{:.0} MB", b as f64/(1u128<<20) as f64) }
        else { format!("{} B", b) }
    };
    let fmt_n = |n: u128| -> String {
        if n >= 1u128<<100 { format!("~2^{}", n.ilog2()) }
        else if n >= 1_000_000_000 { format!("{:.1}G", n as f64/1e9) }
        else if n >= 1_000_000 { format!("{:.1}M", n as f64/1e6) }
        else if n >= 1_000 { format!("{:.1}K", n as f64/1e3) }
        else { format!("{}", n) }
    };

    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║  BSGS-6aut GLV secp256k1 — Faisabilité                    ║");
    println!("╠════════════════════════════════════════════════════════════╣");
    println!("║  range_bits = {range_bits:<5}  k ∈ [0, 2^{range_bits})");
    println!("║  baby_bits  = {baby_bits:<5}  M = 2^{baby_bits} = {}  baby steps", fmt_n(m));
    println!("║  RAM baby   = {}  ({} entrées × 40 oct)", fmt_b(baby_ram), fmt_n(m));
    println!("║  Giant      ≈ {}  (2^{} / M)", fmt_n(giant), range_bits);

    let ok_ram   = baby_ram   < 1u128 << 40;  // < 1 TB
    let ok_giant = giant_exp  < 35;            // < 34G steps

    println!("║");
    if ok_ram && ok_giant {
        let total = m + giant;
        let secs  = total as f64 / 1e7; // ~10M affine add/s CPU
        if secs < 60.0 { println!("║  ✓ Faisable  (~{secs:.1}s CPU)"); }
        else if secs < 3600.0 { println!("║  ✓ Faisable  (~{:.1} min CPU)", secs/60.0); }
        else if secs < 86400.0*2.0 { println!("║  ✓ Faisable  (~{:.1} h CPU)", secs/3600.0); }
        else { println!("║  ✓ Faisable  (~{:.1} j CPU)", secs/86400.0); }
    } else {
        if !ok_ram   { println!("║  ✗ RAM infaisable : baby table = {}", fmt_b(baby_ram)); }
        if !ok_giant { println!("║  ✗ Giant infaisable : ~2^{} steps", giant_exp); }
        println!("║     → Utiliser singraal/ (kangaroo, C≈0.55).");
    }
    println!("╚════════════════════════════════════════════════════════════╝");
}

// ─── Générateur de clé test ───────────────────────────────────────────────────

fn random_key(seed: u64, range_bits: u32) -> Fe {
    let xs = |mut v: u64| -> u64 { v^=v<<13; v^=v>>7; v^=v<<17; v };
    let w0 = xs(seed ^ 0x9e3779b97f4a7c15);
    let w1 = xs(w0  ^ 0x1234567890abcdef);
    let w2 = xs(w1  ^ 0xfedcba9876543210);
    let w3 = xs(w2  ^ 0x0f1e2d3c4b5a6978);
    let mut k = [w0, w1, w2, w3];
    // Masquer au-dessus de range_bits
    let mw = (range_bits / 64) as usize;
    let mb = range_bits % 64;
    if mw < 4 {
        k[mw] &= if mb == 0 { 0 } else { (1u64 << mb) - 1 };
        for i in (mw + 1)..4 { k[i] = 0; }
    }
    // Forcer le bit le plus haut pour être exactement range_bits bits
    if range_bits > 0 {
        let hw = ((range_bits - 1) / 64) as usize;
        let hb = (range_bits - 1) % 64;
        if hw < 4 { k[hw] |= 1u64 << hb; }
    }
    k
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let args = Args::parse();

    // baby_bits optimal : M = sqrt(N) → baby_bits = range_bits / 2
    let baby_bits = args.baby_bits.unwrap_or_else(|| {
        let bb = (args.range_bits / 2).max(4).min(26);
        eprintln!("[auto] baby_bits = {bb}  (M = 2^{bb} = {})", 1u64 << bb);
        bb
    });

    print_feasibility(args.range_bits, baby_bits);
    println!();

    if args.estimate_only { return; }

    // ── Cible ─────────────────────────────────────────────────────────────────
    let (target, expected_k) = if args.selftest {
        let seed_s = args.seed.trim_start_matches("0x");
        let seed   = u64::from_str_radix(seed_s, 16)
            .expect("--seed: entier hex u64 (ex: 0x1a2b3c)");
        let k  = random_key(seed, args.range_bits);
        let pt = scalar_mul(G, k);
        eprintln!("[selftest] k     = 0x{}", fe_to_hex(k));
        eprintln!("[selftest] k·G.x = 0x{}", fe_to_hex(pt.x));
        (pt, Some(k))
    } else {
        let tx = fe_from_hex(&args.target_x).expect("--target-x: 64 hex chars");
        let ty = fe_from_hex(&args.target_y).expect("--target-y: 64 hex chars");
        (Pt { x: tx, y: ty, inf: false }, None)
    };

    // ── Paramètres ────────────────────────────────────────────────────────────
    let m = 1u64 << baby_bits;
    // Giant steps = N / M.  j = floor(k/M) ≤ N/M for any k ∈ [0,N).
    let n_val = if args.range_bits >= 64 { u64::MAX } else { (1u64 << args.range_bits).saturating_sub(1) };
    let giant_max = n_val / m + 2;

    eprintln!("[bsgs2d] M={m}  giant_max≈{giant_max}  range_bits={}", args.range_bits);

    // Garde-fou RAM
    let baby_ram = m as u128 * 40;
    if baby_ram > 1u128 << 37 { // > 128 GB
        eprintln!("[WARN] Baby table = {} GB — peut dépasser la RAM.",
            baby_ram >> 30);
    }

    // ── Phase baby ────────────────────────────────────────────────────────────
    eprintln!("[baby] Construction ({m} entrées)...");
    let t_baby = Instant::now();
    let table  = build_baby_table(m);
    eprintln!("[baby] {} entrées ({:.1}%) en {:.2}s",
        table.len(),
        table.len() as f64 / m as f64 * 100.0,
        t_baby.elapsed().as_secs_f64());
    println!();

    // ── Phase giant ───────────────────────────────────────────────────────────
    eprintln!("[giant] Recherche (max {giant_max} steps)...");
    let t_giant = Instant::now();
    let result  = giant_search(target, &table, m, giant_max, args.range_bits);
    let elapsed = t_giant.elapsed().as_secs_f64();

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
                println!("[VÉRIFIÉ] k·G == cible ✓   ({elapsed:.3}s)");
            } else {
                eprintln!("[ERREUR] k·G ≠ cible — bug dans recover_k_6aut");
            }
            if let Some(ek) = expected_k {
                if ek == k {
                    println!("[SELFTEST ✓] k trouvé == k attendu");
                } else {
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
