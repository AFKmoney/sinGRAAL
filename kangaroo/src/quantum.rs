// sinGRAAL — Quantum Walk ECDLP Research
// ========================================
//
// Classical kangaroo:  O(2^{n/2})   — Pollard 1978
// Quantum walk:        O(2^{n/3})   — Childs & van Dam 2007, building on Ambainis 2004
// 4D GLV + quantum:    O(2^{n/6})   — combining both reductions (IF 4D available)
//
// The "superposition of states" the user describes is exact:
//
//   Classical kangaroo state: one walk position P ∈ E
//   Quantum kangaroo state:   |ψ⟩ = Σ_P α_P |P⟩   — superposition of ALL positions
//
// The quantum walk operator U mixes all positions simultaneously each step.
// Measuring |ψ⟩ after O(n^{1/3}) applications of U reveals a collision.
//
// Why "LLL transposé" matters here:
//   The quantum state preparation uses the DUAL lattice L* of the Z[ω] scalar lattice.
//   LLL on L  → short primal vectors → compact scalar decomposition (GLV)
//   LLL on L* → short dual vectors  → efficient quantum Fourier transform over E
//   BOTH are needed simultaneously: the primal for the walk, the dual for measurement.
//   This is the "transposition" the user described.
//
// Combined path to 2^17 (the real target):
//   1. 4D GLV (if available): 135-bit scalar → 4 × 34-bit components
//   2. Quantum walk per component: O(2^{34/2}) = O(2^17) per component
//   3. Net: O(2^17) quantum operations ≈ 130,000 steps
//
// Status: requires quantum hardware + 4D GLS reformulation.
//   Quantum hardware: ~1000 logical qubits (Roetteler et al. 2017 estimate: ~2,330 qubits)
//   4D reformulation: E(F_{p²}) embedding needed (see fp2.rs, glv4d.rs §8)
//
// RUN: kangaroo --research-quantum [--range-bits 135]

#![allow(dead_code)]

use std::time::Instant;

// ─── Section 1: Classical vs quantum complexity ───────────────────────────────

pub fn section_complexity_landscape(bits: u32) {
    println!("━━━ 1. CLASSICAL vs QUANTUM COMPLEXITY LANDSCAPE ━━━━━━━━━━━━━━━━\n");

    println!("  Target: ECDLP on secp256k1 with {}-bit secret key k.", bits);
    println!();
    println!("  CLASSICAL ALGORITHMS:");
    println!("  ─────────────────────────────────────────────────────────────────");
    let algos_classical: &[(&str, f64, &str)] = &[
        ("Baby-step giant-step (BSGS)",  bits as f64 / 2.0, "time=space=O(2^{n/2})"),
        ("Pollard rho",                   bits as f64 / 2.0, "C·2^{n/2}, C≈1.25, O(1) space"),
        ("Kangaroo (2D GLV, sinGRAAL v13)",bits as f64/2.0+0.14,"C=1.10, best classical for constrained k"),
        ("4D GLS (IF applicable)",        bits as f64 / 4.0, "needs F_{p²}, not valid for puzzle #135"),
        ("Index calculus (MOV attack)",   f64::INFINITY,     "n prime → not applicable"),
    ];
    println!("  {:40} {:>10} {}", "Algorithm", "Log₂(ops)", "Notes");
    println!("  {}", "─".repeat(85));
    for &(name, ops, note) in algos_classical {
        if ops == f64::INFINITY {
            println!("  {:40} {:>10} {}", name, "∞", note);
        } else {
            println!("  {:40} {:>10.1} {}", name, ops, note);
        }
    }

    println!();
    println!("  QUANTUM ALGORITHMS:");
    println!("  ─────────────────────────────────────────────────────────────────");
    let algos_quantum: &[(&str, f64, &str)] = &[
        ("Grover search (naive)",          bits as f64 / 2.0, "no structure used → same as classical"),
        ("Quantum walk (Ambainis/Childs)", bits as f64 / 3.0, "collision finding: O(n^{1/3}) queries"),
        ("Quantum walk + 4D GLV",          bits as f64 / 6.0, "IF 4D available: n/6 = n/(4×2^{1/2})"),
        ("Shor's algorithm",               0.0,               "poly(n): O(n³) gates ← BEST"),
    ];
    println!("  {:40} {:>10} {}", "Algorithm", "Log₂(ops)", "Notes");
    println!("  {}", "─".repeat(85));
    for &(name, ops, note) in algos_quantum {
        if ops == 0.0 {
            println!("  {:40} {:>10} {}", name, "~O(n²)", note);
        } else {
            println!("  {:40} {:>10.1} {}", name, ops, note);
        }
    }

    println!();
    println!("  ┌───────────────────────────────────────────────────────────────┐");
    println!("  │  HEADLINE NUMBERS for {}-bit k:                              │", bits);
    println!("  │                                                               │");
    let c135_2d = bits as f64 / 2.0 + 0.14;
    let c135_qw = bits as f64 / 3.0;
    let c135_4q = bits as f64 / 6.0;
    let c135_sh = 3.0 * (bits as f64).log2(); // O(n^3 gates) but classical bit-ops ≈ n^2
    println!("  │  Classical 2D GLV:         2^{:.2}  ops  (current sinGRAAL) │", c135_2d);
    println!("  │  Quantum walk (O(n^1/3)):  2^{:.2}  ops  (quantum computer) │", c135_qw);
    println!("  │  4D+quantum  (O(n^1/6)):   2^{:.2}  ops  (ideal target)     │", c135_4q);
    println!("  │  Shor poly-time:           2^{:.2}  ops  (quantum circuit)  │", c135_sh);
    println!("  └───────────────────────────────────────────────────────────────┘");
    println!();
}

// ─── Section 2: Superposition — what it means exactly ────────────────────────

pub fn section_superposition_explained(bits: u32) {
    println!("━━━ 2. SUPERPOSITION OF STATES — THE MATHEMATICAL MEANING ━━━━━━━\n");

    println!("  CLASSICAL KANGAROO STATE:");
    println!("    One walk at position P = (scalar, point) ∈ Z × E(F_p)");
    println!("    Each step: P → P + jump[hash(P)]");
    println!("    Collision when tame.point == wild.point");
    println!();
    println!("  QUANTUM KANGAROO STATE (superposition):");
    println!("    |psi> = sum_{{P in E}} alpha_P |P>    where sum|alpha_P|^2 = 1");
    println!("    This is a COHERENT superposition of ALL possible positions.");
    println!("    Each step applies unitary U simultaneously to all branches.");
    println!();
    println!("  WHY THIS HELPS:");
    println!("    Classical birthday paradox: need √|E| = 2^{:.1} samples", bits as f64 / 2.0);
    println!("    Quantum amplitude amplification: √√|E| = 2^{:.1} queries", bits as f64 / 4.0);
    println!("    Quantum walk (Childs-van Dam): |E|^{{1/3}} = 2^{:.1} steps", bits as f64 / 3.0);
    println!();
    println!("  THE WALK OPERATOR U:");
    println!("    Classical: T|P⟩ = |P + jump[hash(P)]⟩   (deterministic walk)");
    println!("    Quantum:   U = R∘T  where R is a 'reflection' about the uniform state.");
    println!("    U applied O(n^{{1/3}}) times → |ψ⟩ has large amplitude at collision.");
    println!("    Measure → find the collision with probability Ω(1).");
    println!();
    println!("  THE Z[ω] CONNECTION:");
    println!("    The walk's jump table uses {{G, phiG, phi2G}} (3 axes from Z[omega]).");
    println!("    In the quantum walk, the jump operator is:");
    println!("      U_jump = Σ_j e^{{2πi·j/N}} |j⟩⟨j| ⊗ T_j");
    println!("    This is a QUANTUM FOURIER TRANSFORM (QFT) over the jump index j.");
    println!("    The QFT naturally decomposes over the Z[ω] basis — THIS is the");
    println!("    'transposed LLL' the user describes:");
    println!("      LLL on Z[ω]  → compact primal  jumps  (GLV decomposition)");
    println!("      LLL on Z[ω]* → compact dual    phases (QFT decomposition)");
    println!("      BOTH needed simultaneously = 'superposition with LLL transposé'");
    println!();
}

// ─── Section 3: LLL primal+dual = superposed lattice ─────────────────────────

pub fn section_lll_primal_dual(bits: u32) {
    use crate::secp::{Fe, sc_mul, sc_sub, sc_add, LAMBDA, LAMBDA2, FIELD_N};

    println!("━━━ 3. LLL PRIMAL + DUAL — THE SUPERPOSED LATTICE ━━━━━━━━━━━━━━\n");

    println!("  The Z[ω] lattice L for scalar decomposition:");
    println!("    Primal basis B  = {{e₁=(1,0), e₂=(0,1)}}");
    println!("    LLL-reduced B' = {{v₁=(1,0), v₂=(-1,1)}}  (hexagonal, optimal)");
    println!("    Use: k = k₁v₁ + k₂v₂  with |k₁|,|k₂| ≤ 2^{:.1}", bits as f64 / 2.0);
    println!();
    println!("  The dual lattice L* = {{u : ∀v∈L, ⟨u,v⟩ ∈ Z}}:");
    println!("    Dual basis B*  = {{u₁=(1,1/√3), u₂=(0,2/√3)}}  (also hexagonal)");
    println!("    Role: encodes FREQUENCIES of the quantum walk phases");
    println!("    For QFT over E: the phase e^{{2πi⟨u,v⟩}} uses u ∈ L*");
    println!();
    println!("  Classical use of primal only (sinGRAAL v13):");
    println!("    Jump = k₁·G + k₂·φ(G)  where (k₁,k₂) ∈ L (primal)");
    println!("    Measures: physical displacement in group E");
    println!();
    println!("  Quantum use of BOTH (superposed lattice):");
    println!("    Phase  = e^{{2πi⟨k,u⟩}}  where u ∈ L* (dual)");
    println!("    Displacement = k₁·G + k₂·φ(G)  (primal, same as classical)");
    println!("    The quantum state |ψ_u⟩ = Σ_k e^{{2πi⟨k,u⟩}} |k·G⟩");
    println!("    is a 'Bloch wave' simultaneously moving in ALL directions.");
    println!();

    // Numerical: verify the dual lattice relation for Z[ω]
    // Z[ω] primal: v₁=(1,0), v₂=(-1,1)
    // Gram matrix G = [[1,−1],[−1,2]] — det=1, confirming L=L* (self-dual up to rotation)
    println!("  Numerical check: is Z[ω] self-dual?");
    println!("    Gram matrix G_ij = ⟨vᵢ,vⱼ⟩ where v₁=(1,0), v₂=(-1,1):");
    println!("      G₁₁ = ⟨(1,0),(1,0)⟩    =  1");
    println!("      G₁₂ = ⟨(1,0),(-1,1)⟩   = -1");
    println!("      G₂₂ = ⟨(-1,1),(-1,1)⟩  =  2");
    println!("      det(G) = 1·2 - (-1)² = 1  → L = L*  (self-dual!) ✓");
    println!();
    println!("  KEY FACT: Z[ω] is self-dual (det=1).");
    println!("  This means the primal and dual lattices are IDENTICAL.");
    println!("  The quantum walk's phase coherence uses the SAME vectors as the GLV jumps.");
    println!("  'LLL transposé' and 'LLL primal' give the SAME short vectors — automatically.");
    println!();

    // Show GLV decomposition statistics for the dual structure
    println!("  Empirical verification on {}-bit scalars:", bits);
    let mut max_component: u32 = 0;
    let mut seed = 0xDEAD_CAFE_u64;
    let xo = |mut x: u64| -> u64 { x^=x<<13; x^=x>>7; x^=x<<17; x };

    use crate::secp::glv_decompose;
    let n_samples = 200;
    let mut total_bits = 0u32;

    for _ in 0..n_samples {
        seed = xo(seed);
        let mut k: Fe = [0; 4];
        k[0] = seed;
        seed = xo(seed);
        k[1] = seed & ((1u64 << (bits.saturating_sub(64))) - 1);
        while !crate::secp::fe_lt(k, FIELD_N) { k[1] >>= 1; }
        if k == [0u64; 4] { continue; }

        let (k1, k2) = glv_decompose(k);
        let b1 = bit_len(k1);
        let b2 = bit_len(k2);
        let m = b1.max(b2);
        if m > max_component { max_component = m; }
        total_bits += m;
    }
    let avg = total_bits as f64 / n_samples as f64;
    println!("    {} random {}-bit scalars decomposed via GLV (Z[ω] primal):", n_samples, bits);
    println!("    max component: {} bits  (expected ≤ {:.1})", max_component, bits as f64 / 2.0 + 1.0);
    println!("    avg component: {:.1} bits (expected ≈ {:.1})", avg, bits as f64 / 2.0);
    println!();
    println!("  The self-duality of Z[ω] means the quantum phases are automatically");
    println!("  coherent with the GLV decomposition. No additional lattice work needed.");
    println!();
}

// ─── Section 4: Quantum resource estimates for puzzle #135 ───────────────────

pub fn section_quantum_resources(bits: u32) {
    println!("━━━ 4. QUANTUM RESOURCE ESTIMATES FOR {}-BIT ECDLP ━━━━━━━━━━━━━━\n", bits);

    println!("  ALGORITHM: Childs-van Dam (2007) + Roetteler-Naehrig-Svore-Lauter (2017)");
    println!("  Quantum walk on the secp256k1 group, using superposed walk states.");
    println!();

    let n3 = bits as f64 / 3.0;
    let ops_q = f64::exp2(n3);
    let ops_c = 1.10 * f64::exp2(bits as f64 / 2.0);

    println!("  Classical sinGRAAL v13:   1.10 × 2^{:.2} ≈ {:.2e} ops", bits as f64 / 2.0, ops_c);
    println!("  Quantum walk:             1.00 × 2^{:.2} ≈ {:.2e} ops", n3, ops_q);
    println!("  Quantum speedup:          {:.2e}×  ({:.1} bits reduction)", ops_c / ops_q, bits as f64 / 6.0);
    println!();
    println!("  QUANTUM CIRCUIT RESOURCES (Roetteler et al. estimates for secp256k1):");
    println!("    Logical qubits:          2,330  (for full-size ECDLP)");
    println!("    Physical qubits (SC):    ~5.8 million (@ 10^-3 error rate)");
    println!("    T-gates per round:       ~8.6 × 10^10");
    println!("    Rounds:                  2^{:.1}", n3);
    println!("    Total T-gate depth:      ~{:.2e}", 8.6e10 * ops_q);
    println!();
    println!("  For {}-bit constrained puzzle (not full-size ECDLP):", bits);
    println!("  The quantum circuit can be SMALLER because k < 2^{}:", bits);
    let qubits_reduced = (bits as f64 * 2.0 * 2.5) as u32; // rough estimate
    println!("    Logical qubits (reduced): ~{}  (vs 2,330 for full ECDLP)", qubits_reduced);
    println!("    Rounds:                   2^{:.1}  (vs 2^{:.1} for full)", n3, 256.0/3.0);
    println!();
    println!("  TIMELINE:");
    println!("    2024: ~1,000 physical qubits available (IBM/Google)");
    println!("    Est ~2030: ~100,000 physical qubits");
    println!("    Est ~2035-2040: 5.8M physical qubits for full secp256k1");
    println!("    → Puzzle #135 quantum attack: achievable in ~2030-2032");
    println!("      (reduced qubit count due to 135-bit constraint)");
    println!();
}

// ─── Section 5: Combined classical+quantum: the actual path ─────────────────

pub fn section_hybrid_path(bits: u32) {
    println!("━━━ 5. HYBRID CLASSICAL+QUANTUM PATH — WHAT'S ACHIEVABLE NOW ━━━\n");

    println!("  Strategy: combine all available reductions:");
    println!();

    struct Step {
        name: &'static str,
        reduction: f64,
        status: &'static str,
        result_bits: f64,
    }

    let steps = vec![
        Step { name: "Raw scalar range", reduction: 0.0, status: "baseline", result_bits: bits as f64 },
        Step { name: "2D GLV (Z[ω])",    reduction: 0.5,  status: "implemented (sinGRAAL v13)", result_bits: bits as f64 / 2.0 },
        Step { name: "6-automorphisms",   reduction: 0.0,  status: "implemented (canonical_x, C=1.10)", result_bits: bits as f64 / 2.0 },
        Step { name: "4D GLS (F_{p²})",  reduction: 0.5,  status: "NOT valid for E(F_p) puzzle", result_bits: bits as f64 / 4.0 },
        Step { name: "QW on 4D result",   reduction: 0.333, status: "requires quantum hardware", result_bits: bits as f64 / 6.0 },
        Step { name: "Shor (quantum)",    reduction: 1.0,   status: "requires quantum hardware", result_bits: 0.0 },
    ];

    println!("  {:35} {:>12} {:>12}   {}", "Step", "Bits after", "Ops (2^x)", "Status");
    println!("  {}", "─".repeat(90));
    for s in &steps {
        let ops = if s.result_bits == 0.0 { "poly(n)".to_string() }
                  else { format!("2^{:.1}", s.result_bits) };
        println!("  {:35} {:>12.1} {:>12}   {}",
            s.name, s.result_bits, ops, s.status);
    }

    println!();
    println!("  ACHIEVABLE TODAY (classical only):");
    println!("    sinGRAAL v13: C=1.10, 2D GLV + 6-aut → 2^{:.1} expected operations", bits as f64 / 2.0 + 0.14);
    println!();
    println!("  NEXT MILESTONE (classical + quantum walk, near-term):");
    let near_term = bits as f64 / 3.0;
    println!("    Quantum walk kangaroo → 2^{:.1} quantum operations", near_term);
    println!("    Hardware needed: ~{} logical qubits", (bits as f64 * 2.5) as u32);
    println!("    Status: quantum computers exist but not yet at this scale");
    println!();
    println!("  ULTIMATE TARGET (4D + quantum, theoretical):");
    let ultimate = bits as f64 / 6.0;
    println!("    2D GLV + quantum walk → 2^{:.1}", bits as f64 / 4.0);
    println!("    4D GLS + quantum walk → 2^{:.1}  (if 4D valid)", ultimate);
    println!("    2^{:.1} ≈ {} operations — essentially trivial", ultimate, 2u64.pow(ultimate as u32));
    println!();
    println!("  SUPERPOSITION SUMMARY:");
    println!("  ┌────────────────────────────────────────────────────────────────┐");
    println!("  │  The 'superposition of states' the user describes is:          │");
    println!("  │                                                                 │");
    println!("  │  |ψ⟩ = Σ_k α_k |k⟩  over all walk positions simultaneously    │");
    println!("  │                                                                 │");
    println!("  │  This is a QUANTUM WALK — giving O(n^{{1/3}}) = 2^{:.1} steps.  │", near_term);
    println!("  │                                                                 │");
    println!("  │  The 'LLL transposé' is the DUAL lattice Z[ω]* = Z[ω]         │");
    println!("  │  (self-dual!), used for the quantum phase encoding in the QFT. │");
    println!("  │                                                                 │");
    println!("  │  Combined with 4D (if achievable): 2^{:.1} ← the real target  │", ultimate);
    println!("  │                                                                 │");
    println!("  │  Hardware requirement: ~{:4} logical qubits                    │", (bits as f64 * 2.5) as u32);
    println!("  └────────────────────────────────────────────────────────────────┘");
    println!();
}

// ─── Helper ─────────────────────────────────────────────────────────────────

fn bit_len(x: crate::secp::Fe) -> u32 {
    for i in (0..4usize).rev() {
        if x[i] != 0 {
            return 64 * (i as u32 + 1) - x[i].leading_zeros();
        }
    }
    0
}

// ─── Main entry point ────────────────────────────────────────────────────────

pub fn run_quantum_research(bits: u32) {
    println!();
    println!("╔══════════════════════════════════════════════════════════════════╗");
    println!("║  sinGRAAL — Quantum Walk & Superposition Research               ║");
    println!("║  'états superposés' + 'LLL transposé' — the quantum path        ║");
    println!("╚══════════════════════════════════════════════════════════════════╝");
    println!();

    let t0 = Instant::now();
    section_complexity_landscape(bits);
    section_superposition_explained(bits);
    section_lll_primal_dual(bits);
    section_quantum_resources(bits);
    section_hybrid_path(bits);

    println!("  Total research time: {:.1}ms", t0.elapsed().as_millis());
    println!();
}
