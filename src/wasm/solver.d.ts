/* tslint:disable */
/* eslint-disable */

export class WasmKangaroo {
    free(): void;
    [Symbol.dispose](): void;
    dps(): number;
    /**
     * Create a new Kangaroo for ECDLP: find k such that k·G = (target_x, target_y)
     * range_low, target_x, target_y: 64-char hex strings (no 0x prefix)
     * range_bits: 2^range_bits is the search range size
     * dp_bits: number of leading zero bits for distinguished points (higher = fewer DPs, faster per DP)
     */
    constructor(target_x: string, target_y: string, range_low: string, range_bits: number, num_tames: number, dp_bits: number);
    /**
     * Show the 6-orbit for a given point (educational / visualization)
     * Returns JSON array: [{name, x, y}, ...]
     */
    static orbit_demo(x: string, y: string): string;
    /**
     * Compute k·G and return hex (x, y) – useful for testing
     */
    static scalar_mul_hex(k_hex: string): string;
    solved(): boolean;
    /**
     * Step the engine by `iterations` per animal per call.
     * Returns JSON: {status, steps, dps, key, wildX, tameX}
     */
    step(iterations: number): string;
    steps(): number;
    /**
     * Verify that k·G = target (sanity check after solve)
     */
    static verify(k_hex: string, target_x: string, target_y: string): boolean;
}

export class WasmSolver {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Add a clause. Pass literals as a flat JS array, e.g. [1, -2, 3]
     */
    add_clause(lits: Int32Array): void;
    /**
     * Add multiple clauses packed as flat array with 0 separators.
     * e.g. [1, -2, 0, 3, 4, 0] = two clauses [1,-2] and [3,4]
     */
    add_clauses_packed(data: Int32Array): void;
    /**
     * Number of assigned variables
     */
    assigned_count(): number;
    /**
     * Number of clauses (including learned)
     */
    clause_count(): number;
    /**
     * After SAT: get assignment as flat i8 array indexed by var-1.
     * val[i] = 1 (true), -1 (false), 0 (unset)
     */
    get_assignment(): Int8Array;
    /**
     * Return all unit-forced literals as JSON array [{var, val}].
     * These are bits of k that CDCL has proven must have a specific value.
     */
    get_forced_bits(k_lits: Int32Array): string;
    constructor(num_vars: number);
    /**
     * Number of bits proven at level-0 (invariant constraints from CDCL)
     */
    proven_bit_count(k_lits: Int32Array): number;
    /**
     * Run up to max_conflicts CDCL steps. Returns JSON status string.
     */
    step(max_conflicts: number): string;
}
