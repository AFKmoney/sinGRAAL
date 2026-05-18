/* tslint:disable */
/* eslint-disable */

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
    constructor(num_vars: number);
    /**
     * Run up to max_conflicts CDCL steps. Returns JSON status string.
     */
    step(max_conflicts: number): string;
}
