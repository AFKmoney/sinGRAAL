export class WasmKangaroo {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmKangarooFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmkangaroo_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    dps() {
        const ret = wasm.wasmkangaroo_dps(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Create a new Kangaroo for ECDLP: find k such that k·G = (target_x, target_y)
     * range_low, target_x, target_y: 64-char hex strings (no 0x prefix)
     * range_bits: 2^range_bits is the search range size
     * dp_bits: number of leading zero bits for distinguished points (higher = fewer DPs, faster per DP)
     * @param {string} target_x
     * @param {string} target_y
     * @param {string} range_low
     * @param {number} range_bits
     * @param {number} num_tames
     * @param {number} dp_bits
     */
    constructor(target_x, target_y, range_low, range_bits, num_tames, dp_bits) {
        const ptr0 = passStringToWasm0(target_x, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(target_y, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passStringToWasm0(range_low, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len2 = WASM_VECTOR_LEN;
        const ret = wasm.wasmkangaroo_new(ptr0, len0, ptr1, len1, ptr2, len2, range_bits, num_tames, dp_bits);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        WasmKangarooFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Show the 6-orbit for a given point (educational / visualization)
     * Returns JSON array: [{name, x, y}, ...]
     * @param {string} x
     * @param {string} y
     * @returns {string}
     */
    static orbit_demo(x, y) {
        let deferred4_0;
        let deferred4_1;
        try {
            const ptr0 = passStringToWasm0(x, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ptr1 = passStringToWasm0(y, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            const ret = wasm.wasmkangaroo_orbit_demo(ptr0, len0, ptr1, len1);
            var ptr3 = ret[0];
            var len3 = ret[1];
            if (ret[3]) {
                ptr3 = 0; len3 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred4_0 = ptr3;
            deferred4_1 = len3;
            return getStringFromWasm0(ptr3, len3);
        } finally {
            wasm.__wbindgen_free(deferred4_0, deferred4_1, 1);
        }
    }
    /**
     * Compute k·G and return hex (x, y) – useful for testing
     * @param {string} k_hex
     * @returns {string}
     */
    static scalar_mul_hex(k_hex) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(k_hex, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.wasmkangaroo_scalar_mul_hex(ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * @returns {boolean}
     */
    solved() {
        const ret = wasm.wasmkangaroo_solved(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Step the engine by `iterations` per animal per call.
     * Returns JSON: {status, steps, dps, key, wildX, tameX}
     * @param {number} iterations
     * @returns {string}
     */
    step(iterations) {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.wasmkangaroo_step(this.__wbg_ptr, iterations);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @returns {number}
     */
    steps() {
        const ret = wasm.wasmkangaroo_steps(this.__wbg_ptr);
        return ret;
    }
    /**
     * Verify that k·G = target (sanity check after solve)
     * @param {string} k_hex
     * @param {string} target_x
     * @param {string} target_y
     * @returns {boolean}
     */
    static verify(k_hex, target_x, target_y) {
        const ptr0 = passStringToWasm0(k_hex, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(target_x, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passStringToWasm0(target_y, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len2 = WASM_VECTOR_LEN;
        const ret = wasm.wasmkangaroo_verify(ptr0, len0, ptr1, len1, ptr2, len2);
        return ret !== 0;
    }
}
if (Symbol.dispose) WasmKangaroo.prototype[Symbol.dispose] = WasmKangaroo.prototype.free;

export class WasmSolver {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmSolverFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmsolver_free(ptr, 0);
    }
    /**
     * Add a clause. Pass literals as a flat JS array, e.g. [1, -2, 3]
     * @param {Int32Array} lits
     */
    add_clause(lits) {
        const ptr0 = passArray32ToWasm0(lits, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.wasmsolver_add_clause(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * Add multiple clauses packed as flat array with 0 separators.
     * e.g. [1, -2, 0, 3, 4, 0] = two clauses [1,-2] and [3,4]
     * @param {Int32Array} data
     */
    add_clauses_packed(data) {
        const ptr0 = passArray32ToWasm0(data, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.wasmsolver_add_clauses_packed(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * Number of assigned variables
     * @returns {number}
     */
    assigned_count() {
        const ret = wasm.wasmsolver_assigned_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of clauses (including learned)
     * @returns {number}
     */
    clause_count() {
        const ret = wasm.wasmsolver_clause_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * After SAT: get assignment as flat i8 array indexed by var-1.
     * val[i] = 1 (true), -1 (false), 0 (unset)
     * @returns {Int8Array}
     */
    get_assignment() {
        const ret = wasm.wasmsolver_get_assignment(this.__wbg_ptr);
        var v1 = getArrayI8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
    }
    /**
     * Return all unit-forced literals as JSON array [{var, val}].
     * These are bits of k that CDCL has proven must have a specific value.
     * @param {Int32Array} k_lits
     * @returns {string}
     */
    get_forced_bits(k_lits) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ptr0 = passArray32ToWasm0(k_lits, wasm.__wbindgen_malloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.wasmsolver_get_forced_bits(this.__wbg_ptr, ptr0, len0);
            deferred2_0 = ret[0];
            deferred2_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * @param {number} num_vars
     */
    constructor(num_vars) {
        const ret = wasm.wasmsolver_new(num_vars);
        this.__wbg_ptr = ret;
        WasmSolverFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Number of bits proven at level-0 (invariant constraints from CDCL)
     * @param {Int32Array} k_lits
     * @returns {number}
     */
    proven_bit_count(k_lits) {
        const ptr0 = passArray32ToWasm0(k_lits, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmsolver_proven_bit_count(this.__wbg_ptr, ptr0, len0);
        return ret >>> 0;
    }
    /**
     * Run up to max_conflicts CDCL steps. Returns JSON status string.
     * @param {number} max_conflicts
     * @returns {string}
     */
    step(max_conflicts) {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.wasmsolver_step(this.__wbg_ptr, max_conflicts);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
}
if (Symbol.dispose) WasmSolver.prototype[Symbol.dispose] = WasmSolver.prototype.free;
export function __wbg___wbindgen_throw_9c31b086c2b26051(arg0, arg1) {
    throw new Error(getStringFromWasm0(arg0, arg1));
}
export function __wbindgen_cast_0000000000000001(arg0, arg1) {
    // Cast intrinsic for `Ref(String) -> Externref`.
    const ret = getStringFromWasm0(arg0, arg1);
    return ret;
}
export function __wbindgen_init_externref_table() {
    const table = wasm.__wbindgen_externrefs;
    const offset = table.grow(4);
    table.set(0, undefined);
    table.set(offset + 0, undefined);
    table.set(offset + 1, null);
    table.set(offset + 2, true);
    table.set(offset + 3, false);
}
const WasmKangarooFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmkangaroo_free(ptr, 1));
const WasmSolverFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmsolver_free(ptr, 1));

function getArrayI8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getInt8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

let cachedInt8ArrayMemory0 = null;
function getInt8ArrayMemory0() {
    if (cachedInt8ArrayMemory0 === null || cachedInt8ArrayMemory0.byteLength === 0) {
        cachedInt8ArrayMemory0 = new Int8Array(wasm.memory.buffer);
    }
    return cachedInt8ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    return decodeText(ptr >>> 0, len);
}

let cachedUint32ArrayMemory0 = null;
function getUint32ArrayMemory0() {
    if (cachedUint32ArrayMemory0 === null || cachedUint32ArrayMemory0.byteLength === 0) {
        cachedUint32ArrayMemory0 = new Uint32Array(wasm.memory.buffer);
    }
    return cachedUint32ArrayMemory0;
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function passArray32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getUint32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_externrefs.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();

if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function (arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {
            read: arg.length,
            written: buf.length
        };
    };
}

let WASM_VECTOR_LEN = 0;


let wasm;
export function __wbg_set_wasm(val) {
    wasm = val;
}
