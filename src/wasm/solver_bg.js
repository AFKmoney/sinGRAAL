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
     * @param {number} num_vars
     */
    constructor(num_vars) {
        const ret = wasm.wasmsolver_new(num_vars);
        this.__wbg_ptr = ret;
        WasmSolverFinalization.register(this, this.__wbg_ptr, this);
        return this;
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
export function __wbindgen_init_externref_table() {
    const table = wasm.__wbindgen_externrefs;
    const offset = table.grow(4);
    table.set(0, undefined);
    table.set(offset + 0, undefined);
    table.set(offset + 1, null);
    table.set(offset + 2, true);
    table.set(offset + 3, false);
}
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

let WASM_VECTOR_LEN = 0;


let wasm;
export function __wbg_set_wasm(val) {
    wasm = val;
}
