// ─── Types ────────────────────────────────────────────────────────────────────

export type Lit = number;       // positive = variable, negative = negation
type Clause = Lit[];
type Word = Lit[];              // 32 bits, LSB-first (index 0 = bit 0 of uint32)

// ─── Variable pool ────────────────────────────────────────────────────────────

let _varCount = 0;
let _TRUE_LIT: Lit = 0;
let _FALSE_LIT: Lit = 0;

function initPool() {
  _varCount = 0;
  const tv = freshVar();
  _TRUE_LIT = tv;
  _FALSE_LIT = -tv;
}

function freshVar(): number { return ++_varCount; }
function totalVars(): number { return _varCount; }

const TRUE_LIT = (): Lit => _TRUE_LIT;
const FALSE_LIT = (): Lit => _FALSE_LIT;

// ─── CNF builder ──────────────────────────────────────────────────────────────

const _clauses: Clause[] = [];

function initCNF() { _clauses.length = 0; }

function addClause(...lits: Lit[]): void { _clauses.push(lits); }

function fix(v: number, val: boolean): void { addClause(val ? v : -v); }

// c = a XOR b  (4 Tseitin clauses)
function xorGate(a: Lit, b: Lit, c: Lit): void {
  addClause(-a, -b, -c);
  addClause(-a,  b,  c);
  addClause( a, -b,  c);
  addClause( a,  b, -c);
}

// c = a AND b  (3 Tseitin clauses)
function andGate(a: Lit, b: Lit, c: Lit): void {
  addClause(-a, -b, c);
  addClause( a, -c);
  addClause( b, -c);
}

// c = NOT a  (2 Tseitin clauses)
function notGate(a: Lit, c: Lit): void {
  addClause( a,  c);
  addClause(-a, -c);
}

// ─── Constant literals ────────────────────────────────────────────────────────

function constBit(val: number): Lit {
  return val ? TRUE_LIT() : FALSE_LIT();
}

function constWord(n: number): Word {
  return Array.from({ length: 32 }, (_, i) => constBit((n >>> i) & 1));
}

// ─── Word operations ──────────────────────────────────────────────────────────

// ROTR: pure permutation, no new variables
function rotr(w: Word, n: number): Word {
  return Array.from({ length: 32 }, (_, i) => w[(i + n) % 32]);
}

// SHR: permutation + fill with 0
function shr(w: Word, n: number): Word {
  return Array.from({ length: 32 }, (_, i) =>
    i + n < 32 ? w[i + n] : FALSE_LIT()
  );
}

// XOR two words, bitwise
function xorWords(a: Word, b: Word): Word {
  return Array.from({ length: 32 }, (_, i) => {
    const c = freshVar();
    xorGate(a[i], b[i], c);
    return c;
  });
}

// XOR three words
function xor3(a: Word, b: Word, c: Word): Word {
  return xorWords(xorWords(a, b), c);
}

// AND two words
function andWords(a: Word, b: Word): Word {
  return Array.from({ length: 32 }, (_, i) => {
    const c = freshVar();
    andGate(a[i], b[i], c);
    return c;
  });
}

// NOT a word
function notWord(a: Word): Word {
  return Array.from({ length: 32 }, (_, i) => {
    const c = freshVar();
    notGate(a[i], c);
    return c;
  });
}

// Full adder → [sum_bit, carry_out]
function fullAdder(a: Lit, b: Lit, cin: Lit): [Lit, Lit] {
  const t = freshVar();
  xorGate(a, b, t);
  const s = freshVar();
  xorGate(t, cin, s);

  const cout = freshVar();
  addClause(-a, -b, cout);
  addClause(-a, -cin, cout);
  addClause(-b, -cin, cout);
  addClause(a, b, -cout);
  addClause(a, cin, -cout);
  addClause(b, cin, -cout);

  return [s, cout];
}

// 32-bit ripple-carry adder (mod 2^32)
function add32(a: Word, b: Word): Word {
  const result: Lit[] = [];
  let carry = FALSE_LIT();
  for (let i = 0; i < 32; i++) {
    const [s, cout] = fullAdder(a[i], b[i], carry);
    result.push(s);
    carry = cout;
  }
  return result; // carry out discarded
}

function addMany(words: Word[]): Word {
  return words.reduce((acc, w) => add32(acc, w));
}

// ─── SHA-256 functions ────────────────────────────────────────────────────────

const bigSigma0 = (a: Word): Word => xor3(rotr(a, 2), rotr(a, 13), rotr(a, 22));
const bigSigma1 = (e: Word): Word => xor3(rotr(e, 6), rotr(e, 11), rotr(e, 25));
const smallSigma0 = (w: Word): Word => xor3(rotr(w, 7), rotr(w, 18), shr(w, 3));
const smallSigma1 = (w: Word): Word => xor3(rotr(w, 17), rotr(w, 19), shr(w, 10));

function ch(e: Word, f: Word, g: Word): Word {
  return xorWords(andWords(e, f), andWords(notWord(e), g));
}

function maj(a: Word, b: Word, c: Word): Word {
  return xor3(andWords(a, b), andWords(a, c), andWords(b, c));
}

// ─── SHA-256 constants ────────────────────────────────────────────────────────

const SHA256_K = [
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

const SHA256_IV = [
  0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
  0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

// ─── Encoder ─────────────────────────────────────────────────────────────────

export interface EncodedSHA256 {
  numVars: number;
  clauses: Clause[];
  inputLits: Lit[];    // free input bits (LSB-first per byte)
  hashLits: Word[];    // 8 output words (32 bits each, LSB-first)
}

/**
 * Encode reduced SHA-256 (numRounds rounds) for a message of inputBytes bytes
 * as a CNF satisfiability instance.
 *
 * The message is a single 512-bit block with standard SHA-256 padding.
 * Supports messages of 1–3 bytes (fits in one block with room for length).
 */
export function encodeSHA256(inputBytes: number, numRounds: number): EncodedSHA256 {
  initPool();
  initCNF();

  // Reserve var 1 as TRUE
  fix(_varCount + 1, true); // will be done by initPool
  // Actually initPool already called freshVar() for the TRUE variable.
  // The TRUE_LIT is already set. Just add its constraint:
  addClause(TRUE_LIT()); // TRUE_LIT is always true

  // ── Free input bits ─────────────────────────────────────────────────────────
  // inputLits[i*8 + b] = bit b of input byte i (b=0 is LSB of the byte)
  const inputLits: Lit[] = Array.from({ length: inputBytes * 8 }, () => freshVar());

  // ── Build message words W[0..15] ─────────────────────────────────────────
  // SHA-256 uses big-endian 32-bit words.
  // W[wi] = (byte[wi*4] << 24) | (byte[wi*4+1] << 16) | (byte[wi*4+2] << 8) | byte[wi*4+3]
  // In our LSB-first bit array:
  //   W[wi][bit k] = bit k of the 32-bit word
  //   byte j of W[wi] occupies bits [(3-j)*8 .. (3-j)*8+7] of W[wi]
  //   So W[wi][(3-j)*8 + b] = bit b of byte (wi*4+j)

  const W: Word[] = [];

  // Length field: goes in last 8 bytes of block (words 14-15)
  // For short messages: word 14 = 0, word 15 = inputBytes * 8
  const msgBitLen = inputBytes * 8;

  for (let wi = 0; wi < 16; wi++) {
    const word: Lit[] = new Array(32);
    for (let j = 0; j < 4; j++) {
      const byteIdx = wi * 4 + j;
      const bitBase = (3 - j) * 8; // bit position in word of byte j's LSB

      for (let b = 0; b < 8; b++) {
        const bitPos = bitBase + b;

        if (byteIdx < inputBytes) {
          // Free input bit
          word[bitPos] = inputLits[byteIdx * 8 + b];
        } else if (byteIdx === inputBytes) {
          // Padding byte 0x80: only MSB (b=7) is 1
          word[bitPos] = b === 7 ? TRUE_LIT() : FALSE_LIT();
        } else if (byteIdx >= 56) {
          // 64-bit big-endian length field (bytes 56–63)
          const lengthByteIdx = byteIdx - 56; // 0..7
          // High 4 bytes (56-59) = 0 for short messages
          const byteVal = lengthByteIdx >= 4
            ? (msgBitLen >>> ((7 - lengthByteIdx) * 8)) & 0xff
            : 0;
          word[bitPos] = constBit((byteVal >>> b) & 1);
        } else {
          // Zero padding
          word[bitPos] = FALSE_LIT();
        }
      }
    }
    W.push(word);
  }

  // ── Expand message schedule to 64 words ──────────────────────────────────
  for (let i = 16; i < 64; i++) {
    const s0 = smallSigma0(W[i - 15]);
    const s1 = smallSigma1(W[i - 2]);
    W.push(addMany([W[i - 16], s0, W[i - 7], s1]));
  }

  // ── Compression ──────────────────────────────────────────────────────────
  let [a, b, c, d, e, f, g, h] = SHA256_IV.map(constWord);

  for (let i = 0; i < numRounds; i++) {
    const S1 = bigSigma1(e);
    const CH = ch(e, f, g);
    const K = constWord(SHA256_K[i]);
    const T1 = addMany([h, S1, CH, K, W[i]]);

    const S0 = bigSigma0(a);
    const MAJ = maj(a, b, c);
    const T2 = add32(S0, MAJ);

    h = g; g = f; f = e;
    e = add32(d, T1);
    d = c; c = b; b = a;
    a = add32(T1, T2);
  }

  // ── Final hash = IV + state ───────────────────────────────────────────────
  const finalState = [a, b, c, d, e, f, g, h];
  const hashLits: Word[] = SHA256_IV.map((iv, i) =>
    add32(constWord(iv), finalState[i])
  );

  return {
    numVars: totalVars(),
    clauses: [..._clauses],
    inputLits,
    hashLits,
  };
}

/**
 * Add unit clauses constraining the hash output to match targetHash (64-char hex).
 */
export function constrainHash(enc: EncodedSHA256, targetHash: string): Clause[] {
  const extra: Clause[] = [];
  for (let wi = 0; wi < 8; wi++) {
    const word = parseInt(targetHash.slice(wi * 8, wi * 8 + 8), 16);
    for (let i = 0; i < 32; i++) {
      const bit = (word >>> i) & 1;
      extra.push(bit ? [enc.hashLits[wi][i]] : [-enc.hashLits[wi][i]]);
    }
  }
  return extra;
}

// ─── DPLL Solver ──────────────────────────────────────────────────────────────

export interface SolveStats {
  decisions: number;
  propagations: number;
  backtracks: number;
  nodesExplored: number;
}

export type SolveStatus = 'sat' | 'unsat' | 'running' | 'cancelled';

export interface SolveResult {
  status: SolveStatus;
  assignment: Map<number, boolean>;
  stats: SolveStats;
}

function evalClause(
  clause: Clause,
  assignment: Map<number, boolean>
): boolean | undefined {
  let hasUndef = false;
  for (const lit of clause) {
    const v = Math.abs(lit);
    const val = assignment.get(v);
    if (val === undefined) { hasUndef = true; continue; }
    if (val === (lit > 0)) return true;
  }
  return hasUndef ? undefined : false;
}

function unitPropagate(
  clauses: Clause[],
  assignment: Map<number, boolean>,
  stats: SolveStats,
  added: number[]
): boolean {
  let changed = true;
  while (changed) {
    changed = false;
    for (const clause of clauses) {
      const e = evalClause(clause, assignment);
      if (e === true) continue;
      if (e === false) return true; // conflict

      // Find unset literals
      let forceLit: Lit | null = null;
      let unsetCount = 0;
      for (const lit of clause) {
        if (!assignment.has(Math.abs(lit))) {
          unsetCount++;
          forceLit = lit;
        }
      }
      if (unsetCount === 0) return true; // conflict (all false)
      if (unsetCount === 1 && forceLit !== null) {
        const v = Math.abs(forceLit);
        if (!assignment.has(v)) {
          assignment.set(v, forceLit > 0);
          added.push(v);
          stats.propagations++;
          changed = true;
        }
      }
    }
  }
  return false;
}

async function dpllRecurse(
  clauses: Clause[],
  assignment: Map<number, boolean>,
  stats: SolveStats,
  onProgress: ((s: SolveStats, a: Map<number, boolean>) => void) | undefined,
  cancelRef: { cancelled: boolean },
  progressCounter: { n: number }
): Promise<boolean> {
  if (cancelRef.cancelled) return false;

  stats.nodesExplored++;
  progressCounter.n++;

  if (progressCounter.n % 200 === 0) {
    onProgress?.(stats, assignment);
    await new Promise(r => setTimeout(r, 0)); // yield to UI
    if (cancelRef.cancelled) return false;
  }

  const propagated: number[] = [];
  const conflict = unitPropagate(clauses, assignment, stats, propagated);

  if (conflict) {
    for (const v of propagated) assignment.delete(v);
    return false;
  }

  // Check all clauses
  let allSat = true;
  for (const clause of clauses) {
    const e = evalClause(clause, assignment);
    if (e === false) {
      for (const v of propagated) assignment.delete(v);
      return false;
    }
    if (e === undefined) allSat = false;
  }
  if (allSat) return true;

  // Choose first unassigned variable from first undetermined clause
  let chosen: number | null = null;
  outer: for (const clause of clauses) {
    if (evalClause(clause, assignment) !== undefined) continue;
    for (const lit of clause) {
      const v = Math.abs(lit);
      if (!assignment.has(v)) { chosen = v; break outer; }
    }
  }

  if (chosen === null) {
    for (const v of propagated) assignment.delete(v);
    return false;
  }

  stats.decisions++;

  // Branch true
  assignment.set(chosen, true);
  if (await dpllRecurse(clauses, assignment, stats, onProgress, cancelRef, progressCounter)) return true;
  assignment.delete(chosen);

  stats.backtracks++;

  // Branch false
  assignment.set(chosen, false);
  if (await dpllRecurse(clauses, assignment, stats, onProgress, cancelRef, progressCounter)) return true;
  assignment.delete(chosen);

  for (const v of propagated) assignment.delete(v);
  return false;
}

export async function solveDPLL(
  allClauses: Clause[],
  onProgress?: (stats: SolveStats, assignment: Map<number, boolean>) => void,
  cancelRef?: { cancelled: boolean }
): Promise<SolveResult> {
  const stats: SolveStats = { decisions: 0, propagations: 0, backtracks: 0, nodesExplored: 0 };
  const assignment = new Map<number, boolean>();
  const cancel = cancelRef ?? { cancelled: false };
  const counter = { n: 0 };

  const sat = await dpllRecurse(allClauses, assignment, stats, onProgress, cancel, counter);

  return {
    status: cancel.cancelled ? 'cancelled' : sat ? 'sat' : 'unsat',
    assignment,
    stats,
  };
}

// ─── Extract input bytes from a satisfying assignment ────────────────────────

export function extractInput(
  assignment: Map<number, boolean>,
  inputLits: Lit[]
): Uint8Array {
  const inputBytes = inputLits.length / 8;
  const result = new Uint8Array(inputBytes);
  for (let i = 0; i < inputBytes; i++) {
    let byte = 0;
    for (let b = 0; b < 8; b++) {
      const lit = inputLits[i * 8 + b];
      const v = Math.abs(lit);
      const val = assignment.get(v);
      const bitVal = val !== undefined ? (val === (lit > 0) ? 1 : 0) : 0;
      byte |= bitVal << b;
    }
    result[i] = byte;
  }
  return result;
}
