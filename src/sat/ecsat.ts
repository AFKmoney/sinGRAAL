// ─── secp256k1 Runtime Arithmetic (BigInt) ───────────────────────────────────

const SECP_P  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2Fn;
const SECP_GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n;
const SECP_GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8n;

function modpow(base: bigint, exp: bigint, mod: bigint): bigint {
  let result = 1n;
  base = ((base % mod) + mod) % mod;
  while (exp > 0n) {
    if (exp & 1n) result = result * base % mod;
    base = base * base % mod;
    exp >>= 1n;
  }
  return result;
}

interface Pt { x: bigint; y: bigint; inf: boolean; }

function ptAdd(A: Pt, B: Pt, p: bigint): Pt {
  if (A.inf) return B;
  if (B.inf) return A;
  if (A.x === B.x) {
    if (A.y !== B.y) return { x: 0n, y: 0n, inf: true };
    const lam = 3n * A.x * A.x % p * modpow(2n * A.y % p, p - 2n, p) % p;
    const x3  = (lam * lam % p - 2n * A.x % p + 2n * p) % p;
    const y3  = (lam * ((A.x - x3 + p) % p) % p - A.y + p) % p;
    return { x: x3, y: y3, inf: false };
  }
  const lam = ((B.y - A.y + p) % p) * modpow((B.x - A.x + p) % p, p - 2n, p) % p;
  const x3  = (lam * lam % p - A.x - B.x + 2n * p) % p;
  const y3  = (lam * ((A.x - x3 + p) % p) % p - A.y + p) % p;
  return { x: x3, y: y3, inf: false };
}

function scalarMul(G: Pt, k: bigint, p: bigint): Pt {
  let result: Pt = { x: 0n, y: 0n, inf: true };
  let addend = G;
  while (k > 0n) {
    if (k & 1n) result = ptAdd(result, addend, p);
    addend = ptAdd(addend, addend, p);
    k >>= 1n;
  }
  return result;
}

// Find all points on y² = x³ + 7 mod p (small p only)
function findCurvePoints(p: number): Array<{x: number; y: number}> {
  const bp = BigInt(p);
  const pts: Array<{x: number; y: number}> = [];
  for (let x = 0; x < p; x++) {
    const rhs = (BigInt(x) ** 3n + 7n) % bp;
    for (let y = 0; y < p; y++) {
      if (BigInt(y) * BigInt(y) % bp === rhs) pts.push({ x, y });
    }
  }
  return pts;
}

// Compute group order via repeated addition (for small p)
function groupOrder(G: Pt, p: bigint, bound: number): number {
  let Q = G;
  for (let n = 1; n <= bound; n++) {
    if (Q.inf) return n;
    Q = ptAdd(Q, G, p);
  }
  return bound;
}

// Public API: generate a secp256k1 mini-puzzle
export function generateMiniPuzzle(keyBits: number): { k: bigint; Px: bigint; Py: bigint } {
  const G: Pt = { x: SECP_GX, y: SECP_GY, inf: false };
  const min = 1n << BigInt(keyBits - 1);
  const max = (1n << BigInt(keyBits)) - 1n;
  const range = max - min + 1n;
  const rand = BigInt(Math.floor(Math.random() * Number(range < 2n ** 32n ? range : 2n ** 32n)));
  const k = min + rand % range;
  const Q = scalarMul(G, k, SECP_P);
  return { k, Px: Q.x, Py: Q.y };
}

export function verifyKey(k: bigint, Px: bigint, Py: bigint): boolean {
  const G: Pt = { x: SECP_GX, y: SECP_GY, inf: false };
  const Q = scalarMul(G, k, SECP_P);
  return !Q.inf && Q.x === Px && Q.y === Py;
}

// Clause count estimate for secp256k1 ECDLP (n unknown bits)
export function estimateSecp256k1Clauses(unknownBits: number) {
  const fieldBits = 256;
  const clausesPerFieldMul = fieldBits * fieldBits * 3 + fieldBits * 14 * 6 + fieldBits * 15;
  const mulsForInverse = 2 * fieldBits;
  const fieldMulsPerAdd = mulsForInverse + 5;
  const clausesPerAdd = fieldMulsPerAdd * clausesPerFieldMul;
  const totalClauses = unknownBits * clausesPerAdd;
  const varsPerFieldMul = fieldBits * fieldBits + fieldBits * 4;
  const totalVars = unknownBits * fieldMulsPerAdd * varsPerFieldMul;
  return { fieldBits, fieldMulsPerAdd, clausesPerAdd, totalClauses, totalVars };
}

// ─── CNF State ────────────────────────────────────────────────────────────────

export type Lit = number;
type FieldLits = Lit[];  // LSB-first bit array

let _vc = 0;
const _cls: number[][] = [];

export function ecInitPool() { _vc = 1; _cls.length = 0; _cls.push([1]); } // var 1 = TRUE
export function ecGetClauses(): number[][] { return _cls; }
export function ecGetVarCount(): number { return _vc; }

function fv(): number { return ++_vc; }
function cl(...lits: number[]): void { _cls.push(lits); }

const TRUE_LIT  = (): Lit =>  1;
const FALSE_LIT = (): Lit => -1;

// ─── Tseitin Gates ───────────────────────────────────────────────────────────

function andG(a: Lit, b: Lit): Lit {
  const c = fv();
  cl(-a, -b, c); cl(a, -c); cl(b, -c);
  return c;
}

function xorG(a: Lit, b: Lit): Lit {
  const c = fv();
  cl(-a, -b, -c); cl(-a, b, c); cl(a, -b, c); cl(a, b, -c);
  return c;
}

// MUX: out = sel ? a : b
function muxG(sel: Lit, a: Lit, b: Lit): Lit {
  if (sel === TRUE_LIT())  return a;
  if (sel === FALSE_LIT()) return b;
  if (a === b) return a;
  const c = fv();
  cl( sel, -b, c); cl( sel, b, -c);
  cl(-sel, -a, c); cl(-sel, a, -c);
  return c;
}

function constLit(v: boolean): Lit  { return v ? TRUE_LIT() : FALSE_LIT(); }
function constVec(v: bigint, n: number): FieldLits {
  return Array.from({ length: n }, (_, i) => constLit(Boolean((v >> BigInt(i)) & 1n)));
}
function newVec(n: number): FieldLits { return Array.from({ length: n }, fv); }

function forceEq(a: FieldLits, b: FieldLits): void {
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) { cl(a[i], -b[i]); cl(-a[i], b[i]); }
}

function forceConst(a: FieldLits, v: bigint): void {
  for (let i = 0; i < a.length; i++) cl(((v >> BigInt(i)) & 1n) ? a[i] : -a[i]);
}

// ─── Ripple-Carry Adder ───────────────────────────────────────────────────────

// Returns sum bits + carry as last element
function addVecs(a: FieldLits, b: FieldLits): FieldLits {
  const n = Math.max(a.length, b.length);
  const sum: FieldLits = [];
  let carry: Lit = FALSE_LIT();
  for (let i = 0; i < n; i++) {
    const ai = i < a.length ? a[i] : FALSE_LIT();
    const bi = i < b.length ? b[i] : FALSE_LIT();
    const ab = xorG(ai, bi);
    sum.push(xorG(ab, carry));
    // carry = MAJ(ai, bi, carry) via OR of three ANDs
    const c0 = andG(ai, bi);
    const c1 = andG(ai, carry);
    const c2 = andG(bi, carry);
    // out = c0 OR c1 OR c2
    const t = fv();
    cl(-c0, t); cl(-c1, t); cl(-c2, t);
    cl(c0, c1, c2, -t);
    carry = t;
  }
  sum.push(carry);
  return sum;
}

// Ripple-borrow subtractor: returns diff bits + noBorrow as last element
// noBorrow = 1 when a >= b (no borrow out)
function subVecs(a: FieldLits, b: FieldLits): FieldLits {
  const n = Math.max(a.length, b.length);
  const diff: FieldLits = [];
  let borrow: Lit = FALSE_LIT();
  for (let i = 0; i < n; i++) {
    const ai = i < a.length ? a[i] : FALSE_LIT();
    const bi = i < b.length ? b[i] : FALSE_LIT();
    const ab = xorG(ai, bi);
    diff.push(xorG(ab, borrow));
    // borrow_out = (NOT ai AND bi) OR (NOT ai AND borrow) OR (bi AND borrow)
    const notAi = -ai;
    const c0 = andG(notAi, bi);
    const c1 = andG(notAi, borrow);
    const c2 = andG(bi, borrow);
    const t = fv();
    cl(-c0, t); cl(-c1, t); cl(-c2, t);
    cl(c0, c1, c2, -t);
    borrow = t;
  }
  diff.push(-borrow); // last element = noBorrow (1 when a >= b)
  return diff;
}

// ─── Field Arithmetic (Mersenne prime p = 2^m - 1) ───────────────────────────
// n = m = bit length of p.  All field elements use n bits.

// Add mod p (Mersenne): result = a + b mod p
// If a + b < p: return a+b. If a + b >= p: return a + b - p.
// For Mersenne p=2^m-1: if carry out, wrap around (add 1).
function fieldAdd(a: FieldLits, b: FieldLits, p: bigint, n: number): FieldLits {
  const rawSum = addVecs(a, b); // n+1 bits, last bit = carry
  const carry = rawSum[n];
  const sumBits = rawSum.slice(0, n);
  // If carry=1, result = sum - p = sum - (2^n - 1) = sum + 1 mod 2^n
  // (since p = 2^n - 1, subtracting p when there's carry is same as wrapping)
  // Simple version: result = carry ? (sumBits + 1) mod 2^n : sumBits
  // But: sumBits + 1 could itself overflow (e.g. if sumBits = 11...1)
  // For Mersenne: a,b < p = 2^n-1, so a+b < 2*(2^n-1) = 2^(n+1) - 2.
  // If carry=1: a+b >= 2^n, so rawSum-p = a+b-2^n+1 = sumBits+1 < 2^n.
  // So sumBits+1 won't overflow.
  // Compute sumPlusOne = sumBits + 1:
  const one: FieldLits = Array.from({ length: n }, (_, i) => i === 0 ? TRUE_LIT() : FALSE_LIT());
  const spo = addVecs(sumBits, one).slice(0, n);
  return Array.from({ length: n }, (_, i) => muxG(carry, spo[i], sumBits[i]));
}

// Negate mod p: neg(a) = p - a = (2^n - 1) - a = NOT a (for Mersenne!)
function fieldNeg(a: FieldLits): FieldLits {
  return a.map(l => -l); // bitwise NOT = p - a for Mersenne p = 2^n - 1
}

// Subtract mod p: a - b = a + (p - b) = a + NOT(b) (for Mersenne)
function fieldSub(a: FieldLits, b: FieldLits, p: bigint, n: number): FieldLits {
  return fieldAdd(a, fieldNeg(b), p, n);
}

// Multiply: a * b mod p (Mersenne, schoolbook + Mersenne reduction)
// For p = 2^n - 1: (x mod 2^n) + (x >> n) gives reduction (possibly one more step)
function fieldMul(a: FieldLits, b: FieldLits, p: bigint, n: number): FieldLits {
  // Step 1: schoolbook multiplication to get 2n-bit product
  // product[k] = XOR over all i+j=k of AND(a[i], b[j])... no, addition not XOR.
  // Use shift-and-add: for each bit b[i], conditionally add a<<i to accumulator.
  let acc: FieldLits = Array(2 * n).fill(FALSE_LIT());
  for (let i = 0; i < n; i++) {
    // Contribution: b[i] ? (a << i) : 0
    const shifted: FieldLits = Array.from({ length: 2 * n }, (_, k) => {
      if (k < i || k >= i + n) return FALSE_LIT();
      return andG(b[i], a[k - i]);
    });
    acc = addVecs(acc, shifted).slice(0, 2 * n);
  }
  // Step 2: Mersenne reduction: p = 2^n - 1
  // a * b mod (2^n - 1) = (lo + hi) mod (2^n - 1)
  // where lo = acc[0..n-1], hi = acc[n..2n-1]
  // If lo + hi < 2^n: result = lo + hi (if < p, done; if = 2^n-1 = p, reduce to 0)
  // If lo + hi >= 2^n: result = lo + hi - p = lo + hi - 2^n + 1
  const lo = acc.slice(0, n);
  const hi = acc.slice(n, 2 * n);
  return fieldAdd(lo, hi, p, n);
}

// Modular exponentiation (base is CNF, exp is constant bigint)
function fieldPow(base: FieldLits, exp: bigint, p: bigint, n: number): FieldLits {
  let result: FieldLits = constVec(1n, n);
  let cur = base;
  let e = exp;
  while (e > 0n) {
    if (e & 1n) result = fieldMul(result, cur, p, n);
    cur = fieldMul(cur, cur, p, n);
    e >>= 1n;
  }
  return result;
}

// Modular inverse via Fermat: a^(p-2) mod p
function fieldInv(a: FieldLits, p: bigint, n: number): FieldLits {
  return fieldPow(a, p - 2n, p, n);
}

// ─── EC Affine Point Addition ─────────────────────────────────────────────────

interface ECPt { x: FieldLits; y: FieldLits; }

// Add S + (Cx, Cy) where Cx, Cy are runtime constants (precomputed G_i values)
function ecAddConst(S: ECPt, Cx: bigint, Cy: bigint, p: bigint, n: number): ECPt {
  const Cy_lit = constVec(Cy, n);
  const Cx_lit = constVec(Cx, n);
  const dy = fieldSub(Cy_lit, S.y, p, n);
  const dx = fieldSub(Cx_lit, S.x, p, n);
  const dxInv = fieldInv(dx, p, n);
  const lam = fieldMul(dy, dxInv, p, n);
  const lamSq = fieldMul(lam, lam, p, n);
  const x3 = fieldSub(fieldSub(lamSq, Cx_lit, p, n), S.x, p, n);
  const sx3 = fieldSub(S.x, x3, p, n);
  const y3 = fieldSub(fieldMul(lam, sx3, p, n), S.y, p, n);
  return { x: x3, y: y3 };
}

// Conditional add: if bit=1, return S + C; else return S
function condEcAdd(bit: Lit, S: ECPt, Cx: bigint, Cy: bigint, p: bigint, n: number): ECPt {
  const added = ecAddConst(S, Cx, Cy, p, n);
  return {
    x: added.x.map((l, i) => muxG(bit, l, S.x[i])),
    y: added.y.map((l, i) => muxG(bit, l, S.y[i])),
  };
}

// ─── ECDLP Puzzle Types ───────────────────────────────────────────────────────

export interface TinyCurveParams {
  p: number;
  Gx: number;
  Gy: number;
  order: number;
  Px: number;
  Py: number;
  k: number;
  unknownBits: number;
}

export interface ECDLPEncoding {
  clauses: number[][];
  numVars: number;
  kBits: Lit[];
  fieldBits: number;
  curveP: bigint;
}

// Only Mersenne primes for the demo (ensures correct reduction)
export const DEMO_PRIMES = [31, 127];

// Build a tiny ECDLP puzzle on y² = x³ + 7 mod p (Mersenne p)
export function buildTinyPuzzle(p: number, unknownBits: number): TinyCurveParams | null {
  const bp = BigInt(p);
  const pts = findCurvePoints(p);
  if (pts.length < 2) return null;
  // Find a point with odd-prime order (avoid 2-torsion)
  let gen = pts.find(pt => pt.y !== 0);
  if (!gen) return null;
  const G0pt: Pt = { x: BigInt(gen.x), y: BigInt(gen.y), inf: false };
  const ord = groupOrder(G0pt, bp, p + 10);
  const min = 1 << (unknownBits - 1);
  const max = (1 << unknownBits) - 1;
  const k = min + Math.floor(Math.random() * (max - min + 1));
  const Qpt = scalarMul(G0pt, BigInt(k), bp);
  if (Qpt.inf) return null;
  return {
    p, Gx: gen.x, Gy: gen.y, order: ord,
    Px: Number(Qpt.x), Py: Number(Qpt.y),
    k, unknownBits,
  };
}

// Encode ECDLP as CNF: find k (unknownBits decision vars) such that k×G = P
export function encodeECDLP(params: TinyCurveParams): ECDLPEncoding {
  const { p, Gx, Gy, Px, Py, unknownBits } = params;
  const bp = BigInt(p);

  // n = bit length of p (must be Mersenne: p = 2^n - 1)
  let n = 0;
  let tmp = bp;
  while (tmp > 0n) { n++; tmp >>= 1n; }

  ecInitPool();

  // Decision variables: k_bits[i] = bit i of k (LSB first)
  const kBits: Lit[] = Array.from({ length: unknownBits }, fv);
  // MSB is always 1 (k in [2^(unknownBits-1), 2^unknownBits))
  cl(kBits[unknownBits - 1]);

  // Precompute G_i = 2^i × G for i = 0..unknownBits-1 (runtime constants)
  const G: Pt = { x: BigInt(Gx), y: BigInt(Gy), inf: false };
  const Gi: Array<{x: bigint; y: bigint}> = [];
  let cur = G;
  for (let i = 0; i < unknownBits; i++) {
    Gi.push({ x: cur.x, y: cur.y });
    cur = ptAdd(cur, cur, bp);
  }

  // Accumulator starts at G_{n-1} (MSB is fixed to 1, so this point is always added)
  let acc: ECPt = {
    x: constVec(Gi[unknownBits - 1].x, n),
    y: constVec(Gi[unknownBits - 1].y, n),
  };

  // Conditionally add G_i for i = 0..unknownBits-2 based on kBits[i]
  for (let i = 0; i < unknownBits - 1; i++) {
    acc = condEcAdd(kBits[i], acc, Gi[i].x, Gi[i].y, bp, n);
  }

  // Constrain accumulator to target public key P
  forceConst(acc.x, BigInt(Px));
  forceConst(acc.y, BigInt(Py));

  return {
    clauses: [..._cls],
    numVars: _vc,
    kBits,
    fieldBits: n,
    curveP: bp,
  };
}

// Extract k from CDCL assignment
export function extractKey(assignment: Map<number, boolean>, kBits: Lit[]): bigint {
  let k = 0n;
  for (let i = 0; i < kBits.length; i++) {
    const lit = kBits[i];
    const v = Math.abs(lit);
    const val = assignment.get(v);
    const bit = val !== undefined ? (val === (lit > 0) ? 1n : 0n) : 0n;
    k |= bit << BigInt(i);
  }
  return k;
}
