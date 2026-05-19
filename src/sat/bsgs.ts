// ─── Baby-Step Giant-Step (BSGS) for secp256k1 ECDLP ────────────────────────
//
// Solves: find k in [low, low + 2^babyBits × 2^giantBits) such that k×G = P
//
// Strategy (2×n split):
//   k = g × 2^babyBits + b
//   k×G = P  ⟺  g×(2^babyBits × G) = P - b×G
//
//   Baby table: {x-coord of b×G → b} for b in [0, 2^babyBits)
//   Giant loop: for g in [0, 2^giantBits): check if P - g×Gstep is in table
//
// Complexity: O(2^babyBits + 2^giantBits) time, O(2^babyBits) memory

const P  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2Fn;
const GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n;
const GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8n;

export interface BsgsPt { x: bigint; y: bigint; inf: boolean; }

function mpow(base: bigint, exp: bigint, mod: bigint): bigint {
  let r = 1n; base = ((base % mod) + mod) % mod;
  while (exp > 0n) { if (exp & 1n) r = r * base % mod; base = base * base % mod; exp >>= 1n; }
  return r;
}

export function ptAdd(A: BsgsPt, B: BsgsPt): BsgsPt {
  if (A.inf) return B; if (B.inf) return A;
  if (A.x === B.x) {
    if (A.y !== B.y) return { x:0n,y:0n,inf:true };
    const lam = 3n*A.x*A.x%P * mpow(2n*A.y%P, P-2n, P) % P;
    const x3  = (lam*lam%P - 2n*A.x%P + 2n*P) % P;
    return { x:x3, y:(lam*((A.x-x3+P)%P)%P - A.y+P)%P, inf:false };
  }
  const lam = ((B.y-A.y+P)%P) * mpow((B.x-A.x+P)%P, P-2n, P) % P;
  const x3  = (lam*lam%P - A.x - B.x + 2n*P) % P;
  return { x:x3, y:(lam*((A.x-x3+P)%P)%P - A.y+P)%P, inf:false };
}

export function ptNeg(A: BsgsPt): BsgsPt {
  return A.inf ? A : { x:A.x, y:(P-A.y)%P, inf:false };
}

export function scalarMul(G: BsgsPt, k: bigint): BsgsPt {
  let r: BsgsPt = {x:0n,y:0n,inf:true}, a = G;
  while (k>0n) { if(k&1n) r=ptAdd(r,a); a=ptAdd(a,a); k>>=1n; }
  return r;
}

export const SECP_G: BsgsPt = { x:GX, y:GY, inf:false };

export interface BsgsStats {
  babySteps: number;
  giantSteps: number;
  tableSize: number;
  found: boolean;
  k: bigint | null;
  elapsedMs: number;
}

// Core BSGS engine. Yields control every `yieldEvery` steps for UI updates.
export async function bsgsSolve(opts: {
  G: BsgsPt;             // generator
  P: BsgsPt;             // target public key
  low: bigint;           // lower bound for k
  babyBits: number;      // baby step count = 2^babyBits
  giantBits: number;     // giant step count = 2^giantBits
  yieldEvery?: number;
  onProgress?: (phase: 'building'|'searching', step: number, total: number) => void;
  cancelRef?: { cancelled: boolean };
}): Promise<BsgsStats> {
  const { G, P: target, low, babyBits, giantBits, yieldEvery = 50_000 } = opts;
  const nBaby  = 1 << babyBits;   // 2^babyBits
  const nGiant = 1 << giantBits;  // 2^giantBits
  const t0 = performance.now();

  // ── Baby steps: build table {x-coord → baby index} ──────────────────────
  const table = new Map<bigint, number>();
  let cur: BsgsPt = { x:0n,y:0n,inf:true }; // start at infinity (b=0 → identity)
  const G_baby = G; // each step adds G once

  for (let b = 0; b < nBaby; b++) {
    // cur = b × G
    if (!cur.inf) table.set(cur.x, b);
    cur = ptAdd(cur, G_baby);
    if (b % yieldEvery === 0) {
      opts.onProgress?.('building', b, nBaby);
      await new Promise<void>(r => requestAnimationFrame(() => r()));
      if (opts.cancelRef?.cancelled) return { babySteps:b, giantSteps:0, tableSize:table.size, found:false, k:null, elapsedMs:performance.now()-t0 };
    }
  }

  // Giant step generator: Gstep = 2^babyBits × G
  const Gstep = scalarMul(G, BigInt(nBaby));
  const negGstep = ptNeg(Gstep);

  // Adjust target for the lower bound: P' = target - low × G
  let searchPt = ptAdd(target, ptNeg(scalarMul(G, low)));

  // ── Giant steps: for g in [0, nGiant): check if searchPt is in table ────
  for (let g = 0; g < nGiant; g++) {
    if (!searchPt.inf && table.has(searchPt.x)) {
      const b = table.get(searchPt.x)!;
      // Recover y sign: baby = b×G, giant shifts by g×nBaby
      const babyPt = scalarMul(G, BigInt(b));
      if (!babyPt.inf && babyPt.y === searchPt.y) {
        const k = low + BigInt(g) * BigInt(nBaby) + BigInt(b);
        return { babySteps:nBaby, giantSteps:g+1, tableSize:table.size, found:true, k, elapsedMs:performance.now()-t0 };
      }
      // Check negation (k might use the negative baby)
      if (!babyPt.inf && (P - babyPt.y) === searchPt.y) {
        const k = low + BigInt(g) * BigInt(nBaby) - BigInt(b);
        if (k >= low) return { babySteps:nBaby, giantSteps:g+1, tableSize:table.size, found:true, k, elapsedMs:performance.now()-t0 };
      }
    }
    searchPt = ptAdd(searchPt, negGstep);
    if (g % yieldEvery === 0) {
      opts.onProgress?.('searching', g, nGiant);
      await new Promise<void>(r => requestAnimationFrame(() => r()));
      if (opts.cancelRef?.cancelled) return { babySteps:nBaby, giantSteps:g, tableSize:table.size, found:false, k:null, elapsedMs:performance.now()-t0 };
    }
  }

  return { babySteps:nBaby, giantSteps:nGiant, tableSize:table.size, found:false, k:null, elapsedMs:performance.now()-t0 };
}

// Compute how many bits can BSGS cover given memory/time budget
export function bsgsCapacity(babyBits: number, giantBits: number) {
  return {
    totalBits: babyBits + giantBits,
    memoryMB: Math.pow(2, babyBits) * 64 / 1_000_000,
    babyOps:  Math.pow(2, babyBits),
    giantOps: Math.pow(2, giantBits),
  };
}

// Estimate projection: for secp256k1 puzzle N, how much SAT needs to fix
export function hybridProjection(puzzleBits: number, babyBits: number, giantBits: number) {
  const bsgsBits = babyBits + giantBits;
  const satMustFix = Math.max(0, puzzleBits - bsgsBits);
  const kangarooOps = Math.pow(2, puzzleBits / 2);
  const hybridOps   = Math.pow(2, babyBits) + Math.pow(2, giantBits);
  const speedup     = kangarooOps / hybridOps;
  return { satMustFix, bsgsBits, kangarooOps, hybridOps, speedup };
}
