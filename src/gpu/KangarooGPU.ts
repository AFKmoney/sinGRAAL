// WebGPU Kangaroo Engine
// 6-automorphism distinguished-point Kangaroo on secp256k1
// Runs thousands of animals in parallel on the GPU.

import shaderSrc from './kangaroo.wgsl?raw';

// ─── secp256k1 constants ──────────────────────────────────────────────────────

const P  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2Fn;
const GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n;
const GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8n;
const N  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141n;
const BETA  = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EEn;
const BETA2 = 0x851695D49A83F8EF919BB86153CBCB16630FB68AED0A766A3EC693D68E6AFA40n;

// ─── BigInt ↔ buffer helpers ──────────────────────────────────────────────────

function bigintToU256LE(n: bigint): Uint32Array {
  const r = new Uint32Array(8);
  let v = ((n % P) + P) % P;
  for (let i = 0; i < 8; i++) {
    r[i] = Number(v & 0xFFFFFFFFn);
    v >>= 32n;
  }
  return r;
}

function u256LEtoBigint(arr: Uint32Array, offset = 0): bigint {
  let r = 0n;
  for (let i = 7; i >= 0; i--) r = (r << 32n) | BigInt(arr[offset + i]);
  return r;
}

// ─── Point arithmetic (JS, for setup) ────────────────────────────────────────

function mpow(b: bigint, e: bigint, m: bigint): bigint {
  let r = 1n; b = ((b % m) + m) % m;
  while (e > 0n) { if (e & 1n) r = r * b % m; b = b * b % m; e >>= 1n; }
  return r;
}

interface Pt { x: bigint; y: bigint; inf: boolean }

function ptAdd(A: Pt, B: Pt): Pt {
  if (A.inf) return B; if (B.inf) return A;
  if (A.x === B.x) {
    if (A.y !== B.y) return { x:0n, y:0n, inf:true };
    const lam = 3n*A.x*A.x%P * mpow(2n*A.y, P-2n, P) % P;
    const x3  = (lam*lam - 2n*A.x%P + 2n*P) % P;
    return { x:x3, y:(lam*(A.x-x3+P)%P - A.y+P)%P, inf:false };
  }
  const lam = ((B.y-A.y+P)%P) * mpow((B.x-A.x+P), P-2n, P) % P;
  const x3  = (lam*lam - A.x - B.x + 2n*P) % P;
  return { x:x3, y:(lam*(A.x-x3+P)%P - A.y+P)%P, inf:false };
}

function scMul(k: bigint): Pt {
  let r: Pt = {x:0n,y:0n,inf:true};
  let a: Pt = {x:GX, y:GY, inf:false};
  while (k > 0n) { if (k&1n) r = ptAdd(r,a); a = ptAdd(a,a); k >>= 1n; }
  return r;
}

function canonX(x: bigint): bigint {
  const x1 = BETA  * x % P;
  const x2 = BETA2 * x % P;
  return [x, x1, x2].reduce((m, v) => v < m ? v : m);
}

// ─── Buffer layout ────────────────────────────────────────────────────────────

// Animal: 3 Points (each 8×u32 x3 for x,y,z) + scalar(8) + is_wild(1) + pad(1) = 3*24 + 9 u32 = 81 u32
// Simplified: store x(8), y(8), z(8), scalar(8), is_wild(1), pad(3) = 36 u32 = 144 bytes
const ANIMAL_U32S = 36;

// JumpPoint: x(8) + y(8) + scalar(8) = 24 u32 = 96 bytes
const JUMP_U32S = 24;

// DP entry: canon_x(8) + scalar(8) + is_wild(1) + pad(3) = 20 u32 = 80 bytes
const DP_U32S = 20;

// Params: 4 u32 = 16 bytes
const PARAMS_U32S = 4;

// ─── Build jump table ─────────────────────────────────────────────────────────

function buildJumps(rangeBits: number, numJumps: number): { pts: Pt[]; scalars: bigint[] } {
  const base = Math.max(0, Math.floor(rangeBits / 2) - 1);
  const pts: Pt[] = [];
  const scalars: bigint[] = [];
  for (let i = 0; i < numJumps; i++) {
    const bit = base + Math.floor(i * rangeBits / (numJumps * 2));
    const sc = 1n << BigInt(bit);
    pts.push(scMul(sc));
    scalars.push(sc % N);
  }
  return { pts, scalars };
}

// ─── GPU Engine ───────────────────────────────────────────────────────────────

export interface GPUKangarooOptions {
  targetX: bigint;
  targetY: bigint;
  rangeLow: bigint;
  rangeBits: number;
  numAnimals?: number;   // default 4096
  numJumps?: number;     // default 32
  dpBits?: number;       // default 28
  stepsPerDispatch?: number; // default 256
}

export interface GPUKangarooResult {
  steps: number;
  dps: number;
  solution: bigint | null;
  gpuAvailable: boolean;
}

interface DPEntry { canonX: bigint; scalar: bigint; isWild: boolean }

export class KangarooGPU {
  private device: GPUDevice | null = null;
  private pipeline: GPUComputePipeline | null = null;

  private animalBuf!: GPUBuffer;
  private jumpBuf!: GPUBuffer;
  private dpBuf!: GPUBuffer;
  private dpCountBuf!: GPUBuffer;
  private dpReadBuf!: GPUBuffer;
  private dpCountReadBuf!: GPUBuffer;
  private paramsBuf!: GPUBuffer;
  private bindGroup!: GPUBindGroup;

  private opts!: Required<GPUKangarooOptions>;
  private dpTable = new Map<string, DPEntry>();

  public steps = 0;
  public dps = 0;
  public solution: bigint | null = null;
  public gpuAvailable = false;

  async init(opts: GPUKangarooOptions): Promise<boolean> {
    this.opts = {
      numAnimals: 4096,
      numJumps: 32,
      dpBits: 20, // start conservative for demo (28 for real)
      stepsPerDispatch: 64,
      ...opts,
    };

    // Try to get WebGPU device
    if (!navigator.gpu) return false;
    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) return false;
    this.device = await adapter.requestDevice();
    this.gpuAvailable = true;

    await this.buildPipeline();
    this.initBuffers();
    return true;
  }

  private async buildPipeline() {
    const dev = this.device!;
    const module = dev.createShaderModule({ code: shaderSrc });

    this.pipeline = await dev.createComputePipelineAsync({
      layout: 'auto',
      compute: { module, entryPoint: 'main' },
    });
  }

  private initBuffers() {
    const dev = this.device!;
    const { numAnimals, numJumps, rangeLow, rangeBits, targetX, targetY } = this.opts;
    const MAX_DPS = 8192;

    // ── Animal buffer ──
    const animalData = new Uint32Array(numAnimals * ANIMAL_U32S);
    const G: Pt = { x: GX, y: GY, inf: false };

    // Tame kangaroos: spread across range
    // Wild kangaroo: starts at target (index 0)
    for (let i = 0; i < numAnimals; i++) {
      const isWild = i === 0;
      const base = i * ANIMAL_U32S;

      let pt: Pt;
      let sc: bigint;

      if (isWild) {
        pt = { x: targetX, y: targetY, inf: false };
        sc = 0n;
      } else {
        // Tame: start at rangeLow + i * (2^rangeBits / numAnimals)
        const stepSize = 1n << BigInt(Math.max(0, rangeBits - 12));
        sc = (rangeLow + BigInt(i) * stepSize) % N;
        pt = scMul(sc);
      }

      // x, y in affine, z = 1 (Jacobian representation of affine point)
      const px = bigintToU256LE(pt.x);
      const py = bigintToU256LE(pt.y);
      const pz = bigintToU256LE(1n);
      const ps = bigintToU256LE(sc);

      for (let j = 0; j < 8; j++) animalData[base + j]      = px[j];
      for (let j = 0; j < 8; j++) animalData[base + 8 + j]  = py[j];
      for (let j = 0; j < 8; j++) animalData[base + 16 + j] = pz[j];
      for (let j = 0; j < 8; j++) animalData[base + 24 + j] = ps[j];
      animalData[base + 32] = isWild ? 1 : 0;
    }

    this.animalBuf = dev.createBuffer({
      size: animalData.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
      mappedAtCreation: true,
    });
    new Uint32Array(this.animalBuf.getMappedRange()).set(animalData);
    this.animalBuf.unmap();

    // ── Jump table ──
    const { pts: jumpPts, scalars: jumpScalars } = buildJumps(rangeBits, numJumps);
    const jumpData = new Uint32Array(numJumps * JUMP_U32S);
    for (let i = 0; i < numJumps; i++) {
      const base = i * JUMP_U32S;
      const jx = bigintToU256LE(jumpPts[i].x);
      const jy = bigintToU256LE(jumpPts[i].y);
      const js = bigintToU256LE(jumpScalars[i]);
      for (let j = 0; j < 8; j++) jumpData[base + j]      = jx[j];
      for (let j = 0; j < 8; j++) jumpData[base + 8 + j]  = jy[j];
      for (let j = 0; j < 8; j++) jumpData[base + 16 + j] = js[j];
    }

    this.jumpBuf = dev.createBuffer({
      size: jumpData.byteLength,
      usage: GPUBufferUsage.STORAGE,
      mappedAtCreation: true,
    });
    new Uint32Array(this.jumpBuf.getMappedRange()).set(jumpData);
    this.jumpBuf.unmap();

    // ── DP buffer ──
    this.dpBuf = dev.createBuffer({
      size: MAX_DPS * DP_U32S * 4,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
    });

    this.dpReadBuf = dev.createBuffer({
      size: MAX_DPS * DP_U32S * 4,
      usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    });

    // ── DP count (atomic) ──
    this.dpCountBuf = dev.createBuffer({
      size: 4,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
    });

    this.dpCountReadBuf = dev.createBuffer({
      size: 4,
      usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    });

    // ── Params ──
    const dpMask = 0xFFFFFFFF >>> (32 - this.opts.dpBits); // top dpBits bits masked
    const paramsData = new Uint32Array([
      this.opts.stepsPerDispatch,
      dpMask,
      numAnimals,
      MAX_DPS,
    ]);
    this.paramsBuf = dev.createBuffer({
      size: 16,
      usage: GPUBufferUsage.UNIFORM,
      mappedAtCreation: true,
    });
    new Uint32Array(this.paramsBuf.getMappedRange()).set(paramsData);
    this.paramsBuf.unmap();

    // ── Bind group ──
    this.bindGroup = dev.createBindGroup({
      layout: this.pipeline!.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: this.animalBuf } },
        { binding: 1, resource: { buffer: this.jumpBuf } },
        { binding: 2, resource: { buffer: this.dpBuf } },
        { binding: 3, resource: { buffer: this.dpCountBuf } },
        { binding: 4, resource: { buffer: this.paramsBuf } },
      ],
    });
  }

  // Dispatch one GPU pass and read back DPs
  async dispatch(): Promise<bigint | null> {
    if (!this.device || !this.pipeline) return null;
    const dev = this.device;

    const workgroups = Math.ceil(this.opts.numAnimals / 64);

    const enc = dev.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, this.bindGroup);
    pass.dispatchWorkgroups(workgroups);
    pass.end();

    // Copy DP count and DP buffer for readback
    enc.copyBufferToBuffer(this.dpCountBuf, 0, this.dpCountReadBuf, 0, 4);
    enc.copyBufferToBuffer(this.dpBuf, 0, this.dpReadBuf, 0, this.dpReadBuf.size);
    dev.queue.submit([enc.finish()]);

    // Read DP count
    await this.dpCountReadBuf.mapAsync(GPUMapMode.READ);
    const countArr = new Uint32Array(this.dpCountReadBuf.getMappedRange());
    const dpCount = Math.min(countArr[0], 8192);
    this.dpCountReadBuf.unmap();

    if (dpCount === 0) {
      this.steps += this.opts.numAnimals * this.opts.stepsPerDispatch;
      return null;
    }

    // Read DPs
    await this.dpReadBuf.mapAsync(GPUMapMode.READ);
    const dpArr = new Uint32Array(this.dpReadBuf.getMappedRange());

    let found: bigint | null = null;
    for (let i = 0; i < dpCount && !found; i++) {
      const base = i * DP_U32S;
      const cx = u256LEtoBigint(dpArr, base);
      const sc = u256LEtoBigint(dpArr, base + 8);
      const isWild = dpArr[base + 16] === 1;

      const key = cx.toString(16);
      const entry = this.dpTable.get(key);

      if (entry && entry.isWild !== isWild) {
        // Collision! Recover k
        const tSc = isWild ? entry.scalar : sc;
        const wSc = isWild ? sc : entry.scalar;
        found = this.recoverK(tSc, wSc);
      } else {
        this.dpTable.set(key, { canonX: cx, scalar: sc, isWild });
      }
      this.dps++;
    }

    this.dpReadBuf.unmap();
    this.steps += this.opts.numAnimals * this.opts.stepsPerDispatch;
    return found;
  }

  private recoverK(tameSc: bigint, wildSc: bigint): bigint | null {
    const lambdas = [1n, /* λ */ 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72n];
    for (const lam of lambdas) {
      for (const sign of [1n, -1n]) {
        const candidate = ((sign * (lam * tameSc % N) - wildSc) % N + N) % N;
        const check = scMul(candidate);
        if (!check.inf && check.x === this.opts.targetX && check.y === this.opts.targetY) {
          return candidate;
        }
      }
    }
    return null;
  }

  destroy() {
    this.animalBuf?.destroy();
    this.jumpBuf?.destroy();
    this.dpBuf?.destroy();
    this.dpReadBuf?.destroy();
    this.dpCountBuf?.destroy();
    this.dpCountReadBuf?.destroy();
    this.paramsBuf?.destroy();
    this.device?.destroy();
  }
}

// ─── Capability detection ─────────────────────────────────────────────────────

export async function detectWebGPU(): Promise<{
  available: boolean;
  adapterName: string;
  maxWorkgroups: number;
}> {
  if (!navigator.gpu) return { available: false, adapterName: 'none', maxWorkgroups: 0 };
  const adapter = await navigator.gpu.requestAdapter();
  if (!adapter) return { available: false, adapterName: 'none', maxWorkgroups: 0 };
  const info = await adapter.requestAdapterInfo();
  const dev = await adapter.requestDevice();
  const maxWG = dev.limits.maxComputeWorkgroupsPerDimension;
  dev.destroy();
  return {
    available: true,
    adapterName: info.description || info.vendor || 'GPU',
    maxWorkgroups: maxWG,
  };
}
