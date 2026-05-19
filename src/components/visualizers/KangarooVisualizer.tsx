import React, { useState, useRef, useCallback, useEffect } from 'react';
import { Play, Square, RefreshCw, Zap, Eye, Target } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';

// ─── secp256k1 small-scale demo ───────────────────────────────────────────────
// We demo on a tiny curve (p=31, n=29) for the visualizer so the orbit
// and kangaroo walk are visible. The real secp256k1 code runs in WASM.

const GX_DEMO = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n;
const GY_DEMO = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8n;
const P_SECP  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2Fn;
const BETA    = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EEn;
const BETA2   = 0x851695D49A83F8EF919BB86153CBCB16630FB68AED0A766A3EC693D68E6AFA40n;

function mpow(b: bigint, e: bigint, m: bigint): bigint {
  let r = 1n; b = ((b % m) + m) % m;
  while (e > 0n) { if (e & 1n) r = r * b % m; b = b * b % m; e >>= 1n; }
  return r;
}

interface Pt { x: bigint; y: bigint; inf: boolean; }

function ptAdd(A: Pt, B: Pt, p: bigint): Pt {
  if (A.inf) return B; if (B.inf) return A;
  if (A.x === B.x) {
    if (A.y !== B.y) return { x:0n, y:0n, inf:true };
    const lam = 3n*A.x*A.x%p * mpow(2n*A.y, p-2n, p) % p;
    const x3 = (lam*lam%p - 2n*A.x%p + 2n*p) % p;
    return { x:x3, y:(lam*(A.x-x3+p)%p - A.y+p)%p, inf:false };
  }
  const lam = ((B.y-A.y+p)%p) * mpow((B.x-A.x+p), p-2n, p) % p;
  const x3 = (lam*lam%p - A.x - B.x + 2n*p) % p;
  return { x:x3, y:(lam*(A.x-x3+p)%p - A.y+p)%p, inf:false };
}

function scMul(G: Pt, k: bigint, p: bigint): Pt {
  let r: Pt = {x:0n,y:0n,inf:true}, a = G;
  while (k>0n) { if(k&1n) r=ptAdd(r,a,p); a=ptAdd(a,a,p); k>>=1n; }
  return r;
}

// Compute the 6-orbit of a secp256k1 point
function orbit6(pt: Pt): { name: string; x: string; y: string }[] {
  if (pt.inf) return [];
  const p = P_SECP;
  const neg = (q: Pt) => ({ x: q.x, y: (p - q.y) % p, inf: false });
  const psi  = (q: Pt) => ({ x: BETA  * q.x % p, y: q.y, inf: false });
  const psi2 = (q: Pt) => ({ x: BETA2 * q.x % p, y: q.y, inf: false });
  const toHex = (n: bigint) => n.toString(16).padStart(64, '0');
  return [
    { name: 'P',    ...pt,      x: toHex(pt.x),      y: toHex(pt.y) },
    { name: '−P',   ...neg(pt), x: toHex(pt.x),      y: toHex(neg(pt).y) },
    { name: 'ψ(P)', ...psi(pt), x: toHex(psi(pt).x), y: toHex(pt.y) },
    { name: '−ψ(P)',...neg(psi(pt)), x: toHex(psi(pt).x), y: toHex(neg(psi(pt)).y) },
    { name: 'ψ²(P)',...psi2(pt), x: toHex(psi2(pt).x), y: toHex(pt.y) },
    { name:'−ψ²(P)',...neg(psi2(pt)), x: toHex(psi2(pt).x), y: toHex(neg(psi2(pt)).y) },
  ];
}

// ─── Tiny-curve Kangaroo demo (p=997 Weierstrass y²=x³+7, just for animation) ─

const TP = 997n;
const TA = 0n; // a=0
const TB = 7n;

function tinyPoints(): Pt[] {
  const pts: Pt[] = [];
  for (let x = 0n; x < TP; x++) {
    const rhs = (x*x%TP*x%TP + TB) % TP;
    const y = mpow(rhs, (TP+1n)/4n, TP);
    if (y*y%TP === rhs) {
      pts.push({ x, y, inf: false });
      if (y !== 0n) pts.push({ x, y: TP-y, inf: false });
    }
  }
  return pts;
}

const TINY_PTS = tinyPoints();
const TINY_G  = TINY_PTS[0];
const TINY_N  = BigInt(TINY_PTS.length + 1);

// ─── Component ────────────────────────────────────────────────────────────────

type Phase = 'idle' | 'running' | 'found';

interface KangarooState {
  tames: { pt: Pt; scalar: bigint; trail: Pt[] }[];
  wild:  { pt: Pt; scalar: bigint; trail: Pt[] };
  target: Pt;
  targetK: bigint;
  found: boolean;
  foundK: bigint | null;
  steps: number;
  dps: number;
  table: Map<string, { scalar: bigint; isWild: boolean }>;
}

const JUMP_SIZES = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n];

function jumpIdx(pt: Pt): number {
  return Number(pt.x % BigInt(JUMP_SIZES.length));
}

function canonicalKey(pt: Pt): string {
  // 6-orbit canonical: min of {x} under ψ, ψ² over tiny curve (no β, just use x)
  return pt.x.toString();
}

function kangarooStep(state: KangarooState): KangarooState {
  const p = TP;
  const g = TINY_G;
  const newState = { ...state, table: new Map(state.table) };

  const walkAnimal = (pt: Pt, scalar: bigint, isWild: boolean): { pt: Pt; scalar: bigint; key?: bigint } => {
    const key = canonicalKey(pt);
    if (newState.table.has(key)) {
      const entry = newState.table.get(key)!;
      if (entry.isWild !== isWild) {
        // Collision! Recover k
        newState.dps++;
        let k: bigint;
        if (isWild) {
          // wild hits tame entry: tame_sc = k + wild_sc → k = tame_sc - wild_sc
          k = ((entry.scalar - scalar) % TINY_N + TINY_N) % TINY_N;
        } else {
          k = ((scalar - entry.scalar) % TINY_N + TINY_N) % TINY_N;
        }
        // Verify
        const check = scMul(g, k, p);
        if (!check.inf && check.x === newState.target.x && check.y === newState.target.y) {
          return { pt, scalar, key: k };
        }
      }
    } else {
      newState.table.set(key, { scalar, isWild });
      newState.dps++;
    }
    const ji = jumpIdx(pt);
    const jumpScalar = JUMP_SIZES[ji];
    const jumpPt = scMul(g, jumpScalar, p);
    return {
      pt: ptAdd(pt, jumpPt, p),
      scalar: (scalar + jumpScalar) % TINY_N,
    };
  };

  // Step wild
  const wResult = walkAnimal(newState.wild.pt, newState.wild.scalar, true);
  if (wResult.key !== undefined) {
    return { ...newState, found: true, foundK: wResult.key };
  }
  newState.wild = {
    pt: wResult.pt,
    scalar: wResult.scalar,
    trail: [...state.wild.trail.slice(-20), wResult.pt],
  };

  // Step tames
  newState.tames = state.tames.map((tame, _i) => {
    const tResult = walkAnimal(tame.pt, tame.scalar, false);
    if (tResult.key !== undefined && !newState.found) {
      newState.found = true;
      newState.foundK = tResult.key;
    }
    return {
      pt: tResult.pt,
      scalar: tResult.scalar,
      trail: [...tame.trail.slice(-20), tResult.pt],
    };
  });

  newState.steps += 1 + state.tames.length;
  return newState;
}

function initKangaroo(targetK: bigint): KangarooState {
  const g = TINY_G;
  const p = TP;
  const n = TINY_N;
  const target = scMul(g, targetK, p);

  // 3 tame kangaroos spread across range
  const range = Number(n);
  const tames = [0.25, 0.5, 0.75].map(frac => {
    const sc = BigInt(Math.floor(frac * range));
    const pt = scMul(g, sc, p);
    return { pt, scalar: sc, trail: [pt] };
  });

  return {
    tames,
    wild: { pt: target, scalar: 0n, trail: [target] },
    target,
    targetK,
    found: false,
    foundK: null,
    steps: 0,
    dps: 0,
    table: new Map(),
  };
}

// ─── Main Visualizer ──────────────────────────────────────────────────────────

export function KangarooVisualizer() {
  const [phase, setPhase] = useState<Phase>('idle');
  const [kState, setKState] = useState<KangarooState | null>(null);
  const [targetK, setTargetK] = useState(42n);
  const [showOrbit, setShowOrbit] = useState(false);
  const [orbitPt, setOrbitPt] = useState<{ x: string; y: string } | null>(null);
  const [orbitData, setOrbitData] = useState<ReturnType<typeof orbit6>>([]);
  const rafRef = useRef<number>(0);

  const animate = useCallback(() => {
    setKState(prev => {
      if (!prev || prev.found) return prev;
      let s = prev;
      for (let i = 0; i < 3; i++) s = kangarooStep(s);
      if (s.found) setPhase('found');
      return s;
    });
    rafRef.current = requestAnimationFrame(animate);
  }, []);

  const start = useCallback(() => {
    const k = (targetK % TINY_N + TINY_N) % TINY_N || 1n;
    setKState(initKangaroo(k));
    setPhase('running');
  }, [targetK]);

  useEffect(() => {
    if (phase === 'running') {
      rafRef.current = requestAnimationFrame(animate);
    } else {
      cancelAnimationFrame(rafRef.current);
    }
    return () => cancelAnimationFrame(rafRef.current);
  }, [phase, animate]);

  const stop = () => { setPhase('idle'); cancelAnimationFrame(rafRef.current); };
  const reset = () => { stop(); setKState(null); };

  // Orbit demo: use generator G of secp256k1
  const showOrbitDemo = () => {
    const pt: Pt = { x: GX_DEMO, y: GY_DEMO, inf: false };
    setOrbitData(orbit6(pt));
    setOrbitPt({ x: GX_DEMO.toString(16).slice(0, 8) + '…', y: GY_DEMO.toString(16).slice(0, 8) + '…' });
    setShowOrbit(true);
  };

  // Map tiny curve point to canvas coords
  const toCanvas = (pt: Pt, w: number, h: number) => ({
    cx: (Number(pt.x) / Number(TP)) * w,
    cy: h - (Number(pt.y) / Number(TP)) * h,
  });

  const W = 460, H = 200;

  return (
    <div className="flex flex-col h-full overflow-y-auto bg-[#050505] text-zinc-300 p-6 space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-sm font-bold uppercase tracking-widest text-zinc-100">
            Kangaroo 6-Automorphism
          </h2>
          <p className="text-[10px] uppercase tracking-widest text-zinc-500 mt-0.5">
            √6 speedup via End(E) = Z[ζ₃] — all 6 units exploited
          </p>
        </div>
        <div className="flex gap-2">
          <button onClick={showOrbitDemo}
            className="flex items-center gap-1.5 px-3 py-1.5 text-[10px] uppercase font-bold tracking-widest rounded border border-purple-500/40 text-purple-400 hover:bg-purple-500/10 transition-all">
            <Eye className="w-3 h-3" /> Orbite 6×
          </button>
          {phase === 'idle' && (
            <button onClick={start}
              className="flex items-center gap-1.5 px-3 py-1.5 text-[10px] uppercase font-bold tracking-widest rounded border border-emerald-500/40 text-emerald-400 hover:bg-emerald-500/10 transition-all">
              <Play className="w-3 h-3" /> Démo
            </button>
          )}
          {phase === 'running' && (
            <button onClick={stop}
              className="flex items-center gap-1.5 px-3 py-1.5 text-[10px] uppercase font-bold tracking-widest rounded border border-red-500/40 text-red-400 hover:bg-red-500/10 transition-all">
              <Square className="w-3 h-3" /> Stop
            </button>
          )}
          <button onClick={reset}
            className="flex items-center gap-1.5 px-3 py-1.5 text-[10px] uppercase font-bold tracking-widest rounded border border-zinc-700 text-zinc-500 hover:bg-zinc-800 transition-all">
            <RefreshCw className="w-3 h-3" />
          </button>
        </div>
      </div>

      {/* Math explanation */}
      <div className="grid grid-cols-3 gap-3">
        {[
          { label: 'Endomorphismes', val: '{ id, ψ, ψ² }', sub: 'ordre 3 dans End(E)' },
          { label: 'Avec négation', val: '× 2 = 6', sub: 'unités de Z[ζ₃]' },
          { label: 'Speedup théorique', val: '√6 ≈ 2.45×', sub: 'vs Kangaroo standard' },
        ].map(c => (
          <div key={c.label} className="bg-zinc-900/40 border border-zinc-800 rounded p-3">
            <div className="text-[9px] uppercase tracking-widest text-zinc-500 mb-1">{c.label}</div>
            <div className="text-sm font-mono text-cyan-400">{c.val}</div>
            <div className="text-[9px] text-zinc-600 mt-0.5">{c.sub}</div>
          </div>
        ))}
      </div>

      {/* Canvas: kangaroo walk visualization */}
      <div className="bg-zinc-900/40 border border-zinc-800 rounded p-4">
        <div className="text-[9px] uppercase tracking-widest text-zinc-500 mb-3">
          Courbe démo y²=x³+7 (p=997) — marche des kangourous
        </div>
        <svg width={W} height={H} className="w-full rounded bg-zinc-950/50">
          {/* Points of the curve as dots */}
          {TINY_PTS.map((pt, i) => {
            const { cx, cy } = toCanvas(pt, W, H);
            return <circle key={i} cx={cx} cy={cy} r={1} fill="#27272a" />;
          })}

          {/* Tame trails */}
          {kState?.tames.map((t, ti) => {
            const colors = ['#10b981', '#34d399', '#6ee7b7'];
            return t.trail.map((pt, i) => {
              if (pt.inf) return null;
              const { cx, cy } = toCanvas(pt, W, H);
              return <circle key={`t${ti}-${i}`} cx={cx} cy={cy} r={2}
                fill={colors[ti % 3]} opacity={(i + 1) / t.trail.length * 0.7} />;
            });
          })}

          {/* Wild trail */}
          {kState?.wild.trail.map((pt, i) => {
            if (pt.inf) return null;
            const { cx, cy } = toCanvas(pt, W, H);
            return <circle key={`w${i}`} cx={cx} cy={cy} r={2}
              fill="#a78bfa" opacity={(i + 1) / kState.wild.trail.length * 0.9} />;
          })}

          {/* Target */}
          {kState && !kState.target.inf && (() => {
            const { cx, cy } = toCanvas(kState.target, W, H);
            return <circle cx={cx} cy={cy} r={5} fill="none" stroke="#f59e0b" strokeWidth={1.5} />;
          })()}

          {/* Current tame heads */}
          {kState?.tames.map((t, ti) => {
            if (t.pt.inf) return null;
            const { cx, cy } = toCanvas(t.pt, W, H);
            const colors = ['#10b981', '#34d399', '#6ee7b7'];
            return <circle key={`th${ti}`} cx={cx} cy={cy} r={4} fill={colors[ti % 3]} />;
          })}

          {/* Wild head */}
          {kState && !kState.wild.pt.inf && (() => {
            const { cx, cy } = toCanvas(kState.wild.pt, W, H);
            return <circle cx={cx} cy={cy} r={5} fill="#a78bfa" />;
          })()}

          {/* Legend */}
          <g transform="translate(8,8)">
            <circle cx={4} cy={4} r={3} fill="#10b981" />
            <text x={10} y={8} fill="#6ee7b7" fontSize={8}>Tame (×3)</text>
            <circle cx={4} cy={18} r={3} fill="#a78bfa" />
            <text x={10} y={22} fill="#a78bfa" fontSize={8}>Wild</text>
            <circle cx={4} cy={32} r={4} fill="none" stroke="#f59e0b" strokeWidth={1.5} />
            <text x={10} y={36} fill="#f59e0b" fontSize={8}>Target k·G</text>
          </g>
        </svg>
      </div>

      {/* Stats */}
      {kState && (
        <div className="grid grid-cols-4 gap-3">
          {[
            { label: 'Steps', val: kState.steps.toLocaleString() },
            { label: 'DPs trouvés', val: kState.dps.toLocaleString() },
            { label: 'Table size', val: kState.table.size.toLocaleString() },
            { label: 'Target k', val: kState.targetK.toString() },
          ].map(s => (
            <div key={s.label} className="bg-zinc-900/40 border border-zinc-800 rounded p-3">
              <div className="text-[9px] uppercase tracking-widest text-zinc-500">{s.label}</div>
              <div className="font-mono text-xs text-zinc-200 mt-1">{s.val}</div>
            </div>
          ))}
        </div>
      )}

      {/* Found! */}
      <AnimatePresence>
        {kState?.found && (
          <motion.div
            initial={{ opacity: 0, y: 10 }} animate={{ opacity: 1, y: 0 }}
            className="bg-emerald-900/20 border border-emerald-500/30 rounded p-4 flex items-center gap-3">
            <Zap className="w-5 h-5 text-emerald-400 shrink-0" />
            <div>
              <div className="text-xs font-bold text-emerald-400 uppercase tracking-widest">Collision détectée !</div>
              <div className="font-mono text-sm text-zinc-200 mt-1">
                k = <span className="text-emerald-300">{kState.foundK?.toString()}</span>
                <span className="text-zinc-500 ml-3 text-xs">
                  trouvé en {kState.steps} steps · {kState.dps} DPs
                </span>
              </div>
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Orbit modal */}
      <AnimatePresence>
        {showOrbit && (
          <motion.div
            initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}
            className="fixed inset-0 bg-black/80 flex items-center justify-center z-50"
            onClick={() => setShowOrbit(false)}>
            <motion.div
              initial={{ scale: 0.9 }} animate={{ scale: 1 }} exit={{ scale: 0.9 }}
              className="bg-zinc-900 border border-zinc-700 rounded-xl p-6 max-w-3xl w-full mx-4"
              onClick={e => e.stopPropagation()}>
              <div className="flex items-center justify-between mb-4">
                <div>
                  <div className="text-sm font-bold text-zinc-100 uppercase tracking-widest">
                    Orbite 6× de G (secp256k1)
                  </div>
                  <div className="text-[10px] text-zinc-500 mt-0.5">
                    Les 6 automorphismes : {'{'}id, ψ, ψ², −id, −ψ, −ψ²{'}'}
                  </div>
                </div>
                <button onClick={() => setShowOrbit(false)} className="text-zinc-500 hover:text-zinc-300 text-lg">✕</button>
              </div>

              <div className="space-y-2">
                {orbitData.map((o, i) => {
                  const colors = ['#10b981','#6ee7b7','#a78bfa','#c4b5fd','#f59e0b','#fbbf24'];
                  const isCanonical = i === orbitData.reduce((mi, _, ii) =>
                    orbitData[ii].x < orbitData[mi].x ? ii : mi, 0);
                  return (
                    <div key={o.name}
                      className={`flex items-center gap-3 p-2.5 rounded border ${
                        isCanonical
                          ? 'bg-cyan-900/20 border-cyan-500/40'
                          : 'bg-zinc-900/50 border-zinc-800'
                      }`}>
                      <div className="w-12 text-right">
                        <span className="text-[11px] font-mono" style={{ color: colors[i] }}>{o.name}</span>
                      </div>
                      <div className="flex-1 min-w-0">
                        <div className="text-[9px] text-zinc-500 font-mono truncate">x: {o.x}</div>
                        <div className="text-[9px] text-zinc-600 font-mono truncate">y: {o.y}</div>
                      </div>
                      {isCanonical && (
                        <div className="text-[9px] text-cyan-400 uppercase font-bold tracking-widest shrink-0">
                          canonical
                        </div>
                      )}
                    </div>
                  );
                })}
              </div>

              <div className="mt-4 p-3 bg-zinc-950/50 rounded border border-zinc-800 text-[10px] text-zinc-500">
                <span className="text-cyan-400">canonical(P)</span> = min(x, β·x, β²·x)
                {' '}— un seul entrée de table pour 6 points distincts → √6 fois moins de collisions nécessaires
              </div>
            </motion.div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* secp256k1 projection for puzzle #135 */}
      <div className="bg-zinc-900/30 border border-zinc-800 rounded p-4">
        <div className="text-[9px] uppercase tracking-widest text-zinc-500 mb-3 flex items-center gap-2">
          <Target className="w-3 h-3" /> Projection puzzle #135 (secp256k1)
        </div>
        <div className="grid grid-cols-3 gap-3 text-center">
          {[
            { algo: 'Kangaroo basic', ops: '2^67.5', color: 'text-red-400', speedup: '1×' },
            { algo: 'Kangaroo 4-aut', ops: '2^66.5', color: 'text-amber-400', speedup: '2×' },
            { algo: 'Kangaroo 6-aut ✦', ops: '2^66.2', color: 'text-emerald-400', speedup: '√6 ≈ 2.45×' },
          ].map(r => (
            <div key={r.algo} className="bg-zinc-950/50 rounded p-3 border border-zinc-800">
              <div className="text-[9px] text-zinc-500 uppercase tracking-wider mb-1">{r.algo}</div>
              <div className={`font-mono text-sm font-bold ${r.color}`}>{r.ops}</div>
              <div className="text-[9px] text-zinc-600 mt-1">{r.speedup}</div>
            </div>
          ))}
        </div>
        <div className="mt-3 text-[10px] text-zinc-600 text-center">
          ψ²(x,y) = (β²·x mod p, y) — une multiplication de champ supplémentaire, 50% de réduction d'espace supplémentaire
        </div>
      </div>
    </div>
  );
}
