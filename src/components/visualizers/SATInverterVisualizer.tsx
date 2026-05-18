import React, { useState, useEffect, useRef } from 'react';
import { Play, Square, RotateCcw, Zap, Hash } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';
import {
  encodeSHA256,
  constrainHash,
  solveCDCL,
  extractInput,
  type SolveStats,
} from '../../sat/sha256sat';

// ─── Reduced SHA-256 (for the hash generator) ─────────────────────────────────

const SHA256_K_GEN = [
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];
const IV_GEN = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
const rotr = (n: number, x: number) => ((x >>> n) | (x << (32 - n))) >>> 0;
const add32 = (...ns: number[]) => ns.reduce((a, b) => (a + b) >>> 0, 0);
const hex8 = (n: number) => (n >>> 0).toString(16).padStart(8, '0');

function reducedSHA256(msg: string, rounds: number): string {
  const bytes = new TextEncoder().encode(msg);
  const l = bytes.length;
  const padLen = ((l + 9 + 63) & ~63);
  const buf = new Uint8Array(padLen);
  buf.set(bytes); buf[l] = 0x80;
  new DataView(buf.buffer).setUint32(padLen - 4, (l * 8) >>> 0, false);
  let [H0, H1, H2, H3, H4, H5, H6, H7] = IV_GEN;
  for (let bs = 0; bs < padLen; bs += 64) {
    const W = new Uint32Array(64);
    const bv = new DataView(buf.buffer, bs, 64);
    for (let i = 0; i < 16; i++) W[i] = bv.getUint32(i * 4, false);
    for (let i = 16; i < 64; i++) {
      const s0 = rotr(7, W[i-15]) ^ rotr(18, W[i-15]) ^ (W[i-15] >>> 3);
      const s1 = rotr(17, W[i-2]) ^ rotr(19, W[i-2]) ^ (W[i-2] >>> 10);
      W[i] = add32(W[i-16], s0, W[i-7], s1);
    }
    let a = H0, b = H1, c = H2, d = H3, e = H4, f = H5, g = H6, h = H7;
    for (let i = 0; i < rounds; i++) {
      const S1 = rotr(6, e) ^ rotr(11, e) ^ rotr(25, e);
      const ch2 = (e & f) ^ (~e & g);
      const t1 = add32(h, S1, ch2, SHA256_K_GEN[i], W[i]);
      const S0 = rotr(2, a) ^ rotr(13, a) ^ rotr(22, a);
      const maj2 = (a & b) ^ (a & c) ^ (b & c);
      const t2 = add32(S0, maj2);
      h = g; g = f; f = e; e = add32(d, t1);
      d = c; c = b; b = a; a = add32(t1, t2);
    }
    H0=add32(H0,a);H1=add32(H1,b);H2=add32(H2,c);H3=add32(H3,d);
    H4=add32(H4,e);H5=add32(H5,f);H6=add32(H6,g);H7=add32(H7,h);
  }
  return [H0,H1,H2,H3,H4,H5,H6,H7].map(hex8).join('');
}

// ─── Component ────────────────────────────────────────────────────────────────

interface BuildInfo { numVars: number; numClauses: number; }
interface SolverState {
  phase: 'idle' | 'building' | 'solving' | 'done';
  stats: SolveStats;
  numAssigned: number;
  result: string | null;
  status: 'sat' | 'unsat' | 'cancelled' | null;
}

const EMPTY_STATS: SolveStats = {
  decisions: 0, propagations: 0, conflicts: 0,
  learnedClauses: 0, backjumps: 0, maxBackjump: 0,
};

export function SATInverterVisualizer() {
  const [rounds, setRounds] = useState(4);
  const [inputBytes, setInputBytes] = useState(1);
  const [genMsg, setGenMsg] = useState('a');
  const [genHash, setGenHash] = useState('');
  const [targetHash, setTargetHash] = useState('');
  const [buildInfo, setBuildInfo] = useState<BuildInfo | null>(null);
  const [solver, setSolver] = useState<SolverState>({
    phase: 'idle', stats: EMPTY_STATS, numAssigned: 0, result: null, status: null,
  });
  const cancelRef = useRef({ cancelled: false });

  useEffect(() => {
    if (genMsg && genMsg.length <= 3) {
      setGenHash(reducedSHA256(genMsg, rounds));
    } else {
      setGenHash('');
    }
  }, [genMsg, rounds]);

  const reset = () => {
    cancelRef.current.cancelled = true;
    setBuildInfo(null);
    setSolver({ phase: 'idle', stats: EMPTY_STATS, numAssigned: 0, result: null, status: null });
  };

  const stop = () => { cancelRef.current.cancelled = true; };

  const runSolver = async () => {
    if (!targetHash || targetHash.length !== 64) return;
    cancelRef.current = { cancelled: false };

    setSolver({ phase: 'building', stats: EMPTY_STATS, numAssigned: 0, result: null, status: null });
    setBuildInfo(null);

    // Yield to let UI render
    await new Promise(r => setTimeout(r, 30));
    if (cancelRef.current.cancelled) return;

    // Build CNF
    const enc = encodeSHA256(inputBytes, rounds);
    const hashClauses = constrainHash(enc, targetHash);
    const allClauses = [...enc.clauses, ...hashClauses];
    const totalVars = enc.numVars;

    setBuildInfo({ numVars: totalVars, numClauses: allClauses.length });
    setSolver(s => ({ ...s, phase: 'solving' }));

    await new Promise(r => setTimeout(r, 30));
    if (cancelRef.current.cancelled) return;

    const result = await solveCDCL(
      allClauses,
      enc.numVars,
      (stats, assigned) => {
        setSolver(s => ({ ...s, stats: { ...stats }, numAssigned: assigned }));
      },
      cancelRef.current
    );

    if (cancelRef.current.cancelled) {
      setSolver(s => ({ ...s, phase: 'done', status: 'cancelled' }));
      return;
    }

    if (result.status === 'sat') {
      const inputBytesArr = extractInput(result.assignment, enc.inputLits);
      const decoded = Array.from(inputBytesArr).map(b => String.fromCharCode(b)).join('');
      setSolver({
        phase: 'done',
        stats: result.stats,
        numAssigned: result.assignment.size,
        result: decoded,
        status: 'sat',
      });
    } else {
      setSolver(s => ({
        ...s,
        phase: 'done',
        stats: result.stats,
        numAssigned: result.assignment.size,
        result: null,
        status: result.status as 'unsat' | 'cancelled',
      }));
    }
  };

  const roundWarning = rounds > 8
    ? `⚠ ${rounds} rounds: le solveur pourrait être très lent`
    : rounds > 6
    ? `Faisable en quelques secondes (${rounds}R)`
    : `Rapide (${rounds}R)`;

  const assignProgress = buildInfo
    ? Math.min((solver.numAssigned / buildInfo.numVars) * 100, 100)
    : 0;

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-4xl mx-auto space-y-8">

          {/* Header */}
          <header className="space-y-2">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase">
              CNF / CDCL
              <span className="text-[10px] text-zinc-500 font-mono align-top ml-3">SHA-256 RÉDUIT — ENCODAGE BOOLÉEN</span>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-tighter max-w-3xl leading-relaxed">
              SHA-256 réduit encodé comme système de clauses CNF (Tseitin).
              Le solveur CDCL encode SHA-256 en clauses CNF (Tseitin), apprend des clauses aux conflits (First UIP), et backjumpe de façon non-chronologique — state of the art.
              Contrairement au brute-force, chaque contradiction élimine un sous-espace entier.
            </p>
          </header>

          {/* Config */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-5">
            <div className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold">Paramètres</div>
            <div className="grid grid-cols-2 gap-6">
              <div className="space-y-2">
                <div className="flex justify-between">
                  <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Rounds actifs</label>
                  <span className="font-mono text-cyan-400 text-[10px]">{rounds}</span>
                </div>
                <input type="range" min={1} max={12} value={rounds}
                  onChange={e => { setRounds(+e.target.value); reset(); }}
                  className="w-full accent-cyan-500" />
                <div className="text-[9px] font-mono text-zinc-600">{roundWarning}</div>
              </div>
              <div className="space-y-2">
                <div className="flex justify-between">
                  <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Octets d'entrée libres</label>
                  <span className="font-mono text-purple-400 text-[10px]">{inputBytes}</span>
                </div>
                <input type="range" min={1} max={3} value={inputBytes}
                  onChange={e => { setInputBytes(+e.target.value); reset(); }}
                  className="w-full accent-purple-500" />
                <div className="text-[9px] font-mono text-zinc-600">{inputBytes * 8} bits libres dans le CNF</div>
              </div>
            </div>
          </div>

          {/* Hash generator */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-4">
            <div className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold flex items-center gap-2">
              <Hash className="w-3 h-3" /> Générateur de hash cible (SHA-256 à {rounds} rounds)
            </div>
            <div className="grid grid-cols-2 gap-4">
              <div className="bg-zinc-900/50 border border-zinc-800 p-3 rounded-lg">
                <label className="text-[9px] text-zinc-500 uppercase tracking-widest block mb-1">
                  Message (max {inputBytes} octet{inputBytes > 1 ? 's' : ''})
                </label>
                <input type="text" value={genMsg}
                  onChange={e => setGenMsg(e.target.value.slice(0, inputBytes))}
                  maxLength={inputBytes}
                  className="w-full bg-transparent text-orange-400 font-mono text-sm focus:outline-none" />
              </div>
              <div className="bg-zinc-900/50 border border-zinc-800 p-3 rounded-lg">
                <label className="text-[9px] text-zinc-500 uppercase tracking-widest block mb-1">Hash réduit</label>
                <div className="font-mono text-[10px] text-zinc-400 break-all">{genHash || '—'}</div>
              </div>
            </div>
            {genHash && (
              <button
                onClick={() => { setTargetHash(genHash); setGenMsg(''); reset(); }}
                className="text-[9px] uppercase tracking-widest font-bold border border-zinc-700 px-3 py-1.5 rounded text-zinc-400 hover:text-white hover:bg-zinc-800 transition"
              >
                Utiliser comme cible →
              </button>
            )}
          </div>

          {/* Target + solve */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-5">
            <div className="space-y-2">
              <label className="text-[10px] text-zinc-300 uppercase tracking-widest font-bold block">Hash cible (64 hex)</label>
              <input type="text" value={targetHash}
                onChange={e => { setTargetHash(e.target.value.toLowerCase().trim()); reset(); }}
                placeholder="Colle un hash SHA-256 réduit ici..."
                className="w-full bg-zinc-900/50 border border-zinc-800 rounded-lg p-3 font-mono text-[10px] text-zinc-300 focus:outline-none focus:border-zinc-600"
              />
            </div>

            <div className="flex items-center gap-3">
              <button onClick={runSolver}
                disabled={targetHash.length !== 64 || solver.phase === 'solving' || solver.phase === 'building'}
                className="px-5 py-2.5 bg-zinc-100 text-black text-[10px] font-bold uppercase tracking-widest hover:bg-white rounded transition disabled:opacity-40 disabled:cursor-not-allowed flex items-center gap-2"
              >
                <Play className="w-3 h-3" /> Résoudre (CDCL)
              </button>
              {(solver.phase === 'solving' || solver.phase === 'building') && (
                <button onClick={stop}
                  className="px-4 py-2 border border-red-700 text-red-400 text-[10px] font-bold uppercase tracking-widest rounded hover:bg-red-950/30 transition flex items-center gap-2"
                >
                  <Square className="w-3 h-3" /> Stop
                </button>
              )}
              <button onClick={reset}
                className="px-4 py-2 border border-zinc-800 text-zinc-500 text-[10px] font-bold uppercase tracking-widest rounded hover:text-zinc-300 transition flex items-center gap-2"
              >
                <RotateCcw className="w-3 h-3" /> Reset
              </button>
            </div>

            {/* Build info */}
            {buildInfo && (
              <div className="grid grid-cols-2 gap-3 text-[9px] font-mono">
                <div className="bg-zinc-900/50 border border-zinc-800 p-3 rounded-lg">
                  <div className="text-zinc-600 uppercase tracking-widest mb-1">Variables CNF</div>
                  <div className="text-cyan-400 text-sm font-bold">{buildInfo.numVars.toLocaleString()}</div>
                </div>
                <div className="bg-zinc-900/50 border border-zinc-800 p-3 rounded-lg">
                  <div className="text-zinc-600 uppercase tracking-widest mb-1">Clauses CNF</div>
                  <div className="text-purple-400 text-sm font-bold">{buildInfo.numClauses.toLocaleString()}</div>
                </div>
              </div>
            )}

            {/* Solver stats */}
            {solver.phase !== 'idle' && solver.phase !== 'building' && (
              <div className="space-y-3">
                <div className="text-[9px] text-zinc-500 uppercase tracking-widest flex items-center gap-2">
                  {solver.phase === 'solving' && <Zap className="w-3 h-3 animate-pulse text-yellow-500" />}
                  {solver.phase === 'solving' ? 'CDCL en cours...' :
                   solver.status === 'sat' ? 'SAT — préimage trouvée' :
                   solver.status === 'unsat' ? 'UNSAT — aucun préimage dans cet espace' :
                   'Annulé'}
                </div>

                {/* Variable assignment progress */}
                {buildInfo && (
                  <div className="space-y-1">
                    <div className="flex justify-between text-[9px] font-mono text-zinc-600">
                      <span>Variables assignées</span>
                      <span>{solver.numAssigned.toLocaleString()} / {buildInfo.numVars.toLocaleString()}</span>
                    </div>
                    <div className="h-px bg-zinc-800 rounded-full overflow-hidden">
                      <div
                        className={`h-full transition-all duration-100 ${
                          solver.status === 'sat' ? 'bg-emerald-500' : 'bg-cyan-500'
                        }`}
                        style={{ width: `${assignProgress}%` }}
                      />
                    </div>
                  </div>
                )}

                {/* Stats grid */}
                <div className="grid grid-cols-3 gap-2 text-[9px] font-mono">
                  {[
                    { label: 'Décisions', val: solver.stats.decisions, color: 'text-cyan-400' },
                    { label: 'Conflits', val: solver.stats.conflicts, color: 'text-red-400' },
                    { label: 'Propagations', val: solver.stats.propagations, color: 'text-purple-400' },
                    { label: 'Clauses apprises', val: solver.stats.learnedClauses, color: 'text-yellow-400' },
                    { label: 'Backjumps', val: solver.stats.backjumps, color: 'text-orange-400' },
                    { label: 'Max backjump', val: solver.stats.maxBackjump, color: 'text-emerald-400', suffix: ' niveaux' },
                  ].map(({ label, val, color, suffix }) => (
                    <div key={label} className="bg-zinc-900/50 border border-zinc-800 p-2 rounded">
                      <div className="text-zinc-600 uppercase tracking-widest text-[8px] mb-1">{label}</div>
                      <div className={`font-bold ${color}`}>{val.toLocaleString()}{suffix ?? ''}</div>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {/* Result */}
            <AnimatePresence>
              {solver.phase === 'done' && solver.status === 'sat' && solver.result !== null && (
                <motion.div
                  initial={{ opacity: 0, y: 8 }}
                  animate={{ opacity: 1, y: 0 }}
                  className="bg-emerald-500/10 border border-emerald-500/30 rounded-lg p-6 space-y-3"
                >
                  <div className="text-[9px] text-emerald-500/70 uppercase tracking-widest font-bold">
                    Préimage trouvée — SHA-256 ({rounds} rounds), {inputBytes} octet{inputBytes > 1 ? 's' : ''}
                  </div>
                  <div className="font-mono text-2xl text-emerald-400">"{solver.result}"</div>
                  <div className="text-[9px] text-zinc-500 font-mono pt-2 border-t border-emerald-500/10 space-y-0.5">
                    <div>
                      {solver.stats.conflicts} conflits · {solver.stats.learnedClauses} clauses apprises · backjump max {solver.stats.maxBackjump} niveaux
                    </div>
                    <div className="text-zinc-600">
                      CDCL: {solver.stats.decisions} décisions vs 2^{inputBytes * 8} = {Math.pow(2, inputBytes * 8).toLocaleString()} pour un brute-force. Chaque conflit élimine un sous-espace entier.
                    </div>
                  </div>
                </motion.div>
              )}
              {solver.phase === 'done' && solver.status === 'unsat' && (
                <motion.div
                  initial={{ opacity: 0, y: 8 }}
                  animate={{ opacity: 1, y: 0 }}
                  className="bg-red-500/10 border border-red-500/20 rounded-lg p-4"
                >
                  <div className="text-[9px] text-red-400 uppercase tracking-widest font-bold mb-1">UNSAT</div>
                  <div className="text-[9px] text-zinc-500 font-mono">
                    Aucun préimage de {inputBytes} octet{inputBytes > 1 ? 's' : ''} n'existe pour ce hash avec {rounds} rounds.
                    Essaie avec plus d'octets ou vérifie que le hash a été généré avec les mêmes paramètres.
                  </div>
                </motion.div>
              )}
            </AnimatePresence>
          </div>

          {/* Roadmap */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-3">
            <div className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold">Roadmap</div>
            {[
              { label: 'Encodage Tseitin + DPLL', done: true, desc: 'XOR/AND/NOT/ADD → clauses CNF · propagation de contraintes' },
              { label: 'CDCL + First UIP + VSIDS', done: true, desc: 'Conflict-Driven Clause Learning · backjump non-chronologique · activité des variables' },
              { label: 'Analyse différentielle', done: false, desc: 'Propager les différences pour contraindre le schedule W' },
              { label: 'Extension à SHA-256 complet (64R)', done: false, desc: 'Horizon ouvert — nécessite des avancées en cryptanalyse' },
            ].map((step, i) => (
              <div key={i} className="flex items-start gap-3">
                <div className={`w-1.5 h-1.5 rounded-full mt-1.5 shrink-0 ${step.done ? 'bg-emerald-500' : 'bg-zinc-700'}`} />
                <div>
                  <div className={`text-[10px] font-bold uppercase tracking-widest ${step.done ? 'text-emerald-400' : 'text-zinc-400'}`}>{step.label}</div>
                  <div className="text-[9px] text-zinc-600 font-mono">{step.desc}</div>
                </div>
              </div>
            ))}
          </div>

        </div>
      </div>
    </div>
  );
}
