import React, { useState, useRef, useCallback } from 'react';
import { Play, Square, RefreshCw, Lock, Unlock, ChevronRight } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';
import {
  buildTinyPuzzle,
  encodeECDLP,
  extractKey,
  estimateSecp256k1Clauses,
  verifyKey,
  DEMO_PRIMES,
  type TinyCurveParams,
  type ECDLPEncoding,
} from '../../sat/ecsat';
import { solveRust as solveCDCL, type RustSolveStats as SolveStats } from '../../sat/rustSolver';

// Mersenne primes only (required for our efficient modular reduction encoding)

// ─── Helpers ─────────────────────────────────────────────────────────────────
const fmt = (n: number) =>
  n >= 1_000_000_000 ? `${(n / 1e9).toFixed(1)}B`
  : n >= 1_000_000   ? `${(n / 1e6).toFixed(1)}M`
  : n >= 1_000       ? `${(n / 1e3).toFixed(1)}K`
  : `${n}`;

const hex = (n: bigint, len = 8) => n.toString(16).padStart(len, '0').toUpperCase();

// ─── Component ───────────────────────────────────────────────────────────────
export function ECPuzzleVisualizer() {
  const [prime, setPrime] = useState(31);
  const [unknownBits, setUnknownBits] = useState(4);
  const [puzzle, setPuzzle] = useState<TinyCurveParams | null>(null);
  const [encoding, setEncoding] = useState<ECDLPEncoding | null>(null);
  const [phase, setPhase] = useState<'idle' | 'encoding' | 'solving' | 'solved' | 'unsat' | 'error'>('idle');
  const [stats, setStats] = useState<SolveStats | null>(null);
  const [assigned, setAssigned] = useState(0);
  const [foundK, setFoundK] = useState<bigint | null>(null);
  const [log, setLog] = useState<string[]>([]);
  const cancelRef = useRef({ cancelled: false });

  const pushLog = (msg: string) => setLog(prev => [...prev.slice(-11), msg]);

  const handleGenerate = useCallback(() => {
    setPhase('encoding');
    setFoundK(null);
    setStats(null);
    setEncoding(null);
    setAssigned(0);
    setLog([]);

    setTimeout(() => {
      let p: TinyCurveParams | null = null;
      // Try this prime, then fall back
      for (const tryP of [prime, ...DEMO_PRIMES.filter(x => x !== prime)]) {
        p = buildTinyPuzzle(tryP, unknownBits);
        if (p) break;
      }
      if (!p) { setPhase('error'); return; }
      setPuzzle(p);
      pushLog(`Courbe: y² = x³ + 7 mod ${p.p}  |  G=(${p.Gx},${p.Gy})`);
      pushLog(`Clé privée secrète: k = ${p.k}  (${unknownBits} bits)`);
      pushLog(`Clé publique: P = (${p.Px}, ${p.Py})`);
      pushLog(`Encodage CNF en cours…`);
      setTimeout(() => {
        try {
          const enc = encodeECDLP(p!);
          setEncoding(enc);
          pushLog(`✓ ${fmt(enc.numVars)} variables, ${fmt(enc.clauses.length)} clauses`);
          pushLog(`Corps: F_${p!.p} (${enc.fieldBits} bits/élément)`);
          pushLog(`Prêt — CDCL va chercher k…`);
          setPhase('idle');
        } catch (e) {
          pushLog(`Erreur encodage: ${e}`);
          setPhase('error');
        }
      }, 50);
    }, 20);
  }, [prime, unknownBits]);

  const handleSolve = useCallback(async () => {
    if (!encoding || !puzzle) return;
    cancelRef.current = { cancelled: false };
    setPhase('solving');
    setStats(null);
    setAssigned(0);
    pushLog(`CDCL démarré — ${fmt(encoding.clauses.length)} clauses, ${fmt(encoding.numVars)} vars`);

    const result = await solveCDCL(
      encoding.clauses,
      encoding.numVars,
      (s, asgn) => {
        setStats({ ...s });
        setAssigned(asgn);
        if (s.conflicts % 500 === 0 && s.conflicts > 0) {
          pushLog(`↺ ${fmt(s.conflicts)} conflits, ${fmt(s.learnedClauses)} clauses apprises`);
        }
      },
      cancelRef.current,
    );

    if (result.status === 'sat') {
      const k = extractKey(result.assignment!, encoding.kBits);
      setFoundK(k);
      const valid = verifyKey(k, BigInt(puzzle.Px), BigInt(puzzle.Py));
      pushLog(`✓ SAT! k = ${k} ${valid ? '✓ vérifié' : '✗ ERREUR'}`);
      pushLog(`Clé attendue: ${puzzle.k}  |  Trouvée: ${k}`);
      setPhase('solved');
    } else if (result.status === 'unsat') {
      pushLog('UNSAT — aucune solution (vérifiez les paramètres)');
      setPhase('unsat');
    } else {
      pushLog('Interrompu');
      setPhase('idle');
    }
  }, [encoding, puzzle]);

  const handleStop = () => { cancelRef.current.cancelled = true; setPhase('idle'); };

  const secp256kStats = estimateSecp256k1Clauses(unknownBits);

  // ─── Render ───────────────────────────────────────────────────────────────
  return (
    <div className="flex flex-col h-full overflow-hidden">
      {/* Header */}
      <div className="border-b border-zinc-800 px-6 py-4 flex items-center justify-between shrink-0">
        <div>
          <h2 className="text-xs uppercase font-bold tracking-widest text-zinc-100">
            EC Puzzle Solver — CDCL/SAT
          </h2>
          <p className="text-[10px] text-zinc-500 mt-0.5 uppercase tracking-wider">
            Encode k×G = P comme CNF · Résout via CDCL
          </p>
        </div>
        <div className="flex items-center gap-2">
          {phase === 'solved' && (
            <div className="flex items-center gap-2 bg-emerald-900/30 border border-emerald-500/30 px-3 py-1.5 rounded">
              <Unlock className="w-3 h-3 text-emerald-400" />
              <span className="text-[10px] text-emerald-300 font-mono uppercase tracking-wider">Clé trouvée</span>
            </div>
          )}
        </div>
      </div>

      <div className="flex flex-1 overflow-hidden">
        {/* Left panel: config + stats */}
        <div className="w-72 border-r border-zinc-800 flex flex-col overflow-y-auto shrink-0">
          {/* Config */}
          <div className="p-4 border-b border-zinc-800">
            <div className="text-[9px] uppercase tracking-widest text-zinc-500 font-bold mb-3">Configuration</div>

            <label className="block mb-3">
              <span className="text-[9px] uppercase tracking-widest text-zinc-600">Prime p (corps fini)</span>
              <select
                value={prime}
                onChange={e => setPrime(Number(e.target.value))}
                className="mt-1 w-full bg-zinc-900 border border-zinc-700 text-zinc-300 text-xs rounded px-2 py-1.5 font-mono"
                disabled={phase === 'solving' || phase === 'encoding'}
              >
                {DEMO_PRIMES.map(p => (
                  <option key={p} value={p}>p = {p} (Mersenne 2^{Math.round(Math.log2(p+1))}-1)</option>
                ))}
              </select>
            </label>

            <label className="block mb-4">
              <span className="text-[9px] uppercase tracking-widest text-zinc-600">
                Bits inconnus de k: {unknownBits}
              </span>
              <input
                type="range" min={2} max={7} value={unknownBits}
                onChange={e => setUnknownBits(Number(e.target.value))}
                className="mt-1 w-full accent-cyan-500"
                disabled={phase === 'solving' || phase === 'encoding'}
              />
              <div className="flex justify-between text-[8px] text-zinc-600 font-mono">
                <span>2</span><span>7</span>
              </div>
            </label>

            <div className="flex gap-2">
              <button
                onClick={handleGenerate}
                disabled={phase === 'solving' || phase === 'encoding'}
                className="flex-1 flex items-center justify-center gap-1.5 px-3 py-2 text-[10px] uppercase tracking-wider font-bold rounded border border-cyan-500/30 text-cyan-400 hover:bg-cyan-900/20 disabled:opacity-40 transition-all"
              >
                <RefreshCw className="w-3 h-3" />
                Générer
              </button>
              {phase === 'solving' ? (
                <button
                  onClick={handleStop}
                  className="flex-1 flex items-center justify-center gap-1.5 px-3 py-2 text-[10px] uppercase tracking-wider font-bold rounded border border-red-500/30 text-red-400 hover:bg-red-900/20 transition-all"
                >
                  <Square className="w-3 h-3" />Stop
                </button>
              ) : (
                <button
                  onClick={handleSolve}
                  disabled={!encoding || phase === 'encoding'}
                  className="flex-1 flex items-center justify-center gap-1.5 px-3 py-2 text-[10px] uppercase tracking-wider font-bold rounded border border-emerald-500/30 text-emerald-400 hover:bg-emerald-900/20 disabled:opacity-40 transition-all"
                >
                  <Play className="w-3 h-3" />Solve
                </button>
              )}
            </div>
          </div>

          {/* Encoding stats */}
          {encoding && (
            <div className="p-4 border-b border-zinc-800">
              <div className="text-[9px] uppercase tracking-widest text-zinc-500 font-bold mb-2">Encodage CNF</div>
              {[
                ['Variables', fmt(encoding.numVars)],
                ['Clauses', fmt(encoding.clauses.length)],
                ['Corps F_p', `${encoding.fieldBits} bits`],
                ['Bits inconnus', `${unknownBits}`],
              ].map(([k, v]) => (
                <div key={k} className="flex justify-between py-0.5">
                  <span className="text-[9px] text-zinc-500 uppercase tracking-wider">{k}</span>
                  <span className="text-[9px] text-cyan-400 font-mono">{v}</span>
                </div>
              ))}
            </div>
          )}

          {/* CDCL stats */}
          {stats && (
            <div className="p-4 border-b border-zinc-800">
              <div className="text-[9px] uppercase tracking-widest text-zinc-500 font-bold mb-2">CDCL Live</div>
              {[
                ['Décisions', fmt(stats.decisions)],
                ['Conflits', fmt(stats.conflicts)],
                ['Clauses apprises', fmt(stats.learnedClauses)],
                ['Backjumps', fmt(stats.backjumps)],
                ['Max backjump', `${stats.maxBackjump}`],
              ].map(([k, v]) => (
                <div key={k} className="flex justify-between py-0.5">
                  <span className="text-[9px] text-zinc-500 uppercase tracking-wider">{k}</span>
                  <span className="text-[9px] text-amber-400 font-mono">{v}</span>
                </div>
              ))}
              {encoding && (
                <div className="mt-2">
                  <div className="h-1 bg-zinc-800 rounded overflow-hidden">
                    <motion.div
                      className="h-full bg-cyan-500"
                      animate={{ width: `${Math.min(100, (assigned / encoding.numVars) * 100).toFixed(1)}%` }}
                      transition={{ duration: 0.3 }}
                    />
                  </div>
                  <div className="text-[8px] text-zinc-600 mt-1 font-mono">
                    {assigned}/{encoding.numVars} vars assignées
                  </div>
                </div>
              )}
            </div>
          )}

          {/* secp256k1 projection */}
          <div className="p-4">
            <div className="text-[9px] uppercase tracking-widest text-zinc-500 font-bold mb-2">
              Projection secp256k1
            </div>
            <div className="text-[8px] text-zinc-600 mb-2 uppercase tracking-wider">
              Pour {unknownBits} bits inconnus sur la vraie courbe Bitcoin
            </div>
            {[
              ['Corps', '256 bits'],
              ['Muls/add EC', fmt(secp256kStats.fieldMulsPerAdd)],
              ['Clauses/add', fmt(secp256kStats.clausesPerAdd)],
              ['Total clauses', fmt(secp256kStats.totalClauses)],
              ['Total vars', fmt(secp256kStats.totalVars)],
            ].map(([k, v]) => (
              <div key={k} className="flex justify-between py-0.5">
                <span className="text-[9px] text-zinc-600 uppercase tracking-wider">{k}</span>
                <span className="text-[9px] text-purple-400 font-mono">{v}</span>
              </div>
            ))}
            <div className="mt-2 p-2 bg-zinc-900/50 rounded border border-zinc-800">
              <p className="text-[8px] text-zinc-500 leading-relaxed">
                L'optimisation clé: avec Kangaroo → réduction de l'espace à ~40 bits inconnus → CNF tractable avec CDCL spécialisé (CryptoMiniSat).
              </p>
            </div>
          </div>
        </div>

        {/* Main panel */}
        <div className="flex-1 flex flex-col overflow-hidden">
          {/* Puzzle info */}
          {puzzle && (
            <div className="border-b border-zinc-800 p-4 shrink-0">
              <div className="text-[9px] uppercase tracking-widest text-zinc-500 font-bold mb-2">
                Puzzle — Courbe y² = x³ + 7 mod {puzzle.p}
              </div>
              <div className="grid grid-cols-3 gap-3">
                <div className="bg-zinc-900/50 rounded p-3 border border-zinc-800">
                  <div className="text-[8px] uppercase tracking-wider text-zinc-500 mb-1">Générateur G</div>
                  <div className="font-mono text-[10px] text-cyan-400">
                    ({puzzle.Gx}, {puzzle.Gy})
                  </div>
                </div>
                <div className="bg-zinc-900/50 rounded p-3 border border-zinc-800">
                  <div className="text-[8px] uppercase tracking-wider text-zinc-500 mb-1">Clé publique P</div>
                  <div className="font-mono text-[10px] text-amber-400">
                    ({puzzle.Px}, {puzzle.Py})
                  </div>
                </div>
                <div className={`bg-zinc-900/50 rounded p-3 border ${phase === 'solved' ? 'border-emerald-500/40' : 'border-zinc-800'}`}>
                  <div className="text-[8px] uppercase tracking-wider text-zinc-500 mb-1">Clé privée k</div>
                  {phase === 'solved' && foundK !== null ? (
                    <div className="font-mono text-[10px] text-emerald-400">
                      {foundK.toString()} ✓
                    </div>
                  ) : (
                    <div className="font-mono text-[10px] text-zinc-600">
                      {'?'.repeat(unknownBits)} ({unknownBits} bits)
                    </div>
                  )}
                </div>
              </div>
            </div>
          )}

          {/* Encoding visualizer — clause density heatmap */}
          {encoding && (
            <div className="border-b border-zinc-800 p-4 shrink-0">
              <div className="text-[9px] uppercase tracking-widest text-zinc-500 font-bold mb-2">
                Structure CNF — densité des clauses
              </div>
              <div className="flex gap-0.5 flex-wrap">
                {Array.from({ length: Math.min(200, encoding.clauses.length) }, (_, i) => {
                  const clause = encoding.clauses[Math.floor(i * encoding.clauses.length / 200)];
                  const len = clause?.length ?? 1;
                  const intensity = len === 1 ? 'bg-red-500' : len === 2 ? 'bg-amber-500' : len === 3 ? 'bg-cyan-500' : 'bg-zinc-600';
                  return <div key={i} className={`w-2 h-2 rounded-sm ${intensity} opacity-70`} />;
                })}
              </div>
              <div className="flex gap-4 mt-2">
                {[['rouge', 'Unitaires (fixes)', 'bg-red-500'], ['ambre', 'Binaires', 'bg-amber-500'], ['cyan', 'Ternaires', 'bg-cyan-500'], ['zinc', 'Longues', 'bg-zinc-600']].map(([, lbl, cls]) => (
                  <div key={lbl} className="flex items-center gap-1">
                    <div className={`w-2 h-2 rounded-sm ${cls} opacity-70`} />
                    <span className="text-[8px] text-zinc-600 uppercase tracking-wider">{lbl}</span>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Log */}
          <div className="flex-1 overflow-y-auto p-4 font-mono">
            {!puzzle && phase === 'idle' && (
              <div className="flex flex-col items-center justify-center h-full text-center">
                <Lock className="w-8 h-8 text-zinc-700 mb-4" />
                <div className="text-xs text-zinc-500 uppercase tracking-wider mb-2">
                  EC Puzzle Solver
                </div>
                <div className="text-[10px] text-zinc-600 max-w-sm leading-relaxed">
                  Encode le problème du logarithme discret elliptique comme instance CNF.
                  Le solveur CDCL trouve k tel que k×G = P — la même logique que pour les puzzles Bitcoin.
                </div>
                <div className="mt-6 flex items-center gap-2 text-[9px] text-zinc-700 uppercase tracking-wider">
                  <ChevronRight className="w-3 h-3" />
                  Choisir un prime et générer un puzzle
                </div>
              </div>
            )}

            <AnimatePresence>
              {log.map((line, i) => (
                <motion.div
                  key={`${i}-${line}`}
                  initial={{ opacity: 0, x: -8 }}
                  animate={{ opacity: 1, x: 0 }}
                  className={`text-[10px] py-0.5 ${
                    line.startsWith('✓') ? 'text-emerald-400'
                    : line.startsWith('↺') ? 'text-amber-400'
                    : line.startsWith('✗') ? 'text-red-400'
                    : 'text-zinc-400'
                  }`}
                >
                  {line}
                </motion.div>
              ))}
            </AnimatePresence>

            {phase === 'solving' && (
              <motion.div
                animate={{ opacity: [0.4, 1, 0.4] }}
                transition={{ repeat: Infinity, duration: 1.2 }}
                className="text-[10px] text-cyan-400 mt-1"
              >
                ▸ CDCL en cours — recherche de k…
              </motion.div>
            )}

            {phase === 'solved' && foundK !== null && puzzle && (
              <motion.div
                initial={{ opacity: 0, scale: 0.95 }}
                animate={{ opacity: 1, scale: 1 }}
                className="mt-4 p-4 bg-emerald-900/20 border border-emerald-500/30 rounded"
              >
                <div className="text-[9px] uppercase tracking-widest text-emerald-500 font-bold mb-2">
                  Logarithme Discret Résolu
                </div>
                <div className="text-xs text-emerald-300 font-mono">
                  k = {foundK.toString()}
                </div>
                <div className="text-[9px] text-emerald-600 mt-1">
                  {puzzle.k === Number(foundK) ? '✓ Identique à la clé secrète' : `Vrai: ${puzzle.k}`}
                </div>
                <div className="mt-2 text-[8px] text-zinc-500 leading-relaxed">
                  C'est la même opération fondamentale que les puzzles Bitcoin — trouver k tel que k×G = P sur secp256k1.
                  La tracabilité algorithmique via CNF/CDCL est démontrée.
                </div>
              </motion.div>
            )}

            {phase === 'error' && (
              <div className="text-[10px] text-red-400 mt-2">
                Erreur — aucun point valide sur cette courbe. Essayer un autre prime.
              </div>
            )}
          </div>

          {/* Architecture roadmap */}
          <div className="border-t border-zinc-800 p-4 shrink-0">
            <div className="text-[9px] uppercase tracking-widest text-zinc-500 font-bold mb-2">
              Feuille de Route — Au-delà de Kangaroo
            </div>
            <div className="flex gap-3 flex-wrap">
              {[
                ['✓', 'Tseitin CNF', 'text-emerald-400'],
                ['✓', 'CDCL + VSIDS', 'text-emerald-400'],
                ['✓', 'Corps fini F_p', 'text-emerald-400'],
                ['✓', 'EC point add CNF', 'text-emerald-400'],
                ['◐', 'Two-watched lits', 'text-amber-400'],
                ['◐', 'Courbe complète 256-bit', 'text-amber-400'],
                ['○', 'Hybrid Kangaroo+SAT', 'text-zinc-500'],
                ['○', 'Clauses algébriques EC', 'text-zinc-500'],
              ].map(([icon, label, color]) => (
                <div key={label} className="flex items-center gap-1">
                  <span className={`text-[9px] ${color}`}>{icon}</span>
                  <span className="text-[9px] text-zinc-500 uppercase tracking-wider">{label}</span>
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
