import React, { useState, useRef, useEffect } from 'react';
import { Play, Square, RotateCcw, Cpu } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';

// ─── Reduced SHA-256 (configurable rounds) ───────────────────────────────────

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

const IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];

const rotr = (n: number, x: number) => ((x >>> n) | (x << (32 - n))) >>> 0;
const add32 = (...ns: number[]) => ns.reduce((a, b) => (a + b) >>> 0, 0);
const hex8 = (n: number) => (n >>> 0).toString(16).padStart(8, '0');

function reducedSHA256(msgBytes: Uint8Array, rounds: number): string {
  const l = msgBytes.length;
  const padLen = ((l + 9 + 63) & ~63);
  const buf = new Uint8Array(padLen);
  buf.set(msgBytes);
  buf[l] = 0x80;
  const dv = new DataView(buf.buffer);
  dv.setUint32(padLen - 4, (l * 8) >>> 0, false);

  let [H0, H1, H2, H3, H4, H5, H6, H7] = IV;

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
      const ch = (e & f) ^ (~e & g);
      const t1 = add32(h, S1, ch, SHA256_K[i], W[i]);
      const S0 = rotr(2, a) ^ rotr(13, a) ^ rotr(22, a);
      const maj = (a & b) ^ (a & c) ^ (b & c);
      const t2 = add32(S0, maj);
      h = g; g = f; f = e; e = add32(d, t1);
      d = c; c = b; b = a; a = add32(t1, t2);
    }
    H0 = add32(H0, a); H1 = add32(H1, b); H2 = add32(H2, c); H3 = add32(H3, d);
    H4 = add32(H4, e); H5 = add32(H5, f); H6 = add32(H6, g); H7 = add32(H7, h);
  }
  return [H0, H1, H2, H3, H4, H5, H6, H7].map(hex8).join('');
}

// ─── Search engine (brute-force over alphabet, async/cancelable) ──────────────

const ALPHABET = 'abcdefghijklmnopqrstuvwxyz0123456789';

interface SearchState {
  status: 'idle' | 'running' | 'found' | 'exhausted';
  candidate: string;
  checked: number;
  total: number;
  found: string | null;
  rate: number; // hashes/sec
}

function buildCandidates(maxLen: number): string[] {
  const results: string[] = [];
  function recurse(cur: string) {
    if (cur.length > 0) results.push(cur);
    if (cur.length < maxLen) {
      for (const c of ALPHABET) recurse(cur + c);
    }
  }
  recurse('');
  return results;
}

// Space sizes for reference
function spaceSize(maxLen: number): number {
  let total = 0;
  for (let l = 1; l <= maxLen; l++) total += Math.pow(ALPHABET.length, l);
  return total;
}

// ─── Component ────────────────────────────────────────────────────────────────

export function SATInverterVisualizer() {
  const [targetHash, setTargetHash] = useState('');
  const [rounds, setRounds] = useState(4);
  const [maxLen, setMaxLen] = useState(2);
  const [search, setSearch] = useState<SearchState>({
    status: 'idle', candidate: '', checked: 0, total: 0, found: null, rate: 0,
  });
  const [previewHash, setPreviewHash] = useState('');
  const [inputMsg, setInputMsg] = useState('');

  const abortRef = useRef(false);
  const rateRef = useRef({ count: 0, ts: Date.now() });

  // Preview: compute reduced hash of typed message
  useEffect(() => {
    if (inputMsg) {
      const h = reducedSHA256(new TextEncoder().encode(inputMsg), rounds);
      setPreviewHash(h);
    } else {
      setPreviewHash('');
    }
  }, [inputMsg, rounds]);

  const stop = () => { abortRef.current = true; };

  const reset = () => {
    abortRef.current = true;
    setSearch({ status: 'idle', candidate: '', checked: 0, total: 0, found: null, rate: 0 });
  };

  const startSearch = async () => {
    if (!targetHash || search.status === 'running') return;
    abortRef.current = false;
    rateRef.current = { count: 0, ts: Date.now() };

    const total = spaceSize(maxLen);
    setSearch({ status: 'running', candidate: '', checked: 0, total, found: null, rate: 0 });

    const encoder = new TextEncoder();
    let checked = 0;
    let lastUpdate = Date.now();
    let lastCount = 0;

    // Iterate through all strings up to maxLen
    const stack: string[] = [''];
    const toCheck: string[] = [];

    // Build list lazily via generator-like loop
    async function run() {
      for (let len = 1; len <= maxLen; len++) {
        await searchLength(len);
        if (abortRef.current) return;
      }
      if (!abortRef.current) {
        setSearch(s => ({ ...s, status: 'exhausted' }));
      }
    }

    async function searchLength(len: number) {
      // Generate all strings of exactly `len` chars and test them
      const indices = new Array(len).fill(0);
      while (true) {
        if (abortRef.current) return;

        const candidate = indices.map(i => ALPHABET[i]).join('');
        const hash = reducedSHA256(encoder.encode(candidate), rounds);
        checked++;

        if (hash === targetHash) {
          setSearch(s => ({ ...s, status: 'found', candidate, checked, found: candidate }));
          abortRef.current = true;
          return;
        }

        // Update UI every ~200ms
        const now = Date.now();
        if (now - lastUpdate > 200) {
          const elapsed = (now - lastUpdate) / 1000;
          const rate = Math.round((checked - lastCount) / elapsed);
          lastUpdate = now;
          lastCount = checked;
          const cap = candidate;
          setSearch(s => ({ ...s, candidate: cap, checked, rate }));
          await new Promise(r => setTimeout(r, 0)); // yield to UI
        }

        // Increment indices (odometer)
        let pos = len - 1;
        while (pos >= 0) {
          indices[pos]++;
          if (indices[pos] < ALPHABET.length) break;
          indices[pos] = 0;
          pos--;
        }
        if (pos < 0) return; // exhausted this length
      }
    }

    run();
  };

  const totalSpace = spaceSize(maxLen);
  const progress = search.total > 0 ? search.checked / search.total : 0;

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-4xl mx-auto space-y-10">

          {/* Header */}
          <header className="space-y-2">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase">
              SHA-256 Réduit
              <span className="text-[10px] text-zinc-500 font-mono align-top ml-3">RECHERCHE DE PRÉIMAGE</span>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-tighter max-w-3xl leading-relaxed">
              SHA-256 avec N rounds au lieu de 64. La réduction affaiblit la diffusion —
              le search space reste exponentiel mais les collisions émergent plus tôt.
              Approche honnête : brute-force sur alphabet connu, espace borné.
            </p>
          </header>

          {/* Config */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-6">
            <div className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold">Configuration</div>
            <div className="grid grid-cols-2 gap-6">
              <div className="space-y-2">
                <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Rounds SHA-256 actifs</label>
                <div className="flex items-center gap-3">
                  <input
                    type="range" min={1} max={64} value={rounds}
                    onChange={e => { setRounds(+e.target.value); reset(); }}
                    className="flex-1 accent-cyan-500"
                  />
                  <span className="font-mono text-cyan-400 text-sm w-8 text-right">{rounds}</span>
                </div>
                <div className="text-[9px] text-zinc-600 font-mono">
                  {rounds <= 4 && 'très faible — collisions fréquentes'}
                  {rounds > 4 && rounds <= 12 && 'faible — recherche faisable'}
                  {rounds > 12 && rounds <= 24 && 'moyen — lent mais possible'}
                  {rounds > 24 && 'fort — brute-force impraticable'}
                </div>
              </div>
              <div className="space-y-2">
                <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Longueur max (alphabet a-z 0-9)</label>
                <div className="flex items-center gap-3">
                  <input
                    type="range" min={1} max={4} value={maxLen}
                    onChange={e => { setMaxLen(+e.target.value); reset(); }}
                    className="flex-1 accent-purple-500"
                  />
                  <span className="font-mono text-purple-400 text-sm w-8 text-right">{maxLen}</span>
                </div>
                <div className="text-[9px] text-zinc-600 font-mono">
                  {totalSpace.toLocaleString()} candidats à tester
                </div>
              </div>
            </div>
          </div>

          {/* Hash generator */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-4">
            <div className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold">
              Générer un hash cible (SHA-256 à {rounds} rounds)
            </div>
            <div className="flex gap-4">
              <div className="flex-1 bg-zinc-900/50 border border-zinc-800 p-3 rounded-lg">
                <label className="text-[9px] text-zinc-500 uppercase tracking-widest block mb-1">Message</label>
                <input
                  type="text"
                  value={inputMsg}
                  onChange={e => setInputMsg(e.target.value)}
                  placeholder="ex: ab, z3, hello..."
                  className="w-full bg-transparent text-orange-400 font-mono text-sm focus:outline-none"
                />
              </div>
              <div className="flex-1 bg-zinc-900/50 border border-zinc-800 p-3 rounded-lg">
                <label className="text-[9px] text-zinc-500 uppercase tracking-widest block mb-1">Hash réduit ({rounds}R)</label>
                <div className="font-mono text-xs text-zinc-400 break-all">{previewHash || '—'}</div>
              </div>
            </div>
            {previewHash && (
              <button
                onClick={() => { setTargetHash(previewHash); setInputMsg(''); reset(); }}
                className="text-[9px] uppercase tracking-widest font-bold border border-zinc-700 px-3 py-1.5 rounded text-zinc-400 hover:text-white hover:bg-zinc-800 transition"
              >
                Utiliser ce hash comme cible →
              </button>
            )}
          </div>

          {/* Target + search */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-6">
            <div className="space-y-2">
              <label className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold block">Hash cible</label>
              <input
                type="text"
                value={targetHash}
                onChange={e => { setTargetHash(e.target.value.toLowerCase().trim()); reset(); }}
                placeholder="Colle un hash SHA-256 réduit ici..."
                className="w-full bg-zinc-900/50 border border-zinc-800 rounded-lg p-3 font-mono text-xs text-zinc-300 focus:outline-none focus:border-zinc-600"
              />
            </div>

            <div className="flex items-center gap-4">
              <button
                onClick={startSearch}
                disabled={!targetHash || search.status === 'running'}
                className="px-5 py-2.5 bg-zinc-100 text-black text-[10px] font-bold uppercase tracking-widest hover:bg-white rounded transition-colors disabled:opacity-40 disabled:cursor-not-allowed flex items-center gap-2"
              >
                <Play className="w-3 h-3" /> Lancer la recherche
              </button>
              {search.status === 'running' && (
                <button onClick={stop} className="px-4 py-2 border border-red-700 text-red-400 text-[10px] font-bold uppercase tracking-widest rounded hover:bg-red-950/30 transition flex items-center gap-2">
                  <Square className="w-3 h-3" /> Stop
                </button>
              )}
              <button onClick={reset} className="px-4 py-2 border border-zinc-800 text-zinc-500 text-[10px] font-bold uppercase tracking-widest rounded hover:text-zinc-300 transition flex items-center gap-2">
                <RotateCcw className="w-3 h-3" /> Reset
              </button>
            </div>

            {/* Progress */}
            {search.status !== 'idle' && (
              <div className="space-y-3">
                <div className="flex justify-between text-[9px] font-mono text-zinc-500 uppercase tracking-widest">
                  <span>
                    {search.status === 'running' && <span className="text-cyan-500 flex items-center gap-1"><Cpu className="w-2.5 h-2.5 animate-spin inline" /> Recherche en cours</span>}
                    {search.status === 'found' && <span className="text-emerald-400">Préimage trouvée</span>}
                    {search.status === 'exhausted' && <span className="text-red-400">Espace épuisé — aucun résultat</span>}
                  </span>
                  <span>{search.checked.toLocaleString()} / {search.total.toLocaleString()}</span>
                </div>
                <div className="h-1 bg-zinc-800 rounded-full overflow-hidden">
                  <motion.div
                    className={`h-full rounded-full ${search.status === 'found' ? 'bg-emerald-500' : search.status === 'exhausted' ? 'bg-red-500' : 'bg-cyan-500'}`}
                    style={{ width: `${Math.min(progress * 100, 100)}%` }}
                    transition={{ duration: 0.1 }}
                  />
                </div>
                {search.status === 'running' && (
                  <div className="flex justify-between text-[9px] font-mono text-zinc-600">
                    <span>Candidat actuel: <span className="text-zinc-400">{search.candidate || '...'}</span></span>
                    <span>{search.rate.toLocaleString()} H/s</span>
                  </div>
                )}
              </div>
            )}

            {/* Result */}
            <AnimatePresence>
              {search.status === 'found' && search.found && (
                <motion.div
                  initial={{ opacity: 0, y: 8 }}
                  animate={{ opacity: 1, y: 0 }}
                  className="bg-emerald-500/10 border border-emerald-500/30 rounded-lg p-6 space-y-3"
                >
                  <div className="text-[9px] text-emerald-500/70 uppercase tracking-widest font-bold">
                    Préimage trouvée — SHA-256 ({rounds} rounds)
                  </div>
                  <div className="font-mono text-2xl text-emerald-400">"{search.found}"</div>
                  <div className="text-[9px] text-zinc-500 font-mono">
                    {search.checked.toLocaleString()} candidats testés · {search.rate.toLocaleString()} H/s
                  </div>
                  <div className="text-[9px] text-zinc-600 uppercase tracking-widest pt-2 border-t border-emerald-500/10">
                    Note: résultat valide pour SHA-256 à {rounds} rounds. SHA-256 complet (64R) reste hors de portée.
                  </div>
                </motion.div>
              )}
              {search.status === 'exhausted' && (
                <motion.div
                  initial={{ opacity: 0, y: 8 }}
                  animate={{ opacity: 1, y: 0 }}
                  className="bg-red-500/10 border border-red-500/20 rounded-lg p-4"
                >
                  <div className="text-[9px] text-red-400 uppercase tracking-widest font-bold mb-1">
                    Espace épuisé sans résultat
                  </div>
                  <div className="text-[9px] text-zinc-500 font-mono">
                    Le préimage n'est pas dans l'alphabet a-z0-9 de longueur ≤ {maxLen}, ou le hash ne correspond pas à ce nombre de rounds.
                  </div>
                </motion.div>
              )}
            </AnimatePresence>
          </div>

          {/* Roadmap */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-4">
            <div className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold">Prochaines étapes vers l'inversion réelle</div>
            <div className="space-y-3">
              {[
                { label: 'Brute-force borné', desc: 'Espace connu, alphabet fixe — implémenté ici', done: true },
                { label: 'Encodage CNF/SAT', desc: 'Transformer SHA-256 réduit en système de clauses booléennes', done: false },
                { label: 'Propagation de contraintes', desc: 'Réduire l\'espace de recherche par unit propagation (DPLL)', done: false },
                { label: 'Analyse différentielle', desc: 'Étudier comment les différences d\'entrée se propagent à travers les rounds', done: false },
                { label: 'MILP sur rounds réduits', desc: 'Modèle linéaire mixte pour des rounds ≤ 24', done: false },
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
    </div>
  );
}
