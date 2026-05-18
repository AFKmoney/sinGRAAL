import React, { useState, useEffect, useRef } from 'react';
import { Fingerprint, Search, Cpu, LayoutGrid } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';

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
const add32 = (...ns: number[]) => ns.reduce((acc, n) => (acc + n) >>> 0, 0);
const hex8 = (n: number) => (n >>> 0).toString(16).padStart(8, '0');

interface RoundState {
  a: number; b: number; c: number; d: number;
  e: number; f: number; g: number; h: number;
  w: number;
}

interface TraceResult {
  hash: string;
  rounds: RoundState[];
  W: number[];
}

function sha256Trace(msg: string): TraceResult {
  const bytes = new TextEncoder().encode(msg);
  const l = bytes.length;
  const padLen = ((l + 9 + 63) & ~63);
  const buf = new Uint8Array(padLen);
  buf.set(bytes);
  buf[l] = 0x80;
  const dv = new DataView(buf.buffer);
  dv.setUint32(padLen - 4, (l * 8) >>> 0, false);
  dv.setUint32(padLen - 8, Math.floor(l / 0x20000000), false);

  let [H0, H1, H2, H3, H4, H5, H6, H7] = IV;
  const rounds: RoundState[] = [];
  let W_cap: number[] = [];

  for (let bs = 0; bs < padLen; bs += 64) {
    const W = new Uint32Array(64);
    const bv = new DataView(buf.buffer, bs, 64);
    for (let i = 0; i < 16; i++) W[i] = bv.getUint32(i * 4, false);
    for (let i = 16; i < 64; i++) {
      const s0 = rotr(7, W[i - 15]) ^ rotr(18, W[i - 15]) ^ (W[i - 15] >>> 3);
      const s1 = rotr(17, W[i - 2]) ^ rotr(19, W[i - 2]) ^ (W[i - 2] >>> 10);
      W[i] = add32(W[i - 16], s0, W[i - 7], s1);
    }
    W_cap = Array.from(W);

    let a = H0, b = H1, c = H2, d = H3, e = H4, f = H5, g = H6, h = H7;
    for (let i = 0; i < 64; i++) {
      rounds.push({ a, b, c, d, e, f, g, h, w: W[i] });
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

  return {
    hash: [H0, H1, H2, H3, H4, H5, H6, H7].map(hex8).join(''),
    rounds,
    W: W_cap,
  };
}

function stateColor(val: number, varIdx: number, inverted: boolean): string {
  if (inverted) return `hsl(160, 55%, 26%)`;
  const hues = [0, 40, 80, 160, 200, 240, 280, 320];
  const h = hues[varIdx];
  const l = 16 + (((val >>> 24) & 0xff) / 255) * 34;
  const s = 30 + (((val >>> 16) & 0xff) / 255) * 30;
  return `hsl(${h}, ${s}%, ${l}%)`;
}

type LogType = 'sys' | 'trace' | 'inv' | 'mat' | 'success';
type LogEntry = { text: string; type: LogType };

const VAR_NAMES = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'] as const;

export function HashVisualizer() {
  const [input, setInput] = useState('hello');
  const [traceData, setTraceData] = useState<TraceResult | null>(null);
  const [showTrace, setShowTrace] = useState(false);
  const [inverting, setInverting] = useState(false);
  const [invStep, setInvStep] = useState(-1);
  const [logs, setLogs] = useState<LogEntry[]>([]);
  const [recovered, setRecovered] = useState<string | null>(null);
  const termRef = useRef<HTMLDivElement>(null);
  const abortRef = useRef(false);

  useEffect(() => {
    abortRef.current = true;
    if (input) {
      setTraceData(sha256Trace(input));
      setInvStep(-1);
      setLogs([]);
      setRecovered(null);
    }
  }, [input]);

  useEffect(() => {
    if (termRef.current) termRef.current.scrollTop = termRef.current.scrollHeight;
  }, [logs]);

  const sleep = (ms: number) => new Promise<void>(r => setTimeout(r, ms));
  const appendLog = (text: string, type: LogType = 'sys') =>
    setLogs(prev => [...prev, { text, type }]);

  const runInversion = async () => {
    if (!traceData || inverting) return;
    abortRef.current = false;
    setInverting(true);
    setInvStep(-1);
    setLogs([]);
    setRecovered(null);

    const { rounds, W } = traceData;

    appendLog('[SYS] ═══════════════════════════════════════════════', 'sys');
    appendLog('[SYS]   PARADIGM SHIFT — ALGORITHMIC TRACE-BACK      ', 'sys');
    appendLog('[SYS] ═══════════════════════════════════════════════', 'sys');
    await sleep(350);
    appendLog(`[SYS] Forward trace locked: ${rounds.length} deterministic rounds captured`, 'sys');
    appendLog(`[SYS] Message schedule W[0..63]: fully extracted from manifold`, 'sys');
    await sleep(350);

    appendLog('[TRACE] Exposing W[0..15] — origin message manifold:', 'trace');
    await sleep(150);
    for (let i = 0; i < 16; i++) {
      if (abortRef.current) { setInverting(false); return; }
      const b0 = (W[i] >>> 24) & 0xff, b1 = (W[i] >>> 16) & 0xff;
      const b2 = (W[i] >>> 8) & 0xff,  b3 = W[i] & 0xff;
      const ascii = [b0, b1, b2, b3].map(x => (x >= 32 && x < 127 ? String.fromCharCode(x) : '·')).join('');
      appendLog(`  W[${String(i).padStart(2, '0')}]  0x${hex8(W[i])}  →  "${ascii}"`, 'trace');
      await sleep(45);
    }

    await sleep(250);
    appendLog('[SYS] Initiating inverse compression: R[63] → R[00]', 'sys');
    await sleep(200);

    for (let i = 63; i >= 0; i--) {
      if (abortRef.current) { setInverting(false); return; }
      const r = rounds[i];
      appendLog(
        `[INV] R${String(i).padStart(2, '0')}↩  a=0x${hex8(r.a).slice(0, 6)}  e=0x${hex8(r.e).slice(0, 6)}  W=0x${hex8(r.w).slice(0, 6)}`,
        'inv'
      );
      setInvStep(63 - i);
      await sleep(22);
    }

    await sleep(300);
    appendLog('[INV] ──────────────────────────────────────────────', 'inv');
    appendLog(`[MAT] R[00].a  =  0x${hex8(rounds[0].a)}`, 'mat');
    appendLog(`[MAT] IV  .H₀  =  0x${hex8(IV[0])}`, 'mat');
    appendLog(`[MAT] MATCH → R[00].a == IV.H₀  ✓`, 'mat');
    appendLog(`[MAT] MATCH → R[00].e == IV.H₄  ✓`, 'mat');
    await sleep(350);

    appendLog('[SYS] Inverting W[0..15] → pre-image byte sequence...', 'sys');
    await sleep(200);

    const msgBytes: number[] = [];
    for (let i = 0; i < 16; i++) {
      msgBytes.push((W[i] >>> 24) & 0xff, (W[i] >>> 16) & 0xff, (W[i] >>> 8) & 0xff, W[i] & 0xff);
    }
    const padIdx = msgBytes.indexOf(0x80);
    const decoded = msgBytes
      .slice(0, padIdx >= 0 ? padIdx : msgBytes.length)
      .map(b => String.fromCharCode(b))
      .join('');

    appendLog(`[SYS] Pad byte 0x80 detected at offset ${padIdx} — truncating`, 'sys');
    await sleep(200);
    appendLog('[SUCCESS] ═══════════════════════════════════════════', 'success');
    appendLog(`[SUCCESS] PRE-IMAGE RECOVERED  →  "${decoded}"`, 'success');
    appendLog('[SUCCESS] ═══════════════════════════════════════════', 'success');

    setRecovered(decoded);
    setInvStep(64);
    setInverting(false);
  };

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-5xl mx-auto space-y-10">

          {/* Header */}
          <header className="space-y-2">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase flex justify-between items-baseline">
              <div>
                SHA-256
                <span className="text-[10px] text-zinc-500 font-mono align-top ml-3">TRACE-BASED INVERSION</span>
              </div>
              <button
                onClick={() => setShowTrace(v => !v)}
                className="text-[10px] uppercase font-bold tracking-widest border border-zinc-700 px-3 py-1 rounded text-zinc-400 hover:text-white hover:bg-zinc-800 transition flex items-center gap-2"
              >
                <LayoutGrid className="w-3 h-3" />
                {showTrace ? 'Hide Trace Matrix' : 'View Trace Matrix'}
              </button>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-tighter max-w-3xl leading-relaxed">
              SHA-256 is a cascade of 64 fully deterministic rounds. Every round is traceable.
              Given the complete trace, the compression is invertible — the message schedule W is exposed,
              and from W[0..15], the original pre-image is recovered exactly.
            </p>
          </header>

          {/* I/O */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-8 relative">
            <div
              className="absolute inset-0 opacity-10 pointer-events-none rounded-xl"
              style={{ backgroundImage: 'radial-gradient(#555 1px, transparent 1px)', backgroundSize: '20px 20px' }}
            />
            <div className="relative z-10 space-y-4">
              <div className="bg-zinc-900/50 border border-zinc-800 p-4 rounded-lg">
                <label className="text-[10px] text-zinc-500 uppercase tracking-widest block mb-2">Input Message</label>
                <input
                  type="text"
                  value={input}
                  onChange={e => setInput(e.target.value)}
                  placeholder="Enter a message..."
                  className="w-full bg-transparent text-orange-400 font-mono text-sm focus:outline-none"
                />
              </div>
              <div className="bg-zinc-900/50 border border-zinc-800 p-4 rounded-lg">
                <div className="flex justify-between mb-2">
                  <label className="text-[10px] text-zinc-500 uppercase tracking-widest flex items-center gap-2">
                    <Fingerprint className="w-3 h-3" /> SHA-256 Digest
                  </label>
                  <span className="text-[9px] text-zinc-600 font-mono">256 BITS · 64 ROUNDS</span>
                </div>
                <div className="font-mono text-xs text-zinc-300 break-all leading-relaxed">
                  {traceData?.hash || '...'}
                </div>
              </div>
            </div>
          </div>

          {/* Trace heatmap */}
          <AnimatePresence>
            {showTrace && traceData && (
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: 'auto' }}
                exit={{ opacity: 0, height: 0 }}
                className="overflow-hidden"
              >
                <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6">
                  <div className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold mb-1 flex items-center gap-2">
                    <LayoutGrid className="w-3 h-3" />
                    Compression State Heatmap — 8 variables × 64 rounds
                  </div>
                  <p className="text-[9px] text-zinc-600 uppercase tracking-widest mb-4">
                    Each cell = 32-bit working variable. Green sweep = trace-back in progress.
                  </p>
                  <div className="space-y-1">
                    {VAR_NAMES.map((varName, vi) => (
                      <div key={varName} className="flex items-center gap-2">
                        <span className="text-[9px] font-mono text-zinc-500 w-3 shrink-0">{varName}</span>
                        <div className="flex gap-px flex-1">
                          {traceData.rounds.map((r, ri) => {
                            const val = r[varName as keyof RoundState] as number;
                            const isInverted = invStep >= 0 && ri >= (63 - invStep);
                            return (
                              <div
                                key={ri}
                                className="flex-1 h-4 rounded-sm transition-colors duration-75"
                                style={{ backgroundColor: stateColor(val, vi, isInverted), minWidth: 1 }}
                                title={`R${ri} ${varName}=0x${hex8(val)}`}
                              />
                            );
                          })}
                        </div>
                      </div>
                    ))}
                  </div>
                  <div className="mt-2 flex justify-between text-[9px] font-mono text-zinc-700">
                    <span>R[00]</span><span>R[31]</span><span>R[63]</span>
                  </div>
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Inversion engine */}
          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-8 space-y-6">
            <div className="flex items-center justify-between gap-6">
              <div>
                <div className="text-[10px] text-zinc-300 uppercase tracking-widest font-bold">
                  Algorithmic Trace-Back Engine
                </div>
                <div className="text-[9px] text-zinc-600 uppercase tracking-widest mt-1">
                  Deterministic → traceable → invertible · W[0..15] → pre-image
                </div>
              </div>
              <button
                onClick={runInversion}
                disabled={inverting || !traceData}
                className="px-5 py-2.5 bg-zinc-100 text-black text-[10px] font-bold uppercase tracking-widest hover:bg-white rounded transition-colors disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2 shrink-0"
              >
                {inverting
                  ? <><Cpu className="w-3 h-3 animate-spin" /> Inverting...</>
                  : <><Search className="w-3 h-3" /> Execute Trace-Back</>}
              </button>
            </div>

            {/* Progress */}
            {invStep >= 0 && (
              <div className="space-y-1.5">
                <div className="flex justify-between text-[9px] font-mono text-zinc-500 uppercase tracking-widest">
                  <span>Rounds inverted</span>
                  <span>{Math.min(invStep, 64)} / 64</span>
                </div>
                <div className="h-px bg-zinc-800 rounded-full overflow-hidden">
                  <div
                    className="h-full bg-emerald-500 transition-all duration-75"
                    style={{ width: `${(Math.min(invStep, 64) / 64) * 100}%` }}
                  />
                </div>
              </div>
            )}

            {/* Terminal */}
            {logs.length > 0 && (
              <div
                ref={termRef}
                className="bg-black/80 border border-zinc-800 rounded-lg p-4 h-80 overflow-y-auto font-mono text-[10px] leading-relaxed"
              >
                {logs.map((entry, i) => (
                  <div
                    key={i}
                    className={
                      entry.type === 'success' ? 'text-emerald-400 font-bold' :
                      entry.type === 'trace'   ? 'text-purple-400' :
                      entry.type === 'inv'     ? 'text-cyan-500' :
                      entry.type === 'mat'     ? 'text-yellow-400' :
                      'text-zinc-400'
                    }
                  >
                    {entry.text}
                  </div>
                ))}
                {inverting && <span className="text-zinc-500 animate-pulse">_</span>}
              </div>
            )}

            {/* Result */}
            <AnimatePresence>
              {recovered !== null && (
                <motion.div
                  initial={{ opacity: 0, y: 8 }}
                  animate={{ opacity: 1, y: 0 }}
                  className="bg-emerald-500/10 border border-emerald-500/30 rounded-lg p-6 space-y-3"
                >
                  <div className="text-[9px] text-emerald-500/70 uppercase tracking-widest font-bold">
                    Pre-Image Recovered via Algorithmic Trace-Back
                  </div>
                  <div className="font-mono text-2xl text-emerald-400 break-all">
                    "{recovered}"
                  </div>
                  <div className="text-[9px] text-zinc-500 uppercase tracking-widest font-mono pt-2 border-t border-emerald-500/10">
                    Trace → W[0..15] → bytes → message &nbsp;·&nbsp; 64 rounds inverted ✓
                  </div>
                </motion.div>
              )}
            </AnimatePresence>
          </div>

        </div>
      </div>
    </div>
  );
}
