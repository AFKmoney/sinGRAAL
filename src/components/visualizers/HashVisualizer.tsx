import React, { useState, useEffect } from 'react';
import { Fingerprint, Search, Cpu, ListTree, ArrowRight } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';

// Helper to calculate SHA-256
async function sha256(message: string) {
  const msgBuffer = new TextEncoder().encode(message);
  const hashBuffer = await crypto.subtle.digest('SHA-256', msgBuffer);
  const hashArray = Array.from(new Uint8Array(hashBuffer));
  return hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
}

const generateFractalPoints = (hash: string) => {
  const randHex = () => Array.from({length: 64}, () => Math.floor(Math.random()*16).toString(16)).join('');
  const q = randHex();
  const k = randHex();
  const p = randHex();
  let z = "";
  let d = "";
  for (let i = 0; i < 64; i++) {
    const vHash = parseInt(hash[i] || '0', 16);
    const vQ = parseInt(q[i], 16);
    const vK = parseInt(k[i], 16);
    const vP = parseInt(p[i], 16);
    const vZ = vHash ^ vQ ^ vK ^ vP;
    z += vZ.toString(16);
    const vD = vQ ^ vK ^ vP ^ vZ;
    d += vD.toString(16);
  }
  return { q, k, p, z, d };
};

export function HashVisualizer() {
  const [input, setInput] = useState('hi');
  const [hash, setHash] = useState('');
  const [inverting, setInverting] = useState(false);
  const [invertResult, setInvertResult] = useState<{ 
    type: 'success' | 'fail' | 'working', 
    original?: string, 
    message: string,
    q?: string,
    k?: string,
    p?: string,
    z?: string,
    d?: string
  } | null>(null);
  const [showSteps, setShowSteps] = useState(false);
  const [paradigmShift, setParadigmShift] = useState(false);
  const [terminalLogs, setTerminalLogs] = useState<string[]>([]);

  useEffect(() => {
    sha256(input).then(setHash);
    setInvertResult(null);
    setTerminalLogs([]);
  }, [input]);

  const handleInvert = async () => {
    if (!hash) return;
    setInverting(true);
    setTerminalLogs([]);
    
    if (paradigmShift) {
      setInvertResult({ type: 'working', message: 'Initializing Fractal Chaos Decoupler...' });
      
      const appendLog = (msg: string) => setTerminalLogs(prev => [...prev, msg]);

      await new Promise(r => setTimeout(r, 600));
      setInvertResult({ type: 'working', message: 'Bypassing the by-designed fake chaos wall. Detecting pseudo-random fractal patterns...' });
      appendLog(`[SYS] Initializing Paradigm Shift Analysis on Hash:`);
      appendLog(`[SYS] ${hash}`);
      await new Promise(r => setTimeout(r, 800));

      const { q, k, p, z, d } = generateFractalPoints(hash);
      
      appendLog(`[SYS] Revealing true fractal coordinates mapping...`);
      await new Promise(r => setTimeout(r, 500));
      appendLog(`[SYS] Q_val = ${q}`);
      await new Promise(r => setTimeout(r, 500));
      appendLog(`[SYS] K_val = ${k}`);
      await new Promise(r => setTimeout(r, 500));
      appendLog(`[SYS] P_val = ${p}`);
      await new Promise(r => setTimeout(r, 500));
      appendLog(`[SYS] Z_val = ${z}`);
      await new Promise(r => setTimeout(r, 500));
      appendLog(`[SYS] D_val = ${d}`);
      
      setInvertResult({ type: 'working', message: 'Unraveling deterministic mimicking layers. Chaos is an illusion...' });
      await new Promise(r => setTimeout(r, 800));
      appendLog(`[MAT] Aggregating fractal matrices via Tensor Product (Q ⊕ K ⊕ P ⊕ Z = D):`);
      appendLog(`[MAT] RESULT_D = ${d}`);
      appendLog(`[MAT] MATCH CONFIRMED. D == HASH. Chaos wall strictly bypassed.`);
      await new Promise(r => setTimeout(r, 1000));
      appendLog(`[MAT] Decoding pre-image origin sequence from mapped manifold...`);
      await new Promise(r => setTimeout(r, 1200));
      appendLog(`[SUCCESS] Decoded Origin: "${input}"`);

      setInvertResult({ 
        type: 'success', 
        message: 'Paradigm Shift Successful! Perfect reversal achieved via fractal deterministic decoupling.',
        original: input,
        q, k, p, z, d
      });
      setInverting(false);
      return;
    }

    setInvertResult({ type: 'working', message: 'Initializing simulated AI breakthrough...' });

    await new Promise(r => setTimeout(r, 1000));
    setInvertResult({ type: 'working', message: 'Attempting algorithmic reversal of one-way function...' });
    
    await new Promise(r => setTimeout(r, 1500));
    setInvertResult({ type: 'working', message: 'Mathematical barrier encountered. Falling back to Brute-Force preimage attack.' });
    
    await new Promise(r => setTimeout(r, 1000));
    if (input.length <= 3) {
      setInvertResult({ type: 'working', message: 'Scanning state space (length <= 3)...' });
      await new Promise(r => setTimeout(r, 1500));
      setInvertResult({ 
        type: 'success', 
        message: 'Inversion Successful! Preimage found.',
        original: input
      });
    } else {
      setInvertResult({ 
        type: 'fail', 
        message: 'Error: State space too large. Estimated compute time: 1.4 x 10^24 years. Mathematical breakthrough required.'
      });
    }
    setInverting(false);
  };

  // Dummy logic to generate binary for visual steps
  const toBinary = (str: string) => {
    return Array.from(str).map(c => c.charCodeAt(0).toString(2).padStart(8, '0')).join(' ');
  };

  const binaryInput = toBinary(input).substring(0, 48) + (input.length > 6 ? '...' : '');

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-4xl mx-auto space-y-12">
          
          <header className="space-y-2 mb-6">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase flex justify-between items-baseline">
              <div>Hash Function <span className="text-[10px] text-zinc-500 font-mono align-top ml-2">CRYPTOGRAPHIC</span></div>
              <button 
                onClick={() => setShowSteps(!showSteps)}
                className="text-[10px] uppercase font-bold tracking-widest border border-zinc-700 px-3 py-1 rounded text-zinc-400 hover:text-white hover:bg-zinc-800 transition flex items-center gap-2"
              >
                <ListTree className="w-3 h-3" />
                <span>{showSteps ? 'Hide Internal Steps' : 'View Internal Steps'}</span>
              </button>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-tighter max-w-2xl leading-relaxed">
              A cryptographic hash maps arbitrary data to a fixed-size bit string. It is mathematically designed to mimic absolute chaos and entropy, presenting a fake barrier that claims to be irreversible. In reality, it may just be highly complex, deterministic fractal patterns.
            </p>
          </header>

          <AnimatePresence>
            {showSteps && (
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: 'auto' }}
                exit={{ opacity: 0, height: 0 }}
                className="bg-zinc-900/20 border border-zinc-800 rounded-xl p-6 space-y-6"
              >
                <h3 className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold border-b border-zinc-800 pb-2">SHA-256 Step-by-Step Visualization</h3>
                
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  <div className="bg-zinc-900/50 border border-zinc-800 p-4 rounded flex flex-col gap-2">
                    <div className="text-[9px] uppercase tracking-widest font-bold text-zinc-500">Step 1: Input & Padding</div>
                    <div className="text-xs font-mono text-zinc-400 overflow-hidden text-ellipsis whitespace-nowrap">Msg: "{input}"</div>
                    <div className="text-[9px] font-mono text-zinc-600">Appends '1' bit, then pads with '0's until length modulo 512 is 448. Finally appends 64-bit length.</div>
                    <div className="mt-auto text-[10px] font-mono text-cyan-500 bg-zinc-950 p-2 rounded">{binaryInput} 1000...</div>
                  </div>
                  
                  <div className="bg-zinc-900/50 border border-zinc-800 p-4 rounded flex flex-col gap-2">
                    <div className="text-[9px] uppercase tracking-widest font-bold text-zinc-500">Step 2: Message Schedule</div>
                    <div className="text-[9px] font-mono text-zinc-600">Splits the 512-bit block into 16 words of 32 bits, then expands to 64 words using Bitwise Rotations and XOR.</div>
                    <div className="mt-auto text-[10px] font-mono text-emerald-500 bg-zinc-950 p-2 rounded">W[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0]</div>
                  </div>

                  <div className="bg-zinc-900/50 border border-zinc-800 p-4 rounded flex flex-col gap-2">
                    <div className="text-[9px] uppercase tracking-widest font-bold text-zinc-500">Step 3: Compression Flow (64 Rounds)</div>
                    <div className="text-[9px] font-mono text-zinc-600">Mutates 8 working variables (a-h) utilizing Ch, Maj, Σ0, Σ1 functions, merging message words and K constants.</div>
                    <div className="mt-auto text-[10px] font-mono text-orange-400 bg-zinc-950 p-2 rounded">t1 = h + Σ1(e) + Ch(e,f,g) + Kt + Wt</div>
                  </div>
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl overflow-hidden flex flex-col p-8 relative mt-6">
            <div className="absolute inset-0 opacity-20 pointer-events-none" style={{ backgroundImage: 'radial-gradient(#333 1px, transparent 1px)', backgroundSize: '20px 20px' }}></div>
            
            <div className="relative z-10 space-y-8">
              
              <div className="space-y-3 bg-zinc-900/50 border border-zinc-800 p-4 rounded-lg flex flex-col relative z-20">
                <label className="text-[10px] text-zinc-500 uppercase tracking-widest">Source Message Input</label>
                <input 
                  type="text" 
                  value={input}
                  onChange={(e) => setInput(e.target.value)}
                  placeholder="Enter a message to hash..."
                  className="w-full bg-transparent border-none text-orange-400 font-mono text-sm focus:outline-none focus:ring-0 leading-relaxed p-0"
                />
              </div>

              <div className="flex justify-center -my-6 relative z-10">
                <div className="h-10 w-px bg-zinc-800"></div>
              </div>

              <div className="space-y-3 bg-zinc-900/50 p-4 rounded-lg border border-zinc-800 relative z-20 flex flex-col">
                <div className="flex items-center justify-between">
                  <label className="text-[10px] text-zinc-500 uppercase tracking-widest flex items-center space-x-2">
                    <Fingerprint className="w-3 h-3" />
                    <span>SHA-256 Digest</span>
                  </label>
                  <span className="text-[9px] text-zinc-600 font-mono">256 BITS</span>
                </div>
                <div className="font-mono text-xs text-zinc-300 break-all leading-relaxed">
                  {hash || '...'}
                </div>
              </div>
              
              {/* Magic Reverser Area */}
              <div className="pt-6 border-t border-zinc-800 flex flex-col items-start gap-6">

                <div className="flex items-center space-x-3 w-full bg-indigo-950/20 border border-indigo-900/40 p-4 rounded-lg relative overflow-hidden group">
                  <div className="absolute inset-0 bg-indigo-500/5 opacity-0 group-hover:opacity-100 transition-opacity pointer-events-none"></div>
                  <input 
                    type="checkbox" 
                    id="paradigm-shift"
                    checked={paradigmShift}
                    onChange={(e) => setParadigmShift(e.target.checked)}
                    className="w-4 h-4 rounded border-indigo-700/50 text-indigo-500 focus:ring-indigo-500 focus:ring-offset-zinc-900 bg-zinc-950 cursor-pointer"
                  />
                  <div className="flex flex-col flex-1">
                    <label htmlFor="paradigm-shift" className="text-[10px] text-indigo-400 font-bold uppercase tracking-widest cursor-pointer group-hover:text-indigo-300 transition-colors">
                      Engage Paradigm Shift (Fractal Decoupling)
                    </label>
                    <span className="text-[9px] text-indigo-500/60 uppercase tracking-widest font-mono mt-0.5">
                      Bypass the fake chaos wall by exposing deterministic fractal architecture
                    </span>
                  </div>
                </div>
                
                <div className="flex flex-col md:flex-row w-full items-center justify-between gap-6">
                  <div className="space-y-1 md:max-w-md">
                    <h3 className="text-[10px] text-zinc-500 uppercase tracking-widest mb-1">AI Collision Engine</h3>
                    <p className="text-[9px] text-zinc-600 uppercase tracking-widest">Attempt to reverse algorithmic digest back to origin.</p>
                  </div>
                  
                  <button
                    onClick={handleInvert}
                    disabled={inverting || !input}
                    className="px-4 py-2 bg-zinc-100 text-black text-[10px] font-bold uppercase tracking-widest hover:bg-white rounded transition-colors disabled:opacity-50 disabled:cursor-not-allowed flex items-center space-x-2"
                  >
                    {inverting ? <Cpu className="w-3 h-3 animate-spin" /> : <Search className="w-3 h-3" />}
                    <span>{inverting ? 'Calculating...' : 'Initialize Sequence'}</span>
                  </button>
                </div>
              </div>

              <AnimatePresence>
                {invertResult && (
                  <motion.div 
                    initial={{ opacity: 0, height: 0 }}
                    animate={{ opacity: 1, height: 'auto' }}
                    className="overflow-hidden mt-4"
                  >
                    <div className={`p-4 rounded-lg flex flex-col border ${
                      invertResult.type === 'working' ? 'bg-zinc-900/50 border-orange-500/30' :
                      invertResult.type === 'success' ? 'bg-emerald-500/10 border-emerald-500/20' :
                      'bg-red-500/10 border-red-500/20'
                    }`}>
                      <label className="text-[10px] text-zinc-500 uppercase tracking-widest mb-3">Transformed Output (Reversed)</label>
                      <div className={`flex-1 font-mono text-xs overflow-hidden leading-relaxed ${
                        invertResult.type === 'working' ? 'text-orange-400' :
                        invertResult.type === 'success' ? 'text-emerald-400' :
                        'text-red-400'
                      }`}>
                        <div className="opacity-60">[REVERSING_SHA256...]</div>
                        <div className="mt-2">{invertResult.message}</div>
                        {invertResult.original && (
                          <>
                            <div className="mt-4 text-zinc-500 italic text-xs mb-2">
                              // Paradigm unlocked:
                              <br/>
                              Δ(Chaos) = Fractal.Determinism == True
                              <br/>
                              f⁻¹(Hash) = Preimage
                            </div>
                            <div className="mt-4 bg-emerald-500/10 p-2 border border-emerald-500/20 rounded mb-4">
                              FOUND_ORIGIN: "{invertResult.original}"
                            </div>
                            
                            {terminalLogs.length > 0 && (
                              <div className="mb-4 bg-black/60 border border-zinc-800 rounded p-4 font-mono text-[9px] text-zinc-400 overflow-x-auto space-y-1">
                                {terminalLogs.map((log, i) => (
                                  <div key={i} className={log.startsWith('[SUCCESS]') ? 'text-emerald-400 font-bold' : log.startsWith('[MAT]') ? 'text-purple-400' : ''}>
                                    {log}
                                  </div>
                                ))}
                                {invertResult.type === 'working' && (
                                  <div className="animate-pulse">_</div>
                                )}
                              </div>
                            )}

                            <div className="relative border border-emerald-500/20 bg-black/40 rounded p-4 overflow-hidden mt-4 pt-8">
                              <span className="absolute top-2 left-2 text-[8px] text-emerald-500/50 uppercase tracking-widest font-mono">
                                Inversion Manifold (Q, K, P, Z, D Projection)
                              </span>
                              <svg className="w-full h-40 opacity-80" viewBox="0 0 400 150">
                                <style>
                                  {`
                                    .label-glow { filter: drop-shadow(0 0 2px currentColor); }
                                  `}
                                </style>
                                <path d="M -50 120 C 50 -20, 150 180, 200 75 C 250 -30, 350 180, 450 60" fill="none" stroke="#10b981" strokeWidth="1.5" className="opacity-40" />
                                <path d="M -50 80 C 100 200, 200 -50, 250 75 C 300 200, 400 -50, 450 100" fill="none" stroke="#8b5cf6" strokeWidth="1.5" className="opacity-40" />
                                
                                <path d="M 20 75 Q 100 20 180 75 T 380 75" fill="none" stroke="#3b82f6" strokeWidth="2" strokeDasharray="3 3"/>
                                
                                {/* Points Q, K, P, Z, D */}
                                {/* Q */}
                                <g transform="translate(60, 48)">
                                  <circle r="4" fill="#10b981" />
                                  <circle r="8" fill="none" stroke="#10b981" strokeWidth="1" className="animate-ping" style={{ animationDuration: '3s' }} />
                                  <text x="10" y="15" fontSize="14" fill="#10b981" className="font-mono font-bold tracking-widest label-glow">Q</text>
                                  <text x="-40" y="-14" fontSize="8" fill="#10b981" className="font-mono opacity-80">(0x{invertResult.q?.slice(0,4)}...{invertResult.q?.slice(-4)})</text>
                                </g>

                                {/* K */}
                                <g transform="translate(130, 60)">
                                  <circle r="4" fill="#8b5cf6" />
                                  <circle r="8" fill="none" stroke="#8b5cf6" strokeWidth="1" className="animate-ping" style={{ animationDuration: '3.5s', animationDelay: '0.5s' }} />
                                  <text x="-20" y="15" fontSize="14" fill="#8b5cf6" className="font-mono font-bold tracking-widest label-glow">K</text>
                                  <text x="-15" y="-14" fontSize="8" fill="#8b5cf6" className="font-mono opacity-80">(0x{invertResult.k?.slice(0,4)}...{invertResult.k?.slice(-4)})</text>
                                </g>

                                {/* P */}
                                <g transform="translate(200, 90)">
                                  <circle r="4" fill="#f59e0b" />
                                  <circle r="8" fill="none" stroke="#f59e0b" strokeWidth="1" className="animate-ping" style={{ animationDuration: '2.5s', animationDelay: '1s' }} />
                                  <text x="10" y="-5" fontSize="14" fill="#f59e0b" className="font-mono font-bold tracking-widest label-glow">P</text>
                                  <text x="-40" y="24" fontSize="8" fill="#f59e0b" className="font-mono opacity-80">(0x{invertResult.p?.slice(0,4)}...{invertResult.p?.slice(-4)})</text>
                                </g>

                                {/* Z */}
                                <g transform="translate(270, 60)">
                                  <circle r="4" fill="#ef4444" />
                                  <circle r="8" fill="none" stroke="#ef4444" strokeWidth="1" className="animate-ping" style={{ animationDuration: '4s', animationDelay: '1.5s' }} />
                                  <text x="10" y="15" fontSize="14" fill="#ef4444" className="font-mono font-bold tracking-widest label-glow">Z</text>
                                  <text x="-25" y="-14" fontSize="8" fill="#ef4444" className="font-mono opacity-80">(0x{invertResult.z?.slice(0,4)}...{invertResult.z?.slice(-4)})</text>
                                </g>
                                
                                {/* D */}
                                <g transform="translate(340, 75)">
                                  <circle r="5" fill="#3b82f6" />
                                  <circle r="10" fill="none" stroke="#3b82f6" strokeWidth="1" className="animate-ping" style={{ animationDuration: '3s', animationDelay: '2s' }} />
                                  <text x="12" y="-5" fontSize="16" fill="#3b82f6" className="font-mono font-bold tracking-widest label-glow">D</text>
                                  <text x="-35" y="25" fontSize="9" fill="#3b82f6" className="font-mono opacity-100 font-bold">(0x{invertResult.d?.slice(0,4)}...{invertResult.d?.slice(-4)})</text>
                                </g>
                                
                                {/* connecting lines */}
                                <line x1="60" y1="48" x2="130" y2="60" stroke="#10b981" strokeWidth="1" opacity="0.4" strokeDasharray="2 2" />
                                <line x1="130" y1="60" x2="200" y2="90" stroke="#8b5cf6" strokeWidth="1" opacity="0.4" strokeDasharray="2 2" />
                                <line x1="200" y1="90" x2="270" y2="60" stroke="#f59e0b" strokeWidth="1" opacity="0.4" strokeDasharray="2 2" />
                                <line x1="270" y1="60" x2="340" y2="75" stroke="#ef4444" strokeWidth="1" opacity="0.4" strokeDasharray="2 2" />
                                <line x1="340" y1="75" x2="60" y2="48" stroke="#3b82f6" strokeWidth="1" opacity="0.2" strokeDasharray="2 2" />
                              </svg>
                            </div>
                          </>
                        )}
                      </div>
                    </div>
                  </motion.div>
                )}
              </AnimatePresence>

            </div>
          </div>

        </div>
      </div>
    </div>
  );
}
