import React, { useState } from 'react';
import { RefreshCw, AlertTriangle } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';

export function OrVisualizer() {
  const [inputA, setInputA] = useState('1101');
  const [inputB, setInputB] = useState('1011');
  const [attemptInvert, setAttemptInvert] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleInput = (val: string, setter: (v: string) => void) => {
    if (/[^01]/.test(val)) {
      setError('Only binary values (0, 1) are permitted.');
    } else {
      setError(null);
    }
    setter(val.replace(/[^01]/g, ''));
    setAttemptInvert(false);
  };

  const formatBits = (str: string) => str.replace(/[^01]/g, '').slice(0, 8);
  const aBits = formatBits(inputA).padEnd(4, '0');
  const bBits = formatBits(inputB).padEnd(4, '0');

  const outBits = Array.from({ length: Math.max(aBits.length, bBits.length) }, (_, i) => {
    return (aBits[i] === '1' || bBits[i] === '1') ? '1' : '0';
  }).join('');

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-4xl mx-auto space-y-12">
          
          <header className="space-y-4">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase flex items-baseline">
              OR_Transform <span className="text-[10px] text-zinc-500 font-mono align-top ml-3 tracking-widest font-bold">LOSSY OPERATION</span>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-widest max-w-2xl leading-relaxed">
              The OR gate destroys information. If the output is "1", we don't know whether A was 1, B was 1, or both were 1. Reversing this operation yields multiple possibilities.
            </p>
          </header>

          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl overflow-hidden flex flex-col p-8 relative">
            <div className="absolute inset-0 opacity-20 pointer-events-none" style={{ backgroundImage: 'radial-gradient(#333 1px, transparent 1px)', backgroundSize: '20px 20px' }}></div>
            
            <div className="relative z-10 space-y-8">
              
              <div className="grid grid-cols-1 md:grid-cols-3 gap-6 items-center">
                
                {/* Inputs block */}
                <div className="space-y-4 relative">
                  <div className="bg-zinc-900/50 border border-zinc-800 rounded-lg p-4">
                    <label className="text-[10px] text-zinc-500 uppercase tracking-widest block mb-2 font-bold">Source Element (A)</label>
                    <input 
                      type="text" 
                      value={inputA}
                      onChange={(e) => handleInput(e.target.value, setInputA)}
                      className="w-full bg-transparent border-none font-mono text-sm text-rose-400 tracking-[0.5em] focus:outline-none focus:ring-0 p-0"
                      maxLength={8}
                    />
                  </div>
                  <div className="bg-zinc-900/50 border border-zinc-800 rounded-lg p-4">
                    <label className="text-[10px] text-zinc-500 uppercase tracking-widest block mb-2 font-bold">Source Element (B)</label>
                    <input 
                      type="text" 
                      value={inputB}
                      onChange={(e) => handleInput(e.target.value, setInputB)}
                      className="w-full bg-transparent border-none font-mono text-sm text-rose-400 tracking-[0.5em] focus:outline-none focus:ring-0 p-0"
                      maxLength={8}
                    />
                  </div>
                  {error && (
                    <motion.div 
                      initial={{ opacity: 0, y: -10 }}
                      animate={{ opacity: 1, y: 0 }}
                      className="absolute -top-12 left-0 right-0 bg-red-950/80 border border-red-900/50 rounded p-2 z-20 backdrop-blur-sm"
                    >
                      <p className="text-[9px] text-red-500 font-bold uppercase tracking-widest text-center">{error}</p>
                    </motion.div>
                  )}
                </div>

                {/* Operation block */}
                <div className="flex flex-col items-center justify-center space-y-4">
                  <div className="w-12 h-12 bg-zinc-900/50 border border-zinc-800 rounded-lg flex items-center justify-center text-zinc-400 font-mono shadow-[0_0_10px_rgba(255,255,255,0.05)]">
                    OR
                  </div>
                  <p className="text-[9px] font-bold uppercase tracking-widest text-zinc-600">
                    MATRIX TRANSLATION
                  </p>
                </div>

                {/* Output block */}
                <div className="relative h-full flex flex-col justify-center">
                  <div className="bg-zinc-900/50 border border-zinc-800 rounded-lg p-4 flex flex-col justify-center min-h-[104px]">
                    <label className="text-[10px] text-zinc-500 uppercase tracking-widest block mb-2 font-bold">
                      Transformed Element
                    </label>
                    <div className="w-full bg-transparent font-mono text-sm text-amber-400 tracking-[0.5em]">
                      {outBits}
                    </div>
                  </div>

                  <AnimatePresence>
                    {attemptInvert && (
                      <motion.div 
                        initial={{ opacity: 0, scale: 0.95 }}
                        animate={{ opacity: 1, scale: 1 }}
                        className="absolute top-full mt-4 left-0 right-0 bg-red-950/20 border border-red-900/50 rounded-lg p-4 backdrop-blur-sm z-20"
                      >
                        <div className="flex items-start space-x-3 text-red-500">
                          <AlertTriangle className="w-4 h-4 shrink-0 mt-0.5" />
                          <div>
                            <h4 className="text-[10px] font-bold uppercase tracking-widest">Inversion Failed</h4>
                            <p className="text-[10px] font-mono mt-1 opacity-80 leading-relaxed uppercase">
                              Entropy was lost. Exact reconstruction is mathematically impossible.
                            </p>
                          </div>
                        </div>
                      </motion.div>
                    )}
                  </AnimatePresence>
                </div>

              </div>
              
              <div className="pt-6 border-t border-zinc-800 flex items-center justify-between">
                 <div className="space-y-1">
                  <h3 className="text-[10px] text-zinc-500 uppercase tracking-widest font-bold">Information Theory</h3>
                  <p className="text-[9px] text-zinc-600 uppercase tracking-widest font-mono">Operations that reduce state space cannot be reversed.</p>
                </div>
                
                <button
                  onClick={() => setAttemptInvert(true)}
                  className="px-4 py-2 border border-amber-500/50 text-amber-500 text-[10px] uppercase font-bold rounded flex items-center space-x-2 transition-colors hover:bg-amber-500/10"
                >
                  <RefreshCw className="w-3 h-3" />
                  <span>Attempt Inversion</span>
                </button>
              </div>

            </div>
          </div>

        </div>
      </div>
    </div>
  );
}
