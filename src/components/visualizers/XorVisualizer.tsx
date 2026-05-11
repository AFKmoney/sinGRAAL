import React, { useState } from 'react';
import { RefreshCw, ArrowRight, ArrowLeft } from 'lucide-react';
import { motion } from 'motion/react';

export function XorVisualizer() {
  const [inputA, setInputA] = useState('1101');
  const [inputB, setInputB] = useState('1011'); // This acts as the Key
  const [isReversing, setIsReversing] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleInput = (val: string, setter: (v: string) => void) => {
    if (/[^01]/.test(val)) {
      setError('Only binary values (0, 1) are permitted.');
    } else {
      setError(null);
    }
    setter(val.replace(/[^01]/g, ''));
  };

  // Pad or trim strings to 8 bits for visual consistency
  const formatBits = (str: string) => {
    const clean = str.replace(/[^01]/g, '').slice(0, 8);
    return clean;
  };

  const aBits = formatBits(inputA).padEnd(4, '0');
  const bBits = formatBits(inputB).padEnd(4, '0');

  // Compute XOR output
  const outBits = Array.from({ length: Math.max(aBits.length, bBits.length) }, (_, i) => {
    const a = aBits[i] === '1' ? 1 : 0;
    const b = bBits[i] === '1' ? 1 : 0;
    return String(a ^ b);
  }).join('');

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-4xl mx-auto space-y-12">
          
          <header className="space-y-4">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase flex items-baseline">
              XOR_Transform <span className="text-[10px] text-zinc-500 font-mono align-top ml-3 tracking-widest font-bold">REVERSIBLE COUPLING</span>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-widest max-w-2xl leading-relaxed">
              The Exclusive-OR (XOR) gate is a fundamental building block of cryptography. 
              Because it preserves entropy perfectly, if you know the Key, you can always reverse the Output to get the original Input.
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
                      className="w-full bg-transparent border-none font-mono text-sm text-emerald-400 tracking-[0.5em] focus:outline-none focus:ring-0 p-0"
                      maxLength={8}
                    />
                  </div>
                  <div className="bg-zinc-900/50 border border-zinc-800 rounded-lg p-4">
                    <label className="text-[10px] text-zinc-500 uppercase tracking-widest block mb-2 font-bold">XOR_Key (K)</label>
                    <input 
                      type="text" 
                      value={inputB}
                      onChange={(e) => handleInput(e.target.value, setInputB)}
                      className="w-full bg-transparent border-none font-mono text-sm text-cyan-400 tracking-[0.5em] focus:outline-none focus:ring-0 p-0"
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
                  <motion.div 
                    animate={{ rotate: isReversing ? 180 : 0 }}
                    className="w-12 h-12 bg-zinc-900/50 border border-zinc-800 rounded-lg flex items-center justify-center text-zinc-400 font-mono shadow-[0_0_10px_rgba(255,255,255,0.05)]"
                  >
                    XOR
                  </motion.div>
                  <p className="text-[9px] font-bold uppercase tracking-widest text-zinc-600">
                    {isReversing ? 'INVERTING...' : 'MATRIX TRANSLATION'}
                  </p>
                </div>

                {/* Output block */}
                <div className="bg-zinc-900/50 border border-zinc-800 rounded-lg p-4 flex flex-col justify-center h-full min-h-[104px]">
                  <label className="text-[10px] text-zinc-500 uppercase tracking-widest block mb-2 font-bold">
                    {isReversing ? 'Recovered Origin' : 'Transformed Element'}
                  </label>
                  <div className="w-full bg-transparent font-mono text-sm text-indigo-400 tracking-[0.5em]">
                    {outBits}
                  </div>
                </div>

              </div>
              
              {/* Magic Reverser */}
              <div className="pt-6 border-t border-zinc-800 flex items-center justify-between">
                <div className="space-y-1">
                  <h3 className="text-[10px] text-zinc-500 uppercase tracking-widest font-bold">Mathematical Proof</h3>
                  <p className="text-[9px] text-zinc-600 uppercase tracking-widest font-mono">A ^ K = C | C ^ K = A</p>
                </div>
                
                <button
                  onClick={() => setIsReversing(!isReversing)}
                  className={`
                    px-4 py-2 text-[10px] uppercase font-bold tracking-widest rounded flex items-center space-x-2 transition-colors
                    ${isReversing 
                      ? 'border border-zinc-700 text-zinc-500 bg-transparent hover:bg-zinc-900/50' 
                      : 'bg-zinc-100 text-black hover:bg-zinc-200'}
                  `}
                >
                  <RefreshCw className={`w-3 h-3 ${isReversing ? 'animate-spin' : ''}`} />
                  <span>{isReversing ? 'Reversing Matrix' : 'Apply Inverter'}</span>
                </button>
              </div>

            </div>
          </div>

        </div>
      </div>
    </div>
  );
}
