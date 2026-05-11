import React, { useState } from 'react';
import { Network, Zap, CornerDownRight } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';

export function TransformerVisualizer() {
  const [inputVal, setInputVal] = useState<number>(-5);
  const [attemptInvert, setAttemptInvert] = useState(false);

  // Simple Feed-Forward + ReLU
  const weight = 2;
  const bias = 3;
  const preRelu = inputVal * weight + bias;
  const output = Math.max(0, preRelu);

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-4xl mx-auto space-y-12">
          
          <header className="space-y-4">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase flex items-baseline">
              Transformer Matrix <span className="text-[10px] text-zinc-500 font-mono align-top ml-3 tracking-widest font-bold">NON-LINEAR LOSS</span>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-widest max-w-2xl leading-relaxed">
              Modern AI models like Transformers rely heavily on non-linear activation functions (like ReLU) to learn complex patterns. These functions intentionally destroy information to create abstraction, making perfect mathematical inversion backward through the network impossible.
            </p>
          </header>

          <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl overflow-hidden flex flex-col p-8 relative">
            <div className="absolute inset-0 opacity-20 pointer-events-none" style={{ backgroundImage: 'radial-gradient(#333 1px, transparent 1px)', backgroundSize: '20px 20px' }}></div>
            
            <div className="relative z-10 space-y-12">
              
              <div className="flex flex-col md:flex-row items-center justify-between gap-8 py-8">
                
                {/* Input node */}
                <div className="flex flex-col items-center">
                  <div className="w-20 h-20 rounded bg-zinc-900/50 border border-zinc-800 flex items-center justify-center shadow-lg relative group">
                    <input 
                      type="number"
                      value={inputVal}
                      onChange={(e) => {
                        setInputVal(Number(e.target.value) || 0);
                        setAttemptInvert(false);
                      }}
                      className="w-16 bg-transparent text-center font-mono text-xl text-blue-400 focus:outline-none"
                    />
                    <div className="absolute -top-6 text-[9px] font-mono font-bold tracking-widest text-zinc-500 uppercase">Input (x)</div>
                  </div>
                </div>

                {/* Linear transform */}
                <div className="flex flex-col items-center flex-1 w-full relative">
                  <div className="h-px w-full bg-zinc-800 absolute top-1/2 -translate-y-1/2 -z-10"></div>
                  <div className="flex flex-col items-center space-y-2 bg-zinc-900/50 border border-zinc-800 p-3 rounded shadow-xl">
                    <div className="text-[9px] uppercase tracking-widest font-mono text-zinc-500">Linear: x * {weight} + {bias}</div>
                    <div className="font-mono text-xs text-zinc-300">{preRelu}</div>
                  </div>
                </div>

                {/* ReLU Node */}
                <div className="flex flex-col items-center flex-1 w-full relative">
                  <div className="h-px w-full bg-zinc-800 absolute top-1/2 -translate-y-1/2 -z-10"></div>
                  <div className="w-16 h-16 rounded bg-zinc-900/50 border border-zinc-800 flex items-center justify-center shadow-xl rotate-45 group hover:border-blue-500/50 transition-colors">
                    <div className="-rotate-45 flex flex-col items-center">
                      <Zap className="w-4 h-4 text-yellow-500 mb-1" />
                      <span className="text-[9px] font-mono font-bold tracking-widest text-zinc-500">ReLU</span>
                    </div>
                  </div>
                </div>

                {/* Arrow to output */}
                <div className="h-px w-16 bg-zinc-800 -z-10 hidden md:block"></div>

                {/* Output Node */}
                <div className="flex flex-col items-center">
                  <div className="w-20 h-20 rounded bg-blue-900/20 border border-blue-500/50 flex items-center justify-center shadow-[0_0_15px_rgba(59,130,246,0.1)] relative">
                    <span className="font-mono text-xl text-blue-300">{output}</span>
                    <div className="absolute -top-6 text-[9px] font-mono font-bold tracking-widest text-blue-400 uppercase">Activation</div>
                  </div>
                </div>

              </div>
              
              <div className="pt-6 border-t border-zinc-800 flex flex-col md:flex-row items-center justify-between gap-6">
                 <div className="space-y-1 md:max-w-md">
                  <h3 className="text-[10px] text-zinc-500 uppercase tracking-widest font-bold">The "Dead Neuron" Problem</h3>
                  <p className="text-[9px] text-zinc-600 uppercase tracking-widest font-mono">If pre-activation is negative, ReLU outputs 0.</p>
                </div>
                
                <button
                  onClick={() => setAttemptInvert(true)}
                  className="px-4 py-2 border border-zinc-700 text-zinc-300 hover:bg-zinc-900/50 text-[10px] uppercase font-bold tracking-widest rounded flex justify-center items-center space-x-2 transition-colors w-full md:w-auto"
                >
                  <Network className="w-3 h-3" />
                  <span>Reverse Propagate</span>
                </button>
              </div>

              <AnimatePresence>
                {attemptInvert && (
                  <motion.div 
                    initial={{ opacity: 0, y: 10 }}
                    animate={{ opacity: 1, y: 0 }}
                    className="bg-zinc-900/50 p-6 rounded-lg border border-zinc-800 shadow-inner flex flex-col"
                  >
                    <div className="space-y-4">
                      
                      <div className="flex items-center space-x-2 font-mono text-[10px] text-zinc-500 uppercase tracking-widest">
                        <span>Output: {output}</span>
                        <CornerDownRight className="w-3 h-3 text-zinc-600" />
                        <span>Inverse Linear: (y - {bias}) / {weight}</span>
                      </div>

                      {output === 0 ? (
                        <div className="p-4 bg-red-950/20 border border-red-900/50 rounded-lg text-red-500 text-[10px] uppercase tracking-widest leading-relaxed">
                          <strong className="block mb-2 font-bold">Inversion Failed: Information Nullified</strong>
                          Because the output is 0, the pre-activation value could have been <em>any number</em> less than or equal to 0. 
                          The exact state ({inputVal}) is unrecoverable.
                        </div>
                      ) : (
                        <div className="p-4 bg-emerald-950/20 border border-emerald-900/50 rounded-lg text-emerald-500 text-[10px] uppercase tracking-widest leading-relaxed">
                          <strong className="block mb-2 font-bold">Inversion Successful (Conditional)</strong>
                          Because the output was strictly greater than 0, it bypassed the destructive nature of ReLU. 
                          <div className="mt-3 text-emerald-400 font-mono text-xs bg-emerald-950/50 p-3 rounded border border-emerald-900/30 inline-block">
                            Original Input = ({output} - {bias}) / {weight} = {inputVal}
                          </div>
                        </div>
                      )}

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
