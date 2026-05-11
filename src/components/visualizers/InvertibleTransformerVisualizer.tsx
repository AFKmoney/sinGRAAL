import React, { useState } from 'react';
import { RotateCcw, FunctionSquare, ShieldAlert, Cpu, Settings2 } from 'lucide-react';
import { motion, AnimatePresence } from 'motion/react';

export function InvertibleTransformerVisualizer() {
  const [inputText, setInputText] = useState('Data');
  const [numLayers, setNumLayers] = useState<number>(3);
  const [keyMode, setKeyMode] = useState<'fixed' | 'custom' | 'auto'>('auto');
  const [fixedKey, setFixedKey] = useState('Bw9s');
  const [customKeys, setCustomKeys] = useState('Key1,Key2,Key3');
  const [autoKeyLen, setAutoKeyLen] = useState<number>(4);

  const [layers, setLayers] = useState<Array<{ id: number, outputText: string, keyUsed: string }>>([]);
  const [reversed, setReversed] = useState<string | null>(null);
  
  const generateKey = (len: number) => {
    let res = '';
    while (res.length < len) res += Math.random().toString(36).substring(2);
    return res.substring(0, len);
  };

  const processForwardNetwork = () => {
    let currentInput = inputText;
    let newLayers = [];
    const splits = customKeys.split(',').map(s => s.trim()).filter(Boolean);
    
    for (let i = 0; i < numLayers; i++) {
        let keyToUse = '';
        if (keyMode === 'fixed') keyToUse = fixedKey || 'A';
        else if (keyMode === 'custom') keyToUse = splits[i % splits.length] || 'A';
        else keyToUse = generateKey(autoKeyLen);
        
        let output = '';
        for(let j=0; j < currentInput.length; j++) {
            let charCode = currentInput.charCodeAt(j);
            let keyChar = keyToUse.charCodeAt(j % keyToUse.length);
            output += String.fromCharCode(charCode ^ keyChar);
        }
        currentInput = output;
        newLayers.push({ id: i + 1, outputText: output, keyUsed: keyToUse });
    }
    
    setLayers(newLayers);
    setReversed(null);
  };

  const applyReversePass = () => {
    let currentData = layers[layers.length - 1].outputText;
    
    // Go backwards
    for (let i = layers.length - 1; i >= 0; i--) {
        let key = layers[i].keyUsed;
        let output = '';
        for(let j=0; j < currentData.length; j++) {
            let charCode = currentData.charCodeAt(j);
            let keyChar = key.charCodeAt(j % key.length);
            output += String.fromCharCode(charCode ^ keyChar);
        }
        currentData = output;
    }
    
    setReversed(currentData);
  };
  
  const reset = () => {
    setLayers([]);
    setReversed(null);
  };

  return (
    <div className="flex flex-col h-full">
      <div className="flex-1 p-8 overflow-y-auto">
        <div className="max-w-4xl mx-auto space-y-12">
          
          <header className="space-y-4">
            <h2 className="text-2xl font-light tracking-widest text-white uppercase flex items-baseline">
              Invertible Transformer <span className="text-[10px] text-zinc-500 font-mono align-top ml-3 tracking-widest font-bold">XOR-BASED MATRICES</span>
            </h2>
            <p className="text-[10px] text-zinc-500 uppercase tracking-widest max-w-2xl leading-relaxed">
              Conceptual framework for an Invertible Transformer. By replacing destructive non-linearities (like ReLU) with inherently reversible operations (like XOR combinations), we guarantee that the entire forward pass can be perfectly reversed to recover the original input. This preserves 100% entropy across infinite layers.
            </p>
          </header>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl overflow-hidden flex flex-col p-6 relative">
               <div className="absolute inset-0 opacity-20 pointer-events-none" style={{ backgroundImage: 'radial-gradient(#333 1px, transparent 1px)', backgroundSize: '10px 10px' }}></div>
               
               <div className="relative z-10 space-y-6">
                 <div className="space-y-2">
                   <label className="text-[10px] text-zinc-500 uppercase tracking-widest font-bold">Input Sequence</label>
                   <input 
                     type="text" 
                     value={inputText}
                     onChange={(e) => { setInputText(e.target.value); reset(); }}
                     className="w-full bg-zinc-900/50 border border-zinc-800 rounded p-3 text-emerald-400 font-mono text-sm focus:outline-none focus:border-zinc-600"
                   />
                 </div>
                 
                 <div className="bg-zinc-900/50 p-4 rounded-lg border border-zinc-800 space-y-4">
                    <div className="flex items-center space-x-2 text-[10px] text-zinc-500 uppercase tracking-widest font-bold border-b border-zinc-800 pb-2">
                       <Settings2 className="w-3 h-3" />
                       <span>Network Configuration</span>
                    </div>
                    
                    <div className="grid grid-cols-2 gap-4">
                        <div className="space-y-2">
                          <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Layers (Num)</label>
                          <input 
                            type="number" min="1" max="20"
                            value={numLayers}
                            onChange={(e) => { setNumLayers(Math.max(1, parseInt(e.target.value) || 1)); reset(); }}
                            className="w-full bg-zinc-950 border border-zinc-800 rounded p-2 text-zinc-300 font-mono text-sm focus:outline-none"
                          />
                        </div>
                        <div className="space-y-2">
                          <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Key Gen Mode</label>
                          <select 
                            value={keyMode} 
                            onChange={(e) => { setKeyMode(e.target.value as 'fixed' | 'custom' | 'auto'); reset(); }}
                            className="w-full bg-zinc-950 border border-zinc-800 rounded p-2 text-zinc-300 font-mono text-xs focus:outline-none h-[38px]"
                          >
                            <option value="auto">Auto (Random)</option>
                            <option value="fixed">Fixed Key</option>
                            <option value="custom">Custom Seq.</option>
                          </select>
                        </div>
                    </div>

                    {keyMode === 'auto' && (
                        <div className="space-y-2">
                          <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Auto Key Length</label>
                          <input 
                            type="number" min="1" max="32"
                            value={autoKeyLen}
                            onChange={(e) => { setAutoKeyLen(Math.max(1, parseInt(e.target.value) || 4)); reset(); }}
                            className="w-full bg-zinc-950 border border-zinc-800 rounded p-2 text-zinc-300 font-mono text-sm focus:outline-none"
                          />
                        </div>
                    )}

                    {keyMode === 'fixed' && (
                        <div className="space-y-2">
                          <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Fixed XOR Key</label>
                          <input 
                            type="text" 
                            value={fixedKey}
                            onChange={(e) => { setFixedKey(e.target.value); reset(); }}
                            className="w-full bg-zinc-950 border border-zinc-800 rounded p-2 text-cyan-400 font-mono text-sm focus:outline-none"
                          />
                        </div>
                    )}

                    {keyMode === 'custom' && (
                        <div className="space-y-2">
                          <label className="text-[9px] text-zinc-500 uppercase tracking-widest">Keys (Comma-separated)</label>
                          <input 
                            type="text" 
                            value={customKeys}
                            onChange={(e) => { setCustomKeys(e.target.value); reset(); }}
                            className="w-full bg-zinc-950 border border-zinc-800 rounded p-2 text-cyan-400 font-mono text-sm focus:outline-none"
                            placeholder="e.g. A, B1, C2"
                          />
                        </div>
                    )}
                 </div>

                 <div className="flex gap-4">
                     <button
                       onClick={processForwardNetwork}
                       className="px-4 py-2 border border-zinc-700 bg-zinc-800 text-white hover:bg-zinc-700 text-[10px] uppercase font-bold tracking-widest rounded flex justify-center items-center space-x-2 transition-colors flex-1"
                     >
                       <FunctionSquare className="w-3 h-3" />
                       <span>Process Full Network</span>
                     </button>
                     <button
                       onClick={reset}
                       className="px-4 py-2 bg-transparent text-zinc-500 hover:text-zinc-300 text-[10px] uppercase font-bold tracking-widest rounded transition-colors border border-zinc-800"
                     >
                       Clear
                     </button>
                 </div>

                 <div className="space-y-2 mt-4 max-h-[300px] overflow-y-auto pr-2 custom-scrollbar">
                   {layers.length > 0 && <label className="text-[10px] text-zinc-500 uppercase tracking-widest font-bold block mb-3">Transformation Layers: {layers.length}</label>}
                   <AnimatePresence>
                     {layers.map((layer, idx) => (
                        <motion.div 
                          key={layer.id}
                          initial={{ opacity: 0, x: -10 }}
                          animate={{ opacity: 1, x: 0 }}
                          className="bg-zinc-950 border border-zinc-800 p-3 rounded flex flex-col space-y-2 mb-2"
                        >
                          <div className="flex justify-between items-center">
                            <span className="text-[9px] text-zinc-500 uppercase tracking-widest">Layer {layer.id}</span>
                            <span className="text-[9px] text-cyan-600 uppercase tracking-widest">Key: {layer.keyUsed}</span>
                          </div>
                          <div className="text-xs font-mono text-indigo-400 break-all">{btoa(unescape(encodeURIComponent(layer.outputText)))}</div>
                        </motion.div>
                     ))}
                   </AnimatePresence>
                 </div>

                 {layers.length > 0 && (
                    <button
                      onClick={applyReversePass}
                      className="w-full px-4 py-3 bg-zinc-100 text-black hover:bg-white text-[10px] uppercase font-bold tracking-widest rounded flex justify-center items-center space-x-2 transition-colors mt-4"
                    >
                      <RotateCcw className="w-3 h-3" />
                      <span>Execute Reverse Pass</span>
                    </button>
                 )}

                 <AnimatePresence>
                    {reversed && (
                      <motion.div 
                        initial={{ opacity: 0, y: 10 }}
                        animate={{ opacity: 1, y: 0 }}
                        className="bg-emerald-500/10 border border-emerald-500/20 p-4 rounded mt-4"
                      >
                        <label className="text-[10px] text-emerald-500 uppercase tracking-widest font-bold block mb-2">Recovered Sequence</label>
                        <div className="text-sm font-mono text-emerald-400 break-all">{reversed}</div>
                      </motion.div>
                    )}
                 </AnimatePresence>
               </div>
            </div>

            <div className="bg-zinc-900/20 border border-zinc-800 rounded-xl flex flex-col p-6 space-y-6">
                <div>
                   <h3 className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold mb-3 border-b border-zinc-800 pb-2">Conceptual Framework</h3>
                   <div className="space-y-3 text-[10px] text-zinc-500 uppercase tracking-widest leading-relaxed">
                     <p><strong className="text-zinc-300">1. Linear Mixing via XOR:</strong> Instead of traditional matrix multiplication, elements interact via bitwise XOR. This prevents dimensional loss.</p>
                     <p><strong className="text-zinc-300">2. Key Scheduling:</strong> Each layer receives a deterministic but pseudo-random key derived from a master sequence, ensuring transformation complexity without destroying state.</p>
                     <p><strong className="text-zinc-300">3. Non-Destructive Activations:</strong> If non-linearity is needed, it must be bijective (1-to-1 map) so its inverse function can perfectly restore the pre-activation state.</p>
                   </div>
                </div>

                <div className="flex-1">
                   <h3 className="text-[10px] text-zinc-400 uppercase tracking-widest font-bold mb-3 border-b border-zinc-800 pb-2">Algorithm Pseudocode</h3>
                   <pre className="bg-zinc-950 border border-zinc-800 p-4 rounded text-[10px] font-mono text-zinc-400 overflow-x-auto">
{`function forward_pass(input_X, keys):
  state = input_X
  for i in range(num_layers):
    K = keys[i]
    state = state ^ K  // Bitwise XOR
    # Optional: Apply Bijective S-Box here
  return state

function backward_pass(output_Y, keys):
  state = output_Y
  for i in reverse(range(num_layers)):
    K = keys[i]
    # Optional: Apply Inverse Bijective S-Box
    state = state ^ K  // XOR is its own inverse!
  return state`}
                   </pre>
                </div>
            </div>
          </div>

        </div>
      </div>
    </div>
  );
}
